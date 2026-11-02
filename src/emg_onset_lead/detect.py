"""Causal streaming onset detectors.

EMG and angle/IMU share one pipeline: causal IIR filtering (``sosfilt`` with carried state),
a feature stream, then a baseline mean + k*SD threshold that must be held for ``hold_s``.
Samples are consumed chunk by chunk; results never depend on the chunking or on future samples.

The angle onset is not a straight-line extrapolation of velocity. A minimum-jerk elbow move
(Flash and Hogan 1985) leaves rest with zero acceleration, so velocity starts as t^2 and a
line fit over a short window fires late. At the trigger the IMU detector instead waits a fixed
horizon (still only past samples) and least-squares matches two rest-to-move families —
minimum-jerk and critically damped — on that same window.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import signal


@dataclass(frozen=True)
class OnsetParams:
    baseline_s: float = 1.0  # fixed calibration window at trial start; the limb must be at rest
    settle_s: float = 0.2  # skip filter start-up transient inside the baseline window
    k: float = 5.0  # threshold = mean + k * SD of baseline feature
    hold_s: float = 0.025  # feature must stay beyond threshold this long
    fit_s: float = 0.07  # live window for the line-fit fallback (EMG, and IMU if the kinematic fit fails)
    # Not at rest when the peak excursion exceeds this many early-baseline SDs.
    # 12 clears the EMG envelope at rest (observed up to ~10) and still flags a
    # 60 deg move that starts at 0.6 s or 0.9 s.
    rest_ratio: float = 12.0


@dataclass(frozen=True)
class EmgParams(OnsetParams):
    k: float = 6.0  # k=5 gave ~2/1000 baseline false alarms on synthetic noise
    hold_s: float = 0.03
    band_hz: tuple[float, float] = (20.0, 450.0)
    envelope_hz: float = 8.0
    order: int = 2


@dataclass(frozen=True)
class ImuParams(OnsetParams):
    smooth_hz: float = 10.0
    order: int = 2
    mode: str = "velocity"  # "velocity" (|d angle/dt|) or "angle" (|angle - baseline mean|)
    k: float = 5.0
    hold_s: float = 0.03
    # After the threshold hold, keep this much more of the move (past samples only) before
    # committing to an onset. A short line fit cannot see a slow curved rise; the likelihood can.
    refine_s: float = 0.55


@dataclass(frozen=True)
class OnsetResult:
    onset_s: float | None  # signal onset (back-dated); None if missed or the baseline was not rest
    detected_s: float | None  # time the detector could first commit (earliest a live system can know)
    baseline_rest_ok: bool = True


class BaselineRest:
    """Reject a trial whose calibration window already contains the movement."""

    @staticmethod
    def excursion_ratio(x, fs: float, settle_s: float, baseline_s: float) -> float:
        """Peak |x - early median| / SD of the early post-settle baseline.

        The early third supplies the noise scale so a move that lands inside the
        baseline cannot inflate the SD and hide itself.
        """
        x = np.asarray(x, dtype=float)
        n0, n1 = round(settle_s * fs), round(baseline_s * fs)
        seg = x[n0:n1]
        if seg.size < 6:
            return 0.0
        early = seg[: seg.size // 3]
        sd = float(early.std(ddof=1))
        if not np.isfinite(sd) or sd < 1e-12:
            sd = 1e-12
        return float(np.max(np.abs(seg - np.median(early))) / sd)

    @staticmethod
    def at_rest(x, fs: float, settle_s: float, baseline_s: float, max_ratio: float = 6.0) -> bool:
        return BaselineRest.excursion_ratio(x, fs, settle_s, baseline_s) <= max_ratio


class KinematicOnset:
    """Sub-sample onset of a rest-to-move angle, using only samples up to ``commit_idx``.

    Two one-scale families are fit on one fixed window (a late guess cannot drop the early
    samples from the residual): minimum-jerk of free duration, and a critically damped step
    of free time constant. The better residual wins. Nothing after the commit sample is read.
    """

    @staticmethod
    def estimate(angle, fs: float, commit_idx: int, baseline_s: float, settle_s: float,
                 lookback_s: float = 1.2) -> float | None:
        angle = np.asarray(angle, dtype=float)
        commit_idx = int(min(commit_idx, angle.size - 1))
        base_n = round(baseline_s * fs)
        if commit_idx <= base_n + 4:
            return None
        w0 = max(base_n, commit_idx - round(lookback_s * fs))
        rest_lo = round(settle_s * fs)
        rest = float(np.median(angle[rest_lo:base_n]))
        y = angle[w0:commit_idx + 1] - rest
        if y.size < 8:
            return None
        mj = KinematicOnset._fit_family(y, fs, KinematicOnset._min_jerk, np.linspace(0.20, 1.60, 15))
        mj = KinematicOnset._zoom(y, fs, KinematicOnset._min_jerk, mj, 0.08, 17)
        mj = KinematicOnset._zoom(y, fs, KinematicOnset._min_jerk, mj, 0.015, 13)
        cd = KinematicOnset._fit_family(y, fs, KinematicOnset._critically_damped, np.linspace(0.05, 0.45, 12))
        cd = KinematicOnset._zoom(y, fs, KinematicOnset._critically_damped, cd, 0.03, 13)
        cd = KinematicOnset._zoom(y, fs, KinematicOnset._critically_damped, cd, 0.008, 9)
        best = mj if mj[0] <= cd[0] else cd
        sse = best[0]
        yy = float(np.dot(y, y))
        if not np.isfinite(sse) or yy <= 0 or 1.0 - sse / yy < 0.5:
            return None
        return (w0 + best[1]) / fs

    @staticmethod
    def _min_jerk(n: int, duration: float, fs: float) -> np.ndarray:
        t = np.arange(n) / fs
        u = np.clip(t / duration, 0.0, 1.0)
        g = 10 * u**3 - 15 * u**4 + 6 * u**5
        g[t >= duration] = 1.0
        return g

    @staticmethod
    def _critically_damped(n: int, tau: float, fs: float) -> np.ndarray:
        t = np.arange(n) / fs
        return 1.0 - (1.0 + t / tau) * np.exp(-t / tau)

    @staticmethod
    def _shape_sse(y: np.ndarray, g: np.ndarray) -> np.ndarray:
        """Residual of y[i:] ≈ c * g[:n-i] at every onset index i (c free, either sign)."""
        n = y.size
        gy = np.correlate(y, g, mode="full")[n - 1:n - 1 + n]
        zz = np.cumsum(g * g)[np.clip(np.arange(n, 0, -1) - 1, 0, None)]
        yy = float(np.dot(y, y))
        with np.errstate(divide="ignore", invalid="ignore"):
            sse = yy - (gy * gy) / np.where(zz > 1e-12, zz, np.inf)
        sse[:2] = np.inf
        sse[-3:] = np.inf
        return sse

    @staticmethod
    def _fit_family(y, fs, builder, params) -> tuple[float, float, float]:
        best = (np.inf, 0.0, float(params[0]))
        for p in params:
            if p <= 0.03:
                continue
            sse = KinematicOnset._shape_sse(y, builder(y.size, float(p), fs))
            j = int(np.argmin(sse))
            if not np.isfinite(sse[j]) or sse[j] >= best[0]:
                continue
            best = (float(sse[j]), KinematicOnset._interp(sse, j), float(p))
        return best

    @staticmethod
    def _zoom(y, fs, builder, current, half: float, n: int):
        if not np.isfinite(current[0]):
            return current
        grid = np.linspace(current[2] - half, current[2] + half, n)
        return KinematicOnset._fit_family(y, fs, builder, grid)

    @staticmethod
    def _interp(sse: np.ndarray, j: int) -> float:
        if j <= 0 or j >= sse.size - 1:
            return float(j)
        s0, s1, s2 = float(sse[j - 1]), float(sse[j]), float(sse[j + 1])
        den = s0 - 2 * s1 + s2
        if den <= 1e-18:
            return float(j)
        return float(j + np.clip(0.5 * (s0 - s2) / den, -1.0, 1.0))


@dataclass
class _StreamDetector:
    fs: float
    params: OnsetParams
    two_sided: bool = False
    _n: int = field(default=0, init=False)
    _base: list = field(default_factory=list, init=False)
    _mean: float | None = field(default=None, init=False)
    _thr: float = field(default=0.0, init=False)
    _run: int = field(default=0, init=False)
    _tail: np.ndarray = field(default_factory=lambda: np.zeros(0), init=False)
    _hold_n: int = field(default=1, init=False)
    _refine_n: int = field(default=0, init=False)
    _pending: int | None = field(default=None, init=False)
    _stopped: bool = field(default=False, init=False)
    _trig_window: np.ndarray | None = field(default=None, init=False)
    _trig_idx: int = field(default=0, init=False)
    result: OnsetResult = field(default_factory=lambda: OnsetResult(None, None), init=False)

    def __post_init__(self):
        self._hold_n = max(1, round(self.params.hold_s * self.fs))
        self._refine_n = max(0, round(getattr(self.params, "refine_s", 0.0) * self.fs))

    def _feature(self, x: np.ndarray) -> np.ndarray:  # pragma: no cover - abstract
        raise NotImplementedError

    def _observed(self) -> np.ndarray | None:
        """Signal the rest check and the kinematic fit read. None skips the rest check."""
        return None

    @staticmethod
    def _sos_stream(sos: np.ndarray):
        """Return a stateful causal filter closure (initialised on the first sample)."""
        state = {"zi": None}

        def apply(x: np.ndarray) -> np.ndarray:
            if state["zi"] is None:
                state["zi"] = signal.sosfilt_zi(sos) * x[0]
            y, state["zi"] = signal.sosfilt(sos, x, zi=state["zi"])
            return y

        return apply

    def process(self, chunk) -> OnsetResult:
        x = np.asarray(chunk, dtype=float)
        if x.size == 0 or self._stopped or self.result.onset_s is not None:
            return self.result
        f = self._feature(x)
        start = self._n
        self._n += f.size
        p = self.params
        base_n = round(p.baseline_s * self.fs)
        settle_n = round(p.settle_s * self.fs)
        if self._mean is None:
            lo = max(settle_n - start, 0)
            hi = min(base_n - start, f.size)
            if hi > lo:
                self._base.append(f[lo:hi])
            if start + f.size < base_n:
                return self.result
            b = np.concatenate(self._base) if self._base else np.zeros(1)
            self._mean = float(b.mean())
            self._thr = p.k * max(float(b.std(ddof=1)) if b.size > 1 else 0.0, 1e-12)
            if not self._rest_ok():
                self._stopped = True
                self.result = OnsetResult(None, None, False)
                return self.result
            f, start = f[max(base_n - start, 0):], max(start, base_n)
            if f.size == 0 and self._pending is None:
                return self.result
        if self._pending is not None:
            return self._maybe_commit()
        if f.size == 0:
            return self.result
        dev = f - self._mean
        keep = max(3, round(self.params.fit_s * self.fs))
        signed = np.abs(dev) if self.two_sided else dev
        buf = np.concatenate([self._tail, signed])
        n_tail = self._tail.size
        self._tail = buf[-keep:]
        above = signed > self._thr
        idx = np.arange(f.size)
        last_false = np.maximum.accumulate(np.where(~above, idx, -1))
        run = np.where(last_false >= 0, idx - last_false, idx + 1 + self._run)
        run = np.where(above, run, 0)
        hits = np.flatnonzero(run >= self._hold_n)
        if hits.size:
            j = int(hits[0])
            window = buf[: n_tail + j + 1][-keep:]  # ends at the trigger sample, never later
            trig = start + j
            self._trig_window = window
            self._trig_idx = trig
            if self._refine_n <= 0:
                self._finish(self._backdate(window, trig), trig)
            else:
                self._pending = trig
                return self._maybe_commit()
        elif run.size:
            self._run = int(run[-1])
        return self.result

    def _rest_ok(self) -> bool:
        obs = self._observed()
        if obs is None:
            return True
        return BaselineRest.at_rest(
            obs, self.fs, self.params.settle_s, self.params.baseline_s, self.params.rest_ratio)

    def _maybe_commit(self) -> OnsetResult:
        commit = self._pending + self._refine_n
        if self._n - 1 < commit:
            return self.result
        self._finish(self._onset_samples(commit), commit)
        return self.result

    def _onset_samples(self, commit_idx: int) -> float:
        return self._backdate(self._trig_window, self._trig_idx)

    def _finish(self, onset_sample: float, detected_idx: int) -> None:
        onset_s = float(onset_sample) / self.fs
        # An onset inside the calibration window means the move started during baseline.
        if not np.isfinite(onset_s) or onset_s < self.params.baseline_s:
            self._stopped = True
            self.result = OnsetResult(None, None, False)
            return
        self.result = OnsetResult(onset_s, detected_idx / self.fs, True)

    def _backdate(self, window: np.ndarray, confirm_idx: int) -> float:
        """Line-fit onset in samples. Used for EMG, and as the IMU fallback.

        Least-squares line through the last ``fit_s`` of the baseline-standardised feature,
        extrapolated to the baseline level. Only samples at or before ``confirm_idx`` are in
        ``window``. A non-positive slope cannot be extrapolated; the threshold-crossing sample
        (the first sample of the hold) is reported instead.
        """
        w = window / (self._thr / self.params.k)
        if w.size < 3:
            return float(confirm_idx - self._hold_n + 1)
        slope, icpt = np.polyfit(np.arange(w.size), w, 1)
        if slope <= 0:
            return float(confirm_idx - self._hold_n + 1)
        return confirm_idx - (w.size - 1) - icpt / slope

    @classmethod
    def run(cls, x, fs: float, params=None, chunk: int | None = None) -> OnsetResult:
        """Process ``x`` in chunks of ``chunk`` samples (None = one shot)."""
        det = cls(fs, params) if params is not None else cls(fs)
        x = np.asarray(x, dtype=float)
        step = chunk or max(x.size, 1)
        for i in range(0, x.size, step):
            det.process(x[i:i + step])
        return det.result


class EmgOnsetDetector(_StreamDetector):
    """Band-pass -> rectify -> causal low-pass envelope -> baseline+k*SD threshold."""

    def __init__(self, fs: float, params: EmgParams = EmgParams()):
        super().__init__(fs, params, two_sided=False)
        lo, hi = params.band_hz
        hi = min(hi, 0.45 * fs)
        self._bp = self._sos_stream(
            signal.butter(params.order, [lo, hi], btype="bandpass", fs=fs, output="sos"))
        self._lp = self._sos_stream(
            signal.butter(params.order, params.envelope_hz, btype="lowpass", fs=fs, output="sos"))
        self._env: list[np.ndarray] = []

    def _feature(self, x):
        f = self._lp(np.abs(self._bp(x)))
        self._env.append(np.asarray(f, dtype=float))
        return f

    def _observed(self):
        return np.concatenate(self._env) if self._env else None


class ImuOnsetDetector(_StreamDetector):
    """Causal low-pass -> (causal first difference) -> baseline+k*SD, then a kinematic onset.

    The threshold only decides that a move is underway. The reported onset is the minimum-jerk
    or critically damped fit on the raw angle up to the commit sample (``refine_s`` after the
    hold). A failed fit falls back to the shared line extrapolation.
    """

    def __init__(self, fs: float, params: ImuParams = ImuParams()):
        if params.mode not in ("velocity", "angle"):
            raise ValueError(f"unknown mode {params.mode!r}")
        super().__init__(fs, params, two_sided=params.mode == "angle")
        self._lp = self._sos_stream(
            signal.butter(params.order, params.smooth_hz, btype="lowpass", fs=fs, output="sos"))
        self._prev: float | None = None
        self._ang: list[np.ndarray] = []

    def _feature(self, x):
        x = np.asarray(x, dtype=float)
        self._ang.append(x.copy())
        y = self._lp(x)
        if self.params.mode == "angle":
            return y
        prev = y[0] if self._prev is None else self._prev
        v = np.diff(np.concatenate(([prev], y))) * self.fs
        self._prev = float(y[-1])
        return np.abs(v)

    def _observed(self):
        return np.concatenate(self._ang) if self._ang else None

    def _onset_samples(self, commit_idx: int) -> float:
        ang = np.concatenate(self._ang)
        est = KinematicOnset.estimate(
            ang, self.fs, commit_idx, self.params.baseline_s, self.params.settle_s)
        if est is None:
            return self._backdate(self._trig_window, self._trig_idx)
        return est * self.fs
