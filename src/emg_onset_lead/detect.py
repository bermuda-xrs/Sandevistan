"""Causal streaming onset detectors.

EMG and angle/IMU share one pipeline: causal IIR filtering (``sosfilt`` with carried state),
a feature stream, then a baseline mean + k*SD threshold that must be held for ``hold_s``.
Samples are consumed chunk by chunk; results never depend on the chunking or on future samples.

The angle onset is not a straight-line extrapolation of velocity. A minimum-jerk elbow move
(Flash and Hogan 1985) leaves rest with zero acceleration, so velocity starts as t^2 and a
line fit over a short window fires late. At the trigger the IMU detector instead waits a fixed
horizon (still only past samples) and least-squares matches two rest-to-move families —
minimum-jerk and critically damped — on that same window.

The EMG onset is not the threshold crossing or a line through the 8 Hz envelope either: both
move with burst rise time and amplitude, because the envelope filter's delay and the k*SD
threshold sit at a fixed absolute level. After the hold the EMG detector waits ``refine_s``
(past samples only) and computes a change-point posterior on the band-passed signal itself
(``RampOnset``): Gaussian noise at the baseline SD before the onset, an amplitude ramp of free
rise time and free gain after it. The onset is the posterior mean over that bounded window.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import fft, signal


@dataclass(frozen=True)
class OnsetParams:
    baseline_s: float = 1.0  # fixed calibration window at trial start; the limb must be at rest
    settle_s: float = 0.2  # skip filter start-up transient inside the baseline window
    k: float = 5.0  # threshold = mean + k * SD of baseline feature
    hold_s: float = 0.03  # feature must stay beyond threshold this long (30 ms, outside the latency budget)
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
    # Onset estimate (``RampOnset``). The envelope threshold only decides that a burst is underway.
    refine_s: float = 0.1  # wait after the hold, then commit; the onset uses samples up to the commit
    lookback_s: float = 0.4  # window start before the hold sample (never before the baseline end)
    rise_max_s: float = 0.5  # rise-time grid 0, step, ..., max (uniform prior)
    rise_step_s: float = 0.02
    gain_range: tuple[float, float] = (1.2, 60.0)  # burst/baseline amplitude grid, log spaced
    n_gain: int = 20
    block_s: float = 0.002  # energy is pooled over blocks this long before the fit (speed)
    # A burst that has already ended by the commit is kept if its envelope peaked at least this many
    # thresholds (baseline-removed). Tradeoff (FC, synthetic white noise, 499 h): false onsets fall
    # from 0.016/h at 1.5x to 0.008/h at 1.6x (noise triggers peaked up to ~1.73x); the cost is
    # 40 ms / SNR 3 short-burst detection 246 -> 234 of 300. The lowest 40-80 ms burst peak at
    # SNR >= 5 was 2.696x (one stream: 60 trials, default_rng(1); a measured minimum, not a bound;
    # a property of the burst, not of this factor). A noise floor 10-20% above baseline dominates
    # over this factor (see README).
    short_burst_peak: float = 1.6


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
    # Time the threshold hold was met: a live controller knows a move is underway here. For a
    # detector with ``refine_s`` > 0 it is ``refine_s`` before ``detected_s``.
    trigger_s: float | None = None


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
    def at_rest(x, fs: float, settle_s: float, baseline_s: float,
                max_ratio: float = OnsetParams.rest_ratio) -> bool:
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


class RampOnset:
    """Onset of an EMG burst from the band-passed signal, using only a bounded past window.

    Model, in units of the baseline SD: before the onset ``j`` the samples are N(0, 1); from
    ``j`` on they are N(0, a(t - j)^2) with a(u) = 1 + (G - 1) * min(u / R, 1), a linear amplitude
    ramp of rise time R to a gain G (R = 0 is a step). For every (R, G) on a grid the Gaussian
    log-likelihood of each ``j`` is one FFT correlation. The onset is the posterior mean of ``j``
    with flat priors on ``j`` over the window, on R (uniform), and on log G.

    The posterior mean, not the arg-max, is used: when the rise is slow and the SNR low, onset
    and rise time trade off, and the arg-max snaps late (toward short R). A log-spaced R grid
    did the same (about +2.5 ms on 80-150 ms rises at SNR 3, synthetic); a uniform R prior
    out to 0.5 s removed it.
    """

    @staticmethod
    def estimate(z, fs: float, rise_grid_s, gain_grid, block: int = 1) -> float | None:
        """Posterior-mean onset, in samples from the start of ``z``. None if ``z`` is too short.

        ``block`` samples are pooled into one energy sum (the likelihood is exact for a block
        whose amplitude is constant; the ramp is evaluated at the block centre). Candidate
        onsets are block boundaries; the posterior mean is still continuous.
        """
        z = np.asarray(z, dtype=float)
        m = max(1, int(block))
        n = z.size // m
        if n < 16 or not np.all(np.isfinite(z)):
            return None
        e = (z[: n * m] ** 2).reshape(n, m).sum(axis=1)  # energy per block
        rise = np.asarray(rise_grid_s, dtype=float)
        gain = np.asarray(gain_grid, dtype=float)
        k = np.arange(n)
        u = (k + 0.5) * m / fs  # seconds from a candidate onset to each block centre
        h = np.ones((rise.size, n))
        slow = rise > 0
        h[slow] = np.clip(u[None, :] / rise[slow, None], 0.0, 1.0)
        a = (1.0 + (gain[None, :, None] - 1.0) * h[:, None, :]).reshape(-1, n)
        # NLL relative to "no onset": sum over blocks b >= j of m log a + e_b / (2 a^2) - e_b / 2.
        log_a = np.cumsum(m * np.log(a), axis=1)  # first n - j terms
        q = 0.5 * (1.0 / a**2 - 1.0)
        nfft = fft.next_fast_len(2 * n - 1, real=True)
        ee = fft.rfft(e[::-1], nfft)
        corr = fft.irfft(fft.rfft(q, nfft, axis=1) * ee[None, :], nfft, axis=1)[:, :n][:, ::-1]
        loglik = -(log_a[:, n - 1 - k] + corr)
        loglik[:, :2] = -np.inf  # need some baseline before and some burst after the onset
        loglik[:, -3:] = -np.inf
        top = float(np.max(loglik))
        if not np.isfinite(top):
            return None
        w = np.exp(loglik - top).sum(axis=0)
        return float(np.dot(w, k) / w.sum() * m)

    @staticmethod
    def energy_delay_s(sos: np.ndarray, fs: float, length_s: float = 0.5) -> float:
        """Energy centroid of the band-pass impulse response (seconds).

        The band-pass smears a change in noise power by its squared impulse response, which
        delays the fitted onset by about this much (about 1.15 ms for 20-450 Hz, order 2,
        2000 Hz). It depends on the filter only, so it is subtracted as a constant.
        """
        imp = np.zeros(max(2, round(length_s * fs)))
        imp[0] = 1.0
        e = signal.sosfilt(sos, imp) ** 2
        return float(np.dot(np.arange(e.size), e) / e.sum() / fs)


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
    _sig_hist: list = field(default_factory=list, init=False)
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
        if f.size:
            dev = f - self._mean
            signed = np.abs(dev) if self.two_sided else dev
            self._sig_hist.append(signed)
        if self._pending is not None:
            return self._maybe_commit()
        if f.size == 0:
            return self.result
        return self._scan(signed, start)

    def _scan(self, signed: np.ndarray, start: int) -> OnsetResult:
        """Threshold-and-hold scan of post-baseline samples ``signed`` starting at ``start``."""
        keep = max(3, round(self.params.fit_s * self.fs))
        buf = np.concatenate([self._tail, signed])
        n_tail = self._tail.size
        self._tail = buf[-keep:]
        above = signed > self._thr
        idx = np.arange(signed.size)
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
        if not self._confirm(commit):
            return self._rearm()
        self._finish(self._onset_samples(commit), commit)
        return self.result

    def _confirm(self, commit_idx: int) -> bool:
        """Whether the trigger still stands at the commit sample. The base detector keeps it."""
        return True

    def _signed_history(self) -> np.ndarray:
        """Baseline-removed feature from the baseline end onward (index 0 = baseline end)."""
        return np.concatenate(self._sig_hist) if self._sig_hist else np.zeros(0)

    def _rearm(self) -> OnsetResult:
        """Drop the pending trigger and rescan, from the sample after it, what is already in."""
        trig = self._pending
        self._pending = None
        self._run = 0
        sig = self._signed_history()
        rel = trig - round(self.params.baseline_s * self.fs)
        keep = max(3, round(self.params.fit_s * self.fs))
        self._tail = sig[max(0, rel + 1 - keep):rel + 1]
        rest = sig[rel + 1:]
        if rest.size == 0:
            return self.result
        return self._scan(rest, trig + 1)

    def _onset_samples(self, commit_idx: int) -> float:
        return self._backdate(self._trig_window, self._trig_idx)

    def _finish(self, onset_sample: float, detected_idx: int) -> None:
        onset_s = float(onset_sample) / self.fs
        # An onset inside the calibration window means the move started during baseline.
        if not np.isfinite(onset_s) or onset_s < self.params.baseline_s:
            self._stopped = True
            self.result = OnsetResult(None, None, False)
            return
        self.result = OnsetResult(onset_s, detected_idx / self.fs, True, self._trig_idx / self.fs)

    # A fitted rise this small, relative to the feature level and the threshold, is
    # numerical noise. ``polyfit`` on a flat window can return about +1e-17 instead of 0.
    _MIN_RISE_FRAC = 1e-6

    def _backdate(self, window: np.ndarray, confirm_idx: int) -> float:
        """Line-fit onset in samples. The fallback for both detectors.

        Least-squares line through the last ``fit_s`` of the baseline-standardised feature,
        extrapolated to the baseline level. Only samples at or before ``confirm_idx`` are in
        ``window``. A slope whose rise across that window is below ``_MIN_RISE_FRAC`` of the
        feature level and of the threshold — including a non-positive slope and a ~1 ulp fit
        on a flat window — is not extrapolated. The threshold-crossing sample (the first
        sample of the hold) is reported instead. The result stays inside the fit window, at
        or after the baseline end, and at or before the confirm sample.
        """
        fallback = float(confirm_idx - self._hold_n + 1)
        w = np.asarray(window, dtype=float)
        n = int(w.size)
        scale = self._thr / self.params.k if self.params.k else self._thr
        if n < 3 or not np.isfinite(scale) or scale <= 0:
            return self._bound_onset(fallback, confirm_idx, max(n, 1))
        w = w / scale
        slope, icpt = np.polyfit(np.arange(n, dtype=float), w, 1)
        level = max(float(np.max(np.abs(w))), float(self.params.k), 1e-12)
        rise = slope * (n - 1)
        if not np.isfinite(rise) or not np.isfinite(icpt) or rise <= self._MIN_RISE_FRAC * level:
            return self._bound_onset(fallback, confirm_idx, n)
        return self._bound_onset(float(confirm_idx - (n - 1) - icpt / slope), confirm_idx, n)

    def _bound_onset(self, onset: float, confirm_idx: int, n: int) -> float:
        """Clamp a back-dated sample so it cannot precede the fit window or the baseline.

        The baseline floor is applied when the baseline ends at or before ``confirm_idx``.
        When the confirm sample itself is still inside the baseline, the floor would sit
        after the confirm sample, so only the fit window and the confirm sample bind.
        """
        window_start = float(confirm_idx - (n - 1))
        lo = window_start
        baseline_end = float(round(self.params.baseline_s * self.fs))
        if baseline_end <= float(confirm_idx):
            lo = max(lo, baseline_end)
        hi = float(confirm_idx)
        if lo > hi:
            lo = hi
        if not np.isfinite(onset):
            return lo
        return float(min(max(float(onset), lo), hi))

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
    """Band-pass -> rectify -> causal low-pass envelope -> baseline+k*SD threshold and hold.

    The hold says a burst is underway. ``refine_s`` later the onset is ``RampOnset`` on the
    band-passed signal from ``lookback_s`` before the hold to the commit sample, scaled by the
    baseline SD, minus the band-pass energy delay. A failed fit falls back to the line fit.
    """

    def __init__(self, fs: float, params: EmgParams = EmgParams()):
        super().__init__(fs, params, two_sided=False)
        lo, hi = params.band_hz
        hi = min(hi, 0.45 * fs)  # 450 Hz at 1000 Hz; the 20-450 Hz band is unchanged from 1000 Hz up
        if not 0 < lo < hi:
            raise ValueError(f"EMG band {params.band_hz} Hz does not fit fs = {fs:g} Hz "
                             f"(upper edge capped at 0.45 * fs = {0.45 * fs:g} Hz)")
        bp_sos = signal.butter(params.order, [lo, hi], btype="bandpass", fs=fs, output="sos")
        self._bp = self._sos_stream(bp_sos)
        self._lp = self._sos_stream(
            signal.butter(params.order, params.envelope_hz, btype="lowpass", fs=fs, output="sos"))
        self._env: list[np.ndarray] = []
        self._bp_hist: list[np.ndarray] = []
        self._bp_delay_n = RampOnset.energy_delay_s(bp_sos, fs) * fs
        self._rise_grid = np.arange(0.0, params.rise_max_s + 1e-9, params.rise_step_s)
        self._gain_grid = np.geomspace(params.gain_range[0], params.gain_range[1], params.n_gain)
        self._block = max(1, round(params.block_s * fs))
        self._fit_end: int | None = None  # set by _confirm

    def _feature(self, x):
        b = self._bp(x)
        self._bp_hist.append(np.asarray(b, dtype=float))
        f = self._lp(np.abs(b))
        self._env.append(np.asarray(f, dtype=float))
        return f

    def _observed(self):
        return np.concatenate(self._env) if self._env else None

    def _confirm(self, commit_idx: int) -> bool:
        """Keep the trigger if a burst is still on at the commit, or was a clear short burst.

        - Envelope above threshold at the commit: a burst is on. The onset fit runs to the commit.
        - Envelope back below threshold, but its peak since the hold reached ``short_burst_peak``
          (1.6) thresholds: an isolated short burst (e.g. 40-80 ms) that has ended. The onset fit
          window ends at that envelope peak, so the post-burst baseline does not pull the
          ramp model.
        - Otherwise a noise excursion that met the 30 ms hold and decayed. The trigger is
          dropped and the detector re-arms (``_rearm``).
        """
        sig = self._signed_history()
        base_n = round(self.params.baseline_s * self.fs)
        rel = commit_idx - base_n
        self._fit_end = commit_idx
        if not 0 <= rel < sig.size:
            return False
        if sig[rel] > self._thr:
            return True
        seg = sig[self._pending - base_n:rel + 1]
        if seg.size == 0 or float(seg.max()) < self.params.short_burst_peak * self._thr:
            return False
        self._fit_end = self._pending + int(np.argmax(seg))
        return True

    def _onset_samples(self, commit_idx: int) -> float:
        p = self.params
        if self._fit_end is not None:
            commit_idx = min(commit_idx, self._fit_end)  # a short burst's fit ends at its peak
        bp = np.concatenate(self._bp_hist)
        base_n, settle_n = round(p.baseline_s * self.fs), round(p.settle_s * self.fs)
        sd = float(np.std(bp[settle_n:base_n], ddof=1)) if base_n - settle_n > 1 else 0.0
        w0 = max(base_n, self._trig_idx - round(p.lookback_s * self.fs))
        est = None
        if np.isfinite(sd) and sd > 0 and commit_idx < bp.size:
            est = RampOnset.estimate(bp[w0:commit_idx + 1] / sd, self.fs, self._rise_grid, self._gain_grid,
                                     self._block)
        if est is None:
            return self._backdate(self._trig_window, self._trig_idx)
        return float(max(w0 + est - self._bp_delay_n, base_n))


class ImuOnsetDetector(_StreamDetector):
    """Causal low-pass -> (causal first difference) -> baseline+k*SD, then a kinematic onset.

    The threshold only decides that a move is underway. The reported onset is the minimum-jerk
    or critically damped fit on the raw angle up to the commit sample (``refine_s`` after the
    hold). A failed fit falls back to the shared line extrapolation.
    """

    def __init__(self, fs: float, params: ImuParams = ImuParams()):
        if params.mode not in ("velocity", "angle"):
            raise ValueError(f"unknown mode {params.mode!r}")
        if not 0 < params.smooth_hz < 0.5 * fs:
            raise ValueError(f"IMU smooth_hz {params.smooth_hz:g} Hz must be below Nyquist ({0.5 * fs:g} Hz)")
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
