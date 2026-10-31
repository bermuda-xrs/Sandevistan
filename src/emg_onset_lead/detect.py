"""Causal streaming onset detectors.

EMG and angle/IMU share one pipeline: causal IIR filtering (``sosfilt`` with carried state),
a feature stream, then a baseline mean + k*SD threshold that must be held for ``hold_s``.
Samples are consumed chunk by chunk; results never depend on the chunking or on future samples.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import signal


@dataclass(frozen=True)
class OnsetParams:
    baseline_s: float = 1.0  # fixed calibration window at trial start
    settle_s: float = 0.2  # skip filter start-up transient inside the baseline window
    k: float = 5.0  # threshold = mean + k * SD of baseline feature
    hold_s: float = 0.025  # feature must stay beyond threshold this long
    fit_s: float = 0.07  # live window for the onset back-dating line fit (same for EMG and angle)


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


@dataclass(frozen=True)
class OnsetResult:
    onset_s: float | None  # first sample of the sustained run
    detected_s: float | None  # time the hold was satisfied (earliest a live system can know)


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
    result: OnsetResult = field(default_factory=lambda: OnsetResult(None, None), init=False)

    def __post_init__(self):
        self._hold_n = max(1, round(self.params.hold_s * self.fs))

    def _feature(self, x: np.ndarray) -> np.ndarray:  # pragma: no cover - abstract
        raise NotImplementedError

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
        if x.size == 0 or self.result.onset_s is not None:
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
            f, start = f[max(base_n - start, 0):], max(start, base_n)
            if f.size == 0:
                return self.result
        dev = f - self._mean
        keep = max(3, round(self.params.fit_s * self.fs))
        buf = np.concatenate([self._tail, np.abs(dev) if self.two_sided else dev])
        n_tail = self._tail.size
        self._tail = buf[-keep:]
        above = (np.abs(dev) if self.two_sided else dev) > self._thr
        idx = np.arange(f.size)
        last_false = np.maximum.accumulate(np.where(~above, idx, -1))
        run = np.where(last_false >= 0, idx - last_false, idx + 1 + self._run)
        run = np.where(above, run, 0)
        hits = np.flatnonzero(run >= self._hold_n)
        if hits.size:
            j = int(hits[0])
            window = buf[: n_tail + j + 1][-keep:]  # ends at the confirmation sample, never later
            self.result = OnsetResult(self._backdate(window, start + j) / self.fs, (start + j) / self.fs)
        else:
            self._run = int(run[-1])
        return self.result

    def _backdate(self, window: np.ndarray, confirm_idx: int) -> float:
        """Signal-onset estimate at confirmation (uses only samples <= confirm_idx).

        A threshold crossing lags the true onset by an amount that depends on rise speed and SNR,
        and that differs between EMG and angle. Both paths therefore apply the same rule over the
        same live window (``fit_s``): least-squares line through the last ``fit_s`` of the
        baseline-standardised feature, extrapolated back to the baseline level. The reported
        onset is that zero-crossing, not the threshold-crossing or the confirmation time.
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

    def _feature(self, x):
        return self._lp(np.abs(self._bp(x)))


class ImuOnsetDetector(_StreamDetector):
    """Causal low-pass -> (causal first difference) -> baseline+k*SD threshold."""

    def __init__(self, fs: float, params: ImuParams = ImuParams()):
        if params.mode not in ("velocity", "angle"):
            raise ValueError(f"unknown mode {params.mode!r}")
        super().__init__(fs, params, two_sided=params.mode == "angle")
        self._lp = self._sos_stream(
            signal.butter(params.order, params.smooth_hz, btype="lowpass", fs=fs, output="sos"))
        self._prev: float | None = None

    def _feature(self, x):
        y = self._lp(x)
        if self.params.mode == "angle":
            return y
        prev = y[0] if self._prev is None else self._prev
        v = np.diff(np.concatenate(([prev], y))) * self.fs
        self._prev = float(y[-1])
        return np.abs(v)
