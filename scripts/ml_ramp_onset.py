"""Maximum-likelihood ramp-model EMG onset fit. Research code; not used by the detectors in src/.

Model (adapted; see docs/ml-ramp-onset-results.md for what is and is not verified in the sources):
  1. The band-passed EMG is whitened with an AR(``AR_ORDER``) filter estimated on the calibration
     baseline (Yule-Walker), and scaled to unit baseline SD. The whitening filter plus likelihood
     approach follows Staude 2001 (IEEE TBME 48(11):1292-1305, doi:10.1109/10.959325), which is
     about kinematic signals.
  2. Whitened samples are N(0, 1) before the ramp start t0, then N(0, a(t)^2) with
     a(t) = 1 + (G - 1) * clip((t - t0) / R, 0, 1).
  3. (t0, R, G) is the arg-max of the Gaussian likelihood: coarse grid (FFT correlation over every
     t0 for each (R, G)), then a local continuous search. It is the arg-max, not the posterior mean
     that src/emg_onset_lead/detect.py::RampOnset reports.
  4. The onset reported is where the FITTED ramp crosses ``ONSET_FRACTION`` of the plateau excursion,
     t0 + ONSET_FRACTION * R, minus the energy delay of band-pass + whitener.

Energy is pooled over ``BLOCK_S`` blocks; (t0, R, G) are continuous in the refinement.
"""

from __future__ import annotations

import numpy as np
from scipy import fft, linalg, signal

ONSET_FRACTION = 0.10
AR_ORDER = 10
BLOCK_S = 0.002
RISE_GRID_S = np.array([0, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 130, 160, 200]) / 1e3
RISE_MAX_S = 0.25
GAIN_GRID = np.geomspace(1.3, 14.0, 12)
GAIN_MIN, GAIN_MAX = 1.05, 30.0
N_REFINE = 6


def ar_whitener(base: np.ndarray, order: int = AR_ORDER) -> tuple[np.ndarray, float]:
    """(FIR prediction-error filter [1, -a1, ...], residual SD) from a baseline segment."""
    x = np.asarray(base, float) - np.mean(base)
    n = x.size
    r = np.array([np.dot(x[: n - k], x[k:]) for k in range(order + 1)]) / n
    a = linalg.solve_toeplitz(r[:-1], r[1:])
    A = np.concatenate(([1.0], -a))
    res = signal.lfilter(A, [1.0], x)[order:]
    return A, float(np.std(res, ddof=1))


def filter_delay_s(bp_impulse: np.ndarray, A: np.ndarray, fs: float) -> float:
    """Energy centroid (s) of the band-pass impulse response followed by the whitener."""
    e = np.convolve(bp_impulse, A) ** 2
    return float(np.dot(np.arange(e.size), e) / e.sum() / fs)


class MlRamp:
    """Vectorised over rows: ``fit(z)`` takes (rows, n) whitened unit-baseline samples."""

    def __init__(self, fs: float, n: int, block_s: float = BLOCK_S):
        self.fs, self.m = fs, max(1, round(block_s * fs))
        self.B = n // self.m
        self.h = self.m / fs  # block length, s
        B, m = self.B, self.m
        # time of a block's mean sample: samples sit at i/fs, so this is half a sample before (b+0.5)*h
        self.c = (np.arange(B) * m + 0.5 * (m - 1)) / fs
        u = self.c
        rows = []
        self.tmpl = [(R, G) for R in RISE_GRID_S for G in GAIN_GRID]
        for R, G in self.tmpl:
            s = np.ones(B) if R <= 0 else np.clip(u / R, 0, 1)
            rows.append(1 + (G - 1) * s)
        a = np.array(rows)
        self.cum_log = np.cumsum(m * np.log(a), axis=1)  # (T, B): first B-j terms
        q = 0.5 * (1 / a**2 - 1)
        self.nfft = fft.next_fast_len(2 * B - 1, real=True)
        self.qf = fft.rfft(q, self.nfft, axis=1)
        self.tR = np.array([t[0] for t in self.tmpl])
        self.tG = np.array([t[1] for t in self.tmpl])

    def block_energy(self, z: np.ndarray) -> np.ndarray:
        z = np.atleast_2d(z)[:, : self.B * self.m]
        return (z**2).reshape(z.shape[0], self.B, self.m).sum(axis=2)

    def coarse_loglik(self, e: np.ndarray) -> np.ndarray:
        """(rows, T, B) log-likelihood relative to 'no onset', for a ramp start at block j."""
        ee = fft.rfft(e[:, ::-1], self.nfft, axis=1)
        corr = fft.irfft(self.qf[None] * ee[:, None, :], self.nfft, axis=2)[:, :, : self.B][:, :, ::-1]
        j = np.arange(self.B)
        ll = -(self.cum_log[:, self.B - 1 - j][None] + corr)
        ll[:, :, :2] = -np.inf
        ll[:, :, -3:] = -np.inf
        return ll

    def loglik(self, e, t0, R, G) -> np.ndarray:
        """Exact block log-likelihood; e (rows, B); t0/R/G (rows, K) -> (rows, K). t0 in s."""
        s = np.clip((self.c - t0[..., None]) / np.maximum(R[..., None], 1e-5), 0, 1)
        a = 1 + (G[..., None] - 1) * s
        return -(self.m * np.log(a) + e[:, None, :] * 0.5 * (1 / (a * a) - 1)).sum(axis=2)

    def fit(self, z: np.ndarray):
        """-> (t0_s from window start, R_s, G, loglik) arrays, one per row."""
        e = self.block_energy(z)
        n = e.shape[0]
        ll = self.coarse_loglik(e).reshape(n, -1)
        k = np.argmax(ll, axis=1)
        T, j = np.unravel_index(k, (len(self.tmpl), self.B))
        t0, R, G = j * self.h, self.tR[T].copy(), self.tG[T].copy()
        st, sr, sg = self.h, 0.005, 0.10
        d = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij")).reshape(3, -1)
        best = None
        for _ in range(N_REFINE):
            T0 = t0[:, None] + st * d[0]
            RR = np.clip(R[:, None] + (sr + 0.15 * R[:, None]) * (d[1] * (st / self.h)), 0, RISE_MAX_S)
            GG = np.clip(G[:, None] * np.exp(sg * d[2]), GAIN_MIN, GAIN_MAX)
            L = self.loglik(e, T0, RR, GG)
            i = np.argmax(L, axis=1)
            r = np.arange(n)
            t0, R, G, best = T0[r, i], RR[r, i], GG[r, i], L[r, i]
            st *= 0.7  # shrink all three search steps together
            sg *= 0.75
            sr *= 0.8
        return t0, R, G, best
