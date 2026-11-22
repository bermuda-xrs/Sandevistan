"""Candidate causal EMG onset alarms for the live path. SYNTHETIC-ONLY RESEARCH CODE.

Nothing here is imported by ``src/``; the production detector (``EmgOnsetDetector``: 8 Hz envelope,
k = 6 SD, 30 ms hold, 0.1 s refine) is not changed. Every function works on a 2-D array
``x[segment, sample]`` of raw EMG at ``fs`` and returns an alarm statistic ``stat[segment, sample]``.
``stat[:, n]`` depends only on ``x[:, :n + 1]`` (causal; ``check_causal`` tests this by truncation).
An alarm is the first post-baseline sample where ``stat > h``. Thresholds ``h`` are not fixed here:
``run_sim.py`` calibrates each one to a false-alarm rate per minute of rest.

Shared front end (same as production): 2nd-order Butterworth band-pass 20-450 Hz, ``sosfilt`` with
the state initialised on the first sample, scaled by the SD of the 0.2-1.0 s baseline (the
production ``baseline_s`` = 1.0 s and ``settle_s`` = 0.2 s; not changed).

Detectors (``DETECTORS`` at the end lists the configured variants):
- ``env``:   |band-pass| -> causal low-pass (fc) -> (env - baseline mean) / baseline SD, then held
             ``hold_s`` (the statistic is the trailing minimum over the hold). fc = 8 Hz, hold 30 ms,
             h = 6 reproduces the production trigger (``trigger_s``) sample for sample.
- ``tkeo``:  causal Teager-Kaiser energy psi[n] = y[n-1]^2 - y[n-2] y[n] (one sample late, 0.5 ms),
             then |psi| -> low-pass -> baseline-standardised, held like ``env``.
- ``rms``:   RMS of the baseline-scaled band-pass over the last 25 or 50 ms (a shorter window).
- ``cusum``: Page's CUSUM on the per-sample log-likelihood ratio of a variance step from 1 to rho^2
             on the baseline-scaled band-pass z: l = 0.5 (1 - 1/rho^2) z^2 - ln rho.
             S_n = C_n - min(0, min_{k<=n} C_k), C = cumulative sum of l from the baseline end.
- ``glr``:   window-limited GLR for an unknown variance increase: max over the last L samples
             (L = 2 ms, 4 ms, ..., W) of 0.5 L (s - 1 - ln s), s = mean z^2 over them, s > 1.
             Only post-baseline samples enter (the baseline sets the scale, not the statistic).
- ``bocpd``: Adams-MacKay run-length recursion on 2 ms block energies, Gamma-conjugate on the
             segment variance, constant hazard, run lengths truncated at R blocks (a run that
             reaches R is merged back into the known-variance rest state). Statistic: log posterior
             odds of "changed 2..R blocks ago to a variance with posterior mean > ``v_min``" against
             "rest". (Odds against everything else saturate near ln(1/hazard): the just-started run
             always holds about the hazard's mass.)
"""

from __future__ import annotations

import numpy as np
from scipy import signal
from scipy.ndimage import minimum_filter1d
from scipy.special import gammaln

BASELINE_S = 1.0  # production baseline_s (protected)
SETTLE_S = 0.2  # production settle_s


def _sos_causal(sos, x):
    """Causal sosfilt along the last axis, state initialised on the first sample (as production)."""
    zi = signal.sosfilt_zi(sos)[:, None, :] * x[:, 0][None, :, None]
    y, _ = signal.sosfilt(sos, x, axis=-1, zi=zi)
    return y


def bandpass(x, fs, band=(20.0, 450.0), order=2):
    lo, hi = band
    sos = signal.butter(order, [lo, min(hi, 0.45 * fs)], btype="bandpass", fs=fs, output="sos")
    return _sos_causal(sos, np.atleast_2d(np.asarray(x, dtype=float)))


def _base_slice(fs):
    return slice(round(SETTLE_S * fs), round(BASELINE_S * fs))


def standardise(f, fs):
    """(f - baseline mean) / baseline SD (ddof=1), per segment, baseline = 0.2-1.0 s."""
    b = f[:, _base_slice(fs)]
    mu = b.mean(axis=1, keepdims=True)
    sd = np.maximum(b.std(axis=1, ddof=1, keepdims=True), 1e-12)
    return (f - mu) / sd


def scaled_bp(x, fs):
    """Band-pass scaled by the baseline SD (no mean removal: the band-pass has zero mean)."""
    y = bandpass(x, fs)
    sd = np.maximum(y[:, _base_slice(fs)].std(axis=1, ddof=1, keepdims=True), 1e-12)
    return y / sd


def trailing_min(d, n_hold, start=0):
    """min(d[n - n_hold + 1 .. n]) (causal). -inf until the window is full of samples >= start.

    ``minimum_filter1d`` windows are [i - size // 2 - origin, ...]; origin = (n_hold - 1) // 2 puts the
    last window sample at i (checked against a brute-force loop in ``test_detectors.py``).
    """
    out = d.copy() if n_hold <= 1 else minimum_filter1d(d, size=n_hold, axis=-1, origin=(n_hold - 1) // 2,
                                                        mode="nearest")
    out[:, : start + n_hold - 1] = -np.inf
    return out


def _lowpass(x, fs, fc, order=2):
    sos = signal.butter(order, fc, btype="lowpass", fs=fs, output="sos")
    return _sos_causal(sos, x)


def stat_env(x, fs, fc=8.0, hold_s=0.03):
    """Production feature (|band-pass| -> LP fc) in baseline-SD units, held hold_s."""
    env = _lowpass(np.abs(bandpass(x, fs)), fs, fc)
    return trailing_min(standardise(env, fs), max(1, round(hold_s * fs)), round(BASELINE_S * fs))


def stat_tkeo(x, fs, fc=20.0, hold_s=0.01):
    y = bandpass(x, fs)
    psi = np.zeros_like(y)
    psi[:, 2:] = y[:, 1:-1] ** 2 - y[:, :-2] * y[:, 2:]  # psi at n uses n-2, n-1, n (causal)
    env = _lowpass(np.abs(psi), fs, fc)
    return trailing_min(standardise(env, fs), max(1, round(hold_s * fs)), round(BASELINE_S * fs))


def stat_rms(x, fs, window_s=0.025, hold_s=0.0):
    """RMS of the baseline-scaled band-pass over the last ``window_s`` (no IIR envelope)."""
    z = scaled_bp(x, fs)
    w = round(window_s * fs)
    q = np.concatenate([np.zeros((z.shape[0], 1)), np.cumsum(z**2, axis=1)], axis=1)
    r = np.zeros(z.shape)
    r[:, w - 1:] = np.sqrt(np.maximum(q[:, w:] - q[:, :-w], 0.0) / w)
    return trailing_min(r, max(1, round(hold_s * fs)), round(BASELINE_S * fs) + w - 1)


def stat_cusum(x, fs, rho=2.0):
    z = scaled_bp(x, fs)
    llr = 0.5 * (1.0 - 1.0 / rho**2) * z**2 - np.log(rho)
    n0 = round(BASELINE_S * fs)
    out = np.full(z.shape, -np.inf)
    c = np.cumsum(llr[:, n0:], axis=1)
    out[:, n0:] = c - np.minimum(np.minimum.accumulate(c, axis=1), 0.0)
    return out


def stat_glr(x, fs, window_s=0.1, step=4):
    n0 = round(BASELINE_S * fs)
    z = scaled_bp(x, fs)
    zp = z[:, n0:]  # post-baseline samples only
    q = np.concatenate([np.zeros((z.shape[0], 1)), np.cumsum(zp**2, axis=1)], axis=1)
    best = np.zeros(zp.shape)
    for L in range(step, min(round(window_s * fs), zp.shape[1]) + 1, step):
        s = q[:, L:] - q[:, :-L]
        s /= L
        np.maximum(s, 1.0, out=s)  # s <= 1 gives 0
        g = np.log(s)
        np.subtract(s, g, out=g)
        g -= 1.0
        g *= 0.5 * L
        np.maximum(best[:, L - 1:], g, out=best[:, L - 1:])
    out = np.full(z.shape, -np.inf)
    out[:, n0:] = best
    return out


def stat_bocpd(x, fs, block_s=0.002, hazard_s=2.0, r_max_s=0.2, alpha0=1.5, beta0=2.0, v_min=1.5):
    """Truncated BOCPD on block energies; returns a per-sample statistic (held between block ends).

    State per segment: posterior log mass of 'rest' (variance 1, known from the baseline scaling) and
    of runs r = 0..R-1, where run r has absorbed r + 1 blocks since a change to an unknown variance
    v ~ InvGamma(alpha0, beta0). A block of m samples has energy e; (2 pi)^(-m/2) is common to every
    hypothesis and cancels. Unlike Adams-MacKay Algorithm 1 step 5, the block at a change is scored
    by the new segment's prior predictive, so a run is never empty. Starts at the baseline end in
    'rest' with probability 1. A run that would exceed R blocks is merged into rest.
    """
    z = scaled_bp(x, fs)
    m = max(1, round(block_s * fs))
    n0 = round(BASELINE_S * fs)
    nb = (z.shape[1] - n0) // m
    half_e = 0.5 * (z[:, n0: n0 + nb * m] ** 2).reshape(z.shape[0], nb, m).sum(axis=2)
    R = max(3, round(r_max_s / block_s))
    h = 1.0 / max(1.0, hazard_s / block_s)
    lh, l1h = np.log(h), np.log1p(-h)
    al = alpha0 + (np.arange(R) + 1) * m / 2.0  # alpha of run r: it holds r + 1 blocks
    lg = gammaln(al + m / 2.0) - gammaln(al)
    lg0 = gammaln(alpha0 + m / 2.0) - gammaln(alpha0) + alpha0 * np.log(beta0)
    ns = z.shape[0]
    log_rest = np.zeros(ns)
    log_run = np.full((ns, R), -np.inf)
    beta = np.full((ns, R), beta0)
    log_beta = np.log(beta)
    out_b = np.zeros((ns, nb))
    for t in range(nb):
        he = half_e[:, t]
        beta_new = beta + he[:, None]
        lb_new = np.log(beta_new)
        pred = lg + al * log_beta - (al + m / 2.0) * lb_new  # log predictive of this block, per run
        grown = log_run + l1h + pred
        merged = log_run[:, -1] + l1h - he  # the oldest run leaves the window: scored as rest
        new_run = np.empty_like(log_run)
        new_run[:, 1:] = grown[:, :-1]
        new_run[:, 0] = lh + lg0 - (alpha0 + m / 2.0) * np.log(beta0 + he)  # total mass is 1 (log 0)
        new_rest = np.logaddexp(log_rest + l1h - he, merged)
        top = np.maximum(new_run.max(axis=1), new_rest)
        p_run = np.exp(new_run - top[:, None])
        p_rest = np.exp(new_rest - top)
        norm = top + np.log(p_run.sum(axis=1) + p_rest)
        log_run, log_rest = new_run - norm[:, None], new_rest - norm
        beta[:, 1:] = beta_new[:, :-1]
        beta[:, 0] = beta0 + he
        log_beta[:, 1:] = lb_new[:, :-1]
        log_beta[:, 0] = np.log(beta[:, 0])
        sel = beta > v_min * (al - 1.0)  # posterior-mean variance beta / (alpha - 1) > v_min
        sel[:, 0] = False  # at least 2 blocks in the run
        p_sel = (p_run * sel).sum(axis=1)
        out_b[:, t] = np.log(np.maximum(p_sel, 1e-300)) + top - norm - log_rest
    out = np.full(z.shape, -np.inf)
    # block t covers samples n0 + t m .. n0 + (t + 1) m - 1; its statistic is known at the last one
    held = np.repeat(out_b, m, axis=1)
    s0 = n0 + m - 1
    out[:, s0: s0 + nb * m] = held[:, : z.shape[1] - s0]
    return out


DETECTORS = {
    # name: (function, kwargs, description)
    "env8_h30": (stat_env, dict(fc=8.0, hold_s=0.03), "production trigger: 8 Hz envelope, 30 ms hold"),
    "env8_h10": (stat_env, dict(fc=8.0, hold_s=0.01), "8 Hz envelope, 10 ms hold"),
    "env20_h10": (stat_env, dict(fc=20.0, hold_s=0.01), "20 Hz envelope, 10 ms hold"),
    "env40_h5": (stat_env, dict(fc=40.0, hold_s=0.005), "40 Hz envelope, 5 ms hold"),
    "tkeo20_h10": (stat_tkeo, dict(fc=20.0, hold_s=0.01), "TKEO, 20 Hz envelope, 10 ms hold"),
    "tkeo40_h5": (stat_tkeo, dict(fc=40.0, hold_s=0.005), "TKEO, 40 Hz envelope, 5 ms hold"),
    "rms25": (stat_rms, dict(window_s=0.025), "25 ms moving RMS, no hold"),
    "rms50": (stat_rms, dict(window_s=0.05), "50 ms moving RMS, no hold"),
    "cusum_r2": (stat_cusum, dict(rho=2.0), "CUSUM, variance step design ratio rho = 2"),
    "cusum_r3": (stat_cusum, dict(rho=3.0), "CUSUM, design ratio rho = 3"),
    "glr_w50": (stat_glr, dict(window_s=0.05), "window-limited GLR, W = 50 ms"),
    "glr_w100": (stat_glr, dict(window_s=0.1), "window-limited GLR, W = 100 ms"),
    "bocpd": (stat_bocpd, dict(), "BOCPD, 2 ms blocks, R = 200 ms, hazard 1/2 s"),
}


def statistic(name, x, fs):
    fn, kw, _ = DETECTORS[name]
    return fn(np.atleast_2d(np.asarray(x, dtype=float)), fs, **kw)


def first_alarm(stat, h, start):
    """Index of the first sample >= start with stat > h, per row; -1 if none."""
    a = stat[:, start:] > h
    idx = np.argmax(a, axis=1)
    return np.where(a.any(axis=1), idx + start, -1)


def check_causal(name, x, fs, cut):
    """(max |stat(x)[:cut] - stat(x[:cut])|, same -inf mask): (0, True) for a causal statistic."""
    full = statistic(name, x, fs)[:, :cut]
    trunc = statistic(name, x[..., :cut], fs)
    fin = np.isfinite(full) & np.isfinite(trunc)
    same_inf = np.array_equal(np.isfinite(full), np.isfinite(trunc))
    return float(np.max(np.abs(full[fin] - trunc[fin]))) if fin.any() else 0.0, same_inf
