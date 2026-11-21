"""Fact Checker's synthetic generators, ported for the tests. Synthetic data only.

``fc_trial`` is ``trial()`` from /workspace/fc_check/bias11/rerun_bias11.py (the source of the
rise-time/SNR EMG bias report). ``MOTION_SHAPES`` and ``fc_motion_trial`` come from
/workspace/fc_check/check_pr2.py (round 2, out-of-family motion). Only names and layout
changed; the signal definitions and constants are FC's.
"""

import numpy as np


def mj(s):
    s = np.clip(s, 0, 1)
    return 10 * s**3 - 15 * s**4 + 6 * s**5


def fc_trial(rng, lead_ms, rise_s, snr, move_T=0.5, span=60.0, fs_e=2000.0, fs_i=150.0, noise=0.3):
    """EMG: white noise whose amplitude ramps linearly from 1 to ``snr`` over ``rise_s``.
    Angle: minimum-jerk move of ``span`` deg over ``move_T`` s, ``lead_ms`` after the EMG onset."""
    te_on = rng.uniform(1.5, 2.5)
    tm = te_on + lead_ms / 1e3
    te = np.arange(int(4 * fs_e)) / fs_e
    r = np.clip((te - te_on) / rise_s, 0, 1)
    emg = (1 + (snr - 1) * r) * rng.normal(0, 1, te.size)
    ti = np.arange(int(4 * fs_i)) / fs_i
    ang = 20 + span * mj((ti - tm) / move_T) + rng.normal(0, noise, ti.size)
    return te_on, tm, emg, ang, fs_e, fs_i


# angle(t - t_move) above rest, deg, for t >= t_move (check_pr2.py SHAPES, out-of-family subset)
MOTION_SHAPES = {
    "skewed_bell_late_peak(0.5s)": lambda d: 60 * mj(np.clip(d / .5, 0, 1) ** 1.8),
    "slow_creep_then_fast": lambda d: 3 * mj(d / .3) + 57 * mj((d - .3) / .4),
    "cosine_ramp_0.5s": lambda d: 60 * (1 - np.cos(np.pi * np.clip(d / .5, 0, 1))) / 2,
    "linear_ramp_0.5s": lambda d: 60 * np.clip(d / .5, 0, 1),
}


def fc_motion_trial(rng, shape, lead_ms=75.0, noise=0.3, fs_e=2000.0, fs_i=150.0):
    """check_pr2.py ``bias()`` trial: 5 s, EMG 40 ms linear rise to SNR 15, angle from ``shape``."""
    t_on = rng.uniform(1.5, 2.5)
    tm = t_on + lead_ms / 1e3
    ti = np.arange(int(5 * fs_i)) / fs_i
    ang = 20 + shape(np.clip(ti - tm, 0, None)) + rng.normal(0, noise, ti.size)
    te = np.arange(int(5 * fs_e)) / fs_e
    r = np.clip((te - t_on) / .04, 0, 1)
    emg = (1 + 14 * r) * rng.normal(0, 1, te.size)
    return t_on, tm, emg, ang, fs_e, fs_i
