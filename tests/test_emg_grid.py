"""EMG onset error across burst rise time x SNR, on Fact Checker's synthetic generator.

Synthetic data only. Every cell of rise 10/20/40/80/150 ms x SNR 3/5/15/30 must have its mean
EMG onset error and the whole 95% t-interval of that mean inside +/-6 ms (the protected margin).
The error is detected EMG onset minus the generator's true onset (FC's ``emg_onset_err``). On
main the same generator gave -10.7 ms (10 ms rise) to +35.8 ms (150 ms rise, SNR 3).

The IMU side is not in this test: it does not depend on the EMG cell (FC: -0.30 ms), and the
lead-error version of this grid is in /workspace/emgfix-grid.md.

Trials per cell: the per-trial SD is ~10-18 ms at SNR 3-5 and under ~3.5 ms at SNR 15-30. At
SNR 3-5 the cell needs ~200 trials before its CI half-width (~2.5 ms) leaves room for seed-to-
seed scatter of the mean; fewer trials pass or fail on the seed, not on the detector.
"""

import numpy as np
import pytest
from fc_synthetic import fc_trial
from scipy import stats

from emg_onset_lead.detect import EmgOnsetDetector, EmgParams

RISES_S = (0.01, 0.02, 0.04, 0.08, 0.15)
SNRS = (3.0, 5.0, 15.0, 30.0)
SEED = 1  # FC's seed. Each cell gets its own stream, default_rng([SEED, rise index, SNR index]).
LEAD_MS = 75.0
MARGIN_MS = 6.0
MAX_MISS = 0.10


def n_per_cell(snr):
    return 200 if snr < 10 else 80


def grid_cell(rise_s, snr, n=None, seed=SEED):
    """Run one cell. Returns (EMG onset errors in ms, number missed)."""
    n = n_per_cell(snr) if n is None else n
    rng = np.random.default_rng([seed, RISES_S.index(rise_s), SNRS.index(snr)])
    errs, missed = [], 0
    for _ in range(n):
        te_on, _, emg, _, fe, _ = fc_trial(rng, LEAD_MS, rise_s, snr)
        a = EmgOnsetDetector.run(emg, fe, EmgParams(), chunk=100)
        if a.onset_s is None:
            missed += 1
            continue
        errs.append((a.onset_s - te_on) * 1e3)
    return np.array(errs), missed


def mean_ci(x, conf=0.95):
    m = float(x.mean())
    h = float(stats.t.ppf(0.5 + conf / 2, x.size - 1) * x.std(ddof=1) / np.sqrt(x.size))
    return m, m - h, m + h


@pytest.mark.parametrize("snr", SNRS)
@pytest.mark.parametrize("rise_s", RISES_S)
def test_emg_onset_error_within_margin(rise_s, snr):
    errs, missed = grid_cell(rise_s, snr)
    assert missed <= MAX_MISS * n_per_cell(snr)
    m, lo, hi = mean_ci(errs)
    assert -MARGIN_MS < lo and hi < MARGIN_MS, f"mean {m:.2f} ms, 95% CI [{lo:.2f}, {hi:.2f}]"


def smoothstep_cell(rise_s, snr, n=80, seed=SEED):
    """Same as ``grid_cell`` but the amplitude follows SyntheticSource's smoothstep 3r^2 - 2r^3."""
    rng = np.random.default_rng([seed, 99, RISES_S.index(rise_s), SNRS.index(snr)])
    errs = []
    for _ in range(n):
        te_on = rng.uniform(1.5, 2.5)
        t = np.arange(int(4 * 2000.0)) / 2000.0
        r = np.clip((t - te_on) / rise_s, 0, 1)
        emg = (1 + (snr - 1) * (3 * r**2 - 2 * r**3)) * rng.normal(0, 1, t.size)
        a = EmgOnsetDetector.run(emg, 2000.0, EmgParams(), chunk=100)
        if a.onset_s is not None:
            errs.append((a.onset_s - te_on) * 1e3)
    return np.array(errs)


@pytest.mark.parametrize("rise_s", [
    pytest.param(r, marks=pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
        f"Known gap: the EMG onset model is a linear amplitude ramp. A smoothstep (S-curve) "
        f"burst of {r * 1e3:.0f} ms reads late by ~10% of its rise (synthetic).")))
    for r in (0.08, 0.15)
])
def test_emg_onset_smoothstep_burst(rise_s):
    m, lo, hi = mean_ci(smoothstep_cell(rise_s, 15.0))
    assert -MARGIN_MS < lo and hi < MARGIN_MS, f"mean {m:.2f} ms, 95% CI [{lo:.2f}, {hi:.2f}]"
