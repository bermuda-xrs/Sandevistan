"""Isolated short EMG bursts: detection rate and onset error. Synthetic data only.

Fact Checker's PR #4 review (row 17/25) found that head 3cbc1d7 detected 0/60 isolated 40-60 ms
bursts (20 ms rise) because the commit check required the envelope to still be above threshold
0.1 s after the hold. Base 608a151 detected 60/60 but with its rise/SNR-dependent bias. The burst
generator below has the same form as FC's ``fc4_harness.py`` ``shortburst`` job: white noise,
amplitude ramping linearly from 1 to SNR over 20 ms, back to 1 after the burst duration.

Trials per cell: 150 at SNR 3 and 5 (per-trial onset SD up to ~10 ms, so 60 trials leave the CI
to the seed), 60 at SNR 15 and 30. Streams: ``default_rng([SEED, SNR, duration ms])``.
Bounds:
- detection: at least 96% in every cell (58/60, 144/150), except 40 ms at SNR 3: at least 65%.
  A 40 ms SNR-3 burst often never meets the 30 ms hold, or peaks below short_burst_peak = 1.6x
  threshold. There this branch detects 121/150 on the test's seed and 107-121/150 (71-81%) over
  seeds 0-9; 243/300 at 300 trials (FC: 234/300). Base detects ~97% (58/60). The 65% bound
  (98/150) is ~4.6 binomial SD below the test seed's 121, so it is not flaky but still catches a
  real loss of short-burst detection. Every other cell measured >= 149/150 (59/60 at worst at
  SNR 15-30) over seeds 0-9. All synthetic.
- onset error of the detected trials: mean and 95% t-interval inside +/-6 ms.
"""

import numpy as np
import pytest
from scipy import stats

from emg_onset_lead.detect import EmgOnsetDetector, EmgParams

FS = 2000.0
SNRS = (3.0, 5.0, 15.0, 30.0)
DURS_S = (0.04, 0.06, 0.08)
RISE_S = 0.02
SEED = 9
MARGIN_MS = 6.0


def n_trials(snr):
    return 150 if snr < 10 else 60


def min_detected_frac(snr, dur):
    return 0.65 if (snr < 4 and dur < 0.05) else 0.96


def short_burst_cell(snr, dur, n=None, seed=SEED):
    """Onset errors (ms) of detected bursts; a detection more than 200 ms away is not the burst."""
    n = n_trials(snr) if n is None else n
    rng = np.random.default_rng([seed, int(snr), int(round(dur * 1000))])
    errs = []
    for _ in range(n):
        t_on = rng.uniform(1.5, 2.5)
        t = np.arange(int(4 * FS)) / FS
        a = 1 + (snr - 1) * np.clip((t - t_on) / RISE_S, 0, 1)
        a = np.where(t >= t_on + dur, 1.0, a)
        r = EmgOnsetDetector.run(a * rng.normal(0, 1, t.size), FS, EmgParams(), chunk=100)
        if r.onset_s is not None and abs(r.onset_s - t_on) < 0.2:
            errs.append((r.onset_s - t_on) * 1e3)
    return np.array(errs)


@pytest.mark.parametrize("dur", DURS_S)
@pytest.mark.parametrize("snr", SNRS)
def test_short_isolated_burst(snr, dur):
    errs = short_burst_cell(snr, dur)
    n = n_trials(snr)
    assert errs.size >= min_detected_frac(snr, dur) * n, f"detected {errs.size}/{n}"
    m = float(errs.mean())
    h = float(stats.t.ppf(0.975, errs.size - 1) * errs.std(ddof=1) / np.sqrt(errs.size))
    assert -MARGIN_MS < m - h and m + h < MARGIN_MS, f"mean {m:.2f}, CI [{m - h:.2f}, {m + h:.2f}]"
