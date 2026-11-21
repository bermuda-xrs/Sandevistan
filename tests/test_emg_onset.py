"""EMG onset detector: refine wait, trigger time, re-arm after a noise trigger, short bursts.

Synthetic only.
"""

import numpy as np
import pytest

from emg_onset_lead.detect import EmgOnsetDetector, EmgParams, RampOnset

FS = 2000.0


def _noise_trigger_then_burst(burst_at=3.2):
    """Pure noise with a real false trigger, then a real burst.

    ``default_rng([7, 259])`` white noise meets the k*SD + 30 ms hold at ~2.46 s by chance (found
    by search; the envelope peaks at only ~1.1x threshold, below short_burst_peak = 1.6x).
    A 40 ms-rise burst to x15 starts at ``burst_at``. The noise before ``burst_at`` is untouched.
    """
    rng = np.random.default_rng([7, 259])
    x = rng.normal(0, 1, round(4.0 * FS))
    t = np.arange(x.size) / FS
    x *= 1.0 + 14.0 * np.clip((t - burst_at) / 0.04, 0, 1)
    return x, burst_at


def _short_burst(dur=0.06, snr=5.0, burst_at=2.0, seed=0):
    """An isolated burst: 20 ms linear rise to ``snr``, back to baseline after ``dur``."""
    rng = np.random.default_rng(seed)
    t = np.arange(round(4.0 * FS)) / FS
    a = 1.0 + (snr - 1.0) * np.clip((t - burst_at) / 0.02, 0, 1)
    a = np.where(t >= burst_at + dur, 1.0, a)
    return a * rng.normal(0, 1, t.size), burst_at


def test_trigger_and_commit_times():
    x, true_on = _noise_trigger_then_burst()
    x[round(2.3 * FS):round(2.7 * FS)] = np.random.default_rng(1).normal(0, 1, round(0.4 * FS))
    r = EmgOnsetDetector.run(x, FS)
    assert r.trigger_s is not None and r.detected_s is not None
    assert r.detected_s - r.trigger_s == pytest.approx(EmgParams().refine_s, abs=1.0 / FS)
    assert r.onset_s < r.trigger_s
    assert abs(r.onset_s - true_on) < 0.006


def test_noise_trigger_is_dropped_and_detector_rearms():
    x, true_on = _noise_trigger_then_burst()
    # The noise alone does trigger the hold (otherwise this test shows nothing).
    det = EmgOnsetDetector(FS)
    det.process(x[:round(2.5 * FS)])
    assert det._pending is not None and abs(det._pending / FS - 2.46) < 0.01
    r = EmgOnsetDetector.run(x, FS)
    assert abs(r.onset_s - true_on) < 0.006  # the burst, not the noise trigger


def test_short_burst_is_kept_and_fit_ends_at_its_peak():
    """A 60 ms burst has decayed by the commit (hold + 0.1 s) but is a burst, not noise."""
    x, true_on = _short_burst()
    det = EmgOnsetDetector(FS)
    r = EmgOnsetDetector.run(x, FS)
    assert r.onset_s is not None and abs(r.onset_s - true_on) < 0.006
    det.process(x)
    commit = round(r.detected_s * FS)
    assert det._fit_end is not None and det._trig_idx <= det._fit_end < commit


@pytest.mark.parametrize("chunk", [1, 37, 100, 999])
@pytest.mark.parametrize("make", [_noise_trigger_then_burst, _short_burst])
def test_rearm_and_short_burst_are_chunk_invariant_and_causal(chunk, make):
    x, _ = make()
    full = EmgOnsetDetector.run(x, FS)
    assert EmgOnsetDetector.run(x, FS, chunk=chunk) == full
    commit = round(full.detected_s * FS)
    assert EmgOnsetDetector.run(x[:commit + 1], FS) == full
    assert EmgOnsetDetector.run(x[:commit], FS).onset_s is None
    y = x.copy()
    y[commit + 1:] = 1e3
    assert EmgOnsetDetector.run(y, FS) == full


def test_ramp_onset_on_clean_step():
    rng = np.random.default_rng(1)
    z = rng.normal(0, 1, 800)
    z[500:] *= 8.0
    det = EmgOnsetDetector(FS)
    est = RampOnset.estimate(z, FS, det._rise_grid, det._gain_grid)
    assert abs(est - 500) < 4  # 2 ms
    assert RampOnset.estimate(z[:10], FS, det._rise_grid, det._gain_grid) is None
