"""EMG onset detector: refine wait, trigger time, re-arm after a decayed trigger. Synthetic only."""

import numpy as np
import pytest

from emg_onset_lead.detect import EmgOnsetDetector, EmgParams, RampOnset

FS = 2000.0


def _blip_then_burst(seed=0, blip_at=1.5, burst_at=2.5):
    """Baseline noise, a 60 ms x5 noise blip (meets the 30 ms hold, then decays), a real burst."""
    rng = np.random.default_rng(seed)
    t = np.arange(round(4.0 * FS)) / FS
    amp = np.ones_like(t)
    amp[(t >= blip_at) & (t < blip_at + 0.06)] = 5.0
    amp += 14.0 * np.clip((t - burst_at) / 0.04, 0, 1)
    return amp * rng.normal(0, 1, t.size), burst_at


def test_trigger_and_commit_times():
    x, true_on = _blip_then_burst()
    x[round(1.4 * FS):round(1.7 * FS)] /= 5.0  # no blip in this one
    r = EmgOnsetDetector.run(x, FS)
    assert r.trigger_s is not None and r.detected_s is not None
    assert r.detected_s - r.trigger_s == pytest.approx(EmgParams().refine_s, abs=1.0 / FS)
    assert r.onset_s < r.trigger_s
    assert abs(r.onset_s - true_on) < 0.006


def test_decayed_trigger_is_dropped_and_detector_rearms():
    x, true_on = _blip_then_burst()
    r = EmgOnsetDetector.run(x, FS)
    assert abs(r.onset_s - true_on) < 0.006  # not the blip at 1.5 s
    # The blip alone does trigger the hold (otherwise this test shows nothing).
    det = EmgOnsetDetector(FS)
    det.process(x[:round(1.62 * FS)])
    assert det._pending is not None and abs(det._pending / FS - 1.5) < 0.1


@pytest.mark.parametrize("chunk", [1, 37, 100, 999])
def test_rearm_is_chunk_invariant_and_causal(chunk):
    x, _ = _blip_then_burst(seed=3)
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
