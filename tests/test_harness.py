import json
import math

import numpy as np
import pytest

from emg_onset_lead import (
    EmgOnsetDetector,
    Harness,
    ImuOnsetDetector,
    ImuParams,
    LatencyConfig,
    LeadStats,
    SyntheticSource,
    Trial,
    TrialSource,
)
from emg_onset_lead.cli import main

SRC = dict(n_subjects=4, n_sessions=3, n_trials=6, seed=3)


@pytest.fixture(scope="module")
def trials():
    return list(SyntheticSource(**SRC))


def _detectors(t):
    yield EmgOnsetDetector, t.emg, t.fs_emg, None
    yield ImuOnsetDetector, t.angle_or_imu, t.fs_imu, None
    yield ImuOnsetDetector, t.angle_or_imu, t.fs_imu, ImuParams(mode="angle")


def test_truncation_causality(trials):
    for t in trials[:12]:
        for cls, x, fs, p in _detectors(t):
            full = cls.run(x, fs, p)
            assert full.onset_s is not None
            confirm = round(full.detected_s * fs)  # sample at which the hold is satisfied
            cut = cls.run(x[:confirm + 1], fs, p)
            assert cut.onset_s == full.onset_s  # nothing after confirmation is needed
            # and nothing earlier is enough: one sample short of confirmation -> no onset
            assert cls.run(x[:confirm], fs, p).onset_s is None
            # corrupting the future must not change the answer
            y = x.copy()
            y[confirm + 1:] = 1e3
            assert cls.run(y, fs, p).onset_s == full.onset_s


@pytest.mark.parametrize("chunk", [1, 7, 64, 333])
def test_chunked_equals_oneshot(trials, chunk):
    for t in trials[:6]:
        for cls, x, fs, p in _detectors(t):
            assert cls.run(x, fs, p, chunk=chunk) == cls.run(x, fs, p)


def _errors(seed):
    res = Harness.run(SyntheticSource(n_subjects=8, n_sessions=3, n_trials=5, seed=seed))
    lead = np.array([r.lead_ms for r in res])
    assert not np.isnan(lead).any()
    return lead - np.array([r.true_lead_ms for r in res])


def test_recovers_known_lead():
    errs = [_errors(seed) for seed in range(6)]
    means = np.array([e.mean() for e in errs])
    pooled = np.concatenate(errs)
    assert abs(pooled.mean()) < 3.0  # pooled mean lead within +/-3 ms of ground truth
    assert np.all(np.abs(means) < 3.0)  # and every seed individually
    # no systematic sign bias: mean error is not significantly different from 0 ...
    assert abs(pooled.mean()) / (pooled.std(ddof=1) / np.sqrt(pooled.size)) < 3.0
    # ... and seeds are not all on one side
    assert 0 < int((means > 0).sum()) < len(means)


def _pooled(true_lead_ms):
    res = Harness.run(SyntheticSource(n_subjects=8, n_sessions=3, n_trials=4, seed=0,
                                      true_lead_ms=true_lead_ms))
    return LeadStats.summarize(res, n_boot=200)


@pytest.mark.parametrize(("true_lead", "expected"), [(40.0, False), (62.0, False), (120.0, True)])
def test_pass_fail_flips(true_lead, expected):
    # 62 ms is just below the 68.67 ms latency: a regression guard against measurement bias
    # inflating the lead (it used to be ~+10 ms, enough to pass).
    s = _pooled(true_lead)
    assert s["pass"] is expected
    if true_lead == 62.0:
        assert s["pooled"]["lower_bound_ms"] < s["pooled"]["required_ms"]


def test_latency_defaults():
    lat = LatencyConfig()
    assert lat.total_ms == 68.67 and lat.required_ms("imu") == 18.67
    res = Harness.run(SyntheticSource(n_subjects=1, n_sessions=1, n_trials=2))
    for r in res:
        assert math.isclose(r.margin_ms, r.lead_ms - 68.67)
        assert math.isclose(r.lead_ms, (r.imu_onset_s - r.emg_onset_s) * 1000)


def test_lower_bound_sanity():
    rng = np.random.default_rng(0)
    v = rng.normal(100, 10, 50)
    m, lb = LeadStats.t_lower_bound(v)
    assert lb < m == pytest.approx(v.mean())
    assert LeadStats.t_lower_bound(v, 0.99)[1] < lb  # higher confidence -> lower bound
    assert LeadStats.t_lower_bound([5.0])[1] is None
    # coverage: lower bound should sit below the true mean ~95% of the time
    hit = np.mean([LeadStats.t_lower_bound(rng.normal(0, 1, 8))[1] < 0 for _ in range(2000)])
    assert 0.92 < hit < 0.98


def test_bootstrap_reproducible_and_nested():
    nested = {f"s{i}": {"a": np.array([i + 100.0]), "b": np.array([i + 101.0])} for i in range(8)}
    a = LeadStats.boot_lower_bound(nested, seed=5, n_boot=300)
    assert a == LeadStats.boot_lower_bound(nested, seed=5, n_boot=300)
    assert a < np.mean([i + 100.5 for i in range(8)])


def test_summary_shape_and_levels():
    s = _pooled(80.0)
    assert len(s["sessions"]) == 24 and len(s["subjects"]) == 8
    assert s["study"]["matches_planned_shape"] is True
    assert s["pooled"]["bootstrap_lower_bound_ms"] is not None
    assert s["pooled"]["lower_bound_ms"] < s["pooled"]["mean_lead_ms"]


def test_seeded_reproducible():
    a = [t.emg[:50].tolist() for t in SyntheticSource(n_subjects=2, n_sessions=1, n_trials=2, seed=9)]
    b = [t.emg[:50].tolist() for t in SyntheticSource(n_subjects=2, n_sessions=1, n_trials=2, seed=9)]
    assert a == b


class ListSource:
    """Custom loader: shows the Protocol is all a dataset adapter needs."""

    def __init__(self, trials):
        self.trials = trials

    def __iter__(self):
        return iter(self.trials)


def test_custom_source_protocol(trials):
    src = ListSource([Trial("X", "X-1", "t1", t.fs_emg, t.emg, t.fs_imu, t.angle_or_imu)
                      for t in trials[:3]])
    assert isinstance(src, TrialSource)
    res = Harness.run(src)
    assert [r.subject_id for r in res] == ["X"] * 3
    assert all(math.isnan(r.true_lead_ms) for r in res)  # no truth supplied
    assert LeadStats.summarize(res, n_boot=50)["pooled"]["n_trials"] == 3


def test_cli_outputs(tmp_path):
    assert main(["run", "--synthetic", "--subjects", "2", "--sessions", "2", "--trials", "3",
                 "--seed", "0", "--out", str(tmp_path)]) == 0
    header = (tmp_path / "trials.csv").read_text().splitlines()[0]
    assert header.startswith("subject_id,session_id,trial_id,emg_onset_s,imu_onset_s,lead_ms,margin_ms")
    s = json.loads((tmp_path / "summary.json").read_text())
    assert {"sessions", "subjects", "pooled", "pass"} <= s.keys()
    assert s["study"]["matches_planned_shape"] is False
