import functools
import json
import math
from dataclasses import replace

import numpy as np
import pytest

from emg_onset_lead import (
    PLANNED_SHAPE,
    EmgOnsetDetector,
    Harness,
    ImuOnsetDetector,
    ImuParams,
    LatencyConfig,
    LeadStats,
    MinJerkSource,
    SampleSize,
    SyntheticSource,
    Trial,
    TrialResult,
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


@functools.lru_cache(maxsize=None)
def _errors(seed):
    res = Harness.run(SyntheticSource(n_subjects=8, n_sessions=3, n_trials=5, seed=seed))
    lead = np.array([r.lead_ms for r in res])
    assert not np.isnan(lead).any()
    return lead - np.array([r.true_lead_ms for r in res])


# SyntheticSource's EMG burst is a 40 ms smoothstep (S-curve) in amplitude. The EMG onset model
# (RampOnset) is a linear amplitude ramp, so it reads this burst about 4.00 ms late and the lead
# about 4.00 ms short (synthetic; scripts/synthetic_source_spread.py). The bound is the 6 ms
# margin; the sign bias is pinned below and kept visible by
# test_no_sign_bias_on_smoothstep_burst (strict xfail).
SMOOTHSTEP_LEAD_BIAS_MS = (-6.0, -1.5)


def test_recovers_known_lead():
    errs = [_errors(seed) for seed in range(6)]
    means = np.array([e.mean() for e in errs])
    pooled = np.concatenate(errs)
    assert abs(pooled.mean()) < 6.0  # pooled mean lead within the 6 ms margin of ground truth
    assert np.all(np.abs(means) < 6.0)  # and every seed individually
    lo, hi = SMOOTHSTEP_LEAD_BIAS_MS
    assert lo < pooled.mean() < hi  # the known shape bias, not something new


@pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
    "Known gap: the linear-ramp EMG onset model reads a 40 ms smoothstep burst ~4.00 ms late "
    "(synthetic), so every seed's lead error is negative. See README 'EMG onset'."))
def test_no_sign_bias_on_smoothstep_burst():
    errs = [_errors(seed) for seed in range(6)]
    means = np.array([e.mean() for e in errs])
    pooled = np.concatenate(errs)
    assert abs(pooled.mean()) < 3.0
    assert abs(pooled.mean()) / (pooled.std(ddof=1) / np.sqrt(pooled.size)) < 3.0
    assert 0 < int((means > 0).sum()) < len(means)


def _pooled(true_lead_ms):
    n_subj, n_sess = PLANNED_SHAPE
    res = Harness.run(SyntheticSource(n_subjects=n_subj, n_sessions=n_sess, n_trials=4, seed=0,
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
    assert len(s["sessions"]) == 48 and len(s["subjects"]) == 16
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
    assert s["study"]["plan_status"] == "unverified"
    assert s["demo"]["single_seed"] is True
    assert "confirmation_delay_ms" in s


# Mean |IMU onset - true motion start| must stay inside this on every sweep cell.
# On main's line-fit detector the same cells read about +15 to +86 ms, all late.
IMU_ONSET_TOL_MS = 5.0


def _imu_onset_errors(src):
    """IMU onset minus the true motion start (ms). This sweep is about the motion families,
    so the EMG side (and its known smoothstep-burst bias) is kept out of it."""
    res = Harness.run(src)
    t_move = np.array([t.meta["t_move_s"] for t in src])
    err = np.array([(r.imu_onset_s - tm) * 1000.0 for r, tm in zip(res, t_move)], dtype=float)
    assert np.mean(~np.isfinite(err)) < 0.15
    return err[np.isfinite(err)]


@pytest.mark.parametrize(("kind", "move_s", "amp", "noise"), [
    ("mj", 0.3, 60.0, 0.3),
    ("mj", 0.5, 60.0, 0.3),
    ("mj", 1.0, 60.0, 0.3),
    ("mj", 0.5, 10.0, 0.3),
    ("mj", 0.5, 120.0, 0.3),
    ("mj", 0.5, 60.0, 1.0),
    ("cd", 0.5, 60.0, 0.3),
    ("cd", 0.5, 60.0, 1.0),
])
def test_imu_onset_unbiased_across_motion(kind, move_s, amp, noise):
    """Duration 0.3/0.5/1.0 s, 10 deg to large, IMU noise to 1 deg, and the old shape.

    Measures the IMU onset against the true motion start. It measured the lead until the EMG
    onset model changed; the EMG side is now covered by tests/test_emg_grid.py.

    Covers only the in-family minimum-jerk and critically damped synthetic sweep. It does not
    show unbiasedness on out-of-family or real motion (Fact Checker round 2 found onset errors
    from 49 ms early to 271 ms late out of family).
    """
    means = []
    for seed in range(4):
        common = dict(n_subjects=3, n_sessions=2, n_trials=4, seed=seed,
                      angle_noise_deg=noise, true_lead_ms=80.0)
        src = (MinJerkSource(move_s=move_s, amplitude_deg=amp, **common) if kind == "mj"
               else SyntheticSource(**common))
        means.append(float(_imu_onset_errors(src).mean()))
    means = np.array(means)
    assert abs(float(means.mean())) < IMU_ONSET_TOL_MS
    # A real bias shows up as every seed on the same side of zero. A sub-ms mean may not.
    if abs(float(means.mean())) >= 1.0:
        assert np.any(means > 0) and np.any(means < 0)


def test_slope_nonpositive_fallback():
    """The line fit must take the threshold-crossing sample when the slope is not rising.

    This branch did not fire on the synthetic flexion runs; call it directly.
    """
    det = ImuOnsetDetector(150.0, ImuParams())
    det._thr = 5.0  # k = 5, so the standardised window equals the raw window
    det._hold_n = 5
    falling = np.array([4.0, 3.5, 3.0, 2.5, 2.0, 1.5])
    assert det._backdate(falling, 100) == 100 - 5 + 1
    flat = np.ones(6)
    assert det._backdate(flat, 80) == 80 - 5 + 1
    rising = np.array([0.2, 0.6, 1.2, 2.0, 3.2, 4.8])
    assert det._backdate(rising, 100) < 100 - 5 + 1


def _manual_rows(n_subj, n_sess, n_trials, lead=120.0, n_miss=0):
    """Spread misses across sessions so every session keeps at least one trial."""
    n_sess_total = n_subj * n_sess
    base, extra = divmod(n_miss, n_sess_total)
    rows = []
    session_i = 0
    for s in range(n_subj):
        for ss in range(n_sess):
            n_this = base + (1 if session_i < extra else 0)
            session_i += 1
            for t in range(n_trials):
                missed = t < n_this
                rows.append(TrialResult(
                    f"S{s:02d}", f"S{s:02d}-D{ss + 1}", f"T{t:02d}",
                    float("nan") if missed else 1.0,
                    float("nan") if missed else 1.0 + lead / 1000.0,
                    float("nan") if missed else lead,
                    float("nan") if missed else lead - 68.67,
                    baseline_rest_ok=not missed,
                ))
    return rows


def test_three_by_one_fails_shape_gate():
    s = LeadStats.summarize(_manual_rows(3, 1, 8, lead=120.0), n_boot=20)
    assert s["study"]["matches_planned_shape"] is False
    assert s["gates"]["lower_bound"] is True
    assert s["pass"] is False


def test_max_missed_trials_fails():
    # 16 x 3 x 10 = 480 trials, 268 missed (0.56 missed).
    # The trials that detected would clear 68.67 ms.
    s = LeadStats.summarize(_manual_rows(16, 3, 10, lead=120.0, n_miss=268), n_boot=20)
    assert s["pooled"]["n_total"] == 480
    assert s["pooled"]["n_missed"] == 268
    assert s["gates"]["lower_bound"] is True
    assert s["gates"]["missed_trials"] is False
    assert s["gates"]["planned_shape"] is True
    assert s["pass"] is False


def test_planned_shape_is_16_by_3():
    assert PLANNED_SHAPE == (16, 3)
    assert SyntheticSource().n_subjects == 16 and SyntheticSource().n_sessions == 3


def test_gate_accepts_16x3():
    s = LeadStats.summarize(_manual_rows(16, 3, 4, lead=120.0), n_boot=20)
    assert s["study"]["matches_planned_shape"] is True
    assert s["study"]["planned"]["subjects"] == 16
    assert s["study"]["planned"]["sessions_per_subject"] == 3
    assert s["gates"]["planned_shape"] is True
    assert s["pass"] is True


@pytest.mark.parametrize(("n_subj", "n_sess"), [(8, 3), (16, 2), (16, 4), (15, 3), (17, 3), (8, 6), (3, 1)])
def test_gate_rejects_wrong_shapes(n_subj, n_sess):
    s = LeadStats.summarize(_manual_rows(n_subj, n_sess, 4, lead=120.0), n_boot=20)
    assert s["gates"]["lower_bound"] is True
    assert s["gates"]["planned_shape"] is False
    assert s["study"]["matches_planned_shape"] is False
    assert s["pass"] is False


def _all_missed(rows, pred):
    """Mark every trial matching ``pred`` as missed (no detection, NaN lead)."""
    nan = float("nan")
    return [replace(r, emg_onset_s=nan, imu_onset_s=nan, lead_ms=nan, margin_ms=nan)
            if pred(r) else r for r in rows]


def test_shape_counts_fully_missed_subject():
    # 17x3 where one whole subject missed: still 17 subjects, not 16.
    rows = _all_missed(_manual_rows(17, 3, 10, lead=120.0), lambda r: r.subject_id == "S16")
    s = LeadStats.summarize(rows, n_boot=20)
    assert s["study"]["actual"]["subjects"] == 17
    assert s["gates"]["planned_shape"] is False and s["pass"] is False


def test_shape_counts_extra_subject_with_missed_session():
    # 16x3 plus a 17th subject whose only session missed every trial.
    extra = _all_missed(_manual_rows(17, 1, 10, lead=120.0)[160:], lambda r: True)
    assert {r.subject_id for r in extra} == {"S16"}
    rows = _manual_rows(16, 3, 10, lead=120.0) + extra
    s = LeadStats.summarize(rows, n_boot=20)
    assert s["gates"]["missed_trials"] is True  # 10/490 is under the cap; shape alone must fail
    assert s["gates"]["planned_shape"] is False and s["pass"] is False


def test_shape_counts_fully_missed_sessions():
    # 16x4 where every subject's 4th session missed: still 4 sessions per subject.
    rows = _all_missed(_manual_rows(16, 4, 10, lead=120.0), lambda r: r.session_id.endswith("D4"))
    s = LeadStats.summarize(rows, n_boot=20)
    assert s["study"]["actual"]["sessions_per_subject_max"] == 4
    assert s["gates"]["planned_shape"] is False and s["pass"] is False


def test_gate_rejects_uneven_sessions():
    rows = _manual_rows(16, 3, 4, lead=120.0)
    rows = [r for r in rows if not (r.subject_id == "S00" and r.session_id.endswith("D3"))]
    s = LeadStats.summarize(rows, n_boot=20)
    assert s["gates"]["planned_shape"] is False and s["pass"] is False


@pytest.mark.parametrize("t_move", [0.6, 0.9])
def test_baseline_movement_is_rejected(t_move):
    # Zero the lead jitter so the flexion itself starts at t_move, inside the 1.0 s baseline.
    src = SyntheticSource(n_subjects=1, n_sessions=1, n_trials=4, seed=0,
                          emg_onset_range_s=(t_move, t_move), true_lead_ms=0.0,
                          subject_sd_ms=0.0, session_sd_ms=0.0, trial_sd_ms=0.0)
    res = Harness.run(src)
    assert res
    assert all(r.baseline_rest_ok is False for r in res)
    assert all(math.isnan(r.lead_ms) and math.isnan(r.imu_onset_s) for r in res)


def test_triggered_mode_adds_confirmation_delay():
    lat = LatencyConfig()
    assert lat.triggered is False and lat.total_ms == 68.67
    # Proportional mode keeps 68.67 even if a delay is supplied.
    assert lat.required_ms("total", 51.5) == 68.67
    src = SyntheticSource(n_subjects=4, n_sessions=2, n_trials=5, seed=0)
    plain = LeadStats.summarize(Harness.run(src), lat, n_boot=30)
    delay = plain["confirmation_delay_ms"]
    # Hold time minus the back-dated onset (synthetic smoothstep burst). The 0.1 s refine wait
    # after the hold is not in it: a triggered controller fires at the hold.
    assert delay["mean"] == pytest.approx(47.0, abs=1.5)
    assert delay["min"] > 40.0  # the 30 ms hold is not the confirmation delay
    assert plain["pooled"]["required_ms"] == 68.67
    assert plain["latency"]["triggered"] is False
    trig = LatencyConfig(triggered=True)
    summary = LeadStats.summarize(Harness.run(src, trig), trig, n_boot=30)
    assert summary["latency"]["triggered"] is True
    assert summary["pooled"]["required_ms"] == pytest.approx(68.67 + delay["mean"], abs=0.05)
    assert summary["confirmation_delay_ms"]["mean"] == pytest.approx(delay["mean"], abs=1e-6)


def test_sample_size_helper_and_unverified_plan():
    # Synthetic variances (subject/session/trial SD 8/6/10 ms), not pilot data. They were chosen
    # so the rejected 8x3 design comes out near 50% power at budget + 6 ms.
    # between: subject sd 8 ms. within: a session mean, sd 6 ms plus 10 trials of sd 10 ms.
    power = SampleSize.power(between_var=64.0, within_var=46.0, effect_margin=6.0,
                             alpha=0.05, n_subjects=8, n_per_subject=3)
    assert power == pytest.approx(0.50, abs=0.03)
    needed = SampleSize.required_subjects(64.0, 46.0, 6.0, 0.05, 0.80, n_per_subject=3)
    assert needed is not None and needed > 8
    # At the same synthetic variances 16 x 3 clears 80%. Illustration only, not evidence for 16.
    assert SampleSize.power(64.0, 46.0, 6.0, 0.05, 16, 3) >= 0.80
    assert SampleSize.design_subject_sd(8.0, 6.0, 10.0, 3, 10) == pytest.approx(8.91, abs=0.01)
    plan = _pooled(80.0)
    assert plan["study"]["plan_status"] == "unverified"
    assert plan["study"]["matches_planned_shape"] is True
