"""T10 `assist_arrived()` validity guard (FC PR #13 items #8, #9). ALL SYNTHETIC torque traces.

A run must not PASS when no assist arrives (zero, noise only, wrong sign, too small, buried in noise), when
the crossing is not after T0, or when the trace has not settled by T0 + 350 ms. Normal steps still pass,
and the committed docs/t10 outputs are reproduced byte for byte.
"""

import hashlib
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("t10_protocol", ROOT / "scripts" / "t10_protocol.py")
t10 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(t10)  # module level only defines; the simulation runs in main()

tm, T0, REPS = t10.tm, t10.T0, 10


def ramp(delay_ms, rise_ms, size=1.0):
    """Clean torque: 0 -> size, linear ramp starting delay_ms after T0."""
    return size * np.clip((tm - T0 - delay_ms / 1e3) / (rise_ms / 1e3), 0, 1)


def run(trace, noise_sd=0.0, seed=1):
    """REPS repeats of one trace plus independent sensor noise (SD noise_sd per repeat)."""
    g = np.random.default_rng(seed)
    return np.tile(trace, (REPS, 1)) + g.standard_normal((REPS, tm.size)) * noise_sd


def check_invalid(r, reason):
    assert r["passed"] is False
    assert r["status"] == "INVALID"
    assert any(reason in x for x in r["reasons"]), r["reasons"]


# ---- no assist: must not pass ----------------------------------------------------------------------------
def test_zero_torque_not_passed():
    r = t10.assist_arrived(np.zeros((REPS, tm.size)))
    check_invalid(r, "baseline SD")


def test_noise_only_not_passed():
    r = t10.assist_arrived(run(np.zeros(tm.size), noise_sd=0.02))  # FC #8 probe: SD 0.02, seed 1
    check_invalid(r, "baseline SD")


def test_reversed_step_not_passed():
    r = t10.assist_arrived(-run((tm >= T0 + 0.02) * 1.0, noise_sd=0.01))  # FC #8 probe: step to -1 at 20 ms
    check_invalid(r, "baseline SD")


def test_step_below_floor_not_passed():
    # clean, fast, settled, far above noise, but only 0.2 (N*m on the bench) < MIN_STEP 0.5
    r = t10.assist_arrived(run(ramp(20, 50, size=0.2), noise_sd=0.001))
    check_invalid(r, "floor")
    assert r["step"] > t10.K_NOISE * r["baseline_sd"] and r["t50"] > 0


def test_step_buried_in_noise_not_passed():
    # step 1.0 (above the floor) but the ensemble-mean baseline SD is about 0.2, so step/SD is about 5 < 10
    r = t10.assist_arrived(run(ramp(20, 1), noise_sd=0.2 * np.sqrt(REPS)))
    check_invalid(r, "baseline SD")
    assert r["step"] > t10.MIN_STEP and 0 < r["t50"] <= t10.T10_PASS_MS


def test_crossing_at_T0_not_passed():
    # torque already at its new level on the T0 sample: t50 = 0.00 ms, which no causal chain can give
    r = t10.assist_arrived(run((tm >= T0) * 1.0, noise_sd=0.001))
    check_invalid(r, "not after T0")
    assert r["t50"] == 0.0


# ---- late settling (FC #9) ---------------------------------------------------------------------------------
def test_fc_slow_ramp_invalid():
    # FC #9: 300 ms delay, 400 ms rise, true t50 500 ms; the 350-550 ms plateau window reads it as 375 ms
    r = t10.assist_arrived(np.tile(ramp(300, 400), (REPS, 1)))
    check_invalid(r, "not settled")
    assert r["t50"] == pytest.approx(375.0)


def test_two_stage_late_response_not_passed():
    # 40 % arrives at 20 ms, the other 60 % ramps from 300 to 700 ms: true t50 (half the final step) is
    # 366.7 ms, but against the unsettled plateau window t50 reads about 20 ms and would PASS unchecked
    trace = 0.4 * (tm >= T0 + 0.02) + ramp(300, 400, size=0.6)
    r = t10.assist_arrived(run(trace, noise_sd=0.01))
    check_invalid(r, "not settled")
    assert 0 < r["t50"] < 25


# ---- FC round 2 (#10-#13): what still slipped past the first guard -------------------------------------
def step_at(at_ms, size=1.0):
    return size * (tm >= T0 + at_ms / 1e3)


def test_late_change_after_plateau_window_invalid():
    # FC #10: 40 % at 20 ms, the other 60 % at 560 ms (after the 350-550 ms plateau window, in the record).
    # True t50 is 560 ms; the plateau sees only the first 40 %, so t50 read 20.5 ms and PASSED.
    r = t10.assist_arrived(run(step_at(20, 0.8) + step_at(560, 1.2), noise_sd=0.02))
    check_invalid(r, "after T0+550")
    assert 0 < r["t50"] < 25


def test_late_change_just_over_tolerance_invalid():
    # 20 % of the total arriving after the plateau window: tail 0.2 of the step, over SETTLE_TOL 0.15
    r = t10.assist_arrived(run(step_at(20, 0.8) + step_at(560, 0.2)))
    check_invalid(r, "after T0+550")


def test_late_pulse_inside_tail_window_invalid():
    # a +0.5 pulse from 560 to 590 ms that is gone again by 590 ms: only the whole 550-600 ms window sees it
    pulse = 0.5 * ((tm >= T0 + 0.56) & (tm < T0 + 0.59))
    r = t10.assist_arrived(run(step_at(20) + pulse))
    check_invalid(r, "after T0+550")


def test_dip_before_T0_invalid():
    # a dip of -0.3 of the step ending on the T0 sample: the offset counts in both directions
    dip = -0.3 * ((tm >= T0 - 0.03) & (tm < T0 + 0.01))
    r = t10.assist_arrived(run(step_at(20) + dip, noise_sd=0.005))
    check_invalid(r, "T0 sample")
    assert r["start_offset"] < 0


def test_small_late_change_still_valid():
    # 10 % of the total arriving after the plateau window is inside SETTLE_TOL (pins it from above)
    r = t10.assist_arrived(run(step_at(20, 0.9) + step_at(560, 0.1)))
    assert r["status"] == "PASS" and r["reasons"] == []


def test_rise_before_T0_invalid():
    # FC #13: torque starts 30 ms before T0 (ramp -30/100). It is 0.28 of the step on the T0 sample.
    r = t10.assist_arrived(run(ramp(-30, 100), noise_sd=0.01))
    check_invalid(r, "T0 sample")
    assert r["t50"] == pytest.approx(21.0, abs=1.0)


def test_small_rise_before_T0_still_valid():
    r = t10.assist_arrived(run(ramp(-5, 100), noise_sd=0.01))  # 0.05 of the step on the T0 sample
    assert r["status"] == "PASS" and r["reasons"] == []


def test_single_sample_spike_invalid():
    # FC #12: slow ramp (150 ms delay, 150 ms rise, true t50 225 ms) + a spike of 8x the step on 1 sample of
    # 1 repeat at 10 ms. The mean is 0.8 of the step there, so the first crossing was 10.0 ms and PASSED.
    x = run(ramp(150, 150), noise_sd=0.01)
    x[0, int(round((T0 + 0.01) * t10.FS))] = 8.0
    r = t10.assist_arrived(x)
    check_invalid(r, "does not hold")
    assert r["t50"] == pytest.approx(10.0)


def test_noisy_valid_step_not_rejected_by_crossing_check():
    # per-trial sensor noise 0.1 of the step, seed 15: the first crossing is early, noise only (pins
    # CONFIRM_LEVEL from above: a level of 1.0 makes this run INVALID)
    r = t10.assist_arrived(run(ramp(20, 100), noise_sd=0.1, seed=15))
    assert r["status"] == "PASS" and r["reasons"] == []


@pytest.mark.parametrize("kind, seed", [("independent", 218), ("common", 760)])
def test_random_walk_noise_not_passed(kind, seed):
    # FC #11: no assist, noise is a random walk (outside the Gaussian assumption). These two draws PASSED the
    # first guard (t50 49.0 and 77.0 ms); both are already off their baseline on the T0 sample.
    if kind == "independent":  # FC section E: cumulative sum of N(0,1) in every repeat, seed 1000 + i
        g = np.random.default_rng(1000 + seed)
        x = np.cumsum(g.standard_normal((REPS, tm.size)), axis=1)
    else:  # FC: the same random walk (x 0.05) in every repeat plus 0.001 noise, seed 5000 + i
        g = np.random.default_rng(5000 + seed)
        w = np.cumsum(g.standard_normal(tm.size)) * 0.05
        x = np.tile(w, (REPS, 1)) + g.standard_normal((REPS, tm.size)) * 0.001
    r = t10.assist_arrived(x)
    assert r["passed"] is False and r["status"] == "INVALID"
    check_invalid(r, "T0 sample")


# ---- normal behaviour unchanged ----------------------------------------------------------------------------
@pytest.mark.parametrize(
    "delay, rise, t50, status",
    [(20, 100, 70.5, "PASS"), (49, 100, 99.5, "PASS"), (52, 100, 102.5, "FAIL")],  # FC #7 probes
)
def test_normal_steps(delay, rise, t50, status):
    r = t10.assist_arrived(np.tile(ramp(delay, rise), (REPS, 1)))
    assert r["status"] == status and r["passed"] is (status == "PASS") and r["reasons"] == []
    assert r["t50"] == pytest.approx(t50)


def test_noisy_normal_step_passes():
    r = t10.assist_arrived(run(ramp(20, 100), noise_sd=0.05))
    assert r["status"] == "PASS" and 60 < r["t50"] < 80


def test_input_checks():
    with pytest.raises(ValueError):
        t10.assist_arrived(np.tile(ramp(20, 100), (9, 1)))  # fewer than 10 repeats
    with pytest.raises(ValueError):
        t10.assist_arrived(np.tile(ramp(20, 100)[:-1], (REPS, 1)))  # wrong time base
    bad = np.tile(ramp(20, 100), (REPS, 1))
    bad[0, 900] = np.nan
    with pytest.raises(ValueError):
        t10.assist_arrived(bad)


def test_committed_outputs_unchanged(tmp_path, capsys):
    """Regenerate docs/t10 (about 30 s): the guard must not change any committed number."""
    t10.rng = np.random.default_rng(7)  # same seed and draw order as a fresh `python scripts/t10_protocol.py`
    t10.main(str(tmp_path))
    for f in ("t10_protocol.md", "t10_protocol.json", "t10_criterion.md", "t10_criterion.json"):
        new = hashlib.sha256((tmp_path / f).read_bytes()).hexdigest()
        old = hashlib.sha256((ROOT / "docs" / "t10" / f).read_bytes()).hexdigest()
        assert new == old, f
    assert capsys.readouterr().out.count("invalid 0/1000") == 4  # no synthetic run is INVALID
