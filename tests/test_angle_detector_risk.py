"""Fast checks for scripts/angle_detector_risk.py (SYNTHETIC model; no recorded data).

The full Monte Carlo (about 100 s) is not run here; ``python scripts/angle_detector_risk.py --check``
reruns it and compares with the committed results and doc.
"""

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "angle_detector_risk", ROOT / "scripts" / "angle_detector_risk.py"
)
adr = importlib.util.module_from_spec(_spec)
sys.modules["angle_detector_risk"] = adr  # dataclasses need the module registered
_spec.loader.exec_module(adr)


@pytest.fixture(scope="module")
def ff():
    return adr.free_fall_points()


@pytest.fixture(scope="module")
def tr10():
    """True 10.000 deg peak with a 30 N*m arrest."""
    return adr.simulate(adr.t_on_for_peak(10.0, 30.0), 30.0)


def test_model_reproduces_pr9(ff):
    ref = adr.PR9_REF
    assert ff["J"] == pytest.approx(ref["J"], rel=1e-12)
    assert ff["G"] == pytest.approx(ref["G"], rel=1e-12)
    assert abs(ff["t10_ms"] - ref["t10_ms"]) < 1e-6
    assert abs(ff["drop200_deg"] - ref["drop200_deg"]) < 1e-8
    assert abs(ff["w10_dps"] - ref["w10_dps"]) < 1e-3


def test_catch_model_hits_target_peak(tr10):
    assert tr10.peak == pytest.approx(10.0, abs=1e-9)
    assert tr10.t_c > tr10.t_on
    # held at the peak after arrest
    assert tr10.drop(0.39) == pytest.approx(tr10.peak, abs=1e-9)


def test_time_rule_trips_below_10_deg():
    """An arm whose arrest completes at the window edge has a true peak below 10 deg (finite arrest)."""
    t_edge = adr.WINDOW_LO_S - adr.U_T_S
    for th in adr.T_ARREST:
        t_on = adr.brentq(lambda to: adr._solve(to, th, dense=False)[1] - t_edge, 0.02, 0.2019, xtol=1e-12)
        assert adr._solve(t_on, th, dense=False)[2] < 10.0


@pytest.mark.parametrize("ex,ey", [(0.0, 0.0), (0.03, -0.02), (-0.04, 0.01), (0.02, 0.04)])
def test_k_equals_dpsi_dphi(ex, ey):
    k, dpsi, _ = adr.k_ratio(ex, ey)
    assert dpsi == pytest.approx(k, rel=1e-7)


def test_disc_sweep_reproduces_revc_k_hi():
    r = adr.pr_sweep(adr.E_SRC)
    assert r["k_max"] == pytest.approx(adr.K_HI, abs=1e-9)
    assert round(r["k_min"], 2) == 0.41 and round(r["k_max"], 2) == 1.59  # Rev C 3.7 "k 0.41 to 1.59"
    assert round(adr.pr_sweep(adr.E_CAP, n_r=11, n_a=181)["k_max"], 2) == 1.27  # Rev C "(k 1.27)" at 46.3 mm


def test_ideal_chain_is_exact(tr10):
    ch = adr.Chain(q=0.0, lat=0.0)
    m = adr.measure(tr10, ch, np.random.default_rng(np.random.SeedSequence(1)), n=50)
    assert np.max(np.abs(m["peak"] - tr10.peak)) < 1e-6


def test_crossing_error_bounds(ff):
    """No quantization or skew: crossing error is latency plus at most one sample period."""
    tr = adr.simulate()
    ch = adr.Chain(q=0.0)
    m = adr.measure(tr, ch, np.random.default_rng(np.random.SeedSequence(2)), n=300)
    e = m["t10"] - ff["t10_ms"] / 1e3
    assert e.min() >= ch.lat - 1e-9
    assert e.max() <= ch.lat + 1 / ch.fs + 1e-9


def test_quantized_drop_is_whole_counts(tr10):
    ch = adr.Chain()
    m = adr.measure(tr10, ch, np.random.default_rng(np.random.SeedSequence(3)), n=100)
    counts = m["peak"] / ch.q
    assert np.allclose(counts, np.round(counts), atol=1e-6)
    assert np.max(np.abs(m["peak"] - tr10.peak)) < ch.q


def test_round1_chain_phase_mean(tr10):
    """Round-1 chain: realised sampling-phase mean only (full design check: test_r2_realised_components)."""
    ch = adr.Chain(q=0.0, lat=0.0, sigma=0.01)
    rng = np.random.default_rng(np.random.SeedSequence(4))
    m = adr.measure(adr.simulate(), ch, rng, n=400)
    assert np.mean(m["phase"]) * ch.fs == pytest.approx(0.5, abs=0.05)
    pre = adr.measure(adr.simulate(), adr.Chain(q=0.0, lat=0.0, sigma=0.01), np.random.default_rng(5), n=1)
    assert pre["peak"][0] > 0  # max of noise after onset is positive


def test_seed_children_distinct():
    ss = np.random.SeedSequence(adr.SEED)
    a = np.random.default_rng(ss.spawn(1)[0]).random(4)
    b = np.random.default_rng(ss.spawn(1)[0]).random(4)
    assert not np.allclose(a, b)


def test_score_threshold():
    g = adr.GUARD_DEG
    assert not adr.score_angle(np.array([10.0 - g - 1e-9]))[0]
    assert adr.score_angle(np.array([10.0 - g + 1e-9]))[0]


def test_stuck_from_onset_flagged_and_gap_flagged(tr10):
    rng = np.random.default_rng(np.random.SeedSequence(6))
    m = adr.measure(tr10, adr.Chain(stuck_at=0.0), rng, n=20)
    assert m["motion_bad"].all() and not adr.score_angle(m["peak"]).any()
    m = adr.measure(tr10, adr.Chain(gap=("peak", 2e-3)), rng, n=20)
    assert m["gap_bad"].all()


def test_device_joint_only_underreads_at_low_k(tr10):
    m = adr.measure(tr10, adr.Chain(k_dev=0.75), np.random.default_rng(7), n=20)
    assert np.all(m["peak"] < 7.6)
    m = adr.measure(tr10, adr.Chain(k_dev=0.75, pivot=True), np.random.default_rng(8), n=20)
    assert np.max(np.abs(m["peak"] - 10.0)) < 2 * adr.Q_REQ + 1e-9


def test_committed_results_render_to_committed_doc():
    R = json.loads((ROOT / "results" / "angle_detector_risk.json").read_text())
    doc = adr.render(R, (ROOT / "scripts" / "angle_detector_risk.template.md").read_text())
    assert doc == (ROOT / "docs" / "angle_detector_risk.md").read_text()
    assert "{{" not in doc


def test_committed_results_key_values():
    R = json.loads((ROOT / "results" / "angle_detector_risk.json").read_text())
    assert abs(R["model"]["d_t10_ms"]) < 1e-6
    assert R["crossing"]["worst_t10_abs_max_ms"] <= R["inputs"]["u_t_ms"]
    assert R["crossing"]["worst_d200_abs_max_deg"] <= R["inputs"]["guard"]
    assert R["peak_oc"]["req"]["max_missed"] == 0.0
    assert all(R["catch_time"][t]["peak_at_time_edge_deg"] < 10 for t in ("15", "20", "30", "50"))
    assert R["plan"]["p_for_5pct"]["n59"] == pytest.approx(1 - 0.95 ** (1 / 59))
    assert math.isclose(R["inputs"]["guard"], adr.GUARD_DEG)


def test_no_repo_detector_import():
    src = (ROOT / "scripts" / "angle_detector_risk.py").read_text()
    assert "emg_onset_lead" not in src.replace("# ", "")


# ---------------------------------------------------------------- fix round 2


def _rng(i):
    return np.random.default_rng(np.random.SeedSequence([2026, i + 1]))


def test_r2_realised_components():
    """Realised noise SD, skew, sampling phase and mount phase match the design (failures log 2026-11-22)."""
    ch = adr.Chain2(pivot=False, q=1e-9, sigma_counts=1e7)  # noise SD = 0.01 deg, quantization negligible
    m = adr.measure2(adr.simulate(), ch, _rng(1), 1000, ref=False)
    assert np.mean(m["base_sd"]) == pytest.approx(0.01, rel=0.01)
    e = ch.skew_uniform
    assert m["skew"].min() >= -e and m["skew"].max() <= e
    assert np.std(m["skew"]) == pytest.approx(e / math.sqrt(3), rel=0.08)
    assert np.mean(m["phase"]) * ch.fs == pytest.approx(0.5, abs=0.03)
    assert np.std(m["phase"]) * ch.fs == pytest.approx(1 / math.sqrt(12), rel=0.08)
    m = adr.measure2(adr.simulate(), adr.Chain2(pivot=False), _rng(2), 1000, ref=False)
    assert np.mean(m["mount_j"]) / adr.Q_REQ == pytest.approx(0.5, abs=0.03)


def test_r2_ideal_two_channel_is_exact(tr10):
    ch = adr.Chain2(k=adr.K_LO_CAP, q=1e-12, sigma_counts=0.0, lat=0.0)
    m = adr.measure2(tr10, ch, _rng(3), 50)
    assert np.max(np.abs(m["peak"] - tr10.peak)) < 1e-6
    assert np.max(np.abs(m["ref_peak"] - tr10.peak)) < 1e-6


def test_r2_joint_minus_pivot_within_two_counts(tr10):
    m = adr.measure2(tr10, adr.Chain2(k=adr.K_LO_CAP, sigma_counts=0.0), _rng(4), 500, ref=False)
    assert np.max(np.abs(m["peak"] - tr10.peak)) < 2 * adr.Q_REQ


def test_r2_pivot_stuck_reads_joint_only(tr10):
    """Pivot frozen from onset: forearm = joint channel = k * true drop."""
    ch = adr.Chain2(k=adr.K_LO_CAP, sigma_counts=0.0, piv=adr.Fault("stuck", t=0.0))
    m = adr.measure2(tr10, ch, _rng(5), 50, ref=False)
    assert np.max(np.abs(m["peak"] - adr.K_LO_CAP * tr10.peak)) < 2 * adr.Q_REQ


@pytest.mark.parametrize("chan", ["joint", "piv"])
def test_r2_slip_underreads_by_step(tr10, chan):
    f = adr.Fault("slip", size=0.2, t=0.05)
    ch = adr.Chain2(k=adr.K_LO_CAP, sigma_counts=0.0, **{chan: f})
    m = adr.measure2(tr10, ch, _rng(6), 50, ref=False)
    assert np.max(np.abs(m["peak"] - (tr10.peak - 0.2))) < 2 * adr.Q_REQ


def test_r2_backlash_underread_within_play(tr10):
    b = 0.2
    ch = adr.Chain2(pivot=False, q=1e-12, sigma_counts=0.0, joint=adr.Fault("backlash", size=b))
    m = adr.measure2(tr10, ch, _rng(7), 300, ref=False)
    under = tr10.peak - m["peak"]
    assert under.min() >= -1e-9 and under.max() <= b + 1e-9
    assert np.mean(under) == pytest.approx(b / 2, rel=0.15)


def test_r2_drift_underread(tr10):
    r = 0.2
    ch = adr.Chain2(pivot=False, q=1e-12, sigma_counts=0.0, joint=adr.Fault("drift", size=r))
    m = adr.measure2(tr10, ch, _rng(8), 20, ref=False)
    # drift counted from the start of the baseline window; baseline is its mean
    assert np.all(m["peak"] < tr10.peak - r * tr10.t_c * 0.9)


def test_r2_longest_run():
    c = np.array([[1, 1, 1, 2, 2, 3, 3, 3, 3, 4]])
    post = np.ones_like(c, bool)
    assert adr._longest_run(c, post)[0] == 3  # four equal values = three repeats
    post[0, :6] = False  # only samples 6-9 count: 3, 3, 3, 4 -> two repeats
    assert adr._longest_run(c, post)[0] == 2


def test_r2_play_operator():
    x = np.array([[0.0, 0.0, -1.0, -2.0, -1.0]])
    rng = np.random.default_rng(0)
    y = adr._play(x, 1.0, rng)
    assert np.all(np.abs(y - x) <= 0.5 + 1e-12)
    assert y[0, 3] == pytest.approx(-1.5)


def test_r2_budget_one_channel_within_exact_sup():
    r = adr.budget_two_channel(1.0, _rng(9), n=20000, pivot=False)
    ex = adr.exact_single_channel(n_phase=201, n_mount=101)
    assert r["t10_err_max_ms"] <= 0.14467 + 1e-4
    assert ex["t10_abs_max_ms"] <= 0.14467 + 1e-4 and ex["t10_abs_max_ms"] > 0.14
    assert ex["lattice_free_t10_ms"] == pytest.approx(
        (adr.ENC_LAT_S + adr.TRIG_SKEW_S + adr.REL_SKEW_S + 1 / adr.FS_REQ) * 1e3
        + adr.Q_REQ / adr.free_fall_points()["w10_dps"] * 1e3
    )


def test_r2_cp_interval():
    lo, hi = adr._cp(105, 400)
    assert lo < 105 / 400 < hi
    assert adr._cp(0, 10)[0] == 0.0 and adr._cp(10, 10)[1] == 1.0


def test_r2_committed_values():
    R = json.loads((ROOT / "results" / "angle_detector_risk.json").read_text())
    r2 = R["r2"]
    b = r2["budget"]
    assert b["one_channel"]["t10_ok"] and b["one_channel"]["d200_ok"]
    assert b["two_ch_d200_fails"]  # joint - pivot fails PR #9's guard (FC round 1 item 12)
    assert b["two_ch_worst_d200_deg"] < -adr.GUARD_DEG
    assert r2["exact"]["t10_abs_max_ms"] <= adr.U_T_S * 1e3
    assert r2["jp_peak_err"]["max"] < 2 * adr.Q_REQ and r2["jp_peak_err"]["min"] > -2 * adr.Q_REQ
    assert r2["faults"]["pivot_stuck_onset"]["missed_base"]["p"] > 0.9
    assert r2["pr9_head"].startswith("335d4e5")


# ---------------------------------------------------------------- fix round 3


def test_r3_one_channel_grid_matches_exact_sup():
    r = adr.exact_worst(1.0, pivot=False, nph=401, nm=128)
    assert 0.14 < r["t10_err_max_ms"] <= 0.14467 + 1e-4
    assert r["d200_err_min_deg"] == pytest.approx(-0.02355480210437655, abs=1e-6)


@pytest.mark.parametrize("k", [0.9, 1.3, 1.5862])
def test_r3_two_channel_d200_is_lattice_set_and_independent_of_k(k):
    r = adr.exact_worst(k, nph=201, nm=48)
    assert r["d200_err_min_deg"] == pytest.approx(-0.02630138413562655, abs=1e-6)
    assert not r["d200_ok"]


def test_r3_crossing_is_below_lattice_bound_and_above_one_channel():
    lb = adr.lattice_bound()
    assert lb == pytest.approx(0.16135, abs=2e-5) and lb > adr.U_T_S * 1e3
    for k in (0.7344681196147431, 1.2257, adr.K_HI):
        r = adr.exact_worst(k, nph=301, nm=48)
        assert r["t10_err_max_ms"] <= lb + 1e-9
        assert r["t10_err_max_ms"] > 0.14


def test_r3_lattice_bound_formula():
    w10 = adr.free_fall_points()["w10_dps"]
    n_up = math.ceil(adr.DROP_LIM / adr.Q_REQ)
    tpart = (1 / adr.FS_REQ + adr.ENC_LAT_S + adr.TRIG_SKEW_S + adr.REL_SKEW_S) * 1e3
    assert adr.lattice_bound() == pytest.approx(tpart + ((n_up + 1) * adr.Q_REQ - 10) / w10 * 1e3)
    assert adr.lattice_bound(pivot=False) == pytest.approx(tpart + (n_up * adr.Q_REQ - 10) / w10 * 1e3)
    assert adr.lattice_bound(pivot=False) < adr.U_T_S * 1e3  # the one-channel chain meets U_t by construction


def test_r3_lattice_bound_options():
    q = adr.Q_REQ
    # equal steps: (n+1) q; mixed 17/18 bit: finer lattice with slack qj+qp, which equals the 17/17 value
    assert adr.lattice_bound(q_p=q / 2) == pytest.approx(adr.lattice_bound())
    assert adr.lattice_bound(q_j=q / 2, q_p=q / 2) < adr.U_T_S * 1e3
    assert adr.lattice_bound(fs=20000.0) < adr.U_T_S * 1e3
    assert adr.lattice_bound(lat=0.0) < adr.lattice_bound()


def test_r3_k_rounding_can_flip_a_row():
    """k_lo vs its 4-digit rounding: same drop, different supremum (coarse grids)."""
    a = adr.exact_worst(adr.K_LO_CAP, nph=401, nm=64)["t10_err_max_ms"]
    b = adr.exact_worst(adr.K_ROUNDED[0], nph=401, nm=64)["t10_err_max_ms"]
    assert abs(a - b) > 1e-3


def test_r3_deceleration_grid_unit():
    """Measured deceleration is an integer number of grid units (q fs^2 / W^2), so a limit 'just above the
    maximum' is a draw-dependent grid value."""
    unit = adr.Q_REQ * adr.FS_REQ**2 / 15**2
    tr = adr.simulate(adr.t_on_for_peak(10.0, 30.0), 30.0)
    m = adr.measure2(tr, adr.Chain2(pivot=False, sigma_counts=1.0), _rng(20), 200, ref=False)
    u = m["max_dec"] / unit
    assert np.allclose(u, np.round(u), atol=1e-6)
    assert unit == pytest.approx(2746.58, abs=0.01)


def test_r3_loosening_bases():
    w10 = adr.free_fall_points()["w10_dps"]
    assert adr.U_T_S * 1e3 + adr.Q_REQ / adr.W10_MIN_PR9 * 1e3 == pytest.approx(0.17711, abs=2e-5)
    assert adr.U_T_S * 1e3 + adr.Q_REQ / w10 * 1e3 == pytest.approx(0.17663, abs=2e-5)


def test_r3_committed_values():
    R = json.loads((ROOT / "results" / "angle_detector_risk.json").read_text())
    r3 = R["r3"]
    ex = r3["exact"]
    assert not ex["k_hi_100mm"]["t10_ok"] and not ex["k_lo_rounded"]["t10_ok"] and ex["k_lo"]["t10_ok"]
    assert ex["one_channel"]["t10_ok"] and ex["one_channel"]["d200_ok"]
    assert (
        r3["sweep"]["n_t10_fail"] >= 1 and r3["sweep"]["t10_max_ms"] <= r3["lattice_bound_ms"]["two_channel"]
    )
    assert r3["t10_worst_ms"] <= r3["lattice_bound_ms"]["two_channel"] + 1e-9
    assert r3["loosen"]["u_t_pr9_basis_ms"] > r3["loosen"]["u_t_gov_ms"]
    for pool in ("clean_two", "clean_one"):
        c = r3[pool]
        assert c["max_units"] == pytest.approx(round(c["max_units"]), abs=1e-6)
        assert c["n_over_limit"] == 0 or c["over_limit"]["p"] > 0
    assert r3["revc"]["ratio_to_guard"] == pytest.approx(0.1 / adr.GUARD_DEG)
    ks = {round(v["dk"], 8): v for v in r3["k_sens"]}
    assert (
        ks[0.0]["t10_ok"] and not ks[round(adr.K_ROUNDED[0] - adr.K_LO_CAP, 8)]["t10_ok"]
    )  # a 3e-5 change in k flips it
    assert r3["limit_sens"]["range_41500"][0] >= r3["limit_sens"]["range_38500"][0]


def test_r3_doc_says_pr9_is_merged():
    doc = (ROOT / "docs" / "angle_detector_risk.md").read_text()
    assert "not merged" not in doc
    assert "merged to main at `335d4e5`" in doc
