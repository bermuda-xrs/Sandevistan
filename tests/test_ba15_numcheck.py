"""scripts/ba15_numcheck.py and docs/ba15-testplan.md. Drop-model numbers are SYNTHETIC."""

import importlib.util
import math
from fractions import Fraction as Fr
from itertools import product
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.stats import beta, binom

from emg_onset_lead.detect import EmgParams
from emg_onset_lead.latency import LatencyConfig
from emg_onset_lead.sources import SyntheticSource

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "ba15_numcheck.py"
_spec = importlib.util.spec_from_file_location("ba15_numcheck", _PATH)
nc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nc)

TEXT = nc.TEMPLATE.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def R():
    return nc.compute()


def test_committed_plan_and_numbers_reproduce(R):
    assert nc.check(R) == []


def test_check_catches_edited_plan_and_numbers(R, tmp_path):
    plan = tmp_path / "plan.md"
    plan.write_text(nc.PLAN.read_text(encoding="utf-8").replace("201.90", "201.91", 1), encoding="utf-8")
    assert any("rendered template" in b for b in nc.check(R, plan=plan))
    numbers = tmp_path / "numbers.json"
    numbers.write_text(nc.NUMBERS.read_text(encoding="utf-8").replace('"n0": 59', '"n0": 60'),
                       encoding="utf-8")
    assert any("fresh run" in b for b in nc.check(R, numbers=numbers))


def test_template_guard_catches_literals_and_missing_labels():
    assert nc.template_problems(TEXT) == []
    assert nc.template_problems("need 61 trials\n")
    assert nc.template_problems("edge {{syn.lo|.2f}} ms\n") and nc.template_problems("{{in.win_lo}}\n")
    assert not nc.template_problems("edge {{syn.lo|.2f}} ms (SYNTHETIC)\n")


def test_synthetic_tables_are_labelled(R):
    for name in ("feas", "feas_icc", "pru", "budget", "absorber", "flip", "matrix"):
        assert "SYNTHETIC" in nc.tables(R)[name].splitlines()[0]


def _upper(x, n):
    """Clopper-Pearson upper bound by bisection on a plain binomial sum."""
    lo, hi = 0.0, 1.0
    for _ in range(80):
        p = (lo + hi) / 2
        tail = sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(x + 1))
        lo, hi = (p, hi) if tail > 0.05 else (lo, p)
    return lo


@pytest.mark.parametrize("c,n", [(0, 59), (1, 93), (2, 124)])
def test_allowed_drop_n_is_smallest(c, n):
    assert nc.n_allowed(c) == n
    assert _upper(c, n) <= 0.05 < _upper(c, n - 1)
    assert nc.cp_upper(c, n) == pytest.approx(_upper(c, n), rel=1e-9)


def _rk4(rate, sto, td, hold, t_end=None, J=None, G=None, k=None, T0=None, ang=None, dt=1e-5):
    """Independent fixed-step RK4 in the absolute elbow angle (90 deg at start).
    Returns the time (ms) at which the drop reaches ``ang`` deg, or the drop (deg) at ``t_end``."""
    A = nc.arm()
    J, G = J or A["J"], G or A["G"]
    k, T0 = k or nc.K_HI, T0 or nc.T0_DEV
    brk = round(sto / dt)

    def acc(th, t, off):
        T = 0.0 if off else (T0 if t < td else max(T0 - rate * (t - td), 0.0))
        return (k * T + hold - G * math.sin(th)) / J

    th, w, prev = math.pi / 2, 0.0, 0.0
    for i in range(round((t_end or 0.4) / dt)):
        t, off = i * dt, i >= brk
        k1 = (w, acc(th, t, off))
        k2 = (w + dt / 2 * k1[1], acc(th + dt / 2 * k1[0], t + dt / 2, off))
        k3 = (w + dt / 2 * k2[1], acc(th + dt / 2 * k2[0], t + dt / 2, off))
        k4 = (w + dt * k3[1], acc(th + dt * k3[0], t + dt, off))
        th += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        w += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        d = 90 - math.degrees(th)
        if ang is not None and d >= ang:
            return 1e3 * (t + dt * (ang - prev) / (d - prev))
        prev = d
    assert ang is None
    return 90 - math.degrees(th)


def _rk4_crossing(rate, sto, td, hold):
    return _rk4(rate, sto, td, hold, ang=nc.DROP_LIM)


def test_window_matches_decision_and_independent_integrator(R):
    syn = R["syn"]
    assert (round(syn["lo"], 2), round(syn["hi"], 2)) == nc.WINDOW_MS
    assert _rk4_crossing(21.0, 0.148, 0.0, 0.0) == pytest.approx(syn["lo"], abs=1e-3)
    assert _rk4_crossing(19.95, 0.142, 0.003, 0.0) == pytest.approx(syn["hi"], abs=1e-3)
    assert syn["kT0"] == pytest.approx(10.0)


def test_earlier_subject_estimate_and_fc_range_reproduce(R):
    assert R["stat"]["range_he"] == [6, 17, 33]
    assert (R["syn"]["fc_p_lo"], R["syn"]["fc_p_hi"]) == pytest.approx(nc.REF["fc_p"], abs=1e-7)
    assert [round(R["syn"][k], 1) for k in ("fc_p_lo", "fc_p_hi")] == [0.2, 0.3]


def test_subjects_ceiling_is_exact():
    assert nc.subjects(59, 10, Fr(1)) == 59
    assert nc.subjects(60, 10, Fr(0)) == 6
    assert nc.subjects(59, 1, Fr("0.5")) == 59


@pytest.mark.parametrize("icc", [0.05, 0.2, 0.5])
def test_beta_binomial_matches_integral(icc):
    p, m, f = 0.05, 10, (1 - icc) / icc
    d = beta(p * f, (1 - p) * f)
    # x = u**8 removes the integrable singularity at 0
    want = quad(lambda u: 8 * u**7 * d.pdf(u**8) * (1 - u**8) ** m if u > 0 else 0.0, 0, 1, limit=400)[0]
    assert nc.pmf_bb(p, icc, m, 0)[0] == pytest.approx(want, rel=1e-5)


@pytest.mark.parametrize("pmf", [nc.pmf_bb, nc.pmf_two])
@pytest.mark.parametrize("icc", [0, 0.1, 0.5, 1])
def test_cluster_models_have_the_requested_mean_and_icc(pmf, icc):
    m, p = 10, 0.05
    q, ks = pmf(p, icc, m, m), np.arange(m + 1)
    mean = float(q @ ks)
    assert q.sum() == pytest.approx(1)
    assert mean == pytest.approx(m * p)
    var = float(q @ ks**2) - mean**2
    assert var == pytest.approx(m * p * (1 - p) * (1 + (m - 1) * icc), rel=1e-9)


@pytest.mark.parametrize("icc", [Fr("0.05"), Fr("0.2"), Fr("0.5"), Fr("0.8")])
def test_worst_case_is_the_never_or_sometimes_mix(icc):
    assert nc.worst_zero(0.05, icc, 10) == pytest.approx(nc.pmf_two(0.05, icc, 10, 0)[0], abs=1e-7)


def test_design_effect_keeps_false_pass_below_5_percent(R):
    for v in R["stat"]["grid"].values():
        assert v["false_pass_worst"] <= 0.05 and v["subjects"] >= max(v["bb"], v["worst"])
    assert R["stat"]["allowed_short"] == []


def test_staged_pass_matches_enumeration():
    p, ns = 0.2, (3, 5, 7)
    total = 0.0
    for seq in product((0, 1), repeat=7):
        if any(sum(seq[:n]) <= i for i, n in enumerate(ns)):
            total += p ** sum(seq) * (1 - p) ** (7 - sum(seq))
    assert nc.staged_pass(ns, p) == pytest.approx(total)
    assert nc.staged_pass_subjects(ns, binom.pmf(range(3), 1, p)) == pytest.approx(total)


def test_extension_controls_false_pass_and_naive_extension_does_not(R):
    s = R["stat"]
    assert s["ext"][0] == 59 and s["ext_size"] <= 0.05 and s["ext_clustered_max"] <= 0.05
    assert nc.staged_pass([59, s["ext"][1] - 1], 0.05) > 0.05
    assert s["naive1"] > 0.07 and s["naive2"] > 0.09


def test_feasibility_is_near_zero(R):
    f = R["syn"]["feas"]
    assert f["indep"]["FC low"][0] == pytest.approx(0.8**59)
    assert f["indep"]["FC high"][0] == pytest.approx(0.67**59)
    assert f["max_pass"] < 1e-4
    assert nc.staged_pass([59], R["stat"]["p80"]) == pytest.approx(0.8)


def test_protected_values_are_untouched():
    assert LatencyConfig().total_ms == pytest.approx(68.67)
    assert EmgParams().hold_s == 0.03 and EmgParams().baseline_s == 1.0 and EmgParams().rest_ratio == 12.0
    src = SyntheticSource()
    assert (src.fs_emg, src.fs_imu) == (nc.EMG_HZ, nc.IMU_HZ)


def test_plan_states_ethics_and_assumptions():
    plan = nc.PLAN.read_text(encoding="utf-8")
    assert "Human testing needs ethics (IRB) approval before any trial" in plan
    assert "| A1 | The 10° drop limit" in plan and "| A2 | The 200 ms catch" in plan


def test_grid_has_every_cell_and_matches_independent_rk4(R):
    grid = nc.build_grid()
    assert len(grid) == 117 and len({k[:4] for k in grid}) == 39
    assert min(v["t10"] for v in grid.values()) == pytest.approx(R["syn"]["lo"])
    A = nc.arm()
    for key in [("b", "estop", 21.0, 148.0, 0.0), ("iv", "comm", 19.95, 142.0, 0.0),
                ("b", "emg", 50.0, None, 1.0), ("iv", "estop", 18.05, 160.0, 0.5)]:
        o, path, rate, sto, hf = key
        v = next(x for kk, x in grid.items() if kk[0] == o and kk[1] == path and kk[3] == sto
                 and kk[4] == hf and kk[2] == pytest.approx(rate))
        hold, td = hf * A["hold_nom"], nc.T_DET[path] / 1e3
        sto_s = 1e9 if sto is None else sto / 1e3
        assert _rk4(rate, sto_s, td, hold, ang=nc.DROP_LIM) == pytest.approx(v["t10"], abs=2e-3)
        assert _rk4(rate, sto_s, td, hold, t_end=0.2) == pytest.approx(v["drop200"], abs=2e-3)


def test_matrix_counts_and_rows(R):
    m = R["syn"]["hw"]["matrix"]
    assert (m["n_cells"], m["n_types"], m["n_b"], m["n_iv"]) == (117, 39, 21, 18)
    assert (m["n_cell"], m["drops"], m["hours"]) == (10, 1170, 19.5)
    assert m["n_formula"] == pytest.approx((norm_z() * 0.15 / 0.1) ** 2)
    assert len(m["rows"]) == 39 and all(a < b < c for a, b, c in (r["t10"] for r in m["rows"]))
    assert min(r["t10"][0] for r in m["rows"]) == pytest.approx(R["syn"]["lo"])
    assert m["upper_edge_t10"] == pytest.approx(R["syn"]["hi"])


def norm_z():
    from scipy.stats import norm
    return float(norm.ppf(0.975))


def test_timing_budget_is_consistent(R):
    b = R["syn"]["hw"]["budget"]
    assert b["margin_ms"] == pytest.approx(R["syn"]["margin_ms"])
    assert b["budget_ms"] == pytest.approx(0.1 * b["margin_ms"])
    assert b["fs_khz"] % 5 == 0 and b["fs_khz"] >= b["fs_min_khz"] > b["fs_khz"] - 5
    assert b["fs_khz"] == 15 and round(b["fs_min_khz"], 1) == 13.1
    assert b["counts"] == 2 ** b["bits"] >= b["counts_min"] > b["counts"] / 2
    assert b["enc_res_deg"] <= b["ang_max_deg"]
    parts = b["period_ms"] + b["trigger_ms"] + b["release_ms"] + b["enc_lat_ms"] + b["ang_time_ms"]
    assert b["u_t_ms"] == pytest.approx(parts) and b["ok"] and b["u_t_ms"] < b["budget_ms"]
    assert b["trigger_ms"] + b["release_ms"] + b["enc_lat_ms"] + b["spare_ms"] + b["period_max_ms"] + \
        b["ang_max_deg"] / b["w10_min"] * 1e3 == pytest.approx(b["budget_ms"])
    assert b["t10_min_pass_ms"] == pytest.approx(200 + b["u_t_ms"])
    assert b["drop200_max_pass_deg"] == pytest.approx(10 - b["ang_unc_deg"])
    assert b["worst_drop200_deg"] < b["drop200_max_pass_deg"]  # every model cell passes the per-trial rule


def test_rig_energies_match_the_ode_and_closed_forms(R):
    rg, A = R["syn"]["hw"]["rig"], nc.arm()
    assert rg["ff_stop_w"] == pytest.approx(rg["w_ff_energy"], rel=1e-9)  # ODE speed = energy speed
    assert rg["ff_stop_w"] == pytest.approx(229.18, abs=0.01)
    assert rg["ff_stop_ms"] == pytest.approx(130.30, abs=0.01)
    th = math.radians
    assert rg["E_ff_total"] == pytest.approx(A["G"] * math.sin(th(25)))
    assert rg["E_cap_total"] == pytest.approx(A["G"] * math.sin(th(25)) + nc.K_HI * nc.T_CAP_DEV * th(25))
    assert rg["E_cap_total"] > rg["E_ff_total"]
    assert rg["T_abs_ideal"] * th(10) == pytest.approx(rg["E_cap_total"])
    assert rg["T_abs_peak"] == pytest.approx(1.5 * rg["T_abs_ideal"])
    assert rg["F_abs_peak"] * 0.20 == pytest.approx(rg["T_abs_peak"])
    assert rg["E_up"] == pytest.approx((rg["kcap"] + A["hold_nom"]) * th(15) - A["G"] * math.sin(th(15)))
    assert rg["trip_Nm"] == pytest.approx(11.0) and rg["r_guard"] == pytest.approx(0.425)
    assert rg["sweep_mm"] == pytest.approx(400 * th(50), rel=1e-9)


def test_flip_values_reach_ten_degrees_in_an_independent_integrator(R):
    f, A = R["syn"]["hw"]["flip"], nc.arm()
    g = R["syn"]["gov_key"]
    rate, sto, td = g[1], g[2] / 1e3, nc.T_DET[g[0]] / 1e3
    base = _rk4(rate, sto, td, 0.0, t_end=0.2)
    assert base < 10
    ms = nc.arm(f["mass"])
    for kw in (dict(J=A["J"] * f["J"]), dict(J=ms["J"], G=ms["G"]), dict(T0=nc.T0_DEV * f["T0"]),
               dict(k=nc.K_HI * f["k"]), dict(hold=-f["offset"])):
        assert _rk4(rate, sto, td, kw.pop("hold", 0.0), t_end=0.2, **kw) == pytest.approx(10.0, abs=2e-3)
    assert f["lag_ms"] == pytest.approx(R["syn"]["margin_ms"])
    assert [round(x, 2) for x in R["syn"]["hw"]["flip_pct"].values()] == [-3.12, 1.29, -0.98, -1.22]


def test_calibration_limits_follow_the_slopes(R):
    c = R["syn"]["hw"]["cal"]
    assert c["share_ms"] == pytest.approx(0.125)
    assert c["net_torque_Nm"] == pytest.approx(0.125 / abs(c["slope_Th"]))
    # raising G by 1 % is the same as lowering the net pre-onset torque by 0.01 G
    assert abs(c["slope_G"]) == pytest.approx(abs(c["slope_Th"]) * nc.arm()["G"] / 100, rel=2e-3)
    assert (round(c["net_torque_Nm"], 4), round(c["G_pct"], 3), round(c["J_pct"], 2)) == (0.0059, 0.057, 0.21)


def test_plan_has_the_new_items():
    plan = nc.PLAN.read_text(encoding="utf-8")
    for need in ("independent of the device e-stop", "Interlocked guard", "Upper runaway stop",
                 "must be verified", "raw encoder counts", "stop contact force and stroke",
                 "calibration records", "measured skew", "1,170 drops"):
        assert need in plan, need
