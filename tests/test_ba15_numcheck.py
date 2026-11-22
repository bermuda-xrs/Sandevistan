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
    for name in ("feas", "feas_icc", "pru"):
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


def _rk4_crossing(rate, sto, td, hold, dt=1e-5):
    """Independent fixed-step RK4 in the absolute elbow angle (90 deg at start)."""
    A = nc.arm()
    brk = round(sto / dt)

    def acc(th, t, off):
        T = 0.0 if off else (nc.T0_DEV if t < td else max(nc.T0_DEV - rate * (t - td), 0.0))
        return (nc.K_HI * T + hold - A["G"] * math.sin(th)) / A["J"]

    th, w, prev = math.pi / 2, 0.0, 0.0
    for i in range(30000):
        t, off = i * dt, i >= brk
        k1 = (w, acc(th, t, off))
        k2 = (w + dt / 2 * k1[1], acc(th + dt / 2 * k1[0], t + dt / 2, off))
        k3 = (w + dt / 2 * k2[1], acc(th + dt / 2 * k2[0], t + dt / 2, off))
        k4 = (w + dt * k3[1], acc(th + dt * k3[0], t + dt, off))
        th += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        w += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        d = 90 - math.degrees(th)
        if d >= nc.DROP_LIM:
            return 1e3 * (t + dt * (nc.DROP_LIM - prev) / (d - prev))
        prev = d
    raise AssertionError


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
