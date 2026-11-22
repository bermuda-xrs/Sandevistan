"""scripts/sd_stress.py: scaled variance components, gate equality, subject-count fit. Synthetic data only."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from emg_onset_lead import PLANNED_SHAPE, SampleSize, SyntheticSource

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "sd_stress.py"
_spec = importlib.util.spec_from_file_location("sd_stress", _PATH)
ss = importlib.util.module_from_spec(_spec)
sys.modules["sd_stress"] = ss
_spec.loader.exec_module(ss)
sim = ss.sim


def test_protected_values_are_read_not_set():
    assert sim.GATE_MS == pytest.approx(68.67)
    assert PLANNED_SHAPE == (16, 3) and ss.N_SUB == 16
    assert ss.LEADS == (74.67, 68.67)
    assert (sim.SD_SUBJECT, sim.SD_SESSION, sim.SD_TRIAL) == ss.DESIGN == (8, 6, 10)


def test_plan_is_fixed():
    assert (ss.N_STUDIES, ss.STUDY0) == (1000, 1000)  # row12_power_sim.py used ids 0-249
    assert ss.CONFIGS[0] == (12, 4, 10)
    assert [c for c in ss.CONFIGS[1:]] == [(s, 6, 10) for s in (8, 10, 12, 14, 16)]
    assert ss.RISES == (10, 20, 30, 40, 80)


@pytest.mark.parametrize("cfg", [(12, 4, 10), (16, 6, 10)])
def test_scaled_components_match_config(cfg):
    ds = [ss.scaled_design(s, cfg, 48) for s in range(64)]  # 3072 subjects
    so, se, eps = (np.concatenate([d[k] for d in ds]) for k in ("so", "se", "eps"))
    # about 4 SE each
    assert so.std(ddof=1) == pytest.approx(cfg[0], rel=0.055)
    assert se.std(ddof=1) == pytest.approx(cfg[1], rel=0.035)
    assert eps.std(ddof=1) == pytest.approx(cfg[2], rel=0.012)
    assert max(abs(np.corrcoef(so, se[:, j])[0, 1]) for j in range(3)) < 0.075
    off = so[:, None, None] + se[:, :, None] + eps
    assert off.mean((1, 2)).std(ddof=1) == pytest.approx(SampleSize.design_subject_sd(*cfg, 3, 10), rel=0.055)


def test_merged_synthetic_source_has_the_stress_spread_too():
    src = SyntheticSource(n_subjects=400, n_sessions=3, n_trials=2, subject_sd_ms=12, session_sd_ms=4, seed=5,
                          duration_s=0.01, fs_emg=200.0)
    m = np.array([t.meta["true_lead_ms"] for t in src]).reshape(400, 6).mean(1)
    assert m.std(ddof=1) == pytest.approx(SampleSize.design_subject_sd(12, 4, 10, 3, 2), rel=0.15)


def test_design_config_reproduces_row12_signals():
    """At 8/6/10 the onsets give exactly row12_power_sim's leads, and the configs share EMG onsets."""
    r = ss.simulate_subject((3, 1))
    d = sim.draw_design(3, 2)
    off = sim.offsets(d)[1]
    ci = ss.CONFIGS.index(ss.DESIGN)
    np.testing.assert_array_equal(r["off"][ci], off)
    for j, k in ((0, 0), (2, 9)):
        for ri, rise in enumerate(ss.RISES[:2]):
            e = sim.detect_emg(sim.emg_signal(3, 1, j, k, rise, d["ton"][1, j, k], d["gain"][1, j]))
            for li, L in enumerate(ss.LEADS):
                im = sim.detect_imu(sim.angle_signal(3, 1, j, k, d["ton"][1, j, k], L + off[j, k]))
                got = (r["imu"][ci, li, j, k] - r["emg"][ri, j, k]) * 1000.0
                np.testing.assert_equal(got, sim.combine(e, im)[0])


def test_gate_blocks_equals_summarize():
    rng = np.random.default_rng(1)
    seen = set()
    for q in range(40):
        lead = (69 + rng.uniform(0, 10) + rng.normal(0, 8, (64, 1, 1)) + rng.normal(0, 6, (64, 3, 1))
                + rng.normal(0, 10, (64, 3, 10)))
        lead[rng.random(lead.shape) < (0.12 if q % 4 == 0 else 0.01)] = np.nan  # miss cap on both sides
        if q % 5 == 0:
            lead[3] = np.nan  # a fully missed subject
            lead[20, 1] = np.nan  # a fully missed session
        m, miss = ss.subject_stats(lead)
        for n in (16, 32):
            got = ss.gate_blocks(m, miss, n)
            want = [ss.gate_summarize(lead[b * n:(b + 1) * n]) for b in range(64 // n)]
            assert got.tolist() == want
            seen |= set(want)
    assert seen == {True, False}


def test_miss_cap_fails_a_clear_pass():
    lead = np.full((16, 3, 10), 90.0) + np.random.default_rng(0).normal(0, 1, (16, 3, 10))
    m, q = ss.subject_stats(lead)
    assert ss.gate_blocks(m, q, 16).tolist() == [True]
    lead[:, :, :2] = np.nan  # 20% missed
    m, q = ss.subject_stats(lead)
    assert ss.gate_blocks(m, q, 16).tolist() == [False]


def test_fit_matches_row12_fit():
    ks, ms = [209, 223, 235, 239, 242, 247, 248, 250, 250], [250] * 9  # docs/row12_power_sim.md, 10 ms
    x = ss.crossing(ss.fit_probit(sim.N_GRID, ks, ms))
    assert x == pytest.approx(sim.probit_sqrt_fit(sim.N_GRID, ks, 250)[0], abs=0.01)
    assert x == pytest.approx(14.7, abs=0.06)


def test_fit_edge_cases():
    n = list(ss.N_GRID)
    assert ss.crossing(ss.fit_probit(n, [0] * len(n), [50] * len(n))) == float("inf")
    assert ss.crossing(ss.fit_probit(n, [20 - i for i in range(len(n))], [50] * len(n))) == float("inf")


def test_subjects_needed_on_normal_subject_means(monkeypatch):
    """Known normal subject means: the crossing lands on the noncentral-t count."""
    monkeypatch.setattr(ss, "N_BOOT", 40)
    sd = 12.35
    m = sim.GATE_MS + 6.0 + sd * np.random.default_rng(2).normal(size=16000)
    d = ss.subjects_needed(m, np.zeros(16000, int))
    want = ss.nct_n(6.0, sd)
    assert want == SampleSize.required_subjects(sd**2, 0.0, 6.0, 0.05, 0.8, 1)
    assert d["n80"] == pytest.approx(want, rel=0.10)
    assert d["lo"] < d["n80"] < d["hi"]
    k = d["check"]
    assert stats.binomtest(k["k"], k["m"], 0.8).pvalue > 1e-3
    assert not ss.subjects_needed(m - 8.0, np.zeros(16000, int))["reachable"]
