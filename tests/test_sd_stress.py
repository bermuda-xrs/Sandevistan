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


# ------------------------------------------------ round 2: tests aimed at the mutants FC found alive
def _fake_data(n_sub=64, shift=-3.0, noise=0.0, seed=0, miss_frac=0.0):
    """Small stand-in for the raw chunks: EMG onset 2 s; IMU onset = EMG + (true lead + shift + noise)."""
    rng = np.random.default_rng(seed)
    nc, nr = len(ss.CONFIGS), len(ss.RISES)
    off = rng.normal(0, 8, (n_sub, nc, 3, 10)) + rng.normal(0, 6, (n_sub, nc, 3, 1))
    emg = np.full((n_sub, nr, 3, 10), 2.0)
    imu = np.empty((n_sub, nc, 2, 3, 10))
    for li, lead in enumerate(ss.LEADS):
        imu[:, :, li] = 2.0 + (lead + off + shift + rng.normal(0, noise, off.shape)) / 1000.0
    if miss_frac:
        imu[rng.random(imu.shape) < miss_frac] = np.nan
    return dict(emg=emg, imu=imu, off=off)


def test_simulate_subject_uses_each_configs_own_offsets():
    r = ss.simulate_subject((3, 1))
    d = [ss.scaled_design(3, c, 2) for c in ss.CONFIGS]
    for ci in range(len(ss.CONFIGS)):
        np.testing.assert_array_equal(r["off"][ci], sim.offsets(d[ci])[1])
    assert not np.array_equal(r["off"][0], r["off"][1])  # 12/4/10 is not the design draw
    ton, gain = d[0]["ton"][1], d[0]["gain"][1]
    for ci in range(len(ss.CONFIGS)):  # the IMU onsets come from that config's offsets, not the design's
        for li, lead in enumerate(ss.LEADS):
            for j, k in ((0, 0), (2, 9)):
                a = sim.angle_signal(3, 1, j, k, ton[j, k], lead + r["off"][ci][j, k])
                assert r["imu"][ci, li, j, k] == ss._onset(sim.detect_imu(a))
    for ri, rise in enumerate(ss.RISES):
        x = sim.emg_signal(3, 1, 1, 4, rise, ton[1, 4], gain[1])
        assert r["emg"][ri, 1, 4] == ss._onset(sim.detect_emg(x))


def test_scaled_design_rescales_all_three_levels():
    base = sim.draw_design(5, 4)
    d = ss.scaled_design(5, (12, 3, 20), 4)
    np.testing.assert_allclose(d["so"], base["so"] * 1.5)
    np.testing.assert_allclose(d["se"], base["se"] * 0.5)
    np.testing.assert_allclose(d["eps"], base["eps"] * 2.0)
    np.testing.assert_array_equal(d["ton"], base["ton"])
    np.testing.assert_array_equal(d["gain"], base["gain"])


def test_onset_drops_unusable_trials():
    from types import SimpleNamespace as NS
    assert ss._onset(NS(baseline_rest_ok=True, onset_s=1.25)) == 1.25
    assert np.isnan(ss._onset(NS(baseline_rest_ok=False, onset_s=1.25)))
    assert np.isnan(ss._onset(NS(baseline_rest_ok=True, onset_s=None)))


def test_subject_mean_is_mean_of_session_means():
    lead = np.full((1, 3, 10), 70.0)
    lead[0, 0, 2:] = np.nan  # session 0 keeps 2 trials, so its mean weighs as much as the full ones
    lead[0, 0, :2] = [60.0, 62.0]
    lead[0, 1] = 70.0
    lead[0, 2] = 80.0
    m, miss = ss.subject_stats(lead)
    assert m[0] == pytest.approx((61.0 + 70.0 + 80.0) / 3)  # all-trial mean would be 73.5
    assert miss[0] == 8


def test_gate_needs_the_bound_strictly_above_the_gate():
    at = np.full(16, sim.GATE_MS)  # zero spread: the lower bound equals the mean exactly
    assert ss.gate_blocks(at, np.zeros(16, int), 16).tolist() == [False]
    assert ss.gate_blocks(at + 1e-6, np.zeros(16, int), 16).tolist() == [True]


def test_gate_ignores_nan_subject_and_needs_two_subjects():
    m = np.array([75.0, 80.0, np.nan, 90.0])
    assert ss.gate_blocks(m, np.zeros(4, int), 4).tolist() == [True]
    assert ss.gate_blocks(np.array([90.0, np.nan]), np.zeros(2, int), 2).tolist() == [False]


def test_score_cell_bias_and_removal_signs(monkeypatch):
    d = _fake_data(shift=-3.0, noise=2.0, seed=1)
    monkeypatch.setattr(ss, "_D", d)
    monkeypatch.setattr(ss, "N_BOOT", 5)
    for ci, ri, li in ((0, 0, 0), (1, 3, 1)):
        out = ss.score_cell((ci, ri, li))
        true = ss.LEADS[li] + d["off"][:, ci]
        assert out["bias"] == pytest.approx(-3.0, abs=0.1)  # detector reads short: negative
        m_true = ss.subject_stats(true)[0]
        # bias removed -> the true mean
        assert out["unb"]["mean"] == pytest.approx(float(np.nanmean(m_true)), abs=0.1)
        assert out["raw"]["mean"] == pytest.approx(float(np.nanmean(m_true)) + out["bias"], abs=1e-9)
        assert out["raw"]["vs_summarize"] == out["unb"]["vs_summarize"] == 0
    out = ss.score_cell((0, 0, 0))
    assert out["lead"] == 74.67 and out["n_trials"] == 64 * 30 and out["n_miss"] == 0
    # the effect size is the mean of the subject means minus the gate, with no extra offset
    assert out["raw"]["need"]["effect"] == pytest.approx(out["raw"]["mean"] - sim.GATE_MS, abs=1e-12)
    assert out["unb"]["need"]["effect"] == pytest.approx(out["unb"]["mean"] - sim.GATE_MS, abs=1e-12)


def test_score_cell_counts_missed_and_early_late_trials(monkeypatch):
    d = _fake_data(shift=0.0, noise=0.0, seed=2)
    d["imu"][0, 0, 0, 0, 0] += 0.5  # +500 ms late
    d["imu"][1, 0, 0, 0, 1] -= 0.5  # -500 ms early
    d["imu"][2, 0, 0, 0, 2] = np.nan
    monkeypatch.setattr(ss, "_D", d)
    monkeypatch.setattr(ss, "N_BOOT", 5)
    out = ss.score_cell((0, 0, 0))
    assert (out["n_late"], out["n_early"], out["n_miss"]) == (1, 1, 1)


def test_score_true_uses_true_leads(monkeypatch):
    d = _fake_data(seed=3)
    monkeypatch.setattr(ss, "_D", d)
    monkeypatch.setattr(ss, "N_BOOT", 5)
    p, q = ss.score_true((1, 0)), ss.score_true((1, 1))
    mp = ss.subject_stats(74.67 + d["off"][:, 1])[0]
    assert p["mean"] == pytest.approx(mp.mean()) and q["mean"] == pytest.approx(mp.mean() - 6.0)
    assert p["need"]["effect"] == pytest.approx(mp.mean() - sim.GATE_MS, abs=1e-12)
    assert "need" not in q


def test_subjects_needed_bootstrap_follows_its_documented_procedure(monkeypatch):
    """95% CI = 2.5/97.5 percentiles of N_BOOT resamples of subjects, as many as there are subjects."""
    monkeypatch.setattr(ss, "N_BOOT", 100)  # enough resamples that 2.5/97.5 and 5/95 pick different ones
    m = sim.GATE_MS + 6.0 + 12.0 * np.random.default_rng(4).normal(size=2000)
    miss = np.zeros(2000, int)
    got = ss.subjects_needed(m, miss)
    rng = np.random.default_rng(ss.BOOT_SEED)
    boot = []
    for _ in range(100):
        pick = rng.integers(0, m.size, m.size)
        boot.append(ss.crossing(ss.fit_probit(ss.N_GRID, *ss.sweep(m[pick], miss[pick]))))
    assert (got["lo"], got["hi"]) == pytest.approx(tuple(np.percentile(boot, [2.5, 97.5], method="nearest")))
    assert got["effect"] == pytest.approx(m.mean() - sim.GATE_MS)
    assert got["sd"] == pytest.approx(m.std(ddof=1))
    assert got["nct"] == ss.nct_n(got["effect"], got["sd"])
    assert ss.nct_n(-1.0, 10.0) is None


def test_realised_components_recover_the_design(monkeypatch):
    ds = [sim.draw_design(s, 48) for s in range(64)]
    off = np.concatenate([sim.offsets(d) for d in ds])  # (3072, 3, 10)
    monkeypatch.setattr(ss, "_D", dict(off=off[:, None]))
    monkeypatch.setattr(ss, "CONFIGS", ((8, 6, 10),))
    v = ss.realised(0)
    assert v["subject"] == pytest.approx(8.0, abs=0.45)
    assert v["session"] == pytest.approx(6.0, abs=0.25)  # without the trial-variance correction it reads 6.8
    assert v["trial"] == pytest.approx(10.0, abs=0.12)
    assert v["design"] == pytest.approx(8.907, abs=1e-3) and v["lo"] < v["sd_mean"] < v["hi"]


def _write_chunk(path, ids):
    n = len(ids)
    np.savez(path, studies=np.array(ids), head="abc", numpy="x", emg=np.zeros((n, 16, 1)),
             imu=np.zeros((n, 16, 1)), off=np.zeros((n, 16, 1)))


def test_load_requires_every_study_exactly_once(tmp_path):
    _write_chunk(tmp_path / "raw_001000.npz", [1000, 1001])
    _write_chunk(tmp_path / "raw_001002.npz", [1002, 1003])
    out = ss.load(tmp_path, 4, 1000)
    assert out["emg"].shape == (64, 1) and out["heads"] == ["abc"]
    with pytest.raises(AssertionError):
        ss.load(tmp_path, 5, 1000)  # a chunk is missing
    with pytest.raises(AssertionError):
        ss.load(tmp_path, 4, 999)  # wrong first id
    _write_chunk(tmp_path / "raw_001004.npz", [1003, 1004])
    with pytest.raises(AssertionError):
        ss.load(tmp_path, 6, 1000)  # study 1003 twice


def test_chunking_covers_the_plan_without_overlap():
    assert ss.CHUNK == 25 and ss.N_STUDIES % ss.CHUNK == 0
    ids = [i for lo in range(ss.STUDY0, ss.STUDY0 + ss.N_STUDIES, ss.CHUNK) for i in range(lo, lo + ss.CHUNK)]
    assert ids == list(range(1000, 2000))
    assert ss.N_GRID[0] == 8 and ss.N_GRID[-1] == 512 and 16 in ss.N_GRID
    assert list(ss.N_GRID) == sorted(set(ss.N_GRID))


def test_row12_reference_counts_match_the_row12_doc():
    text = (Path(__file__).resolve().parents[1] / "docs" / "row12_power_sim.md").read_text()
    block = text.split("### Baseline five-row table")[1].split("###")[0]
    rows = {int(m[0]): (int(m[1]), int(m[2])) for m in
            __import__("re").findall(r"\| (\d+) ms \| (\d+)/250 = [\d.]+ \([^)]*\) \| (\d+)/250", block)}
    assert {r: v[0] for r, v in rows.items()} == ss.ROW12_DOC["power"]
    assert {r: v[1] for r, v in rows.items()} == ss.ROW12_DOC["false_pass"]


def test_committed_tables_are_what_report_prints_from_committed_scores(tmp_path):
    """report is deterministic and the committed tables come from the committed scored.json."""
    root = Path(__file__).resolve().parents[1] / "results" / "sd_stress"
    out = tmp_path / "t.md"
    ss.main(["report", str(root / "scored.json"), "--out", str(out)])
    assert out.read_text() == (root / "tables.md").read_text()
    D = __import__("json").loads((root / "scored.json").read_text())
    assert D["n_grid"] == list(ss.N_GRID) and D["n_studies"] == ss.N_STUDIES and D["first"] == ss.STUDY0
    assert [tuple(c) for c in D["configs"]] == list(ss.CONFIGS) and D["gate_ms"] == sim.GATE_MS


def test_sensitivity_matches_the_main_gate_at_the_protected_margin(tmp_path, monkeypatch):
    d = _fake_data(n_sub=32, shift=-3.0, noise=2.0, seed=5)
    raw = tmp_path / "raw"
    raw.mkdir()
    by_study = {k: v.reshape((2, 16) + v.shape[1:]) for k, v in d.items()}
    np.savez(raw / "raw_000000.npz", studies=np.array([0, 1]), head="h", numpy="x", **by_study)
    out = tmp_path / "s.md"
    ss.main(["sensitivity", "--raw", str(raw), "--out", str(out), "--studies", "2", "--first", "0"])
    rows = [ln.split("|") for ln in out.read_text().splitlines() if ln.startswith("| 12/4/10 | 10 ms rise")]
    lead = (d["imu"][:, 0, 0] - d["emg"][:, 0]) * 1000.0
    for margin, cell in zip(ss.SENS_MARGINS, rows[0][3:-1]):
        m, q = ss.subject_stats(lead + (margin - 6.0))
        p = ss.gate_blocks(m, q, ss.N_SUB)
        assert cell.strip().startswith(f"{int(p.sum())}/{p.size} =")
