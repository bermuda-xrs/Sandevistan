"""Seeding and variance components of scripts/row12_power_sim.py. Synthetic data only.

The old row-12 generator seeded sessions with ``default_rng([seed, i, j])``. SeedSequence pads a seed
list that converts to fewer 32-bit words (not list items) than its pool, 4 by default, with zeros. So
session 0 reused the subject's stream (session-0 offset = 0.75 x subject offset, subject-mean SD
10.55 ms in theory instead of 8.91).
``test_legacy_scheme_fails_the_check`` keeps that scheme as a negative control: the same check must
reject it.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from emg_onset_lead.harness import Harness
from emg_onset_lead.sources import Trial

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "row12_power_sim.py"
_spec = importlib.util.spec_from_file_location("row12_power_sim", _PATH)
sim = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sim)

N_STUDIES, N_SUB = 64, 48  # 3072 subjects, 9216 sessions, 92160 trials


def violations(vc: dict) -> list[str]:
    """Design 8/6/10 ms checks; tolerances are about 4 SE at 3072 subjects."""
    bad = []
    for key, want, tol in (("sd_subject", 8.0, 0.45), ("sd_session", 6.0, 0.2), ("sd_trial", 10.0, 0.12),
                           ("sd_subject_mean", vc["sd_subject_mean_design"], 0.45)):
        if abs(vc[key] - want) > tol:
            bad.append(f"{key} {vc[key]:.3f} vs {want:.3f}")
    for j, c in enumerate(vc["corr_subject_session"]):
        if abs(c) > 0.075:
            bad.append(f"corr(subject, session {j}) = {c:+.3f}")
    for c in vc["corr_session_pairs"]:
        if abs(c) > 0.075:
            bad.append(f"corr(session, session) = {c:+.3f}")
    return bad


def _pooled(draw):
    ds = [draw(s) for s in range(N_STUDIES)]
    return sim.variance_check(*(np.concatenate([d[k] for d in ds]) for k in ("so", "se", "eps")))


@pytest.fixture(scope="module")
def new_vc():
    return _pooled(lambda s: sim.draw_design(s, N_SUB))


def _legacy_design(seed: int) -> dict:
    """The old ps_study.py / fc4_study.py scheme (bug kept on purpose)."""
    so = np.empty(N_SUB)
    se = np.empty((N_SUB, 3))
    eps = np.empty((N_SUB, 3, 10))
    for i in range(N_SUB):
        so[i] = np.random.default_rng([seed, i]).normal(0, 8)
        for j in range(3):
            se[i, j] = np.random.default_rng([seed, i, j]).normal(0, 6)
            for k in range(10):
                eps[i, j, k] = np.random.default_rng([seed, i, j, k + 1]).normal(0, 10)
    return dict(so=so, se=se, eps=eps)


def test_design_sd_is_8_91():
    r = np.random.default_rng(0)
    vc = sim.variance_check(r.normal(size=4), r.normal(size=(4, 3)), r.normal(size=(4, 3, 10)))
    assert vc["sd_subject_mean_design"] == pytest.approx(8.907, abs=1e-3)


def test_variance_components_match_design(new_vc):
    assert violations(new_vc) == [], new_vc


def test_no_session_offset_correlated_with_subject(new_vc):
    assert max(abs(c) for c in new_vc["corr_subject_session"]) < 0.075, new_vc["corr_subject_session"]


def test_legacy_scheme_fails_the_check():
    vc = _pooled(lambda s: _legacy_design(1000 + s))
    bad = violations(vc)
    assert vc["corr_subject_session"][0] > 0.99  # session 0 reuses the subject stream
    assert any("corr(subject, session 0)" in b for b in bad)
    assert any(b.startswith("sd_subject_mean") for b in bad)  # 10.55 ms in theory, not 8.91


def test_seed_nodes_are_distinct():
    keys = [(0,), (0, 0), (0, 0, 0), (0, 0, 1), (0, 0, 1, 0), (0, 0, 1, 1), (0, 0, 1, 1, 0),
            (0, 0, 1, 1, 1), (0, 0, 1, 1, 2), (0, 1), (0, 1, 0), (1,), (1, 0)]
    states = {tuple(sim.node(*k).generate_state(4)) for k in keys}
    assert len(states) == len(keys)


def test_draws_are_nested_across_study_sizes():
    small, big = sim.draw_design(7, 16), sim.draw_design(7, 48)
    for k in ("so", "se", "eps", "ton", "gain"):
        np.testing.assert_array_equal(small[k], big[k][:16])


def test_combine_matches_harness_run():
    d = sim.draw_design(3, 1)
    off = sim.offsets(d)
    trials, mine = [], []
    for k in range(3):
        ton, g = d["ton"][0, 0, k], d["gain"][0, 0]
        emg = sim.emg_signal(3, 0, 0, k, 20, ton, g)
        ang = sim.angle_signal(3, 0, 0, k, ton, 74.67 + off[0, 0, k])
        trials.append(Trial("S0", "S0-D0", f"T{k}", sim.FS_EMG, emg, sim.FS_IMU, ang, {}))
        mine.append(sim.combine(sim.detect_emg(emg), sim.detect_imu(ang)))
    for r, (lead, delay) in zip(Harness.run(trials), mine):
        np.testing.assert_equal([r.lead_ms, r.emg_confirmation_delay_ms], [lead, delay])
