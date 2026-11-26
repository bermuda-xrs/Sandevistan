"""ML ramp-model onset research code (scripts/). SYNTHETIC data only."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import ml_ramp_row12_sim as sim  # noqa: E402
from ml_ramp_onset import MlRamp, ar_whitener  # noqa: E402


@pytest.mark.parametrize("rise", [0.02, 0.08])
def test_fit_recovers_a_linear_ramp(rise):
    fs, n, rng = 2000.0, 1000, np.random.default_rng(1)
    t = np.arange(n) / fs
    z = rng.normal(size=(200, n)) * (1 + 4 * np.clip((t - 0.35) / rise, 0, 1))
    t0, rhat, gain, _ = MlRamp(fs, n).fit(z)
    assert abs(np.mean(t0 + 0.1 * rhat) - (0.35 + 0.1 * rise)) < 1.5e-3
    assert abs(gain.mean() - 5) < 0.25


def test_whitener_leaves_unit_variance_and_less_correlation_than_the_bandpassed_baseline():
    b = sim.stream_filter(sim.BP, np.random.default_rng(2).normal(size=(1, sim.N)))[0]
    A, sd = ar_whitener(b[sim.SETTLE_N:sim.BASE_N])
    w = np.convolve(b, A)[sim.BASE_N:sim.N] / sd

    def lag_corr(v):
        return np.array([np.mean(v[:-k] * v[k:]) for k in range(1, 6)]) / np.mean(v * v)

    assert abs(w.std() - 1) < 0.1
    assert np.abs(lag_corr(w)).max() < 0.15
    assert np.abs(lag_corr(b[sim.BASE_N:])).max() > 0.3


@pytest.mark.parametrize("shape", sim.SHAPES)
@pytest.mark.parametrize("rise", [0.01, 0.08])
def test_true_onset_is_one_plateau_fraction_for_every_shape(shape, rise):
    t0 = sim.ramp_start(shape, rise, 2.0)
    v = float(sim.burst(shape, np.array([2.0]), t0, rise)[0])
    assert abs(v - sim.P) < 1e-9


def test_batch_trigger_equals_the_repo_detector():
    x = sim.emg_rows(np.random.default_rng(3).normal(size=sim.N), 2.0, 0.8)
    b = sim.stream_filter(sim.BP, x)
    f = sim.stream_filter(sim.LP, np.abs(b))
    for r, tr in enumerate(sim.triggers(x, b, f)):
        assert tr == sim.slow_trigger(x[r])


def test_seed_tree_is_fixed_and_session_offset_is_not_tied_to_the_subject_offset():
    a, b = sim.draw_design(0, 2), sim.draw_design(0, 2)
    assert all(np.array_equal(a[k], b[k]) for k in a)
    # the old [seed, i] / [seed, i, 0] collision made session 0 exactly 0.75 x the subject offset
    so = np.array([sim.rng_at(0, i, 0).normal() for i in range(200)])
    se = np.array([sim.rng_at(0, i, 1, 0).normal() for i in range(200)])
    assert abs(np.corrcoef(so, se)[0, 1]) < 0.3


def test_gate_values_are_the_protected_ones():
    assert sim.GATE_MS == 68.67 and sim.N_SUB == 16 and sim.N_SESSIONS == 3
