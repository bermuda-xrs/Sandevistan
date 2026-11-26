"""Checks of the research detectors against brute-force references. SYNTHETIC DATA ONLY.

Run: PYTHONPATH=src python -m pytest research/causal_detectors/test_detectors.py
"""

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import integrate, stats

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "src")]

import detectors as D  # noqa: E402

from emg_onset_lead import EmgOnsetDetector, SyntheticSource  # noqa: E402

FS = 2000.0
N0 = 2000


def _burst(n=3, seed=0, snr=3.0, dur_s=1.6, at_s=1.3):
    rng = np.random.default_rng(seed)
    t = np.arange(round(dur_s * FS)) / FS
    return (1 + (snr - 1) * np.clip((t - at_s) / 0.04, 0, 1)) * rng.normal(0, 1, (n, t.size))


@pytest.mark.parametrize("hold", [1, 2, 3, 10, 11, 60])
def test_trailing_min(hold):
    d = np.random.default_rng(0).normal(size=(2, 200))
    out = D.trailing_min(d, hold, start=50)
    ref = np.full(d.shape, -np.inf)
    for i in range(50 + hold - 1, 200):
        ref[:, i] = d[:, i - hold + 1:i + 1].min(axis=1)
    np.testing.assert_array_equal(out, ref)


@pytest.mark.parametrize("rho", [2.0, 3.0])
def test_cusum_is_pages_recursion(rho):
    x = _burst()
    s = D.stat_cusum(x, FS, rho)
    z = D.scaled_bp(x, FS)
    g = np.zeros(3)
    for i in range(N0, x.shape[1]):
        g = np.maximum(0.0, g + 0.5 * (1 - 1 / rho**2) * z[:, i] ** 2 - np.log(rho))
        np.testing.assert_allclose(s[:, i], g, atol=1e-8)


def test_glr_brute_force():
    x = _burst(n=1)
    s = D.stat_glr(x, FS, window_s=0.05)
    z2 = D.scaled_bp(x, FS)[0] ** 2
    for i in (N0 + 2, N0 + 3, N0 + 40, N0 + 600, N0 + 700, N0 + 1000):
        best = 0.0
        for L in range(4, 101, 4):
            if i - L + 1 < N0:
                break
            v = z2[i - L + 1:i + 1].mean()
            if v > 1:
                best = max(best, 0.5 * L * (v - 1 - np.log(v)))
        assert s[0, i] == pytest.approx(best, abs=1e-8)


def _pred(e, m, blocks, a0, b0):
    """Density of a block energy e given earlier block energies, by numerical integration over v."""
    def joint(v, es):
        out = stats.invgamma.pdf(v, a0, scale=b0)
        for q in es:
            out = out * stats.chi2.pdf(q / v, m) / v
        return out
    lim = dict(limit=400, epsabs=0, epsrel=1e-9)
    f = lambda es: sum(  # noqa: E731
        integrate.quad(lambda u: joint(np.exp(u), es) * np.exp(u), lo, lo + 4, **lim)[0]
        for lo in range(-8, 24, 4))
    return f(blocks + [e]) / f(blocks)


def test_bocpd_matches_numerical_reference():
    """Reference recursion: predictives by quadrature, no conjugate formulas, no log tricks."""
    x = _burst(n=1, dur_s=1.05, at_s=1.01, snr=4.0)[:, : N0 + 4 * 12]
    m, R, a0, b0, hz, vmin = 4, 5, 1.5, 2.0, 1.0 / 1000, 1.5
    s = D.stat_bocpd(x, FS, r_max_s=R * 0.002)
    z = D.scaled_bp(x, FS)[0]
    e = (z[N0: N0 + 48] ** 2).reshape(12, 4).sum(axis=1)
    rest, runs = 1.0, {}  # runs: r -> (mass, blocks)
    for t in range(12):
        p0 = stats.chi2.pdf(e[t], m)
        new = {0: (hz * _pred(e[t], m, [], a0, b0), [e[t]])}
        rest_new = rest * (1 - hz) * p0
        for r, (p, bl) in runs.items():
            if r == R - 1:
                rest_new += p * (1 - hz) * p0
            else:
                new[r + 1] = (p * (1 - hz) * _pred(e[t], m, bl, a0, b0), bl + [e[t]])
        tot = rest_new + sum(p for p, _ in new.values())
        rest, runs = rest_new / tot, {r: (p / tot, bl) for r, (p, bl) in new.items()}
        sel = sum(p for r, (p, bl) in runs.items()
                  if r >= 1 and (b0 + sum(bl) / 2) / (a0 + len(bl) * m / 2 - 1) > vmin)
        got = s[0, N0 + 3 + t * m]
        if sel > 0:
            assert got == pytest.approx(np.log(sel / rest), rel=1e-6, abs=1e-6)
        else:
            assert got < -100


def test_tkeo_known_tone():
    """TKEO of A cos(w n) is A^2 sin^2(w) (Kaiser 1990); here delayed one sample."""
    n = np.arange(400)
    y = 2.0 * np.cos(0.3 * n)[None, :]
    psi = y[:, 1:-1] ** 2 - y[:, :-2] * y[:, 2:]
    np.testing.assert_allclose(psi, 4.0 * np.sin(0.3) ** 2, atol=1e-9)


@pytest.mark.parametrize("name", list(D.DETECTORS))
def test_causal(name):
    x = _burst(n=4, seed=3)
    for cut in (N0 + 103, N0 + 601, N0 + 1002):
        err, same = D.check_causal(name, x, FS, cut)
        assert same and err == 0.0


@pytest.mark.parametrize("snr", [3.0, 5.0, 15.0])
def test_env8_h30_is_the_production_trigger(snr):
    tr = list(SyntheticSource(n_subjects=1, emg_snr=snr, seed=99))
    x = np.stack([t.emg for t in tr])
    a = D.first_alarm(D.statistic("env8_h30", x, FS), 6.0, N0)
    checked = 0
    for i, t in enumerate(tr):
        r = EmgOnsetDetector.run(t.emg, FS, chunk=100)
        if r.trigger_s is not None:
            checked += 1
            assert round(r.trigger_s * FS) == a[i]
    assert checked == len(tr)
