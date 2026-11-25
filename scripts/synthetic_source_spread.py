"""Numbers quoted for ``SyntheticSource`` in README, PLAN.md and test comments. SYNTHETIC DATA ONLY.

No recorded EMG or motion is read. Prints, at this checkout and for the old seeding (kept here as
``LegacySource`` for comparison only):
  1. true-lead spread over ``--pool-seeds`` seeds x 16 subjects x 3 sessions x 10 trials: the SD of
     subject means, its design value (``SampleSize.design_subject_sd``), a chi-square 95% CI, and
     corr(subject offset, session j offset). No detector runs; signals are cut short.
  2. the same subject-mean SD for single seeds 0..5 (what one ``run --synthetic`` shows), and for the
     detected lead at this checkout.
Then, at this checkout only, with the detectors at default settings:
  3. the lead error of ``tests/test_harness.py::test_recovers_known_lead`` (8x3x5, seeds 0-5).
  4. the EMG confirmation delay of ``test_triggered_mode_adds_confirmation_delay`` (4x2x5, seed 0)
     and of the planned 16x3x10 study, seed 0.
The in-family bias and guard band come from ``scripts/bias_guard_band.py``.

Usage (repo root):  python scripts/synthetic_source_spread.py [--pool-seeds 40] [--out out.md]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from emg_onset_lead import (  # noqa: E402
    PLANNED_SHAPE,
    Harness,
    LatencyConfig,
    LeadStats,
    SampleSize,
    SyntheticSource,
)

N_TRIALS = 10
SHORT = dict(duration_s=0.01, fs_emg=200.0, fs_imu=150.0)  # signals are not read in parts 1-2


class LegacySource(SyntheticSource):
    """Seeding before the fix (session 0 reuses the subject stream). Comparison only."""

    def __iter__(self):
        for s in range(self.n_subjects):
            subj_off = np.random.default_rng([self.seed, s]).normal(0, self.subject_sd_ms)
            for ss in range(self.n_sessions):
                rng_s = np.random.default_rng([self.seed, s, ss])
                sess_off = rng_s.normal(0, self.session_sd_ms)
                gain = rng_s.uniform(0.7, 1.3)
                for t in range(self.n_trials):
                    rng = np.random.default_rng([self.seed, s, ss, t + 1])
                    yield self._trial(rng, s, ss, t, subj_off + sess_off, gain)


def offsets(cls, seed: int, **sd) -> np.ndarray:
    """True lead minus the 80 ms mean, shape (subjects, sessions, trials)."""
    n_subj, n_sess = PLANNED_SHAPE
    src = cls(n_subjects=n_subj, n_sessions=n_sess, n_trials=N_TRIALS, seed=seed, **SHORT, **sd)
    lead = np.array([t.meta["true_lead_ms"] for t in src]) - src.true_lead_ms
    return lead.reshape(n_subj, n_sess, N_TRIALS)


def spread(cls, seeds) -> dict:
    means = np.concatenate([offsets(cls, k).mean(axis=(1, 2)) for k in seeds])
    subj = np.concatenate([offsets(cls, k, session_sd_ms=0.0, trial_sd_ms=0.0)[:, 0, 0] for k in seeds])
    sess = np.concatenate([offsets(cls, k, subject_sd_ms=0.0, trial_sd_ms=0.0)[:, :, 0] for k in seeds])
    n, sd = means.size, float(means.std(ddof=1))
    lo, hi = (sd * np.sqrt((n - 1) / stats.chi2.ppf(q, n - 1)) for q in (0.975, 0.025))
    corr = [float(np.corrcoef(subj, sess[:, j])[0, 1]) for j in range(sess.shape[1])]
    return dict(n=n, sd=sd, lo=lo, hi=hi, corr=corr, ratio=float(np.median(sess[:, 0] / subj)))


def detected_sd(seed: int) -> float:
    """SD over subject means of the detected lead, 16x3x10 (the CLI defaults), as LeadStats uses them."""
    summ = LeadStats.summarize(Harness.run(SyntheticSource(seed=seed)), LatencyConfig(), n_boot=1)
    return float(np.std([x["mean_lead_ms"] for x in summ["subjects"]], ddof=1))


def recovers_known_lead() -> tuple[np.ndarray, float]:
    errs = []
    for seed in range(6):
        res = Harness.run(SyntheticSource(n_subjects=8, n_sessions=3, n_trials=5, seed=seed))
        errs.append(np.array([r.lead_ms - r.true_lead_ms for r in res]))
    return np.array([e.mean() for e in errs]), float(np.concatenate(errs).mean())


def delay(n_subj: int, n_sess: int, n_trials: int) -> dict:
    res = Harness.run(SyntheticSource(n_subjects=n_subj, n_sessions=n_sess, n_trials=n_trials, seed=0))
    return LeadStats.summarize(res, LatencyConfig(), n_boot=30)["confirmation_delay_ms"]


def report(pool_seeds: int) -> str:
    n_subj, n_sess = PLANNED_SHAPE
    design = SampleSize.design_subject_sd(8.0, 6.0, 10.0, n_sess, N_TRIALS)
    legacy_theory = float(np.sqrt((8.0 + 6.0 / n_sess) ** 2 + 36.0 * (n_sess - 1) / n_sess**2
                                  + 100.0 / (n_sess * N_TRIALS)))
    L = [f"## SyntheticSource spread (SYNTHETIC; {n_subj}x{n_sess}x{N_TRIALS}, numpy {np.__version__})", "",
         f"Design subject-mean SD {design:.3f} ms. Old seeding in theory: {legacy_theory:.3f} ms.", "",
         f"### 1. True-lead subject-mean SD, seeds 0..{pool_seeds - 1}", "",
         "| seeding | subjects | SD (ms) | 95% CI | corr(subject, session 0/1/2) | median sess. 0 / subj. |",
         "|---|---:|---:|---|---|---:|"]
    for name, cls in (("old list seeding", LegacySource), ("this checkout", SyntheticSource)):
        r = spread(cls, range(pool_seeds))
        L.append(f"| {name} | {r['n']} | {r['sd']:.3f} | [{r['lo']:.3f}, {r['hi']:.3f}] | "
                 + " / ".join(f"{c:+.3f}" for c in r["corr"]) + f" | {r['ratio']:.3f} |")
    L += ["", "### 2. True-lead subject-mean SD of one seed (16 subjects)", "",
          "| seeding | " + " | ".join(f"seed {k}" for k in range(6)) + " |", "|---|" + "---:|" * 6]
    for name, cls in (("old list seeding", LegacySource), ("this checkout", SyntheticSource)):
        sds = [offsets(cls, k).mean(axis=(1, 2)).std(ddof=1) for k in range(6)]
        L.append(f"| {name}, true lead | " + " | ".join(f"{v:.2f}" for v in sds) + " |")
    detected = [detected_sd(k) for k in range(6)]
    L.append("| this checkout, detected lead | " + " | ".join(f"{v:.2f}" for v in detected) + " |")
    means, pooled = recovers_known_lead()
    L += ["", "### 3. Lead error, 8x3x5, seeds 0-5 (test_recovers_known_lead)", "",
          f"Pooled {pooled:+.3f} ms; per seed " + ", ".join(f"{m:+.3f}" for m in means) + " ms.", "",
          "### 4. EMG confirmation delay (mean, min-max), seed 0", ""]
    for shape in ((4, 2, 5), (n_subj, n_sess, N_TRIALS)):
        d = delay(*shape)
        L.append(f"- {shape[0]}x{shape[1]}x{shape[2]}: {d['mean']:.2f} ms ({d['min']:.2f}-{d['max']:.2f})")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pool-seeds", type=int, default=40, help="pool seeds 0..N-1 in part 1 (default 40)")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    text = report(a.pool_seeds)
    print(text)
    if a.out:
        a.out.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
