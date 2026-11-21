"""In-family lead bias and the option C guard band at this checkout. SYNTHETIC DATA ONLY.

No recorded EMG or motion is read. Each seed generates the planned 16 subjects x 3 sessions x 10
trials from ``SyntheticSource``: a critically damped 60 deg step, and EMG with a 40 ms smoothstep
burst at SNR 15 x U(0.7, 1.3) per session. The synthetic subject/session/trial SDs are 8/6/10 ms
and the true lead is 80 ms. Detectors and harness run at default settings. Per-trial error =
detected lead - true lead (ms). The 8/6/10 ms spread cancels out of it, because every trial is
compared with its own true lead. What remains is detector error on synthetic signals.

Per seed, the error is split by a nested ANOVA (method of moments), the same as in
``scripts/bias_split.py`` on analysis/emg-polarity-bias (3e9c27b):
  bias = mean of subject means
  u    = subject-cluster SE of the bias (SD of subject means / sqrt 16)
  u_new_subject = sqrt(s_subj^2 + s_sess^2/3 + s_e^2/30 + u^2)
Across seeds the script reports the mean bias, its seed-to-seed SD, and a t 95% CI over seeds.
Option C's guard band is w = k*u and gate = 68.67 + w; the gate itself is not changed here.

Usage (repo root):  python scripts/bias_guard_band.py [--seeds 6] [--out results.md]
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from emg_onset_lead import PLANNED_SHAPE, Harness, LatencyConfig, SyntheticSource  # noqa: E402

GATE_MS = LatencyConfig().total_ms  # 68.67; reported, never changed
N_TRIALS = 10


def run(seed: int) -> list[dict]:
    n_subj, n_sess = PLANNED_SHAPE
    src = SyntheticSource(n_subjects=n_subj, n_sessions=n_sess, n_trials=N_TRIALS, seed=seed)
    return [dict(subject=r.subject_id, session=r.session_id, err=r.lead_ms - r.true_lead_ms)
            for r in Harness.run(src)]


def components(rows: list[dict]) -> dict:
    """Nested ANOVA for subject > session > trial on non-missed trials (harmonic-mean n)."""
    by = defaultdict(lambda: defaultdict(list))
    missed = 0
    for r in rows:
        if np.isfinite(r["err"]):
            by[r["subject"]][r["session"]].append(r["err"])
        else:
            missed += 1
    subjects = sorted(by)
    a = len(subjects)
    b = len(next(iter(by.values())))
    n_t = [len(v) for s in subjects for v in by[s].values()]
    n = float(stats.hmean(n_t))
    sess = {(s, k): float(np.mean(v)) for s in subjects for k, v in by[s].items()}
    subj = {s: float(np.mean([sess[(s, k)] for k in by[s]])) for s in subjects}
    grand = float(np.mean(list(subj.values())))
    ss_e = sum(float(((np.array(v) - sess[(s, k)]) ** 2).sum()) for s in subjects for k, v in by[s].items())
    ss_b = n * sum((sess[(s, k)] - subj[s]) ** 2 for s in subjects for k in by[s])
    ss_a = n * b * sum((m - grand) ** 2 for m in subj.values())
    ms_e, ms_b, ms_a = ss_e / (sum(n_t) - a * b), ss_b / (a * (b - 1)), ss_a / (a - 1)
    var_e, var_sess, var_subj = ms_e, max((ms_b - ms_e) / n, 0.0), max((ms_a - ms_b) / (n * b), 0.0)
    subj_arr = np.array(list(subj.values()))
    u = float(subj_arr.std(ddof=1) / np.sqrt(a))
    return dict(n_trials=len(rows), missed=missed, bias=grand, u=u,
                sd_subj=var_subj ** 0.5, sd_sess=var_sess ** 0.5, sd_trial=var_e ** 0.5,
                u_new_subject=float(np.sqrt(var_subj + var_sess / b + var_e / (b * n) + u ** 2)))


def report(seeds: list[int]) -> str:
    per = {s: components(run(s)) for s in seeds}
    bias = np.array([c["bias"] for c in per.values()])
    us = np.array([c["u"] for c in per.values()])
    un = np.array([c["u_new_subject"] for c in per.values()])
    m, sd = float(bias.mean()), float(bias.std(ddof=1)) if bias.size > 1 else float("nan")
    h = float(stats.t.ppf(0.975, bias.size - 1) * sd / np.sqrt(bias.size)) if bias.size > 1 else float("nan")
    L = [f"## In-family lead bias and option C guard band (SYNTHETIC; SyntheticSource "
         f"{PLANNED_SHAPE[0]}x{PLANNED_SHAPE[1]}x{N_TRIALS}, seeds {seeds[0]}..{seeds[-1]})", "",
         "| seed | trials | missed | bias (ms) | u = cluster SE (ms) | subject SD | session SD | "
         "trial SD | u_new_subject (ms) |",
         "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for s, c in per.items():
        L.append(f"| {s} | {c['n_trials']} | {c['missed']} | {c['bias']:+.3f} | {c['u']:.3f} | "
                 f"{c['sd_subj']:.3f} | {c['sd_sess']:.3f} | {c['sd_trial']:.3f} | "
                 f"{c['u_new_subject']:.3f} |")
    L += ["", f"- Bias over seeds: mean {m:+.3f} ms, seed-to-seed SD {sd:.3f} ms, "
          f"95% CI over seeds [{m - h:+.3f}, {m + h:+.3f}] ms (t, {bias.size - 1} df).",
          f"- u (cluster SE of one study's bias): mean {us.mean():.3f} ms "
          f"(range {us.min():.3f}-{us.max():.3f}). u_new_subject: mean {un.mean():.3f} ms.",
          "- Negative bias = the lead reads short (EMG onset read late). The gate is not changed here.", "",
          "### Guard band w = k*u, gate = 68.67 + w (u = mean over seeds)", "",
          "| k | u (ms) | w (ms) | gate (ms) | w with u_new_subject (ms) | gate (ms) |",
          "|---:|---:|---:|---:|---:|---:|"]
    u, uns = float(us.mean()), float(un.mean())
    for k in (1, 2, 3):
        L.append(f"| {k} | {u:.3f} | {k * u:.3f} | {GATE_MS + k * u:.3f} | {k * uns:.3f} | "
                 f"{GATE_MS + k * uns:.3f} |")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seeds", type=int, default=6, help="use seeds 0..N-1 (default 6)")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    if a.seeds < 2:
        ap.error("--seeds must be >= 2 (the seed-to-seed SD needs two)")
    text = report(list(range(a.seeds)))
    print(text)
    if a.out:
        a.out.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
