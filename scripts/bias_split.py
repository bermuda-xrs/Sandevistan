"""Split the synthetic onset-lead error by EMG polarity, subject and session.

SYNTHETIC ONLY. No recorded EMG or motion is read. Trials come from ``SyntheticSource``
(critically damped 60 deg step, Gaussian EMG burst) with its synthetic variance components
(subject / session / trial SD 8 / 6 / 10 ms, true lead 80 ms). Detectors and harness run at their
default settings. The planned 16 x 3 x 10 study is generated once with a fixed seed, then run
twice: EMG as generated (polarity +1) and EMG multiplied by -1 (polarity -1). Polarity is not a
source parameter, so the sign flip is the only change between the two runs.

Per-trial error = detected lead - true lead (ms). A trial with a NaN lead is counted as missed.

Usage (from the repo root):
    python scripts/bias_split.py [--seed 0] [--out results.md]
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy import stats

# Use this checkout's src/, not whatever copy happens to be installed.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from emg_onset_lead import PLANNED_SHAPE, Harness, LatencyConfig, SyntheticSource  # noqa: E402

GATE_MS = LatencyConfig().total_ms  # 68.67, not changed here
N_TRIALS = 10


class FlippedEmg:
    """Wrap a source and multiply every EMG sample by ``sign``. Angle and meta are unchanged."""

    def __init__(self, source, sign: float):
        self.source, self.sign = source, sign

    def __iter__(self):
        for t in self.source:
            yield replace(t, emg=self.sign * np.asarray(t.emg))


def run(seed: int) -> list[dict]:
    n_subj, n_sess = PLANNED_SHAPE
    src = SyntheticSource(n_subjects=n_subj, n_sessions=n_sess, n_trials=N_TRIALS, seed=seed)
    rows = []
    for pol in (+1, -1):
        for r in Harness.run(FlippedEmg(src, float(pol))):
            err = r.lead_ms - r.true_lead_ms
            rows.append(dict(polarity=pol, subject=r.subject_id, session=r.session_id,
                             trial=r.trial_id, err=err, emg_onset_s=r.emg_onset_s,
                             missed=not np.isfinite(err)))
    return rows


def describe(v: np.ndarray) -> dict:
    v = v[np.isfinite(v)]
    n = v.size
    if n == 0:
        return dict(n=0, mean=np.nan, sd=np.nan, lo=np.nan, hi=np.nan)
    m = float(v.mean())
    sd = float(v.std(ddof=1)) if n > 1 else np.nan
    if n > 1:
        h = float(stats.t.ppf(0.975, n - 1)) * sd / np.sqrt(n)
        lo, hi = m - h, m + h
    else:
        lo = hi = np.nan
    return dict(n=n, mean=m, sd=sd, lo=lo, hi=hi)


def group(rows, key) -> dict:
    g = defaultdict(list)
    miss = defaultdict(int)
    for r in rows:
        k = key(r)
        g[k].append(r["err"])
        miss[k] += r["missed"]
    return {k: {**describe(np.array(v)), "missed": miss[k], "total": len(v)} for k, v in g.items()}


def table(title: str, groups: dict) -> list[str]:
    out = [f"### {title}", "",
           "| group | n | missed | mean bias (ms) | SD (ms) | 95% CI (ms) |",
           "|---|---:|---:|---:|---:|---|"]
    for k in sorted(groups):
        s = groups[k]
        name = " x ".join(f"{x:+d}" if isinstance(x, int) else str(x)
                          for x in (k if isinstance(k, tuple) else (k,)))
        out.append(f"| {name} | {s['n']} | {s['missed']}/{s['total']} | {s['mean']:+.3f} | "
                   f"{s['sd']:.3f} | [{s['lo']:+.3f}, {s['hi']:+.3f}] |")
    return out + [""]


def nested_components(rows) -> dict:
    """Method-of-moments (ANOVA) variance split for subject > session > trial.

    The design is balanced when nothing is missed. Missed trials are dropped, and the session
    trial count is then the harmonic mean, which is an approximation and is reported as such.
    """
    by = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if not r["missed"]:
            by[r["subject"]][r["session"]].append(r["err"])
    subjects = sorted(by)
    a = len(subjects)
    b = len(next(iter(by.values())))
    n_t = [len(v) for s in subjects for v in by[s].values()]
    n = float(stats.hmean(n_t))
    sess_means = {(s, k): float(np.mean(v)) for s in subjects for k, v in by[s].items()}
    subj_means = {s: float(np.mean([sess_means[(s, k)] for k in by[s]])) for s in subjects}
    grand = float(np.mean(list(subj_means.values())))
    ss_e = sum(float(((np.array(v) - sess_means[(s, k)]) ** 2).sum())
               for s in subjects for k, v in by[s].items())
    df_e = sum(n_t) - a * b
    ss_b = n * sum((sess_means[(s, k)] - subj_means[s]) ** 2 for s in subjects for k in by[s])
    df_b = a * (b - 1)
    ss_a = n * b * sum((m - grand) ** 2 for m in subj_means.values())
    df_a = a - 1
    ms_e, ms_b, ms_a = ss_e / df_e, ss_b / df_b, ss_a / df_a
    var_e = ms_e
    var_sess = max((ms_b - ms_e) / n, 0.0)
    var_subj = max((ms_a - ms_b) / (n * b), 0.0)
    subj_arr = np.array(list(subj_means.values()))
    se = float(subj_arr.std(ddof=1) / np.sqrt(a))  # cluster SE: includes subject + session + trial
    return dict(a=a, b=b, n=n, grand=grand, ms_a=ms_a, ms_b=ms_b, ms_e=ms_e,
                df_a=df_a, df_b=df_b, df_e=df_e,
                F_subj=ms_a / ms_b, p_subj=float(stats.f.sf(ms_a / ms_b, df_a, df_b)),
                F_sess=ms_b / ms_e, p_sess=float(stats.f.sf(ms_b / ms_e, df_b, df_e)),
                sd_subj=var_subj ** 0.5, sd_sess=var_sess ** 0.5, sd_e=var_e ** 0.5,
                sd_subj_means=float(subj_arr.std(ddof=1)),
                sd_sess_means=float(np.std(list(sess_means.values()), ddof=1)),
                se=se, se_formula=float(np.sqrt(var_subj / a + var_sess / (a * b) + var_e / (a * b * n))),
                u_new_subject=float(np.sqrt(var_subj + var_sess / b + var_e / (b * n) + se ** 2)))


def paired_polarity(rows) -> dict:
    """Per-trial difference error(-1) - error(+1), and its subject-cluster SE."""
    pos = {(r["subject"], r["session"], r["trial"]): r for r in rows if r["polarity"] == +1}
    d, onset_d, by_subj = [], [], defaultdict(list)
    for r in rows:
        if r["polarity"] != -1:
            continue
        p = pos[(r["subject"], r["session"], r["trial"])]
        if r["missed"] or p["missed"]:
            continue
        d.append(r["err"] - p["err"])
        onset_d.append((r["emg_onset_s"] - p["emg_onset_s"]) * 1e3)
        by_subj[r["subject"]].append(d[-1])
    d = np.array(d)
    sm = np.array([np.mean(v) for v in by_subj.values()])
    return dict(n=d.size, mean=float(d.mean()), max_abs=float(np.max(np.abs(d))),
                max_abs_onset_ms=float(np.max(np.abs(onset_d))), n_nonzero=int(np.count_nonzero(d)),
                se=float(sm.std(ddof=1) / np.sqrt(sm.size)))


def report(rows, seed: int) -> str:
    L = [f"## Synthetic bias split (seed {seed}, {PLANNED_SHAPE[0]}x{PLANNED_SHAPE[1]}x{N_TRIALS}, "
         f"both EMG polarities)", ""]
    n_tot = len(rows)
    n_miss = sum(r["missed"] for r in rows)
    L.append(f"Trials: {n_tot} ({n_tot // 2} per polarity); missed (NaN lead): {n_miss}.")
    L.append("")
    pp = paired_polarity(rows)
    L += ["### Paired polarity check (same trial, EMG x -1 vs x +1)", "",
          f"- paired trials: {pp['n']}; trials whose error differs between polarities: {pp['n_nonzero']}",
          f"- mean difference error(-1) - error(+1): {pp['mean']:+.6f} ms; max |difference|: "
          f"{pp['max_abs']:.3g} ms; max |EMG onset difference|: {pp['max_abs_onset_ms']:.3g} ms",
          f"- subject-cluster SE of the polarity difference: {pp['se']:.3g} ms", ""]
    L += table("(a) By polarity", group(rows, lambda r: r["polarity"]))
    pos = [r for r in rows if r["polarity"] == +1]
    L += table("(b) By subject (polarity +1; polarity -1 shown in (d))", group(pos, lambda r: r["subject"]))
    L += table("(c) By session (polarity +1)", group(pos, lambda r: r["session"]))
    L += table("(d) Polarity x subject", group(rows, lambda r: (r["polarity"], r["subject"])))
    L += ["### Variance split (nested ANOVA, method of moments) and correction uncertainty", "",
          "| polarity | grand mean (ms) | SD of subject means | SD of session means | "
          "subject SD | session SD | trial (within) SD | F subj (p) | F sess (p) | "
          "u = cluster SE (ms) | u_new_subject (ms) |",
          "|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---:|"]
    comps = {}
    for pol in (+1, -1):
        c = nested_components([r for r in rows if r["polarity"] == pol])
        comps[pol] = c
        L.append(f"| {pol:+d} | {c['grand']:+.3f} | {c['sd_subj_means']:.3f} | {c['sd_sess_means']:.3f} | "
                 f"{c['sd_subj']:.3f} | {c['sd_sess']:.3f} | {c['sd_e']:.3f} | "
                 f"{c['F_subj']:.2f} ({c['p_subj']:.3g}) | {c['F_sess']:.2f} ({c['p_sess']:.3g}) | "
                 f"{c['se']:.3f} | {c['u_new_subject']:.3f} |")
    c = comps[+1]
    L += ["", f"Session trial count used in the split: {c['n']:.2f} (harmonic mean; 10 if nothing missed). "
          f"Formula SE sqrt(s_subj^2/a + s_sess^2/(ab) + s_e^2/(abn)) = {c['se_formula']:.3f} ms "
          f"(polarity +1), vs cluster SE {c['se']:.3f} ms.", ""]
    L += ["### Guard band w = k*u, gate = 68.67 + w ms (u = polarity +1 cluster SE)", "",
          "| k | u (ms) | w (ms) | gate (ms) | w with u_new_subject (ms) | gate (ms) |",
          "|---:|---:|---:|---:|---:|---:|"]
    u, un = c["se"], c["u_new_subject"]
    for k in (1, 2, 3):
        L.append(f"| {k} | {u:.3f} | {k * u:.3f} | {GATE_MS:.2f} + {k * u:.3f} = {GATE_MS + k * u:.3f} | "
                 f"{k * un:.3f} | {GATE_MS + k * un:.3f} |")
    L.append("")
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    text = report(run(a.seed), a.seed)
    print(text)
    if a.out:
        a.out.write_text(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
