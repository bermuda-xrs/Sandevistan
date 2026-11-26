"""Checks on the ML ramp kill-test outputs. SYNTHETIC DATA ONLY.

1. Oracle floor: score the same studies with zero EMG error (true onset known exactly) and the IMU error
   as simulated, and with zero error in both. The zero-error false pass is the gate's own rate; the kill
   limit 0.05 should be read against it.
2. One-sided binomial p against 0.05 per cell and Holm-adjusted across the cells, for context. The kill
   rule itself is the stated one (any cell's one-sided 95% Clopper-Pearson lower bound > 0.05), unadjusted.
3. Where the bias comes from: the fitted ramp start, rise and gain against the truth.
4. What the rare IMU noise triggers do: their count, and oracle pass rates with those trials dropped.

  python scripts/ml_ramp_row12_checks.py --raw 'out/raw/raw_*.npz' \
      --scored out/scored.json --out out/checks.md
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / "src"), str(HERE)]

import ml_ramp_row12_sim as sim  # noqa: E402


def load(files):
    Z = [dict(np.load(f)) for f in sorted(files)]
    return {k: np.concatenate([z[k] for z in Z]) for k in ("off", "imu", "studies", "fit")}


def _gate_job(a):
    off, err, L = a
    tr = L + off
    return sim.gate(tr + err, tr)[0]


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.maximum.accumulate(p[o] * (p.size - np.arange(p.size)))
    out = np.empty_like(p)
    out[o] = np.minimum(adj, 1)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--raw", required=True)
    ap.add_argument("--scored", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    D = load(glob.glob(a.raw))
    S = json.load(open(a.scored))
    n = len(D["studies"])
    assert n == len(S)
    out = [f"Studies: {n}. ALL SYNTHETIC.", ""]
    variants = {}
    with Pool(a.procs) as p:
        for name, errs in (("no EMG or IMU error", lambda li: np.zeros((n, 16, 3, 10))),
                           ("no EMG error; simulated IMU error",
                            lambda li: np.where(np.isnan(D["imu"][:, li]), np.nan, D["imu"][:, li])),
                           ("no EMG error; IMU error, |IMU error| > 100 ms trials dropped",
                            lambda li: np.where(np.abs(D["imu"][:, li]) > 100, np.nan, D["imu"][:, li]))):
            for li, L in enumerate(sim.LEADS):
                E = errs(li)
                k = sum(p.map(_gate_job, [(D["off"][s], E[s], L) for s in range(n)], chunksize=32))
                variants[(name, L)] = k
    out += ["### Oracle floor (gate pass rates with perfect EMG onsets)", "",
            "| Case | Power @ 74.67 | False pass @ 68.67 |", "|---|---|---|"]
    for name in dict.fromkeys(k[0] for k in variants):
        out.append(f"| {name} | {sim.fmt(variants[(name, sim.LEADS[0])], n)} | "
                   f"{sim.fmt(variants[(name, sim.LEADS[1])], n)} |")
    im = D["imu"]
    big = int((np.abs(im) > 100).sum())
    out += ["", f"IMU errors beyond 100 ms: {big} of {im.size} IMU onsets (2 leads share noise); "
            f"trial miss fraction {np.isnan(im).mean():.5f}; IMU error SD excluding those "
            f"{np.nanstd(np.where(np.abs(im) > 100, np.nan, im)):.2f} ms.", ""]
    # multiplicity
    names = [f"{sh} {r} ms SNR {q}" for sh, r, q in sim.CELLS]
    kf = [int(sum(s[f"{sh}_{r}_{q}_{sim.LEADS[1]}_{sim.P}"][0] for s in S)) for sh, r, q in sim.CELLS]
    p1 = [stats.binomtest(k, n, 0.05, alternative="greater").pvalue for k in kf]
    adj = holm(p1)
    out += ["### False pass vs 0.05: one-sided binomial p, unadjusted and Holm over 20 cells", "",
            "| Cell | False pass | CP 1-s LB | p | Holm p | Kill rule (LB > 0.05) |",
            "|---|---|---|---|---|---|"]
    for nm, k, p, h in zip(names, kf, p1, adj):
        lb = sim.cp_lower_1s(k, n)
        out.append(f"| {nm} | {sim.fmt(k, n)} | {lb:.4f} | {p:.2g} | {h:.2g} | "
                   f"{'KILL' if lb > 0.05 else 'no'} |")
    out += ["", f"Cells with Holm p < 0.05: {int((adj < 0.05).sum())}.", ""]
    # fitted parameters
    out += ["### Fitted ramp vs truth, all trials (EMG only)", "",
            "| Cell | start error mean ± SD (ms) | fitted rise mean (ms) | true rise to plateau | "
            "fitted gain mean | P=0.10 crossing error: mean / SD / 1st-99th pct (ms) |",
            "|---|---|---|---|---|---|"]
    for c, (sh, r, q) in enumerate(sim.CELLS):
        f = D["fit"][:, c].astype(float).reshape(-1, 3)
        e = sim.emg_error_ms(f, sh, r, sim.P)
        rise = f"{r}" if sh == "smooth" else f"asymptotic, tau {sim.TAU_PER_R * r:.1f}"
        out.append(f"| {sh} {r} SNR {q} | {np.nanmean(f[:, 0]):+.2f} ± {np.nanstd(f[:, 0]):.2f} | "
                   f"{np.nanmean(f[:, 1]):.1f} | {rise} | {np.nanmean(f[:, 2]):.2f} | "
                   f"{np.nanmean(e):+.2f} / {np.nanstd(e):.2f} / "
                   f"{np.nanpercentile(e, 1):+.1f}..{np.nanpercentile(e, 99):+.1f} |")
    Path(a.out).write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
