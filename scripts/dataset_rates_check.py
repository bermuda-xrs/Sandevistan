"""Lead error at dataset sample rates, paired against 2000/150 Hz. SYNTHETIC DATA ONLY.

README "Sample rates are parameters". No recorded EMG or motion is read. Three parts, at this checkout:

1. Paired rate check. ``SyntheticSource`` and ``MinJerkSource``, 4 subjects x 3 sessions x 6 trials,
   seeds 100-115, SNR 15 (the source default). Each rate runs on the same seeds as 2000/150 Hz, so the
   true lead and the onset draws are shared trial by trial. Per seed, the shift is the mean of
   (error at the rate - error at 2000/150) over trials where both detected; the 95% CI is a t interval
   over the 16 seed means. Also: mean error, per-trial SD and misses at each rate.
2. CLI at rates too low for the 20-450 Hz band (the band is capped at 0.45 * fs, so below 1000 Hz it
   narrows): ``run --synthetic --subjects 2 --sessions 1 --trials 3`` (seed 0, design true lead 80 ms) at
   ``--fs-emg`` 2000 / 100 / 50 / 45. These run with exit 0; the pooled mean lead and misses are printed.
3. The draws of ``tests/test_dataset_rates.py::test_lead_accuracy_at_dataset_rates`` (seeds 0 and 1,
   4x3x6 trials, both sources) at 2000/150, 1000/125 and 2048/125 Hz: mean error, per-trial SD, misses.

4. CLI at low IMU rates: the same run (``--subjects 2 --sessions 1 --trials 3``, seed 0, EMG left at
   2000 Hz) at ``--fs-imu`` 150 / 125 / 100 / 50 / 25 / 20.5. The IMU smoothing needs ``fs_imu`` above
   20 Hz, so these all run with exit 0; the pooled mean lead and misses are printed next to the true lead.

Usage (repo root):  python scripts/dataset_rates_check.py [--procs 2] [--parts 1 2 3] [--out results.md]
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from emg_onset_lead import Harness, MinJerkSource, SyntheticSource  # noqa: E402
from emg_onset_lead.cli import main as cli_main  # noqa: E402

SEEDS = range(100, 116)
SHAPE = dict(n_subjects=4, n_sessions=3, n_trials=6)
BASE = (2000.0, 150.0)
RATES = [(1000.0, 125.0), (2048.0, 125.0), (2000.0, 125.0), (1000.0, 150.0)]
SOURCES = {"SyntheticSource": SyntheticSource, "MinJerkSource": MinJerkSource}
LOW_EMG = [2000.0, 100.0, 50.0, 45.0]
LOW_IMU = [150.0, 125.0, 100.0, 50.0, 25.0, 20.5]


def _errors(job):
    """Lead error (measured - true, ms) per seed x trial; NaN where a trial was missed."""
    name, fs_emg, fs_imu = job
    cls = SOURCES[name]
    rows = [[r.lead_ms - r.true_lead_ms
             for r in Harness.run(cls(**SHAPE, seed=s, fs_emg=fs_emg, fs_imu=fs_imu))] for s in SEEDS]
    return job, np.asarray(rows, float)


def paired_table(procs: int) -> list[str]:
    jobs = [(n, *fs) for n in SOURCES for fs in [BASE, *RATES]]
    if procs > 1:
        with Pool(procs) as pool:
            res = dict(pool.map(_errors, jobs))
    else:
        res = dict(map(_errors, jobs))
    out = ["| Source | Rates (Hz) | Mean error, 2000/150 -> rate (ms) | Shift (95% CI) ms | "
           "Per-trial SD, 2000/150 -> rate (ms) | SD ratio | Misses, 2000/150 / rate |",
           "|---|---|---|---|---|---|---|"]
    shifts, ratios, ci0, misses = [], [], [], []
    for name in SOURCES:
        b = res[(name, *BASE)]
        for fe, fi in RATES:
            x = res[(name, fe, fi)]
            d = np.nanmean(x - b, axis=1)  # per-seed paired shift
            m, se = d.mean(), d.std(ddof=1) / np.sqrt(d.size)
            half = stats.t.ppf(0.975, d.size - 1) * se
            sd_b, sd_x = np.nanstd(b, ddof=1), np.nanstd(x, ddof=1)
            nb, nx = int(np.isnan(b).sum()), int(np.isnan(x).sum())
            shifts.append(abs(m))
            ratios.append(sd_x / sd_b)
            ci0.append(m - half <= 0 <= m + half)
            misses += [nb, nx]
            out.append(f"| {name} | {fe:g}/{fi:g} | {np.nanmean(b):+.2f} -> {np.nanmean(x):+.2f} | "
                       f"{m:+.2f} ({m - half:+.2f}..{m + half:+.2f}) | {sd_b:.2f} -> {sd_x:.2f} | "
                       f"{sd_x / sd_b:.3f} | {nb} / {nx} of {b.size} |")
    out += ["", f"Max |shift| {max(shifts):.3f} ms; every 95% CI contains 0: {all(ci0)}; "
            f"per-trial SD change {100 * (min(ratios) - 1):+.1f}% to {100 * (max(ratios) - 1):+.1f}%; "
            f"misses per cell {min(misses)}-{max(misses)} of {len(SEEDS) * 72}."]
    return out


def low_rate_cli(flag: str = "--fs-emg", rates: list[float] = LOW_EMG) -> list[str]:
    out = [f"| {flag} (Hz) | exit | pooled mean lead (ms) | mean true lead, all trials (ms) | "
           "missed / trials |",
           "|---|---|---|---|---|"]
    for fe in rates:
        with tempfile.TemporaryDirectory() as d, contextlib.redirect_stdout(io.StringIO()):
            rc = cli_main(["run", "--synthetic", "--subjects", "2", "--sessions", "1", "--trials", "3",
                           flag, f"{fe:g}", "--out", d])
            p = json.loads((Path(d) / "summary.json").read_text())["pooled"]
            with open(Path(d) / "trials.csv", newline="") as f:
                true = np.mean([float(r["true_lead_ms"]) for r in csv.DictReader(f)])
        out.append(f"| {fe:g} | {rc} | {p['mean_lead_ms']:.1f} | {true:.1f} | "
                   f"{p['n_missed']} / {p['n_total']} |")
    return out


def suite_seed_draws() -> list[str]:
    out = ["| Source | Rates (Hz) | mean error (ms) | per-trial SD (ms) | missed / trials |",
           "|---|---|---|---|---|"]
    for name, cls in SOURCES.items():
        for fe, fi in [BASE, (1000.0, 125.0), (2048.0, 125.0)]:
            e = np.asarray([r.lead_ms - r.true_lead_ms for seed in range(2)
                            for r in Harness.run(cls(**SHAPE, seed=seed, fs_emg=fe, fs_imu=fi))], float)
            ok = e[np.isfinite(e)]
            out.append(f"| {name} | {fe:g}/{fi:g} | {ok.mean():+.2f} | {ok.std(ddof=1):.2f} | "
                       f"{e.size - ok.size} / {e.size} |")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--parts", type=int, nargs="+", default=[1, 2, 3, 4], choices=[1, 2, 3, 4])
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    lines = ["# Dataset sample rates (SYNTHETIC)", ""]
    if 1 in a.parts:
        lines += ["## 1. Paired rate check, seeds 100-115, 4x3x6 trials, SNR 15", "",
                  *paired_table(a.procs), ""]
    if 2 in a.parts:
        lines += ["## 2. CLI at low EMG rates, 2x1x3 trials, seed 0", "", *low_rate_cli(), ""]
    if 3 in a.parts:
        lines += ["## 3. Test-suite draws, seeds 0-1, 4x3x6 trials", "", *suite_seed_draws(), ""]
    if 4 in a.parts:
        lines += ["## 4. CLI at low IMU rates, 2x1x3 trials, seed 0", "",
                  *low_rate_cli("--fs-imu", LOW_IMU), ""]
    text = "\n".join(lines)
    print(text)
    if a.out:
        a.out.write_text(text)


if __name__ == "__main__":
    main()
