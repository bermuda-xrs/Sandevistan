from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .harness import Harness
from .latency import LatencyConfig
from .sources import SyntheticSource
from .stats import LeadStats


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="emg_onset_lead")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run the onset-lead analysis")
    r.add_argument("--synthetic", action="store_true", help="use the synthetic source (only one shipped)")
    r.add_argument("--subjects", type=int, default=8)
    r.add_argument("--sessions", type=int, default=3)
    r.add_argument("--trials", type=int, default=10, help="trials per session")
    r.add_argument("--seed", type=int, default=0)
    r.add_argument("--true-lead-ms", type=float, default=80.0)
    r.add_argument("--confidence", type=float, default=0.95)
    r.add_argument("--total-latency-ms", type=float, default=None,
                   help="override required latency (default: LatencyConfig total, 68.67)")
    r.add_argument("--out", type=Path, default=Path("out"))
    a = p.parse_args(argv)
    if not a.synthetic:
        p.error("only --synthetic is supported (no real dataset loaders are shipped)")

    lat = LatencyConfig()
    if a.total_latency_ms is not None:
        lat = LatencyConfig(imu_path_ms=a.total_latency_ms - lat.emg_extra_ms)
    src = SyntheticSource(a.subjects, a.sessions, a.trials, true_lead_ms=a.true_lead_ms, seed=a.seed)
    results = Harness.run(src, lat)
    summary = LeadStats.summarize(results, lat, a.confidence, seed=a.seed)
    a.out.mkdir(parents=True, exist_ok=True)
    Harness.write_csv(results, a.out / "trials.csv")
    Harness.write_summary(summary, a.out / "summary.json")
    pl = summary["pooled"]
    print(f"pooled lead {pl['mean_lead_ms']:.2f} ms, lower bound {pl['lower_bound_ms']:.2f} ms "
          f"vs required {pl['required_ms']:.2f} ms -> {'PASS' if summary['pass'] else 'FAIL'}"
          + (f" (ground-truth lead {pl['true_lead_ms_mean']:.2f} ms)" if "true_lead_ms_mean" in pl else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
