from __future__ import annotations

import argparse
import sys
from dataclasses import replace
from pathlib import Path

from .harness import Harness
from .latency import LatencyConfig
from .sources import DEFAULT_FS_EMG, DEFAULT_FS_IMU, SyntheticSource, check_rate
from .stats import PLANNED_SESSIONS, PLANNED_SUBJECTS, LeadStats, SampleSize

# --total-latency-ms is compared to imu_path_ms + emg_extra_ms, which LatencyConfig rounds to 1e-6 ms.
TOTAL_TOL_MS = 1e-6


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="emg_onset_lead")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run the onset-lead analysis")
    r.add_argument("--synthetic", action="store_true", help="use the synthetic source (only one shipped)")
    r.add_argument("--subjects", type=int, default=PLANNED_SUBJECTS)
    r.add_argument("--sessions", type=int, default=PLANNED_SESSIONS)
    r.add_argument("--trials", type=int, default=10, help="trials per session")
    r.add_argument("--seed", type=int, default=0)
    r.add_argument("--seeds", type=int, default=None,
                   help="pool seeds 0..N-1. Default is a single seed (--seed); subject-mean SD is one draw")
    r.add_argument("--true-lead-ms", type=float, default=80.0)
    r.add_argument("--confidence", type=float, default=0.95)
    r.add_argument("--imu-path-ms", type=float, default=None,
                   help="IMU-path latency (default 18.67). The only way to change the budget")
    r.add_argument("--total-latency-ms", type=float, default=None,
                   help="check only: must equal imu path + 50 ms EMG extra (68.67 by default). "
                        "It never changes --imu-path-ms; a mismatch is an error")
    r.add_argument("--triggered", action="store_true",
                   help="add measured EMG confirmation delay to the 68.67 ms budget")
    r.add_argument("--fs-emg", type=float, default=DEFAULT_FS_EMG,
                   help=f"synthetic EMG sample rate in Hz (default {DEFAULT_FS_EMG:g}; "
                        "Lucchetti 2025: 1000, Rojas-Martinez 2020: 2048)")
    r.add_argument("--fs-imu", type=float, default=DEFAULT_FS_IMU,
                   help=f"synthetic angle/IMU sample rate in Hz (default {DEFAULT_FS_IMU:g}; "
                        "Lucchetti 2025 markers: 125)")
    r.add_argument("--out", type=Path, default=Path("out"))
    a = p.parse_args(argv)
    if not a.synthetic:
        p.error("only --synthetic is supported (no real dataset loaders are shipped)")

    try:
        fs_emg, fs_imu = check_rate("--fs-emg", a.fs_emg), check_rate("--fs-imu", a.fs_imu)
    except ValueError as exc:
        p.error(str(exc))

    lat = LatencyConfig(triggered=a.triggered)
    if a.imu_path_ms is not None:
        lat = LatencyConfig(imu_path_ms=a.imu_path_ms, triggered=a.triggered)
    if a.total_latency_ms is not None and abs(a.total_latency_ms - lat.total_ms) > TOTAL_TOL_MS:
        p.error(f"--total-latency-ms {a.total_latency_ms:g} does not match imu path "
                f"{lat.imu_path_ms:g} + EMG extra {lat.emg_extra_ms:g} = {lat.total_ms:g} ms. "
                "The total is not used to back out the IMU path; set --imu-path-ms explicitly")
    seeds = [a.seed] if a.seeds is None else list(range(a.seeds))
    if not seeds:
        p.error("--seeds must be >= 1")
    results = []
    for seed in seeds:
        src = SyntheticSource(a.subjects, a.sessions, a.trials, true_lead_ms=a.true_lead_ms, seed=seed,
                              fs_emg=fs_emg, fs_imu=fs_imu)
        part = Harness.run(src, lat)
        if len(seeds) > 1:
            part = [replace(row, subject_id=f"{row.subject_id}-seed{seed}",
                            session_id=f"{row.session_id}-seed{seed}") for row in part]
        results.extend(part)
    summary = LeadStats.summarize(results, lat, a.confidence, seed=seeds[0])
    design_sd = SampleSize.design_subject_sd(8.0, 6.0, 10.0, PLANNED_SESSIONS, 10)
    summary["demo"] = {
        "single_seed": len(seeds) == 1,
        "seeds": seeds,
        "design_subject_sd_ms": design_sd,
        "fs_emg_hz": fs_emg,
        "fs_imu_hz": fs_imu,
        "note": (
            "Single-seed demo. The subject-mean SD is one draw; under the synthetic variances "
            f"it is expected to be {design_sd:.2f} ms (subject/session/trial SDs 8/6/10 ms, "
            f"{PLANNED_SESSIONS} sessions x 10 trials per subject), "
            "not whatever this seed happened to show."
            if len(seeds) == 1 else
            f"Pooled across seeds {seeds[0]}..{seeds[-1]}. Subject ids are suffixed, so this cohort "
            f"is not the planned {PLANNED_SUBJECTS}x{PLANNED_SESSIONS} and the shape gate does not pass."
        ),
    }
    a.out.mkdir(parents=True, exist_ok=True)
    Harness.write_csv(results, a.out / "trials.csv")
    Harness.write_summary(summary, a.out / "summary.json")
    pl = summary["pooled"]
    delay = summary["confirmation_delay_ms"]
    delay_txt = (f", EMG confirmation delay {delay['mean']:.2f} ms "
                 f"(range {delay['min']:.2f}-{delay['max']:.2f})" if delay["mean"] is not None else "")
    seed_txt = f"single-seed demo (seed {seeds[0]})" if len(seeds) == 1 else f"seeds {seeds[0]}..{seeds[-1]}"
    print(f"{seed_txt}: pooled lead {pl['mean_lead_ms']:.2f} ms, "
          f"lower bound {pl['lower_bound_ms']:.2f} ms "
          f"vs required {pl['required_ms']:.2f} ms -> {'PASS' if summary['pass'] else 'FAIL'}"
          + delay_txt
          + (f" (ground-truth lead {pl['true_lead_ms_mean']:.2f} ms)" if "true_lead_ms_mean" in pl else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
