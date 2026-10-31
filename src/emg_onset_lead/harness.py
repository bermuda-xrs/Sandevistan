from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from .detect import EmgOnsetDetector, EmgParams, ImuOnsetDetector, ImuParams
from .latency import LatencyConfig
from .sources import TrialSource

NAN = float("nan")


@dataclass(frozen=True)
class TrialResult:
    subject_id: str
    session_id: str
    trial_id: str
    emg_onset_s: float
    imu_onset_s: float
    lead_ms: float  # imu_onset - emg_onset
    margin_ms: float  # lead_ms - required latency
    true_lead_ms: float = NAN  # only when the source supplies meta["true_lead_ms"]
    detect_lead_ms: float = NAN  # confirmation-time difference; diagnostic, not gated


class Harness:
    @staticmethod
    def run(source: TrialSource, latency: LatencyConfig = LatencyConfig(),
            emg: EmgParams = EmgParams(), imu: ImuParams = ImuParams(),
            path: str = "total", chunk_s: float = 0.05) -> list[TrialResult]:
        """Detect onsets chunk by chunk (live-style) for every trial."""
        req = latency.required_ms(path)
        out = []
        for t in source:
            e = EmgOnsetDetector.run(t.emg, t.fs_emg, emg, max(1, round(chunk_s * t.fs_emg)))
            i = ImuOnsetDetector.run(t.angle_or_imu, t.fs_imu, imu,
                                     max(1, round(chunk_s * t.fs_imu)))
            es = NAN if e.onset_s is None else e.onset_s
            is_ = NAN if i.onset_s is None else i.onset_s
            lead = (is_ - es) * 1000.0
            dl = NAN if e.detected_s is None or i.detected_s is None else (i.detected_s - e.detected_s) * 1e3
            out.append(TrialResult(t.subject_id, t.session_id, t.trial_id, es, is_, lead,
                                   lead - req, float(t.meta.get("true_lead_ms", NAN)), dl))
        return out

    @staticmethod
    def write_csv(results: list[TrialResult], path: Path) -> None:
        fields = list(TrialResult.__dataclass_fields__)
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(asdict(r) for r in results)

    @staticmethod
    def write_summary(summary: dict, path: Path) -> None:
        Path(path).write_text(json.dumps(summary, indent=2, allow_nan=False))
