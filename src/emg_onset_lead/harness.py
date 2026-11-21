from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

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
    # EMG trigger_s (hold met) - onset_s; added to the budget only if triggered
    emg_confirmation_delay_ms: float = NAN
    baseline_rest_ok: bool = True


class Harness:
    @staticmethod
    def run(source: TrialSource, latency: LatencyConfig = LatencyConfig(),
            emg: EmgParams = EmgParams(), imu: ImuParams = ImuParams(),
            path: str = "total", chunk_s: float = 0.05) -> list[TrialResult]:
        """Detect onsets chunk by chunk (live-style) for every trial."""
        found = []
        for t in source:
            e = EmgOnsetDetector.run(t.emg, t.fs_emg, emg, max(1, round(chunk_s * t.fs_emg)))
            i = ImuOnsetDetector.run(t.angle_or_imu, t.fs_imu, imu,
                                     max(1, round(chunk_s * t.fs_imu)))
            found.append((t, e, i))
        delays = [Harness._confirmation_delay_ms(e) for _, e, _ in found]
        delays = [d for d in delays if np.isfinite(d)]
        # One budget for the run: in triggered mode, the mean EMG confirmation delay.
        req = latency.required_ms(path, float(np.mean(delays)) if delays else 0.0)
        out = []
        for t, e, i in found:
            ok = bool(e.baseline_rest_ok and i.baseline_rest_ok)
            delay = Harness._confirmation_delay_ms(e) if ok else NAN
            if not ok or e.onset_s is None or i.onset_s is None:
                es = is_ = lead = NAN
                dl = NAN
            else:
                es, is_ = e.onset_s, i.onset_s
                lead = (is_ - es) * 1000.0
                if e.detected_s is None or i.detected_s is None:
                    dl = NAN
                else:
                    dl = (i.detected_s - e.detected_s) * 1e3
            out.append(TrialResult(t.subject_id, t.session_id, t.trial_id, es, is_, lead,
                                   lead - req, float(t.meta.get("true_lead_ms", NAN)), dl, delay, ok))
        return out

    @staticmethod
    def _confirmation_delay_ms(emg_result) -> float:
        """Hold time minus back-dated onset: what a triggered controller waits.

        A triggered controller fires when the hold is met (``trigger_s``). The EMG onset
        estimate commits ``refine_s`` later, but the controller does not need it to act, so
        that wait is not part of the confirmation delay.
        """
        if emg_result.onset_s is None or not emg_result.baseline_rest_ok:
            return NAN
        t = emg_result.trigger_s if emg_result.trigger_s is not None else emg_result.detected_s
        if t is None:
            return NAN
        return (t - emg_result.onset_s) * 1000.0

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
