from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LatencyConfig:
    """Exo intent-to-torque latency budget (ms). Defaults: current design."""

    imu_path_ms: float = 18.67
    emg_extra_ms: float = 50.0

    @property
    def total_ms(self) -> float:
        return round(self.imu_path_ms + self.emg_extra_ms, 6)

    def required_ms(self, path: str = "total") -> float:
        """Lead the EMG must provide to beat the given path ("total" or "imu")."""
        if path == "total":
            return self.total_ms
        if path == "imu":
            return self.imu_path_ms
        raise ValueError(f"unknown path {path!r}")
