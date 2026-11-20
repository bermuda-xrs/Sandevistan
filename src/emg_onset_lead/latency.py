from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LatencyConfig:
    """Exo intent-to-torque latency budget (ms).

    ``emg_extra_ms`` is t_win/2: half of the controller's 100 ms EMG window, so the
    default total is 18.67 + 50 = 68.67 ms. The controller shipped here is proportional,
    and that budget is the gate.

    ``triggered`` is off by default. Turn it on only for a triggered controller, which
    cannot act until the EMG onset is confirmed: the total budget then grows by the
    measured confirmation delay: the EMG hold time (``OnsetResult.trigger_s``) minus the
    back-dated EMG onset. ``Harness`` falls back to ``detected_s`` only when ``trigger_s`` is
    None. The EMG refine wait after the hold is not part of it.
    """

    imu_path_ms: float = 18.67
    emg_extra_ms: float = 50.0  # t_win/2, half of the 100 ms EMG window
    triggered: bool = False

    @property
    def total_ms(self) -> float:
        return round(self.imu_path_ms + self.emg_extra_ms, 6)

    def required_ms(self, path: str = "total", confirmation_delay_ms: float = 0.0) -> float:
        """Lead the EMG must provide to beat the given path ("total" or "imu").

        In triggered mode the total path adds ``confirmation_delay_ms``. Proportional
        mode ignores it, so the default gate stays 68.67 ms.
        """
        if path == "total":
            base = self.total_ms
        elif path == "imu":
            base = self.imu_path_ms
        else:
            raise ValueError(f"unknown path {path!r}")
        if self.triggered and path == "total":
            return base + float(confirmation_delay_ms)
        return base
