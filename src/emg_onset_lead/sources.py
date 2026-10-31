from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

import numpy as np


@dataclass(frozen=True)
class Trial:
    subject_id: str
    session_id: str
    trial_id: str
    fs_emg: float
    emg: np.ndarray
    fs_imu: float
    angle_or_imu: np.ndarray  # 1-D joint angle (or IMU-derived angle)
    meta: dict = field(default_factory=dict)  # optional "true_lead_ms" for synthetic data


@runtime_checkable
class TrialSource(Protocol):
    def __iter__(self) -> Iterator[Trial]: ...


@dataclass(frozen=True)
class SyntheticSource:
    """Seeded elbow-flexion trials with a known electromechanical delay (EMD).

    True lead of a trial = motion start - EMG activation start
                         = true_lead_ms + subject offset + session offset + trial jitter.
    Session offsets stand in for electrode re-donning between days.
    """

    n_subjects: int = 8
    n_sessions: int = 3
    n_trials: int = 10
    true_lead_ms: float = 80.0
    subject_sd_ms: float = 8.0
    session_sd_ms: float = 6.0
    trial_sd_ms: float = 10.0
    fs_emg: float = 2000.0
    fs_imu: float = 150.0
    duration_s: float = 4.0
    emg_snr: float = 15.0  # burst amplitude / baseline noise amplitude (per session jitter applied)
    angle_noise_deg: float = 0.3
    seed: int = 0

    def __iter__(self) -> Iterator[Trial]:
        for s in range(self.n_subjects):
            subj_off = np.random.default_rng([self.seed, s]).normal(0, self.subject_sd_ms)
            for ss in range(self.n_sessions):
                rng_s = np.random.default_rng([self.seed, s, ss])
                sess_off = rng_s.normal(0, self.session_sd_ms)
                gain = rng_s.uniform(0.7, 1.3)  # electrode re-don changes amplitude
                for t in range(self.n_trials):
                    rng = np.random.default_rng([self.seed, s, ss, t + 1])
                    yield self._trial(rng, s, ss, t, subj_off + sess_off, gain)

    def _trial(self, rng, s, ss, t, offset_ms, gain) -> Trial:
        lead_ms = self.true_lead_ms + offset_ms + rng.normal(0, self.trial_sd_ms)
        t_emg = rng.uniform(1.5, 2.5)
        t_move = t_emg + lead_ms / 1000.0

        te = np.arange(round(self.duration_s * self.fs_emg)) / self.fs_emg
        rise = np.clip((te - t_emg) / 0.04, 0, 1)
        amp = 1.0 + (self.emg_snr * gain - 1.0) * (3 * rise**2 - 2 * rise**3)
        emg = amp * rng.normal(0, 1, te.size)

        ti = np.arange(round(self.duration_s * self.fs_imu)) / self.fs_imu
        tau, span = 0.12, 60.0  # critically damped step response, 60 deg flexion
        d = np.clip(ti - t_move, 0, None)
        angle = 20.0 + span * (1 - (1 + d / tau) * np.exp(-d / tau))
        angle = angle + rng.normal(0, self.angle_noise_deg, ti.size)

        return Trial(f"S{s + 1:02d}", f"S{s + 1:02d}-D{ss + 1}", f"T{t + 1:02d}",
                     self.fs_emg, emg, self.fs_imu, angle,
                     {"true_lead_ms": float(lead_ms), "emg_gain": float(gain)})
