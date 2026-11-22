from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

import numpy as np

from .stats import PLANNED_SESSIONS, PLANNED_SUBJECTS


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

    n_subjects: int = PLANNED_SUBJECTS
    n_sessions: int = PLANNED_SESSIONS
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
    # EMG burst is drawn inside this window. A point range pins the onset (baseline tests).
    emg_onset_range_s: tuple[float, float] = (1.5, 2.5)

    def __iter__(self) -> Iterator[Trial]:
        # SeedSequence spawn tree: seed -> subject -> session -> trial. Child 0 of a subject or a
        # session is that node's own offset stream, children 1.. are its sessions or trials. Spawn
        # keys are mixed in after the entropy pool, so no two nodes share a stream. Seed lists do:
        # numpy pads a list shorter than the 4-word pool with zeros, so [seed, s] == [seed, s, 0].
        for s, subject in enumerate(np.random.SeedSequence(self.seed).spawn(self.n_subjects)):
            subj_seed, *session_nodes = subject.spawn(1 + self.n_sessions)
            subj_off = np.random.default_rng(subj_seed).normal(0, self.subject_sd_ms)
            for ss, session in enumerate(session_nodes):
                sess_seed, *trial_seeds = session.spawn(1 + self.n_trials)
                rng_s = np.random.default_rng(sess_seed)
                sess_off = rng_s.normal(0, self.session_sd_ms)
                gain = rng_s.uniform(0.7, 1.3)  # electrode re-don changes amplitude
                for t, trial_seed in enumerate(trial_seeds):
                    rng = np.random.default_rng(trial_seed)
                    yield self._trial(rng, s, ss, t, subj_off + sess_off, gain)

    def _trial(self, rng, s, ss, t, offset_ms, gain) -> Trial:
        lead_ms = self.true_lead_ms + offset_ms + rng.normal(0, self.trial_sd_ms)
        t_emg = float(rng.uniform(*self.emg_onset_range_s))
        t_move = t_emg + lead_ms / 1000.0

        te = np.arange(round(self.duration_s * self.fs_emg)) / self.fs_emg
        rise = np.clip((te - t_emg) / 0.04, 0, 1)
        amp = 1.0 + (self.emg_snr * gain - 1.0) * (3 * rise**2 - 2 * rise**3)
        emg = amp * rng.normal(0, 1, te.size)

        ti = np.arange(round(self.duration_s * self.fs_imu)) / self.fs_imu
        angle = self._angle(ti, t_move, rng)

        return Trial(f"S{s + 1:02d}", f"S{s + 1:02d}-D{ss + 1}", f"T{t + 1:02d}",
                     self.fs_emg, emg, self.fs_imu, angle,
                     {"true_lead_ms": float(lead_ms), "emg_gain": float(gain),
                      "motion": self._motion_name(), "t_move_s": float(t_move)})

    def _motion_name(self) -> str:
        return "critically_damped"

    def _angle(self, ti, t_move, rng) -> np.ndarray:
        """Critically damped step into flexion. Acceleration is nonzero at onset."""
        tau, span = 0.12, 60.0
        d = np.clip(ti - t_move, 0, None)
        angle = 20.0 + span * (1 - (1 + d / tau) * np.exp(-d / tau))
        return angle + rng.normal(0, self.angle_noise_deg, ti.size)


@dataclass(frozen=True)
class MinJerkSource(SyntheticSource):
    """Same trials as SyntheticSource, but the elbow follows a minimum-jerk path.

    Flash and Hogan 1985: a point-to-point move of duration ``move_s`` and size
    ``amplitude_deg`` has zero velocity and acceleration at both ends. The critically
    damped source does not (its acceleration steps at onset), so a detector tuned only
    on that shape reads the minimum-jerk onset late.
    """

    move_s: float = 0.5
    amplitude_deg: float = 60.0

    def _motion_name(self) -> str:
        return "minimum_jerk"

    def _angle(self, ti, t_move, rng) -> np.ndarray:
        if self.move_s <= 0:
            raise ValueError("move_s must be positive")
        d = ti - t_move
        y = np.zeros_like(ti)
        moving = (d > 0) & (d < self.move_s)
        u = d[moving] / self.move_s
        y[moving] = self.amplitude_deg * (10 * u**3 - 15 * u**4 + 6 * u**5)
        y[d >= self.move_s] = self.amplitude_deg
        return 20.0 + y + rng.normal(0, self.angle_noise_deg, ti.size)
