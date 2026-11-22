"""SyntheticSource random streams: independent, reproducible, design SDs 8/6/10 ms. Synthetic data only.

The source once seeded the subject with ``default_rng([seed, s])`` and the session with
``default_rng([seed, s, ss])``. numpy's SeedSequence pads a seed shorter than its pool (4 words of 32
bits by default; the count is of words after conversion, not of list items) with zeros, so session 0
read the subject's stream. Its offset was 0.75 x (= 6/8) the subject offset, and the subject-mean
true-lead SD was sqrt((8 + 6/3)^2 + 2 * 36/9 + 100/30) = 10.55 ms (theory) instead of 8.91.

Two kinds of check:
- **Streams** (``_record``). Every generator the source makes is replaced by one that remembers where
  its stream started and how many values it handed out. Two streams overlap if any raw 64-bit word
  they consume is the same, so a stream reused with an offset is caught, not only one restarted from
  its first draw. Seeds are pooled, so reuse across seeds is caught too. The number of values handed
  out must equal the number in the trials, so a draw that bypasses the recorded generators fails.
- **Outputs** (``_draws``). The subject, session and trial offsets, the gain, the onset time, the EMG
  noise and the angle noise are read back from the trials.

``scripts/mutate_synthetic_source.py`` breaks ``sources.py`` in many ways and shows which test fails
for which mutant; every test here fails for at least one.
"""

from functools import lru_cache

import numpy as np
import pytest

from emg_onset_lead import MinJerkSource, SampleSize, SyntheticSource

N_SEEDS, N_SUBJ, N_SESS, N_TRIALS = 64, 16, 3, 10
SHORT = dict(duration_s=0.01, fs_emg=200.0, fs_imu=150.0)  # 2 EMG and 2 angle samples
NOISY = dict(duration_s=0.1, fs_emg=200.0, fs_imu=150.0)  # 20 EMG and 15 angle samples, all before the onset
SOURCES = [SyntheticSource, MinJerkSource]
# Consecutive seeds (a tree rooted at seed + subject), and seeds equal modulo 2**32 (a truncated seed).
STREAM_SEEDS = [0, 1, 2, 3, 7, 2**32 - 1, 2**32, 2**32 + 1, 2**40]


class _Recorder(np.random.Generator):
    """A Generator that knows its starting state and how many values it handed out."""

    def __init__(self, bit_generator):
        super().__init__(bit_generator)
        self.start = bit_generator.state
        self.n_values = 0

    def normal(self, *args, **kwargs):
        out = super().normal(*args, **kwargs)
        self.n_values += np.size(out)
        return out

    def uniform(self, *args, **kwargs):
        out = super().uniform(*args, **kwargs)
        self.n_values += np.size(out)
        return out

    def words(self) -> list[int]:
        """The raw 64-bit words between the starting state and the present one."""
        bg = type(self.bit_generator)()
        bg.state = self.start
        end, out = self.bit_generator.state["state"], []
        while bg.state["state"] != end:
            out.append(bg.random_raw())
            assert len(out) < 10_000, "stream did not end where expected"
        return out


def _record(cls, seed, monkeypatch):
    """The recorded generators and the trials of one 16x3x10 study."""
    made, real = [], np.random.default_rng

    def recording_rng(s=None):
        made.append(_Recorder(real(s).bit_generator))
        return made[-1]

    with monkeypatch.context() as m:
        m.setattr(np.random, "default_rng", recording_rng)
        trials = list(cls(n_subjects=N_SUBJ, n_sessions=N_SESS, n_trials=N_TRIALS, seed=seed, **SHORT))
    return made, trials


def _study(cls, seed, shape, sd=(8.0, 6.0, 10.0), **signal):
    src = cls(n_subjects=shape[0], n_sessions=shape[1], n_trials=shape[2], true_lead_ms=0.0, seed=seed,
              subject_sd_ms=sd[0], session_sd_ms=sd[1], trial_sd_ms=sd[2], **signal)
    return src, list(src)


def _offsets(cls, seed, *sd):
    """True lead minus 0 for every trial, shape (subjects, sessions, trials)."""
    shape = (N_SUBJ, N_SESS, N_TRIALS)
    return np.array([t.meta["true_lead_ms"] for t in _study(cls, seed, shape, sd, **SHORT)[1]]).reshape(shape)


def _draws(cls, seed, shape=(N_SUBJ, N_SESS, N_TRIALS)) -> dict:
    """The random draws behind the outputs, standardised. Signals end before the onset, so the EMG is
    the bare noise (amplitude 1) and the angle is 20 deg plus noise. Arrays lead with ``shape``."""
    src, trials = _study(cls, seed, shape, **NOISY)
    lead = np.array([t.meta["true_lead_ms"] for t in trials])
    onset = np.array([t.meta["t_move_s"] for t in trials]) - lead / 1000.0
    lo, hi = src.emg_onset_range_s
    d = dict(lead=lead, gain=np.array([t.meta["emg_gain"] for t in trials]), onset=(onset - lo) / (hi - lo),
             emg=np.array([t.emg for t in trials]),
             angle=(np.array([t.angle_or_imu for t in trials]) - 20.0) / src.angle_noise_deg)
    return {k: v.reshape(*shape, *v.shape[1:]) for k, v in d.items()}


@lru_cache(maxsize=None)
def _components(cls):
    """Subject (S,), session (S, 3) and trial (S, 3, 10) offsets pooled over seeds, design SDs 8/6/10."""
    subj = np.concatenate([_offsets(cls, k, 8.0, 0.0, 0.0)[:, 0, 0] for k in range(N_SEEDS)])
    sess = np.concatenate([_offsets(cls, k, 0.0, 6.0, 0.0)[:, :, 0] for k in range(N_SEEDS)])
    trial = np.concatenate([_offsets(cls, k, 0.0, 0.0, 10.0) for k in range(N_SEEDS)])
    return subj, sess, trial


def test_variance_components_match_design():
    """Design checks on 1024 subjects. Tolerances are about 4 SE of each estimate."""
    subj, sess, trial = _components(SyntheticSource)
    mean_sd = (subj[:, None, None] + sess[:, :, None] + trial).mean(axis=(1, 2)).std(ddof=1)
    want = SampleSize.design_subject_sd(8.0, 6.0, 10.0, N_SESS, N_TRIALS)
    bad = []
    for name, got, expect, tol in (("subject SD", subj.std(ddof=1), 8.0, 0.7),
                                   ("session SD", sess.std(ddof=1), 6.0, 0.31),
                                   ("trial SD", trial.std(ddof=1), 10.0, 0.17),
                                   ("subject-mean SD", mean_sd, want, 0.79)):
        if abs(got - expect) > tol:
            bad.append(f"{name} {got:.3f} vs {expect:.3f}")
    pairs = [(f"subject, session {j}", subj, sess[:, j]) for j in range(N_SESS)]
    pairs += [(f"session {i}, session {j}", sess[:, i], sess[:, j]) for i, j in ((0, 1), (0, 2), (1, 2))]
    for name, a, b in pairs:
        c = float(np.corrcoef(a, b)[0, 1])
        if abs(c) > 0.125:
            bad.append(f"corr({name}) = {c:+.3f}")
    assert bad == []


@pytest.mark.parametrize("cls", SOURCES)
def test_gain_onset_and_noise_match_design(cls):
    """Gain ~ U(0.7, 1.3) per session, onset ~ U(1.5, 2.5) s, unit EMG noise, 0.3 deg angle noise.
    16 seeds, 768 sessions. Tolerances are about 4 SE."""
    ds = [_draws(cls, k) for k in range(16)]
    gain, onset, emg, angle = (np.concatenate([d[k] for d in ds]) for k in ("gain", "onset", "emg", "angle"))
    sess = np.concatenate([_offsets(cls, k, 0.0, 6.0, 0.0)[:, :, 0] for k in range(16)])
    g = gain[:, :, 0]
    assert (gain == g[:, :, None]).all(), "gain is per session"
    u = (g.ravel() - 0.7) / 0.6
    assert 0.0 <= u.min() and u.max() < 1.0
    sd_u = 1 / np.sqrt(12)
    for name, x, mean, sd, tol_mean, tol_sd in (
            ("gain", u, 0.5, sd_u, 4 * sd_u / np.sqrt(u.size), 0.02),
            ("onset", onset.ravel(), 0.5, sd_u, 4 * sd_u / np.sqrt(onset.size), 0.006),
            ("EMG noise", emg.ravel(), 0.0, 1.0, 4 / np.sqrt(emg.size), 4 / np.sqrt(2 * emg.size)),
            ("angle noise", angle.ravel(), 0.0, 1.0, 4 / np.sqrt(angle.size), 4 / np.sqrt(2 * angle.size))):
        assert abs(x.mean() - mean) < tol_mean, f"{name} mean {x.mean():.4f}"
        assert abs(x.std(ddof=1) - sd) < tol_sd, f"{name} SD {x.std(ddof=1):.4f}"
    tol = 4 / np.sqrt(g.shape[0])
    for i, j in ((0, 1), (0, 2), (1, 2)):
        assert abs(np.corrcoef(g[:, i], g[:, j])[0, 1]) < tol, f"corr(gain session {i}, {j})"
    assert abs(np.corrcoef(u, sess.ravel())[0, 1]) < 4 / np.sqrt(u.size), "corr(gain, session offset)"
    n = angle.shape[-1]
    c = np.corrcoef(emg[..., :n].ravel(), angle.ravel())[0, 1]
    assert abs(c) < 4 / np.sqrt(angle.size), f"corr(EMG noise, angle noise) = {c:+.3f}"


@pytest.mark.parametrize("cls", SOURCES)
def test_no_stream_overlaps_another_within_or_across_seeds(cls, monkeypatch):
    """No raw word is consumed twice by any generator of any seed, at any offset."""
    owner, clashes = {}, []
    for seed in STREAM_SEEDS:
        for i, rec in enumerate(_record(cls, seed, monkeypatch)[0]):
            for pos, word in enumerate(rec.words()):
                first = owner.setdefault(word, (seed, i, pos))
                if first != (seed, i, pos):
                    clashes.append(f"seed {seed} stream {i} word {pos} = seed/stream/word {first}")
    assert not clashes, f"{len(clashes)} shared words, e.g. {clashes[:3]}"


@pytest.mark.parametrize("cls", SOURCES)
def test_every_draw_comes_from_a_recorded_stream(cls, monkeypatch):
    """The overlap check sees every draw: 1 per subject, 2 per session (offset, gain), and per trial 2
    (jitter, onset) plus one per EMG and angle sample."""
    made, trials = _record(cls, 0, monkeypatch)
    n_samples = sum(t.emg.size + t.angle_or_imu.size for t in trials)
    want = N_SUBJ + 2 * N_SUBJ * N_SESS + 2 * N_SUBJ * N_SESS * N_TRIALS + n_samples
    assert sum(r.n_values for r in made) == want


@pytest.mark.parametrize("cls", SOURCES)
def test_no_output_value_repeats(cls):
    """Over seeds 0-3, no normal or uniform draw seen in the outputs appears twice. Chance is < 0.001:
    about 7e4 values, compared to 12 decimals."""
    normal, uniform = [], []
    for seed in range(4):
        d = _draws(cls, seed)
        normal += [_offsets(cls, seed, 1.0, 0.0, 0.0)[:, 0, 0], _offsets(cls, seed, 0.0, 1.0, 0.0)[:, :, 0],
                   _offsets(cls, seed, 0.0, 0.0, 1.0), d["emg"], d["angle"]]
        uniform += [(d["gain"][:, :, 0] - 0.7) / 0.6, d["onset"]]
    for name, parts in (("normal", normal), ("uniform", uniform)):
        x = np.round(np.concatenate([p.ravel() for p in parts]), 12)
        assert x.size == np.unique(x).size, f"{x.size - np.unique(x).size} repeated {name} draws"


@pytest.mark.parametrize("cls", SOURCES)
def test_trials_do_not_depend_on_study_size(cls):
    """Subject i, session j, trial k are the same draws (offsets, gain, onset, noise) at any study size."""
    small, big = _draws(cls, 4, (2, 2, 3)), _draws(cls, 4)
    for k in small:
        np.testing.assert_array_equal(small[k], big[k][:2, :2, :3], err_msg=k)
