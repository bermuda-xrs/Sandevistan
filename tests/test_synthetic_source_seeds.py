"""SyntheticSource random streams: design SDs 8/6/10 ms and no shared streams. Synthetic data only.

The old source seeded the subject with ``default_rng([seed, s])`` and the session with
``default_rng([seed, s, ss])``. numpy's SeedSequence pads a seed list shorter than its 4-word pool
with zeros, so ``[seed, s]`` and ``[seed, s, 0]`` are the same stream: session 0's offset was
0.75 x (= 6/8) the subject offset, and the subject-mean true-lead SD was sqrt((8 + 6/3)^2 + 2 * 36/9
+ 100/30) = 10.55 ms instead of 8.91. ``test_legacy_seeding_fails_these_checks`` keeps that scheme as a
negative control: the same checks must reject it.

The checks read only ``Trial.meta["true_lead_ms"]``. Setting two of the three SDs to 0 isolates one
component, because the draws themselves do not depend on the SDs.
"""

from functools import lru_cache

import numpy as np
import pytest

from emg_onset_lead import SampleSize, SyntheticSource

N_SEEDS, N_SUBJ, N_SESS, N_TRIALS = 64, 16, 3, 10
SHORT = dict(duration_s=0.01, fs_emg=200.0, fs_imu=150.0)  # signals are not read here; keep them tiny


class LegacySource(SyntheticSource):
    """The old seeding, bug kept on purpose."""

    def __iter__(self):
        for s in range(self.n_subjects):
            subj_off = np.random.default_rng([self.seed, s]).normal(0, self.subject_sd_ms)
            for ss in range(self.n_sessions):
                rng_s = np.random.default_rng([self.seed, s, ss])
                sess_off = rng_s.normal(0, self.session_sd_ms)
                gain = rng_s.uniform(0.7, 1.3)
                for t in range(self.n_trials):
                    rng = np.random.default_rng([self.seed, s, ss, t + 1])
                    yield self._trial(rng, s, ss, t, subj_off + sess_off, gain)


def _offsets(cls, seed, subject, session, trial, n_subjects=N_SUBJ, n_sessions=N_SESS, n_trials=N_TRIALS):
    """True lead minus 0 for every trial, shape (subjects, sessions, trials)."""
    src = cls(n_subjects=n_subjects, n_sessions=n_sessions, n_trials=n_trials, true_lead_ms=0.0, seed=seed,
              subject_sd_ms=subject, session_sd_ms=session, trial_sd_ms=trial, **SHORT)
    lead = np.array([t.meta["true_lead_ms"] for t in src])
    return lead.reshape(n_subjects, n_sessions, n_trials)


@lru_cache(maxsize=None)
def _components(cls):
    """Subject (S,), session (S, 3) and trial (S, 3, 10) offsets pooled over seeds, design SDs 8/6/10."""
    subj = np.concatenate([_offsets(cls, k, 8.0, 0.0, 0.0)[:, 0, 0] for k in range(N_SEEDS)])
    sess = np.concatenate([_offsets(cls, k, 0.0, 6.0, 0.0)[:, :, 0] for k in range(N_SEEDS)])
    trial = np.concatenate([_offsets(cls, k, 0.0, 0.0, 10.0) for k in range(N_SEEDS)])
    return subj, sess, trial


def _violations(cls) -> list[str]:
    """Design checks on 1024 subjects. Tolerances are about 4 SE of each estimate."""
    subj, sess, trial = _components(cls)
    mean_sd = (subj[:, None, None] + sess[:, :, None] + trial).mean(axis=(1, 2)).std(ddof=1)
    want = SampleSize.design_subject_sd(8.0, 6.0, 10.0, N_SESS, N_TRIALS)
    bad = []
    for name, got, expect, tol in (("subject SD", subj.std(ddof=1), 8.0, 0.7),
                                   ("session SD", sess.std(ddof=1), 6.0, 0.31),
                                   ("trial SD", trial.std(ddof=1), 10.0, 0.17),
                                   ("subject-mean SD", mean_sd, want, 0.79)):
        if abs(got - expect) > tol:
            bad.append(f"{name} {got:.3f} vs {expect:.3f}")
    for j in range(N_SESS):
        c = float(np.corrcoef(subj, sess[:, j])[0, 1])
        if abs(c) > 0.125:
            bad.append(f"corr(subject, session {j}) = {c:+.3f}")
    for i, j in ((0, 1), (0, 2), (1, 2)):
        c = float(np.corrcoef(sess[:, i], sess[:, j])[0, 1])
        if abs(c) > 0.125:
            bad.append(f"corr(session {i}, session {j}) = {c:+.3f}")
    return bad


def _shared_first_draws(cls, seed=0) -> int:
    """How many (subject, session, trial) first normal draws repeat another one, at one seed."""
    subj = _offsets(cls, seed, 1.0, 0.0, 0.0)[:, 0, 0]
    sess = _offsets(cls, seed, 0.0, 1.0, 0.0)[:, :, 0].ravel()
    trial = _offsets(cls, seed, 0.0, 0.0, 1.0).ravel()
    z = np.concatenate([subj, sess, trial])
    return int(z.size - np.unique(np.round(z, 9)).size)


def test_variance_components_match_design():
    assert _violations(SyntheticSource) == []


@pytest.mark.parametrize("seed", [0, 1, 7, 2**32 - 1, 2**40])
def test_no_stream_is_used_twice(seed):
    assert _shared_first_draws(SyntheticSource, seed) == 0


def test_legacy_seeding_fails_these_checks():
    bad = _violations(LegacySource)
    assert any(b.startswith("corr(subject, session 0)") for b in bad), bad
    assert any(b.startswith("subject-mean SD") for b in bad), bad
    assert _shared_first_draws(LegacySource) == N_SUBJ  # every subject stream is session 0's stream
    subj, sess, _ = _components(LegacySource)
    np.testing.assert_allclose(sess[:, 0], 0.75 * subj)


def test_numpy_pads_short_seed_lists_with_zeros():
    """The numpy behaviour behind the old collision: padding to the pool, not dropping zeros."""
    def same(a, b):
        states = [np.random.SeedSequence(x).generate_state(8) for x in (a, b)]
        return np.array_equal(*states)

    assert same([5, 3], [5, 3, 0]) and same([5, 3], [5, 3, 0, 0])  # up to 4 words: padded, equal
    assert not same([5, 3, 1, 2], [5, 3, 1, 2, 0])  # a 5th word is mixed in, even if 0
    assert not same([5, 3], [5, 3, 0, 0, 0])


def test_trials_do_not_depend_on_study_size():
    """Subject i, session j, trial k are the same draws whatever the study's size (as before)."""
    small = _offsets(SyntheticSource, 4, 8.0, 6.0, 10.0, n_subjects=2, n_sessions=2, n_trials=3)
    big = _offsets(SyntheticSource, 4, 8.0, 6.0, 10.0)
    np.testing.assert_array_equal(small, big[:2, :2, :3])
