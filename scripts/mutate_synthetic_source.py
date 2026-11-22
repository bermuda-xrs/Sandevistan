"""Mutation check of tests/test_synthetic_source_seeds.py. SYNTHETIC DATA ONLY.

Copies ``src`` and the test file to a temporary directory, breaks ``SyntheticSource`` in one way per
mutant, runs the tests and prints which fail. Exits 1 if the unmutated copy fails, if a mutant gets
through, or if a test never fails. M0 is the seeding before #7. M5 and M7 reuse a stream with an
offset; M3, M8 and M18 reuse streams across seeds; M9-M15, M17 and M18 touch only the gain, onset or noise
draws.

Usage (repo root):  python scripts/mutate_synthetic_source.py [--out out.md]     (a few minutes, 1 process)
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEST = "test_synthetic_source_seeds.py"

LEGACY = """\
        for s in range(self.n_subjects):
            subj_off = np.random.default_rng([self.seed, s]).normal(0, self.subject_sd_ms)
            for ss in range(self.n_sessions):
                rng_s = np.random.default_rng([self.seed, s, ss])
                sess_off = rng_s.normal(0, self.session_sd_ms)
                gain = rng_s.uniform(0.7, 1.3)
                for t in range(self.n_trials):
                    rng = np.random.default_rng([self.seed, s, ss, t + 1])
                    yield self._trial(rng, s, ss, t, subj_off + sess_off, gain)
"""
FLAT = """\
        leaves = np.random.SeedSequence(self.seed).spawn(
            self.n_subjects * (1 + self.n_sessions * (1 + self.n_trials)))
        k = 0
        for s in range(self.n_subjects):
            subj_off = np.random.default_rng(leaves[k]).normal(0, self.subject_sd_ms)
            k += 1
            for ss in range(self.n_sessions):
                rng_s = np.random.default_rng(leaves[k])
                k += 1
                sess_off = rng_s.normal(0, self.session_sd_ms)
                gain = rng_s.uniform(0.7, 1.3)
                for t in range(self.n_trials):
                    rng = np.random.default_rng(leaves[k])
                    k += 1
                    yield self._trial(rng, s, ss, t, subj_off + sess_off, gain)
"""
# Anchors in sources.py. Every replacement must match, so a changed source breaks the script loudly.
ROOT_SS = "np.random.SeedSequence(self.seed).spawn(self.n_subjects)"
SUBJ = "subj_off = np.random.default_rng(subj_seed).normal(0, self.subject_sd_ms)"
RNG_S = "rng_s = np.random.default_rng(sess_seed)"
GAIN = "gain = rng_s.uniform(0.7, 1.3)"
RNG = "rng = np.random.default_rng(trial_seed)"
ONSET = "t_emg = float(rng.uniform(*self.emg_onset_range_s))"
EMG = "emg = amp * rng.normal(0, 1, te.size)"
ANGLE_CALL = "angle = self._angle(ti, t_move, rng)"
CD_NOISE = "return angle + rng.normal(0, self.angle_noise_deg, ti.size)"
MJ_NOISE = "return 20.0 + y + rng.normal(0, self.angle_noise_deg, ti.size)"
PAD = " " * 20

MUTANTS: dict[str, list[tuple[str, str]] | str] = {
    "M0 seeding before #7 (session 0 reads the subject stream)": LEGACY,
    "M1 subject reads session 0's stream": [
        (SUBJ, "subj_off = np.random.default_rng(np.random.SeedSequence("
               "self.seed, spawn_key=(s, 1, 0))).normal(0, self.subject_sd_ms)")],
    "M2 trial list seed with a trailing zero": LEGACY.replace("[self.seed, s])", "[self.seed, s, 99])")
                                                     .replace("t + 1]", "t]"),
    "M3 subject tree rooted at seed + subject (cross-seed reuse)": [
        (ROOT_SS, "[np.random.SeedSequence(self.seed + k) for k in range(self.n_subjects)]")],
    "M4 trial streams shared across sessions": [
        (RNG, "rng = np.random.default_rng(np.random.SeedSequence(self.seed, spawn_key=(s, 1, 1 + t)))")],
    "M5 session 0 reads the subject stream from its 2nd draw (one-draw offset)": [
        (RNG_S, "rng_s = np.random.default_rng(sess_seed if ss else subj_seed)\n                if ss == 0:\n"
                + PAD + "rng_s.normal()")],
    "M5b as M5, skipping with advance(1) instead of a draw": [
        (RNG_S, "rng_s = np.random.default_rng(sess_seed if ss else subj_seed)\n                if ss == 0:\n"
                + PAD + "rng_s.bit_generator.advance(1)")],
    "M6 trials 1-9 read trial 0's stream": [(RNG, "rng = np.random.default_rng(trial_seeds[0])")],
    "M7 trial t reads trial 0's stream from draw t (one-draw offsets)": [
        (RNG, "rng = np.random.default_rng(trial_seeds[0])\n" + PAD + "rng.bit_generator.advance(t)")],
    "M8 seed truncated to 32 bits (cross-seed reuse)": [
        (ROOT_SS, ROOT_SS.replace("self.seed", "self.seed % 2**32"))],
    "M9 gain drawn once per subject": [
        (SUBJ, SUBJ + "\n            g0 = np.random.default_rng(subj_seed.spawn(1)[0]).uniform(0.7, 1.3)"),
        (GAIN, "gain = g0")],
    "M10 gain reads the session-offset stream again": [
        (GAIN, "gain = np.random.default_rng(sess_seed).uniform(0.7, 1.3)")],
    "M11 gain range 0.8-1.2": [(GAIN, "gain = rng_s.uniform(0.8, 1.2)")],
    "M12 EMG noise the same in every trial": [(EMG, EMG.replace("rng.", "np.random.default_rng(0)."))],
    "M13 angle noise reuses the EMG noise (copied generator)": [
        (EMG, "import copy\n        rng2 = copy.deepcopy(rng)\n        " + EMG),
        (ANGLE_CALL, "angle = self._angle(ti, t_move, rng2)")],
    "M14 angle noise the same in every trial (critically damped only)": [
        (CD_NOISE, CD_NOISE.replace("rng.", "np.random.default_rng(0)."))],
    "M15 angle noise the same in every trial (minimum jerk only)": [
        (MJ_NOISE, MJ_NOISE.replace("rng.", "np.random.default_rng(0)."))],
    "M16 flat spawn: streams depend on the study size": FLAT,
    "M17 EMG noise SD 1.05": [(EMG, EMG.replace("(0, 1,", "(0, 1.05,"))],
    "M18 onset time ignores the seed (cross-seed reuse)": [
        (ONSET, ONSET.replace("rng.", "np.random.default_rng([s, ss, t, 9])."))],
}
HEAD = "    def __iter__(self) -> Iterator[Trial]:\n"
TAIL = "    def _trial(self, rng, s, ss, t, offset_ms, gain) -> Trial:\n"


def mutate(text: str, change) -> str:
    if isinstance(change, str):
        i, j = text.index(HEAD) + len(HEAD), text.index(TAIL)
        return text[:i] + change + "\n" + text[j:]
    for old, new in change:
        assert text.count(old) == 1, f"anchor matches {text.count(old)} times: {old}"
        text = text.replace(old, new)
    return text


def run(change) -> tuple[list[str], list[str]]:
    """(all tests, failing tests) with ``change`` applied to a copy of src."""
    with tempfile.TemporaryDirectory() as d:
        shutil.copytree(ROOT / "src", Path(d) / "src", ignore=shutil.ignore_patterns("__pycache__"))
        shutil.copy(ROOT / "tests" / TEST, d)
        p = Path(d) / "src" / "emg_onset_lead" / "sources.py"
        if change is not None:
            p.write_text(mutate(p.read_text(), change))
        env = dict(os.environ, PYTHONPATH=str(Path(d) / "src"), PYTHONDONTWRITEBYTECODE="1")
        where = subprocess.run([sys.executable, "-c", "import emg_onset_lead as e; print(e.__file__)"], cwd=d,
                               env=env, capture_output=True, text=True, check=True).stdout
        assert Path(where.strip()).resolve().is_relative_to(Path(d).resolve()), f"wrong package: {where}"
        r = subprocess.run([sys.executable, "-m", "pytest", "-v", "--tb=no", "-p", "no:cacheprovider", TEST],
                           cwd=d, env=env, capture_output=True, text=True)
    res = re.findall(r"::(\S+) (PASSED|FAILED)", r.stdout)
    assert res, r.stdout + r.stderr
    return [n for n, _ in res], [n for n, v in res if v == "FAILED"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    tests, base = run(None)
    kills = {name: run(change)[1] for name, change in MUTANTS.items()}
    L = [f"# Mutation check of tests/{TEST} (SYNTHETIC)", "",
         f"Unmutated: {len(tests) - len(base)} of {len(tests)} pass.", "",
         "| Mutant of sources.py | Failing tests |", "|---|---|"]
    for name, f in kills.items():
        L.append(f"| {name} | {len(f)}: " + ", ".join(x.removeprefix("test_") for x in f) + " |")
    L += ["", "| Test | Mutants it fails for |", "|---|---|"]
    by_test = {t: [n.split()[0] for n, f in kills.items() if t in f] for t in tests}
    L += [f"| {t} | {len(m)}: {', '.join(m)} |" for t, m in by_test.items()]
    survivors, idle = [n for n, f in kills.items() if not f], [t for t, m in by_test.items() if not m]
    L += ["", f"Mutants that get through: {survivors or 'none'}. Tests that never fail: {idle or 'none'}."]
    text = "\n".join(L) + "\n"
    print(text)
    if a.out:
        a.out.write_text(text)
    return 1 if base or survivors or idle else 0


if __name__ == "__main__":
    raise SystemExit(main())
