"""EMG onset error against the wait after the hold (``refine_s``). SYNTHETIC DATA ONLY.

README "EMG onset": what the 0.1 s wait buys at 150 ms rise, SNR 3. No recorded EMG is read. Trials are
Fact Checker's linear-ramp generator (``tests/fc_synthetic.py``), 200 per generator seed, seeded
``default_rng([gen, 150, 3, i])`` as in FC's reconcile_24.py. The same trials are used for every wait.

``refine_s`` = 0 is not "RampOnset without the wait": the detector then uses the line-fit fallback at the
hold (``EmgOnsetDetector``, ``if self._refine_n <= 0``). RampOnset without a wait is ``refine_s`` = 0.5 ms,
one sample at 2000 Hz. Columns: onset error (detected minus true) mean and SD per generator seed. Production
is ``refine_s`` = 100 ms; nothing is changed here.

Usage (repo root):  python scripts/refine_wait_sd.py [--gens 1 2] [--out out.md]     (1 process)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src"), str(ROOT / "tests")]

from fc_synthetic import fc_trial  # noqa: E402

from emg_onset_lead import EmgOnsetDetector, EmgParams  # noqa: E402

RISE_MS, SNR, N = 150, 3, 200
WAITS_S = (0.0, 0.0005, 0.01, 0.025, 0.05, 0.1, 0.15, 0.2)


def errors(gen: int, refine_s: float) -> np.ndarray:
    out = []
    for i in range(N):
        t_on, _, emg, _, fs, _ = fc_trial(np.random.default_rng([gen, RISE_MS, SNR, i]), 75.0, RISE_MS / 1e3,
                                          float(SNR))
        r = EmgOnsetDetector.run(emg, fs, EmgParams(refine_s=refine_s), chunk=100)
        out.append(np.nan if r.onset_s is None else (r.onset_s - t_on) * 1e3)
    return np.array(out)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gens", type=int, nargs="+", default=[1, 2])
    ap.add_argument("--out", type=Path)
    a = ap.parse_args(argv)
    assert EmgParams().refine_s == 0.1, "production wait"
    L = [f"## EMG onset error vs refine_s (SYNTHETIC; linear ramp {RISE_MS} ms, SNR {SNR}, {N} trials per "
         f"generator seed, numpy {np.__version__})", "",
         "| refine_s (ms) | estimator | " + " | ".join(f"gen {g}: n, mean (SD) ms" for g in a.gens) + " |",
         "|---:|---|" + "---|" * len(a.gens)]
    for w in WAITS_S:
        cells = []
        for g in a.gens:
            e = errors(g, w)
            ok = e[np.isfinite(e)]
            cells.append(f"{ok.size}, {ok.mean():+.2f} ({ok.std(ddof=1):.2f})")
        est = "line fit at the hold" if w == 0 else "RampOnset"
        L.append(f"| {w * 1e3:g} | {est} | " + " | ".join(cells) + " |")
    text = "\n".join(L) + "\n"
    print(text)
    if a.out:
        a.out.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
