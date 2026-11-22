# Live EMG onset detectors (research only, SYNTHETIC data)

Nothing in `src/` imports this and no production or protected value changes.

- `detectors.py`: causal alarm statistics (envelope variants, TKEO, RMS, CUSUM, GLR, BOCPD).
- `test_detectors.py`: brute-force references, causality, and identity of `env8_h30` at h = 6 with the production
  trigger. Not in the default `testpaths`; run `PYTHONPATH=src python -m pytest research/causal_detectors`.
- `run_sim.py`: rest noise, latency trials, artefact stress test, report (commands in its docstring).
- `results/`: `digest.md` (short), `tables.md` (all cells), `report.json` (every number), `equiv.txt`.

Raw per-batch output (571 MB) is not committed. It is seeded; one batch (400) was regenerated with the committed
code and matched the stored file bit for bit. The other 719 batches were not re-checked.
