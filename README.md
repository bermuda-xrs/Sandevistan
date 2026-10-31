# emg-onset-lead

Dataset-agnostic harness for an EMG-controlled assistive elbow exoskeleton. It measures how far
surface-EMG movement onset **leads** joint-angle/IMU onset, per trial, and tests whether that lead
exceeds the exo's intent-to-torque latency.

- `lead_ms = imu_onset - emg_onset`
- `margin_ms = lead_ms - required latency` (default 68.67 ms total = 18.67 ms IMU path + 50 ms EMG extra)

## Causality guarantee
EMG and angle/IMU go through the same chunked pipeline: causal IIR filters (`scipy.signal.sosfilt`
with carried state), a fixed past-only baseline window (first `baseline_s`), and a
mean + k*SD threshold that must be held for `hold_s`. No `filtfilt`, no zero-phase, no centered
windows. Onset is stamped at the first sample of the sustained run; it is *confirmed* `hold_s`
later. Tests check that truncating at the confirmation sample gives the same onset, that
truncating one sample earlier gives none, that corrupting future samples changes nothing, and that
chunked == one-shot processing.

Known limitation: causal filter delay biases each onset late. On the synthetic data the measured
lead exceeds the true lead by about +10 ms (EMG-envelope and angle-smoothing delays differ), so
a pass near the latency threshold should be read with that bias in mind.

## Run
```
pip install -e .[test]
pytest
python -m emg_onset_lead run --synthetic --subjects 8 --sessions 3 --seed 0 --out out/
```
Writes `out/trials.csv` and `out/summary.json`.

## Output
`trials.csv`: `subject_id, session_id, trial_id, emg_onset_s, imu_onset_s, lead_ms, margin_ms,
true_lead_ms` (`true_lead_ms` is NaN unless the source provides `meta["true_lead_ms"]`).

`summary.json`: `latency`, `study` (planned 8 subjects x 3 sessions over >= 3 days, electrodes
re-donned; actual shape), `sessions[]`, `subjects[]`, `pooled`, `pass`. Each level has
`mean_lead_ms`, one-sided `lower_bound_ms` (default 95%), `required_ms`, `margin_ms`, `pass`
(lower bound > required latency). Pooled is the mean of subject means with a t-bound on subject
means (gating); `bootstrap_lower_bound_ms` is a seeded cluster bootstrap over subjects/sessions
(reported only). See `PLAN.md`.

## Data
Only synthetic data is used. Meta EMG datasets are excluded for licensing; no external data has
been cleared yet. To add a source, implement `TrialSource` (`__iter__` yielding `Trial`).
