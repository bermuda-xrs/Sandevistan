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

## What `lead_ms` is
`lead_ms` is the **signal-onset lead**: the back-dated angle onset minus the back-dated EMG onset.
A bare threshold crossing lags the true onset by a path-dependent amount (about +20 ms for the EMG
envelope, +30 ms for angle velocity on the synthetic data), which inflated the lead by ~10 ms. At
confirmation time both paths therefore apply the same rule over the same live window (`fit_s` = 70 ms,
ending at the confirmation sample): fit a line to the baseline-standardised feature and extrapolate
to the baseline level. This is still causal (nothing after confirmation is used). Signal-onset lead is
the quantity comparable to the exo latency because the 68.67 ms budget already includes the EMG path's
own processing/detection delay; detection-time lead would count that twice. `detect_lead_ms` (difference
of confirmation times) is written as a diagnostic and is not gated.

Limitation: `fit_s` was chosen on the synthetic generator so that both paths' onset errors are near
zero (synthetic recovery: pooled lead error within +/-3 ms, mixed sign across seeds). Real EMG/angle
rise shapes differ, so the residual path bias on real data is unverified; check it before trusting a
pass near the threshold. Absolute onsets are only meaningful to a few ms; the lead is the target.

## Run
```
pip install -e .[test]
pytest
python -m emg_onset_lead run --synthetic --subjects 8 --sessions 3 --seed 0 --out out/
```
Writes `out/trials.csv` and `out/summary.json`.

## Output
`trials.csv`: `subject_id, session_id, trial_id, emg_onset_s, imu_onset_s, lead_ms, margin_ms,
true_lead_ms, detect_lead_ms` (`true_lead_ms` is NaN unless the source provides
`meta["true_lead_ms"]`).

`summary.json`: `latency`, `study` (planned 8 subjects x 3 sessions over >= 3 days, electrodes
re-donned; actual shape), `sessions[]`, `subjects[]`, `pooled`, `pass`. Each level has
`mean_lead_ms`, one-sided `lower_bound_ms` (default 95%), `required_ms`, `margin_ms`, `pass`
(lower bound > required latency). Pooled is the mean of subject means with a t-bound on subject
means (gating); `bootstrap_lower_bound_ms` is a seeded cluster bootstrap over subjects/sessions
(reported only). See `PLAN.md`.

## Data
Only synthetic data is used. Meta EMG datasets are excluded for licensing; no external data has
been cleared yet. To add a source, implement `TrialSource` (`__iter__` yielding `Trial`).
