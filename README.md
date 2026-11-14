# emg-onset-lead

Dataset-agnostic harness for an EMG-controlled assistive elbow exoskeleton. It measures how far
surface-EMG movement onset **leads** joint-angle/IMU onset, per trial, and tests whether that lead
exceeds the exo's intent-to-torque latency.

- `lead_ms = imu_onset - emg_onset`
- `margin_ms = lead_ms - required latency` (default 68.67 ms total = 18.67 ms IMU path + 50 ms EMG extra)
- The 50 ms term is **t_win/2**, half of the controller's 100 ms EMG window.
- The default gate is the **proportional** controller. `--triggered` (off unless you pass it) adds the measured EMG confirmation delay on top of 68.67 ms.

## Causality guarantee
EMG and angle/IMU go through the same chunked pipeline: causal IIR filters (`scipy.signal.sosfilt`
with carried state), a fixed past-only baseline window (first `baseline_s` = **1.0 s**; the limb
is required to be at rest in that window), and a mean + k*SD threshold that must be held for
`hold_s`. No `filtfilt`, no zero-phase, no centered windows. A trial that moves inside the 1.0 s
baseline (the excursion of the post-settle window exceeds 12 early-baseline SDs, or the onset
would fall inside it) is rejected. Tests check that truncating at the commit sample gives the
same onset, that truncating one sample earlier gives none, that corrupting future samples changes
nothing, and that chunked == one-shot processing.

## What `lead_ms` is
`lead_ms` is the **signal-onset lead**: the angle onset minus the EMG onset.

EMG still line-fits the last `fit_s` (70 ms) of the baseline-standardised envelope and extrapolates
to baseline. That line is computed at the hold, from samples at or before the hold only. A
slope whose rise across the window is below one part in a million of the feature level and of
the threshold (a negative slope, or the ~1e-17 slope `polyfit` can return on a flat window)
is not extrapolated; the detector then reports the threshold-crossing sample. The back-dated
index stays inside the fit window and does not move earlier than the end of the 1.0 s baseline.

The angle onset is not that line. A minimum-jerk flexion (Flash and Hogan 1985) leaves rest with
zero acceleration, so velocity starts out curved and a short line reads the onset late (on the order
of +15 to +90 ms as the move gets slower, smaller, or noisier). After the velocity hold, the IMU
detector waits `refine_s` (0.55 s) — still only samples already seen — and least-squares matches two
rest-to-move families on that fixed window: minimum-jerk of free duration, and a critically damped
step of free time constant. The better residual is the onset. The old critically damped generator is
one of the two families, so it is not a separate tuning. If the fit does not beat a flat baseline,
the line-fit fallback is used.

Signal-onset lead is the quantity compared to the exo latency because the 68.67 ms budget already
includes the EMG path's processing (`t_win/2`). `detect_lead_ms` (difference of commit times) is a
diagnostic and is not gated. `emg_confirmation_delay_ms` is the EMG commit time minus the back-dated
EMG onset (about 51 ms on the synthetic burst, not the 30 ms hold). It is always reported. It is
added to the budget only when `LatencyConfig(triggered=True)`.

The pass/fail gate is all of: the one-sided lower bound on subject means above the budget, the run
having the planned 16 subjects x 3 sessions (`PLANNED_SHAPE` in `stats.py`), and no more than 10% of trials
missed. A 3x1 run and the old 8x3 plan both fail. A run that misses 268 of 480 trials fails even if the
detected trials look fine. The pilot is 16x3 since 2026-11-03: 8x3 was rejected at about 53% power, and the
adopted power calculation assumed 3 sessions per subject.

Limitation: the kinematic onset was checked on synthetic minimum-jerk and critically damped moves
only. It is **unverified on real motion** until it is compared to a recorded bench move with a known
true onset. Do not treat a pass near the threshold as settled. The 16x3 plan is likewise unverified
until real pilot variances exist (`SampleSize` takes those variances; it does not invent them).

Status: **Fact Checker round 2 failed (2026-11-11).** On out-of-family motion the onset error ran from
49 ms early to 271 ms late, and the EMG onset showed a sign bias of about +/-11 ms, which is over the
6 ms margin. The detectors here are not validated for real EMG or real motion, and the synthetic
bias tests do not show otherwise.

## Run
```
pip install -e .[test]
pytest
python -m emg_onset_lead run --synthetic --subjects 16 --sessions 3 --seed 0 --out out/
```
That command is a **single-seed demo**. Seed 0's subject-mean spread can be a couple of milliseconds
against a design SD of about 8.91 ms (`SampleSize.design_subject_sd` with subject/session/trial SDs
8/6/10 ms, 3 sessions, 10 trials). Pass `--seeds 6` to pool seeds 0..5. Writes `out/trials.csv` and
`out/summary.json`.

## Output
`trials.csv`: `subject_id, session_id, trial_id, emg_onset_s, imu_onset_s, lead_ms, margin_ms,
true_lead_ms, detect_lead_ms, emg_confirmation_delay_ms, baseline_rest_ok`
(`true_lead_ms` is NaN unless the source provides `meta["true_lead_ms"]`). A rejected baseline
has `baseline_rest_ok` false and a NaN lead.

`summary.json`: `latency` (including `triggered` and whether the 50 ms term is t_win/2),
`confirmation_delay_ms` (mean, min, max), `study` (planned 16x3 over >= 3 days, electrodes
re-donned; `plan_status` is `unverified` until real pilot variance exists; actual shape),
`gates` (lower bound, planned shape, missed-trial fraction), `sessions[]`, `subjects[]`, `pooled`,
`pass`. `pass` requires the lower bound, the 16x3 shape, and a miss fraction at or below 10%.
Pooled is the mean of subject means with a t-bound on subject means; `bootstrap_lower_bound_ms`
is a seeded cluster bootstrap (reported only). See `PLAN.md`. The CLI adds `demo.single_seed`.

## Data
Only synthetic data is used. `SyntheticSource` is a critically damped 60 deg step (nonzero
acceleration at onset). `MinJerkSource` is the Flash and Hogan minimum-jerk path (duration and
amplitude set per trial). Meta EMG datasets are excluded for licensing; no external data has been
cleared yet. To add a source, implement `TrialSource` (`__iter__` yielding `Trial`).
