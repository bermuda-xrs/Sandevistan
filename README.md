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
`hold_s` (30 ms). No `filtfilt`, no zero-phase, no centered windows. A trial that moves inside the 1.0 s
baseline (the excursion of the post-settle window exceeds 12 early-baseline SDs, or the onset
would fall inside it) is rejected. Tests check that truncating at the commit sample gives the
same onset, that truncating one sample earlier gives none, that corrupting future samples changes
nothing, and that chunked == one-shot processing.

## What `lead_ms` is
`lead_ms` is the **signal-onset lead**: the angle onset minus the EMG onset.

### EMG onset
The 8 Hz envelope threshold and 30 ms hold only decide that a burst is underway. The onset is
estimated separately, because the old estimate (a line through the envelope at the hold) moved
with the burst. The envelope filter's delay and a fixed k*SD threshold both depend on how fast
and how high the burst rises. On Fact Checker's synthetic grid the old lead error ran from +12.3 ms
(10 ms rise, SNR 30) to -36.1 ms (150 ms rise, SNR 3).

`RampOnset` (in `detect.py`) works on the band-passed signal, not the envelope, so no envelope
delay enters. After the hold the detector waits `refine_s` = 0.1 s (samples already seen only)
and takes the window from `lookback_s` = 0.4 s before the hold (never before the end of the 1.0 s
baseline) to the commit sample. It scales that window by the baseline SD and models it as:

- before the onset `j`: N(0, 1)
- from `j` on: N(0, a^2), with a = 1 + (G - 1) * min((t - j) / R, 1). This is a linear amplitude
  ramp of rise time R to gain G; R = 0 is a step.

For every R (0 to 0.5 s in 20 ms steps, uniform prior) and G (1.2 to 60, 20 log-spaced values),
the Gaussian log-likelihood of every `j` is one FFT correlation. Energy is pooled in 2 ms blocks
for speed. The onset is the **posterior mean** of `j` (flat prior over the window). The band-pass
filter's energy delay (1.15 ms at 20-450 Hz, order 2, 2000 Hz) is then subtracted; it depends on
the filter only. If the fit cannot run, the old line fit is the fallback.

Why this removes the rise-time and SNR dependence: the likelihood uses the burst's own gain and
rise, so it does not depend on where a fixed absolute threshold sits on the burst. It also
reads the raw band-passed samples, so the 8 Hz filter delay is not involved. Two choices matter
on the slow, low-SNR cells (synthetic):

- The posterior mean, not the arg-max. The arg-max snaps late when onset and rise time trade off.
- A uniform prior on R. A log-spaced grid weights short rises and read about +2.5 ms late at
  80-150 ms rise, SNR 3.

The `refine_s` wait gives the fit more of the rise. At 150 ms / SNR 3 the per-trial SD is 24 ms
without the wait and 14 ms with it (synthetic, 200 trials). A wait beyond 0.1 s did not help further.

At the commit sample the trigger is checked:
- The envelope is still above threshold: a burst is on, and the fit runs to the commit.
- The envelope is back below threshold, but it peaked at >= `short_burst_peak` = 1.5 thresholds
  (baseline-removed) since the hold: an **isolated short burst** (40-80 ms) that has ended. It is
  kept, and the fit window ends at the envelope peak, so post-burst baseline does not pull the
  ramp model.
- Otherwise it is a baseline noise excursion that met the hold and decayed. That trigger is
  dropped and the detector **re-arms**: it rescans the samples it already has, so a later real
  burst is still found.

Why 1.5 (synthetic): in ~490 h of white noise, 229 false triggers peaked at 1.00-1.57x
threshold, and 3 of them (1.3%) reached 1.5x. Isolated bursts with a 20 ms rise peak above
2.7x at SNR >= 5. At SNR 3 / 40 ms the median peak is 2.3x and the minimum 1.2x. Without any commit check, FC's
generator produced about one onset in noise per ~3000 trials, hundreds of ms early. The
previous check (envelope above threshold at the commit) dropped every isolated 40-60 ms burst.
Isolated bursts detected, 20 ms rise, 60 trials (synthetic): 60/60 at SNR 3/5/15/30 x
40/60/80/100/150 ms, except SNR 3 / 40 ms at 53/60 (base 608a151: 58/60). Mean onset error is
within +/-3 ms in every cell (`tests/test_emg_short_burst.py`).

Results (synthetic only; full table in the PR): on FC's generator (linear amplitude ramp), rise
10/20/40/80/150 ms x SNR 3/5/15/30, every cell's mean EMG onset error is within about 3.3 ms (max |mean| 2.23 ms on the test's
seed; Fact Checker measured -3.32 ms at 10 ms / SNR 3, n 200, on another seed). Its 95% CI
is within +/-6 ms at 200 trials per cell (SNR 3-5) and 80 (SNR 15-30); `tests/test_emg_grid.py`.
Per-trial SD is 1-3.5 ms at SNR 15-30 but 5-14 ms at SNR 3, rising with rise time (about 14.2-14.5 ms at
150 ms / SNR 3; one 40-trial draw showed 16.6). At low SNR the
mean needs many trials: FC's 40-trial run leaves the 150 ms / SNR 3 cell's CI at about +/-5 ms.

Limitations (synthetic):
- The model is a linear amplitude ramp. An S-curve (smoothstep) burst reads late by about 10% of
  its rise: +3.7 ms at 40 ms, about +8 ms at 80 ms, +14 ms at 150 ms (SNR 15). `SyntheticSource`
  uses a 40 ms smoothstep, so leads from it read about 3.9 ms short (see below). The 80/150 ms
  cases and the sign bias on `SyntheticSource` are strict xfails.
- A burst that has ended by the commit must have peaked at >= 1.5x threshold, or it is treated
  as noise. Isolated 40 ms bursts at SNR 3 are missed about 13-22% of the time.
- The EMG onset now commits 0.1 s after the hold. `detected_s` is that commit. `trigger_s` is the
  hold time.

### In-family bias and the option C guard band (synthetic)
`scripts/bias_guard_band.py` (about 40 s) runs the planned 16x3x10 `SyntheticSource` study for
seeds 0-5. That is the repo's in-family case: a 40 ms smoothstep burst at SNR 15 x U(0.7, 1.3)
and a critically damped move. For each seed it computes the lead error (detected minus true)
and a nested ANOVA (subject > session > trial). It reports:
- bias = the mean of subject means
- u = the subject-cluster SE of that bias
- u_new_subject

This is the same method as `bias_split.py` on analysis/emg-polarity-bias (3e9c27b).

Values at this checkout (synthetic):
- **bias = -3.93 ms** (lead reads short). The seed-to-seed SD is 0.12 ms, and the 95% CI over
  6 seeds is [-4.05, -3.81] ms.
- **u = 0.139 ms** (range 0.087-0.180 over seeds). u_new_subject = 0.649 ms.
- Guard band w = k*u: 0.14 / 0.28 / 0.42 ms for k = 1/2/3, so a gate of 68.81 / 68.95 /
  69.09 ms. With u_new_subject: 69.32 / 69.97 / 70.62 ms.

These replace the stale base-608a151 values (bias +0.35 ms, u = 0.17 ms).

Uncertainty: u and the seed spread cover only detector noise on this one synthetic burst shape.
They do not cover shape or rise time, which dominate. On smoothstep bursts at SNR 15 the bias
runs from about -0.3 ms (10 ms rise) to -14.5 ms (150 ms rise); Fact Checker PR #4, table B.
The real rise-time distribution is unknown (no recorded EMG). Option C only raises the gate
(68.67 + w). The -3.93 ms lateness already uses about 4 of the 6 ms margin, so any w > 0 cuts
power further. If option C instead added an estimated bias back, a wrong estimate would be
anti-conservative.

The angle onset is not that line. A minimum-jerk flexion (Flash and Hogan 1985) leaves rest with
zero acceleration, so velocity starts out curved and a short line reads the onset late (on the order
of +15 to +90 ms as the move gets slower, smaller, or noisier). After the velocity hold, the IMU
detector waits `refine_s` (0.55 s) — still only samples already seen — and least-squares matches two
rest-to-move families on that fixed window: minimum-jerk of free duration, and a critically damped
step of free time constant. The better residual is the onset. The old critically damped generator is
one of the two families, so it is not a separate tuning. If the fit does not beat a flat baseline,
the line-fit fallback is used.

Signal-onset lead is the quantity compared to the exo latency because the 68.67 ms budget already
includes the EMG path's processing (`t_win/2`). `detect_lead_ms` (difference of commit times; the EMG
commit is now 0.1 s after the hold) is a diagnostic and is not gated. `emg_confirmation_delay_ms`
is the EMG hold time (`trigger_s`) minus the back-dated EMG onset: about 47 ms on the synthetic
burst, not the 30 ms hold. A triggered controller fires at the hold, so the 0.1 s refine wait is
not in it. It is always reported. It is added to the budget only when `LatencyConfig(triggered=True)`.

The pass/fail gate is all of: the one-sided lower bound on subject means above the budget, the run
having the planned 16 subjects x 3 sessions (`PLANNED_SHAPE` in `stats.py`), and no more than 10% of trials
missed. A 3x1 run and the old 8x3 plan both fail. A run that misses 268 of 480 trials fails even if the
detected trials look fine. The pilot is 16x3 since 2026-11-03: 8x3 was rejected at about 53% power (computed
from synthetic variances, subject/session/trial SD 8/6/10 ms, not pilot data), and the
adopted power calculation assumed 3 sessions per subject.

Limitation: the kinematic onset was checked on synthetic minimum-jerk and critically damped moves
only. It is **unverified on real motion** until it is compared to a recorded bench move with a known
true onset. Do not treat a pass near the threshold as settled. The 16x3 plan is likewise unverified
until real pilot variances exist (`SampleSize` takes those variances; it does not invent them).

Status: **Fact Checker round 2 failed (2026-11-11).** On out-of-family motion the onset error ran from
49 ms early to 271 ms late, and the EMG onset showed a sign bias of about +/-11 ms, which is over the
6 ms margin. The EMG side is addressed on synthetic data by the ramp-onset estimator above. The
out-of-family motion cases are unchanged and are kept as strict xfails
(`tests/test_imu_out_of_family.py`): late velocity peak +142.8 ms, slow creep then fast move
+271.3 ms, cosine ramp -42.3 ms, linear ramp -49.3 ms. The detectors here are not validated for real
EMG or real motion, and the synthetic tests do not show otherwise.

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

The latency budget is 18.67 ms IMU path + 50 ms EMG extra = 68.67 ms. `--imu-path-ms` is the only
way to change it. `--total-latency-ms` is a check only. It must equal the IMU path + 50 ms, and a
mismatch is an error. It used to set the IMU path to total - 50 without saying so.

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
