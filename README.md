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

The `refine_s` wait gives the fit more of the rise. It cuts the spread, not the bias: at 150 ms / SNR 3
`RampOnset`'s per-trial SD is about 19-20 ms without it (19.2 and 20.0 ms; `refine_s` = 0.5 ms, one sample
after the hold) and about 14-15 ms with the 0.1 s wait (13.8 and 15.0 ms), with a mean error within 1.6 ms
either way (synthetic, 200 trials on each of 2 seeds; `scripts/refine_wait_sd.py`). `refine_s` = 0 is not
"`RampOnset` without the wait": the detector then falls back to the line fit at the hold, which reads
about 36 ms late there. Waiting 0.15 or 0.2 s did not help further (SD 14.9-16.1 ms, same script).

At the commit sample the trigger is checked:
- The envelope is still above threshold: a burst is on, and the fit runs to the commit.
- The envelope is back below threshold, but it peaked at >= `short_burst_peak` = 1.6 thresholds
  (baseline-removed) since the hold: an **isolated short burst** (40-80 ms) that has ended. It is
  kept, and the fit window ends at the envelope peak, so post-burst baseline does not pull the
  ramp model.
- Otherwise it is a baseline noise excursion that met the hold and decayed. That trigger is
  dropped and the detector **re-arms**: it rescans the samples it already has, so a later real
  burst is still found.

Why 1.6 (synthetic; Fact Checker's factor sweep, PR #4 round 2, 499 h of white noise):
- At the baseline noise level, false onsets fall from 0.016/h at 1.5x to 0.008/h at 1.6x. Noise
  triggers peaked at up to ~1.73x, so no factor is a hard bound.
- The cost is short-burst detection at 40 ms / SNR 3: 246 -> 234 of 300. Every other short-burst
  cell stays at 297-300/300.
- With a noise floor 10-20% above the baseline the floor dominates, not the factor. At +10%:
  0.60/h at 1.5x -> 0.28/h at 1.6x. At +20%: 63/h -> 47/h. About 67% of the 47/h at 1.6x
  (31.5/h) comes from the pre-existing path (envelope still above threshold at the commit), which
  this factor cannot fix.

The envelope peaks of real bursts do not depend on the factor. Re-measured at this head with
confirmation forced (60 trials per cell, `default_rng(1)` restarted per cell, synthetic): the
lowest peak of a 20 ms-rise isolated burst at SNR >= 5 was 2.696x (SNR 5 / 40 ms). That is the
measured minimum of one stream, not a bound. At SNR 3 / 40 ms the median peak is 2.3x and the
minimum 1.2x. Without any commit check, FC's
generator produced about one onset in noise per ~3000 trials, hundreds of ms early. The
previous check (envelope above threshold at the commit) dropped every isolated 40-60 ms burst.
Isolated bursts at 1.6x (synthetic; 20 ms linear rise, test generator and seed,
`tests/test_emg_short_burst.py`), 300 trials per cell: detected / 300, mean onset error (SD):

| SNR | 40 ms | 60 ms | 80 ms |
|---:|---|---|---|
| 3 | 243/300, -1.19 ms (9.8) | 299/300, -1.91 ms (8.4) | 300/300, -2.39 ms (8.5) |
| 5 | 300/300, +0.29 ms (4.2) | 300/300, +0.27 ms (4.0) | 300/300, -0.56 ms (4.0) |
| 15 | 300/300, -0.81 ms (1.6) | 300/300, -0.36 ms (1.5) | 300/300, -0.76 ms (1.8) |
| 30 | 300/300, -0.68 ms (0.9) | 300/300, -0.50 ms (1.0) | 300/300, -1.09 ms (1.2) |

At 60 trials (FC's short-burst form, 40-150 ms) every cell is 60/60 except SNR 3 / 40 ms at
49/60 (base 608a151: 58/60). The mean onset error is within +/-3 ms in every cell.

The 65% detection bound for 40 ms / SNR 3 in `tests/test_emg_short_burst.py` is a loose
regression guard. It catches a drop to about 60% or below, not small changes. Measured over test
seeds 0-9: 107-121/150 (synthetic).

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
  uses a 40 ms smoothstep, so leads from it read about 3.9-4.0 ms short (3.927 ms over 16x3x10, 4.00 ms
  over the smaller 8x3x5; see below). The 80/150 ms
  cases and the sign bias on `SyntheticSource` are strict xfails.
- A burst that has ended by the commit must have peaked at >= 1.6x threshold, or it is treated
  as noise. Isolated 40 ms bursts at SNR 3 are missed about 22% of the time (FC: 234/300 kept).
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
- **bias = -3.927 ms** (lead reads short). The seed-to-seed SD is 0.070 ms, and the 95% CI over
  6 seeds is [-4.000, -3.854] ms.
- **u = 0.135 ms** (range 0.097-0.166 over seeds). u_new_subject = 0.637 ms.
- Guard band w = k*u: 0.135 / 0.270 / 0.405 ms for k = 1/2/3, so a gate of 68.805 / 68.940 /
  69.075 ms. With u_new_subject: 69.307 / 69.944 / 70.580 ms.

The 4.00 ms in `tests/test_harness.py` is the same bias at 8x3x5, seeds 0-5: -3.996 ms, 95% CI
[-4.258, -3.734]. Those 720 trials are a subset of the 2880 in the 16x3x10 study. They differ from
the other 2160 (2159 detected) by -0.092 ms [-0.469, +0.284] (paired over 6 seeds, p = 0.556), so no
difference is shown: 3.9 and 4.00 are one bias read at two sample sizes. The interval of the
difference is 0.75 ms wide, so a difference as small as 0.07 ms is not excluded either.

One of the 2880 trials is missed and left out of every 16x3x10 number here, including those of
`bias_guard_band.py`: seed 1, S11-D3 T02, rejected by the EMG baseline-rest check. It is not one of the
8x3x5 trials. Its subject mean uses 29 of 30 trials. Leaving out its whole session or its whole
subject instead moves the bias over seeds from -3.927 to -3.918 or -3.923 ms (synthetic;
`scripts/synthetic_source_spread.py`, part 5).

These replace the stale base-608a151 values (bias +0.35 ms, u = 0.17 ms), and the values from before
the `SyntheticSource` seed fix (bias -3.93 ms, 95% CI [-4.05, -3.81], u = 0.139 ms).

Uncertainty: u and the seed spread cover only detector noise on this one synthetic burst shape.
They do not cover shape or rise time, which dominate. On smoothstep bursts at SNR 15 the bias
runs from about -0.3 ms (10 ms rise) to -14.5 ms (150 ms rise); Fact Checker PR #4, table B.
The real rise-time distribution is unknown (no recorded EMG). Option C only raises the gate
(68.67 + w). The -3.927 ms lateness already uses about 4 of the 6 ms margin, so any w > 0 cuts
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
That command is a **single-seed demo**. Its subject-mean SD is one draw from 16 subjects. Over seeds
0-5 that SD is 6.89-9.53 ms for the true lead and 6.86-9.39 ms for the detected lead (synthetic), against a
design SD of 8.91 ms (`SampleSize.design_subject_sd` with subject/session/trial SDs 8/6/10 ms,
3 sessions, 10 trials). Pass `--seeds 6` to pool seeds 0..5. Writes `out/trials.csv` and
`out/summary.json`.

`SyntheticSource` seeds each subject, session and trial from its own node of a `SeedSequence` spawn
tree. Until that fix, the subject used `default_rng([seed, s])` and session 0 used
`default_rng([seed, s, 0])`. numpy pads a seed list that converts to fewer 32-bit words (not
list items) than its 4-word pool with zeros, so those were the same stream: session 0's offset was
0.75x the subject offset. Subject-mean true-lead SD, 16x3x10 (synthetic; `scripts/synthetic_source_spread.py`, `--pool-seeds 400` for the second row):

| pooled seeds | subjects | old seeding (95% CI) | now (95% CI) |
|---|---:|---|---|
| 0-39 | 640 | 10.55 ms [10.01, 11.17] | 9.14 ms [8.66, 9.67] |
| 0-399 | 6400 | 10.59 ms [10.41, 10.77] | 8.96 ms [8.81, 9.12] |

In theory it is 10.55 ms with the old seeding and 8.91 ms with independent streams.
`tests/test_synthetic_source_seeds.py` guards this: no raw word of any stream is used twice, at any
offset, within a seed or across the 12 seeds it checks (neighbouring seeds, strides 8 and 64, seeds equal
modulo 2**32 or 10**6; reuse with seeds outside that list is not checked), and the gain, onset and noise
draws are checked too.
`scripts/mutate_synthetic_source.py` breaks the source and shows that each of those tests can fail.

**Seeded outputs from before #7 do not reproduce.** Every `SyntheticSource` and `MinJerkSource` trial is a
new draw at the same seed. Example: the true-lead subject-mean SD at seed 0 (16x3x10) went from 9.40 to
8.76 ms (synthetic; both regenerated by `scripts/synthetic_source_spread.py`, part 2). Do not compare
`trials.csv` or `summary.json` from before and after #7 trial by trial.

The latency budget is 18.67 ms IMU path + 50 ms EMG extra = 68.67 ms. `--imu-path-ms` is the only
way to change it. `--total-latency-ms` is a check only. It must equal the IMU path + 50 ms, and a
mismatch is an error. It used to set the IMU path to total - 50 without saying so.

Sample rates are parameters. The synthetic defaults are EMG 2000 Hz and angle/IMU 150 Hz; `--fs-emg`
and `--fs-imu` override them (written to `summary.json` as `demo.fs_emg_hz` and `demo.fs_imu_hz`).
The public datasets checked so far use other rates: Lucchetti et al. 2025
(doi:10.1038/s41597-025-06174-3) has EMG at 1000 Hz and markers at 125 Hz (`EmgFreq`/`KinFreq` in
HS02.mat), and Rojas-Martinez et al. 2020 (doi:10.1038/s41597-020-00717-6) has HD-EMG at 2048 Hz. On
synthetic data the rates do not move the mean lead error: paired against 2000/150 Hz, the shift at
1000/125, 2048/125, 2000/125 and 1000/150 Hz is at most 0.17 ms and every 95% interval contains 0
(16 seeds x 72 trials, `SyntheticSource` and `MinJerkSource`, no misses). The per-trial SD changes by
-1.0% to +6.6%, and the known smoothstep bias (about -4 ms) is unchanged. This check is at SNR 15
only (synthetic; `scripts/dataset_rates_check.py`, part 1, with the current `SeedSequence` seeding).
`tests/test_dataset_rates.py` guards 1000/125 and 2048/125. The detectors take `fs` from each `Trial`.
The EMG band's upper edge is capped at 0.45 * fs, so it stays 20-450 Hz down to 1000 Hz. Below that
the band narrows and the leads can be badly wrong, with no error: on a 2x1x3 run (seed 0, mean true
lead 91.6 ms) the mean measured lead is 88.3 ms at 2000 Hz, 111.7 ms at 100 Hz, 213.1 ms at 50 Hz
(1 of 6 trials missed) and 328.2 ms at 45 Hz (synthetic; same script, part 2). Nothing below 1000 Hz
is tested for accuracy. A rate too low for the default EMG band (below about 44.4 Hz) or the IMU smoothing
(20 Hz or less) raises a `ValueError`; in the CLI it is an argparse error (exit 2).

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
is a seeded cluster bootstrap (reported only). See `PLAN.md`. The CLI adds `demo.single_seed`
and the sample rates `demo.fs_emg_hz` / `demo.fs_imu_hz`.

## Data
Only synthetic data is used. `SyntheticSource` is a critically damped 60 deg step (nonzero
acceleration at onset). `MinJerkSource` is the Flash and Hogan minimum-jerk path (duration and
amplitude set per trial). Meta EMG datasets are excluded for licensing; no external data has been
cleared yet. To add a source, implement `TrialSource` (`__iter__` yielding `Trial`).
