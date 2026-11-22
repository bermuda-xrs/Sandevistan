# Plan: emg-onset-lead v1
Measures lead = IMU/angle onset - EMG onset per trial; tests pooled one-sided lower bound > exo latency.

## Modules (src/emg_onset_lead/)
- `latency.py` `LatencyConfig(imu_path_ms=18.67, emg_extra_ms=50.0, triggered=False)`.
  `emg_extra_ms` is t_win/2 (half of the 100 ms EMG window). `.total_ms` = 68.67.
  `.required_ms(path, confirmation_delay_ms)` adds the delay only when `triggered` and path is total.
- `detect.py` causal streaming detectors. `OnsetParams(baseline_s=1.0, settle_s, k, hold_s, fit_s, rest_ratio)`,
  `EmgParams`, `ImuParams(mode, refine_s=0.55)`, `OnsetResult(onset_s, detected_s, baseline_rest_ok)`.
  `OnsetResult.trigger_s` = hold time. `EmgParams(refine_s=0.1, lookback_s=0.4, rise_max_s=0.5,
  rise_step_s=0.02, gain_range=(1.2, 60), n_gain=20, block_s=0.002, short_burst_peak=1.6)`.
  EMG: envelope k*SD + 30 ms hold triggers. 0.1 s later the onset is `RampOnset`: the posterior mean of
  a change point on the band-passed signal (baseline-SD units). Before it: N(0,1). After it: a linear
  amplitude ramp, rise R uniform 0-0.5 s, gain log 1.2-60. Minus the band-pass energy delay (1.15 ms).
  At the commit the envelope must still be above threshold, or (an ended short burst) must have
  peaked at >= 1.6x threshold since the hold, in which case the fit ends at that peak. Otherwise
  the trigger is dropped and the detector re-arms on the stored samples. The line fit over `fit_s` is the fallback. IMU: after the hold, fit
  minimum-jerk and critically damped families on the raw angle (`KinematicOnset`); line fit is the fallback.
  `BaselineRest` rejects a move inside the 1.0 s baseline.
- `sources.py` `Trial`, `TrialSource`, `SyntheticSource` (critically damped), `MinJerkSource` (Flash & Hogan 1985).
  Streams are a `SeedSequence` spawn tree (subject > session > trial); spread numbers: `scripts/synthetic_source_spread.py`.
- `stats.py` `LeadStats` (t bound, cluster bootstrap, summarize) and `SampleSize` (power / required N from
  pilot between- and within-subject variance). 16x3 is `plan_status: unverified` until real variances exist.
- `harness.py` `TrialResult` (includes confirmation delay and baseline flag), `Harness.run`.
- `cli.py`: `run --synthetic`. `--imu-path-ms` sets the budget; `--total-latency-ms` must match it or errors. Default output is a single seed. `--seeds N` pools 0..N-1. `--triggered` opts in.
## Reuse (not reimplemented)
scipy.signal `butter(sos)`, `sosfilt`, `sosfilt_zi`; `scipy.stats.t` and `nct`; numpy `default_rng`, `polyfit`, `correlate`;
scipy.fft `rfft`/`irfft`/`next_fast_len` (RampOnset likelihoods).
## Causality
No filtfilt / zero-phase / centered windows. Baseline = first 1.0 s, and it must be at rest. EMG commits 0.1 s
after the hold, IMU 0.55 s after, each using only samples up to that index. A dropped EMG trigger is rescanned
from stored samples only. Truncating at commit matches;
one sample earlier does not.
## Gate
Pass iff lower bound > required (68.67 proportional, or 68.67 + measured EMG confirmation delay if triggered),
the run is 16 subjects x 3 sessions (`PLANNED_SHAPE`; 8x3 was rejected at ~53% power, 2026-11-03), and
missed trials are at most 10% (268/480 fails). Confirmation delay (EMG hold time minus onset) is
reported either way (synthetic bursts land near 47 ms, not the 30 ms hold).
## Statistics
Subject > session > trial. Pooled gate: mean of subject means, one-sided t (df = n_subjects-1). Cluster
bootstrap is reported, not gated. `SampleSize.power` / `required_subjects` take pilot variances, margin, alpha.
## Known limits
Kinematic onset is synthetic-only. Unverified on real motion until a bench move with a known true onset.
16x3 power is unverified until those pilot variances are real.
EMG onset (synthetic): FC's rise 10-150 ms x SNR 3-30 grid holds every cell's mean and 95% CI inside +/-6 ms
at 200/80 trials per cell (`tests/test_emg_grid.py`); at FC's 40 trials the 150 ms / SNR 3 CI reaches -7.3 ms. An S-curve burst reads late by ~10% of its rise (strict xfails). Per-trial SD is
5-14 ms at SNR 3. Isolated 40-80 ms bursts: >= 96% detected except 40 ms / SNR 3 (71-81% over test seeds 0-9, 243/300 at 300 trials; synthetic);
`tests/test_emg_short_burst.py`. Its 65% bound for 40 ms / SNR 3 is a loose regression guard: it catches a drop to about
60% or below, not small changes (measured 107-121/150 over seeds 0-9; synthetic).
In-family bias (synthetic, `scripts/bias_guard_band.py`, 16x3x10 SyntheticSource, seeds 0-5): -3.927 ms
(95% CI over seeds [-4.000, -3.854]); option C u = 0.135 ms, w = k*u -> gate 68.805/68.940/69.075 ms (k = 1/2/3).
u covers detector noise on one burst shape only. Rise/shape dependence (about -0.3 to -14.5 ms) dominates
and is unknown without recorded EMG. Out-of-family motion is unchanged (strict xfails in `tests/test_imu_out_of_family.py`).
Fact Checker round 2 FAILED 2026-11-11: out-of-family motion 49 ms early to 271 ms late, EMG sign bias ~±11 ms (> 6 ms margin). A single demo seed is not the 8.91 ms design SD.
