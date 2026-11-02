# Plan: emg-onset-lead v1
Measures lead = IMU/angle onset - EMG onset per trial; tests pooled one-sided lower bound > exo latency.

## Modules (src/emg_onset_lead/)
- `latency.py` `LatencyConfig(imu_path_ms=18.67, emg_extra_ms=50.0, triggered=False)`.
  `emg_extra_ms` is t_win/2 (half of the 100 ms EMG window). `.total_ms` = 68.67.
  `.required_ms(path, confirmation_delay_ms)` adds the delay only when `triggered` and path is total.
- `detect.py` causal streaming detectors. `OnsetParams(baseline_s=1.0, settle_s, k, hold_s, fit_s, rest_ratio)`,
  `EmgParams`, `ImuParams(mode, refine_s=0.55)`, `OnsetResult(onset_s, detected_s, baseline_rest_ok)`.
  EMG: line fit over `fit_s`. A rise below 1e-6 of the feature level and the threshold
  (so slope<=0, and a ~1e-17 fit on a flat window) falls back to the threshold crossing.
  The onset stays inside the fit window and at or after the 1.0 s baseline. IMU: after the hold, fit
  minimum-jerk and critically damped families on the raw angle (`KinematicOnset`); line fit is the fallback.
  `BaselineRest` rejects a move inside the 1.0 s baseline.
- `sources.py` `Trial`, `TrialSource`, `SyntheticSource` (critically damped), `MinJerkSource` (Flash & Hogan 1985).
- `stats.py` `LeadStats` (t bound, cluster bootstrap, summarize) and `SampleSize` (power / required N from
  pilot between- and within-subject variance). 8x3 is `plan_status: unverified` until real variances exist.
- `harness.py` `TrialResult` (includes confirmation delay and baseline flag), `Harness.run`.
- `cli.py`: `run --synthetic`. Default output is a single seed. `--seeds N` pools 0..N-1. `--triggered` opts in.
## Reuse (not reimplemented)
scipy.signal `butter(sos)`, `sosfilt`, `sosfilt_zi`; `scipy.stats.t` and `nct`; numpy `default_rng`, `polyfit`, `correlate`.
## Causality
No filtfilt / zero-phase / centered windows. Baseline = first 1.0 s, and it must be at rest. EMG commits at the
hold. IMU commits `refine_s` later, still using only samples up to that index. Truncating at commit matches;
one sample earlier does not.
## Gate
Pass iff lower bound > required (68.67 proportional, or 68.67 + measured EMG confirmation delay if triggered),
the run is 8 subjects x 3 sessions, and missed trials are at most 10% (134/240 fails). Confirmation delay is
reported either way (synthetic bursts land near 51 ms, not the 30 ms hold).
## Statistics
Subject > session > trial. Pooled gate: mean of subject means, one-sided t (df = n_subjects-1). Cluster
bootstrap is reported, not gated. `SampleSize.power` / `required_subjects` take pilot variances, margin, alpha.
## Known limits
Kinematic onset is synthetic-only. Unverified on real motion until a bench move with a known true onset.
8x3 power is unverified until those pilot variances are real. A single demo seed is not the 8.91 ms design SD.
