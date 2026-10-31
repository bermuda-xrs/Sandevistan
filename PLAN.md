# Plan: emg-onset-lead v1

Measures lead = IMU/angle onset - EMG onset per trial; tests pooled one-sided lower bound > exo latency.

## Modules (src/emg_onset_lead/)
- `latency.py` `LatencyConfig(imu_path_ms=18.67, emg_extra_ms=50.0)`; `.total_ms` (=68.67),
  `.required_ms(path="total"|"imu")`. Needed: single source of truth for the pass criterion.
- `detect.py` causal streaming detectors. `OnsetParams(baseline_s, settle_s, k, hold_s)`,
  `EmgParams`, `ImuParams(mode="velocity"|"angle")`, `OnsetResult(onset_s, detected_s)`.
  `EmgOnsetDetector` / `ImuOnsetDetector` (subclasses of `_StreamDetector`): `.process(chunk)`;
  `.run(x, fs, params, chunk=None)` static. Both use the same chunked path: `sosfilt` with carried
  `zi`, then one shared baseline mean+k*SD threshold + minimum-hold state machine.
  Needed: EMG and IMU share window/filter/threshold logic, so causality is tested once.
- `sources.py` `Trial` dataclass, `TrialSource` Protocol (`__iter__` -> Trial), `SyntheticSource`
  (seeded, known EMD, subject/session/trial variation, EMG 2000 Hz, angle 150 Hz).
- `stats.py` `LeadStats`: `t_lower_bound`, `boot_lower_bound`, `summarize(results, latency, conf)`.
- `harness.py` `TrialResult`, `Harness.run(source, ...)`, `write_csv`, `write_summary`.
- `cli.py` / `__main__.py`: argparse `run --synthetic`.

## Reuse (not reimplemented)
scipy.signal `butter(output="sos")`, `sosfilt`, `sosfilt_zi`; `scipy.stats.t.ppf`;
`numpy.random.default_rng`; csv/json from the stdlib.
## Causality and reported quantity
No filtfilt/zero-phase/centered windows. Baseline = first `baseline_s` (past only). A run above
mean+k*SD held `hold_s` is *confirmed* at `detected_s`; at that instant only, the onset is back-dated:
line fit over the last `fit_s` (70 ms, same for EMG and angle) of the standardised feature,
extrapolated to baseline. `lead_ms` = signal-onset lead (back-dated imu - emg). Why: raw threshold
crossings lag by path-dependent amounts (~+20 ms EMG, ~+30 ms angle => ~+10 ms lead bias), and the
68.67 ms budget already includes the EMG path's detection time, so detection-time lead would double
count it. `detect_lead_ms` (confirmation difference) is a CSV diagnostic only, not gated.
Truncating at confirmation gives the same onset; one sample earlier gives none.

## Statistics (nesting subject > session > trial)
- Per session: mean of trial leads, t lower bound over trials.
- Per subject: mean of session means, t lower bound over session means (trials if 1 session).
- Pooled **primary**: mean of subject means (equal subject weight), one-sided t lower bound on the
  subject means (df = n_subjects-1). Pass = lower bound > required latency (68.67 ms default).
- Pooled **secondary** (reported, not gating): seeded cluster bootstrap (resample subjects, then
  sessions within subject; 5th percentile). Percentile bootstrap undercovers with ~8 clusters,
  hence the t-bound gates.
- Planned study: 8 subjects x 3 sessions over >= 3 days, electrodes re-donned each session.
  Summary records actual vs planned shape.
