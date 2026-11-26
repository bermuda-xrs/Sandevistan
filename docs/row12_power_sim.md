# Row-12 power simulation with fixed seeds (SYNTHETIC)

**All numbers here are synthetic.** No recorded EMG or motion was used, and nothing here validates the detectors on real signals. `src/` is unchanged (detectors, the 68.67 ms gate and `short_burst_peak` = 1.6).

Script: `scripts/row12_power_sim.py`. Test: `tests/test_row12_power_sim_seeds.py`. This replaces the earlier scratch results (16 subjects: 0.72 at 10 ms rise, 0.54 at 20 ms; "21-23 / 30-31 subjects"), which are **withdrawn**.

## What was wrong
The scratch generator seeded the subject with `default_rng([seed, i])` and the session with `default_rng([seed, i, j])`. numpy's SeedSequence pads a seed list with zeros if it converts to fewer 32-bit words than the pool holds (4 by default) ([`mix_entropy`](https://github.com/numpy/numpy/blob/v2.5.5/numpy/random/bit_generator.pyx#L363-L370)), so `[seed, i]` and `[seed, i, 0]` are the same seed. A list that already fills the pool is different: `[1, 2, 3, 4]` and `[1, 2, 3, 4, 0]` are not the same seed. The count is of 32-bit words, not list items: an entry of 2**32 or more takes two words, so `[2**40, 1, 2]` already fills the pool. So session 0 reused the subject stream: its offset was exactly 0.75 x the subject offset, and the SD of a subject mean was 10.55 ms (theory) instead of 8.91. Fact Checker (FC) found this. The earlier explanation for the gap to the analytic power ("detector noise") was wrong.

## Method
- Seeds are a SeedSequence spawn-key tree: study -> subject -> session -> trial.
- 250 studies per cell (ids 0-249). Each is simulated once at 48 subjects x 3 sessions x 10 trials for 10 and 20 ms rise, and at 16 subjects for 30, 40 and 80 ms. Smaller sizes are the first n subjects.
- **Cells are not independent.** Sizes are nested, and all rises and both true leads share offsets, EMG noise and angle noise. Chance in the offsets moves all cells together.
- Gate: `LeadStats.summarize` with `planned=(n, 3)`. Leads are combined as in `Harness.run`; a test checks they are equal.
- Correction: bias = -0.109 x rise_est + 0.51 ms, lead - bias. rise_est is the 10-90% rise of the trial's own envelope / 0.6084, never using the detected onset. The coefficients are taken as given and were fitted on synthetic smoothstep bursts.
- Burst shape is smoothstep only, SNR 15 x U(0.7, 1.3).
- Percent intervals are 95% Clopper-Pearson.

## Results
Studies scored: 250. ALL SYNTHETIC.

### Variance check (all simulated subjects)

| Component | Realised | Design |
|---|---|---|
| Subject SD (ms) | 7.953 | 8 |
| Session SD (ms) | 6.013 | 6 |
| Trial SD (ms) | 9.978 | 10 |
| Subject-mean true-lead SD, 3x10 (ms) | 8.871 | 8.907 |
| corr(subject offset, session j offset), j = 0, 1, 2 | -0.0011, +0.0077, +0.0006 | 0 |
| corr between session offsets (0-1, 0-2, 1-2) | -0.0064, -0.0021, +0.0016 | 0 |
| Subjects pooled | 12000 | |

### Pass rate vs subjects (x 3 sessions x 10 trials), true lead 74.67 ms, gate 68.67 ms

| Subjects | Trials/study | 10 ms rise | 20 ms rise |
|---|---|---|---|
| 16 | 480 | 209/250 = 0.836 (0.784–0.880) | 167/250 = 0.668 (0.606–0.726) |
| 20 | 600 | 223/250 = 0.892 (0.847–0.928) | 189/250 = 0.756 (0.698–0.808) |
| 24 | 720 | 235/250 = 0.940 (0.903–0.966) | 206/250 = 0.824 (0.771–0.869) |
| 28 | 840 | 239/250 = 0.956 (0.923–0.978) | 217/250 = 0.868 (0.820–0.907) |
| 32 | 960 | 242/250 = 0.968 (0.938–0.986) | 228/250 = 0.912 (0.870–0.944) |
| 36 | 1080 | 247/250 = 0.988 (0.965–0.998) | 228/250 = 0.912 (0.870–0.944) |
| 40 | 1200 | 248/250 = 0.992 (0.971–0.999) | 234/250 = 0.936 (0.898–0.963) |
| 44 | 1320 | 250/250 = 1.000 (0.985–1.000) | 239/250 = 0.956 (0.923–0.978) |
| 48 | 1440 | 250/250 = 1.000 (0.985–1.000) | 245/250 = 0.980 (0.954–0.993) |

### Subjects for 0.8 power: probit-in-sqrt(subjects) binomial ML fit over all grid points; 95% CI = 2000x bootstrap over studies (each resampled study keeps all its nested sizes)

| Rise | n at 0.8 (point) | 95% bootstrap CI | Trials/study at the point | Fit chi2 p |
|---|---|---|---|---|
| 10 ms | 14.7 | 12.4–16.8 | 442 | 0.86 (chi2 3.2, 7 df; nested points are correlated, so this is indicative) |
| 20 ms | 22.8 | 20.2–25.4 | 685 | 0.92 (chi2 2.6, 7 df; nested points are correlated, so this is indicative) |

### Baseline five-row table (no correction), 16x3x10, fixed seeds

| Rise | Power at 74.67 ms | False pass at 68.67 ms | Mean lead err (ms) |
|---|---|---|---|
| 10 ms | 209/250 = 0.836 (0.784–0.880) | 12/250 = 0.048 (0.025–0.082) | -0.25 |
| 20 ms | 167/250 = 0.668 (0.606–0.726) | 6/250 = 0.024 (0.009–0.052) | -1.38 |
| 30 ms | 100/250 = 0.400 (0.339–0.464) | 0/250 = 0.000 (0.000–0.015) | -3.01 |
| 40 ms | 59/250 = 0.236 (0.185–0.294) | 0/250 = 0.000 (0.000–0.015) | -3.95 |
| 80 ms | 3/250 = 0.012 (0.002–0.035) | 0/250 = 0.000 (0.000–0.015) | -7.96 |

### Power, true lead 74.67 ms (row 12), 16x3x10, correction primary = 40 Hz envelope

| Rise | Power uncorrected | Mean lead err (ms) | Power corrected | Mean corrected lead err (ms) | Rise-est err mean ± SD (ms) |
|---|---|---|---|---|---|
| 10 ms | 209/250 = 0.836 (0.784–0.880) | -0.25 | 245/250 = 0.980 (0.954–0.993) | +2.27 | +17.8 ± 15.9 |
| 20 ms | 167/250 = 0.668 (0.606–0.726) | -1.38 | 241/250 = 0.964 (0.933–0.983) | +1.64 | +12.4 ± 16.8 |
| 30 ms | 100/250 = 0.400 (0.339–0.464) | -3.01 | 236/250 = 0.944 (0.908–0.969) | +0.72 | +8.9 ± 17.7 |
| 40 ms | 59/250 = 0.236 (0.185–0.294) | -3.95 | 234/250 = 0.936 (0.898–0.963) | +0.62 | +6.6 ± 18.4 |
| 80 ms | 3/250 = 0.012 (0.002–0.035) | -7.96 | 217/250 = 0.868 (0.820–0.907) | +0.16 | -0.8 ± 21.4 |

### False pass, true lead exactly 68.67 ms (nominal 0.05), 16x3x10, correction primary = 40 Hz envelope

| Rise | False pass uncorrected | Mean lead err (ms) | False pass corrected | Mean corrected lead err (ms) | Rise-est err mean ± SD (ms) | One-sided binomial p (>0.05), uncorr / corr |
|---|---|---|---|---|---|---|
| 10 ms | 12/250 = 0.048 (0.025–0.082) | -0.25 | 69/250 = 0.276 (0.222–0.336) | +2.27 | +17.8 ± 15.9 | 0.6 / 9.5e-32 |
| 20 ms | 6/250 = 0.024 (0.009–0.052) | -1.38 | 43/250 = 0.172 (0.127–0.225) | +1.64 | +12.4 ± 16.8 | 0.99 / 1.7e-12 |
| 30 ms | 0/250 = 0.000 (0.000–0.015) | -3.01 | 23/250 = 0.092 (0.059–0.135) | +0.73 | +8.9 ± 17.7 | 1 / 0.0039 |
| 40 ms | 0/250 = 0.000 (0.000–0.015) | -3.95 | 22/250 = 0.088 (0.056–0.130) | +0.62 | +6.6 ± 18.4 | 1 / 0.0078 |
| 80 ms | 0/250 = 0.000 (0.000–0.015) | -7.96 | 15/250 = 0.060 (0.034–0.097) | +0.16 | -0.8 ± 21.4 | 1 / 0.27 |

### Cutoff sensitivity (16x3x10, same studies)

| Rise | Cutoff | Power corrected | False pass corrected | Rise-est err mean ± SD (ms) |
|---|---|---|---|---|
| 10 ms | 20 Hz | 249/250 = 0.996 (0.978–1.000) | 153/250 = 0.612 (0.549–0.673) | +37.3 ± 24.6 (failed 0) |
| 10 ms | 40 Hz | 245/250 = 0.980 (0.954–0.993) | 69/250 = 0.276 (0.222–0.336) | +17.8 ± 15.9 (failed 0) |
| 10 ms | 80 Hz | 241/250 = 0.964 (0.933–0.983) | 31/250 = 0.124 (0.086–0.171) | +8.3 ± 10.2 (failed 0) |
| 20 ms | 20 Hz | 249/250 = 0.996 (0.978–1.000) | 121/250 = 0.484 (0.421–0.548) | +30.0 ± 25.2 (failed 0) |
| 20 ms | 40 Hz | 241/250 = 0.964 (0.933–0.983) | 43/250 = 0.172 (0.127–0.225) | +12.4 ± 16.8 (failed 0) |
| 20 ms | 80 Hz | 236/250 = 0.944 (0.908–0.969) | 25/250 = 0.100 (0.066–0.144) | +4.3 ± 10.9 (failed 0) |
| 30 ms | 20 Hz | 246/250 = 0.984 (0.960–0.996) | 72/250 = 0.288 (0.233–0.348) | +24.3 ± 26.0 (failed 0) |
| 30 ms | 40 Hz | 236/250 = 0.944 (0.908–0.969) | 23/250 = 0.092 (0.059–0.135) | +8.9 ± 17.7 (failed 0) |
| 30 ms | 80 Hz | 213/250 = 0.852 (0.802–0.894) | 14/250 = 0.056 (0.031–0.092) | +1.4 ± 11.7 (failed 0) |
| 40 ms | 20 Hz | 245/250 = 0.980 (0.954–0.993) | 61/250 = 0.244 (0.192–0.302) | +20.0 ± 26.8 (failed 0) |
| 40 ms | 40 Hz | 234/250 = 0.936 (0.898–0.963) | 22/250 = 0.088 (0.056–0.130) | +6.6 ± 18.4 (failed 0) |
| 40 ms | 80 Hz | 208/250 = 0.832 (0.780–0.876) | 13/250 = 0.052 (0.028–0.087) | -1.6 ± 12.7 (failed 0) |
| 80 ms | 20 Hz | 241/250 = 0.964 (0.933–0.983) | 35/250 = 0.140 (0.099–0.189) | +10.9 ± 30.0 (failed 0) |
| 80 ms | 40 Hz | 217/250 = 0.868 (0.820–0.907) | 15/250 = 0.060 (0.034–0.097) | -0.8 ± 21.4 (failed 0) |
| 80 ms | 80 Hz | 166/250 = 0.664 (0.602–0.722) | 5/250 = 0.020 (0.007–0.046) | -15.1 ± 17.6 (failed 0) |

Max per-study missed-trial fraction (addendum cells): 0.0021.

## Reading the subject counts
| Estimate of subjects for 0.8 power | 10 ms rise | 20 ms rise |
|---|---|---|
| **Primary** (grid 16-48, as simulated) | **14.7 (12.4-16.8)** | **22.8 (20.2-25.4)** |
| Post hoc: grid 8-48, as simulated | 15.4 (13.8-17.1) | 23.2 (20.8-25.6) |
| Post hoc: grid 16-48, offsets centred | 16.0 (14.0-17.9) | 24.6 (22.1-27.1) |
| Post hoc: grid 8-48, offsets centred | 16.2 (14.8-17.8) | 24.8 (22.1-27.4) |
| FC r1 (120 studies, grid 16-48) | 17.0 (14.3-19.4) | 25.3 (21.5-29.1) |
| Analytic nct (measured bias, design SD + detector error) | 17 | 25 |

- **The 10 ms primary value is an extrapolation.** Power is already 0.836 at 16 subjects, the smallest size in the planned grid.
- The primary values are probably a little low. The first-16-subject offsets average +0.20 ms over the 250 studies (expected 0, SE 0.14, z = 1.4). Removing that chance shift gives the "centred" rows. The post-hoc rows were added after seeing the primary result, so they are sensitivity checks and not replacements.
- For planning, about **17 subjects at 10 ms and 25 at 20 ms** is the safer reading. The planned 16 reach 0.8 at neither rise on the analytic value (0.79 and 0.63), and only at 10 ms here.
- Agreement with FC: no uncorrected power or false-pass cell differs (Fisher p >= 0.17). Corrected power is higher here at 30 ms (0.944 vs 0.842, p = 0.003) and 80 ms (0.868 vs 0.783, p = 0.048), out of 20 comparisons. Opposite chance offsets (+0.20 ms here, -0.28 ms in FC's 30 ms pool) go in that direction; this was not tested further.

## Reading the correction
- It only looks like a power gain. False pass is 0.276 at 10 ms and 0.172 at 20 ms against a nominal 0.05, and 0.092 and 0.088 at 30 and 40 ms. **At 80 ms (0.060, p = 0.27) an excess is not shown.** The earlier "at every rise" was an overstatement.
- No cutoff (20/40/80 Hz) keeps false pass at or below 0.05 at 10-20 ms rise.
- One-sided p values are unadjusted. Holm over the five primary rises leaves 10, 20, 30 and 40 ms below 0.05.
- It concerns the offline back-dated lead only, not live `detect_lead_ms`.

## Other finding: rare early IMU onsets
11 of 359,947 detected trials at 10 ms rise have a lead error of -0.5 to -1.4 s. In all 11 the angle detector fired early on angle noise, and the EMG onset was within 3.4 ms. They are kept in every number above. Each one lowers the lead and inflates the study's spread: at 16 subjects, 10 ms rise, 2 of the 4 affected studies fail where they would pass with a typical error. FC's pool has them too (max error 1.39 s).

## Seed lists elsewhere in the repo (504c035)
- `src/emg_onset_lead/sources.py:56-62`: **`SyntheticSource` had the same collision** (session-0 offset = 0.75 x subject offset; subject-mean true-lead SD 10.55 ms over 640 subjects). It affected the subject spread of `run --synthetic`, for which README:191-193 (at 504c035) quoted 8.91 ms. **Fixed in #7**, after this run; see README "Run". The results above don't use `SyntheticSource`.
- `tests/test_emg_grid.py:38,67`: 10 of 40 keys end in 0, with no collision between cells. `[1, 0, 0]` equals `default_rng(1)`. No subject structure.
- `tests/test_emg_short_burst.py:47`: keys end in the duration (40/60/80), never 0.
- The PR #4 bounds do not depend on between-subject spread:

| Bound | Merged seeds | Offsets fixed only | Fully re-seeded |
|---|---|---|---|
| Lead error, 6 ms (pooled mean, ms) | -3.879 | -3.907 | -3.996 |
| Triggered delay 47.0 +/- 1.5 ms (mean) | 46.56 | 46.56 (identical) | 47.36 |
| Guard band bias / u (ms) | -3.930 / 0.139 | -3.943 / 0.147 | -3.927 / 0.135 |
| 40 ms / SNR 3, at least 98/150 | 121 | no subjects | 110-128 over 10 fresh streams |

`test_pass_fail_flips` gives the same verdicts. With `SyntheticSource` re-seeded in a scratch copy, the whole suite still gives 101 passed, 7 xfailed. No bound was changed.

## Reproduce
```
python scripts/row12_power_sim.py simulate --studies 0:25 --out out/row12/raw_000.npz   # ... to 225:250
python scripts/row12_power_sim.py score out/row12/raw_*.npz --out out/row12/scored.json
python scripts/row12_power_sim.py report out/row12/scored.json --out out/row12/tables.md
```
About 45 min on 8 cores for all 250 studies.
