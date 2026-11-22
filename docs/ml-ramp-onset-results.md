# ML ramp-model EMG onset fit: row-12 kill test (ALL DATA SYNTHETIC)

**Verdict: KILL.** 12 of 20 cells have a false-pass rate at 68.67 ms whose one-sided 95% Clopper-Pearson lower bound is above 0.05. The rule needs only one. The worst cell is exponential, 80 ms, SNR 3: 145/1000 = 0.145 (0.124–0.168). All 10 exponential cells fail, and so do both smoothstep 10 ms cells. Do not use this estimator as a row-12 gate method.

No recorded EMG or motion was used. `src/` is unchanged, and so are all protected values: 68.67 ms gate, proportional mode, 30 ms hold, 6 ms margin, 10% miss cap, 1.0 s / 12-SD baseline, 16x3 `PLANNED_SHAPE`.

## What was built
`scripts/ml_ramp_onset.py`:
1. Whiten the band-passed EMG with an AR(10) filter estimated on the baseline (Yule-Walker), scaled to unit baseline SD.
2. Model: amplitude 1 before t0, then 1 + (G − 1)·clip((t − t0)/R, 0, 1). Gaussian likelihood on 2 ms energy blocks.
3. (t0, R, G) is the likelihood arg-max: an FFT grid over every t0 (14 rises, 12 gains), then a continuous local search.
4. Reported onset: where the fitted ramp crosses 10% of its plateau, t0 + 0.10·R, minus the energy delay of band-pass plus whitener.

This differs from `RampOnset` in `src/`, which reports the posterior mean of the ramp start on the unwhitened band-passed signal.

## True onset (one definition, stated before any gate result)
The **first time the noise-free amplitude crosses P = 0.10 of its plateau excursion**, (a − 1)/(G − 1) = 0.10, for every shape. True lead is the true IMU onset minus that time. This is later than the row-12 sim's definition (the rise start), so these numbers are not directly comparable with that sim's uncorrected baseline table.

## Design
- **Shapes:** smoothstep of full duration R, and exponential 1 − exp(−t/tau) with tau = 0.6084·R/ln 9. Both have the same 10–90% span, 0.6084·R. "Rise" is R.
- **Cells:** 2 shapes × rise 10/20/30/40/80 ms × SNR 3/5. SNR is the plateau/baseline amplitude ratio, times the per-session gain U(0.7, 1.3).
- **Studies:** 1000 per cell, each 16x3x10. Subject/session/trial SD 8/6/10 ms. True leads 74.67 and 68.67 ms.
- **Seeds:** the fixed SeedSequence spawn-key tree of `scripts/row12_power_sim.py` (commit eb339a1, branch `sim/row12-power-seedfix`), root 20261122, studies 0–999. The design draws are identical to that script's (study 0, `prechecks.txt`). The drawn EMG time is used as the P-crossing instead of the ramp start.
- **Reuse:** noise and IMU onsets are shared by all cells, and each EMG fit serves both leads. Cells are therefore paired, so 12 kills are not 12 independent pieces of evidence.
- **Detection:** trigger, rest check and IMU onset are the repo's. The batch trigger equals `EmgOnsetDetector` (480 of 480 in `prechecks.txt`, plus a test). The fit reads the repo's causal window: trigger − 0.4 s to trigger + 0.1 s (last 1000 samples).
- **Gate:** `LeadStats.summarize` at its defaults.
- **Interpretation choices (mine, not in the task text):** P = 0.10, AR order 10, 2 ms blocks, the matched 10–90% span for the exponential, and the grid limits. One configuration was run; nothing was tuned on gate results.

## Results
Studies per cell: 1000, each 16x3x10. ALL SYNTHETIC. True onset = first crossing of 0.10 of the plateau excursion, every shape.

### Primary table (P = 0.10)

Lead error = measured minus true lead (positive = reads high = unsafe). Bias = mean over studies of the study-mean error, with its Monte-Carlo SE. Spread = SD across studies of the study-mean error. Both use the 74.67 ms studies. EMG columns are the onset-estimate error alone (estimated minus true onset; an early estimate is negative and makes the lead read high). Rare IMU noise triggers (hundreds of ms early) are left in and enter the lead columns only. Pass rates carry two-sided 95% Clopper-Pearson intervals; 'CP 1-s LB' is the one-sided 95% lower bound the kill rule uses.

| Shape | Rise (ms) | SNR | Lead bias (ms) | Spread across studies (ms) | EMG onset bias (ms) | EMG spread across studies (ms) | EMG per-trial SD (ms) | Missed trials | Power @ 74.67 | False pass @ 68.67 | CP 1-s LB | Kill? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| smooth | 10 | 3 | +0.65 ± 0.01 | 0.31 | -0.71 | 0.21 | 4.52 | 0.0002 | 892/1000 = 0.892 (0.871–0.911) | 79/1000 = 0.079 (0.063–0.097) | 0.0654 | **KILL** |
| smooth | 10 | 5 | +0.42 ± 0.01 | 0.27 | -0.48 | 0.12 | 2.75 | 0.0002 | 875/1000 = 0.875 (0.853–0.895) | 69/1000 = 0.069 (0.054–0.087) | 0.0563 | **KILL** |
| smooth | 20 | 3 | -0.02 ± 0.01 | 0.33 | -0.04 | 0.24 | 5.32 | 0.0002 | 823/1000 = 0.823 (0.798–0.846) | 44/1000 = 0.044 (0.032–0.059) | 0.0338 | no |
| smooth | 20 | 5 | +0.11 ± 0.01 | 0.28 | -0.17 | 0.14 | 3.19 | 0.0002 | 842/1000 = 0.842 (0.818–0.864) | 51/1000 = 0.051 (0.038–0.067) | 0.0401 | no |
| smooth | 30 | 3 | -0.36 ± 0.01 | 0.36 | +0.30 | 0.27 | 6.04 | 0.0002 | 789/1000 = 0.789 (0.762–0.814) | 34/1000 = 0.034 (0.024–0.047) | 0.0251 | no |
| smooth | 30 | 5 | -0.00 ± 0.01 | 0.29 | -0.06 | 0.16 | 3.58 | 0.0002 | 827/1000 = 0.827 (0.802–0.850) | 46/1000 = 0.046 (0.034–0.061) | 0.0356 | no |
| smooth | 40 | 3 | -0.54 ± 0.01 | 0.39 | +0.48 | 0.30 | 6.68 | 0.0002 | 764/1000 = 0.764 (0.736–0.790) | 28/1000 = 0.028 (0.019–0.040) | 0.0200 | no |
| smooth | 40 | 5 | -0.03 ± 0.01 | 0.30 | -0.04 | 0.17 | 3.90 | 0.0002 | 823/1000 = 0.823 (0.798–0.846) | 45/1000 = 0.045 (0.033–0.060) | 0.0347 | no |
| smooth | 80 | 3 | -0.74 ± 0.01 | 0.46 | +0.68 | 0.39 | 8.69 | 0.0002 | 733/1000 = 0.733 (0.704–0.760) | 26/1000 = 0.026 (0.017–0.038) | 0.0183 | no |
| smooth | 80 | 5 | +0.08 ± 0.01 | 0.33 | -0.15 | 0.22 | 4.93 | 0.0002 | 832/1000 = 0.832 (0.807–0.855) | 43/1000 = 0.043 (0.031–0.057) | 0.0330 | no |
| exp | 10 | 3 | +1.32 ± 0.01 | 0.31 | -1.39 | 0.20 | 4.43 | 0.0002 | 939/1000 = 0.939 (0.922–0.953) | 135/1000 = 0.135 (0.114–0.158) | 0.1175 | **KILL** |
| exp | 10 | 5 | +0.93 ± 0.01 | 0.27 | -0.99 | 0.12 | 2.63 | 0.0002 | 921/1000 = 0.921 (0.903–0.937) | 101/1000 = 0.101 (0.083–0.121) | 0.0857 | **KILL** |
| exp | 20 | 3 | +1.05 ± 0.01 | 0.32 | -1.12 | 0.23 | 4.95 | 0.0002 | 923/1000 = 0.923 (0.905–0.939) | 109/1000 = 0.109 (0.090–0.130) | 0.0932 | **KILL** |
| exp | 20 | 5 | +0.75 ± 0.01 | 0.27 | -0.81 | 0.13 | 2.89 | 0.0002 | 902/1000 = 0.902 (0.882–0.920) | 88/1000 = 0.088 (0.071–0.107) | 0.0737 | **KILL** |
| exp | 30 | 3 | +0.98 ± 0.01 | 0.34 | -1.04 | 0.25 | 5.44 | 0.0002 | 920/1000 = 0.920 (0.901–0.936) | 100/1000 = 0.100 (0.082–0.120) | 0.0848 | **KILL** |
| exp | 30 | 5 | +0.70 ± 0.01 | 0.27 | -0.77 | 0.14 | 3.11 | 0.0002 | 901/1000 = 0.901 (0.881–0.919) | 83/1000 = 0.083 (0.067–0.102) | 0.0691 | **KILL** |
| exp | 40 | 3 | +1.00 ± 0.01 | 0.36 | -1.06 | 0.27 | 5.88 | 0.0002 | 918/1000 = 0.918 (0.899–0.934) | 103/1000 = 0.103 (0.085–0.124) | 0.0876 | **KILL** |
| exp | 40 | 5 | +0.73 ± 0.01 | 0.28 | -0.79 | 0.15 | 3.30 | 0.0002 | 902/1000 = 0.902 (0.882–0.920) | 84/1000 = 0.084 (0.068–0.103) | 0.0700 | **KILL** |
| exp | 80 | 3 | +1.39 ± 0.01 | 0.41 | -1.46 | 0.32 | 7.25 | 0.0002 | 944/1000 = 0.944 (0.928–0.957) | 145/1000 = 0.145 (0.124–0.168) | 0.1270 | **KILL** |
| exp | 80 | 5 | +0.97 ± 0.01 | 0.30 | -1.03 | 0.17 | 3.88 | 0.0002 | 920/1000 = 0.920 (0.901–0.936) | 101/1000 = 0.101 (0.083–0.121) | 0.0857 | **KILL** |

**Verdict by the stated rule (any cell with CP one-sided 95% LB > 0.05): KILL.** 12 of 20 cells exceed.
Cells over the line: smooth 10 ms SNR 3: 79/1000, LB 0.0654; smooth 10 ms SNR 5: 69/1000, LB 0.0563; exp 10 ms SNR 3: 135/1000, LB 0.1175; exp 10 ms SNR 5: 101/1000, LB 0.0857; exp 20 ms SNR 3: 109/1000, LB 0.0932; exp 20 ms SNR 5: 88/1000, LB 0.0737; exp 30 ms SNR 3: 100/1000, LB 0.0848; exp 30 ms SNR 5: 83/1000, LB 0.0691; exp 40 ms SNR 3: 103/1000, LB 0.0876; exp 40 ms SNR 5: 84/1000, LB 0.0700; exp 80 ms SNR 3: 145/1000, LB 0.1270; exp 80 ms SNR 5: 101/1000, LB 0.0857.

### Sensitivity to the onset definition (same trials, true lead and IMU onset kept relative to each definition's own crossing; not the kill test)

| Shape | Rise | SNR | P=0: bias / power / false pass | P=0.05: bias / power / false pass | P=0.1: bias / power / false pass | P=0.25: bias / power / false pass | P=0.5: bias / power / false pass |
|---|---|---|---|---|---|---|---|
| smooth | 10 | 3 | -0.26 / 0.807 / 0.034 | +0.57 / 0.884 / 0.074 | +0.65 / 0.892 / 0.079 | +0.38 / 0.871 / 0.063 | -0.50 / 0.772 / 0.028 |
| smooth | 10 | 5 | -0.62 / 0.761 / 0.027 | +0.27 / 0.858 / 0.059 | +0.42 / 0.875 / 0.069 | +0.35 / 0.870 / 0.062 | -0.20 / 0.801 / 0.037 |
| smooth | 20 | 3 | -2.29 / 0.489 / 0.005 | -0.40 / 0.785 / 0.032 | -0.02 / 0.823 / 0.044 | +0.12 / 0.840 / 0.051 | -0.54 / 0.769 / 0.028 |
| smooth | 20 | 5 | -2.19 / 0.500 / 0.005 | -0.29 / 0.797 / 0.032 | +0.11 / 0.842 / 0.051 | +0.30 / 0.866 / 0.061 | -0.26 / 0.799 / 0.037 |
| smooth | 30 | 3 | -3.91 / 0.225 / 0.000 | -1.01 / 0.706 / 0.017 | -0.36 / 0.789 / 0.034 | +0.07 / 0.835 / 0.050 | -0.53 / 0.772 / 0.027 |
| smooth | 30 | 5 | -3.52 / 0.294 / 0.001 | -0.64 / 0.759 / 0.024 | -0.00 / 0.827 / 0.046 | +0.38 / 0.873 / 0.064 | -0.29 / 0.793 / 0.035 |
| smooth | 40 | 3 | -5.34 / 0.076 / 0.000 | -1.45 / 0.629 / 0.011 | -0.54 / 0.764 / 0.028 | +0.12 / 0.843 / 0.048 | -0.52 / 0.773 / 0.029 |
| smooth | 40 | 5 | -4.74 / 0.133 / 0.000 | -0.89 / 0.729 / 0.020 | -0.03 / 0.823 / 0.045 | +0.52 / 0.886 / 0.073 | -0.32 / 0.793 / 0.032 |
| smooth | 80 | 3 | -10.37 / 0.000 / 0.000 | -2.56 / 0.434 / 0.003 | -0.74 / 0.733 / 0.026 | +0.64 / 0.884 / 0.081 | -0.57 / 0.759 / 0.028 |
| smooth | 80 | 5 | -9.37 / 0.000 / 0.000 | -1.65 / 0.599 / 0.010 | +0.08 / 0.832 / 0.043 | +1.21 / 0.934 / 0.125 | -0.43 / 0.787 / 0.030 |
| exp | 10 | 3 | +2.03 / 0.965 / 0.218 | +1.67 / 0.956 / 0.174 | +1.32 / 0.939 / 0.135 | +0.33 / 0.865 / 0.060 | -1.04 / 0.701 / 0.018 |
| exp | 10 | 5 | +1.44 / 0.947 / 0.146 | +1.18 / 0.932 / 0.124 | +0.93 / 0.921 / 0.101 | +0.23 / 0.853 / 0.058 | -0.65 / 0.760 / 0.024 |
| exp | 20 | 3 | +1.92 / 0.961 / 0.207 | +1.48 / 0.944 / 0.146 | +1.05 / 0.923 / 0.109 | -0.11 / 0.817 / 0.042 | -1.48 / 0.630 / 0.011 |
| exp | 20 | 5 | +1.42 / 0.946 / 0.143 | +1.07 / 0.928 / 0.111 | +0.75 / 0.902 / 0.088 | -0.12 / 0.809 / 0.041 | -1.01 / 0.706 / 0.017 |
| exp | 30 | 3 | +2.02 / 0.961 / 0.222 | +1.49 / 0.945 / 0.154 | +0.98 / 0.920 / 0.100 | -0.38 / 0.786 / 0.030 | -1.80 / 0.581 / 0.010 |
| exp | 30 | 5 | +1.54 / 0.954 / 0.155 | +1.11 / 0.930 / 0.119 | +0.70 / 0.901 / 0.083 | -0.35 / 0.790 / 0.029 | -1.27 / 0.670 / 0.012 |
| exp | 40 | 3 | +2.23 / 0.966 / 0.252 | +1.60 / 0.949 / 0.171 | +1.00 / 0.918 / 0.103 | -0.58 / 0.770 / 0.027 | -2.08 / 0.532 / 0.006 |
| exp | 40 | 5 | +1.74 / 0.958 / 0.178 | +1.22 / 0.939 / 0.126 | +0.73 / 0.902 / 0.084 | -0.53 / 0.774 / 0.026 | -1.49 / 0.633 / 0.010 |
| exp | 80 | 3 | +3.38 / 0.988 / 0.428 | +2.35 / 0.967 / 0.267 | +1.39 / 0.944 / 0.145 | -1.05 / 0.698 / 0.017 | -2.86 / 0.393 / 0.002 |
| exp | 80 | 5 | +2.62 / 0.976 / 0.314 | +1.76 / 0.960 / 0.186 | +0.97 / 0.920 / 0.101 | -0.98 / 0.712 / 0.018 | -1.96 / 0.545 / 0.007 |

## Check tables (`ml_ramp_row12_checks.py`)

Studies: 1000. ALL SYNTHETIC.

### Oracle floor (gate pass rates with perfect EMG onsets)

| Case | Power @ 74.67 | False pass @ 68.67 |
|---|---|---|
| no EMG or IMU error | 833/1000 = 0.833 (0.808–0.856) | 44/1000 = 0.044 (0.032–0.059) |
| no EMG error; simulated IMU error | 822/1000 = 0.822 (0.797–0.845) | 42/1000 = 0.042 (0.030–0.056) |
| no EMG error; IMU error, |IMU error| > 100 ms trials dropped | 828/1000 = 0.828 (0.803–0.851) | 42/1000 = 0.042 (0.030–0.056) |

IMU errors beyond 100 ms: 32 of 960000 IMU onsets (2 leads share noise); trial miss fraction 0.00000; IMU error SD excluding those 1.99 ms.

### False pass vs 0.05: one-sided binomial p, unadjusted and Holm over 20 cells

| Cell | False pass | CP 1-s LB | p | Holm p | Kill rule (LB > 0.05) |
|---|---|---|---|---|---|
| smooth 10 ms SNR 3 | 79/1000 = 0.079 (0.063–0.097) | 0.0654 | 5.9e-05 | 0.00059 | KILL |
| smooth 10 ms SNR 5 | 69/1000 = 0.069 (0.054–0.087) | 0.0563 | 0.0051 | 0.046 | KILL |
| smooth 20 ms SNR 3 | 44/1000 = 0.044 (0.032–0.059) | 0.0338 | 0.83 | 1 | no |
| smooth 20 ms SNR 5 | 51/1000 = 0.051 (0.038–0.067) | 0.0401 | 0.46 | 1 | no |
| smooth 30 ms SNR 3 | 34/1000 = 0.034 (0.024–0.047) | 0.0251 | 0.99 | 1 | no |
| smooth 30 ms SNR 5 | 46/1000 = 0.046 (0.034–0.061) | 0.0356 | 0.74 | 1 | no |
| smooth 40 ms SNR 3 | 28/1000 = 0.028 (0.019–0.040) | 0.0200 | 1 | 1 | no |
| smooth 40 ms SNR 5 | 45/1000 = 0.045 (0.033–0.060) | 0.0347 | 0.79 | 1 | no |
| smooth 80 ms SNR 3 | 26/1000 = 0.026 (0.017–0.038) | 0.0183 | 1 | 1 | no |
| smooth 80 ms SNR 5 | 43/1000 = 0.043 (0.031–0.057) | 0.0330 | 0.86 | 1 | no |
| exp 10 ms SNR 3 | 135/1000 = 0.135 (0.114–0.158) | 0.1175 | 5.3e-25 | 1e-23 | KILL |
| exp 10 ms SNR 5 | 101/1000 = 0.101 (0.083–0.121) | 0.0857 | 3.9e-11 | 6.3e-10 | KILL |
| exp 20 ms SNR 3 | 109/1000 = 0.109 (0.090–0.130) | 0.0932 | 5.8e-14 | 1e-12 | KILL |
| exp 20 ms SNR 5 | 88/1000 = 0.088 (0.071–0.107) | 0.0737 | 3.5e-07 | 4.5e-06 | KILL |
| exp 30 ms SNR 3 | 100/1000 = 0.100 (0.082–0.120) | 0.0848 | 8.4e-11 | 1.2e-09 | KILL |
| exp 30 ms SNR 5 | 83/1000 = 0.083 (0.067–0.102) | 0.0691 | 6.8e-06 | 7.5e-05 | KILL |
| exp 40 ms SNR 3 | 103/1000 = 0.103 (0.085–0.124) | 0.0876 | 8.2e-12 | 1.4e-10 | KILL |
| exp 40 ms SNR 5 | 84/1000 = 0.084 (0.068–0.103) | 0.0700 | 3.8e-06 | 4.6e-05 | KILL |
| exp 80 ms SNR 3 | 145/1000 = 0.145 (0.124–0.168) | 0.1270 | 6.1e-30 | 1.2e-28 | KILL |
| exp 80 ms SNR 5 | 101/1000 = 0.101 (0.083–0.121) | 0.0857 | 3.9e-11 | 6.3e-10 | KILL |

Cells with Holm p < 0.05: 12.

### Fitted ramp vs truth, all trials (EMG only)

| Cell | start error mean ± SD (ms) | fitted rise mean (ms) | true rise to plateau | fitted gain mean | P=0.10 crossing error: mean / SD / 1st-99th pct (ms) |
|---|---|---|---|---|---|
| smooth 10 SNR 3 | +0.20 ± 5.50 | 10.5 | 10 | 3.03 | -0.71 / 4.53 / -16.3..+6.2 |
| smooth 10 SNR 5 | +0.56 ± 3.16 | 9.2 | 10 | 5.04 | -0.48 / 2.75 / -6.7..+3.7 |
| smooth 20 SNR 3 | +2.22 ± 6.50 | 16.5 | 20 | 3.03 | -0.04 / 5.32 / -17.0..+9.7 |
| smooth 20 SNR 5 | +2.13 ± 3.77 | 16.1 | 20 | 5.04 | -0.17 / 3.20 / -7.4..+6.0 |
| smooth 30 SNR 3 | +3.85 ± 7.42 | 23.2 | 30 | 3.03 | +0.30 / 6.05 / -18.0..+12.7 |
| smooth 30 SNR 5 | +3.46 ± 4.28 | 23.5 | 30 | 5.04 | -0.06 / 3.58 / -8.1..+7.7 |
| smooth 40 SNR 3 | +5.28 ± 8.23 | 30.4 | 40 | 3.03 | +0.48 / 6.69 / -18.8..+15.3 |
| smooth 40 SNR 5 | +4.68 ± 4.71 | 31.1 | 40 | 5.04 | -0.04 / 3.90 / -8.7..+8.9 |
| smooth 80 SNR 3 | +10.30 ± 10.78 | 60.4 | 80 | 3.03 | +0.68 / 8.70 / -22.5..+23.1 |
| smooth 80 SNR 5 | +9.31 ± 6.03 | 62.1 | 80 | 5.04 | -0.15 / 4.94 / -11.2..+11.7 |
| exp 10 SNR 3 | -2.09 ± 5.42 | 10.0 | asymptotic, tau 2.8 | 3.03 | -1.39 / 4.43 / -17.1..+5.4 |
| exp 10 SNR 5 | -1.50 ± 3.03 | 8.0 | asymptotic, tau 2.8 | 5.03 | -0.99 / 2.63 / -7.3..+2.6 |
| exp 20 SNR 3 | -1.98 ± 6.14 | 14.5 | asymptotic, tau 5.5 | 3.02 | -1.12 / 4.95 / -18.1..+7.4 |
| exp 20 SNR 5 | -1.48 ± 3.44 | 12.5 | asymptotic, tau 5.5 | 5.02 | -0.81 / 2.89 / -8.0..+3.9 |
| exp 30 SNR 3 | -2.08 ± 6.82 | 19.1 | asymptotic, tau 8.3 | 3.01 | -1.04 / 5.45 / -19.3..+9.1 |
| exp 30 SNR 5 | -1.61 ± 3.79 | 17.1 | asymptotic, tau 8.3 | 5.00 | -0.77 / 3.11 / -8.7..+4.8 |
| exp 40 SNR 3 | -2.29 ± 7.45 | 24.0 | asymptotic, tau 11.1 | 3.01 | -1.06 / 5.89 / -20.5..+10.5 |
| exp 40 SNR 5 | -1.80 ± 4.09 | 21.8 | asymptotic, tau 11.1 | 4.99 | -0.79 / 3.30 / -9.3..+5.5 |
| exp 80 SNR 3 | -3.44 ± 9.44 | 43.2 | asymptotic, tau 22.2 | 2.97 | -1.46 / 7.26 / -24.3..+14.6 |
| exp 80 SNR 5 | -2.69 ± 5.02 | 39.9 | asymptotic, tau 22.2 | 4.91 | -1.03 / 3.88 / -11.4..+7.1 |

## Why it fails
- The model is a straight ramp. On an exponential burst the fit starts early (mean start error −1.5 to −3.4 ms) and its 10% crossing is 0.8–1.5 ms early, so the lead reads 0.7–1.4 ms high. At 10 ms smoothstep the crossing is 0.5–0.7 ms early.
- A bias that small is enough. The spread across studies is only 0.27–0.46 ms, and the analytic one-sided t false pass is 0.072 / 0.086 / 0.147 at a bias of +0.42 / +0.65 / +1.39 ms. The simulated values are 0.069 / 0.079 / 0.145 (`out/audit.txt`).
- The gate itself is not the cause. With perfect EMG onsets on the same studies, false pass is 0.044 (no error) or 0.042 (simulated IMU error).
- Smoothstep at 20–80 ms survives: bias −0.74 to +0.11 ms, false pass 0.026–0.051. "Unbiased" holds only for that shape and range.
- On a linear ramp, where the model is right, the 10% crossing is +0.1 to +1.0 ms late (`prechecks.txt`). The early reading comes from the burst shape.
- Power in the killed cells (0.875–0.944) is inflated by the same bias and must not be quoted as a gain.

## Onset-definition sensitivity (post hoc, not the kill test)
Cells over the kill line by P (`out/audit.txt`): 10/20 at P = 0, 11/20 at 0.05, **12/20 at 0.10**, 5/20 at 0.25, 0/20 at 0.50. P = 0.50 is not a rescue. Choosing it after seeing this table would be tuning on the kill test, it reads late in every cell (−0.20 to −2.86 ms), and only 1 of 20 cells reaches 0.8 power.

## Checks
- An independent numpy gate agrees with `LeadStats` on 40,000 of 40,000 decisions. Studies 3 and 617 re-run bit-identically (`audit.txt`).
- Realised SD subject/session/trial 7.961/5.994/9.985 ms, subject mean 8.832 ms (design 8.907). Subject–session correlations are within ±0.013, so the old seed collision is absent.
- 96 of 480,000 trials are missed, the same ones in every cell. All 96 are rejected by the unchanged baseline rest check (`miss_cause.txt`). The largest per-study miss fraction is 0.0042 (cap 0.10).
- 32 of 960,000 IMU onsets are more than 100 ms off, all early (−1233 to −540 ms, `prechecks.txt`), consistent with noise triggers in the unchanged IMU detector. They are kept. They push the lead down, the safe direction.
- All 12 killed cells also have Holm-adjusted p < 0.05 (largest: smoothstep 10 ms SNR 5, p = 0.046).

## Process disclosures
- **Pilot:** 12 studies on seeds 100000–100011, outside the main run. I saw its bias and gate tables before the full run. The estimator did not change afterwards; only EMG-only columns were added to the report. The estimator's one fix (half-sample block-centre time) came earlier, from linear-ramp smoke tests.
- **Aborted start:** the first full-run launch was stopped for lint fixes before any shard finished. The run used a frozen copy of commit b175d36. `ml_ramp_row12_checks.py` (1c8c247) was committed before any full-run result.
- **Runtime:** pilot forecast about 95 min (12 studies in 68 s on 6 processes; console output only, not saved to a file). Actual: 6571 s (110 min), summed over shards in `run_full.log`, plus 4 min scoring (`out/score.log`).

## Sources
- Staude GH. "Precise onset detection of human motor responses using a whitening filter and the log-likelihood-ratio test." *IEEE Trans Biomed Eng* 2001;48(11):1292–1305. https://doi.org/10.1109/10.959325 (PMID 11686628). **Abstract read; full text not read.** The abstract is about kinematic signals (e.g. joint angle), not EMG. It supports the whitening + likelihood idea only.
- Staude G, Flachenecker C, Daumer M, Wolf W. "Onset Detection in Surface Electromyographic Signals: A Systematic Comparison of Methods." *EURASIP J Adv Signal Process* 2001. https://doi.org/10.1155/S1110865701000191. **Abstract read on the Springer page; full text not obtained (sign-in page).** It reports that detection depends on onset rise time, SNR and background activity, and uses a maximum-likelihood estimator as its reference.
- **Not found:** I could not verify in either paper a page or equation for a ramp amplitude profile. The ramp model is this repo's (`RampOnset`) with the whitening and arg-max changes above. Do not cite it as Staude's published equation, and do not read this kill as a verdict on Staude's own methods.

## Limits
Synthetic only: white noise, sustained bursts, minimum-jerk angle. No real EMG, no linear-ramp or other shapes at the gate, SNR 3 and 5 only, one AR order. The test covers the offline back-dated lead and says nothing about live detection delay. The incumbent `RampOnset` was not re-run under the P = 0.10 definition, so there is no paired comparison.

## Reproduce
```
python scripts/ml_ramp_row12_sim.py simulate --studies 0:25 --out out/ml/raw_0000.npz --procs 6   # ... to 975:1000
python scripts/ml_ramp_row12_sim.py score out/ml/raw_*.npz --out out/ml/scored.json
python scripts/ml_ramp_row12_sim.py report out/ml/scored.json --out out/ml/tables.md
python scripts/ml_ramp_row12_checks.py --raw 'out/ml/raw_*.npz' --scored out/ml/scored.json --out out/ml/checks.md
```
Raw data, scored results, audit scripts and logs: `/workspace/claude_sessions/3_ml-ramp-onset/` (`out/`, `audit.py`, `miss_cause.py`, `prechecks.py`, `run_full.sh`, `run_full.log`). They are not committed.
