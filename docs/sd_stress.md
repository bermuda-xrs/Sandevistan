# Row-12 power under larger SD (SD stress runs, SYNTHETIC)

**Every number here is SYNTHETIC.** No recorded EMG or motion was used, and nothing here validates the detectors on real signals. The design values stay 8/6/10 ms (subject/session/trial SD). The gate (68.67 ms, proportional), the 16 x 3 PLANNED_SHAPE and `src/` are unchanged.

Scripts: `scripts/sd_stress.py` (run, score, report) and `tests/test_sd_stress.py`. Output: `results/sd_stress/scored.json` (all counts), `results/sd_stress/tables.md` (every table, with Clopper-Pearson intervals), `results/sd_stress/raw_manifest.sha256` (sha256 of the 40 raw chunks, 79 MB, not committed).

## Question
Row 12 plans 16 subjects x 3 sessions x 10 trials against a gate of 68.67 ms with a true lead of 74.67 ms (gate + 6 ms margin). `docs/row12_power_sim.md` assumed 8/6/10 ms. The session-4 literature report (`claude_sessions/4_published-variance/REPORT.md` section 1, not in this repo) recommends 12/4/10 as a nominal split and 8-16 ms as the plausible range of the subject SD. This run asks what happens to power, to the subject count for 0.8 power and to the false-pass rate if the real spread is larger.

Caveats on those inputs, from the Fact Checker's review of the session-4 report (`/workspace/fc_check/sessions_r1/S4_verdict.md`, items 12, 19 and 20; quoted, not re-derived here):
- 12 ms is in the upper half of the elbow-flexion evidence, not its midpoint (midrange 10.85 ms, median 11.45 ms). One-session between-subject SDs already contain session and trial variance, so using 12 ms as a pure subject SD and adding a session SD double-counts: the pure subject SD is about 11.3 ms, which is the conservative direction.
- Lowering the session SD from 6 to 4 ms is the one anti-conservative change in 12/4/10. It rests on abstract-level evidence (2.3-4.9 ms), and keeping 6 ms costs nothing. In this run 12/4/10 and 12/6/10 differ little (power at 16 subjects, 10 ms rise, with bias: 0.569 against 0.555).
- The one source with paired session and subject SDs (Almosnino) has an SEM/SD ratio of 0.74, which fits 8/6/10 (0.72), not 12/4/10 (0.505).

## Study count (fixed before the 1000-study run)
**1000 independent studies per cell, ids 1000-1999**, written into the script and pushed as commit c818622 before the first chunk was computed. Before that commit an 8-study pipeline smoke run (ids 900000-900007, uncommitted script at afafe8b, outputs in the session folder `scratch/smoke.json` and `scratch/smoke.md`) was scored and its tables were looked at. It found one report crash (division by zero when a subject size has no study) and that is all it was used for; the study count and configs were already in the script, and its ids are disjoint from 1000-1999. The raw data were produced at c818622 (`head` stored in every chunk). Afterwards main was merged into the branch (the detector changes on main are parameter validation only, so the settings used here are unaffected) and, after the run, three things were added that read the stored results without changing them: a compact summary section and a Claim checks section in `report`, the `sensitivity` command, and stronger tests. The simulation and scoring code are unchanged since c818622.

## Method
- **Script:** `scripts/sd_stress.py`, which imports `scripts/row12_power_sim.py` from main unchanged. It uses the same smoothstep EMG bursts (SNR 15 × U(0.7, 1.3), 2 kHz), the same minimum-jerk 60° angle move (150 Hz, 0.3° noise), the same detectors at default settings, and the same SeedSequence spawn-key tree (study → subject → session → trial).
  - It does not use `SyntheticSource`: that class hard-codes a 40 ms rise, and the row-12 study is built on `row12_power_sim.py`. Both use spawn-tree seeding (PR #6 and PR #7, both on main).
  - `tests/test_sd_stress.py` also checks that the merged `SyntheticSource` gives the 12/4/10 spread.
- **Configs (subject/session/trial SD, ms):** 12/4/10, then 8/6/10 (design), 10/6/10, 12/6/10, 14/6/10 and 16/6/10. The sweep keeps session and trial at the design 6/10. The design values in the repo stay 8/6/10; nothing in `src/` changed.
- **Draws:** each config rescales the same standard-normal draws of the row-12 tree (subject × SD/8, session × SD/6, trial × SD/10). So the configs differ only in spread, never in chance. Gains, EMG onsets, EMG noise and angle noise are also shared.
- **Studies:** 1000 independent 16 × 3 × 10 studies per cell, ids 1000–1999. The row-12 doc used ids 0–249, so no draw is reused. That is 16,000 subjects and 480,000 trials per cell.
- **Subject counts other than 16:** the same 16,000 subjects, in id order, are cut into floor(16000/n) disjoint n-subject studies, for n = 8, 12, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 128, 160, 192, 256, 320, 400 and 512.
- **Gate:** the same as `LeadStats.summarize`: a one-sided 95% t lower bound on subject means (each the mean of session means) must exceed 68.67 ms, at most 10% of trials may be missed, and `planned = (n, 3)`. The vectorised gate is checked against `LeadStats.summarize` on every 16-subject study.
- **True leads:** 74.67 ms (gate + 6 ms margin) for power; exactly 68.67 ms for false pass.
- **With bias:** the detected lead, as the gate sees it.
- **Bias removed:** the same leads minus the pooled mean detector error (detected − true) of that config, rise and true lead. Detector noise stays. This needs the true lead, so it is an oracle reference, not a correction you could apply to recorded data. It is not the PR #6 rise-time correction.
- **No detector:** the gate run on the true leads.
- **Subjects for 0.8 power:**
  - Fit: a probit-in-√n binomial ML fit over all grid sizes (the row-12 method). The 95% CI comes from 1000 resamples of subjects with replacement, each cut into studies again.
  - "Not reachable" means the pooled mean lead is at or below the gate.
  - Values above 512 are marked as extrapolated.
  - Cross-check: noncentral-t (`SampleSize.required_subjects`) with the measured mean and SD of subject means.
  - Also reported: the measured power at ceil(n).
- **Intervals:** every proportion carries a 95% Clopper-Pearson interval.
- **Cells are paired, not independent.** All configs, rises and both true leads share draws. Only the 1000 studies within one cell are independent.

## Results (SYNTHETIC; 1000 studies per cell, 16 x 3 x 10, gate 68.67 ms)
The summary tables show point values. Every proportion with its 95% Clopper-Pearson interval, and the subject sweeps, are in `results/sd_stress/tables.md`.

### Realised spread of the true leads (SYNTHETIC)
| SD subject/session/trial | Subject | Session | Trial | Subject-mean SD (95% CI) | Design | Mean offset (ms) |
|---|---|---|---|---|---|---|
| 12/4/10 | 12.05 | 3.98 | 10.03 | 12.40 (12.27–12.54) | 12.36 | +0.077 |
| 8/6/10 | 8.03 | 5.97 | 10.03 | 8.93 (8.84–9.03) | 8.91 | +0.035 |
| 10/6/10 | 10.04 | 5.97 | 10.03 | 10.77 (10.66–10.89) | 10.74 | +0.051 |
| 12/6/10 | 12.05 | 5.97 | 10.03 | 12.66 (12.53–12.80) | 12.62 | +0.068 |
| 14/6/10 | 14.05 | 5.97 | 10.03 | 14.59 (14.43–14.75) | 14.54 | +0.084 |
| 16/6/10 | 16.06 | 5.97 | 10.03 | 16.53 (16.35–16.71) | 16.47 | +0.100 |

### Detector error, detected - true lead (design config; other configs differ from it by at most 0.006 ms in the mean, 0.013 ms in the median and 0.004 ms in the SD; see Claim checks)
| Rise | True lead | Mean (bias) | Median | SD | Missed | < -100 ms | > +100 ms |
|---|---|---|---|---|---|---|---|
| 10 ms | 74.67 | -0.26 | -0.22 | 5.75 | 87/480000 | 13 | 1 |
| 10 ms | 68.67 | -0.26 | -0.22 | 5.72 | 87/480000 | 13 | 1 |
| 20 ms | 74.67 | -1.39 | -1.31 | 5.60 | 87/480000 | 13 | 1 |
| 20 ms | 68.67 | -1.38 | -1.31 | 5.58 | 87/480000 | 13 | 1 |
| 30 ms | 74.67 | -3.01 | -3.02 | 5.74 | 87/480000 | 13 | 1 |
| 30 ms | 68.67 | -3.01 | -3.01 | 5.71 | 87/480000 | 13 | 1 |
| 40 ms | 74.67 | -3.95 | -3.85 | 5.90 | 87/480000 | 13 | 1 |
| 40 ms | 68.67 | -3.95 | -3.84 | 5.88 | 87/480000 | 13 | 1 |
| 80 ms | 74.67 | -7.95 | -7.84 | 6.30 | 87/480000 | 13 | 1 |
| 80 ms | 68.67 | -7.95 | -7.84 | 6.28 | 87/480000 | 13 | 1 |

Range over configs and true leads of the mean error: 10 ms -0.26 to -0.25; 20 ms -1.39 to -1.38; 30 ms -3.01 to -3.00; 40 ms -3.95 to -3.94; 80 ms -7.95 to -7.94 (all cells -7.95 to -0.25).

### Summary: Power at 16 subjects, with bias (true lead 74.67 ms) (SYNTHETIC)
| SD subject/session/trial | 10 ms rise | 20 ms rise | 30 ms rise | 40 ms rise | 80 ms rise | no detector |
|---|---|---|---|---|---|---|
| 12/4/10 | 0.569 | 0.418 | 0.223 | 0.127 | 0.010 | — |
| 8/6/10 (design) | 0.807 | 0.650 | 0.361 | 0.202 | 0.008 | — |
| 10/6/10 | 0.670 | 0.508 | 0.272 | 0.154 | 0.010 | — |
| 12/6/10 | 0.555 | 0.408 | 0.206 | 0.119 | 0.011 | — |
| 14/6/10 | 0.437 | 0.342 | 0.171 | 0.111 | 0.013 | — |
| 16/6/10 | 0.385 | 0.278 | 0.145 | 0.102 | 0.017 | — |

### Summary: Power at 16 subjects, bias removed (oracle) (SYNTHETIC)
| SD subject/session/trial | 10 ms rise | 20 ms rise | 30 ms rise | 40 ms rise | 80 ms rise | no detector |
|---|---|---|---|---|---|---|
| 12/4/10 | 0.606 | 0.598 | 0.599 | 0.597 | 0.601 | 0.610 |
| 8/6/10 (design) | 0.832 | 0.832 | 0.832 | 0.830 | 0.830 | 0.839 |
| 10/6/10 | 0.703 | 0.703 | 0.703 | 0.702 | 0.703 | 0.702 |
| 12/6/10 | 0.582 | 0.582 | 0.585 | 0.585 | 0.584 | 0.588 |
| 14/6/10 | 0.473 | 0.473 | 0.474 | 0.473 | 0.471 | 0.473 |
| 16/6/10 | 0.406 | 0.403 | 0.404 | 0.403 | 0.402 | 0.412 |

### Summary: Subjects for 0.8 power, with bias (95% CI) (SYNTHETIC)
| SD subject/session/trial | 10 ms rise | 20 ms rise | 30 ms rise | 40 ms rise | 80 ms rise | no detector |
|---|---|---|---|---|---|---|
| 12/4/10 | 28.8 (27.4–32.1) | 43.9 (40.5–49.1) | 97.6 (88.6–117.6) | 203.1 (171.9–257.1) | not reachable | — |
| 8/6/10 (design) | 16.4 (15.6–17.5) | 23.9 (22.9–26.4) | 54.9 (50.2–62.1) | 109.5 (99.1–133.3) | not reachable | — |
| 10/6/10 | 22.6 (21.5–24.8) | 34.1 (31.7–37.6) | 77.3 (69.4–89.5) | 152.9 (136.4–195.5) | not reachable | — |
| 12/6/10 | 30.6 (28.5–33.5) | 45.6 (42.3–51.4) | 101.2 (92.2–123.6) | 207.0 (179.6–270.6) | not reachable | — |
| 14/6/10 | 39.4 (36.5–43.8) | 59.7 (54.3–67.8) | 131.9 (117.9–165.1) | 280.7 (227.5–363.7) | not reachable | — |
| 16/6/10 | 50.3 (45.8–56.3) | 75.6 (67.7–87.0) | 165.4 (146.7–211.7) | 335.1 (279.0–475.3) | not reachable | — |

### Summary: Subjects for 0.8 power, bias removed (95% CI) (SYNTHETIC)
| SD subject/session/trial | 10 ms rise | 20 ms rise | 30 ms rise | 40 ms rise | 80 ms rise | no detector |
|---|---|---|---|---|---|---|
| 12/4/10 | 26.6 (25.3–29.6) | 26.6 (25.2–29.5) | 26.5 (25.3–29.5) | 26.5 (25.3–29.5) | 26.6 (25.3–29.5) | 26.3 (25.1–29.3) |
| 8/6/10 (design) | 15.1 (14.4–16.2) | 15.1 (14.4–16.2) | 15.1 (14.4–16.2) | 15.1 (14.4–16.2) | 15.2 (14.4–16.2) | 14.9 (14.2–15.9) |
| 10/6/10 | 20.8 (19.9–22.8) | 20.9 (19.9–22.7) | 20.8 (19.9–22.7) | 20.7 (19.9–22.7) | 20.9 (19.9–22.7) | 20.7 (19.7–22.5) |
| 12/6/10 | 27.7 (26.3–30.8) | 27.8 (26.3–30.9) | 27.7 (26.3–30.8) | 27.6 (26.3–30.8) | 27.7 (26.4–30.8) | 27.6 (26.2–30.6) |
| 14/6/10 | 36.5 (33.8–40.3) | 36.3 (33.7–40.4) | 36.2 (33.8–40.4) | 36.2 (33.8–40.4) | 36.3 (33.8–40.3) | 36.3 (33.6–40.1) |
| 16/6/10 | 46.5 (42.4–51.3) | 46.1 (42.4–51.4) | 46.1 (42.4–51.5) | 46.2 (42.4–51.5) | 46.0 (42.5–51.5) | 46.0 (42.1–51.2) |

### Summary: False pass at 68.67 ms, with bias (nominal 0.05) (SYNTHETIC)
| SD subject/session/trial | 10 ms rise | 20 ms rise | 30 ms rise | 40 ms rise | 80 ms rise | no detector |
|---|---|---|---|---|---|---|
| 12/4/10 | 0.035 | 0.016 | 0.005 | 0.003 | 0.000 | — |
| 8/6/10 (design) | 0.030 | 0.010 | 0.002 | 0.000 | 0.000 | — |
| 10/6/10 | 0.030 | 0.013 | 0.005 | 0.001 | 0.000 | — |
| 12/6/10 | 0.036 | 0.015 | 0.005 | 0.003 | 0.000 | — |
| 14/6/10 | 0.033 | 0.017 | 0.008 | 0.005 | 0.000 | — |
| 16/6/10 | 0.036 | 0.019 | 0.010 | 0.005 | 0.000 | — |

### Summary: False pass at 68.67 ms, bias removed (oracle) (SYNTHETIC)
| SD subject/session/trial | 10 ms rise | 20 ms rise | 30 ms rise | 40 ms rise | 80 ms rise | no detector |
|---|---|---|---|---|---|---|
| 12/4/10 | 0.046 | 0.044 | 0.043 | 0.045 | 0.043 | 0.044 |
| 8/6/10 (design) | 0.041 | 0.041 | 0.043 | 0.041 | 0.044 | 0.039 |
| 10/6/10 | 0.042 | 0.042 | 0.038 | 0.038 | 0.044 | 0.043 |
| 12/6/10 | 0.044 | 0.044 | 0.043 | 0.044 | 0.044 | 0.043 |
| 14/6/10 | 0.047 | 0.048 | 0.047 | 0.048 | 0.043 | 0.043 |
| 16/6/10 | 0.048 | 0.046 | 0.047 | 0.048 | 0.045 | 0.045 |

### 8/6/10 here against docs/row12_power_sim.md (SYNTHETIC; different seeds, two-sided Fisher p)
| Rise | Power here | Power, row 12 doc | p | False pass here | False pass, row 12 doc | p |
|---|---|---|---|---|---|---|
| 10 ms | 807/1000 = 0.807 (0.781–0.831) | 209/250 = 0.836 (0.784–0.880) | 0.32 | 30/1000 = 0.030 (0.020–0.043) | 12/250 = 0.048 (0.025–0.082) | 0.17 |
| 20 ms | 650/1000 = 0.650 (0.620–0.680) | 167/250 = 0.668 (0.606–0.726) | 0.60 | 10/1000 = 0.010 (0.005–0.018) | 6/250 = 0.024 (0.009–0.052) | 0.11 |
| 30 ms | 361/1000 = 0.361 (0.331–0.392) | 100/250 = 0.400 (0.339–0.464) | 0.27 | 2/1000 = 0.002 (0.000–0.007) | 0/250 = 0.000 (0.000–0.015) | 1.00 |
| 40 ms | 202/1000 = 0.202 (0.178–0.228) | 59/250 = 0.236 (0.185–0.294) | 0.26 | 0/1000 = 0.000 (0.000–0.004) | 0/250 = 0.000 (0.000–0.015) | 1.00 |
| 80 ms | 8/1000 = 0.008 (0.003–0.016) | 3/250 = 0.012 (0.002–0.035) | 0.47 | 0/1000 = 0.000 (0.000–0.004) | 0/250 = 0.000 (0.000–0.015) | 1.00 |

### Claim checks (numbers quoted in docs/sd_stress.md; SYNTHETIC)
- False pass at 68.67 ms, with bias: 30 cells, point 0.000-0.036, CP upper bound 0.0037-0.0495; upper bound above 0.05 in 0 of 30 cells.
- False pass at 68.67 ms, bias removed: 30 cells, point 0.038-0.048, CP upper bound 0.0518-0.0631; upper bound above 0.05 in 30 of 30 cells.
- False pass at 68.67 ms, no detector: 6 cells, point 0.039-0.045, CP upper bound 0.0529-0.0598; upper bound above 0.05 in 6 of 6 cells.
- Power at 16, bias removed against no detector: largest gap 0.013 (12/4/10, 40 ms rise).
- Noncentral-t count over fitted count, with bias: 24 cells, +1.5% to +9.8%.
- Noncentral-t count over fitted count, bias removed: 30 cells, +1.1% to +6.3%.
- Noncentral-t count over fitted count, no detector: 6 cells, +0.6% to +6.4%.
- Detector error, largest difference of any config from 8/6/10 over all rises and both true leads: mean error 0.006 ms, median error 0.013 ms, error SD 0.004 ms.

### Checks
- Vectorised gate vs `LeadStats.summarize`: 0 mismatches in 120000 16-subject studies.
- Subjects with no usable trial: 0.
- Largest missed-trial fraction of a cell: 0.00018.
- Bootstrap resamples whose fitted curve did not rise (counted as infinite): 0.
- Probit fit chi2 over 18 df (points share subjects, so indicative): p < 0.01 in 0 of 54 fits.

## Sensitivity to the margin (SYNTHETIC)
The protected margin is 6 ms (true lead 74.67 ms against the 68.67 ms gate). `python scripts/sd_stress.py sensitivity` re-scores the stored chunks with the detected leads shifted by margin - 6 ms (output `results/sd_stress/sensitivity.md`; nothing is re-simulated, and for the detected rows it is an approximation, see the script). Power at 16 subjects, 10 ms rise unless stated, 95% Clopper-Pearson:

| SD | Lead | 5 ms | 5.5 ms | 6 ms | 6.5 ms | 7 ms | 7.5 ms | 8 ms | 9 ms |
|---|---|---|---|---|---|---|---|---|---|
| 12/4/10 | 10 ms rise, with bias | 424/1000 = 0.424 (0.393–0.455) | 495/1000 = 0.495 (0.464–0.526) | 569/1000 = 0.569 (0.538–0.600) | 626/1000 = 0.626 (0.595–0.656) | 686/1000 = 0.686 (0.656–0.715) | 752/1000 = 0.752 (0.724–0.778) | 796/1000 = 0.796 (0.770–0.821) | 871/1000 = 0.871 (0.849–0.891) |
| 12/4/10 | 20 ms rise, with bias | 297/1000 = 0.297 (0.269–0.326) | 365/1000 = 0.365 (0.335–0.396) | 418/1000 = 0.418 (0.387–0.449) | 471/1000 = 0.471 (0.440–0.502) | 549/1000 = 0.549 (0.518–0.580) | 612/1000 = 0.612 (0.581–0.642) | 671/1000 = 0.671 (0.641–0.700) | 789/1000 = 0.789 (0.762–0.814) |
| 12/4/10 | no detector | 451/1000 = 0.451 (0.420–0.482) | 537/1000 = 0.537 (0.506–0.568) | 610/1000 = 0.610 (0.579–0.640) | 659/1000 = 0.659 (0.629–0.688) | 725/1000 = 0.725 (0.696–0.752) | 780/1000 = 0.780 (0.753–0.805) | 811/1000 = 0.811 (0.785–0.835) | 891/1000 = 0.891 (0.870–0.910) |
| 8/6/10 | 10 ms rise, with bias | 667/1000 = 0.667 (0.637–0.696) | 738/1000 = 0.738 (0.710–0.765) | 807/1000 = 0.807 (0.781–0.831) | 858/1000 = 0.858 (0.835–0.879) | 906/1000 = 0.906 (0.886–0.923) | 936/1000 = 0.936 (0.919–0.950) | 953/1000 = 0.953 (0.938–0.965) | 974/1000 = 0.974 (0.962–0.983) |
| 8/6/10 | no detector | 704/1000 = 0.704 (0.675–0.732) | 775/1000 = 0.775 (0.748–0.801) | 839/1000 = 0.839 (0.815–0.861) | 886/1000 = 0.886 (0.865–0.905) | 923/1000 = 0.923 (0.905–0.939) | 953/1000 = 0.953 (0.938–0.965) | 967/1000 = 0.967 (0.954–0.977) | 985/1000 = 0.985 (0.975–0.992) |

- **The design value is close to 0.8.** At 8/6/10, 10 ms rise, with bias, power is 0.807 (0.781-0.831) at the 6 ms margin; its interval includes values below 0.8, and at a 5.5 ms margin it is 0.738 (0.710-0.765). This PR does not claim the design reaches 0.8.
- **At 12/4/10 the conclusion "16 subjects do not reach 0.8" depends on the 6 ms margin.** With bias, 10 ms rise, power is 0.796 (0.770-0.821) at an 8 ms margin and 0.871 at 9 ms; with no detector 0.811 at 8 ms. The point estimate crosses 0.8 between 8 and 9 ms in this grid (margins 6.5, 7, 7.5, 8 and 9 ms were run).

## Reading the results
- **The 8/6/10 design reproduces the row-12 doc.** At 8/6/10 here, power at 16 subjects is 807/1000 (10 ms rise) against 209/250 in the row-12 doc, and no cell of the five-row table differs (two-sided Fisher p >= 0.11 over the 10 comparisons above). Subjects for 0.8 power at 10 ms rise: 16.4 (15.6-17.5), against 14.7 (12.4-16.8) there and 17 analytically.
- **At 12/4/10 the planned 16 subjects do not reach 0.8 at any rise time.** With bias, power is 0.569 (10 ms rise), 0.418, 0.223, 0.127 and 0.010 (80 ms). Even with no detector at all (the gate on the true leads) it is 0.610 (0.579-0.640). Subjects for 0.8: 28.8 (27.4-32.1) at 10 ms, 43.9 (40.5-49.1) at 20 ms, 97.6 (88.6-117.6) at 30 ms, 203.1 (171.9-257.1) at 40 ms; not reachable at 80 ms, where the mean detected lead is 1.87 ms below the gate. With bias removed it is 26.6 (25.3-29.6), close to the 28 of the session-4 analytic value and to the noncentral-t 28 here.
- **Subject-SD sweep (session and trial at 6/10).** Power at 16 subjects, bias removed: 0.832, 0.703, 0.582, 0.473, 0.406 at subject SD 8, 10, 12, 14, 16 ms. Subjects for 0.8, bias removed (10 ms rise): 15.1, 20.8, 27.7, 36.5, 46.5. With bias at 10 ms rise: 16.4, 22.6, 30.6, 39.4, 50.3. The count grows about with the variance of a subject mean (realised 8.93 to 16.53 ms), not with the subject SD alone.
- **Detector bias is what makes rise time matter.** The mean error is -0.26, -1.39, -3.01, -3.95 and -7.95 ms at 10, 20, 30, 40 and 80 ms rise, the same in every config. With the bias removed (an oracle reference, see Method), power and subject counts are the same at every rise time, to within 0.01 in power and 0.5 subjects (largest at 46 subjects), and within 0.013 of the no-detector power (largest gap 12/4/10, 40 ms rise; see Claim checks). So for any SD the rise-time loss is the detector bias, and the SD sets the baseline.
- **False pass: shown below 0.05 with bias; only consistent with 0.05 once the bias is removed.** With bias it is 0.000-0.036 and the Clopper-Pearson upper bound is 0.0495 or less in all 30 cells, so a rate at or below 0.05 is shown (reading short is conservative at the 68.67 ms boundary). With the bias removed it is 0.038-0.048, but the upper bound is 0.0518-0.0631 in **all 30 cells** (no detector: 0.0529-0.0598 in all 6), at every SD including 8/6/10. So those results are consistent with the nominal 0.05 at every SD, and do not show a rate at or below 0.05. Source: Claim checks in `tables.md`. The session-4 argument that the t-gate holds about 5% for any SD is analytic and is only supported, not proved, here.
- **Subject counts that need a model.** The count for 0.8 power comes from a probit-in-sqrt(n) fit (the row-12 method). Where it is above the largest grid size it would be extrapolated; none is (largest 335.1, grid to 512). At the fitted ceil(n) the measured power is 0.75-0.85 in every cell (table "Power measured at the fitted subject count" in `tables.md`), and the noncentral-t counts from the measured mean and SD are +0.6% to +9.8% above the fit (largest: 368 against 335 at 16/6/10, 40 ms rise; see Claim checks).
- **Caveats.** (1) One burst shape (smoothstep, SNR 15 x U(0.7, 1.3)), one motion shape (minimum jerk), one set of detector settings; nothing here says how the real distribution of rise times looks. (2) The SDs are inputs: this run says what each would cost, not which is true. (3) Counts above 100 subjects are far beyond the planned design and rest on the fit's tail. (4) Cells are paired through shared draws, so differences between cells are less noisy than their separate intervals suggest, and the cells cannot be treated as independent replicates. (5) The subject sweeps for n other than 16 use disjoint cuts of the same 16000 subjects, so a larger n means fewer studies (down to 31 at n = 512; 47 at the fitted n of the 16/6/10, 40 ms rise cell).
- **Early IMU onsets.** In every cell (the draws are shared) 13 of 480,000 trials have a detected lead more than 100 ms too short (IMU fired on noise) and 1 has one more than 100 ms too long; 87 trials are missed. All are kept, as in the row-12 doc.

## Reproduce
```
python scripts/sd_stress.py run --raw raw-data --procs 6        # 40 restartable chunks of 25 studies; about 16 h wall on 6 shared cores
python scripts/sd_stress.py score --raw raw-data --out results/sd_stress/scored.json --procs 6    # about 3 min
python scripts/sd_stress.py report results/sd_stress/scored.json --out results/sd_stress/tables.md
```
`run` defaults to the pre-registered 1000 studies, ids 1000-1999; `--studies` and `--first` exist for smoke tests only. numpy 2.5.5, scipy 1.18.1, Python 3.12.16.
