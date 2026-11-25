# Angle-detector risk: the "drop past 10°" score of B-A15 / B9 (SYNTHETIC)

**Status: draft for Fact Checker. Every number here is SYNTHETIC: model and Monte Carlo output, not bench or human data. No trial was run and no recorded data was read.**

- Every number comes from `scripts/angle_detector_risk.py` (output `results/angle_detector_risk.json`). This file is rendered from `scripts/angle_detector_risk.template.md`; do not edit it by hand. `python scripts/angle_detector_risk.py --check` recomputes everything and compares.
- **ASSUMPTION** marks anything without a source; the list is in section 9. Inputs quoted from PR #9 (`docs/ba15-testplan.md`, `scripts/ba15_numcheck.py`, `docs/ba15_numbers.json` at `58cce8d`, not merged) and from the Rev C package (`exo/EMG_Exo_Package_RevC.md` §3.7 and BOM, not in this repo) are marked as quoted.
- Constraints (Charlie): synthetic only, force assist only, no exact parts. Requirements below are stated as numbers, not as parts.

## What the "angle detector" is

The joint-angle chain that turns a physical arm drop into the B-A15 trial score (PR #9 §3, "Trial score"): angle sensor → sampling → zero reference → peak drop, 10° crossing time and catch time → timestamps against fault onset → DROP / NO DROP / INVALID. It verifies B9 (drop ≤ 10°). It does not act on the device.

- A **missed detection** (a true drop past 10° scored NO DROP) is a false pass of a safety test.
- A **false trip** (a spurious DROP) fails the whole test under the zero-drop rule (0 in 59, Charlie 2b).
- Not studied: the repo's `ImuOnsetDetector` (EMG-lead reference onset). It shares the name, not the role.

## Short version (all SYNTHETIC)

1. **The PR #9 requirement set meets its own budget, with thin room.** At 15 kHz, 17 bits, 0.019 ms latency and ±0.057 ms skew, the 10° crossing time is off by at most 0.140 ms (budget U_t 0.160 ms) and the drop at 200 ms by at most 0.0236° (budget 0.0263°). The peak drop is off by -0.0027° to +0.0026°. No missed detection was seen at any true drop above 10° in the grid.
2. **The score must use the forearm angle, not the device joint angle.** Under the P+R interface the device angle moves k times the forearm angle. At the Rev C offset limit of 46.3 mm (k 0.734–1.266), a device reading of 10° is a forearm drop of 7.90–13.60°; at 100 mm it is 6.31–24.00°. Joint plus pivot angle removes this (peak error within ±0.0031°).
3. **The time rule scores arms that never pass 10° as DROP.** "Catch time" is arrest complete (speed zero), but the 201.90 ms edge is the 10° crossing of an arm stopped instantly. With a finite arresting torque (ASSUMPTION, 15–50 N·m), a trial whose arrest completes just inside the edge has a true peak of only 4.69° to 8.39°. A true peak of exactly 10° completes its arrest at 213.3–273.7 ms, after the edge. So the time rule tests a stricter quantity than B9's 10°; PR #9 A9 (how time and angle combine) needs this.
4. **Catch time cannot be measured to U_t from the angle.** Near speed zero the angle is flat. With 17 bits and no noise, the argmax time is early by up to 2.55 ms; a 15-sample speed-zero rule is off by up to 3.31 ms; with 1 count of noise the errors reach 236 ms. Meeting U_t would need 22–26 bits and no noise.
5. **Gain and zero errors are not in the 0.0263° budget.** A gain error of -0.264% alone, or a zero-reference error of -0.0266° alone, makes a true 10.002° drop read NO DROP. A −0.5 % gain error missed 100% of true 10.01° drops, and no plausibility check flags it.
6. **Noise never causes a miss here but widens the false-trip band.** The peak is a maximum, so noise biases it up. With 1 count of noise every arm at 0.03° below the limit or closer is scored DROP; with 4 counts, every arm within 0.05°. A 3-sample median brings the 4-count bias from +0.0292…+0.0570° to +0.0184…+0.0389°.
7. **A stuck encoder is the main missed-detection fault.** Frozen from onset, it is caught by a motion check 100% of the time. Frozen at a random time during the drop, it misses 97.5% of true 10.01° drops; motion plus deceleration checks still leave 26.2% missed. Tightening the deceleration limit flags clean noisy trials first. A second, independent angle reading is needed.
8. **The Rev C joint-encoder requirement (0.1°) and the repo IMU default (150 Hz) both fail the budget.** At 0.1°: crossing time off by 0.609 ms, drop at 200 ms off by 0.089°. At 150 Hz: crossing time off by 6.66 ms (margin 1.905 ms), drop at 200 ms off by 1.045°. Defaults are quoted, not changed.
9. **Plan level.** To keep the chance of at least one false DROP at 5 % or less, the per-trial false-trip rate must be at most 8.69e-04 over 59 trials and 1.55e-04 over 330. A rate of 1e-3 gives 5.7% and 28.1%.

## 1. Drop model and checks

- Re-implementation of HE's drop model as used in PR #9 (`scripts/ba15_numcheck.py`; same inputs), governing cell (b) e-stop, ramp 21.00 N·m/s, STO 148 ms, hold 0. Torque basis: T0 and ramp device side, × k = 1.5862 at the arm; G, hold and arrest torque arm side.
- 10° crossing 201.904848 ms (PR #9: 201.904848; difference -5.2e-09 ms). Drop at 200 ms 9.688777° (PR #9 9.688777). Speed at 10° 165.046 °/s (PR #9 165.046). Drop at 100 ms 1.210°.
- **Catch model (ASSUMPTION):** from a catch-onset time the user adds a constant arresting torque T_h (15, 20, 30, 50 N·m at the arm) until speed is zero; then the arm is held at the peak. The onset is solved so the true peak is 10° + δ. No human data is behind T_h.

| T_h (N·m) | Deceleration near 10° (°/s²) | Catch onset for a 10° peak (ms) | Arrest complete (ms) |
|---|---|---|---|
| 15 | 836 | 161.05 | 273.74 |
| 20 | 1696 | 173.69 | 242.10 |
| 30 | 3416 | 184.26 | 223.59 |
| 50 | 6855 | 191.76 | 213.26 |

## 2. Measurement chain

Per trial: random sampling phase, encoder mounting phase (where the count edges fall) and skew; latency; floor quantization; optional white noise, gain, zero error, 3-sample median and faults. Baseline: mean of the samples 50 ms before onset. Peak drop: maximum after onset. Crossing: first sample at or above 10°, no interpolation (as in PR #9's budget). Drop at 200 ms: last sample logged at or before 200 ms. Score: DROP if peak + 0.0263° > 10°. 400 trials per point unless stated; seeds from one `SeedSequence(20261126)` by spawn.

Requirement chain ("req"): 15000 Hz, 17 bits (0.00275°), latency 0.0190 ms, skew uniform within ±0.0571 ms (trigger 0.0190 + release 0.0381; ASSUMPTION: release skew acts as a timing offset). All quoted from PR #9.

## 3. Crossing channel (Phase A dummy matrix, per-trial pass)

Error = measured − true. A positive crossing error or a negative drop error helps a unit pass.

| Chain | 10° crossing error (ms), min / mean / max | Drop at 200 ms error (°), min / max | Within U_t / guard |
|---|---|---|---|
| req, skew random | -0.0391 / 0.0459 / 0.1354 | -0.0236 / 0.0067 | True / True |
| req, skew at −edge | 0.0624 / 0.1037 / 0.1399 | -0.0236 / -0.0126 | True / True |
| req, skew at +edge | -0.0502 / -0.0112 / 0.0285 | -0.0071 / 0.0067 | True / True |
| 0.1° steps (Rev C joint encoder) | -0.6088 / -0.2536 / 0.0942 | -0.0888 / 0.0112 | False / False |
| 150 Hz (repo IMU default) | -0.0358 / 3.3017 / 6.6618 | -1.0453 / 0.0012 | False / False |

- Room left at the requirement chain: 0.0201 ms of U_t and 0.0027° of the angle guard (worst rows).
- At 150 Hz, linear interpolation between samples brings the crossing error to -0.098 to 0.074 ms on this smooth model trace; that is not a bound for real data.

## 4. Peak-drop channel (Phase B angle score)

P(DROP) at true peak 10° + δ, T_h 30 N·m. δ > 0 is a real drop past 10° (a NO DROP there is a miss).

| δ (°) | req | 1 count noise | 4 counts noise | 4 counts + median-3 | 0.1° steps | 150 Hz |
|---|---|---|---|---|---|---|
| -0.0500 | 0.000 | 0.000 | 1.000 | 0.660 | 0.477 | 0.000 |
| -0.0300 | 0.000 | 1.000 | 1.000 | 1.000 | 0.743 | 0.000 |
| -0.0263 | 0.307 | 1.000 | 1.000 | 1.000 | 0.743 | 0.345 |
| -0.0200 | 1.000 | 1.000 | 1.000 | 1.000 | 0.802 | 1.000 |
| -0.0100 | 1.000 | 1.000 | 1.000 | 1.000 | 0.890 | 1.000 |
| -0.0050 | 1.000 | 1.000 | 1.000 | 1.000 | 0.965 | 1.000 |
| +0.0000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| +0.0020 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| +0.0050 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| +0.0100 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| +0.0200 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| +0.0300 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| +0.0500 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

| Chain | Peak error (°), min / max | Largest miss rate (δ > 0) | Some DROP from δ = | All DROP from δ = |
|---|---|---|---|---|
| req | -0.0027 / +0.0026 | 0.000 | -0.0263 | -0.0200 |
| 1 count noise | +0.0060 / +0.0148 | 0.000 | -0.0300 | -0.0300 |
| 4 counts noise | +0.0292 / +0.0570 | 0.000 | -0.0500 | -0.0500 |
| 4 counts + median-3 | +0.0184 / +0.0389 | 0.000 | -0.0500 | -0.0300 |
| 0.1° steps | -0.0950 / +0.0980 | 0.000 | -0.0500 | +0.0000 |
| 150 Hz | -0.0027 / +0.0026 | 0.000 | -0.0263 | -0.0200 |

- The guard makes every arm within 0.0263° below 10° a possible DROP by design (PR #9: the score is meant to be strict).
- Peak error at the requirement chain over all four arrest torques: -0.0024° to +0.0021°. Sampling rate barely matters for the peak, because the speed is zero there and the arm is held.

### Single-source limits for a missed detection (true peak 10.002°, worst of 200 draws)

| Source | Limit (this chain) | Simple bound (guard + 0.002°) |
|---|---|---|
| Gain (scale) error, reading low | -0.2641% | -0.2827% |
| Zero-reference error, drop read low | -0.02658° | -0.02827° |

Resolution scan (true peaks 10.002–10.05°, worst miss rate over 13 peaks, 100 trials each):

| Bits | Step (°) | Worst miss, no noise | Worst miss, 0.5-count noise |
|---|---|---|---|
| 10 | 0.35156 | 0.59 | 0.00 |
| 11 | 0.17578 | 0.09 | 0.00 |
| 12 | 0.08789 | 0.18 | 0.00 |
| 13 | 0.04395 | 0.00 | 0.00 |
| 14 | 0.02197 | 0.00 | 0.00 |
| 15 | 0.01099 | 0.00 | 0.00 |
| 16 | 0.00549 | 0.00 | 0.00 |
| 17 | 0.00275 | 0.00 | 0.00 |

- Without noise the measured drop is a whole number of counts, so whether a miss shows up depends on where 10° falls on the count grid. Misses were seen up to 12 bits. A step below the guard needs at least 14 bits; the 17-bit requirement has room.

## 5. Catch-time channel (time rule t_c + U_t > 201.90 ms)

| T_h (N·m) | True peak when arrest completes at the edge (°) | Argmax time error, no noise (ms) | Speed-zero W = 15 error (ms) | W = 75 error (ms) | W = 1 error (ms) | Argmax, 1 count noise (ms) | W = 15, 1 count noise (ms) | Flat time for 1 count (ms) | Bits for U_t |
|---|---|---|---|---|---|---|---|---|---|
| 15 | 4.693 | -2.55 to -0.07 | -3.31 to -2.04 | -0.01 to +2.37 | -273.27 to -264.69 | -3.09 to +124.25 | -236.49 to -7.28 | 2.564 | 26 |
| 20 | 6.018 | -1.80 to -0.13 | -1.59 to -0.48 | +0.73 to +2.43 | -241.32 to -233.01 | -1.68 to +149.10 | -219.61 to -2.75 | 1.800 | 24 |
| 30 | 7.339 | -1.27 to -0.03 | -0.81 to +0.34 | +1.23 to +2.51 | -223.13 to -214.55 | -0.63 to +174.71 | -188.62 to -0.99 | 1.268 | 23 |
| 50 | 8.394 | -0.86 to +0.02 | -0.39 to +0.56 | +1.62 to +2.52 | -212.62 to -204.21 | -0.49 to +183.50 | -181.09 to -0.20 | 0.895 | 22 |

- Errors are min to max over 200 trials. "Flat time for 1 count" = √(2q/a): how long the angle stays within one count around the turning point. "Bits for U_t" = bits whose step is ½·a·U_t².
- A one-sample difference (W = 1) reads zero speed as soon as one sample has no count change, so it marks arrest at the start of the motion.
- With noise the argmax picks a noise maximum during the hold after arrest, up to 184 ms late.

## 6. Which angle: device joint vs forearm (P+R interface, Rev C §3.7)

- Geometry quoted from Rev C: pivot height 54 mm, s_mid 175 mm, link 0.183 m; d = √(L² − (h_B − e_y)²), s_B = e_x + d, k = s_B/d. In the script the device link angle is ψ = φ + α with sin α = (h_B − e_y′)/L, which gives dψ/dφ = 1 + e_x′/d = s_B/d = k: the device angle moves k times the forearm angle. Numerical check at one offset: k 1.179076555, dψ/dφ 1.179076555.
- Disc sweep, 41 radii × 721 directions (as Rev C). The 100 mm k_max reproduces Rev C's k_hi 1.5862 (1.5862286755).

| Offset bound | k range | Device change for a 10° forearm drop (°) | Forearm drop when the device reads 10° (°) |
|---|---|---|---|
| 46.3 mm (Rev C pass bound) | 0.7345–1.2655 | 7.348–12.652 | 7.903–13.605 |
| 100 mm, all offsets | 0.4138–1.5862 | 4.143–15.857 | 6.305–24.004 |
| 100 mm, carriage inside the 100–250 mm rail | 0.5473–1.5836 | 5.227–15.856 | 6.308–19.930 |

- Through the chain at a true 10.000° forearm peak (46.3 mm extremes, linear k): joint angle only, peak error -2.655° at k_lo and +2.655° at k_hi. Joint minus pivot angle (two 17-bit channels): -0.0024° to +0.0031° at k_lo.
- On the dummy (Phase A) the offset is set by the build; in people (Phase B) it is not known, which is why Rev C added the carriage sensor and pivot angle encoder for the torque cap.

## 7. Encoder faults and plausibility checks

True peak 10.01° (a real drop), T_h 30 N·m, 400 trials. Checks (ASSUMPTION): **motion**, measured drop at 100 ms at least half the model's 1.210° (needs a relaxed subject who cannot catch before 100 ms); **deceleration**, no measured deceleration above 50000 °/s² (1 ms speed windows; the largest modelled is 6855 °/s²); **sample counter**, any lost sample makes the trial INVALID (PR #9 logs the sample-drop count).

| Case | Missed (no checks) | Flagged INVALID | Missed after checks |
|---|---|---|---|
| Clean | 0.0% | 0.0% | 0.0% |
| Clean, 1 count noise | 0.0% | 0.0% | 0.0% |
| Stuck from onset | 100.0% | 100.0% | 0.0% |
| Stuck at a random time during the drop | 97.5% | 71.2% | 26.2% |
| 2 ms of samples lost at the peak, arm held | 0.0% | 100.0% | 0.0% |
| 2 ms lost at the peak, arm rebounds (ASSUMPTION) | 0.0% | 100.0% | 0.0% |
| Gain −0.5 % | 100.0% | 0.0% | 100.0% |

- Clean trials at all four arrest torques with 1 count of noise: flagged 0.0%; largest measured deceleration 38452 °/s² (noise dominates it).
- Deceleration limit trade-off (stuck during the drop vs clean noisy trials):

| Limit (°/s²) | Stuck faults missed after checks | Clean trials flagged |
|---|---|---|
| 10000 | 0.0% | 100.0% |
| 20000 | 0.0% | 100.0% |
| 50000 | 26.2% | 0.0% |
| 100000 | 50.5% | 0.0% |

- A gap at a rebounding peak under-reads by about ½·a·(D/2)²: 0.00171° for 2 ms; a gap up to 7.84 ms stays inside the guard.
- One-sample glitch of 0.1° on an arm truly at 9.95°: false trip 29.5% of trials (when it lands near the peak); with a 3-sample median 0.0%.

## 8. Requirements that follow (by requirement, no parts; SYNTHETIC basis)

1. Score the **forearm** angle (device joint angle combined with the pivot angle), never the device joint angle alone. Log both raw channels.
2. Keep PR #9's chain (≥ 15 kHz, ≤ 0.0031° step, ≤ 0.019 ms latency, skews) and add to the angle budget: **gain/linearity error over the 0–15° span** and **zero-reference error**. Alone, each uses the whole guard at -0.264% and -0.0266°. Calibrate the angle scale on the dummy against a reference.
3. **Angle noise ≤ 1 count RMS** at the scoring rate, or a stated filter (3-sample median at least) whose upward bias is added to the false-trip analysis.
4. **Do not score the catch time from speed zero of the angle.** Either drop the time rule (the angle rule is the B9 quantity), or redefine it as a quantity the angle can measure at speed (for example the 10° crossing time) and recompute its window. FC and Charlie decide (PR #9 A9 is already open).
5. **A second, independent angle reading** cross-checked every trial, plus the motion and sample-counter checks, to catch a stuck or frozen reading. Motion, deceleration and sample-counter checks together still leave 26.2% of stuck-during-drop faults missed.
6. Do not use the Rev C 0.1° joint-encoder requirement or the repo's 150 Hz IMU path for scoring.

## 9. Assumptions

| # | ASSUMPTION |
|---|---|
| A1 | PR #9's drop model and governing cell are the reference trajectory (HE / PR #9; SYNTHETIC). |
| A2 | Catch model: constant arresting torque 15–50 N·m from a catch onset until speed is zero, then held. Rebound only in the gap case. No source. |
| A3 | Truth for "drop past 10°" is the true forearm peak drop > 10°. |
| A4 | Noise is white Gaussian, in counts; faults are injected one at a time. Quantization is floor with a random mounting phase. |
| A5 | Release skew acts as a timing offset, like trigger skew; the two add to ±0.0571 ms. |
| A6 | P+R relation from Rev C §3.7 geometry; the chain test uses a constant k over the 10° drop (the disc sweep uses the exact geometry). |
| A7 | Plausibility thresholds (motion at 100 ms, deceleration 50000 °/s²) are illustrations, not settings. |
| A8 | Plan-level rates assume independent trials. |

## 10. Reproduce

    python scripts/angle_detector_risk.py --check     # about 100 s on one core; exit 0 if both files reproduce
    pytest tests/test_angle_detector_risk.py          # fast checks (model vs PR #9, geometry, chain, rendering)
