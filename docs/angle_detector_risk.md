# Angle-detector risk: the "drop past 10°" score of B-A15 / B9 (SYNTHETIC)

**Status: draft for Fact Checker, fix round 3 (FC round 1 = PARTIAL, round 2 = PASS with minors; this round fixes the minors). Every number here is SYNTHETIC: model and Monte Carlo output, not bench or human data. No trial was run and no recorded data was read.**

**Headline (rounds 2 and 3): the forearm-angle requirement (joint − pivot, section 8 item 1) FAILS PR #9's per-trial thresholds when scored with PR #9's chain.** The drop-at-200 ms error reaches -0.02630° against the 0.02627° guard (over by 3.1e-05°, at every k ≠ 1). The crossing error reaches 0.16094 ms against U_t 0.15999 ms (over by 9.6e-04 ms), and this crossing failure is **not specific to one k**: 4 of 41 k values across the 46.3 mm range fail (worst k 1.2257), and the lattice-aware bound for two 17-bit channels is 0.16135 ms. The thresholds are not loosened here; options are in section 10.1 and none is chosen.

- Every number comes from `scripts/angle_detector_risk.py` (output `results/angle_detector_risk.json`). This file is rendered from `scripts/angle_detector_risk.template.md`; do not edit it by hand. `python scripts/angle_detector_risk.py --check` recomputes everything and compares.
- **ASSUMPTION** marks anything without a source; the list is in section 9. Inputs quoted from PR #9 (`docs/ba15-testplan.md`, `scripts/ba15_numcheck.py`, `docs/ba15_numbers.json`; merged to main at `335d4e5` on 2026-11-26 10:49 KST, merge commit `a7b41e0` (GitHub API `merged_at` 01:49:42Z); every quoted input is identical to `58cce8d`, which round 1 cited) and from the Rev C package (`exo/EMG_Exo_Package_RevC.md` §3.7 and BOM, not in this repo) are marked as quoted.
- Constraints (Charlie): synthetic only, force assist only, no exact parts. Requirements below are stated as numbers, not as parts.

## What the "angle detector" is

The joint-angle chain that turns a physical arm drop into the B-A15 trial score (PR #9 §3, "Trial score"): angle sensor → sampling → zero reference → peak drop, 10° crossing time and catch time → timestamps against fault onset → DROP / NO DROP / INVALID. It verifies B9 (drop ≤ 10°). It does not act on the device.

- A **missed detection** (a true drop past 10° scored NO DROP) is a false pass of a safety test.
- A **false trip** (a spurious DROP) fails the whole test under the zero-drop rule (0 in 59, Charlie 2b).
- Not studied: the repo's `ImuOnsetDetector` (EMG-lead reference onset). It shares the name, not the role.

## Short version (all SYNTHETIC)

1. **The PR #9 requirement set meets its own budget, with thin room.** At 15 kHz, 17 bits, 0.019 ms latency and ±0.057 ms skew, the 10° crossing time is off by at most 0.14467 ms (exact worst case over 2001 sampling phases × 1001 mount phases × both skew edges; the round-1 figure 0.140 ms was a 400-trial maximum) against U_t 0.15999 ms, room 0.0153 ms on this trajectory. That room depends on where 10° falls on the count grid; without that luck (every timing term at its edge plus one count) the error can reach 0.1595 ms, room 0.0005 ms. The drop at 200 ms is off by at most 0.02355° (exact; guard 0.02627°, room 0.0027°). The peak drop is off by -0.0027° to +0.0026°. No missed detection was seen at any true drop above 10° in the grid.
2. **The score must use the forearm angle, not the device joint angle.** Under the P+R interface the device angle moves k times the forearm angle. At the Rev C offset limit of 46.3 mm (k 0.734–1.266), a device reading of 10° is a forearm drop of 7.90–13.60°; at 100 mm it is 6.31–24.00°. Joint minus pivot angle removes this: its peak error is bounded by two counts, ±0.0055° (observed -0.0024° to +0.0031° over 4000 trials). But scored through PR #9's chain it fails PR #9's thresholds (headline above).
3. **The time rule scores arms that never pass 10° as DROP.** "Catch time" is arrest complete (speed zero), but the 201.90 ms edge is the 10° crossing of an arm stopped instantly. With a finite arresting torque (ASSUMPTION, 15–50 N·m), a trial whose arrest completes just inside the edge has a true peak of only 4.69° to 8.39°. A true peak of exactly 10° completes its arrest at 213.3–273.7 ms, after the edge. So the time rule tests a stricter quantity than B9's 10°; PR #9 A9 (how time and angle combine) needs this.
4. **Catch time cannot be measured to U_t from the angle.** Near speed zero the angle is flat. With 17 bits and no noise, the argmax time is early by up to 2.55 ms; a 15-sample speed-zero rule is off by up to 3.31 ms; with 1 count of noise the errors reach 236 ms. Meeting U_t would need 22–26 bits and no noise.
5. **Gain and zero errors are not in the 0.0263° budget.** A gain error of -0.264% alone, or a zero-reference error of -0.0266° alone, makes a true 10.002° drop read NO DROP. A −0.5 % gain error missed 100% of true 10.01° drops, and no plausibility check flags it.
6. **Noise never causes a miss here but widens the false-trip band.** The peak is a maximum, so noise biases it up. With 1 count of noise every arm at 0.03° below the limit or closer is scored DROP; with 4 counts, every arm within 0.05°. A 3-sample median brings the 4-count bias from +0.0292…+0.0570° to +0.0184…+0.0389°.
7. **A stuck reading is the main missed-detection fault, and a frozen-count check catches it when the reading has noise.** These rates are *conditional*: they assume the freeze happens (no fault rate is known), at a uniform time in [0, t_c], on a true 10.01° peak with T_h 30 N·m (all ASSUMPTION). Under that model one channel frozen during the drop misses 98.0% (95 % CI 97.5%–98.4%, n 4000). Motion plus deceleration (50,000 °/s²) checks leave 26.4% (25.0%–27.8%); over the sensitivity cases 13.9%–35.7%. With the deceleration limit just above the clean maximum seen in round 2 (38500 °/s², a draw-dependent grid value; section 10.4) they leave 9.3%–25.2%. Adding a frozen-count check (no run of identical raw counts longer than 10 ms) leaves 0.0% (CI upper 0.09%) with 1 count of noise. That check needs at least about 0.5 count of noise; with none, a held clean arm is "frozen" too (section 10.3). Round 1's "97.5 % / 26.2 %" was n 400 (±4.4 points) under the same conditions. A second reading is an option, not shown to be needed (section 10.2).
8. **The Rev C joint-encoder requirement (0.1°) and the repo IMU default (150 Hz) both fail the budget.** At 0.1° taken as the quantization step (the Rev C BOM says "absolute, ≤ 0.1 deg" without saying whether that is resolution or accuracy; if it is an accuracy bound, a 0.1° error is 3.8 times the 0.0263° guard, so the requirement fails either way): crossing time off by 0.609 ms, drop at 200 ms off by 0.089°. At 150 Hz: crossing time off by 6.66 ms (margin 1.905 ms), drop at 200 ms off by 1.045°. Defaults are quoted, not changed.
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

- Room left at the requirement chain, 400-trial maximum (round 1): 0.0201 ms of U_t and 0.0027° of the guard. **Exact worst case (round 2):** 0.14467 ms and 0.02355°, room 0.0153 ms and 0.0027°. Both rooms rely on where 10° and the 200 ms drop fall on the count grid; the lattice-free crossing error is 0.1595 ms (room 0.0005 ms), which is PR #9's linear budget by design.
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

- **Model scope (ASSUMPTION):** planar, device and elbow axes parallel, rigid cuffs and links, offset fixed during the drop. Out-of-plane axis tilt, elbow-axis wander and cuff or soft-tissue motion are not modelled. Rev C itself marks its 100 mm offset source as not verified for the elbow.
- Geometry quoted from Rev C: pivot height 54 mm, s_mid 175 mm, link 0.183 m; d = √(L² − (h_B − e_y)²), s_B = e_x + d, k = s_B/d. In the script the device link angle is ψ = φ + α with sin α = (h_B − e_y′)/L, which gives dψ/dφ = 1 + e_x′/d = s_B/d = k: the device angle moves k times the forearm angle. Numerical check at one offset: k 1.179076555, dψ/dφ 1.179076555.
- Disc sweep, 41 radii × 721 directions (as Rev C). The 100 mm k_max reproduces Rev C's k_hi 1.5862 (1.5862286755).

| Offset bound | k range | Device change for a 10° forearm drop (°) | Forearm drop when the device reads 10° (°) |
|---|---|---|---|
| 46.3 mm (Rev C pass bound) | 0.7345–1.2655 | 7.348–12.652 | 7.903–13.605 |
| 100 mm, all offsets | 0.4138–1.5862 | 4.143–15.857 | 6.305–24.004 |
| 100 mm, carriage inside the 100–250 mm rail | 0.5473–1.5836 | 5.227–15.856 | 6.308–19.930 |

- Through the chain at a true 10.000° forearm peak (46.3 mm extremes, linear k): joint angle only, peak error -2.655° at k_lo and +2.655° at k_hi. Joint minus pivot angle (two 17-bit channels): bound ±2 counts = ±0.0055°; observed -0.0024° to +0.0031° over 4000 trials (round 1 quoted a 200-trial range, -0.0024° to +0.0031°, as if it were the bound).
- On the dummy (Phase A) the offset is set by the build; in people (Phase B) it is not known, which is why Rev C added the carriage sensor and pivot angle encoder for the torque cap.

## 7. Encoder faults and plausibility checks

Round 1 (this section) injects faults on one channel only. Faults on both channels of the joint − pivot reading, common-mode faults, a frozen-count check and the stuck-fault sensitivity are in section 10.2–10.4.

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

1. Score the **forearm** angle (device joint angle combined with the pivot angle), never the device joint angle alone. Log both raw channels. **With PR #9's chain this FAILS PR #9's thresholds** (drop-at-200 -0.02630° vs 0.02627°; crossing up to 0.16094 ms vs 0.15999 ms, failing at 4 of 41 k across the 46.3 mm range). Something in section 10.1 has to change; Charlie decides which.
2. Keep PR #9's chain (≥ 15 kHz, ≤ 0.0031° step, ≤ 0.019 ms latency, skews) and add to the angle budget: **gain/linearity error over the 0–15° span** and **zero-reference error**. Alone, each uses the whole guard at -0.264% and -0.0266°. Calibrate the angle scale on the dummy against a reference.
3. **Angle noise ≤ 1 count RMS** at the scoring rate, or a stated filter (3-sample median at least) whose upward bias is added to the false-trip analysis.
4. **Do not score the catch time from speed zero of the angle.** Either drop the time rule (the angle rule is the B9 quantity), or redefine it as a quantity the angle can measure at speed (for example the 10° crossing time) and recompute its window. FC and Charlie decide (PR #9 A9 is already open).
5. **Plausibility checks on every channel (round 2 wording: option set, not a proven need).** A frozen-count check plus a deceleration limit just above the clean maximum caught every stuck reading in this model when the reading has ≥ 0.5 count of noise (section 10.2–10.4). **Option:** a second, independent forearm reading with a cross-check. Among the faults modelled it is the only check that sees slow faults: backlash, drift and, partly, gain. Its own coverage, common-cause faults and rig fit are open (section 10.2). Round 1 called it "needed"; that was not shown.
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
| A9 | Record length 0.40 s after onset (the noise bias of the maximum grows with the hold length) and a glitch window of 0.3 s after onset. |
| A10 | Arresting torque 15–50 N·m: **ASSUMPTION, no source.** Section 5, finding 3 and the stuck sensitivity scale with it. |
| A11 | Round 2 faults: one at a time; sizes (5 ms holds at 20 per s, 0.05° / 0.2° slip, 0.05° / 0.2° backlash, 0.2 °/s drift, −0.5 % common gain) are illustrations. Slip and drift take the sign that makes the drop read low. Backlash starts at a random place in the play. |
| A12 | Frozen-count limit 10 ms; cross-check tolerance 0.05°; the second reading is a 17-bit forearm angle with its own noise and mount, sharing only the injected common fault. |
| A13 | Two-channel budget check: free-fall governing cell, skew fixed at the edge that makes times read late, 200000 random phase and mount draws per row. An MC maximum is a lower bound on the worst case, so a FAIL is a found counterexample, but a pass is not proof. |

## 10. Fix round 2 (FC round 1 on PR #14 @ 1df8a94)

### 10.1 Forearm angle (joint − pivot) against PR #9's thresholds: FAILS

Two independent 17-bit channels, PR #9's chain otherwise unchanged (15 kHz, latency 0.019 ms, skew at the edge that reads late), free-fall governing cell. **Exact columns** are a grid brute force over 2001 sampling phases x 128 joint mounts x 128 pivot mounts; a grid maximum is a lower bound of the true supremum, within the grid spacing. **MC column** is the round-2 random draw maximum, 200000 trials (also a lower bound).

| Chain | Exact worst crossing error (ms) vs U_t 0.15999 | Result | Exact worst drop-at-200 error (°) vs −0.02627 | Result | MC crossing (ms) |
|---|---|---|---|---|---|
| joint − pivot, k 0.7345 (exact k_lo = 0.7344681, 46.3 mm low end) | 0.15713 | pass | -0.02630 | **FAIL** | 0.15662 |
| joint − pivot, k 0.7345 (rounded; round 1 and FC) | 0.16006 | **FAIL** | -0.02630 | **FAIL** | 0.15979 |
| joint − pivot, k 1.0001 | 0.15513 | pass | -0.02630 | **FAIL** | 0.15519 |
| joint − pivot, k 1.2655 (exact k = 1.2655319, 46.3 mm high end) | 0.16036 | **FAIL** | -0.02630 | **FAIL** | 0.16025 |
| joint − pivot, k 1.2655 (rounded) | 0.14873 | pass | -0.02630 | **FAIL** | 0.14871 |
| joint − pivot, k 1.5862 (100 mm high end) | 0.16079 | **FAIL** | -0.02630 | **FAIL** | 0.15034 |
| one channel (forearm read directly) | 0.14457 | pass | -0.02355 | pass | 0.14455 |

- **Round 3 correction.** Round 2 marked the k 1.5862 crossing row "pass" at 0.15034 ms; that was an MC maximum. The grid brute force gives 0.16079 ms, a FAIL. The crossing failure is not specific to the high end of k. Over 41 k values across the 46.3 mm range (coarser grid, 401 phases x 64 x 64 mounts) the worst crossing error is 0.16094 ms at k 1.2257; 4 of 41 k fail (k = 0.7743, 1.0398, 1.2257, 1.2655). The lattice-aware bound for two 17-bit channels (every timing term at its edge, plus the time the arm needs at the 10° crossing speed to cover n+1 counts above 10°) is 0.16135 ms, above U_t; one channel: 0.14471 ms.
- The one-channel exact value here (0.14457 ms) uses 128 mount phases; section 3 quotes 0.14467 ms from 1001 mount phases.
- The drop-at-200 error is -0.02630° at every k ≠ 1: the measured drop is a multiple of q that is more than 2q below the true drop at the last sample before 200 ms, whatever k is. One channel: -0.02355°.
- The overrun is small (drop-at-200 by 3.1e-05°, crossing by 9.6e-04 ms). But the thresholds are hard per-trial pass lines, and the second channel's step is not in PR #9's budget. **This is a failing result.** PR #9's thresholds are not loosened here.
- **Why FC's round-1 run read 0.16004 ms at k_lo while round 2's run read 0.15662 ms there (and why the rows differ).** It is not random draws alone. The crossing supremum depends on where k·d and (k−1)·d fall on the count grid, so a change of k of 3·10⁻⁵ can flip a row. FC's round-1 run used the rounded k 0.7345, where the grid supremum is 0.16006 ms (FAIL; my MC at 0.7345: 0.15979 ms). Round 2 used the exact k_lo = 0.7344681, where the supremum is 0.15713 ms (pass; round-2 MC 0.15662 ms). The 0.16025 ms of round 2 is a different row, k_hi = 1.2655319 (grid supremum 0.16036 ms). All MC values here are lower bounds of their suprema. The k-sensitivity of the supremum near k_lo (grid 2001 sampling phases x 128 joint mounts x 128 pivot mounts):

| k | Exact worst crossing error (ms) | Result |
|---|---|---|
| k_lo +0.0e+00 = 0.7344681 | 0.15713 | pass |
| k_lo +1.0e-06 = 0.7344691 | 0.15736 | pass |
| k_lo +1.0e-05 = 0.7344781 | 0.15939 | pass |
| k_lo +3.2e-05 = 0.7345000 | 0.16006 | **FAIL** |
| k_lo +1.0e-04 = 0.7345681 | 0.15443 | pass |

- Options that pass (exact columns at k_lo / 1.0001 / k_hi, grid 1001 sampling phases x 96 joint mounts x 96 pivot mounts; worst over the three k; the last column is the lattice-aware bound for that option). **None is chosen**; each is for HE / Charlie:

| Option | Worst crossing error (ms) | Worst drop-at-200 error (°) | Result | Lattice-aware bound (ms) |
|---|---|---|---|---|
| pivot channel 18 bits (joint 17) | 0.16029 | -0.02630 | **fails** | 0.16135 |
| both channels 18 bits | 0.15136 | -0.02493 | passes | 0.15303 |
| sampling 20 kHz | 0.14366 | -0.02355 | passes | 0.14468 |
| angle-path latency 0 (bound case) | 0.14124 | -0.02355 | passes | 0.14230 |
| linear interpolation of crossing and drop at 200 ms | 0.10750 | -0.01741 | passes | n/a |

  The MC values for the same options (round 2) were: pivot channel 18 bits (joint 17) 0.16012 ms / -0.02630°; both channels 18 bits 0.15135 ms / -0.02493°; sampling 20 kHz 0.14366 ms / -0.02355°; angle-path latency 0 (bound case) 0.14101 ms / -0.02355°; linear interpolation of crossing and drop at 200 ms 0.10728 ms / -0.01742°.
- Loosening option (Charlie's call only, not proposed): add one pivot count to each budget. The crossing figure depends on the speed used for the count: U_t 0.1771 ms on PR #9's own basis (`w10_min` 160.373 °/s, the basis of PR #9's angle term), or 0.1766 ms at the governing cell's 165.046 °/s. Guard 0.0290° on either basis.
- "Latency 0" is a bound case showing that latency is one of the terms that matters, not a buildable value. Interpolation changes the scoring rule, which is PR #9's. A pass here is a grid result (A13).

### 10.2 Faults on both channels, common-mode faults, and checks

Forearm = joint − pivot at k 0.7345 (the 46.3 mm low end, where a frozen pivot under-reads), true peak 10.01°, T_h 30 N·m, 1 count of noise per channel, 2000 trials per row. Checks are cumulative left to right: base = motion + deceleration limit 58000 °/s² (chosen just above the largest clean two-channel deceleration seen, 57678 °/s², which is a draw-dependent maximum on a grid; see the next bullet); + frozen count = no run of identical raw counts on either channel longer than 10 ms; + cross-check = second reading within 0.05°.

| Fault | Missed, no checks | Missed after base | After + frozen count | After + cross-check | Flagged INVALID (all checks) |
|---|---|---|---|---|---|
| No fault | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| Joint stuck, random time during the drop | 97.6% | 40.6% | 0.0% | 0.0% | 100.0% |
| Pivot stuck from onset | 100.0% | 100.0% | 0.0% | 0.0% | 100.0% |
| Pivot stuck during the drop | 95.2% | 95.2% | 0.0% | 0.0% | 100.0% |
| Joint partial sticking (5 ms holds, 20 per s) | 0.0% | 0.0% | 0.0% | 0.0% | 97.7% |
| Pivot partial sticking (5 ms holds, 20 per s) | 0.0% | 0.0% | 0.0% | 0.0% | 89.3% |
| Joint coupling slip 0.05° | 46.7% | 0.9% | 0.9% | 0.8% | 98.2% |
| Joint coupling slip 0.2° | 97.8% | 0.0% | 0.0% | 0.0% | 100.0% |
| Pivot coupling slip 0.05° | 46.6% | 1.0% | 1.0% | 0.9% | 97.8% |
| Joint backlash 0.05° | 1.3% | 1.3% | 1.3% | 1.3% | 0.0% |
| Joint backlash 0.2° | 75.8% | 75.8% | 75.8% | 1.8% | 74.0% |
| Joint zero drift 0.2 °/s | 96.2% | 96.2% | 96.2% | 57.9% | 38.4% |
| Pivot zero drift 0.2 °/s | 95.5% | 95.5% | 95.5% | 57.3% | 38.2% |
| Common: acquisition freeze, both device channels | 98.2% | 27.0% | 0.0% | 0.0% | 100.0% |
| Common: acquisition freeze incl. the second reading | 97.7% | 26.1% | 0.0% | 0.0% | 100.0% |
| Common: gain −0.5 % on both device channels | 49.9% | 49.9% | 49.9% | 48.5% | 1.5% |
| Pivot stuck from onset, no noise ¹ | 100.0% | 100.0% | 0.0% | 0.0% | 100.0% |
| Joint stuck during the drop, no noise ¹ | 97.8% | 47.0% | 0.0% | 0.0% | 100.0% |

¹ With no noise the frozen-count check also flags every clean trial (the held arm reads constant; 10.3), so the "+ frozen count" columns are not usable for these two rows.

- **The limit is draw-dependent.** The measured deceleration is a second difference of quantized angle over 1 ms windows, so it takes values on a grid of 2746.58 °/s². A larger clean pool (round 3): 24000 clean trials; largest measured deceleration 57678.2 °/s² = 21.00 grid units of 2746.58 °/s²; limit 58000 °/s² = 21.12 units; trials over the limit 0, 95 % upper bound on the false-INVALID rate 0.015%. A maximum over a finite pool is not a bound: another set of draws can reach one grid unit higher (on the one-channel chain a larger pool did reach 41,199 °/s², one unit above the round-2 maximum, as FC's round-2 draws had; section 10.4). The false-INVALID rate at the limit is small but not zero.
- Clean trials, all four arrest torques × k 0.7345 / 1.0001 / 1.2655, 12000 trials: flagged 0.0% (base), 0.0% (+ frozen count), 0.0% (+ cross-check). Longest clean run 1.13 ms.
- Pivot frozen from onset at k 1.2655 over-reads: a true 9.95° arm is scored DROP 100.0% of the time (a false trip); the frozen-count check flags 100.0%.
- What the rows show: a frozen pivot passes the round-1 checks every time (as FC found). The frozen-count check catches stuck readings on either channel, including a common acquisition freeze, but only when the reading has noise. Slips of 0.05° or more trip the deceleration check. Partial sticking is flagged but causes no miss here, because the arm is held at its peak. Backlash, drift and common gain are slow, and no check except the cross-check sees them. The cross-check cuts backlash 0.2° misses from 75.8% to 1.8% and drift from 96.2% to 57.9%. It does not help with a −0.5 % common gain (48.5%), because the under-read is about its 0.05° tolerance.
- **Second reading, open:** its own accuracy on a person (cuff and soft tissue), its tolerance budget, what a disagreement scores (INVALID or DROP), its common-cause links (shaft or coupling, supply, clock, acquisition), and whether it fits on the rig are not modelled. No part is named (Charlie's constraint).
- **Not modelled:** encoder nonlinearity beyond a scale error; electromagnetic bursts; counter roll-over or wrong sign; clock drift between channels; partial supply loss; joint-encoder placement on the motor side of the gearbox (Rev C BOM does not say); out-of-plane misalignment.

### 10.3 Frozen-count check against noise (one channel, k 1, T_h 30, 1000 clean and 1000 stuck trials per row)

| Noise (counts RMS) | Clean longest run (ms) | Stuck shortest run (ms) | Clean flagged at 10 ms | Fully separated |
|---|---|---|---|---|
| 0 | 177.60 | 177.1 | 100.0% | no |
| 0.25 | 15.33 | 176.7 | 0.9% | yes |
| 0.5 | 2.13 | 176.4 | 0.0% | yes |
| 1 | 0.80 | 176.7 | 0.0% | yes |

### 10.4 Stuck reading: assumptions and sensitivity (one channel, as round 1; n 4000 per row, 95 % Clopper-Pearson)

Conditions: the freeze always happens (no fault rate), at a uniform time over the window shown, on one true peak and one arrest torque per row. These are ASSUMPTIONS; the rates are conditional on a freeze.

| Case | Missed, no checks | After motion + decel 50,000 | After motion + decel 38500 | Share of the 50,000-limit misses reading < 5° |
|---|---|---|---|---|
| Base: δ 0.01°, T_h 30, freeze U[0, t_c], no noise | 98.0% (97.5%–98.4%) | 26.4% (25.0%–27.8%) | 18.1% (17.0%–19.4%) | 82% |
| Base with 1 count noise | 98.2% (97.8%–98.6%) | 23.4% (22.1%–24.8%) | 14.1% (13.0%–15.2%) | 84% |
| δ 0.1° | 95.8% (95.2%–96.4%) | 23.2% (21.9%–24.5%) | 14.9% (13.8%–16.0%) | 89% |
| δ 1° | 89.6% (88.6%–90.5%) | 19.7% (18.5%–20.9%) | 12.7% (11.7%–13.7%) | 100% |
| δ 3° | 82.7% (81.5%–83.9%) | 20.3% (19.1%–21.6%) | 13.8% (12.8%–14.9%) | 100% |
| T_h 15 N·m | 97.2% (96.6%–97.7%) | 35.7% (34.2%–37.2%) | 25.2% (23.8%–26.6%) | 47% |
| T_h 20 N·m | 98.1% (97.7%–98.5%) | 30.2% (28.8%–31.7%) | 20.2% (19.0%–21.5%) | 68% |
| T_h 50 N·m | 98.5% (98.1%–98.9%) | 25.2% (23.9%–26.6%) | 17.0% (15.8%–18.2%) | 93% |
| Freeze U[0, 10° crossing] | 99.2% (98.8%–99.4%) | 27.7% (26.3%–29.1%) | 18.4% (17.2%–19.6%) | 83% |
| Freeze U[0, end of record] | 54.3% (52.8%–55.9%) | 13.9% (12.8%–15.0%) | 9.3% (8.4%–10.2%) | 81% |

- 38500 °/s² is just above the largest clean one-channel deceleration seen, 38452 °/s² (1 count noise, four arrest torques, 4000 trials). It is a draw-dependent maximum on a grid of 2746.58 °/s²; larger clean pool (round 3): 16000 clean trials; largest measured deceleration 41198.7 °/s² = 15.00 grid units of 2746.58 °/s²; limit 38500 °/s² = 14.02 units; trials over the limit 2, 95 % upper bound on the false-INVALID rate 0.045%. The larger pool reached 41199 °/s² (15 units), one unit above the round-2 maximum, so 2 clean trials would have been flagged at 38500. **Sensitivity to the limit** (stuck during the drop, one channel with 1 count of noise, n 4000 per row; rows: base T_h 30, T_h 15, T_h 50, freeze window to end of record): after-check misses 7.5%–21.2% at 38,500; 9.1%–23.8% at 41,500 (one grid unit above the larger pool's maximum); 12.4%–32.4% at 50,000. The conclusion does not depend on which of these is used; the exact residual does. With the frozen-count check added and 1 count of noise: 0.0% missed (CI 0.00%–0.09%). Clean trials flagged by the deceleration limit: 0.0% in the round-2 pool (4000 trials), 2 of 16000 in the larger pool; the frozen-count check flagged 0.
- Most after-check misses are gross under-reads (reading < 5°), not near-10° cases.

## 11. Reproduce

    python scripts/angle_detector_risk.py --check     # about 20 min on one core since round 2; exit 0 if both files reproduce
    pytest tests/test_angle_detector_risk.py          # fast checks (model vs PR #9, geometry, chain, rendering)
