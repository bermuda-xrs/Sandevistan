# Angle-detector risk: the "drop past 10°" score of B-A15 / B9 (SYNTHETIC)

**Status: draft for Fact Checker, fix round 2 (FC round 1 = PARTIAL). Every number here is SYNTHETIC: model and Monte Carlo output, not bench or human data. No trial was run and no recorded data was read.**

**Fix round 2, headline: the forearm-angle requirement (joint − pivot, section 8 item 1) FAILS PR #9's per-trial thresholds when scored with PR #9's chain.** The drop-at-200 ms error reaches {{r2.budget.two_ch_worst_d200_deg|.5f}}° against the {{inputs.guard|.5f}}° guard (over by {{r2.budget.two_ch_d200_over_deg|.1e}}°, at every k tested), and the crossing error reaches {{r2.budget.two_ch_worst_t10_ms|.5f}} ms against U_t {{inputs.u_t_ms|.5f}} ms (over by {{r2.budget.two_ch_t10_over_ms|.1e}} ms, at k 1.2655). The thresholds are not loosened here; options are in section 10.1 and none is chosen.

- Every number comes from `scripts/angle_detector_risk.py` (output `results/angle_detector_risk.json`). This file is rendered from `scripts/angle_detector_risk.template.md`; do not edit it by hand. `python scripts/angle_detector_risk.py --check` recomputes everything and compares.
- **ASSUMPTION** marks anything without a source; the list is in section 9. Inputs quoted from PR #9 (`docs/ba15-testplan.md`, `scripts/ba15_numcheck.py`, `docs/ba15_numbers.json`; current head `335d4e5`, not merged; every quoted input is identical to `58cce8d`, which round 1 cited) and from the Rev C package (`exo/EMG_Exo_Package_RevC.md` §3.7 and BOM, not in this repo) are marked as quoted.
- Constraints (Charlie): synthetic only, force assist only, no exact parts. Requirements below are stated as numbers, not as parts.

## What the "angle detector" is

The joint-angle chain that turns a physical arm drop into the B-A15 trial score (PR #9 §3, "Trial score"): angle sensor → sampling → zero reference → peak drop, 10° crossing time and catch time → timestamps against fault onset → DROP / NO DROP / INVALID. It verifies B9 (drop ≤ 10°). It does not act on the device.

- A **missed detection** (a true drop past 10° scored NO DROP) is a false pass of a safety test.
- A **false trip** (a spurious DROP) fails the whole test under the zero-drop rule (0 in 59, Charlie 2b).
- Not studied: the repo's `ImuOnsetDetector` (EMG-lead reference onset). It shares the name, not the role.

## Short version (all SYNTHETIC)

1. **The PR #9 requirement set meets its own budget, with thin room.** At 15 kHz, 17 bits, 0.019 ms latency and ±{{inputs.skew_total_ms|.3f}} ms skew, the 10° crossing time is off by at most {{r2.exact.t10_abs_max_ms|.5f}} ms (exact worst case over {{r2.exact.n_phase}} sampling phases × {{r2.exact.n_mount}} mount phases × both skew edges; the round-1 figure {{crossing.worst_t10_abs_max_ms|.3f}} ms was a 400-trial maximum) against U_t {{inputs.u_t_ms|.5f}} ms, room {{r2.exact.u_t_room_ms|.4f}} ms on this trajectory. That room depends on where 10° falls on the count grid; without that luck (every timing term at its edge plus one count) the error can reach {{r2.exact.lattice_free_t10_ms|.4f}} ms, room {{r2.exact.lattice_free_room_ms|.4f}} ms. The drop at 200 ms is off by at most {{r2.exact.d200_abs_max_deg|.5f}}° (exact; guard {{inputs.guard|.5f}}°, room {{r2.exact.guard_room_deg|.4f}}°). The peak drop is off by {{peak_oc.req.err_min|.4f}}° to +{{peak_oc.req.err_max|.4f}}°. No missed detection was seen at any true drop above 10° in the grid.
2. **The score must use the forearm angle, not the device joint angle.** Under the P+R interface the device angle moves k times the forearm angle. At the Rev C offset limit of 46.3 mm (k {{pr.e_cap.k_min|.3f}}–{{pr.e_cap.k_max|.3f}}), a device reading of 10° is a forearm drop of {{pr.e_cap.hum_at_dev10_min|.2f}}–{{pr.e_cap.hum_at_dev10_max|.2f}}°; at 100 mm it is {{pr.e_src.hum_at_dev10_min|.2f}}–{{pr.e_src.hum_at_dev10_max|.2f}}°. Joint minus pivot angle removes this: its peak error is bounded by two counts, ±{{r2.budget.two_count_bound_deg|.4f}}° (observed {{r2.jp_peak_err.min|+.4f}}° to {{r2.jp_peak_err.max|+.4f}}° over 4000 trials). But scored through PR #9's chain it fails PR #9's thresholds (headline above).
3. **The time rule scores arms that never pass 10° as DROP.** "Catch time" is arrest complete (speed zero), but the 201.90 ms edge is the 10° crossing of an arm stopped instantly. With a finite arresting torque (ASSUMPTION, 15–50 N·m), a trial whose arrest completes just inside the edge has a true peak of only {{catch_time.15.peak_at_time_edge_deg|.2f}}° to {{catch_time.50.peak_at_time_edge_deg|.2f}}°. A true peak of exactly 10° completes its arrest at {{catch_time.50.t_c_true_ms|.1f}}–{{catch_time.15.t_c_true_ms|.1f}} ms, after the edge. So the time rule tests a stricter quantity than B9's 10°; PR #9 A9 (how time and angle combine) needs this.
4. **Catch time cannot be measured to U_t from the angle.** Near speed zero the angle is flat. With 17 bits and no noise, the argmax time is early by up to {{catch_time.worst_arg_ms|.2f}} ms; a 15-sample speed-zero rule is off by up to {{catch_time.worst_v15_ms|.2f}} ms; with 1 count of noise the errors reach {{catch_time.worst_noise_ms|.0f}} ms. Meeting U_t would need {{catch_time.50.bits_for_u_t}}–{{catch_time.15.bits_for_u_t}} bits and no noise.
5. **Gain and zero errors are not in the {{inputs.guard|.4f}}° budget.** A gain error of {{miss_limits.gain|.3%}} alone, or a zero-reference error of {{miss_limits.zero|.4f}}° alone, makes a true 10.002° drop read NO DROP. A −0.5 % gain error missed {{faults.gain_minus_0p5pct.missed|.0%}} of true 10.01° drops, and no plausibility check flags it.
6. **Noise never causes a miss here but widens the false-trip band.** The peak is a maximum, so noise biases it up. With 1 count of noise every arm at {{peak_oc.noise_1lsb.trip_all_delta_abs_deg|.2f}}° below the limit or closer is scored DROP; with 4 counts, every arm within {{peak_oc.noise_4lsb.trip_all_delta_abs_deg|.2f}}°. A 3-sample median brings the 4-count bias from +{{peak_oc.noise_4lsb.err_min|.4f}}…+{{peak_oc.noise_4lsb.err_max|.4f}}° to +{{peak_oc.noise_4lsb_med3.err_min|.4f}}…+{{peak_oc.noise_4lsb_med3.err_max|.4f}}°.
7. **A stuck reading is the main missed-detection fault, and a frozen-count check catches it when the reading has noise.** These rates are *conditional*: they assume the freeze happens (no fault rate is known), at a uniform time in [0, t_c], on a true 10.01° peak with T_h 30 N·m (all ASSUMPTION). Under that model one channel frozen during the drop misses {{r2.stuck_sens.base.missed.p|.1%}} (95 % CI {{r2.stuck_sens.base.missed.ci.0|.1%}}–{{r2.stuck_sens.base.missed.ci.1|.1%}}, n {{r2.stuck_sens.base.missed.n}}). Motion plus deceleration (50,000 °/s²) checks leave {{r2.stuck_sens.base.after_a50k.p|.1%}} ({{r2.stuck_sens.base.after_a50k.ci.0|.1%}}–{{r2.stuck_sens.base.after_a50k.ci.1|.1%}}); over the sensitivity cases {{r2.stuck_sens.range_after_a50k.0|.1%}}–{{r2.stuck_sens.range_after_a50k.1|.1%}}. With the deceleration limit just above the clean maximum ({{r2.clean_1ch.a_lim_tight|.0f}} °/s²) they leave {{r2.stuck_sens.range_after_a_tight.0|.1%}}–{{r2.stuck_sens.range_after_a_tight.1|.1%}}. Adding a frozen-count check (no run of identical raw counts longer than {{r2.assumptions.run_lim_ms|.0f}} ms) leaves {{r2.stuck_sens.base_noise1.after_tight_run.p|.1%}} (CI upper {{r2.stuck_sens.base_noise1.after_tight_run.ci.1|.2%}}) with 1 count of noise. That check needs at least about 0.5 count of noise; with none, a held clean arm is "frozen" too (section 10.3). Round 1's "97.5 % / 26.2 %" was n 400 (±{{r2.stuck_sens.ci_halfwidth_n400_pts|.1f}} points) under the same conditions. A second reading is an option, not shown to be needed (section 10.2).
8. **The Rev C joint-encoder requirement (0.1°) and the repo IMU default (150 Hz) both fail the budget.** At 0.1°: crossing time off by {{crossing.revc_joint.t10_abs_max_ms|.3f}} ms, drop at 200 ms off by {{crossing.revc_joint.d200_abs_max_deg|.3f}}°. At 150 Hz: crossing time off by {{crossing.imu150.t10_abs_max_ms|.2f}} ms (margin {{model.margin_ms|.3f}} ms), drop at 200 ms off by {{crossing.imu150.d200_abs_max_deg|.3f}}°. Defaults are quoted, not changed.
9. **Plan level.** To keep the chance of at least one false DROP at 5 % or less, the per-trial false-trip rate must be at most {{plan.p_for_5pct.n59|.2e}} over 59 trials and {{plan.p_for_5pct.n330|.2e}} over 330. A rate of 1e-3 gives {{plan.p1e-03.n59|.1%}} and {{plan.p1e-03.n330|.1%}}.

## 1. Drop model and checks

- Re-implementation of HE's drop model as used in PR #9 (`scripts/ba15_numcheck.py`; same inputs), governing cell (b) e-stop, ramp 21.00 N·m/s, STO 148 ms, hold 0. Torque basis: T0 and ramp device side, × k = 1.5862 at the arm; G, hold and arrest torque arm side.
- 10° crossing {{model.t10_ms|.6f}} ms (PR #9: {{model.ref.t10_ms|.6f}}; difference {{model.d_t10_ms|.1e}} ms). Drop at 200 ms {{model.drop200_deg|.6f}}° (PR #9 {{model.ref.drop200_deg|.6f}}). Speed at 10° {{model.w10_dps|.3f}} °/s (PR #9 {{model.ref.w10_dps|.3f}}). Drop at 100 ms {{model.drop100_deg|.3f}}°.
- **Catch model (ASSUMPTION):** from a catch-onset time the user adds a constant arresting torque T_h (15, 20, 30, 50 N·m at the arm) until speed is zero; then the arm is held at the peak. The onset is solved so the true peak is 10° + δ. No human data is behind T_h.

| T_h (N·m) | Deceleration near 10° (°/s²) | Catch onset for a 10° peak (ms) | Arrest complete (ms) |
|---|---|---|---|
| 15 | {{catch_time.15.decel_dps2|.0f}} | {{catch_time.15.t_on_for_10_ms|.2f}} | {{catch_time.15.t_c_true_ms|.2f}} |
| 20 | {{catch_time.20.decel_dps2|.0f}} | {{catch_time.20.t_on_for_10_ms|.2f}} | {{catch_time.20.t_c_true_ms|.2f}} |
| 30 | {{catch_time.30.decel_dps2|.0f}} | {{catch_time.30.t_on_for_10_ms|.2f}} | {{catch_time.30.t_c_true_ms|.2f}} |
| 50 | {{catch_time.50.decel_dps2|.0f}} | {{catch_time.50.t_on_for_10_ms|.2f}} | {{catch_time.50.t_c_true_ms|.2f}} |

## 2. Measurement chain

Per trial: random sampling phase, encoder mounting phase (where the count edges fall) and skew; latency; floor quantization; optional white noise, gain, zero error, 3-sample median and faults. Baseline: mean of the samples 50 ms before onset. Peak drop: maximum after onset. Crossing: first sample at or above 10°, no interpolation (as in PR #9's budget). Drop at 200 ms: last sample logged at or before 200 ms. Score: DROP if peak + {{inputs.guard|.4f}}° > 10°. {{inputs.n_mc}} trials per point unless stated; seeds from one `SeedSequence({{inputs.seed}})` by spawn.

Requirement chain ("req"): {{inputs.fs_req|.0f}} Hz, {{inputs.bits_req}} bits ({{inputs.q_req|.5f}}°), latency {{inputs.enc_lat_ms|.4f}} ms, skew uniform within ±{{inputs.skew_total_ms|.4f}} ms (trigger {{inputs.trig_skew_ms|.4f}} + release {{inputs.rel_skew_ms|.4f}}; ASSUMPTION: release skew acts as a timing offset). All quoted from PR #9.

## 3. Crossing channel (Phase A dummy matrix, per-trial pass)

Error = measured − true. A positive crossing error or a negative drop error helps a unit pass.

| Chain | 10° crossing error (ms), min / mean / max | Drop at 200 ms error (°), min / max | Within U_t / guard |
|---|---|---|---|
| req, skew random | {{crossing.req.t10_err_ms.min|.4f}} / {{crossing.req.t10_err_ms.mean|.4f}} / {{crossing.req.t10_err_ms.max|.4f}} | {{crossing.req.d200_err_deg.min|.4f}} / {{crossing.req.d200_err_deg.max|.4f}} | {{crossing.req.within_u_t}} / {{crossing.req.within_guard}} |
| req, skew at −edge | {{crossing.req_skew_neg.t10_err_ms.min|.4f}} / {{crossing.req_skew_neg.t10_err_ms.mean|.4f}} / {{crossing.req_skew_neg.t10_err_ms.max|.4f}} | {{crossing.req_skew_neg.d200_err_deg.min|.4f}} / {{crossing.req_skew_neg.d200_err_deg.max|.4f}} | {{crossing.req_skew_neg.within_u_t}} / {{crossing.req_skew_neg.within_guard}} |
| req, skew at +edge | {{crossing.req_skew_pos.t10_err_ms.min|.4f}} / {{crossing.req_skew_pos.t10_err_ms.mean|.4f}} / {{crossing.req_skew_pos.t10_err_ms.max|.4f}} | {{crossing.req_skew_pos.d200_err_deg.min|.4f}} / {{crossing.req_skew_pos.d200_err_deg.max|.4f}} | {{crossing.req_skew_pos.within_u_t}} / {{crossing.req_skew_pos.within_guard}} |
| 0.1° steps (Rev C joint encoder) | {{crossing.revc_joint.t10_err_ms.min|.4f}} / {{crossing.revc_joint.t10_err_ms.mean|.4f}} / {{crossing.revc_joint.t10_err_ms.max|.4f}} | {{crossing.revc_joint.d200_err_deg.min|.4f}} / {{crossing.revc_joint.d200_err_deg.max|.4f}} | {{crossing.revc_joint.within_u_t}} / {{crossing.revc_joint.within_guard}} |
| 150 Hz (repo IMU default) | {{crossing.imu150.t10_err_ms.min|.4f}} / {{crossing.imu150.t10_err_ms.mean|.4f}} / {{crossing.imu150.t10_err_ms.max|.4f}} | {{crossing.imu150.d200_err_deg.min|.4f}} / {{crossing.imu150.d200_err_deg.max|.4f}} | {{crossing.imu150.within_u_t}} / {{crossing.imu150.within_guard}} |

- Room left at the requirement chain, 400-trial maximum (round 1): {{crossing.u_t_room_ms|.4f}} ms of U_t and {{crossing.guard_room_deg|.4f}}° of the guard. **Exact worst case (round 2):** {{r2.exact.t10_abs_max_ms|.5f}} ms and {{r2.exact.d200_abs_max_deg|.5f}}°, room {{r2.exact.u_t_room_ms|.4f}} ms and {{r2.exact.guard_room_deg|.4f}}°. Both rooms rely on where 10° and the 200 ms drop fall on the count grid; the lattice-free crossing error is {{r2.exact.lattice_free_t10_ms|.4f}} ms (room {{r2.exact.lattice_free_room_ms|.4f}} ms), which is PR #9's linear budget by design.
- At 150 Hz, linear interpolation between samples brings the crossing error to {{crossing.imu150.t10i_err_ms.min|.3f}} to {{crossing.imu150.t10i_err_ms.max|.3f}} ms on this smooth model trace; that is not a bound for real data.

## 4. Peak-drop channel (Phase B angle score)

P(DROP) at true peak 10° + δ, T_h 30 N·m. δ > 0 is a real drop past 10° (a NO DROP there is a miss).

| δ (°) | req | 1 count noise | 4 counts noise | 4 counts + median-3 | 0.1° steps | 150 Hz |
|---|---|---|---|---|---|---|
{{oc_table}}

| Chain | Peak error (°), min / max | Largest miss rate (δ > 0) | Some DROP from δ = | All DROP from δ = |
|---|---|---|---|---|
{{oc_summary}}

- The guard makes every arm within {{inputs.guard|.4f}}° below 10° a possible DROP by design (PR #9: the score is meant to be strict).
- Peak error at the requirement chain over all four arrest torques: {{peak_err_by_Th.worst_min|.4f}}° to +{{peak_err_by_Th.worst_max|.4f}}°. Sampling rate barely matters for the peak, because the speed is zero there and the arm is held.

### Single-source limits for a missed detection (true peak 10.002°, worst of 200 draws)

| Source | Limit (this chain) | Simple bound (guard + 0.002°) |
|---|---|---|
| Gain (scale) error, reading low | {{miss_limits.gain|.4%}} | {{miss_limits.analytic.gain|.4%}} |
| Zero-reference error, drop read low | {{miss_limits.zero|.5f}}° | {{miss_limits.analytic.zero|.5f}}° |

Resolution scan (true peaks 10.002–10.05°, worst miss rate over 13 peaks, 100 trials each):

| Bits | Step (°) | Worst miss, no noise | Worst miss, 0.5-count noise |
|---|---|---|---|
{{q_table}}

- Without noise the measured drop is a whole number of counts, so whether a miss shows up depends on where 10° falls on the count grid. Misses were seen up to {{miss_limits.q_bits_last_miss}} bits. A step below the guard needs at least {{miss_limits.q_bits_guard}} bits; the 17-bit requirement has room.

## 5. Catch-time channel (time rule t_c + U_t > 201.90 ms)

| T_h (N·m) | True peak when arrest completes at the edge (°) | Argmax time error, no noise (ms) | Speed-zero W = 15 error (ms) | W = 75 error (ms) | W = 1 error (ms) | Argmax, 1 count noise (ms) | W = 15, 1 count noise (ms) | Flat time for 1 count (ms) | Bits for U_t |
|---|---|---|---|---|---|---|---|---|---|
{{tc_table}}

- Errors are min to max over 200 trials. "Flat time for 1 count" = √(2q/a): how long the angle stays within one count around the turning point. "Bits for U_t" = bits whose step is ½·a·U_t².
- A one-sample difference (W = 1) reads zero speed as soon as one sample has no count change, so it marks arrest at the start of the motion.
- With noise the argmax picks a noise maximum during the hold after arrest, up to {{catch_time.worst_noise_arg_ms|.0f}} ms late.

## 6. Which angle: device joint vs forearm (P+R interface, Rev C §3.7)

- **Model scope (ASSUMPTION):** planar, device and elbow axes parallel, rigid cuffs and links, offset fixed during the drop. Out-of-plane axis tilt, elbow-axis wander and cuff or soft-tissue motion are not modelled. Rev C itself marks its 100 mm offset source as not verified for the elbow.
- Geometry quoted from Rev C: pivot height 54 mm, s_mid 175 mm, link 0.183 m; d = √(L² − (h_B − e_y)²), s_B = e_x + d, k = s_B/d. In the script the device link angle is ψ = φ + α with sin α = (h_B − e_y′)/L, which gives dψ/dφ = 1 + e_x′/d = s_B/d = k: the device angle moves k times the forearm angle. Numerical check at one offset: k {{pr.k_check.0|.9f}}, dψ/dφ {{pr.k_check.1|.9f}}.
- Disc sweep, 41 radii × 721 directions (as Rev C). The 100 mm k_max reproduces Rev C's k_hi 1.5862 ({{pr.e_src.k_max|.10f}}).

| Offset bound | k range | Device change for a 10° forearm drop (°) | Forearm drop when the device reads 10° (°) |
|---|---|---|---|
| 46.3 mm (Rev C pass bound) | {{pr.e_cap.k_min|.4f}}–{{pr.e_cap.k_max|.4f}} | {{pr.e_cap.dev_for_10_min|.3f}}–{{pr.e_cap.dev_for_10_max|.3f}} | {{pr.e_cap.hum_at_dev10_min|.3f}}–{{pr.e_cap.hum_at_dev10_max|.3f}} |
| 100 mm, all offsets | {{pr.e_src.k_min|.4f}}–{{pr.e_src.k_max|.4f}} | {{pr.e_src.dev_for_10_min|.3f}}–{{pr.e_src.dev_for_10_max|.3f}} | {{pr.e_src.hum_at_dev10_min|.3f}}–{{pr.e_src.hum_at_dev10_max|.3f}} |
| 100 mm, carriage inside the 100–250 mm rail | {{pr.e_src_clip.k_min|.4f}}–{{pr.e_src_clip.k_max|.4f}} | {{pr.e_src_clip.dev_for_10_min|.3f}}–{{pr.e_src_clip.dev_for_10_max|.3f}} | {{pr.e_src_clip.hum_at_dev10_min|.3f}}–{{pr.e_src_clip.hum_at_dev10_max|.3f}} |

- Through the chain at a true 10.000° forearm peak (46.3 mm extremes, linear k): joint angle only, peak error {{pr.joint_only_k_lo_peak_err.mean|.3f}}° at k_lo and +{{pr.joint_only_k_hi_peak_err.mean|.3f}}° at k_hi. Joint minus pivot angle (two 17-bit channels): bound ±2 counts = ±{{r2.budget.two_count_bound_deg|.4f}}°; observed {{r2.jp_peak_err.min|+.4f}}° to {{r2.jp_peak_err.max|+.4f}}° over 4000 trials (round 1 quoted a 200-trial range, {{pr.joint_pivot_k_lo_peak_err.min|.4f}}° to +{{pr.joint_pivot_k_lo_peak_err.max|.4f}}°, as if it were the bound).
- On the dummy (Phase A) the offset is set by the build; in people (Phase B) it is not known, which is why Rev C added the carriage sensor and pivot angle encoder for the torque cap.

## 7. Encoder faults and plausibility checks

Round 1 (this section) injects faults on one channel only. Faults on both channels of the joint − pivot reading, common-mode faults, a frozen-count check and the stuck-fault sensitivity are in section 10.2–10.4.

True peak 10.01° (a real drop), T_h 30 N·m, {{inputs.n_mc}} trials. Checks (ASSUMPTION): **motion**, measured drop at 100 ms at least half the model's {{model.drop100_deg|.3f}}° (needs a relaxed subject who cannot catch before 100 ms); **deceleration**, no measured deceleration above {{inputs.a_lim|.0f}} °/s² (1 ms speed windows; the largest modelled is {{faults.model_max_decel|.0f}} °/s²); **sample counter**, any lost sample makes the trial INVALID (PR #9 logs the sample-drop count).

| Case | Missed (no checks) | Flagged INVALID | Missed after checks |
|---|---|---|---|
| Clean | {{faults.clean.missed|.1%}} | {{faults.clean.flagged|.1%}} | {{faults.clean.missed_after_checks|.1%}} |
| Clean, 1 count noise | {{faults.clean_noise1.missed|.1%}} | {{faults.clean_noise1.flagged|.1%}} | {{faults.clean_noise1.missed_after_checks|.1%}} |
| Stuck from onset | {{faults.stuck_from_onset.missed|.1%}} | {{faults.stuck_from_onset.flagged|.1%}} | {{faults.stuck_from_onset.missed_after_checks|.1%}} |
| Stuck at a random time during the drop | {{faults.stuck_during_drop.missed|.1%}} | {{faults.stuck_during_drop.flagged|.1%}} | {{faults.stuck_during_drop.missed_after_checks|.1%}} |
| 2 ms of samples lost at the peak, arm held | {{faults.gap_2ms_at_peak_held.missed|.1%}} | {{faults.gap_2ms_at_peak_held.flagged|.1%}} | {{faults.gap_2ms_at_peak_held.missed_after_checks|.1%}} |
| 2 ms lost at the peak, arm rebounds (ASSUMPTION) | {{faults.gap_2ms_at_peak_rebound.missed|.1%}} | {{faults.gap_2ms_at_peak_rebound.flagged|.1%}} | {{faults.gap_2ms_at_peak_rebound.missed_after_checks|.1%}} |
| Gain −0.5 % | {{faults.gain_minus_0p5pct.missed|.1%}} | {{faults.gain_minus_0p5pct.flagged|.1%}} | {{faults.gain_minus_0p5pct.missed_after_checks|.1%}} |

- Clean trials at all four arrest torques with 1 count of noise: flagged {{faults.clean_flag_rate_all_Th|.1%}}; largest measured deceleration {{faults.clean_max_decel_all_Th|.0f}} °/s² (noise dominates it).
- Deceleration limit trade-off (stuck during the drop vs clean noisy trials):

| Limit (°/s²) | Stuck faults missed after checks | Clean trials flagged |
|---|---|---|
{{decel_table}}

- A gap at a rebounding peak under-reads by about ½·a·(D/2)²: {{faults.gap_rebound_underread_2ms_analytic|.5f}}° for 2 ms; a gap up to {{faults.gap_rebound_max_ms_for_guard|.2f}} ms stays inside the guard.
- One-sample glitch of 0.1° on an arm truly at 9.95°: false trip {{faults.glitch_0p1.false_trip|.1%}} of trials (when it lands near the peak); with a 3-sample median {{faults.glitch_0p1_med3.false_trip|.1%}}.

## 8. Requirements that follow (by requirement, no parts; SYNTHETIC basis)

1. Score the **forearm** angle (device joint angle combined with the pivot angle), never the device joint angle alone. Log both raw channels. **Round 2: with PR #9's chain this FAILS PR #9's thresholds** (drop-at-200 {{r2.budget.two_ch_worst_d200_deg|.5f}}° vs {{inputs.guard|.5f}}°; crossing {{r2.budget.two_ch_worst_t10_ms|.5f}} ms vs {{inputs.u_t_ms|.5f}} ms). Something in section 10.1 has to change; Charlie decides which.
2. Keep PR #9's chain (≥ 15 kHz, ≤ 0.0031° step, ≤ 0.019 ms latency, skews) and add to the angle budget: **gain/linearity error over the 0–15° span** and **zero-reference error**. Alone, each uses the whole guard at {{miss_limits.gain|.3%}} and {{miss_limits.zero|.4f}}°. Calibrate the angle scale on the dummy against a reference.
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
| A5 | Release skew acts as a timing offset, like trigger skew; the two add to ±{{inputs.skew_total_ms|.4f}} ms. |
| A6 | P+R relation from Rev C §3.7 geometry; the chain test uses a constant k over the 10° drop (the disc sweep uses the exact geometry). |
| A7 | Plausibility thresholds (motion at 100 ms, deceleration {{inputs.a_lim|.0f}} °/s²) are illustrations, not settings. |
| A8 | Plan-level rates assume independent trials. |
| A9 | Record length {{r2.assumptions.t_end_s|.2f}} s after onset (the noise bias of the maximum grows with the hold length) and a glitch window of {{r2.assumptions.glitch_window_s|.1f}} s after onset. |
| A10 | Arresting torque 15–50 N·m: **ASSUMPTION, no source.** Section 5, finding 3 and the stuck sensitivity scale with it. |
| A11 | Round 2 faults: one at a time; sizes (5 ms holds at 20 per s, 0.05° / 0.2° slip, 0.05° / 0.2° backlash, 0.2 °/s drift, −0.5 % common gain) are illustrations. Slip and drift take the sign that makes the drop read low. Backlash starts at a random place in the play. |
| A12 | Frozen-count limit {{r2.assumptions.run_lim_ms|.0f}} ms; cross-check tolerance {{r2.assumptions.xchk_tol_deg|.2f}}°; the second reading is a 17-bit forearm angle with its own noise and mount, sharing only the injected common fault. |
| A13 | Two-channel budget check: free-fall governing cell, skew fixed at the edge that makes times read late, {{r2.assumptions.n_budget}} random phase and mount draws per row. An MC maximum is a lower bound on the worst case, so a FAIL is a found counterexample, but a pass is not proof. |

## 10. Fix round 2 (FC round 1 on PR #14 @ 1df8a94)

### 10.1 Forearm angle (joint − pivot) against PR #9's thresholds: FAILS

Two independent 17-bit channels, PR #9's chain otherwise unchanged (15 kHz, latency 0.019 ms, skew at the edge that reads late), free-fall governing cell, {{r2.assumptions.n_budget}} trials per row:

| Chain | Worst crossing error (ms) vs U_t {{inputs.u_t_ms|.5f}} | Result | Worst drop-at-200 error (°) vs −{{inputs.guard|.5f}} | Result |
|---|---|---|---|---|
{{r2_budget_table}}

- The overrun is small (drop-at-200 by {{r2.budget.two_ch_d200_over_deg|.1e}}°, crossing by {{r2.budget.two_ch_t10_over_ms|.1e}} ms). But the thresholds are hard per-trial pass lines, and the second channel's step is not in PR #9's budget. **This is a failing result.** PR #9's thresholds are not loosened here.
- Options that pass in this check (same MC, k 0.7345 / 1.0001 / 1.2655; worst over the three). **None is chosen**; each is for HE / Charlie:

| Option | Worst crossing error (ms) | Worst drop-at-200 error (°) | Result |
|---|---|---|---|
{{r2_options_table}}

- Loosening option (Charlie's call only, not proposed): add one pivot count to each budget, giving U_t {{r2.budget.loosen.u_t_ms|.4f}} ms and guard {{r2.budget.loosen.guard_deg|.4f}}°.
- "Latency 0" is a bound case showing that latency is one of the terms that matters, not a buildable value. Interpolation changes the scoring rule, which is PR #9's. A pass here is an MC result (A13).

### 10.2 Faults on both channels, common-mode faults, and checks

Forearm = joint − pivot at k 0.7345 (the 46.3 mm low end, where a frozen pivot under-reads), true peak 10.01°, T_h 30 N·m, 1 count of noise per channel, {{r2.assumptions.n_r2}} trials per row. Checks are cumulative left to right: base = motion + deceleration limit {{r2.clean.a_lim_tight|.0f}} °/s² (just above the clean two-channel maximum {{r2.clean.max_dec|.0f}}); + frozen count = no run of identical raw counts on either channel longer than {{r2.assumptions.run_lim_ms|.0f}} ms; + cross-check = second reading within {{r2.assumptions.xchk_tol_deg|.2f}}°.

| Fault | Missed, no checks | Missed after base | After + frozen count | After + cross-check | Flagged INVALID (all checks) |
|---|---|---|---|---|---|
{{r2_fault_table}}

¹ With no noise the frozen-count check also flags every clean trial (the held arm reads constant; 10.3), so the "+ frozen count" columns are not usable for these two rows.

- Clean trials, all four arrest torques × k 0.7345 / 1.0001 / 1.2655, {{r2.clean.n}} trials: flagged {{r2.clean.flagged_base|.1%}} (base), {{r2.clean.flagged_run|.1%}} (+ frozen count), {{r2.clean.flagged_xchk|.1%}} (+ cross-check). Longest clean run {{r2.clean.max_run_ms|.2f}} ms.
- Pivot frozen from onset at k 1.2655 over-reads: a true 9.95° arm is scored DROP {{r2.faults.pivot_stuck_onset_khi_false_drop.false_drop.p|.1%}} of the time (a false trip); the frozen-count check flags {{r2.faults.pivot_stuck_onset_khi_false_drop.flagged_run|.1%}}.
- What the rows show: a frozen pivot passes the round-1 checks every time (as FC found). The frozen-count check catches stuck readings on either channel, including a common acquisition freeze, but only when the reading has noise. Slips of 0.05° or more trip the deceleration check. Partial sticking is flagged but causes no miss here, because the arm is held at its peak. Backlash, drift and common gain are slow, and no check except the cross-check sees them. The cross-check cuts backlash 0.2° misses from {{r2.faults.joint_backlash_0p2.missed_run.p|.1%}} to {{r2.faults.joint_backlash_0p2.missed_xchk.p|.1%}} and drift from {{r2.faults.joint_drift_0p2dps.missed_run.p|.1%}} to {{r2.faults.joint_drift_0p2dps.missed_xchk.p|.1%}}. It does not help with a −0.5 % common gain ({{r2.faults.common_gain_m0p5pct.missed_xchk.p|.1%}}), because the under-read is about its {{r2.assumptions.xchk_tol_deg|.2f}}° tolerance.
- **Second reading, open:** its own accuracy on a person (cuff and soft tissue), its tolerance budget, what a disagreement scores (INVALID or DROP), its common-cause links (shaft or coupling, supply, clock, acquisition), and whether it fits on the rig are not modelled. No part is named (Charlie's constraint).
- **Not modelled:** encoder nonlinearity beyond a scale error; electromagnetic bursts; counter roll-over or wrong sign; clock drift between channels; partial supply loss; joint-encoder placement on the motor side of the gearbox (Rev C BOM does not say); out-of-plane misalignment.

### 10.3 Frozen-count check against noise (one channel, k 1, T_h 30, 1000 clean and 1000 stuck trials per row)

| Noise (counts RMS) | Clean longest run (ms) | Stuck shortest run (ms) | Clean flagged at {{r2.assumptions.run_lim_ms|.0f}} ms | Fully separated |
|---|---|---|---|---|
{{r2_sep_table}}

### 10.4 Stuck reading: assumptions and sensitivity (one channel, as round 1; n {{r2.assumptions.n_sens}} per row, 95 % Clopper-Pearson)

Conditions: the freeze always happens (no fault rate), at a uniform time over the window shown, on one true peak and one arrest torque per row. These are ASSUMPTIONS; the rates are conditional on a freeze.

| Case | Missed, no checks | After motion + decel 50,000 | After motion + decel {{r2.clean_1ch.a_lim_tight|.0f}} | Share of the 50,000-limit misses reading < 5° |
|---|---|---|---|---|
{{r2_sens_table}}

- {{r2.clean_1ch.a_lim_tight|.0f}} °/s² is just above the clean one-channel maximum {{r2.clean_1ch.max_dec|.0f}} °/s² (1 count noise, four arrest torques, {{r2.clean_1ch.n}} trials). With the frozen-count check added and 1 count of noise: {{r2.stuck_sens.base_noise1.after_tight_run.p|.1%}} missed (CI {{r2.stuck_sens.base_noise1.after_tight_run.ci.0|.2%}}–{{r2.stuck_sens.base_noise1.after_tight_run.ci.1|.2%}}), and clean trials flagged {{r2.clean_1ch.flagged_tight_run|.1%}}.
- Most after-check misses are gross under-reads (reading < 5°), not near-10° cases.

## 11. Reproduce

    python scripts/angle_detector_risk.py --check     # about 20 min on one core since round 2; exit 0 if both files reproduce
    pytest tests/test_angle_detector_risk.py          # fast checks (model vs PR #9, geometry, chain, rendering)
