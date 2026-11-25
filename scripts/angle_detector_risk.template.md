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

1. **The PR #9 requirement set meets its own budget, with thin room.** At 15 kHz, 17 bits, 0.019 ms latency and ±{{inputs.skew_total_ms|.3f}} ms skew, the 10° crossing time is off by at most {{crossing.worst_t10_abs_max_ms|.3f}} ms (budget U_t {{inputs.u_t_ms|.3f}} ms) and the drop at 200 ms by at most {{crossing.worst_d200_abs_max_deg|.4f}}° (budget {{inputs.guard|.4f}}°). The peak drop is off by {{peak_oc.req.err_min|.4f}}° to +{{peak_oc.req.err_max|.4f}}°. No missed detection was seen at any true drop above 10° in the grid.
2. **The score must use the forearm angle, not the device joint angle.** Under the P+R interface the device angle moves k times the forearm angle. At the Rev C offset limit of 46.3 mm (k {{pr.e_cap.k_min|.3f}}–{{pr.e_cap.k_max|.3f}}), a device reading of 10° is a forearm drop of {{pr.e_cap.hum_at_dev10_min|.2f}}–{{pr.e_cap.hum_at_dev10_max|.2f}}°; at 100 mm it is {{pr.e_src.hum_at_dev10_min|.2f}}–{{pr.e_src.hum_at_dev10_max|.2f}}°. Joint plus pivot angle removes this (peak error within ±{{pr.joint_pivot_abs_max|.4f}}°).
3. **The time rule scores arms that never pass 10° as DROP.** "Catch time" is arrest complete (speed zero), but the 201.90 ms edge is the 10° crossing of an arm stopped instantly. With a finite arresting torque (ASSUMPTION, 15–50 N·m), a trial whose arrest completes just inside the edge has a true peak of only {{catch_time.15.peak_at_time_edge_deg|.2f}}° to {{catch_time.50.peak_at_time_edge_deg|.2f}}°. A true peak of exactly 10° completes its arrest at {{catch_time.50.t_c_true_ms|.1f}}–{{catch_time.15.t_c_true_ms|.1f}} ms, after the edge. So the time rule tests a stricter quantity than B9's 10°; PR #9 A9 (how time and angle combine) needs this.
4. **Catch time cannot be measured to U_t from the angle.** Near speed zero the angle is flat. With 17 bits and no noise, the argmax time is early by up to {{catch_time.worst_arg_ms|.2f}} ms; a 15-sample speed-zero rule is off by up to {{catch_time.worst_v15_ms|.2f}} ms; with 1 count of noise the errors reach {{catch_time.worst_noise_ms|.0f}} ms. Meeting U_t would need {{catch_time.50.bits_for_u_t}}–{{catch_time.15.bits_for_u_t}} bits and no noise.
5. **Gain and zero errors are not in the {{inputs.guard|.4f}}° budget.** A gain error of {{miss_limits.gain|.3%}} alone, or a zero-reference error of {{miss_limits.zero|.4f}}° alone, makes a true 10.002° drop read NO DROP. A −0.5 % gain error missed {{faults.gain_minus_0p5pct.missed|.0%}} of true 10.01° drops, and no plausibility check flags it.
6. **Noise never causes a miss here but widens the false-trip band.** The peak is a maximum, so noise biases it up. With 1 count of noise every arm at {{peak_oc.noise_1lsb.trip_all_delta_abs_deg|.2f}}° below the limit or closer is scored DROP; with 4 counts, every arm within {{peak_oc.noise_4lsb.trip_all_delta_abs_deg|.2f}}°. A 3-sample median brings the 4-count bias from +{{peak_oc.noise_4lsb.err_min|.4f}}…+{{peak_oc.noise_4lsb.err_max|.4f}}° to +{{peak_oc.noise_4lsb_med3.err_min|.4f}}…+{{peak_oc.noise_4lsb_med3.err_max|.4f}}°.
7. **A stuck encoder is the main missed-detection fault.** Frozen from onset, it is caught by a motion check {{faults.stuck_from_onset.flagged|.0%}} of the time. Frozen at a random time during the drop, it misses {{faults.stuck_during_drop.missed|.1%}} of true 10.01° drops; motion plus deceleration checks still leave {{faults.stuck_during_drop.missed_after_checks|.1%}} missed. Tightening the deceleration limit flags clean noisy trials first. A second, independent angle reading is needed.
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

- Room left at the requirement chain: {{crossing.u_t_room_ms|.4f}} ms of U_t and {{crossing.guard_room_deg|.4f}}° of the angle guard (worst rows).
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

- Geometry quoted from Rev C: pivot height 54 mm, s_mid 175 mm, link 0.183 m; d = √(L² − (h_B − e_y)²), s_B = e_x + d, k = s_B/d. In the script the device link angle is ψ = φ + α with sin α = (h_B − e_y′)/L, which gives dψ/dφ = 1 + e_x′/d = s_B/d = k: the device angle moves k times the forearm angle. Numerical check at one offset: k {{pr.k_check.0|.9f}}, dψ/dφ {{pr.k_check.1|.9f}}.
- Disc sweep, 41 radii × 721 directions (as Rev C). The 100 mm k_max reproduces Rev C's k_hi 1.5862 ({{pr.e_src.k_max|.10f}}).

| Offset bound | k range | Device change for a 10° forearm drop (°) | Forearm drop when the device reads 10° (°) |
|---|---|---|---|
| 46.3 mm (Rev C pass bound) | {{pr.e_cap.k_min|.4f}}–{{pr.e_cap.k_max|.4f}} | {{pr.e_cap.dev_for_10_min|.3f}}–{{pr.e_cap.dev_for_10_max|.3f}} | {{pr.e_cap.hum_at_dev10_min|.3f}}–{{pr.e_cap.hum_at_dev10_max|.3f}} |
| 100 mm, all offsets | {{pr.e_src.k_min|.4f}}–{{pr.e_src.k_max|.4f}} | {{pr.e_src.dev_for_10_min|.3f}}–{{pr.e_src.dev_for_10_max|.3f}} | {{pr.e_src.hum_at_dev10_min|.3f}}–{{pr.e_src.hum_at_dev10_max|.3f}} |
| 100 mm, carriage inside the 100–250 mm rail | {{pr.e_src_clip.k_min|.4f}}–{{pr.e_src_clip.k_max|.4f}} | {{pr.e_src_clip.dev_for_10_min|.3f}}–{{pr.e_src_clip.dev_for_10_max|.3f}} | {{pr.e_src_clip.hum_at_dev10_min|.3f}}–{{pr.e_src_clip.hum_at_dev10_max|.3f}} |

- Through the chain at a true 10.000° forearm peak (46.3 mm extremes, linear k): joint angle only, peak error {{pr.joint_only_k_lo_peak_err.mean|.3f}}° at k_lo and +{{pr.joint_only_k_hi_peak_err.mean|.3f}}° at k_hi. Joint minus pivot angle (two 17-bit channels): {{pr.joint_pivot_k_lo_peak_err.min|.4f}}° to +{{pr.joint_pivot_k_lo_peak_err.max|.4f}}° at k_lo.
- On the dummy (Phase A) the offset is set by the build; in people (Phase B) it is not known, which is why Rev C added the carriage sensor and pivot angle encoder for the torque cap.

## 7. Encoder faults and plausibility checks

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

1. Score the **forearm** angle (device joint angle combined with the pivot angle), never the device joint angle alone. Log both raw channels.
2. Keep PR #9's chain (≥ 15 kHz, ≤ 0.0031° step, ≤ 0.019 ms latency, skews) and add to the angle budget: **gain/linearity error over the 0–15° span** and **zero-reference error**. Alone, each uses the whole guard at {{miss_limits.gain|.3%}} and {{miss_limits.zero|.4f}}°. Calibrate the angle scale on the dummy against a reference.
3. **Angle noise ≤ 1 count RMS** at the scoring rate, or a stated filter (3-sample median at least) whose upward bias is added to the false-trip analysis.
4. **Do not score the catch time from speed zero of the angle.** Either drop the time rule (the angle rule is the B9 quantity), or redefine it as a quantity the angle can measure at speed (for example the 10° crossing time) and recompute its window. FC and Charlie decide (PR #9 A9 is already open).
5. **A second, independent angle reading** cross-checked every trial, plus the motion and sample-counter checks, to catch a stuck or frozen reading. Motion, deceleration and sample-counter checks together still leave {{faults.stuck_during_drop.missed_after_checks|.1%}} of stuck-during-drop faults missed.
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

## 10. Reproduce

    python scripts/angle_detector_risk.py --check     # about 100 s on one core; exit 0 if both files reproduce
    pytest tests/test_angle_detector_risk.py          # fast checks (model vs PR #9, geometry, chain, rendering)
