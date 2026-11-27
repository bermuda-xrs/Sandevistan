# B-A15 arm-drop test plan (written plan, no trials run)

**Status: draft for Fact Checker. No human or bench trial has been run. No human trial may start before ethics (IRB) approval (section 2).**

- Every number comes from `scripts/ba15_numcheck.py` (output: `docs/ba15_numbers.json`). This file is rendered from `scripts/ba15_testplan.template.md`; do not edit it by hand.
- Numbers marked **SYNTHETIC** come from the drop model or from an extrapolation through it. They are not bench or human data. Sample-size numbers are exact binomial calculations, not simulations.
- **ASSUMPTION** marks anything with no source. The full list is in section 8.
- The window edges and FC's drop range are SYNTHETIC wherever they appear.
- "FC" is Fact Checker, "HE" is Hardware Engineer, "STO" is safe torque off. Citations are in the References.

## Short version

1. **Window (Charlie, 1a):** one fixed relaxed-hold catch-time window for every subject, 201.90–203.72 ms (SYNTHETIC). No per-subject recompute.
2. **Pass rule (Charlie, 2b):** show P(drop past 10°) ≤ 5 % at 95 % one-sided confidence with zero drops in 59 independent trials.
3. **Design (Charlie, 2026-11-26 01:25 KST): 59 subjects × one fault trial each.** Its false-pass chance is 4.85 % at any ICC (section 4). The earlier 33 × 10 default is in "Rejected alternatives" (section 4b).
4. **Any drop:** stop, investigate, and the test has FAILED. No extension unless it was registered before trial 1.
5. **Dummy phase first:** 117 cells × 10 drops = 1,170 drops, about 19.5 h (ASSUMPTION: 1 min per reset). Timing uncertainty is a requirement, at most 0.160 ms (SYNTHETIC budget, section 3), until it is measured on the rig.
6. **OPEN SAFETY ITEM THAT BLOCKS PHASE B:** there is no sourced human tolerance limit for an arm hitting the 15° stop (section 2). Phase B cannot start until there is one.
7. **Feasibility: this plan will almost certainly fail on the current design.** If the true drop probability is FC's 20 %–33 %, the chosen design passes with chance 1.9e-06 to 5.5e-11 (SYNTHETIC). Section 6 gives the bound over the other options and says what it means.

## 1. What is tested

B-A15 is the project's human catch model: the user arrests the arm at exactly 200 ms after fault onset. It has no human data behind it (ASSUMPTION). The 10° drop limit (B9) has no human-safety source either (ASSUMPTION): the S7 search found none and FC verified that (`S7_verdict.md` row 14).

**The fixed window.** Arm horizontal (elbow 90°), user hold torque zero ("relaxed"). The drop model (torque basis: ramp on the device side, times k = 1.5862 at the arm) reaches 10° at:

| Edge (SYNTHETIC) | Converged value | Decided value | Governing cell |
|---|---|---|---|
| Lower: design (b), SYNTHETIC | 201.9048 ms | 201.90 ms | e-stop, ramp 21.00 N·m/s, STO 148 ms |
| Upper: backup (iv), SYNTHETIC | 203.7165 ms | 203.72 ms | comm loss, ramp 19.95 N·m/s, STO 142 ms |

- Both decided values (SYNTHETIC) reproduce at two decimals. The lower one is rounded down by 0.0048 ms (safe side); the upper one is rounded up by 0.0035 ms (not the safe side, but far below the 0.1 ms timing resolution).
- The lower edge is the earliest crossing over all 21 (b) cells (path × rate tolerance × STO band). The next paths cross later: comm loss 201.96 ms, EMG loss 205.87 ms (SYNTHETIC).
- Margin over the 200 ms catch: 1.90 ms. Drop at 200 ms: 9.69°. Arm speed at 10°: 165 °/s (all SYNTHETIC).
- This decision replaces, for this test only, the Rev U row 14 rule that recomputes the window from each measured hold torque. One reason it cannot be recomputed well: a ±5 % error in a subject's gravity moment (ASSUMPTION) is 0.49 N·m, 1.66 times the whole nominal hold of 0.297 N·m (SYNTHETIC).
- **Limit of one window.** A subject whose gravity moment is one percent above the model crosses 10° -2.15 ms from the edge (minus is earlier); one percent more inertia moves it +0.60 ms (SYNTHETIC). The window is only 1.81 ms wide (SYNTHETIC). So the trial score below also uses the directly measured angle.

## 2. Ethics and safety before any trial

- **Human testing needs ethics (IRB) approval before any trial**, including pilots and familiarisation. The Declaration of Helsinki, paragraph 23, requires the protocol to be approved by a research ethics committee before the research begins [WMA24]. Where US rules apply, the IRB approves, requires changes to, or disapproves the research [CFR46]. Which law applies here is not known (ASSUMPTION: whichever applies, approval comes first).
- The IRB must see section 6. Exposing people to repeated arm drops in a test that is expected to fail is a risk with little expected benefit.
- Inclusion and exclusion criteria, consent text, medical cover and the adverse-event definition are for the IRB protocol. This plan does not set them (open item).
- **Rig safety (HE scoping; ASSUMPTIONS until built and proven on the dummy).** A person may not use the rig until every item is proven:
  1. **Rig power cut, independent of the device e-stop.** One cut removes power from the rig's latch, hold emulator and drives, and works with the device e-stop bypassed or failed.
  2. **Interlocked guard** over the 0.425 m sweep (payload tip 0.40 m plus 25 mm clearance). The tip sweeps 349 mm of arc between 65° and 115° (SYNTHETIC geometry). Opening the guard cuts rig power.
  3. **Passive lower stop** at 15° drop with a 10° absorber stroke, sized in the table below. In the Phase B test cell an uncaught arm reaches it 26.58 ms after its 10° crossing, at 228.49 ms. The earliest contact of any cell is 228.16 ms, in (b) comm loss 21.00 N·m/s, STO 142 ms, 26.20 ms after that cell's own 10° crossing (SYNTHETIC).
  4. **Upper runaway stop** at +15° above the start, for a device that drives upward.
  5. **Torque trip** at 11.0 N·m at the arm (SYNTHETIC k·T0 basis, ASSUMPTION).
  6. **Stop holder (requirement).** One person, named in the IRB protocol, whose only task during a trial is to hold the rig power cut and use it. Not the operator running the trial.
  7. **Subject stop button (requirement).** A stop the subject holds in the free hand, which cuts rig power. "Subject asks to stop" must not depend on speech alone.

| SYNTHETIC rig case | Energy | Speed or force |
|---|---|---|
| No-hold free fall, device torque lost at t = 0 (stop at 130.3 ms; 10° at 106.3 ms) | 2.665 J at the stop, 4.352 J through the stroke | 229 °/s at the stop |
| Device driving the arm down at its cap (HE's runaway case, k × 7.899 N·m = 12.53 N·m) | 5.945 J at the stop, **9.82 J through the stroke** | 342 °/s at the stop |
| Absorber, constant force over the stroke | ideal 56.3 N·m; peak 84.4 N·m (× 1.5, ASSUMPTION) | 422 N at 0.20 m |
| Upper runaway stop at +15°, device driving up at its cap | 0.693 J | 117 °/s; ideal 16.8 N·m, 126 N peak |
| Rig torque trip (+10 % over k × T0) | 11.0 N·m at the arm | ASSUMPTION |

  The absorber is sized by the larger case. The no-hold free fall sets the lowest case; with the device driving at its cap the stroke must absorb 9.82 J, against 4.35 J for free fall (SYNTHETIC). The stop is hit at 342 °/s in the larger case, not the 229 °/s of free fall (SYNTHETIC).

- **OPEN SAFETY ITEM THAT BLOCKS PHASE B: no human tolerance limit for hitting the 15° stop.**
  - The absorber above is sized for the rig, not for a person.
  - In Phase B an uncaught arm reaches the stop at 211 °/s in the test cell (SYNTHETIC). In the runaway case it is 342 °/s, with 84.4 N·m and 422 N at the contact (SYNTHETIC).
  - With the S7 catch model (normal, SD 16 ms), the share of Phase B trials whose catch comes after the test cell's stop contact is 4.34 % (nominal-hold mean) to 6.69 % (zero-hold mean) (SYNTHETIC). FC's 4.53 %–6.96 % uses the earliest contact of any cell (SYNTHETIC).
  - The stop rules halt at the first DROP, so at most one person hits the stop. But that person still takes it: the chance that the study ends on a stop contact is 9.05 % to 11.86 % when a drop is a catch time later than the window edge, and 8.97 % to 11.77 % when the edge is lowered by U_t, as the DROP rule in section 3 does (SYNTHETIC). Both assume exactly 59 released fault trials: every released INVALID fault trial that step 8 of Phase B replaces with a new subject adds one more exposure that can end on a stop contact.
  - **No sourced limit exists in this plan.** FC notes that ISO/TS 15066 has body-region limits but holds no clause text (`PR9_verdict.md` row 12); the S7 search could not read it either. No other source was found. **No limit is invented here.**
  - **Phase B is blocked until** a sourced human tolerance for this contact (speed, force or energy at the forearm) is in the plan, and the stop is shown on the dummy to stay within it.

## 3. Procedure

**Phase A: dummy arm, no people.**
1. **Build the rig** with a dummy forearm and the device fitted, start latch, hold emulator, angle and torque sensors, and every item of section 2. The dummy plus device must reproduce the model's inertia J = 0.3331 kg·m² and gravity moment G = 10.297 N·m (SYNTHETIC) within the calibration limits below.
2. **Prove the safety items** of section 2 before any drop with the device on, including the rig power cut with the device e-stop bypassed.
3. **Prove the timing budget.** Every term below is measured on the rig, not taken from a data sheet:

| SYNTHETIC timing budget item | Requirement | Derivation |
|---|---|---|
| Smallest model margin | 1.905 ms | governing cell, over all 117 cells |
| Budget | 0.190 ms | 10 % of the margin |
| Sampling | at least 15 kHz | period at most 0.076 ms (40 %) needs at least 13.1 kHz; rounded up to 5 kHz steps |
| Trigger skew | at most 0.019 ms | 10 % |
| Release skew | at most 0.038 ms | 20 % |
| Encoder latency | at most 0.019 ms | 10 % |
| Encoder bandwidth | at least 18.4 kHz | 0.35 / encoder latency (first-order rise time, ASSUMPTION) |
| Angle resolution | at most 0.0031° (17-bit: 131,072 counts, 0.00275°) | 10 % at the slowest 160.4 °/s crossing; 0.0171 ms |
| Spare | 0.019 ms | 10 % |
| **Used, U_t** | **at most 0.160 ms** | 0.0667 + 0.0190 + 0.0381 + 0.0190 + 0.0171 ms, against 0.190 ms |
| Headroom | 0.030 ms | 0.190 − 0.160 ms |
| Angle uncertainty | at most 0.0263° | encoder step plus 164.7 °/s times the timing terms |

   Every term is a **requirement** until it is measured on the rig. If a measured term is over its share, Phase A stops until it is fixed. In Phase B, a trial whose measured release skew is over 0.038 ms is INVALID (SYNTHETIC limit). A late release holds the arm while the torque ramps down, which makes the trial easier. A drop in such a trial still triggers the stop rules. The torque and force channels have no bandwidth requirement yet (open item).
   The sampling period of 0.0667 ms is also within the 0.1 ms timing resolution that Rev U row 14 asks for (a tenth of the 1.81 ms window, rounded down; SYNTHETIC).
4. **Calibrate** the dummy. The flip values show how little the model can be off before a failing unit would read as passing; the calibration limits come from splitting the 0.5 ms model agreement (ASSUMPTION) equally over four inputs (25 % each):

| SYNTHETIC governing cell, hold 0 | Flip value (drop at the catch reaches 10°) | Implied calibration (0.125 ms each) |
|---|---|---|
| Dummy human-side mass | +1.29 % |  |
| Inertia J | -3.12 % | ±0.21 % |
| Gravity moment G |  | ±0.057 % (±0.0059 N·m) |
| Net pre-onset torque, downward offset | 0.091 N·m | ±0.0059 N·m |
| Clamp T0 | -0.98 % |  |
| Torque ratio k | -1.22 % | ±0.080 % |
| Timing lag | more than 1.90 ms | U_t 0.160 ms |

   The limits for G and for net pre-onset torque are one error, counted twice (conservative): 0.057 % of G is 0.0059 N·m (SYNTHETIC).
5. **Firmware test mode.** The governing cell needs ramp 21.00 N·m/s and STO 148 ms, the edge of the tolerance band, not the product settings (20.00 N·m/s, 150 ms). That needs a test mode in the device firmware. **ASSUMPTION and open risk: it must be verified** that such a mode can set these values without changing the fault path, and the ramp and STO actually delivered must be measured and logged on every trial. A mode that changes the fault path makes the whole test invalid.
6. **Dummy matrix.** Every cell is run 10 times: 117 cells × 10 = 1,170 drops, about 19.5 h at 1 min per reset (ASSUMPTION). A cell is one of 21 (b) cells and 18 (iv) cells (path × ramp −5 %, nominal, +5 % × STO band edges and nominal), at three holds: none, half and the full nominal hold (0.000, 0.149, 0.297 N·m; SYNTHETIC). That covers every D-4 corner, and the upper-edge cell of the window: (iv), comm loss, 19.95 N·m/s, STO 142 ms, crossing 203.72 ms at hold 0 (SYNTHETIC). Model crossing times:

| SYNTHETIC cell: option, path, ramp, STO | Crossing at hold 0 (ms) | at 0.149 N·m (ms) | at 0.297 N·m (ms) |
|---|---|---|---|
| (b) estop 19.00 N·m/s, STO 148 ms | 205.69 | 208.98 | 212.36 |
| (b) estop 19.00 N·m/s, STO 150 ms | 206.34 | 209.66 | 213.09 |
| (b) estop 19.00 N·m/s, STO 160 ms | 209.30 | 212.81 | 216.42 |
| (b) estop 20.00 N·m/s, STO 148 ms | 203.78 | 207.02 | 210.35 |
| (b) estop 20.00 N·m/s, STO 150 ms | 204.39 | 207.66 | 211.02 |
| (b) estop 20.00 N·m/s, STO 160 ms | 207.11 | 210.55 | 214.09 |
| (b) estop 21.00 N·m/s, STO 148 ms | 201.90 | 205.09 | 208.37 |
| (b) estop 21.00 N·m/s, STO 150 ms | 202.46 | 205.68 | 208.99 |
| (b) estop 21.00 N·m/s, STO 160 ms | 204.96 | 208.33 | 211.81 |
| (b) comm 19.00 N·m/s, STO 142 ms | 205.33 | 208.55 | 211.86 |
| (b) comm 19.00 N·m/s, STO 150 ms | 208.16 | 211.53 | 214.99 |
| (b) comm 19.00 N·m/s, STO 160 ms | 211.27 | 214.82 | 218.47 |
| (b) comm 20.00 N·m/s, STO 142 ms | 203.63 | 206.81 | 210.08 |
| (b) comm 20.00 N·m/s, STO 150 ms | 206.28 | 209.60 | 213.01 |
| (b) comm 20.00 N·m/s, STO 160 ms | 209.15 | 212.64 | 216.22 |
| (b) comm 21.00 N·m/s, STO 142 ms | 201.96 | 205.09 | 208.31 |
| (b) comm 21.00 N·m/s, STO 150 ms | 204.43 | 207.69 | 211.05 |
| (b) comm 21.00 N·m/s, STO 160 ms | 207.07 | 210.49 | 214.02 |
| (b) emg 47.50 N·m/s, STO none | 210.50 | 213.82 | 217.24 |
| (b) emg 50.00 N·m/s, STO none | 208.09 | 211.31 | 214.63 |
| (b) emg 52.50 N·m/s, STO none | 205.87 | 209.00 | 212.23 |
| (iv) estop 18.05 N·m/s, STO 148 ms | 207.52 | 210.86 | 214.30 |
| (iv) estop 18.05 N·m/s, STO 150 ms | 208.22 | 211.60 | 215.08 |
| (iv) estop 18.05 N·m/s, STO 160 ms | 211.43 | 214.99 | 218.66 |
| (iv) estop 19.00 N·m/s, STO 148 ms | 205.69 | 208.98 | 212.36 |
| (iv) estop 19.00 N·m/s, STO 150 ms | 206.34 | 209.66 | 213.09 |
| (iv) estop 19.00 N·m/s, STO 160 ms | 209.30 | 212.81 | 216.42 |
| (iv) estop 19.95 N·m/s, STO 148 ms | 203.88 | 207.12 | 210.45 |
| (iv) estop 19.95 N·m/s, STO 150 ms | 204.48 | 207.76 | 211.13 |
| (iv) estop 19.95 N·m/s, STO 160 ms | 207.22 | 210.66 | 214.21 |
| (iv) comm 18.05 N·m/s, STO 142 ms | 206.96 | 210.22 | 213.58 |
| (iv) comm 18.05 N·m/s, STO 150 ms | 209.98 | 213.39 | 216.91 |
| (iv) comm 18.05 N·m/s, STO 160 ms | 213.32 | 216.93 | 220.64 |
| (iv) comm 19.00 N·m/s, STO 142 ms | 205.33 | 208.55 | 211.86 |
| (iv) comm 19.00 N·m/s, STO 150 ms | 208.16 | 211.53 | 214.99 |
| (iv) comm 19.00 N·m/s, STO 160 ms | 211.27 | 214.82 | 218.47 |
| (iv) comm 19.95 N·m/s, STO 142 ms | 203.72 | 206.89 | 210.17 |
| (iv) comm 19.95 N·m/s, STO 150 ms | 206.38 | 209.69 | 213.11 |
| (iv) comm 19.95 N·m/s, STO 160 ms | 209.25 | 212.74 | 216.33 |

   18 of the 117 cells (6 types × 3 holds) repeat a setting already in the matrix, because the (b) −5 % rate equals the (iv) nominal rate: e-stop and comm loss at 19.00 N·m/s, all three STO values (SYNTHETIC count). They are kept as replication: 180 drops, 3.0 h. Dropping them leaves 99 distinct cells, 990 drops and 16.5 h. HE decides.
   The 10 per cell is the Rev U row 14 floor; a trial SD of 0.15 ms and a wanted precision of 0.1 ms on the mean ask for 8.64 drops (both ASSUMPTIONS; SYNTHETIC formula).
   **Per-trial pass:** crossing at or after 200.16 ms (200 ms plus U_t), **or** drop at 200 ms at most 9.974° (10° minus the angle uncertainty) (SYNTHETIC thresholds). Either reading may be used; a trial where only one passes is reported. **Per cell:** the mean crossing agrees with the model within ±0.5 ms (ASSUMPTION). Every model cell passes the per-trial rule with room: the worst model drop at 200 ms is 9.689° (SYNTHETIC).
7. Any cell that fails stops Phase A. If the governing cell's crossing is outside the model by more than the agreement, the fixed window is wrong: stop and go to Charlie (via Seren).

**Gate to Phase B:** Phase A passed; **the open safety item of section 2 closed** (a sourced human tolerance for stop contact, met on the dummy); stop holder named and subject stop button fitted; IRB approval in hand; allowed drops registered (default zero); and Charlie's explicit go after reading section 6.

**Phase B: people.**
1. One test cell only: e-stop fault, ramp 21.00 N·m/s, STO 148 ms, relaxed hold. It has the earliest model crossing, so a given catch time is hardest to pass there (ASSUMPTION: a subject's catch time does not depend on which path failed). **A pass covers only the cell tested.** To claim each path separately, each needs its own full sample.
2. Seat the subject, elbow at 90°, payload in hand, forearm resting on the start latch. A relaxed arm would sag before the fault, because the device carries slightly less than the arm's weight, so the latch must hold it and release at fault onset (ASSUMPTION: release adds no usable cue).
3. Tell the subject to relax and to stop the arm only if it falls. Use no countdown and mask sound.
4. Random wait of 1 to 4 s before onset, as in [PRU09]. A sham is a release where the latch opens and no fault follows. Sham scheduling is fixed before trial 1 (ASSUMPTION, a design choice, not a computed number): before the session each subject is randomised, with equal chance, to one of two schedules, **fault only** (the fault is the subject's first release) or **one sham, then the fault**. A sham is always before the fault, never after it. The subject is not told how many releases to expect or in which order. A relaxed arm still sags slowly in a sham, so shams are never scored.
5. Fault onset is time zero. The catch time is from onset to arrest complete (arm speed zero), as in Rev U row 14.
6. **One fault trial per subject** (Charlie's design). The session ends right after the fault trial and the check on the subject, so nothing is released after it: no sham follows the fault trial.
7. No practice with the real fault. The unprepared user is the hazard, and [PRU09] subjects were practised. With one fault trial per subject, every fault is a first exposure.
8. **Re-running an INVALID trial** (this keeps the one-fault rule and the first-exposure basis). A trial found INVALID **before** the release, when nothing was released, may be restarted on the same subject, who is still unexposed. An INVALID **fault** trial where the latch was released is **not repeated on that subject**: a second fault would not be a first exposure and would make that subject two trials, which the ICC-free basis of section 4 does not cover. The subject is closed, the trial is not counted, and a new subject is enrolled in their place, so that the plan ends with 59 valid first-exposure fault trials. An INVALID sham is not repeated and the subject goes on to the fault. Every INVALID trial and every replacement is reported. The replacement is not chosen after seeing the angle trace, because INVALID is decided without it (see the score rules). A drop in an INVALID trial is not removed: the stop rules and the zero-drop failure rule still apply.

**Trial score (fixed before trial 1).** Let t_c be the measured catch time and U_t the measured timing uncertainty.

| Score | Rule |
|---|---|
| **DROP** | t_c + U_t > 201.90 ms (SYNTHETIC edge), where U_t is the measured timing uncertainty, required to be at most 0.160 ms (SYNTHETIC budget requirement), **or** measured peak drop + angle uncertainty (0.0263°) > 10°, **or** the passive stop is touched |
| **NO DROP** | valid trial and none of the above |
| **INVALID** | rig fault (ramp slower than 21.00 N·m/s, STO later than 148 ms, release skew over 0.038 ms (SYNTHETIC limit), lost sensor or clock), or EMG activity before onset showing the subject was not relaxed |

- Using both time and angle makes the score stricter than either alone, so a pass holds for both readings of "drop past 10°" (ASSUMPTION: this is how 1a and 2b combine; for FC and Charlie to confirm).
- Not relaxed means the repo EMG onset detector, at its unchanged defaults, fires during the wait (ASSUMPTION).
- INVALID is decided from pre-onset and rig channels only, by someone who has not seen the angle trace after onset. An INVALID trial is not counted. It is re-run or replaced as in Phase B step 8, never repeated as a second fault on the same subject, and all are reported. **A drop in an invalid trial still triggers the stop rules.**
- The rig may not be easier than the cell. Going past it by one percent in ramp rate moves the crossing -0.39 ms, and STO one ms early moves it -0.29 ms (SYNTHETIC). The allowed overshoot is open.

## 4. Sample size

**Independent trials.** Zero drops in n shows P ≤ 5 % at 95 % one-sided confidence when the chance of zero drops at exactly 5 % is at most 5 % (the exact one-sided Clopper–Pearson bound [CP34]; zero-event case [HL83]). 59 is the smallest such n: at n = 59 the chance is 0.04849; one trial fewer gives 0.05105. **Verified.**

**Chosen design (Charlie, 2026-11-26 01:25 KST): 59 subjects × one fault trial each, 59 fault trials in total.**

**Repeated trials in one subject are correlated.** The design effect is DE = 1 + (m − 1)·ICC, for m trials per subject. Effective sample size = trials / DE [DBB81, RS92]. We checked the formula on the published example in [KMP04] (DE 1.527, effective size 83.8; the paper gives the same DE and a rounded size). Subjects = ceil(n·DE/m). Entries are subjects (total fault trials):

| Trials per subject m | ICC 0 | ICC 0.05 | ICC 0.1 | ICC 0.2 | ICC 0.3 | ICC 0.5 | ICC 0.8 | ICC 1 |
|---|---|---|---|---|---|---|---|---|
| 1 | 59 (59) | 59 (59) | 59 (59) | 59 (59) | 59 (59) | 59 (59) | 59 (59) | 59 (59) |
| 2 | 30 (60) | 31 (62) | 33 (66) | 36 (72) | 39 (78) | 45 (90) | 54 (108) | 59 (118) |
| 5 | 12 (60) | 15 (75) | 17 (85) | 22 (110) | 26 (130) | 36 (180) | 50 (250) | 59 (295) |
| 10 | 6 (60) | 9 (90) | 12 (120) | 17 (170) | 22 (220) | 33 (330) | 49 (490) | 59 (590) |
| 20 | 3 (60) | 6 (120) | 9 (180) | 15 (300) | 20 (400) | 31 (620) | 48 (960) | 59 (1180) |

**Why 59 × one is robust to the ICC.** The m = 1 row is 59 at every ICC. With one fault trial per subject, subjects' drop probabilities may differ in any way: the chance that one subject does not drop is still one minus the mean drop probability. So zero drops in 59 subjects at a true 5 % has chance 4.85 % at every ICC in the table, under both the beta-binomial and the worst-case mix (computed: 0.04849 to 0.04849). No ICC has to be assumed, estimated or registered. The other rows show what the ICC would cost with more trials per subject.

- **The earlier estimate of 6–33 subjects is reproduced**, for m = 10 and ICC 0, 0.2 and 0.5: 6, 17, 33. **It is not an upper bound.** It stops at ICC 0.5 and assumes m = 10. Over the full ICC range 0 to 1 it becomes 6 to 59. No source gives the ICC of arm drops (ASSUMPTION).

## 4b. Rejected alternatives

**33 subjects × 10 fault trials (the earlier planning default, DE 5.5, 330 fault trials). Not chosen: Charlie chose 59 × one (2026-11-26 01:25 KST, relayed by Seren; the rejection of this default is read from that choice, not a separate statement).** Reasons:
1. **It is only valid if the drop ICC is at most 0.5, and that ICC must be fixed before trial 1.** A passing data set cannot tell you its ICC: with zero drops there is no variation to estimate it from. A pilot could measure the ICC of catch time. If catch times are normal (ASSUMPTION), a catch-time ICC of 0.25, 0.5, 0.75 or 0.9 gives a drop ICC of 0.08, 0.20, 0.41 or 0.62 at 5 %.
2. **Above an ICC of 0.5104, its 330 trials are worth fewer than 59 independent ones, and the false-pass chance at a true 5 % goes over 5 %:**

| True drop ICC | False pass, worst-case mix | False pass, beta-binomial | Chosen 59 × one |
|---|---|---|---|
| 0.5 | 3.69 % | 0.75 % | 4.85 % |
| 0.6 | 6.24 % | 1.73 % | 4.85 % |
| 0.8 | 12.21 % | 6.53 % | 4.85 % |
| 1 | 18.40 % | 18.40 % | 4.85 % |

3. **A high ICC is plausible here.** By the mapping in reason 1, a catch-time ICC of 0.75 to 0.9 already gives a drop ICC of 0.41 to 0.62.
4. **Too small an ICC breaks it the other way.** Sizing for ICC 0 (6 subjects) when the truth is 0.5 gives a false-pass chance of 41.04 % (beta-binomial) to 54.87 % (worst case), not 5 %.
5. Repeated faults teach the subject what is coming; the unprepared response is the hazard (section 3).

For the record, the design-effect check for m = 10:

**Is the design effect safe for a zero-drop test?** It is a variance correction, so we checked it exactly at P = 5 %, m = 10, against (a) a beta-binomial model [W75] and (b) the worst case over every possible mix of subject drop rates with that mean and ICC. The worst case (some subjects at one rate, the rest never dropping) is confirmed by a linear program (they agree within 1e-09):

| ICC | Design effect | Subjects, design effect (rejected) | Subjects, beta-binomial | Subjects, worst case | False-pass chance, plan size, worst case |
|---|---|---|---|---|---|
| 0 | 1.00 | 6 | 6 | 6 | 4.61 % |
| 0.05 | 1.45 | 9 | 8 | 8 | 2.76 % |
| 0.1 | 1.90 | 12 | 9 | 10 | 2.19 % |
| 0.2 | 2.80 | 17 | 11 | 14 | 2.51 % |
| 0.3 | 3.70 | 22 | 14 | 19 | 3.05 % |
| 0.5 | 5.50 | 33 | 21 | 30 | 3.69 % |
| 0.8 | 8.20 | 49 | 37 | 48 | 4.41 % |
| 1 | 10.00 | 59 | 59 | 59 | 4.85 % |

- The design effect is enough in every row. On a grid of ICC 0 to 1 in steps of one hundredth it stays enough for every m up to 46, and first falls short at the next m. **Do not use more than 46 trials per subject with this formula.**

**The other trial counts per subject** in the table above carry the same ICC dependence and are rejected for the same reasons.

## 5. Pass, fail and stop rules

**PASS:** every planned valid trial done, zero DROP scores.
**FAIL:** any DROP score.

**If a drop occurs**
1. **Stop** the session and all enrolment at once. Check on the subject. Report as the IRB requires.
2. **Investigate** before anything else: rig, ramp, STO, hold, catch time, angle.
3. **The zero-drop test has FAILED.** A drop is never removed afterwards. Report the exact one-sided upper bound: one drop in 59 independent trials gives P ≤ 7.79 %; two give 10.29 %.
4. **No extending after the fact.** "Add trials and allow one drop" chosen after seeing the drop has a false-pass chance of 7.48 % (59 then 93) or 9.26 % (59, 93, then 124), not 5 %.
5. **Where the drop falls** (from t_c; SYNTHETIC edges; decides who is told what, not a pass): above 201.90 and up to 203.72 ms, (b) has failed and (iv) becomes the candidate. (iv) needs Charlie's switch decision and its own new full test in its own cell. Above 203.72 ms both fail: escalate to Charlie (via Seren).
6. **Retest** only after a documented design change, from zero. Trying again without a change until it passes gives a false-pass chance of 9.46 % over two tries and 13.85 % over three.

**Allowed-drop plans** (must replace the zero-drop plan before trial 1, never after a drop). With one trial per subject, subjects = trials at any ICC:

| Allowed drops c | n (independent) | Size at 5 % | Subjects (one trial each) |
|---|---|---|---|
| 0 | 59 | 4.849 % | 59 |
| 1 | 93 | 4.998 % | 93 |
| 2 | 124 | 4.953 % | 124 |

**Registered extension that keeps the first look at 59** (also before trial 1 only). The overall false-pass chance is 0.0499997; with one trial per subject it does not depend on the ICC. It is costly because the first look already spends 0.0485:

| Look | Pass if total drops at most | Total subjects (one trial each) |
|---|---|---|
| 1 | 0 | 59 |
| 2 | 1 | 149 |
| 3 | 2 | 273 |

**Default: zero drops, no extension.** A registered extension also means more exposure after a known drop, so it needs explicit IRB approval.

**Other stop rules**

| Trigger | Action |
|---|---|
| Passive stop touched, pain, injury, or any adverse event | Stop the study. IRB report. Scored DROP if the stop was touched |
| Torque trip, ramp or STO outside the cell, lost sensor or clock | Stop the session, fix and re-prove on the dummy. Trial INVALID; drop rule still applies |
| Subject presses the stop button | Rig power is cut. Stop that subject; the trial is INVALID and the drop rule still applies |
| Subject asks to stop | Stop that subject. Keep and report their valid trials. Recruit a replacement to reach the planned size |
| Dummy check at start or end of a day outside tolerance | Stop. That day's trials go to FC before they count |
| All planned trials done with no drop | Stop. PASS. No further trials |

## 6. Feasibility before testing (SYNTHETIC)

**In plain words: if the true drop probability is 20 %–33 %, the chance of passing is effectively zero (SYNTHETIC).**
- **Chosen design (59 × one):** 1.9e-06 at 20 % and 5.5e-11 at 33 % (SYNTHETIC), at any ICC.
- **Every option in this plan's tables** (trials per subject, ICC, allowed drops and the registered extension): at most 2.4e-05, about one in 41,267 (SYNTHETIC). This bound holds **only for the tabulated options**.
- **Any trial count per subject up to 46** (where the design-effect check holds), any ICC in steps of 1/2000, zero drops, best-case mix: at most 2.1e-03 at 20 % and 5.7e-04 at 33 %, both at m = 46, ICC 0.012, 2 subjects (SYNTHETIC). All of those options are rejected (section 4b).

Chosen design and the allowed-drop sizes (independent trials, one per subject):

| SYNTHETIC true drop probability | c = 0, n = 59 (chosen) | c = 1, n = 93 | c = 2, n = 124 | Mean subjects to first drop |
|---|---|---|---|---|
| 20 % (FC low) | 1.9e-06 | 2.4e-08 | 4.9e-10 | 5.0 |
| 33 % (FC high) | 5.5e-11 | 3.1e-15 | 5.2e-19 | 3.0 |
| 47.95 % (window (b) edge) | 1.8e-17 | 3.6e-25 | 4.5e-32 | 2.1 |

Rejected clustered designs, subjects inflated for the ICC (m = 10). "Best case" is the mix of subjects most favourable to passing:

| SYNTHETIC: ICC | Subjects (trials) | 20 %: beta-binomial / best case | 33 %: beta-binomial / best case | 47.95 %: beta-binomial / best case |
|---|---|---|---|---|
| 0 | 6 (60) | 1.5e-06 / 1.5e-06 | 3.7e-11 / 3.7e-11 | 9.6e-18 / 9.6e-18 |
| 0.05 | 9 (90) | 9.5e-08 / 1.2e-06 | 4.1e-13 / 1.2e-09 | 2.1e-20 / 2.9e-12 |
| 0.1 | 12 (120) | 1.3e-08 / 8.7e-07 | 1.6e-14 / 7.7e-10 | 2.2e-22 / 8.2e-13 |
| 0.2 | 17 (170) | 3.1e-09 / 1.3e-06 | 1.8e-15 / 7.3e-10 | 1.4e-23 / 1.9e-13 |
| 0.3 | 22 (220) | 1.7e-09 / 1.7e-06 | 7.0e-16 / 5.3e-10 | 4.2e-24 / 3.9e-14 |
| 0.5 | 33 (330) | 1.3e-09 / 1.5e-06 | 4.7e-16 / 1.5e-10 | 2.3e-24 / 1.1e-15 |
| 0.8 | 49 (490) | 3.6e-08 / 1.6e-06 | 9.9e-14 / 6.2e-11 | 3.4e-21 / 4.9e-17 |
| 1 | 59 (590) | 1.9e-06 / 1.9e-06 | 5.5e-11 / 5.5e-11 | 1.8e-17 / 1.8e-17 |

- **More subjects do not help.** No row gives any real chance of passing. Under-sizing raises the chance only by breaking the plan: sized for ICC 0 with a true ICC of 0.5, a design that drops 20 % of the time would pass with chance 2.42 % under the beta-binomial and 8.78 % under the best-case mix (SYNTHETIC).
- **The first drop is expected almost at once:** on average at subject 3.0 to 5.0 (SYNTHETIC).
- **FC's 20 %–33 % is for the nominal hold, not this window.** We reproduce it: 19.86 %–32.57 % at mean 201.1 ms, SD 16 ms (SYNTHETIC; a linear extrapolation of [PRU09] outside its data). With the **relaxed** window Charlie chose and the same catch times it is 43.47 % at the (iv) edge and 47.95 % at the (b) edge (SYNTHETIC). Extending the trend to zero hold (mean 204.5 ms) gives 51.95 % and 56.44 % (SYNTHETIC). So FC's range is the optimistic case here.
- **No [PRU09] condition passes either.** Reading ± as SD and times as normal (ASSUMPTIONS), the best is 7.35 % drops and a 1.11 % chance of passing (SYNTHETIC):

| SYNTHETIC: PRU09 condition (mean ± SD, ms) | P(catch later than the (b) edge) | P(0 drops in 59) |
|---|---|---|
| reversal, 1 N·m background (193 ± 16) | 28.9 % | 1.8e-09 |
| reversal, 3 N·m background (170 ± 22) | 7.3 % | 1.1e-02 |
| maximal excursion, 1.25 N·m step (178 ± 22) | 13.9 % | 1.5e-04 |
| maximal excursion, 2.5 N·m step (184 ± 15) | 11.6 % | 6.8e-04 |

- **What would pass.** Chance of passing by true drop probability (exact, independent trials):

| True drop probability | Pass, c = 0 (n = 59) | Pass, c = 1 (n = 93) | Pass, c = 2 (n = 124) |
|---|---|---|---|
| 0.1 % | 94.27 % | 99.60 % | 99.97 % |
| 1 % | 55.27 % | 76.16 % | 87.16 % |
| 2 % | 30.36 % | 44.27 % | 54.77 % |
| 5 % | 4.85 % | 5.00 % | 4.95 % |

  A passing chance of four in five needs a true drop probability of 0.38 % or less; an even chance needs 1.17 %. **A design exactly at the 5 % target would almost always fail.** With an SD of 16 ms (SYNTHETIC use of [PRU09]), a true 5 % at the (b) edge needs a mean catch of 175.6 ms, and a four-in-five chance needs 159.2 ms (SYNTHETIC). At FC's mean catch the 10° crossing would have to move to 227.4 ms and 243.8 ms (SYNTHETIC).

**What it implies**
1. The test is a valid way to *reject* the design, and would do so quickly. It is not a realistic way to *accept* the current one.
2. Running it on people as planned exposes subjects to drops for a result already expected. **Recommendation: do Phase A only, and do not start Phase B on the current design.** Charlie decides.
3. FC's estimate is an extrapolation and could be wrong. A small IRB-approved pilot that measures catch time (not pass/fail) would test that at lower cost than Phase B. It is not part of this plan.
4. This plan proposes no design fix. Any fix must first be checked against the failures log.

## 7. Data logged per trial

| Group | Channels |
|---|---|
| Identity | coded subject, session, trial number, fault or sham, assigned sham schedule, random wait, test cell, operator, firmware commit, plan version |
| Commands | ramp rate and STO time commanded **and measured as delivered** (firmware test mode), fault-onset trigger, latch-release time |
| Torque | commanded and measured joint torque (device side), arm-side torque at the cuff, k used, torque-trip state |
| Hold | latch load before release (arm side), payload mass, subject arm mass estimate |
| Motion | elbow angle and speed, **raw encoder counts**, peak drop, drop at 200 ms, catch time t_c |
| Rig | **stop contact force and stroke**, **rig power-cut and guard-interlock state**, torque-trip state |
| EMG / IMU | biceps and triceps EMG and IMU at the repo defaults (2000 / 150 Hz, unchanged), pre-onset EMG onset flag |
| Timing | **fault and release timestamps on one clock, with the measured skew between them**, clock-sync error, sample-drop count, U_t in force |
| Dummy | **mass, inertia and gravity-moment calibration records** (before and after each day), the matrix cell and the model crossing it is compared with |
| Outcome | DROP / NO DROP / INVALID, reason, who decided and when, blinding flag, adverse event, comfort report |

All rig channels are logged at 15 kHz or faster (the budget in section 3; one sample every 0.0667 ms, SYNTHETIC). The 150 Hz IMU has one sample every 6.7 ms and cannot score the criterion; it is logged only as a cross-check at the repo defaults. Raw data are kept. Nothing is deleted. Logged data rate (SYNTHETIC; ASSUMPTION A21): 14 rig channels at 15 kHz with 4-byte samples is 0.84 MB/s, or 8.4 MB for a 10 s trial record.

## 8. Assumptions (unsourced unless stated)

| # | ASSUMPTION |
|---|---|
| A1 | The 10° drop limit (B9). No human-safety source |
| A2 | The 200 ms catch. No human data |
| A3 | Drop model: instant catch, zero elbow damping, HE's arm and device parameters, k = 1.5862 |
| A4 | Zero hold torque represents a relaxed user, and one window fits every subject (Charlie 1a) |
| A5 | ICC of drops: not needed by the chosen design; used only to compare designs (0.5 was the rejected default) |
| A6 | Not an assumption: one fault trial per subject is Charlie's decision (2026-11-26). It is below the Rev U row 14 floor of 10 trials per subject (open item) |
| A7 | Subjects are independent and represent the users; subject population not set |
| A8 | Testing only the governing cell is enough; catch time does not depend on the path |
| A9 | Score combines catch time and angle, with uncertainty added against a pass |
| A10 | Subject gravity moment known to ±5 % (HE scoping) |
| A11 | Latch release and shams give no usable cue. Weaker under 59 × one: with the sham-then-fault schedule a subject who works out that the second release is real is prepared, which (ASSUMPTION) makes a drop less likely and so favours a false pass. Not tested; open item 9 |
| A12 | "Relaxed" means the EMG onset detector at unchanged defaults stays silent |
| A13 | (SYNTHETIC) FC's 20 %–33 %: linear extrapolation of [PRU09], normal catch times, SD 16 ms, ± read as SD |
| A14 | [PRU09] transfers to this task (step torque, horizontal plane, practised subjects there; ramp, vertical, unprepared here) |
| A15 | Timing budget: 10 % of the smallest model margin, split 40 % sampling, 10 % trigger, 20 % release, 10 % angle, 10 % encoder latency, 10 % spare; sampling rounded up to 5 kHz steps (HE, reusing the Rev U row-8 rule) |
| A16 | Model agreement of 0.5 ms, split equally over four inputs (net torque, J, G, k) |
| A17 | Dummy trial SD 0.15 ms, wanted precision 0.1 ms, 1 min per reset |
| A18 | Rig: stop at 15° with a 10° stroke, upper stop at +15°, contact radius 0.20 m, clearance 25 mm, absorber peak 1.5 × the ideal force, torque trip 10 % over k·T0, device-at-cap runaway case (HE) |
| A19 | **Firmware test mode** that sets ramp 21.00 N·m/s and STO 148 ms without changing the fault path. Open risk: **must be verified** |
| A20 | Encoder bandwidth from the first-order rise-time rule, bandwidth = 0.35 / rise time, with the encoder-latency share as the rise time (HE's method) |
| A21 | (SYNTHETIC) Data rate: 14 rig channels (FC's count of the channels named in section 7), 4 bytes per sample, a 10 s record. Not HE's figure |

## 9. Open items

1. **OPEN SAFETY ITEM, BLOCKS PHASE B:** a sourced human tolerance limit for stop contact (section 2), and proof on the dummy that the stop stays within it.
2. Charlie: run Phase B at all, given section 6?
3. Charlie: which paths need their own sample? (Design decided: 59 × one, 2026-11-26.)
4. FC / Charlie: confirm the combined score (A9).
5. HE: allowed ramp and STO overshoot on the rig; latch design and whether the release skew can be met; absorber proof on the dummy; bandwidth for the torque and force channels. Rev U row 14 still says the window is recomputed per hold torque and asks for 10 trials per subject; both differ from Charlie's decisions.
6. HE / firmware: **the test mode (A19) must be verified**: ramp 21.00 N·m/s and STO 148 ms delivered as set, fault path unchanged.
7. Before Phase B: name the stop holder (IRB protocol) and fit the subject stop button.
8. IRB: criteria, consent (what the subject is told about releases), adverse events, rest times. The sham schedule is set in Phase B step 4 (ASSUMPTION) and left open in item 9.
9. Charlie / FC: sham schedule (Phase B step 4) and its cue risk (A11). Zero shams (every subject fault-only) is the other option: it removes the sham cue. Outcomes by schedule are reported descriptively; they do not change the pass rule. Not decided here.

## References

Bibliographic data for the DOI entries were checked on Crossref on 2026-11-23. "Text read" says how much was read.

- [PRU09] Pruszynski JA, Kurtzer I, Lillicrap TP, Scott SH (2009). Temporal evolution of "automatic gain-scaling". J Neurophysiol 102(2):992–1003. doi:10.1152/jn.00085.2009, PMC2724331. Text read: Methods (1–4 s hold, about five familiarisation repeats) and Results, "Behavior and kinematics" (times). The ± on the times is not labelled there.
- [CP34] Clopper CJ, Pearson ES (1934). Biometrika 26(4):404–413. doi:10.1093/biomet/26.4.404. Metadata only.
- [HL83] Hanley JA, Lippman-Hand A (1983). If nothing goes wrong, is everything all right? JAMA 249(13):1743. doi:10.1001/jama.1983.03330370053031. Metadata only.
- [DBB81] Donner A, Birkett N, Buck C (1981). Randomization by cluster: sample size requirements and analysis. Am J Epidemiol 114(6):906–914. doi:10.1093/oxfordjournals.aje.a113261. Metadata only.
- [RS92] Rao JNK, Scott AJ (1992). A simple method for the analysis of clustered binary data. Biometrics 48(2):577. doi:10.2307/2532311. Abstract read (design effect and effective sample size).
- [KMP04] Killip S, Mahfoud Z, Pearce K (2004). What is an intracluster correlation coefficient? Ann Fam Med 2(3):204–208. doi:10.1370/afm.141, PMC1466680. Read through a fetch agent: worked example m = 32, ICC 0.017, DE 1.527, 128 subjects, effective size 84. The formula itself is an image there, so we checked it by that example. One sentence, as fetched, gives DE 1 at ICC 1; that conflicts with its own statement that the effective size falls to the number of clusters, which needs DE = m as used here.
- [W75] Williams DA (1975). The analysis of binary responses from toxicological experiments involving reproduction and teratogenicity. Biometrics 31(4):949. doi:10.2307/2529820. Abstract read (beta-distributed response probability within clusters).
- [WMA24] World Medical Association. Declaration of Helsinki, October 2024 revision, paragraph 23. https://www.wma.net/policies-post/wma-declaration-of-helsinki/ (read through a fetch agent 2026-11-23). doi:10.1001/jama.2024.21972.
- [CFR46] 45 CFR 46.109(a), https://www.law.cornell.edu/cfr/text/45/46.109 (read in part through a fetch agent; eCFR and HHS refused the fetch).
- Project files (not in this repo): `exo/calc_BA15_testplan.py` and `exo/BA15_testplan_out.json` (HE scoping; rerun byte-identical), `exo/SN5_BenchC_rules_RevU.md` (row 14, B-A15), `fc_check/sessions_r1/S7_verdict.md` (FC, rows 1, 12, 14, 16), `exo/benchC_rules_out.json` (k, T0).
