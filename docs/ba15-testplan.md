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
3. **Subjects:** with 10 trials each, 6 to 33 subjects for an ICC of 0 to 0.5; 59 for any ICC. Planning default: 33 subjects × 10 trials (ICC 0.5, ASSUMPTION).
4. **Any drop:** stop, investigate, and the test has FAILED. No extension unless it was registered before trial 1.
5. **Feasibility: this plan will almost certainly fail on the current design.** If the true drop probability is FC's 20 %–33 %, the chance of passing is at most 2.4e-05 (SYNTHETIC) (about one in 41,267) over every ICC, trials-per-subject and allowed-drop option in this plan. Section 6 says what that means.

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
- Passive stop below the arm at 15° (ASSUMPTION, HE scoping). In the governing cell an uncaught arm reaches it at 228.5 ms and 211 °/s (SYNTHETIC). The stop's absorber must be proven on the dummy before any person uses the rig.

## 3. Procedure

**Phase A: dummy arm, no people.**
1. Build the rig with a dummy forearm, start latch, angle sensor, torque sensor, passive stop and independent torque trip.
2. Calibrate. Prove timing uncertainty U_t ≤ 0.1 ms and angle uncertainty ≤ 0.0165° (SYNTHETIC: 0.1 ms at the model's 10° speed).
3. Run the governing cell. The dummy's measured 10° crossing must agree with the model within 0.5 ms (ASSUMPTION, HE scoping). If it does not, the fixed window is wrong: stop and go to Charlie (via Seren).
4. Prove the passive stop and the torque trip.

**Gate to Phase B:** Phase A passed, IRB approval in hand, ICC / trials per subject / allowed drops registered, and Charlie's explicit go after reading section 6.

**Phase B: people.**
1. One test cell only: e-stop fault, ramp 21.00 N·m/s, STO 148 ms, relaxed hold. It has the earliest model crossing, so a given catch time is hardest to pass there (ASSUMPTION: a subject's catch time does not depend on which path failed). **A pass covers only the cell tested.** To claim each path separately, each needs its own full sample.
2. Seat the subject, elbow at 90°, payload in hand, forearm resting on the start latch. A relaxed arm would sag before the fault, because the device carries slightly less than the arm's weight, so the latch must hold it and release at fault onset (ASSUMPTION: release adds no usable cue).
3. Tell the subject to relax and to stop the arm only if it falls. Use no countdown and mask sound.
4. Random wait of 1 to 4 s before onset, as in [PRU09]. Mix in sham releases (latch opens, no fault) at random, one sham per fault trial (ASSUMPTION). A relaxed arm still sags slowly in a sham, so shams are never scored.
5. Fault onset is time zero. The catch time is from onset to arrest complete (arm speed zero), as in Rev U row 14.
6. Reset, rest, repeat up to 10 fault trials. Rest length is open.
7. No practice with the real fault. The unprepared user is the hazard, and [PRU09] subjects were practised. Trial number is logged so a learning trend can be seen.

**Trial score (fixed before trial 1).** Let t_c be the measured catch time and U_t the proven timing uncertainty.

| Score | Rule |
|---|---|
| **DROP** | t_c + U_t > 201.90 ms (SYNTHETIC edge), **or** measured peak drop + angle uncertainty > 10°, **or** the passive stop is touched |
| **NO DROP** | valid trial and none of the above |
| **INVALID** | rig fault (ramp slower than 21.00 N·m/s, STO later than 148 ms, lost sensor or clock), or EMG activity before onset showing the subject was not relaxed |

- Using both time and angle makes the score stricter than either alone, so a pass holds for both readings of "drop past 10°" (ASSUMPTION: this is how 1a and 2b combine; for FC and Charlie to confirm).
- Not relaxed means the repo EMG onset detector, at its unchanged defaults, fires during the wait (ASSUMPTION).
- INVALID is decided from pre-onset and rig channels only, by someone who has not seen the angle trace after onset. Invalid trials are repeated and all are reported. **A drop in an invalid trial still triggers the stop rules.**
- The rig may not be easier than the cell. Going past it by one percent in ramp rate moves the crossing -0.39 ms, and STO one ms early moves it -0.29 ms (SYNTHETIC). The allowed overshoot is open.

## 4. Sample size

**Independent trials.** Zero drops in n shows P ≤ 5 % at 95 % one-sided confidence when the chance of zero drops at exactly 5 % is at most 5 % (the exact one-sided Clopper–Pearson bound [CP34]; zero-event case [HL83]). 59 is the smallest such n: at n = 59 the chance is 0.04849; one trial fewer gives 0.05105. **Verified.**

**Repeated trials in one subject are correlated.** The design effect is DE = 1 + (m − 1)·ICC, for m trials per subject. Effective sample size = trials / DE [DBB81, RS92]. We checked the formula on the published example in [KMP04] (DE 1.527, effective size 83.8; the paper gives the same DE and a rounded size). Subjects = ceil(n·DE/m). Entries are subjects (total fault trials):

| Trials per subject m | ICC 0 | ICC 0.05 | ICC 0.1 | ICC 0.2 | ICC 0.3 | ICC 0.5 | ICC 0.8 | ICC 1 |
|---|---|---|---|---|---|---|---|---|
| 1 | 59 (59) | 59 (59) | 59 (59) | 59 (59) | 59 (59) | 59 (59) | 59 (59) | 59 (59) |
| 2 | 30 (60) | 31 (62) | 33 (66) | 36 (72) | 39 (78) | 45 (90) | 54 (108) | 59 (118) |
| 5 | 12 (60) | 15 (75) | 17 (85) | 22 (110) | 26 (130) | 36 (180) | 50 (250) | 59 (295) |
| 10 | 6 (60) | 9 (90) | 12 (120) | 17 (170) | 22 (220) | 33 (330) | 49 (490) | 59 (590) |
| 20 | 3 (60) | 6 (120) | 9 (180) | 15 (300) | 20 (400) | 31 (620) | 48 (960) | 59 (1180) |

- **The earlier estimate of 6–33 subjects is reproduced**, for m = 10 and ICC 0, 0.2 and 0.5: 6, 17, 33. **It is not an upper bound.** It stops at ICC 0.5 and assumes m = 10. Over the full ICC range 0 to 1 it becomes 6 to 59. No source gives the ICC of arm drops (ASSUMPTION).
- **Is the design effect safe for a zero-drop test?** It is a variance correction, so we checked it exactly at P = 5 %, m = 10, against (a) a beta-binomial model [W75] and (b) the worst case over every possible mix of subject drop rates with that mean and ICC. The worst case (some subjects at one rate, the rest never dropping) is confirmed by a linear program (they agree within 1e-09):

| ICC | Design effect | Subjects, design effect (plan) | Subjects, beta-binomial | Subjects, worst case | False-pass chance, plan size, worst case |
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
- **Too small an ICC breaks the test.** Sizing for ICC 0 (6 subjects) when the truth is 0.5 gives a false-pass chance of 41.04 % (beta-binomial) to 54.87 % (worst case), not 5 %.
- **A passing data set cannot tell you its ICC:** with zero drops there is no variation to estimate it from. So the ICC must be fixed before trial 1. A pilot could measure the ICC of catch time. If catch times are normal (ASSUMPTION), a catch-time ICC of 0.25, 0.5, 0.75 or 0.9 gives a drop ICC of 0.08, 0.20, 0.41 or 0.62 at 5 %.

**Planning default:** 33 subjects × 10 trials = 330 fault trials (DE 5.5). The 10 trials come from Rev U row 14 and the ICC of 0.5 is the top of the earlier range; both are ASSUMPTIONS.
**Recommended alternative for Charlie:** 59 subjects × one trial. It needs no ICC assumption, has the fewest fault trials (59) and tests the unprepared first response. It needs more people and does not meet Rev U's ten-trials-per-subject floor.

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

**Allowed-drop plans** (must replace the zero-drop plan before trial 1, never after a drop). m = 10, subjects by design effect:

| Allowed drops c | n (independent) | Size at 5 % | Subjects, ICC 0 | Subjects, ICC 0.05 | Subjects, ICC 0.1 | Subjects, ICC 0.2 | Subjects, ICC 0.3 | Subjects, ICC 0.5 | Subjects, ICC 0.8 | Subjects, ICC 1 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 59 | 4.849 % | 6 | 9 | 12 | 17 | 22 | 33 | 49 | 59 |
| 1 | 93 | 4.998 % | 10 | 14 | 18 | 27 | 35 | 52 | 77 | 93 |
| 2 | 124 | 4.953 % | 13 | 18 | 24 | 35 | 46 | 69 | 102 | 124 |

The design effect also covers the beta-binomial and the never/sometimes-drop mix for c = 1 and 2 at every ICC in the table (rows where it falls short: 0). For c above 0 that mix is not proven to be the worst case.

**Registered extension that keeps the first look at 59** (also before trial 1 only). The overall false-pass chance is 0.0499997 for independent trials, and at most 0.0485 with the subject counts below under either clustered model (m = 10). It is costly because the first look already spends 0.0485:

| Look | Pass if total drops at most | Total trials (independent) | Subjects, ICC 0 | Subjects, ICC 0.05 | Subjects, ICC 0.1 | Subjects, ICC 0.2 | Subjects, ICC 0.3 | Subjects, ICC 0.5 | Subjects, ICC 0.8 | Subjects, ICC 1 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | 59 | 6 | 9 | 12 | 17 | 22 | 33 | 49 | 59 |
| 2 | 1 | 149 | 15 | 22 | 29 | 42 | 56 | 82 | 123 | 149 |
| 3 | 2 | 273 | 28 | 40 | 52 | 77 | 102 | 151 | 224 | 273 |

**Default: zero drops, no extension.** A registered extension also means more exposure after a known drop, so it needs explicit IRB approval.

**Other stop rules**

| Trigger | Action |
|---|---|
| Passive stop touched, pain, injury, or any adverse event | Stop the study. IRB report. Scored DROP if the stop was touched |
| Torque trip, ramp or STO outside the cell, lost sensor or clock | Stop the session, fix and re-prove on the dummy. Trial INVALID; drop rule still applies |
| Subject asks to stop | Stop that subject. Keep and report their valid trials. Recruit a replacement to reach the planned size |
| Dummy check at start or end of a day outside tolerance | Stop. That day's trials go to FC before they count |
| All planned trials done with no drop | Stop. PASS. No further trials |

## 6. Feasibility before testing (SYNTHETIC)

**In plain words: if the true drop probability is 20 %–33 %, the chance of passing is effectively zero.** At most 2.4e-05, about one in 41,267 (SYNTHETIC), over every m, ICC and allowed-drop option in the tables, and the registered extension.

Independent trials:

| SYNTHETIC true drop probability | c = 0, n = 59 | c = 1, n = 93 | c = 2, n = 124 | Mean trials to first drop | Drop within one subject's 10 trials, ICC 0 / 0.5 |
|---|---|---|---|---|---|
| 20 % (FC low) | 1.9e-06 | 2.4e-08 | 4.9e-10 | 5.0 | 89.3 % / 46.2 % |
| 33 % (FC high) | 5.5e-11 | 3.1e-15 | 5.2e-19 | 3.0 | 98.2 % / 65.7 % |
| 48 % (window (b) edge) | 1.8e-17 | 3.6e-25 | 4.5e-32 | 2.1 | 99.9 % / 80.8 % |

With subjects inflated for the ICC (m = 10). "Best case" is the mix of subjects most favourable to passing:

| SYNTHETIC: ICC | Subjects (trials) | 20 %: beta-binomial / best case | 33 %: beta-binomial / best case | 48 %: beta-binomial / best case |
|---|---|---|---|---|
| 0 | 6 (60) | 1.5e-06 / 1.5e-06 | 3.7e-11 / 3.7e-11 | 9.6e-18 / 9.6e-18 |
| 0.05 | 9 (90) | 9.5e-08 / 1.2e-06 | 4.1e-13 / 1.2e-09 | 2.1e-20 / 2.9e-12 |
| 0.1 | 12 (120) | 1.3e-08 / 8.7e-07 | 1.6e-14 / 7.7e-10 | 2.2e-22 / 8.2e-13 |
| 0.2 | 17 (170) | 3.1e-09 / 1.3e-06 | 1.8e-15 / 7.3e-10 | 1.4e-23 / 1.9e-13 |
| 0.3 | 22 (220) | 1.7e-09 / 1.7e-06 | 7.0e-16 / 5.3e-10 | 4.2e-24 / 3.9e-14 |
| 0.5 | 33 (330) | 1.3e-09 / 1.5e-06 | 4.7e-16 / 1.5e-10 | 2.3e-24 / 1.1e-15 |
| 0.8 | 49 (490) | 3.6e-08 / 1.6e-06 | 9.9e-14 / 6.2e-11 | 3.4e-21 / 4.9e-17 |
| 1 | 59 (590) | 1.9e-06 / 1.9e-06 | 5.5e-11 / 5.5e-11 | 1.8e-17 / 1.8e-17 |

- **More subjects do not help.** No row gives any real chance of passing. Under-sizing raises the chance only by breaking the plan: sized for ICC 0 with a true ICC of 0.5, a design that drops 20 % of the time would pass with chance 2.42 % (SYNTHETIC).
- **The first drop is expected almost at once:** on average at trial 3.0 to 5.0 (SYNTHETIC, independent). At ICC 0.5, the chance of a drop within the first two subjects is 71.09 % to 88.23 %.
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
| Identity | coded subject, session, trial number, fault or sham, random wait, test cell, operator, firmware commit, plan version |
| Commands | ramp rate, STO time, fault-onset trigger, latch-release time |
| Torque | commanded and measured joint torque (device side), arm-side torque at the cuff, k used, torque-trip state |
| Hold | latch load before release (arm side), payload mass, subject arm mass estimate |
| Motion | elbow angle and speed, peak drop, drop at 200 ms, catch time t_c, passive-stop contact |
| EMG / IMU | biceps and triceps EMG and IMU at the repo defaults (2000 / 150 Hz, unchanged), pre-onset EMG onset flag |
| Timing | clock-sync error, sample-drop count, U_t in force |
| Outcome | DROP / NO DROP / INVALID, reason, who decided and when, blinding flag, adverse event, comfort report |

The 150 Hz IMU has one sample every 6.7 ms and cannot score a 0.1 ms criterion. The score uses the rig's own angle and torque channels at 10 kHz or faster (0.1 ms per sample; Rev U row 14 asks for a tenth of the window, 0.18 ms SYNTHETIC, rounded down). Raw data are kept. Nothing is deleted.

## 8. Assumptions (unsourced unless stated)

| # | ASSUMPTION |
|---|---|
| A1 | The 10° drop limit (B9). No human-safety source |
| A2 | The 200 ms catch. No human data |
| A3 | Drop model: instant catch, zero elbow damping, HE's arm and device parameters, k = 1.5862 |
| A4 | Zero hold torque represents a relaxed user, and one window fits every subject (Charlie 1a) |
| A5 | ICC of drops (0.5 for planning) and that it is the same for all subjects |
| A6 | 10 fault trials per subject (Rev U row 14) |
| A7 | Subjects are independent and represent the users; subject population not set |
| A8 | Testing only the governing cell is enough; catch time does not depend on the path |
| A9 | Score combines catch time and angle, with uncertainty added against a pass |
| A10 | Passive stop at 15°, model agreement 0.5 ms, gravity-moment error ±5 % (HE scoping) |
| A11 | Latch release and shams give no usable cue; one sham per fault trial |
| A12 | "Relaxed" means the EMG onset detector at unchanged defaults stays silent |
| A13 | (SYNTHETIC) FC's 20 %–33 %: linear extrapolation of [PRU09], normal catch times, SD 16 ms, ± read as SD |
| A14 | [PRU09] transfers to this task (step torque, horizontal plane, practised subjects there; ramp, vertical, unprepared here) |

## 9. Open items

1. Charlie: run Phase B at all, given section 6?
2. Charlie: 33 × 10 or 59 × one? ICC value? Which paths need their own sample?
3. FC / Charlie: confirm the combined score (A9).
4. HE: allowed ramp and STO overshoot on the rig; latch design; absorber proof; Rev U row 14 still says the window is recomputed per hold torque.
5. IRB: criteria, consent, adverse events, rest times, sham ratio.

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
