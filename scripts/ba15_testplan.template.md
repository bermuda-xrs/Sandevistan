# B-A15 arm-drop test plan (written plan, no trials run)

**Status: draft for Fact Checker. No human or bench trial has been run. No human trial may start before ethics (IRB) approval (section 2).**

- Every number comes from `scripts/ba15_numcheck.py` (output: `docs/ba15_numbers.json`). This file is rendered from `scripts/ba15_testplan.template.md`; do not edit it by hand.
- Numbers marked **SYNTHETIC** come from the drop model or from an extrapolation through it. They are not bench or human data. Sample-size numbers are exact binomial calculations, not simulations.
- **ASSUMPTION** marks anything with no source. The full list is in section 8.
- The window edges and FC's drop range are SYNTHETIC wherever they appear.
- "FC" is Fact Checker, "HE" is Hardware Engineer, "STO" is safe torque off. Citations are in the References.

## Short version

1. **Window (Charlie, 1a):** one fixed relaxed-hold catch-time window for every subject, {{in.win_lo|.2f}}–{{in.win_hi|.2f}} ms (SYNTHETIC). No per-subject recompute.
2. **Pass rule (Charlie, 2b):** show P(drop past {{in.drop|.0f}}°) ≤ {{in.p0|pct}} at {{in.conf|pct}} one-sided confidence with zero drops in {{in.n0}} independent trials.
3. **Design (Charlie, 2026-11-26 01:25 KST): {{in.n0}} subjects × one fault trial each.** Its false-pass chance is {{stat.chosen.false_pass_max|pct}} at any ICC (section 4). The earlier {{stat.plan.subjects}} × {{in.m_alt}} default is in "Rejected alternatives" (section 4b).
4. **Any drop:** stop, investigate, and the test has FAILED. No extension unless it was registered before trial 1.
5. **Dummy phase first:** {{syn.hw.matrix.n_cells}} cells × {{syn.hw.matrix.n_cell}} drops = {{syn.hw.matrix.drops|,}} drops, about {{syn.hw.matrix.hours|.1f}} h (ASSUMPTION: {{in.reset_min|g}} min per reset). Timing uncertainty is a requirement, at most {{syn.hw.budget.u_t_ms|.3f}} ms (SYNTHETIC budget, section 3), until it is measured on the rig.
6. **OPEN SAFETY ITEM THAT BLOCKS PHASE B:** there is no sourced human tolerance limit for an arm hitting the {{in.stop|.0f}}° stop (section 2). Phase B cannot start until there is one.
7. **Feasibility: this plan will almost certainly fail on the current design.** If the true drop probability is FC's {{in.fc_lo|pct}}–{{in.fc_hi|pct}}, the chosen design passes with chance {{syn.feas.chosen_pass.FC low|.1e}} to {{syn.feas.chosen_pass.FC high|.1e}} (SYNTHETIC). Section 6 gives the bound over the other options and says what it means.

## 1. What is tested

B-A15 is the project's human catch model: the user arrests the arm at exactly {{in.catch|.0f}} ms after fault onset. It has no human data behind it (ASSUMPTION). The {{in.drop|.0f}}° drop limit (B9) has no human-safety source either (ASSUMPTION): the S7 search found none and FC verified that (`S7_verdict.md` row 14).

**The fixed window.** Arm horizontal (elbow 90°), user hold torque zero ("relaxed"). The drop model (torque basis: ramp on the device side, times k = {{in.k|.4f}} at the arm) reaches {{in.drop|.0f}}° at:

| Edge (SYNTHETIC) | Converged value | Decided value | Governing cell |
|---|---|---|---|
| Lower: design (b), SYNTHETIC | {{syn.lo|.4f}} ms | {{in.win_lo|.2f}} ms | e-stop, ramp {{in.rate|.2f}} N·m/s, STO {{in.sto|.0f}} ms |
| Upper: backup (iv), SYNTHETIC | {{syn.hi|.4f}} ms | {{in.win_hi|.2f}} ms | comm loss, ramp {{in.rate_iv|.2f}} N·m/s, STO {{in.sto_comm|.0f}} ms |

- Both decided values (SYNTHETIC) reproduce at two decimals. The lower one is rounded down by {{ref.lo_round_down_ms|.4f}} ms (safe side); the upper one is rounded up by {{ref.hi_round_up_ms|.4f}} ms (not the safe side, but far below the {{in.res|g}} ms timing resolution).
- The lower edge is the earliest crossing over all {{syn.n_cells_b}} (b) cells (path × rate tolerance × STO band). The next paths cross later: comm loss {{syn.lo_comm|.2f}} ms, EMG loss {{syn.lo_emg|.2f}} ms (SYNTHETIC).
- Margin over the {{in.catch|.0f}} ms catch: {{syn.margin_ms|.2f}} ms. Drop at {{in.catch|.0f}} ms: {{syn.drop200|.2f}}°. Arm speed at {{in.drop|.0f}}°: {{syn.w10|.0f}} °/s (all SYNTHETIC).
- This decision replaces, for this test only, the Rev U row 14 rule that recomputes the window from each measured hold torque. One reason it cannot be recomputed well: a ±{{in.g_unc|pct}} error in a subject's gravity moment (ASSUMPTION) is {{syn.hold_unc|.2f}} N·m, {{syn.hold_unc_ratio|.2f}} times the whole nominal hold of {{syn.hold_nom|.3f}} N·m (SYNTHETIC).
- **Limit of one window.** A subject whose gravity moment is one percent above the model crosses {{in.drop|.0f}}° {{syn.dG_1pct|+.2f}} ms from the edge (minus is earlier); one percent more inertia moves it {{syn.dJ_1pct|+.2f}} ms (SYNTHETIC). The window is only {{syn.width|.2f}} ms wide (SYNTHETIC). So the trial score below also uses the directly measured angle.

## 2. Ethics and safety before any trial

- **Human testing needs ethics (IRB) approval before any trial**, including pilots and familiarisation. The Declaration of Helsinki, paragraph 23, requires the protocol to be approved by a research ethics committee before the research begins [WMA24]. Where US rules apply, the IRB approves, requires changes to, or disapproves the research [CFR46]. Which law applies here is not known (ASSUMPTION: whichever applies, approval comes first).
- The IRB must see section 6. Exposing people to repeated arm drops in a test that is expected to fail is a risk with little expected benefit.
- Inclusion and exclusion criteria, consent text, medical cover and the adverse-event definition are for the IRB protocol. This plan does not set them (open item).
- **Rig safety (HE scoping; ASSUMPTIONS until built and proven on the dummy).** A person may not use the rig until every item is proven:
  1. **Rig power cut, independent of the device e-stop.** One cut removes power from the rig's latch, hold emulator and drives, and works with the device e-stop bypassed or failed.
  2. **Interlocked guard** over the {{syn.hw.rig.r_guard|.3f}} m sweep (payload tip {{syn.hw.rig.r_tip|.2f}} m plus {{in.clr_mm|.0f}} mm clearance). The tip sweeps {{syn.hw.rig.sweep_mm|.0f}} mm of arc between {{syn.hw.rig.ang_bottom|.0f}}° and {{syn.hw.rig.ang_top|.0f}}° (SYNTHETIC geometry). Opening the guard cuts rig power.
  3. **Passive lower stop** at {{in.stop|.0f}}° drop with a {{in.th_stroke|.0f}}° absorber stroke, sized in the table below. In the Phase B test cell an uncaught arm reaches it {{syn.hw.rig.gap_gov_ms|.2f}} ms after its {{in.drop|.0f}}° crossing, at {{syn.hw.rig.t_stop_gov|.2f}} ms. The earliest contact of any cell is {{syn.hw.rig.t_stop_min|.2f}} ms, in (b) comm loss {{in.rate|.2f}} N·m/s, STO {{in.sto_comm|.0f}} ms, {{syn.hw.rig.gap_first15_ms|.2f}} ms after that cell's own {{in.drop|.0f}}° crossing (SYNTHETIC).
  4. **Upper runaway stop** at +{{in.th_up|.0f}}° above the start, for a device that drives upward.
  5. **Torque trip** at {{syn.hw.rig.trip_Nm|.1f}} N·m at the arm (SYNTHETIC k·T0 basis, ASSUMPTION).
  6. **Stop holder (requirement).** One person, named in the IRB protocol, whose only task during a trial is to hold the rig power cut and use it. Not the operator running the trial.
  7. **Subject stop button (requirement).** A stop the subject holds in the free hand, which cuts rig power. "Subject asks to stop" must not depend on speech alone.

{{table:absorber}}

  The absorber is sized by the larger case. The no-hold free fall sets the lowest case; with the device driving at its cap the stroke must absorb {{syn.hw.rig.E_cap_total|.2f}} J, against {{syn.hw.rig.E_ff_total|.2f}} J for free fall (SYNTHETIC). The stop is hit at {{syn.hw.rig.w_cap_stop|.0f}} °/s in the larger case, not the {{syn.hw.rig.ff_stop_w|.0f}} °/s of free fall (SYNTHETIC).

- **OPEN SAFETY ITEM THAT BLOCKS PHASE B: no human tolerance limit for hitting the {{in.stop|.0f}}° stop.**
  - The absorber above is sized for the rig, not for a person.
  - In Phase B an uncaught arm reaches the stop at {{syn.hw.rig.w_stop_gov|.0f}} °/s in the test cell (SYNTHETIC). In the runaway case it is {{syn.hw.rig.w_cap_stop|.0f}} °/s, with {{syn.hw.rig.T_abs_peak|.1f}} N·m and {{syn.hw.rig.F_abs_peak|.0f}} N at the contact (SYNTHETIC).
  - With the S7 catch model (normal, SD {{syn.fc_sd|.0f}} ms), the share of Phase B trials whose catch comes after the test cell's stop contact is {{syn.feas.contact.nominal hold.p_contact_gov|pct}} (nominal-hold mean) to {{syn.feas.contact.zero hold.p_contact_gov|pct}} (zero-hold mean) (SYNTHETIC). FC's {{syn.feas.contact.nominal hold.p_contact_first15|pct}}–{{syn.feas.contact.zero hold.p_contact_first15|pct}} uses the earliest contact of any cell (SYNTHETIC).
  - The stop rules halt at the first DROP, so at most one person hits the stop. But that person still takes it: the chance that the study ends on a stop contact is {{syn.feas.contact.nominal hold.p_study_contact|pct}} to {{syn.feas.contact.zero hold.p_study_contact|pct}} (SYNTHETIC).
  - **No sourced limit exists in this plan.** FC notes that ISO/TS 15066 has body-region limits but holds no clause text (`PR9_verdict.md` row 12); the S7 search could not read it either. No other source was found. **No limit is invented here.**
  - **Phase B is blocked until** a sourced human tolerance for this contact (speed, force or energy at the forearm) is in the plan, and the stop is shown on the dummy to stay within it.

## 3. Procedure

**Phase A: dummy arm, no people.**
1. **Build the rig** with a dummy forearm and the device fitted, start latch, hold emulator, angle and torque sensors, and every item of section 2. The dummy plus device must reproduce the model's inertia J = {{syn.J|.4f}} kg·m² and gravity moment G = {{syn.G|.3f}} N·m (SYNTHETIC) within the calibration limits below.
2. **Prove the safety items** of section 2 before any drop with the device on, including the rig power cut with the device e-stop bypassed.
3. **Prove the timing budget.** Every term below is measured on the rig, not taken from a data sheet:

{{table:budget}}

   Every term is a **requirement** until it is measured on the rig. If a measured term is over its share, Phase A stops until it is fixed. In Phase B, a trial whose measured release skew is over {{syn.hw.budget.release_ms|.3f}} ms is INVALID (SYNTHETIC limit). A late release holds the arm while the torque ramps down, which makes the trial easier. A drop in such a trial still triggers the stop rules. The torque and force channels have no bandwidth requirement yet (open item).
   The sampling period of {{syn.hw.budget.period_ms|.4f}} ms is also within the {{in.res|g}} ms timing resolution that Rev U row 14 asks for (a tenth of the {{syn.width|.2f}} ms window, rounded down; SYNTHETIC).
4. **Calibrate** the dummy. The flip values show how little the model can be off before a failing unit would read as passing; the calibration limits come from splitting the {{in.tol_model|g}} ms model agreement (ASSUMPTION) equally over four inputs ({{in.share|pct}} each):

{{table:flip}}

   The limits for G and for net pre-onset torque are one error, counted twice (conservative): {{syn.hw.cal.G_pct|.3f}} % of G is {{syn.hw.cal.G_Nm|.4f}} N·m (SYNTHETIC).
5. **Firmware test mode.** The governing cell needs ramp {{in.rate|.2f}} N·m/s and STO {{in.sto|.0f}} ms, the edge of the tolerance band, not the product settings ({{in.rate_nom|.2f}} N·m/s, {{in.sto_nom|.0f}} ms). That needs a test mode in the device firmware. **ASSUMPTION and open risk: it must be verified** that such a mode can set these values without changing the fault path, and the ramp and STO actually delivered must be measured and logged on every trial. A mode that changes the fault path makes the whole test invalid.
6. **Dummy matrix.** Every cell is run {{syn.hw.matrix.n_cell}} times: {{syn.hw.matrix.n_cells}} cells × {{syn.hw.matrix.n_cell}} = {{syn.hw.matrix.drops|,}} drops, about {{syn.hw.matrix.hours|.1f}} h at {{in.reset_min|g}} min per reset (ASSUMPTION). A cell is one of {{syn.hw.matrix.n_b}} (b) cells and {{syn.hw.matrix.n_iv}} (iv) cells (path × ramp −5 %, nominal, +5 % × STO band edges and nominal), at three holds: none, half and the full nominal hold ({{syn.hw.matrix.holds.0|.3f}}, {{syn.hw.matrix.holds.1|.3f}}, {{syn.hw.matrix.holds.2|.3f}} N·m; SYNTHETIC). That covers every D-4 corner, and the upper-edge cell of the window: (iv), comm loss, {{in.rate_iv|.2f}} N·m/s, STO {{in.sto_comm|.0f}} ms, crossing {{syn.hw.matrix.upper_edge_t10|.2f}} ms at hold 0 (SYNTHETIC). Model crossing times:

{{table:matrix}}

   {{syn.hw.matrix.dup_cells}} of the {{syn.hw.matrix.n_cells}} cells ({{syn.hw.matrix.dup_types}} types × 3 holds) repeat a setting already in the matrix, because the (b) −5 % rate equals the (iv) nominal rate: e-stop and comm loss at {{in.rate_dup|.2f}} N·m/s, all three STO values (SYNTHETIC count). They are kept as replication: {{syn.hw.matrix.dup_drops}} drops, {{syn.hw.matrix.dup_hours|.1f}} h. Dropping them leaves {{syn.hw.matrix.distinct}} distinct cells, {{syn.hw.matrix.distinct_drops}} drops and {{syn.hw.matrix.distinct_hours|.1f}} h. HE decides.
   The {{syn.hw.matrix.n_cell}} per cell is the Rev U row 14 floor; a trial SD of {{in.s_dummy|g}} ms and a wanted precision of {{in.e_dummy|g}} ms on the mean ask for {{syn.hw.matrix.n_formula|.2f}} drops (both ASSUMPTIONS; SYNTHETIC formula).
   **Per-trial pass:** crossing at or after {{syn.hw.budget.t10_min_pass_ms|.2f}} ms ({{in.catch|.0f}} ms plus U_t), **or** drop at {{in.catch|.0f}} ms at most {{syn.hw.budget.drop200_max_pass_deg|.3f}}° ({{in.drop|.0f}}° minus the angle uncertainty) (SYNTHETIC thresholds). Either reading may be used; a trial where only one passes is reported. **Per cell:** the mean crossing agrees with the model within ±{{in.tol_model|g}} ms (ASSUMPTION). Every model cell passes the per-trial rule with room: the worst model drop at {{in.catch|.0f}} ms is {{syn.hw.budget.worst_drop200_deg|.3f}}° (SYNTHETIC).
7. Any cell that fails stops Phase A. If the governing cell's crossing is outside the model by more than the agreement, the fixed window is wrong: stop and go to Charlie (via Seren).

**Gate to Phase B:** Phase A passed; **the open safety item of section 2 closed** (a sourced human tolerance for stop contact, met on the dummy); stop holder named and subject stop button fitted; IRB approval in hand; allowed drops registered (default zero); and Charlie's explicit go after reading section 6.

**Phase B: people.**
1. One test cell only: e-stop fault, ramp {{in.rate|.2f}} N·m/s, STO {{in.sto|.0f}} ms, relaxed hold. It has the earliest model crossing, so a given catch time is hardest to pass there (ASSUMPTION: a subject's catch time does not depend on which path failed). **A pass covers only the cell tested.** To claim each path separately, each needs its own full sample.
2. Seat the subject, elbow at 90°, payload in hand, forearm resting on the start latch. A relaxed arm would sag before the fault, because the device carries slightly less than the arm's weight, so the latch must hold it and release at fault onset (ASSUMPTION: release adds no usable cue).
3. Tell the subject to relax and to stop the arm only if it falls. Use no countdown and mask sound.
4. Random wait of 1 to 4 s before onset, as in [PRU09]. Mix in sham releases (latch opens, no fault) at random, one sham per fault trial (ASSUMPTION). A relaxed arm still sags slowly in a sham, so shams are never scored.
5. Fault onset is time zero. The catch time is from onset to arrest complete (arm speed zero), as in Rev U row 14.
6. **One fault trial per subject** (Charlie's design), with the shams mixed in before and after it. Then the subject's session ends.
7. No practice with the real fault. The unprepared user is the hazard, and [PRU09] subjects were practised. With one fault trial per subject, every fault is a first exposure.

**Trial score (fixed before trial 1).** Let t_c be the measured catch time and U_t the measured timing uncertainty.

| Score | Rule |
|---|---|
| **DROP** | t_c + U_t > {{in.win_lo|.2f}} ms (SYNTHETIC edge), where U_t is the measured timing uncertainty, required to be at most {{syn.hw.budget.u_t_ms|.3f}} ms (SYNTHETIC budget requirement), **or** measured peak drop + angle uncertainty ({{syn.hw.budget.ang_unc_deg|.4f}}°) > {{in.drop|.0f}}°, **or** the passive stop is touched |
| **NO DROP** | valid trial and none of the above |
| **INVALID** | rig fault (ramp slower than {{in.rate|.2f}} N·m/s, STO later than {{in.sto|.0f}} ms, release skew over {{syn.hw.budget.release_ms|.3f}} ms (SYNTHETIC limit), lost sensor or clock), or EMG activity before onset showing the subject was not relaxed |

- Using both time and angle makes the score stricter than either alone, so a pass holds for both readings of "drop past {{in.drop|.0f}}°" (ASSUMPTION: this is how 1a and 2b combine; for FC and Charlie to confirm).
- Not relaxed means the repo EMG onset detector, at its unchanged defaults, fires during the wait (ASSUMPTION).
- INVALID is decided from pre-onset and rig channels only, by someone who has not seen the angle trace after onset. Invalid trials are repeated and all are reported. **A drop in an invalid trial still triggers the stop rules.**
- The rig may not be easier than the cell. Going past it by one percent in ramp rate moves the crossing {{syn.drate_1pct|+.2f}} ms, and STO one ms early moves it {{syn.dsto_1ms|+.2f}} ms (SYNTHETIC). The allowed overshoot is open.

## 4. Sample size

**Independent trials.** Zero drops in n shows P ≤ {{in.p0|pct}} at {{in.conf|pct}} one-sided confidence when the chance of zero drops at exactly {{in.p0|pct}} is at most {{in.p0|pct}} (the exact one-sided Clopper–Pearson bound [CP34]; zero-event case [HL83]). {{in.n0}} is the smallest such n: at n = {{in.n0}} the chance is {{stat.q59|.5f}}; one trial fewer gives {{stat.q58|.5f}}. **Verified.**

**Chosen design (Charlie, 2026-11-26 01:25 KST): {{in.n0}} subjects × one fault trial each, {{in.n0}} fault trials in total.**

**Repeated trials in one subject are correlated.** The design effect is DE = 1 + (m − 1)·ICC, for m trials per subject. Effective sample size = trials / DE [DBB81, RS92]. We checked the formula on the published example in [KMP04] (DE {{stat.killip.deff|.3f}}, effective size {{stat.killip.ess|.1f}}; the paper gives the same DE and a rounded size). Subjects = ceil(n·DE/m). Entries are subjects (total fault trials):

{{table:subjects}}

**Why {{in.n0}} × one is robust to the ICC.** The m = 1 row is {{in.n0}} at every ICC. With one fault trial per subject, subjects' drop probabilities may differ in any way: the chance that one subject does not drop is still one minus the mean drop probability. So zero drops in {{in.n0}} subjects at a true {{in.p0|pct}} has chance {{stat.chosen.false_pass_max|pct}} at every ICC in the table, under both the beta-binomial and the worst-case mix (computed: {{stat.chosen.false_pass_min|.5f}} to {{stat.chosen.false_pass_max|.5f}}). No ICC has to be assumed, estimated or registered. The other rows show what the ICC would cost with more trials per subject.

- **The earlier estimate of {{stat.range_he.0}}–{{stat.range_he.2}} subjects is reproduced**, for m = {{in.m_alt}} and ICC {{stat.range_icc.0|g}}, {{stat.range_icc.1|g}} and {{stat.range_icc.2|g}}: {{stat.range_he.0}}, {{stat.range_he.1}}, {{stat.range_he.2}}. **It is not an upper bound.** It stops at ICC {{stat.range_icc.2|g}} and assumes m = {{in.m_alt}}. Over the full ICC range 0 to 1 it becomes {{stat.range_he.0}} to {{stat.plan.subjects_any}}. No source gives the ICC of arm drops (ASSUMPTION).

## 4b. Rejected alternatives

**{{stat.plan.subjects}} subjects × {{in.m_alt}} fault trials (the earlier planning default, DE {{stat.plan.deff|g}}, {{stat.plan.trials}} fault trials). Rejected by Charlie, 2026-11-26.** Reasons:
1. **It is only valid if the drop ICC is at most {{in.icc_alt|g}}, and that ICC must be fixed before trial 1.** A passing data set cannot tell you its ICC: with zero drops there is no variation to estimate it from. A pilot could measure the ICC of catch time. If catch times are normal (ASSUMPTION), a catch-time ICC of {{stat.latent_icc.0|g}}, {{stat.latent_icc.1|g}}, {{stat.latent_icc.2|g}} or {{stat.latent_icc.3|g}} gives a drop ICC of {{stat.latent.0|.2f}}, {{stat.latent.1|.2f}}, {{stat.latent.2|.2f}} or {{stat.latent.3|.2f}} at {{in.p0|pct}}.
2. **Above an ICC of {{stat.alt_icc_break|.4f}}, its {{stat.plan.trials}} trials are worth fewer than {{in.n0}} independent ones, and the false-pass chance at a true {{in.p0|pct}} goes over {{in.p0|pct}}:**

{{table:alt_risk}}

3. **A high ICC is plausible here.** By the mapping in reason 1, a catch-time ICC of {{stat.latent_icc.2|g}} to {{stat.latent_icc.3|g}} already gives a drop ICC of {{stat.latent.2|.2f}} to {{stat.latent.3|.2f}}.
4. **Too small an ICC breaks it the other way.** Sizing for ICC 0 ({{stat.range_he.0}} subjects) when the truth is {{in.icc_alt|g}} gives a false-pass chance of {{stat.misspec.0|pct}} (beta-binomial) to {{stat.misspec.1|pct}} (worst case), not {{in.p0|pct}}.
5. Repeated faults teach the subject what is coming; the unprepared response is the hazard (section 3).

For the record, the design-effect check for m = {{in.m_alt}}:

**Is the design effect safe for a zero-drop test?** It is a variance correction, so we checked it exactly at P = {{in.p0|pct}}, m = {{in.m_alt}}, against (a) a beta-binomial model [W75] and (b) the worst case over every possible mix of subject drop rates with that mean and ICC. The worst case (some subjects at one rate, the rest never dropping) is confirmed by a linear program (they agree within {{in.lp_tol|.0e}}):

{{table:check}}

- The design effect is enough in every row. On a grid of ICC 0 to 1 in steps of one hundredth it stays enough for every m up to {{stat.deff_ok_to_m}}, and first falls short at the next m. **Do not use more than {{stat.deff_ok_to_m}} trials per subject with this formula.**

**The other trial counts per subject** in the table above carry the same ICC dependence and are rejected for the same reasons.

## 5. Pass, fail and stop rules

**PASS:** every planned valid trial done, zero DROP scores.
**FAIL:** any DROP score.

**If a drop occurs**
1. **Stop** the session and all enrolment at once. Check on the subject. Report as the IRB requires.
2. **Investigate** before anything else: rig, ramp, STO, hold, catch time, angle.
3. **The zero-drop test has FAILED.** A drop is never removed afterwards. Report the exact one-sided upper bound: one drop in {{in.n0}} independent trials gives P ≤ {{stat.cp.1|pct}}; two give {{stat.cp.2|pct}}.
4. **No extending after the fact.** "Add trials and allow one drop" chosen after seeing the drop has a false-pass chance of {{stat.naive1|pct}} ({{stat.n0}} then {{stat.n1}}) or {{stat.naive2|pct}} ({{stat.n0}}, {{stat.n1}}, then {{stat.n2}}), not {{in.p0|pct}}.
5. **Where the drop falls** (from t_c; SYNTHETIC edges; decides who is told what, not a pass): above {{in.win_lo|.2f}} and up to {{in.win_hi|.2f}} ms, (b) has failed and (iv) becomes the candidate. (iv) needs Charlie's switch decision and its own new full test in its own cell. Above {{in.win_hi|.2f}} ms both fail: escalate to Charlie (via Seren).
6. **Retest** only after a documented design change, from zero. Trying again without a change until it passes gives a false-pass chance of {{stat.retest.0|pct}} over two tries and {{stat.retest.1|pct}} over three.

**Allowed-drop plans** (must replace the zero-drop plan before trial 1, never after a drop). With one trial per subject, subjects = trials at any ICC:

{{table:allowed}}

**Registered extension that keeps the first look at {{in.n0}}** (also before trial 1 only). The overall false-pass chance is {{stat.ext_size|.7f}}; with one trial per subject it does not depend on the ICC. It is costly because the first look already spends {{stat.q59|.4f}}:

{{table:ext}}

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

**In plain words: if the true drop probability is {{in.fc_lo|pct}}–{{in.fc_hi|pct}}, the chance of passing is effectively zero (SYNTHETIC).**
- **Chosen design ({{in.n0}} × one):** {{syn.feas.chosen_pass.FC low|.1e}} at {{in.fc_lo|pct}} and {{syn.feas.chosen_pass.FC high|.1e}} at {{in.fc_hi|pct}} (SYNTHETIC), at any ICC.
- **Every option in this plan's tables** (trials per subject, ICC, allowed drops and the registered extension): at most {{syn.feas.max_pass|.1e}}, about one in {{syn.feas.one_in|,.0f}} (SYNTHETIC). This bound holds **only for the tabulated options**.
- **Any trial count per subject up to {{stat.deff_ok_to_m}}** (where the design-effect check holds), any ICC in steps of 1/{{in.sup_steps}}, zero drops, best-case mix: at most {{syn.feas.sup_any_m.lo.P|.1e}} at {{in.fc_lo|pct}} and {{syn.feas.sup_any_m.hi.P|.1e}} at {{in.fc_hi|pct}}, both at m = {{syn.feas.sup_any_m.lo.m}}, ICC {{syn.feas.sup_any_m.lo.icc|g}}, {{syn.feas.sup_any_m.lo.subjects}} subjects (SYNTHETIC). All of those options are rejected (section 4b).

Chosen design and the allowed-drop sizes (independent trials, one per subject):

{{table:feas}}

Rejected clustered designs, subjects inflated for the ICC (m = {{in.m_alt}}). "Best case" is the mix of subjects most favourable to passing:

{{table:feas_icc}}

- **More subjects do not help.** No row gives any real chance of passing. Under-sizing raises the chance only by breaking the plan: sized for ICC 0 with a true ICC of {{in.icc_alt|g}}, a design that drops {{in.fc_lo|pct}} of the time would pass with chance {{syn.feas.undersized|pct}} under the beta-binomial and {{syn.feas.undersized_best|pct}} under the best-case mix (SYNTHETIC).
- **The first drop is expected almost at once:** on average at subject {{syn.feas.first_drop.FC high|.1f}} to {{syn.feas.first_drop.FC low|.1f}} (SYNTHETIC).
- **FC's {{in.fc_lo|pct}}–{{in.fc_hi|pct}} is for the nominal hold, not this window.** We reproduce it: {{syn.fc_p_lo|pct}}–{{syn.fc_p_hi|pct}} at mean {{syn.fc_mean|.1f}} ms, SD {{syn.fc_sd|.0f}} ms (SYNTHETIC; a linear extrapolation of [PRU09] outside its data). With the **relaxed** window Charlie chose and the same catch times it is {{syn.p_hi_same_mean|pct}} at the (iv) edge and {{syn.p_lo_same_mean|pct}} at the (b) edge (SYNTHETIC). Extending the trend to zero hold (mean {{syn.mean_h0|.1f}} ms) gives {{syn.p_hi_h0_mean|pct}} and {{syn.p_lo_h0_mean|pct}} (SYNTHETIC). So FC's range is the optimistic case here.
- **No [PRU09] condition passes either.** Reading ± as SD and times as normal (ASSUMPTIONS), the best is {{syn.pru_best_p|pct}} drops and a {{syn.pru_best_pass|pct}} chance of passing (SYNTHETIC):

{{table:pru}}

- **What would pass.** Chance of passing by true drop probability (exact, independent trials):

{{table:oc}}

  A passing chance of four in five needs a true drop probability of {{stat.p80|pct}} or less; an even chance needs {{stat.p50|pct}}. **A design exactly at the {{in.p0|pct}} target would almost always fail.** With an SD of {{in.pru_sd|.0f}} ms (SYNTHETIC use of [PRU09]), a true {{in.p0|pct}} at the (b) edge needs a mean catch of {{syn.mean_for_p0|.1f}} ms, and a four-in-five chance needs {{syn.mean_for_80|.1f}} ms (SYNTHETIC). At FC's mean catch the {{in.drop|.0f}}° crossing would have to move to {{syn.cross_for_p0|.1f}} ms and {{syn.cross_for_80|.1f}} ms (SYNTHETIC).

**What it implies**
1. The test is a valid way to *reject* the design, and would do so quickly. It is not a realistic way to *accept* the current one.
2. Running it on people as planned exposes subjects to drops for a result already expected. **Recommendation: do Phase A only, and do not start Phase B on the current design.** Charlie decides.
3. FC's estimate is an extrapolation and could be wrong. A small IRB-approved pilot that measures catch time (not pass/fail) would test that at lower cost than Phase B. It is not part of this plan.
4. This plan proposes no design fix. Any fix must first be checked against the failures log.

## 7. Data logged per trial

| Group | Channels |
|---|---|
| Identity | coded subject, session, trial number, fault or sham, random wait, test cell, operator, firmware commit, plan version |
| Commands | ramp rate and STO time commanded **and measured as delivered** (firmware test mode), fault-onset trigger, latch-release time |
| Torque | commanded and measured joint torque (device side), arm-side torque at the cuff, k used, torque-trip state |
| Hold | latch load before release (arm side), payload mass, subject arm mass estimate |
| Motion | elbow angle and speed, **raw encoder counts**, peak drop, drop at {{in.catch|.0f}} ms, catch time t_c |
| Rig | **stop contact force and stroke**, **rig power-cut and guard-interlock state**, torque-trip state |
| EMG / IMU | biceps and triceps EMG and IMU at the repo defaults ({{in.emg_hz}} / {{in.imu_hz}} Hz, unchanged), pre-onset EMG onset flag |
| Timing | **fault and release timestamps on one clock, with the measured skew between them**, clock-sync error, sample-drop count, U_t in force |
| Dummy | **mass, inertia and gravity-moment calibration records** (before and after each day), the matrix cell and the model crossing it is compared with |
| Outcome | DROP / NO DROP / INVALID, reason, who decided and when, blinding flag, adverse event, comfort report |

All rig channels are logged at {{syn.hw.budget.fs_khz|.0f}} kHz or faster (the budget in section 3; one sample every {{syn.hw.budget.period_ms|.4f}} ms, SYNTHETIC). The {{in.imu_hz}} Hz IMU has one sample every {{in.imu_ms|.1f}} ms and cannot score the criterion; it is logged only as a cross-check at the repo defaults. Raw data are kept. Nothing is deleted.

## 8. Assumptions (unsourced unless stated)

| # | ASSUMPTION |
|---|---|
| A1 | The {{in.drop|.0f}}° drop limit (B9). No human-safety source |
| A2 | The {{in.catch|.0f}} ms catch. No human data |
| A3 | Drop model: instant catch, zero elbow damping, HE's arm and device parameters, k = {{in.k|.4f}} |
| A4 | Zero hold torque represents a relaxed user, and one window fits every subject (Charlie 1a) |
| A5 | ICC of drops: not needed by the chosen design; used only to compare designs ({{in.icc_alt|g}} was the rejected default) |
| A6 | Not an assumption: one fault trial per subject is Charlie's decision (2026-11-26). It is below the Rev U row 14 floor of {{in.m_alt}} trials per subject (open item) |
| A7 | Subjects are independent and represent the users; subject population not set |
| A8 | Testing only the governing cell is enough; catch time does not depend on the path |
| A9 | Score combines catch time and angle, with uncertainty added against a pass |
| A10 | Subject gravity moment known to ±{{in.g_unc|pct}} (HE scoping) |
| A11 | Latch release and shams give no usable cue; one sham per fault trial |
| A12 | "Relaxed" means the EMG onset detector at unchanged defaults stays silent |
| A13 | (SYNTHETIC) FC's {{in.fc_lo|pct}}–{{in.fc_hi|pct}}: linear extrapolation of [PRU09], normal catch times, SD {{in.pru_sd|.0f}} ms, ± read as SD |
| A14 | [PRU09] transfers to this task (step torque, horizontal plane, practised subjects there; ramp, vertical, unprepared here) |
| A15 | Timing budget: {{in.budget_pct|pct}} of the smallest model margin, split {{in.alloc.sample|pct}} sampling, {{in.alloc.trigger|pct}} trigger, {{in.alloc.release|pct}} release, {{in.alloc.angle|pct}} angle, {{in.alloc.enc_lat|pct}} encoder latency, {{in.alloc.spare|pct}} spare; sampling rounded up to {{in.step_khz}} kHz steps (HE, reusing the Rev U row-8 rule) |
| A16 | Model agreement of {{in.tol_model|g}} ms, split equally over four inputs (net torque, J, G, k) |
| A17 | Dummy trial SD {{in.s_dummy|g}} ms, wanted precision {{in.e_dummy|g}} ms, {{in.reset_min|g}} min per reset |
| A18 | Rig: stop at {{in.stop|.0f}}° with a {{in.th_stroke|.0f}}° stroke, upper stop at +{{in.th_up|.0f}}°, contact radius {{in.r_stopc|.2f}} m, clearance {{in.clr_mm|.0f}} mm, absorber peak {{in.abs_profile|g}} × the ideal force, torque trip {{in.trip_pct|pct}} over k·T0, device-at-cap runaway case (HE) |
| A19 | **Firmware test mode** that sets ramp {{in.rate|.2f}} N·m/s and STO {{in.sto|.0f}} ms without changing the fault path. Open risk: **must be verified** |
| A20 | Encoder bandwidth from the first-order rise-time rule, bandwidth = {{in.rise_k|g}} / rise time, with the encoder-latency share as the rise time (HE's method) |

## 9. Open items

1. **OPEN SAFETY ITEM, BLOCKS PHASE B:** a sourced human tolerance limit for stop contact (section 2), and proof on the dummy that the stop stays within it.
2. Charlie: run Phase B at all, given section 6?
3. Charlie: which paths need their own sample? (Design decided: {{in.n0}} × one, 2026-11-26.)
4. FC / Charlie: confirm the combined score (A9).
5. HE: allowed ramp and STO overshoot on the rig; latch design and whether the release skew can be met; absorber proof on the dummy; bandwidth for the torque and force channels. Rev U row 14 still says the window is recomputed per hold torque and asks for {{in.m_alt}} trials per subject; both differ from Charlie's decisions.
6. HE / firmware: **the test mode (A19) must be verified**: ramp {{in.rate|.2f}} N·m/s and STO {{in.sto|.0f}} ms delivered as set, fault path unchanged.
7. Before Phase B: name the stop holder (IRB protocol) and fit the subject stop button.
8. IRB: criteria, consent, adverse events, rest times, sham ratio.

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
