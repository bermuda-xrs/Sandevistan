# T10: when "assist arrived" (decision A)

T10's pass criterion: **assist arrived = torque t50 on the ensemble mean** (at least 10 steps, preferably 20), limit still at most 100 ms. t10 and t90 are reported alongside it but are not the criterion.

- Reference implementation: `scripts/t10_protocol.py` (`assist_arrived()`). Running `python scripts/t10_protocol.py docs/t10` regenerates every file in `docs/t10/`. All values there are SYNTHETIC.
- Source of the definition: the RevC design package (`EMG_Exo_Package_RevC.md`, line 203, plus its change log), which lives outside this repo.
- **What the 100 ms is:** RevC's engineering target (A18, line 43: "Engineering target. Must be verified against user-felt-lag testing"; risk #7, line 228). Torque t50 is not the quantity behind the often-cited 100 ms of Farrell & Weir 2007 (doi:10.1109/TNSRE.2007.891391): theirs is a pure delay from EMG to motor command, actuator excluded. The half-window mapping that links the two comes from Farrell & Weir 2008 and covers state-select controllers only, not RevC's proportional mode. See session 6 `REPORT_v2.md` (`/workspace/claude_sessions/6_latency-budget/REPORT_v2.md`, outside this repo): §(a) rows L1 and L2b (lines 36, 38), Verdict (a) item 1 and its "Caveat (FC N2)" (lines 46–48).
- Status: spec only. No harness or production code changes. Protected values are unchanged.

## Validity guard: a run that is not valid does not pass

`assist_arrived()` returns `status` = PASS, FAIL (t50 > 100 ms) or INVALID. Only PASS counts as passed. A run is INVALID if any of these fails on the ensemble-mean trace (`guard()` in the script):

| Check | Rule | Why |
|---|---|---|
| Step above noise, positive | step > `K_NOISE` = 10 × baseline SD | Without it, zero torque, sensor noise only, or a step the wrong way gave t50 = 0.00 ms and PASS (FC PR #13 item #8). With 10 SDs the t50 threshold is more than 5 SDs above baseline. For Gaussian noise, P(z > 5) = 2.867e-7 per sample, and the union bound over the 1200 samples after T0 is 3.44e-4 (`work/tail.py` → `results/tail_out.txt`; session-14 paths are under `/workspace/claude_sessions/14_t10-fix/`, outside this repo). **ASSUMPTION** (engineering choice, not sourced; assumes Gaussian baseline noise) |
| Step above an absolute floor | step > `MIN_STEP` = 0.5, in the torque's own units (N·m at the joint on the bench) | 0.5 N·m is RevC's T5 felt-torque limit for transparency (line 198). A step smaller than the residual the user may already feel is not an assist. **ASSUMPTION** (using T5 for this is our choice). The synthetic chain's units are arbitrary. Its smallest step is 10.45, far above the floor |
| Crossing after T0 | t50 > 0 ms | A causal chain cannot reach the threshold on the T0 sample itself |
| Settled by T0 + 350 ms | the two halves of the plateau window (T0 + 350–450 and 450–550 ms) differ by at most `SETTLE_TOL` = 15 % of the step | The plateau is the mean over T0 + 350–550 ms and the record ends at T0 + 600 ms, so a trace still moving there reads early (FC item #9: a ramp with 300 ms delay and 400 ms rise, true t50 500 ms, reads 375 ms; its drift is 0.668 of the step). **ASSUMPTION**: 15 % sits above the drift that MAV ripple alone gives on the synthetic chain (largest 0.0810 of 4000 runs; mean 0.0132–0.0193, SD 0.0103–0.0146 per row) |

On the 4000 synthetic runs behind `docs/t10/t10_criterion.*`, the guard flags 0 runs, and all 4 output files are byte-identical to before the guard. The smallest step / baseline SD is 55.7, against the 10 required (stdout of `t10_protocol.py`, saved as session 14 `results/t10_protocol_stdout.txt`). SYNTHETIC.

**Limit that remains (stated, not fixed).** t50 is taken against the plateau over T0 + 350–550 ms, not against the final value. A slow first-order (exponential) torque response is nearly flat in that window but still below its final value, so the drift check cannot see it. The sweep covered linear ramps (delay 0–500 ms, rise 5–800 ms) and first-order responses (delay 0–500 ms, τ 2–300 ms), with true t50 taken against the final value. Linear ramps gave no false pass at any tolerance from 5 % to 30 %. First-order responses gave 30 false passes at every tolerance in that range: true t50 100.11–108.13 ms read as 94.50–100.00 ms, at most 8.63 ms early, with τ 80–156 ms and drift at most 0.0395 of the step. That drift is below the synthetic chain's own ripple, so no tolerance that keeps normal runs valid can catch these. They are torque loops of roughly 1–2 Hz, far below RevC T7 (−3 dB ≥ 30 Hz). Treat a bench t50 within about 9 ms of the limit as marginal when the torque response looks first-order and slow. Source: session 14 `work/settle_sweep.py` → `results/settle_sweep_out.txt`. SYNTHETIC.

Tests: `tests/test_t10_guard.py` (zero, noise only, reversed, below the floor, buried in noise, crossing at T0, FC slow ramp, two-stage late response, normal steps, unchanged outputs).
