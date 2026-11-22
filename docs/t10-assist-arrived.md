# T10: when "assist arrived" (decision A)

T10's pass criterion: **assist arrived = torque t50 on the ensemble mean** (at least 10 steps, preferably 20), limit still at most 100 ms. t10 and t90 are reported alongside it but are not the criterion.

- Reference implementation: `scripts/t10_protocol.py` (`assist_arrived()`). Running `python scripts/t10_protocol.py docs/t10` regenerates every file in `docs/t10/`. All values there are SYNTHETIC.
- Source of the definition: the RevC design package (`EMG_Exo_Package_RevC.md`, line 203, plus its change log), which lives outside this repo.
- Status: spec only. No harness or production code changes. Protected values are unchanged.
