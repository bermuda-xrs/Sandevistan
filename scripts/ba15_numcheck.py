"""Numbers for the B-A15 arm-drop test plan (docs/ba15-testplan.md).

Every number in the plan is computed here. The plan is rendered from ba15_testplan.template.md, so a
number cannot be typed into it by hand: literals left in the template must be listed in LITERALS.

    python scripts/ba15_numcheck.py           check the committed plan and JSON; exit 1 on any mismatch
    python scripts/ba15_numcheck.py --write   rewrite docs/ba15-testplan.md and docs/ba15_numbers.json

ALL DROP-MODEL NUMBERS ARE SYNTHETIC (model output, not bench or human data); they live under "syn".
The sample-size numbers are exact binomial calculations. No random numbers are used anywhere.

Torque basis: T0 and the ramp rates are device side; G and the user's hold are arm side; arm torque
from the device = k * T(t). Nothing here imports or changes the detectors, the gate or PLANNED_SHAPE.
"""

import json
import math
import re
import sys
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, linprog
from scipy.stats import beta, betabinom, binom, norm

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "scripts" / "ba15_testplan.template.md"
PLAN = ROOT / "docs" / "ba15-testplan.md"
NUMBERS = ROOT / "docs" / "ba15_numbers.json"

# ---------------------------------------------------------------- inputs
# Charlie's decisions, 2026-11-23 (fixed).
WINDOW_MS = (201.90, 203.72)  # (1a) one relaxed-hold window for every subject
P0, CONF, N0 = 0.05, 0.95, 59  # (2b)

# Drop model (Hardware Engineer, exo/calc_BA15_testplan.py; exo/ is not in this repo).
ARM = dict(m_body=80.0, f_fa=0.016, f_hand=0.006, L_fa=0.27, com_fa=0.430, rg_fa=0.303, hand_off=0.06,
           m_pay=2.0, pay_off=0.08, m_dev=0.30, r_dev=0.12, m_carr=0.04, s_mid=0.175, J_rot=0.5e-5,
           gear=30, g=9.80665)
K_HI, T0_DEV = 1.5862286754583168, 6.304261267443454  # exo/benchC_rules_out.json: revD.k_hi, chain_nom
DROP_LIM, CATCH_MS, STOP_DEG = 10.0, 200.0, 15.0
T_DET = {"estop": 0.0, "comm": 3.0, "emg": 50.0}  # ms
STO = {"estop": (148.0, 150.0, 160.0), "comm": (142.0, 150.0, 160.0), "emg": (None,)}
# Logged data rate (ASSUMPTION, not HE's): 14 rig channels (FC's count from section 7), 4-byte
# samples, 10 s record.
DATA_CHANNELS, DATA_BYTES, DATA_TRIAL_S = 14, 4, 10
RATE = {"b": {"estop": 20.0, "comm": 20.0, "emg": 50.0}, "iv": {"estop": 19.0, "comm": 19.0}}
TOL = 0.05
G_UNC = 0.05  # subject gravity moment known to +-5 % (ASSUMPTION, HE)
RES_MS = 0.1  # RevU row 14 timing resolution
TOL_MODEL_MS = 0.5  # dummy agreement with the model (ASSUMPTION, HE)
EMG_HZ, IMU_HZ = 2000, 150  # repo defaults, quoted only (protected, not changed)
T_CAP_DEV = 7.898894154818326  # device chain cap, N*m (exo/benchC_rules_out.json chain_nom cap)

# Timing budget and rig (ASSUMPTIONS from HE's scoping, exo/calc_BA15_testplan.py, unless noted).
BUDGET_FRAC = 0.10  # total timing uncertainty <= 1/10 of the smallest margin (Rev U row-8 rule, reused)
ALLOC = dict(sample=0.40, trigger=0.10, release=0.20, angle=0.10, enc_lat=0.10, spare=0.10)
FS_STEP_KHZ = 5  # sampling rate rounded up to this step
RISE_K = 0.35  # first-order rise-time rule t_r * f_3dB (ASSUMPTION, HE's method)
INPUT_SHARE = 0.25  # each of 4 model inputs (net torque, J, G, k) may use this share of TOL_MODEL_MS
HOLD_FRACS = (0.0, 0.5, 1.0)  # dummy hold torque as a fraction of G - k*T0 (0 = relaxed)
# Dummy trial SD (ms), wanted precision of the mean (ms), Rev U row 14 floor of drops per cell, reset (min)
S_DUMMY_MS, E_DUMMY_MS, N_MIN_CELL, T_RESET_MIN = 0.15, 0.1, 10, 1.0
Z_TWO_SIDED = 0.975  # two-sided 95 % normal quantile level for the dummy precision
# th_stroke: absorber stroke beyond engagement; th_up: upper runaway stop above 90 deg (deg);
# r_stopc, r_tip_extra, clr: stop contact radius, payload half-length, guard clearance (m);
# abs_profile: real absorber peak / ideal torque; trip_frac: rig torque trip / (k*T0)
RIG = dict(th_stroke=10.0, th_up=15.0, r_stopc=0.20, r_tip_extra=0.05, clr=0.025, abs_profile=1.5,
           trip_frac=1.10)

# PRU09 (Pruszynski et al. 2009, Results, "Behavior and kinematics"), ms, read as mean and SD.
PRU = {"reversal, 1 N·m background": (193.0, 16.0), "reversal, 3 N·m background": (170.0, 22.0),
       "maximal excursion, 1.25 N·m step": (178.0, 22.0), "maximal excursion, 2.5 N·m step": (184.0, 15.0)}
FC_P = (0.20, 0.33)  # Fact Checker, S7_verdict.md row 12

# Hardware Engineer's and Fact Checker's shipped values, for comparison only.
REF = dict(window=(201.90484840282784, 203.7164677450151), n_allow=(59, 93, 124), subjects=(6, 17, 33),
           fc_p=(0.19864655721181523, 0.3257000686903164), h0_p_worst=0.4795370246007952)

M_CHOSEN = 1  # Charlie's decision 2026-11-26 01:25 KST: 59 subjects x 1 fault trial
# rejected alternative 33 x 10 (Rev U row 14 floor; top of the earlier 6-33 range)
M_ALT, ICC_ALT = 10, Fr("0.5")
# true ICCs at which the rejected design is checked
ALT_TRUE_ICCS = [Fr(s) for s in ("0.5", "0.6", "0.8", "1")]
SUP_ICC_STEPS, SUP_M_MAX = 2000, 46  # continuous-ICC sweep for the bound over every allowed m
ICCS = [Fr(s) for s in ("0", "0.05", "0.1", "0.2", "0.3", "0.5", "0.8", "1")]
MS = (1, 2, 5, 10, 20)
P_ROWS = (0.001, 0.01, 0.02, 0.05)
LP_TOL = 1e-9


# ---------------------------------------------------------------- drop model (SYNTHETIC)
def arm(mass_scale=1.0) -> dict:
    """Arm plus device inertia and gravity moment. ``mass_scale`` scales the human-side masses only."""
    a, ms = ARM, mass_scale
    parts = [(a["f_fa"] * a["m_body"] * ms, a["com_fa"] * a["L_fa"], a["rg_fa"] * a["L_fa"]),
             (a["f_hand"] * a["m_body"] * ms, a["L_fa"] + a["hand_off"], 0.0),
             (a["m_pay"] * ms, a["L_fa"] + a["pay_off"], 0.0)]
    J_h = sum(m * (r**2 + rg**2) for m, r, rg in parts)
    Mr_h = sum(m * r for m, r, _ in parts)
    J = J_h + a["m_dev"] * a["r_dev"] ** 2 + a["m_carr"] * a["s_mid"] ** 2 + a["J_rot"] * a["gear"] ** 2
    G = (Mr_h + a["m_dev"] * a["r_dev"] + a["m_carr"] * a["s_mid"]) * a["g"]
    return dict(J=J, G=G, G_human=Mr_h * a["g"], kT0=K_HI * T0_DEV, hold_nom=G - K_HI * T0_DEV)


def _integrate(path, rate, sto_ms, hold, J=None, G=None, k=None, T0=None, step_loss=False, ang=None,
               at_ms=None):
    """Piecewise DOP853 solve at every torque breakpoint.
    x = [angle from horizontal (negative = drop), speed].
    ``ang``: return (time ms, speed deg/s) when the drop first reaches ``ang`` deg. ``at_ms``: return
    (drop deg, speed deg/s) at that time. ``step_loss``: device torque lost at t = 0 (rig fault case)."""
    A = arm()
    J, G = A["J"] if J is None else J, A["G"] if G is None else G
    k, T0 = K_HI if k is None else k, T0_DEV if T0 is None else T0
    td, sto = T_DET[path] / 1e3, math.inf if sto_ms is None else sto_ms / 1e3
    end = 0.4 if at_ms is None else at_ms / 1e3

    def tq(t):
        return T0 if t < td else max(T0 - rate * (t - td), 0.0)

    edges = sorted({0.0, end} | {e for e in (td, td + T0 / rate, sto) if 0 < e < end})
    y = [0.0, 0.0]

    def hit(t, x):
        return x[0] + math.radians(ang)

    hit.terminal = True
    for a, b in zip(edges[:-1], edges[1:]):
        off = step_loss or a >= sto  # STO opens the torque at its edge; the ramp is smooth inside a segment

        def f(t, x, off=off):
            return [x[1], (k * (0.0 if off else tq(t)) + hold - G * math.cos(x[0])) / J]

        s = solve_ivp(f, (a, b), y, method="DOP853", rtol=1e-11, atol=1e-12, events=hit if ang else None)
        if ang and s.t_events[0].size:
            return float(s.t_events[0][0]) * 1e3, -math.degrees(s.y_events[0][0][1])
        y = s.y[:, -1]
    if at_ms is not None:
        return -math.degrees(y[0]), -math.degrees(y[1])
    raise ValueError("no crossing")


def crossing(path, rate, sto_ms, hold, ang=DROP_LIM, **kw):
    """Time (ms) and speed (deg/s) at which the drop from 90 deg reaches ``ang`` deg."""
    return _integrate(path, rate, sto_ms, hold, ang=ang, **kw)


def drop_at(t_ms, path, rate, sto_ms, hold, **kw):
    """Drop (deg) and speed (deg/s) at ``t_ms``."""
    return _integrate(path, rate, sto_ms, hold, at_ms=t_ms, **kw)


def cells(opt):
    for path, r in RATE[opt].items():
        for tol in (-TOL, 0.0, TOL):
            for s in STO[path]:
                yield path, r * (1 + tol), s


def pru_mean(hold):
    """FC's linear extrapolation of the PRU09 reversal time to ``hold`` N*m (outside the data)."""
    (t1, _), (t3, _) = list(PRU.values())[:2]
    return t1 + (t3 - t1) / 2 * (hold - 1)


def build_grid() -> dict:
    """All 117 dummy-matrix cells: option x path x rate x STO x hold fraction (SYNTHETIC model values)."""
    A, grid = arm(), {}
    for opt in RATE:
        for path, rate, sto in cells(opt):
            for hf in HOLD_FRACS:
                hold = hf * A["hold_nom"]
                t10, w10 = crossing(path, rate, sto, hold)
                t_stop, w_stop = crossing(path, rate, sto, hold, ang=STOP_DEG)
                d200, w200 = drop_at(CATCH_MS, path, rate, sto, hold)
                grid[(opt, path, rate, sto, hf)] = dict(t10=t10, w10=w10, t_stop=t_stop, w_stop=w_stop,
                                                        drop200=d200, w200=w200)
    return grid


def model(grid) -> dict:
    A = arm()
    b = {k[1:4]: v["t10"] for k, v in grid.items() if k[0] == "b" and k[4] == 0.0}
    iv = {k[1:4]: v["t10"] for k, v in grid.items() if k[0] == "iv" and k[4] == 0.0}
    gov, gov_iv = min(b, key=b.get), min(iv, key=iv.get)
    lo, hi = b[gov], iv[gov_iv]
    g0 = grid[("b",) + gov + (0.0,)]
    w10, t_stop, w_stop = g0["w10"], g0["t_stop"], g0["w_stop"]
    out = dict(A, lo=lo, hi=hi, width=hi - lo, gov="|".join(map(str, gov)), gov_key=list(gov),
               gov_iv="|".join(map(str, gov_iv)),
               lo_comm=min(v for (p, _, _), v in b.items() if p == "comm"),
               lo_emg=min(v for (p, _, _), v in b.items() if p == "emg"),
               hi_estop=min(v for (p, _, _), v in iv.items() if p == "estop"),
               n_cells_b=len(b), n_cells_iv=len(iv),
               margin_ms=lo - CATCH_MS, drop200=g0["drop200"], w10=w10,
               res_req=(hi - lo) / 10, ang_res=w10 * RES_MS / 1e3, t_stop=t_stop, w_stop=w_stop)
    # limits of one fixed window: other subjects cross at other times (decision 1a stays fixed)
    path, r, s = gov
    out["dG_1pct"] = crossing(path, r, s, 0.0, G=A["G"] * 1.01)[0] - lo
    out["dJ_1pct"] = crossing(path, r, s, 0.0, J=A["J"] * 1.01)[0] - lo
    out["drate_1pct"] = crossing(path, r * 1.01, s, 0.0)[0] - lo
    out["dsto_1ms"] = crossing(path, r, s - 1, 0.0)[0] - lo
    out["hold_unc"] = G_UNC * A["G_human"]
    out["hold_unc_ratio"] = out["hold_unc"] / A["hold_nom"]
    # Fact Checker's estimate: nominal hold, five cells, PRU09 trend, SD 16 ms
    sd = PRU["reversal, 1 N·m background"][1]
    fc_cells = [("estop", 20.0, 150.0), ("comm", 20.0, 150.0), ("emg", 50.0, None), ("estop", 21.0, 148.0),
                ("comm", 21.0, 142.0)]
    x = [crossing(*c, A["hold_nom"])[0] for c in fc_cells]
    mu = pru_mean(A["hold_nom"])
    out.update(fc_mean=mu, fc_sd=sd, fc_x=x, fc_p_lo=float(norm.sf(max(x), mu, sd)),
               fc_p_hi=float(norm.sf(min(x), mu, sd)), mean_h0=pru_mean(0.0))
    out["p_lo_same_mean"] = float(norm.sf(lo, mu, sd))
    out["p_hi_same_mean"] = float(norm.sf(hi, mu, sd))
    out["p_lo_h0_mean"] = float(norm.sf(lo, out["mean_h0"], sd))
    out["p_hi_h0_mean"] = float(norm.sf(hi, out["mean_h0"], sd))
    out["pru"] = {k: dict(p=float(norm.sf(lo, *v)), P_pass=float(norm.cdf(lo, *v)) ** N0)
                  for k, v in PRU.items()}
    out["pru_best_p"] = min(v["p"] for v in out["pru"].values())
    out["pru_best_pass"] = max(v["P_pass"] for v in out["pru"].values())
    z = float(norm.isf(P0))
    p80 = 1 - 0.8 ** (1 / N0)
    out.update(mean_for_p0=lo - z * sd, mean_for_80=lo - float(norm.isf(p80)) * sd,
               cross_for_p0=mu + z * sd, cross_for_80=mu + float(norm.isf(p80)) * sd)
    return out


# ------------------------------------------------------- rig, timing budget, dummy matrix (SYNTHETIC)
def hardware(grid, mod) -> dict:
    """Timing budget, rig energies, dummy matrix and calibration limits. Same method as HE's scoping."""
    A = arm()
    margin = min(v["t10"] for v in grid.values()) - CATCH_MS
    w10_min, w200_max = min(v["w10"] for v in grid.values()), max(v["w200"] for v in grid.values())
    t_stop_min, w_stop_max = min(v["t_stop"] for v in grid.values()), max(v["w_stop"] for v in grid.values())
    B = BUDGET_FRAC * margin
    fs_min = 1 / (ALLOC["sample"] * B)  # kHz, because B is in ms
    fs = math.ceil(fs_min / FS_STEP_KHZ) * FS_STEP_KHZ
    ang_max = w10_min * ALLOC["angle"] * B / 1e3
    bits = math.ceil(math.log2(360 / ang_max))
    enc_res = 360 / 2**bits
    ang_t = enc_res / w10_min * 1e3
    bud = dict(margin_ms=margin, budget_ms=B, fs_min_khz=fs_min, fs_khz=fs, period_ms=1 / fs,
               period_max_ms=ALLOC["sample"] * B, trigger_ms=ALLOC["trigger"] * B,
               release_ms=ALLOC["release"] * B, enc_lat_ms=ALLOC["enc_lat"] * B, spare_ms=ALLOC["spare"] * B,
               ang_max_deg=ang_max, counts_min=360 / ang_max, bits=bits, counts=2**bits, enc_res_deg=enc_res,
               ang_time_ms=ang_t, w10_min=w10_min, w200_max=w200_max)
    bud["u_t_ms"] = bud["period_ms"] + bud["trigger_ms"] + bud["release_ms"] + bud["enc_lat_ms"] + ang_t
    bud["headroom_ms"] = B - bud["u_t_ms"]
    bud["bw_enc_khz"] = RISE_K / bud["enc_lat_ms"]  # first-order rise time t_r = 0.35 / f (ASSUMPTION, HE)
    bud["ok"] = bool(bud["u_t_ms"] <= B)
    bud["ang_unc_deg"] = enc_res + w200_max * (bud["u_t_ms"] - ang_t) / 1e3
    bud["t10_min_pass_ms"] = CATCH_MS + bud["u_t_ms"]
    # logged data rate (ASSUMPTION: channel count and sample size; FC's count of the channels in section 7)
    bud["data_channels"], bud["data_bytes"], bud["data_trial_s"] = DATA_CHANNELS, DATA_BYTES, DATA_TRIAL_S
    bud["data_mb_s"] = DATA_CHANNELS * fs * 1e3 * DATA_BYTES / 1e6
    bud["data_mb_trial"] = bud["data_mb_s"] * DATA_TRIAL_S
    bud["drop200_max_pass_deg"] = DROP_LIM - bud["ang_unc_deg"]
    bud["worst_drop200_deg"] = max(v["drop200"] for v in grid.values())
    # dummy matrix
    z = float(norm.ppf(Z_TWO_SIDED))
    n_formula = (z * S_DUMMY_MS / E_DUMMY_MS) ** 2
    n_cell = max(N_MIN_CELL, math.ceil(n_formula))
    mat = dict(n_cells=len(grid), n_types=len(grid) // len(HOLD_FRACS), n_cell=n_cell, n_formula=n_formula,
               drops=n_cell * len(grid), hours=n_cell * len(grid) * T_RESET_MIN / 60,
               holds=[hf * A["hold_nom"] for hf in HOLD_FRACS], n_b=sum(k[0] == "b" for k in grid) // 3,
               n_iv=sum(k[0] == "iv" for k in grid) // 3)
    settings = {}
    for (o, p, r, sto, hf) in grid:
        settings.setdefault((p, round(r, 9), sto, hf), []).append(o)
    mat["distinct"] = len(settings)
    mat["dup_cells"] = len(grid) - len(settings)
    mat["dup_types"] = mat["dup_cells"] // len(HOLD_FRACS)
    mat["dup_drops"] = mat["dup_cells"] * n_cell
    mat["dup_hours"] = mat["dup_drops"] * T_RESET_MIN / 60
    mat["distinct_drops"] = len(settings) * n_cell
    mat["distinct_hours"] = mat["distinct_drops"] * T_RESET_MIN / 60
    up = grid[("iv", "comm", RATE["iv"]["comm"] * (1 + TOL), STO["comm"][0], 0.0)]
    mat["upper_edge_t10"] = up["t10"]
    mat["rows"] = [dict(cell=f"{'(b)' if o == 'b' else '(iv)'} {p} {r:.2f} N·m/s, STO "
                        f"{'none' if sto is None else f'{sto:.0f} ms'}",
                        t10=[grid[(o, p, r, sto, hf)]["t10"] for hf in HOLD_FRACS])
                   for (o, p, r, sto, hf0) in grid if hf0 == 0.0]
    # rig energies and speeds (Rig fault: torque lost at t = 0, or the device driving at its cap)
    th_s, th_k, th_u = (math.radians(x) for x in (STOP_DEG, RIG["th_stroke"], RIG["th_up"]))
    J, Gm, kcap = A["J"], A["G"], K_HI * T_CAP_DEV
    ff_t10 = crossing("estop", 20.0, None, 0.0, step_loss=True)
    ff_stop = crossing("estop", 20.0, None, 0.0, step_loss=True, ang=STOP_DEG)
    rig = dict(kcap=kcap, ff_t10_ms=ff_t10[0], ff_stop_ms=ff_stop[0], ff_stop_w=ff_stop[1],
               E_ff_stop=Gm * math.sin(th_s), E_ff_total=Gm * math.sin(th_s + th_k),
               E_cap_total=Gm * math.sin(th_s + th_k) + kcap * (th_s + th_k),
               E_cap_stop=Gm * math.sin(th_s) + kcap * th_s,
               E_up=(kcap + A["hold_nom"]) * th_u - Gm * math.sin(th_u))
    rig["w_ff_energy"] = math.degrees(math.sqrt(2 * rig["E_ff_stop"] / J))
    rig["w_cap_stop"] = math.degrees(math.sqrt(2 * rig["E_cap_stop"] / J))
    rig["w_up"] = math.degrees(math.sqrt(2 * rig["E_up"] / J))
    rig["T_abs_ideal"] = rig["E_cap_total"] / th_k
    rig["T_abs_peak"] = rig["T_abs_ideal"] * RIG["abs_profile"]
    rig["F_abs_peak"] = rig["T_abs_peak"] / RIG["r_stopc"]
    rig["T_up_ideal"] = rig["E_up"] / th_k + kcap + A["hold_nom"]
    rig["F_up_peak"] = rig["T_up_ideal"] * RIG["abs_profile"] / RIG["r_stopc"]
    rig["trip_Nm"] = RIG["trip_frac"] * A["kT0"]
    rig["trip_pct"] = (RIG["trip_frac"] - 1) * 100
    r_pay = ARM["L_fa"] + ARM["pay_off"]
    rig["r_tip"] = r_pay + RIG["r_tip_extra"]
    rig["r_guard"] = rig["r_tip"] + RIG["clr"]
    rig["ang_bottom"] = 90 - STOP_DEG - RIG["th_stroke"]
    rig["ang_top"] = 90 + RIG["th_up"] + RIG["th_stroke"]
    rig["sweep_mm"] = rig["r_tip"] * math.radians(rig["ang_top"] - rig["ang_bottom"]) * 1e3
    first15 = min(grid, key=lambda c: grid[c]["t_stop"])
    gov = grid[("b",) + tuple(mod["gov_key"]) + (0.0,)]
    rig["first15_cell"] = "|".join(map(str, first15))
    rig["gap_first15_ms"] = grid[first15]["t_stop"] - grid[first15]["t10"]  # same cell
    rig["gap_gov_ms"] = gov["t_stop"] - gov["t10"]
    rig["t_stop_gov"], rig["w_stop_gov"] = gov["t_stop"], gov["w_stop"]
    rig.update(t_stop_min=t_stop_min, w_stop_max=w_stop_max)
    # flip values at the governing cell (hold 0): the parameter change that takes the drop at the
    # 200 ms catch to exactly 10 deg, i.e. where the model would pass a unit that really fails
    path, rate, sto = mod["gov_key"][0], mod["gov_key"][1], mod["gov_key"][2]

    def d200(**kw):
        return drop_at(CATCH_MS, path, rate, sto, kw.pop("hold", 0.0), **kw)[0] - DROP_LIM

    flip = dict(
        J=brentq(lambda m: d200(J=A["J"] * m), 0.3, 1.0, xtol=1e-9),
        mass=brentq(lambda m: d200(J=arm(m)["J"], G=arm(m)["G"]), 1.0, 1.3, xtol=1e-9),
        T0=brentq(lambda m: d200(T0=T0_DEV * m), 0.9, 1.0, xtol=1e-9),
        k=brentq(lambda m: d200(k=K_HI * m), 0.9, 1.0, xtol=1e-9),
        offset=brentq(lambda dh: d200(hold=-dh), 0.0, 2.0, xtol=1e-9), lag_ms=margin)
    h = 1e-3

    def t10(**kw):
        return crossing(path, rate, sto, kw.pop("hold", 0.0), **kw)[0]

    sl = dict(Th=(t10(hold=h) - t10(hold=-h)) / (2 * h),
              G=(t10(G=A["G"] * 1.01) - t10(G=A["G"] * 0.99)) / 2,
              J=(t10(J=A["J"] * 1.01) - t10(J=A["J"] * 0.99)) / 2,
              k=(t10(k=K_HI * 1.01) - t10(k=K_HI * 0.99)) / 2)
    ti = TOL_MODEL_MS * INPUT_SHARE
    cal = dict(share_ms=ti, net_torque_Nm=ti / abs(sl["Th"]), G_pct=ti / abs(sl["G"]),
               J_pct=ti / abs(sl["J"]), k_pct=ti / abs(sl["k"]), slope_Th=sl["Th"], slope_G=sl["G"],
               slope_J=sl["J"], slope_k=sl["k"])
    cal["G_Nm"] = cal["G_pct"] / 100 * A["G"]
    flip_pct = {k: (v - 1) * 100 for k, v in flip.items() if k in ("J", "mass", "T0", "k")}
    return dict(budget=bud, matrix=mat, rig=rig, flip=flip, flip_pct=flip_pct, cal=cal,
                gov_cell=mod["gov"], t_stop_min=t_stop_min)


# ---------------------------------------------------------------- sample size (exact)
def n_allowed(c, p=P0, alpha=1 - CONF):
    n = c + 1
    while binom.cdf(c, n, p) > alpha:
        n += 1
    return n


def cp_upper(x, n):
    """One-sided Clopper-Pearson upper bound."""
    return 1.0 if x >= n else float(beta.ppf(CONF, x + 1, n - x))


def deff(m, icc: Fr) -> Fr:
    return 1 + (m - 1) * Fr(icc)


def subjects(n, m, icc: Fr) -> int:
    """ceil(n * design effect / m), in exact arithmetic."""
    return math.ceil(n * deff(m, icc) / m)


def pmf_bb(p, icc, m, c):
    """P(0..c drops in one subject), beta-binomial with mean p and correlation icc."""
    icc, ks = float(icc), np.arange(c + 1)
    if icc == 0:
        return binom.pmf(ks, m, p)
    if icc == 1:
        out = np.zeros(c + 1)
        out[0] = 1 - p
        if c >= m:
            out[m] = p
        return out
    f = (1 - icc) / icc
    return betabinom.pmf(ks, m, p * f, (1 - p) * f)


def pmf_two(p, icc, m, c):
    """Same, when a share w of subjects drops with probability q and the rest never drop."""
    q = p + (1 - p) * float(icc)
    out = p / q * binom.pmf(np.arange(c + 1), m, q)
    out[0] += 1 - p / q
    return out


def p_total_le(pm, k, c):
    """P(total drops over k iid subjects <= c), pm = per-subject pmf on 0..c."""
    acc = np.zeros(c + 1)
    acc[0] = 1.0
    base, n = np.array(pm, float), k
    while n:
        if n & 1:
            acc = np.convolve(acc, base)[: c + 1]
        base = np.convolve(base, base)[: c + 1]
        n >>= 1
    return float(acc.sum())


def subjects_model(pmf, c, m, icc):
    k = 1
    while p_total_le(pmf(P0, icc, m, c), k, c) > 1 - CONF:
        k += 1
    return k


def worst_zero(p, icc, m, n_grid=2001):
    """Largest P(0 drops in m trials) over every distribution of subject drop rates with mean p and
    the given ICC (linear program over a grid of rates)."""
    x = np.linspace(0, 1, n_grid)
    v = p * p + float(icc) * p * (1 - p)
    s = linprog(-((1 - x) ** m), A_eq=np.array([np.ones_like(x), x, x * x]), b_eq=[1, p, v],
                bounds=(0, None))
    return -float(s.fun)


def staged_pass(ns, p):
    """P(pass) of a test that passes at stage i (total n = ns[i]) if total drops <= i, else goes on."""
    alive, prev, out = {0: 1.0}, 0, 0.0
    for i, n in enumerate(ns):
        new = {}
        for d, w in alive.items():
            for j in range(len(ns) - d):
                new[d + j] = new.get(d + j, 0.0) + w * float(binom.pmf(j, n - prev, p))
        out += sum(w for d, w in new.items() if d <= i)
        alive, prev = {d: w for d, w in new.items() if d > i}, n
    return out


def staged_pass_subjects(ks, pm):
    """staged_pass for cumulative subject counts ks and a per-subject drop pmf on 0..len(ks)-1."""
    c = len(ks) - 1
    alive, prev, out = np.eye(c + 1)[0], 0, 0.0
    for i, k in enumerate(ks):
        for _ in range(k - prev):
            alive = np.convolve(alive, pm)[: c + 1]
        out += float(alive[: i + 1].sum())
        alive[: i + 1] = 0.0
        prev = k
    return out


def extension():
    """Smallest totals n1, n2 that keep the whole staged test at <= 5 % with the first look at N0."""
    n1 = N0 + 1
    while staged_pass([N0, n1], P0) > 1 - CONF:
        n1 += 1
    n2 = n1 + 1
    while staged_pass([N0, n1, n2], P0) > 1 - CONF:
        n2 += 1
    return N0, n1, n2


def icc_binary(icc_latent, p=P0):
    """ICC of the drop/no-drop outcome if catch time is normal with subject share icc_latent."""
    z, r = norm.isf(p), icc_latent
    e2 = quad(lambda u: norm.pdf(u) * norm.sf((z - math.sqrt(r) * u) / math.sqrt(1 - r)) ** 2, -10, 10)[0]
    return (e2 - p * p) / (p * (1 - p))


def stats() -> dict:
    alpha = 1 - CONF
    ext = extension()
    nc = [n_allowed(c) for c in range(3)]
    out = dict(n0=nc[0], q59=0.95**N0, q58=0.95 ** (N0 - 1), n1=nc[1], n2=nc[2],
               cp=[cp_upper(x, N0) for x in range(4)],
               size=[float(binom.cdf(c, n, P0)) for c, n in enumerate(nc)],
               naive1=staged_pass(nc[:2], P0), naive2=staged_pass(nc, P0), ext=ext,
               ext_size=staged_pass(ext, P0), retest=[1 - (1 - 0.95**N0) ** r for r in (2, 3)],
               p80=1 - 0.8 ** (1 / N0), p50=1 - 0.5 ** (1 / N0),
               oc=[dict(p=p, P=[float(binom.cdf(c, n, p)) for c, n in enumerate(nc)]) for p in P_ROWS])
    out["killip"] = dict(deff=float(deff(32, Fr("0.017"))), ess=float(128 / deff(32, Fr("0.017"))))
    grid = {}
    for m in MS:
        for r in ICCS:
            grid[f"{m}|{float(r):g}"] = dict(
                deff=float(deff(m, r)), subjects=subjects(N0, m, r), trials=m * subjects(N0, m, r),
                bb=subjects_model(pmf_bb, 0, m, r), worst=subjects_model(pmf_two, 0, m, r),
                false_pass_worst=p_total_le(pmf_two(P0, r, m, 0), subjects(N0, m, r), 0))
    out["grid"] = grid
    out["range_he"] = [grid[f"10|{v:g}"]["subjects"] for v in (0, 0.2, 0.5)]
    out["worst_lp_gap"] = max(abs(worst_zero(P0, r, m) - pmf_two(P0, r, m, 0)[0])
                              for m in (2, M_ALT) for r in ICCS)
    fine = [Fr(i, 100) for i in range(101)]

    def covered(m):
        return all(p_total_le(pmf_two(P0, r, m, 0), subjects(N0, m, r), 0) <= alpha for r in fine)

    out["deff_ok_to_m"] = next(m for m in range(1, 200) if not covered(m + 1))
    out["allowed"] = {f"{c}|{float(r):g}": dict(deff=subjects(nc[c], M_ALT, r),
                                                bb=subjects_model(pmf_bb, c, M_ALT, r),
                                                two=subjects_model(pmf_two, c, M_ALT, r))
                      for c in (0, 1, 2) for r in ICCS}
    out["allowed_short"] = [k for k, v in out["allowed"].items() if v["deff"] < max(v["bb"], v["two"])]
    out["ext_subjects"] = {f"{float(r):g}": [subjects(n, M_ALT, r) for n in ext] for r in ICCS}
    out["ext_clustered_max"] = max(
        staged_pass_subjects(out["ext_subjects"][f"{float(r):g}"], f(P0, r, M_ALT, 2))
        for r in ICCS for f in (pmf_bb, pmf_two))
    out["latent_icc"] = [0.25, 0.5, 0.75, 0.9]
    out["latent"] = [icc_binary(r) for r in out["latent_icc"]]
    out["range_icc"] = [0.0, 0.2, 0.5]
    k = subjects(N0, M_ALT, ICC_ALT)
    # plan sized for ICC 0 when the true ICC is ICC_ALT: false-pass chance at p = P0
    k0 = subjects(N0, M_ALT, Fr(0))
    out["misspec"] = [p_total_le(f(P0, ICC_ALT, M_ALT, 0), k0, 0) for f in (pmf_bb, pmf_two)]
    out["chosen"] = dict(m=M_CHOSEN, subjects=subjects(N0, M_CHOSEN, Fr(0)),
                         false_pass=[p_total_le(f(P0, r, M_CHOSEN, 0), N0, 0)
                                     for r in ICCS for f in (pmf_bb, pmf_two)])
    out["chosen"]["false_pass_max"] = max(out["chosen"]["false_pass"])
    out["chosen"]["false_pass_min"] = min(out["chosen"]["false_pass"])
    out["alt_icc_break"] = float((Fr(k * M_ALT, N0) - 1) / (M_ALT - 1))  # n_eff < N0 above this ICC
    out["alt_risk"] = [dict(icc=float(r), worst=p_total_le(pmf_two(P0, r, M_ALT, 0), k, 0),
                            bb=p_total_le(pmf_bb(P0, r, M_ALT, 0), k, 0)) for r in ALT_TRUE_ICCS]
    out["plan"] = dict(m=M_ALT, icc=float(ICC_ALT), subjects=k, trials=k * M_ALT,
                       deff=float(deff(M_ALT, ICC_ALT)), subjects_any=subjects(N0, M_ALT, Fr(1)),
                       trials_any=subjects(N0, M_ALT, Fr(1)) * M_ALT)
    return out


def _plab(name, p) -> str:
    """Row label for a true drop probability: whole percent, but two decimals for the window-edge row
    (47.95 %, not 48 %; FC PR #9 round-2 nit)."""
    return f"{100 * p:.2f}" if name.startswith("window") else f"{100 * p:.0f}"


# ---------------------------------------------------------------- feasibility (SYNTHETIC)
def feasibility(mod, hw) -> dict:
    ps = {"FC low": FC_P[0], "FC high": FC_P[1], "window (b) edge": mod["p_lo_same_mean"]}
    nc = [n_allowed(c) for c in range(3)]
    out = {"ps": ps, "indep": {k: [float(binom.cdf(c, n, p)) for c, n in enumerate(nc)]
                               for k, p in ps.items()}}
    out["clustered"] = {}
    for k, p in ps.items():
        for r in ICCS:
            n_s = subjects(N0, M_ALT, r)
            out["clustered"][f"{k}|{float(r):g}"] = dict(
                subjects=n_s, bb=p_total_le(pmf_bb(p, r, M_ALT, 0), n_s, 0),
                two=p_total_le(pmf_two(p, r, M_ALT, 0), n_s, 0))
    # largest chance of passing over every option offered, for FC's range
    best = []
    for p in FC_P:
        best += [float(binom.cdf(c, n, p)) for c, n in enumerate(nc)] + [staged_pass(extension(), p)]
        for m in MS:
            for r in ICCS:
                ext = [subjects(n, m, r) for n in extension()]
                for f in (pmf_bb, pmf_two):
                    best += [p_total_le(f(p, r, m, c), subjects(n, m, r), c) for c, n in enumerate(nc)]
                    best.append(staged_pass_subjects(ext, f(p, r, m, 2)))
    out["max_pass"], out["one_in"] = max(best), 1 / max(best)
    # same bound over every m the plan's design-effect check allows (1..46) and a fine ICC sweep, c = 0
    sup = {}
    for p in FC_P:
        top = (0.0, None)
        for m in range(1, SUP_M_MAX + 1):
            for i in range(SUP_ICC_STEPS + 1):
                r = Fr(i, SUP_ICC_STEPS)
                k = subjects(N0, m, r)
                v = p_total_le(pmf_two(p, r, m, 0), k, 0)
                if v > top[0]:
                    top = (v, dict(m=m, icc=float(r), subjects=k))
        sup["lo" if p == FC_P[0] else "hi"] = dict(P=top[0], **top[1])
    out["sup_any_m"] = sup
    out["chosen_pass"] = {k: (1 - p) ** N0 for k, p in ps.items()}
    # stop contact in Phase B: catch later than the 15 deg stop (no-catch arm reaches it), S7 catch model
    sd, means = mod["fc_sd"], {"nominal hold": mod["fc_mean"], "zero hold": mod["mean_h0"]}
    rg = hw["rig"]
    out["contact"] = {}
    for lab, mu in means.items():
        pc_gov = float(norm.sf(rg["t_stop_gov"], mu, sd))
        pc_first = float(norm.sf(rg["t_stop_min"], mu, sd))
        pd = float(norm.sf(mod["lo"], mu, sd))  # any DROP (catch after the (b) edge) halts the study
        pd_ut = float(norm.sf(mod["lo"] - hw["budget"]["u_t_ms"], mu, sd))  # DROP rule: t_c + U_t > edge
        out["contact"][lab] = dict(mean=mu, p_contact_gov=pc_gov, p_contact_first15=pc_first, p_drop=pd,
                                   p_study_contact=pc_gov / pd * (1 - (1 - pd) ** N0),
                                   p_study_contact_ut=pc_gov / pd_ut * (1 - (1 - pd_ut) ** N0))
    out["first_drop"] = {k: 1 / p for k, p in ps.items()}
    out["drop_in_10"] = {k: 1 - (1 - p) ** 10 for k, p in ps.items()}
    out["drop_in_10_icc"] = {k: 1 - float(pmf_bb(p, ICC_ALT, 10, 0)[0]) for k, p in ps.items()}
    out["drop_in_2_subj"] = {k: 1 - (1 - v) ** 2 for k, v in out["drop_in_10_icc"].items()}
    out["undersized"] = p_total_le(pmf_bb(FC_P[0], ICC_ALT, M_ALT, 0), subjects(N0, M_ALT, Fr(0)), 0)
    out["undersized_best"] = p_total_le(pmf_two(FC_P[0], ICC_ALT, M_ALT, 0), subjects(N0, M_ALT, Fr(0)), 0)
    return out


def reference_check(mod, st) -> dict:
    return dict(window_diff_ms=[mod["lo"] - REF["window"][0], mod["hi"] - REF["window"][1]],
                window_rounds=[round(mod["lo"], 2) == WINDOW_MS[0], round(mod["hi"], 2) == WINDOW_MS[1]],
                hi_round_up_ms=WINDOW_MS[1] - mod["hi"], lo_round_down_ms=mod["lo"] - WINDOW_MS[0],
                n=[st["n0"], st["n1"], st["n2"]] == list(REF["n_allow"]),
                subjects=st["range_he"] == list(REF["subjects"]),
                allowed_short_n=len(st["allowed_short"]),
                fc_p_diff=[mod["fc_p_lo"] - REF["fc_p"][0], mod["fc_p_hi"] - REF["fc_p"][1]],
                h0_p_diff=mod["p_lo_same_mean"] - REF["h0_p_worst"])


def compute() -> dict:
    grid = build_grid()
    mod, st = model(grid), stats()
    hw = hardware(grid, mod)
    feas = feasibility(mod, hw)
    inp = dict(win_lo=WINDOW_MS[0], win_hi=WINDOW_MS[1], p0=P0, conf=CONF, n0=N0, drop=DROP_LIM,
               catch=CATCH_MS, stop=STOP_DEG, res=RES_MS, rate=RATE["b"]["estop"] * (1 + TOL),
               rate_iv=RATE["iv"]["estop"] * (1 + TOL), sto=STO["estop"][0], sto_comm=STO["comm"][0],
               tol_model=TOL_MODEL_MS, emg_hz=EMG_HZ, imu_hz=IMU_HZ, imu_ms=1e3 / IMU_HZ,
               daq_khz=1 / RES_MS, lp_tol=LP_TOL, m_alt=M_ALT, icc_alt=float(ICC_ALT), m_chosen=M_CHOSEN,
               rate_dup=RATE["iv"]["estop"], sup_steps=SUP_ICC_STEPS, rise_k=RISE_K,
               k=K_HI, T0=T0_DEV, tol=TOL,
               g_unc=G_UNC, fc_lo=FC_P[0], fc_hi=FC_P[1], pru_sd=mod["fc_sd"],
               rate_nom=RATE["b"]["estop"], sto_nom=STO["estop"][1], budget_pct=BUDGET_FRAC,
               step_khz=FS_STEP_KHZ,
               share=INPUT_SHARE, s_dummy=S_DUMMY_MS, e_dummy=E_DUMMY_MS, n_min_cell=N_MIN_CELL,
               reset_min=T_RESET_MIN, th_stroke=RIG["th_stroke"], th_up=RIG["th_up"], r_stopc=RIG["r_stopc"],
               clr_mm=RIG["clr"] * 1e3, abs_profile=RIG["abs_profile"], trip_pct=(RIG["trip_frac"] - 1),
               alloc=dict(ALLOC), hold_fracs=list(HOLD_FRACS), t_cap=T_CAP_DEV)
    return {"NOTE": "keys under 'syn' are SYNTHETIC (model output, not bench or human data)",
            "in": inp, "stat": st, "syn": dict(mod, feas=feas, hw=hw), "ref": reference_check(mod, st)}


# ---------------------------------------------------------------- rendering
def _pct(x):
    return f"{100 * x:.2f} %"


def tables(R) -> dict:
    g, al, f, st = R["stat"]["grid"], R["stat"]["allowed"], R["syn"]["feas"], R["stat"]
    rs = [f"{float(r):g}" for r in ICCS]
    t = {}
    t["subjects"] = ["| Trials per subject m | " + " | ".join(f"ICC {r}" for r in rs) + " |",
                     "|---" * (len(rs) + 1) + "|"]
    t["subjects"] += [f"| {m} | " + " | ".join(
        f"{g[f'{m}|{r}']['subjects']} ({g[f'{m}|{r}']['trials']})" for r in rs) + " |" for m in MS]
    t["check"] = ["| ICC | Design effect | Subjects, design effect (rejected) | Subjects, beta-binomial "
                  "| Subjects, worst case | False-pass chance, plan size, worst case |", "|---" * 6 + "|"]
    t["check"] += [f"| {r} | {v['deff']:.2f} | {v['subjects']} | {v['bb']} | {v['worst']} | "
                   f"{_pct(v['false_pass_worst'])} |" for r in rs for v in [g[f"{M_ALT}|{r}"]]]
    t["allowed"] = ["| Allowed drops c | n (independent) | Size at 5 % | Subjects (one trial each) |",
                    "|---" * 4 + "|"]
    t["allowed"] += [f"| {c} | {st[f'n{c}']} | {100 * st['size'][c]:.3f} % | {st[f'n{c}']} |"
                     for c in (0, 1, 2)]
    t["ext"] = ["| Look | Pass if total drops at most | Total subjects (one trial each) |", "|---" * 3 + "|"]
    t["ext"] += [f"| {i + 1} | {i} | {n} |" for i, n in enumerate(st["ext"])]
    t["alt_risk"] = ["| True drop ICC | False pass, worst-case mix | False pass, beta-binomial "
                     f"| Chosen {N0} × one |", "|---" * 4 + "|"]
    t["alt_risk"] += [f"| {r['icc']:g} | {_pct(r['worst'])} | {_pct(r['bb'])} | "
                      f"{_pct(st['chosen']['false_pass_max'])} |" for r in st["alt_risk"]]
    t["oc"] = ["| True drop probability | Pass, c = 0 (n = %d) | Pass, c = 1 (n = %d) "
               "| Pass, c = 2 (n = %d) |" % (st["n0"], st["n1"], st["n2"]), "|---" * 4 + "|"]
    t["oc"] += [f"| {100 * r['p']:g} % | " + " | ".join(_pct(x) for x in r["P"]) + " |" for r in st["oc"]]
    t["feas"] = ["| SYNTHETIC true drop probability | c = 0, n = %d (chosen) | c = 1, n = %d | c = 2, n = %d "
                 "| Mean subjects to first drop |" % (st["n0"], st["n1"], st["n2"]), "|---" * 5 + "|"]
    t["feas"] += [f"| {_plab(k, p)} % ({k}) | " + " | ".join(f"{x:.1e}" for x in f["indep"][k])
                  + f" | {f['first_drop'][k]:.1f} |" for k, p in f["ps"].items()]
    t["feas_icc"] = ["| SYNTHETIC: ICC | Subjects (trials) | " + " | ".join(
        f"{_plab(k, p)} %: beta-binomial / best case" for k, p in f["ps"].items()) + " |",
        "|---" * (2 + len(f["ps"])) + "|"]
    t["feas_icc"] += [f"| {r} | {v0['subjects']} ({v0['subjects'] * M_ALT}) | " + " | ".join(
        f"{f['clustered'][f'{k}|{r}']['bb']:.1e} / {f['clustered'][f'{k}|{r}']['two']:.1e}"
        for k in f["ps"]) + " |" for r in rs for v0 in [f["clustered"][f"FC low|{r}"]]]
    hw, ni = R["syn"]["hw"], R["in"]
    b, rg, mt, cl, fl, al = hw["budget"], hw["rig"], hw["matrix"], hw["cal"], hw["flip"], ni["alloc"]

    def tab(head, rows):
        return [head, "|---" * (head.count("|") - 1) + "|"] + [f"| {' | '.join(r)} |" for r in rows]

    pc = lambda x: f"{100 * x:.0f} %"  # noqa: E731
    t["budget"] = tab("| SYNTHETIC timing budget item | Requirement | Derivation |", [
        ("Smallest model margin", f"{b['margin_ms']:.3f} ms",
         f"governing cell, over all {mt['n_cells']} cells"),
        ("Budget", f"{b['budget_ms']:.3f} ms", f"{pc(ni['budget_pct'])} of the margin"),
        ("Sampling", f"at least {b['fs_khz']:.0f} kHz",
         f"period at most {b['period_max_ms']:.3f} ms ({pc(al['sample'])}) needs at least "
         f"{b['fs_min_khz']:.1f} kHz; rounded up to {ni['step_khz']} kHz steps"),
        ("Trigger skew", f"at most {b['trigger_ms']:.3f} ms", pc(al["trigger"])),
        ("Release skew", f"at most {b['release_ms']:.3f} ms", pc(al["release"])),
        ("Encoder latency", f"at most {b['enc_lat_ms']:.3f} ms", pc(al["enc_lat"])),
        ("Encoder bandwidth", f"at least {b['bw_enc_khz']:.1f} kHz",
         f"{ni['rise_k']:g} / encoder latency (first-order rise time, ASSUMPTION)"),
        ("Angle resolution", f"at most {b['ang_max_deg']:.4f}° ({b['bits']}-bit: {b['counts']:,} counts, "
         f"{b['enc_res_deg']:.5f}°)", f"{pc(al['angle'])} at the slowest {b['w10_min']:.1f} °/s crossing; "
         f"{b['ang_time_ms']:.4f} ms"),
        ("Spare", f"{b['spare_ms']:.3f} ms", pc(al["spare"])),
        ("**Used, U_t**", f"**at most {b['u_t_ms']:.3f} ms**",
         f"{b['period_ms']:.4f} + {b['trigger_ms']:.4f} + {b['release_ms']:.4f} + {b['enc_lat_ms']:.4f} + "
         f"{b['ang_time_ms']:.4f} ms, against {b['budget_ms']:.3f} ms"),
        ("Headroom", f"{b['headroom_ms']:.3f} ms", f"{b['budget_ms']:.3f} − {b['u_t_ms']:.3f} ms"),
        ("Angle uncertainty", f"at most {b['ang_unc_deg']:.4f}°",
         f"encoder step plus {b['w200_max']:.1f} °/s times the timing terms")])
    t["absorber"] = tab("| SYNTHETIC rig case | Energy | Speed or force |", [
        (f"No-hold free fall, device torque lost at t = 0 (stop at {rg['ff_stop_ms']:.1f} ms; 10° at "
         f"{rg['ff_t10_ms']:.1f} ms)", f"{rg['E_ff_stop']:.3f} J at the stop, "
         f"{rg['E_ff_total']:.3f} J through the stroke", f"{rg['ff_stop_w']:.0f} °/s at the stop"),
        (f"Device driving the arm down at its cap (HE's runaway case, k × {ni['t_cap']:.3f} N·m = "
         f"{rg['kcap']:.2f} N·m)", f"{rg['E_cap_stop']:.3f} J at the stop, "
         f"**{rg['E_cap_total']:.2f} J through the stroke**", f"{rg['w_cap_stop']:.0f} °/s at the stop"),
        ("Absorber, constant force over the stroke", f"ideal {rg['T_abs_ideal']:.1f} N·m; peak "
         f"{rg['T_abs_peak']:.1f} N·m (× {ni['abs_profile']:g}, ASSUMPTION)",
         f"{rg['F_abs_peak']:.0f} N at {ni['r_stopc']:.2f} m"),
        (f"Upper runaway stop at +{ni['th_up']:.0f}°, device driving up at its cap", f"{rg['E_up']:.3f} J",
         f"{rg['w_up']:.0f} °/s; ideal {rg['T_up_ideal']:.1f} N·m, {rg['F_up_peak']:.0f} N peak"),
        (f"Rig torque trip (+{pc(ni['trip_pct'])} over k × T0)", f"{rg['trip_Nm']:.1f} N·m at the arm",
         "ASSUMPTION")])
    pf = hw["flip_pct"]
    t["flip"] = tab("| SYNTHETIC governing cell, hold 0 | Flip value (drop at the catch reaches "
                    f"{ni['drop']:.0f}°) | Implied calibration ({cl['share_ms']:.3f} ms each) |", [
        ("Dummy human-side mass", f"{pf['mass']:+.2f} %", ""),
        ("Inertia J", f"{pf['J']:+.2f} %", f"±{cl['J_pct']:.2f} %"),
        ("Gravity moment G", "", f"±{cl['G_pct']:.3f} % (±{cl['G_Nm']:.4f} N·m)"),
        ("Net pre-onset torque, downward offset", f"{fl['offset']:.3f} N·m",
         f"±{cl['net_torque_Nm']:.4f} N·m"),
        ("Clamp T0", f"{pf['T0']:+.2f} %", ""),
        ("Torque ratio k", f"{pf['k']:+.2f} %", f"±{cl['k_pct']:.3f} %"),
        ("Timing lag", f"more than {fl['lag_ms']:.2f} ms", f"U_t {b['u_t_ms']:.3f} ms")])
    t["matrix"] = tab("| SYNTHETIC cell: option, path, ramp, STO | Crossing at hold 0 (ms) | "
                      f"at {mt['holds'][1]:.3f} N·m (ms) | at {mt['holds'][2]:.3f} N·m (ms) |",
                      [(r["cell"], *(f"{x:.2f}" for x in r["t10"])) for r in mt["rows"]])
    t["pru"] = ["| SYNTHETIC: PRU09 condition (mean ± SD, ms) | P(catch later than the (b) edge) "
                f"| P(0 drops in {N0}) |", "|---" * 3 + "|"]
    t["pru"] += [f"| {k} ({PRU[k][0]:.0f} ± {PRU[k][1]:.0f}) | {100 * v['p']:.1f} % | {v['P_pass']:.1e} |"
                 for k, v in R["syn"]["pru"].items()]
    return {k: "\n".join(v) for k, v in t.items()}


# Numbers allowed as literals in the template (everything else must be a placeholder).
LITERALS = {"1": "counts, list items, N·m background load (PRU09)", "2": "counts, c", "3": "counts, N·m",
            "0": "zero drops, c = 0", "4": "PRU09 foreperiod upper limit, s", "5": "section number",
            "6": "section number", "7": "section number", "8": "section number", "9": "section number",
            "14": "row number (Rev U, S7 verdict)",
            "23": "Declaration of Helsinki paragraph", "2026-11-26": "date of Charlie's design decision",
            "01": "time of Charlie's decision (01:25 KST)", "25": "time of Charlie's decision (01:25 KST)",
            "15066": "standard number, ISO/TS 15066", "12": "row number in FC's PR9_verdict.md",
            "90": "elbow start angle, deg"}
SYN_KEYS = ("syn.", "in.win_", "in.fc_", "ref.lo_", "ref.hi_")
PLACEHOLDER = re.compile(r"\{\{([^{}|]+)(?:\|([^{}]*))?\}\}")


def lookup(R, key):
    for part in key.split("."):
        R = R[int(part)] if isinstance(R, (list, tuple)) else R[part]
    return R


def render(R, text) -> str:
    tab = tables(R)

    def sub(m):
        key, spec = m.group(1), m.group(2) or ""
        if key.startswith("table:"):
            return tab[key[6:]]
        v = lookup(R, key)
        return _pct(v).replace(".00 %", " %") if spec == "pct" else format(v, spec)

    return PLACEHOLDER.sub(sub, text)


def template_problems(text) -> list[str]:
    """Literal numbers without provenance, and SYNTHETIC values on a line without the label."""
    body = text.split("## References")[0]
    bad = []
    for i, line in enumerate(body.splitlines(), 1):
        keys = [m.group(1) for m in PLACEHOLDER.finditer(line)]
        if any(k.startswith(SYN_KEYS) for k in keys) and "SYNTHETIC" not in line:
            bad.append(f"line {i}: SYNTHETIC label missing")
        rest = re.sub(r"`[^`]*`", " ", PLACEHOLDER.sub(" ", line))
        for tok in re.findall(r"(?<![\w.\-/])\d[\d.\-]*\d|(?<![\w.\-/])\d", rest):
            if tok not in LITERALS:
                bad.append(f"line {i}: literal {tok!r} has no provenance")
    return bad


def close(a, b, rel=1e-6) -> bool:
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(close(a[k], b[k], rel) for k in a)
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(close(x, y, rel) for x, y in zip(a, b))
    if isinstance(a, float) or isinstance(b, float):
        return math.isclose(a, b, rel_tol=rel, abs_tol=1e-8)
    return a == b


def check(R, plan=PLAN, numbers=NUMBERS, template=TEMPLATE) -> list[str]:
    text = template.read_text(encoding="utf-8")
    bad = template_problems(text)
    if render(R, text) != plan.read_text(encoding="utf-8"):
        bad.append("docs/ba15-testplan.md differs from the rendered template")
    if not close(json.loads(json.dumps(R)), json.loads(numbers.read_text(encoding="utf-8"))):
        bad.append("docs/ba15_numbers.json differs from a fresh run")
    ref = R["ref"]
    if R["stat"]["worst_lp_gap"] > LP_TOL or R["stat"]["ext_clustered_max"] > 1 - CONF:
        bad.append("a claim in the plan (worst case / extension under clustering) no longer holds")
    if not (all(ref["window_rounds"]) and ref["n"] and ref["subjects"]):
        bad.append(f"decided values not reproduced: {ref}")
    return bad


def main(argv) -> int:
    R = compute()
    if "--write" in argv:
        NUMBERS.write_text(json.dumps(R, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        PLAN.write_text(render(R, TEMPLATE.read_text(encoding="utf-8")), encoding="utf-8")
    bad = check(R)
    print("\n".join(bad) or "ba15 numcheck OK: every number in docs/ba15-testplan.md reproduces")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
