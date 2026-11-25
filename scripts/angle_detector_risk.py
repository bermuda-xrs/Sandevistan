"""Angle-detector risk for the B-A15 / B9 "drop past 10 deg" score. ALL NUMBERS ARE SYNTHETIC.

The "angle detector" is the joint-angle chain that turns a physical arm drop into the B-A15 trial
score (PR #9, docs/ba15-testplan.md @ 58cce8d, section 3 and "Trial score"):
angle sensor -> sampling -> zero reference -> peak drop / 10 deg crossing / catch time -> timestamps
against fault onset -> DROP / NO DROP. A missed detection is a false pass of the B9 safety test; a
false trip is a spurious DROP, and under the zero-drop rule (0 in 59) one spurious DROP fails the test.

No recorded human or bench data is read. Nothing here imports or changes the detectors, the 68.67 ms
gate, PLANNED_SHAPE or any repo default; the 150 Hz IMU default is only quoted.

    python scripts/angle_detector_risk.py            compute, write results/ and docs/
    python scripts/angle_detector_risk.py --check    recompute and compare with the committed files

Torque basis: T0 and the ramp rate are device side; G, the user's hold and the user's arrest torque
are arm side; arm torque from the device = k * T(t) (same as PR #9 scripts/ba15_numcheck.py).
Seeding: one SeedSequence, children by spawn (no seed lists; see /workspace/exo-failures.md
2026-11-22 seed-list bug).
"""

from __future__ import annotations

import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "angle_detector_risk.json"
TEMPLATE = ROOT / "scripts" / "angle_detector_risk.template.md"
DOC = ROOT / "docs" / "angle_detector_risk.md"

# ---------------------------------------------------------------- inputs (quoted, not re-derived)
# Drop model: PR #9 scripts/ba15_numcheck.py @ 58cce8d (Hardware Engineer's model; exo/ is not in repo).
ARM = dict(m_body=80.0, f_fa=0.016, f_hand=0.006, L_fa=0.27, com_fa=0.430, rg_fa=0.303, hand_off=0.06,
           m_pay=2.0, pay_off=0.08, m_dev=0.30, r_dev=0.12, m_carr=0.04, s_mid=0.175, J_rot=0.5e-5,
           gear=30, g=9.80665)
K_HI, T0_DEV = 1.5862286754583168, 6.304261267443454  # PR #9 (exo/benchC_rules_out.json revD.k_hi)
DROP_LIM = 10.0  # deg, B9
CATCH_MS = 200.0  # B-A15
# Governing cell of the B-A15 window (PR #9 section 1): (b) e-stop, ramp 21.00 N*m/s, STO 148 ms, hold 0
GOV = dict(t_det=0.0, rate=21.0, sto=0.148, hold=0.0)
WINDOW_LO_S = 0.20190  # decided lower edge (Charlie 1a, PR #9), s

# PR #9 measurement requirements (docs/ba15_numbers.json syn.hw.budget @ 58cce8d), quoted as inputs
FS_REQ = 15000.0  # Hz
BITS_REQ = 17
Q_REQ = 360.0 / 2**BITS_REQ  # deg
ENC_LAT_S = 0.019048484011547375e-3
TRIG_SKEW_S = 0.019048484011547375e-3
REL_SKEW_S = 0.03809696802309475e-3
U_T_S = 0.1599868492362151e-3  # proven timing uncertainty used in the score
GUARD_DEG = 0.026270778655533122  # angle uncertainty added to the measured peak drop
# PR #9 published model values, used only to check this re-implementation
PR9_REF = dict(t10_ms=201.90484840115474, drop200_deg=9.688776970044117, w10_dps=165.04629561071022,
               J=0.33313727020800005, G=10.297060953199999)

# Rev C package (exo/EMG_Exo_Package_RevC.md section 3.7 and BOM), quoted
H_B, S_MID = 0.054, 0.175  # pivot height and mid carriage position, m
RAIL = (0.100, 0.250)  # carriage window, m
E_SRC, E_CAP = 0.100, 0.0463  # offset bound and the bound that passes, m
REVC_JOINT_RES_DEG = 0.1  # Rev C BOM: joint encoder "absolute, <= 0.1 deg"
IMU_HZ_DEFAULT = 150.0  # repo default (protected), quoted only

# ---------------------------------------------------------------- ASSUMPTIONS (this session)
T_ARREST = (15.0, 20.0, 30.0, 50.0)  # user's constant arresting torque at the arm, N*m (swept)
T_ARREST_MAIN = 30.0
SEED = 20261126
N_MC = 400  # Monte Carlo trials per point
BASE_WIN = (-0.050, 0.0)  # baseline window before onset, s
T_END = 0.40  # end of the scored record, s
V_MOVING = 20.0  # deg/s: speed the arm must exceed before "speed zero" can mark arrest
A_LIM = 50000.0  # deg/s^2: plausibility limit on measured deceleration (about 7x the largest modelled)
MOTION_T, MOTION_FRAC = 0.100, 0.5  # motion check: measured drop at 100 ms >= half the model's
N_PLAN = (59, 330)  # PR #9: 0 in 59 (2b); planning default 33 x 10
FINE_DT = 2e-6  # s, grid for the true trajectory (linear interpolation error ~ a*dt^2/8 < 1e-8 deg)


# ---------------------------------------------------------------- drop model (SYNTHETIC)
def arm() -> dict:
    a = ARM
    parts = [(a["f_fa"] * a["m_body"], a["com_fa"] * a["L_fa"], a["rg_fa"] * a["L_fa"]),
             (a["f_hand"] * a["m_body"], a["L_fa"] + a["hand_off"], 0.0),
             (a["m_pay"], a["L_fa"] + a["pay_off"], 0.0)]
    J_h = sum(m * (r**2 + rg**2) for m, r, rg in parts)
    Mr_h = sum(m * r for m, r, _ in parts)
    J = J_h + a["m_dev"] * a["r_dev"] ** 2 + a["m_carr"] * a["s_mid"] ** 2 + a["J_rot"] * a["gear"] ** 2
    G = (Mr_h + a["m_dev"] * a["r_dev"] + a["m_carr"] * a["s_mid"]) * a["g"]
    return dict(J=J, G=G)


@dataclass
class Traj:
    """True trajectory. ``d`` is the drop in deg (positive down) on the fine grid ``t`` (s)."""
    t: np.ndarray
    d: np.ndarray
    t_c: float  # arrest complete (speed zero), s; inf if no catch
    peak: float  # true peak drop, deg
    t_on: float | None
    T_h: float | None
    rebound: bool
    t_peak: float  # time of the true peak, s

    def drop(self, ts):
        return np.interp(ts, self.t, self.d, left=0.0, right=self.d[-1])


def _solve(t_on=None, T_h=None, rebound=False, cell=GOV, end=0.45, dense=True):
    """Piecewise DOP853. x = [angle from horizontal (rad, negative = drop), speed (rad/s)].
    After arrest (speed back to zero) the arm is held at the peak, unless ``rebound``: then the
    arrest torque keeps acting and the arm rises again (worst case for a gap at the peak)."""
    A = arm()
    J, G = A["J"], A["G"]
    td, r, sto, hold = cell["t_det"], cell["rate"], cell["sto"], cell["hold"]

    def dev(t):
        return T0_DEV if t < td else max(T0_DEV - r * (t - td), 0.0)

    cuts = {0.0, end, td, td + T0_DEV / r, sto} | ({t_on} if t_on is not None else set())
    edges = sorted(e for e in cuts if 0.0 <= e <= end)
    y = np.array([0.0, 0.0])
    segs, t_c = [], math.inf

    def arrest(t, x):
        return x[1]

    arrest.terminal, arrest.direction = True, 1
    for a, b in zip(edges[:-1], edges[1:]):
        if b <= a:
            continue
        off = a >= sto
        hum = T_h if (t_on is not None and a >= t_on) else 0.0

        def f(t, x, off=off, hum=hum):
            return [x[1], (K_HI * (0.0 if off else dev(t)) + hold + hum - G * math.cos(x[0])) / J]

        ev = arrest if (hum and not rebound) else None
        s = solve_ivp(f, (a, b), y, method="DOP853", rtol=1e-11, atol=1e-13, dense_output=dense,
                      events=ev)
        segs.append((a, float(s.t[-1]), s.sol))
        y = s.y[:, -1]
        if ev is not None and s.t_events[0].size:
            t_c = float(s.t_events[0][0])
            y = s.y_events[0][0]
            break
    return segs, t_c, -math.degrees(y[0]), y  # peak is final only without rebound


def simulate(t_on=None, T_h=None, rebound=False, cell=GOV, end=0.45) -> Traj:
    segs, t_c, peak, y = _solve(t_on, T_h, rebound, cell, end)
    t = np.arange(0.0, end + FINE_DT / 2, FINE_DT)
    d = np.empty_like(t)
    filled = np.zeros(t.size, bool)
    for a, b, sol in segs:
        m = (t >= a) & (t <= b) & ~filled
        d[m] = -np.degrees(sol(t[m])[0])
        filled |= m
    d[~filled] = -math.degrees(y[0])  # held at the peak after arrest
    if rebound:
        peak = float(d.max())
    t_peak = t_c if math.isfinite(t_c) else float(t[np.argmax(d)])
    return Traj(t, d, t_c, float(peak), t_on, T_h, rebound, t_peak)


def free_fall_points(cell=GOV) -> dict:
    """10 deg crossing time and speed, drop at 200 ms (no catch); compare PR #9."""
    tr = simulate(cell=cell)
    i10 = int(np.argmax(tr.d >= DROP_LIM))
    t10 = float(np.interp(DROP_LIM, tr.d[i10 - 1:i10 + 1], tr.t[i10 - 1:i10 + 1]))
    w10 = float(np.interp(t10, tr.t, np.gradient(tr.d, FINE_DT)))
    return dict(t10_ms=t10 * 1e3, w10_dps=w10, drop200_deg=float(tr.drop(CATCH_MS / 1e3)),
                drop100_deg=float(tr.drop(MOTION_T)), **arm())


_FF_CACHE: dict = {}


def free_fall_drop(t):
    """Model drop (deg) at time t with no catch (governing cell): what a relaxed arm must show."""
    if "tr" not in _FF_CACHE:
        _FF_CACHE["tr"] = simulate()
    return float(_FF_CACHE["tr"].drop(t))


def t_on_for_peak(target, T_h, rebound=False, cell=GOV):
    def g(t_on):
        if rebound:
            return simulate(t_on, T_h, True, cell).peak - target
        return _solve(t_on, T_h, False, cell, dense=False)[2] - target
    return brentq(g, 0.05, 0.30, xtol=1e-12)


# ---------------------------------------------------------------- measurement chain (SYNTHETIC)
@dataclass
class Chain:
    fs: float = FS_REQ
    q: float = Q_REQ  # deg per count; 0 = no quantization
    lat: float = ENC_LAT_S  # angle-path latency, s
    skew: float = 0.0  # logged onset minus true onset, s (positive: times read early)
    skew_uniform: float = 0.0  # if > 0, skew ~ U(-skew_uniform, +skew_uniform) per trial
    sigma: float = 0.0  # white noise on the reading, deg
    gain: float = 0.0  # scale error (reading = (1+gain) * angle)
    zero: float = 0.0  # zero-reference error added to the baseline, deg (positive: drop read high)
    med3: bool = False  # 3-sample median before scoring
    k_dev: float | None = None  # read the device joint angle only (device angle = k * forearm angle)
    pivot: bool = False  # forearm = joint - pivot, two quantized channels (needs k_dev)
    stuck_at: float | None = None  # reading frozen from this time (s) on; nan = uniform over [0, t_c]
    gap: tuple | None = None  # (start s or "peak", duration s): samples lost
    glitch: float = 0.0  # one-sample spike (deg, drop direction) at a random post-onset sample


def _median3(x):
    y = x.copy()
    y[:, 1:-1] = np.median(np.stack([x[:, :-2], x[:, 1:-1], x[:, 2:]]), axis=0)
    return y


def measure(tr: Traj, ch: Chain, rng: np.random.Generator, n: int = N_MC, v_win: int | None = None):
    """Sample, quantize and score ``n`` trials of one true trajectory. Returns per-trial arrays."""
    nper = int(round((T_END - BASE_WIN[0]) * ch.fs))
    phase = rng.uniform(0.0, 1.0 / ch.fs, n)
    ts = BASE_WIN[0] + phase[:, None] + np.arange(nper)[None, :] / ch.fs  # true sample instants
    ang = -tr.drop(ts - ch.lat)  # forearm angle from start, deg (up positive)
    if ch.k_dev is not None:
        joint = ch.k_dev * ang
        reading = (1 + ch.gain) * joint
    else:
        reading = (1 + ch.gain) * ang
    noise = rng.normal(0.0, ch.sigma, ts.shape) if ch.sigma > 0 else None
    if noise is not None:
        reading = reading + noise
    mount = rng.uniform(0.0, ch.q, (n, 1)) if ch.q > 0 else None

    def quant(x, mnt):
        return np.floor((x + mnt) / ch.q) * ch.q if ch.q > 0 else x

    val = quant(reading, mount)
    if ch.pivot:  # pivot angle (forearm relative to link) = (k - 1) * forearm angle, own quantizer
        mount2 = rng.uniform(0.0, ch.q, (n, 1))
        piv = quant((1 + ch.gain) * (ch.k_dev - 1.0) * ang, mount2)
        val = val - piv
    skew = rng.uniform(-ch.skew_uniform, ch.skew_uniform, n) if ch.skew_uniform > 0 else np.full(n, ch.skew)
    tm = ts - skew[:, None]  # logged time from logged onset
    if ch.stuck_at is not None:
        t_s = (rng.uniform(0.0, tr.t_c, n) if math.isnan(ch.stuck_at) else np.full(n, ch.stuck_at))
        frozen = ts >= t_s[:, None]
        first = np.argmax(frozen, axis=1)
        hold_v = val[np.arange(n), np.maximum(first - 1, 0)]
        val = np.where(frozen, hold_v[:, None], val)
    lost = np.zeros(ts.shape, bool)
    if ch.gap is not None:
        g0, dur = ch.gap
        g0 = tr.t_peak - dur / 2 if g0 == "peak" else g0
        lost = (ts >= g0) & (ts < g0 + dur)
    base_m = (tm >= BASE_WIN[0]) & (tm < BASE_WIN[1])
    base = np.nanmean(np.where(base_m, val, np.nan), axis=1) + ch.zero
    dm = base[:, None] - val
    if ch.glitch:
        post_idx = np.argmax(tm >= 0, axis=1)
        j = post_idx + rng.integers(0, int(0.3 * ch.fs), n)
        dm[np.arange(n), j] += ch.glitch
    dm = np.where(lost, np.nan, dm)
    if ch.med3:
        dm = _median3(dm)
    post = tm >= 0
    dpost = np.where(post, dm, -np.inf)
    dpost = np.where(np.isnan(dpost), -np.inf, dpost)
    peak = dpost.max(axis=1)
    # 10 deg crossing: first sample at or above 10 deg (no interpolation), and linear interpolation
    hit = dpost >= DROP_LIM
    has = hit.any(axis=1)
    i10 = np.argmax(hit, axis=1)
    rows = np.arange(n)
    t10 = np.where(has, tm[rows, i10], np.nan)
    i0 = np.maximum(i10 - 1, 0)
    d0, d1 = dm[rows, i0], dm[rows, i10]
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(d1 > d0, (DROP_LIM - d0) / (d1 - d0), 1.0)
    t10i = np.where(has, tm[rows, i0] + frac * (tm[rows, i10] - tm[rows, i0]), np.nan)
    # drop at 200 ms: last sample logged at or before 200 ms
    i200 = np.sum(tm <= CATCH_MS / 1e3, axis=1) - 1
    d200 = dm[rows, i200]
    # catch time, argmax rule: first logged time at which the measured drop reaches its maximum
    iarg = np.argmax(dpost >= peak[:, None], axis=1)
    tc_arg = tm[rows, iarg]
    out = dict(peak=peak, t10=t10, t10i=t10i, d200=d200, tc_arg=tc_arg, skew=skew, phase=phase)
    if v_win is not None:  # catch time, speed-zero rule with a W-sample difference
        W = v_win
        v = (dm[:, W:] - dm[:, :-W]) * ch.fs / W
        tv = 0.5 * (tm[:, W:] + tm[:, :-W])
        moving = np.maximum.accumulate(np.where(np.isnan(v), -np.inf, v), axis=1) > V_MOVING
        stop = moving & (v <= 0) & (tv > 0)
        has_s = stop.any(axis=1)
        out["tc_v"] = np.where(has_s, tv[rows, np.argmax(stop, axis=1)], np.nan)
    # plausibility checks (INVALID flags)
    i100 = np.sum(tm <= MOTION_T, axis=1) - 1
    out["motion_bad"] = dm[rows, i100] < MOTION_FRAC * free_fall_drop(MOTION_T)
    Wc = max(1, int(round(1e-3 * ch.fs)))  # 1 ms speed windows
    vc = (dm[:, Wc:] - dm[:, :-Wc]) * ch.fs / Wc
    dec = -(vc[:, Wc:] - vc[:, :-Wc]) * ch.fs / Wc
    out["max_decel"] = np.nanmax(np.where(np.isnan(dec), -np.inf, dec), axis=1)
    out["decel_bad"] = out["max_decel"] > A_LIM
    out["gap_bad"] = lost.any(axis=1)
    return out


def score_angle(peak):
    return peak + GUARD_DEG > DROP_LIM


# ---------------------------------------------------------------- P+R interface (Rev C 3.7)
L_LINK = math.hypot(S_MID, H_B)


def device_angle(phi, ex, ey):
    """Device link angle (rad) for forearm angle ``phi`` (rad), offset (ex, ey) given in the forearm
    frame at phi = 0 (fixed in the upper-arm frame). Also returns s_B and d."""
    c, s = math.cos(-phi), math.sin(-phi)
    exp, eyp = c * ex - s * ey, s * ex + c * ey  # offset seen in the rotated forearm frame
    d = math.sqrt(L_LINK**2 - (H_B - eyp) ** 2)
    alpha = math.atan2(H_B - eyp, d)
    return phi + alpha, exp + d, d


def k_ratio(ex, ey, h=1e-6):
    """k = s_B / d (Rev C) and the numerical d(psi)/d(phi), for the check k = dpsi/dphi."""
    _, sB, d = device_angle(0.0, ex, ey)
    dpsi = (device_angle(h, ex, ey)[0] - device_angle(-h, ex, ey)[0]) / (2 * h)
    return sB / d, dpsi, sB


def pr_sweep(e_max, n_r=41, n_a=721, clip=False):
    """Disc sweep (Rev C: 41 radii x 721 angles). For each offset: k, device change for a 10 deg
    forearm drop, and the forearm drop when the device reads 10 deg."""
    ks, dev10, hum10 = [], [], []
    for rr in np.linspace(0.0, e_max, n_r):
        for th in np.linspace(0.0, 2 * math.pi, n_a):
            ex, ey = rr * math.cos(th), rr * math.sin(th)
            k, _, sB = k_ratio(ex, ey)
            if clip and not (RAIL[0] <= sB <= RAIL[1]):
                continue
            psi0 = device_angle(0.0, ex, ey)[0]
            dv = math.degrees(psi0 - device_angle(-math.radians(DROP_LIM), ex, ey)[0])
            hd = brentq(lambda ph: math.degrees(psi0 - device_angle(-math.radians(ph), ex, ey)[0])
                        - DROP_LIM, 0.5, 60.0, xtol=1e-12)
            ks.append(k)
            dev10.append(dv)
            hum10.append(hd)
    return dict(k_min=min(ks), k_max=max(ks), dev_for_10_min=min(dev10), dev_for_10_max=max(dev10),
                hum_at_dev10_min=min(hum10), hum_at_dev10_max=max(hum10), n=len(ks))


# ---------------------------------------------------------------- analyses
def _stats(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return dict(min=float(x.min()), max=float(x.max()), mean=float(x.mean()), sd=float(x.std(ddof=1)))


def compute() -> dict:
    ss = np.random.SeedSequence(SEED)

    def rng():  # a new, distinct child stream per call (spawn keys, no seed lists)
        return np.random.default_rng(ss.spawn(1)[0])

    R: dict = {"inputs": dict(fs_req=FS_REQ, bits_req=BITS_REQ, q_req=Q_REQ, enc_lat_ms=ENC_LAT_S * 1e3,
                              trig_skew_ms=TRIG_SKEW_S * 1e3, rel_skew_ms=REL_SKEW_S * 1e3,
                              u_t_ms=U_T_S * 1e3, guard=GUARD_DEG, n_mc=N_MC, t_arrest=list(T_ARREST),
                              t_arrest_main=T_ARREST_MAIN, a_lim=A_LIM, seed=SEED)}
    A = arm()
    # 1. model check against PR #9
    ff = free_fall_points()
    R["model"] = dict(ff, ref=PR9_REF, d_t10_ms=ff["t10_ms"] - PR9_REF["t10_ms"],
                      d_drop200=ff["drop200_deg"] - PR9_REF["drop200_deg"])
    trf = simulate()
    skew_max = TRIG_SKEW_S + REL_SKEW_S  # ASSUMPTION: release skew acts as a timing offset like trigger skew

    # 2. crossing-time channel (Phase A per-trial rule, dummy matrix) at the PR #9 requirements
    cross = {}
    for name, ch in dict(req=Chain(skew_uniform=skew_max),
                         req_skew_pos=Chain(skew=skew_max), req_skew_neg=Chain(skew=-skew_max),
                         revc_joint=Chain(q=REVC_JOINT_RES_DEG, skew_uniform=skew_max),
                         imu150=Chain(fs=IMU_HZ_DEFAULT, skew_uniform=skew_max)).items():
        m = measure(trf, ch, rng())
        e10 = (m["t10"] - ff["t10_ms"] / 1e3) * 1e3
        e10i = (m["t10i"] - ff["t10_ms"] / 1e3) * 1e3
        e200 = m["d200"] - ff["drop200_deg"]
        cross[name] = dict(t10_err_ms=_stats(e10), t10i_err_ms=_stats(e10i), d200_err_deg=_stats(e200),
                           t10_abs_max_ms=float(np.nanmax(np.abs(e10))),
                           d200_abs_max_deg=float(np.nanmax(np.abs(e200))),
                           within_u_t=bool(np.nanmax(np.abs(e10)) <= U_T_S * 1e3),
                           within_guard=bool(np.nanmax(np.abs(e200)) <= GUARD_DEG))
    R["crossing"] = cross

    # 3. peak-drop channel: operating characteristic around 10 deg
    deltas = [-0.05, -0.03, -0.0263, -0.02, -0.01, -0.005, 0.0, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05]
    trajs = {}
    for Th in T_ARREST:
        for dl in deltas:
            t_on = t_on_for_peak(DROP_LIM + dl, Th)
            trajs[(Th, dl)] = simulate(t_on, Th)
    oc = {}
    chains = dict(req=Chain(skew_uniform=skew_max),
                  noise_1lsb=Chain(sigma=Q_REQ, skew_uniform=skew_max),
                  noise_4lsb=Chain(sigma=4 * Q_REQ, skew_uniform=skew_max),
                  noise_4lsb_med3=Chain(sigma=4 * Q_REQ, med3=True, skew_uniform=skew_max),
                  revc_joint=Chain(q=REVC_JOINT_RES_DEG, skew_uniform=skew_max),
                  imu150=Chain(fs=IMU_HZ_DEFAULT, skew_uniform=skew_max))
    for cname, ch in chains.items():
        rows = []
        for dl in deltas:
            m = measure(trajs[(T_ARREST_MAIN, dl)], ch, rng())
            err = m["peak"] - trajs[(T_ARREST_MAIN, dl)].peak
            rows.append(dict(delta=dl, p_drop=float(score_angle(m["peak"]).mean()), err=_stats(err)))
        miss = [r["p_drop"] for r in rows if r["delta"] > 0]
        oc[cname] = dict(rows=rows, max_missed=float(1 - min(miss)),
                         err_min=min(r["err"]["min"] for r in rows),
                         err_max=max(r["err"]["max"] for r in rows))
    for cname, v in oc.items():  # false-trip band: lowest true margin with any DROP, and with all DROP
        v["trip_any_delta"] = min(r["delta"] for r in v["rows"] if r["p_drop"] > 0)
        v["trip_all_delta"] = min(
            r["delta"] for r in v["rows"]
            if all(r2["p_drop"] == 1.0 for r2 in v["rows"] if r2["delta"] >= r["delta"]))
    R["peak_oc"] = oc
    # peak error over arrest torques at the requirement chain
    R["peak_err_by_Th"] = {f"{Th:g}": _stats(np.concatenate([
        measure(trajs[(Th, dl)], Chain(skew_uniform=skew_max), rng(), n=200)["peak"] - trajs[(Th, dl)].peak
        for dl in (-0.01, 0.0, 0.01)])) for Th in T_ARREST}

    # 4. single-source limits for a missed detection (true peak just above 10 deg)
    tr_edge = trajs[(T_ARREST_MAIN, 0.002)]

    crn = ss.spawn(1)[0]  # common random numbers: the same draws at every bisection step

    def worst_peak(ch, n=200):
        return float(measure(tr_edge, ch, np.random.default_rng(crn), n=n)["peak"].min())

    lim = {}
    lim["gain"] = brentq(lambda g: worst_peak(Chain(gain=g)) + GUARD_DEG - DROP_LIM, -0.02, 0.0, xtol=1e-7)
    lim["zero"] = brentq(lambda z: worst_peak(Chain(zero=z)) + GUARD_DEG - DROP_LIM, -0.1, 0.0, xtol=1e-7)
    # resolution. Without noise the measured drop is a whole number of counts, so whether a given true
    # drop is missed depends on where 10 deg falls on the count lattice. To avoid that artefact, the true
    # peak is spread over 10 + (0, 0.05] deg and the worst miss rate is kept, without noise and with
    # 0.5-count noise (dither). The bit counts at which a miss is still seen are reported.
    edge_trajs = [simulate(t_on_for_peak(DROP_LIM + dl, T_ARREST_MAIN), T_ARREST_MAIN)
                  for dl in np.linspace(0.002, 0.05, 13)]
    scan = []
    for bits in (10, 11, 12, 13, 14, 15, 16, 17):
        q = 360.0 / 2**bits
        row = dict(bits=bits, q=q)
        for nm, sg in (("plain", 0.0), ("dither", 0.5 * q)):
            row[f"miss_max_{nm}"] = max(
                float((~score_angle(measure(tr, Chain(q=q, sigma=sg, skew_uniform=skew_max), rng(),
                                            n=100)["peak"])).mean()) for tr in edge_trajs)
        scan.append(row)
    lim["q_scan"] = scan
    lim["q_bits_last_miss"] = max(r["bits"] for r in scan
                                  if r["miss_max_plain"] > 0 or r["miss_max_dither"] > 0)
    lim["q_bits_guard"] = math.ceil(math.log2(360.0 / GUARD_DEG))  # smallest bits with q < guard
    lim["analytic"] = dict(gain=-(GUARD_DEG + 0.002) / (DROP_LIM + 0.002),
                           zero=-(GUARD_DEG + 0.002))
    R["miss_limits"] = lim

    # 5. catch-time channel (speed zero / argmax) vs U_t
    tc = {}
    for Th in T_ARREST:
        tr = trajs[(Th, 0.0)]
        A_dec = (Th - A["G"] * math.cos(math.radians(DROP_LIM))) / A["J"]  # rad/s^2 near the peak
        row = dict(t_c_true_ms=tr.t_c * 1e3, decel_dps2=math.degrees(A_dec),
                   flat_ms_1lsb=math.sqrt(2 * Q_REQ / math.degrees(A_dec)) * 1e3)
        for W in (1, 15, 75):
            m = measure(tr, Chain(skew_uniform=skew_max), rng(), n=200, v_win=W)
            row[f"tc_v_W{W}_err_ms"] = _stats((m["tc_v"] - tr.t_c) * 1e3)
        m = measure(tr, Chain(skew_uniform=skew_max), rng(), n=200)
        row["tc_arg_err_ms"] = _stats((m["tc_arg"] - tr.t_c) * 1e3)
        m = measure(tr, Chain(sigma=Q_REQ, skew_uniform=skew_max), rng(), n=200, v_win=15)
        row["tc_v_W15_noise1_err_ms"] = _stats((m["tc_v"] - tr.t_c) * 1e3)
        row["tc_arg_noise1_err_ms"] = _stats((m["tc_arg"] - tr.t_c) * 1e3)
        # time rule (t_c + U_t > 201.90 ms) vs angle truth: the largest true peak drop whose arrest
        # completes just inside the window edge; every slower-to-arrest trial is a time-rule DROP
        t_edge = WINDOW_LO_S - U_T_S
        t_on_e = brentq(lambda to: _solve(to, Th, dense=False)[1] - t_edge, 0.02, 0.2019, xtol=1e-12)
        row["peak_at_time_edge_deg"] = _solve(t_on_e, Th, dense=False)[2]
        row["t_on_at_time_edge_ms"] = t_on_e * 1e3
        row["t_on_for_10_ms"] = tr.t_on * 1e3
        # angle error needed for the argmax time to stay within U_t: 0.5 * a * U_t^2
        row["ang_for_u_t_deg"] = 0.5 * math.degrees(A_dec) * U_T_S**2
        row["bits_for_u_t"] = math.ceil(math.log2(360.0 / row["ang_for_u_t_deg"]))
        tc[f"{Th:g}"] = row
    R["catch_time"] = tc

    # 6. P+R interface: which angle is measured
    pr = dict(e_src=pr_sweep(E_SRC), e_src_clip=pr_sweep(E_SRC, clip=True), e_cap=pr_sweep(E_CAP),
              k_check=k_ratio(0.03, -0.02)[:2])
    # peak channel when the score reads the device joint only, at the e_cap extremes (linear k)
    for nm, kk in (("k_lo", pr["e_cap"]["k_min"]), ("k_hi", pr["e_cap"]["k_max"])):
        m = measure(trajs[(T_ARREST_MAIN, 0.0)], Chain(k_dev=kk, skew_uniform=skew_max), rng(), n=200)
        pr[f"joint_only_{nm}_peak_err"] = _stats(m["peak"] - DROP_LIM)
        m = measure(trajs[(T_ARREST_MAIN, 0.0)], Chain(k_dev=kk, pivot=True, skew_uniform=skew_max), rng(),
                    n=200)
        pr[f"joint_pivot_{nm}_peak_err"] = _stats(m["peak"] - DROP_LIM)
    R["pr"] = pr

    # 7. encoder faults and plausibility checks (true peak 10.01 deg: a real drop past 10)
    tr_f = trajs[(T_ARREST_MAIN, 0.01)]
    tr_reb = simulate(t_on_for_peak(DROP_LIM + 0.01, T_ARREST_MAIN, rebound=True), T_ARREST_MAIN,
                      rebound=True)
    faults = {}

    def fault(name, ch, tr=tr_f):
        m = measure(tr, ch, rng())
        drop = score_angle(m["peak"])
        inval = m["motion_bad"] | m["decel_bad"] | m["gap_bad"]
        faults[name] = dict(missed=float((~drop).mean()), missed_after_checks=float((~drop & ~inval).mean()),
                            flagged=float(inval.mean()), motion_bad=float(m["motion_bad"].mean()),
                            decel_bad=float(m["decel_bad"].mean()), gap_bad=float(m["gap_bad"].mean()))
        return m

    fault("clean", Chain(skew_uniform=skew_max))
    fault("clean_noise1", Chain(sigma=Q_REQ, skew_uniform=skew_max))
    fault("stuck_from_onset", Chain(stuck_at=0.0, skew_uniform=skew_max))
    st = fault("stuck_during_drop", Chain(stuck_at=float("nan"), skew_uniform=skew_max))
    fault("gap_2ms_at_peak_held", Chain(gap=("peak", 2e-3), skew_uniform=skew_max))
    fault("gap_2ms_at_peak_rebound", Chain(gap=("peak", 2e-3), skew_uniform=skew_max), tr=tr_reb)
    fault("gain_minus_0p5pct", Chain(gain=-0.005, skew_uniform=skew_max))
    # false trips from a one-sample glitch on a true 9.95 deg arm
    tr_ok = trajs[(T_ARREST_MAIN, -0.05)]
    for nm, ch in (("glitch_0p1", Chain(glitch=0.1, skew_uniform=skew_max)),
                   ("glitch_0p1_med3", Chain(glitch=0.1, med3=True, skew_uniform=skew_max))):
        m = measure(tr_ok, ch, rng())
        faults[nm] = dict(false_trip=float(score_angle(m["peak"]).mean()))
    # clean trials must not be flagged: check over all arrest torques at delta 0
    flags = [measure(trajs[(Th, 0.0)], Chain(sigma=Q_REQ, skew_uniform=skew_max), rng(), n=200)
             for Th in T_ARREST]
    faults["clean_flag_rate_all_Th"] = float(np.mean([(f["motion_bad"] | f["decel_bad"]).mean()
                                                      for f in flags]))
    faults["clean_max_decel_all_Th"] = float(max(f["max_decel"].max() for f in flags))
    faults["model_max_decel"] = max(math.degrees((Th - A["G"] * math.cos(math.radians(DROP_LIM))) / A["J"])
                                    for Th in T_ARREST)
    # decel threshold sweep: stuck-during-drop misses after checks vs clean false flags (noise 1 LSB)
    st_drop = score_angle(st["peak"])
    sweep = []
    for a_lim in (10000.0, 20000.0, 50000.0, 100000.0):
        miss = ~st_drop & ~(st["motion_bad"] | (st["max_decel"] > a_lim) | st["gap_bad"])
        clean = np.mean([(f["motion_bad"] | (f["max_decel"] > a_lim)).mean() for f in flags])
        sweep.append(dict(a_lim=a_lim, stuck_missed_after_checks=float(miss.mean()),
                          clean_flagged=float(clean)))
    faults["decel_sweep"] = sweep
    # rebound gap: under-read vs gap length (analytic 0.5*a*(D/2)^2 check)
    a_reb = math.degrees((T_ARREST_MAIN - A["G"] * math.cos(math.radians(DROP_LIM))) / A["J"])
    faults["gap_rebound_underread_2ms_analytic"] = 0.5 * a_reb * (1e-3) ** 2
    faults["gap_rebound_max_ms_for_guard"] = 2 * math.sqrt(2 * GUARD_DEG / a_reb) * 1e3
    R["faults"] = faults

    # 8. plan level: per-trial false-trip p -> P(at least one false DROP) in n trials (exact)
    plan = {}
    for p in (1e-4, 1e-3, 1e-2, 0.05):
        plan[f"p{p:.0e}"] = dict(p=p, **{f"n{n}": 1 - (1 - p) ** n for n in N_PLAN})
    # largest per-trial false-trip p that keeps P(any false DROP) <= 5 % over n trials
    plan["p_for_5pct"] = {f"n{n}": 1 - 0.95 ** (1 / n) for n in N_PLAN}
    R["plan"] = plan
    derive(R)
    return R


def derive(R) -> None:
    """Summary values quoted in the doc (all computed from the results above)."""
    R["inputs"]["skew_total_ms"] = R["inputs"]["trig_skew_ms"] + R["inputs"]["rel_skew_ms"]
    R["model"]["margin_ms"] = R["model"]["t10_ms"] - CATCH_MS
    req = [R["crossing"][k] for k in ("req", "req_skew_pos", "req_skew_neg")]
    c = R["crossing"]
    c["worst_t10_abs_max_ms"] = max(r["t10_abs_max_ms"] for r in req)
    c["worst_d200_abs_max_deg"] = max(r["d200_abs_max_deg"] for r in req)
    c["u_t_room_ms"] = R["inputs"]["u_t_ms"] - c["worst_t10_abs_max_ms"]
    c["guard_room_deg"] = R["inputs"]["guard"] - c["worst_d200_abs_max_deg"]
    for v in R["peak_oc"].values():
        v["trip_all_delta_abs_deg"] = -v["trip_all_delta"]
    pe = R["peak_err_by_Th"]
    pe["worst_min"] = min(v["min"] for v in pe.values())
    pe["worst_max"] = max(v["max"] for v in pe.values() if isinstance(v, dict))
    pr = R["pr"]
    pr["joint_pivot_abs_max"] = max(abs(pr[f"joint_pivot_{k}_peak_err"][m]) for k in ("k_lo", "k_hi")
                                    for m in ("min", "max"))
    rows = [v for v in R["catch_time"].values() if isinstance(v, dict)]

    def ext(key):
        return max(max(abs(r[key]["min"]), abs(r[key]["max"])) for r in rows)

    ct = R["catch_time"]
    ct["worst_arg_ms"] = ext("tc_arg_err_ms")
    ct["worst_v15_ms"] = ext("tc_v_W15_err_ms")
    ct["worst_noise_ms"] = max(ext("tc_arg_noise1_err_ms"), ext("tc_v_W15_noise1_err_ms"))
    ct["worst_noise_arg_ms"] = max(r["tc_arg_noise1_err_ms"]["max"] for r in rows)


def _f(x, fmt):
    return format(x, fmt)


def tables(R) -> dict:
    """Markdown table bodies for the doc, built only from R."""
    oc = R["peak_oc"]
    names = ("req", "noise_1lsb", "noise_4lsb", "noise_4lsb_med3", "revc_joint", "imu150")
    lab = dict(req="req", noise_1lsb="1 count noise", noise_4lsb="4 counts noise",
               noise_4lsb_med3="4 counts + median-3", revc_joint="0.1° steps", imu150="150 Hz")
    deltas = [r["delta"] for r in oc["req"]["rows"]]
    oc_rows = []
    for i, dl in enumerate(deltas):
        cells = [_f(oc[n]["rows"][i]["p_drop"], ".3f") for n in names]
        oc_rows.append(f"| {dl:+.4f} | " + " | ".join(cells) + " |")
    summ = [f"| {lab[n]} | {_f(oc[n]['err_min'], '+.4f')} / {_f(oc[n]['err_max'], '+.4f')} | "
            f"{_f(oc[n]['max_missed'], '.3f')} | {oc[n]['trip_any_delta']:+.4f} | "
            f"{oc[n]['trip_all_delta']:+.4f} |"
            for n in names]
    q_rows = [f"| {r['bits']} | {r['q']:.5f} | {r['miss_max_plain']:.2f} | {r['miss_max_dither']:.2f} |"
              for r in R["miss_limits"]["q_scan"]]

    def rng_(d):
        return f"{d['min']:+.2f} to {d['max']:+.2f}"

    tc_rows = []
    for th in ("15", "20", "30", "50"):
        r = R["catch_time"][th]
        tc_rows.append(f"| {th} | {r['peak_at_time_edge_deg']:.3f} | {rng_(r['tc_arg_err_ms'])} | "
                       f"{rng_(r['tc_v_W15_err_ms'])} | {rng_(r['tc_v_W75_err_ms'])} | "
                       f"{rng_(r['tc_v_W1_err_ms'])} | "
                       f"{rng_(r['tc_arg_noise1_err_ms'])} | {rng_(r['tc_v_W15_noise1_err_ms'])} | "
                       f"{r['flat_ms_1lsb']:.3f} | {r['bits_for_u_t']} |")
    dec_rows = [f"| {r['a_lim']:.0f} | {r['stuck_missed_after_checks']:.1%} | {r['clean_flagged']:.1%} |"
                for r in R["faults"]["decel_sweep"]]
    return dict(oc_table="\n".join(oc_rows), oc_summary="\n".join(summ), q_table="\n".join(q_rows),
                tc_table="\n".join(tc_rows), decel_table="\n".join(dec_rows))


# ---------------------------------------------------------------- doc rendering
PH = re.compile(r"\{\{([^{}|]+)(?:\|([^{}]*))?\}\}")


def lookup(R, key):
    v = R
    for part in key.split("."):
        v = v[int(part)] if isinstance(v, list) else v[part]
    return v


def render(R, text):
    ctx = dict(R, **tables(R))

    def sub(m):
        v = lookup(ctx, m.group(1).strip())
        fmt = m.group(2)
        return format(v, fmt) if fmt else str(v)
    return PH.sub(sub, text)


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    return x


def main(argv) -> int:
    R = _jsonable(compute())
    text = json.dumps(R, indent=1, sort_keys=True) + "\n"
    doc = render(R, TEMPLATE.read_text())
    if "--check" in argv:
        bad = []
        if RESULTS.read_text() != text:
            bad.append(str(RESULTS))
        if DOC.read_text() != doc:
            bad.append(str(DOC))
        print("MISMATCH: " + ", ".join(bad) if bad else "OK: results and doc reproduce")
        return 1 if bad else 0
    RESULTS.parent.mkdir(exist_ok=True)
    RESULTS.write_text(text)
    DOC.write_text(doc)
    print(f"wrote {RESULTS.relative_to(ROOT)} and {DOC.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
