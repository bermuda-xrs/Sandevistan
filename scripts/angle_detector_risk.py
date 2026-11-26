"""Angle-detector risk for the B-A15 / B9 "drop past 10 deg" score. ALL NUMBERS ARE SYNTHETIC.

The "angle detector" is the joint-angle chain that turns a physical arm drop into the B-A15 trial
score (PR #9, docs/ba15-testplan.md @ 335d4e5, section 3 and "Trial score"):
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
# Drop model: PR #9 scripts/ba15_numcheck.py @ 335d4e5 (inputs identical at 58cce8d; round 1 cited that).
# It is Hardware Engineer's model; exo/ is not in this repo.
ARM = dict(
    m_body=80.0,
    f_fa=0.016,
    f_hand=0.006,
    L_fa=0.27,
    com_fa=0.430,
    rg_fa=0.303,
    hand_off=0.06,
    m_pay=2.0,
    pay_off=0.08,
    m_dev=0.30,
    r_dev=0.12,
    m_carr=0.04,
    s_mid=0.175,
    J_rot=0.5e-5,
    gear=30,
    g=9.80665,
)
K_HI, T0_DEV = 1.5862286754583168, 6.304261267443454  # PR #9 (exo/benchC_rules_out.json revD.k_hi)
DROP_LIM = 10.0  # deg, B9
CATCH_MS = 200.0  # B-A15
# Governing cell of the B-A15 window (PR #9 section 1): (b) e-stop, ramp 21.00 N*m/s, STO 148 ms, hold 0
GOV = dict(t_det=0.0, rate=21.0, sto=0.148, hold=0.0)
WINDOW_LO_S = 0.20190  # decided lower edge (Charlie 1a, PR #9), s

# PR #9 measurement requirements (docs/ba15_numbers.json syn.hw.budget @ 335d4e5, identical at 58cce8d),
# quoted as inputs
FS_REQ = 15000.0  # Hz
BITS_REQ = 17
Q_REQ = 360.0 / 2**BITS_REQ  # deg
ENC_LAT_S = 0.019048484011547375e-3
TRIG_SKEW_S = 0.019048484011547375e-3
REL_SKEW_S = 0.03809696802309475e-3
U_T_S = 0.1599868492362151e-3  # proven timing uncertainty used in the score
GUARD_DEG = 0.026270778655533122  # angle uncertainty added to the measured peak drop
# PR #9 published model values, used only to check this re-implementation
PR9_REF = dict(
    t10_ms=201.90484840115474,
    drop200_deg=9.688776970044117,
    w10_dps=165.04629561071022,
    J=0.33313727020800005,
    G=10.297060953199999,
)

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
    parts = [
        (a["f_fa"] * a["m_body"], a["com_fa"] * a["L_fa"], a["rg_fa"] * a["L_fa"]),
        (a["f_hand"] * a["m_body"], a["L_fa"] + a["hand_off"], 0.0),
        (a["m_pay"], a["L_fa"] + a["pay_off"], 0.0),
    ]
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
        s = solve_ivp(f, (a, b), y, method="DOP853", rtol=1e-11, atol=1e-13, dense_output=dense, events=ev)
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
    t10 = float(np.interp(DROP_LIM, tr.d[i10 - 1 : i10 + 1], tr.t[i10 - 1 : i10 + 1]))
    w10 = float(np.interp(t10, tr.t, np.gradient(tr.d, FINE_DT)))
    return dict(
        t10_ms=t10 * 1e3,
        w10_dps=w10,
        drop200_deg=float(tr.drop(CATCH_MS / 1e3)),
        drop100_deg=float(tr.drop(MOTION_T)),
        **arm(),
    )


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
        t_s = rng.uniform(0.0, tr.t_c, n) if math.isnan(ch.stuck_at) else np.full(n, ch.stuck_at)
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
            hd = brentq(
                lambda ph: math.degrees(psi0 - device_angle(-math.radians(ph), ex, ey)[0]) - DROP_LIM,
                0.5,
                60.0,
                xtol=1e-12,
            )
            ks.append(k)
            dev10.append(dv)
            hum10.append(hd)
    return dict(
        k_min=min(ks),
        k_max=max(ks),
        dev_for_10_min=min(dev10),
        dev_for_10_max=max(dev10),
        hum_at_dev10_min=min(hum10),
        hum_at_dev10_max=max(hum10),
        n=len(ks),
    )


# ---------------------------------------------------------------- analyses
def _stats(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return dict(min=float(x.min()), max=float(x.max()), mean=float(x.mean()), sd=float(x.std(ddof=1)))


def compute() -> dict:
    ss = np.random.SeedSequence(SEED)

    def rng():  # a new, distinct child stream per call (spawn keys, no seed lists)
        return np.random.default_rng(ss.spawn(1)[0])

    R: dict = {
        "inputs": dict(
            fs_req=FS_REQ,
            bits_req=BITS_REQ,
            q_req=Q_REQ,
            enc_lat_ms=ENC_LAT_S * 1e3,
            trig_skew_ms=TRIG_SKEW_S * 1e3,
            rel_skew_ms=REL_SKEW_S * 1e3,
            u_t_ms=U_T_S * 1e3,
            guard=GUARD_DEG,
            n_mc=N_MC,
            t_arrest=list(T_ARREST),
            t_arrest_main=T_ARREST_MAIN,
            a_lim=A_LIM,
            seed=SEED,
        )
    }
    A = arm()
    # 1. model check against PR #9
    ff = free_fall_points()
    R["model"] = dict(
        ff,
        ref=PR9_REF,
        d_t10_ms=ff["t10_ms"] - PR9_REF["t10_ms"],
        d_drop200=ff["drop200_deg"] - PR9_REF["drop200_deg"],
    )
    trf = simulate()
    skew_max = TRIG_SKEW_S + REL_SKEW_S  # ASSUMPTION: release skew acts as a timing offset like trigger skew

    # 2. crossing-time channel (Phase A per-trial rule, dummy matrix) at the PR #9 requirements
    cross = {}
    for name, ch in dict(
        req=Chain(skew_uniform=skew_max),
        req_skew_pos=Chain(skew=skew_max),
        req_skew_neg=Chain(skew=-skew_max),
        revc_joint=Chain(q=REVC_JOINT_RES_DEG, skew_uniform=skew_max),
        imu150=Chain(fs=IMU_HZ_DEFAULT, skew_uniform=skew_max),
    ).items():
        m = measure(trf, ch, rng())
        e10 = (m["t10"] - ff["t10_ms"] / 1e3) * 1e3
        e10i = (m["t10i"] - ff["t10_ms"] / 1e3) * 1e3
        e200 = m["d200"] - ff["drop200_deg"]
        cross[name] = dict(
            t10_err_ms=_stats(e10),
            t10i_err_ms=_stats(e10i),
            d200_err_deg=_stats(e200),
            t10_abs_max_ms=float(np.nanmax(np.abs(e10))),
            d200_abs_max_deg=float(np.nanmax(np.abs(e200))),
            within_u_t=bool(np.nanmax(np.abs(e10)) <= U_T_S * 1e3),
            within_guard=bool(np.nanmax(np.abs(e200)) <= GUARD_DEG),
        )
    R["crossing"] = cross

    # 3. peak-drop channel: operating characteristic around 10 deg
    deltas = [-0.05, -0.03, -0.0263, -0.02, -0.01, -0.005, 0.0, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05]
    trajs = {}
    for Th in T_ARREST:
        for dl in deltas:
            t_on = t_on_for_peak(DROP_LIM + dl, Th)
            trajs[(Th, dl)] = simulate(t_on, Th)
    oc = {}
    chains = dict(
        req=Chain(skew_uniform=skew_max),
        noise_1lsb=Chain(sigma=Q_REQ, skew_uniform=skew_max),
        noise_4lsb=Chain(sigma=4 * Q_REQ, skew_uniform=skew_max),
        noise_4lsb_med3=Chain(sigma=4 * Q_REQ, med3=True, skew_uniform=skew_max),
        revc_joint=Chain(q=REVC_JOINT_RES_DEG, skew_uniform=skew_max),
        imu150=Chain(fs=IMU_HZ_DEFAULT, skew_uniform=skew_max),
    )
    for cname, ch in chains.items():
        rows = []
        for dl in deltas:
            m = measure(trajs[(T_ARREST_MAIN, dl)], ch, rng())
            err = m["peak"] - trajs[(T_ARREST_MAIN, dl)].peak
            rows.append(dict(delta=dl, p_drop=float(score_angle(m["peak"]).mean()), err=_stats(err)))
        miss = [r["p_drop"] for r in rows if r["delta"] > 0]
        oc[cname] = dict(
            rows=rows,
            max_missed=float(1 - min(miss)),
            err_min=min(r["err"]["min"] for r in rows),
            err_max=max(r["err"]["max"] for r in rows),
        )
    for cname, v in oc.items():  # false-trip band: lowest true margin with any DROP, and with all DROP
        v["trip_any_delta"] = min(r["delta"] for r in v["rows"] if r["p_drop"] > 0)
        v["trip_all_delta"] = min(
            r["delta"]
            for r in v["rows"]
            if all(r2["p_drop"] == 1.0 for r2 in v["rows"] if r2["delta"] >= r["delta"])
        )
    R["peak_oc"] = oc
    # peak error over arrest torques at the requirement chain
    R["peak_err_by_Th"] = {
        f"{Th:g}": _stats(
            np.concatenate(
                [
                    measure(trajs[(Th, dl)], Chain(skew_uniform=skew_max), rng(), n=200)["peak"]
                    - trajs[(Th, dl)].peak
                    for dl in (-0.01, 0.0, 0.01)
                ]
            )
        )
        for Th in T_ARREST
    }

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
    edge_trajs = [
        simulate(t_on_for_peak(DROP_LIM + dl, T_ARREST_MAIN), T_ARREST_MAIN)
        for dl in np.linspace(0.002, 0.05, 13)
    ]
    scan = []
    for bits in (10, 11, 12, 13, 14, 15, 16, 17):
        q = 360.0 / 2**bits
        row = dict(bits=bits, q=q)
        for nm, sg in (("plain", 0.0), ("dither", 0.5 * q)):
            row[f"miss_max_{nm}"] = max(
                float(
                    (
                        ~score_angle(
                            measure(tr, Chain(q=q, sigma=sg, skew_uniform=skew_max), rng(), n=100)["peak"]
                        )
                    ).mean()
                )
                for tr in edge_trajs
            )
        scan.append(row)
    lim["q_scan"] = scan
    lim["q_bits_last_miss"] = max(
        r["bits"] for r in scan if r["miss_max_plain"] > 0 or r["miss_max_dither"] > 0
    )
    lim["q_bits_guard"] = math.ceil(math.log2(360.0 / GUARD_DEG))  # smallest bits with q < guard
    lim["analytic"] = dict(gain=-(GUARD_DEG + 0.002) / (DROP_LIM + 0.002), zero=-(GUARD_DEG + 0.002))
    R["miss_limits"] = lim

    # 5. catch-time channel (speed zero / argmax) vs U_t
    tc = {}
    for Th in T_ARREST:
        tr = trajs[(Th, 0.0)]
        A_dec = (Th - A["G"] * math.cos(math.radians(DROP_LIM))) / A["J"]  # rad/s^2 near the peak
        row = dict(
            t_c_true_ms=tr.t_c * 1e3,
            decel_dps2=math.degrees(A_dec),
            flat_ms_1lsb=math.sqrt(2 * Q_REQ / math.degrees(A_dec)) * 1e3,
        )
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
    pr = dict(
        e_src=pr_sweep(E_SRC),
        e_src_clip=pr_sweep(E_SRC, clip=True),
        e_cap=pr_sweep(E_CAP),
        k_check=k_ratio(0.03, -0.02)[:2],
    )
    # peak channel when the score reads the device joint only, at the e_cap extremes (linear k)
    for nm, kk in (("k_lo", pr["e_cap"]["k_min"]), ("k_hi", pr["e_cap"]["k_max"])):
        m = measure(trajs[(T_ARREST_MAIN, 0.0)], Chain(k_dev=kk, skew_uniform=skew_max), rng(), n=200)
        pr[f"joint_only_{nm}_peak_err"] = _stats(m["peak"] - DROP_LIM)
        m = measure(
            trajs[(T_ARREST_MAIN, 0.0)], Chain(k_dev=kk, pivot=True, skew_uniform=skew_max), rng(), n=200
        )
        pr[f"joint_pivot_{nm}_peak_err"] = _stats(m["peak"] - DROP_LIM)
    R["pr"] = pr

    # 7. encoder faults and plausibility checks (true peak 10.01 deg: a real drop past 10)
    tr_f = trajs[(T_ARREST_MAIN, 0.01)]
    tr_reb = simulate(
        t_on_for_peak(DROP_LIM + 0.01, T_ARREST_MAIN, rebound=True), T_ARREST_MAIN, rebound=True
    )
    faults = {}

    def fault(name, ch, tr=tr_f):
        m = measure(tr, ch, rng())
        drop = score_angle(m["peak"])
        inval = m["motion_bad"] | m["decel_bad"] | m["gap_bad"]
        faults[name] = dict(
            missed=float((~drop).mean()),
            missed_after_checks=float((~drop & ~inval).mean()),
            flagged=float(inval.mean()),
            motion_bad=float(m["motion_bad"].mean()),
            decel_bad=float(m["decel_bad"].mean()),
            gap_bad=float(m["gap_bad"].mean()),
        )
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
    for nm, ch in (
        ("glitch_0p1", Chain(glitch=0.1, skew_uniform=skew_max)),
        ("glitch_0p1_med3", Chain(glitch=0.1, med3=True, skew_uniform=skew_max)),
    ):
        m = measure(tr_ok, ch, rng())
        faults[nm] = dict(false_trip=float(score_angle(m["peak"]).mean()))
    # clean trials must not be flagged: check over all arrest torques at delta 0
    flags = [
        measure(trajs[(Th, 0.0)], Chain(sigma=Q_REQ, skew_uniform=skew_max), rng(), n=200) for Th in T_ARREST
    ]
    faults["clean_flag_rate_all_Th"] = float(
        np.mean([(f["motion_bad"] | f["decel_bad"]).mean() for f in flags])
    )
    faults["clean_max_decel_all_Th"] = float(max(f["max_decel"].max() for f in flags))
    faults["model_max_decel"] = max(
        math.degrees((Th - A["G"] * math.cos(math.radians(DROP_LIM))) / A["J"]) for Th in T_ARREST
    )
    # decel threshold sweep: stuck-during-drop misses after checks vs clean false flags (noise 1 LSB)
    st_drop = score_angle(st["peak"])
    sweep = []
    for a_lim in (10000.0, 20000.0, 50000.0, 100000.0):
        miss = ~st_drop & ~(st["motion_bad"] | (st["max_decel"] > a_lim) | st["gap_bad"])
        clean = np.mean([(f["motion_bad"] | (f["max_decel"] > a_lim)).mean() for f in flags])
        sweep.append(
            dict(a_lim=a_lim, stuck_missed_after_checks=float(miss.mean()), clean_flagged=float(clean))
        )
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
    R["r2"] = round2(rng, trajs)  # fix round 2: appended so round-1 draws are unchanged
    derive(R)
    return R


# ---------------------------------------------------------------- fix round 2 (FC round 1 on PR #14)
# Two-channel chain (forearm = joint - pivot), per-channel and common-mode faults, frozen-count and
# second-reading checks. ALL SYNTHETIC. ASSUMPTIONS for this part are listed in R2_ASSUMPTIONS.
K_LO_CAP, K_HI_CAP = 0.7344681196147431, 1.2655318803852569  # 46.3 mm disc extremes (pr.e_cap, round 1)
N_R2 = 2000  # trials per fault row
N_SENS = 4000  # trials per stuck-sensitivity row
N_BUDGET = 200000  # trials for the two-channel budget check
RUN_LIM = 0.010  # s: frozen-count check, longest run of identical raw counts after onset (ASSUMPTION)
XCHK_TOL = 0.05  # deg: second-reading cross-check tolerance on the peak drop (ASSUMPTION)
R2_ASSUMPTIONS = dict(
    run_lim_ms=RUN_LIM * 1e3,
    xchk_tol_deg=XCHK_TOL,
    n_r2=N_R2,
    n_sens=N_SENS,
    n_budget=N_BUDGET,
    t_end_s=T_END,
    glitch_window_s=0.3,
)


@dataclass
class Fault:
    """kind: stuck (hold from t), intermittent (hold for ``size`` s episodes at ``rate`` per s), slip (step of
    ``size`` deg at t, sign that makes the forearm drop read low), backlash (play of width ``size`` deg,
    random start inside the play), drift (``size`` deg/s, forearm drop read low), gain (scale error ``size``).
    t = nan: uniform over [0, t_c] per trial."""

    kind: str
    size: float = 0.0
    t: float = float("nan")
    rate: float = 0.0


@dataclass
class Chain2:
    k: float = 1.0
    pivot: bool = True
    fs: float = FS_REQ
    q: float = Q_REQ
    q_pivot: float | None = None
    lat: float = ENC_LAT_S
    skew_uniform: float = TRIG_SKEW_S + REL_SKEW_S
    sigma_counts: float = 1.0  # noise per channel, in counts
    joint: Fault | None = None
    piv: Fault | None = None
    common: Fault | None = None  # applied to both device channels (and to the reference if common_ref)
    common_ref: bool = False
    a_lim: float = A_LIM


def _play(x, b, rng):
    """Backlash (play) operator, width b, random initial position inside the play."""
    y = np.empty_like(x)
    cur = x[:, 0] + rng.uniform(-b / 2, b / 2, x.shape[0])
    for j in range(x.shape[1]):
        cur = np.clip(cur, x[:, j] - b / 2, x[:, j] + b / 2)
        y[:, j] = cur
    return y


def _hold(v, mask):
    """Hold the last value before each True run of ``mask`` (row-wise)."""
    idx = np.where(mask, 0, np.arange(v.shape[1])[None, :])
    idx = np.maximum.accumulate(idx, axis=1)
    return np.take_along_axis(v, idx, axis=1)


def _longest_run(c, post):
    """Longest run (in samples) of identical values in c where post is True."""
    same = np.zeros(c.shape, bool)
    same[:, 1:] = (c[:, 1:] == c[:, :-1]) & post[:, 1:] & post[:, :-1]
    ar = np.arange(c.shape[1])[None, :]
    last = np.maximum.accumulate(np.where(~same, ar, 0), axis=1)
    return (ar - last).max(axis=1)


def measure2(tr: Traj, ch: Chain2, rng: np.random.Generator, n: int, ref: bool = True, chunk: int = 250):
    outs = [_measure2(tr, ch, rng, min(chunk, n - i), ref) for i in range(0, n, chunk)]
    return {k: np.concatenate([o[k] for o in outs]) for k in outs[0]}


def _measure2(tr, ch, rng, n, ref):
    nper = int(round((T_END - BASE_WIN[0]) * ch.fs))
    phase = rng.uniform(0.0, 1.0 / ch.fs, n)
    ts = BASE_WIN[0] + phase[:, None] + np.arange(nper)[None, :] / ch.fs
    skew = rng.uniform(-ch.skew_uniform, ch.skew_uniform, n)
    tm = ts - skew[:, None]
    ang = -tr.drop(ts - ch.lat)
    t_c = tr.t_c if math.isfinite(tr.t_c) else tr.t_peak

    def ftime(f):
        return rng.uniform(0.0, t_c, n) if math.isnan(f.t) else np.full(n, f.t)

    mounts: list = []  # realised mount phases (for the design checks in the tests)

    def channel(x, q, f, sign, common):
        """x: true channel angle; sign +1 if a positive reading error makes the forearm drop read low."""
        for flt in (f, common):
            if flt is None:
                continue
            if flt.kind == "gain":
                x = (1 + flt.size) * x
            elif flt.kind == "drift":
                x = x + sign * flt.size * (ts - BASE_WIN[0])
            elif flt.kind == "backlash":
                x = _play(x, flt.size, rng)
        x = x + rng.normal(0.0, ch.sigma_counts * q, x.shape) if ch.sigma_counts > 0 else x
        mount = rng.uniform(0.0, q, (n, 1))
        mounts.append(mount[:, 0])
        v = np.floor((x + mount) / q) * q
        for flt, tt in ((f, None), (common, common_t)):
            if flt is None:
                continue
            if flt.kind == "stuck":
                t_s = ftime(flt) if tt is None else tt
                v = _hold(v, ts >= t_s[:, None])
            elif flt.kind == "intermittent":
                starts = rng.uniform(
                    BASE_WIN[0], T_END, (n, max(1, int(rng.poisson(flt.rate * (T_END - BASE_WIN[0]))) + 1))
                )
                ep = np.zeros(ts.shape, bool)
                for s0 in starts.T:
                    ep |= (ts >= s0[:, None]) & (ts < s0[:, None] + flt.size)
                v = _hold(v, ep)
            elif flt.kind == "slip":
                t_s = ftime(flt) if tt is None else tt
                v = v + sign * flt.size * (ts >= t_s[:, None])
        return v

    common_t = ftime(ch.common) if ch.common is not None and ch.common.kind in ("stuck", "slip") else None
    vj = channel(ch.k * ang, ch.q, ch.joint, +1.0, ch.common)
    if ch.pivot:
        vp = channel((ch.k - 1.0) * ang, ch.q_pivot or ch.q, ch.piv, -1.0, ch.common)
        val = vj - vp
    else:
        vp, val = None, vj
    base_m = (tm >= BASE_WIN[0]) & (tm < BASE_WIN[1])
    base = np.where(base_m, val, 0).sum(1) / base_m.sum(1)
    dm = base[:, None] - val
    post = tm >= 0
    peak = np.where(post, dm, -np.inf).max(1)
    rows = np.arange(n)
    hit = (dm >= DROP_LIM) & post
    i10 = np.argmax(hit, 1)
    t10 = np.where(hit.any(1), tm[rows, i10], np.nan)
    i200 = (tm <= CATCH_MS / 1e3).sum(1) - 1
    d200 = dm[rows, i200]
    i100 = (tm <= MOTION_T).sum(1) - 1
    motion_bad = dm[rows, i100] < MOTION_FRAC * free_fall_drop(MOTION_T)
    Wc = max(1, int(round(1e-3 * ch.fs)))
    vc = (dm[:, Wc:] - dm[:, :-Wc]) * ch.fs / Wc
    max_dec = (-(vc[:, Wc:] - vc[:, :-Wc]) * ch.fs / Wc).max(1)
    run_j = _longest_run(vj, post) / ch.fs
    run_p = _longest_run(vp, post) / ch.fs if vp is not None else np.zeros(n)
    out = dict(
        peak=peak,
        t10=t10,
        d200=d200,
        motion_bad=motion_bad,
        max_dec=max_dec,
        run=np.maximum(run_j, run_p),
        skew=skew,
        phase=phase,
        mount_j=mounts[0],
        base_sd=np.sqrt(np.where(base_m, (val - base[:, None]) ** 2, 0).sum(1) / (base_m.sum(1) - 1)),
    )
    if ref:  # second, independent forearm reading: own quantizer, mount, noise; shares only the common fault
        x = ang + rng.normal(0.0, ch.sigma_counts * ch.q, ang.shape) if ch.sigma_counts > 0 else ang
        vr = np.floor((x + rng.uniform(0.0, ch.q, (n, 1))) / ch.q) * ch.q
        if ch.common_ref and ch.common is not None and ch.common.kind == "stuck":
            vr = _hold(vr, ts >= common_t[:, None])
        br = np.where(base_m, vr, 0).sum(1) / base_m.sum(1)
        out["ref_peak"] = np.where(post, br[:, None] - vr, -np.inf).max(1)
    return out


def _cp(x, n, conf=0.95):
    from scipy.stats import beta

    a = (1 - conf) / 2
    lo = 0.0 if x == 0 else float(beta.ppf(a, x, n - x + 1))
    hi = 1.0 if x == n else float(beta.ppf(1 - a, x + 1, n - x))
    return [lo, hi]


def _rate(mask):
    return dict(p=float(mask.mean()), ci=_cp(int(mask.sum()), mask.size), n=int(mask.size))


def fault_row(m, a_lim, run_lim=RUN_LIM, tol=XCHK_TOL):
    drop = score_angle(m["peak"])
    base = m["motion_bad"] | (m["max_dec"] > a_lim)
    run = base | (m["run"] > run_lim)
    xchk = run | (np.abs(m["peak"] - m["ref_peak"]) > tol)
    return dict(
        drop=_rate(drop),
        missed=_rate(~drop),
        missed_base=_rate(~drop & ~base),
        missed_run=_rate(~drop & ~run),
        missed_xchk=_rate(~drop & ~xchk),
        flagged_base=float(base.mean()),
        flagged_run=float(run.mean()),
        flagged_xchk=float(xchk.mean()),
    )


def budget_two_channel(
    k,
    rng,
    n=N_BUDGET,
    fs=FS_REQ,
    q=Q_REQ,
    q_pivot=None,
    lat=ENC_LAT_S,
    skew=-(TRIG_SKEW_S + REL_SKEW_S),
    interp=False,
    pivot=True,
    chunk=5000,
):
    """Free-fall governing cell, skew at the edge that makes times read late: worst crossing error (ms) and
    worst drop-at-200 ms error (deg) over n random phase/mount draws. Joint - pivot when ``pivot``."""
    tr = free_fall_traj()
    ff = free_fall_points()
    t10, d200 = ff["t10_ms"] / 1e3, ff["drop200_deg"]
    js = np.arange(int(0.235 * fs), int(0.262 * fs))
    qp = q_pivot or q
    wt, wd = -np.inf, np.inf
    for _ in range(0, n, chunk):
        ph = rng.uniform(0, 1 / fs, chunk)
        ts = BASE_WIN[0] + ph[:, None] + js[None, :] / fs
        tm = ts - skew
        ang = -tr.drop(ts - lat)
        m1 = rng.uniform(0, q, (chunk, 1))
        v = np.floor((k * ang + m1) / q) * q
        base = np.floor(m1 / q) * q
        if pivot:
            m2 = rng.uniform(0, qp, (chunk, 1))
            v = v - np.floor(((k - 1) * ang + m2) / qp) * qp
            base = base - np.floor(m2 / qp) * qp
        dm = (base - v) / (k if not pivot else 1.0)
        rows = np.arange(chunk)
        i1 = np.argmax(dm >= DROP_LIM, 1)
        if interp:
            d0, d1 = dm[rows, i1 - 1], dm[rows, i1]
            fr = np.where(d1 > d0, (DROP_LIM - d0) / np.where(d1 > d0, d1 - d0, 1), 1.0)
            tcr = tm[rows, i1 - 1] + fr / fs
        else:
            tcr = tm[rows, i1]
        i2 = (tm <= CATCH_MS / 1e3).sum(1) - 1
        if interp:
            fr2 = (CATCH_MS / 1e3 - tm[rows, i2]) * fs
            dd = dm[rows, i2] + fr2 * (dm[rows, i2 + 1] - dm[rows, i2])
        else:
            dd = dm[rows, i2]
        wt = max(wt, float(((tcr - t10) * 1e3).max()))
        wd = min(wd, float((dd - d200).min()))
    return dict(
        t10_err_max_ms=wt, d200_err_min_deg=wd, t10_ok=bool(wt <= U_T_S * 1e3), d200_ok=bool(-wd <= GUARD_DEG)
    )


def free_fall_traj():
    if "tr" not in _FF_CACHE:
        _FF_CACHE["tr"] = simulate()
    return _FF_CACHE["tr"]


def exact_single_channel(n_phase=2001, n_mount=1001):
    """Brute-force worst case of the one-channel requirement chain over sampling phase x mount phase, at
    both skew edges (the sup of a 400-trial MC is a lower bound only)."""
    tr = free_fall_traj()
    ff = free_fall_points()
    t10, d200 = ff["t10_ms"] / 1e3, ff["drop200_deg"]
    fs, q, lat = FS_REQ, Q_REQ, ENC_LAT_S
    js = np.arange(int(0.235 * fs), int(0.262 * fs))
    mn = np.linspace(0, q, n_mount, endpoint=False)[:, None]
    res = {}
    for sk in (-(TRIG_SKEW_S + REL_SKEW_S), TRIG_SKEW_S + REL_SKEW_S):
        wt_hi, wt_lo, wd_lo, wd_hi = -np.inf, np.inf, np.inf, -np.inf
        for p in np.linspace(0, 1 / fs, n_phase, endpoint=False):
            ts = BASE_WIN[0] + p + js / fs
            tm = ts - sk
            dm = np.floor(mn / q) * q - np.floor((-tr.drop(ts - lat)[None, :] + mn) / q) * q
            i1 = np.argmax(dm >= DROP_LIM, 1)
            e = (tm[i1] - t10) * 1e3
            i2 = int((tm <= CATCH_MS / 1e3).sum()) - 1
            ed = dm[:, i2] - d200
            wt_hi, wt_lo = max(wt_hi, e.max()), min(wt_lo, e.min())
            wd_lo, wd_hi = min(wd_lo, ed.min()), max(wd_hi, ed.max())
        res["neg" if sk < 0 else "pos"] = dict(
            t10_lo=float(wt_lo), t10_hi=float(wt_hi), d200_lo=float(wd_lo), d200_hi=float(wd_hi)
        )
    wt = max(max(abs(v["t10_lo"]), abs(v["t10_hi"])) for v in res.values())
    wd = max(max(abs(v["d200_lo"]), abs(v["d200_hi"])) for v in res.values())
    ffp = free_fall_points()
    w10 = ffp["w10_dps"]
    # lattice-free bound: every timing term at its edge plus one count at the crossing speed (PR #9's sum)
    lf_t = (ENC_LAT_S + TRIG_SKEW_S + REL_SKEW_S + 1 / FS_REQ) * 1e3 + Q_REQ / w10 * 1e3
    return dict(
        edges=res,
        t10_abs_max_ms=float(wt),
        d200_abs_max_deg=float(wd),
        u_t_room_ms=U_T_S * 1e3 - float(wt),
        guard_room_deg=GUARD_DEG - float(wd),
        lattice_free_t10_ms=lf_t,
        lattice_free_room_ms=U_T_S * 1e3 - lf_t,
        n_phase=n_phase,
        n_mount=n_mount,
    )


def round2(rng, trajs) -> dict:
    out: dict = dict(assumptions=R2_ASSUMPTIONS, pr9_head="335d4e586dd51cf818e0153672a31d0e37d90ea0")
    # item 6: exact worst case of the one-channel chain
    out["exact"] = exact_single_channel()
    # items 1 and 5: joint - pivot through PR #9's chain; options that would pass (none chosen)
    bud = {}
    for nm, kk in (("k_lo", K_LO_CAP), ("k_1", 1.0001), ("k_hi", K_HI_CAP), ("k_hi_100mm", K_HI)):
        bud[nm] = budget_two_channel(kk, rng())
    bud["one_channel"] = budget_two_channel(1.0, rng(), pivot=False)
    opts = dict(
        pivot_18bit=dict(q_pivot=Q_REQ / 2),
        both_18bit=dict(q=Q_REQ / 2),
        fs_20k=dict(fs=20000.0),
        lat_0=dict(lat=0.0),
        interp=dict(interp=True),
    )
    bud["options"] = {}
    for nm, kw in opts.items():
        rows = {
            kn: budget_two_channel(kk, rng(), **kw)
            for kn, kk in (("k_lo", K_LO_CAP), ("k_1", 1.0001), ("k_hi", K_HI_CAP))
        }
        rows["passes_all"] = all(r["t10_ok"] and r["d200_ok"] for r in rows.values())
        bud["options"][nm] = rows
    # option 6 (Charlie's call, loosens thresholds): add one pivot count to the budget terms
    bud["loosen"] = dict(
        u_t_ms=U_T_S * 1e3 + Q_REQ / free_fall_points()["w10_dps"] * 1e3, guard_deg=GUARD_DEG + Q_REQ
    )
    bud["two_count_bound_deg"] = 2 * Q_REQ
    out["budget"] = bud
    # healthy joint - pivot peak error (true 10.000 deg), many trials: the bound is 2 counts
    m = measure2(trajs[(T_ARREST_MAIN, 0.0)], Chain2(k=K_LO_CAP, sigma_counts=0.0), rng(), 4000, ref=False)
    out["jp_peak_err"] = _stats(m["peak"] - DROP_LIM)
    # items 2 and 4: fault rows (k_lo, true peak 10.01, T_h 30, 1 count noise per channel)
    tr_f = trajs[(T_ARREST_MAIN, 0.01)]
    clean = [
        measure2(trajs[(Th, 0.0)], Chain2(k=kk), rng(), 1000)
        for Th in T_ARREST
        for kk in (K_LO_CAP, 1.0001, K_HI_CAP)
    ]
    cmax = float(max(c["max_dec"].max() for c in clean))
    a_tight = math.ceil(cmax / 500.0) * 500.0
    out["clean"] = dict(
        max_dec=cmax,
        a_lim_tight=a_tight,
        max_run_ms=float(max(c["run"].max() for c in clean)) * 1e3,
        flagged_base=float(np.mean([((c["motion_bad"]) | (c["max_dec"] > a_tight)).mean() for c in clean])),
        flagged_run=float(
            np.mean(
                [((c["motion_bad"]) | (c["max_dec"] > a_tight) | (c["run"] > RUN_LIM)).mean() for c in clean]
            )
        ),
        flagged_xchk=float(
            np.mean(
                [
                    (
                        (c["motion_bad"])
                        | (c["max_dec"] > a_tight)
                        | (c["run"] > RUN_LIM)
                        | (np.abs(c["peak"] - c["ref_peak"]) > XCHK_TOL)
                    ).mean()
                    for c in clean
                ]
            )
        ),
        n=int(sum(c["peak"].size for c in clean)),
    )
    U = float("nan")
    rows_def = {
        "clean": Chain2(k=K_LO_CAP),
        "joint_stuck_during": Chain2(k=K_LO_CAP, joint=Fault("stuck", t=U)),
        "pivot_stuck_onset": Chain2(k=K_LO_CAP, piv=Fault("stuck", t=0.0)),
        "pivot_stuck_during": Chain2(k=K_LO_CAP, piv=Fault("stuck", t=U)),
        "joint_intermittent_5ms": Chain2(k=K_LO_CAP, joint=Fault("intermittent", size=0.005, rate=20.0)),
        "pivot_intermittent_5ms": Chain2(k=K_LO_CAP, piv=Fault("intermittent", size=0.005, rate=20.0)),
        "joint_slip_0p05": Chain2(k=K_LO_CAP, joint=Fault("slip", size=0.05, t=U)),
        "joint_slip_0p2": Chain2(k=K_LO_CAP, joint=Fault("slip", size=0.2, t=U)),
        "pivot_slip_0p05": Chain2(k=K_LO_CAP, piv=Fault("slip", size=0.05, t=U)),
        "joint_backlash_0p05": Chain2(k=K_LO_CAP, joint=Fault("backlash", size=0.05)),
        "joint_backlash_0p2": Chain2(k=K_LO_CAP, joint=Fault("backlash", size=0.2)),
        "joint_drift_0p2dps": Chain2(k=K_LO_CAP, joint=Fault("drift", size=0.2)),
        "pivot_drift_0p2dps": Chain2(k=K_LO_CAP, piv=Fault("drift", size=0.2)),
        "common_daq_freeze": Chain2(k=K_LO_CAP, common=Fault("stuck", t=U)),
        "common_daq_freeze_incl_ref": Chain2(k=K_LO_CAP, common=Fault("stuck", t=U), common_ref=True),
        "common_gain_m0p5pct": Chain2(k=K_LO_CAP, common=Fault("gain", size=-0.005)),
        "pivot_stuck_onset_noise0": Chain2(k=K_LO_CAP, piv=Fault("stuck", t=0.0), sigma_counts=0.0),
        "joint_stuck_during_noise0": Chain2(k=K_LO_CAP, joint=Fault("stuck", t=U), sigma_counts=0.0),
    }
    out["faults"] = {nm: fault_row(measure2(tr_f, ch, rng(), N_R2), a_tight) for nm, ch in rows_def.items()}
    # pivot stuck at k_hi over-reads: false DROP on an arm truly at 9.95 deg
    m = measure2(trajs[(T_ARREST_MAIN, -0.05)], Chain2(k=K_HI_CAP, piv=Fault("stuck", t=0.0)), rng(), N_R2)
    fr = fault_row(m, a_tight)
    out["faults"]["pivot_stuck_onset_khi_false_drop"] = dict(
        false_drop=fr["drop"], flagged_run=fr["flagged_run"], flagged_xchk=fr["flagged_xchk"]
    )
    # frozen-count check vs noise: clean max run and stuck min run (one channel, k 1, T_h 30)
    sep = []
    for sc in (0.0, 0.25, 0.5, 1.0):
        c = measure2(
            trajs[(T_ARREST_MAIN, 0.0)], Chain2(pivot=False, sigma_counts=sc), rng(), 1000, ref=False
        )
        s = measure2(
            tr_f, Chain2(pivot=False, sigma_counts=sc, joint=Fault("stuck", t=U)), rng(), 1000, ref=False
        )
        sep.append(
            dict(
                sigma_counts=sc,
                clean_max_run_ms=float(c["run"].max() * 1e3),
                stuck_min_run_ms=float(s["run"].min() * 1e3),
                clean_flag_rate=float((c["run"] > RUN_LIM).mean()),
                separated=bool(c["run"].max() < s["run"].min()),
            )
        )
    out["run_sep"] = sep
    # item 3: stuck-during-drop sensitivity, one channel (as round 1), n = N_SENS, 95 % CP intervals.
    # Tight deceleration limit for one channel: just above the clean one-channel max (1 count noise)
    clean1 = [measure2(trajs[(Th, 0.0)], Chain2(pivot=False), rng(), 1000, ref=False) for Th in T_ARREST]
    cmax1 = float(max(c["max_dec"].max() for c in clean1))
    a_tight1 = math.ceil(cmax1 / 500.0) * 500.0
    out["clean_1ch"] = dict(
        max_dec=cmax1,
        a_lim_tight=a_tight1,
        max_run_ms=float(max(c["run"].max() for c in clean1)) * 1e3,
        flagged_tight_run=float(
            np.mean(
                [(c["motion_bad"] | (c["max_dec"] > a_tight1) | (c["run"] > RUN_LIM)).mean() for c in clean1]
            )
        ),
        n=int(sum(c["peak"].size for c in clean1)),
    )
    sens = {}

    def stuck_case(tr, t_hi, sigma_counts=0.0):
        ch = Chain2(pivot=False, sigma_counts=sigma_counts, joint=Fault("stuck", t=float("nan")))
        tr2 = Traj(
            tr.t, tr.d, t_hi, tr.peak, tr.t_on, tr.T_h, tr.rebound, tr.t_peak
        )  # freeze window U[0, t_hi]
        m = measure2(tr2, ch, rng(), N_SENS, ref=False)
        drop = score_angle(m["peak"])
        res = dict(missed=_rate(~drop))
        for nm, al in (("a50k", A_LIM), ("a_tight", a_tight1)):
            res[f"after_{nm}"] = _rate(~drop & ~(m["motion_bad"] | (m["max_dec"] > al)))
        res["after_tight_run"] = _rate(
            ~drop & ~(m["motion_bad"] | (m["max_dec"] > a_tight1) | (m["run"] > RUN_LIM))
        )
        miss = ~drop & ~(m["motion_bad"] | (m["max_dec"] > A_LIM))
        res["share_gross_lt5"] = float((m["peak"][miss] < 5.0).mean()) if miss.any() else 0.0
        return res

    base_tr = trajs[(T_ARREST_MAIN, 0.01)]
    sens["base"] = stuck_case(base_tr, base_tr.t_c)
    sens["base_noise1"] = stuck_case(base_tr, base_tr.t_c, 1.0)
    for dl in (0.1, 1.0, 3.0):
        tr = simulate(t_on_for_peak(DROP_LIM + dl, T_ARREST_MAIN), T_ARREST_MAIN)
        sens[f"delta_{dl:g}"] = stuck_case(tr, tr.t_c)
    for Th in (15.0, 20.0, 50.0):
        tr = trajs[(Th, 0.01)]
        sens[f"Th_{Th:g}"] = stuck_case(tr, tr.t_c)
    sens["window_t10"] = stuck_case(base_tr, float(base_tr.t[np.argmax(base_tr.d >= DROP_LIM)]))
    sens["window_tend"] = stuck_case(base_tr, T_END)
    a50 = [v["after_a50k"]["p"] for v in sens.values()]
    sens["range_after_a50k"] = [min(a50), max(a50)]
    sens["round1_ci_n400"] = _cp(105, 400)  # round 1: 26.2 % = 105 / 400
    out["stuck_sens"] = sens
    return out


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
    pe["worst_min"] = min(v["min"] for v in pe.values() if isinstance(v, dict))
    pe["worst_max"] = max(v["max"] for v in pe.values() if isinstance(v, dict))
    pr = R["pr"]
    pr["joint_pivot_abs_max"] = max(
        abs(pr[f"joint_pivot_{k}_peak_err"][m]) for k in ("k_lo", "k_hi") for m in ("min", "max")
    )
    rows = [v for v in R["catch_time"].values() if isinstance(v, dict)]

    def ext(key):
        return max(max(abs(r[key]["min"]), abs(r[key]["max"])) for r in rows)

    ct = R["catch_time"]
    ct["worst_arg_ms"] = ext("tc_arg_err_ms")
    ct["worst_v15_ms"] = ext("tc_v_W15_err_ms")
    ct["worst_noise_ms"] = max(ext("tc_arg_noise1_err_ms"), ext("tc_v_W15_noise1_err_ms"))
    ct["worst_noise_arg_ms"] = max(r["tc_arg_noise1_err_ms"]["max"] for r in rows)
    # fix round 2 summaries
    r2 = R["r2"]
    b = r2["budget"]
    ks = ("k_lo", "k_1", "k_hi", "k_hi_100mm")
    b["two_ch_worst_t10_ms"] = max(b[k]["t10_err_max_ms"] for k in ks)
    b["two_ch_worst_d200_deg"] = min(b[k]["d200_err_min_deg"] for k in ks)
    b["two_ch_t10_fails"] = [k for k in ks if not b[k]["t10_ok"]]
    b["two_ch_d200_fails"] = [k for k in ks if not b[k]["d200_ok"]]
    b["two_ch_t10_over_ms"] = b["two_ch_worst_t10_ms"] - U_T_S * 1e3
    b["two_ch_d200_over_deg"] = -b["two_ch_worst_d200_deg"] - GUARD_DEG
    sens = r2["stuck_sens"]
    cases = [v for v in sens.values() if isinstance(v, dict) and "missed" in v]
    for key in ("missed", "after_a50k", "after_a_tight"):
        vals = [c[key]["p"] for c in cases]
        sens[f"range_{key}"] = [min(vals), max(vals)]
    # the frozen-count check is only usable with noise (a held, noiseless arm is frozen too: r2.run_sep),
    # so its result is quoted only for the 1-count-noise row
    sens.pop("range_after_tight_run", None)
    sens["ci_halfwidth_n400_pts"] = (sens["round1_ci_n400"][1] - sens["round1_ci_n400"][0]) / 2 * 100


def _f(x, fmt):
    return format(x, fmt)


def tables(R) -> dict:
    """Markdown table bodies for the doc, built only from R."""
    oc = R["peak_oc"]
    names = ("req", "noise_1lsb", "noise_4lsb", "noise_4lsb_med3", "revc_joint", "imu150")
    lab = dict(
        req="req",
        noise_1lsb="1 count noise",
        noise_4lsb="4 counts noise",
        noise_4lsb_med3="4 counts + median-3",
        revc_joint="0.1° steps",
        imu150="150 Hz",
    )
    deltas = [r["delta"] for r in oc["req"]["rows"]]
    oc_rows = []
    for i, dl in enumerate(deltas):
        cells = [_f(oc[n]["rows"][i]["p_drop"], ".3f") for n in names]
        oc_rows.append(f"| {dl:+.4f} | " + " | ".join(cells) + " |")
    summ = [
        f"| {lab[n]} | {_f(oc[n]['err_min'], '+.4f')} / {_f(oc[n]['err_max'], '+.4f')} | "
        f"{_f(oc[n]['max_missed'], '.3f')} | {oc[n]['trip_any_delta']:+.4f} | "
        f"{oc[n]['trip_all_delta']:+.4f} |"
        for n in names
    ]
    q_rows = [
        f"| {r['bits']} | {r['q']:.5f} | {r['miss_max_plain']:.2f} | {r['miss_max_dither']:.2f} |"
        for r in R["miss_limits"]["q_scan"]
    ]

    def rng_(d):
        return f"{d['min']:+.2f} to {d['max']:+.2f}"

    tc_rows = []
    for th in ("15", "20", "30", "50"):
        r = R["catch_time"][th]
        tc_rows.append(
            f"| {th} | {r['peak_at_time_edge_deg']:.3f} | {rng_(r['tc_arg_err_ms'])} | "
            f"{rng_(r['tc_v_W15_err_ms'])} | {rng_(r['tc_v_W75_err_ms'])} | "
            f"{rng_(r['tc_v_W1_err_ms'])} | "
            f"{rng_(r['tc_arg_noise1_err_ms'])} | {rng_(r['tc_v_W15_noise1_err_ms'])} | "
            f"{r['flat_ms_1lsb']:.3f} | {r['bits_for_u_t']} |"
        )
    dec_rows = [
        f"| {r['a_lim']:.0f} | {r['stuck_missed_after_checks']:.1%} | {r['clean_flagged']:.1%} |"
        for r in R["faults"]["decel_sweep"]
    ]
    t = dict(
        oc_table="\n".join(oc_rows),
        oc_summary="\n".join(summ),
        q_table="\n".join(q_rows),
        tc_table="\n".join(tc_rows),
        decel_table="\n".join(dec_rows),
    )
    if "r2" in R:
        t.update(r2_tables(R["r2"]))
    return t


def _pc(r):
    """Rate with its 95 % Clopper-Pearson interval, in percent."""
    return f"{r['p']:.1%} ({r['ci'][0]:.1%}–{r['ci'][1]:.1%})"


def r2_tables(r2) -> dict:
    b = r2["budget"]
    lab = dict(
        k_lo="0.7345 (46.3 mm low)",
        k_1="1.0001",
        k_hi="1.2655 (46.3 mm high)",
        k_hi_100mm="1.5862 (100 mm high)",
    )
    bud = []
    for k, nm in lab.items():
        r = b[k]
        bud.append(
            f"| joint − pivot, k {nm} | {r['t10_err_max_ms']:.5f} | "
            f"{'pass' if r['t10_ok'] else '**FAIL**'} | "
            f"{r['d200_err_min_deg']:.5f} | {'pass' if r['d200_ok'] else '**FAIL**'} |"
        )
    r = b["one_channel"]
    bud.append(
        f"| one channel (forearm read directly) | {r['t10_err_max_ms']:.5f} | "
        f"{'pass' if r['t10_ok'] else '**FAIL**'} | {r['d200_err_min_deg']:.5f} | "
        f"{'pass' if r['d200_ok'] else '**FAIL**'} |"
    )
    olab = dict(
        pivot_18bit="pivot channel 18 bits (joint 17)",
        both_18bit="both channels 18 bits",
        fs_20k="sampling 20 kHz",
        lat_0="angle-path latency 0 (bound case)",
        interp="linear interpolation of crossing and drop at 200 ms",
    )
    opt = []
    for k, nm in olab.items():
        v = b["options"][k]
        wt = max(v[x]["t10_err_max_ms"] for x in ("k_lo", "k_1", "k_hi"))
        wd = min(v[x]["d200_err_min_deg"] for x in ("k_lo", "k_1", "k_hi"))
        opt.append(f"| {nm} | {wt:.5f} | {wd:.5f} | {'passes' if v['passes_all'] else '**fails**'} |")
    flab = [
        ("clean", "No fault"),
        ("joint_stuck_during", "Joint stuck, random time during the drop"),
        ("pivot_stuck_onset", "Pivot stuck from onset"),
        ("pivot_stuck_during", "Pivot stuck during the drop"),
        ("joint_intermittent_5ms", "Joint partial sticking (5 ms holds, 20 per s)"),
        ("pivot_intermittent_5ms", "Pivot partial sticking (5 ms holds, 20 per s)"),
        ("joint_slip_0p05", "Joint coupling slip 0.05°"),
        ("joint_slip_0p2", "Joint coupling slip 0.2°"),
        ("pivot_slip_0p05", "Pivot coupling slip 0.05°"),
        ("joint_backlash_0p05", "Joint backlash 0.05°"),
        ("joint_backlash_0p2", "Joint backlash 0.2°"),
        ("joint_drift_0p2dps", "Joint zero drift 0.2 °/s"),
        ("pivot_drift_0p2dps", "Pivot zero drift 0.2 °/s"),
        ("common_daq_freeze", "Common: acquisition freeze, both device channels"),
        ("common_daq_freeze_incl_ref", "Common: acquisition freeze incl. the second reading"),
        ("common_gain_m0p5pct", "Common: gain −0.5 % on both device channels"),
        ("pivot_stuck_onset_noise0", "Pivot stuck from onset, no noise ¹"),
        ("joint_stuck_during_noise0", "Joint stuck during the drop, no noise ¹"),
    ]
    frows = []
    for k, nm in flab:
        f = r2["faults"][k]
        frows.append(
            f"| {nm} | {f['missed']['p']:.1%} | {f['missed_base']['p']:.1%} | {f['missed_run']['p']:.1%} | "
            f"{f['missed_xchk']['p']:.1%} | {f['flagged_xchk']:.1%} |"
        )
    sep = [
        f"| {r['sigma_counts']:g} | {r['clean_max_run_ms']:.2f} | {r['stuck_min_run_ms']:.1f} | "
        f"{r['clean_flag_rate']:.1%} | {'yes' if r['separated'] else 'no'} |"
        for r in r2["run_sep"]
    ]
    slab = [
        ("base", "Base: δ 0.01°, T_h 30, freeze U[0, t_c], no noise"),
        ("base_noise1", "Base with 1 count noise"),
        ("delta_0.1", "δ 0.1°"),
        ("delta_1", "δ 1°"),
        ("delta_3", "δ 3°"),
        ("Th_15", "T_h 15 N·m"),
        ("Th_20", "T_h 20 N·m"),
        ("Th_50", "T_h 50 N·m"),
        ("window_t10", "Freeze U[0, 10° crossing]"),
        ("window_tend", "Freeze U[0, end of record]"),
    ]
    srows = []
    for k, nm in slab:
        v = r2["stuck_sens"][k]
        srows.append(
            f"| {nm} | {_pc(v['missed'])} | {_pc(v['after_a50k'])} | {_pc(v['after_a_tight'])} | "
            f"{v['share_gross_lt5']:.0%} |"
        )
    return dict(
        r2_budget_table="\n".join(bud),
        r2_options_table="\n".join(opt),
        r2_fault_table="\n".join(frows),
        r2_sep_table="\n".join(sep),
        r2_sens_table="\n".join(srows),
    )


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
