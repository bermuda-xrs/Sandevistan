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
RATE = {"b": {"estop": 20.0, "comm": 20.0, "emg": 50.0}, "iv": {"estop": 19.0, "comm": 19.0}}
TOL = 0.05
G_UNC = 0.05  # subject gravity moment known to +-5 % (ASSUMPTION, HE)
RES_MS = 0.1  # RevU row 14 timing resolution
TOL_MODEL_MS = 0.5  # dummy agreement with the model (ASSUMPTION, HE)
EMG_HZ, IMU_HZ = 2000, 150  # repo defaults, quoted only (protected, not changed)

# PRU09 (Pruszynski et al. 2009, Results, "Behavior and kinematics"), ms, read as mean and SD.
PRU = {"reversal, 1 N·m background": (193.0, 16.0), "reversal, 3 N·m background": (170.0, 22.0),
       "maximal excursion, 1.25 N·m step": (178.0, 22.0), "maximal excursion, 2.5 N·m step": (184.0, 15.0)}
FC_P = (0.20, 0.33)  # Fact Checker, S7_verdict.md row 12

# Hardware Engineer's and Fact Checker's shipped values, for comparison only.
REF = dict(window=(201.90484840282784, 203.7164677450151), n_allow=(59, 93, 124), subjects=(6, 17, 33),
           fc_p=(0.19864655721181523, 0.3257000686903164), h0_p_worst=0.4795370246007952)

M_PLAN, ICC_PLAN = 10, Fr("0.5")  # ASSUMPTIONS (RevU row 14; upper end of the earlier 6-33 range)
ICCS = [Fr(s) for s in ("0", "0.05", "0.1", "0.2", "0.3", "0.5", "0.8", "1")]
MS = (1, 2, 5, 10, 20)
P_ROWS = (0.001, 0.01, 0.02, 0.05)
LP_TOL = 1e-9


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
    return dict(J=J, G=G, G_human=Mr_h * a["g"], kT0=K_HI * T0_DEV, hold_nom=G - K_HI * T0_DEV)


def crossing(path, rate, sto_ms, hold, J=None, G=None, ang=DROP_LIM):
    """Time (ms) and speed (deg/s) at which the drop from 90 deg reaches ``ang`` deg."""
    A = arm()
    J, G = J or A["J"], G or A["G"]
    td, sto = T_DET[path] / 1e3, math.inf if sto_ms is None else sto_ms / 1e3

    def rhs(tq):
        return lambda t, x: [x[1], (K_HI * tq(t) + hold - G * math.cos(x[0])) / J]

    # x = [drop angle measured up from the start (negative = falling), speed]; gravity lever = cos(drop)
    def tq(t):
        return T0_DEV if t < td else max(T0_DEV - rate * (t - td), 0.0)

    edges = sorted({0.0, 0.4} | {e for e in (td, td + T0_DEV / rate, sto) if 0 < e < 0.4})
    y = [0.0, 0.0]

    def hit(t, x):
        return x[0] + math.radians(ang)

    hit.terminal = True
    for a, b in zip(edges[:-1], edges[1:]):
        f = rhs(tq if a < sto else (lambda t: 0.0))
        s = solve_ivp(f, (a, b), y, method="DOP853", rtol=1e-11, atol=1e-12, events=hit)
        if s.t_events[0].size:
            return float(s.t_events[0][0]) * 1e3, -math.degrees(s.y_events[0][0][1])
        y = s.y[:, -1]
    raise ValueError("no crossing")


def drop_at(t_ms, path, rate, sto_ms, hold):
    return brentq(lambda d: crossing(path, rate, sto_ms, hold, ang=d)[0] - t_ms, 0.5, 60.0, xtol=1e-10)


def cells(opt):
    for path, r in RATE[opt].items():
        for tol in (-TOL, 0.0, TOL):
            for s in STO[path]:
                yield path, r * (1 + tol), s


def pru_mean(hold):
    """FC's linear extrapolation of the PRU09 reversal time to ``hold`` N*m (outside the data)."""
    (t1, _), (t3, _) = list(PRU.values())[:2]
    return t1 + (t3 - t1) / 2 * (hold - 1)


def model() -> dict:
    A = arm()
    b = {c: crossing(*c, 0.0)[0] for c in cells("b")}
    iv = {c: crossing(*c, 0.0)[0] for c in cells("iv")}
    gov, gov_iv = min(b, key=b.get), min(iv, key=iv.get)
    lo, hi = b[gov], iv[gov_iv]
    w10 = crossing(*gov, 0.0)[1]
    t_stop, w_stop = crossing(*gov, 0.0, ang=STOP_DEG)
    out = dict(A, lo=lo, hi=hi, width=hi - lo, gov="|".join(map(str, gov)),
               gov_iv="|".join(map(str, gov_iv)),
               lo_comm=min(v for (p, _, _), v in b.items() if p == "comm"),
               lo_emg=min(v for (p, _, _), v in b.items() if p == "emg"),
               hi_estop=min(v for (p, _, _), v in iv.items() if p == "estop"),
               n_cells_b=len(b), n_cells_iv=len(iv),
               margin_ms=lo - CATCH_MS, drop200=drop_at(CATCH_MS, *gov, 0.0), w10=w10,
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
                              for m in (2, M_PLAN) for r in ICCS)
    fine = [Fr(i, 100) for i in range(101)]

    def covered(m):
        return all(p_total_le(pmf_two(P0, r, m, 0), subjects(N0, m, r), 0) <= alpha for r in fine)

    out["deff_ok_to_m"] = next(m for m in range(1, 200) if not covered(m + 1))
    out["allowed"] = {f"{c}|{float(r):g}": dict(deff=subjects(nc[c], M_PLAN, r),
                                                bb=subjects_model(pmf_bb, c, M_PLAN, r),
                                                two=subjects_model(pmf_two, c, M_PLAN, r))
                      for c in (0, 1, 2) for r in ICCS}
    out["allowed_short"] = [k for k, v in out["allowed"].items() if v["deff"] < max(v["bb"], v["two"])]
    out["ext_subjects"] = {f"{float(r):g}": [subjects(n, M_PLAN, r) for n in ext] for r in ICCS}
    out["ext_clustered_max"] = max(
        staged_pass_subjects(out["ext_subjects"][f"{float(r):g}"], f(P0, r, M_PLAN, 2))
        for r in ICCS for f in (pmf_bb, pmf_two))
    out["latent_icc"] = [0.25, 0.5, 0.75, 0.9]
    out["latent"] = [icc_binary(r) for r in out["latent_icc"]]
    out["range_icc"] = [0.0, 0.2, 0.5]
    k = subjects(N0, M_PLAN, ICC_PLAN)
    # plan sized for ICC 0 when the true ICC is ICC_PLAN: false-pass chance at p = P0
    k0 = subjects(N0, M_PLAN, Fr(0))
    out["misspec"] = [p_total_le(f(P0, ICC_PLAN, M_PLAN, 0), k0, 0) for f in (pmf_bb, pmf_two)]
    out["plan"] = dict(m=M_PLAN, icc=float(ICC_PLAN), subjects=k, trials=k * M_PLAN,
                       deff=float(deff(M_PLAN, ICC_PLAN)), subjects_any=subjects(N0, M_PLAN, Fr(1)),
                       trials_any=subjects(N0, M_PLAN, Fr(1)) * M_PLAN)
    return out


# ---------------------------------------------------------------- feasibility (SYNTHETIC)
def feasibility(mod) -> dict:
    ps = {"FC low": FC_P[0], "FC high": FC_P[1], "window (b) edge": mod["p_lo_same_mean"]}
    nc = [n_allowed(c) for c in range(3)]
    out = {"ps": ps, "indep": {k: [float(binom.cdf(c, n, p)) for c, n in enumerate(nc)]
                               for k, p in ps.items()}}
    out["clustered"] = {}
    for k, p in ps.items():
        for r in ICCS:
            n_s = subjects(N0, M_PLAN, r)
            out["clustered"][f"{k}|{float(r):g}"] = dict(
                subjects=n_s, bb=p_total_le(pmf_bb(p, r, M_PLAN, 0), n_s, 0),
                two=p_total_le(pmf_two(p, r, M_PLAN, 0), n_s, 0))
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
    out["first_drop"] = {k: 1 / p for k, p in ps.items()}
    out["drop_in_10"] = {k: 1 - (1 - p) ** 10 for k, p in ps.items()}
    out["drop_in_10_icc"] = {k: 1 - float(pmf_bb(p, ICC_PLAN, 10, 0)[0]) for k, p in ps.items()}
    out["drop_in_2_subj"] = {k: 1 - (1 - v) ** 2 for k, v in out["drop_in_10_icc"].items()}
    out["undersized"] = p_total_le(pmf_bb(FC_P[0], ICC_PLAN, M_PLAN, 0), subjects(N0, M_PLAN, Fr(0)), 0)
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
    mod, st = model(), stats()
    feas = feasibility(mod)
    inp = dict(win_lo=WINDOW_MS[0], win_hi=WINDOW_MS[1], p0=P0, conf=CONF, n0=N0, drop=DROP_LIM,
               catch=CATCH_MS, stop=STOP_DEG, res=RES_MS, rate=RATE["b"]["estop"] * (1 + TOL),
               rate_iv=RATE["iv"]["estop"] * (1 + TOL), sto=STO["estop"][0], sto_comm=STO["comm"][0],
               tol_model=TOL_MODEL_MS, emg_hz=EMG_HZ, imu_hz=IMU_HZ, imu_ms=1e3 / IMU_HZ,
               daq_khz=1 / RES_MS, lp_tol=LP_TOL, m=M_PLAN, icc=float(ICC_PLAN), k=K_HI, T0=T0_DEV, tol=TOL,
               g_unc=G_UNC, fc_lo=FC_P[0], fc_hi=FC_P[1], pru_sd=mod["fc_sd"])
    return {"NOTE": "keys under 'syn' are SYNTHETIC (model output, not bench or human data)",
            "in": inp, "stat": st, "syn": dict(mod, feas=feas), "ref": reference_check(mod, st)}


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
    t["check"] = ["| ICC | Design effect | Subjects, design effect (plan) | Subjects, beta-binomial "
                  "| Subjects, worst case | False-pass chance, plan size, worst case |", "|---" * 6 + "|"]
    t["check"] += [f"| {r} | {v['deff']:.2f} | {v['subjects']} | {v['bb']} | {v['worst']} | "
                   f"{_pct(v['false_pass_worst'])} |" for r in rs for v in [g[f"{M_PLAN}|{r}"]]]
    t["allowed"] = ["| Allowed drops c | n (independent) | Size at 5 % | "
                    + " | ".join(f"Subjects, ICC {r}" for r in rs) + " |", "|---" * (len(rs) + 3) + "|"]
    t["allowed"] += [f"| {c} | {st[f'n{c}']} | {100 * st['size'][c]:.3f} % | "
                     + " | ".join(str(al[f"{c}|{r}"]["deff"]) for r in rs) + " |" for c in (0, 1, 2)]
    t["ext"] = ["| Look | Pass if total drops at most | Total trials (independent) | "
                + " | ".join(f"Subjects, ICC {r}" for r in rs) + " |", "|---" * (len(rs) + 3) + "|"]
    t["ext"] += [f"| {i + 1} | {i} | {n} | "
                 + " | ".join(str(st["ext_subjects"][r][i]) for r in rs) + " |"
                 for i, n in enumerate(st["ext"])]
    t["oc"] = ["| True drop probability | Pass, c = 0 (n = %d) | Pass, c = 1 (n = %d) "
               "| Pass, c = 2 (n = %d) |" % (st["n0"], st["n1"], st["n2"]), "|---" * 4 + "|"]
    t["oc"] += [f"| {100 * r['p']:g} % | " + " | ".join(_pct(x) for x in r["P"]) + " |" for r in st["oc"]]
    t["feas"] = ["| SYNTHETIC true drop probability | c = 0, n = %d | c = 1, n = %d | c = 2, n = %d "
                 "| Mean trials to first drop | Drop within one subject's 10 trials, ICC 0 / %g |"
                 % (st["n0"], st["n1"], st["n2"], st["plan"]["icc"]), "|---" * 6 + "|"]
    t["feas"] += [f"| {100 * p:.0f} % ({k}) | " + " | ".join(f"{x:.1e}" for x in f["indep"][k])
                  + f" | {f['first_drop'][k]:.1f} | {100 * f['drop_in_10'][k]:.1f} % / "
                  f"{100 * f['drop_in_10_icc'][k]:.1f} % |" for k, p in f["ps"].items()]
    t["feas_icc"] = ["| SYNTHETIC: ICC | Subjects (trials) | " + " | ".join(
        f"{100 * p:.0f} %: beta-binomial / best case" for p in f["ps"].values()) + " |",
        "|---" * (2 + len(f["ps"])) + "|"]
    t["feas_icc"] += [f"| {r} | {v0['subjects']} ({v0['subjects'] * M_PLAN}) | " + " | ".join(
        f"{f['clustered'][f'{k}|{r}']['bb']:.1e} / {f['clustered'][f'{k}|{r}']['two']:.1e}"
        for k in f["ps"]) + " |" for r in rs for v0 in [f["clustered"][f"FC low|{r}"]]]
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
            "23": "Declaration of Helsinki paragraph", "90": "elbow start angle, deg"}
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
