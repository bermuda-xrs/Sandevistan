"""T10 bench-protocol check (ANALYSIS ONLY, SYNTHETIC). Same S0 chain as latency_chain.py
(MAV 100 ms -> 7 ms delay -> 1st-order 30 Hz; no AAF term), 2 kHz.
Compares estimators of t50 (time to 50 % of the step change in torque):
  (a) single trial, band-limited Gaussian carrier (EMG-like), first crossing
  (b) ensemble mean of M trials (M = 10, 20), then first crossing of the averaged torque
  (c) deterministic 150 Hz sine-burst carrier (MAV ripple-free because 100 ms holds 30 half-periods),
      single trial
Torque-sensor noise SD = 2 % or 5 % of the step. Reports mean, SD, and fraction of estimates above 68.67 ms
(the protected gate value, used here only as a reference line).
T10 CRITERION (Charlie, decision A, 2026-11-23): "assist arrived" = torque t50 on the ENSEMBLE MEAN. PASS if
t50 <= 100 ms (RevC T10 limit, unchanged). t10 and t90 of the same ensemble-mean trace are reported and never
decide pass/fail. `assist_arrived` is that rule: trials (rows) are averaged first, then the mean trace is
crossed; it needs at least 10 repeats (20 preferred). Its results on the synthetic chain go to
t10_criterion.md/.json.
VALIDITY GUARD (FC PR #13 items #8/#9): a run is INVALID, so not passed, unless the ensemble-mean step is
positive, larger than K_NOISE x the baseline SD and larger than MIN_STEP, the t50 crossing is strictly after
T0, and the trace has settled by T0 + 350 ms (plateau drift within SETTLE_TOL of the step). See ``guard``.
The estimator comparison above and its outputs t10_protocol.md/.json are unchanged (same draws, same numbers).
Usage: python3 t10_protocol.py <out_dir>
"""

import json
import os
import sys

import numpy as np
from scipy import signal

FS = 2000.0
T0 = 0.4
TL = 1.0
SNR = 15.0
tm = np.arange(int(TL * FS)) / FS
BP = signal.butter(2, [20, 450], btype="bandpass", fs=FS, output="sos")
tau = 1 / (2 * np.pi * 30.0)
alpha = np.exp(-1 / (FS * tau))
n = int(0.1 * FS)
d = int(round(0.007 * FS))
rng = np.random.default_rng(7)
GATE = 68.67


def chain(emg):
    cs = np.cumsum(np.abs(emg), axis=1)
    env = np.empty_like(cs)
    env[:, :n] = cs[:, :n] / np.arange(1, n + 1)
    env[:, n:] = (cs[:, n:] - cs[:, :-n]) / n
    env = np.concatenate([np.repeat(env[:, :1], d, axis=1), env[:, :-d]], axis=1)
    return signal.lfilter([1 - alpha], [1, -alpha], env, axis=1, zi=env[:, :1] * alpha)[0]


T10_PASS_MS = 100.0  # RevC T10 limit (line 203), not changed here
MIN_REPEATS = 10
FRACTIONS = {"t10": 0.1, "t50": 0.5, "t90": 0.9}
BASE = slice(int(0.15 * FS), int(T0 * FS))  # baseline window, 0.15 s .. T0
PLAT = slice(int((T0 + 0.35) * FS), int((T0 + 0.55) * FS))  # plateau window, T0 + 350 .. 550 ms
PLAT_A = slice(PLAT.start, (PLAT.start + PLAT.stop) // 2)  # first and second half of the plateau window
PLAT_B = slice((PLAT.start + PLAT.stop) // 2, PLAT.stop)
# Validity guard. ASSUMPTION (not sourced): the three constants are engineering choices, justified here.
# K_NOISE: the step must exceed 10 baseline SDs, so the t50 threshold (base + step/2) is more than 5 SDs
#   above the baseline mean. For Gaussian noise the one-sided tail beyond 5 SD is 2.9e-7 per sample, so
#   noise alone reaches the t50 threshold with probability of at most 3.4e-4 over the 1200 samples after
#   T0 (union bound).
# MIN_STEP: absolute floor in the torque's own units, meant as N*m at the joint on the bench. 0.5 N*m is
#   the RevC T5 felt-torque limit for transparency (line 198): a step smaller than the residual torque the
#   user may already feel is not an assist. Synthetic chain units are arbitrary; its steps are far above it.
# SETTLE_TOL: the means of the two halves of the plateau window may differ by at most 15 % of the step;
#   otherwise the trace is still moving at T0 + 350 ms and t10/t50/t90 read early (FC #9 slow ramp:
#   drift 0.67 of the step). 15 % sits above the drift that MAV ripple alone gives on the synthetic chain
#   (max 0.081 over 4000 runs), so normal runs stay valid. It cannot catch slow first-order tails (drift
#   < 0.04): see the spec (docs/t10-assist-arrived.md) for the measured bias that remains.
K_NOISE = 10.0
MIN_STEP = 0.5
SETTLE_TOL = 0.15


def crossing(y, frac):
    """First time after T0 (ms) that each row reaches base + frac * step. NaN if it never does."""
    base = y[:, BASE].mean(axis=1)
    ss = y[:, PLAT].mean(axis=1)
    thr = base + frac * (ss - base)
    hit = (y >= thr[:, None]) & (tm >= T0)[None, :]
    return np.where(hit.any(axis=1), (tm[hit.argmax(axis=1)] - T0) * 1e3, np.nan)


def t50(y):
    return crossing(y, 0.5)


def guard(y, min_step=MIN_STEP):
    """Validity of each row (an ensemble-mean trace). Returns step, baseline SD, plateau drift and reasons."""
    base = y[:, BASE].mean(axis=1)
    step = y[:, PLAT].mean(axis=1) - base
    base_sd = y[:, BASE].std(axis=1, ddof=1)
    drift = y[:, PLAT_B].mean(axis=1) - y[:, PLAT_A].mean(axis=1)
    tt = {k: crossing(y, f) for k, f in FRACTIONS.items()}
    reasons = [[] for _ in range(len(y))]
    for i in range(len(y)):
        if not step[i] > K_NOISE * base_sd[i]:
            reasons[i].append(f"step {step[i]:.4g} not > {K_NOISE:g} x baseline SD {base_sd[i]:.4g}")
        if not step[i] > min_step:
            reasons[i].append(f"step {step[i]:.4g} not > floor {min_step:g}")
        if not tt["t50"][i] > 0:  # also catches NaN (no crossing)
            reasons[i].append(f"t50 crossing {tt['t50'][i]:.2f} ms is not after T0")
        if not abs(drift[i]) <= SETTLE_TOL * abs(step[i]):
            reasons[i].append(f"not settled by T0+350 ms: drift {drift[i]:.4g} > {SETTLE_TOL:g} x step")
    valid = np.array([not r for r in reasons])
    passed = valid & (tt["t50"] <= T10_PASS_MS)
    return dict(tt, step=step, base_sd=base_sd, drift=drift, reasons=reasons, valid=valid, passed=passed)


def assist_arrived(torque, min_step=MIN_STEP):
    """T10 result for one bench run: torque (repeats, samples), same time base as ``tm``.
    status is PASS, FAIL (t50 > 100 ms) or INVALID (guard failed); only PASS counts as passed."""
    torque = np.asarray(torque, float)
    if torque.ndim != 2 or len(torque) < MIN_REPEATS:
        raise ValueError(f"T10 needs at least {MIN_REPEATS} repeats")
    if torque.shape[1] != tm.size or not np.isfinite(torque).all():
        raise ValueError(f"T10 needs finite torque on the time base tm ({tm.size} samples, T0 = {T0} s)")
    g = guard(torque.mean(axis=0, keepdims=True), min_step)
    r = {k: float(g[k][0]) for k in FRACTIONS}
    status = "INVALID" if not g["valid"][0] else ("PASS" if g["passed"][0] else "FAIL")
    return dict(
        r,
        n=len(torque),
        assist_arrived_ms=r["t50"],
        passed=bool(g["passed"][0]),
        status=status,
        reasons=g["reasons"][0],
        step=float(g["step"][0]),
        baseline_sd=float(g["base_sd"][0]),
        plateau_drift=float(g["drift"][0]),
    )


amp = 1 + (SNR - 1) * (tm >= T0)


def noise_trials(N):
    car = signal.sosfilt(BP, rng.standard_normal((N, tm.size + 400)), axis=1)[:, 400:]
    car /= car[:, : int(0.3 * FS)].std(axis=1, keepdims=True)
    return car * amp


def sine_trials(N):
    ph = rng.uniform(0, 2 * np.pi, (N, 1))
    car = np.sqrt(2) * np.sin(2 * np.pi * 150 * tm[None, :] + ph)  # unit RMS like the noise carrier
    return car * amp


def main(OUT):
    os.makedirs(OUT, exist_ok=True)
    crit = {}
    lines = [
        "# t10_protocol.py output (ALL SYNTHETIC)",
        f"Step EMG amplitude x{SNR:.0f}; reference line {GATE} ms",
        "| estimator | sensor noise | n estimates | t50 mean ms | t50 SD ms | frac > 68.67 |",
        "|---|---|---|---|---|---|",
    ]
    res = {}
    for sn in (0.02, 0.05):
        for name, M, gen in (
            ("single trial, noise carrier", 1, noise_trials),
            ("mean of 10, noise carrier", 10, noise_trials),
            ("mean of 20, noise carrier", 20, noise_trials),
            ("single trial, 150 Hz sine carrier", 1, sine_trials),
        ):
            K = 1000
            ys = []
            for _ in range(M):
                y = chain(gen(K))
                step = (y[:, int((T0 + 0.35) * FS) :].mean(axis=1) - y[:, : int(T0 * FS)].mean(axis=1)).mean()
                ys.append(y + rng.standard_normal(y.shape) * sn * step)
            mean = np.mean(ys, axis=0)
            est = t50(mean)
            if M >= MIN_REPEATS:
                g = guard(mean)
                crit[f"{name}|{sn}"] = g
                # the rule on one run equals the vectorised one
                one = assist_arrived(np.array([y[0] for y in ys]))
                assert one["assist_arrived_ms"] == est[0] and one["passed"] == g["passed"][0]
                # guard margins on the synthetic runs (stdout only; the output files are unchanged)
                snr, drel = g["step"] / g["base_sd"], np.abs(g["drift"]) / g["step"]
                print(
                    f"guard {name} {sn}: invalid {int((~g['valid']).sum())}/{K}, "
                    f"min step {g['step'].min():.4f}, "
                    f"min step/baseline SD {snr.min():.1f} (need > {K_NOISE:g}), "
                    f"|drift|/step mean {drel.mean():.4f} SD {drel.std(ddof=1):.4f} "
                    f"max {drel.max():.4f} (need <= {SETTLE_TOL:g}), "
                    f"min t50 {np.nanmin(g['t50']):.2f} ms"
                )
            res[f"{name}|{sn}"] = dict(
                mean=float(est.mean()), sd=float(est.std(ddof=1)), frac_gt=float((est > GATE).mean())
            )
            lines.append(
                f"| {name} | {sn * 100:.0f} % | {K} | {est.mean():.2f} | {est.std(ddof=1):.2f} | {(est > GATE).mean():.3f} |"  # noqa: E501
            )
    print("\n".join(lines))
    open(os.path.join(OUT, "t10_protocol.md"), "w").write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(OUT, "t10_protocol.json"), "w"), indent=1)

    # T10 criterion on the synthetic chain: ensemble mean, t50 decides, t10/t90 reported, guard applied.
    L = [
        "# T10 criterion: assist arrived = t50 of the ensemble mean, pass if <= 100 ms (ALL SYNTHETIC)",
        "1000 simulated bench runs per row, noise carrier, S0 chain (RevC assumptions, no AAF). Mean (SD) in ms.",  # noqa: E501
        "| run | sensor noise | t10 (reported) | t50 (criterion) | t90 (reported) | runs that pass (t50 <= 100) | t90 > 100 | t50 > 68.67 (reference only) |",  # noqa: E501
        "|---|---|---|---|---|---|---|---|",
    ]
    js = {}
    for key, g in crit.items():
        name, sn = key.split("|")
        js[key] = {
            k: dict(
                mean=float(np.nanmean(g[k])),
                sd=float(np.nanstd(g[k], ddof=1)),
                n_nan=int(np.isnan(g[k]).sum()),
            )
            for k in FRACTIONS
        }
        js[key].update(
            n=int(g["t50"].size),
            n_pass=int(g["passed"].sum()),
            n_t90_gt_100=int((g["t90"] > T10_PASS_MS).sum()),
            n_t50_gt_gate=int((g["t50"] > GATE).sum()),
        )
        j = js[key]
        L.append(
            f"| {name} | {float(sn) * 100:.0f} % | "
            + " | ".join(f"{j[k]['mean']:.2f} ({j[k]['sd']:.2f})" for k in FRACTIONS)
            + f" | {j['n_pass']}/{j['n']} | {j['n_t90_gt_100']}/{j['n']} | {j['n_t50_gt_gate']}/{j['n']} |"
        )
    print("\n".join(L))
    open(os.path.join(OUT, "t10_criterion.md"), "w").write("\n".join(L) + "\n")
    json.dump(js, open(os.path.join(OUT, "t10_criterion.json"), "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "out")
