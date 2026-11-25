# ruff: noqa: E501  (long docstring and markdown table literals; content unchanged)
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
t50 <= 100 ms (RevC T10 limit, unchanged). t10 and t90 of the same ensemble-mean trace are reported and never decide
pass/fail. `assist_arrived` is that rule: trials (rows) are averaged first, then the mean trace is crossed; it needs
at least 10 repeats (20 preferred). Its results on the synthetic chain go to t10_criterion.md/.json.
The estimator comparison above and its outputs t10_protocol.md/.json are unchanged (same draws, same numbers).
Usage: python3 t10_protocol.py <out_dir>
"""

import json
import os
import sys

import numpy as np
from scipy import signal

OUT = sys.argv[1] if len(sys.argv) > 1 else "out"
os.makedirs(OUT, exist_ok=True)
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


def crossing(y, frac):
    """First time after T0 (ms) that each row reaches base + frac * step. NaN if it never does."""
    base = y[:, int(0.15 * FS) : int(T0 * FS)].mean(axis=1)
    ss = y[:, int((T0 + 0.35) * FS) : int((T0 + 0.55) * FS)].mean(axis=1)
    thr = base + frac * (ss - base)
    hit = (y >= thr[:, None]) & (tm >= T0)[None, :]
    return np.where(hit.any(axis=1), (tm[hit.argmax(axis=1)] - T0) * 1e3, np.nan)


def t50(y):
    return crossing(y, 0.5)


def assist_arrived(torque):
    """T10 result for one bench run: torque (repeats, samples), same time base as ``tm``."""
    torque = np.asarray(torque, float)
    if torque.ndim != 2 or len(torque) < MIN_REPEATS:
        raise ValueError(f"T10 needs at least {MIN_REPEATS} repeats")
    mean = torque.mean(axis=0, keepdims=True)
    r = {k: float(crossing(mean, f)[0]) for k, f in FRACTIONS.items()}
    return dict(r, n=len(torque), assist_arrived_ms=r["t50"], passed=bool(r["t50"] <= T10_PASS_MS))


amp = 1 + (SNR - 1) * (tm >= T0)


def noise_trials(N):
    car = signal.sosfilt(BP, rng.standard_normal((N, tm.size + 400)), axis=1)[:, 400:]
    car /= car[:, : int(0.3 * FS)].std(axis=1, keepdims=True)
    return car * amp


def sine_trials(N):
    ph = rng.uniform(0, 2 * np.pi, (N, 1))
    car = np.sqrt(2) * np.sin(2 * np.pi * 150 * tm[None, :] + ph)  # unit RMS like the noise carrier
    return car * amp


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
            crit[f"{name}|{sn}"] = {k: crossing(mean, f) for k, f in FRACTIONS.items()}
            one = assist_arrived(
                np.array([y[0] for y in ys])
            )  # the rule on one run equals the vectorised estimate
            assert one["assist_arrived_ms"] == est[0]
        res[f"{name}|{sn}"] = dict(
            mean=float(est.mean()), sd=float(est.std(ddof=1)), frac_gt=float((est > GATE).mean())
        )
        lines.append(
            f"| {name} | {sn * 100:.0f} % | {K} | {est.mean():.2f} | {est.std(ddof=1):.2f} | {(est > GATE).mean():.3f} |"
        )
print("\n".join(lines))
open(os.path.join(OUT, "t10_protocol.md"), "w").write("\n".join(lines) + "\n")
json.dump(res, open(os.path.join(OUT, "t10_protocol.json"), "w"), indent=1)

# T10 criterion on the synthetic chain: ensemble mean, t50 decides, t10/t90 reported.
L = [
    "# T10 criterion: assist arrived = t50 of the ensemble mean, pass if <= 100 ms (ALL SYNTHETIC)",
    "1000 simulated bench runs per row, noise carrier, S0 chain (RevC assumptions, no AAF). Mean (SD) in ms.",
    "| run | sensor noise | t10 (reported) | t50 (criterion) | t90 (reported) | runs that pass (t50 <= 100) | t90 > 100 | t50 > 68.67 (reference only) |",
    "|---|---|---|---|---|---|---|---|",
]
js = {}
for key, c in crit.items():
    name, sn = key.split("|")
    js[key] = {
        k: dict(mean=float(np.nanmean(x)), sd=float(np.nanstd(x, ddof=1)), n_nan=int(np.isnan(x).sum()))
        for k, x in c.items()
    }
    js[key].update(
        n=int(c["t50"].size),
        n_pass=int((c["t50"] <= T10_PASS_MS).sum()),
        n_t90_gt_100=int((c["t90"] > T10_PASS_MS).sum()),
        n_t50_gt_gate=int((c["t50"] > GATE).sum()),
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
