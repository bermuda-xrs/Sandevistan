"""Row-12 kill test of the ML ramp-model EMG onset fit. SYNTHETIC DATA ONLY.

No recorded EMG or motion is read. Estimator: scripts/ml_ramp_onset.py. Nothing in src/ is changed and
no protected value is set here: the gate is LeadStats.summarize at its defaults (68.67 ms, 10% miss cap,
16x3 shape, proportional mode).

True onset (one definition, every shape): the FIRST time the noise-free amplitude crosses P = 0.10 of its
plateau excursion, (a - 1)/(G - 1) = 0.10. Every row-12 lead is true IMU onset minus THAT time. P was set
before any gate result; other P are scored post hoc from the same trials (``P_SENS``).

Shapes (x = (t - t0)/R for the nominal rise R; both have the same 10-90% span, 0.6084 R):
  smooth  a smoothstep of full duration R:  x^2 (3 - 2x)
  exp     1 - exp(-(t - t0)/tau), tau = 0.6084 R / ln 9

Design and seeds are the fixed scheme of scripts/row12_power_sim.py (commit eb339a1, branch
sim/row12-power-seedfix): ROOT_ENTROPY and spawn keys study -> subject -> 1+session -> 1+trial, with
children 0 (design draws), 1 (EMG noise), 2 (angle noise). The one change: the drawn EMG time U(1.5, 2.5)
is the P-crossing (the true onset), and the ramp start is placed before it. All cells share the noise and
the IMU signals (common random numbers), so cells are paired, not independent.

EMG path per trial and cell: the repo's band-pass, 8 Hz envelope, 6-SD threshold and 30 ms hold (batch
re-implementation, checked against EmgOnsetDetector, which also handles every non-trivial trial) decide
the trigger; the fit reads the same causal window as the repo detector (trigger - 0.4 s .. trigger + 0.1 s,
last 1000 samples). A trial with a rest failure, no trigger, or an onset before 1.0 s is missed.
IMU path: ImuOnsetDetector unchanged, once per (trial, true lead).

  python scripts/ml_ramp_row12_sim.py simulate --studies 0:25 --out out/ml/raw_000.npz [--procs 6]
  python scripts/ml_ramp_row12_sim.py score out/ml/raw_*.npz --out out/ml/scored.json
  python scripts/ml_ramp_row12_sim.py report out/ml/scored.json --out docs/ml-ramp-onset-tables.md
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402
from scipy import optimize, signal, stats  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / "src"), str(HERE)]

from ml_ramp_onset import ONSET_FRACTION, MlRamp, ar_whitener, filter_delay_s  # noqa: E402

from emg_onset_lead.detect import (  # noqa: E402
    BaselineRest,
    EmgOnsetDetector,
    EmgParams,
    ImuOnsetDetector,
    ImuParams,
)
from emg_onset_lead.harness import TrialResult  # noqa: E402
from emg_onset_lead.latency import LatencyConfig  # noqa: E402
from emg_onset_lead.stats import LeadStats  # noqa: E402

ROOT_ENTROPY = 20261122
SD_SUBJECT, SD_SESSION, SD_TRIAL = 8.0, 6.0, 10.0
N_SUB, N_SESSIONS, N_TRIALS = 16, 3, 10
FS, FS_IMU, T_LEN = 2000.0, 150.0, 4.0
N = int(T_LEN * FS)
GATE_MS = LatencyConfig().total_ms
LEADS = (74.67, 68.67)
P = ONSET_FRACTION
P_SENS = (0.0, 0.05, 0.10, 0.25, 0.50)
SHAPES, RISES, SNRS = ("smooth", "exp"), (10, 20, 30, 40, 80), (3, 5)
CELLS = [(s, r, q) for s in SHAPES for r in RISES for q in SNRS]
EP = EmgParams()
assert abs(GATE_MS - 68.67) < 1e-9 and EP.short_burst_peak == 1.6 and EP.hold_s == 0.03


def _smooth_x(p):
    return 0.0 if p <= 0 else optimize.brentq(lambda x: 3 * x * x - 2 * x**3 - p, 0, 1)


FRAC_1090 = _smooth_x(0.9) - _smooth_x(0.1)
TAU_PER_R = FRAC_1090 / np.log(9)


def q_cross(shape: str, p: float) -> float:
    """Time from ramp start to the first crossing of plateau fraction p, in units of R."""
    return _smooth_x(p) if shape == "smooth" else -TAU_PER_R * np.log1p(-p)


# ------------------------------------------------------------------ seeds (as row12_power_sim.py)
def node(*key):
    return np.random.SeedSequence(ROOT_ENTROPY, spawn_key=tuple(int(k) for k in key))


def rng_at(*key):
    return np.random.default_rng(node(*key))


def draw_design(study, n_sub=N_SUB):
    so, se, gain = np.empty(n_sub), np.empty((n_sub, N_SESSIONS)), np.empty((n_sub, N_SESSIONS))
    eps, ton = (np.empty((n_sub, N_SESSIONS, N_TRIALS)) for _ in range(2))
    for i in range(n_sub):
        so[i] = rng_at(study, i, 0).normal(0, SD_SUBJECT)
        for j in range(N_SESSIONS):
            r = rng_at(study, i, 1 + j, 0)
            se[i, j], gain[i, j] = r.normal(0, SD_SESSION), r.uniform(0.7, 1.3)
            for k in range(N_TRIALS):
                r = rng_at(study, i, 1 + j, 1 + k, 0)
                eps[i, j, k], ton[i, j, k] = r.normal(0, SD_TRIAL), r.uniform(1.5, 2.5)
    return dict(so=so, se=se, gain=gain, eps=eps, ton=ton)


# ------------------------------------------------------------------ EMG path
_BAND = [EP.band_hz[0], min(EP.band_hz[1], 0.45 * FS)]
BP = signal.butter(EP.order, _BAND, btype="bandpass", fs=FS, output="sos")
LP = signal.butter(EP.order, EP.envelope_hz, btype="lowpass", fs=FS, output="sos")
BASE_N, SETTLE_N = round(EP.baseline_s * FS), round(EP.settle_s * FS)
HOLD_N, REFINE_N, LOOK_N = round(EP.hold_s * FS), round(EP.refine_s * FS), round(EP.lookback_s * FS)
WIN = LOOK_N + REFINE_N
_imp = np.zeros(1000)
_imp[0] = 1.0
BP_IMP = signal.sosfilt(BP, _imp)
TE = np.arange(N) / FS


def stream_filter(sos, x):
    """Causal filter along the last axis, started like the repo's _sos_stream."""
    zi = signal.sosfilt_zi(sos)[:, None, :] * x[None, :, :1]
    return signal.sosfilt(sos, x, axis=-1, zi=zi)[0]


class _Probe(EmgOnsetDetector):
    """The repo detector, recording which trigger and fit end it settles on."""

    probe = None

    def _onset_samples(self, commit_idx):
        end = commit_idx if self._fit_end is None else min(commit_idx, self._fit_end)
        self.probe = (self._trig_idx, end)
        return float(commit_idx)  # any valid onset: the caller does its own fit


def slow_trigger(x):
    """(trigger, fit end) from the repo detector, or None when it misses / rejects the baseline."""
    det = _Probe(FS, EP)
    for i in range(0, x.size, 100):
        det.process(x[i:i + 100])
    return det.probe if det.result.onset_s is not None else None


def triggers(x, b, f):
    """Per row: (trigger, fit end) or None. x raw, b band-passed, f envelope, all (rows, N)."""
    out = []
    seg = f[:, SETTLE_N:BASE_N]
    mean, thr = seg.mean(axis=1), EP.k * np.maximum(seg.std(axis=1, ddof=1), 1e-12)
    above = (f[:, BASE_N:] - mean[:, None]) > thr[:, None]
    c = np.cumsum(above, axis=1, dtype=np.int32)
    run = c[:, HOLD_N - 1:].copy()
    run[:, 1:] -= c[:, :-HOLD_N]
    full = run == HOLD_N  # HOLD_N consecutive samples above, ending at BASE_N + HOLD_N - 1 + idx
    first = np.argmax(full, axis=1)
    for r in range(x.shape[0]):
        t = BASE_N + HOLD_N - 1 + first[r]
        commit = t + REFINE_N
        if not BaselineRest.at_rest(f[r], FS, EP.settle_s, EP.baseline_s, EP.rest_ratio):
            out.append(None)
        elif full[r, first[r]] and commit < N and above[r, commit - BASE_N]:
            out.append((int(t), int(commit)))
        else:
            out.append(slow_trigger(x[r]))
    return out


_ML = {}


def ml_fit(z):
    n = z.shape[1]
    if n not in _ML:
        _ML[n] = MlRamp(FS, n)
    return _ML[n].fit(z)


def emg_fit(x):
    """(t0_s, R_s, G) per row of raw EMG (rows, N); NaN rows are missed. Rows share one baseline."""
    b = stream_filter(BP, x)
    f = stream_filter(LP, np.abs(b))
    tr = triggers(x, b, f)
    out = np.full((x.shape[0], 3), np.nan)
    # every row carries the same baseline noise up to the earliest ramp start, so one whitener serves all
    A, sd = ar_whitener(b[0, SETTLE_N:BASE_N])
    w = signal.lfilter(A, [1.0], b, axis=1) / sd
    delay = filter_delay_s(BP_IMP, A, FS)
    groups = {}
    for r, t in enumerate(tr):
        if t is None:
            continue
        end = t[1] + 1
        w0 = max(BASE_N, t[0] - LOOK_N + 1, end - WIN)
        if end - w0 >= 64:
            groups.setdefault(end - w0, []).append((r, w0, end))
    for rows in groups.values():
        z = np.stack([w[r, a:e] for r, a, e in rows])
        t0, R, G, _ = ml_fit(z)
        for q, (r, a, _e) in enumerate(rows):
            ts = a / FS + t0[q] - delay
            if ts + P * R[q] >= EP.baseline_s:
                out[r] = ts, R[q], G[q]
    return out


def burst(shape, t, t0, R):
    if shape == "smooth":
        u = np.clip((t - t0) / R, 0, 1)
        return u * u * (3 - 2 * u)
    return np.where(t >= t0, -np.expm1(-(np.maximum(t - t0, 0)) / (TAU_PER_R * R)), 0.0)


def ramp_start(shape, R_s, ton):
    """Ramp start for a burst whose P-crossing is at ``ton``."""
    return ton - q_cross(shape, P) * R_s


def emg_rows(noise, ton, gain):
    rows = np.empty((len(CELLS), N))
    for c, (sh, r, snr) in enumerate(CELLS):
        R = r / 1e3
        rows[c] = (1 + (snr * gain - 1) * burst(sh, TE, ramp_start(sh, R, ton), R)) * noise
    return rows


def _minjerk(x):
    x = np.clip(x, 0, 1)
    return x**3 * (10 - 15 * x + 6 * x * x)


def imu_error_ms(noise, tm):
    """IMU onset error (ms) vs the true move start ``tm``; NaN if missed or rest-rejected."""
    ti = np.arange(int(T_LEN * FS_IMU)) / FS_IMU
    ang = 20 + 60 * _minjerk((ti - tm) / 0.5) + noise
    res = ImuOnsetDetector.run(ang, FS_IMU, ImuParams(), max(1, round(0.05 * FS_IMU)))
    ok = res.baseline_rest_ok and res.onset_s is not None
    return (res.onset_s - tm) * 1e3 if ok else np.nan


def simulate_study(study):
    d = draw_design(study)
    off = d["so"][:, None, None] + d["se"][:, :, None] + d["eps"]
    nc = len(CELLS)
    shp = (N_SUB, N_SESSIONS, N_TRIALS)
    fit = np.full((nc,) + shp + (3,), np.nan, np.float32)  # t0 error vs ramp start (ms), R (ms), G
    imu = np.full((len(LEADS),) + shp, np.nan, np.float32)
    R0 = np.array([r / 1e3 for _, r, _ in CELLS])
    for i in range(N_SUB):
        for j in range(N_SESSIONS):
            for k in range(N_TRIALS):
                ton, g = d["ton"][i, j, k], d["gain"][i, j]
                noise = rng_at(study, i, 1 + j, 1 + k, 1).normal(0, 1, N)
                est = emg_fit(emg_rows(noise, ton, g))
                start = np.array([ramp_start(sh, R, ton) for (sh, _r, _s), R in zip(CELLS, R0)])
                fit[:, i, j, k, 0] = (est[:, 0] - start) * 1e3
                fit[:, i, j, k, 1] = est[:, 1] * 1e3
                fit[:, i, j, k, 2] = est[:, 2]
                an = rng_at(study, i, 1 + j, 1 + k, 2).normal(0, 0.3, int(T_LEN * FS_IMU))
                for li, L in enumerate(LEADS):
                    imu[li, i, j, k] = imu_error_ms(an, ton + (L + off[i, j, k]) / 1e3)
    return dict(study=study, off=off, gain=d["gain"], fit=fit, imu=imu)


def cmd_simulate(a):
    lo, hi = (int(v) for v in a.studies.split(":"))
    t = time.time()
    with Pool(a.procs) as p:
        res = p.map(simulate_study, range(lo, hi), chunksize=1)
    out = {k: np.stack([r[k] for r in res]) for k in ("off", "gain", "fit", "imu")}
    out["studies"] = np.arange(lo, hi)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, **out)
    print(f"done {a.out} {hi - lo} studies in {time.time() - t:.0f} s", flush=True)


# ------------------------------------------------------------------ scoring
def emg_error_ms(fit, shape, r_ms, p):
    """EMG onset error (ms) at plateau fraction p: fitted crossing minus true first crossing."""
    return fit[..., 0] + p * fit[..., 1] - q_cross(shape, p) * r_ms


def gate(lead, true_ms):
    """Project gate on one study; lead (16, 3, 10) -> (pass, lower bound, missed fraction)."""
    res = [TrialResult(f"S{i}", f"S{i}-D{j}", f"T{k}", np.nan, np.nan, float(lead[i, j, k]),
                       float(lead[i, j, k]) - GATE_MS, float(true_ms[i, j, k]))
           for i in range(N_SUB) for j in range(N_SESSIONS) for k in range(N_TRIALS)]
    s = LeadStats.summarize(res, n_boot=1)
    lb = s["pooled"]["lower_bound_ms"]
    return bool(s["pass"]), float("nan") if lb is None else float(lb), s["pooled"]["miss_fraction"]


def score_study(Z):
    out = dict(study=int(Z["studies"]))
    for c, (sh, r, snr) in enumerate(CELLS):
        fit = Z["fit"][c].astype(float)
        for p in P_SENS:
            ee = emg_error_ms(fit, sh, r, p)
            for li, L in enumerate(LEADS):
                tr = L + Z["off"]
                lead = tr + Z["imu"][li] - ee
                ok, lb, miss = gate(lead, tr)
                e = lead - tr
                out[f"{sh}_{r}_{snr}_{L}_{p}"] = (ok, float(np.nanmean(e)), float(np.nanstd(e)), lb, miss,
                                                  float(np.nanmean(ee)), float(np.nanstd(ee)))
    return out


def cmd_score(a):
    rows = []
    for f in sorted(a.raw):
        Z = dict(np.load(f))
        per = [{k: v[i] for k, v in Z.items()} for i in range(len(Z["studies"]))]
        with Pool(a.procs) as p:
            rows += p.map(score_study, per, chunksize=1)
        print(f, len(rows), flush=True)
    ids = [r["study"] for r in rows]
    assert len(ids) == len(set(ids)), "duplicate studies"
    json.dump(sorted(rows, key=lambda r: r["study"]), open(a.out, "w"))


# ------------------------------------------------------------------ report
def cp(k, n, alpha=0.05):
    lo = 0.0 if k == 0 else stats.beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    return float(lo), float(hi)


def cp_lower_1s(k, n):
    return 0.0 if k == 0 else float(stats.beta.ppf(0.05, k, n - k + 1))


def fmt(k, n):
    lo, hi = cp(k, n)
    return f"{k}/{n} = {k / n:.3f} ({lo:.3f}–{hi:.3f})"


def cmd_report(a):
    S = json.load(open(a.scored))
    n = len(S)

    def col(sh, r, snr, L, p, i):
        return np.array([s[f"{sh}_{r}_{snr}_{L}_{p}"][i] for s in S], float)

    out = [f"Studies per cell: {n}, each {N_SUB}x{N_SESSIONS}x{N_TRIALS}. ALL SYNTHETIC. "
           f"True onset = first crossing of {P:.2f} of the plateau excursion, every shape.", "",
           "### Primary table (P = 0.10)", "",
           "Lead error = measured minus true lead (positive = reads high = unsafe). Bias = mean over studies "
           "of the study-mean error, with its Monte-Carlo SE. Spread = SD across studies of the study-mean "
           "error. Both use the 74.67 ms studies. EMG columns are the onset-estimate error alone (estimated "
           "minus true onset; an early estimate is negative and makes the lead read high). Rare IMU noise "
           "triggers (hundreds of ms early) are left in and enter the lead columns only. Pass rates carry "
           "two-sided 95% Clopper-Pearson intervals; 'CP 1-s LB' is the one-sided 95% lower bound the kill "
           "rule uses.", "",
           "| Shape | Rise (ms) | SNR | Lead bias (ms) | Spread across studies (ms) | EMG onset bias (ms) | "
           "EMG spread across studies (ms) | EMG per-trial SD (ms) | Missed trials | "
           "Power @ 74.67 | False pass @ 68.67 | CP 1-s LB | Kill? |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    killed, worst = [], []
    for sh, r, snr in CELLS:
        e = col(sh, r, snr, LEADS[0], P, 1)
        ee = col(sh, r, snr, LEADS[0], P, 5)
        sd = col(sh, r, snr, LEADS[0], P, 6)
        m = col(sh, r, snr, LEADS[0], P, 4)
        kp = int(col(sh, r, snr, LEADS[0], P, 0).sum())
        kf = int(col(sh, r, snr, LEADS[1], P, 0).sum())
        lb = cp_lower_1s(kf, n)
        if lb > 0.05:
            killed.append((sh, r, snr, kf, lb))
        worst.append((kf, sh, r, snr))
        out.append(f"| {sh} | {r} | {snr} | {e.mean():+.2f} ± {e.std(ddof=1) / np.sqrt(n):.2f} | "
                   f"{e.std(ddof=1):.2f} | {ee.mean():+.2f} | {ee.std(ddof=1):.2f} | "
                   f"{np.sqrt(np.mean(sd**2)):.2f} | {m.mean():.4f} | {fmt(kp, n)} | "
                   f"{fmt(kf, n)} | {lb:.4f} | {'**KILL**' if lb > 0.05 else 'no'} |")
    out += ["", f"**Verdict by the stated rule (any cell with CP one-sided 95% LB > 0.05): "
            f"{'KILL' if killed else 'SURVIVES'}.** {len(killed)} of {len(CELLS)} cells exceed."]
    if killed:
        out.append("Cells over the line: " + "; ".join(
            f"{sh} {r} ms SNR {q}: {k}/{n}, LB {b:.4f}" for sh, r, q, k, b in killed) + ".")
    out += ["", "### Sensitivity to the onset definition (same trials, true lead and IMU onset kept "
            "relative to each definition's own crossing; not the kill test)", "",
            "| Shape | Rise | SNR | "
            + " | ".join(f"P={p:g}: bias / power / false pass" for p in P_SENS) + " |",
            "|---|---|---|" + "---|" * len(P_SENS)]
    for sh, r, snr in CELLS:
        cells = []
        for p in P_SENS:
            b = col(sh, r, snr, LEADS[0], p, 1).mean()
            kp = int(col(sh, r, snr, LEADS[0], p, 0).sum())
            kf = int(col(sh, r, snr, LEADS[1], p, 0).sum())
            cells.append(f"{b:+.2f} / {kp / n:.3f} / {kf / n:.3f}")
        out.append(f"| {sh} | {r} | {snr} | " + " | ".join(cells) + " |")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text("\n".join(out) + "\n")
    print("\n".join(out))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("simulate")
    s.add_argument("--studies", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--procs", type=int, default=6)
    s = sub.add_parser("score")
    s.add_argument("raw", nargs="+")
    s.add_argument("--out", required=True)
    s.add_argument("--procs", type=int, default=4)
    s = sub.add_parser("report")
    s.add_argument("scored")
    s.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"simulate": cmd_simulate, "score": cmd_score, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
