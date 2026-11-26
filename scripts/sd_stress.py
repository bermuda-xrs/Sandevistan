"""Row-12 power under other variance components (SD stress runs). SYNTHETIC DATA ONLY.

No recorded EMG or motion is read. The design values stay 8/6/10 ms (subject/session/trial SD); this
script only asks what happens to the row-12 gate if the real spread is larger. Signals, seed tree and
detector calls come from ``row12_power_sim.py`` unchanged. The gate (68.67 ms, proportional), the 10%
miss cap and PLANNED_SHAPE are read from the package and never set here; as in row12_power_sim.py the
shape gate follows the simulated size, ``planned=(n, 3)``.

Fixed before any result was computed (the commit that adds this file precedes the run):
  - CONFIGS: 12/4/10, then subject SD 8, 10, 12, 14, 16 with session/trial at the design 6/10.
  - N_STUDIES = 1000 independent 16x3x10 studies per cell, study ids 1000-1999 (row12_power_sim.py
    used 0-249, so nothing is reused). 16000 subjects and 480000 trials per cell.
  - Subject counts other than 16 come from cutting the same 16000 subjects, in id order, into
    floor(16000/n) disjoint studies of n subjects, for every n in N_GRID.
  - True lead 74.67 ms for power (gate + 6 ms margin), exactly 68.67 ms for false pass.
  - "with bias": the detected lead as the gate sees it. "bias removed": the same leads minus the pooled
    mean detector error (detected - true, over detected trials) of that config, rise and true lead. The
    detector noise stays. That needs the true lead, so it is an oracle reference and not a correction
    anyone could apply to recorded data. "true lead": the gate on the true leads (no detector).
  - Subjects for 0.8 power: probit-in-sqrt(n) binomial ML fit over all of N_GRID (the row12 method).
    95% CI = 2.5/97.5 percentiles of N_BOOT resamples of subjects (each resample is cut into studies
    again). "not reachable" when the mean lead is at or below the gate. A value above max(N_GRID) is
    an extrapolation and is marked. The cross-check is noncentral t with the measured mean and SD of
    subject means.
  - 95% Clopper-Pearson intervals on every proportion. No interim look; chunks are scored only when
    all of them exist.

Cells are not independent: every config uses the same standard-normal draws (scaled), EMG noise and
angle noise, and every rise uses the same angle detections. Only the 1000 studies within a cell are
independent of each other.

Usage (repo root):
  python scripts/sd_stress.py run --raw out/sd_stress/raw --procs 6     # restartable, skips done chunks
  python scripts/sd_stress.py score --raw out/sd_stress/raw --out out/sd_stress/scored.json --procs 6
  python scripts/sd_stress.py report out/sd_stress/scored.json --out out/sd_stress/tables.md
  python scripts/sd_stress.py sensitivity --raw out/sd_stress/raw --out out/sd_stress/sensitivity.md
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import row12_power_sim as sim  # noqa: E402  (also puts src/ on the path)

from emg_onset_lead.detect import EmgParams  # noqa: E402
from emg_onset_lead.harness import TrialResult  # noqa: E402
from emg_onset_lead.stats import PLANNED_SHAPE, LeadStats, SampleSize  # noqa: E402

CONFIGS = ((12, 4, 10), (8, 6, 10), (10, 6, 10), (12, 6, 10), (14, 6, 10), (16, 6, 10))
DESIGN = (8, 6, 10)
RISES = sim.ADDENDUM_RISES  # 10, 20, 30, 40, 80 ms
LEADS = (sim.TRUE_POWER, sim.TRUE_NULL)
N_STUDIES, STUDY0, CHUNK = 1000, 1000, 25
N_SUB = PLANNED_SHAPE[0]
N_SESSIONS, N_TRIALS = sim.N_SESSIONS, sim.N_TRIALS
N_GRID = (8, 12, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 128, 160, 192, 256, 320, 400, 512)
MAX_MISSED = 0.10  # LeadStats.summarize default; tests check the two gates agree
N_BOOT, BOOT_SEED = 1000, 20261123
POWER_TARGET = 0.8
# docs/row12_power_sim.md, "Baseline five-row table", 250 studies at 8/6/10, study ids 0-249
ROW12_DOC = {"power": {10: 209, 20: 167, 30: 100, 40: 59, 80: 3},
             "false_pass": {10: 12, 20: 6, 30: 0, 40: 0, 80: 0}, "n": 250}


def cfg_name(c) -> str:
    return "/".join(str(v) for v in c)


# ---------------------------------------------------------------- simulation
def scaled_design(study: int, cfg, n_sub: int = N_SUB) -> dict:
    """row12 draws rescaled from 8/6/10 to ``cfg`` (same standard normals, same gain and EMG onset)."""
    d = sim.draw_design(study, n_sub)
    return dict(d, so=d["so"] * (cfg[0] / sim.SD_SUBJECT), se=d["se"] * (cfg[1] / sim.SD_SESSION),
                eps=d["eps"] * (cfg[2] / sim.SD_TRIAL))


def _onset(res) -> float:
    return res.onset_s if res.baseline_rest_ok and res.onset_s is not None else np.nan


def simulate_subject(key) -> dict:
    """Onset times (s) of one subject: EMG per rise, IMU per config and true lead. NaN = not usable."""
    study, i = key
    ds = [{k: v[i:i + 1] for k, v in scaled_design(study, c, i + 1).items()} for c in CONFIGS]
    offs = [sim.offsets(d)[0] for d in ds]
    ton, gain = ds[0]["ton"][0], ds[0]["gain"][0]
    emg = np.full((len(RISES), N_SESSIONS, N_TRIALS), np.nan)
    imu = np.full((len(CONFIGS), len(LEADS), N_SESSIONS, N_TRIALS), np.nan)
    for j in range(N_SESSIONS):
        for k in range(N_TRIALS):
            for ri, r in enumerate(RISES):
                x = sim.emg_signal(study, i, j, k, r, ton[j, k], gain[j])
                emg[ri, j, k] = _onset(sim.detect_emg(x))
            for ci in range(len(CONFIGS)):
                for li, L in enumerate(LEADS):
                    a = sim.angle_signal(study, i, j, k, ton[j, k], L + offs[ci][j, k])
                    imu[ci, li, j, k] = _onset(sim.detect_imu(a))
    return dict(emg=emg, imu=imu, off=np.stack(offs))


def _head() -> str:
    p = subprocess.run(["git", "-C", str(Path(__file__).resolve().parent), "describe", "--always", "--dirty",
                        "--abbrev=40"], capture_output=True, text=True)
    return p.stdout.strip() or "unknown"


def cmd_run(a):
    assert EmgParams().short_burst_peak == 1.6 and abs(sim.GATE_MS - 68.67) < 1e-9
    assert (sim.N_SESSIONS, N_SUB) == (PLANNED_SHAPE[1], 16)
    a.raw.mkdir(parents=True, exist_ok=True)
    head = _head()
    with Pool(a.procs) as p:
        for lo in range(a.first, a.first + a.studies, CHUNK):
            ids = list(range(lo, min(lo + CHUNK, a.first + a.studies)))
            f = a.raw / f"raw_{lo:06d}.npz"
            if f.exists():
                continue
            res = p.map(simulate_subject, [(s, i) for s in ids for i in range(N_SUB)], chunksize=1)
            shape = (len(ids), N_SUB)
            arrs = {k: np.stack([r[k] for r in res]).reshape(shape + res[0][k].shape) for k in res[0]}
            tmp = f.with_suffix(".tmp.npz")
            np.savez_compressed(tmp, studies=np.array(ids), head=head, numpy=np.__version__, **arrs)
            os.replace(tmp, f)
            print("chunk", f.name, len(ids), "studies", flush=True)
    print("ALLDONE", flush=True)


# ---------------------------------------------------------------- gate on subject means
def subject_stats(lead: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """lead (S, sessions, trials) -> subject means (mean of session means) and missed trials per subject."""
    ok = np.isfinite(lead)
    cnt = ok.sum(-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        sess = np.where(cnt > 0, np.where(ok, lead, 0).sum(-1) / cnt, np.nan)
        n_sess = np.isfinite(sess).sum(-1)
        mean = np.where(n_sess > 0, np.nansum(sess, -1) / n_sess, np.nan)
    return mean, (~ok).sum(axis=(1, 2))


def gate_blocks(mean: np.ndarray, miss: np.ndarray, n: int) -> np.ndarray:
    """Pass/fail of floor(S/n) disjoint n-subject studies: the LeadStats.summarize gate, vectorised."""
    k = mean.size // n
    m, q = mean[:k * n].reshape(k, n), miss[:k * n].reshape(k, n)
    c = np.isfinite(m).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = np.nansum(m, 1) / c
        sd = np.sqrt(np.nansum((m - mu[:, None]) ** 2, 1) / (c - 1))
        lb = mu - stats.t.ppf(0.95, c - 1) * sd / np.sqrt(c)
    return (c >= 2) & (lb > sim.GATE_MS) & (q.sum(1) / (n * N_SESSIONS * N_TRIALS) <= MAX_MISSED)


def gate_summarize(lead: np.ndarray) -> bool:
    """LeadStats.summarize on one study, lead (n, sessions, trials). The reference for ``gate_blocks``."""
    n = lead.shape[0]
    res = [TrialResult(f"S{i}", f"S{i}-D{j}", f"T{k}", np.nan, np.nan, float(lead[i, j, k]),
                       float(lead[i, j, k]) - sim.GATE_MS)
           for i in range(n) for j in range(N_SESSIONS) for k in range(N_TRIALS)]
    return bool(LeadStats.summarize(res, n_boot=1, planned=(n, N_SESSIONS))["pass"])


# ---------------------------------------------------------------- subjects for 0.8 power
def fit_probit(ns, ks, ms) -> np.ndarray | None:
    """Binomial ML of P(pass) = Phi(b0 + b1*sqrt(n)) by Fisher scoring; None if it has no finite fit."""
    ks, ms = np.asarray(ks, float), np.asarray(ms, float)
    if ks.sum() == 0 or ks.sum() == ms.sum():
        return None
    X = np.column_stack([np.ones(len(ns)), np.sqrt(np.asarray(ns, float))])
    p0 = np.clip((ks + 0.5) / (ms + 1), 1e-3, 1 - 1e-3)
    w = np.sqrt(ms * p0 * (1 - p0))
    b = np.linalg.lstsq(X * w[:, None], stats.norm.ppf(p0) * w, rcond=None)[0]
    ll = -np.inf
    for _ in range(100):
        eta = np.clip(X @ b, -8, 8)
        p, phi = stats.norm.cdf(eta), stats.norm.pdf(eta)
        v = np.maximum(p * (1 - p), 1e-15)
        W = ms * phi**2 / v
        try:
            step = np.linalg.solve(X.T @ (X * W[:, None]), X.T @ (phi * (ks - ms * p) / v))
        except np.linalg.LinAlgError:
            return None
        new, scale = None, 1.0
        for _ in range(30):  # step halving
            e = np.clip(X @ (b + scale * step), -8, 8)
            q = np.clip(stats.norm.cdf(e), 1e-15, 1 - 1e-15)
            new = float(np.sum(ks * np.log(q) + (ms - ks) * np.log1p(-q)))
            if new >= ll - 1e-12:
                break
            scale /= 2
        b = b + scale * step
        if abs(new - ll) < 1e-10:
            break
        ll = new
    return b


def crossing(b) -> float:
    """n at POWER_TARGET; inf if the fitted curve does not rise with n."""
    if b is None or b[1] <= 0:
        return float("inf")
    root = (stats.norm.ppf(POWER_TARGET) - b[0]) / b[1]
    return float(root**2) if root > 0 else 0.0


def nct_n(effect: float, sd: float, max_n: int = 5000) -> int | None:
    """Smallest n whose noncentral-t power is at least POWER_TARGET."""
    if effect <= 0:
        return None
    return SampleSize.required_subjects(sd**2, 0.0, effect, 0.05, POWER_TARGET, 1, max_n)


def sweep(mean, miss) -> tuple[list[int], list[int]]:
    ps = [gate_blocks(mean, miss, n) for n in N_GRID]
    return [int(p.sum()) for p in ps], [int(p.size) for p in ps]


def subjects_needed(mean, miss) -> dict:
    """Crossing, subject-bootstrap CI and the noncentral-t cross-check for one power cell."""
    effect, sd = float(np.nanmean(mean) - sim.GATE_MS), float(np.nanstd(mean, ddof=1))
    ks, ms = sweep(mean, miss)
    out = dict(ks=ks, ms=ms, effect=effect, sd=sd, nct=nct_n(effect, sd), reachable=effect > 0)
    if effect <= 0:
        return out
    b = fit_probit(N_GRID, ks, ms)
    x = crossing(b)
    rng = np.random.default_rng(BOOT_SEED)
    boot = np.empty(N_BOOT)
    for q in range(N_BOOT):
        pick = rng.integers(0, mean.size, mean.size)
        boot[q] = crossing(fit_probit(N_GRID, *sweep(mean[pick], miss[pick])))
    out["n80"], out["b"] = x, None if b is None else [float(v) for v in b]
    out["lo"], out["hi"] = (float(v) for v in np.percentile(boot, [2.5, 97.5], method="nearest"))
    out["boot_not_rising"] = int(np.isinf(boot).sum())
    if b is not None and min(ms) > 0:
        pr = np.clip(stats.norm.cdf(b[0] + b[1] * np.sqrt(N_GRID)), 1e-12, 1 - 1e-12)
        out["chi2"] = float(np.sum((np.array(ks) - np.array(ms) * pr) ** 2 / (np.array(ms) * pr * (1 - pr))))
    if np.isfinite(x) and x <= max(N_GRID):
        n = max(2, int(np.ceil(x)))
        p = gate_blocks(mean, miss, n)
        out["check"] = dict(n=n, k=int(p.sum()), m=int(p.size))
    return out


# ---------------------------------------------------------------- scoring
def load(raw: Path, studies: int, first: int) -> dict:
    files = sorted(raw.glob("raw_*.npz"))
    files = [f for f in files if not f.name.endswith(".tmp.npz")]
    Z = [dict(np.load(f)) for f in files]
    ids = np.concatenate([z["studies"] for z in Z])
    assert ids.tolist() == list(range(first, first + studies)), "missing, duplicate or unexpected studies"
    out = {k: np.concatenate([z[k] for z in Z]) for k in ("emg", "imu", "off")}
    out = {k: v.reshape((-1,) + v.shape[2:]) for k, v in out.items()}  # (subject, ...), id order
    out["heads"] = sorted({str(z["head"]) for z in Z})
    out["numpy"] = sorted({str(z["numpy"]) for z in Z})
    return out


_D: dict = {}


def score_cell(key) -> dict:
    ci, ri, li = key
    L = LEADS[li]
    true = L + _D["off"][:, ci]
    lead = (_D["imu"][:, ci, li] - _D["emg"][:, ri]) * 1000.0
    err = lead - true
    ok = np.isfinite(err)
    bias = float(err[ok].mean())
    med = float(np.median(err[ok]))
    out = dict(cfg=ci, rise=RISES[ri], lead=L, n_trials=int(err.size), n_miss=int((~ok).sum()), bias=bias,
               median_err=med, err_sd=float(err[ok].std(ddof=1)),
               n_early=int((err[ok] < -100).sum()), n_late=int((err[ok] > 100).sum()))
    for tag, x in (("raw", lead), ("unb", lead - bias)):
        m, q = subject_stats(x)
        p = gate_blocks(m, q, N_SUB)
        xs = x.reshape((-1, N_SUB, N_SESSIONS, N_TRIALS))
        ref = [gate_summarize(v) for v in xs]
        out[tag] = dict(k=int(p.sum()), m=int(p.size), vs_summarize=int(np.sum(p != np.array(ref))),
                        mean=float(np.nanmean(m)), sd=float(np.nanstd(m, ddof=1)),
                        n_subj_missing=int((~np.isfinite(m)).sum()))
        if L == sim.TRUE_POWER:
            out[tag]["need"] = subjects_needed(m, q)
    return out


def score_true(key) -> dict:
    """The gate on the true leads (no detector, no rise)."""
    ci, li = key
    m, q = subject_stats(LEADS[li] + _D["off"][:, ci])
    p = gate_blocks(m, q, N_SUB)
    out = dict(cfg=ci, lead=LEADS[li], k=int(p.sum()), m=int(p.size), sd=float(m.std(ddof=1)),
               mean=float(m.mean()))
    if LEADS[li] == sim.TRUE_POWER:
        out["need"] = subjects_needed(m, q)
    return out


def realised(ci: int) -> dict:
    """Realised spread of the true leads of one config (method-of-moments from the nested layout)."""
    off = _D["off"][:, ci]
    s = off.mean(-1)  # session means
    v_tr = off.var(-1, ddof=1).mean()
    v_se = s.var(-1, ddof=1).mean() - v_tr / N_TRIALS
    m = off.mean((1, 2))
    v_su = m.var(ddof=1) - v_se / N_SESSIONS - v_tr / (N_SESSIONS * N_TRIALS)
    sd, n = float(m.std(ddof=1)), m.size
    lo, hi = (sd * np.sqrt((n - 1) / stats.chi2.ppf(q, n - 1)) for q in (0.975, 0.025))
    return dict(subject=float(np.sqrt(v_su)), session=float(np.sqrt(v_se)), trial=float(np.sqrt(v_tr)),
                sd_mean=sd, lo=float(lo), hi=float(hi), mean=float(m.mean()), n=int(n),
                design=SampleSize.design_subject_sd(*CONFIGS[ci], N_SESSIONS, N_TRIALS))


def cmd_score(a):
    _D.update(load(a.raw, a.studies, a.first))
    keys = [(c, r, i) for c in range(len(CONFIGS)) for r in range(len(RISES)) for i in range(len(LEADS))]
    with Pool(a.procs) as p:
        cells = p.map(score_cell, keys, chunksize=1)
        true = p.map(score_true, [(c, i) for c in range(len(CONFIGS)) for i in range(len(LEADS))],
                     chunksize=1)
    out = dict(synthetic=True, n_studies=a.studies, first=a.first, configs=CONFIGS, n_grid=N_GRID,
               n_boot=N_BOOT, gate_ms=sim.GATE_MS, sim_head=_D["heads"], score_head=_head(),
               numpy=_D["numpy"],
               realised=[realised(c) for c in range(len(CONFIGS))], cells=cells, true=true)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=1))
    print("done", a.out)


# ---------------------------------------------------------------- report
def _num(x: float) -> str:
    return f"> {max(N_GRID)}" if not np.isfinite(x) or x > max(N_GRID) else f"{x:.1f}"


def _fmt(k: int, n: int) -> str:
    return sim.fmt(k, n) if n else "—"


def _need(d: dict) -> str:
    if not d["reachable"]:
        return f"not reachable (mean lead {d['effect']:+.2f} ms vs gate)"
    if np.isfinite(d["n80"]) and d["n80"] <= max(N_GRID):
        s = f"{d['n80']:.1f} ({_num(d['lo'])}–{_num(d['hi'])})"
    else:
        s = f"> {max(N_GRID)} (fit {d['n80']:.0f}, extrapolated)"
    return s + f"; nct {d['nct'] if d['nct'] is not None else 'none'}"


def _summary(D: dict, C: dict, T: dict) -> list[str]:
    """Compact config x rise tables of the headline numbers (the long tables follow)."""
    P, Q = LEADS
    head = "| SD subject/session/trial | " + " | ".join(f"{r} ms rise" for r in RISES) + " | no detector |"
    rule = "|---|" + "---|" * (len(RISES) + 1)

    def short(k, m):
        return f"{k / m:.3f}"

    def need(d):
        if not d["reachable"]:
            return "not reachable"
        if not (np.isfinite(d["n80"]) and d["n80"] <= max(N_GRID)):
            return f"> {max(N_GRID)}"
        return f"{d['n80']:.1f} ({_num(d['lo'])}–{_num(d['hi'])})"

    L = []
    for title, fn in (
        (f"Power at {N_SUB} subjects, with bias (true lead {P} ms)",
         lambda c, r: short(C[(c, r, P)]["raw"]["k"], C[(c, r, P)]["raw"]["m"])),
        (f"Power at {N_SUB} subjects, bias removed (oracle)",
         lambda c, r: short(C[(c, r, P)]["unb"]["k"], C[(c, r, P)]["unb"]["m"])),
        (f"Subjects for {POWER_TARGET} power, with bias (95% CI)",
         lambda c, r: need(C[(c, r, P)]["raw"]["need"])),
        (f"Subjects for {POWER_TARGET} power, bias removed (95% CI)",
         lambda c, r: need(C[(c, r, P)]["unb"]["need"])),
        (f"False pass at {Q} ms, with bias (nominal 0.05)",
         lambda c, r: short(C[(c, r, Q)]["raw"]["k"], C[(c, r, Q)]["raw"]["m"])),
        (f"False pass at {Q} ms, bias removed (oracle)",
         lambda c, r: short(C[(c, r, Q)]["unb"]["k"], C[(c, r, Q)]["unb"]["m"])),
    ):
        L += ["", f"### Summary: {title} (SYNTHETIC)", "", head, rule]
        for c in D["configs"]:
            c = tuple(c)
            row = [fn(c, r) for r in RISES]
            if "Power at" in title and "removed" in title:
                tt = T[(c, P)]
                row.append(short(tt["k"], tt["m"]))
            elif "Subjects" in title and "removed" in title:
                row.append(need(T[(c, P)]["need"]))
            elif "False pass" in title and "removed" in title:
                tt = T[(c, Q)]
                row.append(short(tt["k"], tt["m"]))
            else:
                row.append("—")
            name = cfg_name(c) + (" (design)" if c == DESIGN else "")
            L.append(f"| {name} | " + " | ".join(row) + " |")
    return L + [""]


def _claims(D: dict, C: dict, T: dict) -> list[str]:
    """Numbers quoted in docs/sd_stress.md prose, computed here so each has a script and an output."""
    P, Q = LEADS
    cfgs = [tuple(c) for c in D["configs"]]
    L = ["", "### Claim checks (numbers quoted in docs/sd_stress.md; SYNTHETIC)", ""]
    fp = {"with bias": [C[(c, r, Q)]["raw"] for c in cfgs for r in RISES],
          "bias removed": [C[(c, r, Q)]["unb"] for c in cfgs for r in RISES],
          "no detector": [T[(c, Q)] for c in cfgs]}
    for name, cells in fp.items():
        ups = [sim.clopper_pearson(x["k"], x["m"])[1] for x in cells]
        ps = [x["k"] / x["m"] for x in cells]
        L.append(f"- False pass at {Q} ms, {name}: {len(cells)} cells, point {min(ps):.3f}-{max(ps):.3f}, "
                 f"CP upper bound {min(ups):.4f}-{max(ups):.4f}; upper bound above 0.05 in "
                 f"{sum(u > 0.05 for u in ups)} of {len(ups)} cells.")
    gaps = [(abs(T[(c, P)]["k"] / T[(c, P)]["m"] - C[(c, r, P)]["unb"]["k"] / C[(c, r, P)]["unb"]["m"]), c, r)
            for c in cfgs for r in RISES]
    g = max(gaps)
    L.append(f"- Power at {N_SUB}, bias removed against no detector: largest gap {g[0]:.3f} "
             f"({cfg_name(g[1])}, {g[2]} ms rise).")
    for name, cells in (("with bias", [C[(c, r, P)]["raw"]["need"] for c in cfgs for r in RISES]),
                        ("bias removed", [C[(c, r, P)]["unb"]["need"] for c in cfgs for r in RISES]),
                        ("no detector", [T[(c, P)]["need"] for c in cfgs])):
        q = [d["nct"] / d["n80"] for d in cells if d["reachable"] and d.get("nct") and np.isfinite(d["n80"])]
        L.append(f"- Noncentral-t count over fitted count, {name}: {len(q)} cells, "
                 f"{(min(q) - 1) * 100:+.1f}% to {(max(q) - 1) * 100:+.1f}%.")
    dev = {}
    for key, lab in (("bias", "mean error"), ("median_err", "median error"), ("err_sd", "error SD")):
        dev[lab] = max(abs(C[(c, r, ld)][key] - C[(DESIGN, r, ld)][key]) for c in cfgs for r in RISES
                       for ld in LEADS)
    L.append("- Detector error, largest difference of any config from 8/6/10 over all rises and both true "
             "leads: " + ", ".join(f"{k} {v:.3f} ms" for k, v in dev.items()) + ".")
    return L


def cmd_report(a):
    D = json.loads(a.scored.read_text())
    C = {(tuple(D["configs"][c["cfg"]]), c["rise"], c["lead"]): c for c in D["cells"]}
    T = {(tuple(D["configs"][t["cfg"]]), t["lead"]): t for t in D["true"]}
    P, Q = LEADS
    last = D["first"] + D["n_studies"] - 1
    L = [f"ALL SYNTHETIC. {D['n_studies']} studies per cell (ids {D['first']}-{last}), "
         f"16x3x10, gate {D['gate_ms']} ms. Simulated at commit {', '.join(D['sim_head'])}, "
         f"numpy {', '.join(D['numpy'])}. Proportions are k/n = p (95% Clopper-Pearson).", ""]
    L += ["### Realised spread of the true leads (SYNTHETIC)", "",
          "| SD subject/session/trial | Subject | Session | Trial | Subject-mean SD (95% CI) | Design "
          "| Mean offset (ms) |",
          "|---|---|---|---|---|---|---|"]
    for c, v in zip(D["configs"], D["realised"]):
        L.append(f"| {cfg_name(c)} | {v['subject']:.2f} | {v['session']:.2f} | {v['trial']:.2f} | "
                 f"{v['sd_mean']:.2f} ({v['lo']:.2f}–{v['hi']:.2f}) | {v['design']:.2f} | {v['mean']:+.3f} |")
    L += ["", "### Detector error, detected - true lead (SYNTHETIC; design config, both true leads)", "",
          "| Rise | True lead | Mean (bias) | Median | SD | Missed | < -100 ms | > +100 ms |",
          "|---|---|---|---|---|---|---|---|"]
    for r in RISES:
        for lead in LEADS:
            c = C[(DESIGN, r, lead)]
            L.append(f"| {r} ms | {lead} | {c['bias']:+.2f} | {c['median_err']:+.2f} | {c['err_sd']:.2f} | "
                     f"{c['n_miss']}/{c['n_trials']} | {c['n_early']} | {c['n_late']} |")
    b = [c["bias"] for c in D["cells"]]
    L += ["", "Range over configs and true leads of the mean error: " + "; ".join(
        f"{r} ms {min(c['bias'] for c in D['cells'] if c['rise'] == r):+.2f} to "
        f"{max(c['bias'] for c in D['cells'] if c['rise'] == r):+.2f}" for r in RISES)
        + f" (all cells {min(b):+.2f} to {max(b):+.2f})."]
    L += _summary(D, C, T)
    for c in D["configs"]:
        c = tuple(c)
        t, tn = T[(c, P)], T[(c, Q)]
        L += ["", f"### SD {cfg_name(c)} ms" + (" (design)" if c == DESIGN else "") + " (SYNTHETIC)", "",
              f"| Rise | Power at {N_SUB}, with bias | Power at {N_SUB}, bias removed "
              f"| Subjects for {POWER_TARGET}, with bias | Subjects for {POWER_TARGET}, bias removed "
              f"| False pass at {Q}, with bias | False pass at {Q}, bias removed |",
              "|---|---|---|---|---|---|---|"]
        for r in RISES:
            x, y = C[(c, r, P)], C[(c, r, Q)]
            L.append(f"| {r} ms | {_fmt(x['raw']['k'], x['raw']['m'])} | "
                     f"{_fmt(x['unb']['k'], x['unb']['m'])} | "
                     f"{_need(x['raw']['need'])} | {_need(x['unb']['need'])} | "
                     f"{_fmt(y['raw']['k'], y['raw']['m'])} | {_fmt(y['unb']['k'], y['unb']['m'])} |")
        L.append(f"| no detector (true lead) | — | {_fmt(t['k'], t['m'])} | — | {_need(t['need'])} | — | "
                 f"{_fmt(tn['k'], tn['m'])} |")
    L += ["", "### Power measured at the fitted subject count, ceil(n) (SYNTHETIC)", "",
          "Disjoint studies cut from the same 16000 subjects.", "",
          "| SD | Rise | With bias | Bias removed |", "|---|---|---|---|"]
    for c in D["configs"]:
        for r in RISES:
            row = []
            for tag in ("raw", "unb"):
                k = C[(tuple(c), r, P)][tag]["need"].get("check")
                row.append(f"n={k['n']}: {_fmt(k['k'], k['m'])}" if k else "—")
            L.append(f"| {cfg_name(c)} | {r} ms | " + " | ".join(row) + " |")
    L += ["", "### Full subject sweep, with bias, passes / studies (95% CI) (SYNTHETIC)", "",
          "| SD | Subjects | Studies | " + " | ".join(f"{r} ms rise" for r in RISES) + " |",
          "|---|---|---|" + "---|" * len(RISES)]
    for c in D["configs"]:
        for gi, n in enumerate(D["n_grid"]):
            nd = [C[(tuple(c), r, P)]["raw"]["need"] for r in RISES]
            L.append(f"| {cfg_name(c)} | {n} | {nd[0]['ms'][gi]} | "
                     + " | ".join(_fmt(d["ks"][gi], d["ms"][gi]) for d in nd) + " |")
    L += ["", "### Full subject sweep, bias removed (SYNTHETIC)", "",
          "| SD | Subjects | Studies | " + " | ".join(f"{r} ms rise" for r in RISES) + " |",
          "|---|---|---|" + "---|" * len(RISES)]
    for c in D["configs"]:
        for gi, n in enumerate(D["n_grid"]):
            nd = [C[(tuple(c), r, P)]["unb"]["need"] for r in RISES]
            L.append(f"| {cfg_name(c)} | {n} | {nd[0]['ms'][gi]} | "
                     + " | ".join(_fmt(d["ks"][gi], d["ms"][gi]) for d in nd) + " |")
    L += ["", "### 8/6/10 here against docs/row12_power_sim.md "
          "(SYNTHETIC; different seeds, two-sided Fisher p)", "",
          "| Rise | Power here | Power, row 12 doc | p | False pass here | False pass, row 12 doc | p |",
          "|---|---|---|---|---|---|---|"]
    for r in RISES:
        row = [f"{r} ms"]
        for lead, lab in ((P, "power"), (Q, "false_pass")):
            c, k, n = C[(DESIGN, r, lead)]["raw"], ROW12_DOC[lab][r], ROW12_DOC["n"]
            pv = stats.fisher_exact([[c["k"], c["m"] - c["k"]], [k, n - k]])[1]
            row += [_fmt(c["k"], c["m"]), _fmt(k, n), f"{pv:.2f}"]
        L.append("| " + " | ".join(row) + " |")
    L += _claims(D, C, T)
    L += ["", "### Checks", ""]
    mm = sum(c[t]["vs_summarize"] for c in D["cells"] for t in ("raw", "unb"))
    nc = sum(c[t]["m"] for c in D["cells"] for t in ("raw", "unb"))
    flat = sum(c[t]["need"].get("boot_not_rising", 0) for c in D["cells"] if c["lead"] == P
               for t in ("raw", "unb"))
    L += [f"- Vectorised gate vs `LeadStats.summarize`: {mm} mismatches in {nc} 16-subject studies.",
          f"- Subjects with no usable trial: {sum(c['raw']['n_subj_missing'] for c in D['cells'])}.",
          "- Largest missed-trial fraction of a cell: "
          f"{max(c['n_miss'] / c['n_trials'] for c in D['cells']):.5f}.",
          f"- Bootstrap resamples whose fitted curve did not rise (counted as infinite): {flat}."]
    chi = [(c[t]["need"]["chi2"], D["configs"][c["cfg"]], c["rise"], t) for c in D["cells"] if c["lead"] == P
           for t in ("raw", "unb") if "chi2" in c[t]["need"]]
    df = len(N_GRID) - 2
    bad = [x for x in chi if stats.chi2.sf(x[0], df) < 0.01]
    names = ", ".join(f"{cfg_name(x[1])} {x[2]} ms {x[3]} ({x[0]:.0f})" for x in bad)
    L += [f"- Probit fit chi2 over {df} df (points share subjects, so indicative): p < 0.01 in {len(bad)} of "
          f"{len(chi)} fits" + (f": {names}" if bad else "") + "."]
    text = "\n".join(L) + "\n"
    a.out.write_text(text)
    print(text)


SENS_MARGINS = (5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0)  # ms; 6.0 is the protected margin


def cmd_sensitivity(a):
    """Gate pass rate when the margin between the true lead and the gate is not 6 ms.

    Analysis of the stored chunks only (nothing is simulated). The detected leads of the 74.67 ms
    trials are shifted by a constant, margin - 6 ms. That is exact for the no-detector rows; for the
    detected rows it ignores that a different true lead moves the angle samples relative to the 150 Hz
    grid, so it is an approximation. At margin 6.0 it reproduces the main tables.
    """
    _D.update(load(a.raw, a.studies, a.first))
    L = ["ALL SYNTHETIC. Pass rate of the 16 x 3 x 10 gate vs the margin (true lead - 68.67 ms); detected "
         f"leads shifted by margin - 6 ms. {a.studies} studies per cell, 95% Clopper-Pearson.", "",
         "| SD | Lead | " + " | ".join(f"{m:g} ms" for m in SENS_MARGINS) + " |",
         "|---|---|" + "---|" * len(SENS_MARGINS)]
    for ci, c in enumerate(CONFIGS):
        rows = [(f"{r} ms rise, with bias", (_D["imu"][:, ci, 0] - _D["emg"][:, ri]) * 1000.0)
                for ri, r in enumerate(RISES)]
        rows.append(("no detector", LEADS[0] + _D["off"][:, ci]))
        for name, lead in rows:
            cells = []
            for m in SENS_MARGINS:
                mean, miss = subject_stats(lead + (m - 6.0))
                p = gate_blocks(mean, miss, N_SUB)
                cells.append(sim.fmt(int(p.sum()), int(p.size)))
            L.append(f"| {cfg_name(c)} | {name} | " + " | ".join(cells) + " |")
    text = "\n".join(L) + "\n"
    a.out.write_text(text)
    print(text)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "score", "sensitivity"):
        s = sub.add_parser(name)
        s.add_argument("--raw", type=Path, required=True)
        s.add_argument("--studies", type=int, default=N_STUDIES, help="smoke tests only")
        s.add_argument("--first", type=int, default=STUDY0)
        s.add_argument("--procs", type=int, default=6)
        if name != "run":
            s.add_argument("--out", type=Path, required=True)
    s = sub.add_parser("report")
    s.add_argument("scored", type=Path)
    s.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    {"run": cmd_run, "score": cmd_score, "report": cmd_report, "sensitivity": cmd_sensitivity}[a.cmd](a)


if __name__ == "__main__":
    main()
