"""Row-12 power simulation with independent variance components. SYNTHETIC DATA ONLY.

No recorded EMG or motion is read. This replaces the scratch ``ps_study.py`` / ``fc4_study.py``
generator. Those seeded the subject offset with ``default_rng([seed, i])`` and the session with
``default_rng([seed, i, j])``. numpy's SeedSequence drops trailing zeros, so session 0 drew from the
subject's stream and its offset was exactly 0.75 x the subject offset. That made the SD of a subject
mean 10.57 ms instead of the 8.91 ms design value.

Seeding here is a SeedSequence tree, study -> subject -> session -> trial, built from spawn keys
(``node``). Spawn keys keep trailing zeros, so no two nodes share a stream. Each node draws only from
its own children:

    study (s,)
      subject (s, i)           child 0: subject offset ~ N(0, 8 ms)
        session (s, i, 1+j)    child 0: session offset ~ N(0, 6 ms), gain ~ U(0.7, 1.3)
          trial (..., 1+k)     child 0: trial jitter ~ N(0, 10 ms), EMG onset ~ U(1.5, 2.5) s
                               child 1: EMG noise, child 2: angle noise

Signals follow the row-12 setup (fc4_study.py): a 4 s trial with EMG at 2 kHz, white noise whose
amplitude ramps 1 -> 15*gain over a smoothstep of the given rise, and a 60 deg minimum-jerk angle
move over 0.5 s at 150 Hz (0.3 deg noise). The move starts true_lead ms after the EMG onset.

Efficiency:
  - Detectors run once per trial signal. EMG depends on (study, rise), not on the true lead; the
    angle depends on (study, true lead), not on the rise. So the same EMG detections serve both
    leads and the same angle detections serve every rise. Noise is shared across rises and leads
    (common random numbers), so cells are paired, not independent.
  - The subject sweep is nested: each study is simulated at the largest N once and scored on its
    first n subjects.
Per trial, lead and confirmation delay are combined exactly as ``Harness.run`` does
(``tests/test_row12_power_sim_seeds.py`` checks this against Harness.run). The pass/fail gate is
``LeadStats.summarize`` (68.67 ms, short_burst_peak 1.6), with planned=(n, 3) so the shape gate
follows the simulated size.

Correction addendum: rise_est = 10-90% span of the trial's own envelope / 0.6084 (smoothstep). The
envelope is zero-phase band-pass 20-450 Hz, rectified, then zero-phase 2nd-order low-pass at 20, 40
or 80 Hz (40 is primary). Baseline and plateau levels use fixed windows, never the detected onset.
Then bias = -0.109*rise_est + 0.51 ms and lead_corrected = lead - bias.

Usage (repo root):
  python scripts/row12_power_sim.py simulate --studies 0:250 --out out/row12/raw_000.npz [--procs 8]
  python scripts/row12_power_sim.py score out/row12/raw_*.npz --out out/row12/scored.json
  python scripts/row12_power_sim.py report out/row12/scored.json --out out/row12/tables.md
"""

from __future__ import annotations

import argparse
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import optimize, signal, stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from emg_onset_lead.detect import EmgOnsetDetector, EmgParams, ImuOnsetDetector, ImuParams  # noqa: E402
from emg_onset_lead.harness import Harness, TrialResult  # noqa: E402
from emg_onset_lead.latency import LatencyConfig  # noqa: E402
from emg_onset_lead.stats import LeadStats  # noqa: E402

ROOT_ENTROPY = 20261122  # one root; studies are spawn keys under it
SD_SUBJECT, SD_SESSION, SD_TRIAL = 8.0, 6.0, 10.0  # ms, design values
N_SESSIONS, N_TRIALS = 3, 10
FS_EMG, FS_IMU, T_LEN = 2000.0, 150.0, 4.0
SNR = 15.0
GATE_MS = LatencyConfig().total_ms  # 68.67, read only
TRUE_POWER, TRUE_NULL = 74.67, 68.67
N_MAX = 48
POWER_RISES = (10, 20)  # simulated at N_MAX subjects for the sweep
ADDENDUM_RISES = (10, 20, 30, 40, 80)  # scored at 16 subjects, both true leads
N_ADDENDUM = 16
N_GRID = tuple(range(16, N_MAX + 1, 4))
CUTS = (20, 40, 80)
PRIMARY_CUT = 40
BIAS_SLOPE, BIAS_ICPT = -0.109, 0.51


def _smoothstep_1090() -> float:
    f = lambda r, y: 3 * r * r - 2 * r**3 - y  # noqa: E731
    return optimize.brentq(f, 0, 1, args=(0.9,)) - optimize.brentq(f, 0, 1, args=(0.1,))


FRAC_1090 = _smoothstep_1090()  # 0.60840
_BP = signal.butter(2, [20, 450], btype="bandpass", fs=FS_EMG, output="sos")
_LP = {c: signal.butter(2, c, btype="lowpass", fs=FS_EMG, output="sos") for c in CUTS}


# ---------------------------------------------------------------- seeding / design draws
def node(*key: int) -> np.random.SeedSequence:
    """SeedSequence at spawn key ``key`` under the root (same as repeated ``.spawn`` indexing)."""
    return np.random.SeedSequence(ROOT_ENTROPY, spawn_key=tuple(int(k) for k in key))


def rng_at(*key: int) -> np.random.Generator:
    return np.random.default_rng(node(*key))


def draw_design(study: int, n_sub: int = N_MAX) -> dict:
    """Variance-component draws for one study (no signals). Arrays are (subject, session, trial)."""
    so = np.empty(n_sub)
    se = np.empty((n_sub, N_SESSIONS))
    gain = np.empty((n_sub, N_SESSIONS))
    eps = np.empty((n_sub, N_SESSIONS, N_TRIALS))
    ton = np.empty((n_sub, N_SESSIONS, N_TRIALS))
    for i in range(n_sub):
        so[i] = rng_at(study, i, 0).normal(0, SD_SUBJECT)
        for j in range(N_SESSIONS):
            r = rng_at(study, i, 1 + j, 0)
            se[i, j] = r.normal(0, SD_SESSION)
            gain[i, j] = r.uniform(0.7, 1.3)
            for k in range(N_TRIALS):
                r = rng_at(study, i, 1 + j, 1 + k, 0)
                eps[i, j, k] = r.normal(0, SD_TRIAL)
                ton[i, j, k] = r.uniform(1.5, 2.5)
    return dict(so=so, se=se, gain=gain, eps=eps, ton=ton)


def offsets(d: dict) -> np.ndarray:
    """Lead offset (ms) of every trial from the design draws: subject + session + trial."""
    return d["so"][:, None, None] + d["se"][:, :, None] + d["eps"]


def variance_check(so: np.ndarray, se: np.ndarray, eps: np.ndarray) -> dict:
    """Realised components. so (S,), se (S, n_sess), eps (S, n_sess, n_tr); S pooled over studies."""
    lead = so[:, None, None] + se[:, :, None] + eps
    subj_mean = lead.mean(axis=(1, 2))
    n_sess, n_tr = eps.shape[1], eps.shape[2]
    return dict(
        n_subjects=int(so.size),
        sd_subject=float(so.std(ddof=1)),
        sd_session=float(se.std(ddof=1)),
        sd_trial=float(eps.std(ddof=1)),
        sd_subject_mean=float(subj_mean.std(ddof=1)),
        sd_subject_mean_design=float(np.sqrt(SD_SUBJECT**2 + SD_SESSION**2 / n_sess
                                             + SD_TRIAL**2 / (n_sess * n_tr))),
        corr_subject_session=[float(np.corrcoef(so, se[:, j])[0, 1]) for j in range(n_sess)],
        corr_session_pairs=[float(np.corrcoef(se[:, a], se[:, b])[0, 1])
                            for a in range(n_sess) for b in range(a + 1, n_sess)],
    )


# ---------------------------------------------------------------- signals and detection
def _minjerk(x):
    x = np.clip(x, 0, 1)
    return x**3 * (10 - 15 * x + 6 * x * x)


def emg_signal(study, i, j, k, rise_ms, ton, gain):
    te = np.arange(int(T_LEN * FS_EMG)) / FS_EMG
    r = np.clip((te - ton) / (rise_ms / 1e3), 0, 1)
    r = r * r * (3 - 2 * r)
    return (1 + (SNR * gain - 1) * r) * rng_at(study, i, 1 + j, 1 + k, 1).normal(0, 1, te.size)


def angle_signal(study, i, j, k, ton, lead_ms):
    ti = np.arange(int(T_LEN * FS_IMU)) / FS_IMU
    tm = ton + lead_ms / 1e3
    return 20 + 60 * _minjerk((ti - tm) / 0.5) + rng_at(study, i, 1 + j, 1 + k, 2).normal(0, 0.3, ti.size)


def rise_estimate(emg: np.ndarray) -> list[float]:
    """10-90% envelope rise / FRAC_1090 (ms) per cutoff in CUTS; fixed windows, onset-independent."""
    rect = np.abs(signal.sosfiltfilt(_BP, emg))
    i0, i1, itop = int(0.3 * FS_EMG), int(1.3 * FS_EMG), emg.size - int(1.0 * FS_EMG)
    out = []
    for c in CUTS:
        e = signal.sosfiltfilt(_LP[c], rect)
        base, top = np.median(e[i0:i1]), np.median(e[itop:])
        lo, hi = base + 0.1 * (top - base), base + 0.9 * (top - base)
        above = np.flatnonzero(e >= hi)
        if top <= base or above.size == 0:
            out.append(np.nan)
            continue
        below = np.flatnonzero(e[: above[0]] <= lo)
        out.append(np.nan if below.size == 0 else (above[0] - below[-1]) / FS_EMG * 1e3 / FRAC_1090)
    return out


def detect_emg(emg):
    return EmgOnsetDetector.run(emg, FS_EMG, EmgParams(), max(1, round(0.05 * FS_EMG)))


def detect_imu(angle):
    return ImuOnsetDetector.run(angle, FS_IMU, ImuParams(), max(1, round(0.05 * FS_IMU)))


def combine(e, im) -> tuple[float, float]:
    """(lead_ms, confirmation_delay_ms) for one trial, exactly as Harness.run builds them."""
    ok = bool(e.baseline_rest_ok and im.baseline_rest_ok)
    delay = Harness._confirmation_delay_ms(e) if ok else np.nan
    if not ok or e.onset_s is None or im.onset_s is None:
        return np.nan, delay
    return (im.onset_s - e.onset_s) * 1000.0, delay


def simulate_study(study: int) -> dict:
    """Per-trial outputs for one study. Unsimulated cells are NaN."""
    d = draw_design(study, N_MAX)
    off = offsets(d)
    shape = (N_MAX, N_SESSIONS, N_TRIALS)
    rises = tuple(sorted(set(POWER_RISES) | set(ADDENDUM_RISES)))
    leads = (TRUE_POWER, TRUE_NULL)
    n_for_rise = {r: (N_MAX if r in POWER_RISES else N_ADDENDUM) for r in rises}
    n_for_lead = {TRUE_POWER: N_MAX, TRUE_NULL: N_ADDENDUM}
    rest = {r: np.full(shape + (len(CUTS),), np.nan) for r in rises}
    lead = {(r, L): np.full(shape, np.nan) for r in rises for L in leads}
    delay = {(r, L): np.full(shape, np.nan) for r in rises for L in leads}
    for i in range(N_MAX):
        for j in range(N_SESSIONS):
            for k in range(N_TRIALS):
                ton, g = d["ton"][i, j, k], d["gain"][i, j]
                imu = {L: detect_imu(angle_signal(study, i, j, k, ton, L + off[i, j, k]))
                       for L in leads if i < n_for_lead[L]}
                for r in rises:
                    if i >= n_for_rise[r]:
                        continue
                    x = emg_signal(study, i, j, k, r, ton, g)
                    rest[r][i, j, k] = rise_estimate(x)
                    e = detect_emg(x)
                    for L, im in imu.items():
                        lead[(r, L)][i, j, k], delay[(r, L)][i, j, k] = combine(e, im)
    out = dict(study=study, so=d["so"], se=d["se"], eps=d["eps"])
    for r in rises:
        out[f"rest_r{r}"] = rest[r]
        for L in leads:
            out[f"lead_r{r}_L{L}"] = lead[(r, L)]
            out[f"delay_r{r}_L{L}"] = delay[(r, L)]
    return out


def cmd_simulate(a):
    lo, hi = (int(x) for x in a.studies.split(":"))
    assert EmgParams().short_burst_peak == 1.6 and abs(GATE_MS - 68.67) < 1e-9
    with Pool(a.procs) as p:
        res = p.map(simulate_study, range(lo, hi), chunksize=1)
    keys = [k for k in res[0] if k != "study"]
    arrs = {k: np.stack([r[k] for r in res]) for k in keys}
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, studies=np.array([r["study"] for r in res]), **arrs)
    print("done", a.out, hi - lo, "studies")


# ---------------------------------------------------------------- scoring
def gate(lead: np.ndarray, delay: np.ndarray, true_ms: np.ndarray, n: int) -> tuple[bool, float]:
    """LeadStats.summarize on the first n subjects -> (pass, pooled lower bound)."""
    res = [TrialResult(f"S{i}", f"S{i}-D{j}", f"T{k}", np.nan, np.nan, float(lead[i, j, k]),
                       float(lead[i, j, k]) - GATE_MS, float(true_ms[i, j, k]), np.nan,
                       float(delay[i, j, k]))
           for i in range(n) for j in range(N_SESSIONS) for k in range(N_TRIALS)]
    s = LeadStats.summarize(res, n_boot=1, planned=(n, N_SESSIONS))
    lb = s["pooled"]["lower_bound_ms"]
    return bool(s["pass"]), float("nan") if lb is None else float(lb)


def score_study(Z: dict) -> dict:
    """Score one study; ``Z`` holds that study's arrays (as written by ``simulate``)."""
    off = Z["so"][:, None, None] + Z["se"][:, :, None] + Z["eps"]
    out = dict(study=int(Z["studies"]))
    for r in POWER_RISES:
        L = TRUE_POWER
        ld, dl = Z[f"lead_r{r}_L{L}"], Z[f"delay_r{r}_L{L}"]
        for n in N_GRID:
            p, lb = gate(ld, dl, L + off, n)
            out[f"pw_r{r}_n{n}"] = p
            out[f"lb_r{r}_n{n}"] = lb
    n = N_ADDENDUM
    for r in ADDENDUM_RISES:
        rest = Z[f"rest_r{r}"][:n]
        for L in (TRUE_POWER, TRUE_NULL):
            ld, dl, tr = Z[f"lead_r{r}_L{L}"][:n], Z[f"delay_r{r}_L{L}"][:n], (L + off)[:n]
            tag = f"r{r}_L{L}"
            out[f"raw_{tag}"] = gate(ld, dl, tr, n)[0]
            out[f"err_{tag}"] = float(np.nanmean(ld - tr))
            out[f"miss_{tag}"] = float(np.isnan(ld).mean())
            for ci, c in enumerate(CUTS):
                re = rest[..., ci]
                re = np.where(np.isfinite(re), re, np.nanmedian(re))  # study median for failed estimates
                lc = ld - (BIAS_SLOPE * re + BIAS_ICPT)
                out[f"c{c}_{tag}"] = gate(lc, dl, tr, n)[0]
                out[f"cerr{c}_{tag}"] = float(np.nanmean(lc - tr))
        for ci, c in enumerate(CUTS):
            e = rest[..., ci] - r
            out[f"rerr{c}_r{r}"] = [float(np.nansum(e)), float(np.nansum(e**2)), int(np.isfinite(e).sum()),
                                    int(np.isnan(e).sum())]
    return out


def cmd_score(a):
    rows, var = [], dict(so=[], se=[], eps=[])
    for f in a.raw:
        Z = dict(np.load(f))
        for k in var:
            var[k].append(Z[k].reshape((-1,) + Z[k].shape[2:]))
        with Pool(a.procs) as p:
            per = [{k: v[i] for k, v in Z.items()} for i in range(len(Z["studies"]))]
            rows += p.map(score_study, per, chunksize=2)
    seen = [r["study"] for r in rows]
    assert len(seen) == len(set(seen)), "duplicate study ids across raw files"
    vc = variance_check(*(np.concatenate(var[k]) for k in ("so", "se", "eps")))
    json.dump(dict(variance=vc, studies=sorted(rows, key=lambda r: r["study"])), open(a.out, "w"))
    print("done", a.out, len(rows), "studies")


# ---------------------------------------------------------------- report
def clopper_pearson(k: int, n: int) -> tuple[float, float]:
    lo = 0.0 if k == 0 else stats.beta.ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(0.975, k + 1, n - k)
    return float(lo), float(hi)


def fmt(k: int, n: int) -> str:
    lo, hi = clopper_pearson(k, n)
    return f"{k}/{n} = {k / n:.3f} ({lo:.3f}–{hi:.3f})"


def probit_sqrt_fit(ns, ks, n_studies) -> tuple[float, np.ndarray]:
    """Binomial ML fit P(pass) = Phi(b0 + b1*sqrt(n)); returns (n at P = 0.8, b)."""
    x = np.sqrt(np.asarray(ns, float))
    ks = np.asarray(ks, float)

    def nll(b):
        q = np.clip(stats.norm.cdf(b[0] + b[1] * x), 1e-12, 1 - 1e-12)
        return -np.sum(ks * np.log(q) + (n_studies - ks) * np.log(1 - q))

    b = optimize.minimize(nll, [-3.0, 1.0], method="Nelder-Mead",
                          options=dict(xatol=1e-9, fatol=1e-11, maxiter=8000)).x
    return float(((stats.norm.ppf(0.8) - b[0]) / b[1]) ** 2), b


def cmd_report(a):
    D = json.load(open(a.scored))
    S, vc = D["studies"], D["variance"]
    N = len(S)
    rng = np.random.default_rng(2026)
    L = [f"Studies scored: {N}. ALL SYNTHETIC.", "", "### Variance check (all simulated subjects)", "",
         "| Component | Realised | Design |", "|---|---|---|",
         f"| Subject SD (ms) | {vc['sd_subject']:.3f} | 8 |",
         f"| Session SD (ms) | {vc['sd_session']:.3f} | 6 |",
         f"| Trial SD (ms) | {vc['sd_trial']:.3f} | 10 |",
         f"| Subject-mean true-lead SD, 3x10 (ms) | {vc['sd_subject_mean']:.3f} | "
         f"{vc['sd_subject_mean_design']:.3f} |",
         f"| corr(subject offset, session j offset), j = 0, 1, 2 | "
         f"{', '.join(f'{c:+.4f}' for c in vc['corr_subject_session'])} | 0 |",
         f"| corr between session offsets (0-1, 0-2, 1-2) | "
         f"{', '.join(f'{c:+.4f}' for c in vc['corr_session_pairs'])} | 0 |",
         f"| Subjects pooled | {vc['n_subjects']} | |", ""]
    L += ["### Pass rate vs subjects (x 3 sessions x 10 trials), true lead 74.67 ms, gate 68.67 ms", "",
          "| Subjects | Trials/study | " + " | ".join(f"{r} ms rise" for r in POWER_RISES) + " |",
          "|---|---|" + "---|" * len(POWER_RISES)]
    M = {r: np.array([[s[f"pw_r{r}_n{n}"] for n in N_GRID] for s in S], float) for r in POWER_RISES}
    for ni, n in enumerate(N_GRID):
        cells = " | ".join(fmt(int(M[r][:, ni].sum()), N) for r in POWER_RISES)
        L.append(f"| {n} | {n * 30} | {cells} |")
    L += ["", "### Subjects for 0.8 power: probit-in-sqrt(subjects) binomial ML fit over all grid points; "
          "95% CI = 2000x bootstrap over studies (each resampled study keeps all its nested sizes)", "",
          "| Rise | n at 0.8 (point) | 95% bootstrap CI | Trials/study at the point | Fit chi2 p |",
          "|---|---|---|---|---|"]
    cross = {}
    for r in POWER_RISES:
        ks = M[r].sum(0)
        x, b = probit_sqrt_fit(N_GRID, ks, N)
        bs = [probit_sqrt_fit(N_GRID, M[r][rng.integers(0, N, N)].sum(0), N)[0] for _ in range(2000)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        pred = stats.norm.cdf(b[0] + b[1] * np.sqrt(N_GRID))
        chi = float(np.sum((ks - N * pred) ** 2 / np.maximum(N * pred * (1 - pred), 1e-9)))
        pchi = stats.chi2.sf(chi, len(N_GRID) - 2)
        cross[r] = (x, lo, hi)
        L.append(f"| {r} ms | {x:.1f} | {lo:.1f}–{hi:.1f} | {x * 30:.0f} | {pchi:.2f} (chi2 {chi:.1f}, "
                 f"{len(N_GRID) - 2} df; nested points are correlated, so this is indicative) |")

    def rerr(r, c):
        t = np.array([s[f"rerr{c}_r{r}"] for s in S], float).sum(0)
        m = t[0] / t[2]
        return m, np.sqrt((t[1] - t[2] * m * m) / (t[2] - 1)), int(t[3])

    L += ["", f"### Baseline five-row table (no correction), {N_ADDENDUM}x3x10, fixed seeds", "",
          f"| Rise | Power at {TRUE_POWER} ms | False pass at {TRUE_NULL} ms | Mean lead err (ms) |",
          "|---|---|---|---|"]
    for r in ADDENDUM_RISES:
        kp = sum(s[f"raw_r{r}_L{TRUE_POWER}"] for s in S)
        kf = sum(s[f"raw_r{r}_L{TRUE_NULL}"] for s in S)
        L.append(f"| {r} ms | {fmt(kp, N)} | {fmt(kf, N)} | "
                 f"{np.mean([s[f'err_r{r}_L{TRUE_POWER}'] for s in S]):+.2f} |")
    for title, Lv, lab in ((f"Power, true lead {TRUE_POWER} ms (row 12)", TRUE_POWER, "Power"),
                           (f"False pass, true lead exactly {TRUE_NULL} ms (nominal 0.05)", TRUE_NULL,
                            "False pass")):
        L += ["", f"### {title}, {N_ADDENDUM}x3x10, correction primary = {PRIMARY_CUT} Hz envelope", "",
              f"| Rise | {lab} uncorrected | Mean lead err (ms) | {lab} corrected "
              f"| Mean corrected lead err (ms) | Rise-est err mean ± SD (ms) |"
              + (" One-sided binomial p (>0.05), uncorr / corr |" if Lv == TRUE_NULL else ""),
              "|---|---|---|---|---|---|" + ("---|" if Lv == TRUE_NULL else "")]
        for r in ADDENDUM_RISES:
            t = f"r{r}_L{Lv}"
            k0 = sum(s[f"raw_{t}"] for s in S)
            k1 = sum(s[f"c{PRIMARY_CUT}_{t}"] for s in S)
            m, sd, nf = rerr(r, PRIMARY_CUT)
            row = (f"| {r} ms | {fmt(k0, N)} | {np.mean([s[f'err_{t}'] for s in S]):+.2f} | {fmt(k1, N)} | "
                   f"{np.mean([s[f'cerr{PRIMARY_CUT}_{t}'] for s in S]):+.2f} | {m:+.1f} ± {sd:.1f} |")
            if Lv == TRUE_NULL:
                p0 = stats.binomtest(k0, N, 0.05, alternative="greater").pvalue
                p1 = stats.binomtest(k1, N, 0.05, alternative="greater").pvalue
                row += f" {p0:.2g} / {p1:.2g} |"
            L.append(row)
    L += ["", f"### Cutoff sensitivity ({N_ADDENDUM}x3x10, same studies)", "",
          "| Rise | Cutoff | Power corrected | False pass corrected | Rise-est err mean ± SD (ms) |",
          "|---|---|---|---|---|"]
    for r in ADDENDUM_RISES:
        for c in CUTS:
            kp = sum(s[f"c{c}_r{r}_L{TRUE_POWER}"] for s in S)
            kf = sum(s[f"c{c}_r{r}_L{TRUE_NULL}"] for s in S)
            m, sd, nf = rerr(r, c)
            L.append(f"| {r} ms | {c} Hz | {fmt(kp, N)} | {fmt(kf, N)} | {m:+.1f} ± {sd:.1f} (failed {nf}) |")
    miss = max(s[k] for s in S for k in s if k.startswith("miss_"))
    L += ["", f"Max per-study missed-trial fraction (addendum cells): {miss:.4f}."]
    Path(a.out).write_text("\n".join(L) + "\n")
    json.dump({str(r): v for r, v in cross.items()}, open(Path(a.out).with_suffix(".cross.json"), "w"))
    print("\n".join(L))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("simulate")
    s.add_argument("--studies", required=True, help="lo:hi study ids")
    s.add_argument("--out", required=True)
    s.add_argument("--procs", type=int, default=8)
    s = sub.add_parser("score")
    s.add_argument("raw", nargs="+")
    s.add_argument("--out", required=True)
    s.add_argument("--procs", type=int, default=8)
    s = sub.add_parser("report")
    s.add_argument("scored")
    s.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    {"simulate": cmd_simulate, "score": cmd_score, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
