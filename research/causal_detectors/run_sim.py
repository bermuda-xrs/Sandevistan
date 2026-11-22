"""Live EMG onset alarms: false alarms per minute of rest, and latency. ALL DATA SYNTHETIC.

  python3 research/causal_detectors/run_sim.py noise  --batches 0:720 --out RAW   (resumable; splittable)
  python3 research/causal_detectors/run_sim.py equiv  --batches 0:4   --out RAW
  python3 research/causal_detectors/run_sim.py trials --out RAW
  python3 research/causal_detectors/run_sim.py report --out RAW --results research/causal_detectors/results

noise   Rest from the repo generator: ``SyntheticSource`` with the burst pushed past the record, so every
        sample is its N(0, 1) baseline noise. One batch = 100 segments x (1 s baseline + 30 s rest) = 50
        rest-min. Stores each statistic's post-baseline running-max records, which give the first alarm for
        any threshold. Batches < FLOOR_BATCHES also run with the post-baseline noise 10% and 20% louder than
        the baseline (the README's raised-floor conditions). The unchanged production ``EmgOnsetDetector``
        (defaults, 50 ms chunks) runs on every segment whose ``env8_h30`` statistic exceeds 6, i.e. where
        its hold is met; elsewhere it cannot report an onset (``equiv`` checks that on full batches).
equiv   Production on ALL segments of some batches: are the trigger / no-trigger calls the same?
trials  ``SyntheticSource`` defaults (16 x 3 x 10, 40 ms smoothstep), seeds 0 and 1, SNR 3 / 5 / 15 (times
        the generator's U(0.7, 1.3) session gain); FC's linear-ramp ``tests/fc_synthetic.fc_trial`` (rise 10
        and 150 ms, 300 trials per cell, lead 75 ms). Production EMG and IMU detectors unchanged.
artefact  480 rest trials with one 5 ms transient of 5, 10 or 30 baseline SDs (not in the repo generator).
report  Even batches calibrate thresholds, odd batches validate. False-alarm targets:
        ``trigger``  white rest, the production trigger's own measured rate;
        ``mid``      white rest, 0.1 per hour;
        ``commit``   white rest, 0.008 per hour = the production confirmed-onset rate FC measured (4 in 499 h;
                     README, "Why 1.6"). The threshold is an exponential tail extrapolation from the top
                     100 calibration segments, so it is checked, not fitted, on the odd batches;
        ``floor10``  rest 10% louder than the baseline, the production trigger's rate in that condition;
        ``floor20c`` rest 20% louder, the production COMMIT rate in that condition.

A segment contributes rest time up to its first alarm. Rate = alarms / rest time; CIs are exact Poisson
(Garwood).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "src"), str(HERE.parents[1] / "tests")]

import detectors as D  # noqa: E402
from fc_synthetic import fc_trial  # noqa: E402

from emg_onset_lead import EmgOnsetDetector, EmgParams, ImuOnsetDetector, SyntheticSource  # noqa: E402

FS = 2000.0
N0 = round(D.BASELINE_S * FS)
SEG = 100
REST_S = 30.0
FLOORS = {"floor10": 1.10, "floor20": 1.20}
FLOOR_BATCHES = 40
SNRS = (3.0, 5.0, 15.0)
TRIAL_SEEDS = (0, 1)
FC_RISES = (0.010, 0.150)
FC_N = 300
PROD_K = EmgParams().k
MID_PER_MIN = 0.1 / 60.0
COMMIT_PER_MIN = 0.008 / 60.0  # FC, README 'Why 1.6': 4 false onsets in 499 h of white noise
ART_AMPS = (5.0, 10.0, 30.0)
N_ART = 480
KEEP = 40  # records kept per segment (the last ones: they decide every threshold a target can pick)
CHUNK = round(0.05 * FS)


def rest_batch(b: int) -> np.ndarray:
    src = SyntheticSource(n_subjects=1, n_sessions=1, n_trials=SEG, emg_snr=1.0, seed=10_000 + b,
                          duration_s=D.BASELINE_S + REST_S, emg_onset_range_s=(1e6, 1e6))
    return np.stack([t.emg for t in src])


def records(stat, keep=None):
    """Running-max records after the baseline: (row, sample, value), at most the last ``keep`` per row."""
    s = stat[:, N0:]
    cm = np.maximum.accumulate(s, axis=1)
    prev = np.concatenate([np.full((s.shape[0], 1), -np.inf), cm[:, :-1]], axis=1)
    rr, cc = np.nonzero(cm > prev)
    if keep is not None:
        # records left in the row, this one included
        n_after = np.cumsum(np.bincount(rr, minlength=s.shape[0]))[rr] - np.arange(rr.size)
        m = n_after <= keep
        rr, cc = rr[m], cc[m]
    return rr.astype(np.int32), (cc + N0).astype(np.int32), cm[rr, cc]


def first_alarm(row, idx, val, n_rows, h):
    first = np.full(n_rows, -1, dtype=np.int64)
    m = val > h
    r, pos = np.unique(row[m], return_index=True)  # records are in time order within a row
    first[r] = idx[m][pos]
    return first


def prod(x, scale=1.0):
    xi = x.copy()
    xi[N0:] *= scale
    r = EmgOnsetDetector.run(xi, FS, chunk=CHUNK)
    t = [np.nan if v is None else v for v in (r.trigger_s, r.detected_s, r.onset_s)]
    return t + [float(r.baseline_rest_ok)]


def cmd_noise(a):
    out = Path(a.out) / "noise"
    out.mkdir(parents=True, exist_ok=True)
    for b in range(*map(int, a.batches.split(":"))):
        f = out / f"b{b:04d}.npz"
        if f.exists():
            continue
        x = rest_batch(b)
        res = {"n": np.int64(x.shape[1])}
        tags = {"white": 1.0, **(FLOORS if b < FLOOR_BATCHES else {})}
        for tag, sc in tags.items():
            xx = x.copy()
            xx[:, N0:] *= sc
            for name in D.DETECTORS:
                r, c, v = records(D.statistic(name, xx, FS), KEEP)
                res[f"{name}|{tag}|row"], res[f"{name}|{tag}|idx"], res[f"{name}|{tag}|val"] = r, c, v
            hit = np.unique(res[f"env8_h30|{tag}|row"][res[f"env8_h30|{tag}|val"] > PROD_K])
            res[f"prod|{tag}|seg"] = hit
            res[f"prod|{tag}|res"] = np.array([prod(x[i], sc) for i in hit]).reshape(-1, 4)
        tmp = f.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, **res)
        tmp.rename(f)
        print("batch", b, flush=True)


def cmd_equiv(a):
    n = same = trig_off = 0
    for b in range(*map(int, a.batches.split(":"))):
        x = rest_batch(b)
        for tag, sc in {"white": 1.0, "floor20": 1.2}.items():
            xx = x.copy()
            xx[:, N0:] *= sc
            fa = D.first_alarm(D.statistic("env8_h30", xx, FS), PROD_K, N0)
            for i in range(SEG if tag == "white" else 25):
                # _pending / _trig_idx record the hold even when the trigger is later dropped
                det = EmgOnsetDetector(FS)
                first = -1
                for j in range(0, x.shape[1], CHUNK):
                    det.process(xx[i, j:j + CHUNK])
                    if first < 0 and det._trig_window is not None:
                        first = det._trig_idx
                        break
                n += 1
                same += int(first == fa[i])
                trig_off += int((first < 0) != (fa[i] < 0))
    line = (f"SYNTHETIC. production first hold sample vs env8_h30 > {PROD_K:g}: {same}/{n} rest segments "
            f"identical, {trig_off} differ in trigger / no trigger (batches {a.batches}, "
            "white 100 + floor20 25 per batch)")
    print(line)
    (Path(a.out) / "equiv.txt").write_text(line + "\n")


def trial_sets():
    for snr in SNRS:
        yield f"repo40_snr{snr:g}", [
            (t.emg, t.angle_or_imu, t.fs_emg, t.fs_imu, t.meta["t_move_s"] - t.meta["true_lead_ms"] / 1e3,
             f"{seed}-{t.session_id}")
            for seed in TRIAL_SEEDS for t in SyntheticSource(emg_snr=snr, seed=seed)]
    for rise in FC_RISES:
        for snr in SNRS:
            tr = []
            for i in range(FC_N):
                rng = np.random.default_rng([555, round(rise * 1e3), round(snr), i])
                t_on, _, emg, ang, fe, fi = fc_trial(rng, 75.0, rise, snr)
                tr.append((emg, ang, fe, fi, t_on, str(i)))
            yield f"fc{round(rise * 1e3)}_snr{snr:g}", tr


def cmd_trials(a):
    out = Path(a.out) / "trials"
    out.mkdir(parents=True, exist_ok=True)
    for cell, tr in trial_sets():
        f = out / f"{cell}.npz"
        if f.exists():
            continue
        x = np.stack([t[0] for t in tr])
        res = dict(t_on=np.array([t[4] for t in tr]), group=np.array([t[5] for t in tr]),
                   prod=np.array([prod(t[0]) for t in tr]))
        imu = []
        for t in tr:
            r = ImuOnsetDetector.run(t[1], t[3], chunk=max(1, round(0.05 * t[3])))
            imu.append([np.nan if v is None else v for v in (r.trigger_s, r.detected_s)])
        res["imu"] = np.array(imu)
        for name in D.DETECTORS:
            res[f"{name}|row"], res[f"{name}|idx"], res[f"{name}|val"] = records(D.statistic(name, x, FS))
        tmp = f.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, **res)
        tmp.rename(f)
        print("trials", cell, len(tr), flush=True)


def cmd_artefact(a):
    """Rest plus one 5 ms transient (white noise at amp x the baseline SD) at 2.0 s. NOT part of the repo
    generator: a stress test of detectors without a hold. Same rest trials for every amplitude."""
    out = Path(a.out) / "artefact"
    out.mkdir(parents=True, exist_ok=True)
    x0 = np.stack([t.emg for t in SyntheticSource(seed=777, emg_snr=1.0, duration_s=3.0,
                                                  emg_onset_range_s=(1e6, 1e6))])
    assert x0.shape[0] == N_ART
    i0, w = round(2.0 * FS), round(0.005 * FS)
    tr = np.random.default_rng(778).normal(0, 1, (N_ART, w))
    for amp in ART_AMPS:
        x = x0.copy()
        x[:, i0:i0 + w] += amp * tr
        res = dict(i0=i0, prod=np.array([prod(r) for r in x]))
        for name in D.DETECTORS:
            res[f"{name}|row"], res[f"{name}|idx"], res[f"{name}|val"] = records(D.statistic(name, x, FS))
        np.savez_compressed(out / f"amp{amp:g}.npz", **res)
        print("artefact", amp, flush=True)


# ---------------------------------------------------------------- report

def pois(k, t):
    lo = float(stats.chi2.ppf(0.025, 2 * k) / 2 / t) if k else 0.0
    return [lo, float(stats.chi2.ppf(0.975, 2 * k + 2) / 2 / t)]


def rate(alarms, minutes):
    return dict(alarms=int(alarms), minutes=round(float(minutes), 1), per_min=alarms / minutes,
                ci95=pois(alarms, minutes))


class Noise:
    """All batches of one tag, concatenated per detector."""

    def __init__(self, files, tag):
        self.tag, self.nseg, self.n, self.data, self.prod = tag, SEG * len(files), 0, {}, []
        parts = {k: ([], [], []) for k in D.DETECTORS}
        for j, f in enumerate(files):
            with np.load(f) as z:
                self.n = int(z["n"])
                for k in D.DETECTORS:
                    parts[k][0].append(z[f"{k}|{tag}|row"] + SEG * j)
                    parts[k][1].append(z[f"{k}|{tag}|idx"])
                    parts[k][2].append(z[f"{k}|{tag}|val"])
                self.prod.append(z[f"prod|{tag}|res"])
        self.data = {k: tuple(np.concatenate(q) for q in v) for k, v in parts.items()}
        self.prod = np.concatenate(self.prod)
        self.max = {k: np.maximum.reduceat(v[2], np.flatnonzero(np.diff(v[0], prepend=-1)))
                    for k, v in self.data.items()}

    def rate(self, name, h):
        fa = first_alarm(*self.data[name], self.nseg, h)
        mins = np.where(fa >= 0, fa - N0 + 1, self.n - N0).sum() / FS / 60
        return rate(int((fa >= 0).sum()), mins)

    def prod_confirmed(self):
        on = np.isfinite(self.prod[:, 2])
        full = self.nseg * (self.n - N0) / FS / 60
        mins = full - ((self.n / FS) - self.prod[on, 1]).sum() / 60
        return dict(rate(int(on.sum()), mins), triggers_checked=int(self.prod.shape[0]),
                    rejected_rest=int((self.prod[:, 3] == 0).sum()))


def threshold(noise, name, target, tail=False):
    """Threshold for ``target`` alarms/min. tail: exponential fit to the top 100 segment maxima."""
    mx = np.sort(noise.max[name])[::-1]
    if tail:
        k = target * noise.nseg * (noise.n - N0) / FS / 60
        u = mx[100]
        return float(u + (mx[:100] - u).mean() * np.log(100.0 / k))
    lo, hi = float(mx[-1]), float(mx[0]) + 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if noise.rate(name, mid)["per_min"] > target else (lo, mid)
    return hi


def ms(v):
    v = np.asarray(v, float)
    return [float(v.mean()), float(v.std(ddof=1))] if v.size > 1 else [None, None]


def summary(z, alarm):
    """Latency of ``alarm`` (s) and paired differences against production trigger and commit."""
    t_on = z["t_on"]
    hit = np.isfinite(alarm) & (alarm >= t_on)
    lat = (alarm[hit] - t_on[hit]) * 1e3
    out = dict(n=int(t_on.size), detected=int(hit.sum()), early=int((alarm < t_on).sum()),
               missed=int((~np.isfinite(alarm)).sum()), lat=ms(lat), median=float(np.median(lat)),
               p95=float(np.percentile(lat, 95)), over_300ms=int((lat > 300).sum()),
               detect_lead=ms((z["imu"][:, 1] - alarm)[hit & np.isfinite(z["imu"][:, 1])] * 1e3),
               trigger_lead=ms((z["imu"][:, 0] - alarm)[hit & np.isfinite(z["imu"][:, 0])] * 1e3))
    for lab, col in (("vs_prod_trigger", 0), ("vs_prod_commit", 1)):
        p = z["prod"][:, col]
        m = hit & np.isfinite(p) & (p >= t_on)
        d = (alarm[m] - p[m]) * 1e3
        g = np.array([d[z["group"][m] == s].mean() for s in np.unique(z["group"][m])])  # session means
        se = float(g.std(ddof=1) / np.sqrt(g.size))
        out[lab] = dict(n=int(m.sum()), mean=float(d.mean()), se_cluster=se)
    return out


def onset_err(z):
    e = (z["prod"][:, 2] - z["t_on"]) * 1e3
    return e[np.isfinite(e)]


def cmd_report(a):
    raw, res = Path(a.out), Path(a.results)
    res.mkdir(parents=True, exist_ok=True)
    files = sorted((raw / "noise").glob("b[0-9][0-9][0-9][0-9].npz"))
    ids = [int(f.stem[1:]) for f in files]
    cal = Noise([f for f, i in zip(files, ids) if i % 2 == 0], "white")
    val = Noise([f for f, i in zip(files, ids) if i % 2 == 1], "white")
    fl = {t: Noise([f for f, i in zip(files, ids) if i < FLOOR_BATCHES], t) for t in FLOORS}
    flc = {t: Noise([f for f, i in zip(files, ids) if i < FLOOR_BATCHES and i % 2 == 0], t) for t in FLOORS}
    flv = {t: Noise([f for f, i in zip(files, ids) if i < FLOOR_BATCHES and i % 2 == 1], t) for t in FLOORS}
    trig = {"cal": cal.rate("env8_h30", PROD_K), "val": val.rate("env8_h30", PROD_K),
            **{t: n.rate("env8_h30", PROD_K) for t, n in fl.items()}}
    c_cal, c_val = cal.prod_confirmed(), val.prod_confirmed()
    commit = {"all_white": rate(c_cal["alarms"] + c_val["alarms"], c_cal["minutes"] + c_val["minutes"]),
              "cal": c_cal, "val": c_val, **{t: n.prod_confirmed() for t, n in fl.items()}}
    targets = {"trigger": trig["cal"]["per_min"], "mid": MID_PER_MIN, "commit": COMMIT_PER_MIN,
               "floor10": flc["floor10"].rate("env8_h30", PROD_K)["per_min"],
               "floor20c": flc["floor20"].prod_confirmed()["per_min"]}
    on = {"floor10": (flc["floor10"], flv["floor10"]), "floor20c": (flc["floor20"], flv["floor20"])}
    art = {}
    for amp in ART_AMPS:
        with np.load(raw / "artefact" / f"amp{amp:g}.npz") as z:
            art[f"artefact x{amp:g}"] = {k: z[k] for k in z.files}
    trials = {}
    for f in sorted((raw / "trials").glob("*snr*.npz")):
        if ".tmp" not in f.name:
            with np.load(f) as z:
                trials[f.stem] = {k: z[k] for k in z.files}
    rep = dict(note="ALL DATA SYNTHETIC", batches=len(files),
               rest_hours_white=len(files) * SEG * REST_S / 3600,
               prod_trigger=trig, prod_commit=commit, prod_artefact={
                   k: int(np.isfinite(z["prod"][:, 2]).sum()) for k, z in art.items()},
               targets_per_min=targets,
               prod_latency={c: dict(trigger=summary(z, z["prod"][:, 0]), commit=summary(z, z["prod"][:, 1]),
                                     onset_err=ms(onset_err(z)))
                             for c, z in trials.items()}, cand={})
    for tn, tv in targets.items():
        rep["cand"][tn] = {}
        for name, (_, kw, desc) in D.DETECTORS.items():
            c_set, v_set = on.get(tn, (cal, val))
            fixed = name == "env8_h30" and tn in ("trigger", "floor10")  # the production threshold itself
            h = PROD_K if fixed else threshold(c_set, name, tv, tail=tn == "commit")
            lat = {}
            for c, z in trials.items():
                fa = first_alarm(z[f"{name}|row"], z[f"{name}|idx"], z[f"{name}|val"], z["t_on"].size, h)
                lat[c] = summary(z, np.where(fa >= 0, fa / FS, np.nan))
            art_n = {k: int((first_alarm(z[f"{name}|row"], z[f"{name}|idx"], z[f"{name}|val"], N_ART, h)
                             >= 0).sum()) for k, z in art.items()}
            rep["cand"][tn][name] = dict(artefact=art_n, desc=desc, kw=kw, h=h, cal=c_set.rate(name, h),
                                         val=v_set.rate(name, h),
                                         white=val.rate(name, h),
                                         **{t: n.rate(name, h) for t, n in fl.items()}, latency=lat)
    (res / "report.json").write_text(json.dumps(rep))
    tables(rep, res / "tables.md")
    digest(rep, res / "digest.md")


def f1(p, d=1):
    return "n/a" if p[0] is None else f"{p[0]:.{d}f} ({p[1]:.{d}f})"


def tables(rep, path):
    r = lambda q: f"{q['per_min']:.3g} [{q['ci95'][0]:.2g}, {q['ci95'][1]:.2g}] ({q['alarms']} in {q['minutes']:.0f} min)"  # noqa: E731,E501
    L = ["# Tables. ALL DATA SYNTHETIC. Latency in ms after the true onset: mean (SD).", "",
         f"`run_sim.py report`; {rep['batches']} batches = {rep['rest_hours_white']:.0f} h of white rest.",
         "",
         "## Production false alarms per minute of rest [95% CI]", "", 
         "| rest | trigger (hold met) | commit (reported onset) |",
         "|---|---|---|"]
    for t in ("cal", "val", "floor10", "floor20"):
        L.append(f"| {t} | {r(rep['prod_trigger'][t])} | {r(rep['prod_commit'][t])} |")
    L += [f"| all white | | {r(rep['prod_commit']['all_white'])} |", "", "## Production latency", "",
          "| cell | detected / n | early | trigger | commit | detect_lead (commit - commit) | onset error |",
          "|---|---|---|---|---|---|---|"]
    for c, p in rep["prod_latency"].items():
        k = p["commit"]
        L.append(f"| {c} | {k['detected']}/{k['n']} | {k['early']} | {f1(p['trigger']['lat'])} | "
                 f"{f1(k['lat'])} | {f1(k['detect_lead'])} | {f1(p['onset_err'], 2)} |")
    for tn, tv in rep["targets_per_min"].items():
        cand = rep["cand"][tn]
        L += ["", f"## False-alarm target `{tn}` = {tv:.3g} per min", "",
              "| detector | h | held-out rate at the target condition | white (odd batches) | +10% floor "
              "| +20% floor |",
              "|---|---|---|---|---|---|"]
        for n, q in cand.items():
            L.append(f"| {n} | {q['h']:.4g} | {r(q['val'])} | {q['white']['per_min']:.3g} | "
                     f"{q['floor10']['per_min']:.3g} | {q['floor20']['per_min']:.3g} |")
        cells = list(rep["prod_latency"])
        for group in ("repo40", "fc10", "fc150"):
            cs = [c for c in cells if c.startswith(group + "_")]
            cs.sort(key=lambda c: float(c.split("snr")[1]))
            L += ["", f"### {group}: latency mean (SD) [detected/n, early]; "
                  "paired difference vs production trigger", "",
                  "| detector | " + " | ".join(cs) + " |", "|---|" + "---|" * len(cs)]
            for n, q in cand.items():
                cols = []
                for c in cs:
                    s = q["latency"][c]
                    d = s["vs_prod_trigger"]
                    cols.append(f"{f1(s['lat'])} [{s['detected']}/{s['n']}, {s['early']}]; "
                                f"{d['mean']:+.1f} ± {d['se_cluster']:.1f}")
                L.append(f"| {n} | " + " | ".join(cols) + " |")
        cs = sorted((c for c in cells if c.startswith("repo40_")), key=lambda c: float(c.split("snr")[1]))
        L += ["", "repo40 live detect_lead_ms (IMU detected_s - alarm) mean (SD); IMU trigger_s - alarm; "
              "95th percentile latency; paired difference vs production commit (= gain in detect_lead)", "",
              "| detector | " + " | ".join(cs) + " |", "|---|" + "---|" * len(cs)]
        for n, q in cand.items():
            cols = []
            for c in cs:
                s = q["latency"][c]
                cols.append(f"{f1(s['detect_lead'])}; {f1(s['trigger_lead'])}; p95 {s['p95']:.1f}; "
                            f"{s['vs_prod_commit']['mean']:+.1f}")
            L.append(f"| {n} | " + " | ".join(cols) + " |")
    path.write_text("\n".join(L) + "\n")


def digest(rep, path):
    """One short table per target: repo40 latency, gain over the production commit, false alarms."""
    cs = [f"repo40_snr{x:g}" for x in SNRS]
    L = ["# Digest. ALL DATA SYNTHETIC (repo generator, 40 ms smoothstep; 960 trials per SNR).", "",
         "Latency = alarm - true onset, ms, mean (SD), at SNR 3 / 5 / 15. Gain = production commit - alarm,",
         "ms, paired mean = the gain in live detect_lead_ms. FA = false alarms per rest-minute, odd batches.",
         "Artefact: see cmd_artefact; share of 480 rest trials with a 5 ms transient that raise an alarm."]
    L += ["", "Production COMMITTED onsets on the same artefact trials (its trigger is env8_h30 at h = 6): "
          + "; ".join(f"{k}: {v}/{N_ART}" for k, v in rep["prod_artefact"].items())]
    for tn, tv in rep["targets_per_min"].items():
        L += ["", f"## target `{tn}` = {tv:.3g}/min", "",
              "| detector | SNR 3 | SNR 5 | SNR 15 | gain 3 / 5 / 15 | detected | early | FA at target "
              "| FA white | FA +10% | FA +20% | " + " | ".join(rep["prod_artefact"]) + " |",
              "|---|" + "---|" * (10 + len(rep["prod_artefact"]))]
        for n, q in rep["cand"][tn].items():
            lat = [q["latency"][c] for c in cs]
            gain = " / ".join(f"{-x['vs_prod_commit']['mean']:.1f}" for x in lat)
            c = q["val"]["ci95"]
            ok = "" if c[0] <= tv <= c[1] else " (target outside CI)"
            det, n_tr, early = (sum(x[k] for x in lat) for k in ("detected", "n", "early"))
            art = " | ".join(f"{q['artefact'][k]}/{N_ART}" for k in rep["prod_artefact"])
            L.append(f"| {n} | " + " | ".join(f1(x["lat"]) for x in lat)
                     + f" | {gain} | {det}/{n_tr} | {early} | {q['val']['per_min']:.2g} "
                     f"[{c[0]:.2g}, {c[1]:.2g}]{ok} | {q['white']['per_min']:.2g} | "
                     f"{q['floor10']['per_min']:.2g} | {q['floor20']['per_min']:.2g} | {art} |")
    path.write_text("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser(description="SYNTHETIC-only live onset alarm study")
    ap.add_argument("cmd", choices=["noise", "equiv", "trials", "artefact", "report"])
    ap.add_argument("--batches", default="0:720")
    ap.add_argument("--out", required=True)
    ap.add_argument("--results", default=str(HERE / "results"))
    a = ap.parse_args()
    {"noise": cmd_noise, "equiv": cmd_equiv, "trials": cmd_trials, "artefact": cmd_artefact,
     "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
