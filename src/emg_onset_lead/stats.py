from __future__ import annotations

import numpy as np
from scipy import stats as sps

from .latency import LatencyConfig


class LeadStats:
    """Lower confidence bounds respecting subject > session > trial nesting."""

    @staticmethod
    def t_lower_bound(values, conf: float = 0.95) -> tuple[float, float | None]:
        """(mean, one-sided lower bound). Bound is None when n < 2."""
        v = np.asarray(values, dtype=float)
        if v.size == 0:
            return float("nan"), None
        m = float(v.mean())
        if v.size < 2:
            return m, None
        se = float(v.std(ddof=1)) / np.sqrt(v.size)
        return m, m - float(sps.t.ppf(conf, v.size - 1)) * se

    @staticmethod
    def boot_lower_bound(nested: dict[str, dict[str, np.ndarray]], conf: float = 0.95,
                         n_boot: int = 2000, seed: int = 0) -> float | None:
        """Cluster bootstrap: resample subjects, then sessions within subject.

        Statistic = mean over subjects of (mean over sessions of session-mean lead).
        """
        subj = [[float(np.mean(v)) for v in sess.values()] for sess in nested.values()]
        if len(subj) < 2:
            return None
        rng = np.random.default_rng(seed)
        stat = np.empty(n_boot)
        for b in range(n_boot):
            pick = rng.integers(0, len(subj), len(subj))
            means = []
            for i in pick:
                sm = subj[i]
                means.append(np.mean(rng.choice(sm, size=len(sm), replace=True)))
            stat[b] = np.mean(means)
        return float(np.quantile(stat, 1 - conf))

    @staticmethod
    def summarize(results, latency: LatencyConfig = LatencyConfig(), conf: float = 0.95,
                  path: str = "total", n_boot: int = 2000, seed: int = 0,
                  planned: tuple[int, int] = (8, 3)) -> dict:
        req = latency.required_ms(path)
        nested: dict[str, dict[str, np.ndarray]] = {}
        missed: dict[tuple[str, str], int] = {}
        trues: dict[tuple[str, str], list[float]] = {}
        for r in results:
            key = (r.subject_id, r.session_id)
            nested.setdefault(r.subject_id, {}).setdefault(r.session_id, [])
            if np.isfinite(r.lead_ms):
                nested[r.subject_id][r.session_id].append(r.lead_ms)
            else:
                missed[key] = missed.get(key, 0) + 1
            if np.isfinite(r.true_lead_ms):
                trues.setdefault(key, []).append(r.true_lead_ms)
        nested = {s: {k: np.asarray(v) for k, v in ss.items() if len(v)} for s, ss in nested.items()}
        nested = {s: ss for s, ss in nested.items() if ss}

        def row(mean, lb):
            return {"mean_lead_ms": mean, "lower_bound_ms": lb, "required_ms": req,
                    "margin_ms": None if lb is None else lb - req,
                    "pass": bool(lb is not None and lb > req)}

        sessions = []
        for s, ss in nested.items():
            for k, v in ss.items():
                m, lb = LeadStats.t_lower_bound(v, conf)
                sessions.append({"subject_id": s, "session_id": k, "n_trials": int(v.size),
                                 "n_missed": missed.get((s, k), 0),
                                 "sd_lead_ms": float(v.std(ddof=1)) if v.size > 1 else None,
                                 **row(m, lb)})
        subjects, sub_means = [], []
        for s, ss in nested.items():
            smeans = np.array([v.mean() for v in ss.values()])
            vals = smeans if smeans.size > 1 else next(iter(ss.values()))
            _, lb = LeadStats.t_lower_bound(vals, conf)
            m = float(smeans.mean())
            sub_means.append(m)
            subjects.append({"subject_id": s, "n_sessions": len(ss), **row(m, lb),
                             "bound_basis": "session_means" if smeans.size > 1 else "trials"})
        m, lb = LeadStats.t_lower_bound(sub_means, conf)
        pooled = {"n_subjects": len(nested), "n_sessions": len(sessions),
                  "n_trials": int(sum(s["n_trials"] for s in sessions)),
                  "n_missed": int(sum(missed.values())),
                  "method": "t_on_subject_means",
                  "bootstrap_lower_bound_ms": LeadStats.boot_lower_bound(nested, conf, n_boot, seed),
                  **row(m, lb)}
        all_true = [x for v in trues.values() for x in v]
        if all_true:
            pooled["true_lead_ms_mean"] = float(np.mean(all_true))
        shape = (len(nested), max((len(v) for v in nested.values()), default=0))
        return {"confidence": conf, "criterion_path": path,
                "latency": {"imu_path_ms": latency.imu_path_ms, "emg_extra_ms": latency.emg_extra_ms,
                            "total_ms": latency.total_ms},
                "study": {"planned": {"subjects": planned[0], "sessions_per_subject": planned[1],
                                      "min_days": 3, "electrodes_redonned_between_sessions": True},
                          "actual": {"subjects": shape[0], "sessions_per_subject_max": shape[1]},
                          "matches_planned_shape": shape == planned},
                "sessions": sessions, "subjects": subjects, "pooled": pooled,
                "pass": pooled["pass"]}
