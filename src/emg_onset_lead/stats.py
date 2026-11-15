from __future__ import annotations

import numpy as np
from scipy import stats as sps

from .latency import LatencyConfig

# Planned onset pilot: 16 subjects x 3 sessions per subject (decided 2026-11-03; 8x3 was
# rejected at about 53% power, computed from synthetic variances, subject/session/trial SD
# 8/6/10 ms, not pilot data). The adopted power calculation assumed 3 sessions per subject.
PLANNED_SUBJECTS = 16
PLANNED_SESSIONS = 3
PLANNED_SHAPE = (PLANNED_SUBJECTS, PLANNED_SESSIONS)


class SampleSize:
    """One-sided t power for a subject-mean lead.

    ``between_var`` is the variance of the subject effect. ``within_var`` is the variance
    of one within-subject replicate (a session mean, if that is what the pilot reports).
    The variance of a subject mean is ``between_var + within_var / n_per_subject``.
    The test is the same one the harness gates on: subject means versus a margin above
    the latency budget. Inputs are pilot variances; nothing here is estimated from the
    synthetic generator, and the planned 16x3 is unverified until those come from real data.
    """

    @staticmethod
    def subject_mean_variance(between_var: float, within_var: float, n_per_subject: int) -> float:
        if n_per_subject < 1:
            raise ValueError("n_per_subject must be >= 1")
        return float(between_var + within_var / n_per_subject)

    @staticmethod
    def design_subject_sd(subject_sd: float, session_sd: float, trial_sd: float,
                          n_sessions: int, n_trials: int) -> float:
        """SD of a subject mean under the synthetic variance components.

        With the generator defaults (8, 6, 10 ms) and 3 sessions x 10 trials this is 8.91 ms.
        A single seed can land far from it.
        """
        if n_sessions < 1 or n_trials < 1:
            raise ValueError("n_sessions and n_trials must be >= 1")
        var = subject_sd**2 + session_sd**2 / n_sessions + trial_sd**2 / (n_sessions * n_trials)
        return float(np.sqrt(var))

    @staticmethod
    def power(between_var: float, within_var: float, effect_margin: float, alpha: float,
              n_subjects: int, n_per_subject: int) -> float:
        """Power of a one-sided one-sample t-test on ``n_subjects`` subject means."""
        if n_subjects < 2:
            raise ValueError("n_subjects must be >= 2")
        if not 0 < alpha < 1:
            raise ValueError("alpha must be in (0, 1)")
        sigma = np.sqrt(SampleSize.subject_mean_variance(between_var, within_var, n_per_subject))
        if sigma == 0:
            return 1.0 if effect_margin > 0 else float(alpha)
        df = n_subjects - 1
        ncp = effect_margin / (sigma / np.sqrt(n_subjects))
        tcrit = float(sps.t.ppf(1 - alpha, df))
        return float(sps.nct.sf(tcrit, df, ncp))

    @staticmethod
    def required_subjects(between_var: float, within_var: float, effect_margin: float, alpha: float,
                          power: float, n_per_subject: int, max_n: int = 500) -> int | None:
        """Smallest subject count with power at least ``power``, or None if none up to ``max_n``."""
        if not 0 < power < 1:
            raise ValueError("power must be in (0, 1)")
        for n in range(2, max_n + 1):
            if SampleSize.power(between_var, within_var, effect_margin, alpha, n, n_per_subject) >= power:
                return n
        return None


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
                  planned: tuple[int, int] = PLANNED_SHAPE, max_missed_fraction: float = 0.10) -> dict:
        delays = [float(r.emg_confirmation_delay_ms) for r in results
                  if np.isfinite(getattr(r, "emg_confirmation_delay_ms", np.nan))]
        delay_mean = float(np.mean(delays)) if delays else 0.0
        req = latency.required_ms(path, delay_mean)
        nested: dict[str, dict[str, np.ndarray]] = {}
        missed: dict[tuple[str, str], int] = {}
        trues: dict[tuple[str, str], list[float]] = {}
        # Study shape comes from the input trials (detected or not), so a fully missed
        # subject or session still counts against the planned shape.
        input_sessions: dict[str, set[str]] = {}
        for r in results:
            key = (r.subject_id, r.session_id)
            input_sessions.setdefault(r.subject_id, set()).add(r.session_id)
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
        n_missed = int(sum(missed.values()))
        n_total = len(results)
        miss_fraction = (n_missed / n_total) if n_total else 1.0
        pooled = {"n_subjects": len(nested), "n_sessions": len(sessions),
                  "n_trials": int(sum(s["n_trials"] for s in sessions)),
                  "n_missed": n_missed,
                  "n_total": n_total,
                  "miss_fraction": miss_fraction,
                  "method": "t_on_subject_means",
                  "bootstrap_lower_bound_ms": LeadStats.boot_lower_bound(nested, conf, n_boot, seed),
                  **row(m, lb)}
        all_true = [x for v in trues.values() for x in v]
        if all_true:
            pooled["true_lead_ms_mean"] = float(np.mean(all_true))
        sess_counts = [len(v) for v in input_sessions.values()]
        shape_ok = (len(input_sessions) == planned[0] and bool(sess_counts)
                    and all(c == planned[1] for c in sess_counts))
        # 134/240 missed is 0.56 and must fail. A lower-bound pass on the trials that
        # happened to detect is not enough, and neither is a study that is not the planned shape.
        miss_ok = bool(n_total) and miss_fraction <= max_missed_fraction
        bound_ok = bool(pooled["pass"])
        passed = bool(bound_ok and shape_ok and miss_ok)
        if delays:
            delay_report = {"mean": delay_mean, "min": float(np.min(delays)),
                            "max": float(np.max(delays)), "n": len(delays)}
        else:
            delay_report = {"mean": None, "min": None, "max": None, "n": 0}
        return {"confidence": conf, "criterion_path": path,
                "latency": {"imu_path_ms": latency.imu_path_ms, "emg_extra_ms": latency.emg_extra_ms,
                            "total_ms": latency.total_ms, "triggered": bool(latency.triggered),
                            "emg_extra_is": "t_win/2",
                            "confirmation_delay_ms": delay_mean if latency.triggered else 0.0},
                "confirmation_delay_ms": delay_report,
                "study": {"planned": {"subjects": planned[0], "sessions_per_subject": planned[1],
                                      "min_days": 3, "electrodes_redonned_between_sessions": True},
                          "actual": {"subjects": len(input_sessions),
                                     "sessions_per_subject_max": max(sess_counts, default=0)},
                          "matches_planned_shape": shape_ok,
                          "plan_status": "unverified",
                          "plan_note": (
                              f"{planned[0]} subjects x {planned[1]} sessions is unverified until real data "
                              "supply between-subject and within-subject variance. Use SampleSize with "
                              "those pilot variances; synthetic power is not a substitute."
                          )},
                "gates": {"lower_bound": bound_ok, "planned_shape": shape_ok,
                          "missed_trials": miss_ok, "max_missed_fraction": max_missed_fraction,
                          "miss_fraction": miss_fraction},
                "sessions": sessions, "subjects": subjects, "pooled": pooled,
                "pass": passed}
