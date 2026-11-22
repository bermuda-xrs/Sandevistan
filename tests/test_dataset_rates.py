"""Sample rates are parameters, not assumptions (synthetic data only).

The synthetic defaults stay EMG 2000 Hz / angle 150 Hz. The public datasets differ:
Lucchetti et al. 2025 records EMG at 1000 Hz and markers at 125 Hz (HS02.mat ``EmgFreq``/``KinFreq``),
and Rojas-Martinez et al. 2020 records HD-EMG at 2048 Hz. The harness must run at those rates with
the same lead accuracy, and the defaults must not move.

Measured (SYNTHETIC, 16 seeds x 72 trials, paired against 2000/150 Hz, both sources): mean lead-error
shift at most 0.14 ms at 1000/125, 2048/125, 2000/125 and 1000/150 Hz, every 95% interval containing 0;
per-trial SD 3-17% larger; one missed trial of 1152 at every rate. The common offset of about -4 ms is the
known smoothstep-burst bias (``SMOOTHSTEP_LEAD_BIAS_MS``), not a rate effect.
"""

import json
import math

import numpy as np
import pytest
from test_harness import SMOOTHSTEP_LEAD_BIAS_MS

from emg_onset_lead import (
    EmgOnsetDetector,
    EmgParams,
    Harness,
    ImuOnsetDetector,
    ImuParams,
    MinJerkSource,
    SyntheticSource,
)
from emg_onset_lead.cli import main
from emg_onset_lead.sources import DEFAULT_FS_EMG, DEFAULT_FS_IMU, check_rate

DATASET_RATES = [
    pytest.param(1000.0, 125.0, id="lucchetti-1000-125"),
    pytest.param(2048.0, 125.0, id="emg2048-imu125"),
]


def test_default_rates_unchanged():
    assert (DEFAULT_FS_EMG, DEFAULT_FS_IMU) == (2000.0, 150.0)
    src = SyntheticSource(n_subjects=1, n_sessions=1, n_trials=1)
    assert (src.fs_emg, src.fs_imu) == (2000.0, 150.0)
    t = next(iter(src))
    assert (t.fs_emg, t.fs_imu) == (2000.0, 150.0)
    assert (t.emg.size, t.angle_or_imu.size) == (8000, 600)  # 4 s at each rate


@pytest.mark.parametrize(("fs_emg", "fs_imu"), DATASET_RATES)
def test_trial_carries_rates(fs_emg, fs_imu):
    t = next(iter(MinJerkSource(n_subjects=1, n_sessions=1, n_trials=1, fs_emg=fs_emg, fs_imu=fs_imu)))
    assert (t.fs_emg, t.fs_imu) == (fs_emg, fs_imu)
    assert t.emg.size == round(4.0 * fs_emg)
    assert t.angle_or_imu.size == round(4.0 * fs_imu)


@pytest.mark.parametrize("cls", [SyntheticSource, MinJerkSource])
@pytest.mark.parametrize(("fs_emg", "fs_imu"), DATASET_RATES)
def test_lead_accuracy_at_dataset_rates(cls, fs_emg, fs_imu):
    err, missed = [], 0
    for seed in range(2):
        src = cls(n_subjects=4, n_sessions=3, n_trials=6, seed=seed, fs_emg=fs_emg, fs_imu=fs_imu)
        for r in Harness.run(src):
            if math.isfinite(r.lead_ms):
                err.append(r.lead_ms - r.true_lead_ms)
            else:
                missed += 1
    err = np.asarray(err)
    assert missed == 0
    lo, hi = SMOOTHSTEP_LEAD_BIAS_MS
    assert lo < err.mean() < hi  # same known shape bias as at 2000/150 Hz
    assert err.std(ddof=1) < 6.0  # 3.4-3.7 ms measured; 2000/150 Hz gives 3.1-3.2


@pytest.mark.parametrize("fs", [0.0, -1000.0, float("nan"), float("inf")])
def test_bad_rates_rejected(fs):
    with pytest.raises(ValueError):
        check_rate("fs", fs)
    with pytest.raises(ValueError, match="fs_emg"):
        SyntheticSource(fs_emg=fs)
    with pytest.raises(ValueError, match="fs_imu"):
        MinJerkSource(fs_imu=fs)


def test_detectors_reject_rates_their_filters_cannot_use():
    # 20 Hz band edge vs 0.45 * 40 Hz = 18 Hz cap; scipy would otherwise fail with an opaque message.
    with pytest.raises(ValueError, match="EMG band"):
        EmgOnsetDetector(40.0, EmgParams())
    with pytest.raises(ValueError, match="Nyquist"):
        ImuOnsetDetector(15.0, ImuParams())  # 10 Hz smoothing above the 7.5 Hz Nyquist
    EmgOnsetDetector(1000.0, EmgParams())  # 20-450 Hz still fits at 1000 Hz
    ImuOnsetDetector(125.0, ImuParams())


def test_cli_rates_default_and_override(tmp_path):
    assert main(["run", "--synthetic", "--subjects", "2", "--sessions", "1", "--trials", "2",
                 "--out", str(tmp_path / "d")]) == 0
    demo = json.loads((tmp_path / "d" / "summary.json").read_text())["demo"]
    assert (demo["fs_emg_hz"], demo["fs_imu_hz"]) == (2000.0, 150.0)

    assert main(["run", "--synthetic", "--subjects", "2", "--sessions", "1", "--trials", "2",
                 "--fs-emg", "1000", "--fs-imu", "125", "--out", str(tmp_path / "l")]) == 0
    demo = json.loads((tmp_path / "l" / "summary.json").read_text())["demo"]
    assert (demo["fs_emg_hz"], demo["fs_imu_hz"]) == (1000.0, 125.0)
    # The rates must not touch the gate: still 68.67 ms in the default proportional mode.
    s = json.loads((tmp_path / "l" / "summary.json").read_text())
    assert s["pooled"]["required_ms"] == pytest.approx(68.67)


@pytest.mark.parametrize("flag", ["--fs-emg", "--fs-imu"])
@pytest.mark.parametrize("value", ["0", "-125", "nan"])
def test_cli_rejects_bad_rates(tmp_path, flag, value):
    with pytest.raises(SystemExit) as exc:
        main(["run", "--synthetic", "--subjects", "1", "--sessions", "1", "--trials", "1",
              flag, value, "--out", str(tmp_path)])
    assert exc.value.code == 2
