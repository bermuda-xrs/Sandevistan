"""Defaults and CLI behaviour fixed in the emg-onset-robust cleanup. Synthetic data only."""

import inspect
import json

import numpy as np
import pytest

from emg_onset_lead.cli import main
from emg_onset_lead.detect import BaselineRest, EmgParams, ImuParams, OnsetParams


def test_hold_default_is_30ms():
    # The 30 ms hold is the real value everywhere; the base default used to say 25 ms.
    assert OnsetParams().hold_s == 0.03
    assert EmgParams().hold_s == 0.03 and ImuParams().hold_s == 0.03


def test_at_rest_default_matches_rest_ratio():
    default = inspect.signature(BaselineRest.at_rest).parameters["max_ratio"].default
    assert default == 12.0 == OnsetParams().rest_ratio
    # An excursion of ~9 early-baseline SDs is rest under 12 (it was rejected under the old 6).
    rng = np.random.default_rng(0)
    fs, x = 100.0, rng.normal(0, 1, 100)
    early_sd = x[20:46].std(ddof=1)
    x[80] = np.median(x[20:46]) + 9.0 * early_sd
    ratio = BaselineRest.excursion_ratio(x, fs, 0.2, 1.0)
    assert 6.0 < ratio <= 12.0
    assert BaselineRest.at_rest(x, fs, 0.2, 1.0) is True
    assert BaselineRest.at_rest(x, fs, 0.2, 1.0, max_ratio=6.0) is False


def _summary(tmp_path, *extra):
    assert main(["run", "--synthetic", "--subjects", "2", "--sessions", "1", "--trials", "2",
                 "--seed", "0", "--out", str(tmp_path), *extra]) == 0
    return json.loads((tmp_path / "summary.json").read_text())


def test_cli_default_budget(tmp_path):
    lat = _summary(tmp_path)["latency"]
    assert lat["imu_path_ms"] == 18.67 and lat["total_ms"] == 68.67


def test_cli_matching_total_is_accepted(tmp_path):
    lat = _summary(tmp_path, "--total-latency-ms", "68.67")["latency"]
    assert lat["imu_path_ms"] == 18.67 and lat["total_ms"] == 68.67


def test_cli_imu_path_is_explicit(tmp_path):
    lat = _summary(tmp_path, "--imu-path-ms", "30", "--total-latency-ms", "80")["latency"]
    assert lat["imu_path_ms"] == 30.0 and lat["total_ms"] == 80.0


def test_cli_conflicting_total_is_an_error(tmp_path, capsys):
    for extra in (["--total-latency-ms", "80"],  # used to set imu_path_ms to 30 silently
                  ["--imu-path-ms", "18.67", "--total-latency-ms", "80"]):
        with pytest.raises(SystemExit) as e:
            main(["run", "--synthetic", "--subjects", "1", "--sessions", "1", "--trials", "1",
                  "--out", str(tmp_path), *extra])
        assert e.value.code == 2
        assert "--imu-path-ms" in capsys.readouterr().err
    assert not (tmp_path / "summary.json").exists()
