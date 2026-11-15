"""Defaults and CLI behaviour fixed in the emg-onset-robust cleanup. Synthetic data only."""

import inspect

import numpy as np

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
