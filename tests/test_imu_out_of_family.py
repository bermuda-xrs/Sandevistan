"""Out-of-family motion: known IMU onset gap, kept visible as strict xfails. Synthetic only.

The IMU onset fits minimum-jerk and critically damped families. Fact Checker round 2
(/workspace/fc_check/check_pr2.py, 40 trials, seed 11) found these shapes far outside the 6 ms
margin. Each test asserts the IMU onset is within the margin, so it fails today (XFAIL). If a
change fixes a shape, the strict xfail turns into a failure and the marker must come off.
"""

import numpy as np
import pytest
from fc_synthetic import MOTION_SHAPES, fc_motion_trial

from emg_onset_lead.detect import ImuOnsetDetector, ImuParams

MARGIN_MS = 6.0
N_TRIALS = 40
SEED = 11

# FC round 2 mean IMU onset error (ms) on main's IMU detector, which this branch does not change.
FC_ROUND2_MS = {
    "skewed_bell_late_peak(0.5s)": 142.8,
    "slow_creep_then_fast": 271.3,
    "cosine_ramp_0.5s": -42.3,
    "linear_ramp_0.5s": -49.3,
}


def imu_onset_errors(shape):
    rng = np.random.default_rng(SEED)
    errs = []
    for _ in range(N_TRIALS):
        _, tm, _, ang, _, fi = fc_motion_trial(rng, MOTION_SHAPES[shape])
        b = ImuOnsetDetector.run(ang, fi, ImuParams(), chunk=8)
        if b.onset_s is not None and np.isfinite(b.onset_s):
            errs.append((b.onset_s - tm) * 1e3)
    return np.array(errs)


@pytest.mark.parametrize("shape", [
    pytest.param(name, marks=pytest.mark.xfail(strict=True, raises=AssertionError, reason=(
        f"Known gap: {name} is outside the minimum-jerk / critically damped families; "
        f"FC round 2 measured a mean IMU onset error of {ms:+.1f} ms (synthetic).")))
    for name, ms in FC_ROUND2_MS.items()
])
def test_imu_onset_out_of_family(shape):
    errs = imu_onset_errors(shape)
    assert errs.size >= 0.9 * N_TRIALS
    assert abs(errs.mean()) < MARGIN_MS
