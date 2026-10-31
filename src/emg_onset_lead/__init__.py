"""Causal EMG-vs-IMU onset-lead analysis harness."""

from .detect import EmgOnsetDetector, EmgParams, ImuOnsetDetector, ImuParams, OnsetResult
from .harness import Harness, TrialResult
from .latency import LatencyConfig
from .sources import SyntheticSource, Trial, TrialSource
from .stats import LeadStats

__all__ = [
    "EmgOnsetDetector", "EmgParams", "ImuOnsetDetector", "ImuParams", "OnsetResult",
    "Harness", "TrialResult", "LatencyConfig", "SyntheticSource", "Trial", "TrialSource",
    "LeadStats",
]
