"""Causal EMG-vs-IMU onset-lead analysis harness."""

from .detect import EmgOnsetDetector, EmgParams, ImuOnsetDetector, ImuParams, OnsetResult
from .harness import Harness, TrialResult
from .latency import LatencyConfig
from .sources import MinJerkSource, SyntheticSource, Trial, TrialSource
from .stats import LeadStats, SampleSize

__all__ = [
    "EmgOnsetDetector", "EmgParams", "ImuOnsetDetector", "ImuParams", "OnsetResult",
    "Harness", "TrialResult", "LatencyConfig", "SyntheticSource", "MinJerkSource",
    "Trial", "TrialSource", "LeadStats", "SampleSize",
]
