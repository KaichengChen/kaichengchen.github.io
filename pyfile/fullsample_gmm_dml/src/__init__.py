"""Full-sample debiased GMM with multiway clustered inference."""

from .core import GMMResult, OptimizerDiagnostics, fit_gmm
from .covariance import ClusterMomentResult, compute_cluster_moment_matrix
from .model_selection import MultiwayGroupKFold
from .nuisance import (
    NuisanceResult,
    NuisanceSpec,
    SklearnNuisanceFitter,
    fit_sklearn_nuisances,
    sklearn_nuisance_fit,
)

__all__ = [
    "ClusterMomentResult",
    "GMMResult",
    "MultiwayGroupKFold",
    "NuisanceResult",
    "NuisanceSpec",
    "OptimizerDiagnostics",
    "SklearnNuisanceFitter",
    "compute_cluster_moment_matrix",
    "fit_gmm",
    "fit_sklearn_nuisances",
    "sklearn_nuisance_fit",
]

__version__ = "0.1.0"
