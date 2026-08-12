"""Generic full-sample one- and two-step GMM estimation."""

from __future__ import annotations

import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import pandas as pd
from scipy.optimize import OptimizeResult, minimize
from scipy.stats import norm

from .covariance import VarianceType, compute_cluster_moment_matrix
from .nuisance import SklearnNuisanceFitter


ScoreFunction = Callable[[pd.DataFrame, np.ndarray, Any], np.ndarray]
JacobianFunction = Callable[[pd.DataFrame, np.ndarray, Any], np.ndarray]


@dataclass(frozen=True)
class OptimizerDiagnostics:
    """Compact diagnostics from a GMM optimization step."""

    success: bool
    status: int
    message: str
    objective: float
    n_iterations: int | None
    n_function_evaluations: int | None


@dataclass(frozen=True)
class GMMResult:
    """Results from full-sample clustered GMM."""

    theta: np.ndarray
    se: np.ndarray
    ci: np.ndarray
    covariance: np.ndarray
    moment: np.ndarray
    jacobian: np.ndarray
    middle_matrix: np.ndarray
    preliminary_middle_matrix: np.ndarray
    weight_matrix: np.ndarray
    first_step_theta: np.ndarray
    n_eff: int
    cluster_counts: dict[str, int]
    n_cells: int
    n_observations: int
    nuisance: Any
    n_steps: Literal[1, 2]
    weight_type: VarianceType
    covariance_type: VarianceType
    first_step_optimizer: OptimizerDiagnostics
    second_step_optimizer: OptimizerDiagnostics | None

    @property
    def coef(self) -> np.ndarray:
        """Alias for the parameter estimate."""

        return self.theta

    @property
    def standard_errors(self) -> np.ndarray:
        """Alias for ``se``."""

        return self.se

    @property
    def confidence_interval(self) -> np.ndarray:
        """Alias for the two-column 95% confidence interval."""

        return self.ci


def _validate_variance_type(value: str, argument: str) -> VarianceType:
    if value not in {"PSD", "CGM"}:
        raise ValueError(f"{argument} must be exactly 'PSD' or 'CGM'")
    return value  # type: ignore[return-value]


def _evaluate_score(
    score: ScoreFunction,
    data: pd.DataFrame,
    theta: np.ndarray,
    nuisance: Any,
    expected_q: int | None = None,
) -> np.ndarray:
    values = np.asarray(score(data, theta, nuisance), dtype=float)
    if values.ndim != 2:
        raise ValueError("score must return a two-dimensional n-by-q array")
    if values.shape[0] != len(data):
        raise ValueError("score must return one row per observation in data")
    if values.shape[1] == 0:
        raise ValueError("score must return at least one moment")
    if expected_q is not None and values.shape[1] != expected_q:
        raise ValueError("score returned a different number of moments across calls")
    if not np.isfinite(values).all():
        raise ValueError("score returned non-finite values")
    return values


def _diagnostics(result: OptimizeResult) -> OptimizerDiagnostics:
    return OptimizerDiagnostics(
        success=bool(result.success),
        status=int(result.status),
        message=str(result.message),
        objective=float(result.fun),
        n_iterations=(None if result.get("nit") is None else int(result.nit)),
        n_function_evaluations=(
            None if result.get("nfev") is None else int(result.nfev)
        ),
    )


def _optimize(
    score: ScoreFunction,
    data: pd.DataFrame,
    nuisance: Any,
    start: np.ndarray,
    weight: np.ndarray,
    q: int,
    method: str,
    options: Mapping[str, Any] | None,
) -> tuple[np.ndarray, OptimizerDiagnostics]:
    def objective(theta: np.ndarray) -> float:
        moment = _evaluate_score(
            score, data, np.asarray(theta, dtype=float), nuisance, q
        ).mean(axis=0)
        return float(moment @ weight @ moment)

    try:
        result = minimize(
            objective,
            start,
            method=method,
            options=None if options is None else dict(options),
        )
    except Exception as exc:
        raise RuntimeError(f"GMM optimization could not be run: {exc}") from exc
    diagnostics = _diagnostics(result)
    if not diagnostics.success:
        raise RuntimeError(f"GMM optimization failed: {diagnostics.message}")
    theta = np.asarray(result.x, dtype=float)
    if theta.shape != start.shape or not np.isfinite(theta).all():
        raise RuntimeError("GMM optimization returned an invalid parameter vector")
    return theta, diagnostics


def _numerical_jacobian(
    score: ScoreFunction,
    data: pd.DataFrame,
    theta: np.ndarray,
    nuisance: Any,
    q: int,
) -> np.ndarray:
    result = np.empty((q, len(theta)), dtype=float)
    for index in range(len(theta)):
        step = 1e-6 * max(1.0, abs(theta[index]))
        upper = theta.copy()
        lower = theta.copy()
        upper[index] += step
        lower[index] -= step
        upper_moment = _evaluate_score(
            score, data, upper, nuisance, q
        ).mean(axis=0)
        lower_moment = _evaluate_score(
            score, data, lower, nuisance, q
        ).mean(axis=0)
        result[:, index] = -(upper_moment - lower_moment) / (2.0 * step)
    return result


def _jacobian(
    jacobian: JacobianFunction | None,
    score: ScoreFunction,
    data: pd.DataFrame,
    theta: np.ndarray,
    nuisance: Any,
    q: int,
    p: int,
) -> np.ndarray:
    if jacobian is None:
        result = _numerical_jacobian(score, data, theta, nuisance, q)
    else:
        result = np.asarray(jacobian(data, theta, nuisance), dtype=float)
    if result.shape != (q, p):
        raise ValueError(f"jacobian must return an array with shape {(q, p)}")
    if not np.isfinite(result).all():
        raise ValueError("jacobian returned non-finite values")
    return result


def _positive_definite_inverse(
    matrix: np.ndarray,
    method: VarianceType,
) -> np.ndarray:
    symmetric = (matrix + matrix.T) / 2.0
    eigenvalues = np.linalg.eigvalsh(symmetric)
    tolerance = (
        np.finfo(float).eps
        * max(1, len(symmetric))
        * max(1.0, float(np.max(np.abs(eigenvalues))))
        * 100.0
    )
    if float(eigenvalues.min()) <= tolerance:
        if method == "CGM":
            raise RuntimeError(
                "the preliminary CGM moment matrix is not positive definite; "
                "use weight_type='PSD'"
            )
        raise RuntimeError(
            "the preliminary PSD moment matrix is not positive definite; "
            "a two-step weight cannot be constructed"
        )
    return np.linalg.inv(symmetric)


def _coefficient_covariance(
    jacobian: np.ndarray,
    weight: np.ndarray,
    middle: np.ndarray,
    n_eff: int,
) -> np.ndarray:
    gram = jacobian.T @ weight @ jacobian
    try:
        bread = np.linalg.inv(gram)
    except np.linalg.LinAlgError as exc:
        raise RuntimeError(
            "the weighted Jacobian is singular; the parameter is not identified"
        ) from exc
    projection = bread @ jacobian.T @ weight
    covariance = projection @ middle @ projection.T / n_eff
    covariance = (covariance + covariance.T) / 2.0
    if not np.isfinite(covariance).all():
        raise RuntimeError("the coefficient covariance contains non-finite values")
    return covariance


def _standard_errors_and_ci(
    theta: np.ndarray,
    covariance: np.ndarray,
    covariance_type: VarianceType,
) -> tuple[np.ndarray, np.ndarray]:
    variances = np.diag(covariance).copy()
    scale = max(1.0, float(np.max(np.abs(variances))))
    tolerance = np.finfo(float).eps * len(theta) * scale * 100.0
    materially_negative = variances < -tolerance
    if materially_negative.any():
        warnings.warn(
            f"{covariance_type} produced a negative coefficient variance; "
            "the corresponding standard error and confidence interval are undefined",
            RuntimeWarning,
            stacklevel=3,
        )
    variances[(variances < 0.0) & ~materially_negative] = 0.0
    se = np.full_like(variances, np.nan)
    valid = variances >= 0.0
    se[valid] = np.sqrt(variances[valid])
    critical_value = float(norm.ppf(0.975))
    ci = np.column_stack((theta - critical_value * se, theta + critical_value * se))
    return se, ci


def fit_gmm(
    data: pd.DataFrame,
    score: ScoreFunction,
    theta_start: np.ndarray | Sequence[float],
    cluster_cols: Sequence[str],
    nuisance: Any = None,
    nuisance_fit: Callable[[pd.DataFrame], Any] | None = None,
    jacobian: JacobianFunction | None = None,
    weight_type: VarianceType = "PSD",
    covariance_type: VarianceType = "PSD",
    optimizer_method: str = "BFGS",
    optimizer_options: Mapping[str, Any] | None = None,
    n_steps: Literal[1, 2] = 2,
) -> GMMResult:
    """Estimate generic full-sample GMM with arbitrary-way clustered inference.

    The analytic Jacobian, when supplied, must equal the negative derivative
    of the average score with respect to ``theta`` and have shape q-by-p.
    """

    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")
    if len(data) == 0:
        raise ValueError("data must contain at least one observation")
    if (nuisance is None) == (nuisance_fit is None):
        raise ValueError("provide exactly one of nuisance or nuisance_fit")
    selected_weight = _validate_variance_type(weight_type, "weight_type")
    selected_covariance = _validate_variance_type(
        covariance_type, "covariance_type"
    )
    if n_steps not in {1, 2}:
        raise ValueError("n_steps must be 1 or 2")

    start = np.asarray(theta_start, dtype=float)
    if start.ndim != 1 or len(start) == 0 or not np.isfinite(start).all():
        raise ValueError("theta_start must be a finite, nonempty one-dimensional array")

    if isinstance(nuisance_fit, SklearnNuisanceFitter):
        fitted_nuisance = nuisance_fit(data, cluster_cols=cluster_cols)
    elif nuisance_fit is not None:
        fitted_nuisance = nuisance_fit(data)
    else:
        fitted_nuisance = nuisance
    initial_scores = _evaluate_score(score, data, start, fitted_nuisance)
    q = initial_scores.shape[1]
    p = len(start)
    if q < p:
        raise ValueError(
            "the number of moments must be at least the number of parameters"
        )

    identity_weight = np.eye(q)
    first_theta, first_diagnostics = _optimize(
        score,
        data,
        fitted_nuisance,
        start,
        identity_weight,
        q,
        optimizer_method,
        optimizer_options,
    )
    first_scores = _evaluate_score(
        score, data, first_theta, fitted_nuisance, q
    )
    preliminary = compute_cluster_moment_matrix(
        first_scores, data, cluster_cols, selected_weight
    )

    if n_steps == 2:
        final_weight = _positive_definite_inverse(
            preliminary.matrix, selected_weight
        )
        theta, second_diagnostics = _optimize(
            score,
            data,
            fitted_nuisance,
            first_theta,
            final_weight,
            q,
            optimizer_method,
            optimizer_options,
        )
    else:
        final_weight = identity_weight
        theta = first_theta.copy()
        second_diagnostics = None

    final_scores = _evaluate_score(score, data, theta, fitted_nuisance, q)
    final_middle = compute_cluster_moment_matrix(
        final_scores, data, cluster_cols, selected_covariance
    )
    estimated_jacobian = _jacobian(
        jacobian,
        score,
        data,
        theta,
        fitted_nuisance,
        q,
        p,
    )
    covariance = _coefficient_covariance(
        estimated_jacobian,
        final_weight,
        final_middle.matrix,
        final_middle.n_eff,
    )
    if selected_covariance == "PSD":
        covariance_eigenvalues = np.linalg.eigvalsh(covariance)
        covariance_tolerance = (
            np.finfo(float).eps
            * max(1, p)
            * max(1.0, float(np.max(np.abs(covariance_eigenvalues))))
            * 100.0
        )
        if float(covariance_eigenvalues.min()) < -covariance_tolerance:
            raise RuntimeError("PSD produced a non-positive-semidefinite covariance")

    se, ci = _standard_errors_and_ci(theta, covariance, selected_covariance)
    return GMMResult(
        theta=theta,
        se=se,
        ci=ci,
        covariance=covariance,
        moment=final_scores.mean(axis=0),
        jacobian=estimated_jacobian,
        middle_matrix=final_middle.matrix,
        preliminary_middle_matrix=preliminary.matrix,
        weight_matrix=final_weight,
        first_step_theta=first_theta,
        n_eff=final_middle.n_eff,
        cluster_counts=final_middle.cluster_counts,
        n_cells=final_middle.n_cells,
        n_observations=final_middle.n_observations,
        nuisance=fitted_nuisance,
        n_steps=n_steps,
        weight_type=selected_weight,
        covariance_type=selected_covariance,
        first_step_optimizer=first_diagnostics,
        second_step_optimizer=second_diagnostics,
    )
