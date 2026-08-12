"""Arbitrary-way clustered score middle matrices."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Literal, Sequence

import numpy as np
import pandas as pd


VarianceType = Literal["PSD", "CGM"]


@dataclass(frozen=True)
class ClusterMomentResult:
    """A clustered score middle matrix and its sampling metadata."""

    matrix: np.ndarray
    method: VarianceType
    n_eff: int
    cluster_counts: dict[str, int]
    n_cells: int
    n_observations: int


def _validate_method(method: str) -> VarianceType:
    if method not in {"PSD", "CGM"}:
        raise ValueError("method must be exactly 'PSD' or 'CGM'")
    return method  # type: ignore[return-value]


def _validate_inputs(
    scores: np.ndarray,
    data: pd.DataFrame,
    cluster_cols: Sequence[str],
) -> tuple[np.ndarray, tuple[str, ...], dict[str, int]]:
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")

    columns = tuple(cluster_cols)
    if not columns:
        raise ValueError("cluster_cols must contain at least one column")
    if len(set(columns)) != len(columns):
        raise ValueError("cluster_cols must not contain repeated column names")
    missing = [column for column in columns if column not in data.columns]
    if missing:
        raise ValueError(f"cluster columns not found in data: {missing}")
    if data.loc[:, list(columns)].isna().any(axis=None):
        raise ValueError("cluster columns must not contain missing values")

    array = np.asarray(scores, dtype=float)
    if array.ndim != 2:
        raise ValueError("scores must be a two-dimensional n-by-q array")
    if len(data) == 0:
        raise ValueError("data and scores must contain at least one observation")
    if array.shape[0] != len(data):
        raise ValueError("scores and data must have the same number of rows")
    if array.shape[1] == 0:
        raise ValueError("scores must contain at least one moment")
    if not np.isfinite(array).all():
        raise ValueError("scores must contain only finite values")

    counts = {
        column: int(data[column].nunique(dropna=False)) for column in columns
    }
    too_small = [column for column, count in counts.items() if count < 2]
    if too_small:
        raise ValueError(
            "each clustering dimension must contain at least two clusters; "
            f"invalid columns: {too_small}"
        )
    return array, columns, counts


def _cell_scores(
    scores: np.ndarray,
    data: pd.DataFrame,
    cluster_cols: tuple[str, ...],
) -> pd.DataFrame:
    index = pd.MultiIndex.from_frame(
        data.loc[:, list(cluster_cols)].reset_index(drop=True)
    )
    score_frame = pd.DataFrame(scores, index=index)
    return score_frame.groupby(
        level=list(range(len(cluster_cols))), sort=False, observed=True
    ).sum()


def _subset_outer_sum(
    cell_scores: pd.DataFrame,
    levels: tuple[int, ...],
) -> np.ndarray:
    grouped = cell_scores.groupby(
        level=list(levels), sort=False, observed=True
    ).sum()
    values = grouped.to_numpy(dtype=float)
    return values.T @ values


def compute_cluster_moment_matrix(
    scores: np.ndarray,
    data: pd.DataFrame,
    cluster_cols: Sequence[str],
    method: VarianceType = "PSD",
) -> ClusterMomentResult:
    """Compute the paper's ``PSD`` or inclusion--exclusion ``CGM`` matrix.

    Unit scores are first summed within each observed full cluster cell. The
    returned matrix is normalized by ``n_eff / n_observations**2``, where
    ``n_eff`` is the smallest number of clusters among the dimensions.
    """

    selected_method = _validate_method(method)
    score_array, columns, counts = _validate_inputs(scores, data, cluster_cols)
    cells = _cell_scores(score_array, data, columns)
    n_dimensions = len(columns)
    q = score_array.shape[1]
    unscaled = np.zeros((q, q), dtype=float)

    if selected_method == "PSD":
        for level in range(n_dimensions):
            unscaled += _subset_outer_sum(cells, (level,))
    else:
        for subset_size in range(1, n_dimensions + 1):
            sign = 1.0 if subset_size % 2 == 1 else -1.0
            for levels in combinations(range(n_dimensions), subset_size):
                unscaled += sign * _subset_outer_sum(cells, levels)

    n_eff = min(counts.values())
    scale = n_eff / len(data) ** 2
    matrix = scale * unscaled
    matrix = (matrix + matrix.T) / 2.0
    return ClusterMomentResult(
        matrix=matrix,
        method=selected_method,
        n_eff=n_eff,
        cluster_counts=counts,
        n_cells=len(cells),
        n_observations=len(data),
    )

