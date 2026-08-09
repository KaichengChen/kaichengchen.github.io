"""Cluster-aware cross-validation splitters for nuisance tuning."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from itertools import product
from math import prod

import numpy as np
import pandas as pd


class MultiwayGroupKFold:
    """Strict arbitrary-way grouped cross-validation.

    Validation observations lie in an intersection of held-out cluster bins.
    Training observations must lie outside the held-out bin in every cluster
    dimension. Observations in the intervening border regions are omitted from
    that split, ensuring that training and validation share no cluster label.
    """

    def __init__(
        self,
        n_splits_per_dimension: int | Sequence[int] = 3,
        *,
        shuffle: bool = False,
        random_state: int | None = None,
    ) -> None:
        self.n_splits_per_dimension = n_splits_per_dimension
        self.shuffle = shuffle
        self.random_state = random_state

    def _validate_groups(
        self,
        groups: object,
        n_samples: int | None = None,
    ) -> np.ndarray:
        if groups is None:
            raise ValueError("groups are required for MultiwayGroupKFold")
        values = np.asarray(groups, dtype=object)
        if values.ndim == 1:
            values = values[:, None]
        if values.ndim != 2 or values.shape[1] == 0:
            raise ValueError("groups must be an n-by-K array with K >= 1")
        if n_samples is not None and values.shape[0] != n_samples:
            raise ValueError("groups and X must contain the same number of rows")
        if values.shape[0] == 0:
            raise ValueError("groups must contain at least one observation")
        if pd.isna(values).any():
            raise ValueError("groups must not contain missing values")
        return values

    def _fold_counts(self, n_dimensions: int) -> tuple[int, ...]:
        requested = self.n_splits_per_dimension
        if isinstance(requested, (int, np.integer)):
            counts = (int(requested),) * n_dimensions
        else:
            counts = tuple(int(value) for value in requested)
            if len(counts) != n_dimensions:
                raise ValueError(
                    "n_splits_per_dimension must have one value per group dimension"
                )
        if any(count < 2 for count in counts):
            raise ValueError("every dimension must use at least two folds")
        return counts

    def _observation_bins(
        self,
        groups: np.ndarray,
        fold_counts: tuple[int, ...],
    ) -> np.ndarray:
        rng = np.random.default_rng(self.random_state)
        bins = np.empty(groups.shape, dtype=int)
        for dimension, n_folds in enumerate(fold_counts):
            codes, unique_groups = pd.factorize(groups[:, dimension], sort=False)
            if len(unique_groups) < n_folds:
                raise ValueError(
                    f"group dimension {dimension} has {len(unique_groups)} distinct "
                    f"clusters but requires at least {n_folds}"
                )
            cluster_order = np.arange(len(unique_groups))
            if self.shuffle:
                rng.shuffle(cluster_order)
            cluster_bins = np.empty(len(unique_groups), dtype=int)
            for fold, members in enumerate(np.array_split(cluster_order, n_folds)):
                cluster_bins[members] = fold
            bins[:, dimension] = cluster_bins[codes]
        return bins

    def split(
        self,
        X: object,
        y: object = None,
        groups: object = None,
    ) -> Iterator[tuple[np.ndarray, np.ndarray]]:
        """Yield strict cluster-disjoint training and validation indices."""

        del y
        try:
            n_samples = len(X)  # type: ignore[arg-type]
        except TypeError as exc:
            raise ValueError("X must have a defined number of rows") from exc
        group_values = self._validate_groups(groups, n_samples)
        fold_counts = self._fold_counts(group_values.shape[1])
        observation_bins = self._observation_bins(group_values, fold_counts)

        for held_out in product(*(range(count) for count in fold_counts)):
            held_out_array = np.asarray(held_out)
            validation = np.flatnonzero(
                np.all(observation_bins == held_out_array, axis=1)
            )
            if len(validation) == 0:
                raise ValueError(
                    "a multiway validation intersection is empty; reduce the number "
                    "of folds per dimension"
                )
            training = np.flatnonzero(
                np.all(observation_bins != held_out_array, axis=1)
            )
            if len(training) == 0:
                raise ValueError(
                    "a multiway training set is empty; adjust the folds or clusters"
                )
            yield training, validation

    def get_n_splits(
        self,
        X: object = None,
        y: object = None,
        groups: object = None,
    ) -> int:
        """Return the Cartesian product of dimension-specific fold counts."""

        del X, y
        group_values = self._validate_groups(groups)
        fold_counts = self._fold_counts(group_values.shape[1])
        for dimension, n_folds in enumerate(fold_counts):
            distinct = len(pd.unique(group_values[:, dimension]))
            if distinct < n_folds:
                raise ValueError(
                    f"group dimension {dimension} has {distinct} distinct clusters "
                    f"but requires at least {n_folds}"
                )
        return prod(fold_counts)

