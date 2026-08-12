from __future__ import annotations

from itertools import product

import numpy as np
import pytest

from fullsampleDML import MultiwayGroupKFold


def complete_groups(n_dimensions: int, n_clusters: int = 6) -> np.ndarray:
    return np.asarray(list(product(range(n_clusters), repeat=n_dimensions)))


@pytest.mark.parametrize("n_dimensions", [1, 2, 3])
def test_multiway_splits_have_no_cluster_overlap_and_full_validation_coverage(
    n_dimensions,
):
    groups = complete_groups(n_dimensions)
    splitter = MultiwayGroupKFold(n_splits_per_dimension=3)
    validation_counts = np.zeros(len(groups), dtype=int)
    splits = list(splitter.split(np.zeros((len(groups), 1)), groups=groups))

    assert len(splits) == 3**n_dimensions
    assert splitter.get_n_splits(groups=groups) == 3**n_dimensions
    for training, validation in splits:
        validation_counts[validation] += 1
        assert not np.intersect1d(training, validation).size
        for dimension in range(n_dimensions):
            assert not set(groups[training, dimension]) & set(
                groups[validation, dimension]
            )
    np.testing.assert_array_equal(validation_counts, np.ones(len(groups), dtype=int))


def test_multiway_splits_support_different_fold_counts():
    groups = complete_groups(2)
    splitter = MultiwayGroupKFold(n_splits_per_dimension=(2, 3))
    assert len(list(splitter.split(groups, groups=groups))) == 6
    assert splitter.get_n_splits(groups=groups) == 6


def test_multiway_splits_support_heterogeneous_and_missing_cells():
    groups = complete_groups(2, n_clusters=4)
    keep = ~(
        ((groups[:, 0] == 0) & (groups[:, 1] == 0))
        | ((groups[:, 0] == 2) & (groups[:, 1] == 2))
    )
    sparse = groups[keep]
    heterogeneous = np.vstack((sparse, sparse[:5], sparse[:2]))
    splitter = MultiwayGroupKFold(n_splits_per_dimension=2)
    splits = list(splitter.split(heterogeneous, groups=heterogeneous))

    assert len(splits) == 4
    validation_counts = np.zeros(len(heterogeneous), dtype=int)
    for training, validation in splits:
        validation_counts[validation] += 1
        for dimension in range(2):
            assert not set(heterogeneous[training, dimension]) & set(
                heterogeneous[validation, dimension]
            )
    np.testing.assert_array_equal(
        validation_counts, np.ones(len(heterogeneous), dtype=int)
    )


def test_shuffled_splits_are_reproducible():
    groups = complete_groups(2)
    first = MultiwayGroupKFold(3, shuffle=True, random_state=19)
    second = MultiwayGroupKFold(3, shuffle=True, random_state=19)
    first_splits = list(first.split(groups, groups=groups))
    second_splits = list(second.split(groups, groups=groups))
    for left, right in zip(first_splits, second_splits, strict=True):
        np.testing.assert_array_equal(left[0], right[0])
        np.testing.assert_array_equal(left[1], right[1])


def test_multiway_validation_rejects_invalid_inputs():
    groups = complete_groups(2, n_clusters=3).astype(float)
    with pytest.raises(ValueError, match="required"):
        list(MultiwayGroupKFold(2).split(groups))
    with pytest.raises(ValueError, match="same number"):
        list(MultiwayGroupKFold(2).split(groups[:-1], groups=groups))
    invalid = groups.copy()
    invalid[0, 0] = np.nan
    with pytest.raises(ValueError, match="missing"):
        list(MultiwayGroupKFold(2).split(groups, groups=invalid))
    with pytest.raises(ValueError, match="requires at least 4"):
        list(MultiwayGroupKFold(4).split(groups, groups=groups))
    with pytest.raises(ValueError, match="one value per"):
        list(MultiwayGroupKFold((2, 2, 2)).split(groups, groups=groups))


def test_empty_validation_and_training_splits_are_rejected():
    diagonal = np.array([[0, 0], [1, 1]])
    with pytest.raises(ValueError, match="validation intersection is empty"):
        list(MultiwayGroupKFold(2).split(diagonal, groups=diagonal))

    missing_training_corner = np.array([[0, 0], [0, 1], [1, 0]])
    with pytest.raises(ValueError, match="training set is empty"):
        list(
            MultiwayGroupKFold(2).split(
                missing_training_corner, groups=missing_training_corner
            )
        )


def test_string_groups_are_supported():
    groups = complete_groups(2, n_clusters=4).astype(str)
    splits = list(MultiwayGroupKFold(2).split(groups, groups=groups))
    assert len(splits) == 4
