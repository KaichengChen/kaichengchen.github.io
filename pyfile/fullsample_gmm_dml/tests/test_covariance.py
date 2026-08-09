from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
import pytest

from fullsample_gmm import compute_cluster_moment_matrix


CLUSTERS = pd.DataFrame(
    {
        "one": [0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 3],
        "two": [0, 0, 1, 0, 1, 0, 1, 2, 0, 1, 2, 2],
        "three": ["a", "b", "a", "b", "c", "a", "c", "b", "a", "c", "b", "c"],
        "four": [True, False, True, True, False, False, True, False, True, False, True, False],
    }
)
SCORES = np.random.default_rng(2026).normal(size=(len(CLUSTERS), 3))


def brute_force(scores, data, columns, method):
    labels = data.loc[:, columns].to_numpy()
    q = scores.shape[1]
    total = np.zeros((q, q))
    if method == "PSD":
        subsets = [((index,), 1.0) for index in range(len(columns))]
    else:
        subsets = []
        for size in range(1, len(columns) + 1):
            sign = 1.0 if size % 2 else -1.0
            subsets.extend(
                (subset, sign)
                for subset in combinations(range(len(columns)), size)
            )

    for subset, sign in subsets:
        for left in range(len(data)):
            for right in range(len(data)):
                if all(labels[left, level] == labels[right, level] for level in subset):
                    total += sign * np.outer(scores[left], scores[right])
    counts = [data[column].nunique() for column in columns]
    return min(counts) / len(data) ** 2 * total


@pytest.mark.parametrize("n_dimensions", [1, 2, 3, 4])
@pytest.mark.parametrize("method", ["PSD", "CGM"])
def test_k_way_matrix_matches_pairwise_definition(n_dimensions, method):
    columns = list(CLUSTERS.columns[:n_dimensions])
    result = compute_cluster_moment_matrix(SCORES, CLUSTERS, columns, method)
    expected = brute_force(SCORES, CLUSTERS, columns, method)

    np.testing.assert_allclose(result.matrix, expected, atol=1e-12)
    assert result.n_eff == min(CLUSTERS[column].nunique() for column in columns)
    assert result.n_cells == len(CLUSTERS.groupby(columns))


@pytest.mark.parametrize("n_dimensions", [1, 2, 3, 4])
def test_psd_is_positive_semidefinite(n_dimensions):
    columns = list(CLUSTERS.columns[:n_dimensions])
    matrix = compute_cluster_moment_matrix(
        SCORES, CLUSTERS, columns, "PSD"
    ).matrix
    assert np.linalg.eigvalsh(matrix).min() >= -1e-12


def test_permutation_and_cluster_relabeling_do_not_change_matrix():
    columns = list(CLUSTERS.columns)
    original = compute_cluster_moment_matrix(
        SCORES, CLUSTERS, columns, "PSD"
    ).matrix
    order = np.random.default_rng(4).permutation(len(CLUSTERS))
    permuted_data = CLUSTERS.iloc[order].reset_index(drop=True).copy()
    permuted_scores = SCORES[order]
    permuted_data["one"] = permuted_data["one"].map(
        {0: "north", 1: "south", 2: "east", 3: "west"}
    )
    changed = compute_cluster_moment_matrix(
        permuted_scores, permuted_data, columns, "PSD"
    ).matrix
    np.testing.assert_allclose(changed, original, atol=1e-12)


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda frame: frame.assign(one=np.nan), "missing"),
        (lambda frame: frame.assign(one=1), "at least two"),
    ],
)
def test_invalid_cluster_identifiers_are_rejected(change, message):
    with pytest.raises(ValueError, match=message):
        compute_cluster_moment_matrix(SCORES, change(CLUSTERS), ["one"], "PSD")


def test_invalid_cluster_columns_and_method_are_rejected():
    with pytest.raises(ValueError, match="repeated"):
        compute_cluster_moment_matrix(SCORES, CLUSTERS, ["one", "one"], "PSD")
    with pytest.raises(ValueError, match="not found"):
        compute_cluster_moment_matrix(SCORES, CLUSTERS, ["absent"], "PSD")
    with pytest.raises(ValueError, match="exactly"):
        compute_cluster_moment_matrix(SCORES, CLUSTERS, ["one"], "psd")


def test_invalid_scores_are_rejected():
    with pytest.raises(ValueError, match="two-dimensional"):
        compute_cluster_moment_matrix(SCORES[:, 0], CLUSTERS, ["one"], "PSD")
    invalid = SCORES.copy()
    invalid[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        compute_cluster_moment_matrix(invalid, CLUSTERS, ["one"], "PSD")

