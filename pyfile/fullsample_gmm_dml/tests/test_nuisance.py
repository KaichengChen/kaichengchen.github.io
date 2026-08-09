from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from fullsample_gmm import (
    MultiwayGroupKFold,
    NuisanceSpec,
    fit_gmm,
    fit_sklearn_nuisances,
    sklearn_nuisance_fit,
)


class RecordingRegressor(RegressorMixin, BaseEstimator):
    fit_sizes: list[int] = []

    def __init__(self, offset=0.0):
        self.offset = offset

    def fit(self, features, target):
        type(self).fit_sizes.append(len(features))
        self.mean_ = float(np.mean(target)) + self.offset
        return self

    def predict(self, features):
        return np.full(len(features), self.mean_)


class NonfiniteRegressor(RegressorMixin, BaseEstimator):
    def fit(self, features, target):
        return self

    def predict(self, features):
        return np.full(len(features), np.nan)


@pytest.fixture
def data():
    x = np.linspace(-1.0, 1.0, 18)
    return pd.DataFrame(
        {
            "x": x,
            "x2": x**2,
            "y": 1.0 + x,
            "d": x - x**2,
            "large_group": np.tile(np.arange(6), 3),
            "small_group": np.repeat(np.arange(3), 6),
            "small_group_tie": np.tile(np.arange(3), 6),
        }
    )


def test_separate_searches_refit_on_the_complete_sample(data):
    RecordingRegressor.fit_sizes = []
    specs = {
        "ell": NuisanceSpec(
            estimator=GridSearchCV(
                RecordingRegressor(), {"offset": [-0.1, 0.0, 0.1]}, cv=3, refit=True
            ),
            target="y",
            features=["x", "x2"],
        ),
        "m": NuisanceSpec(
            estimator=GridSearchCV(
                RecordingRegressor(), {"offset": [-0.2, 0.0]}, cv=2, refit=True
            ),
            target="d",
            features=lambda frame: frame[["x"]],
        ),
    }
    result = fit_sklearn_nuisances(data, specs)

    assert set(result.models) == {"ell", "m"}
    assert set(result.predictions) == {"ell", "m"}
    assert result.predictions["ell"].shape == (len(data),)
    assert result.predictions["m"].shape == (len(data),)
    assert RecordingRegressor.fit_sizes.count(len(data)) == 2
    assert result.models["ell"].cv == 3
    assert result.models["m"].cv == 2


def test_callback_wrapper_returns_full_sample_results(data):
    callback = sklearn_nuisance_fit(
        {
            "ell": NuisanceSpec(
                estimator=RecordingRegressor(), target="y", features="x"
            )
        }
    )
    result = callback(data)
    assert result.predictions["ell"].shape == (len(data),)


def grouped_spec(n_splits=3, groups=None):
    return NuisanceSpec(
        estimator=GridSearchCV(
            RecordingRegressor(),
            {"offset": [-0.1, 0.0, 0.1]},
            cv=GroupKFold(n_splits=n_splits),
            refit=True,
        ),
        target="y",
        features=["x", "x2"],
        groups=groups,
    )


def test_smallest_cluster_dimension_is_selected_automatically(data):
    result = fit_sklearn_nuisances(
        data,
        {"ell": grouped_spec()},
        cluster_cols=["large_group", "small_group"],
    )
    assert result.group_columns["ell"] == ("small_group",)
    assert result.group_counts["ell"] == {"small_group": 3}

    groups = data["small_group"].to_numpy()
    for training, validation in result.models["ell"].cv.split(data, groups=groups):
        assert not set(groups[training]) & set(groups[validation])


def test_cluster_count_ties_use_first_listed_dimension(data):
    result = fit_sklearn_nuisances(
        data,
        {"ell": grouped_spec()},
        cluster_cols=["small_group_tie", "small_group"],
    )
    assert result.group_columns["ell"] == ("small_group_tie",)


def test_explicit_groups_override_automatic_selection(data):
    result = fit_sklearn_nuisances(
        data,
        {"ell": grouped_spec(groups="large_group")},
        cluster_cols=["large_group", "small_group"],
    )
    assert result.group_columns["ell"] == ("large_group",)
    assert result.group_counts["ell"] == {"large_group": 6}


def test_explicit_multiway_groups_reach_grid_search(data):
    spec = NuisanceSpec(
        estimator=GridSearchCV(
            Ridge(),
            {"alpha": [0.1, 1.0]},
            cv=MultiwayGroupKFold(n_splits_per_dimension=(2, 3)),
            refit=True,
        ),
        target="y",
        features=["x", "x2"],
        groups=["large_group", "small_group"],
    )
    result = fit_sklearn_nuisances(
        data, {"ell": spec}, cluster_cols=["large_group", "small_group"]
    )
    assert result.group_columns["ell"] == ("large_group", "small_group")
    assert result.models["ell"].best_params_["alpha"] in {0.1, 1.0}


def test_plain_pipeline_does_not_receive_tuning_groups(data):
    spec = NuisanceSpec(
        estimator=make_pipeline(StandardScaler(), Ridge()),
        target="y",
        features=["x", "x2"],
    )
    result = fit_sklearn_nuisances(
        data,
        {"ell": spec},
        cluster_cols=["large_group", "small_group"],
    )
    assert result.predictions["ell"].shape == (len(data),)
    assert result.group_columns["ell"] == ()


def test_fit_gmm_passes_cluster_columns_to_sklearn_adapter(data):
    result = fit_gmm(
        data=data,
        score=lambda frame, theta, nuisance: (
            frame["y"].to_numpy() - theta[0]
        )[:, None],
        theta_start=[0.0],
        cluster_cols=["large_group", "small_group"],
        nuisance_fit=sklearn_nuisance_fit({"ell": grouped_spec()}),
        jacobian=lambda frame, theta, nuisance: np.ones((1, 1)),
        n_steps=1,
    )
    assert result.nuisance.group_columns["ell"] == ("small_group",)
    assert result.nuisance.group_counts["ell"] == {"small_group": 3}


def test_group_column_validation(data):
    missing = data.copy()
    missing.loc[0, "small_group"] = np.nan
    with pytest.raises(ValueError, match="must not contain missing"):
        fit_sklearn_nuisances(
            missing,
            {"ell": grouped_spec()},
            cluster_cols=["large_group", "small_group"],
        )
    with pytest.raises(ValueError, match="not found"):
        fit_sklearn_nuisances(
            data,
            {"ell": grouped_spec(groups="absent")},
            cluster_cols=["large_group", "small_group"],
        )


def test_search_with_refit_false_is_rejected(data):
    spec = NuisanceSpec(
        estimator=GridSearchCV(
            RecordingRegressor(), {"offset": [0.0, 0.1]}, cv=2, refit=False
        ),
        target="y",
        features=["x"],
    )
    with pytest.raises(ValueError, match="refit=False"):
        fit_sklearn_nuisances(data, {"ell": spec})


def test_nonfinite_predictions_are_rejected(data):
    spec = NuisanceSpec(
        estimator=NonfiniteRegressor(), target="y", features=["x"]
    )
    with pytest.raises(ValueError, match="non-finite"):
        fit_sklearn_nuisances(data, {"ell": spec})


def test_invalid_nuisance_specs_are_rejected(data):
    with pytest.raises(ValueError, match="at least one"):
        fit_sklearn_nuisances(data, {})
    with pytest.raises(ValueError, match="not found"):
        fit_sklearn_nuisances(
            data,
            {
                "ell": NuisanceSpec(
                    estimator=RecordingRegressor(), target="absent", features=["x"]
                )
            },
        )
