#!/usr/bin/env python3
"""Overidentified full-sample GMM on the project's synthetic data."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

from fullsample_gmm import (
    MultiwayGroupKFold,
    NuisanceSpec,
    fit_gmm,
    sklearn_nuisance_fit,
)


DATA_PATH = Path(__file__).resolve().parent / "data" / "synthetic_gmm.csv"
FEATURES = [f"x{index}" for index in range(1, 6)]
CLUSTER_COLS = ["row_id", "col_id"]

# Change only this line to switch the nuisance-tuning scheme.
CV_MODE = "oneway"  # available values: "oneway", "multiway"
MULTIWAY_FOLDS_PER_DIMENSION = 3


def tuned_ridge(
    alphas: list[float],
    one_way_folds: int,
    cv_mode: str,
) -> GridSearchCV:
    """Return an independently configured grouped nuisance search."""

    estimator = make_pipeline(
        PolynomialFeatures(degree=2, include_bias=False),
        StandardScaler(),
        Ridge(),
    )
    if cv_mode == "oneway":
        cv = GroupKFold(n_splits=one_way_folds)
    elif cv_mode == "multiway":
        cv = MultiwayGroupKFold(
            n_splits_per_dimension=MULTIWAY_FOLDS_PER_DIMENSION,
            shuffle=True,
            random_state=42,
        )
    else:
        raise ValueError("CV_MODE must be 'oneway' or 'multiway'")

    return GridSearchCV(
        estimator,
        {"ridge__alpha": alphas},
        cv=cv,
        scoring="neg_mean_squared_error",
        refit=True,
    )


# The tuning parameters demonstrate that each nuisance learner can have its
# own tuning configuration.
def main(cv_mode: str = CV_MODE) -> None:
    data = pd.read_csv(DATA_PATH)
    tuning_groups = None if cv_mode == "oneway" else CLUSTER_COLS
    specs = {
        "ell": NuisanceSpec(
            estimator=tuned_ridge([0.1, 1.0, 10.0], 5, cv_mode),
            target="y",
            features=FEATURES,
            groups=tuning_groups,
        ),
        "r": NuisanceSpec(
            estimator=tuned_ridge([0.01, 0.1, 1.0], 4, cv_mode),
            target="d",
            features=FEATURES,
            groups=tuning_groups,
        ),
        "m1": NuisanceSpec(
            estimator=tuned_ridge([0.1, 1.0, 5.0], 3, cv_mode),
            target="z1",
            features=FEATURES,
            groups=tuning_groups,
        ),
        "m2": NuisanceSpec(
            estimator=tuned_ridge([0.5, 2.0, 10.0], 6, cv_mode),
            target="z2",
            features=FEATURES,
            groups=tuning_groups,
        ),
    }

    def score(frame, theta, fitted):
        predictions = fitted.predictions
        y_residual = frame["y"].to_numpy() - predictions["ell"]
        d_residual = frame["d"].to_numpy() - predictions["r"]
        z_residual = np.column_stack(
            (
                frame["z1"].to_numpy() - predictions["m1"],
                frame["z2"].to_numpy() - predictions["m2"],
            )
        )
        return z_residual * (y_residual - theta[0] * d_residual)[:, None]

    # this particular model has a Jacobian that does not depend on theta
    # but the function must accept theta because the package expects it
    def jacobian(frame, theta, fitted):
        del theta
        predictions = fitted.predictions
        d_residual = frame["d"].to_numpy() - predictions["r"]
        z_residual = np.column_stack(
            (
                frame["z1"].to_numpy() - predictions["m1"],
                frame["z2"].to_numpy() - predictions["m2"],
            )
        )
        return np.mean(z_residual * d_residual[:, None], axis=0)[:, None]

    result = fit_gmm(
        data=data,
        score=score,
        theta_start=[0.0],
        cluster_cols=CLUSTER_COLS,
        nuisance_fit=sklearn_nuisance_fit(specs),
        jacobian=jacobian,
        weight_type="PSD",
        covariance_type="PSD",
    )

    print("Full-sample overidentified GMM")
    print(f"Nuisance tuning CV: {cv_mode}")
    print(f"Observations: {result.n_observations}; cells: {result.n_cells}")
    print(f"Clusters: {result.cluster_counts}; n_eff: {result.n_eff}")
    tuning_columns = result.nuisance.group_columns["ell"]
    tuning_counts = result.nuisance.group_counts["ell"]
    print(f"Tuning groups: {tuning_columns}; counts: {tuning_counts}")
    print(f"First step: {result.first_step_theta[0]:.6f}")
    print(f"Two step:   {result.theta[0]:.6f}")
    print(f"PSD SE:     {result.se[0]:.6f}")
    print(f"95% CI:     [{result.ci[0, 0]:.6f}, {result.ci[0, 1]:.6f}]")
    print("Selected nuisance tuning:")
    for name, model in result.nuisance.models.items():
        print(f"  {name}: {model.best_params_}")


if __name__ == "__main__":
    main()
