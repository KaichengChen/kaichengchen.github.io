"""Convenience utilities for full-sample scikit-learn nuisance fits."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from inspect import Parameter, signature
from typing import Any

import numpy as np
import pandas as pd
from sklearn.base import clone


ColumnSelector = str | Sequence[str] | Callable[[pd.DataFrame], Any]
GroupSelector = str | Sequence[str]


@dataclass(frozen=True)
class NuisanceSpec:
    """Specification for one independently configured nuisance learner."""

    estimator: Any
    target: ColumnSelector
    features: ColumnSelector
    groups: GroupSelector | None = None
    predict_method: str = "predict"
    fit_params: Mapping[str, Any] = field(default_factory=dict)
    predict_params: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class NuisanceResult:
    """Fitted nuisance models and their complete-sample predictions."""

    predictions: dict[str, np.ndarray]
    models: dict[str, Any]
    group_columns: dict[str, tuple[str, ...]]
    group_counts: dict[str, dict[str, int]]


@dataclass(frozen=True)
class SklearnNuisanceFitter:
    """Callable adapter that receives inference clusters from ``fit_gmm``."""

    specs: Mapping[str, NuisanceSpec]

    def __call__(
        self,
        data: pd.DataFrame,
        cluster_cols: Sequence[str] | None = None,
    ) -> NuisanceResult:
        return fit_sklearn_nuisances(data, self.specs, cluster_cols=cluster_cols)


def _select(data: pd.DataFrame, selector: ColumnSelector, role: str) -> Any:
    if callable(selector):
        return selector(data)
    if isinstance(selector, str):
        if selector not in data.columns:
            raise ValueError(f"{role} column not found in data: {selector!r}")
        if role == "target":
            return data[selector]
        return data[[selector]]

    columns = list(selector)
    if not columns:
        raise ValueError(f"{role} columns must not be empty")
    missing = [column for column in columns if column not in data.columns]
    if missing:
        raise ValueError(f"{role} columns not found in data: {missing}")
    selected = data.loc[:, columns]
    if role == "target" and len(columns) == 1:
        return selected.iloc[:, 0]
    return selected


def _validate_predictions(name: str, predictions: Any, n_rows: int) -> np.ndarray:
    array = np.asarray(predictions)
    if array.ndim == 0 or array.shape[0] != n_rows:
        raise ValueError(
            f"nuisance {name!r} returned predictions with an invalid first dimension"
        )
    try:
        finite = np.isfinite(array.astype(float, copy=False)).all()
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"nuisance {name!r} predictions must be numeric"
        ) from exc
    if not finite:
        raise ValueError(f"nuisance {name!r} returned non-finite predictions")
    return array


def _validate_group_columns(
    data: pd.DataFrame,
    columns: Sequence[str],
    role: str,
) -> tuple[tuple[str, ...], dict[str, int]]:
    selected = tuple(columns)
    if not selected:
        raise ValueError(f"{role} must contain at least one column")
    if len(set(selected)) != len(selected):
        raise ValueError(f"{role} must not contain repeated column names")
    missing = [column for column in selected if column not in data.columns]
    if missing:
        raise ValueError(f"{role} not found in data: {missing}")
    if data.loc[:, list(selected)].isna().any(axis=None):
        raise ValueError(f"{role} must not contain missing values")
    counts = {
        column: int(data[column].nunique(dropna=False)) for column in selected
    }
    return selected, counts


def _resolve_groups(
    data: pd.DataFrame,
    groups: GroupSelector | None,
    cluster_cols: Sequence[str] | None,
) -> tuple[np.ndarray | None, tuple[str, ...], dict[str, int]]:
    if groups is None:
        if cluster_cols is None:
            return None, (), {}
        candidates, counts = _validate_group_columns(
            data, cluster_cols, "cluster_cols"
        )
        # ``min`` retains the first candidate when cluster counts tie.
        selected = min(candidates, key=counts.__getitem__)
        columns = (selected,)
        selected_counts = {selected: counts[selected]}
    else:
        requested = (groups,) if isinstance(groups, str) else tuple(groups)
        columns, selected_counts = _validate_group_columns(
            data, requested, "group columns"
        )

    values = data.loc[:, list(columns)].to_numpy()
    if len(columns) == 1:
        values = values[:, 0]
    return values, columns, selected_counts


def _fit_accepts_groups(model: Any) -> bool:
    parameters = tuple(signature(model.fit).parameters.values())
    if any(parameter.name == "groups" for parameter in parameters):
        return True
    # Search estimators accept group metadata through **params. Do not send
    # groups to ordinary pipelines, which also expose **params but reject it.
    accepts_keywords = any(
        parameter.kind == Parameter.VAR_KEYWORD for parameter in parameters
    )
    return accepts_keywords and hasattr(model, "cv")


def fit_sklearn_nuisances(
    data: pd.DataFrame,
    specs: Mapping[str, NuisanceSpec],
    cluster_cols: Sequence[str] | None = None,
) -> NuisanceResult:
    """Tune, refit, and predict each nuisance learner on the complete sample."""

    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")
    if not specs:
        raise ValueError("specs must contain at least one nuisance learner")

    predictions: dict[str, np.ndarray] = {}
    models: dict[str, Any] = {}
    group_columns: dict[str, tuple[str, ...]] = {}
    group_counts: dict[str, dict[str, int]] = {}
    for name, spec in specs.items():
        if not isinstance(name, str) or not name:
            raise ValueError("every nuisance name must be a nonempty string")
        if getattr(spec.estimator, "refit", True) is False:
            raise ValueError(
                f"nuisance {name!r} has refit=False; search estimators must "
                "refit the selected learner on the complete sample"
            )

        features = _select(data, spec.features, "feature")
        target = _select(data, spec.target, "target")
        model = clone(spec.estimator)
        groups, selected_columns, selected_counts = _resolve_groups(
            data, spec.groups, cluster_cols
        )
        fit_params = dict(spec.fit_params)
        if "groups" in fit_params:
            raise ValueError(
                f"nuisance {name!r} must specify groups through NuisanceSpec.groups"
            )
        if groups is not None and _fit_accepts_groups(model):
            fit_params["groups"] = groups
            group_columns[name] = selected_columns
            group_counts[name] = selected_counts
        else:
            group_columns[name] = ()
            group_counts[name] = {}
        model.fit(features, target, **fit_params)
        if not hasattr(model, spec.predict_method):
            raise ValueError(
                f"nuisance {name!r} has no {spec.predict_method!r} method"
            )
        predictor = getattr(model, spec.predict_method)
        predicted = predictor(features, **dict(spec.predict_params))
        predictions[name] = _validate_predictions(name, predicted, len(data))
        models[name] = model

    return NuisanceResult(
        predictions=predictions,
        models=models,
        group_columns=group_columns,
        group_counts=group_counts,
    )


def sklearn_nuisance_fit(
    specs: Mapping[str, NuisanceSpec],
) -> SklearnNuisanceFitter:
    """Build a ``nuisance_fit`` callback from scikit-learn specifications."""

    return SklearnNuisanceFitter(specs=dict(specs))
