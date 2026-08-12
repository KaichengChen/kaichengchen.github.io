from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from pathlib import Path

from fullsampleDML import fit_gmm


def test_cross_language_gmm_fixture():
    fixture_dir = Path(__file__).parent / "fixtures"
    data = pd.read_csv(fixture_dir / "golden_gmm_data.csv")
    expected_frame = pd.read_csv(fixture_dir / "golden_gmm_expected.csv")
    expected = dict(zip(expected_frame["quantity"], expected_frame["value"]))

    def score(frame, theta, nuisance):
        instruments = frame[["z1", "z2"]].to_numpy()
        residual = frame["y"].to_numpy() - theta[0] * frame["d"].to_numpy()
        return instruments * residual[:, None]

    def jacobian(frame, theta, nuisance):
        instruments = frame[["z1", "z2"]].to_numpy()
        return np.mean(
            instruments * frame["d"].to_numpy()[:, None], axis=0
        )[:, None]

    result = fit_gmm(
        data, score, [0.0], ["row", "col"], nuisance={}, jacobian=jacobian
    )
    assert result.first_step_theta[0] == pytest.approx(expected["first_step_theta"])
    assert result.theta[0] == pytest.approx(expected["theta"])
    assert result.se[0] == pytest.approx(expected["se"])
    assert result.covariance[0, 0] == pytest.approx(expected["covariance"])
    np.testing.assert_allclose(
        result.moment, [expected["moment_1"], expected["moment_2"]]
    )


def example_data():
    rng = np.random.default_rng(881)
    n_side = 20
    row = np.repeat(np.arange(n_side), n_side)
    column = np.tile(np.arange(n_side), n_side)
    z1 = rng.normal(size=n_side**2)
    z2 = 0.25 * z1 + rng.normal(size=n_side**2)
    d = 0.8 * z1 + 0.5 * z2 + rng.normal(size=n_side**2)
    error = (
        rng.normal(scale=0.4, size=n_side)[row]
        + rng.normal(scale=0.4, size=n_side)[column]
        + rng.normal(size=n_side**2)
    )
    y = 1.5 * d + error
    return pd.DataFrame(
        {"row": row, "column": column, "z1": z1, "z2": z2, "d": d, "y": y}
    )


def iv_score(data, theta, nuisance):
    del nuisance
    instruments = data[["z1", "z2"]].to_numpy()
    residual = data["y"].to_numpy() - theta[0] * data["d"].to_numpy()
    return instruments * residual[:, None]


def iv_jacobian(data, theta, nuisance):
    del theta, nuisance
    instruments = data[["z1", "z2"]].to_numpy()
    return np.mean(instruments * data["d"].to_numpy()[:, None], axis=0)[:, None]


def test_two_step_overidentified_gmm_and_result_metadata():
    data = example_data()
    calls = 0

    def fit_nuisance(frame):
        nonlocal calls
        calls += 1
        assert frame is data
        return {"fitted_on": len(frame)}

    result = fit_gmm(
        data,
        iv_score,
        [0.0],
        ["row", "column"],
        nuisance_fit=fit_nuisance,
        jacobian=iv_jacobian,
    )

    assert calls == 1
    assert result.theta[0] == pytest.approx(1.5, abs=0.15)
    assert result.se[0] > 0
    assert result.ci.shape == (1, 2)
    assert result.n_eff == 20
    assert result.cluster_counts == {"row": 20, "column": 20}
    assert result.n_cells == 400
    assert result.n_observations == 400
    assert result.weight_type == "PSD"
    assert result.covariance_type == "PSD"
    assert result.second_step_optimizer is not None
    assert np.linalg.eigvalsh(result.covariance).min() >= -1e-12
    np.testing.assert_allclose(result.coef, result.theta)
    np.testing.assert_allclose(result.standard_errors, result.se)


def test_numerical_and_analytic_jacobians_agree():
    data = example_data()
    common = dict(
        data=data,
        score=iv_score,
        theta_start=[0.0],
        cluster_cols=["row", "column"],
        nuisance={},
    )
    numeric = fit_gmm(**common)
    analytic = fit_gmm(**common, jacobian=iv_jacobian)
    np.testing.assert_allclose(numeric.theta, analytic.theta, atol=1e-10)
    np.testing.assert_allclose(numeric.jacobian, analytic.jacobian, rtol=1e-7)
    np.testing.assert_allclose(numeric.covariance, analytic.covariance, rtol=1e-7)


def test_one_step_uses_identity_weight():
    data = example_data()
    result = fit_gmm(
        data,
        iv_score,
        [0.0],
        ["row", "column"],
        nuisance={},
        jacobian=iv_jacobian,
        n_steps=1,
    )
    np.testing.assert_allclose(result.theta, result.first_step_theta)
    np.testing.assert_allclose(result.weight_matrix, np.eye(2))
    assert result.second_step_optimizer is None


def test_psd_weight_can_be_combined_with_cgm_covariance():
    data = example_data()
    result = fit_gmm(
        data,
        iv_score,
        [0.0],
        ["row", "column"],
        nuisance={},
        jacobian=iv_jacobian,
        weight_type="PSD",
        covariance_type="CGM",
    )
    assert result.weight_type == "PSD"
    assert result.covariance_type == "CGM"
    assert result.middle_matrix.shape == (2, 2)
    assert np.isfinite(result.covariance).all()


def test_exactly_one_nuisance_input_is_required():
    data = example_data()
    with pytest.raises(ValueError, match="exactly one"):
        fit_gmm(data, iv_score, [0.0], ["row", "column"])
    with pytest.raises(ValueError, match="exactly one"):
        fit_gmm(
            data,
            iv_score,
            [0.0],
            ["row", "column"],
            nuisance={},
            nuisance_fit=lambda frame: {},
        )


def test_score_and_parameter_validation():
    data = example_data()
    with pytest.raises(ValueError, match="at least the number"):
        fit_gmm(
            data,
            lambda frame, theta, nuisance: np.ones((len(frame), 1)),
            [0.0, 0.0],
            ["row"],
            nuisance={},
        )
    with pytest.raises(ValueError, match="two-dimensional"):
        fit_gmm(
            data,
            lambda frame, theta, nuisance: np.ones(len(frame)),
            [0.0],
            ["row"],
            nuisance={},
        )
    with pytest.raises(ValueError, match="non-finite"):
        fit_gmm(
            data,
            lambda frame, theta, nuisance: np.full((len(frame), 1), np.nan),
            [0.0],
            ["row"],
            nuisance={},
        )
    with pytest.raises(ValueError, match="n_steps"):
        fit_gmm(
            data,
            iv_score,
            [0.0],
            ["row"],
            nuisance={},
            n_steps=3,
        )


def test_singular_psd_weight_is_rejected():
    data = example_data()

    def duplicated_score(frame, theta, nuisance):
        del nuisance
        residual = frame["y"].to_numpy() - theta[0]
        return np.column_stack((residual, residual))

    with pytest.raises(RuntimeError, match="PSD moment matrix"):
        fit_gmm(
            data,
            duplicated_score,
            [0.0],
            ["row", "column"],
            nuisance={},
        )


def test_indefinite_cgm_weight_is_rejected_with_psd_recommendation():
    data = pd.DataFrame({"row": [0, 0, 1, 1], "column": [0, 1, 0, 1]})
    fixed_scores = np.array(
        [
            [0.345584, 0.821618],
            [0.330437, -1.303157],
            [0.905356, 0.446375],
            [-0.536953, 0.581118],
        ]
    )

    with pytest.raises(RuntimeError, match="use weight_type='PSD'"):
        fit_gmm(
            data,
            lambda frame, theta, nuisance: fixed_scores,
            [0.0],
            ["row", "column"],
            nuisance={},
            weight_type="CGM",
        )


def test_optimizer_failure_is_reported():
    data = example_data()
    with pytest.raises(RuntimeError, match="optimization failed"):
        fit_gmm(
            data,
            iv_score,
            [100.0],
            ["row", "column"],
            nuisance={},
            optimizer_options={"maxiter": 0},
        )
