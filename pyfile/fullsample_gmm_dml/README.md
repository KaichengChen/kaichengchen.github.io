# `fullsample-gmm-dml`

`fullsample-gmm-dml` implements full-sample debiased GMM with user-supplied
moments and arbitrary multiway clustered inference. It accompanies
*Cross-Fitting-Free Debiased Machine Learning with Multiway Dependence*,
[arXiv:2602.11333](https://arxiv.org/abs/2602.11333).

The package estimates nuisance once on the full sample (allowing cluster-aware
cross-validation tuning), accepts any finite unit-level score matrix, performs
identity-weighted one-step or feasible two-step GMM, and calculates `PSD` or
`CGM` inference.

## Installation

Python 3.10 or newer is required. There are two installation options.

### Option 1: install directly from GitHub

This is the simplest option for users who only need to import and use the
package. It requires Git to be installed:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install \
  "fullsample-gmm-dml @ git+https://github.com/kaichengchen/kaichengchen.github.io.git#subdirectory=pyfile/fullsample_gmm_dml"
```

This installs the package from the `pyfile/fullsample_gmm_dml` subdirectory
of the repository. It does not create a local copy of the examples, data, or
tests. 

### Option 2: clone the repository

Use this option to obtain the examples, synthetic data, and tests as well as
the package source:

```sh
git clone https://github.com/kaichengchen/kaichengchen.github.io.git
cd kaichengchen.github.io/pyfile/fullsample_gmm_dml
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
```

The example can then be run from the package folder:

```sh
python examples/synthetic_overidentified.py
```

On Windows, use `.venv\Scripts\activate` instead of
`source .venv/bin/activate`. The package's runtime dependencies are NumPy,
pandas, SciPy, and scikit-learn.

The downloaded folder and Python import package are named
`fullsample_gmm_dml`. The distribution installed by pip is named
`fullsample-gmm-dml`.:

```python
from fullsample_gmm_dml import fit_gmm
```

For development after cloning, install the source in editable mode together
with the test dependency and run the tests:

```sh
python -m pip install -e ".[test]"
python -m pytest
```

## Generic GMM interface

```python
from fullsample_gmm_dml import fit_gmm

result = fit_gmm(
    data=data,
    score=score,
    theta_start=[0.0],
    cluster_cols=["firm", "market", "year"],
    nuisance_fit=fit_nuisance,
    weight_type="PSD",
    covariance_type="PSD",
)
```

The score has signature

```python
def score(data, theta, nuisance):
    # Return one row per observation and one column per moment.
    return unit_scores
```

For $p$ parameters it may return any $q \ge p$ moments. Supply either
`nuisance_fit`, which is called exactly once with the complete data, or a
precomputed `nuisance` object. A nuisance-free model can pass `nuisance={}`.
There is no sample-splitting or cross-fitting stage.

By default, `fit_gmm` first minimizes the identity-weighted criterion. It then
inverts the preliminary clustered moment matrix and minimizes the feasible
second-step criterion. Set `n_steps=1` to retain the identity-weighted
estimate. An optional analytic Jacobian must return the $q \times p$ matrix

$$
J(\theta)=-\partial_\theta E[\psi(X,\theta,\eta)].
$$

Otherwise, the package uses central numerical derivatives. The result stores
the coefficient, standard error, 95% confidence interval, covariance,
Jacobian, GMM weight, first-step estimate, sample metadata, nuisance object,
and optimizer diagnostics.

## Full-sample nuisance fitting and tuning

The core accepts any nuisance callback, so learners from XGBoost, LightGBM,
PyTorch, statsmodels, or custom code can be used. The convenience adapter
fits scikit-learn-compatible learners:

```python
from fullsample_gmm_dml import NuisanceSpec, sklearn_nuisance_fit
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Lasso
from sklearn.model_selection import GridSearchCV, GroupKFold, RandomizedSearchCV

specs = {
    "ell": NuisanceSpec(
        estimator=GridSearchCV(
            Lasso(max_iter=20_000),
            {"alpha": [0.001, 0.01, 0.1]},
            cv=GroupKFold(n_splits=5),
            refit=True,
        ),
        target="y",
        features=x_cols,
    ),
    "m": NuisanceSpec(
        estimator=RandomizedSearchCV(
            RandomForestRegressor(random_state=1),
            {"max_depth": [3, 5, 10, None]},
            n_iter=4,
            cv=GroupKFold(n_splits=4),
            refit=True,
        ),
        target="d",
        features=x_cols,
    ),
}
nuisance_fit = sklearn_nuisance_fit(specs)
```

Every nuisance component may use its own estimator, target, features, scoring
rule, search space, CV splitter, prediction method, and fit parameters. The
tuning groups default to the clustering dimension with the fewest distinct
clusters among the `cluster_cols` passed to `fit_gmm`. Ties are resolved by
the order of `cluster_cols`. It is a practical default based on the
effective sample size, not a tuning rule uniquely required by the paper.

Each `GroupKFold` search above automatically receives that selected group
vector. Set `NuisanceSpec(groups="firm")` to override the automatic choice
for one learner. Direct calls to `fit_sklearn_nuisances` can provide the same
candidate dimensions through `cluster_cols`. 

Use `MultiwayGroupKFold` when tuning training and validation samples must not
share identifiers in any clustering dimension:

```python
from fullsample_gmm_dml import MultiwayGroupKFold

search = GridSearchCV(
    estimator,
    parameter_grid,
    cv=MultiwayGroupKFold(
        n_splits_per_dimension=3,
        shuffle=True,
        random_state=42,
    ),
    refit=True,
)
spec = NuisanceSpec(
    estimator=search,
    target="y",
    features=x_cols,
    groups=["firm", "market", "year"],
)
```

The splitter assigns the clusters in each dimension to fold bins. Validation
uses the intersection of one held-out bin from every dimension. Training uses
only observations outside every held-out bin. This guarantees no
train--validation cluster overlap in any requested dimension. With $F$ bins
in each of $K$ dimensions, it evaluates $F^K$ splits, so it is substantially
more expensive than the one-way default.

For example, suppose row and column clusters are each assigned to three bins:

```
             Column bin
              C1   C2   C3
Row bin R1    11   12   13
        R2    21   22   23
        R3    31   32   33
```

For the split where block 11 is validation:

```
              C1      C2      C3
R1         Validate   Omit    Omit
R2           Omit    Train   Train
R3           Omit    Train   Train
```
The next split could validate on block ```12```. It continues until each of the 
nine blocks has served as validation exactly once.


The tuning folds choose hyperparameters; with `refit=True`, each selected learner
is subsequently fitted on all observations and predicts those same
observations. Search estimators with `refit=False` are rejected (refit belongs 
to scikit-learn’s GridSearchCV, not to fullsample_gmm_dml).


## Arbitrary multiway inference

The implementation first sums unit scores inside each observed full
$K$-way cell. For `PSD`, it sums the outer products of the cluster totals
separately over every dimension and multiplies the result by
$n_{\mathrm{eff}}/\bar M^2$. This is the default variance estimator in
[arXiv:2602.11333](https://arxiv.org/abs/2602.11333). It is positive
semidefinite and requires approximately $K$ aggregation passes.

The `CGM`-type variance estimator, with the double-counting term removed, 
is also supported in the package. `CGM` can be indefinite and
becomes more expensive as $K$ grows. It can be used for the reported
covariance independently of the GMM weight:

```python
result = fit_gmm(
    ...,
    weight_type="PSD",
    covariance_type="CGM",
)
```

Using `weight_type="CGM"` requires its preliminary matrix to be positive
definite. Otherwise, estimation stops with an error recommending `PSD` as the
weight. Heterogeneous cell sizes, missing cells, string cluster identifiers,
and nonrectangular arrays are supported. Cluster identifiers themselves must
be observed, and every dimension must contain at least two clusters.

## Synthetic overidentified example

The worked example uses heterogeneous-cell synthetic data ```synthetic_gmm.csv```
generated by the script ```generate_synthetic.py```. There are two
instruments, one target parameter, and four independently tuned nuisance
regressions:

```sh
python examples/synthetic_overidentified.py
```

At the top of the script, change one setting to choose the tuning scheme:

```python
CV_MODE = "oneway"   # automatic smallest clustering dimension
CV_MODE = "multiway" # strict row-and-column-disjoint tuning
```

The multiway setting uses three fold bins per clustering dimension, producing
$3^2=9$ tuning splits. For each split, one block is validation, the four
blocks sharing neither held-out dimension are training, and the four border
blocks are omitted.

Its two-component orthogonal score is

$$
\psi_i(\theta)=
\begin{pmatrix}\tilde Z_{1i}\\\tilde Z_{2i}\end{pmatrix}
(\tilde Y_i-\theta\tilde D_i),
$$

so the scalar coefficient is overidentified. The example prints the
identity-weighted estimate, feasible two-step estimate, `PSD` standard error,
confidence interval, automatically selected tuning group, and selected tuning
parameters. In the supplied data, the automatic rule selects `col_id` because
it has 35 clusters while `row_id` has 40.

## Tests

Run the package tests with

```sh
python -m pip install -e ".[test]"
python -m pytest
```
