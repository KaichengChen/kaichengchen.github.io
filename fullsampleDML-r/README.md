# fullsampleDML for R

`fullsampleDML` implements full-sample debiased GMM with user-supplied moments
and arbitrary multiway clustered inference. It accompanies *Cross-Fitting-Free
Debiased Machine Learning with Multiway Dependence*
([arXiv:2602.11333](https://arxiv.org/abs/2602.11333)).

The package is native R: it does not require Python or `reticulate`. Its core
uses base/recommended R only. An optional `mlr3` adapter independently tunes
nuisance learners with cluster-aware folds and then refits each selected
learner on every observation. Those tuning folds are not DML cross-fitting.

## Installation

Install the current development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("kaichengchen/fullsampleDML-r")
```

After the first release, pin version 0.1.0 with:

```r
remotes::install_github("kaichengchen/fullsampleDML-r@v0.1.0")
```

## Generic GMM interface

```r
library(fullsampleDML)

result <- fit_gmm(
  data = data,
  score = score,
  theta_start = 0,
  cluster_cols = c("firm", "market", "year"),
  nuisance_fit = fit_nuisance,
  weight_type = "PSD",
  covariance_type = "PSD"
)

coef(result)
vcov(result)
confint(result)
summary(result)
```

The score callback has the signature below and must return one row per
observation and one column per moment:

```r
score <- function(data, theta, nuisance) {
  unit_scores
}
```

For `p` parameters the score may provide any `q >= p` moments. Supply exactly
one of `nuisance_fit`, which is called once on the complete data, or a
precomputed `nuisance` object. A nuisance-free model uses `nuisance = list()`.

The default is feasible two-step GMM: an identity-weighted first step followed
by a second step using the inverse clustered moment matrix. Set `n_steps = 1`
for the identity-weighted estimate. An optional analytic Jacobian returns the
`q`-by-`p` negative derivative of the average score; otherwise central
numerical derivatives are used.

## Clustered inference

`compute_cluster_moment_matrix()` first sums scores inside observed full cells.
`method = "PSD"` adds the outer products of cluster totals separately across
all dimensions. `method = "CGM"` applies inclusion-exclusion over all nonempty
subsets. Both use `n_eff = min(cluster_counts)` and support heterogeneous or
missing cells and nonnumeric cluster identifiers.

`PSD` is positive semidefinite and is the recommended GMM weight. A CGM weight
is rejected when its preliminary matrix is not positive definite. The weight
and reported covariance choices can differ.

## Optional mlr3 nuisance fitting

```r
library(mlr3)
library(mlr3tuning)
library(paradox)

learner <- auto_tuner(
  tuner = tnr("grid_search", resolution = 10),
  learner = lrn("regr.rpart"),
  resampling = rsmp("cv", folds = 5),
  measure = msr("regr.mse"),
  term_evals = 10,
  search_space = ps(cp = p_dbl(0.001, 0.1))
)

specs <- list(
  ell = nuisance_spec(learner, target = "y", features = x_cols),
  m = nuisance_spec(learner, target = "d", features = x_cols)
)

result <- fit_gmm(
  data, score, theta_start = 0,
  cluster_cols = c("firm", "year"),
  nuisance_fit = mlr3_nuisance_fit(specs)
)
```

Unspecified tuning groups default to the clustering dimension with the fewest
distinct clusters; ties follow `cluster_cols` order. Request strict multiway
tuning with:

```r
strict_groups <- multiway_group_kfold(
  group_cols = c("firm", "year"),
  n_splits_per_dimension = 3,
  shuffle = TRUE,
  random_state = 42
)

spec <- nuisance_spec(
  learner, target = "y", features = x_cols, groups = strict_groups
)
```

Strict validation sets are intersections of held-out cluster bins. Training
uses observations outside the held-out bin in every dimension; border regions
are omitted. With `F` bins in each of `K` dimensions this produces `F^K`
tuning splits.

## Reproducible example and tests

The canonical 3,521-observation synthetic data are installed at:

```r
system.file("extdata", "synthetic_gmm.csv", package = "fullsampleDML")
```

See `vignette("synthetic-overidentified", package = "fullsampleDML")` for the
worked overidentified example. Development tests run with:

```r
testthat::test_local()
```

Users remain responsible for verifying the nuisance-rate, complexity,
orthogonality, and dependence conditions in the paper.
