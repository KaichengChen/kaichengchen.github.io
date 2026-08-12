mlr_test_data <- function() {
  x <- seq(-1, 1, length.out = 36)
  data.frame(
    x = x,
    x2 = x^2,
    y = 1 + x,
    d = x - x^2,
    large_group = rep(1:6, 6),
    small_group = rep(1:3, each = 12),
    small_group_tie = rep(1:3, 12)
  )
}

make_tuned_learner <- function(folds = 3L) {
  mlr3tuning::auto_tuner(
    tuner = mlr3tuning::tnr("grid_search", resolution = 3L),
    learner = mlr3::lrn("regr.rpart"),
    resampling = mlr3::rsmp("cv", folds = folds),
    measure = mlr3::msr("regr.mse"),
    term_evals = 3L,
    search_space = paradox::ps(cp = paradox::p_dbl(0.001, 0.1))
  )
}

test_that("mlr3 components fit and predict on the complete sample", {
  skip_if_not_installed("mlr3")
  data <- mlr_test_data()
  specs <- list(
    ell = nuisance_spec(mlr3::lrn("regr.featureless"), "y", c("x", "x2")),
    m = nuisance_spec(mlr3::lrn("regr.featureless"), "d", function(x) x["x"])
  )
  result <- fit_mlr3_nuisances(data, specs)
  expect_named(result$models, c("ell", "m"))
  expect_named(result$predictions, c("ell", "m"))
  expect_length(result$predictions$ell, nrow(data))
  expect_true(all(is.finite(result$predictions$m)))
  expect_equal(result$group_columns$ell, character())
})

test_that("tuned components select the smallest cluster with stable ties", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3tuning")
  skip_if_not_installed("paradox")
  skip_if_not_installed("mlr3learners")
  data <- mlr_test_data()
  spec <- nuisance_spec(make_tuned_learner(), "y", c("x", "x2"))
  result <- suppressWarnings(fit_mlr3_nuisances(
    data, list(ell = spec), cluster_cols = c("large_group", "small_group")
  ))
  expect_equal(result$group_columns$ell, "small_group")
  expect_equal(unname(result$group_counts$ell), 3L)
  tied <- fit_mlr3_nuisances(
    data, list(ell = spec), cluster_cols = c("small_group_tie", "small_group")
  )
  expect_equal(tied$group_columns$ell, "small_group_tie")
})

test_that("explicit one-way and strict multiway tuning groups work", {
  skip_if_not_installed("mlr3")
  skip_if_not_installed("mlr3tuning")
  skip_if_not_installed("paradox")
  skip_if_not_installed("mlr3learners")
  data <- mlr_test_data()
  one <- nuisance_spec(
    make_tuned_learner(3), "y", c("x", "x2"), groups = "large_group"
  )
  one_result <- fit_mlr3_nuisances(
    data, list(ell = one), cluster_cols = c("large_group", "small_group")
  )
  expect_equal(one_result$group_columns$ell, "large_group")
  strict <- nuisance_spec(
    make_tuned_learner(2), "y", c("x", "x2"),
    groups = multiway_group_kfold(
      c("large_group", "small_group"), 2, shuffle = TRUE, random_state = 7
    )
  )
  strict_result <- fit_mlr3_nuisances(data, list(ell = strict))
  expect_equal(
    strict_result$group_columns$ell, c("large_group", "small_group")
  )
  expect_true(all(is.finite(strict_result$predictions$ell)))
})

test_that("the mlr3 callback receives fit_gmm clusters", {
  skip_if_not_installed("mlr3")
  data <- mlr_test_data()
  callback <- mlr3_nuisance_fit(list(
    ell = nuisance_spec(mlr3::lrn("regr.featureless"), "y", "x")
  ))
  result <- fit_gmm(
    data,
    score = function(data, theta, nuisance) {
      matrix(data$y - theta[1], ncol = 1)
    },
    theta_start = 0,
    cluster_cols = c("large_group", "small_group"),
    nuisance_fit = callback,
    jacobian = function(...) matrix(1, 1, 1),
    n_steps = 1
  )
  expect_length(result$nuisance$predictions$ell, nrow(data))
})

test_that("invalid nuisance specifications and group columns fail clearly", {
  skip_if_not_installed("mlr3")
  data <- mlr_test_data()
  expect_error(fit_mlr3_nuisances(data, list()), "at least one")
  expect_error(
    nuisance_spec(mlr3::lrn("regr.featureless"), "y", "x", groups = 1),
    "groups"
  )
  bad <- nuisance_spec(mlr3::lrn("regr.featureless"), "absent", "x")
  expect_error(fit_mlr3_nuisances(data, list(ell = bad)), "not found")
  missing <- data; missing$small_group[1] <- NA
  expect_error(
    fit_mlr3_nuisances(
      missing,
      list(ell = nuisance_spec(mlr3::lrn("regr.featureless"), "y", "x")),
      cluster_cols = "small_group"
    ),
    "missing"
  )
})
