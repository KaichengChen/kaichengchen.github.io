test_that("R agrees with the Python-generated GMM fixture", {
  fixture <- test_path("fixtures")
  data <- utils::read.csv(file.path(fixture, "golden_gmm_data.csv"))
  expected_data <- utils::read.csv(file.path(fixture, "golden_gmm_expected.csv"))
  expected <- stats::setNames(expected_data$value, expected_data$quantity)
  score <- function(data, theta, nuisance) {
    as.matrix(data[c("z1", "z2")]) * (data$y - theta[1] * data$d)
  }
  jacobian <- function(data, theta, nuisance) {
    matrix(colMeans(as.matrix(data[c("z1", "z2")]) * data$d), ncol = 1)
  }
  result <- fit_gmm(
    data, score, 0, c("row", "col"), nuisance = list(), jacobian = jacobian
  )
  expect_equal(unname(result$first_step_theta), expected[["first_step_theta"]], tolerance = 1e-6)
  expect_equal(unname(result$theta), expected[["theta"]], tolerance = 1e-6)
  expect_equal(unname(result$se), expected[["se"]], tolerance = 1e-6)
  expect_equal(result$covariance[1, 1], expected[["covariance"]], tolerance = 1e-6)
  expect_equal(
    unname(result$moment),
    c(expected[["moment_1"]], expected[["moment_2"]]),
    tolerance = 1e-6
  )
})

test_that("two-step overidentified GMM returns estimates and metadata", {
  data <- example_gmm_data()
  calls <- 0L
  fit_nuisance <- function(frame) {
    calls <<- calls + 1L
    list(fitted_on = nrow(frame))
  }
  result <- fit_gmm(
    data, iv_score, 0, c("row", "column"),
    nuisance_fit = fit_nuisance, jacobian = iv_jacobian
  )
  expect_equal(calls, 1L)
  expect_equal(unname(result$theta[1]), 1.5, tolerance = 0.15)
  expect_gt(result$se[1], 0)
  expect_equal(dim(result$ci), c(1L, 2L))
  expect_equal(result$n_eff, 20L)
  expect_equal(unname(result$cluster_counts), c(20L, 20L))
  expect_equal(result$n_cells, 400L)
  expect_equal(result$n_observations, 400L)
  expect_s3_class(result, "fullsampleDML_fit")
  expect_equal(coef(result), result$theta)
  expect_equal(vcov(result), result$covariance)
  expect_equal(confint(result), result$ci)
  expect_output(print(result), "Full-sample debiased GMM fit")
  expect_output(print(summary(result)), "Full-sample debiased GMM summary")
})

test_that("analytic and numerical Jacobians agree", {
  data <- example_gmm_data()
  numeric <- fit_gmm(data, iv_score, 0, c("row", "column"), nuisance = list())
  analytic <- fit_gmm(
    data, iv_score, 0, c("row", "column"), nuisance = list(),
    jacobian = iv_jacobian
  )
  expect_equal(numeric$theta, analytic$theta, tolerance = 1e-7)
  expect_equal(numeric$jacobian, analytic$jacobian, tolerance = 1e-6)
  expect_equal(numeric$covariance, analytic$covariance, tolerance = 1e-6)
})

test_that("one-step and mixed covariance choices work", {
  data <- example_gmm_data()
  one <- fit_gmm(
    data, iv_score, 0, c("row", "column"), nuisance = list(),
    jacobian = iv_jacobian, n_steps = 1
  )
  expect_equal(one$theta, one$first_step_theta)
  expect_equal(one$weight_matrix, diag(2))
  expect_null(one$second_step_optimizer)
  mixed <- fit_gmm(
    data, iv_score, 0, c("row", "column"), nuisance = list(),
    jacobian = iv_jacobian, weight_type = "PSD", covariance_type = "CGM"
  )
  expect_equal(mixed$covariance_type, "CGM")
  expect_true(all(is.finite(mixed$covariance)))
})

test_that("GMM validates nuisance, score, and parameter inputs", {
  data <- example_gmm_data()
  expect_error(fit_gmm(data, iv_score, 0, "row"), "exactly one")
  expect_error(
    fit_gmm(data, iv_score, 0, "row", nuisance = list(), nuisance_fit = identity),
    "exactly one"
  )
  expect_error(
    fit_gmm(data, function(...) matrix(1, nrow(data), 1), c(0, 0), "row", nuisance = list()),
    "at least the number"
  )
  expect_error(
    fit_gmm(data, function(...) rep(1, nrow(data)), 0, "row", nuisance = list()),
    "two-dimensional"
  )
  expect_error(
    fit_gmm(data, function(...) matrix(NA_real_, nrow(data), 1), 0, "row", nuisance = list()),
    "non-finite"
  )
  expect_error(fit_gmm(data, iv_score, 0, "row", nuisance = list(), n_steps = 3), "n_steps")
})

test_that("singular and indefinite two-step weights are rejected", {
  data <- example_gmm_data()
  duplicate_score <- function(data, theta, nuisance) {
    residual <- data$y - theta[1]
    cbind(residual, residual)
  }
  expect_error(
    fit_gmm(data, duplicate_score, 0, c("row", "column"), nuisance = list()),
    "PSD moment matrix"
  )
  small <- data.frame(row = c(0, 0, 1, 1), column = c(0, 1, 0, 1))
  fixed <- matrix(c(
    0.345584, 0.821618, 0.330437, -1.303157,
    0.905356, 0.446375, -0.536953, 0.581118
  ), byrow = TRUE, ncol = 2)
  expect_error(
    fit_gmm(
      small, function(...) fixed, 0, c("row", "column"),
      nuisance = list(), weight_type = "CGM"
    ),
    "use weight_type='PSD'"
  )
})

test_that("optimizer failures are reported", {
  data <- example_gmm_data()
  expect_error(
    fit_gmm(
      data, iv_score, 100, c("row", "column"), nuisance = list(),
      optimizer_options = list(maxit = 1)
    ),
    "optimization failed"
  )
})
