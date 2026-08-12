#' Fit full-sample debiased GMM with multiway clustered inference
#'
#' `score` must return an n-by-q numeric matrix. The optional analytic
#' `jacobian` must return the q-by-p negative derivative of the average score.
#' Supply exactly one of a precomputed `nuisance` object or a `nuisance_fit`
#' callback. Use `nuisance = list()` for nuisance-free models.
#'
#' @param data A data frame with one row per unit observation.
#' @param score Function with arguments `(data, theta, nuisance)`.
#' @param theta_start Finite starting parameter vector.
#' @param cluster_cols Character vector of clustering variables.
#' @param nuisance Precomputed nuisance object.
#' @param nuisance_fit Callback fitted exactly once on the complete data.
#' @param jacobian Optional analytic negative Jacobian callback.
#' @param weight_type `"PSD"` or `"CGM"` for the feasible GMM weight.
#' @param covariance_type `"PSD"` or `"CGM"` for reported inference.
#' @param optimizer_method Method passed to [stats::optim()].
#' @param optimizer_options Control list passed to [stats::optim()].
#' @param n_steps Either one or two GMM steps.
#' @return A `fullsampleDML_fit` object.
#' @export
fit_gmm <- function(data, score, theta_start, cluster_cols,
                    nuisance = NULL, nuisance_fit = NULL, jacobian = NULL,
                    weight_type = "PSD", covariance_type = "PSD",
                    optimizer_method = "BFGS", optimizer_options = NULL,
                    n_steps = 2L) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame", call. = FALSE)
  }
  if (!nrow(data)) {
    stop("data must contain at least one observation", call. = FALSE)
  }
  if (!is.function(score)) {
    stop("score must be a function", call. = FALSE)
  }
  if (is.null(nuisance) == is.null(nuisance_fit)) {
    stop("provide exactly one of nuisance or nuisance_fit", call. = FALSE)
  }
  if (!is.null(nuisance_fit) && !is.function(nuisance_fit)) {
    stop("nuisance_fit must be a function", call. = FALSE)
  }
  if (!is.null(jacobian) && !is.function(jacobian)) {
    stop("jacobian must be NULL or a function", call. = FALSE)
  }
  selected_weight <- .validate_variance_type(weight_type, "weight_type")
  selected_covariance <- .validate_variance_type(
    covariance_type, "covariance_type"
  )
  if (!is.numeric(n_steps) || length(n_steps) != 1L ||
      is.na(n_steps) || !n_steps %in% c(1, 2)) {
    stop("n_steps must be 1 or 2", call. = FALSE)
  }
  n_steps <- as.integer(n_steps)
  start <- as.numeric(theta_start)
  if (!length(start) || is.null(theta_start) ||
      is.matrix(theta_start) || any(!is.finite(start))) {
    stop("theta_start must be a finite, nonempty one-dimensional vector", call. = FALSE)
  }
  if (!is.character(optimizer_method) || length(optimizer_method) != 1L) {
    stop("optimizer_method must be one character value", call. = FALSE)
  }
  if (!is.null(optimizer_options) && !is.list(optimizer_options)) {
    stop("optimizer_options must be NULL or a list", call. = FALSE)
  }

  fitted_nuisance <- if (!is.null(nuisance_fit)) {
    if (inherits(nuisance_fit, "fullsampleDML_mlr3_nuisance_fitter")) {
      nuisance_fit(data, cluster_cols = cluster_cols)
    } else {
      nuisance_fit(data)
    }
  } else {
    nuisance
  }

  initial_scores <- .evaluate_score(score, data, start, fitted_nuisance)
  q <- ncol(initial_scores)
  p <- length(start)
  if (q < p) {
    stop(
      "the number of moments must be at least the number of parameters",
      call. = FALSE
    )
  }
  identity_weight <- diag(q)
  first <- .optimize_gmm(
    score, data, fitted_nuisance, start, identity_weight, q,
    optimizer_method, optimizer_options
  )
  first_scores <- .evaluate_score(
    score, data, first$theta, fitted_nuisance, expected_q = q
  )
  preliminary <- compute_cluster_moment_matrix(
    first_scores, data, cluster_cols, selected_weight
  )

  if (n_steps == 2L) {
    final_weight <- .positive_definite_inverse(
      preliminary$matrix, selected_weight
    )
    second <- .optimize_gmm(
      score, data, fitted_nuisance, first$theta, final_weight, q,
      optimizer_method, optimizer_options
    )
    theta <- second$theta
    second_diagnostics <- second$diagnostics
  } else {
    final_weight <- identity_weight
    theta <- first$theta
    second_diagnostics <- NULL
  }

  final_scores <- .evaluate_score(
    score, data, theta, fitted_nuisance, expected_q = q
  )
  final_middle <- compute_cluster_moment_matrix(
    final_scores, data, cluster_cols, selected_covariance
  )
  estimated_jacobian <- .estimate_jacobian(
    jacobian, score, data, theta, fitted_nuisance, q, p
  )
  covariance <- .coefficient_covariance(
    estimated_jacobian, final_weight, final_middle$matrix, final_middle$n_eff
  )
  if (selected_covariance == "PSD") {
    eigenvalues <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
    tolerance <- .matrix_tolerance(eigenvalues, p)
    if (min(eigenvalues) < -tolerance) {
      stop("PSD produced a non-positive-semidefinite covariance", call. = FALSE)
    }
  }
  inference <- .standard_errors_and_ci(theta, covariance, selected_covariance)
  if (is.null(names(theta))) names(theta) <- paste0("theta", seq_along(theta))
  names(first$theta) <- names(theta)
  names(inference$se) <- names(theta)
  rownames(inference$ci) <- names(theta)

  result <- list(
    theta = theta,
    se = inference$se,
    ci = inference$ci,
    covariance = covariance,
    moment = colMeans(final_scores),
    jacobian = estimated_jacobian,
    middle_matrix = final_middle$matrix,
    preliminary_middle_matrix = preliminary$matrix,
    weight_matrix = final_weight,
    first_step_theta = first$theta,
    n_eff = final_middle$n_eff,
    cluster_counts = final_middle$cluster_counts,
    n_cells = final_middle$n_cells,
    n_observations = final_middle$n_observations,
    nuisance = fitted_nuisance,
    n_steps = n_steps,
    weight_type = selected_weight,
    covariance_type = selected_covariance,
    first_step_optimizer = first$diagnostics,
    second_step_optimizer = second_diagnostics,
    call = match.call()
  )
  class(result) <- "fullsampleDML_fit"
  result
}

.evaluate_score <- function(score, data, theta, nuisance, expected_q = NULL) {
  values <- tryCatch(
    score(data, theta, nuisance),
    error = function(error) {
      stop("score could not be evaluated: ", conditionMessage(error), call. = FALSE)
    }
  )
  if (!is.matrix(values) || !is.numeric(values) || length(dim(values)) != 2L) {
    stop("score must return a two-dimensional numeric n-by-q matrix", call. = FALSE)
  }
  storage.mode(values) <- "double"
  if (nrow(values) != nrow(data)) {
    stop("score must return one row per observation in data", call. = FALSE)
  }
  if (!ncol(values)) {
    stop("score must return at least one moment", call. = FALSE)
  }
  if (!is.null(expected_q) && ncol(values) != expected_q) {
    stop("score returned a different number of moments across calls", call. = FALSE)
  }
  if (any(!is.finite(values))) {
    stop("score returned non-finite values", call. = FALSE)
  }
  values
}

.optimize_gmm <- function(score, data, nuisance, start, weight, q,
                          method, control) {
  objective <- function(theta) {
    moment <- colMeans(.evaluate_score(
      score, data, as.numeric(theta), nuisance, expected_q = q
    ))
    as.numeric(crossprod(moment, weight %*% moment))
  }
  result <- tryCatch(
    stats::optim(start, objective, method = method, control = control),
    error = function(error) {
      stop(
        "GMM optimization could not be run: ", conditionMessage(error),
        call. = FALSE
      )
    }
  )
  message <- if (is.null(result$message)) "" else as.character(result$message)
  diagnostics <- list(
    success = identical(result$convergence, 0L),
    status = as.integer(result$convergence),
    message = message,
    objective = as.numeric(result$value),
    n_iterations = NA_integer_,
    n_function_evaluations = unname(as.integer(result$counts[["function"]]))
  )
  class(diagnostics) <- "fullsampleDML_optimizer_diagnostics"
  if (!diagnostics$success) {
    detail <- if (nzchar(message)) message else paste("status", diagnostics$status)
    stop("GMM optimization failed: ", detail, call. = FALSE)
  }
  theta <- as.numeric(result$par)
  if (length(theta) != length(start) || any(!is.finite(theta))) {
    stop("GMM optimization returned an invalid parameter vector", call. = FALSE)
  }
  list(theta = theta, diagnostics = diagnostics)
}

.estimate_jacobian <- function(jacobian, score, data, theta, nuisance, q, p) {
  if (is.null(jacobian)) {
    result <- matrix(NA_real_, q, p)
    for (index in seq_len(p)) {
      step <- 1e-6 * max(1, abs(theta[index]))
      upper <- lower <- theta
      upper[index] <- upper[index] + step
      lower[index] <- lower[index] - step
      upper_moment <- colMeans(.evaluate_score(
        score, data, upper, nuisance, expected_q = q
      ))
      lower_moment <- colMeans(.evaluate_score(
        score, data, lower, nuisance, expected_q = q
      ))
      result[, index] <- -(upper_moment - lower_moment) / (2 * step)
    }
  } else {
    result <- jacobian(data, theta, nuisance)
  }
  if (!is.matrix(result) || !is.numeric(result) ||
      !identical(dim(result), c(q, p))) {
    stop(
      "jacobian must return a numeric matrix with shape ", q, "-by-", p,
      call. = FALSE
    )
  }
  storage.mode(result) <- "double"
  if (any(!is.finite(result))) {
    stop("jacobian returned non-finite values", call. = FALSE)
  }
  result
}

.matrix_tolerance <- function(eigenvalues, dimension) {
  .Machine$double.eps * max(1, dimension) *
    max(1, max(abs(eigenvalues))) * 100
}

.positive_definite_inverse <- function(matrix, method) {
  symmetric <- (matrix + t(matrix)) / 2
  eigenvalues <- eigen(symmetric, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= .matrix_tolerance(eigenvalues, nrow(symmetric))) {
    if (method == "CGM") {
      stop(
        "the preliminary CGM moment matrix is not positive definite; use weight_type='PSD'",
        call. = FALSE
      )
    }
    stop(
      "the preliminary PSD moment matrix is not positive definite; a two-step weight cannot be constructed",
      call. = FALSE
    )
  }
  solve(symmetric)
}

.coefficient_covariance <- function(jacobian, weight, middle, n_eff) {
  gram <- crossprod(jacobian, weight %*% jacobian)
  bread <- tryCatch(
    solve(gram),
    error = function(error) {
      stop(
        "the weighted Jacobian is singular; the parameter is not identified",
        call. = FALSE
      )
    }
  )
  projection <- bread %*% t(jacobian) %*% weight
  covariance <- projection %*% middle %*% t(projection) / n_eff
  covariance <- (covariance + t(covariance)) / 2
  if (any(!is.finite(covariance))) {
    stop("the coefficient covariance contains non-finite values", call. = FALSE)
  }
  covariance
}

.standard_errors_and_ci <- function(theta, covariance, covariance_type) {
  variances <- diag(covariance)
  tolerance <- .Machine$double.eps * length(theta) *
    max(1, max(abs(variances))) * 100
  materially_negative <- variances < -tolerance
  if (any(materially_negative)) {
    warning(
      covariance_type,
      " produced a negative coefficient variance; the corresponding standard error and confidence interval are undefined",
      call. = FALSE
    )
  }
  variances[variances < 0 & !materially_negative] <- 0
  se <- rep(NA_real_, length(variances))
  valid <- variances >= 0
  se[valid] <- sqrt(variances[valid])
  critical <- stats::qnorm(0.975)
  ci <- cbind(lower = theta - critical * se, upper = theta + critical * se)
  list(se = se, ci = ci)
}

#' @export
coef.fullsampleDML_fit <- function(object, ...) object$theta

#' @export
vcov.fullsampleDML_fit <- function(object, ...) object$covariance

#' @export
confint.fullsampleDML_fit <- function(object, parm, level = 0.95, ...) {
  if (!is.numeric(level) || length(level) != 1L ||
      is.na(level) || level <= 0 || level >= 1) {
    stop("level must be strictly between zero and one", call. = FALSE)
  }
  estimate <- object$theta
  se <- object$se
  if (!missing(parm)) {
    estimate <- estimate[parm]
    se <- se[parm]
  }
  critical <- stats::qnorm(1 - (1 - level) / 2)
  cbind(lower = estimate - critical * se, upper = estimate + critical * se)
}

#' @export
print.fullsampleDML_fit <- function(x, ...) {
  cat("Full-sample debiased GMM fit\n")
  cat("Observations:", x$n_observations, "; cells:", x$n_cells,
      "; n_eff:", x$n_eff, "\n")
  cat("Steps:", x$n_steps, "; weight:", x$weight_type,
      "; covariance:", x$covariance_type, "\n")
  print(cbind(estimate = x$theta, std.error = x$se), ...)
  invisible(x)
}

#' @export
summary.fullsampleDML_fit <- function(object, ...) {
  result <- object
  result$coefficients <- cbind(
    estimate = object$theta,
    std.error = object$se,
    conf.low = object$ci[, 1L],
    conf.high = object$ci[, 2L]
  )
  class(result) <- c("summary.fullsampleDML_fit", class(object))
  result
}

#' @export
print.summary.fullsampleDML_fit <- function(x, ...) {
  cat("Full-sample debiased GMM summary\n")
  cat("Observations:", x$n_observations, "; cells:", x$n_cells,
      "; n_eff:", x$n_eff, "\n")
  cat("Clusters:",
      paste(names(x$cluster_counts), x$cluster_counts, sep = "=", collapse = ", "),
      "\n")
  stats::printCoefmat(x$coefficients, ...)
  invisible(x)
}
