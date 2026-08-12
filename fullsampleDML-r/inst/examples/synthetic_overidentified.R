# Overidentified full-sample GMM on the canonical synthetic data.
# Set CV_MODE to "multiway" for strict row-and-column-disjoint tuning.

library(fullsampleDML)

required <- c("mlr3", "mlr3tuning", "mlr3learners", "paradox", "glmnet")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install the example dependencies: ", paste(missing, collapse = ", "))
}

CV_MODE <- "oneway"
CLUSTER_COLS <- c("row_id", "col_id")
data <- utils::read.csv(system.file(
  "extdata", "synthetic_gmm.csv", package = "fullsampleDML"
))

features <- function(frame) {
  stats::model.matrix(
    ~ (x1 + x2 + x3 + x4 + x5)^2 +
      I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2) + I(x5^2) - 1,
    data = frame
  )
}

make_learner <- function() {
  mlr3tuning::auto_tuner(
    tuner = mlr3tuning::tnr("grid_search", resolution = 10L),
    learner = mlr3::lrn("regr.glmnet", alpha = 0),
    resampling = mlr3::rsmp("cv", folds = 5L),
    measure = mlr3::msr("regr.mse"),
    term_evals = 10L,
    search_space = paradox::ps(
      lambda = paradox::p_dbl(log(1e-3), log(1e3), trafo = exp)
    )
  )
}

groups <- if (CV_MODE == "multiway") {
  multiway_group_kfold(
    CLUSTER_COLS, 3L, shuffle = TRUE, random_state = 2026L
  )
} else {
  NULL
}
specs <- list(
  ell = nuisance_spec(make_learner(), "y", features, groups),
  r = nuisance_spec(make_learner(), "d", features, groups),
  m1 = nuisance_spec(make_learner(), "z1", features, groups),
  m2 = nuisance_spec(make_learner(), "z2", features, groups)
)

score <- function(data, theta, nuisance) {
  p <- nuisance$predictions
  y_residual <- data$y - p$ell
  d_residual <- data$d - p$r
  z_residual <- cbind(data$z1 - p$m1, data$z2 - p$m2)
  z_residual * (y_residual - theta[1] * d_residual)
}

jacobian <- function(data, theta, nuisance) {
  p <- nuisance$predictions
  d_residual <- data$d - p$r
  z_residual <- cbind(data$z1 - p$m1, data$z2 - p$m2)
  matrix(colMeans(z_residual * d_residual), ncol = 1L)
}

result <- fit_gmm(
  data, score, 0, CLUSTER_COLS,
  nuisance_fit = mlr3_nuisance_fit(specs), jacobian = jacobian
)

cat("Full-sample overidentified GMM (R)\n")
cat("Nuisance tuning CV:", CV_MODE, "\n")
cat("Observations:", result$n_observations, "; cells:", result$n_cells, "\n")
cat("Clusters:", paste(result$cluster_counts, collapse = ", "),
    "; n_eff:", result$n_eff, "\n")
cat("First step:", result$first_step_theta[1], "\n")
cat("Two step:", result$theta[1], "\n")
cat("PSD SE:", result$se[1], "\n")
print(confint(result))
