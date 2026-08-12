example_gmm_data <- function() {
  set.seed(881)
  n_side <- 20L
  row <- rep(seq_len(n_side), each = n_side)
  column <- rep(seq_len(n_side), times = n_side)
  z1 <- rnorm(n_side^2)
  z2 <- 0.25 * z1 + rnorm(n_side^2)
  d <- 0.8 * z1 + 0.5 * z2 + rnorm(n_side^2)
  error <- rnorm(n_side, sd = 0.4)[row] +
    rnorm(n_side, sd = 0.4)[column] + rnorm(n_side^2)
  y <- 1.5 * d + error
  data.frame(row, column, z1, z2, d, y)
}

iv_score <- function(data, theta, nuisance) {
  instruments <- as.matrix(data[c("z1", "z2")])
  residual <- data$y - theta[1L] * data$d
  instruments * residual
}

iv_jacobian <- function(data, theta, nuisance) {
  matrix(colMeans(as.matrix(data[c("z1", "z2")]) * data$d), ncol = 1L)
}
