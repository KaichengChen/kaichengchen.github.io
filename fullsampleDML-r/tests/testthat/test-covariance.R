test_that("cluster matrices match cross-language golden values", {
  data <- data.frame(
    one = c(0, 0, 1, 1),
    two = c(0, 1, 0, 1)
  )
  scores <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), byrow = TRUE, ncol = 2)
  fixture <- utils::read.csv(test_path("fixtures", "cluster_moments.csv"))
  for (method in c("PSD", "CGM")) {
    expected_rows <- fixture[fixture$method == method, ]
    expected <- matrix(expected_rows$value, 2, 2, byrow = TRUE)
    result <- compute_cluster_moment_matrix(scores, data, c("one", "two"), method)
    expect_equal(result$matrix, expected, tolerance = 1e-12)
    expect_equal(result$n_eff, 2L)
    expect_equal(result$n_cells, 4L)
  }
})

brute_force_middle <- function(scores, data, columns, method) {
  q <- ncol(scores)
  total <- matrix(0, q, q)
  subsets <- unlist(lapply(seq_along(columns), function(size) {
    utils::combn(seq_along(columns), size, simplify = FALSE)
  }), recursive = FALSE)
  if (method == "PSD") subsets <- lapply(seq_along(columns), identity)
  for (subset in subsets) {
    sign <- if (method == "PSD" || length(subset) %% 2L) 1 else -1
    for (left in seq_len(nrow(data))) {
      for (right in seq_len(nrow(data))) {
        same <- vapply(subset, function(index) {
          data[[columns[index]]][left] == data[[columns[index]]][right]
        }, logical(1))
        if (all(same)) total <- total + sign * tcrossprod(scores[left, ], scores[right, ])
      }
    }
  }
  counts <- vapply(data[columns], function(x) length(unique(x)), integer(1))
  min(counts) / nrow(data)^2 * total
}

test_that("one- through four-way matrices match the pairwise definition", {
  clusters <- data.frame(
    one = c(0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 3),
    two = c(0, 0, 1, 0, 1, 0, 1, 2, 0, 1, 2, 2),
    three = c("a", "b", "a", "b", "c", "a", "c", "b", "a", "c", "b", "c"),
    four = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )
  set.seed(2026)
  scores <- matrix(rnorm(nrow(clusters) * 3L), ncol = 3L)
  for (k in 1:4) {
    columns <- names(clusters)[seq_len(k)]
    for (method in c("PSD", "CGM")) {
      actual <- compute_cluster_moment_matrix(scores, clusters, columns, method)$matrix
      expected <- brute_force_middle(scores, clusters, columns, method)
      expect_equal(actual, expected, tolerance = 1e-12)
    }
    psd <- compute_cluster_moment_matrix(scores, clusters, columns, "PSD")$matrix
    expect_gte(min(eigen(psd, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
  }
})

test_that("cluster matrices are invariant to order and relabeling", {
  data <- data.frame(a = rep(1:3, each = 4), b = rep(letters[1:4], 3))
  set.seed(4)
  scores <- matrix(rnorm(24), ncol = 2)
  original <- compute_cluster_moment_matrix(scores, data, c("a", "b"))$matrix
  order <- sample(seq_len(nrow(data)))
  changed <- data[order, ]
  changed$a <- c("north", "south", "east")[changed$a]
  relabeled <- compute_cluster_moment_matrix(scores[order, ], changed, c("a", "b"))$matrix
  expect_equal(relabeled, original, tolerance = 1e-12)
})

test_that("invalid clustered covariance inputs are rejected", {
  data <- data.frame(a = c(1, 1, 2, 2), b = c(1, 2, 1, 2))
  scores <- matrix(1:8, ncol = 2)
  expect_error(compute_cluster_moment_matrix(scores, data, c("a", "a")), "repeated")
  expect_error(compute_cluster_moment_matrix(scores, data, "absent"), "not found")
  expect_error(compute_cluster_moment_matrix(scores, data, "a", "psd"), "exactly")
  expect_error(compute_cluster_moment_matrix(scores[, 1], data, "a"), "two-dimensional")
  bad <- scores; bad[1, 1] <- NA_real_
  expect_error(compute_cluster_moment_matrix(bad, data, "a"), "finite")
  expect_error(compute_cluster_moment_matrix(scores, transform(data, a = NA), "a"), "missing")
  expect_error(compute_cluster_moment_matrix(scores, transform(data, a = 1), "a"), "at least two")
})
