complete_groups <- function(k, n_clusters = 6L) {
  grid <- do.call(expand.grid, rep(list(seq_len(n_clusters)), k))
  as.data.frame(grid)
}

test_that("strict folds are disjoint and cover validation exactly once", {
  for (k in 1:3) {
    groups <- complete_groups(k)
    result <- multiway_group_splits(groups, 3)
    expect_equal(result$n_splits, 3^k)
    counts <- integer(nrow(groups))
    for (index in seq_len(result$n_splits)) {
      train <- result$train[[index]]
      test <- result$test[[index]]
      counts[test] <- counts[test] + 1L
      expect_length(intersect(train, test), 0L)
      for (dimension in seq_len(k)) {
        expect_length(
          intersect(groups[[dimension]][train], groups[[dimension]][test]), 0L
        )
      }
    }
    expect_equal(counts, rep(1L, nrow(groups)))
  }
})

test_that("fold counts, sparse cells, and strings are supported", {
  groups <- complete_groups(2)
  expect_equal(multiway_group_splits(groups, c(2, 3))$n_splits, 6L)
  sparse <- groups[!(groups[[1]] == 1 & groups[[2]] == 1), ]
  heterogeneous <- rbind(sparse, sparse[1:5, ], sparse[1:2, ])
  result <- multiway_group_splits(heterogeneous, 2)
  expect_equal(result$n_splits, 4L)
  strings <- as.data.frame(lapply(groups, as.character))
  expect_equal(multiway_group_splits(strings, 2)$n_splits, 4L)
})

test_that("shuffled folds are reproducible without changing global RNG", {
  groups <- complete_groups(2)
  set.seed(99)
  before <- .Random.seed
  first <- multiway_group_splits(groups, 3, shuffle = TRUE, random_state = 19)
  after <- .Random.seed
  second <- multiway_group_splits(groups, 3, shuffle = TRUE, random_state = 19)
  expect_equal(first$train, second$train)
  expect_equal(first$test, second$test)
  expect_equal(after, before)
})

test_that("invalid and empty multiway splits are rejected", {
  groups <- complete_groups(2, 3)
  expect_error(multiway_group_splits(NULL, 2), "required")
  invalid <- groups; invalid[1, 1] <- NA
  expect_error(multiway_group_splits(invalid, 2), "missing")
  expect_error(multiway_group_splits(groups, 4), "requires at least 4")
  expect_error(multiway_group_splits(groups, c(2, 2, 2)), "one value per")
  diagonal <- data.frame(a = 1:2, b = 1:2)
  expect_error(multiway_group_splits(diagonal, 2), "validation intersection is empty")
  missing_corner <- data.frame(a = c(1, 1, 2), b = c(1, 2, 1))
  expect_error(multiway_group_splits(missing_corner, 2), "training set is empty")
})

test_that("multiway specifications retain their configuration", {
  spec <- multiway_group_kfold(
    c("firm", "year"), c(2, 3), shuffle = TRUE, random_state = 42
  )
  expect_s3_class(spec, "fullsampleDML_multiway_kfold")
  expect_equal(spec$group_cols, c("firm", "year"))
  expect_equal(spec$n_splits_per_dimension, c(2, 3))
})
