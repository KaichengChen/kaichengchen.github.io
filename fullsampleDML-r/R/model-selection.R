#' Specify strict multiway grouped cross-validation
#'
#' @param group_cols Optional clustering columns to attach to the specification.
#' @param n_splits_per_dimension One integer or one integer per group dimension.
#' @param shuffle Whether to randomize the cluster-to-bin assignment.
#' @param random_state Optional integer seed.
#' @return A `fullsampleDML_multiway_kfold` specification.
#' @export
multiway_group_kfold <- function(group_cols = NULL,
                                 n_splits_per_dimension = 3L,
                                 shuffle = FALSE,
                                 random_state = NULL) {
  if (!is.null(group_cols) && (!is.character(group_cols) || !length(group_cols))) {
    stop("group_cols must be NULL or a nonempty character vector", call. = FALSE)
  }
  if (!is.logical(shuffle) || length(shuffle) != 1L || is.na(shuffle)) {
    stop("shuffle must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(random_state) &&
      (!is.numeric(random_state) || length(random_state) != 1L ||
       is.na(random_state))) {
    stop("random_state must be NULL or one finite number", call. = FALSE)
  }
  result <- list(
    group_cols = group_cols,
    n_splits_per_dimension = n_splits_per_dimension,
    shuffle = shuffle,
    random_state = if (is.null(random_state)) NULL else as.integer(random_state)
  )
  class(result) <- "fullsampleDML_multiway_kfold"
  result
}

#' Materialize strict multiway grouped folds
#'
#' Validation observations lie in one intersection of held-out cluster bins.
#' Training observations must be outside the held-out bin in every dimension;
#' border observations are omitted from that split.
#'
#' @param groups A vector, matrix, or data frame of group identifiers.
#' @param n_splits_per_dimension One integer or one integer per dimension.
#' @param shuffle Whether to randomize cluster-to-bin assignment.
#' @param random_state Optional integer seed.
#' @return A list containing one-based `train` and `test` index lists.
#' @export
multiway_group_splits <- function(groups, n_splits_per_dimension = 3L,
                                  shuffle = FALSE, random_state = NULL) {
  groups <- .validate_groups(groups)
  k <- ncol(groups)
  fold_counts <- .fold_counts(n_splits_per_dimension, k)
  if (!is.logical(shuffle) || length(shuffle) != 1L || is.na(shuffle)) {
    stop("shuffle must be TRUE or FALSE", call. = FALSE)
  }

  restore_seed <- FALSE
  if (shuffle && !is.null(random_state)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(random_state))
    restore_seed <- TRUE
  }
  invisible(restore_seed)

  bins <- matrix(NA_integer_, nrow(groups), k)
  for (dimension in seq_len(k)) {
    labels <- unique(groups[[dimension]])
    n_folds <- fold_counts[dimension]
    if (length(labels) < n_folds) {
      stop(
        "group dimension ", dimension, " has ", length(labels),
        " distinct clusters but requires at least ", n_folds, call. = FALSE
      )
    }
    order <- seq_along(labels)
    if (shuffle) order <- sample(order)
    sizes <- rep(length(labels) %/% n_folds, n_folds)
    if (length(labels) %% n_folds) {
      sizes[seq_len(length(labels) %% n_folds)] <-
        sizes[seq_len(length(labels) %% n_folds)] + 1L
    }
    cluster_bins <- integer(length(labels))
    cursor <- 1L
    for (fold in seq_len(n_folds)) {
      members <- order[cursor:(cursor + sizes[fold] - 1L)]
      cluster_bins[members] <- fold
      cursor <- cursor + sizes[fold]
    }
    bins[, dimension] <- cluster_bins[match(groups[[dimension]], labels)]
  }

  values <- lapply(fold_counts, seq_len)
  grid <- do.call(expand.grid, c(rev(values), stringsAsFactors = FALSE))
  grid <- grid[, rev(seq_len(k)), drop = FALSE]
  train <- vector("list", nrow(grid))
  test <- vector("list", nrow(grid))
  for (split_index in seq_len(nrow(grid))) {
    held_out <- as.integer(grid[split_index, ])
    validation <- which(rowSums(sweep(bins, 2L, held_out, `==`)) == k)
    if (!length(validation)) {
      stop(
        "a multiway validation intersection is empty; reduce the number of folds per dimension",
        call. = FALSE
      )
    }
    training <- which(rowSums(sweep(bins, 2L, held_out, `!=`)) == k)
    if (!length(training)) {
      stop(
        "a multiway training set is empty; adjust the folds or clusters",
        call. = FALSE
      )
    }
    train[[split_index]] <- training
    test[[split_index]] <- validation
  }
  result <- list(
    train = train,
    test = test,
    n_splits = length(train),
    fold_counts = fold_counts,
    bins = bins
  )
  class(result) <- "fullsampleDML_multiway_splits"
  result
}

.validate_groups <- function(groups) {
  if (is.null(groups)) {
    stop("groups are required", call. = FALSE)
  }
  if (is.atomic(groups) && is.null(dim(groups))) {
    groups <- data.frame(group = groups, check.names = FALSE)
  } else {
    groups <- as.data.frame(groups, stringsAsFactors = FALSE)
  }
  if (!nrow(groups) || !ncol(groups)) {
    stop("groups must be an n-by-K object with K >= 1", call. = FALSE)
  }
  if (anyNA(groups)) {
    stop("groups must not contain missing values", call. = FALSE)
  }
  groups
}

.fold_counts <- function(requested, k) {
  if (!is.numeric(requested) || anyNA(requested) ||
      any(!is.finite(requested)) || any(requested != as.integer(requested))) {
    stop("n_splits_per_dimension must contain integers", call. = FALSE)
  }
  if (length(requested) == 1L) {
    counts <- rep(as.integer(requested), k)
  } else {
    counts <- as.integer(requested)
    if (length(counts) != k) {
      stop(
        "n_splits_per_dimension must have one value per group dimension",
        call. = FALSE
      )
    }
  }
  if (any(counts < 2L)) {
    stop("every dimension must use at least two folds", call. = FALSE)
  }
  counts
}
