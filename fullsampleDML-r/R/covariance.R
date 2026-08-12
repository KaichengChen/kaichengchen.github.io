#' Compute an arbitrary-way clustered score middle matrix
#'
#' Unit scores are first aggregated within each observed full cluster cell.
#' The result is normalized by `n_eff / nrow(data)^2`, where `n_eff` is the
#' smallest number of clusters across the requested dimensions.
#'
#' @param scores A finite numeric matrix with one row per observation.
#' @param data A data frame containing the clustering variables.
#' @param cluster_cols A nonempty character vector of clustering columns.
#' @param method Either `"PSD"` or `"CGM"`.
#' @return An object of class `fullsampleDML_cluster_moment`.
#' @export
compute_cluster_moment_matrix <- function(scores, data, cluster_cols,
                                          method = "PSD") {
  method <- .validate_variance_type(method, "method")
  validated <- .validate_cluster_inputs(scores, data, cluster_cols)
  scores <- validated$scores
  cluster_cols <- validated$cluster_cols
  counts <- validated$counts

  q <- ncol(scores)
  unscaled <- matrix(0, nrow = q, ncol = q)
  k <- length(cluster_cols)
  if (identical(method, "PSD")) {
    for (dimension in seq_len(k)) {
      unscaled <- unscaled + .subset_outer_sum(
        scores, data, cluster_cols[dimension]
      )
    }
  } else {
    for (subset_size in seq_len(k)) {
      sign <- if (subset_size %% 2L == 1L) 1 else -1
      subsets <- utils::combn(seq_len(k), subset_size, simplify = FALSE)
      for (subset in subsets) {
        unscaled <- unscaled + sign * .subset_outer_sum(
          scores, data, cluster_cols[subset]
        )
      }
    }
  }

  n_eff <- min(counts)
  matrix_value <- (n_eff / nrow(data)^2) * unscaled
  matrix_value <- (matrix_value + t(matrix_value)) / 2
  result <- list(
    matrix = matrix_value,
    method = method,
    n_eff = as.integer(n_eff),
    cluster_counts = counts,
    n_cells = .number_of_groups(data, cluster_cols),
    n_observations = nrow(data)
  )
  class(result) <- "fullsampleDML_cluster_moment"
  result
}

.validate_cluster_inputs <- function(scores, data, cluster_cols) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame", call. = FALSE)
  }
  if (!is.character(cluster_cols) || !length(cluster_cols)) {
    stop("cluster_cols must contain at least one column", call. = FALSE)
  }
  if (anyDuplicated(cluster_cols)) {
    stop("cluster_cols must not contain repeated column names", call. = FALSE)
  }
  missing_cols <- setdiff(cluster_cols, names(data))
  if (length(missing_cols)) {
    stop(
      "cluster columns not found in data: ",
      paste(missing_cols, collapse = ", "), call. = FALSE
    )
  }
  if (anyNA(data[cluster_cols])) {
    stop("cluster columns must not contain missing values", call. = FALSE)
  }
  if (!is.matrix(scores) || !is.numeric(scores) || length(dim(scores)) != 2L) {
    stop("scores must be a two-dimensional numeric n-by-q matrix", call. = FALSE)
  }
  if (!nrow(data)) {
    stop("data and scores must contain at least one observation", call. = FALSE)
  }
  if (nrow(scores) != nrow(data)) {
    stop("scores and data must have the same number of rows", call. = FALSE)
  }
  if (!ncol(scores)) {
    stop("scores must contain at least one moment", call. = FALSE)
  }
  if (any(!is.finite(scores))) {
    stop("scores must contain only finite values", call. = FALSE)
  }
  counts <- vapply(data[cluster_cols], function(x) length(unique(x)), integer(1))
  invalid <- names(counts)[counts < 2L]
  if (length(invalid)) {
    stop(
      "each clustering dimension must contain at least two clusters; invalid columns: ",
      paste(invalid, collapse = ", "), call. = FALSE
    )
  }
  list(scores = scores, cluster_cols = cluster_cols, counts = counts)
}

.group_factor <- function(data, columns) {
  factors <- lapply(data[columns], function(x) factor(x, exclude = NULL))
  if (length(factors) == 1L) {
    return(factors[[1L]])
  }
  do.call(interaction, c(factors, list(drop = TRUE, lex.order = TRUE)))
}

.number_of_groups <- function(data, columns) {
  length(unique(.group_factor(data, columns)))
}

.subset_outer_sum <- function(scores, data, columns) {
  totals <- rowsum(scores, .group_factor(data, columns), reorder = FALSE)
  crossprod(totals)
}

.validate_variance_type <- function(value, argument) {
  if (!is.character(value) || length(value) != 1L ||
      is.na(value) || !value %in% c("PSD", "CGM")) {
    stop(argument, " must be exactly 'PSD' or 'CGM'", call. = FALSE)
  }
  value
}
