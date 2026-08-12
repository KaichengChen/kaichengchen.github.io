#' Specify one mlr3 nuisance learner
#'
#' @param learner An `mlr3` learner or `mlr3tuning` auto-tuner.
#' @param target One target column name or a selector function.
#' @param features Feature column names or a selector function.
#' @param groups Optional one or more clustering column names used for tuning.
#' @param predict_type The mlr3 prediction type, normally `"response"`.
#' @return A `fullsampleDML_nuisance_spec` object.
#' @export
nuisance_spec <- function(learner, target, features, groups = NULL,
                          predict_type = "response") {
  if (is.null(learner)) stop("learner must be supplied", call. = FALSE)
  if (!(is.character(target) || is.function(target))) {
    stop("target must be a column name or selector function", call. = FALSE)
  }
  if (!(is.character(features) || is.function(features))) {
    stop("features must be column names or a selector function", call. = FALSE)
  }
  if (!is.null(groups) && !inherits(groups, "fullsampleDML_multiway_kfold") &&
      (!is.character(groups) || !length(groups))) {
    stop(
      "groups must be NULL, a nonempty character vector, or multiway_group_kfold()",
      call. = FALSE
    )
  }
  if (inherits(groups, "fullsampleDML_multiway_kfold") &&
      is.null(groups$group_cols)) {
    stop("multiway_group_kfold() must specify group_cols here", call. = FALSE)
  }
  if (!is.character(predict_type) || length(predict_type) != 1L) {
    stop("predict_type must be one character value", call. = FALSE)
  }
  result <- list(
    learner = learner,
    target = target,
    features = features,
    groups = groups,
    predict_type = predict_type
  )
  class(result) <- "fullsampleDML_nuisance_spec"
  result
}

#' Fit mlr3 nuisance components on the complete sample
#'
#' Tuning folds select hyperparameters. The selected learner is then refitted
#' on every observation and predicts those same observations; this is not DML
#' cross-fitting.
#'
#' @param data A data frame.
#' @param specs A named list created with [nuisance_spec()].
#' @param cluster_cols Optional candidate clustering columns. When a tuned
#'   component does not set `groups`, the dimension with the fewest clusters is
#'   used, with ties resolved by the order given here.
#' @return A `fullsampleDML_nuisance_result` object.
#' @export
fit_mlr3_nuisances <- function(data, specs, cluster_cols = NULL) {
  if (!requireNamespace("mlr3", quietly = TRUE)) {
    stop("fit_mlr3_nuisances() requires the suggested package 'mlr3'", call. = FALSE)
  }
  if (!is.data.frame(data)) stop("data must be a data frame", call. = FALSE)
  if (!is.list(specs) || !length(specs)) {
    stop("specs must contain at least one nuisance learner", call. = FALSE)
  }
  if (is.null(names(specs)) || any(!nzchar(names(specs))) || anyDuplicated(names(specs))) {
    stop("every nuisance name must be a unique nonempty string", call. = FALSE)
  }
  if (!is.null(cluster_cols)) {
    .validate_group_columns(data, cluster_cols, "cluster_cols")
  }

  predictions <- models <- group_columns <- group_counts <- list()
  for (name in names(specs)) {
    spec <- specs[[name]]
    if (!inherits(spec, "fullsampleDML_nuisance_spec")) {
      stop("every component in specs must be created by nuisance_spec()", call. = FALSE)
    }
    learner <- tryCatch(
      spec$learner$clone(deep = TRUE),
      error = function(error) {
        stop("nuisance '", name, "' learner must be an mlr3 learner", call. = FALSE)
      }
    )
    is_tuned <- inherits(learner, "AutoTuner")
    selected_groups <- if (is_tuned) {
      .resolve_groups(data, spec$groups, cluster_cols)
    } else {
      list(columns = character(), counts = integer())
    }
    task_data <- .make_task_data(data, spec$target, spec$features, name)

    if (length(selected_groups$columns) == 1L) {
      task_data$data$.fullsampleDML_group <- data[[selected_groups$columns]]
    }
    task <- mlr3::TaskRegr$new(
      id = paste0("fullsampleDML_", name),
      backend = task_data$data,
      target = task_data$target
    )
    if (length(selected_groups$columns) == 1L) {
      task$set_col_roles(".fullsampleDML_group", roles = "group")
    } else if (length(selected_groups$columns) > 1L) {
      if (!requireNamespace("mlr3tuning", quietly = TRUE)) {
        stop(
          "strict multiway tuning requires the suggested package 'mlr3tuning'",
          call. = FALSE
        )
      }
      template <- if (inherits(spec$groups, "fullsampleDML_multiway_kfold")) {
        spec$groups
      } else {
        multiway_group_kfold()
      }
      splits <- multiway_group_splits(
        data[selected_groups$columns],
        n_splits_per_dimension = template$n_splits_per_dimension,
        shuffle = template$shuffle,
        random_state = template$random_state
      )
      custom <- mlr3::rsmp("custom")
      custom$instantiate(task, train_sets = splits$train, test_sets = splits$test)
      learner$instance_args$resampling <- custom
    }

    learner$predict_type <- spec$predict_type
    tryCatch(
      learner$train(task),
      error = function(error) {
        stop(
          "nuisance '", name, "' could not be trained: ", conditionMessage(error),
          call. = FALSE
        )
      }
    )
    prediction <- tryCatch(
      learner$predict(task),
      error = function(error) {
        stop(
          "nuisance '", name, "' could not predict: ", conditionMessage(error),
          call. = FALSE
        )
      }
    )
    values <- prediction$response
    if (!is.numeric(values) || length(values) != nrow(data) ||
        any(!is.finite(values))) {
      stop(
        "nuisance '", name,
        "' returned non-finite predictions or an invalid first dimension",
        call. = FALSE
      )
    }
    predictions[[name]] <- as.numeric(values)
    models[[name]] <- learner
    group_columns[[name]] <- selected_groups$columns
    group_counts[[name]] <- selected_groups$counts
  }

  result <- list(
    predictions = predictions,
    models = models,
    group_columns = group_columns,
    group_counts = group_counts
  )
  class(result) <- "fullsampleDML_nuisance_result"
  result
}

#' Build an mlr3 nuisance callback for fit_gmm
#'
#' @param specs A named list of nuisance specifications.
#' @return A callback accepted by [fit_gmm()].
#' @export
mlr3_nuisance_fit <- function(specs) {
  callback <- function(data, cluster_cols = NULL) {
    fit_mlr3_nuisances(data, specs, cluster_cols = cluster_cols)
  }
  class(callback) <- c("fullsampleDML_mlr3_nuisance_fitter", class(callback))
  callback
}

.validate_group_columns <- function(data, columns, role) {
  if (!is.character(columns) || !length(columns)) {
    stop(role, " must contain at least one column", call. = FALSE)
  }
  if (anyDuplicated(columns)) {
    stop(role, " must not contain repeated column names", call. = FALSE)
  }
  missing_cols <- setdiff(columns, names(data))
  if (length(missing_cols)) {
    stop(role, " not found in data: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  if (anyNA(data[columns])) {
    stop(role, " must not contain missing values", call. = FALSE)
  }
  stats::setNames(
    vapply(data[columns], function(x) length(unique(x)), integer(1)),
    columns
  )
}

.resolve_groups <- function(data, groups, cluster_cols) {
  if (is.null(groups)) {
    if (is.null(cluster_cols)) {
      return(list(columns = character(), counts = integer()))
    }
    counts <- .validate_group_columns(data, cluster_cols, "cluster_cols")
    selected <- cluster_cols[which.min(counts)]
    return(list(columns = selected, counts = counts[selected]))
  }
  if (inherits(groups, "fullsampleDML_multiway_kfold")) {
    columns <- groups$group_cols
    counts <- .validate_group_columns(data, columns, "group columns")
    return(list(columns = columns, counts = counts))
  }
  counts <- .validate_group_columns(data, groups, "group columns")
  list(columns = groups, counts = counts)
}

.make_task_data <- function(data, target, features, name) {
  feature_values <- if (is.function(features)) features(data) else data[features]
  if (is.null(feature_values)) {
    stop("nuisance '", name, "' feature selector returned NULL", call. = FALSE)
  }
  feature_values <- as.data.frame(feature_values, check.names = FALSE)
  if (!ncol(feature_values) || nrow(feature_values) != nrow(data)) {
    stop("nuisance '", name, "' has invalid features", call. = FALSE)
  }
  if (is.function(target)) {
    target_values <- target(data)
  } else {
    if (length(target) != 1L || !target %in% names(data)) {
      stop("target column not found in data: ", paste(target, collapse = ", "),
           call. = FALSE)
    }
    target_values <- data[[target]]
  }
  if (!is.numeric(target_values) || length(target_values) != nrow(data) ||
      any(!is.finite(target_values))) {
    stop("nuisance '", name, "' target must be one finite numeric vector",
         call. = FALSE)
  }
  names(feature_values) <- make.unique(paste0("feature_", seq_len(ncol(feature_values))))
  target_name <- ".fullsampleDML_target"
  feature_values[[target_name]] <- as.numeric(target_values)
  list(data = feature_values, target = target_name)
}
