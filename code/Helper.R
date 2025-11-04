#' Complete SHAP Analysis Package for PatientLevelPrediction
#'
#' This file contains all functions needed for SHAP analysis of PLP models


#' Calculate SHAP values for PatientLevelPrediction models
#'
#' @param plpResult A plpResult object from PatientLevelPrediction
#' @param plpData The plpData object used for training
#' @param n_samples Number of test samples to explain (default: 100)
#' @param background_size Number of background samples for SHAP (default: 100)
#' @param nsim Number of Monte Carlo simulations for SHAP (default: 50)
#' @param parallel Use parallel processing (default: TRUE, auto-disabled for XGBoost)
#' @param n_cores Number of cores to use for parallel processing (default: detectCores() - 1)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A plpShap object containing SHAP values and metadata
#'
#' @details
#' This function works with multiple model types:
#' - Lasso Logistic Regression (supports parallel)
#' - Random Forest (supports parallel)
#' - XGBoost/Gradient Boosting (parallel automatically disabled due to serialization issues)
#'
calculate_plp_shap <- function(plpResult,
                               plpData,
                               n_samples = 100,
                               background_size = 100,
                               nsim = 50,
                               parallel = TRUE,
                               n_cores = NULL,
                               verbose = TRUE) {
  
  if (verbose) cat("=== SHAP Calculation ===\n\n")
  
  model <- plpResult$model
  population <- plpResult$prediction
  covariateData <- plpData$covariateData
  
  # Detect model type
  model_type <- attr(model, "modelType")
  is_xgboost <- !is.null(model$model) && "xgb.Booster" %in% class(model$model)
  
  if (verbose) {
    cat("Model type:", if(!is.null(model_type)) model_type else "unknown", "\n")
  }
  
  # Get train/test split
  test_pop <- population[population$evaluationType == "Test", ]
  train_pop <- population[population$evaluationType == "Train", ]
  
  if (nrow(test_pop) > n_samples) {
    test_pop <- test_pop[sample(nrow(test_pop), n_samples), ]
  }
  
  if (nrow(train_pop) > background_size) {
    train_pop <- train_pop[sample(nrow(train_pop), background_size), ]
  }
  
  if (verbose) {
    cat("Test samples:", nrow(test_pop), "\n")
    cat("Background samples:", nrow(train_pop), "\n")
  }
  
  # Get covariates and metadata BEFORE parallel setup
  selected_covariates <- unique(model$covariateImportance$covariateId)
  if (verbose) cat("Features:", length(selected_covariates), "\n")
  
  covariate_ref <- covariateData$covariateRef |> dplyr::collect()
  original_metadata <- plpData$metaData
  
  # Auto-disable parallel for XGBoost
  if (is_xgboost && parallel) {
    if (verbose) {
      cat("\nXGBoost model detected - parallel processing disabled\n")
      cat("(XGBoost models cannot be serialized to parallel workers)\n")
    }
    parallel <- FALSE
  }
  
  # Setup parallel processing if enabled (AFTER all variables are defined)
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    n_cores <- max(1, min(n_cores, parallel::detectCores() - 1))
    
    if (verbose) cat("\nUsing parallel processing with", n_cores, "cores\n")
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    # Export necessary objects to workers
    parallel::clusterExport(cl, 
                           c("covariate_ref", "original_metadata"), 
                           envir = environment())
    
    # Load required packages on workers
    parallel::clusterEvalQ(cl, {
      library(PatientLevelPrediction)
      library(dplyr)
    })
    
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  
  if (verbose) cat("\n")
  
  # Create matrices
  create_matrix <- function(rowIds, covariateIds) {
    subset_cov <- covariateData$covariates |>
      dplyr::filter(rowId %in% rowIds & covariateId %in% covariateIds) |>
      dplyr::collect()
    
    n_rows <- length(rowIds)
    n_cols <- length(covariateIds)
    dense_matrix <- matrix(0, nrow = n_rows, ncol = n_cols)
    colnames(dense_matrix) <- as.character(covariateIds)
    
    if (nrow(subset_cov) > 0) {
      row_mapping <- setNames(1:n_rows, as.character(rowIds))
      col_mapping <- setNames(1:n_cols, as.character(covariateIds))
      row_indices <- row_mapping[as.character(subset_cov$rowId)]
      col_indices <- col_mapping[as.character(subset_cov$covariateId)]
      valid_idx <- !is.na(row_indices) & !is.na(col_indices)
      if (any(valid_idx)) {
        dense_matrix[cbind(row_indices[valid_idx], col_indices[valid_idx])] <- 
          subset_cov$covariateValue[valid_idx]
      }
    }
    return(dense_matrix)
  }
  
  background_matrix <- create_matrix(train_pop$rowId, selected_covariates)
  test_matrix <- create_matrix(test_pop$rowId, selected_covariates)
  
  if (verbose) {
    cat("Background matrix:", nrow(background_matrix), "x", ncol(background_matrix), "\n")
    cat("Test matrix:", nrow(test_matrix), "x", ncol(test_matrix), "\n\n")
  }
  
  # Prediction wrapper
  predict_wrapper <- function(object, newdata) {
    
    if (is.matrix(newdata)) {
      newdata <- as.data.frame(newdata)
    }
    
    covariate_ids <- as.numeric(colnames(newdata))
    n_obs <- nrow(newdata)
    
    # Build covariates
    batch_covariates_list <- lapply(1:n_obs, function(idx) {
      row_data <- as.numeric(newdata[idx, ])
      nonzero_idx <- which(row_data != 0)
      if (length(nonzero_idx) > 0) {
        data.frame(
          rowId = idx,
          covariateId = covariate_ids[nonzero_idx],
          covariateValue = row_data[nonzero_idx]
        )
      } else {
        NULL
      }
    })
    
    batch_covariates_list <- Filter(Negate(is.null), batch_covariates_list)
    
    if (length(batch_covariates_list) > 0) {
      batch_covariates <- do.call(rbind, batch_covariates_list)
    } else {
      batch_covariates <- data.frame(
        rowId = integer(0),
        covariateId = numeric(0),
        covariateValue = numeric(0)
      )
    }
    
    batch_population <- data.frame(
      rowId = 1:n_obs,
      subjectId = 1:n_obs,
      cohortStartDate = as.Date("2020-01-01"),
      outcomeCount = 0,
      timeAtRisk = 365,
      survivalTime = 365
    )
    
    batch_covariate_data <- list(
      covariates = batch_covariates,
      covariateRef = covariate_ref,
      analysisRef = data.frame(
        analysisId = 1,
        analysisName = "Demographics",
        domainId = "Demographics",
        startDay = 0,
        endDay = 0,
        isBinary = "Y",
        missingMeansZero = "Y"
      ),
      metaData = list(populationSize = n_obs)
    )
    
    class(batch_covariate_data) <- "CovariateData"
    attr(batch_covariate_data, "metaData") <- batch_covariate_data$metaData
    
    batch_plpData <- list(
      covariateData = batch_covariate_data,
      cohorts = batch_population,
      outcomes = data.frame(
        rowId = integer(0),
        outcomeId = integer(0),
        daysToEvent = integer(0)
      ),
      metaData = original_metadata
    )
    
    class(batch_plpData) <- "plpData"
    
    # Predict with output suppression
    pred <- suppressMessages(suppressWarnings(
      PatientLevelPrediction::predictPlp(
        plpModel = object,
        plpData = batch_plpData,
        population = batch_population
      )
    ))
    
    pred_values <- pred$value[match(1:n_obs, pred$rowId)]
    
    if (any(is.na(pred_values))) {
      mean_pred <- mean(pred_values, na.rm = TRUE)
      if (is.na(mean_pred)) mean_pred <- 0.5
      pred_values[is.na(pred_values)] <- mean_pred
    }
    
    return(as.numeric(pred_values))
  }
  
  # Test wrapper
  if (verbose) cat("Testing prediction wrapper...\n")
  test_pred <- invisible(predict_wrapper(model, test_matrix[1:min(3, nrow(test_matrix)), , drop = FALSE]))
  
  if (verbose) {
    cat("Test predictions:", paste(round(test_pred, 4), collapse = ", "), "\n")
    if (all(test_pred == 0.5)) {
      warning("Test predictions all 0.5 - wrapper may be failing")
    } else {
      cat("✓ Wrapper working\n\n")
    }
  }
  
  # Calculate SHAP
  if (verbose) {
    cat("Calculating SHAP values (this may take several minutes)...\n")
    if (parallel) {
      cat("Using parallel processing with", n_cores, "cores\n")
    } else {
      cat("Running sequentially (no parallel processing)\n")
    }
    cat("\n")
  }
  
  # Calculate SHAP values
  shap_values <- fastshap::explain(
    object = model,
    X = background_matrix,
    pred_wrapper = function(object, newdata) {
      invisible(predict_wrapper(object, newdata))
    },
    nsim = nsim,
    newdata = test_matrix,
    adjust = TRUE,
    parallel = parallel
  )
  
  if (verbose) {
    cat("\n\nSHAP calculation complete!\n")
    cat("Dimensions:", nrow(shap_values), "x", ncol(shap_values), "\n")
    
    non_na <- sum(!is.na(shap_values))
    total <- length(shap_values)
    na_pct <- 100 * (1 - non_na / total)
    
    cat("NA percentage:", round(na_pct, 2), "%\n")
    
    if (non_na > 0) {
      cat("Range: [", round(min(shap_values, na.rm = TRUE), 4), ", ",
          round(max(shap_values, na.rm = TRUE), 4), "]\n", sep = "")
    }
    cat("\n")
  }
  
  # Add feature names
  if (nrow(covariate_ref) > 0) {
    cov_names <- setNames(
      covariate_ref$covariateName,
      as.character(covariate_ref$covariateId)
    )
    
    matched_names <- cov_names[colnames(shap_values)]
    matched_names[is.na(matched_names)] <- colnames(shap_values)[is.na(matched_names)]
    
    # Clean feature names: remove prefix before ":" and capitalize first letter
    matched_names <- sapply(matched_names, function(name) {
      # Split on colon and take the part after it
      parts <- strsplit(name, ":", fixed = TRUE)[[1]]
      if (length(parts) > 1) {
        # Take everything after the colon and trim whitespace
        clean_name <- trimws(parts[2])
      } else {
        # No colon found, use original name
        clean_name <- name
      }
      # Capitalize first letter
      if (nchar(clean_name) > 0) {
        substr(clean_name, 1, 1) <- toupper(substr(clean_name, 1, 1))
      }
      return(clean_name)
    }, USE.NAMES = FALSE)
    
    colnames(shap_values) <- matched_names
  }
  
  # Get predictions
  predictions <- invisible(predict_wrapper(model, test_matrix))
  
  result <- list(
    shap_values = shap_values,
    rowIds = test_pop$rowId,
    feature_names = colnames(shap_values),
    covariateIds = selected_covariates,
    method = "fastshap",
    background_size = nrow(background_matrix),
    predictions = predictions,
    test_matrix = test_matrix,
    background_matrix = background_matrix,
    model_type = model_type,
    parallel = parallel,
    n_cores = if(parallel) n_cores else NA
  )
  
  class(result) <- c("plpShap", "list")
  return(result)
}


#' Summary method for plpShap objects
#'
#' @param object A plpShap object
#' @param ... Additional arguments (unused)
#'
summary.plpShap <- function(object, ...) {
  cat("SHAP Values Summary\n")
  cat("===================\n\n")
  cat("Model type:", if(!is.null(object$model_type)) object$model_type else "unknown", "\n")
  cat("Method:", object$method, "\n")
  cat("Parallel processing:", object$parallel, "\n")
  if (object$parallel && !is.na(object$n_cores)) {
    cat("Cores used:", object$n_cores, "\n")
  }
  cat("\nNumber of observations:", nrow(object$shap_values), "\n")
  cat("Number of features:", ncol(object$shap_values), "\n")
  cat("Background sample size:", object$background_size, "\n\n")
  
  cat("Prediction range:", round(min(object$predictions), 4), "to", 
      round(max(object$predictions), 4), "\n\n")
  
  cat("Top 10 Most Important Features (Mean |SHAP|):\n")
  importance <- colMeans(abs(object$shap_values))
  top_features <- head(sort(importance, decreasing = TRUE), 10)
  print(round(top_features, 6))
  
  invisible(object)
}


#' Plot SHAP values and save to file
#'
#' @param shap_result A plpShap object
#' @param plot_type Type of plot: "importance", "waterfall", or "beeswarm"
#' @param observation Which observation to plot for waterfall (default: 1)
#' @param max_features Maximum number of features to display (default: 20)
#' @param save_path File path to save the plot
#' @param title The plot's title
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @param dpi Resolution in dots per inch (default: 300)
#'
#'
plot_plp_shap <- function(shap_result,
                          plot_type = "importance",
                          observation = 1,
                          max_features = 20,
                          save_path,
                          title = NULL,
                          width = 10,
                          height = 8,
                          dpi = 600) {
  
  if (missing(save_path) || is.null(save_path)) {
    stop("save_path is required. Please provide a file path to save the plot.")
  }
  
  shap_values <- shap_result$shap_values
  
  if (plot_type == "importance") {
    # Feature importance plot
    importance_df <- data.frame(
      feature = colnames(shap_values),
      importance = colMeans(abs(shap_values))
    )
    importance_df <- importance_df[order(-importance_df$importance), ]
    importance_df <- head(importance_df, max_features)
    
    # Set default title if not provided
    if (is.null(title)) {
      title <- "Feature Importance (Mean |SHAP|)"
    }
    
    p <- ggplot2::ggplot(
      importance_df,
      ggplot2::aes(x = reorder(feature, importance), y = importance)
    ) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = title,
        x = "Feature",
        y = "Mean Absolute SHAP Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(size = 14, face = "bold")
      )
    
  } else if (plot_type == "waterfall") {
    if (observation > nrow(shap_values)) {
      stop("Observation index out of bounds")
    }
    
    obs_shap <- shap_values[observation, ]
    obs_df <- data.frame(
      feature = names(obs_shap),
      shap_value = as.numeric(obs_shap)
    )
    obs_df <- obs_df[order(-abs(obs_df$shap_value)), ]
    obs_df <- head(obs_df, 15)
    obs_df$feature <- factor(obs_df$feature, levels = rev(obs_df$feature))
    
    # Set default title if not provided
    if (is.null(title)) {
      title <- paste("SHAP Values for Observation", observation)
    }
    
    p <- ggplot2::ggplot(
      obs_df,
      ggplot2::aes(x = feature, y = shap_value, fill = shap_value > 0)
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(
        values = c("TRUE" = "#FF0D57", "FALSE" = "#1E88E5"),
        labels = c("Increases risk", "Decreases risk")
      ) +
      ggplot2::labs(
        title = title,
        subtitle = paste("Prediction:", round(shap_result$predictions[observation], 4)),
        x = "Feature",
        y = "SHAP Value",
        fill = "Effect"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 16, face = "bold"),
        axis.text = ggplot2::element_text(size = 14)
      )
    
  } else if (plot_type == "beeswarm") {
    
    # Get feature values from test_matrix
    test_matrix <- shap_result$test_matrix
    
    if (is.null(test_matrix)) {
      stop("test_matrix not found in shap_result. Cannot create beeswarm plot with feature values.")
    }
    
    # Convert SHAP values to long format
    shap_long <- as.data.frame(shap_values)
    shap_long$observation <- 1:nrow(shap_long)
    shap_long <- tidyr::pivot_longer(
      shap_long,
      cols = -observation,
      names_to = "feature",
      values_to = "shap_value"
    )
    
    # Convert feature values to long format
    feature_long <- as.data.frame(test_matrix)
    feature_long$observation <- 1:nrow(feature_long)
    
    # Match column names from SHAP (which have feature names) to test_matrix (which has covariate IDs)
    covariate_ids <- shap_result$covariateIds
    if (!is.null(covariate_ids)) {
      colnames(test_matrix) <- as.character(covariate_ids)
      feature_long <- as.data.frame(test_matrix)
      feature_long$observation <- 1:nrow(feature_long)
    }
    
    feature_long <- tidyr::pivot_longer(
      feature_long,
      cols = -observation,
      names_to = "feature_id",
      values_to = "feature_value"
    )
    
    # Create mapping from feature names back to IDs
    feature_name_to_id <- setNames(
      as.character(shap_result$covariateIds),
      shap_result$feature_names
    )
    
    # Add feature_id to shap_long
    shap_long$feature_id <- feature_name_to_id[shap_long$feature]
    
    # Merge SHAP values with feature values
    shap_long <- shap_long |>
      dplyr::left_join(
        feature_long,
        by = c("observation", "feature_id")
      )
    
    # Calculate importance and get top features
    importance <- shap_long |>
      dplyr::group_by(feature) |>
      dplyr::summarise(importance = mean(abs(shap_value)), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(importance))
    
    top_features <- head(importance$feature, max_features)
    shap_long <- shap_long |>
      dplyr::filter(feature %in% top_features) |>
      dplyr::mutate(feature = factor(feature, levels = rev(top_features)))
    
    # Normalize feature values per feature (0 to 1 scale within each feature)
    # This makes binary and continuous variables equally visible
    shap_long <- shap_long |>
      dplyr::group_by(feature) |>
      dplyr::mutate(
        feature_value_normalized = {
          min_val <- min(feature_value, na.rm = TRUE)
          max_val <- max(feature_value, na.rm = TRUE)
          # Handle case where all values are the same
          if (max_val == min_val) {
            rep(0.5, dplyr::n())
          } else {
            (feature_value - min_val) / (max_val - min_val)
          }
        }
      ) |>
      dplyr::ungroup()
    
    # Set default title if not provided
    if (is.null(title)) {
      title <- "SHAP Summary Plot"
    }
    
    # Create the plot colored by normalized feature value
    p <- ggplot2::ggplot(
      shap_long,
      ggplot2::aes(x = shap_value, y = feature)
    ) +
      ggplot2::geom_jitter(
        ggplot2::aes(color = feature_value_normalized),
        alpha = 0.6,
        height = 0.2,
        size = 1.5
      ) +
      ggplot2::scale_color_gradient(
        low = "#FDE725", # Yellow for low values (within feature)
        high = "#440154", # Purple for high values (within feature)
        name = "Feature\nValue",
        labels = c("Low", "High"),
        breaks = c(0, 1)
      ) +
      ggplot2::labs(
        title = title,
        x = "SHAP Value (impact on model output)",
        y = "Feature"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 26, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 10),
        axis.text = ggplot2::element_text(size = 24),
        axis.title.y = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(size = 24)
      )
    
  } else {
    stop("Invalid plot_type. Choose 'importance', 'waterfall', or 'beeswarm'")
  }
  
  # Create directory if it doesn't exist
  dir_path <- dirname(save_path)
  if (!dir.exists(dir_path) && dir_path != ".") {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Save using ragg (works without X11)
  ggplot2::ggsave(
    filename = save_path,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    device = ragg::agg_png
  )
  
  cat("Plot saved to:", save_path, "\n")
  
  return(invisible(save_path))
}


#' Export SHAP values to data frame
#'
#' @param shap_result Result from calculate_plp_shap()
#' @return A data frame with SHAP values in long format
#'
export_shap_values <- function(shap_result) {
  
  shap_df <- as.data.frame(shap_result$shap_values)
  shap_df$rowId <- shap_result$rowIds
  shap_df$prediction <- shap_result$predictions
  
  # Convert to long format
  shap_long <- shap_df |>
    tidyr::pivot_longer(
      cols = -c(rowId, prediction),
      names_to = "feature",
      values_to = "shap_value"
    )
  
  return(shap_long)
}

' Calculate feature importance ranks
#'
#' @param shap_result Result from calculate_plp_shap()
#' @return A data frame with feature rankings
#'
#' @details
#' For each patient, features are ranked by absolute SHAP value (1 = most important).
#' The function returns:
#' - mean_rank: Average rank across all patients (lower = more consistently important)
#' - median_rank: Median rank across all patients
#' - mean_abs_shap: Average absolute SHAP value (overall importance magnitude)
#' - times_top_5: How many patients had this feature in their top 5
#' - times_top_10: How many patients had this feature in their top 10
#'
get_feature_ranks <- function(shap_result) {
  
  shap_values <- shap_result$shap_values
  n_patients <- nrow(shap_values)
  
  # Calculate absolute SHAP values
  abs_shap <- abs(shap_values)
  
  # For each patient (row), rank features by absolute SHAP value
  # Rank 1 = highest absolute SHAP (most important)
  ranks <- t(apply(abs_shap, 1, function(row) {
    rank(-row, ties.method = "average")
  }))
  
  # Calculate statistics for each feature
  feature_stats <- data.frame(
    feature = colnames(shap_values),
    mean_rank = colMeans(ranks),
    median_rank = apply(ranks, 2, median),
    sd_rank = apply(ranks, 2, sd),
    min_rank = apply(ranks, 2, min),
    max_rank = apply(ranks, 2, max),
    mean_abs_shap = colMeans(abs_shap),
    median_abs_shap = apply(abs_shap, 2, median),
    times_top_5 = colSums(ranks <= 5),
    times_top_10 = colSums(ranks <= 10),
    pct_top_5 = 100 * colSums(ranks <= 5) / n_patients,
    pct_top_10 = 100 * colSums(ranks <= 10) / n_patients
  )
  
  # Sort by mean rank (most consistently important first)
  feature_stats <- feature_stats[order(feature_stats$mean_rank), ]
  rownames(feature_stats) <- NULL
  
  return(feature_stats)
}


#' Print feature rank summary
#'
#' @param shap_result Result from calculate_plp_shap()
#' @param top_n Number of top features to display (default: 20)
#'
print_feature_ranks <- function(shap_result, top_n = 20) {
  
  ranks <- get_feature_ranks(shap_result)
  ranks_display <- head(ranks, top_n)
  
  cat("Feature Importance Rankings\n")
  cat("============================\n\n")
  cat("Based on", nrow(shap_result$shap_values), "patients\n")
  cat("Lower mean rank = more consistently important across patients\n\n")
  
  # Format for display
  ranks_display$mean_rank <- round(ranks_display$mean_rank, 2)
  ranks_display$mean_abs_shap <- round(ranks_display$mean_abs_shap, 4)
  
  print(ranks_display[, c("feature", "mean_rank", "median_rank", "mean_abs_shap", 
                          "times_top_5", "pct_top_5")], 
        row.names = FALSE)
  
  cat("\n")
  cat("Columns:\n")
  cat("  mean_rank: Average rank across patients (1 = most important)\n")
  cat("  median_rank: Median rank across patients\n")
  cat("  mean_abs_shap: Average absolute SHAP value\n")
  cat("  times_top_5: Number of patients where feature was in top 5\n")
  cat("  pct_top_5: Percentage of patients where feature was in top 5\n")
  
  invisible(ranks)
}
