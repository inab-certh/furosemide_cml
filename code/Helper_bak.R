#' Calculate SHAP values for PatientLevelPrediction models (Parallelized)
#'
#' @param plpResult A plpResult object from PatientLevelPrediction::runPlp()
#' @param plpData The plpData object used for training/testing
#' @param population The population data frame (optional, defaults to plpResult$prediction)
#' @param indices Indices of patients to explain (default: uses test set or all data)
#' @param subset Which subset to use: "test", "train", "all", or NULL (default: "test")
#' @param n_samples Number of samples to explain (if NULL, uses all from subset)
#' @param background_size Number of background samples for SHAP (default: 100)
#' @param method Method for SHAP calculation: "fastshap" or "python"
#' @param nsim Number of simulations for fastshap (default: 50)
#' @param parallel Use parallel processing (default: TRUE)
#' @param n_cores Number of cores to use (default: detectCores() - 1)
#' @return A list containing SHAP values matrix and metadata
#'
#' Calculate SHAP values for PatientLevelPrediction models (with suppressed messages)
#'
calculate_plp_shap <- function(plpResult,
                               plpData,
                               population = NULL,
                               indices = NULL,
                               subset = "test",
                               n_samples = NULL,
                               background_size = 100,
                               method = "fastshap",
                               nsim = 50,
                               parallel = TRUE,
                               n_cores = NULL,
                               verbose = TRUE) {
  
  if (!requireNamespace("PatientLevelPrediction", quietly = TRUE)) {
    stop("PatientLevelPrediction package required")
  }
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package required")
  }
  
  if (parallel) {
    if (!requireNamespace("foreach", quietly = TRUE) || 
        !requireNamespace("doParallel", quietly = TRUE)) {
      warning("foreach and doParallel packages required for parallel processing.")
      parallel <- FALSE
    }
  }
  
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    n_cores <- max(1, min(n_cores, parallel::detectCores() - 1))
    if (verbose) cat("Setting up parallel processing with", n_cores, "cores...\n")
    
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    parallel::clusterEvalQ(cl, {
      library(PatientLevelPrediction)
    })
    
    on.exit(parallel::stopCluster(cl), add = TRUE)
  }
  
  if (verbose) cat("Extracting model and data...\n")
  
  model <- plpResult$model
  model_type <- attr(model, "modelType")
  if (verbose) cat("Model type:", if(!is.null(model_type)) model_type else "unknown", "\n")
  
  covariateData <- plpData$covariateData
  
  if (is.null(population)) {
    if (!is.null(plpResult$prediction)) {
      population <- plpResult$prediction
    } else {
      stop("No population found.")
    }
  }
  
  if (verbose) cat("Checking population structure...\n")
  
  if (!is.null(indices)) {
    selected_pop <- population[population$rowId %in% indices, ]
  } else if (!is.null(subset)) {
    if ("indexes" %in% names(population)) {
      if (subset == "test") {
        selected_pop <- population[population$indexes > 0, ]
      } else if (subset == "train") {
        selected_pop <- population[population$indexes < 0, ]
      } else if (subset == "all") {
        selected_pop <- population
      } else {
        stop("Invalid subset.")
      }
    } else {
      selected_pop <- population
    }
  } else {
    selected_pop <- population
  }
  
  if (nrow(selected_pop) == 0) {
    stop("No data found for subset '", subset, "'. Try subset='all'")
  }
  
  if (!is.null(n_samples) && nrow(selected_pop) > n_samples) {
    if (verbose) cat("Sampling", n_samples, "from", nrow(selected_pop), "observations\n")
    selected_pop <- selected_pop[sample(nrow(selected_pop), n_samples), ]
  }
  
  indices <- selected_pop$rowId
  if (verbose) cat("Found", length(indices), "patients to explain\n")
  
  if (!is.null(model$covariateImportance)) {
    selected_covariates <- unique(model$covariateImportance$covariateId)
    if (verbose) cat("Using", length(selected_covariates), "covariates from model\n")
  } else {
    if (verbose) cat("Collecting all covariates from database...\n")
    all_covariates <- covariateData$covariates |> 
      dplyr::select(covariateId) |>
      dplyr::distinct() |>
      dplyr::collect()
    selected_covariates <- unique(all_covariates$covariateId)
    if (verbose) cat("Using all", length(selected_covariates), "covariates\n")
  }
  
  if (verbose) cat("Preparing background data...\n")
  
  if ("indexes" %in% names(population)) {
    train_pop <- population[population$indexes < 0, ]
    if (nrow(train_pop) == 0) {
      train_pop <- population[!population$rowId %in% indices, ]
      if (nrow(train_pop) == 0) {
        train_pop <- population
      }
    }
  } else {
    train_pop <- population[!population$rowId %in% indices, ]
    if (nrow(train_pop) == 0) {
      train_pop <- population
    }
  }
  
  if (nrow(train_pop) > background_size) {
    background_ids <- sample(train_pop$rowId, background_size)
  } else {
    background_ids <- train_pop$rowId
    if (verbose) cat("Using all", length(background_ids), "samples as background\n")
  }
  
  create_dense_matrix <- function(rowIds, covariateIds) {
    if (verbose) cat("  Querying database for", length(rowIds), "rows...\n")
    
    subset_cov <- covariateData$covariates |>
      dplyr::filter(rowId %in% rowIds & covariateId %in% covariateIds) |>
      dplyr::collect()
    
    if (verbose) cat("  Retrieved", nrow(subset_cov), "non-zero covariate values\n")
    
    n_rows <- length(rowIds)
    n_cols <- length(covariateIds)
    
    dense_matrix <- matrix(0, nrow = n_rows, ncol = n_cols)
    colnames(dense_matrix) <- as.character(covariateIds)
    rownames(dense_matrix) <- as.character(rowIds)
    
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
    
    if (verbose) cat("  Created dense matrix:", nrow(dense_matrix), "x", ncol(dense_matrix), "\n")
    return(dense_matrix)
  }
  
  background_matrix <- create_dense_matrix(background_ids, selected_covariates)
  test_matrix <- create_dense_matrix(indices, selected_covariates)
  
  if (verbose) {
    cat("Background matrix:", nrow(background_matrix), "x", ncol(background_matrix), "\n")
    cat("Test matrix:", nrow(test_matrix), "x", ncol(test_matrix), "\n")
  }
  
  if (verbose) cat("Collecting covariate reference...\n")
  covariate_ref <- covariateData$covariateRef |>
    dplyr::collect()
  
  # Get metadata from original plpData
  original_metadata <- plpData$metaData
  if (is.null(original_metadata)) {
    original_metadata <- list()
  }
  
  if (verbose) cat("Creating prediction wrapper (messages will be suppressed)...\n")
  
  predict_wrapper <- function(object, newdata) {
    
    if (is.matrix(newdata)) {
      newdata <- as.data.frame(newdata)
    }
    
    covariate_ids <- as.numeric(colnames(newdata))
    n_obs <- nrow(newdata)
    
    batch_covariates_list <- list()
    batch_population <- data.frame(
      rowId = 1:n_obs,
      subjectId = 1:n_obs,
      cohortStartDate = as.Date("2020-01-01"),
      outcomeCount = 0,
      timeAtRisk = 365,
      survivalTime = 365,
      indexes = 1:n_obs
    )
    
    for (idx in 1:n_obs) {
      row_data <- newdata[idx, , drop = FALSE]
      nonzero_cols <- which(as.numeric(row_data) != 0)
      
      if (length(nonzero_cols) > 0) {
        batch_covariates_list[[idx]] <- data.frame(
          rowId = idx,
          covariateId = covariate_ids[nonzero_cols],
          covariateValue = as.numeric(row_data[nonzero_cols])
        )
      }
    }
    
    if (length(batch_covariates_list) > 0) {
      batch_covariates <- do.call(rbind, batch_covariates_list)
    } else {
      batch_covariates <- data.frame(
        rowId = integer(0),
        covariateId = numeric(0),
        covariateValue = numeric(0)
      )
    }
    
    # Create properly structured CovariateData with correct class
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
      metaData = list(
        populationSize = n_obs
      )
    )
    
    # Set proper class attributes
    class(batch_covariate_data) <- "CovariateData"
    attr(batch_covariate_data, "metaData") <- batch_covariate_data$metaData
    
    # Create properly structured plpData
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
    
    # Set proper class
    class(batch_plpData) <- "plpData"
    
    # Suppress all output and predict
    preds <- suppressMessages(suppressWarnings({
      tryCatch({
        pred <- PatientLevelPrediction::predictPlp(
          plpModel = object,
          plpData = batch_plpData,
          population = batch_population
        )
        
        result <- pred$value[match(1:n_obs, pred$rowId)]
        
        if (any(is.na(result))) {
          mean_pred <- mean(result, na.rm = TRUE)
          if (is.na(mean_pred)) mean_pred <- 0.5
          result[is.na(result)] <- mean_pred
        }
        
        result
        
      }, error = function(e) {
        if (verbose) cat("Prediction error:", e$message, "\n")
        rep(0.5, n_obs)
      })
    }))
    
    return(preds)
  }
  
  if (verbose) cat("Testing prediction wrapper...\n")
  test_pred <- predict_wrapper(model, background_matrix[1:min(5, nrow(background_matrix)), , drop = FALSE])
  if (verbose) cat("Test predictions:", paste(round(test_pred, 4), collapse = ", "), "\n")
  
  if (all(abs(test_pred - 0.5) < 0.001)) {
    warning("All test predictions are 0.5. The prediction wrapper may not be working correctly.")
  }
  
  if (method == "fastshap") {
    if (!requireNamespace("fastshap", quietly = TRUE)) {
      stop("fastshap package required")
    }
    
    if (verbose) {
      cat("Calculating SHAP values using fastshap", 
          if(parallel) paste("(", n_cores, " cores)") else "(sequential)", "...\n")
      cat("Estimated time:", 
          round(nsim * nrow(test_matrix) * ncol(test_matrix) / (if(parallel) n_cores * 10 else 100)), 
          "minutes\n")
    }
    
    shap_values <- fastshap::explain(
      object = model,
      X = background_matrix,
      pred_wrapper = predict_wrapper,
      nsim = nsim,
      newdata = test_matrix,
      adjust = TRUE,
      parallel = parallel
    )
    
  } else if (method == "python") {
    stop("Python method not fully implemented yet")
  } else {
    stop("Invalid method")
  }
  
  if (verbose) cat("SHAP calculation complete!\n")
  
  if (nrow(covariate_ref) > 0) {
    cov_names <- setNames(
      covariate_ref$covariateName,
      as.character(covariate_ref$covariateId)
    )
    
    matched_names <- cov_names[colnames(shap_values)]
    matched_names[is.na(matched_names)] <- colnames(shap_values)[is.na(matched_names)]
    colnames(shap_values) <- matched_names
  }
  
  if (verbose) cat("Calculating final predictions...\n")
  predictions <- predict_wrapper(model, test_matrix)
  
  if (verbose) {
    cat("Prediction summary:\n")
    cat("  Min:", round(min(predictions), 4), "\n")
    cat("  Mean:", round(mean(predictions), 4), "\n")
    cat("  Max:", round(max(predictions), 4), "\n")
  }
  
  result <- list(
    shap_values = shap_values,
    rowIds = indices,
    feature_names = colnames(shap_values),
    covariateIds = selected_covariates,
    method = method,
    background_size = nrow(background_matrix),
    predictions = predictions,
    test_matrix = test_matrix,
    background_matrix = background_matrix,
    parallel = parallel,
    n_cores = if(parallel) n_cores else NA
  )
  
  class(result) <- c("plpShap", "list")
  return(result)
}


#' Plot SHAP values and save to file (using ragg for headless servers)
#'
plot_plp_shap <- function(shap_result, 
                          plot_type = "importance",
                          observation = 1,
                          max_features = 20,
                          save_path,
                          width = 10,
                          height = 8,
                          dpi = 300) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required")
  }
  
  if (!requireNamespace("ragg", quietly = TRUE)) {
    stop("ragg package required for headless plotting. Install with: install.packages('ragg')")
  }
  
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
    
    p <- ggplot2::ggplot(
      importance_df,
      ggplot2::aes(x = reorder(feature, importance), y = importance)
    ) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Feature Importance (Mean |SHAP|)",
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
        title = paste("SHAP Values for Observation", observation),
        subtitle = paste("Prediction:", round(shap_result$predictions[observation], 4)),
        x = "Feature",
        y = "SHAP Value",
        fill = "Effect"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold")
      )
    
  } else if (plot_type == "beeswarm") {
    if (!requireNamespace("tidyr", quietly = TRUE)) {
      stop("tidyr package required for beeswarm plot")
    }
    
    shap_long <- as.data.frame(shap_values) |>
      tibble::rowid_to_column("observation") |>
      tidyr::pivot_longer(
        -observation,
        names_to = "feature",
        values_to = "shap_value"
      )
    
    importance <- shap_long |>
      dplyr::group_by(feature) |>
      dplyr::summarise(importance = mean(abs(shap_value)), .groups = "drop") |>
      dplyr::arrange(dplyr::desc(importance))
    
    top_features <- head(importance$feature, max_features)
    shap_long <- shap_long |>
      dplyr::filter(feature %in% top_features) |>
      dplyr::mutate(feature = factor(feature, levels = rev(top_features)))
    
    p <- ggplot2::ggplot(
      shap_long,
      ggplot2::aes(x = shap_value, y = feature)
    ) +
      ggplot2::geom_jitter(
        ggplot2::aes(color = shap_value),
        alpha = 0.6,
        height = 0.2,
        size = 1.5
      ) +
      ggplot2::scale_color_gradient2(
        low = "#1E88E5",
        mid = "white",
        high = "#FF0D57",
        midpoint = 0
      ) +
      ggplot2::labs(
        title = "SHAP Summary Plot",
        x = "SHAP Value",
        y = "Feature",
        color = "SHAP Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold")
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


#' Summary method for plpShap objects
#'
#' @param object A plpShap object
#'
summary.plpShap <- function(object, ...) {
  cat("SHAP Values Summary\n")
  cat("===================\n\n")
  cat("Method:", object$method, "\n")
  if (!is.na(object$parallel)) {
    cat("Parallel:", object$parallel, "\n")
    if (object$parallel && !is.na(object$n_cores)) {
      cat("Cores used:", object$n_cores, "\n")
    }
  }
  cat("Number of observations:", nrow(object$shap_values), "\n")
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

#' Export SHAP values to data frame
#'
#' @param shap_result Result from calculate_plp_shap()
#' @return A data frame with SHAP values in long format
#'
export_shap_values <- function(shap_result) {
  
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("tidyr package required")
  }
  
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
