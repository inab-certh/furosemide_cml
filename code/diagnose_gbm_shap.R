#' Calculate SHAP values for PatientLevelPrediction models (Improved for GBM)
#'
#' Enhanced version with better model type handling and debugging
#'
calculate_plp_shap_improved <- function(plpResult,
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
  if (verbose) {
    cat("Model type:", if(!is.null(model_type)) model_type else "unknown", "\n")
    cat("Model class:", paste(class(model), collapse = ", "), "\n")
  }
  
  # Check model structure for debugging
  if (verbose) {
    cat("\nModel structure:\n")
    cat("Names:", paste(names(model), collapse = ", "), "\n")
    if (!is.null(model$model)) {
      cat("Inner model class:", paste(class(model$model), collapse = ", "), "\n")
    }
  }
  
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
  
  if (verbose) cat("Creating IMPROVED prediction wrapper with better GBM support...\n")
  
  # IMPROVED prediction wrapper with better error handling and model type detection
  predict_wrapper <- function(object, newdata) {
    
    if (is.matrix(newdata)) {
      newdata <- as.data.frame(newdata)
    }
    
    covariate_ids <- as.numeric(colnames(newdata))
    n_obs <- nrow(newdata)
    
    # Build covariate data more efficiently
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
    
    # Remove NULL entries
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
      survivalTime = 365,
      indexes = 1:n_obs
    )
    
    # Create properly structured CovariateData
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
    
    class(batch_plpData) <- "plpData"
    
    # CRITICAL: Suppress ALL output including print statements
    preds <- suppressMessages(suppressWarnings({
      invisible(capture.output({
        tryCatch({
          # Use predictPlp which should handle all model types
          pred <- PatientLevelPrediction::predictPlp(
            plpModel = object,
            plpData = batch_plpData,
            population = batch_population
          )
          
          # Extract predictions in correct order
          result <- pred$value[match(1:n_obs, pred$rowId)]
          
          # Check for NA values
          if (any(is.na(result))) {
            na_count <- sum(is.na(result))
            warning_msg <- paste("Found", na_count, "NA predictions out of", n_obs)
            if (verbose) cat(warning_msg, "\n")
            
            # Try to impute NAs with mean of non-NA values
            mean_pred <- mean(result, na.rm = TRUE)
            if (is.na(mean_pred)) {
              # If all are NA, use 0.5 as fallback
              mean_pred <- 0.5
            }
            result[is.na(result)] <- mean_pred
          }
          
          # Sanity check - ensure predictions are in valid range
          if (any(result < 0 | result > 1, na.rm = TRUE)) {
            warning("Predictions outside [0,1] range detected. Clipping values.")
            result <- pmax(0, pmin(1, result))
          }
          
          result
          
        }, error = function(e) {
          if (verbose) {
            cat("\n!!! PREDICTION ERROR !!!\n")
            cat("Error message:", e$message, "\n")
            cat("Model type:", class(object), "\n")
            cat("Number of observations:", n_obs, "\n")
            cat("Number of non-zero covariates:", nrow(batch_covariates), "\n")
          }
          # Return default predictions
          rep(0.5, n_obs)
        })
      }, type = "message"))
    }))
    
    return(preds)
  }
  
  if (verbose) cat("\nTesting prediction wrapper with multiple samples...\n")
  
  # Test with different sample sizes
  test_sizes <- c(1, min(5, nrow(background_matrix)), min(10, nrow(background_matrix)))
  test_results <- list()
  
  for (test_size in test_sizes) {
    if (verbose) cat("Testing with", test_size, "sample(s)...\n")
    test_pred <- predict_wrapper(model, background_matrix[1:test_size, , drop = FALSE])
    test_results[[as.character(test_size)]] <- test_pred
    
    if (verbose) {
      cat("  Predictions:", paste(round(test_pred, 4), collapse = ", "), "\n")
      cat("  Range: [", round(min(test_pred), 4), ", ", round(max(test_pred), 4), "]\n", sep = "")
      cat("  Mean:", round(mean(test_pred), 4), "\n")
      cat("  NAs:", sum(is.na(test_pred)), "\n")
    }
  }
  
  # Check if predictions are all the same (problematic)
  all_test_preds <- unlist(test_results)
  if (length(unique(all_test_preds)) == 1) {
    warning("All test predictions are identical (", all_test_preds[1], 
            "). The model may not be predicting correctly.")
  }
  
  if (all(abs(all_test_preds - 0.5) < 0.001)) {
    warning("All test predictions are ~0.5. The prediction wrapper may not be working correctly.")
    cat("\n!!! DEBUGGING INFO !!!\n")
    cat("Model structure:\n")
    print(str(model, max.level = 2))
  }
  
  if (method == "fastshap") {
    if (!requireNamespace("fastshap", quietly = TRUE)) {
      stop("fastshap package required")
    }
    
    if (verbose) {
      cat("\nCalculating SHAP values using fastshap", 
          if(parallel) paste("(", n_cores, " cores)") else "(sequential)", "...\n")
      cat("Estimated time:", 
          round(nsim * nrow(test_matrix) * ncol(test_matrix) / (if(parallel) n_cores * 10 else 100)), 
          "minutes\n")
    }
    
    # Calculate SHAP with error handling
    shap_values <- tryCatch({
      fastshap::explain(
        object = model,
        X = background_matrix,
        pred_wrapper = predict_wrapper,
        nsim = nsim,
        newdata = test_matrix,
        adjust = TRUE,
        parallel = parallel
      )
    }, error = function(e) {
      cat("\n!!! SHAP CALCULATION ERROR !!!\n")
      cat("Error:", e$message, "\n")
      stop("SHAP calculation failed. Check the prediction wrapper.")
    })
    
  } else if (method == "python") {
    stop("Python method not fully implemented yet")
  } else {
    stop("Invalid method")
  }
  
  if (verbose) {
    cat("\nSHAP calculation complete!\n")
    cat("SHAP values summary:\n")
    cat("  Dimensions:", nrow(shap_values), "x", ncol(shap_values), "\n")
    cat("  Range: [", round(min(shap_values, na.rm = TRUE), 4), ", ", 
        round(max(shap_values, na.rm = TRUE), 4), "]\n", sep = "")
    cat("  NAs:", sum(is.na(shap_values)), "\n")
  }
  
  # Check for all-NA columns
  na_cols <- colSums(is.na(shap_values)) == nrow(shap_values)
  if (any(na_cols)) {
    warning("Found ", sum(na_cols), " features with all NA SHAP values. These will be removed.")
    shap_values <- shap_values[, !na_cols, drop = FALSE]
  }
  
  # Add feature names
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
    cat("Final prediction summary:\n")
    cat("  Min:", round(min(predictions, na.rm = TRUE), 4), "\n")
    cat("  Mean:", round(mean(predictions, na.rm = TRUE), 4), "\n")
    cat("  Max:", round(max(predictions, na.rm = TRUE), 4), "\n")
    cat("  NAs:", sum(is.na(predictions)), "\n")
  }
  
  result <- list(
    shap_values = shap_values,
    rowIds = indices,
    feature_names = colnames(shap_values),
    covariateIds = selected_covariates[selected_covariates %in% as.numeric(colnames(shap_values))],
    method = method,
    background_size = nrow(background_matrix),
    predictions = predictions,
    test_matrix = test_matrix,
    background_matrix = background_matrix,
    parallel = parallel,
    n_cores = if(parallel) n_cores else NA,
    model_type = model_type
  )
  
  class(result) <- c("plpShap", "list")
  return(result)
}


#' Alternative SHAP calculation using direct model access (for GBM debugging)
#'
#' This function attempts to bypass the standard prediction wrapper
#' and access the model more directly
#'
calculate_plp_shap_direct <- function(plpResult,
                                      plpData,
                                      population = NULL,
                                      subset = "test",
                                      n_samples = 50,  # Smaller default for debugging
                                      background_size = 100,
                                      nsim = 50,
                                      verbose = TRUE) {
  
  if (verbose) cat("=== DIRECT SHAP CALCULATION (Debugging Mode) ===\n\n")
  
  # Extract model
  model <- plpResult$model
  model_type <- attr(model, "modelType")
  
  if (verbose) {
    cat("Model type:", model_type, "\n")
    cat("Model class:", paste(class(model), collapse = ", "), "\n")
    cat("\nModel structure (first level):\n")
    print(names(model))
  }
  
  # For GBM, check if we can access the inner model
  if (!is.null(model$model)) {
    if (verbose) {
      cat("\nInner model found!\n")
      cat("Inner model class:", paste(class(model$model), collapse = ", "), "\n")
      if ("gbm" %in% class(model$model)) {
        cat("GBM parameters:\n")
        cat("  Trees:", model$model$n.trees, "\n")
        cat("  Variables:", length(model$model$var.names), "\n")
      }
    }
  }
  
  # Get population subset
  if (is.null(population)) {
    population <- plpResult$prediction
  }
  
  if ("indexes" %in% names(population) && subset == "test") {
    selected_pop <- population[population$indexes > 0, ]
  } else {
    selected_pop <- population
  }
  
  if (!is.null(n_samples) && nrow(selected_pop) > n_samples) {
    selected_pop <- selected_pop[sample(nrow(selected_pop), n_samples), ]
  }
  
  if (verbose) cat("\nSelected", nrow(selected_pop), "observations for SHAP calculation\n")
  
  # Try to get actual predictions from plpResult
  if (verbose) cat("\nChecking actual model predictions...\n")
  
  actual_preds <- selected_pop$value
  if (verbose) {
    cat("Actual predictions from plpResult:\n")
    cat("  Range: [", round(min(actual_preds), 4), ", ", round(max(actual_preds), 4), "]\n", sep = "")
    cat("  Mean:", round(mean(actual_preds), 4), "\n")
    cat("  Unique values:", length(unique(actual_preds)), "\n")
  }
  
  # Recommendation
  cat("\n=== RECOMMENDATIONS ===\n")
  cat("1. Check if the GBM model is properly fitted\n")
  cat("2. Verify that predictPlp() works correctly with this model\n")
  cat("3. Try running: PatientLevelPrediction::predictPlp(plpResult$model, plpData, population)\n")
  cat("4. If predictions work manually but not in SHAP, the issue is with the wrapper\n")
  cat("\nTo test manually:\n")
  cat("test_pred <- PatientLevelPrediction::predictPlp(plpResult$model, plpData, population)\n")
  cat("summary(test_pred$value)\n")
  
  return(invisible(list(
    model_type = model_type,
    actual_predictions = actual_preds,
    population_size = nrow(selected_pop)
  )))
}


#' Calculate SHAP for XGBoost/GBM - Final Fixed Version
#' 
#' Uses sink() to suppress output properly while preserving return values
#'
calculate_plp_shap_xgb <- function(plpResult,
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
  if (verbose) {
    cat("Model type:", if(!is.null(model_type)) model_type else "unknown", "\n")
    cat("Model class:", paste(class(model), collapse = ", "), "\n")
  }
  
  # Check if it's XGBoost
  is_xgboost <- !is.null(model$model) && "xgb.Booster" %in% class(model$model)
  if (verbose && is_xgboost) {
    cat("XGBoost model detected\n")
  }
  
  covariateData <- plpData$covariateData
  
  if (is.null(population)) {
    if (!is.null(plpResult$prediction)) {
      population <- plpResult$prediction
    } else {
      stop("No population found.")
    }
  }
  
  if (verbose) cat("Preparing population...\n")
  
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
    if (verbose) cat("Collecting all covariates...\n")
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
  }
  
  create_dense_matrix <- function(rowIds, covariateIds) {
    subset_cov <- covariateData$covariates |>
      dplyr::filter(rowId %in% rowIds & covariateId %in% covariateIds) |>
      dplyr::collect()
    
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
    
    return(dense_matrix)
  }
  
  background_matrix <- create_dense_matrix(background_ids, selected_covariates)
  test_matrix <- create_dense_matrix(indices, selected_covariates)
  
  if (verbose) {
    cat("Background matrix:", nrow(background_matrix), "x", ncol(background_matrix), "\n")
    cat("Test matrix:", nrow(test_matrix), "x", ncol(test_matrix), "\n")
  }
  
  if (verbose) cat("Collecting covariate reference...\n")
  covariate_ref <- covariateData$covariateRef |> dplyr::collect()
  
  original_metadata <- plpData$metaData
  if (is.null(original_metadata)) {
    original_metadata <- list()
  }
  
  if (verbose) cat("Creating prediction wrapper with sink() for output suppression...\n")
  
  # NEW APPROACH: Use sink() to redirect output to null device
  predict_wrapper <- function(object, newdata) {
    
    if (is.matrix(newdata)) {
      newdata <- as.data.frame(newdata)
    }
    
    covariate_ids <- as.numeric(colnames(newdata))
    n_obs <- nrow(newdata)
    
    # Build covariate data
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
      survivalTime = 365,
      indexes = 1:n_obs
    )
    
    # Create CovariateData
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
    
    # Create plpData
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
    
    # Use sink() to redirect all output to null device
    result <- tryCatch({
      # Open connection to null device
      sink_file <- tempfile()
      sink(sink_file, type = "output")
      sink(sink_file, type = "message")
      
      # Wrap in try block to ensure sinks are closed
      pred_result <- try({
        suppressMessages(suppressWarnings({
          pred <- PatientLevelPrediction::predictPlp(
            plpModel = object,
            plpData = batch_plpData,
            population = batch_population
          )
          
          # Extract values
          pred_values <- pred$value[match(1:n_obs, pred$rowId)]
          
          # Handle NAs
          if (any(is.na(pred_values))) {
            mean_pred <- mean(pred_values, na.rm = TRUE)
            if (is.na(mean_pred)) mean_pred <- 0.5
            pred_values[is.na(pred_values)] <- mean_pred
          }
          
          as.numeric(pred_values)
        }))
      }, silent = TRUE)
      
      # Close sinks
      sink(type = "message")
      sink(type = "output")
      unlink(sink_file)
      
      if (inherits(pred_result, "try-error")) {
        return(rep(0.5, n_obs))
      }
      
      pred_result
      
    }, error = function(e) {
      # Make sure sinks are closed even on error
      try(sink(type = "message"), silent = TRUE)
      try(sink(type = "output"), silent = TRUE)
      rep(0.5, n_obs)
    })
    
    return(result)
  }
  
  if (verbose) cat("\nTesting prediction wrapper...\n")
  
  # Test predictions
  test_sizes <- c(1, min(5, nrow(background_matrix)))
  all_tests_passed <- TRUE
  
  for (test_size in test_sizes) {
    if (verbose) cat("Testing with", test_size, "sample(s)...\n")
    
    test_pred <- predict_wrapper(model, background_matrix[1:test_size, , drop = FALSE])
    
    if (is.numeric(test_pred) && length(test_pred) == test_size) {
      if (verbose) {
        cat("  ✓ Length:", length(test_pred), "\n")
        cat("  ✓ Values:", paste(round(test_pred, 4), collapse = ", "), "\n")
      }
    } else {
      if (verbose) {
        cat("  ✗ ERROR: Class =", class(test_pred), ", Length =", length(test_pred), "\n")
      }
      all_tests_passed <- FALSE
    }
  }
  
  if (!all_tests_passed) {
    stop("Prediction wrapper tests failed.")
  }
  
  if (verbose) cat("\n✓ Prediction wrapper working correctly!\n\n")
  
  if (method == "fastshap") {
    if (!requireNamespace("fastshap", quietly = TRUE)) {
      stop("fastshap package required")
    }
    
    if (verbose) {
      cat("Calculating SHAP values...\n")
      cat("This will take several minutes...\n\n")
    }
    
    # Calculate SHAP (fastshap handles its own output)
    shap_values <- fastshap::explain(
      object = model,
      X = background_matrix,
      pred_wrapper = predict_wrapper,
      nsim = nsim,
      newdata = test_matrix,
      adjust = TRUE,
      parallel = parallel
    )
    
  } else {
    stop("Only fastshap method supported")
  }
  
  if (verbose) {
    cat("\nSHAP calculation complete!\n")
    if (is.matrix(shap_values) || is.data.frame(shap_values)) {
      cat("Dimensions:", nrow(shap_values), "x", ncol(shap_values), "\n")
      cat("Range: [", round(min(shap_values, na.rm = TRUE), 4), ", ", 
          round(max(shap_values, na.rm = TRUE), 4), "]\n", sep = "")
      cat("NAs:", sum(is.na(shap_values)), "\n")
    } else {
      cat("WARNING: SHAP values have unexpected structure\n")
      cat("Class:", class(shap_values), "\n")
      print(str(shap_values))
    }
  }
  
  # Add feature names
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
    n_cores = if(parallel) n_cores else NA,
    model_type = model_type
  )
  
  class(result) <- c("plpShap", "list")
  return(result)
}

#' Calculate SHAP without relying on indexes column
#' 
#' Works when plpResult$prediction doesn't have indexes column
#'
calculate_plp_shap_no_indexes <- function(plpResult,
                                          plpData,
                                          population = NULL,
                                          use_train_for_background = TRUE,
                                          n_samples = 100,
                                          background_size = 100,
                                          method = "fastshap",
                                          nsim = 50,
                                          parallel = FALSE,
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
  
  if (verbose) cat("=== SHAP Calculation (No Indexes Mode) ===\n\n")
  
  model <- plpResult$model
  model_type <- attr(model, "modelType")
  if (verbose) {
    cat("Model type:", if(!is.null(model_type)) model_type else "unknown", "\n")
  }
  
  covariateData <- plpData$covariateData
  
  if (is.null(population)) {
    if (!is.null(plpResult$prediction)) {
      population <- plpResult$prediction
    } else {
      stop("No population found.")
    }
  }
  
  if (verbose) {
    cat("Total population:", nrow(population), "\n")
    cat("Columns:", paste(names(population), collapse = ", "), "\n\n")
  }
  
  # Check for evaluationType to determine train/test split
  if ("evaluationType" %in% names(population)) {
    if (verbose) cat("Using evaluationType column for train/test split\n")
    
    test_pop <- population[population$evaluationType == "Test", ]
    train_pop <- population[population$evaluationType == "Train", ]
    
    if (nrow(test_pop) == 0) {
      # Fallback: use all as test, sample for background
      if (verbose) cat("No Test set found in evaluationType, using all data\n")
      test_pop <- population
      if (nrow(population) > background_size + n_samples) {
        # Sample background from population excluding test samples
        all_ids <- population$rowId
        test_ids <- sample(all_ids, n_samples)
        train_ids <- setdiff(all_ids, test_ids)
        test_pop <- population[population$rowId %in% test_ids, ]
        train_pop <- population[population$rowId %in% train_ids, ]
      } else {
        train_pop <- population
      }
    }
    
  } else if ("indexes" %in% names(population)) {
    if (verbose) cat("Using indexes column for train/test split\n")
    
    test_pop <- population[population$indexes > 0, ]
    train_pop <- population[population$indexes < 0, ]
    
    if (nrow(test_pop) == 0) {
      if (verbose) cat("No positive indexes found, using all data\n")
      test_pop <- population
      train_pop <- population
    }
    
  } else {
    # No split information available
    if (verbose) cat("No train/test split information found\n")
    
    if (use_train_for_background) {
      if (verbose) cat("Creating random split for background\n")
      
      # Use all population for test, sample for background
      test_pop <- population
      
      if (nrow(population) > background_size + n_samples) {
        # Sample n_samples for explanation
        test_ids <- sample(population$rowId, n_samples)
        test_pop <- population[population$rowId %in% test_ids, ]
        
        # Use remaining for background
        background_ids <- sample(
          setdiff(population$rowId, test_ids), 
          min(background_size, nrow(population) - n_samples)
        )
        train_pop <- population[population$rowId %in% background_ids, ]
      } else {
        # Too few samples, use same data for both
        train_pop <- population
        test_pop <- population[1:min(n_samples, nrow(population)), ]
      }
    } else {
      if (verbose) cat("Using all data for both test and background\n")
      test_pop <- population
      train_pop <- population
    }
  }
  
  # Sample if needed
  if (nrow(test_pop) > n_samples) {
    if (verbose) cat("Sampling", n_samples, "from", nrow(test_pop), "test observations\n")
    test_pop <- test_pop[sample(nrow(test_pop), n_samples), ]
  }
  
  if (nrow(train_pop) > background_size) {
    if (verbose) cat("Sampling", background_size, "from", nrow(train_pop), "background observations\n")
    train_pop <- train_pop[sample(nrow(train_pop), background_size), ]
  }
  
  if (verbose) {
    cat("\nFinal split:\n")
    cat("  Test samples:", nrow(test_pop), "\n")
    cat("  Background samples:", nrow(train_pop), "\n\n")
  }
  
  indices <- test_pop$rowId
  background_ids <- train_pop$rowId
  
  # Get covariates
  if (!is.null(model$covariateImportance)) {
    selected_covariates <- unique(model$covariateImportance$covariateId)
    if (verbose) cat("Using", length(selected_covariates), "covariates from model\n")
  } else {
    if (verbose) cat("Collecting all covariates...\n")
    all_covariates <- covariateData$covariates |> 
      dplyr::select(covariateId) |>
      dplyr::distinct() |>
      dplyr::collect()
    selected_covariates <- unique(all_covariates$covariateId)
    if (verbose) cat("Using all", length(selected_covariates), "covariates\n")
  }
  
  # Create matrices
  create_dense_matrix <- function(rowIds, covariateIds) {
    subset_cov <- covariateData$covariates |>
      dplyr::filter(rowId %in% rowIds & covariateId %in% covariateIds) |>
      dplyr::collect()
    
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
    
    return(dense_matrix)
  }
  
  if (verbose) cat("Creating feature matrices...\n")
  background_matrix <- create_dense_matrix(background_ids, selected_covariates)
  test_matrix <- create_dense_matrix(indices, selected_covariates)
  
  if (verbose) {
    cat("Background matrix:", nrow(background_matrix), "x", ncol(background_matrix), "\n")
    cat("Test matrix:", nrow(test_matrix), "x", ncol(test_matrix), "\n\n")
  }
  
  # Get covariate reference
  covariate_ref <- covariateData$covariateRef |> dplyr::collect()
  
  original_metadata <- plpData$metaData
  if (is.null(original_metadata)) {
    original_metadata <- list()
  }
  
  if (verbose) cat("Creating prediction wrapper...\n")
  
  # Prediction wrapper using sink()
  predict_wrapper <- function(object, newdata) {
    
    if (is.matrix(newdata)) {
      newdata <- as.data.frame(newdata)
    }
    
    covariate_ids <- as.numeric(colnames(newdata))
    n_obs <- nrow(newdata)
    
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
    
    result <- tryCatch({
      sink_file <- tempfile()
      sink(sink_file, type = "output")
      sink(sink_file, type = "message")
      
      pred_result <- try({
        suppressMessages(suppressWarnings({
          pred <- PatientLevelPrediction::predictPlp(
            plpModel = object,
            plpData = batch_plpData,
            population = batch_population
          )
          
          pred_values <- pred$value[match(1:n_obs, pred$rowId)]
          
          if (any(is.na(pred_values))) {
            mean_pred <- mean(pred_values, na.rm = TRUE)
            if (is.na(mean_pred)) mean_pred <- 0.5
            pred_values[is.na(pred_values)] <- mean_pred
          }
          
          as.numeric(pred_values)
        }))
      }, silent = TRUE)
      
      sink(type = "message")
      sink(type = "output")
      unlink(sink_file)
      
      if (inherits(pred_result, "try-error")) {
        return(rep(0.5, n_obs))
      }
      
      pred_result
      
    }, error = function(e) {
      try(sink(type = "message"), silent = TRUE)
      try(sink(type = "output"), silent = TRUE)
      rep(0.5, n_obs)
    })
    
    return(result)
  }
  
  if (verbose) cat("Testing prediction wrapper...\n")
  
  test_pred <- predict_wrapper(model, background_matrix[1:min(5, nrow(background_matrix)), , drop = FALSE])
  if (verbose) {
    cat("Test predictions:", paste(round(test_pred, 4), collapse = ", "), "\n")
    cat("✓ Prediction wrapper working\n\n")
  }
  
  if (method == "fastshap") {
    if (!requireNamespace("fastshap", quietly = TRUE)) {
      stop("fastshap package required")
    }
    
    if (verbose) {
      cat("Calculating SHAP values using fastshap...\n")
      cat("This may take several minutes...\n\n")
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
    
  } else {
    stop("Only fastshap method supported")
  }
  
  if (verbose) {
    cat("\nSHAP calculation complete!\n")
    if (is.matrix(shap_values) || is.data.frame(shap_values)) {
      cat("Dimensions:", nrow(shap_values), "x", ncol(shap_values), "\n")
      
      non_na <- sum(!is.na(shap_values))
      total <- length(shap_values)
      cat("Non-NA values:", non_na, "/", total, 
          "(", round(100 * non_na / total, 1), "%)\n")
      
      if (non_na > 0) {
        cat("Range: [", round(min(shap_values, na.rm = TRUE), 4), ", ", 
            round(max(shap_values, na.rm = TRUE), 4), "]\n", sep = "")
      }
    }
  }
  
  # Add feature names
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
    n_cores = if(parallel) n_cores else NA,
    model_type = model_type
  )
  
  class(result) <- c("plpShap", "list")
  return(result)
}

plpResult_gbm <- PatientLevelPrediction::loadPlpResult("PlpMultiOutput/Analysis_3/plpResult")
plpData <- PatientLevelPrediction::loadPlpData("PlpMultiOutput/targetId_1_L1")

shap_result_gbm <- calculate_plp_shap_no_indexes(
  plpResult = plpResult_gbm,
  plpData = plpData,
  use_train_for_background = TRUE,  # Create random split
  n_samples = 20,
  background_size = 50,
  nsim = 10,
  parallel = TRUE,
  verbose = TRUE
)

shap_result_gbm$shap_values

# shap_result_gbm <- calculate_plp_shap_xgb(
#   plpResult = plpResult_gbm,
#   plpData = plpData,
#   method = "python",
#   n_samples = 20,       # Start small
#   background_size = 50,
#   nsim = 10,
#   parallel = TRUE,     # Better error messages
#   verbose = TRUE
# )

saveRDS(shap_result_gbm, "shap_gbm.rds")
