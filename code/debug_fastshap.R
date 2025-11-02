# Diagnostic script to find out why predictions return 0.5
source("Helper.R")

library(PatientLevelPrediction)
library(dplyr)

plpData <- PatientLevelPrediction::loadPlpData("PlpMultiOutput/targetId_1_L1")
plpResult_gbm <- PatientLevelPrediction::loadPlpResult("PlpMultiOutput/Analysis_3/plpResult")

#' Final Working SHAP Calculation
#'
#' Based on diagnostic findings - predictions work, just need proper output handling
#'
calculate_plp_shap <- function(plpResult,
                                       plpData,
                                       n_samples = 100,
                                       background_size = 100,
                                       nsim = 50,
                                       parallel = FALSE,
                                       n_cores = NULL,
                                       verbose = TRUE) {
  
  library(PatientLevelPrediction)
  library(dplyr)
  library(fastshap)
  
  if (parallel && !is.null(n_cores)) {
    library(doParallel)
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }
  
  if (verbose) cat("=== SHAP Calculation (Xgboost) ===\n\n")
  
  model <- plpResult$model
  population <- plpResult$prediction
  covariateData <- plpData$covariateData
  
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
  
  # Get covariates
  selected_covariates <- unique(model$covariateImportance$covariateId)
  if (verbose) cat("Features:", length(selected_covariates), "\n\n")
  
  covariate_ref <- covariateData$covariateRef |> collect()
  original_metadata <- plpData$metaData
  
  # Create matrices
  create_matrix <- function(rowIds, covariateIds) {
    subset_cov <- covariateData$covariates |>
      filter(rowId %in% rowIds & covariateId %in% covariateIds) |>
      collect()
    
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
  
  # WORKING prediction wrapper - capture output but preserve return value
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
    
    # KEY FIX: Use capture.output but assign result FIRST
    pred <- suppressMessages(suppressWarnings(
      PatientLevelPrediction::predictPlp(
        plpModel = object,
        plpData = batch_plpData,
        population = batch_population
      )
    ))
    
    # Discard the verbose output that was printed
    invisible(capture.output(type = "output"))
    
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
    if (parallel) cat("Using", n_cores, "cores\n")
    cat("\n")
  }
  
  # Suppress the verbose output from fastshap iterations
  if (verbose) {
    # Show progress but suppress individual predictions
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
  } else {
    # Suppress everything
    shap_values <- suppressMessages(
      fastshap::explain(
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
    )
  }
  
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
    parallel = parallel,
    n_cores = if(parallel) n_cores else NA
  )
  
  class(result) <- c("plpShap", "list")
  return(result)
}

shap_result_gbm <- calculate_plp_shap_working(
  plpResult = plpResult_gbm,
  plpData = plpData,
  n_samples = 100,
  background_size = 100,
  nsim = 50,
  parallel = FALSE,
  n_cores = 10,
  verbose = TRUE
)

# Check results
summary(shap_result_gbm)

# Save if successful
na_pct <- 100 * sum(is.na(shap_result_gbm$shap_values)) / length(shap_result_gbm$shap_values)
cat("\nNA percentage:", round(na_pct, 2), "%\n")

if (na_pct < 10) {
  saveRDS(shap_result_gbm, "shap_result_gbm_FINAL.rds")
  
  # Create plots
  plot_plp_shap(
    shap_result_gbm,
    plot_type = "importance",
    save_path = "shap_importance_FINAL.png",
    width = 14,
    height = 10,
    dpi = 600
  )
  
  plot_plp_shap(
    shap_result_gbm,
    plot_type = "beeswarm",
    max_features = 15,
    save_path = "shap_beeswarm_FINAL.png",
    width = 15,
    height = 8,
    dpi = 600
  )
  
  cat("\n✓ Success! Plots saved.\n")
}
