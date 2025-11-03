source("code/Helper.R")

library(PatientLevelPrediction)
library(dplyr)

allDirs <- list.dirs(
  path = file.path(resultsDir, "models"),
  recursive = FALSE,
  full.names = TRUE
)

evaluationDir <- file.path(resultsDir, "evaluation")
if (!dir.exists(evaluationDir)) dir.create(evaluationDir)

modelDirs <- allDirs[stringr::str_detect(allDirs, "Analysis")]
plpDataDir <- allDirs[stringr::str_detect(allDirs, "target")]

# if (length(plpDataDir >= 2)) {
#   stop("Only one plpData directory at a time can be processed")
# }

plpData <- PatientLevelPrediction::loadPlpData(plpDataDir)

for (i in seq_along(modelDirs)) {
  analysisName <- basename(modelDirs[i]) |> tolower()
  message(glue::glue("Running model explanation for: { analysisName }..."))
  plpResult <- PatientLevelPrediction::loadPlpResult(
    file.path(modelDirs[i], "plpResult")
  )
  
  message("Calculating  SHAP values. This may take a while...")
  shap_result <- calculate_plp_shap(
    plpResult = plpResult,
    plpData = plpData,
    n_samples = 100,
    parallel = TRUE
  )

  message("Finished calculating SHAP values.")

  saveRDS(
    shap_result,
    file = file.path(
      evaluationDir,
      glue::glue("shap_result_{ analysisName }.rds")
    )
  )

  message("Saved the shap values calculations.")

  ranks <- get_feature_ranks(shap_result)
  message("==== SHAP value ranks ====")
  ranks |> dplyr::filter(times_top_10 > 10)
  readr::write_csv(
    ranks,
    file.path(
      evaluationDir,
      glue::glue("shap_ranks_{ analysisName }.csv")
    )
  )
  message("Saved the shap values ranks.")

  plot_plp_shap(
    shap_result,
    "beeswarm",
    save_path = file.path(
      evaluationDir,
      glue::glue("shap_beeswarm_{ analysisName }.png")
    ),
    dpi = 600
  )
  message("Saved the beeswarm plot")
}
