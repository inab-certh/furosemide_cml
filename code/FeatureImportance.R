# source("helper.R")
# 
# library(PatientLevelPrediction)
# 
# library(doParallel)
# library(foreach)
# 
# plpData <- PatientLevelPrediction::loadPlpData("PlpMultiOutput/targetId_1_L1")
# plpResult <- PatientLevelPrediction::loadPlpResult("PlpMultiOutput/Analysis_3/plpResult")
# 
# shap_result <- calculate_plp_shap(
#   plpResult = plpResult,
#   plpData = plpData,
#   population = plpResult$prediction,
#   background_size = 100,  # Number of background samples
#   method = "fastshap",     # or "python" if you have Python SHAP installed
#   nsim = 50                # Number of simulations
# )
# 
# saveRDS(shap_result, "shap_result.rds")
# 
# summary(shap_result)
# 
# plot_plp_shap(
#   shap_result, 
#   plot_type = "importance",
#   save_path = "shap_importance.png",
#   width = 14,
#   height = 10,
#   dpi = 600  # High resolution
# )
# 
# plot_plp_shap(
#   shap_result, 
#   plot_type = "beeswarm",
#   max_features = 15,
#   save_path = "shap_beeswarm.png",
#   width = 15,
#   height = 8,
#   dpi = 600  # High resolution
# )
