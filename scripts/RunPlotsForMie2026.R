source("code/Helper.R")

shap_result_rf_full <- readRDS(
  "results_full/evaluation/shap_result_analysis_2.rds"
)
shap_result_rf_no_ckd <- readRDS(
  "results_no_ckd/evaluation/shap_result_analysis_2.rds"
)

plot_plp_shap(
  shap_result = shap_result_rf_full,
  plot_type = "beeswarm",
  save_path = "shap_beeswarm_rf_full.png",
  title = "SHAP Summary Plot: Full population",
  max_features = 10,
  width = 14,
  dpi = 1200
)

plot_plp_shap(
  shap_result = shap_result_rf_no_ckd,
  plot_type = "beeswarm",
  save_path = "shap_beeswarm_rf_no_ckd.png",
  title = "SHAP Summary Plot: Population without CKD",
  max_features = 10,
  width = 14,
  dpi = 1200
)
