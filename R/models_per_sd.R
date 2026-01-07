#' Models per PRS Standard Deviation
#' @description Functions to obtain models with PRS + covariates, with only covariates or only PRS
#' @param train_ctrl_stats Training set of the PRS data.
#' @param test_ctrl_stats Testing set of the PRS data.
#' @param log_reg A parsnip logistic regression model.
#' @param norm_prs name of the normalized PRS column.
#' @param recipe_var A recipe from recipes::recipe() object, to be provided optionally, in case the recipe originally included in the model (status ~ prs + scaled_centered_pc(1:10) + age_analysis) is not what you need for your PRS.
#'
#' @return A list with OR, AUC, delta AUC, and ROC curve for AUC with and without PRS.
#' @importFrom generics fit
#' @importFrom stats predict
#' @importFrom stats as.formula
#' @importFrom broom tidy
#' @importFrom recipes recipe
#' @importFrom recipes step_rm
#' @importFrom recipes step_normalize
#' @import rsample
#' @import dplyr
#' @import tidyr
#' @import parsnip
#' @import workflows
#' @import yardstick
#' @import pROC
#' @import tune
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' # TBD

# Modeling with PRS function -------------------------------
model_with_prs <- function(
  train_ctrl_stats,
  test_ctrl_stats,
  log_reg,
  norm_prs,
  recipe_var = NULL
) {
  # If external recipe is not provided, use the basic one with PCs 1:10 and age_analysis as covariates
  if (is.null(recipe_var)) {
    form_full_model <- as.formula(paste(
      "status ~",
      norm_prs,
      "+",
      paste0("pc", 1:10, sep = " + ", collapse = " "),
      "age_analysis"
    ))

    rec_prs <- recipes::recipe(
      form_full_model,
      data = train_ctrl_stats
    ) |>
      recipes::step_normalize(dplyr::starts_with("pc"))
  } else {
    rec_prs <- recipe_var
  }

  wflow_prs <- workflows::workflow() |>
    workflows::add_model(log_reg) |>
    workflows::add_recipe(rec_prs)

  # Odds Ratio
  fit_prs <- wflow_prs |>
    generics::fit(data = train_ctrl_stats)

  or <- broom::tidy(fit_prs, exponentiate = T, conf.int = T) |>
    dplyr::select(!dplyr::all_of(c("std.error", "statistic")))

  # AUC training set
  data_folds <- rsample::vfold_cv(train_ctrl_stats, v = 10, strata = status)

  fit_resamples_full_model <- tune::fit_resamples(
    wflow_prs,
    resamples = data_folds,
    control = tune::control_resamples(save_pred = TRUE, save_workflow = TRUE)
  )

  auc_training <- fit_resamples_full_model |>
    tune::collect_metrics(summarize = F) |>
    dplyr::filter(.metric == "roc_auc") |>
    dplyr::select(id, .metric, .estimate)

  # AUC testing set
  pred_prob_prs_vars <- predict(
    fit_prs,
    new_data = test_ctrl_stats,
    type = "prob"
  )

  test_results_with_prs <- test_ctrl_stats |>
    dplyr::bind_cols(pred_prob_prs_vars)

  roc_auc_test <- pROC::roc(
    response = test_results_with_prs$status,
    predictor = test_results_with_prs$.pred_0
  )

  ci_delong <- pROC::ci.auc(roc_auc_test, method = "delong")

  df_auc_test <- data.frame(
    auc = roc_auc_test$auc[[1]],
    lower = ci_delong[[1]],
    upper = ci_delong[[3]]
  )

  # Get data for ROC curve
  roc_curve_test_data <- data.frame(
    model = "with_prs",
    sensitivity = roc_auc_test$sensitivities,
    specificity = roc_auc_test$specificities
  )

  res <- list(
    or = or,
    auc_train_with_prs = auc_training,
    auc_test_with_prs = df_auc_test,
    roc_curve_test_with_prs_data = roc_curve_test_data,
    full_roc_test_with_prs = roc_auc_test,
    fit_with_prs = fit_prs,
    pred_with_prs = test_results_with_prs
  )

  return(res)
}


# Modeling without PRS function ---------------------------------
model_without_prs <- function(
  train_ctrl_stats,
  test_ctrl_stats,
  log_reg,
  norm_prs,
  recipe_var
) {
  # If external recipe is not provided, use the basic one with PCs 1:10 and age_analysis as covariates
  if (is.null(recipe_var)) {
    # Formula with predictors wanted
    form_model_wo_prs <- as.formula(paste(
      "status ~",
      paste0("pc", 1:10, sep = " + ", collapse = " "),
      "age_analysis"
    ))

    rec_wo_prs <- recipes::recipe(
      form_model_wo_prs,
      data = train_ctrl_stats
    ) |>
      recipes::step_normalize(dplyr::starts_with("pc"))
  } else {
    rec_wo_prs <- recipe_var |>
      step_rm(norm_prs)
  }

  wflow_wo_prs <- workflows::workflow() |>
    workflows::add_model(log_reg) |>
    workflows::add_recipe(rec_wo_prs)

  # Fit for test AUC
  fit_wo_prs <- wflow_wo_prs |>
    generics::fit(data = train_ctrl_stats)

  # AUC training set
  data_folds <- rsample::vfold_cv(train_ctrl_stats, v = 10, strata = status)

  fit_resamples_model_wo_prs <- tune::fit_resamples(
    wflow_wo_prs,
    resamples = data_folds,
    control = tune::control_resamples(save_pred = TRUE, save_workflow = TRUE)
  )

  auc_training <- fit_resamples_model_wo_prs |>
    tune::collect_metrics(summarize = F) |>
    dplyr::filter(.metric == "roc_auc") |>
    dplyr::select(id, .metric, .estimate)

  # AUC testing set
  pred_prob_wo_prs <- predict(
    fit_wo_prs,
    new_data = test_ctrl_stats,
    type = "prob"
  )

  test_results_wo_prs <- test_ctrl_stats |>
    dplyr::bind_cols(pred_prob_wo_prs)

  roc_auc_test <- pROC::roc(
    response = test_results_wo_prs$status,
    predictor = test_results_wo_prs$.pred_0
  )

  ci_delong <- pROC::ci.auc(roc_auc_test, method = "delong")

  df_auc_test <- data.frame(
    auc = roc_auc_test$auc[[1]],
    lower = ci_delong[[1]],
    upper = ci_delong[[3]]
  )

  # Get data for ROC curve
  roc_curve_test_data <- data.frame(
    model = "without_prs",
    sensitivity = roc_auc_test$sensitivities,
    specificity = roc_auc_test$specificities
  )

  res <- list(
    auc_train_wo_prs = auc_training,
    auc_test_wo_prs = df_auc_test,
    roc_curve_test_wo_prs_data = roc_curve_test_data,
    full_roc_test_wo_prs = roc_auc_test,
    fit_wo_prs = fit_wo_prs
  )

  return(res)
}


# Model PRS only -----------------------------------------------
model_prs_only <- function(
  train_ctrl_stats,
  test_ctrl_stats,
  log_reg,
  norm_prs
) {
  form_prs_only <- as.formula(paste("status ~", norm_prs))

  rec_prs_only <- recipes::recipe(
    form_prs_only,
    data = train_ctrl_stats
  )

  wflow_rec_prs_only <- workflows::workflow() |>
    workflows::add_model(log_reg) |>
    workflows::add_recipe(rec_prs_only)

  # Model to obtain test set AUC
  fit_prs_only <- wflow_rec_prs_only |>
    generics::fit(data = train_ctrl_stats)

  # AUC training set
  data_folds <- rsample::vfold_cv(train_ctrl_stats, v = 10, strata = status)

  fit_resamples_prs_only <- tune::fit_resamples(
    wflow_rec_prs_only,
    resamples = data_folds,
    control = tune::control_resamples(save_pred = TRUE, save_workflow = TRUE)
  )

  auc_training <- fit_resamples_prs_only |>
    tune::collect_metrics(summarize = F) |>
    dplyr::filter(.metric == "roc_auc") |>
    dplyr::select(id, .metric, .estimate)

  # Predict probabilities
  pred_prob_prs_only <- predict(
    fit_prs_only,
    new_data = test_ctrl_stats,
    type = "prob"
  )

  # Add to testing data
  test_results_prs_only <- test_ctrl_stats |>
    dplyr::bind_cols(pred_prob_prs_only)

  # ROC AUC with confidence interval using DeLong method
  roc_auc_test <- pROC::roc(
    response = test_results_prs_only$status,
    predictor = test_results_prs_only$.pred_0
  )

  ci_delong <- pROC::ci.auc(roc_auc_test, method = "delong")

  df_auc_test <- data.frame(
    auc = roc_auc_test$auc[[1]],
    lower = ci_delong[[1]],
    upper = ci_delong[[3]]
  )

  # Plot ROC curve
  p <- pROC::ggroc(roc_auc_test, color = "steelblue", linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste0(
      "ROC Curve (AUC = ",
      round(df_auc_test$auc, 3),
      ")"
    )) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 1,
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::labs(x = "Specificity", y = "Sensitivity")

  res <- list(
    auc_train_prs_only = auc_training,
    auc_ci_test_prs_only = df_auc_test,
    roc_curve_test_prs_only = p
  )

  return(res)
}
