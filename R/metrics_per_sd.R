#' Metrics per Standard Deviation
#' @description Get OR and AUC metrics for a certain PRS
#' @param dataset A dataframe with columns status, age at analysis or diagnosis, impersonate code and id.
#' @param prs_col Name of PRS column (character).
#' @param seed A random number to be set as a seed for the training and testing sampling to be reproducible.
#'
#' @return A list with OR, AUC, delta AUC, and ROC curve for AUC with and without PRS.
#' @importFrom generics fit
#' @importFrom stats predict
#' @importFrom stats as.formula
#' @importFrom stats anova
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
#' # PCs 1 to 10
#' n_cols <- 10
#' n_rows <- 100
#' data_matrix <- replicate(
#'   n_cols,
#'   rnorm(n_rows, mean = runif(1, 0, 10), sd = runif(1, 1, 5))
#' )
#' df_pcs <- as.data.frame(data_matrix)
#' colnames(df_pcs) <- paste0("pc", 1:n_cols)
#' # Data mock
#' data_mock <- data.frame(
#'   impersonate_code = paste0(
#'     sample(LETTERS, 100, replace = TRUE),
#'     sample(LETTERS, 100, replace = TRUE),
#'     sample(LETTERS, 100, replace = TRUE),
#'     sample(1:9, 100, replace = TRUE),
#'     sample(1:9, 100, replace = TRUE),
#'     sample(1:9, 100, replace = TRUE)
#'   ),
#'   status = as.factor(sample(c(0, 1), size = 100, replace = TRUE)),
#'   age_analysis = round(runif(100, 19, 80)),
#'   tier1 = as.factor(sample(
#'     c(0, 1),
#'     size = 100,
#'     prob = c(0.8, 0.2),
#'     replace = TRUE
#'   )),
#'   prs_test = rnorm(n = 100),
#'   version = sample(c("v1", "v2"), size = 100, replace = TRUE)
#' ) |>
#'   cbind(df_pcs)
#'
#' # Arguments
#' seed <- 28
#' prs_col_mock <- "prs_test"
#'
#' per_sd_metrics(data_mock, prs_col_mock, seed)
#'

per_sd_metrics <- function(dataset, prs_col, seed) {
  stopifnot(
    is.double(dataset |> dplyr::pull({{ prs_col }})),
    is.factor(dataset |> dplyr::pull(status))
  )

  set.seed(seed)

  # Set PRS name with "norm" because we use the normalized version for the analyses
  norm_prs <- paste0("norm_", prs_col)

  split <- rsample::initial_split(dataset, strata = status, prop = 0.75)

  train <- rsample::training(split)
  test <- rsample::testing(split)

  control_stats_train <- get_control_stats(train) |> as.data.frame()

  control_train_mean <- control_stats_train[, paste0(prs_col, "_mean")]
  control_train_sd <- control_stats_train[, paste0(prs_col, "_sd")]

  train_ctrl_stats <- train |>
    dplyr::mutate(
      "{norm_prs}" := (get(prs_col) - control_train_mean) /
        control_train_sd
    ) |>
    dplyr::select(
      !dplyr::starts_with("prs")
    )

  control_stats_test <- get_control_stats(test) |> as.data.frame()

  control_test_mean <- control_stats_test[, paste0(prs_col, "_mean")]
  control_test_sd <- control_stats_test[, paste0(prs_col, "_sd")]

  test_ctrl_stats <- test |>
    dplyr::mutate(
      "{norm_prs}" := (get(prs_col) - control_test_mean) /
        control_test_sd
    ) |>
    dplyr::select(
      !starts_with("prs")
    )

  log_reg <- parsnip::logistic_reg() |>
    parsnip::set_engine("glm") |>
    parsnip::set_mode("classification")

  # Modeling with PRS -------------------
  res_model_with_prs <- model_with_prs(
    train_ctrl_stats,
    test_ctrl_stats,
    log_reg,
    norm_prs
  )

  # Modeling without PRS ------------------------------------
  res_model_without_prs <- model_without_prs(
    train_ctrl_stats,
    test_ctrl_stats,
    log_reg,
    norm_prs
  )

  # Modeling only PRS-----------------------
  res_model_prs_only <- model_prs_only(
    train_ctrl_stats,
    test_ctrl_stats,
    log_reg,
    norm_prs
  )

  # Results---------------

  # AUC training set:
  auc_training_prs_only <- res_model_prs_only[["auc_train_prs_only"]]

  # AUC testing set:
  auc_with_prs <- round(res_model_with_prs[["auc_test_with_prs"]], 3)
  auc_wo_prs <- round(res_model_without_prs[["auc_test_wo_prs"]], 3)
  auc_prs_only <- round(res_model_prs_only[["auc_ci_test_prs_only"]], 3)

  # Likelihood Ratio Test between model with and without PRS
  fit_with_prs <- extract_fit_engine(res_model_with_prs[["fit_with_prs"]])
  fit_wo_prs <- extract_fit_engine(res_model_without_prs[["fit_wo_prs"]])

  lrt_res <- anova(fit_wo_prs, fit_with_prs, test = "Chisq")

  # Delta AUC
  delta_auc <- auc_with_prs$auc_test - auc_wo_prs$auc_test

  # Get roc_auc together #TODO: get data back here
  all_roc_auc <- rbind(
    res_model_with_prs[["roc_curve_test_with_prs_data"]],
    res_model_without_prs[["roc_curve_test_without_prs_data"]]
  )

  p <- ggplot2::ggplot(
    all_roc_auc,
    ggplot2::aes(x = 1 - specificity, y = sensitivity, color = model)
  ) +
    ggplot2::geom_line() +
    ggplot2::theme_bw()

  return(list(
    or = res_model_with_prs[["or"]],
    lrt_res = lrt_res,
    auc_with_prs = auc_with_prs,
    auc_wo_prs = auc_wo_prs,
    delta_auc = delta_auc,
    roc_comparative_curve = p,
    roc_with_prs = res_model_with_prs[["full_roc_test_with_prs"]],
    roc_wo_prs = res_model_with_prs[["full_roc_test_wo_prs"]],
    roc_prs_only = res_model_prs_only[["roc_auc_prs_only"]],
    auc_prs_only_training_set = auc_training_prs_only,
    auc_prs_only_testing_set = auc_prs_only
  ))
}

# Modeling with PRS function -------------------------------
model_with_prs <- function(
  train_ctrl_stats,
  test_ctrl_stats,
  log_reg,
  norm_prs
) {
  # Formula with predictors wanted
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

  wflow_prs <- workflows::workflow() |>
    workflows::add_model(log_reg) |>
    workflows::add_recipe(rec_prs)

  # Odds Ratio
  fit_prs <- wflow_prs |>
    generics::fit(data = train_ctrl_stats)

  or <- broom::tidy(fit_prs, exponentiate = T, conf.int = T) |>
    dplyr::select(!dplyr::all_of(c("std.error", "statistic")))

  # AUC training set
  data_folds <- vfold_cv(train_ctrl_stats, v = 10, strata = status)

  fit_resamples_full_model <- tune::fit_resamples(
    wflow_prs,
    resamples = data_folds,
    control = tune::control_resamples(save_pred = TRUE, save_workflow = TRUE)
  )

  auc_training <- fit_resamples_full_model |>
    tune::collect_metrics(summarize = F) |>
    filter(.metric == "roc_auc") |>
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
    auc_test = roc_auc_test$auc[[1]],
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
    fit_with_prs = fit_prs
  )

  return(res)
}


# Modeling without PRS function ---------------------------------
model_without_prs <- function(
  train_ctrl_stats,
  test_ctrl_stats,
  log_reg,
  norm_prs
) {
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

  wflow_wo_prs <- workflows::workflow() |>
    workflows::add_model(log_reg) |>
    workflows::add_recipe(rec_wo_prs)

  # Fit for test AUC
  fit_wo_prs <- wflow_wo_prs |>
    generics::fit(data = train_ctrl_stats)

  # AUC training set
  data_folds <- vfold_cv(train_ctrl_stats, v = 10, strata = status)

  fit_resamples_model_wo_prs <- tune::fit_resamples(
    wflow_wo_prs,
    resamples = data_folds,
    control = tune::control_resamples(save_pred = TRUE, save_workflow = TRUE)
  )

  auc_training <- fit_resamples_model_wo_prs |>
    tune::collect_metrics(summarize = F) |>
    filter(.metric == "roc_auc") |>
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
    auc_test = roc_auc_test$auc[[1]],
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
    auc_test = roc_auc_test$auc[[1]],
    lower = ci_delong[[1]],
    upper = ci_delong[[3]]
  )

  # Plot ROC curve
  p <- pROC::ggroc(roc_auc_test, color = "steelblue", size = 1) +
    theme_minimal() +
    ggtitle(paste0("ROC Curve (AUC = ", round(auc(roc_auc_test), 3), ")")) +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey") +
    labs(x = "Specificity", y = "Sensitivity")

  res <- list(
    auc_train_prs_only = auc_training,
    auc_ci_test_prs_only = df_auc_test,
    roc_curve_test_prs_only = p
  )

  return(res)
}
