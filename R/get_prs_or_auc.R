#' Get OR and AUC metrics for a certain PRS
#'
#' @param dataset A dataframe with columns status, age, tier1(0/1 for presence of pathogenic or likely-pathogenic variants in tier1 list for the disease evaluated).
#' @param prs_col Name of PRS column (character).
#' @param seed A random number to be set as a seed for the training and testing sampling to be reproducible.
#'
#' @return A list with OR, AUC, delta AUC, and ROC curve for AUC with and without PRS.
#' @importFrom generics fit
#' @importFrom stats predict
#' @import rsample
#' @import dplyr
#' @import tidyr
#' @importFrom broom tidy
#' @import parsnip
#' @import recipes
#' @import workflows
#' @import yardstick
#' @import ggplot2
#' @export
#'
#' @examples
#' # get_prs_or_auc(breast_v2_v3_mut, x, seed = 28)
#' # TBD

# TODO: add ensurance that status is a factor variable in all model functions or in the get_pgs... functions

# Get control statistics ------------------------------------
get_control_stats <- function(df) {
  df |>
    dplyr::filter(status == 0) |>
    dplyr::summarise(dplyr::across(
      dplyr::starts_with("prs"),
      list(mean = mean, sd = sd),
      .names = "{.col}_{.fn}"
    ))
}

# Modeling with PRS function -------------------------------
model_with_prs <- function(train_ctrl_stats, test_ctrl_stats, log_reg) {
  rec_prs <- recipes::recipe(
    status ~ .,
    data = train_ctrl_stats
  ) |>
    recipes::step_rm(
      version,
      impersonate_code
    ) |>
    recipes::step_normalize(starts_with("pc"))

  wflow_prs <- workflows::workflow() |>
    workflows::add_model(log_reg) |>
    workflows::add_recipe(rec_prs)

  fit_prs <- wflow_prs |>
    generics::fit(data = train_ctrl_stats)

  or <- broom::tidy(fit_prs, exponentiate = T, conf.int = T) |>
    # filter(!grepl("pc", term), term != "(Intercept)") |>
    dplyr::select(!dplyr::all_of(c("std.error", "statistic")))

  # Predict probabilities
  pred_prob_prs_vars <- stats::predict(
    fit_prs,
    new_data = test_ctrl_stats,
    type = "prob"
  )

  # Add to testing data
  results_prs_vars <- test_ctrl_stats |>
    dplyr::bind_cols(pred_prob_prs_vars)

  # Obtain ROC AUC **with PRS**
  auc_vars <- yardstick::roc_auc(results_prs_vars, truth = status, .pred_0) |>
    dplyr::select(!dplyr::any_of(".estimator"))

  # Get ROC AUC values
  roc_auc_vars <- results_prs_vars |>
    yardstick::roc_curve(truth = status, .pred_0) |>
    dplyr::mutate(model = "with PRS")

  res <- list(or = or, auc_vars = auc_vars, roc_auc_vars = roc_auc_vars)

  return(res)
}


# Modeling without PRS function ---------------------------------
model_without_prs <- function(
  train_ctrl_stats,
  test_ctrl_stats,
  log_reg,
  norm_prs
) {
  rec_no_prs <- recipes::recipe(
    status ~ .,
    data = train_ctrl_stats
  ) |>
    recipes::step_rm(
      version,
      all_of(norm_prs),
      impersonate_code
    ) |>
    recipes::step_normalize(starts_with("pc"))

  wflow_rec_no_prs <- workflows::workflow() |>
    workflows::add_model(log_reg) |>
    workflows::add_recipe(rec_no_prs)

  fit_no_prs <- wflow_rec_no_prs |>
    generics::fit(data = train_ctrl_stats)

  # Predict probabilities
  pred_prob_no_prs <- stats::predict(
    fit_no_prs,
    new_data = test_ctrl_stats,
    type = "prob"
  )

  # Add to testing data
  results_no_prs <- test_ctrl_stats |>
    dplyr::bind_cols(pred_prob_no_prs)

  # AUC value
  auc_no_prs <- yardstick::roc_auc(results_no_prs, truth = status, .pred_0) |>
    dplyr::select(!dplyr::any_of(".estimator"))

  # Get ROC AUC values
  roc_auc_no_prs <- results_no_prs |>
    yardstick::roc_curve(truth = status, .pred_0) |>
    dplyr::mutate(model = "without PRS")

  res <- list(auc_no_prs = auc_no_prs, roc_auc_no_prs = roc_auc_no_prs)

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

  fit_prs_only <- wflow_rec_prs_only |>
    generics::fit(data = train_ctrl_stats)

  # Predict probabilities
  pred_prob_prs_only <- stats::predict(
    fit_prs_only,
    new_data = test_ctrl_stats,
    type = "prob"
  )

  # Add to testing data
  results_prs_only <- test_ctrl_stats |>
    dplyr::bind_cols(pred_prob_prs_only)

  # Obtain ROC AUC
  auc_prs_only <- yardstick::roc_auc(
    results_prs_only,
    truth = status,
    .pred_0
  ) |>
    dplyr::select(!dplyr::any_of(".estimator"))

  # Plot ROC AUC
  roc_auc_prs_only <- results_prs_only |>
    yardstick::roc_curve(truth = status, .pred_0) |>
    ggplot2::autoplot()

  res <- list(auc_prs_only = auc_prs_only, roc_auc_prs_only = roc_auc_prs_only)

  return(res)
}


# Get PRS OR, AUC and ROC curves --------------------------------
get_prs_or_auc <- function(dataset, prs_col, seed) {
  stopifnot(is.double(dataset[, prs_col]), is.factor(dataset[, "status"]))

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
      !starts_with("prs")
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

  # Modeling with PRS ------------------------------------------------------------------
  res_model_with_prs <- model_with_prs(
    train_ctrl_stats,
    test_ctrl_stats,
    log_reg
  )

  # Modeling without PRS ----------------------------------------------------------------------------------------------
  res_model_without_prs <- model_without_prs(
    train_ctrl_stats,
    test_ctrl_stats,
    log_reg,
    norm_prs
  )

  # Modeling only PRS------------------------------------------------------------------------------------------------
  res_model_prs_only <- model_prs_only(
    train_ctrl_stats,
    test_ctrl_stats,
    log_reg,
    norm_prs
  )

  # Results-------------------------------------------------------------------------------

  # AUC:
  auc_with_prs <- round(res_model_with_prs[["auc_vars"]]$.estimate, 3)
  auc_wo_prs <- round(res_model_without_prs[["auc_no_prs"]]$.estimate, 3)
  auc_prs_only <- round(res_model_prs_only[["auc_prs_only"]]$.estimate, 3)

  # Delta AUC
  delta_auc <- auc_with_prs - auc_wo_prs

  # get roc_auc together
  all_roc_auc <- rbind(
    res_model_with_prs[["roc_auc_vars"]],
    res_model_without_prs[["roc_auc_no_prs"]]
  )

  p <- ggplot2::ggplot(
    all_roc_auc,
    ggplot2::aes(x = 1 - specificity, y = sensitivity, color = model)
  ) +
    ggplot2::geom_line() +
    ggplot2::theme_bw()

  return(list(
    or = res_model_with_prs[["or"]],
    auc_prs_only = auc_prs_only,
    auc_with_prs = auc_with_prs,
    auc_wo_prs = auc_wo_prs,
    delta_auc = delta_auc,
    roc_comparative_curve = p,
    roc_prs_only = res_model_prs_only[["roc_auc_prs_only"]]
  ))
}
