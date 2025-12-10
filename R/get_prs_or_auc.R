#!/usr/bin/env Rscript

# Separate prs auc or in small functions so it is more readable

# Modeling with PRS function -------------------------------
model_with_prs <- function(train_ctrl_stats, test_ctrl_stats, log_reg) {
  rec_prs <- recipe(
    status ~ .,
    data = train_ctrl_stats
  ) |>
    step_rm(
      version,
      impersonate_code
    ) |>
    step_normalize(starts_with("pc"))

  wflow_prs <- workflow() |>
    add_model(log_reg) |>
    add_recipe(rec_prs)

  fit_prs <- wflow_prs |>
    fit(data = train_ctrl_stats)

  or <- tidy(fit_prs, exponentiate = T, conf.int = T) |>
    # filter(!grepl("pc", term), term != "(Intercept)") |>
    select(!all_of(c("std.error", "statistic")))

  # Predict probabilities
  pred_prob_prs_vars <- predict(
    fit_prs,
    new_data = test_ctrl_stats,
    type = "prob"
  )

  # Add to testing data
  results_prs_vars <- test_ctrl_stats |>
    bind_cols(pred_prob_prs_vars)

  # Obtain ROC AUC **with PRS**
  auc_vars <- roc_auc(results_prs_vars, truth = status, .pred_0) |>
    select(!any_of(".estimator"))

  # Get ROC AUC values
  roc_auc_vars <- results_prs_vars |>
    roc_curve(truth = status, .pred_0) |>
    mutate(model = "with PRS")

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
  rec_no_prs <- recipe(
    status ~ .,
    data = train_ctrl_stats
  ) |>
    step_rm(
      version,
      all_of(norm_prs),
      impersonate_code
    ) |>
    step_normalize(starts_with("pc"))

  wflow_rec_no_prs <- workflow() |> add_model(log_reg) |> add_recipe(rec_no_prs)

  fit_no_prs <- wflow_rec_no_prs |>
    fit(data = train_ctrl_stats)

  # Predict probabilities
  pred_prob_no_prs <- predict(
    fit_no_prs,
    new_data = test_ctrl_stats,
    type = "prob"
  )

  # Add to testing data
  results_no_prs <- test_ctrl_stats |>
    bind_cols(pred_prob_no_prs)

  # AUC value
  auc_no_prs <- roc_auc(results_no_prs, truth = status, .pred_0) |>
    select(!any_of(".estimator"))

  # Get ROC AUC values
  roc_auc_no_prs <- results_no_prs |>
    roc_curve(truth = status, .pred_0) |>
    mutate(model = "without PRS")

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

  rec_prs_only <- recipe(
    form_prs_only,
    data = train_ctrl_stats
  )

  wflow_rec_prs_only <- workflow() |>
    add_model(log_reg) |>
    add_recipe(rec_prs_only)

  fit_prs_only <- wflow_rec_prs_only |>
    fit(data = train_ctrl_stats)

  # Predict probabilities
  pred_prob_prs_only <- predict(
    fit_prs_only,
    new_data = test_ctrl_stats,
    type = "prob"
  )

  # Add to testing data
  results_prs_only <- test_ctrl_stats |>
    bind_cols(pred_prob_prs_only)

  # Obtain ROC AUC
  auc_prs_only <- roc_auc(results_prs_only, truth = status, .pred_0) |>
    select(!any_of(".estimator"))

  # Plot ROC AUC
  roc_auc_prs_only <- results_prs_only |>
    roc_curve(truth = status, .pred_0) |>
    autoplot()

  res <- list(auc_prs_only = auc_prs_only, roc_auc_prs_only = roc_auc_prs_only)

  return(res)
}

# Get PRS OR, AUC and ROC curves --------------------------------
get_prs_or_auc <- function(dataset, prs_col, seed) {
  set.seed(seed)

  # Set PRS name with "norm" because we use the normalized version for the analyses
  norm_prs <- paste0("norm_", prs_col)

  split <- initial_split(dataset, strata = status, prop = 0.75)

  train <- training(split)
  test <- testing(split)

  control_stats_train <- get_control_stats(train) |> as.data.frame()

  control_train_mean <- control_stats_train[, paste0(prs_col, "_mean")]
  control_train_sd <- control_stats_train[, paste0(prs_col, "_sd")]

  train_ctrl_stats <- train |>
    mutate(
      "{norm_prs}" := (get(prs_col) - control_train_mean) /
        control_train_sd
    ) |>
    select(
      !starts_with("prs")
    )

  control_stats_test <- get_control_stats(test) |> as.data.frame()

  control_test_mean <- control_stats_test[, paste0(prs_col, "_mean")]
  control_test_sd <- control_stats_test[, paste0(prs_col, "_sd")]

  test_ctrl_stats <- test |>
    mutate(
      "{norm_prs}" := (get(prs_col) - control_test_mean) /
        control_test_sd
    ) |>
    select(
      !starts_with("prs")
    )

  log_reg <- logistic_reg() |>
    set_engine("glm") |>
    set_mode("classification")

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

  p <- ggplot(
    all_roc_auc,
    aes(x = 1 - specificity, y = sensitivity, color = model)
  ) +
    geom_line() +
    theme_bw()

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
