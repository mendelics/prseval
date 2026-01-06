#' Metrics per Standard Deviation
#' @description Get OR and AUC metrics for a certain PRS
#' @param dataset A dataframe with columns status, age at analysis or diagnosis, impersonate code and id.
#' @param prs_col Name of PRS column (character).
#' @param seed A random number to be set as a seed for the training and testing sampling to be reproducible.
#' @param recipe_var A recipe from recipes::recipe() object, to be provided optionally, in case the recipe originally included in the model (status ~ prs + scaled_centered_pc(1:10) + age_analysis) is not what you need for your PRS.
#'
#' @return A list with OR, AUC, delta AUC, and ROC curve for AUC with and without PRS.
#' @importFrom stats anova
#' @import rsample
#' @import dplyr
#' @import parsnip
#' @import workflows
#' @import pROC
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

per_sd_metrics <- function(dataset, prs_col, seed, recipe_var = NULL) {
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
    norm_prs,
    recipe_var
  )

  # Modeling without PRS ------------------------------------
  res_model_without_prs <- model_without_prs(
    train_ctrl_stats,
    test_ctrl_stats,
    log_reg,
    norm_prs,
    recipe_var
  )

  # Modeling only PRS-----------------------
  res_model_prs_only <- model_prs_only(
    train_ctrl_stats,
    test_ctrl_stats,
    log_reg,
    norm_prs,
    recipe_var
  )

  # Results---------------

  # AUC training set:
  auc_training_prs_only <- res_model_prs_only[["auc_train_prs_only"]]

  # AUC testing set:
  auc_with_prs <- round(res_model_with_prs[["auc_test_with_prs"]], 3)
  auc_wo_prs <- round(res_model_without_prs[["auc_test_wo_prs"]], 3)
  auc_prs_only <- round(res_model_prs_only[["auc_ci_test_prs_only"]], 3)

  # Likelihood Ratio Test between model with and without PRS
  fit_with_prs <- workflows::extract_fit_engine(res_model_with_prs[[
    "fit_with_prs"
  ]])
  fit_wo_prs <- workflows::extract_fit_engine(res_model_without_prs[[
    "fit_wo_prs"
  ]])

  lrt_res <- anova(fit_wo_prs, fit_with_prs, test = "LRT")
  pval_lrt <- lrt_res$`Pr(>Chi)`[2]

  # Delta AUC with DeLong test
  delong_nested_test <- pROC::roc.test(
    res_model_with_prs[["full_roc_test_with_prs"]],
    res_model_without_prs[["full_roc_test_wo_prs"]],
    method = "delong"
  )

  res_delong_nested <- data.frame(
    auc_with = auc_with_prs,
    auc_wo = auc_wo_prs,
    delta = auc_with_prs - auc_wo_prs,
    p_value = delong_nested_test$p.value
  )

  # Get roc_auc together
  all_roc_auc <- rbind(
    res_model_with_prs[["roc_curve_test_with_prs_data"]],
    res_model_without_prs[["roc_curve_test_wo_prs_data"]]
  ) |>
    dplyr::mutate(model = as.factor(model))

  p <- ggplot2::ggplot(
    all_roc_auc,
    ggplot2::aes(x = 1 - specificity, y = sensitivity, color = model)
  ) +
    ggplot2::geom_line() +
    ggplot2::theme_bw()

  return(list(
    or = res_model_with_prs[["or"]],
    lrt_res = pval_lrt,
    auc_with_prs = auc_with_prs,
    auc_wo_prs = auc_wo_prs,
    delta_auc = res_delong_nested,
    roc_comparative_curve = p,
    roc_with_prs = res_model_with_prs[["full_roc_test_with_prs"]],
    roc_wo_prs = res_model_without_prs[["full_roc_test_wo_prs"]],
    roc_curve_plot_prs_only = res_model_prs_only[["roc_curve_test_prs_only"]],
    pred_with_prs = res_model_with_prs[["pred_with_prs"]],
    auc_prs_only_training_set = auc_training_prs_only,
    auc_prs_only_testing_set = auc_prs_only
  ))
}
