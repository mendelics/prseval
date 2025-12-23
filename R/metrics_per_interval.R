#' Metrics per Interval of PRS Values
#' @description Get OR metrics for a certain PRS per defined intervals according to PRS values.
#' @param dataset A dataframe with columns status, age_analysis (age at analysis for controls or age at diagnosis for cases), and impersonate_code (individuals id).
#' @param prs_col Name of PRS column (character).
#' @param intervals_mode Either "decile" or "mavaddat" to define the intervals that should be created in the dataset. See classify_decile() and classify_mavaddat() to understand the difference.
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
#' @import parsnip
#' @import workflows
#' @import yardstick
#'
#' @export
#'
#' @examples
#' #TBD

per_interval_metrics <- function(
  dataset,
  prs_col,
  intervals_mode
) {
  stopifnot(
    is.double(dataset |> pull({{ prs_col }})),
    is.factor(dataset |> pull(status)),
    is.character(intervals_mode)
  )

  # Set function to classify intervals based on intervals_mode
  # 1. classificar PRS em decis (interval)
  if (intervals_mode == "decile") {
    calc_intervals <- classify_decile(dataset, prs_col)
  } else if (intervals_mode == "mavaddat") {
    calc_intervals <- classify_mavaddat(dataset, prs_col)
  } else {
    stop("Unknown definition of intervals_mode.")
  }

  df_intervals <- dataset |>
    dplyr::mutate(interval = calc_intervals)

  # Model base to be executed below
  log_reg <- parsnip::logistic_reg() |>
    parsnip::set_engine("glm") |>
    parsnip::set_mode("classification")

  # All intervals to generate metrics
  all_intervals <- setdiff(
    levels(df_intervals$interval),
    c("40\u201350%", "50%\u201360%")
  )

  list_of_results <- lapply(all_intervals, function(intv) {
    model_by_interval(
      df_intervals = df_intervals,
      prs_col = prs_col,
      model = log_reg,
      interval = intv
    )
  })

  res <- do.call(rbind, list_of_results) |>
    mutate(
      interval = factor(interval, levels = levels(calc_intervals), ordered = T)
    )

  return(res)
}


# Apply logistic regression model to a specific interval ----
model_by_interval <- function(df_intervals, prs_col, model, interval) {
  df_spec_interval <- df_intervals |>
    # drop everyone that is not the selected interval or the middle deciles
    filter(interval %in% c("40%–50%", "50%–60%", {{ interval }})) |>
    mutate(
      is_interval = if_else(interval == {{ interval }}, 1, 0) # important! all categorical predictors need to be numeric for logistic regression with glm.
    ) |>
    dplyr::select(
      impersonate_code,
      status,
      is_interval,
      age_analysis,
      starts_with("pc")
    )

  # Recipe for logistic regression
  rec <- recipes::recipe(
    status ~ .,
    data = df_spec_interval
  ) |>
    recipes::step_rm(
      impersonate_code
    )

  # Create workflow
  wflow <- workflows::workflow() |>
    add_model(model) |>
    add_recipe(rec)

  # Fit
  fit_prs <- wflow |>
    generics::fit(data = df_spec_interval)

  # Extract metrics
  logreg_stats <- broom::tidy(fit_prs, exponentiate = T, conf.int = T) |>
    dplyr::filter(term == "is_interval") |>
    dplyr::mutate(
      prs = {{ prs_col }},
      interval = {{ interval }},
      .before = 1
    ) |>
    dplyr::select(!dplyr::all_of(c("std.error", "statistic", "term"))) |>
    dplyr::relocate(p.value, .after = dplyr::everything())

  return(logreg_stats)
}
