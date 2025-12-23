#' Obtain mean and standard deviation for only controls of all PRSss
#'
#' @param df A dataframe with columns status, and "prs" in the start of all PRS columns.
#' @return A tibble with columns named "prs\\[something\\]_mean", "prs\\[something\\]_sd" for all PRS columns, where mean and sd were calculated based on controls only (status == 0).
#'
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @importFrom dplyr across
#' @importFrom dplyr starts_with
#' @importFrom stats sd
#'
#' @export
#'
#' @examples
#' #TBD

# Get control statistics ------------------------------------
get_control_stats <- function(df) {
  df |>
    dplyr::filter(status == 0) |>
    dplyr::summarise(dplyr::across(
      dplyr::starts_with("prs"),
      base::list(mean = base::mean, sd = stats::sd),
      .names = "{.col}_{.fn}"
    ))
}
