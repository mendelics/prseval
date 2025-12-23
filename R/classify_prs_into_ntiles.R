#' Classify Subjects Into Intervals According To PRS Values
#' @description Classification of individuals into intervals according to deciles and to breaks present in Mavaddat et al., 2019.
#'
#' @param prs_values A numeric vector with PRS values.
#' @return A character vector with intervals in which each PRS value was classified.
#'
#' @importFrom stats quantile
#' @importFrom stats sd
#' @export
#'
#' @examples
#' #TBD

# Decile cutoffs ---------------------------------------------------------
decile_breaks <- function(prs_values) {
  breaks <- quantile(
    prs_values,
    probs = base::seq(0, 1.0, by = 0.1),
    na.rm = T
  )

  return(breaks)
}


# Mavaddat paper cutoffs ------------------------------------------
mavaddat_breaks <- function(prs_values) {
  breaks <- quantile(
    prs_values,
    probs = base::c(
      0,
      0.01,
      0.05,
      0.10,
      0.20,
      0.40,
      0.50,
      0.60,
      0.80,
      0.90,
      0.95,
      0.99,
      1.0
    ),
    na.rm = T
  )

  return(breaks)
}
