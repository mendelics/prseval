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

# Categorize samples into deciles -------------------------------------------------------
classify_decile <- function(df, prs_col) {
  stopifnot(is.character(prs_col), length(prs_col) <= 1)
  prs <- df[, prs_col] |> unlist()

  breaks <- decile_breaks(prs)

  # Adjust breaks for min and max values to be 2 units less and 2 units more than
  # This is because if we need to classify new samples into the breaks obtained here,
  breaks_adj <- breaks
  breaks_adj[1] <- breaks[1] - 2
  breaks_adj[length(breaks_adj)] <- breaks[length(breaks)] + 2

  # Create names for the intervals
  names_left <- names(breaks_adj)[-length(breaks_adj)]
  names_right <- names(breaks_adj)[-1]
  intervals <- paste0(names_left, sep = "\u2013", names_right)

  decile_intervals <- cut(
    prs,
    breaks = breaks_adj,
    labels = intervals,
    right = T, # value on the left of the interval is not included (only bigger than that value), and value on the right of the interval is included
    ordered_result = T
  )
  return(decile_intervals)
}


# Categorize samples into Mavaddat percentiles ----------------------------------------------------
classify_mavaddat <- function(df, prs_col) {
  stopifnot(is.character(prs_col), length(prs_col) <= 1)
  prs <- df[, prs_col] |> unlist()

  breaks <- mavaddat_breaks(prs)

  # Adjust breaks for min and max values to be 2 units less and 2 units more than
  # This is because if we need to classify new samples into the breaks obtained here,
  breaks_adj <- breaks
  breaks_adj[1] <- breaks[1] - 2
  breaks_adj[length(breaks_adj)] <- breaks[length(breaks)] + 2

  # Create names for the intervals
  names_left <- names(breaks_adj)[-length(breaks_adj)]
  names_right <- names(breaks_adj)[-1]
  intervals <- paste0(names_left, sep = "\u2013", names_right)

  mavaddat_intervals <- cut(
    prs,
    breaks = breaks_adj,
    labels = intervals,
    right = T,
    ordered_result = T
  )
  return(mavaddat_intervals)
}
