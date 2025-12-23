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
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' #TBD
