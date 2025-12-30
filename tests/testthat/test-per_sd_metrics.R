setup_mock_df <- function(n = 100) {
  # Mock PCs 1 to 10
  n_cols <- 10
  n_rows <- n
  data_matrix <- replicate(
    n_cols,
    rnorm(n_rows, mean = runif(1, 0, 10), sd = runif(1, 1, 5))
  )
  df_pcs <- as.data.frame(data_matrix)
  colnames(df_pcs) <- paste0("pc", 1:n_cols)

  # Mock df
  data_mock <- data.frame(
    impersonate_code = paste0(
      sample(LETTERS, n, replace = T),
      sample(LETTERS, n, replace = T),
      sample(LETTERS, n, replace = T),
      sample(1:9, n, replace = T),
      sample(1:9, n, replace = T),
      sample(1:9, n, replace = T)
    ),
    status = as.factor(sample(c(0, 1), size = n, replace = T)),
    age_analysis = round(runif(n, 19, 80)),
    tier1 = as.factor(sample(
      c(0, 1),
      size = n,
      prob = c(0.8, 0.2),
      replace = T
    )),
    prs_test = rnorm(n = n),
    version = sample(c("v1", "v2"), size = n, replace = T)
  ) |>
    cbind(df_pcs)
}

test_that("per_sd_metrics returns a list", {
  data_mock <- setup_mock_df()

  # Run
  res <- per_sd_metrics(dataset = data_mock, prs_col = "prs_test", seed = 28)

  # Test
  expect_type(res, "list")
})

test_that("per_sd_metrics returns the correct list structure", {
  data_mock <- setup_mock_df()

  # Run
  res <- per_sd_metrics(dataset = data_mock, prs_col = "prs_test", seed = 28)

  # Test
  expect_named(
    res,
    c(
      "or",
      "lrt_res",
      "auc_with_prs",
      "auc_wo_prs",
      "delta_auc",
      "roc_comparative_curve",
      "roc_with_prs",
      "roc_wo_prs",
      "roc_prs_only",
      "auc_prs_only_training_set",
      "auc_prs_only_testing_set"
    )
  )
})

test_that("per_sd_metrics returns the correct auc_prs_only structure", {
  data_mock <- setup_mock_df()

  # Run
  res <- per_sd_metrics(dataset = data_mock, prs_col = "prs_test", seed = 82)

  # Test
  expect_s3_class(res[["auc_prs_only_training_set"]], "data.frame")
})
