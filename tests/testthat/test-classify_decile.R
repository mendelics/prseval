test_that("classify_decile intervals have the en-dash ASCII-compatible", {
  set.seed(28)
  prs <- rnorm(100)
  df <- data.frame(prs)

  deciles <- classify_decile(df, "prs")

  expect_true(all(grepl("\u2013", deciles)))
})
