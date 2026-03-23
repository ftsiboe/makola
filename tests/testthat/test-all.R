test_that("all functions work correctly", {

  data <- simulate_market_prices(
    start_date           = "2000-01-01",
    end_date             = Sys.Date(),
    frequency            = "month",
    n_markets            = 5,
    n_integrated_markets = 3,
    seed                 = 123
  )

  expect_equal(max(nrow(data) > 1), 1)

  res <- estimate_pair_integration_model(
    pair    = c("market_1","market_2"),
    data    = data,
    nboot   = 100,
    lag.max = 2,
    lag.opt = 2,
    trim    = 0.10)

  expect_equal(max(length(res) > 1), 1)

  for(target_model in c("TVECM3", "TVECM2", "AVECM")) {
    res$model <- target_model
    res$est   <- estimate_threshold_vecm(res)
    expect_equal(max(length(res) > 1), 1)
  }

  for(target_model in c("VECM")) {
    res$model <- target_model
    res$est <- estimate_linear_vecm(res)
    expect_equal(max(length(res) > 1), 1)
  }

  for(target_model in c("DTVARDL2", "DTVARDL3", "TVARDL3", "TVARDL2")) {
    res$model <- target_model
    res$est <- estimate_threshold_vardl(res)
    expect_equal(max(length(res) > 1), 1)
  }

  for(target_model in c("VARDL", "DVARDL")) {
    res$model <- target_model
    res$est <- estimate_linear_vardl(res)
    expect_equal(max(length(res) > 1), 1)
  }

})
