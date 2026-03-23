# Hard reset of workspace
rm(list = ls(all = TRUE));gc()
library(data.table)

devtools::document()

data <- simulate_market_prices(
  start_date           = "2000-01-01",
  end_date             = Sys.Date(),
  frequency            = "month",
  n_markets            = 5,
  n_integrated_markets = 3,
  seed                 = 123
)

res <- estimate_pair_integration_model(
    pair    = c("market_1","market_2"),
    data    = data,
    nboot   = 100,
    lag.max = 2,
    lag.opt = 2,
    trim    = 0.10)





length(res) > 1






res <- select_pair_integration_model(
  pair = c("market_1","market_2"),
  data    = data,
  nboot   = 100,
  lag.max = 2,
  lag.opt = 2,
  trim    = 0.10)



fallback_models <- switch(
  EXPR = res$model,
  "TVECM3"   = c("TVECM3", "TVECM2", "AVECM", "VECM"),
  "TVECM2"   = c("TVECM2", "AVECM", "VECM"),
  "AVECM"    = c("AVECM", "VECM"),
  "VECM"     = c("VECM"),
  "TVARDL3"  = c("TVARDL3", "TVARDL2", "VARDL"),
  "TVARDL2"  = c("TVARDL2", "VARDL"),
  "VARDL"    = c("VARDL"),
  "DTVARDL3" = c("DTVARDL3", "DTVARDL2", "DVARDL"),
  "DTVARDL2" = c("DTVARDL2", "DVARDL"),
  "DVARDL"   = c("DVARDL"),
  res$model
)







