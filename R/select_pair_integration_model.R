#' Select a bivariate linear or threshold time-series model
#'
#' Applies a sequence of lag selection, Granger-ordering, unit root,
#' cointegration, and threshold diagnostic tests to a pair of time series, then
#' assigns a model class based on the combined test outcomes.
#'
#' @param pair Character or integer vector of length 2 identifying the two
#'   columns in `data` to use.
#' @param data A data frame or object coercible to a data frame containing the
#'   time-series variables.
#' @param nboot Positive integer giving the number of bootstrap replications
#'   used in the Hansen-Seo and Hansen threshold tests. Default is `100`.
#' @param lag.max Positive integer giving the maximum lag order considered in
#'   `select_lag_order()`. Default is `2`.
#' @param lag.opt Positive integer giving a user-imposed lower bound on the lag
#'   order used in the final testing sequence. Default is `2`.
#' @param trim Numeric scalar in `(0, 0.5)` giving the trimming proportion used
#'   in `hansen_seo_threshold_test()`. Default is `0.10`.
#'
#' @return A list containing:
#' \describe{
#'   \item{names}{Character vector giving the reordered series names returned by
#'   `order_pair_by_granger()`.}
#'   \item{P1}{First reordered series.}
#'   \item{P2}{Second reordered series.}
#'   \item{granger}{Data frame of Granger-causality test results.}
#'   \item{URoot}{Data frame of unit root screening results from
#'   `identify_i1_specification()`.}
#'   \item{TRACE}{Data frame of Johansen trace test results from
#'   `trace_pair_test()`.}
#'   \item{HStesT}{Data frame of Hansen-Seo threshold cointegration test
#'   results from `hansen_seo_threshold_test()`.}
#'   \item{HtesT}{Data frame of Hansen linearity test results from
#'   `hansen_linearity_test()`.}
#'   \item{model}{Selected model class.}
#'   \item{lag}{Selected lag order used in the final testing sequence.}
#'   \item{cor}{Sample correlation matrix of the selected pair.}
#'   \item{cov}{Sample covariance matrix of the selected pair.}
#' }
#'
#' @details
#' The function proceeds as follows:
#' \enumerate{
#'   \item selects a lag order using `select_lag_order()`,
#'   \item reorders the pair using `order_pair_by_granger()`,
#'   \item screens the pair for I(1)-consistent behavior using
#'   `identify_i1_specification()`,
#'   \item updates the lag order using the unit root output and user-supplied
#'   lower bound,
#'   \item applies a Johansen trace test using `trace_pair_test()`,
#'   \item applies the Hansen-Seo threshold cointegration test using
#'   `hansen_seo_threshold_test()`,
#'   \item applies Hansen's linearity test using `hansen_linearity_test()`,
#'   \item selects a model class based on the resulting decision rules.
#' }
#'
#' The candidate model classes are:
#' \describe{
#'   \item{`"TVECM3"`}{three-regime threshold VECM}
#'   \item{`"TVECM2"`}{two-regime threshold VECM}
#'   \item{`"AVECM"`}{asymmetric VECM}
#'   \item{`"VECM"`}{linear VECM}
#'   \item{`"DTVARDL3"`}{threshold differenced VARDL with three regimes}
#'   \item{`"DTVARDL2"`}{threshold differenced VARDL with two regimes}
#'   \item{`"DVARDL"`}{linear differenced VARDL}
#'   \item{`"TVARDL3"`}{threshold VARDL with three regimes}
#'   \item{`"TVARDL2"`}{threshold VARDL with two regimes}
#'   \item{`"VARDL"`}{linear VARDL}
#' }
#' @export
select_pair_integration_model  <- function(
    pair,
    data,
    nboot = 100,
    lag.max = 2,
    lag.opt = 2,
    trim = 0.10
) {
  selected_lag <- select_lag_order(
    pair = pair,
    data = data,
    lag.max = lag.max
  )

  Granger <- order_pair_by_granger(
    pair = pair,
    data = data,
    lag = selected_lag
  )

  URoot <- identify_i1_specification(
    pair = pair,
    data = data
  )

  selected_lag <- max(
    c(
      lag.opt,
      min(c(URoot[grepl("lags", URoot$name), "Estimate"], selected_lag))
    )
  )

  TRACE <- trace_pair_test(
    pair = pair,
    data = data,
    lag = selected_lag,
    model = unique(URoot$model)
  )

  HStesT <- hansen_seo_threshold_test(
    pair = pair,
    data = data,
    lag = selected_lag,
    nboot = nboot,
    intercept = unique(URoot$model) %in% c("trend", "const", "both"),
    trim = trim
  )

  HtesT <- hansen_linearity_test(
    pair = pair,
    data = data,
    lag = selected_lag,
    nboot = nboot
  )

  O1    <- max(URoot$O1) > 0
  Coint <- TRACE[TRACE$name %in% "Trace0", "Pval"] <= 0.10
  TVECM <- HStesT$Pval <= 0.10
  H1vs2 <- HtesT[HtesT$name %in% "H1vs2", "Pval"]
  H1vs3 <- HtesT[HtesT$name %in% "H1vs3", "Pval"]

  model <- NULL

  if (is.null(model) && O1 %in% TRUE && Coint %in% TRUE && TVECM %in% TRUE && H1vs3 <= 0.10) {
    model <- "TVECM3"
  }
  if (is.null(model) && O1 %in% TRUE && Coint %in% TRUE && TVECM %in% TRUE && H1vs2 <= 0.10) {
    model <- "TVECM2"
  }
  if (is.null(model) && O1 %in% TRUE && Coint %in% TRUE && H1vs2 <= 0.10) {
    model <- "AVECM"
  }
  if (is.null(model) && O1 %in% TRUE && Coint %in% TRUE) {
    model <- "VECM"
  }

  if (is.null(model) && O1 %in% TRUE && Coint %in% FALSE && H1vs3 <= 0.10) {
    model <- "DTVARDL3"
  }
  if (is.null(model) && O1 %in% TRUE && Coint %in% FALSE && H1vs2 <= 0.10) {
    model <- "DTVARDL2"
  }
  if (is.null(model) && O1 %in% TRUE && Coint %in% FALSE) {
    model <- "DVARDL"
  }

  if (is.null(model) && O1 %in% FALSE && Coint %in% FALSE && H1vs3 <= 0.10) {
    model <- "TVARDL3"
  }
  if (is.null(model) && O1 %in% FALSE && Coint %in% FALSE && H1vs2 <= 0.10) {
    model <- "TVARDL2"
  }
  if (is.null(model) && O1 %in% FALSE && Coint %in% FALSE) {
    model <- "VARDL"
  }

  res <- Granger
  res[["URoot"]]  <- URoot
  res[["TRACE"]]  <- TRACE
  res[["HStesT"]] <- HStesT
  res[["HtesT"]]  <- HtesT
  res[["model"]]  <- model
  res[["lag"]]    <- selected_lag
  res[["cor"]]    <- cor(data[pair])
  res[["cov"]]    <- cov(data[pair])

  return(res)
}
