#' Estimate the selected integration-based model for a pair of time series
#'
#' Selects an integration-based model for a pair of time series using
#' `select_pair_integration_model()`, then attempts to estimate the selected
#' model. If estimation fails, the function progressively falls back to simpler
#' nested model classes until estimation succeeds.
#'
#' @param pair Character or integer vector of length 2 identifying the two
#'   columns in `data` to use.
#' @param data A data frame or object coercible to a data frame containing the
#'   time-series variables.
#' @param nboot Positive integer giving the number of bootstrap replications
#'   used in threshold-based model selection tests. Default is `100`.
#' @param lag.max Positive integer giving the maximum lag order considered in
#'   `select_pair_integration_model()`. Default is `2`.
#' @param lag.opt Positive integer giving the user-imposed lower bound on the
#'   lag order used in the final testing sequence. Default is `2`.
#' @param trim Numeric scalar in `(0, 0.5)` giving the trimming proportion used
#'   in threshold cointegration testing. Default is `0.10`.
#' @param marketPairs Optional data frame used to relabel `X1` and `X2` in the
#'   final output. If supplied, it must contain columns `"P1"` and `"P2"`.
#' @param mpair Optional row index or row name identifying the relevant row of
#'   `marketPairs`.
#'
#' @return A list returned by `select_pair_integration_model()` with added
#'   components:
#' \describe{
#'   \item{est}{A data frame of estimation results from the successfully fitted
#'   model, when estimation succeeds.}
#'   \item{estimation_error}{Character vector of captured estimation error
#'   messages, if any estimation attempts fail.}
#' }
#'
#' The returned list may also contain an updated `model` value if the function
#' falls back to a simpler specification during estimation.
#'
#' @details
#' The function first applies `select_pair_integration_model()` to determine the
#' preferred model class. It then attempts estimation using the corresponding
#' estimator:
#' \describe{
#'   \item{`"TVECM3"`, `"TVECM2"`, `"AVECM"`}{`estimate_threshold_vecm()`}
#'   \item{`"VECM"`}{`estimate_linear_vecm()`}
#'   \item{`"DTVARDL2"`, `"DTVARDL3"`, `"TVARDL2"`, `"TVARDL3"`}{`estimate_threshold_vardl()`}
#'   \item{`"VARDL"`, `"DVARDL"`}{`estimate_linear_vardl()`}
#' }
#'
#' If estimation fails, the function attempts fallback models in the following
#' order:
#' \enumerate{
#'   \item `TVECM3 -> TVECM2 -> AVECM -> VECM`
#'   \item `TVARDL3 -> TVARDL2 -> VARDL`
#'   \item `DTVARDL3 -> DTVARDL2 -> DVARDL`
#' }
#'
#' Error messages from failed estimation attempts are stored in
#' `res$estimation_error` rather than being silently discarded.
#'
#' If `marketPairs` and `mpair` are supplied, the function relabels `X1` and
#' `X2` in the estimation output using the corresponding market names. This
#' relabeling is applied only if `res$est` exists and contains `X1` and `X2`.
#' @export
estimate_pair_integration_model <- function(
    pair,
    data,
    nboot = 100,
    lag.max = 2,
    lag.opt = 2,
    trim = 0.10,
    marketPairs = NULL,
    mpair = NULL
) {

  #---------------------------------------------------------------------------
  # Step 1: Select the preferred integration-based model for the series pair
  #---------------------------------------------------------------------------
  res <- select_pair_integration_model(
    pair = pair,
    data = data,
    nboot = nboot,
    lag.max = lag.max,
    lag.opt = lag.opt,
    trim = trim
  )

  # Store error messages from failed estimation attempts
  res$estimation_error <- character(0)

  #---------------------------------------------------------------------------
  # Internal helper to estimate a given model and capture any error message
  #---------------------------------------------------------------------------
  attempt_estimation <- function(res, target_model) {
    out <- res
    out$model <- target_model

    tryCatch({
      if (target_model %in% c("TVECM3", "TVECM2", "AVECM")) {
        out$est <- estimate_threshold_vecm(out)
      }

      if (target_model %in% c("VECM")) {
        out$est <- estimate_linear_vecm(out)
      }

      if (target_model %in% c("DTVARDL2", "DTVARDL3", "TVARDL3", "TVARDL2")) {
        out$est <- estimate_threshold_vardl(out)
      }

      if (target_model %in% c("VARDL", "DVARDL")) {
        out$est <- estimate_linear_vardl(out)
      }

      return(list(success = TRUE, res = out, error = NULL))
    }, error = function(e) {
      return(list(
        success = FALSE,
        res = out,
        error = paste0(target_model, ": ", conditionMessage(e))
      ))
    })
  }

  #---------------------------------------------------------------------------
  # Step 2: Build fallback path based on the initially selected model
  #---------------------------------------------------------------------------
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

  #---------------------------------------------------------------------------
  # Step 3: Try selected model, then fallback models if needed
  #---------------------------------------------------------------------------
  CHECK <- FALSE

  for (candidate_model in fallback_models) {
    attempt <- attempt_estimation(res = res, target_model = candidate_model)

    if (attempt$success) {
      res <- attempt$res
      CHECK <- TRUE
      break
    } else {
      res$estimation_error <- c(res$estimation_error, attempt$error)
    }
  }

  if (!CHECK) {
    warning("All estimation attempts failed. See `res$estimation_error` for details.")
  }

  #---------------------------------------------------------------------------
  # Step 4: Optional relabeling of X1 and X2 using external market-pair lookup
  # Apply only if estimation output exists and contains X1 and X2.
  #---------------------------------------------------------------------------
  # if (!is.null(marketPairs) &&
  #     !is.null(mpair) &&
  #     !is.null(res$est) &&
  #     all(c("X1", "X2") %in% names(res$est))) {
  #
  #   res$est$X1 <- ifelse(
  #     res$est$X1 %in% "P1",
  #     marketPairs[mpair, "P1"],
  #     marketPairs[mpair, "P2"]
  #   )
  #
  #   res$est$X2 <- ifelse(
  #     res$est$X2 %in% "P1",
  #     marketPairs[mpair, "P1"],
  #     marketPairs[mpair, "P2"]
  #   )
  # }

  return(res)
}
