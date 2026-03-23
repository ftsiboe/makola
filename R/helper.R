#' Select lag order for a pair of time series
#'
#' Selects the lag order for a pair of variables using
#' `tsDyn::lags.select()`. The function subsets the requested columns,
#' removes incomplete observations, and returns the AIC-based lag choice.
#' A minimum lag order of 2 is enforced.
#'
#' @param pair Character or integer vector of length 2 identifying the two
#'   columns in `data` to use.
#' @param data A data frame or object coercible to a data frame containing the
#'   time-series variables.
#' @param lag.max Positive integer giving the maximum lag order considered.
#'   Default is `4`.
#'
#' @return An integer scalar giving the selected lag order. The returned value
#'   is at least `2`.
#'
#' @details
#' The function applies `tsDyn::lags.select()` to the selected pair of series
#' with a constant term included in the candidate models, and extracts the
#' AIC-minimizing lag order from the returned object. After selection, the lag
#' order is constrained to be no smaller than 2.
#'
#' @export
select_lag_order <- function(pair, data, lag.max = 4) {
  x <- as.data.frame(data)[, pair, drop = FALSE]
  x <- x[stats::complete.cases(x), , drop = FALSE]

  selected_lag <- tsDyn::lags.select(
    x,
    lag.max = lag.max,
    include = "const",
    fitMeasure = "LL",
    sameSample = TRUE
  )$AIC_min[2]

  selected_lag <- max(2L, as.integer(selected_lag))

  return(selected_lag)
}


#' Order a pair of series by Granger-causality strength
#'
#' Compares bidirectional Granger-causality tests for a pair of time series and
#' orders the series so that the direction with the stronger Granger-causality
#' evidence is assigned as `P2 -> P1`. The function returns the reordered series
#' along with the Granger test results in both directions.
#'
#' @param pair Character or integer vector of length 2 identifying the two
#'   columns in `data` to use.
#' @param data A data frame or object coercible to a data frame containing the
#'   time-series variables.
#' @param lag Positive integer giving the lag order used in the Granger tests.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{names}{Character vector of length 2 giving the reordered series
#'   names.}
#'   \item{P1}{Numeric vector containing the first reordered series.}
#'   \item{P2}{Numeric vector containing the second reordered series.}
#'   \item{granger}{A data frame summarizing the Granger test statistics and
#'   p-values for the directions `P2 -> P1` and `P1 -> P2`.}
#' }
#'
#' @details
#' The function first extracts the selected pair of variables and removes
#' incomplete observations. It then runs two preliminary Granger-causality
#' tests:
#' \describe{
#'   \item{Series 2 -> Series 1}{`grangertest(y = x[, 1], x = x[, 2])`}
#'   \item{Series 1 -> Series 2}{`grangertest(y = x[, 2], x = x[, 1])`}
#' }
#'
#' If the p-value for `Series 2 -> Series 1` is smaller, the original ordering
#' is retained. Otherwise, the pair is swapped. After ordering, the function
#' reruns the Granger tests and reports results in the reordered system.
#'
#' Note that this function orders series using relative statistical evidence
#' from the two p-values; it does not impose a significance threshold.
#'
#' @export
order_pair_by_granger <- function(pair, data, lag) {
  x <- as.data.frame(data)[, pair, drop = FALSE]
  x <- x[stats::complete.cases(x), , drop = FALSE]

  if (length(pair) != 2) {
    stop("`pair` must identify exactly two columns.")
  }

  if (nrow(x) == 0) {
    stop("No complete observations remain after removing missing values.")
  }

  if (!is.numeric(lag) || length(lag) != 1 || lag < 1 || lag %% 1 != 0) {
    stop("`lag` must be a positive integer.")
  }

  g1 <- lmtest::grangertest(y = x[, 1], x = x[, 2], order = lag)
  g2 <- lmtest::grangertest(y = x[, 2], x = x[, 1], order = lag)

  if (g1$`Pr(>F)`[2] < g2$`Pr(>F)`[2]) {
    P1 <- x[, 1]
    P2 <- x[, 2]
    selected_names <- c(pair[1], pair[2])
  } else {
    P1 <- x[, 2]
    P2 <- x[, 1]
    selected_names <- c(pair[2], pair[1])
  }

  granger12 <- lmtest::grangertest(y = P1, x = P2, order = lag)
  granger21 <- lmtest::grangertest(y = P2, x = P1, order = lag)

  granger <- data.frame(
    Estimate = c(granger12$F[2], granger21$F[2]),
    Pval = c(granger12$`Pr(>F)`[2], granger21$`Pr(>F)`[2]),
    name = c("granger12", "granger21")
  )

  return(list(
    names = selected_names,
    P1 = P1,
    P2 = P2,
    granger = granger
  ))
}









