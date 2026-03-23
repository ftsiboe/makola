#' Simulate synthetic market price series with a tunable integrated subgroup
#'
#' Generates synthetic market price series over a user-defined date range and
#' temporal resolution. A selected subset of markets shares a common stochastic
#' trend, producing price series that are integrated through a common long-run
#' component by construction, while the remaining markets evolve independently.
#'
#' The function is intended for tool development, testing, demonstrations, and
#' examples involving market integration and related time-series workflows.
#'
#' @param start_date A valid start date coercible to `Date`.
#' @param end_date A valid end date coercible to `Date`.
#' @param frequency Character scalar giving the temporal resolution of the
#'   simulated data. Supported values are `"day"`, `"daily"`, `"week"`,
#'   `"weekly"`, `"month"`, `"monthly"`, `"quarter"`, `"quarterly"`, `"year"`,
#'   and `"yearly"`.
#' @param n_markets Positive integer giving the total number of market price
#'   series to simulate.
#' @param n_integrated_markets Nonnegative integer giving the number of markets
#'   that share a common stochastic trend. This value cannot exceed
#'   `n_markets`.
#' @param seed Optional integer seed for reproducibility. If `NULL`, no seed is
#'   set.
#' @param base_price Numeric scalar giving the base price level used as the
#'   starting intercept for the first market. Subsequent markets are shifted
#'   upward by a fixed increment.
#' @param trend_drift Numeric scalar giving the drift term for the stochastic
#'   trend process.
#' @param trend_sd Numeric scalar giving the standard deviation of innovations
#'   in the stochastic trend process.
#' @param stationary_ar Numeric scalar giving the autoregressive coefficient
#'   used to generate stationary deviations around the common trend for
#'   integrated markets.
#' @param stationary_sd Numeric scalar giving the standard deviation of the
#'   stationary deviation process for integrated markets.
#' @param noise_sd Numeric scalar giving the standard deviation of additional
#'   idiosyncratic noise added to each market series.
#' @param output Character scalar specifying the output format. Must be one of
#'   `"wide"` or `"long"`.
#'
#' @details
#' Integrated markets are simulated using a shared cumulative stochastic trend,
#' market-specific loadings on that trend, a stationary autoregressive
#' deviation, and idiosyncratic noise. Non-integrated markets are simulated
#' using independent stochastic trends plus idiosyncratic noise.
#'
#' In `"wide"` output, the returned object includes one row per date and one
#' column per market. In `"long"` output, the data are reshaped to three
#' columns: `date`, `market`, and `price`.
#'
#' The returned object includes the following attributes:
#' \describe{
#'   \item{integrated_markets}{Character vector naming the markets simulated
#'   with a shared stochastic trend.}
#'   \item{non_integrated_markets}{Character vector naming the markets simulated
#'   independently.}
#'   \item{frequency}{Character scalar giving the normalized temporal increment
#'   used in `seq.Date()`.}
#' }
#'
#' @return
#' A data frame of simulated market prices. If `output = "wide"`, the result
#' contains a `date` column and one price column per market. If
#' `output = "long"`, the result contains the columns:
#' \describe{
#'   \item{date}{Date index for the simulated observation.}
#'   \item{market}{Market identifier.}
#'   \item{price}{Simulated market price.}
#' }
#' @export
simulate_market_prices <- function(
    start_date,
    end_date,
    frequency = "month",
    n_markets = 2,
    n_integrated_markets = 2,
    seed = NULL,
    base_price = 100,
    trend_drift = 0.2,
    trend_sd = 1,
    stationary_ar = 0.6,
    stationary_sd = 0.8,
    noise_sd = 0.25,
    output = c("wide", "long")
) {
  output <- match.arg(output)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)

  if (is.na(start_date) || is.na(end_date)) {
    stop("`start_date` and `end_date` must be valid dates.")
  }

  if (end_date < start_date) {
    stop("`end_date` must be on or after `start_date`.")
  }

  if (!is.numeric(n_markets) || n_markets < 1 || n_markets %% 1 != 0) {
    stop("`n_markets` must be a positive integer.")
  }

  if (!is.numeric(n_integrated_markets) ||
      n_integrated_markets < 0 ||
      n_integrated_markets %% 1 != 0) {
    stop("`n_integrated_markets` must be a nonnegative integer.")
  }

  if (n_integrated_markets > n_markets) {
    stop("`n_integrated_markets` cannot exceed `n_markets`.")
  }

  # Map common labels to seq.Date() increments
  by_map <- c(
    day = "day",
    daily = "day",
    week = "week",
    weekly = "week",
    month = "month",
    monthly = "month",
    quarter = "quarter",
    quarterly = "quarter",
    year = "year",
    yearly = "year"
  )

  frequency_key <- tolower(frequency)

  if (!frequency_key %in% names(by_map)) {
    stop("Unsupported `frequency`. Use one of: day, week, month, quarter, year.")
  }

  dates <- seq.Date(
    from = start_date,
    to   = end_date,
    by   = by_map[[frequency_key]]
  )

  n_periods <- length(dates)

  if (n_periods < 2) {
    stop("The chosen date range and frequency produce fewer than 2 observations.")
  }

  market_names <- paste0("market_", seq_len(n_markets))

  out <- data.frame(date = dates)

  # Shared stochastic trend for integrated markets
  if (n_integrated_markets > 0) {
    common_trend <- cumsum(stats::rnorm(n_periods, mean = trend_drift, sd = trend_sd))

    for (i in seq_len(n_integrated_markets)) {
      stationary_dev <- as.numeric(
        stats::arima.sim(
          model = list(ar = stationary_ar),
          n = n_periods,
          sd = stationary_sd
        )
      )

      idiosyncratic_noise <- stats::rnorm(n_periods, mean = 0, sd = noise_sd)

      intercept_i <- base_price + (i - 1) * 5
      loading_i   <- 0.8 + 0.2 * i

      out[[market_names[i]]] <-
        intercept_i +
        loading_i * common_trend +
        stationary_dev +
        idiosyncratic_noise
    }
  }

  # Independent non-cointegrated markets
  if (n_integrated_markets < n_markets) {
    for (i in (n_integrated_markets + 1):n_markets) {
      own_trend <- cumsum(stats::rnorm(n_periods, mean = trend_drift, sd = trend_sd))
      idiosyncratic_noise <- stats::rnorm(n_periods, mean = 0, sd = noise_sd)

      intercept_i <- base_price + (i - 1) * 5

      out[[market_names[i]]] <-
        intercept_i +
        own_trend +
        idiosyncratic_noise
    }
  }

  attr(out, "integrated_markets") <- market_names[seq_len(n_integrated_markets)]
  attr(out, "non_integrated_markets") <- market_names[seq.int(n_integrated_markets + 1, n_markets)]
  attr(out, "frequency") <- by_map[[frequency_key]]

  if (output == "long") {
    out <- tidyr::pivot_longer(
      out,
      cols = -date,
      names_to = "market",
      values_to = "price"
    )
  }

  out
}
