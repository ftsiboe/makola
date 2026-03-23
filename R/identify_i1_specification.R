# ==============================================================================
# Unit root screening for paired time series
#
# Purpose:
# Evaluate a pair of time series using multiple unit root tests in levels and
# first differences, then retain deterministic specifications that are
# consistent with both series being integrated of order one, I(1).
#
# Main functions:
# - identify_i1_specification()
# - root_test_specification()
#
# Notes:
# - This version stays close to your original two-function design.
# - It evaluates ADF, ERS, KPSS, PP, and ZA tests.
# - For ADF/ERS/PP/ZA, I(1)-consistent evidence means:
#     fail to reject in levels, reject in first differences.
# - For KPSS, where the null is stationarity, the logic is reversed:
#     reject in levels, fail to reject in first differences.
# ==============================================================================

#' Identify an I(1)-consistent deterministic specification for a series pair
#'
#' Applies a collection of unit root tests to a pair of time series in levels
#' and first differences, then retains test/specification combinations that are
#' consistent with both series being integrated of order one, I(1). The
#' function then selects a preferred deterministic specification among `"none"`,
#' `"const"`, and `"both"` based on the surviving evidence.
#'
#' @param pair Character or integer vector of length 2 identifying the two
#'   columns in `data` to test.
#' @param data A data frame or object coercible to a data frame containing the
#'   time-series variables.
#'
#' @return A data frame containing the retained unit root test results for the
#'   selected specification. The output includes:
#' \describe{
#'   \item{ur}{Unit root test family.}
#'   \item{series}{Series identifier (`"I1"` or `"I2"`).}
#'   \item{AdfL}{Test statistic for the level series.}
#'   \item{AdfD}{Test statistic for the differenced series.}
#'   \item{AdfLp}{Approximate p-value category for the level series.}
#'   \item{AdfDp}{Approximate p-value category for the differenced series.}
#'   \item{model}{Selected deterministic specification.}
#'   \item{lag}{Lag used by the test when applicable.}
#'   \item{O1}{Indicator that the evidence is consistent with I(1).}
#' }
#'
#' @details
#' The function evaluates multiple test families, including:
#' \describe{
#'   \item{ADF}{Augmented Dickey-Fuller test.}
#'   \item{ERS}{Elliott-Rothenberg-Stock test.}
#'   \item{KPSS}{Kwiatkowski-Phillips-Schmidt-Shin stationarity test.}
#'   \item{PP}{Phillips-Perron test.}
#'   \item{ZA}{Zivot-Andrews structural break unit root test.}
#' }
#'
#' For ADF-, ERS-, PP-, and ZA-type tests, evidence is treated as I(1)-consistent
#' when the level series appears nonstationary and the differenced series appears
#' stationary. For KPSS, where the null is stationarity, the logic is reversed.
#'
#' After pooling across tests, the function keeps only specifications for which
#' both series are classified as I(1), then prioritizes deterministic terms in
#' the order `"none"`, `"const"`, and `"both"`.
#'
#' @examples
#' \dontrun{
#' toy_data <- data.frame(
#'   market_1 = cumsum(rnorm(200)),
#'   market_2 = cumsum(rnorm(200))
#' )
#'
#' identify_i1_specification(
#'   pair = c("market_1", "market_2"),
#'   data = toy_data
#' )
#' }
#'
#' @export
identify_i1_specification <- function(pair, data) {

  x <- as.data.frame(data)[, pair, drop = FALSE]
  x <- x[stats::complete.cases(x), , drop = FALSE]

  #-----------------------------------------------
  # Augmented Dickey-Fuller Unit Root Test
  # Null hypothesis: non-stationarity
  ur_df <- rbind(
    data.frame(type = c("none", "drift", "trend"), selectlags = "AIC"),
    data.frame(type = c("none", "drift", "trend"), selectlags = "BIC")
  )
  ur_df <- as.data.frame(data.table::rbindlist(
    lapply(
      seq_len(nrow(ur_df)),
      root_test_specification,
      alt = ur_df,
      ur_test = "df",
      x = x
    ),
    fill = TRUE
  ))
  ur_df$model <- ifelse(ur_df$model %in% "drift", "const", ur_df$model)
  ur_df$O1 <- as.numeric(ur_df$AdfLp > 0.1 & ur_df$AdfDp <= 0.1)
  ur_df$ur <- "ur_df"

  #-----------------------------------------------
  # Elliott, Rothenberg & Stock test
  ur_ers <- rbind(
    data.frame(type = c("DF-GLS", "P-test"), model = "constant"),
    data.frame(type = c("DF-GLS", "P-test"), model = "trend")
  )
  ur_ers <- as.data.frame(data.table::rbindlist(
    lapply(
      seq_len(nrow(ur_ers)),
      root_test_specification,
      alt = ur_ers,
      ur_test = "ers",
      x = x
    ),
    fill = TRUE
  ))
  ur_ers$model <- ifelse(ur_ers$model %in% "constant", "const", ur_ers$model)
  ur_ers$O1 <- as.numeric(ur_ers$AdfLp > 0.1 & ur_ers$AdfDp <= 0.1)
  ur_ers$ur <- "ur_ers"

  #-----------------------------------------------
  # KPSS test
  # Null hypothesis: stationarity
  ur_kpss <- rbind(
    data.frame(lags = c("short", "long", "nil"), type = "mu"),
    data.frame(lags = c("short", "long", "nil"), type = "tau")
  )
  ur_kpss <- as.data.frame(data.table::rbindlist(
    lapply(
      seq_len(nrow(ur_kpss)),
      root_test_specification,
      alt = ur_kpss,
      ur_test = "kpss",
      x = x
    ),
    fill = TRUE
  ))
  ur_kpss$O1 <- as.numeric(ur_kpss$AdfLp <= 0.1 & ur_kpss$AdfDp > 0.1)
  ur_kpss$model <- ifelse(ur_kpss$model %in% "mu", "const", ur_kpss$model)
  ur_kpss$model <- ifelse(ur_kpss$model %in% "tau", "trend", ur_kpss$model)
  ur_kpss$ur <- "ur_kpss"

  #-----------------------------------------------
  # Phillips-Perron test
  ur_pp <- rbind(
    data.frame(type = c("Z-alpha", "Z-tau"), model = "constant"),
    data.frame(type = c("Z-alpha", "Z-tau"), model = "trend")
  )
  ur_pp <- rbind(
    data.frame(ur_pp, lags = "short"),
    data.frame(ur_pp, lags = "long")
  )
  ur_pp <- as.data.frame(data.table::rbindlist(
    lapply(
      seq_len(nrow(ur_pp)),
      root_test_specification,
      alt = ur_pp,
      ur_test = "pp",
      x = x
    ),
    fill = TRUE
  ))
  ur_pp$O1 <- as.numeric(ur_pp$AdfLp > 0.1 & ur_pp$AdfDp <= 0.1)
  ur_pp$ur <- "ur_pp"

  #-----------------------------------------------
  # Zivot-Andrews test
  ur_za <- data.frame(model = c("intercept", "trend", "both"), lag = 6)
  ur_za <- as.data.frame(data.table::rbindlist(
    lapply(
      seq_len(nrow(ur_za)),
      root_test_specification,
      alt = ur_za,
      ur_test = "za",
      x = x
    ),
    fill = TRUE
  ))
  ur_za$O1 <- as.numeric(ur_za$AdfLp > 0.1 & ur_za$AdfDp <= 0.1)
  ur_za$ur <- "ur_za"

  #-----------------------------------------------
  # Cleanup and choose preferred deterministic specification
  ur_all <- rbind(
    ur_df[c("ur", "id", "series", "AdfL", "AdfD", "AdfLp", "AdfDp", "model", "lag", "O1")],
    ur_ers[c("ur", "id", "series", "AdfL", "AdfD", "AdfLp", "AdfDp", "model", "lag", "O1")],
    ur_kpss[c("ur", "id", "series", "AdfL", "AdfD", "AdfLp", "AdfDp", "model", "lag", "O1")],
    ur_pp[c("ur", "id", "series", "AdfL", "AdfD", "AdfLp", "AdfDp", "model", "lag", "O1")],
    ur_za[c("ur", "id", "series", "AdfL", "AdfD", "AdfLp", "AdfDp", "model", "lag", "O1")]
  )

  O1 <- ur_all[c("ur", "id", "model", "lag", "series", "O1")] |>
    tidyr::spread(series, O1)

  ur_all <- dplyr::inner_join(ur_all, O1, by = c("ur", "id", "model", "lag"))
  ur_all$O1 <- ur_all$I1 %in% 1 & ur_all$I2 %in% 1
  ur_all <- ur_all[ur_all$O1 %in% 1, ]

  if (nrow(ur_all) == 0) {
    return(ur_all)
  }

  if (any(ur_all$model %in% "none" & ur_all$O1 %in% 1)) {
    ur_all <- ur_all[ur_all$model %in% "none", ]
    ur_all$model <- "none"
  } else if (any(ur_all$model %in% c("const", "constant", "intercept") & ur_all$O1 %in% 1)) {
    ur_all <- ur_all[ur_all$model %in% c("const", "constant", "intercept"), ]
    ur_all$model <- "const"
  } else if (any(ur_all$model %in% c("trend", "both") & ur_all$O1 %in% 1)) {
    ur_all <- ur_all[ur_all$model %in% c("trend", "both"), ]
    ur_all$model <- "both"
  } else {
    ur_all <- ur_all[ur_all$model %in% "none", ]
    ur_all$model <- "none"
  }

  ur_all <- unique(ur_all[names(ur_all)[!names(ur_all) %in% "id"]])

  return(ur_all)
}

#' Run one unit root test specification for both series in a pair
#'
#' Internal helper that applies a single unit root test configuration to each
#' series in a two-series dataset, both in levels and in first differences, and
#' returns a pooled summary of test statistics, approximate p-value categories,
#' deterministic specification, and lag selection.
#'
#' @param j Integer index identifying the row of `alt` to use.
#' @param alt Data frame of alternative test specifications.
#' @param ur_test Character scalar identifying the test family. Supported values
#'   are `"df"`, `"ers"`, `"kpss"`, `"pp"`, `"sp"`, and `"za"`.
#' @param x A two-column data frame or matrix containing the paired series after
#'   removing incomplete observations.
#'
#' @return A data frame with one row per series summarizing the selected test
#'   specification and the corresponding results for levels and first
#'   differences. Returns `NULL` if the test call fails.
#'
#' @details
#' This function is designed for internal use by
#' `identify_i1_specification()`. It applies the requested test family under one
#' candidate deterministic specification, computes the level and differenced
#' test statistics, and maps the results to coarse p-value categories.
#'
#' The returned columns include:
#' \describe{
#'   \item{id}{Row index of the specification in `alt`.}
#'   \item{series}{Series identifier (`"I1"` or `"I2"`).}
#'   \item{AdfL}{Test statistic for the level series.}
#'   \item{AdfD}{Test statistic for the differenced series.}
#'   \item{AdfLp}{Approximate p-value category for the level series.}
#'   \item{AdfDp}{Approximate p-value category for the differenced series.}
#'   \item{model}{Deterministic specification label.}
#'   \item{lag}{Lag used by the test when applicable.}
#'   \item{...}{Columns carried over from the specification row in `alt`.}
#' }
#'
#' @keywords internal
root_test_specification <- function(j, alt, ur_test = "df", x) {
  DONE <- NULL

  tryCatch({
    if (ur_test %in% "df") {
      selection <- alt[j, c("type", "selectlags")]
      model <- alt$type[j]
    }
    if (ur_test %in% "ers") {
      selection <- alt[j, c("type")]
      model <- alt$model[j]
    }
    if (ur_test %in% "kpss") {
      selection <- alt[j, c("type", "lags")]
      model <- alt$type[j]
    }
    if (ur_test %in% "pp") {
      selection <- alt[j, c("type", "lags")]
      model <- alt$model[j]
    }
    if (ur_test %in% "sp") {
      selection <- alt[j, c("type", "pol.deg", "signif")]
      model <- alt$type[j]
    }
    if (ur_test %in% "za") {
      selection <- data.frame(x = NA)
      model <- alt$model[j]
    }

    for (i in 1:2) {
      ADF <- list()

      if (ur_test %in% "df") {
        ADF[[1]] <- urca::ur.df(x[, i], type = alt$type[j], selectlags = alt$selectlags[j])
        lag <- ADF[[1]]@lags
        ADF[[2]] <- urca::ur.df(diff(x[, i]), type = alt$type[j], lags = ADF[[1]]@lags)
      }
      if (ur_test %in% "ers") {
        ADF[[1]] <- urca::ur.ers(x[, i], type = alt$type[j], model = alt$model[j], lag.max = 6)
        lag <- ADF[[1]]@lag
        ADF[[2]] <- urca::ur.ers(diff(x[, i]), type = alt$type[j], model = alt$model[j], lag.max = 6)
      }
      if (ur_test %in% "kpss") {
        ADF[[1]] <- urca::ur.kpss(x[, i], type = alt$type[j], lags = alt$lags[j])
        lag <- ADF[[1]]@lag[1]
        ADF[[2]] <- urca::ur.kpss(diff(x[, i]), type = alt$type[j], lags = alt$lags[j])
      }
      if (ur_test %in% "pp") {
        ADF[[1]] <- urca::ur.pp(x[, i], type = alt$type[j], model = alt$model[j], lags = alt$lags[j])
        lag <- ADF[[1]]@lag[1]
        ADF[[2]] <- urca::ur.pp(diff(x[, i]), type = alt$type[j], model = alt$model[j], use.lag = ADF[[1]]@lag[1])
      }
      if (ur_test %in% "sp") {
        ADF[[1]] <- urca::ur.sp(x[, i], type = alt$type[j], pol.deg = alt$pol.deg[j], signif = alt$signif[j])
        lag <- NA
        ADF[[2]] <- urca::ur.sp(diff(x[, i]), type = alt$type[j], pol.deg = alt$pol.deg[j], signif = alt$signif[j])
      }
      if (ur_test %in% "za") {
        ADF[[1]] <- urca::ur.za(x[, i], model = alt$model[j], lag = alt$lag[j])
        lag <- ADF[[1]]@lag[1]
        ADF[[2]] <- urca::ur.za(diff(x[, i]), model = alt$model[j], lag = ADF[[1]]@lag[1])
      }

      ADFp <- list()
      for (k in 1:2) {
        cval <- abs(ifelse(
          ur_test %in% "sp",
          ADF[[k]]@cval[1],
          ifelse(ur_test %in% "za", ADF[[k]]@cval[3], ADF[[k]]@cval[1, 3])
        ))

        ADFp[[k]] <- ifelse(
          abs(ADF[[k]]@teststat[1]) > cval, 0.01,
          ifelse(
            abs(ADF[[k]]@teststat[1]) > cval, 0.05,
            ifelse(abs(ADF[[k]]@teststat[1]) > cval, 0.10, 1)
          )
        )
      }

      DONE[[length(DONE) + 1]] <- data.frame(
        id = j,
        series = paste0("I", i),
        AdfL = ADF[[1]]@teststat[1],
        AdfD = ADF[[2]]@teststat[1],
        AdfLp = ADFp[[1]],
        AdfDp = ADFp[[2]],
        model = model,
        lag = lag,
        selection
      )
    }

    DONE <- as.data.frame(data.table::rbindlist(DONE, fill = TRUE))
  }, error = function(e) {})

  return(DONE)
}
