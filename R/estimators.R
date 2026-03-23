#' Estimate a threshold or asymmetric vector error correction model
#'
#' Estimates a threshold vector error correction model (TVECM) or asymmetric
#' vector error correction model (AVECM) from a pre-screened model-selection
#' result object. The function fits the selected model, extracts coefficient and
#' variance-covariance information, computes nonlinear error-correction and
#' adjustment summaries, and returns a formatted result table.
#'
#' @param res A list-like object produced by the model-selection workflow. It is
#'   expected to contain at least:
#'   \describe{
#'     \item{`P1`, `P2`}{The ordered series used in estimation.}
#'     \item{`URoot`}{Unit root screening results containing a `model` field.}
#'     \item{`lag`}{Lag order used in estimation.}
#'     \item{`model`}{Selected model class, such as `"TVECM2"`, `"TVECM3"`, or
#'     `"AVECM"`.}
#'     \item{`names`}{Identifiers for the ordered pair of series.}
#'   }
#'
#' @return A data frame containing estimated coefficients, nonlinear
#'   combinations, regime summaries, and model metadata. The returned columns
#'   include:
#'   \describe{
#'     \item{`Model`}{Selected model class.}
#'     \item{`name`}{Name of the estimated coefficient or derived statistic.}
#'     \item{`Estimate`}{Point estimate.}
#'     \item{`SE`}{Standard error when available.}
#'     \item{`Tvalue`}{t-statistic when available.}
#'     \item{`Pval`}{Two-sided p-value when available.}
#'   }
#'
#' @details
#' The deterministic specification is chosen from `res$URoot$model`, with
#' `"trend"` recoded to `"both"` for `tsDyn::TVECM()`. The function then:
#' \enumerate{
#'   \item fits a two- or three-regime TVECM, or an AVECM with threshold fixed
#'   at zero,
#'   \item extracts coefficients and their variance-covariance matrix using
#'   `extract_tsdyn_vcov()`,
#'   \item computes regime-specific error-correction summaries and symmetry
#'   measures using `compute_nlcom_delta()`,
#'   \item computes regime-specific adjustment-length summaries across
#'   percentages from 1 to 99,
#'   \item appends regime shares, thresholds, transaction-cost proxies, and
#'   log-likelihood.
#' }
#' @export
estimate_threshold_vecm <- function(res){
  cat(crayon::black("Threshold Vector Error Correction model (TVECM)"), fill = TRUE)

  include <- ifelse(unique(res$URoot$model) %in% "trend", "both", unique(res$URoot$model))
  x <- data.frame(P1 = res$P1, P2 = res$P2)

  Pmean     <- mean(c(x[, 1], x[, 2]), na.rm = TRUE)
  Pmean_log <- mean(c(exp(x[, 1]), exp(x[, 2])), na.rm = TRUE)

  nthresh <- 1
  if (res$model %in% c("TVECM3")) {
    nthresh <- 2
  }

  if (res$model %in% c("TVECM3", "TVECM2")) {
    capture.output(
      TVECM <- tsDyn::TVECM(
        x,
        nthresh = nthresh,
        lag = res$lag,
        common = "All",
        include = include,
        plot = FALSE,
        ngridBeta = 200,
        ngridTh = 200,
        trim = 0.10
      )
    )
  }

  if (res$model %in% "AVECM") {
    capture.output(
      TVECM <- tsDyn::TVECM(
        x,
        nthresh = nthresh,
        lag = res$lag,
        common = "All",
        include = include,
        plot = FALSE,
        th1 = list(exact = 0),
        trim = 0.10
      )
    )
  }

  list2env(extract_tsdyn_vcov(TVECM), environment())

  Vht <- vcov[grepl("ECT", colnames(vcov)), grepl("ECT", colnames(vcov))]
  Tht <- coef[grepl("ECT", colnames(vcov))]

  ECTnlcom <- c(
    as.formula("~eq2_L_ECT-eq1_L_ECT"),
    as.formula("~eq2_H_ECT-eq1_H_ECT"),
    as.formula("~(eq2_L_ECT-eq1_L_ECT)/(eq2_H_ECT-eq1_H_ECT)")
  )
  ECTnlcom <- compute_nlcom_delta(func = ECTnlcom, vcMat = Vht, coefs = Tht, df = df)
  ECTnlcom$name <- c("rhoL", "rhoH", "SymmSize")

  Adjnlcom <- as.data.frame(
    data.table::rbindlist(
      lapply(
        1:99,
        function(adj, vcMat, coefs, df){
          DONE <- NULL
          Adjnlcom <- NULL
          AdjRegH  <- NULL
          AdjRegL  <- NULL
          ESTNAME  <- NULL

          if (coef["eq2_L_ECT"] - coef["eq1_L_ECT"] > -1 &
              coef["eq2_L_ECT"] - coef["eq1_L_ECT"] < 0) {
            Adjnlcom <- c(
              Adjnlcom,
              as.formula(paste0("~log(1.0001-", adj/100, ")/log(1.0001+(eq2_L_ECT-eq1_L_ECT))"))
            )
            AdjRegL <- paste0("(log(1.0001-", adj/100, ")/log(1.0001+(eq2_L_ECT-eq1_L_ECT)))")
            ESTNAME <- c(ESTNAME, paste0("AdjRegL_", adj))
          }

          if (coef["eq2_L_ECT"] - coef["eq1_L_ECT"] > 0 &
              coef["eq2_L_ECT"] - coef["eq1_L_ECT"] < 1) {
            Adjnlcom <- c(
              Adjnlcom,
              as.formula(paste0("~log(1.0001-", adj/100, ")/log(1.0001-(eq2_L_ECT-eq1_L_ECT))"))
            )
            AdjRegL <- paste0("(log(1.0001-", adj/100, ")/log(1.0001-(eq2_L_ECT-eq1_L_ECT)))")
            ESTNAME <- c(ESTNAME, paste0("AdjRegL_", adj))
          }

          if (coef["eq2_H_ECT"] - coef["eq1_H_ECT"] > -1 &
              coef["eq2_H_ECT"] - coef["eq1_H_ECT"] < 0) {
            Adjnlcom <- c(
              Adjnlcom,
              as.formula(paste0("~log(1.0001-", adj/100, ")/log(1.0001+(eq2_H_ECT-eq1_H_ECT))"))
            )
            AdjRegH <- paste0("(log(1.0001-", adj/100, ")/log(1.0001+(eq2_H_ECT-eq1_H_ECT)))")
            ESTNAME <- c(ESTNAME, paste0("AdjRegH_", adj))
          }

          if (coef["eq2_H_ECT"] - coef["eq1_H_ECT"] > 0 &
              coef["eq2_H_ECT"] - coef["eq1_H_ECT"] > 1) {
            Adjnlcom <- c(
              Adjnlcom,
              as.formula(paste0("~log(1.0001-", adj/100, ")/log(1.0001-(eq2_H_ECT-eq1_H_ECT))"))
            )
            AdjRegH <- paste0("(log(1.0001-", adj/100, ")/log(1.0001-(eq2_H_ECT-eq1_H_ECT)))")
            ESTNAME <- c(ESTNAME, paste0("AdjRegH_", adj))
          }

          if (is.null(AdjRegL) & is.null(AdjRegL)) {
            Adjnlcom <- c(Adjnlcom, as.formula(paste0("~", AdjRegL, "/", AdjRegH)))
            ESTNAME  <- c(ESTNAME, paste0("SymmLenght_", adj))
          }

          Adjnlcom <- compute_nlcom_delta(func = Adjnlcom, vcMat = Vht, coefs = Tht, df = df)
          Adjnlcom$name <- ESTNAME

          if (abs(coef["eq2_L_ECT"] - coef["eq1_L_ECT"]) > 1) {
            Adjnlcom$Estimate <- ifelse(grepl("AdjRegL", Adjnlcom$name), 0, Adjnlcom$Estimate)
            Adjnlcom$Estimate <- ifelse(
              grepl("SymmLenght", Adjnlcom$name),
              Adjnlcom[grepl("AdjRegH", Adjnlcom$name), "Estimate"] / 1,
              Adjnlcom$Estimate
            )
          }

          if (abs(coef["eq2_H_ECT"] - coef["eq1_H_ECT"]) > 1) {
            Adjnlcom$Estimate <- ifelse(grepl("AdjRegH", Adjnlcom$name), 0, Adjnlcom$Estimate)
            Adjnlcom$Estimate <- ifelse(
              grepl("SymmLenght", Adjnlcom$name),
              Adjnlcom[grepl("AdjRegL", Adjnlcom$name), "Estimate"] / 1,
              Adjnlcom$Estimate
            )
          }

          DONE <- Adjnlcom
          return(DONE)
        },
        vcMat = Vht,
        coefs = Tht,
        df = df
      ),
      fill = TRUE
    )
  )

  Obs <- sum(TVECM$nobs_regimes, na.rm = TRUE)
  Scalers <- t(data.frame(
    Beta      = -1 * TVECM$model.specific$coint[2],
    ObsLo     = TVECM$nobs_regimes[1] / Obs,
    ObsMi     = TVECM$nobs_regimes[2] / Obs,
    ObsUp     = TVECM$nobs_regimes[3] / Obs,
    TrshUp    = TVECM$model.specific$Thresh[2],
    TrshLo    = TVECM$model.specific$Thresh[1],
    Pmean     = Pmean,
    TCost     = (mean(abs(TVECM$model.specific$Thresh), na.rm = TRUE) / Pmean) * 100,
    TCost_log = (mean(abs(exp(TVECM$model.specific$Thresh)), na.rm = TRUE) / Pmean_log) * 100,
    logLik    = logLik(TVECM)
  ))

  TVECMRes <- rbind(
    data.frame(Estimate = coef, SE = diag(vcov)^0.5, name = NAMES$name),
    data.frame(Estimate = Scalers[, 1], name = rownames(Scalers), SE = NA)
  )

  TVECMRes$Tvalue <- TVECMRes$Estimate / TVECMRes$SE
  TVECMRes$Pval <- pt(abs(TVECMRes$Tvalue), df = df, lower.tail = FALSE) +
    pt(-abs(TVECMRes$Tvalue), df = df, lower.tail = TRUE)

  TVECMRes <- rbind(TVECMRes, ECTnlcom, Adjnlcom)
  TVECMRes$Model <- res$model
  TVECMRes <- data.frame(
    data.frame(t(res$names)),
    TVECMRes[c("Model", "name", "Estimate", "SE", "Tvalue", "Pval")]
  )
  rownames(TVECMRes) <- 1:nrow(TVECMRes)

  return(TVECMRes)
}


#' Estimate a linear vector error correction model
#'
#' Estimates a linear vector error correction model (VECM) from a pre-screened
#' model-selection result object. The function fits the VECM, extracts
#' coefficients and covariance information, computes nonlinear error-correction
#' and adjustment summaries, and returns a formatted result table.
#'
#' @param res A list-like object produced by the model-selection workflow. It is
#'   expected to contain at least `P1`, `P2`, `URoot`, `lag`, `model`, and
#'   `names`.
#'
#' @return A data frame of coefficient estimates and derived summaries with the
#'   columns `Model`, `name`, `Estimate`, `SE`, `Tvalue`, and `Pval`.
#'
#' @details
#' The function fits a rank-one VECM using `tsDyn::VECM()` under the selected
#' deterministic specification. It then:
#' \enumerate{
#'   \item extracts coefficients and the associated covariance matrix,
#'   \item computes nonlinear summaries for the combined error-correction term,
#'   \item computes adjustment-length summaries for percentages from 1 to 99,
#'   \item computes a symmetry measure for transmission,
#'   \item appends beta and log-likelihood summaries.
#' }
#'
#' @export
estimate_linear_vecm <- function(res){
  cat(crayon::black("Linear Vector Error Correction Model (VECM)"), fill = TRUE)

  x <- data.frame(P1 = res$P1, P2 = res$P2)
  include <- ifelse(unique(res$URoot$model) %in% "trend", "both", unique(res$URoot$model))

  capture.output(
    VECM <- tsDyn::VECM(
      x,
      lag = res$lag,
      r = 1,
      include = include,
      beta = NULL,
      estim = c("ML"),
      exogen = NULL
    )
  )

  coef <- data.frame(summary(VECM)$coefMat)[, 1]
  vcov <- summary(VECM)$sigma %x% summary(VECM)$cov.unscaled

  NAMES <- rownames(data.frame(summary(VECM)$coefMat))
  NAMES <- gsub("2-", "2_t", gsub("1-", "1_t", NAMES))
  NAMES <- gsub("P1:", "eq1_", gsub("P2:", "eq2_", NAMES))

  colnames(vcov) <- rownames(vcov) <- names(coef) <- NAMES

  Vht <- vcov[grepl("ECT", colnames(vcov)), grepl("ECT", colnames(vcov))]
  Tht <- coef[grepl("ECT", colnames(vcov))]
  df  <- VECM$df.residual

  ECTnlcom <- compute_nlcom_delta(
    func = c(as.formula("~eq2_ECT-eq1_ECT")),
    vcMat = Vht,
    coefs = Tht,
    df = df
  )
  ECTnlcom$name <- "ECT"

  Adjnlcom <- NULL
  Adjnames <- NULL

  for (adj in seq(1, 99, 1)) {
    Adjnames <- c(Adjnames, paste0("AdjRegL0_", adj))

    Adjnlcom <- c(
      Adjnlcom,
      as.formula(
        ifelse(
          coef["eq2_ECT"] - coef["eq1_ECT"] > 0,
          paste0("~log(1.0001-", adj/100, ")/log(1.0001-(eq2_ECT-eq1_ECT))"),
          paste0("~log(1.0001-", adj/100, ")/log(1.0001+(eq2_ECT-eq1_ECT))")
        )
      )
    )
  }

  Adjnlcom <- compute_nlcom_delta(func = Adjnlcom, vcMat = Vht, coefs = Tht, df = df)
  Adjnlcom$name <- Adjnames

  Symmnlcom <- compute_nlcom_delta(
    func = c(as.formula("~eq1_ECT-eq2_ECT")),
    vcMat = Vht,
    coefs = Tht,
    df = df
  )
  Symmnlcom$name <- c("Symm0")

  VECMRes <- data.frame(summary(VECM)$coefMat)[c("Estimate", "Std..Error", "t.value", "Pr...t..")]
  names(VECMRes) <- c("Estimate", "SE", "Tvalue", "Pval")
  VECMRes$name <- NAMES

  VECMRes <- rbind(
    VECMRes,
    ECTnlcom,
    Adjnlcom,
    Symmnlcom,
    data.frame(Estimate = -1 * VECM$model.specific$beta[2, 1], name = "Beta", SE = NA, Pval = NA, Tvalue = NA),
    data.frame(Estimate = logLik(VECM), name = "logLik", SE = NA, Pval = NA, Tvalue = NA)
  )

  VECMRes$Model <- res$model
  VECMRes <- VECMRes[c("Model", "name", "Estimate", "SE", "Tvalue", "Pval")]
  rownames(VECMRes) <- 1:nrow(VECMRes)
  VECMRes <- data.frame(data.frame(t(res$names)), VECMRes[c("Model", "name", "Estimate", "SE", "Tvalue", "Pval")])

  return(VECMRes)
}


#' Estimate a threshold vector autoregressive distributed lag model
#'
#' Estimates a threshold vector autoregressive distributed lag model (TVARDL)
#' from a pre-screened model-selection result object. Because the model is not
#' estimated directly as a system in the original workflow, the function first
#' fits a threshold VAR and then uses its generated regressors to estimate the
#' associated TVARDL system with `systemfit`.
#'
#' @param res A list-like object produced by the model-selection workflow. It is
#'   expected to contain at least `P1`, `P2`, `lag`, `model`, and `names`.
#'
#' @return A data frame of coefficient estimates and derived summaries with the
#'   columns `Model`, `name`, `Estimate`, `SE`, `Tvalue`, and `Pval`.
#'
#' @details
#' Depending on `res$model`, the function fits either a two-regime or
#' three-regime threshold model, optionally using differenced data for
#' `"DTVARDL2"` and `"DTVARDL3"`. It then:
#' \enumerate{
#'   \item fits a threshold VAR using `tsDyn::TVAR()`,
#'   \item extracts the generated regressors and reparameterizes them into a
#'   TVARDL system estimated with `systemfit::systemfit()`,
#'   \item computes regime-specific error-correction and long-run summaries,
#'   \item computes adjustment-length summaries for percentages from 1 to 99,
#'   \item appends regime shares, thresholds, transaction-cost proxies, and
#'   log-likelihood.
#' }
#'
#' @export
estimate_threshold_vardl <- function(res){
  cat(crayon::black("Threshold Vector Autoregressive distributed lag (TVARDL)"), fill = TRUE)

  x <- data.frame(P1 = res$P1, P2 = res$P2)
  Pmean     <- mean(c(x[, 1], x[, 2]), na.rm = TRUE)
  Pmean_log <- mean(c(exp(x[, 1]), exp(x[, 2])), na.rm = TRUE)

  if (res$model %in% c("TVARDL3", "DTVARDL3")) nthresh <- 2
  if (res$model %in% c("TVARDL2", "DTVARDL2")) nthresh <- 1
  if (res$model %in% c("DTVARDL2", "DTVARDL3")) x <- x - dplyr::lag(x); x <- x[complete.cases(x), ]

  capture.output(
    TVAR <- tsDyn::TVAR(x, lag = res$lag, nthresh = nthresh, thDelay = 1, plot = FALSE, mTh = 1)
  )

  list2env(extract_tsdyn_vcov(TVAR), environment())
  df <- df - (nthresh + 1)

  xx    <- TVAR$model
  NAMES <- names(TVAR$model)
  NAMES <- gsub("1 -", "1_t", gsub("2 -", "2_t", NAMES))
  NAMES <- gsub("L ", "L_", gsub("H ", "H_", NAMES))
  names(xx) <- NAMES

  xvars <- names(xx)
  xvars <- xvars[!xvars %in% c("P1", "P2")]

  if (res$model %in% c("TVARDL3", "DTVARDL3")) {
    TVARDL <- list(
      eq1 = as.formula(paste0("P1~-1+L_Intercept*P2+Intercept*P2+H_Intercept*P2-P2+", paste0(xvars, collapse = "+"))),
      eq2 = as.formula(paste0("P2~-1+L_Intercept*P1+Intercept*P1+H_Intercept*P1-P1+", paste0(xvars, collapse = "+")))
    )
  }

  if (res$model %in% c("TVARDL2", "DTVARDL2")) {
    TVARDL <- list(
      eq1 = as.formula(paste0("P1~-1+L_Intercept*P2+H_Intercept*P2-P2+", paste0(xvars, collapse = "+"))),
      eq2 = as.formula(paste0("P2~-1+L_Intercept*P1+H_Intercept*P1-P1+", paste0(xvars, collapse = "+")))
    )
  }

  TVARDL <- systemfit::systemfit(TVARDL, data = xx)

  vcov <- vcov(TVARDL)
  coef <- coef(TVARDL)

  NAMES <- names(coef)
  NAMES <- gsub("eq1_L_Intercept:P2", "eq1_L_P2_t0", NAMES)
  NAMES <- gsub("eq1_P2:Intercept",   "eq1_P2_t0", NAMES)
  NAMES <- gsub("eq1_P2:H_Intercept", "eq1_H_P2_t0", NAMES)
  NAMES <- gsub("eq2_L_Intercept:P1", "eq2_L_P1_t0", NAMES)
  NAMES <- gsub("eq2_P1:Intercept",   "eq2_P1_t0", NAMES)
  NAMES <- gsub("eq2_P1:H_Intercept", "eq2_H_P1_t0", NAMES)

  names(coef) <- colnames(vcov) <- rownames(vcov) <- NAMES

  rhoL <- paste0(
    "((-1*(1-", paste0(names(coef)[grepl("eq2_L_P2", names(coef))], collapse = "-"), "))",
    "-",
    "(-1*(1-", paste0(names(coef)[grepl("eq1_L_P1", names(coef))], collapse = "-"), ")))"
  )

  rhoH <- paste0(
    "((-1*(1-", paste0(names(coef)[grepl("eq2_H_P2", names(coef))], collapse = "-"), "))",
    "-",
    "(-1*(1-", paste0(names(coef)[grepl("eq1_H_P1", names(coef))], collapse = "-"), ")))"
  )

  ECTnlcom <- c(
    as.formula(paste0("rhoL~", rhoL)),
    as.formula(paste0("rhoH~", rhoH)),
    as.formula(paste0("SymmSize~", rhoL, "/", rhoH))
  )

  for (i in 1:2) {
    j <- ifelse(i == 1, 2, 1)
    ECTnlcom <- c(
      ECTnlcom,
      as.formula(paste0("eq", i, "_L_ECT~-1*(1-", paste0(names(coef)[grepl(paste0("eq", i, "_L_P", i), names(coef))], collapse = "-"), ")")),
      as.formula(paste0("eq", i, "_H_ECT~-1*(1-", paste0(names(coef)[grepl(paste0("eq", i, "_H_P", i), names(coef))], collapse = "-"), ")")),
      as.formula(paste0("eq", i, "_L_LR~(", paste0(names(coef)[grepl(paste0("eq", i, "_L_P", j), names(coef))], collapse = "+"), ")/",
                        paste0("(1-", paste0(names(coef)[grepl(paste0("eq", i, "_L_P", i), names(coef))], collapse = "-"), ")"))),
      as.formula(paste0("eq", i, "_H_LR~(", paste0(names(coef)[grepl(paste0("eq", i, "_H_P", j), names(coef))], collapse = "+"), ")/",
                        paste0("(1-", paste0(names(coef)[grepl(paste0("eq", i, "_H_P", i), names(coef))], collapse = "-"), ")")))
    )
  }

  eqns_lhs <- lapply(X = ECTnlcom, FUN = function(x) x[[2L]])
  eqns_rhs <- lapply(X = ECTnlcom, FUN = function(x) x[[3L]])
  vnam     <- sapply(X = eqns_lhs, FUN = as.character)

  Adjustment <- data.frame(t(coef))
  Adjustment <- as.data.frame(lapply(X = eqns_rhs, FUN = eval, envir = Adjustment))
  names(Adjustment) <- vnam
  Adjustment <- data.frame(Adjustment, t(coef))

  ECTnlcom <- NULL
  for (i in 1:length(eqns_rhs)) {
    ECTnlcom <- c(ECTnlcom, as.formula(paste0("~", eqns_rhs[i])))
  }

  ECTnlcom <- compute_nlcom_delta(func = ECTnlcom, vcMat = vcov, coefs = coef, df = TVARDL$df.residual)
  ECTnlcom$name <- vnam

  Adjnlcom <- as.data.frame(
    data.table::rbindlist(
      lapply(
        1:99,
        function(adj){
          DONE <- NULL
          tryCatch({
            if (Adjustment$rhoL > -1 & Adjustment$rhoL < 0) AdjRegL <- log(1.0001 - adj/100) / log(1.0001 + Adjustment$rhoL)
            if (Adjustment$rhoL >  0 & Adjustment$rhoL < 1) AdjRegL <- log(1.0001 - adj/100) / log(1.0001 - Adjustment$rhoL)
            if (abs(Adjustment$rhoL) > 1) AdjRegL <- 0

            if (Adjustment$rhoH > -1 & Adjustment$rhoH < 0) AdjRegH <- log(1.0001 - adj/100) / log(1.0001 + Adjustment$rhoH)
            if (Adjustment$rhoH >  0 & Adjustment$rhoH < 1) AdjRegH <- log(1.0001 - adj/100) / log(1.0001 - Adjustment$rhoH)
            if (abs(Adjustment$rhoH) > 1) AdjRegH <- 0

            DONE <- data.frame(
              Estimate = c(AdjRegL, AdjRegH, AdjRegL / AdjRegH),
              name = c(paste0("AdjRegL_", adj), paste0("AdjRegH_", adj), paste0("SymmLenght_", adj))
            )
            DONE$SE <- DONE$Tvalue <- DONE$Pval <- NA
          }, error = function(e){})
          return(DONE)
        }
      ),
      fill = TRUE
    )
  )

  Obs <- sum(TVAR$nobs_regimes, na.rm = TRUE)
  Scalers <- t(data.frame(
    Beta      = NA,
    ObsLo     = TVAR$nobs_regimes[1] / Obs,
    ObsMi     = TVAR$nobs_regimes[2] / Obs,
    ObsUp     = TVAR$nobs_regimes[3] / Obs,
    TrshUp    = TVAR$model.specific$Thresh[2],
    TrshLo    = TVAR$model.specific$Thresh[1],
    Pmean     = Pmean,
    TCost     = (mean(abs(TVAR$model.specific$Thresh), na.rm = TRUE) / Pmean) * 100,
    TCost_log = (mean(abs(exp(TVAR$model.specific$Thresh)), na.rm = TRUE) / Pmean_log) * 100,
    logLik    = as.numeric(logLik(TVARDL))
  ))

  TVARDL <- rbind(
    data.frame(Estimate = coef, SE = diag(vcov)^0.5, name = NAMES),
    data.frame(Estimate = Scalers[, 1], name = rownames(Scalers), SE = NA)
  )
  TVARDL$Tvalue <- TVARDL$Estimate / TVARDL$SE
  TVARDL$Pval <- pt(abs(TVARDL$Tvalue), df = df, lower.tail = FALSE) +
    pt(-abs(TVARDL$Tvalue), df = df, lower.tail = TRUE)

  TVARDL <- rbind(TVARDL, ECTnlcom, Adjnlcom)
  TVARDL$Model <- res$model
  TVARDL <- data.frame(
    data.frame(t(res$names)),
    TVARDL[c("Model", "name", "Estimate", "SE", "Tvalue", "Pval")]
  )
  rownames(TVARDL) <- 1:nrow(TVARDL)

  return(TVARDL)
}


#' Estimate a linear vector autoregressive distributed lag model
#'
#' Estimates a linear vector autoregressive distributed lag model (VARDL) from
#' a pre-screened model-selection result object. Because the model is not
#' estimated directly as a system in the original workflow, the function first
#' fits a linear VAR and then uses its generated regressors to estimate the
#' associated VARDL system with `systemfit`.
#'
#' @param res A list-like object produced by the model-selection workflow. It is
#'   expected to contain at least `P1`, `P2`, `lag`, `model`, and `names`.
#'
#' @return A data frame of coefficient estimates and derived summaries with the
#'   columns `Model`, `name`, `Estimate`, `SE`, `Tvalue`, and `Pval`.
#'
#' @details
#' For `"DVARDL"`, the input series are first differenced. The function then:
#' \enumerate{
#'   \item fits a linear VAR using `tsDyn::lineVar()`,
#'   \item extracts generated regressors and reparameterizes them into a VARDL
#'   system estimated with `systemfit::systemfit()`,
#'   \item computes error-correction and long-run summaries using
#'   `compute_nlcom_delta()`,
#'   \item computes symmetry and adjustment-length summaries,
#'   \item appends log-likelihood and formats the output table.
#' }
#' @export
estimate_linear_vardl <- function(res){
  cat(crayon::black("Vector Autoregressive distributed lag (VARDL)"), fill = TRUE)

  x <- data.frame(P1 = res$P1, P2 = res$P2)
  if (res$model %in% c("DVARDL")) x <- x - dplyr::lag(x); x <- x[complete.cases(x), ]

  capture.output(
    VAR <- tsDyn::lineVar(x, lag = res$lag, model = "VAR", estim = "ML", I = c("level"))
  )

  xx    <- VAR$model
  NAMES <- names(VAR$model)
  NAMES <- gsub("1 -", "1_t", gsub("2 -", "2_t", NAMES))
  names(xx) <- NAMES
  xvars <- names(xx)
  xvars <- xvars[!xvars %in% c("P1", "P2")]

  VARDL <- systemfit::systemfit(
    list(
      eq1 = as.formula(paste0("P1~-1+Intercept*P2-P2+", paste0(xvars, collapse = "+"))),
      eq2 = as.formula(paste0("P2~-1+Intercept*P1-P1+", paste0(xvars, collapse = "+")))
    ),
    data = xx
  )

  vcov <- vcov(VARDL)
  coef <- coef(VARDL)

  NAMES <- names(coef)
  NAMES <- gsub("eq1_Intercept:P2", "eq1_P2_t0", NAMES)
  NAMES <- gsub("eq2_Intercept:P1", "eq2_P1_t0", NAMES)
  names(coef) <- colnames(vcov) <- rownames(vcov) <- NAMES

  ECT <- LR <- list()
  for (i in 1:2) {
    j <- ifelse(i == 1, 2, 1)
    ECT[[length(ECT) + 1]] <- paste0("(-1*(1-", paste0(names(coef)[grepl(paste0("eq", i, "_P", i), names(coef))], collapse = "-"), "))")
    LR[[length(LR) + 1]] <- paste0(
      "(",
      paste0(names(coef)[grepl(paste0("eq", i, "_P", j), names(coef))], collapse = "+"),
      ")/",
      paste0("(1-", paste0(names(coef)[grepl(paste0("eq", i, "_P", i), names(coef))], collapse = "-"), ")")
    )
  }

  ECTNAmes <- ECTnlcom <- NULL
  for (i in 1:2) {
    ECTNAmes <- c(ECTNAmes, paste0("eq", i, "_ECT"), paste0("eq", i, "_LR"))
    ECTnlcom <- c(
      ECTnlcom,
      as.formula(paste0("~", ECT[[i]])),
      as.formula(paste0("~", LR[[i]]))
    )
  }

  ECTnlcom[[length(ECTnlcom) + 1]] <- as.formula(paste0("~(", ECT[[2]], ")-(", ECT[[1]], ")"))
  ECTNAmes <- c(ECTNAmes, "ECT")

  ECTnlcom <- compute_nlcom_delta(func = ECTnlcom, vcMat = vcov, coefs = coef, df = VAR$df.residual - 2)
  ECTnlcom$name <- ECTNAmes
  row.names(ECTnlcom) <- ECTNAmes

  Symmnlcom <- compute_nlcom_delta(
    func = c(as.formula(paste0("~", ECT[[1]], "-", ECT[[2]]))),
    vcMat = vcov,
    coefs = coef,
    df = VAR$df.residual - 2
  )
  Symmnlcom$name <- c("Symm0")

  Adjnlcom <- Adjnames <- NULL
  for (adj in seq(1, 99, 1)) {
    Adjnames <- c(Adjnames, paste0("AdjRegL0_", adj))

    Adjnlcom <- c(
      Adjnlcom,
      as.formula(
        ifelse(
          ECTnlcom["eq2_ECT", 1] - ECTnlcom["eq1_ECT", 1] > 0,
          paste0("~log(1-", adj/100, ")/log(1-(", ECT[[2]], "-", ECT[[1]], "))"),
          paste0("~log(1-", adj/100, ")/log(1+(", ECT[[2]], "-", ECT[[1]], "))")
        )
      )
    )
  }

  Adjnlcom <- compute_nlcom_delta(func = Adjnlcom, vcMat = vcov, coefs = coef, df = VAR$df.residual - 2)
  Adjnlcom$name <- Adjnames

  logLik <- logLik(VARDL)
  VARDL <- data.frame(summary(VARDL)$coefficients)
  names(VARDL) <- c("Estimate", "SE", "Tvalue", "Pval")
  VARDL <- rbind(VARDL, data.frame(Estimate = logLik, SE = NA, Pval = NA, Tvalue = NA))
  VARDL$name <- c(NAMES, "logLik")

  VARDL <- rbind(VARDL, ECTnlcom, Adjnlcom, Symmnlcom)
  VARDL$Model <- res$model
  VARDL <- VARDL[c("Model", "name", "Estimate", "SE", "Tvalue", "Pval")]

  rownames(VARDL) <- 1:nrow(VARDL)
  VARDL <- data.frame(data.frame(t(res$names)), VARDL[c("Model", "name", "Estimate", "SE", "Tvalue", "Pval")])

  return(VARDL)
}
