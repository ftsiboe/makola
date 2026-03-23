#' Compute nonlinear combinations of estimators using the delta method
#'
#' Evaluates nonlinear functions of model coefficients using
#' `car::deltaMethod()`, then returns estimates, standard errors, t-statistics,
#' and p-values.
#'
#' @param func A named list of expressions defining the nonlinear combinations
#'   to evaluate. Each element should be an expression compatible with
#'   `car::deltaMethod()`.
#' @param vcMat Variance-covariance matrix of the coefficient estimates.
#' @param coefs Named numeric vector of coefficient estimates.
#' @param df Degrees of freedom used to compute t-based p-values.
#'
#' @return A data frame with one row per nonlinear combination and the
#'   following columns:
#' \describe{
#'   \item{Estimate}{Estimated value of the nonlinear combination.}
#'   \item{SE}{Delta-method standard error.}
#'   \item{Tvalue}{t-statistic computed as `Estimate / SE`.}
#'   \item{Pval}{Two-sided p-value based on the t distribution with `df`
#'   degrees of freedom.}
#' }
#'
#' @details
#' The function converts each supplied expression in `func` to a character
#' representation, evaluates it with `car::deltaMethod()`, and then combines
#' the resulting estimates and standard errors into a single summary table.
#'
#' P-values are computed manually using the t distribution rather than being
#' taken directly from `car::deltaMethod()`.
#'
#' @param func A named list of expressions. The names are currently not used in
#'   the returned object, but they are useful for organizing the input.
#'
#' @export
compute_nlcom_delta <- function(func, vcMat, coefs, df) {
  func_names <- names(func)

  temp <- lapply(func, function(x) {
    as.character(as.expression(x[[length(x)]]))
  })

  func <- data.frame(form = unlist(temp))

  lisRes <- apply(
    func,
    1,
    function(x) {
      car::deltaMethod(
        object = coefs,
        g = x,
        vcov. = vcMat,
        level = 0.95
      )
    }
  )

  val <- plyr::ldply(lisRes)
  val <- cbind(func, val)

  lenVal <- length(val[1, ])
  retDF <- val[, c(-(lenVal - 1), -lenVal)]
  retDF <- retDF[c("Estimate", "SE")]
  retDF$Tvalue <- retDF$Estimate / retDF$SE
  retDF$Pval <- pt(abs(retDF$Tvalue), df = df, lower.tail = FALSE) +
    pt(-abs(retDF$Tvalue), df = df, lower.tail = TRUE)

  return(retDF)
}


#' Extract coefficients and variance-covariance matrix from a tsDyn estimator
#'
#' Builds a coefficient vector and variance-covariance matrix for a fitted
#' `tsDyn` model object, along with associated coefficient names and residual
#' degrees of freedom.
#'
#' @param est A fitted `tsDyn` model object containing residuals, model matrix,
#'   coefficient matrices, and summary information.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{coef}{Named numeric vector of extracted coefficients.}
#'   \item{df}{Residual degrees of freedom used for inference.}
#'   \item{vcov}{Variance-covariance matrix of the extracted coefficients.}
#'   \item{NAMES}{Data frame linking row indices, coefficient names, and
#'   coefficient values.}
#' }
#'
#' @details
#' The function computes the residual covariance matrix from the model
#' residuals, applies a degrees-of-freedom adjustment, constructs the regressor
#' matrix used in estimation, and then forms the variance-covariance matrix
#' using a Kronecker product.
#'
#' Coefficient names are taken from `summary(est)$model`, cleaned, and expanded
#' into equation-specific names with prefixes `"eq1_"` and `"eq2_"`.
#'
#' This function assumes a two-equation `tsDyn` model structure with coefficient
#' blocks stored in `est$coefficients$Bdown`, `est$coefficients$Bmiddle`, and
#' `est$coefficients$Bup`.
#' @export
extract_tsdyn_vcov <- function(est) {
  Sigmabest <- matrix(
    1 / est$t * crossprod(est$residuals),
    ncol = est$k
  )

  SigmabestOls <- Sigmabest * (est$t / (est$t - ncol(est$coeffmat)))

  Z <- t(as.matrix(tail.matrix(est$model[, -c(1:est$k)], est$t)))
  vcov <- solve(tcrossprod(Z)) %x% SigmabestOls
  coef <- est$coefficients
  df <- (ncol(Z) - nrow(Z))

  NAMES <- names(summary(est)$model)
  NAMES <- NAMES[3:length(NAMES)]
  NAMES <- gsub("t_-", "t", gsub(" ", "_", NAMES))

  Idex <- 1:nrow(vcov)

  NAMES <- rbind(
    data.frame(
      indx = Idex[!Idex %in% seq(2, nrow(vcov), 2)],
      name = paste0("eq1_", NAMES),
      coef = c(coef$Bdown[1, ], coef$Bmiddle[1, ], coef$Bup[1, ])
    ),
    data.frame(
      indx = Idex[Idex %in% seq(2, nrow(vcov), 2)],
      name = paste0("eq2_", NAMES),
      coef = c(coef$Bdown[2, ], coef$Bmiddle[2, ], coef$Bup[2, ])
    )
  )

  NAMES <- NAMES[order(NAMES$indx), ]

  coef <- NAMES$coef

  names(coef) <- colnames(vcov) <- rownames(vcov) <- NAMES$name

  return(list(
    coef = coef,
    df = df,
    vcov = vcov,
    NAMES = NAMES
  ))
}
