#' @name PACKAGE_GLOBALVARIABLES
#' @title global_names
#' @description A combined dataset for global_names
#' @format A list of global names.
#' @source Internal innovation
PACKAGE_GLOBALVARIABLES <-  strsplit(
  "NAMES as.formula capture.output coef complete.cases cor cov df logLik
    pt series tail.matrix vcov",
  "\\s+"
)[[1]]

