

#' @title Perform Linear Regression-Based Deconvolution for Cell-Type Estimation
#'
#' @param cc A numeric matrix or vector containing the independent variable(s)
#'           (e.g., cell-type-specific reference profiles).
#' @param G A numeric matrix of observed mixed signals (features × samples).
#' @param logRm Logical indicating whether to apply inverse logarithmic transformation
#'              to the input data \code{G} (default: \code{FALSE}).
#' @param logBase Base of the logarithm used for transformation (default: 2).
#'
#'
#' @return A list containing two components:
#' \itemize{
#'   \item \code{ghat}: Estimated coefficients (cell-type proportions or reference profiles).
#'         If \code{logRm = TRUE}, coefficients are log-transformed.
#'   \item \code{residuals}: Residuals from the linear regression fit.
#' }
#'
#' @importFrom stats lsfit
#' @export
#'

csfit1 <- function (cc, G, logRm = FALSE, logBase = 2)
{
  if (logRm == TRUE) {
    G = logBase^G
  }
  fit1 = lsfit(cc, G, intercept = FALSE)
  #se1 = ls.diag(fit1)$std.err
  if (logRm == TRUE) {
    ghat = log(fit1$coefficients, logBase)
    ghat[is.nan(ghat)] = 0
    #se = log(se1, logBase)
    return(list(ghat = ghat, residuals = fit1$residuals))
  }
  else {
    return(list(ghat = fit1$coefficients, residuals = fit1$residuals))
  }
}
