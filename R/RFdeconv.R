#' @title Feature-selection-based reference-free deconvolution method
#'
#' @description
#' Estimates cell-type proportions from bulk omics data through iterative refinement of
#' cell-type-specific (CTS) feature selection. Supports both gene expression and DNA
#' methylation data with automatic marker optimization.
#'
#' @param DataType either 'Gene expression' or 'DNA methylation'
#' @param Y_raw Numeric matrix. Raw input data (features × samples) where rows
#'   represent genes/probes and columns represent biological samples.
#' @param K Integer (≥2). Number of cell types to deconvolve.
#' @param CTSoption feature selection options.
#' DEVarSelect_CV – Coefficient of variation; DEVarSelect_VAR – Variance;
#' DEVarSelect_1VSother – Single vs. Composite (SvC); DEVarSelect_2VSother – Dual vs. Composite (DvC);
#' DEVarSelect_pairwise – Pairwise Direct (PwD); DEVarSelect_RFdecd – RFdecd.
#' Default: `DEVarSelect_1VSother`
#' @param nMarker number of cell type specific markers. Default: 1000.
#' @param InitMarker Numeric/Logical vector. Initial feature indices for first
#'   iteration. If `NULL` (default), automatically identifies markers via
#'   cross-validation (`findRefinx.CV`).
#' @param TotalIter number of iterations. Default: 30.
#'
#' @return List containing:
#'   - `allRMSE`: Numeric vector of RMSE values across all iterations
#'   - `allProp`: List of proportion matrices from each iteration
#'   - `estProp`: Optimal proportion matrix (minimum RMSE iteration)
#' @export
#'

RFdeconv <- function(DataType = "Gene expression",
                             Y_raw,
                             K,
                             CTSoption = DEVarSelect_1VSother,
                             nMarker = 1000,
                             InitMarker = NULL,
                             TotalIter = 30) {

  if (is.null(InitMarker)) {
    InitMarker <- findRefinx.CV(Y_raw, nMarker = nMarker)
  }

  Prop0 = computeProp(DataType = DataType, Y_raw, InitMarker, K)

  allProp <- list()
  allRMSE <- rep(0, TotalIter + 1)

  out_all <- csSAM::csfit(Prop0, t(Y_raw))
  prof <- t(out_all$ghat)
  prof[prof < 0] = 0
  tmpmat <- prof %*% t(Prop0)

  allProp[[1]] <- Prop0
  allRMSE[1] <- sqrt(mean((t(Y_raw) - t(tmpmat)) ^ 2))

  message("+========================================+")
  message("+======= Total iterations = ",
          TotalIter, " ==========+")

  for (i in seq_len(TotalIter)) {
    message("Current iter = ", i)

    updatedInx <- CTSoption(Y_raw, Prop0, nMarker)
    Prop0 = computeProp(DataType = DataType, Y_raw, updatedInx, K)
    allProp[[i + 1]] <- Prop0

    out_all <- csfit1(Prop0, t(Y_raw))
    prof <- t(out_all$ghat)
    prof[prof < 0] = 0
    tmpmat <- prof %*% t(Prop0)
    allRMSE[i + 1] <- sqrt(mean((t(Y_raw) - t(tmpmat)) ^ 2))
  }

  min_idx <- which.min(allRMSE)
  Prop0 <- allProp[[min_idx]]

  return(list(allRMSE = allRMSE,
              allProp = allProp,
              estProp = Prop0
  ))
}
