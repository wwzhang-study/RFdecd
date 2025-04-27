#' @title Feature-selection-based reference-free deconvolution method
#'
#' @param DataType either 'Gene expression' or 'DNA methylation'
#' @param Y_raw A raw data matrix of complex samples
#' @param K Number of cell type
#' @param CTSoption feature selection options
#' @param nMarker number of cell type specific markers
#' @param InitMarker initial cell type specific markers
#' @param TotalIter number of iterations
#'
#' @return A list containing rmse,proportions with 30 iterations,
#' the estimated proportion matrix corresponding to the iteration with the smallest RMSE
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
