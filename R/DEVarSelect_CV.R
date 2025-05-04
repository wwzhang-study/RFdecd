#' @title Coefficient of variation-based Cell-Type Marker Selection (SvC)
#'
#' @param Y.raw Raw data matrix (features × samples) with features as rows and
#'              samples as columns.
#' @param Prop0 Matrix (samples × K) of estimated cell type proportions from
#'             `RFdeconv` or `computeProp`. Each row should sum to ~1.
#' @param nMarker number of cell type specific markers. Default: 1000.
#'
#' @return estimated cell type specific markers
#' @export
#'
#'
DEVarSelect_CV <- function(Y.raw, Prop0, nMarker = 1000){

  K = dim(Prop0)[2]
  N_sample = dim(Prop0)[1]

  out_all = csSAM::csfit(Prop0, t(Y.raw))
  prof <- t(out_all$ghat)
  cvidx <- findRefinx.CV(prof, nMarker = nMarker)

  return(cvidx)
}
