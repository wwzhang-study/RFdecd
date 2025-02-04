#' @title Select the top 1000 features with the largest coefficient of variation in the estimated cell-type profiles (named CV)
#'
#' @param Y.raw A raw data matrix of complex samples
#' @param Prop0 estimated cell type proportion matrix
#' @param nMarker number of cell type specific markers
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
