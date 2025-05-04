#' @title Single-vs-Composite Cell-Type Marker Selection (SvC)
#'
#' @param Y.raw Raw data matrix (features × samples) with features as rows and
#'              samples as columns.
#' @param Prop0 Matrix (samples × K) of estimated cell type proportions from
#'             `RFdeconv` or `computeProp`. Each row should sum to ~1.
#' @param nMarker number of cell type specific markers. Default: 1000.
#'
#' @return estimated cell-type specific markers
#' @export
#'

DEVarSelect_1VSother <- function(Y.raw, Prop0, nMarker = 1000){

  K = dim(Prop0)[2]
  N_sample = dim(Prop0)[1]

  ## find tissue specific genes
  idx = NULL
  for(k in 1:K) {
    cvec = rep(-1/(K-1),K)
    cvec[k] = 1
    design = rep(0,N_sample)
    tmp = DEKTissue(K, Y=Y.raw, Prop=Prop0, design=design, contrast.vec=cvec)
    idx[[k]] = sort(abs(tmp$t.stat), decreasing =TRUE, index=TRUE)$ix
  }
  nMarker.tissue = nMarker/K * 1.2 ## number of markers per tissue. Consider overlaps
  idxMarker = NULL
  for(k in 1:K) {
    idxMarker = c(idxMarker, idx[[k]][1:nMarker.tissue])
  }
  idxMarker = unique(idxMarker)

  return(idxMarker)
}
