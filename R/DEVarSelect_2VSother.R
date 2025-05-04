#' @title Dual-vs-Composite Cell-Type Marker Selection (DvC)
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

DEVarSelect_2VSother <- function(Y.raw, Prop0, nMarker = 1000){

  K = dim(Prop0)[2]
  N_sample = dim(Prop0)[1]

  ## find tissue specific genes
  idx = NULL
  ncount <- 1
  for(k in 1:K) {
    for(q in 2:K) {
      if(k < q) {
        cvec = rep(-2/(K-2),K)
        cvec[k] = 1
        cvec[q] = 1
        design = rep(0,N_sample)
        tmp = DEKTissue(K, Y=Y.raw, Prop=Prop0, design=design, contrast.vec=cvec)
        idx[[ncount]] = sort(abs(tmp$t.stat), decreasing =TRUE, index=TRUE)$ix
        ncount <- ncount + 1
      }
    }
  }

  tmpii <- K*(K-1)/4
  nMarker.tissue = nMarker/tmpii * 1.2 ## number of markers per tissue. Consider overlaps
  idxMarker = NULL
  for(ii in 1:tmpii) {
    idxMarker = c(idxMarker, idx[[ii]][1:nMarker.tissue])
  }
  idxMarker = unique(idxMarker)

  return(idxMarker)
}
