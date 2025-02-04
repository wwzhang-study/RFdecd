#' @title Cross-cell type differential analysis between two cell types and the other cell types (named as 2_O)
#'
#' @param Y.raw A raw data matrix of complex samples
#' @param Prop0 estimated cell type proportion matrix
#' @param nMarker the number of cell type specific markers
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
        idx[[ncount]] = sort(abs(tmp$t.stat), dec=TRUE, index=TRUE)$ix
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
