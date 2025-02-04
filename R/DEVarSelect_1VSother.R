#' @title Cross-cell type differential analysis between one cell type and the other cell types (named as 1_O)
#'
#' @param Y.raw a raw data matrix of complex samples
#' @param Prop0 estimated cell-type proportion matrix
#' @param nMarker the number of cell type specific markers
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
    idx[[k]] = sort(abs(tmp$t.stat), dec=TRUE, index=TRUE)$ix
  }
  nMarker.tissue = nMarker/K * 1.2 ## number of markers per tissue. Consider overlaps
  idxMarker = NULL
  for(k in 1:K) {
    idxMarker = c(idxMarker, idx[[k]][1:nMarker.tissue])
  }
  idxMarker = unique(idxMarker)

  return(idxMarker)
}
