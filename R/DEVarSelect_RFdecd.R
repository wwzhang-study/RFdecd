#' @title a hybrid approach integrating SvC and DvC Cell-Type Marker Selection (RFdecd)
#'
#' @param Y.raw Raw data matrix (features × samples) with features as rows and
#'              samples as columns.
#' @param Prop0 Matrix (samples × K) of estimated cell type proportions from
#'             `RFdeconv` or `computeProp`. Each row should sum to ~1.
#' @param nMarker number of cell type specific markers. Default: 1000.
#'
#' @return estimated cell type specific markers
#' @export

DEVarSelect_RFdecd <- function(Y.raw,Prop0,nMarker = 1000){
  K = dim(Prop0)[2]
  N_sample = dim(Prop0)[1]

  ## 1-vs-others
  p_1 = matrix(NA,nrow = nrow(Y.raw),ncol = K)
  for(k in 1:K){
    cvec = rep(-1/(K-1),K)
    cvec[k] = 1
    design = rep(0,N_sample)
    tmp = DEKTissue(K, Y=Y.raw, Prop=Prop0, design=design, contrast.vec=cvec)
    muEst = tmp[,4:(3+K)]
    p_1[,k] = tmp$t.pval
  }

  ## 2-vs-others
  p_2 = matrix(NA,nrow = nrow(Y.raw),ncol = K*(K-1)/2)
  ncount <- 1
  for(k in 1:K) {
    for(q in 2:K) {
      if(k < q) {
        cvec = rep(-2/(K-2),K)
        cvec[k] = 1
        cvec[q] = 1
        design = rep(0,N_sample)
        tmp = DEKTissue(K, Y=Y.raw, Prop=Prop0, design=design, contrast.vec=cvec)
        muEst = tmp[,4:(3+K)]
        p_2[,ncount] = tmp$t.pval
        ncount <- ncount + 1
      }
    }
  }
  rownames(p_1) = rownames(p_2) = rownames(Y.raw)

  ## select CTS features
  nMarker.tissue_1 = 100
  markers <- NULL
  for(k in 1:K){
    p_1 = p_1[setdiff(rownames(p_1),markers),]
    markers = c(markers,names(p_1[,k])[sort(p_1[,k],decreasing = F,index = TRUE)$ix][1:nMarker.tissue_1])
  }

  nMarker.tissue_2 = 100
  for(k in 1:(K*(K-1)/4)){
    p_2 = p_2[setdiff(rownames(p_2),markers),]
    markers = c(markers,names(p_2[,k])[sort(p_2[,k],decreasing = F,index = TRUE)$ix][1:nMarker.tissue_2])
  }

  p_1 = p_1[setdiff(rownames(p_1),markers),]
  p_2 = p_2[setdiff(rownames(p_2),markers),]

  p = c(p_1[,1],p_1[,2],p_1[,3],p_1[,4],p_2[,1],p_2[,2],p_2[,3])
  markers = c(markers,names(p)[sort(p,decreasing = FALSE,index = TRUE)$ix][1:300])

  markers = unique(markers)
  idxMarker = match(markers,rownames(Y.raw))
  return(idxMarker)
}
