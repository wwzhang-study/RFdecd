
#' @title Perform reference-free deconvolution to estimate cell-type proportion
#'
#' @param DataType either 'Gene expression' or 'DNA methylation'
#' @param Y.raw A raw data matrix of complex samples
#' @param refinx estimated cell-type specific features
#' @param K the number of cell type
#'
#' @return estimated cell-type proportion matrix
#' @importFrom deconf deconfounding
#' @importFrom RefFreeEWAS RefFreeCellMix
#' @importFrom RefFreeEWAS RefFreeCellMixInitialize
#' @export
#'

computeProp <- function(DataType = "Gene expression", Y.raw, refinx, K){
  if (is.null(DataType)) {
    stop("Please specify 'DataType', should be either 'Gene expression' or 'DNA methylation'.")
  }
  Y = as.matrix(Y.raw[refinx,])

  if(DataType == DataType){
    outY = deconfounding(Y,K)
    Prop = t(outY$C$Matrix)
  }else{
    outY = RefFreeCellMix(Y, mu0 = RefFreeCellMixInitialize(Y,K = K))
    Prop = outY$Omega
  }
  return(Prop)
}
