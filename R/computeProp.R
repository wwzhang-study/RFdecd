
#' @title Perform reference-free deconvolution to estimate cell-type proportion
#'
#' @param DataType Character specifying data type. Must be either:
#'   - "Gene expression" (default): Uses deconfounding via `deconf` package
#'   - "DNA methylation": Uses Reference-Free EWAS via `RefFreeEWAS` package
#' @param Y.raw A raw data matrix of complex samples where rows represent features (genes/CpG sites)
#'  and columns represent samples.
#' @param refinx estimated cell-type specific features
#' @param K the number of cell type
#'
#' @return estimated cell-type proportion matrix (K x samples), where rows
#' correspond to cell types and columns to samples.
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
