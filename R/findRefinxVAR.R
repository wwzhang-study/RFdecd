#' @title find index for cell-type specific features based on the largest variation in the raw data
#' @param rawdata A raw data matrix from complex samples.
#' @param nMarker Number of cell-type specific features.
#'
#' @return The index of top features with the largest variation.
#' @import matrixStats
#' @export

findRefinx.VAR <- function(rawdata, nMarker=1000) {
  vv = rowVars(log(rawdata+1))
  vv[is.na(vv)] = 0
  ix = sort(vv, dec=TRUE, index=TRUE)$ix
  ix[1:nMarker]
}
