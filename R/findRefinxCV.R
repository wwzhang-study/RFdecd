
#' @title Find index for cell-type specific features based on the largest coefficient
## of variation in the raw data
#'
#' @param rawdata A raw data matrix from complex samples.
#' @param nMarker Number of cell-type specific features.
#'
#' @return The index of top features with the largest coefficient of variation.
#' @import matrixStats
#' @export
#'

findRefinx.CV <- function(rawdata, nMarker=1000){
  mm = rowMeans(rawdata)
  vv = rowVars(rawdata)
  cv = sqrt(vv) / mm
  ##cv = vv
  cv[is.na(cv)] = 0
  ix = sort(cv, dec=TRUE, index=TRUE)$ix
  ix[1:nMarker]
}
