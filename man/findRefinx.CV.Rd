% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findRefinxCV.R
\name{findRefinx.CV}
\alias{findRefinx.CV}
\title{Find index for cell-type specific features based on the largest coefficient}
\usage{
findRefinx.CV(rawdata, nMarker = 1000)
}
\arguments{
\item{rawdata}{A raw data matrix from complex samples.}

\item{nMarker}{Number of cell-type specific features.}
}
\value{
The index of top features with the largest coefficient of variation.
}
\description{
Find index for cell-type specific features based on the largest coefficient
}
