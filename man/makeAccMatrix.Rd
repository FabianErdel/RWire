% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeAccMatrix.R
\name{makeAccMatrix}
\alias{makeAccMatrix}
\title{This function makes accessibility matrices from a set of BED files}
\usage{
makeAccMatrix(path, rois, nmax = 0)
}
\arguments{
\item{path}{Path to BED files}

\item{rois}{GRanges object containing ROIs}

\item{nmax}{Number of cells used for the analysis. When set to 0, all cells will be considered.}
}
\value{
Accessibility matrix
}
\description{
This function makes accessibility matrices from a set of BED files
}
