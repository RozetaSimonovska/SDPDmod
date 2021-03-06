% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isrownor.R
\name{isrownor}
\alias{isrownor}
\title{Is the matrix row normalised}
\usage{
isrownor(W)
}
\arguments{
\item{W}{spatial weights matrix}
}
\value{
Logical value. If the weights matrix is row normalised
such that all rows sum up to 1, the value is TRUE.
}
\description{
Check if a spatial weights matrix is row normalised.
}
\seealso{
\code{\link{rownor}}
}
\author{
Rozeta Simonovska
}
