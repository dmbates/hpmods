\name{msubs}
\alias{msubs}
\title{Extract model subsets}
\usage{
  msubs(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"hpmods"}}

  \item{\dots}{optional, additional arguments.  At present,
  none are used.}
}
\value{
  a list of model subsets
}
\description{
  Extract model subsets from an hpmods object
}
\details{
  After fitting a model with several terms the residual
  sums of squares (and the coefficient estimates, if
  needed) for any leading subset of those terms can be
  determined directly.  This function extracts distinct
  subsets of the rank-preserving models that can where each
  subset can be analyzed from a single model fit.
}
\author{
  Douglas Bates
}

