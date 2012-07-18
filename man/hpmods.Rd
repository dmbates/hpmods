\name{hpmods}
\alias{hpmods}
\title{Hierarchy-preserving models}
\usage{
  hpmods(formula, data, subset, weights, na.action, offset,
    ...)
}
\arguments{
  \item{formula}{as in \code{"\link{lm}"}}

  \item{data}{as in \code{"\link{lm}"}}

  \item{subset}{as in \code{"\link{lm}"}}

  \item{weights}{as in \code{"\link{lm}"}}

  \item{na.action}{as in \code{"\link{lm}"}}

  \item{offset}{as in \code{"\link{lm}"}}

  \item{\dots}{optional, additional arguments.  At present,
  none are used.}
}
\value{
  an object of class \code{"hpmods"}
}
\description{
  Determine hierarchy-preserving models
}
\details{
  Determine the hierarchy-preserving models from a formula
}
\examples{
set.seed(12321)
fr <- data.frame(y = rnorm(900), x1 = rnorm(900), x2 = rnorm(900),
                 x3 = gl(3,10,900), x4 = gl(10,3,900))
(hpm <- hpmods(y ~ (x1 + x2 + x3 + x4)^4, fr))
}
\author{
  Douglas Bates
}
\keyword{models}

