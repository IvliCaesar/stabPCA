\name{bbp_population}
\alias{bbp_population}
\title{Population spectrum estimated from a sample spectrum}
\description{
Replaces the eigenvalues above the Marchenko-Pastur edge by population estimates. It inverts the map of Baik,
Ben Arous and Peche over a flat bulk whose level is fixed by \eqn{\mathrm{tr}\,C = p}{tr C = p}.
}
\usage{
bbp_population(lambda, n, p)
}
\arguments{
\item{lambda}{sample eigenvalues of a correlation matrix, in decreasing order.}
\item{n}{number of observations.}
\item{p}{number of variables.}
}
\value{A numeric vector of the same length as \code{lambda}.}
\references{
Baik, J., Ben Arous, G. and Peche, S. (2005). Phase transition of the largest eigenvalue for nonnull complex sample
covariance matrices. \emph{Ann. Probab.} 33, 1643--1697.
}
\examples{
bbp_population(c(20, 8, 2, 1.5, 1, 0.5), n = 100, p = 6)
}
