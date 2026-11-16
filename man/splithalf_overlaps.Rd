\name{splithalf_overlaps}
\alias{splithalf_overlaps}
\title{Empirical split-half overlaps of principal components}
\description{
Splits the observations into two random halves repeatedly and returns the mean absolute overlap between the
corresponding eigenvectors, and optionally the overlap of subspaces.
}
\usage{
splithalf_overlaps(data, k_max = 10, n_splits = 20, groups = NULL, center = TRUE, scale. = TRUE)
}
\arguments{
\item{data}{numeric matrix or data frame, observations in rows.}
\item{k_max}{number of leading components.}
\item{n_splits}{number of random splits into halves.}
\item{groups}{optional list of integer vectors; for each group \eqn{G} the subspace overlap
\eqn{\mathrm{tr}(P_1P_2)/|G|}{tr(P1 P2)/|G|} is returned.}
\item{center, scale.}{passed to \code{scale}.}
}
\value{A list with \code{individual}, the mean absolute overlap per component, and \code{subspace}.}
\examples{
set.seed(2)
X <- matrix(rnorm(400 * 50), 400) \%*\% diag(c(4, 3, rep(1, 48)))
splithalf_overlaps(X, k_max = 3, n_splits = 5)
}
