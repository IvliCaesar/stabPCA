\name{stabPCA}
\alias{stabPCA}
\alias{print.stabPCA}
\alias{plot.stabPCA}
\title{Reproducibility of principal components}
\description{
Predicts, from the eigenvalue spectrum alone, the overlap between the principal components computed on two random
halves of the data. Components that mix are grouped into subspaces, and the number of observations needed for a
target reproducibility is given.
}
\usage{
stabPCA(data, k_max = NULL, tau = 0.1, eps = 0.1, center = TRUE, scale. = TRUE,
        debias = c("none", "bbp"), sig = NULL)
\method{print}{stabPCA}(x, ...)
\method{plot}{stabPCA}(x, observed = NULL, ...)
}
\arguments{
\item{data}{numeric matrix or data frame, observations in rows and variables in columns, as in \code{sigPCA}.}
\item{k_max}{number of leading components to analyse; by default those above the Marchenko-Pastur edge.}
\item{tau}{threshold on the mutual term of two consecutive components above which they are grouped.}
\item{eps}{target: a split-half overlap of at least \code{1 - eps}.}
\item{center, scale.}{passed to \code{scale}; the defaults give the correlation matrix, as in \code{sigPCA}.}
\item{debias}{\code{"none"} uses the sample eigenvalues and is conservative; \code{"bbp"} uses
\code{\link{bbp_population}}.}
\item{sig}{optional result of \code{sigPCA::sigPCA} or \code{sigPCA::sigPCA_mp} on the same data; its significant
components are used.}
\item{x}{a \code{stabPCA} object.}
\item{observed}{optional numeric vector of observed overlaps, e.g. from \code{\link{splithalf_overlaps}}.}
\item{...}{passed to \code{plot}; ignored by \code{print}.}
}
\details{
With \eqn{m} observations, the \eqn{k}-th sample eigenvector fluctuates with variance
\deqn{A_k(m) = \frac{1}{m}\sum_{l \neq k}\frac{\lambda_k\lambda_l}{(\lambda_k - \lambda_l)^2}}{A_k(m) = (1/m) sum_{l != k} lambda_k lambda_l / (lambda_k - lambda_l)^2}
(Anderson, 1963). The predicted overlap between eigenvectors computed on two halves is \eqn{1 - A_k(n/2)}, and its
square \eqn{1 - 2A_k(n/2)}. For a group \eqn{G} of mixing components the predicted subspace overlap
\eqn{\mathrm{tr}(P_1P_2)/|G|}{tr(P1 P2)/|G|} keeps only the gaps to components outside \eqn{G}. An overlap of at least
\eqn{1 - \epsilon}{1 - eps} needs \eqn{n \ge 2a_k/\epsilon}{n >= 2 a_k / eps} with \eqn{a_k = mA_k(m)}{a_k = m A_k(m)}.
The approximation is first order: when \eqn{A_k} is of order one the prediction is negative and the component
is not reproducible.
}
\value{
An object of class \code{stabPCA}, a list with
\item{table}{one row per component: eigenvalue, significance, \code{A_half}, predicted overlap and its square,
group, observations required, and whether the component is reproducible on its own.}
\item{subspaces}{one row per group of components, with its predicted subspace overlap and the observations required.}
\item{eigenvalues, vectors}{the spectrum and the leading eigenvectors.}
\item{mp_bounds, significant_components}{as in \code{sigPCA_mp}.}
\item{n, p, tau, eps, debias}{the inputs.}
}
\references{
Anderson, T. W. (1963). Asymptotic theory for principal component analysis. \emph{Ann. Math. Statist.} 34, 122--148.
}
\seealso{\code{\link{splithalf_overlaps}}, \code{\link{bbp_population}}}
\examples{
set.seed(1)
n <- 600; p <- 200
ell <- c(30, 15, 9, 8.5, rep(1, p - 4))
Q <- qr.Q(qr(matrix(rnorm(p * p), p)))
X <- matrix(rnorm(n * p), n) \%*\% (t(Q) * sqrt(ell))
res <- stabPCA(X)
res
res$table
}
