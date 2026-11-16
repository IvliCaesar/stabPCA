# stabPCA: reproducibility of principal components from the eigenvalue spectrum.
#
# For n observations, the fluctuation of the k-th sample eigenvector has variance (Anderson, 1963)
#   A_k(m) = (1/m) sum_{l != k} lambda_k lambda_l / (lambda_k - lambda_l)^2,   m = sample size.
# Between eigenvectors computed on two independent halves (m = n/2) the expected overlap is 1 - A_k(n/2).
# Close eigenvalues mix; the subspace G they span has expected tr(P1 P2)/|G| = 1 - (2/|G|) sum_{k in G} A_k^out,
# with A_k^out restricted to l outside G. An overlap of at least 1 - eps needs n >= 2 a_k / eps, a_k = m A_k(m).

.mp_edge <- function(p, n) (1 + sqrt(p / n))^2

.prepare <- function(data, center, scale.) {
  stopifnot(is.matrix(data) || is.data.frame(data))
  X <- as.matrix(data)
  keep <- apply(X, 2, function(x) all(is.finite(x)) && stats::var(x) > 0)
  if (!all(keep)) warning(sum(!keep), " column(s) with non-finite values or zero variance removed")
  X <- X[, keep, drop = FALSE]
  if (center || scale.) X <- scale(X, center = center, scale = scale.)
  X
}

.pair_var <- function(lambda, kmax, m) {
  V <- matrix(0, kmax, length(lambda))
  for (k in seq_len(kmax)) {
    v <- lambda[k] * lambda / (m * (lambda[k] - lambda)^2)
    v[k] <- 0
    V[k, ] <- v
  }
  V
}

.group <- function(V, kmax, tau) {
  groups <- list(); g <- 1L
  if (kmax > 1) for (k in 2:kmax) {
    if (V[k - 1, k] > tau) g <- c(g, k) else { groups[[length(groups) + 1]] <- g; g <- k }
  }
  groups[[length(groups) + 1]] <- g
  groups
}

#' Population spectrum estimated from a sample spectrum
#'
#' Replaces the eigenvalues above the Marchenko-Pastur edge by population estimates, inverting the
#' Baik-Ben Arous-Peche map over a flat bulk whose level is fixed by trace = p.
#'
#' @param lambda sample eigenvalues of a correlation matrix, in decreasing order.
#' @param n number of observations.
#' @param p number of variables.
#' @return numeric vector of the same length as \code{lambda}.
#' @export
bbp_population <- function(lambda, n, p) {
  g <- p / n
  s <- sum(lambda > .mp_edge(p, n))
  if (s == 0) return(rep(1, length(lambda)))
  b <- 1
  for (it in 1:50) {
    cc <- lambda[1:s] + b - g * b
    ell <- (cc + sqrt(pmax(cc^2 - 4 * lambda[1:s] * b, 0))) / 2
    b <- (p - sum(ell)) / (p - s)
  }
  c(ell, rep(b, length(lambda) - s))
}

#' Reproducibility of principal components
#'
#' Predicts, from the eigenvalue spectrum alone, the overlap between the principal components computed on two
#' random halves of the data, groups components that mix into subspaces, and gives the number of observations
#' needed for a target reproducibility.
#'
#' @param data numeric matrix or data frame, observations in rows and variables in columns (as in sigPCA).
#' @param k_max number of leading components to analyse; by default those above the Marchenko-Pastur edge.
#' @param tau threshold on the mutual term of two consecutive components above which they are grouped.
#' @param eps target: a split-half overlap of at least \code{1 - eps}.
#' @param center,scale. passed to \code{scale}; the defaults give the correlation matrix, as in sigPCA.
#' @param debias \code{"none"} uses the sample eigenvalues (conservative); \code{"bbp"} uses \code{bbp_population}.
#' @param sig optional result of \code{sigPCA::sigPCA} or \code{sigPCA::sigPCA_mp} on the same data; its
#'   significant components are reported alongside.
#' @return an object of class \code{stabPCA}: a list with \code{table} (one row per component),
#'   \code{subspaces}, \code{eigenvalues}, \code{vectors}, \code{mp_bounds}, \code{significant_components},
#'   \code{n}, \code{p}.
#' @references Anderson, T. W. (1963). Asymptotic theory for principal component analysis.
#'   Ann. Math. Statist. 34, 122-148.
#' @examples
#' set.seed(1)
#' n <- 600; p <- 200
#' ell <- c(30, 15, 9, 8.5, rep(1, p - 4))
#' Q <- qr.Q(qr(matrix(rnorm(p * p), p)))
#' X <- matrix(rnorm(n * p), n) %*% (t(Q) * sqrt(ell))
#' res <- stabPCA(X)
#' res
#' res$table
#' @export
stabPCA <- function(data, k_max = NULL, tau = 0.1, eps = 0.1, center = TRUE, scale. = TRUE,
                    debias = c("none", "bbp"), sig = NULL) {
  debias <- match.arg(debias)
  X <- .prepare(data, center, scale.)
  n <- nrow(X); p <- ncol(X)
  e <- eigen(stats::cov(X), symmetric = TRUE)
  lambda <- e$values
  edge <- .mp_edge(p, n)
  sigc <- which(lambda > edge)
  if (!is.null(sig)) {
    s <- if (!is.null(sig$mp)) sig$mp else sig
    if (!is.null(s$significant_components)) sigc <- s$significant_components
  }
  if (is.null(k_max)) k_max <- max(length(sigc), 1)
  k_max <- max(1, min(k_max, length(lambda) - 1))
  m <- n / 2
  lam_used <- if (debias == "bbp") bbp_population(lambda, n, p) else lambda
  V <- .pair_var(lam_used, k_max, m)
  A <- rowSums(V)
  groups <- .group(V, k_max, tau)
  gid <- integer(k_max)
  for (i in seq_along(groups)) gid[groups[[i]]] <- i
  sub <- do.call(rbind, lapply(seq_along(groups), function(i) {
    G <- groups[[i]]
    Aout <- vapply(G, function(k) sum(V[k, -G]), numeric(1))
    data.frame(group = i, components = paste(G, collapse = "-"), size = length(G),
               predicted_overlap_sq = 1 - 2 * sum(Aout) / length(G),
               n_required = ceiling(4 * m * sum(Aout) / length(G) / (2 * eps - eps^2)))
  }))
  tab <- data.frame(k = seq_len(k_max), eigenvalue = lambda[seq_len(k_max)],
                    significant = seq_len(k_max) %in% sigc,
                    A_half = A, predicted_overlap = 1 - A, predicted_overlap_sq = 1 - 2 * A,
                    group = gid, n_required = ceiling(2 * m * A / eps))
  tab$reproducible <- tab$significant & tab$predicted_overlap >= 1 - eps & sub$size[tab$group] == 1
  structure(list(table = tab, subspaces = sub, eigenvalues = lambda,
                 vectors = e$vectors[, seq_len(k_max), drop = FALSE],
                 mp_bounds = list(lambda_min = (1 - sqrt(p / n))^2, lambda_max = edge),
                 significant_components = sigc, n = n, p = p, tau = tau, eps = eps, debias = debias),
            class = "stabPCA")
}

#' @export
print.stabPCA <- function(x, ...) {
  cat(sprintf("stabPCA: n = %d, p = %d, Marchenko-Pastur edge %.3f, %d significant components\n",
              x$n, x$p, x$mp_bounds$lambda_max, length(x$significant_components)))
  rep_k <- x$table$k[x$table$reproducible]
  cat(sprintf("Reproducible components (split-half overlap >= %.2f, isolated): %s\n", 1 - x$eps,
              if (length(rep_k)) paste(rep_k, collapse = ", ") else "none"))
  multi <- x$subspaces[x$subspaces$size > 1, ]
  if (nrow(multi)) {
    cat("Components that mix; report them as subspaces:\n")
    print(multi, row.names = FALSE, digits = 3)
  }
  invisible(x)
}

#' Empirical split-half overlaps of principal components
#'
#' @param data numeric matrix or data frame, observations in rows.
#' @param k_max number of leading components.
#' @param n_splits number of random splits into halves.
#' @param groups optional list of integer vectors; for each, the subspace overlap tr(P1 P2)/|G| is returned.
#' @param center,scale. passed to \code{scale}.
#' @return a list with \code{individual} (mean absolute overlap per component) and \code{subspace}.
#' @examples
#' set.seed(2)
#' X <- matrix(rnorm(400 * 50), 400) %*% diag(c(4, 3, rep(1, 48)))
#' splithalf_overlaps(X, k_max = 3, n_splits = 5)
#' @export
splithalf_overlaps <- function(data, k_max = 10, n_splits = 20, groups = NULL, center = TRUE, scale. = TRUE) {
  X <- .prepare(data, center, scale.)
  n <- nrow(X); h <- n %/% 2
  vec <- function(Y) eigen(stats::cov(Y), symmetric = TRUE)$vectors[, seq_len(k_max), drop = FALSE]
  ind <- matrix(NA_real_, n_splits, k_max)
  subs <- if (is.null(groups)) NULL else matrix(NA_real_, n_splits, length(groups))
  for (s in seq_len(n_splits)) {
    idx <- sample.int(n)
    V1 <- vec(X[idx[1:h], , drop = FALSE]); V2 <- vec(X[idx[(h + 1):(2 * h)], , drop = FALSE])
    ind[s, ] <- abs(colSums(V1 * V2))
    if (!is.null(groups)) for (i in seq_along(groups)) {
      G <- groups[[i]]
      subs[s, i] <- sum(crossprod(V1[, G, drop = FALSE], V2[, G, drop = FALSE])^2) / length(G)
    }
  }
  list(individual = colMeans(ind), subspace = if (is.null(subs)) NULL else colMeans(subs))
}

#' Plot predicted against observed split-half overlaps
#'
#' @param x a \code{stabPCA} object.
#' @param observed optional numeric vector of observed overlaps, e.g. from \code{splithalf_overlaps}.
#' @param ... passed to \code{plot}.
#' @export
plot.stabPCA <- function(x, observed = NULL, ...) {
  t <- x$table
  graphics::plot(t$k, pmax(t$predicted_overlap, 0), type = "b", pch = 1, lty = 2, col = "#D55E00",
                 ylim = c(0, 1), xlab = "component k", ylab = "overlap between halves", ...)
  if (!is.null(observed)) graphics::points(seq_along(observed), observed, pch = 19, col = "#0072B2")
  graphics::abline(h = 1 - x$eps, lty = 3, col = "grey40")
  graphics::legend("bottomleft", bty = "n", pch = c(1, 19), lty = c(2, NA), col = c("#D55E00", "#0072B2"),
                   legend = c("predicted 1 - A_k(n/2)", "observed"))
  invisible(x)
}
