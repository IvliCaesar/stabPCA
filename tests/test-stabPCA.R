library(stabPCA)

# spiked model with known spectrum: predictions within 0.05 of the empirical split-half overlaps
set.seed(1)
n <- 1000; p <- 300
ell <- c(40, 20, 12, 11.6, rep(1, p - 4))
Q <- qr.Q(qr(matrix(rnorm(p * p), p)))
X <- matrix(rnorm(n * p), n) %*% (t(Q) * sqrt(ell))

res <- stabPCA(X, k_max = 4)
stopifnot(inherits(res, "stabPCA"), nrow(res$table) == 4)
stopifnot(all(res$table$significant))
stopifnot(identical(res$subspaces$components[res$subspaces$size > 1], "3-4"))   # 12 and 11.6 mix

obs <- splithalf_overlaps(X, k_max = 4, n_splits = 10, groups = list(3:4))
stopifnot(all(abs(obs$individual[1:2] - res$table$predicted_overlap[1:2]) < 0.05))
plane <- res$subspaces$predicted_overlap_sq[res$subspaces$components == "3-4"]
stopifnot(abs(obs$subspace - plane) < 0.08)
stopifnot(obs$subspace > mean(obs$individual[3:4]^2))                         # the plane is more stable

# sample-size rule is monotone in eps
r1 <- stabPCA(X, k_max = 2, eps = 0.1); r2 <- stabPCA(X, k_max = 2, eps = 0.05)
stopifnot(all(r2$table$n_required >= r1$table$n_required))

# population correction moves spikes down and keeps the trace
pop <- bbp_population(res$eigenvalues, n, p)
stopifnot(abs(sum(pop) - p) < 1e-6, pop[1] < res$eigenvalues[1])

# sigPCA-compatible input: list with $mp$significant_components
fake_sig <- list(mp = list(significant_components = 1:3))
stopifnot(identical(stabPCA(X, sig = fake_sig)$significant_components, 1:3))

# zero-variance columns are dropped with a warning
Xz <- cbind(X[, 1:20], 1)
w <- tryCatch(stabPCA(Xz, k_max = 2), warning = function(w) "warned")
stopifnot(identical(w, "warned"))
cat("stabPCA tests passed\n")
