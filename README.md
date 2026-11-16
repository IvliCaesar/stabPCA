# stabPCA

**Which principal components are reproducible?** An eigenvalue above the Marchenko–Pastur edge is significant, but
its eigenvector need not survive a change of sample. `stabPCA` predicts, from the eigenvalue spectrum alone, the
overlap between the components computed on two random halves of the data. It groups components that mix into
subspaces and gives the number of observations needed for a target reproducibility.

It complements [sigPCA](https://github.com/guillermodeandajauregui/sigPCA) (CRAN), which tests which eigenvalues
are significant, and takes the same input.

**Significant is not reproducible.** In three large breast-cancer cohorts, 45 to 69 components pass the
Marchenko–Pastur test, but only 3 to 7 are recovered with overlap ≥ 0.9 from two halves of the same cohort.
`stabPCA` predicts which ones without resampling.

| | significance of eigenvalues | per-component reproducibility | mixing subspaces | sample size | no resampling |
|---|---|---|---|---|---|
| sigPCA, PCAtools (`chooseMarchenkoPastur`) | yes | no | no | no | yes |
| bootstrap of loadings (bootSVD, multivarious) | no | yes | no | no | no |
| split-half matching (SHEM) | no | yes | no | no | no |
| **stabPCA** | yes | yes | yes | yes | yes |

## Installation

```r
# install.packages("remotes")
remotes::install_github("IvliCaesar/stabPCA")
```

Base R only (`stats`, `graphics`).

## Use

```r
library(stabPCA)
res <- stabPCA(X)            # rows = samples, columns = variables (genes), as in sigPCA
res                          # reproducible components, and components to report as subspaces
res$table                    # per component: A_k(n/2), predicted overlap, group, observations required
res$subspaces                # groups of mixing components and their predicted subspace overlap

# with sigPCA
sig <- sigPCA::sigPCA(X, method = "mp")
stabPCA(X, sig = sig)

# population correction of the spikes, and an empirical check
stabPCA(X, debias = "bbp")
obs <- splithalf_overlaps(X, k_max = 10, n_splits = 20)
plot(res, observed = obs$individual)
```

## What it computes

With $m$ observations, the $k$-th sample eigenvector fluctuates with variance (Anderson, 1963)

$$A_k(m) = \frac1m\sum_{l\neq k}\frac{\lambda_k\lambda_l}{(\lambda_k-\lambda_l)^2}.$$

- **Single components.** The expected overlap between eigenvectors computed on two halves is $1 - A_k(n/2)$.
- **Subspaces.** Two close eigenvalues mix their eigenvectors, but the plane they span is controlled only by the
  gaps to the rest of the spectrum. Components whose mutual term exceeds `tau` are reported together, with the
  predicted overlap $\mathrm{tr}(P_1P_2)/|G|$.
- **Sample size.** An overlap of at least $1-\varepsilon$ needs $n \ge 2a_k/\varepsilon$, where $a_k = mA_k(m)$. The
  plug-in estimate is conservative.
- **Bulk.** In the bulk the gaps repel linearly, the fluctuation has a tail of index 1 and no mean, and
  eigenvectors are not reproducible.

## Validation

Four breast-cancer cohorts (2,000 most variable genes, 20 random splits into halves):

| Cohort | Tumours | Platform | Mean error of the predicted split-half overlap, isolated components among the first 14 |
|---|---|---|---|
| SCAN-B (GSE96058) | 3,273 | RNA-seq | 0.009 |
| TCGA-BRCA | 1,097 | RNA-seq | 0.024 (0.019 with `debias = "bbp"`) |
| METABRIC | 1,980 | Illumina HT-12 | 0.026 (0.023) |
| GSE87049 (Mexico) | 124 | Affymetrix Gene 1.0 ST | one reproducible component at p = 2,000 |

In SCAN-B, components 10 and 11 have a mean squared overlap of 0.46 each, and their plane 0.85, against a
predicted 0.84.

## Citation

C. Galindo and E. Hernández-Lemus, *Which eigenvectors of a co-expression matrix are reproducible?*, manuscript, 2026.

## License

MIT.
