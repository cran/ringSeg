# ringSeg

Rank-based, asymptotic **distribution-free change-point detection** (the RING method,
Zhou & Chen 2025, IEEE TIT). Detects a single change-point or a changed interval in a
sequence of (possibly high-dimensional / non-Euclidean) observations, with **analytic
distribution-free p-values** (no permutation needed) plus an optional skewness
correction and optional permutation p-values.

The method ranks each observation's neighbours within a sparse k-nearest-neighbour
graph and reduces the scan to three rank statistics -- weighted (`WR`), max-type
(`MR`), and generalized (`TR`) -- whose null distribution is asymptotically
distribution-free, so p-values are available in closed form. Inputs may be a data
matrix, a distance matrix, or a precomputed rank matrix; the data need not be
Euclidean.

## Install
```r
# from the built tarball
install.packages("ringSeg_0.1.0.tar.gz", repos = NULL, type = "source")
# or from source dir:  R CMD INSTALL ringSeg
```

## Use
```r
library(ringSeg)

set.seed(1)
n <- 200; d <- 10; tau <- 100
X <- matrix(rnorm(n * d), n, d)
X[(tau + 1):n, ] <- X[(tau + 1):n, ] * 1.6        # a dispersion change at t = 100

res <- ring_cpd(X, skew.corr = TRUE)              # data in (rows = time order)
res$scanZ$TR$tau        # estimated change-point
res$pval.appr$TR        # generalized-statistic analytic p-value
res$pval.appr$WR.cor    # skewness-corrected weighted-statistic p-value
```

By default `ring_cpd()` builds RING's **row-rank k-nearest-neighbour graph**
(`k = round(n^0.65)`) -- the sparse "new ranking scheme" the method is built on.
It also accepts a distance matrix (`is.distance = TRUE`) or a precomputed rank matrix
(`is.rank = TRUE`). For full control, build the graph with `ring_graph()` and call
`rcpd()`:
```r
D   <- as.matrix(dist(X))
R   <- ring_graph(D, k = 13)                   # RING k-NN rank graph
res <- rcpd(R, skew.corr = TRUE, B = 0)        # B > 0 adds permutation p-values
```

## Functions
The package exposes four functions, layered from high-level (data in) to low-level:

| Function | Input | Does |
|---|---|---|
| `ring_cpd(x, ...)` | data matrix (rows = time order), or a distance / rank matrix | **Top-level wrapper.** Builds the RING graph from `x` (via `ring_graph`) and runs the scan (via `rcpd`). Start here. |
| `rcpd(R, ...)` | an `n x n` rank/weight matrix `R` | **Core engine.** Scans for a change-point / interval, returns the WR/MR/TR statistics, the estimated change-point(s), analytic p-values, and (if `B > 0`) permutation p-values. |
| `ring_graph(D, k)` | an `n x n` distance matrix `D` | Builds RING's sparse **row-rank k-NN rank graph** -- the "new ranking scheme" `rcpd` operates on. Returns `R`. |
| `Rise_Rank(S, method)` | an `n x n` similarity matrix `S` | Low-level utility: ranks the entries of `S` (`"row"` = per-row ranks, used by `ring_graph`; `"overall"` = all pairs jointly). Rarely called directly. |

Relationship: `ring_cpd(x)` is essentially `rcpd(ring_graph(as.matrix(dist(x))))`. Use
`ring_cpd` for "data in, result out"; drop to `ring_graph` + `rcpd` when you want full
control over the similarity / graph (e.g. non-Euclidean data, a custom kernel, or
reusing one `R`).

## Statistics
- **WR** -- weighted rank statistic (`Z_w`)
- **MR** -- max-type, `max(|Z_diff|, Z_w)`
- **TR** -- generalized omnibus, `Z_w^2 + Z_diff^2` (recommended default)

## Input handling (warn vs. error)
Inputs are checked with a consistent policy: coercible quirks are accepted with a
**warning**, genuinely broken inputs are a hard **error**.

| Input | Behaviour |
|---|---|
| `D` (distance) with **negative** entries | **warning**, proceeds -- `S = max(D) - D` still ranks, so a general (dis)similarity is allowed |
| `n0` / `n1` / `k` **non-integer** but in range | **warning**, rounded to the nearest integer |
| `R` **not symmetric** | **warning**, symmetrized as `(R + t(R)) / 2` |
| `n0` / `n1` outside `[0, n]`, or `k` outside `1..n-1` | **error** |
| `D`, `x`, or `R` with `NA` / `NaN` / `Inf` | **error** (not coercible) |
| `R` not square / negative weights / non-zero diagonal | **error** (invalid rank/weight matrix) |
| fewer than 4 observations | **error** (the scan statistic is undefined) |
| `B` not a single non-negative integer; `pval.appr` / `skew.corr` not a single `TRUE`/`FALSE` | **error** |

Diagnostic notes (e.g. the scan-range adjustment) are emitted via `message()`, so they
can be silenced with `suppressMessages()`.

## Notes on calibration
The skewness-corrected analytic p-value (`*.cor`, the default) is well-calibrated at
finite *n*; the uncorrected Gaussian/field approximation can be mildly
anti-conservative for small *n*. For small samples, `B > 0` (permutation) is always
available as an exact alternative. Permutation p-values use the add-one estimator
`(#{perm >= obs} + 1) / (B + 1)`, so they are always strictly positive.

## Maintainer / license
Maintainer Hao Chen <hxchen@ucdavis.edu>. License `GPL (>= 2)`, matching the sibling
packages `gSeg` and `kerSeg`. Citation (`citation("ringSeg")`): Zhou, D. & Chen, H.
(2025), *Asymptotic Distribution-Free Change-Point Detection for Modern Data Based on a
New Ranking Scheme*, IEEE Transactions on Information Theory, **71**(8), 6183-6197.
[doi:10.1109/TIT.2025.3575858](https://doi.org/10.1109/TIT.2025.3575858).
