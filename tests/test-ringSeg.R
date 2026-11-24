## Base-R tests (run by R CMD check). No testthat dependency.
library(ringSeg)

## 1. Analytic critical-value self-test: the b-values from the paper give p ~ 0.05.
n <- 1000; n0 <- 25; n1 <- n - n0; l0 <- 50; l1 <- n - l0
stopifnot(abs(ringSeg:::pTR1(14.11, n, n0, n1) - 0.05) < 0.01)
stopifnot(abs(ringSeg:::pWR1(3.14,  n, n0, n1) - 0.05) < 0.01)
stopifnot(abs(ringSeg:::pMR1(3.38,  n, n0, n1) - 0.05) < 0.01)
stopifnot(abs(ringSeg:::pTR2(23.96, n, l0, l1) - 0.05) < 0.01)
stopifnot(abs(ringSeg:::pWR2(4.22,  n, l0, l1) - 0.05) < 0.01)
stopifnot(abs(ringSeg:::pMR2(4.34,  n, l0, l1) - 0.05) < 0.01)

## 2. ring_cpd detects a clear change and localizes it near the truth.
set.seed(1)
nn <- 200; d <- 10; tau <- 100
X <- matrix(rnorm(nn * d), nn, d)
X[(tau + 1):nn, ] <- X[(tau + 1):nn, ] * 1.8       # dispersion change
res <- suppressMessages(ring_cpd(X, skew.corr = TRUE, B = 0))
stopifnot(res$pval.appr$TR < 0.01)                 # detected
stopifnot(abs(res$scanZ$TR$tau - tau) <= 10)       # localized

## 3. ring_graph builds the RING k-NN rank graph (sparse) and rcpd on it matches
##    ring_cpd's own pipeline (is.rank path).
R <- ring_graph(as.matrix(dist(X)), k = round(nn^0.65))
stopifnot(dim(R)[1] == nn, dim(R)[2] == nn, all(diag(R) == 0))
stopifnot(mean(R > 0) < 0.5)                       # k-NN graph is sparse
res2 <- suppressMessages(rcpd(R, B = 0))
stopifnot(isTRUE(all.equal(res2$scanZ$TR$tau, res$scanZ$TR$tau)))   # same graph as ring_cpd default

## 4. Interval alternative runs and returns a (start, end) pair.
resI <- suppressMessages(ring_cpd(X, alternative = "interval", B = 0))
stopifnot(length(resI$scanZ$TR$tau) == 2)

## 5. Xuming Guo's scale-change example. The max-type and weighted scans
##    have different maxima, so using WR.max for MR.cor is detectable.
set.seed(1234567)
n <- 100L
generate_ar1 <- function() {
  x <- matrix(rnorm(n * n), nrow = n)
  for (j in 2:n) {
    x[, j] <- 0.3 * x[, j - 1L] + sqrt(1 - 0.3^2) * x[, j]
  }
  x
}
base_data <- replicate(5L, generate_ar1(), simplify = FALSE)
for (sigma in c(2, 2.5)) {
  for (x in base_data) {
    x[51:100, ] <- sqrt(sigma) * x[51:100, ]
    R <- ring_graph(as.matrix(dist(x)), k = floor(n^0.65))
    fit <- suppressMessages(rcpd(R, n0 = 10, n1 = 90, B = 0))
    stopifnot(fit$scanZ$MR$tau == 50L,
              fit$scanZ$MR$MR.max > fit$scanZ$WR$WR.max,
              fit$pval.appr$MR.cor > 0, fit$pval.appr$MR.cor < 1e-6)
  }
}

## 6. Plant the scale change within an interval. Check both corrected tails
##    at the MR maximum and the standalone weighted tail at its own maximum.
x <- base_data[[1L]]
x[26:75, ] <- sqrt(2.5) * x[26:75, ]
R <- ring_graph(as.matrix(dist(x)), k = floor(n^0.65))
fit <- suppressMessages(rcpd(R, n0 = 10, n1 = 90, B = 0,
                             alternative = "interval"))
b.mr <- fit$scanZ$MR$MR.max
b.wr <- fit$scanZ$WR$WR.max
stopifnot(b.mr > b.wr)
edges_by_node <- lapply(seq_len(n), function(i) {
  neighbors <- which(R[i, ] > 0)
  list(node = neighbors, weight = R[i, neighbors])
})
moments <- ringSeg:::ThirdMoment(R, edges_by_node)
t <- seq_len(n - 1L)
x.w <- ringSeg:::rho_one_Rw(n, t)
x.w[is.infinite(x.w)] <- 0
x.diff <- n / (2 * t * (n - t))
p.w.mr <- suppressMessages(ringSeg:::pI_correct_w(n, b.mr, moments$EZw3, x.w, 10, 90))
p.w.wr <- suppressMessages(ringSeg:::pI_correct_w(n, b.wr, moments$EZw3, x.w, 10, 90))
p.diff <- suppressMessages(ringSeg:::pI_correct_diff(n, b.mr, moments$EZdiff3, x.diff, 10, 90))
stopifnot(abs(p.w.mr - p.w.wr) > 0.01,
          isTRUE(all.equal(fit$pval.appr$MR.cor, p.w.mr + (1 - p.w.mr) * p.diff)),
          isTRUE(all.equal(fit$pval.appr$WR.cor, p.w.wr)))

## On a regular cycle, the weighted component still contributes to MR and TR.
n <- 12L
cycle <- matrix(0, n, n)
cycle[cbind(seq_len(n), c(2:n, 1L))] <- 1
cycle <- cycle + t(cycle)
flat <- matrix(1, n, n)
diag(flat) <- 0
for (alt in c("single", "interval")) {
  set.seed(23)
  fit <- suppressWarnings(suppressMessages(rcpd(cycle, B = 19, alternative = alt)))
  stopifnot(all(is.finite(unlist(fit$pval.perm))),
            all(is.na(unlist(fit$pval.appr))),
            all(is.finite(fit$scanZ$MR$tau)),
            isTRUE(all.equal(fit$perm.dat$MR, pmax(0, fit$perm.dat$WR))))
  fit <- suppressMessages(rcpd(flat, pval.appr = FALSE, B = 9, alternative = alt))
  stopifnot(all(unlist(fit$pval.perm) == 1),
            all(vapply(fit$scanZ, function(z) all(is.na(z$tau)), logical(1))))
}

## All eligible WR statistics are negative. The earlier implementation
## selected a length-one interval with a zero score from outside the search.
R <- as.matrix(read.table(text = "
0 9 18 17 4 6 10 10
9 0 12 5 6 10 10 10
18 12 0 9 14 2 10 5
17 5 9 0 11 15 12 15
4 6 14 11 0 12 13 11
6 10 2 15 12 0 3 11
10 10 10 12 13 3 0 6
10 10 5 15 11 11 6 0"))
set.seed(32)
fit <- suppressWarnings(rcpd(R, n0 = 6, n1 = 6, B = 79, alternative = "interval"))
stopifnot(diff(fit$scanZ$WR$tau) == 6, fit$scanZ$WR$WR.max < 0,
          fit$pval.perm$WR == 0.875, all(is.na(unlist(fit$pval.appr))))
fit <- suppressWarnings(rcpd(R, n0 = 4, n1 = 4, B = 19))
stopifnot(fit$scanZ$MR$tau == 4, all(is.na(unlist(fit$pval.appr))),
          all(is.finite(unlist(fit$pval.perm))))

## Ranking raw distances preserves distinct near neighbours and treats ties symmetrically.
D <- matrix(c(0, 1, 2, 1, 0, 1e18, 2, 1e18, 0), 3, 3)
stopifnot(identical(ring_graph(D, 1), matrix(c(0, 2, 1, 2, 0, 0, 1, 0, 0), 3, 3)))
D <- matrix(1, 6, 6)
diag(D) <- 0
ids <- c(6, 3, 1, 4, 2, 5)
stopifnot(identical(ring_graph(D[ids, ids], 3), ring_graph(D, 3)[ids, ids]))
for (scale in c(1e-200, 1e200)) {
  ref <- rcpd(R, n0 = 2, n1 = 6, pval.appr = FALSE)
  fit <- rcpd(scale * R, n0 = 2, n1 = 6, pval.appr = FALSE)
  stopifnot(isTRUE(all.equal(ref$scanZ, fit$scanZ)))
}
## A tie across the cutoff shares the retained weights without inventing
## a preference based on observation order.
group <- rep(1:2, each = 4L)
D <- outer(group, group, function(a, b) 1 + (a != b))
diag(D) <- 0
ids <- c(2, 5, 1, 8, 4, 7, 3, 6)
inside <- outer(group, group, "==") & row(D) != col(D)
for (k in c(1, 2, 3, 4)) {
  R <- ring_graph(D, k)
  stopifnot(all(R[inside] > 0),
            isTRUE(all.equal(rowSums(R), rep(k * (k + 1), 8L))),
            identical(R[ids, ids], ring_graph(D[ids, ids], k)))
}
stopifnot(isTRUE(all.equal(ring_graph(D, 1)[inside], rep(2/3, sum(inside)))),
          all(ring_graph(D, 1)[!inside] == 0))

## n = 6 can leave no valid right-tail extrapolation anchors. Analytic
## fallback and permutation-only calls use the same permutation draws.
set.seed(4)
X <- matrix(rnorm(18), 6, 3)
for (alt in c("single", "interval")) {
  set.seed(91)
  fit <- suppressMessages(ring_cpd(X, B = 19, alternative = alt))
  set.seed(91)
  perm <- suppressMessages(ring_cpd(X, pval.appr = FALSE, B = 19, alternative = alt))
  stopifnot(all(is.finite(unlist(fit$pval.appr))),
            identical(fit$pval.perm, perm$pval.perm),
            identical(fit$perm.dat, perm$perm.dat))
}

## Preserve the public order used by 0.1.0, including analytic fallback.
p_fields <- c("TR", "MR", "WR", "MR.cor", "WR.cor")
for (alt in c("single", "interval")) {
  fit <- suppressMessages(ring_cpd(X, alternative = alt))
  fixed <- suppressWarnings(ring_cpd(X, n0 = 2, n1 = 2, alternative = alt))
  stopifnot(identical(names(fit$pval.appr), p_fields),
            identical(names(fixed$pval.appr), p_fields))
}

## The skewness correction tends to the Gaussian factor as skewness tends
## to zero. Each scan integrates over the same grid up to quadrature error.
n <- 50L
b <- 3
t <- seq_len(n - 1L)
x.w <- ringSeg:::rho_one_Rw(n, t)
x.w[!is.finite(x.w)] <- 0
x.diff <- n / (2 * t * (n - t))
bounds <- c(5L, 45L)
case <- list(
  list(ringSeg:::pS_correct_w, ringSeg:::pWR1, x.w),
  list(ringSeg:::pI_correct_w, ringSeg:::pWR2, x.w),
  list(ringSeg:::pS_correct_diff, ringSeg:::pDR1, x.diff),
  list(ringSeg:::pI_correct_diff, ringSeg:::pDR2, x.diff)
)
for (parts in case) {
  gaussian <- parts[[2]](b, n, bounds[1], bounds[2])
  corrected <- vapply(c(0, -1e-17, 1e-17), function(skewness) {
    parts[[1]](n, b, rep(skewness, n - 1L), parts[[3]], bounds[1], bounds[2])
  }, numeric(1))
  stopifnot(all(is.finite(corrected)), all(corrected > 0),
            max(abs(corrected - corrected[1])) < 1e-10,
            abs(corrected[1] / gaussian - 1) < 0.05)
}

## Preserve scans and permutation counts when a large constant edge weight
## is added. The expected counts come from enumerating subsets for the
## permutation moments and rescanning each permuted graph directly.
set.seed(92)
n <- 8L
weights <- matrix(sample(0:10, n*n, replace = TRUE), n)
weights <- weights + t(weights)
for (spread in c(1e-10, 1e-13)) {
  weights_shifted <- 1 + spread * weights
  diag(weights_shifted) <- 0
  for (alt in c("single", "interval")) {
    set.seed(123)
    fit <- rcpd(weights_shifted, n0 = 2, n1 = 6, B = 49,
                pval.appr = FALSE, alternative = alt)
    expected <- if (alt == "single") c(.32, .36, .22) else c(.88, .62, .56)
    stopifnot(identical(unlist(fit$pval.perm, use.names = FALSE), expected))
  }
}

## An isolated invalid saddlepoint may not trigger extrapolation. It must
## still be excluded, with a fallback to the uncorrected component tail.
n <- 50L
t <- seq_len(n - 1L)
w <- ringSeg:::rho_one_Rw(n, t)
w[!is.finite(w)] <- 0
d <- n / (2*t*(n-t))
r <- replace(numeric(n - 1L), 2L, -.2)
b <- 3
cases <- list(
  list(ringSeg:::pS_correct_w, ringSeg:::pWR1, w),
  list(ringSeg:::pI_correct_w, ringSeg:::pWR2, w),
  list(ringSeg:::pS_correct_diff, ringSeg:::pDR1, d),
  list(ringSeg:::pI_correct_diff, ringSeg:::pDR2, d)
)
for (parts in cases) {
  observed <- suppressMessages(parts[[1]](n, b, r, parts[[3]], 2, 45))
  stopifnot(identical(observed, parts[[2]](b, n, 2, 45)))
}

## Arrange an almost-zero WR maximum at two adjacent split sizes.
## The old small-threshold approximation reported a spurious small p-value.
set.seed(338)
n <- 8L
edges <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
for (alt in c("single", "interval")) {
  splits <- if (alt == "single") cbind(0L, 2:3) else edges[diff(t(edges)) %in% 2:3, ]
  A <- apply(splits, 1, function(split) {
    k <- diff(split)
    i <- edges[, 1] > split[1] & edges[, 1] <= split[2]
    j <- edges[, 2] > split[1] & edges[, 2] <= split[2]
    2*((n-k-1)*(i & j)+(k-1)*(!i & !j))/(n-2) - 2*(k-1)*(n-k-1)/((n-1)*(n-2))
  })
  weights <- 4 + qr.resid(qr(A), rnorm(nrow(edges))) + 1e-6*A[, 1]
  R <- matrix(0, n, n)
  R[edges] <- weights
  R <- R + t(R)
  set.seed(19)
  fit <- suppressWarnings(suppressMessages(rcpd(R, n0 = 2, n1 = 3, B = 19,
                                               alternative = alt)))
  set.seed(19)
  perm <- rcpd(R, n0 = 2, n1 = 3, B = 19, pval.appr = FALSE,
               alternative = alt)
  stopifnot(fit$scanZ$WR$WR.max > 0, fit$scanZ$WR$WR.max < 1e-3,
            is.na(fit$pval.appr$WR), is.na(fit$pval.appr$WR.cor),
            identical(names(fit$pval.appr), p_fields),
            identical(fit$scanZ, perm$scanZ),
            identical(fit$pval.perm, perm$pval.perm),
            identical(fit$perm.dat, perm$perm.dat))
}

## The product of the two survival complements rounds to one when
## both component tails are very small, but their union is still positive.
b <- 10
n <- 150L
for (parts in list(list(ringSeg:::pMR1, ringSeg:::pWR1, ringSeg:::pDR1),
                   list(ringSeg:::pMR2, ringSeg:::pWR2, ringSeg:::pDR2))) {
  weighted <- parts[[2]](b, n, 15, 135)
  difference <- parts[[3]](b, n, 15, 135)
  actual <- parts[[1]](b, n, 15, 135)
  expected <- -expm1(log1p(-weighted) + log1p(-difference))
  stopifnot(is.finite(actual), actual > 0, actual < .Machine$double.eps,
            abs(actual / expected - 1) < 1e-14)
}

cat("ringSeg: all tests passed.\n")
