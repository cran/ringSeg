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
R <- ring_graph(as.matrix(dist(X)), k = round(sqrt(nn)))
stopifnot(dim(R)[1] == nn, dim(R)[2] == nn, all(diag(R) == 0))
stopifnot(mean(R > 0) < 0.5)                       # k-NN graph is sparse
res2 <- suppressMessages(rcpd(R, B = 0))
stopifnot(isTRUE(all.equal(res2$scanZ$TR$tau, res$scanZ$TR$tau)))   # same graph as ring_cpd default

## 4. Interval alternative runs and returns a (start, end) pair.
resI <- suppressMessages(ring_cpd(X, alternative = "interval", B = 0))
stopifnot(length(resI$scanZ$TR$tau) == 2)

cat("ringSeg: all tests passed.\n")
