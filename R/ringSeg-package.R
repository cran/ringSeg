#' ringSeg: rank-based asymptotic distribution-free change-point detection
#'
#' Implements the RING change-point detector (Zhou & Chen, 2025, IEEE TIT): a
#' rank-based, asymptotic distribution-free scan for a single change-point or a
#' changed interval in a sequence of (possibly high-dimensional / non-Euclidean)
#' observations. The detector operates on a rank matrix built from a sparse
#' similarity graph (the row-rank k-nearest-neighbour graph) and returns
#' the weighted (WR), max-type (MR), and generalized (TR) statistics with analytic
#' distribution-free p-value approximations (optionally skewness-corrected) and
#' optional permutation p-values.
#'
#' Main entry points: \code{\link{ring_cpd}} (data / distance in),
#' \code{\link{rcpd}} (rank matrix in), \code{\link{ring_graph}} (build the rank
#' graph), \code{\link{Rise_Rank}} (low-level ranking of a similarity matrix).
#'
#' @keywords internal
"_PACKAGE"

#' Build a RING rank graph from a distance matrix
#'
#' Constructs the sparse, rank-based k-nearest-neighbour graph that the RING
#' statistics operate on -- RING's ranking scheme: rank each observation's
#' neighbours within its row of the similarity \eqn{S = \max(D) - D}, keep the
#' \code{k} nearest (weights \code{k}..1 from nearest to k-th), and symmetrize.
#'
#' @param D an \eqn{n \times n} distance matrix.
#' @param k number of nearest neighbours (default \code{round(n^0.65)}, a value in
#'   the k-range used by Zhou & Chen).
#' @return an \eqn{n \times n} rank matrix with zero diagonal, suitable as the
#'   \code{R} argument of \code{\link{rcpd}}.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 5), 60, 5)
#' R <- ring_graph(as.matrix(dist(X)), k = 8)
#' @seealso \code{\link{ring_cpd}}, \code{\link{rcpd}}, \code{\link{Rise_Rank}}
#' @export
ring_graph <- function(D, k = NULL) {
  D <- as.matrix(D)
  n <- nrow(D)
  if (ncol(D) != n) stop("'D' must be a square n x n distance matrix.")
  if (n < 3) stop("'D' must have at least 3 observations.")
  if (!is.numeric(D) || any(!is.finite(D)))
    stop("'D' must be a finite numeric distance matrix (no NA/NaN/Inf).")
  if (any(D < 0))
    warning("'D' has negative entries; treating it as a general (dis)similarity. ",
            "RING ranks neighbours within each row, so this is allowed but unusual ",
            "for a distance matrix.")
  if (is.null(k)) k <- round(n^0.65)              # RING / RISE recommendation for the k-NNG
  if (length(k) != 1L || !is.numeric(k) || is.na(k))
    stop("'k' must be a single number.")
  if (k < 1 || k > n - 1)
    stop("'k' must be between 1 and n-1 (", n - 1, ").")
  if (k != round(k)) {
    warning("'k' rounded to the nearest integer (", round(k), ").")
    k <- round(k)
  }
  S <- max(D) - D                                  # similarity (rank-equiv to -D)
  Rrow <- Rise_Rank(S, method = "row")             # 0..n-1 row ranks (self lowest)
  R <- Rrow - n + 1 + k                             # keep top-k: weights k..1, else <=0
  R[R < 0] <- 0
  R <- R + t(R)
  diag(R) <- 0
  R
}

#' RING change-point detection from data or a distance matrix
#'
#' Convenience wrapper: builds the RING rank graph from the input (unless a rank
#' matrix is supplied) and calls \code{\link{rcpd}}. Observations must be in time
#' order.
#'
#' @param x One of: an \eqn{n \times d} data matrix (rows = time-ordered
#'   observations; the default), an \eqn{n \times n} distance matrix
#'   (\code{is.distance = TRUE}), or an \eqn{n \times n} rank matrix
#'   (\code{is.rank = TRUE}).
#' @param is.distance logical; if \code{TRUE}, \code{x} is a distance matrix.
#' @param is.rank logical; if \code{TRUE}, \code{x} is already a rank matrix (e.g.
#'   from \code{\link{ring_graph}}) and is used as-is.
#' @param k number of nearest neighbours for the RING k-NN graph (default
#'   \code{round(n^0.65)}).
#' @param dist.method distance passed to \code{\link[stats]{dist}} when \code{x} is a
#'   data matrix (default \code{"euclidean"}).
#' @param n0,n1 scan range (defaults \code{ceiling(0.05*n)} / \code{floor(0.95*n)};
#'   for \code{alternative = "interval"} these are the minimum / maximum interval
#'   lengths).
#' @param pval.appr logical; compute analytic p-value approximations.
#' @param skew.corr logical; also compute the skewness-corrected approximation.
#' @param B number of permutations for permutation p-values (\code{0} = none).
#' @param alternative \code{"single"} (a single change-point) or \code{"interval"}
#'   (a changed interval).
#'
#' @return A list with \code{scanZ} (the WR/MR/TR statistics and estimated
#'   change-point(s)), and, as requested, \code{pval.appr} and \code{pval.perm}.
#'   See \code{\link{rcpd}}.
#'
#' @references Zhou, D. and Chen, H. (2025). Asymptotic Distribution-Free
#'   Change-Point Detection for Modern Data Based on a New Ranking Scheme.
#'   \emph{IEEE Transactions on Information Theory}, \bold{71}(8), 6183--6197.
#'   \doi{10.1109/TIT.2025.3575858}.
#'
#' @examples
#' set.seed(1)
#' n <- 200; d <- 10; tau <- 100
#' X <- matrix(rnorm(n * d), n, d)
#' X[(tau + 1):n, ] <- X[(tau + 1):n, ] * 1.6   # a dispersion change at t = 100
#' res <- ring_cpd(X, skew.corr = TRUE)
#' res$scanZ$TR$tau          # estimated change-point
#' res$pval.appr$TR          # generalized-statistic p-value
#'
#' @seealso \code{\link{ring_graph}}, \code{\link{rcpd}}, \code{\link{Rise_Rank}}
#' @export
ring_cpd <- function(x, is.distance = FALSE, is.rank = FALSE, k = NULL,
                     dist.method = "euclidean",
                     n0 = NULL, n1 = NULL, pval.appr = TRUE, skew.corr = TRUE,
                     B = 0, alternative = c("single", "interval")) {
  alternative <- match.arg(alternative)
  if (is.rank) {
    R <- as.matrix(x)
  } else {
    x <- as.matrix(x)
    if (!is.numeric(x) || any(!is.finite(x)))
      stop("'x' must be finite and numeric (no NA/NaN/Inf).")
    D <- if (is.distance) x else as.matrix(stats::dist(x, method = dist.method))
    R <- ring_graph(D, k = k)
  }
  rcpd(R, n0 = n0, n1 = n1, pval.appr = pval.appr, skew.corr = skew.corr,
       B = B, alternative = alternative)
}
