\name{ring_cpd}
\alias{ring_cpd}
\title{RING change-point detection from data or a distance matrix}
\usage{
ring_cpd(x, is.distance = FALSE, is.rank = FALSE, k = NULL,
         dist.method = "euclidean",
         n0 = NULL, n1 = NULL, pval.appr = TRUE, skew.corr = TRUE,
         B = 0, alternative = c("single", "interval"))
}
\arguments{
\item{x}{One of: an \eqn{n \times d} data matrix (rows = time-ordered
observations; the default), an \eqn{n \times n} distance matrix
(\code{is.distance = TRUE}), or an \eqn{n \times n} rank matrix
(\code{is.rank = TRUE}).}
\item{is.distance}{logical; if \code{TRUE}, \code{x} is a distance matrix.}
\item{is.rank}{logical; if \code{TRUE}, \code{x} is already a rank matrix (e.g. from
\code{\link{ring_graph}}) and is used as-is.}
\item{k}{number of nearest neighbours for the RING k-NN graph
(default \code{round(n^0.65)}).}
\item{dist.method}{distance passed to \code{\link[stats]{dist}} when \code{x} is a
data matrix (default \code{"euclidean"}).}
\item{n0, n1}{scan range (defaults \code{ceiling(0.05*n)} / \code{floor(0.95*n)};
for \code{alternative = "interval"} these are the minimum / maximum interval
lengths).}
\item{pval.appr}{logical; compute analytic p-value approximations.
Unavailable approximations, including those below the tail-range cutoff,
are returned as \code{NA}. See \code{\link{rcpd}}.}
\item{skew.corr}{logical; also compute the skewness-corrected approximation.}
\item{B}{number of permutations for permutation p-values (\code{0} = none).}
\item{alternative}{\code{"single"} (a single change-point) or \code{"interval"}
(a changed interval).}
}
\value{
A list with \code{scanZ} (the WR / MR / TR statistics and estimated
change-point(s)) and, as requested, \code{pval.appr} and \code{pval.perm}.
For a changed interval, \code{tau = c(i, j)} identifies rows
\code{(i + 1):j}. See \code{\link{rcpd}}.
}
\description{
Convenience wrapper: builds the RING rank graph from the input (unless a rank
matrix is supplied) and calls \code{\link{rcpd}}. Observations must be in time order.
}
\examples{
set.seed(1)
n <- 200; d <- 10; tau <- 100
X <- matrix(rnorm(n * d), n, d)
X[(tau + 1):n, ] <- X[(tau + 1):n, ] * 1.6   # a dispersion change at t = 100
res <- ring_cpd(X, skew.corr = TRUE)
res$scanZ$TR$tau
res$pval.appr$TR
}
\references{
Zhou, D. and Chen, H. (2025). Asymptotic Distribution-Free Change-Point
Detection for Modern Data Based on a New Ranking Scheme.
\emph{IEEE Transactions on Information Theory}, \bold{71}(8), 6183--6197.
\doi{10.1109/TIT.2025.3575858}.
}
\seealso{\code{\link{ring_graph}}, \code{\link{rcpd}}, \code{\link{Rise_Rank}}}
