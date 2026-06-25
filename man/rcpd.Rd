\name{rcpd}
\alias{rcpd}
\title{RING rank-based change-point detection (rank matrix in)}
\usage{
rcpd(R, n0 = NULL, n1 = NULL, pval.appr = TRUE, skew.corr = TRUE,
     B = 0, alternative = c("single", "interval"))
}
\arguments{
\item{R}{an \eqn{n \times n} rank matrix with zero diagonal, as produced by
\code{\link{ring_graph}} (the RING k-NN rank graph). Observations (nodes) are
indexed in time order.}
\item{n0, n1}{scan range; defaults \code{ceiling(0.05*n)} and
\code{floor(0.95*n)}, clamped to \code{[2, n-2]}. For
\code{alternative = "interval"} these are the minimum (\code{l0}) and maximum
(\code{l1}) interval lengths.}
\item{pval.appr}{logical; compute analytic distribution-free p-value
approximations for the WR, MR and TR statistics.}
\item{skew.corr}{logical; additionally compute the skewness-corrected
approximation (\code{WR.cor}, \code{MR.cor}). The skewness-corrected p-values
are better calibrated at finite \eqn{n}; the uncorrected Gaussian/field
approximation can be anti-conservative for small \eqn{n}.}
\item{B}{number of permutations for permutation p-values (\code{0} = none).}
\item{alternative}{\code{"single"} for a single change-point, or
\code{"interval"} for a changed interval.}
}
\value{
A list containing
\item{scanZ}{the scan statistics. For \code{alternative = "single"}, each of
\code{TR} (generalized, \eqn{Z_w^2 + Z_{\mathrm{diff}}^2}), \code{MR} (max-type,
\eqn{\max(|Z_{\mathrm{diff}}|, Z_w)}) and \code{WR} (weighted, \eqn{Z_w}) holds
the full scan curve, the estimated change-point \code{tau}, and the maximum.}
\item{pval.appr}{(if \code{pval.appr}) analytic p-values \code{TR}, \code{MR},
\code{WR}, and (if \code{skew.corr}) the skewness-corrected \code{WR.cor},
\code{MR.cor}.}
\item{pval.perm}{(if \code{B > 0}) permutation p-values for TR, MR, WR.}
}
\description{
Core RING detector operating directly on a rank matrix. The analytic
p-value approximations are distribution-free and require no permutation;
they assume the sparse RING k-NN graph, so build \code{R} with
\code{\link{ring_graph}}, or call \code{\link{ring_cpd}} to build it for you.
}
\examples{
set.seed(1)
n <- 150; d <- 8
X <- matrix(rnorm(n * d), n, d)
X[76:n, ] <- X[76:n, ] + 0.8
R <- ring_graph(as.matrix(dist(X)), k = 13)
res <- rcpd(R, B = 0)
res$scanZ$TR$tau
res$pval.appr$TR
}
\references{
Zhou, D. and Chen, H. (2025). Asymptotic Distribution-Free Change-Point
Detection for Modern Data Based on a New Ranking Scheme.
\emph{IEEE Transactions on Information Theory}, \bold{71}(8), 6183--6197.
\doi{10.1109/TIT.2025.3575858}.
}
\seealso{\code{\link{ring_cpd}}, \code{\link{Rise_Rank}}}
