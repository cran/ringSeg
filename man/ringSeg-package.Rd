\name{ringSeg-package}
\alias{ringSeg-package}
\alias{ringSeg}
\docType{package}
\title{ringSeg: rank-based asymptotic distribution-free change-point detection}
\description{
Implements the RING change-point detector (Zhou & Chen, 2025, IEEE TIT): a
rank-based, asymptotic distribution-free scan for a single change-point or a
changed interval in a sequence of (possibly high-dimensional / non-Euclidean)
observations. The method operates on a rank matrix built from a pairwise
similarity and returns the weighted (WR), max-type (MR), and generalized (TR)
statistics with analytic distribution-free p-value approximations (optionally
skewness-corrected) and optional permutation p-values.
}
\details{
Main entry points: \code{\link{ring_cpd}} (data / distance / rank-matrix in),
\code{\link{rcpd}} (rank matrix in), \code{\link{ring_graph}} (build the RING
k-NN rank graph), \code{\link{Rise_Rank}} (low-level ranker).
}
\references{
Zhou, D. and Chen, H. (2025). Asymptotic Distribution-Free Change-Point
Detection for Modern Data Based on a New Ranking Scheme.
\emph{IEEE Transactions on Information Theory}, \bold{71}(8), 6183--6197.
\doi{10.1109/TIT.2025.3575858}.
}
\author{Doudou Zhou, Hao Chen (maintainer)}
