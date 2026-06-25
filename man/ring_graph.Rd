\name{ring_graph}
\alias{ring_graph}
\title{Build a RING rank graph from a distance matrix}
\usage{
ring_graph(D, k = NULL)
}
\arguments{
\item{D}{an \eqn{n \times n} distance matrix.}
\item{k}{number of nearest neighbours (default \code{round(n^0.65)}, a value in the
k-range used by Zhou & Chen).}
}
\value{
An \eqn{n \times n} rank matrix with zero diagonal, suitable as the \code{R}
argument of \code{\link{rcpd}}.
}
\description{
Constructs the sparse, rank-based k-nearest-neighbour graph that the RING
statistics operate on -- RING's ranking scheme: rank each observation's
neighbours within its row of the similarity \eqn{S = \max(D) - D}, keep the \code{k}
nearest (weights \code{k}..1 from nearest to k-th), and symmetrize.
}
\examples{
set.seed(1)
X <- matrix(rnorm(60 * 5), 60, 5)
R <- ring_graph(as.matrix(dist(X)), k = 8)
dim(R)
}
\seealso{\code{\link{ring_cpd}}, \code{\link{rcpd}}, \code{\link{Rise_Rank}}}
