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
statistics operate on. Rank each observation's neighbours by increasing
distance, give the \code{k} nearest weights \code{k}..1 when the distances
are distinct, and symmetrize. Tied neighbours share the mean of the
weights at their positions, including zero weight beyond \code{k}.
A tie at the cutoff can retain more than \code{k} neighbours. The sum of
outgoing weights before symmetrization is \eqn{k(k+1)/2} for each node.
}
\examples{
set.seed(1)
X <- matrix(rnorm(60 * 5), 60, 5)
R <- ring_graph(as.matrix(dist(X)), k = 8)
dim(R)
}
\seealso{\code{\link{ring_cpd}}, \code{\link{rcpd}}, \code{\link{Rise_Rank}}}
