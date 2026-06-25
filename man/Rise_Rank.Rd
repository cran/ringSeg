\name{Rise_Rank}
\alias{Rise_Rank}
\title{Build a rank matrix from a similarity matrix}
\usage{
Rise_Rank(S, method = c("overall", "row"))
}
\arguments{
\item{S}{an \eqn{n \times n} similarity matrix (larger = more similar; e.g.
the negative of a distance matrix).}
\item{method}{\code{"overall"} ranks all distinct pairs jointly (the symmetric
rank matrix is returned); \code{"row"} ranks within each row.}
}
\value{
An \eqn{n \times n} rank matrix with zero diagonal.
}
\description{
Low-level ranker for the RING statistics: turns a similarity matrix into a
rank matrix. \code{method = "row"} is the per-row ranking used to build the
sparse RING k-NN graph. For change-point detection, build the \code{R} argument
of \code{\link{rcpd}} with \code{\link{ring_graph}} (the validated sparse graph);
the analytic p-value theory assumes that graph rather than a dense
\code{"overall"} ranking.
}
\examples{
set.seed(1)
X <- matrix(rnorm(50 * 5), 50, 5)
R <- Rise_Rank(-as.matrix(dist(X)), method = "overall")
dim(R)
}
\seealso{\code{\link{rcpd}}, \code{\link{ring_cpd}}}
