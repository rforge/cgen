\name{ccross}
\alias{ccross}
\title{
ccross
}
\description{
Computation of the following matrix-product: \eqn{\mathbf{XDX}'}
Where \eqn{\mathbf{D}} is a diagonal matrix, which is being passed to the function as a vector.
}
\usage{
W <- ccross(X,D=NULL)	
}

\arguments{
  \item{X}{matrix}
  \item{D}{numeric vector, will be used as a weighting diagonal matrix of dimension ncol(X). 
If omitted an identity matrix will be assigned.}
}

\value{
Square matrix of dimension nrow(X)
}

\examples{

# Computing the matrix-square-root of a genomic relationship matrix:
## Not run:

data(QM16)

A <- cgrm.A(M,lambda=0.01)

# this is the implementation of the matrix power-operator '\%**\%'
A_sqrt <- with(eigen(A), ccross(vectors,values**0.5))
 
## End(Not run)

}

\keyword{Tools}


