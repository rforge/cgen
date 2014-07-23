\name{\%c\%}
\alias{\%c\%}
\title{
(Parallel) crossproduct operator
}
\description{
This operator computes the crossproduct between two matrices.
It can be used as a replacement for \code{\%*\%} in most cases. The operator
only accepts matrices of types: \code{matrix} or \code{dgCMatrix}. 
In the case of two dense matrices the operator will compute the crossproduct 
in parallel (Eigen + OpenMP)
}
\usage{
W <- X\%c\%Z	
}

\arguments{
\item{}{Two matrices or vectors (treated as column-vectors) of type: \code{matrix} or \code{dgCMatrix}}
}

\value{
Matrix of type: \code{matrix} or \code{dgCMatrix}
}

\examples{
## Not run:

# Least Squares Solving

# Generate random data 

n = 1000
p = 500

M <- matrix(rnorm(n*p),n,p)
y <- rnorm(n)

# least squares solution:

b <- csolve(t(M) \%c\% M, t(M) \%c\% y)

## End(Not run)
}

\keyword{Tools}


