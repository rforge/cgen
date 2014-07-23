\name{\%**\%}
\alias{\%**\%}
\title{
Square matrix power operator
}
\description{
This operator computes an arbitrary power of a positive definite square matrix
using an Eigen-decomposition: \eqn{\mathbf{X}^p = \mathbf{UD}^{p}\mathbf{U}'}
}
\usage{
X_p <- X \%**\% p
}

\arguments{
\item{X}{Positive definite square matrix}
\item{p}{numeric scalar - desired power of X}
}

\value{
Matrix \code{X} to the power \code{p}
}

\examples{
## Not run:

# Inverse Square Root of a positive definite square matrix
X <- matrix(rnorm(100*5000),100,1000)

XX <- ccross(X)

XX_InvSqrt <- XX \%**\% -0.5

# check result: ((XX')^-0.5 (XX')^-0.5)^-1 = XX'
table(round(csolve(XX_InvSqrt \%c\% XX_InvSqrt),dig=2) == round(XX,dig=2) )

## End(Not run)
}

\keyword{Tools}


