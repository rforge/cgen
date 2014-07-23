\name{rand_data}
\alias{rand_data}
\title{
Generate random data for test purposes
}
\description{
Generates a random marker-matrix in -1,0,1 coding and a phenotype vector.
Both objects are independent and the only purpose is to test 
functions and show examples
}
\usage{
rand_data(n,p_marker)
}

\arguments{
  \item{n}{Number of oberservations}
  \item{p_marker}{Number of markers}
}

\value{
No return value. Generates two objects globally (\code{M} and \code{y}) that can be used after the execution of the function.
\code{M} is the marker matrix and \code{y} the phenotype vector
}

\examples{

# Generate random data with 100 observations and 500 markers
rand_data(100,500)

# check that objects have been created
str(M)
str(y)

}

\keyword{Genomic Prediction}
