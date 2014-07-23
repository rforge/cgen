\name{cmaf}
\alias{cmaf}
\title{
cmaf
}
\description{
Computes the minor allele frequencies of a marker-matrix.
}
\usage{
maf <- cmaf(X)	
}

\arguments{
  \item{X}{Marker matrix with \{-1,0,1\} coding}
}

\value{
Numeric Vector of minor allele frequencies for every column in X
}

\examples{
## Not run:

data(QM16)
mafs <- cmaf(M)

## End(Not run)
}

\keyword{Tools}
