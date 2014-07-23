\name{cscanx}
\alias{cscanx}
\title{
Read in a matrix from a file
}
\description{
Reads in a matrix from file (no header, no row-names, no NA's, space or tab-delimiter)
and returns the according R-matrix. No Need to specify dimensions.
}
\usage{
X <- cscanx("/path/to/file")	
}

\arguments{
\item{"/path/to/file"}{character - location of the file to be read}
}

\value{
Matrix shaped as in the file 
}

\examples{

# random matrix
X <- matrix(rnorm(10,5),10,5)

# write that matrix to a file
write.table(X,file="X",col.names=F,row.names=F,quote=F)

# read in the matrix to object Z
Z <- cscanx("X")

}

\keyword{Tools}
