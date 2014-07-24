\name{get_num_threads}
\alias{get_num_threads}
\title{
Get the number of threads for \code{cgen2}
}
\description{
Check the variable that specifies the number of threads being used by \code{cgen2}-functions
}
\usage{
get_num_threads()
}


\value{
Returns the value of the global variable \code{cgen2.threads}
}

\seealso{\code{\link{set_num_threads}, \link{get_max_threads}, \link{check_openmp}}}

\examples{

# set the number of threads to 1
set_num_threads(1)

# check
get_num_threads()

# set number of threads to the value reported by get_max_threads()
n_threads <- get_max_threads()
set_num_threads(n_threads)

# check
get_num_threads()


}

\keyword{Parallelization}


