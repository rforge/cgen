\name{set_num_threads}
\alias{set_num_threads}
\title{
Set the number of OpenMP threads used by the functions of package \code{cgen}
}
\description{
This function sets the value of the global variable stored in \code{options()$cgen.threads} to the assigned integer.
Note: The assigned value may exceed the number of physical cores present but that might lead to dramatical
decrease in performance.
}
\usage{
set_num_threads(x)
}

\arguments{
\item{x}{Integer scalar that specifies the number of threads to be used by \code{cgen}-functions}
}

\value{
Changes the global variable \code{cgen.threads} to the value in \code{x}
}

\seealso{\code{\link{get_num_threads}, \link{get_max_threads}, \link{check_openmp}}}

\examples{

# Control the number of threads being used in an R-session:

# set the number of threads to 1
set_num_threads(1)

# Use a parallelized cgen-function

data(QM16)

# check single-threaded performance
system.time(W <- M\%c\%t(M))

# set number of threads to the value reported by get_max_threads()

n_threads <- get_max_threads()
set_num_threads(n_threads)

# check multi-threaded performance
system.time(W <- M\%c\%t(M))


}

\keyword{Parallelization}


