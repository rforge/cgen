.onAttach <- function( libname, pkgname ){
  options( cgenpp.threads = get_max_threads() )
}


