#' Migration and Range Change Analysis in R
#' 
#'  A collection of functions for performing a migration and range change analysis (MRSA) as described in by Gurarie et al. (2017).  The key features are estimation of precise times, distances, and locations of a one or two step range shift in movement data. 
#' 
#' @details Some key functions for using \code{marcher} are:
#' 
#' 1. \code{\link{estimate_shift}} {Estimate a range shift process.}
#' 
#' 2. \code{\link{simulate_shift}} {Simulate a range shift process.}
#' 
#' 3. \code{\link{plot.shiftfit}} {Visualize a range shift process.}
#' 
#' 4. \code{\link{test_rangeshift}} {Test whether a range shift occurred.}
#' 
#' 5. \code{\link{test_return}} {Test whether a migration was a return migration.}
#' 
#' 6. \code{\link{test_stopover}} {Test whether a stopover occurred during a migration.}
#' 
#' Several simulated datasets are in the \code{\link{SimulatedTracks}} data object.  
#' 
#' One roe deer (\emph{Capreolus capreolus}) track is in the \code{\link{Michela}} object.
#' 
#' See the respective help files and \code{vignette("marcher")} for more details and examples. 
#' 
#' @references Gurarie, E., F. Cagnacci, W. Peters, C. Fleming, J. Calabrese, T. Mueller and W. Fagan (2017)  A framework for modeling range shifts and migrations: asking whether, whither, when, and will it return. \emph{Journal of Animal Ecology}, 86(4):943-59. DOI: 10.1111/1365-2656.12674 
#' @seealso 
#' Useful links:
#' \itemize{
#'  \item Report bugs at \url{https://github.com/EliGurarie/marcher/issues}
#' }
#' @author \strong{Maintainer}: Eliezer Gurarie \email{egurarie@umd.edu}
#' 
#' Other contributors:
#' \itemize{\item Faridedin Cheraghi [contributor]}
"_PACKAGE"
