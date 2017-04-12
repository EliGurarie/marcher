#' Movement track of Michela, a roe deer
#'
#' GPS tracks of one roe deer (\emph{Capreolus capreolus}) in the Italian alps. 
#' This deer performs two seasonal migrations, from a wintering ground to a summering ground, 
#' back its wintering ground.  For several ways to analyze these data, see examples in the 
#' \link{marcher} vignette.
#'
#' @usage 
#' data(Michela)
#' 
#' @format Data frame containing movements of roe deer with the following columns:
#' \describe{
#'   \item{id}{ID of animal}
#'   \item{name}{Names - for mnemonic convenience - of Italian authors.}
#'   \item{x,y}{In Easting Westing}
#'   \item{latitude, longitude}{}
#'   \item{time}{POSIXct object}
#'   \item{day}{Day of year, counting from January 1 of the first year of observations 
#'   (thus day 367 is January 2 or the following year).}}
#' 
#' @examples 
#' data(Michela)
#' with(Michela, scan_track(time = time, x = x, y = y))
#' 
#' @keywords data
#' @references 
#' For more details, see: \url{Eurodeer.org}
#'
"Michela"

#' Simulated range shift tracks
#' 
#' Five simulated tracks: \code{MWN.sim},  \code{MOU.sim}, \code{MOUF.sim} are simulated two-range shifts with different levels of position and velocity autocorrelation, \code{MOUF.sim.random} which has 100 observations random times, and \code{MOU.3range} which is a MOU process with two range shifts (and 200 observations).
#' 
#' @format Each of these is a data frame with 100 observations of 
#' three numeric variables (except for MOU.3range, which has 200 observations).  
#' The columns are: \code{T}, \code{X}, \code{Y}.
#' 
#' @details The data frames are also \code{track} class object frame.
#' \describe{
#'   \item{MOU.3range}{Simulated migratory Ornstein-Uhlenbeck with 3 range}
#'   \item{MOU.sim}{Simulated migratory Ornstein-Uhlenbeck}
#'   \item{MOUF.sim}{Simulated migratory Ornstein-Uhlenbeck Flemming}
#'   \item{MOUF.sim.random}{Simulated migratory Ornstein-Uhlenbeck Flemming at random or arbitrary times of observation}
#'   \item{MWN.sim}{Simulated migratory white noise ranging model}
#'   }
#' @source Code to simulate tracks like these are provided in the marcher vignette. 
#' @keywords data
#' @usage data("SimulatedTracks")
#' @examples
#' data(SimulatedTracks)
#' scan_track(MWN.sim)
#' scan_track(MOU.sim)
#' scan_track(MOUF.sim)
#' scan_track(MOUF.sim.random)
#' scan_track(MOU.3range)
#' @docType data
#' @name SimulatedTracks
#' @aliases MOU.3range MOU.sim MOUF.sim MOUF.sim.random MWN.sim
NULL