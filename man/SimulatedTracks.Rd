\name{SimulatedTracks}
\alias{SimulatedTracks}
\docType{data}
\title{
Simulated range shift tracks
}
\description{
Five simulated tracks: \code{MWN.sim},  \code{MOU.sim}, \code{MOUF.sim} are simulated two-range shifts with different levels of position and velocity autocorrelation, \code{MOUF.sim.random} which has 100 observations random times, and \code{MOU.3range} which is a MOU process with two range shifts (and 200 observations).
}
\usage{data("SimulatedTracks")}
\format{
  Each of these is a data frame with 100 observations of three numeric variables (except for MOU.3range, which has 200 observations).  The columns are: \code{T}, \code{X}, \code{Y}.
}
\details{
The data frames are also \code{track} class object frame.
}
\source{
Code to simulate tracks like these are provided in the marcher vignette. 
}
\references{
}
\examples{
data(SimulatedTracks)
scan.track(MWN.sim)
scan.track(MOU.sim)
scan.track(MOUF.sim)
scan.track(MOUF.sim.random)
scan.track(MOU.3range)
}
\keyword{datasets}
  