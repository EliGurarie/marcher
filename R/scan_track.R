#' scan_track
#' 
#' Plotting x-y, time-x, time-y scan of a track.  This function will take x, y, and time coordinates or a \code{track} class object
#' 
#' @param track a \code{track} class object, or any data-frame that contains (at least) three columns labeled "T", "X" and "Y"
#' @param time time (can be a \code{\link{POSIXt}})
#' @param x x Coordinate. x,y coordiantes an be two separate vectors OR a complex "x" OR a two-column matrix/date-frame. 
#' @param y y coordinate. 
#' @param layout the default layout places the x-y plot on the left and the respective 1-d time series on the right. 
#' @param auto.par by default, uses a decent looking default layout.  Otherwise can be a \code{\link{par}} list, or, e.g. FALSE to keep externally defined settings.
#' @param ... options to be passed to plot functions
#' 
#' @examples
#' 
#' ## Roe deer data
#' 
#' data(Michela)
#' par(bty="l", mar = c(0,4,0,2), oma=c(4,0,4,0), xpd=NA) 
#' with(Michela, scan_track(time = time, x = x, y = y, main="Michela"))
#' 
#' ## Simulated track
#' 
#' time <- 1:200
#' Mean <- getMu(T = time, p.m = c(x1 = 0, y1 = 0, x2 = 10, y2 = 10, t1 = 90, dt = 20))
#' SimTrack <- simulate_shift(T = time, tau = c(tau.z = 5), mu = Mean, A = 40)
#' with(SimTrack, scan_track(time = T, x = X, y = Y))
#' 
#' # OR (because SimTrack is a "track")
#' scan_track(SimTrack)
#' @export

scan_track <- function(track=NULL, time, x, y=NULL, layout = NULL, auto.par = NULL, ...)
{
  if(inherits(track, "track") | all(c("T","X","Y") %in% names(track))){
    time <- track[['T']]
    x <- track[['X']]
    y <- track[['Y']] 
  } else  if(is.null(y)) if(is.complex(x)){y <- Im(x); x <- Re(x)} else if(ncol(x) == 2){y <- x[,2]; x <- x[,1]}
  
  if(is.null(layout)) layout(rbind(c(1,2), c(1,3))) 
  if(is.null(auto.par)) par(mar = c(0,4,0,0), oma = c(4,0,4,4), xpd = NA)  else 
    par(auto.par)
  
  plot(x,y,asp=1, type="o", pch=19, col=rgb(0,0,0,.5), cex=0.5, ...)
  plot(time,x, type="o", pch=19, col=rgb(0,0,0,.5), xaxt="n", xlab="", cex=0.5, ...)
  plot(time,y, type="o", pch=19, col=rgb(0,0,0,.5), cex=0.5, ...)
}
