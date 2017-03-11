#' Interactive locating of range shifting
#' 
#' Plots an x-y, time-x, time-y track of a potential migration process and prompts the user to click on the figure to obtain initial estimates of range centroids and timing of start and end of migrations. 
#'
#' @param time time (can be a \code{\link[base]{POSIXt}})
#' @param x x and y coordinates.  Can be two separate vectors OR a complex "x" OR a two-column matrix/date-frame. 
#' @param y see x
#' @param n.clust number of ranges (either 2 or 3)
#' @param ... additional parameters to pass to plot functions
#' @return a named vector of initial estimates: 
#' if \code{n.clust = 2}, c(x1, x2, y1, y2, t1, dt) 
#' if \code{n.clust = 3}, c(x1, x2, x3, y1, y2, y3, t1, t2, dt1, dt2)
#' @seealso \code{\link{quickfit}}, code{\link{locator}}
#' @export
locate_shift <- function(time, x, y, n.clust = 2, ...){
  
  layout(rbind(c(1,2), c(1,3)))

  if(n.clust == 2){
		plot(x,y,asp=1, type="o", pch=19, col=rgb(0,0,0,.5), cex=0.5, main = "click on two centroids ", ...)
		points(x[c(1, length(x))], y[c(1, length(x))], col= c("green", "red"), pch = c(19, 20), cex=2)
		xy <- locator(2)
		
		plot(time,x, type="o", pch=19, col=rgb(0,0,0,.5), xaxt="n", xlab="", cex=0.5, main = "click on migration start and end times", ...)
		plot(time,y, type="o", pch=19, col=rgb(0,0,0,.5), cex=0.5, ...)
		ts <- locator(2)
		
		p.m0 <- vector()
		p.m0['x1'] <- xy$x[1]
		p.m0['x2'] <- xy$x[2]
		
		p.m0['y1'] <- xy$y[1]
		p.m0['y2'] <- xy$y[2]
		
		p.m0['t1'] <- ts$x[1]
		p.m0['dt'] <- ts$x[2]-ts$x[1]	
	}  else if(n.clust == 3){
			
		plot(x,y,asp=1, type="o", pch=19, col=rgb(0,0,0,.5), cex=0.5, main = "click on three centroids ", ...)
		points(x[c(1, length(x))], y[c(1, length(x))], col= c("green", "red"), pch = c(19, 20), cex=2)
		xy <- locator(3)
		
		plot(time,x, type="o", pch=19, col=rgb(0,0,0,.5), xaxt="n", xlab="", cex=0.5, main = "click on start and end of migration 1 and migration 2", ...)
		plot(time,y, type="o", pch=19, col=rgb(0,0,0,.5), cex=0.5, ...)
		ts <- locator(4)
		
		p.m0 <- vector()
		p.m0['x1'] <- xy$x[1]
		p.m0['x2'] <- xy$x[2]
		p.m0['x3'] <- xy$x[3]
		
		p.m0['y1'] <- xy$y[1]
		p.m0['y2'] <- xy$y[2]
		p.m0['y3'] <- xy$y[3]
		
		p.m0['t1'] <- ts$x[1]
		p.m0['t2'] <- ts$x[3]
		p.m0['dt1'] <- ts$x[2] - ts$x[1]
		p.m0['dt2'] <- ts$x[4] - ts$x[3]
	}
  return(p.m0)
}

