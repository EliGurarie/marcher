#' Plot results of an range-shift fit
#' 
#' Plotting functions for illustrating the results of a range-shift fit. 
#' 
#' @param x a fitted range shift object, i.e. output of the \code{\link{estimate_shift}}
#' @param ns a vector of 3 simulation values, useful for smoothing the bars in the dumbbell plot.  For smoothing, it might be recommended to increase the first value, \code{n.sims} - the number of draws from the fitted migation process.
#' @param plot.ts whether or not to plot the time series as well
#' @param stretch an extra parameter to extend the bars on the dumbbells (in real distance units).
#' @param CI.cols three shading colors, from lightest to darkest. The default is a sequence of blues.
#' @param layout the default layout places the x-y plot on the left and - if \code{plot.ts==TRUE} - the respective 1-d time series on the right. 
#' @param par graphics window parameters that, by default, look nice with the default layout.
#' @param pt.cex point character expansion.
#' @param pt.col points color.
#' @param ... additional parameters to pass to plot function (e.g. labels, title, etc.)
#'
#' @example ./demo/estimate_shift_example.r
#' @export	
plot.shiftfit <- function(x, ns = c(n.sims = 1e3, n.times = 1e2, n.bins = 10), 
                          plot.ts = TRUE, stretch = 0, pt.cex = 0.8, 
                          pt.col = "antiquewhite", CI.cols = NULL, 
                          layout = NULL, par = NULL, ...){	
  #x is FIT
	X <- x$X
	Y <- x$Y
	Z <- X + 1i*Y
	T <- x$T
	
	p.CI <- x$p.CI  
	
	
	# SOMETIMES - there are infinities in the time estimates (usually with no real range shift). We need to tidy those with the following ugly code. 
	
	ts <- which(grepl("t", row.names(p.CI)) & !(grepl("tau", row.names(p.CI))))
	if( max(abs(p.CI[ts,2:3])) == Inf) cat("Some confidence intervals around the timing parameters are infinitely wide.  I'll plot something fat but not quite THAT wide.  The range shift model might not be appropriate.\n")
	
	p.CI[ts,]$CI.high <- pmin(p.CI[ts,]$CI.high, p.CI[ts,]$p.hat*10)
	ts <-  which(row.names(p.CI) %in% c("t1","t2"))
	p.CI[ts,]$CI.low <- pmax(p.CI[ts,]$CI.low, min(x$T))

	# carrying on ...
	
	p.hat <- t(p.CI)[1,] %>% t %>% data.frame

	z.95 <- sqrt(-2*log(1-0.95))
	z.50 <- sqrt(-2*log(1-0.5))
	r95 <- sqrt(as.numeric(p.hat$A)/pi) 
	r50 <- r95 * (z.50/z.95)

	if(is.null(CI.cols)){
		blues <- brewer.pal(9, "Blues")[c(3,6,9)]
		g1 <- blues[1] #rgb(0.7,0.7,1)
		g2 <- blues[2] #g2 <- rgb(0.5,0.5,.8)
		g3 <- blues[3] #g3 <- rgb(0.05,0.05, 0.4)
	} else {g1 <- CI.cols[1]; g2 <- CI.cols[2]; g3 <- CI.cols[3]}

	which.mu <- which(!names(p.hat) %in% c("A","tau.z","tau.v"))
	p.mu.hat <- p.hat[which.mu]
	p.mu.se <- p.mu.hat*0 + with(p.CI[which.mu,], (CI.high - p.hat)/2)
	n.clust <- length(grep("x", names(p.mu.hat)))	

	if(is.null(layout) & plot.ts) layout(rbind(c(1,2), c(1,3))) 
	if(is.null(par)) par(mar = c(0,4,0,0), oma = c(4,0,4,4), xpd = NA) 
	
	# simulate and draw bars

	plot(X,Y, asp=1, type = "n", ...)
	if(n.clust == 2){
		Z.bars <- getBars(p.mu.hat, p.mu.se, n.sims = ns[1], n.time = ns[2], n.bins = ns[3], stretch = stretch)

		# draw big bar
		polygon.Z(Z.bars[,1], Z.bars[,4], col=g1)

		# draw bells
		with(p.mu.hat, {
				polygon.circle(x1, y1, r95, col=g1, border=NA)
				polygon.circle(x2, y2, r95, col=g1, border=NA)
				polygon.circle(x1, y1, r50, col=g2, border=NA)
				polygon.circle(x2, y2, r50, col=g2, border=NA)})
			
		# draw small bar
		polygon.Z(Z.bars[,2], Z.bars[,3], col=g2)
	}
	
	if(n.clust == 3){
		p.mu.hat1 <- with(p.mu.hat, data.frame(x1, x2, y1, y2, t1, dt=dt1))
		p.mu.se1 <- with(p.mu.se, data.frame(x1, x2, y1, y2, t1, dt=dt1))
		
		p.mu.hat2 <- with(p.mu.hat, data.frame(x1=x3, x2, y1=y3, y2, t1 = t2, dt = dt2))
		p.mu.se2  <- with(p.mu.se,  data.frame(x1=x3, x2, y1=y3, y2, t1 = t2, dt = dt2))
		
		Z.bars1 <- getBars(p.mu.hat1, p.mu.se1, n.sims = ns[1], n.time = ns[2], n.bins = ns[3])
		Z.bars2 <- getBars(p.mu.hat2, p.mu.se2, n.sims = ns[1], n.time = ns[2], n.bins = ns[3])

		# plot wide bars
		polygon.Z(Z.bars1[,1], Z.bars1[,4], col=g1)
		polygon.Z(Z.bars2[,1], Z.bars2[,4], col=g1)
		
		# plot big circles
		with(p.mu.hat, {
				polygon.circle(x1, y1, r95, col=g1, border=NA)
				polygon.circle(x2, y2, r95, col=g1, border=NA)
				polygon.circle(x3, y3, r95, col=g1, border=NA)}
				)
		
		# plot narrow bars
		polygon.Z(Z.bars1[,2], Z.bars1[,3], col=g2)
		polygon.Z(Z.bars2[,2], Z.bars2[,3], col=g2)
		
		# plot small circles
		with(p.mu.hat, {
				polygon.circle(x1, y1, r50, col=g2, border=NA)
				polygon.circle(x2, y2, r50, col=g2, border=NA)
				polygon.circle(x3, y3, r50, col=g2, border=NA)
				
				points(x1, y1, pch=19, cex=2, col=g3)
				points(x2, y2, pch=19, cex=2, col=g3)
				points(x3, y3, pch=19, cex=2, col=g3)
				})
	}
	# draw points
	points(X,Y, pch=21, bg = pt.col, type="o", cex=pt.cex)
	
	####################
	# Plot Time Series
	####################
	x$p.CI <- p.CI
	if(plot.ts)		plotFit.ts(x, CI.cols, pt.cex = pt.cex)
}


plotFit.ts <- function(FIT, pt.col = "antiquewhite", CI.cols = NULL, n.sims = 1e3, xlab = "Time", pt.cex = 1, ...){

  gr1 <- grey(.2)
 
	if(is.null(CI.cols)){
			blues <- brewer.pal(9, "Blues")[c(3,6,9)]
			g1 <- blues[1] #rgb(0.7,0.7,1)
			g2 <- blues[2] #g2 <- rgb(0.5,0.5,.8)
			g3 <- blues[3] #g3 <- rgb(0.05,0.05, 0.4)
		} else {g1 <- CI.cols[1]; g2 <- CI.cols[2]; g3 <- CI.cols[3]}

	X <- FIT$X
	Y <- FIT$Y
	T <- FIT$T
	
	p.CI <- FIT$p.CI
	which.mu <- which(!(row.names(p.CI) %in% c("A","tau.z","tau.v")))
	p.CI <- p.CI[which.mu,]
	
	p.hat <- t(p.CI)[1,]
	p.se <- p.hat*0 + with(p.CI, (p.hat-CI.low)/2)
	
	n.clust <- length(grep("x", names(p.hat)))	
	if(n.clust == 3)  GETMU <- getMu_multi else GETMU <- getMu
	XY.hat <-GETMU(T, p.hat)
	X.hat <- XY.hat[,1]
	Y.hat <- XY.hat[,2]
	
	p.sim <- rmvnorm(n.sims, mean = p.hat, sigma = diag(p.se^2))

	T.sim <- seq(min(T), max(T), length=500)
	XY.sim <- aaply(p.sim, 1, function(p) GETMU(T.sim, p))
	
	X.high <- apply(XY.sim[,,1], 2, quantile, p=c(0.75, 0.975))
	X.low  <- apply(XY.sim[,,1], 2, quantile, p=c(0.25, 0.025))
	Y.high <- apply(XY.sim[,,2], 2, quantile, p=c(0.75, 0.975))
	Y.low  <- apply(XY.sim[,,2], 2, quantile, p=c(0.25, 0.025))

	plot(T ,X, type="l", xaxt="n", xlab="")
	polygon.CI(T.sim, X.low[2,], X.high[2,], col=g1)
	polygon.CI(T.sim, X.low[1,], X.high[1,], col=g2)
	lines(T, X.hat, col=g3, lwd=3)
	lines(T, X, pch=21, bg=pt.col, col=gr1, type="o", cex=pt.cex)
	
	plot(T, Y, type="l", xlab=xlab)
	polygon.CI(T.sim, Y.low[2,], Y.high[2,], col=g1)
	polygon.CI(T.sim, Y.low[1,], Y.high[1,], col=g2)
	lines(T, Y.hat, col=g3, lwd=3)
	lines(T, Y, pch=21, bg=pt.col, col=gr1, type="o", cex=pt.cex)
}

	
getBars <- function(p.mu.hat, p.mu.se, n.sims = 1e3, n.time = 1e2, n.bins = 25, stretch = 0){	

	p.mu.sim <- rmvnorm(n.sims, mean = as.matrix(p.mu.hat)[1,], 
														 sigma = diag(as.matrix(p.mu.se^2)[1,]))
	
	T.sim <- with(p.mu.hat, seq(t1, t1 + dt, length=n.time))
	X.hat <- getMu(T.sim, p.mu.hat)[,1]
  Y.hat <- getMu(T.sim, p.mu.hat)[,2]
	XY.sim <- aaply(p.mu.sim, 1, function(p) getMu(T.sim, p))
	Z.sim <- XY.sim[,,1] + 1i*XY.sim[,,2]
	theta <- with(p.mu.hat, Arg(x2 + 1i*y2 - (x1 + 1i*y1)))	
	
	M.rotate <- complex(argument = -theta, modulus = 1)
	Zr <- Z.sim * M.rotate
	
	Z.center <- with(p.mu.hat, c(x1+1i*y1, x2 + 1i*y2)) 
	Zr.center <- Z.center * M.rotate
	
	Xr <- Re(Zr)
	dX <- diff(range(Re(Zr.center))) / n.bins
	Xr.trim <- seq(min(Re(Zr.center)) - dX, max(Re(Zr.center)) + dX, dX)
	Xr.bins <- cut(Xr, Xr.trim)
	
	# rotate	
	Xr.envelope <- (Xr.trim[-1] + Xr.trim[-length(Xr.trim)])/2
	Yr.envelope <- tapply(Im(Zr), Xr.bins, quantile, prob = c(0.025,0.25, 0.75, 0.975)) %>% ldply %>% mutate(.id = NULL)
	Zr.envelope <- Xr.envelope + 1i*Yr.envelope
	
	# snip at limits of x1,y1,x2,y2

	Zr.envelope <- subset(Zr.envelope, Xr.envelope > min(Re(Zr.center)) - stretch &  Xr.envelope < max(Re(Zr.center)) + stretch)
	Zr.envelope / M.rotate
}	
	
# helper functions for plots

polygon.circle <- function(x, y, r, n=100, ...){
  z.circle <- complex(argument = seq(0, 2*pi, length=n), modulus = r) + x + 1i*y
  polygon(Re(z.circle), Im(z.circle), ...)
}
polygon.XY <- function(X1, X2, Y1, Y2, ...)
  polygon(c(X1, X2[length(X2):1]), c(Y1, Y2[length(Y2):1]), border=NA, ...)
polygon.Z <- function(Z1, Z2, ...)
  polygon(c(Re(Z1), Re(Z2[length(Z2):1])), c(Im(Z1), Im(Z2[length(Z2):1])), border=NA, ...)
polygon.CI <- function(T, low, high, ...)
  polygon(c(T, T[length(T):1]), c(low, high[length(high):1]), border=NA, ...)

