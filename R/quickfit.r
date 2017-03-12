#' Quick fit of one-step migration
#' 
#' Using k-means clustering to get quick fits of 2 or 3 cluster centers in X-Y coordinates.  
#' 
#' @details This function does estimates the locations and times of migration, but not the duration (dt).  It is most useful for obtaining a "null" estimate for seeding the likelihood estimation. 
#' 
#' @param T time
#' @param X x coordinate of movement
#' @param Y y coordinate of movement
#' @param dt duration of migration (arbitrarily = 1)
#' @param n.clust number of clusters (2 or 3)
#' @param plotme whether or not to plot the result
#' @return a named vector of initial estimates: 
#' \itemize{
#' \item{if \code{n.clust = 2}} {returns \code{t1, dt, x1, y1, x2, y2}} 
#' \item{if \code{n.clust = 3}} {returns \code{t1, dt1, t2, dt2, x1, y1, x2, y2, x3, y3} }}
#' @example ./demo/quickfit_example.r
#' @export
quickfit <-  function(T,X,Y, dt=1, n.clust = 2, plotme=TRUE){
  clusterme <- cbind(X, Y)
  r.cl <- kmeans(clusterme, n.clust)
  
  if(n.clust == 2){
    if(lm(r.cl$cluster~T)$coef[2] < 0){
      r.cl$cluster <- 3 - r.cl$cluster
      r.cl$centers <- r.cl$centers[2:1,]
    }
    x0y0 <- r.cl$center[rbind(c(1,1),c(1,2),c(2,1),c(2,2))]
    names(x0y0) <- c("x1", "y1", "x2", "y2")
    clustering <- zoo(r.cl$cluster, order.by = T)
    t0 <- T[which(diff(smooth(clustering))!=0)][1]
    p0.hat <- c(t1 = max(0, t0-dt/2), dt = dt, x0y0)
  }
  
  if(n.clust == 3){
    n <- length(X)	
    z.centers <- r.cl$center[,1] + 1i*r.cl$center[,2]
    
    first <- names(which.min(Mod(z.centers - (X[1] + 1i * Y[1]))))
    last <- names(which.min(Mod(z.centers[-as.numeric(first)] - (X[n] + 1i * Y[n]))))
    middle <- which(!(1:3 %in% c(first,last)))
    
    c.raw <- r.cl$cluster
    c.new <- rep(2, length(c.raw))
    c.new[c.raw == first] <- 1
    c.new[c.raw == last] <- 3
    
    xs <- Re(z.centers)[c(first,middle,last)]
    ys <- Im(z.centers)[c(first,middle,last)]
    
    c.smooth <- zoo(c.new, order.by = T) %>% smooth %>% as.numeric
    t.starts.all <- T[which(diff(c.smooth) == 1)]
    t.starts <- t.starts.all[c(1,length(t.starts.all))]
    
    names(xs) <- paste0("x",1:3)
    names(ys) <- paste0("y",1:3)
    names(t.starts) <- paste0("t",1:2)
    
    p0.hat <- c(xs, ys, t.starts, dt1 = dt, dt2 = dt)
  }
	
  if(plotme){
    plot(clusterme, type="l",  col="grey", asp=1)
    points(clusterme,col=r.cl$cluster, asp=1)
  }
	
  return(p0.hat)
}