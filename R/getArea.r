#' Compute area 
#'
#' Compute predicted area at given alpha level (e.g. 50\% or 90\%) of a migration model fit
#' 
#' @details For sufficient data (i.e. where the range in the times is much greater than the ) This function estimates the (symmetric) 95\% area of use from a bivariate Gaussian 
#' 
#' @param p estimated mouf parameter vector (tau.z, tau.v, t1, dt, x1, y1, x2, y2)
#' @param {T,X,Y} time and location data
#' @param alpha proportion of area used to be computed 
#' @param model one of "wn", "ou", "ouf" - whether or not the velocity autocorrelation needs to be taken into account. 
#' @param n.clust number of ranges
#' @param method one of "nls" or "like"
#' @param direct whether or not to compute the area directly (i.e. fitting a symmetric bivariate normal to the residuals) or to account for the autocorrelation.  The default behavior (NULL) computes directly for the "wn" model, and uses the autocorrelation (which is slower) only if the estmated spatial time scale is greater that 1/30 of the total time range. 

getArea <- function(p, T, X, Y, alpha=0.95, model=c("wn", "ou", "ouf")[1], direct = NULL) {
  directA <- function(Z.res){
    s2 <- (var(Re(Z.res)) + var(Im(Z.res)))/2
    z.p <- sqrt(-2 * log(0.05))
    A.hat <- z.p^2 * pi * s2
    A.se <- A.hat/sqrt(length(T))
    list(A.hat=A.hat, A.se=A.se)
  }
  
  timescaleA <- function(Z.res, model){
    n <- length(T)
    S <-  outer(T, T, getCov, model=model, p=tau)
    S.inv <- Matrix::solve(S)
    
    Xt <- Re(Z.res)
    Yt <- Im(Z.res)
    
    s.xx <- as.numeric(t(Xt) %*% S.inv %*% Xt)/n
    s.yy <- as.numeric(t(Yt) %*% S.inv %*% Yt)/n
    
    s.hat <- (s.xx + s.yy)/2
    A.hat <- -2 * pi * s.hat * log(1-alpha)
    A.se <- A.hat/sqrt(n)
    list(A.hat=A.hat, A.se=A.se)
  }
  
  tau <- p[c("tau.z","tau.v")]
	p.m <- p[which(! names(p) %in% c("tau.z","tau.v"))]

	n.clust <- length(grep("x", names(p.m)))
	
	if(n.clust == 2) XY.hat <- getMu(T, p.m) else XY.hat <- getMu.multi(T, p.m) 
  Z.hat <- XY.hat[,1] + 1i*XY.hat[,2]
  Z <- X + 1i*Y 
  Z.res <- Z  - Z.hat
	
  if(model == "wn") direct <- TRUE else 
    if(is.null(direct)){
      if(diff(range(T))/tau["tau.z"] > 30){
        cat("Calculating area directly (time range much greater than time scale).\n")
        direct <- TRUE
      } else{ 
        cat("Calculating area from time scale constants due to insufficient time range.\n")
        direct <- FALSE}}
  
  if(direct) return(directA(Z.res)) else 
    return(timescaleA(Z.res, model)) 
  }