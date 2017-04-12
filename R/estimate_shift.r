#' Estimating range shifts
#' 
#' Estimation and helper functions for nls fit of migration model
#' 
#' @details This algorithm minimizes the square of the distance of the locations from a double-headed hockeystick curve, then estimates the times scale using the ARMA/AR models.  Confidence intervals are obtained by bootstrapping the data and reestimating. See example and vignette for implementation. 
#'
#' @param T time
#' @param n.clust the number of ranges to estimate.  Two is relatively easy and robust, and three works fairly will (with good initial guesses).  More can be prohibitively slow. 
#' @param p.m0 initial parameter guesses - a named vector with (e.g.) elements x1, x2, y1, y2, t1, dt.  It helps if this is close - the output of \code{\link{quickfit}} can be helpful, as can plotting the curve and using \code{\link{locator}}. If left as NULL, the function will make some guesses for you - starting with \code{quickfit}. 
#' @param dt0 initial guess for duration of migration
#' @param method one of `ar` or `like` (case insenstive), whether or not to use the AR equivalence method (faster, needs regular data - with some tolerance for gaps) or Likelihood method, which is slower but robust for irregular data.
#' @param CI whether or not to estimate confidence intervals
#' @param nboot number of bootstraps 
#' @param X x coordinate
#' @param Y y coordinate
#' @param area.direct passed as direct argument to getArea
#' @param model one of "MWN", "MOU" or "MOUF" (case insensitive).  By default, the algorithm selects the best one according to AIC using the \code{\link{selectModel}} function.
#'
#' @return a list with the following elements
#' \item{T,X,Y}{Longitude coordinate with NA at prediction times} 
#' \item{p.hat}{Point estimates of parameters}
#' \item{p.CI}{Data frame of parameter estimates with (approximate) confidence intervals.}
#' \item{model}{One of "wn", "ou" or "ouf" - the selected model for the residuals.}
#' \item{hessian}{The hessian of the mean parameters.}
#' 
#' @aliases estimate.mouf.nls getP.nls geMu.nls
#' @example ./demo/estimate_shift_example.r
#' @export
estimate_shift <- function(T, X, Y, n.clust = 2,  
                           p.m0 = NULL, dt0=min(5, diff(range(T))/20),  
                           method = c("ar","like")[1], CI=TRUE, 
                           nboot = 100, model = NULL, 
                           area.direct = NULL)
{
  method <- tolower(method)
  
  n <- length(T)
  hessian <- NULL
  use.quickfit <- FALSE
  
  if(is.null(p.m0)){ 
    p.m0 <- try(quickfit(T,X,Y, dt=dt0, n.clust = n.clust, plotme=FALSE))
    use.quickfit <- TRUE
  }
  
  if(inherits(p.m0, "try-error")){
    warning("Having a hard time guessing initial centroids using k-means clustering. Going to take a wild guess, but you should consider providing guesses.\n")
    if(n.clust == 2)
      p.m0 <- c(t1 = round((T[n]-T[1])*0.5+T[1]), dt = mean(diff(T))*2, x1 = X[1], y1 = Y[1], x2 = X[n], y2 = Y[n]) 
    if(n.clust == 3)
      p.m0 <- c(t1 = round((T[n]-T[1])*0.3+T[1]),
                t2 = round((T[n]-T[1])*0.6+T[1]), 
                dt1 = mean(diff(T))*2, 
                dt2 = mean(diff(T))*2, 
                x1 = X[1], y1 = Y[1], 
                x2 = X[n], y2 = Y[n])
    if(n.clust > 3) stop("Scratch that. You're going to have to provide some decent guesses if you really want to fit more than 3 ranges at once.")
  }
  
  ##################################################
  # obtain nls fit and residuals
  ##################################################
  
  # mean parameters
  
  if(n.clust == 2) FUN <- getMu else if(n.clust > 2) FUN <- getMu_multi
  p.mu.fit <- getP.mu(T, X, Y, p.m0, FUN = FUN, hessian = TRUE)
  p.mu.hat <- p.mu.fit$par
  
  dts <- grep("dt", names(p.mu.hat))
  if(any(p.mu.hat[dts] < 0)){
    p.mu.fit <- getP.mu(T, X, Y, p.m0, FUN = FUN, hessian = TRUE, bounds=TRUE)
    p.mu.hat <- p.mu.fit$par
  }
  
  XY.hat <- FUN(T, p.mu.hat)
  Z.res <- (X + 1i*Y) - (XY.hat[,1] + 1i*XY.hat[,2])
  
  rangemodel.given <- TRUE
  if(is.null(model)){
    rangemodel.given <- FALSE
    model <- selectModel(Z.res, T, method = method)
  }

  model <- tolower(model)
  if(grepl("m", model)) model <- substring(model, 2)
  tau.fit	 <- try(getTau(Z.res, T=T, model=model, method = method, CI=CI))
  
  if(inherits(tau.fit, "try-error") & model == "ouf"){
    warning("Couldn't fit OUF to the residuals. Switching to OU model.\n")
    model <- "ou"
    tau.fit <- try(getTau(Z.res, T=T, model=model, method = method, CI=CI))
  }
  
  if(inherits(tau.fit, "try-error") & model == "ou"){
    warning("Couldn't fit OU to the residuals. Switching to WN model.\n")
    model <- "wn"
    tau.fit <- try(getTau(Z.res, T=T, model=model, method = method, CI=CI))
  }
  
  if(method == "ar" & model == "ouf") model <- "ou"
  tau.hat <- tau.fit$tau.hat
  ll <- tau.fit$ll
  aic <- computeAIC(ll, model, n.clust)
  
  # compile results
  p.hat <- c(tau.hat, p.mu.hat)
  A <- getArea(p.hat,T,X,Y, model=model, direct = area.direct)
  
  mrsa.fit <- list(p.hat = c(A = A$A.hat, p.hat), 
                    model = model, 
                    ll=c(ll=ll), aic = aic['aic'], df = aic['df'],
                    method = method,
                    n.clust = n.clust,
                    rangemodel.given = rangemodel.given,
                    p.m0 = p.m0, 
                    use.quickfit = use.quickfit,
                    X=X, Y=Y, T=T, Z.res=Z.res)	
  
  if(CI) {
    # obtain mean parameters
    hessian <- p.mu.fit$hessian
    
    # hideous, hideous kludge to fix occasinal, inexplicable, enormous CI's for time related estimates
    if(n.clust == 2) 
      p.mu.sd <- try(sqrt(1/diag(hessian))) else {	
			#try(sqrt(diag(solve(hessian)))) else {	
        which.xy <-	c(grep("x", names(p.mu.hat)), grep("y", names(p.mu.hat)))
        which.t <-	grep("t", names(p.mu.hat))
        p.mu.sd <- try(sqrt( c(diag(solve(hessian[which.xy, which.xy])), 1/diag(hessian[which.t,which.t]))))
      }
    
    if(inherits(p.mu.sd, "try-error")){
      warning(cat("Can't obtain some confidence intervals (probably the migration time is too fast).\n"))
      p.mu.sd <- p.mu.hat * NA
      zeros <- which(diag(hessian) == 0)		
      p.mu.sd[-zeros] <- sqrt(diag(solve(hessian[-zeros,-zeros])))
    }
    if(any(is.na(p.mu.sd))){
      warning(cat("Can't obtain some confidence intervals - most likely because the migration time is too fast.\n"))
      p.mu.sd[is.na(p.mu.sd)] <- 0
    }
    
    # obtain variance parameters
		
		tau.bs <- matrix(ncol = length(tau.hat), nrow = nboot)
		
    if(model == "wn") tau.withCI <- NULL
    if(model == "ou" | model == "ouf"){
      if(method == "ar"){
			
				count.errors <- 0

        cat("Using bootstrap for confidence intervals around the variance parameters.\n")
				for(i in 1:nboot){
				sample <- sample(1:n, replace=TRUE) %>% sort %>% unique
        tau.fit <- try(getTau(Z.res = Z.res[sample], T = T[sample], 
															 tau0 = tau.hat, model = model, method="ar"), silent = TRUE)
															 
				 if(inherits(tau.fit, "try-error")){
							tau.bs[i,] <- rep(0, length(tau.hat))
							count.errors <- count.errors + 1
							}
					else tau.bs[i,] <- tau.fit$tau.hat
				}
								
        tau.withCI <- apply(t(t(tau.bs)), 2, quantile, prob = c(0.5,0.025, 0.975), na.rm=TRUE) %>% data.frame
        names(tau.withCI) <- names(tau.hat)
        tau.withCI <- t(tau.withCI)
				
				if(count.errors >= 10)
					cat("Several bootstrapped samples failed to estimate time-scale parameters.  Consider trying a simpler model.\n")
      } 
      if(method == "like"){ 
        tau.CI <- tau.fit$tau.CI
        if(model == "ou") {tau.withCI <- data.frame(tau.z = c(tau.hat, tau.CI$tau.z)) %>% t} else
          tau.withCI <- cbind(tau.hat, t(tau.CI))
      }	
    }
    
    p.mu.CI <-  cbind(p.mu.hat, p.mu.hat - 1.96*p.mu.sd, p.mu.hat + 1.96*p.mu.sd)
    
    p.CI <- rbind(A = c(A$A.hat, A$A.hat + c(-1.96,1.96) * A$A.se), 
                  tau.withCI, 
                  p.mu.CI) %>% data.frame 
    names(p.CI) <- c("p.hat", "CI.low", "CI.high")
    
    # truncate negative dt's
    dts <- grep("dt", row.names(p.CI))
    p.CI[dts,]$CI.low <- pmax(0,p.CI[dts,]$CI.low)
    mrsa.fit$p.CI <- p.CI	
  }
  class(mrsa.fit) <- "shiftfit"
  return(mrsa.fit)
}

#' @export
print.shiftfit <- function(x, ...){
  # x is the fit
  with(x, {
    cat("Range shift estimation:\n\n")
    n <- length(X)
       cat(paste0(n, " observations over time range ", min(T), " to ", max(T),".\n"))
    if(n.clust == 2)
      cat(paste0(n.clust, " ranges (one shift) fitted using the ", toupper(method), " method.\n\n"))  else
        cat(paste0(n.clust, " ranges (", n.clust-1, " shifts) fitted using the ", toupper(method), " method.\n")) 
   if(use.quickfit){ 
     cat("Initial shift parameters obtained from `quickfit`:\n") 
     print(signif(p.m0, 4),...)}
       
   if(!rangemodel.given) 
     cat(paste0("The AIC-selected ranging model is ", toupper(model), ".\n\n")) else
       cat(paste0("The fitted ranging model is ", toupper(model), ".\n\n"))
    
    cat("Parameter Estimates:\n")
    if(!("p.CI" %in% names(x))) print(t(p.hat),...) else print(p.CI,...)
    cat(paste0("\nlog-likelihood: ", signif(ll,4), "; AIC: ", signif(aic,4), "; degrees of freedom: ", signif(df,4), "\n"))
  })
  invisible(x)
}
#' @export
summary.shiftfit <- function(object,...){print.shiftfit(object,...); invisible(object)}
#' @export
AIC.shiftfit <- function(object,...) object$aic
#' @export
logLik.shiftfit <- function(object,...) object$ll


computeAIC <- function(ll, model, n.clust){
  if(model %in% c("mouf","ouf")) k.s <- 3 else
    if(model %in% c("mou","ou")) k.s <- 2 else 
      if(model %in% c("mwn", "wn")) k.s <- 1
      k.m <-  n.clust * 2 + (n.clust - 1) * 2
      
      return(c(aic = -2*ll + 2*(k.s + k.m), df = k.s + k.m))
}

getP.mu <- function(T,X,Y,p.m0, FUN = getMu, hessian = TRUE, bounds = TRUE){
  getL.nls <- function(p.m, T, X, Y, FUN){
    Mu <- FUN(T, p.m)
    sum(c(X - Mu[,1], Y - Mu[,2])^2)
  }
  
  xs <- grep("x", names(p.m0))
  ys <- grep("y", names(p.m0))
  dts <- grep("dt", names(p.m0))
  ts <- (1:length(p.m0))[-c(xs,ys,dts)]
  
  if(!bounds) fit <- try(optim(p.m0, getL.nls, T=T, X=X, Y=Y, hessian=hessian, FUN = FUN)) else {
    upper <- p.m0*0
    lower <- p.m0*0
    
    lower[xs] <- min(X)
    lower[ys] <- min(Y)
    lower[ts] <- min(T)
    lower[dts] <- 0
    upper[xs] <- max(X)
    upper[ys] <- max(Y)
    upper[ts] <- max(T)
    upper[dts] <- diff(range(T))
    fit <- try(optim(p.m0, getL.nls, T=T, X=X, Y=Y, hessian=hessian, FUN = FUN, method="L-BFGS-B", lower = lower, upper = upper))
  }
  fit
}
