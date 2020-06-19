#' Test range shift using net-squared displacement
#' 
#' @details The test below assumes that the net squared displacement (NSD) for a migrating organism is well characterized 
#' by the logistic formula: E(NSD(t)) =  a / (1 + exp [(b-t)/c]) as described in Boerger and Fryxell (2012).  In practice, 
#' the square root of the NSD, i.e., the linear displacement, is fitted to the square root of the formula assuming Gaussian 
#' residuals with constant variance 's'.  A likelihood ratio test against a null model of no-dispersal is provided at a 95\% 
#' significance level.
#' @param T time
#' @param plotme whether or not to plot the result
#' @param setpar whether or not to run par(mfrow = c(1,2)) before plotting
#' @param X x coordinate
#' @param Y y coordinate
#' @param ... additional parameters to pass to \code{plot}
#' @encoding UTF-8
#' @references Boerger, L. & Fryxell, J. (2012) Quantifying individual differences in dispersal using net squared displacement. 
#' Dispersal Ecology and Evolution (eds. J. Clobert, M. Baguette, T. Benton & J. Bullock), pp. 222-228. Oxford University 
#' Press, Oxford, UK.
#'
#' @return a list with a vector of four parameter estimates, and a vector with test statistics (likelihood, AIC and p.values)
#' @example ./demo/fitNSD_example.r
#' @export
#' 
fitNSD <- function(T,X,Y,plotme=FALSE,setpar=TRUE, ...){
  
  Z <- X + 1i*Y
  NSD <- Mod((Z-Z[1])^2)
  f.NSD <- function(t,a,b,c) a / (1 + exp((b - t)/c))
  getLL <- function(p, x, y){
    y.hat <- f.NSD(x, p["a"], p["b"], p["c"])
    -sum(dnorm(sqrt(y), sqrt(y.hat), exp(p["s"]), log=TRUE))	
  }
  y <- NSD
  x <- T
  p0 <- c(a = max(NSD), b = mean(x), c = 1, s = sd(sqrt(y)))
  NSD.fit <- suppressWarnings(optim(p0, getLL, x = x, y = y))
  p.hat <- NSD.fit$par
  LL <- -NSD.fit$value
  AIC <- -2*LL + 8
  
  # null model
  mu.null <- mean(sqrt(y))
  sigma.null <- sd(sqrt(y))
  LL.null <- sum(dnorm(sqrt(y), mu.null, sigma.null, log=TRUE))	
  AIC.null <- -2*LL.null + 4
  
  lrt <- 2*(LL - LL.null)
  p.value <- 1 - pchisq(lrt, 2)
  
  if(plotme){
    y.hat <- f.NSD(x, p.hat["a"], p.hat["b"], p.hat["c"])
    residual <- sqrt(y) - sqrt(y.hat)
    if(setpar) par(mfrow=c(1,2))
    plot(x, y, type="o", ...)
    lines(x, y.hat, col=2, lwd=2)
    qqnorm(residual)
  }	
  
  list(estimate = p0, 
       test = c(LL=LL, AIC=AIC, LL.null=LL.null, AIC.null=AIC.null, p.value=p.value))
}
