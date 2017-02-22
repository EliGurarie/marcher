#' Compute time scale parameters
#' 
#' A mostly internal function that takes the "residuals" of a range-shift process and estimates \deqn{\tau_z} and, if necessary, \deqn{\tau_v}.
#' 
#' @param Z.res complex vector of isotropic Gaussian, possibly autocorrelated time series of points
#' @param T time vector
#' @param model one of  \code{"wn"} (white noise),  \code{"ou"} or \code{"ouf"} (case insensitive), denoting, respectively, no autocorrelation, position autocorrelation, or velocity and position autocorrelation.  If \code{model = NULL} and \code{method = "ar"}, the algorithm will select a model using AIC comparisons of the three. If the selected model iswhite noise, the function will return 0's for both parameters.
#' @param tau0 initial values of parameter estimates - a named vector: \code{c(tau.z = tau0[1], tau.v = tau0[2])}
#' @param method either \code{"like"} or \code{"ar"}.  The former refers to the likelihood method - it is most general (i.e. works with irregular sampling). The latter refers to the auto-rgressive model equivalence, which is faster but only works with regular sampling. 
#' @param CI whether or not to compute the confidence intervals (temporarily only available for \code{like} method).


getTau <- function(Z.res, T = T, model=c("wn", "ou", "ouf")[1], tau0 = NULL, CI = FALSE, method = c("like", "ar")[1]){
  Z.res <- Z.res-mean(Z.res)
  
  model <- tolower(model)
  method <- tolower(method)
  
  if(model == "wn")  tau.fit <- getTau.WN(Z.res, T) else { 
    if(method == "like") tau.fit <- getTau.OUF.like(Z.res, T, CI=CI, tau0=tau0, model = model)
    if(method == "ar") tau.fit <- getTau.OU.ar(Z.res, T, CI=CI, tau0=tau0, model = model)
  }
  return(tau.fit)	
}


getTau.WN <- function(Z.res, T = T){
  tau.hat = NULL
  tau.CI = NULL
  ll <- getLikelihood.res(p.s=0, T, Z.res, "wn")
  return( list(tau.hat = tau.hat, model = "wn", ll=ll, tau.CI = tau.CI))
}

#########################
# Likelihood Estimation
#########################
 
# logtau.zs <- seq(-2,1,length=30)
# kappas <- seq(1e-13,0.8,length=30)
# 
# LL <- outer(logtau.zs, kappas, 
# 					Vectorize(function(x,y) getLikelihood.res(c(logtau.z = x, kappa = y), 
#  T=T, Z.res = Z.res, model="ouf")))
#   
#  require(rgl)
#  z.len <- round(max(LL) - min(LL))+1
#  cols <- terrain.colors(z.len)[round(LL - min(LL, na.rm=TRUE))]
#  terrain3d(logtau.zs, kappas, LL/100, color = cols)
#  
#  a.max <- which(LL == max(LL, na.rm=TRUE), arr.ind = TRUE)
#  logtau.zs[a.max[1]]
#  kappas[a.max[2]]
#  axes3d()
 

getTau.OUF.like <- function(Z.res, T = T, CI = FALSE, tau0 = NULL, model = model){
  
  if(is.null(tau0)) 
    tau0 <- c(tau.z = diff(range(T))/10, tau.v = diff(range(T))/100)
  
  logtau0 <- c(logtau.z = as.numeric(log(tau0["tau.z"])))
  logtau0.lower <- c(logtau.z = log(mean(diff(T))/10))
  logtau0.upper <- c(logtau.z = log(diff(range(T))))
  
  if(model == "ouf"){ 
    logtau0["kappa"] <- tau0["tau.v"]/tau0["tau.z"]
    logtau0.lower["kappa"] <- 1e-13
    logtau0.upper["kappa"] <- 0.5
  }
  waserror <- FALSE
  FIT <- try(optim(logtau0, getLikelihood.res, T=T, Z.res = Z.res, 
                   model = model, 
                   method="L-BFGS-B",
                   lower = logtau0.lower, 
                   upper = logtau0.upper,
                   control = list(fnscale = -1), 
                   hessian = TRUE))	
  if(inherits(FIT, "try-error")){
    waserror <- TRUE
    cat("Trying to fit the model without a Hessian... \n")
    FIT <- try(optim(logtau0, getLikelihood.res, T=T, Z.res = Z.res, 
                     model = model, 
                     method="L-BFGS-B",
                     # method="Nelder-Mead",
                     lower = logtau0.lower, 
                     upper = logtau0.upper,
                     control = list(fnscale = -1), 
                     hessian = FALSE))
  }
  tau.z.hat <- as.numeric(exp(FIT$par["logtau.z"]))
  tau.hat <- c(tau.z = tau.z.hat)
  if(model == "ouf") tau.hat["tau.v"] <- FIT$par["kappa"] * tau.hat["tau.z"]
  
  if(CI & !waserror){
    if(model == "ou") logtau.se <- diag(1/sqrt(-FIT$hessian) )
    if(model == "ouf") logtau.se <- sqrt(diag(solve(-FIT$hessian)))
  } else logtau.se <- NA
  
  tau.CI <- data.frame(tau.z = with(FIT, exp(par['logtau.z'] + c(low = -2, high = 2)*logtau.se['logtau.z'])))
  if(model == "ouf") 
    tau.CI <- data.frame(tau.CI, tau.v = tau.hat["tau.z"] * pmin(1, pmax(0, FIT$par["kappa"] + c(-2,2)*logtau.se["kappa"])))
  ll <- FIT$value
  
  return( list(tau.hat = tau.hat, model = "ou", ll=ll, tau.CI = tau.CI))
}


#########################
# AR Estimation
#########################

getTau.OU.ar <- function(Z.res, T = T, CI = FALSE, tau0 = NULL, model = model){
  method <- "ar"
  if(is.null(model)){
    model <- selectModel(Z.res)
    cat(paste0("The selected model based on AIC comparisons is ", toupper(model), ".\n"))
  }
  if(model == "ouf"){ 
    model <- "ou"
    warning("AR approximation for OUF model not yet incorporated. I will fit an OU model for now, but the likelihood method is a better bet.\n")
  }
  
  fit.x <- arima(Re(Z.res), order=c(1,0,0), include.mean = FALSE, method="ML")
  fit.y <- arima(Im(Z.res), order=c(1,0,0), include.mean = FALSE, method="ML")
  phi.hat <- pmax((fit.x$coef + fit.y$coef)/2,0) %>% as.numeric
  tau.hat <- c(tau.z = -1/log(phi.hat))
  
  ll <- getLikelihood.res(p.s = c(logtau.z = as.numeric(log(tau.hat))), 
                          T, Z.res, model = "ou")
  
  return(list(tau.hat = tau.hat, model = model, ll=ll))
}

## AR and ARMA equations for Area

# s2.ar <- (fit1$sigma2 + fit2$sigma2)/2
# s2.f <- 2 * s2.ar / (tau.z.hat * (1 - phi.hat)^2)
#A.hat <- -2*pi*s2.f*log(0.05)  

# s2.arma <- (fit1$sigma2 + fit2$sigma2)/2
# s2.f <- 2*( (1 + theta.hat)^2 / (1 - phi.hat)^2) * s2.arma * ((tau.z.hat + tau.v.hat)/(tau.v.hat^2 * tau.z.hat^2))	  

#    if(model == "ouf"){
#        fit.x <- arima(Re(Z.res), order=c(2,0,0), include.mean = FALSE)
#        fit.y <- arima(Im(Z.res), order=c(2,0,0), include.mean = FALSE)  
#        phi.hat <- (fit.x$coef["ar1"] + fit.y$coef["ar1"])/2
#        theta.hat <- (fit.x$coef["ar2"] + fit.y$coef["ar2"])/2
#        tau.z.hat <- -1/log(phi.hat)
#        tau.v.hat <- -1/log(theta.hat)
#      } 

