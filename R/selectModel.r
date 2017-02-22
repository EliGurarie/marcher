#' Select residual model
#' 
#' Given a complex vector of movement residuals, will use AIC to select the order of the autocorrelation, i.e. white noise (WN), position autocorrelation (OU), or position and velocity autocorrelation (OUF)
#' 
#' @param Z.res complex vector of residuals
#' @param T time vector (only needed for method = 'like')
#' @param method One of 'ar' or 'like' - whether to use the AR equivalence (faster, but needs to be regular) or likelihood estimation.
#' @param showtable whether to return the AIC values of the respective models
#' @return A character string - \code{'wn'}, \code{'ou'} or \code{'ouf'}.  Optionally also the AIC table. 
#' @example ./examples/selectModel.example.r

selectModel <- function(Z.res, T=NULL, method = c("ar", "like")[1], showtable = FALSE){
  method <- tolower(method)
  if(!(method %in% c("ar", "like")))
    stop("Please select one of two possible methods: 'ar' and 'like'.")
  
  if(method == "ar"){
    ar1.fit.x <- try(arima(Re(Z.res), order=c(1,0,0), include.mean = FALSE, method="ML"))
    ar1.fit.y <- try(arima(Im(Z.res), order=c(1,0,0), include.mean = FALSE, method="ML"))
    ar2.fit.x <- try(arima(Re(Z.res), order=c(2,0,0), include.mean = FALSE, method="ML"))
    ar2.fit.y <- try(arima(Im(Z.res), order=c(2,0,0), include.mean = FALSE, method="ML"))
		aic.wn <- -2*(sum(dnorm(Re(Z.res), 0, sd = sd(Re(Z.res)), log=TRUE)) + sum(dnorm(Im(Z.res), 0, sd = sd(Im(Z.res)), log=TRUE))) + 4
    
		if(inherits(ar1.fit.x, "try-error") | inherits(ar1.fit.y, "try-error")){
			warning("Can't estimate AR1 model for residuals.\n")
			aic.ou <- NA
		} else aic.ou <- AIC(ar1.fit.x) + AIC(ar1.fit.y)
		
		if(inherits(ar2.fit.x, "try-error") | inherits(ar2.fit.y, "try-error")){
			warning("Can't estimate AR2 model for residuals.\n")
			aic.ouf <- NA
		} else aic.ouf <- AIC(ar2.fit.x) + AIC(ar2.fit.y)
    aic <- c(aic.wn, aic.ou, aic.ouf)
  }
  
  if(tolower(method) == "like"){
    ll.wn <- getLikelihood.res(T = T, Z.res = Z.res, model = "wn")	
    
    # ou
    p.s <- c(logtau.z = log(diff(range(T))/10))
    p.s.lower <- c(logtau.z = -1e5)
    p.s.upper <- c(logtau.z = log(diff(range(T))/2))
    fit.ou <- try(optim(p.s, getLikelihood.res, T = T, Z = Z.res, model = "ou", method="L-BFGS-B",
                        lower = p.s.lower, upper = p.s.upper, control = list(fnscale = -1)))
    ll.ou <- ifelse(!inherits(fit.ou, "try-error"), fit.ou$value, NA)
    
    #ouf
    if(!inherits(fit.ou, "try-error")) 
      p.s <- c(fit.ou$par, kappa = 0.1) else p.s <- c(p.s, kappa = 0.1)
    p.s.lower['kappa'] <- 1e-5
    p.s.upper['kappa'] <- 0.8
    fit.ouf <- try(optim(p.s, getLikelihood.res, T = T, Z = Z.res, model = "ouf", method="L-BFGS-B",
                         lower = p.s.lower, upper = p.s.upper, control = list(fnscale = -1)))
    ll.ouf <- ifelse(!inherits(fit.ouf, "try-error"), fit.ouf$value, NA)
    aic <- -2*c(ll.wn, ll.ou, ll.ouf) + 2*c(2,3,4)
  }
  names(aic) <- c("wn", "ou", "ouf")
  model <- c("wn","ou","ouf")[which.min(aic)]
  if(!showtable) return(model) else
  return(list(model = model, aic = aic))
}

