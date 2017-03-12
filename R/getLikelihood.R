#' Estimate likelihoods and AICs
#' 
#' Estimate likelihoods and AIC for several possible migration models. 
#' 
#' @param p initial parameters: [tau.z, tau.v, t1, t2, x1, x2, y1, y2]
#' @param T,X,Y time,x and y coordinates
#' @param model "wn", "ou", "ouf", "mou" or "mouf",  - whether or not to estimate tau.v
#' @aliases getLikelihood.res getAIC.nls
#' @export
getLikelihood <- function(p, T, X, Y, model = c('mouf', 'mou', 'mwn', 'ouf', 'ou', 'wn')[1])
{
  # parameters of MOUF are: [log(tau.z), k.tau, t1, dt, x1, x2, y1, y2]
  
  if(model %in% c("mouf", "ouf"))
    p.s <- c(tau.z = as.numeric(exp(p['logtau.z'])), 
             tau.v = as.numeric(p['kappa']*exp(p['logtau.z']))) else 
               if(model %in% c("mou", "ou"))
                 p.s <- c(tau.z = as.numeric(exp(p['logtau.z']) )) else 
                   p.s <- NULL

  if(model %in% c("mouf", "mou", 'mwn')){
    p.m <- p[c("t1","dt","x1","y1","x2","y2")]
    M <- getMu(T, p.m=p.m)
    Xt <- X - M[,1]
    Yt <- Y - M[,2]
  } else {
    Xt <- X - mean(X)
    Yt <- Y - mean(Y)
  }
  
  s.model <- gsub("m", "", model)
	
  n <- length(T)
  S <-  outer(T, T, getCov, model=s.model, p=p.s)
  S.inv <- Matrix::solve(S)
  logS.det <- as.numeric(Matrix::determinant(S, logarithm=TRUE)$mod)
  
  s.xx <- as.numeric(t(Xt) %*% S.inv %*% Xt)/n
  s.yy <- as.numeric(t(Yt) %*% S.inv %*% Yt)/n
  
  s.hat <- (s.xx + s.yy)/2
  -logS.det - n*log(s.hat) - n*log(2*pi) - n
}


getLikelihood.res <- function(p.s=NULL, T, Z.res, model = c("wn", "ou", "ouf")[1]){

	if(model == "wn") p <- 0 else
	if(model == "ou" | model == "ouf"){
		p <- c(tau.z = as.numeric(exp(p.s['logtau.z'])), 
           tau.v = as.numeric(p.s['kappa']*exp(p.s['logtau.z'])))
					 }
			
    n <- length(T)
    S <-  outer(T, T, getCov, model=model, p=p)
    S[is.na(S)] <- 0
    S.inv <- Matrix::solve(S)
    Xt <- Re(Z.res)
    Yt <- Im(Z.res)
    logS.det <- as.numeric(Matrix::determinant(S, logarithm=TRUE)$mod)
    s.xx <- as.numeric(t(Xt) %*% S.inv %*% Xt)/n
    s.yy <- as.numeric(t(Yt) %*% S.inv %*% Yt)/n
    s.hat <- (s.xx + s.yy)/2
		ll <- -logS.det - n*log(s.hat) - n*log(2*pi) - n
#		else {
#    residuals <- c(Re(Z.res), Im(Z.res))
#    sigma.hat <- sd(residuals)		
#    ll <- sum(dnorm(residuals, 0, sigma.hat, log=TRUE))
#  }
  return(ll)
}
				
getAIC <- function(FIT, model = NULL, compute.null = TRUE)
{
  if(is.null(model))  model <- FIT$model
  p.hat <- t(FIT$p.CI$p.hat)
  names(p.hat) <- rownames(FIT$p.CI)
  
  p.s <- p.hat[c("tau.z", "tau.v")]
  p.m <- p.hat[c("t1", "dt", "x1", "y1", "x2", "y2")]
  Z.hat <- getMu(FIT$T, p.m)
  Z <- with(FIT, X + 1i*Y)
  Z.res <- Z - Z.hat
  
  ll <- getLikelihood.res(p.s = p.s, T=FIT$T, Z.res=Z.res, model=model)
  
  k <- 7  
  if(model == "ou") k <- 8 else if(model == "ouf") k <- 9
  aic <- -2*ll + 2*k
  aic.table <- data.frame(ll, k, aic)
  
  if(compute.null)
  {
    Z.res <- (FIT$X + 1i*FIT$Y) - (mean(FIT$X) + 1i*mean(FIT$Y))
    p.s <- p.hat[c("tau.z", "tau.v")]
    ll.null <- getLikelihood.res(p.s = p.s, T=FIT$T, Z.res=Z.res, model=model)
    k.null <- k-4
    aic.null <-  -2*ll.null + 2*k.null
    aic.table <- data.frame(ll = c(ll, ll.null), k = c(k, k.null), aic = c(aic, aic.null))
    row.names(aic.table) <- c("with shift", "without shift")
    if(aic.table$aic[1] < aic.table$aic[2])
    cat("According to AIC, a range shift is likely.\n") else cat("According to AIC, a range shift is less likely.\n") 
  }
  return(aic.table)
}
