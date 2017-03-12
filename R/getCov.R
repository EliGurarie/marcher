#' Estimation Helper Functions
#' 
#' functions which provide the theoretical covariance [getCov()] and area [getArea()] for specific models and parameter values
#' 
#' @param t1 time 1
#' @param t2  time 2
#' @param model the model
#' @param p vector of the parameters
#'
#' @details \code{getCov(t1, t2, model, p)} calculates the covariance matrix for different models. \code{mvrnorm2} is a slightly more efficient multivariate normal function. 
#' 
#' 
#' @aliases getCov mvrnorm2
#' @export
getCov <- function(t1, t2, model, p){
  dt <- abs(t1-t2)
  if(model == "wn")  # p = variance
    return(as.numeric(dt == 0)) else
  if(model %in% c("ou","mou")) # p = tau
    return(exp(-dt/p["tau.z"])) else
  if(model %in% c("ouf", "mouf")) # p = c(tau.z, tau.u)
    return((p["tau.v"]*exp(-dt/p["tau.v"])-p["tau.z"]*exp(-dt/p["tau.z"]))/(p["tau.v"]-p["tau.z"])) else
  stop("No such model!")
}

#' @export
mvrnorm2 <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE) 
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  eS <- eigen(Sigma, symmetric = TRUE, EISPACK = FALSE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% 
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1) 
    drop(X)
  else t(X)
}
