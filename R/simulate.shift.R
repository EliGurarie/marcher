#' Simulate MOUF process
#' 
#' Simulate MOUF process
#'
#'@param T sampling times 
#'@param tau variance parameters - named vector with 'tau.z' and 'tau.v'
#'@param A 95\% area parameter
#'@param mu mean vector - typically output of \code{\link{getMu}}. Can also be any complex or a two-column matrix, or a multi-column matrix with some named columns "x" and "y" (case-insensitive)).
#'@return a data frame with Time, X, and Y columns.
#'@example ./examples/simulate.shift.example.r
#'@seealso \code{\link{getMu}}

simulate.shift <- function (T, tau = NULL, mu, A) {
  if(is.complex(mu)) 
    mu <- cbind(x = Re(mu), y = Im(mu))
  if(ncol(mu) == 2 & !all(c("x","y") %in% colnames(mu))) 
    colnames(mu) <- c("x","y")
  
  # lower case "X" and "Y"
  colnames(mu)[c("X","Y") %in% colnames(mu)] <- 
    tolower(colnames(mu)[c("X","Y") %in% colnames(mu)])
  
  if(!all(c("x","y") %in% colnames(mu)))
    if(ncol(mu)==2){  
        colnames(mu) <- c("x","y")
        warning("I'm renaming your 'mu' columns 'x' and 'y'.")} else 
    stop("Apologies for stickling, but if you're not going to give me a two column or complex 'mu', you need to at least name your relevant columns 'x' and 'y'.")
  
  
  z.p <- sqrt(-2 * log(0.05))
  s <- A/(pi * z.p^2)
  
  tau.v <- tau["tau.v"]
  tau.z <- tau["tau.z"]
  if(is.null(tau.v)) tau.v <- 0 else if(is.na(tau.v)) tau.v <- 0 
  if(is.null(tau.z)) tau.z <- 0 else if(is.na(tau.z)) tau.z <- 0 
  
  model <- "mouf"
  if(tau.v == 0) model <- "ou" 
  if(tau.v == 0  & tau.z == 0 ) model <- "wn" 
  if(tau.v > tau.z) warning("You have proposed a condition which is difficult to wrap one's head around in which the position auto-correlation is greater than the pseudo-velocity correlation.  Problems may ensue.")
  
  S <- s * outer(T, T, getCov, p=tau, model=model)
  X <- mvrnorm2(n = 1, mu[,"x"], S)
  Y <- mvrnorm2(n = 1, mu[,"y"], S)
  track <- data.frame(T, X, Y)
  class(track) <- c("track", "data.frame")
  return(track)
}
