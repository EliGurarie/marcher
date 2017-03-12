#' Range shift hypothesis tests
#'
#' Three tests for three hypotheses to test on fitted range shifts: Was the range shift significant?  Did an animal that performed two consecutive seasonal migrations return to the same location it began?  Was there a stopover during a migration? 
#' @param FIT a fitted range shift (output of \code{\link{estimate_shift}})
#' @param verbose whether to print verbose message
#' @return Outputs a summary of the test results and returns a list of test results including: 
#' \itemize{
#' \item{\code{aic.table}} {an AIC table comparing models}
#' \item{\code{lrt}} {a likelihood ratio test statistic}
#' \item{\code{df}} {degrees of freedom for the l.r.t.}
#' \item{\code{p.value}} {a p.value for the l.r.t.}
#' }
#' 
#' @aliases test_rangeshift test_return test_stopover

#' @describeIn test_rangeshift Compare a two range fitted model to a null model of no range shift.
#' @export
test_rangeshift <- function(FIT, verbose = TRUE){
  #method one of "ar" or "like", for AR equivalence or likelihood method
	model <- FIT$model
	method <- FIT$method
  p.hat <- FIT$p.hat
  ll <- FIT$ll
  k <- FIT$df 	
	aic <- FIT$aic
	
  test.table <- data.frame(ll, k, aic)

	# now - compute null model
		Z.res.null <- (FIT$X + 1i*FIT$Y)
		Z.res.null <- Z.res.null - mean(Z.res.null)
		
		# turn tau.hat into appropriate form
		
		p.s.null <- NULL
		if(model != "wn"){
			tau.hat.null <- getTau(Z.res.null, T = FIT$T, model = model, method = method)$tau.hat
			if(model %in% c("ou", "ouf")) p.s.null <- c(logtau.z = log(t(tau.hat.null['tau.z'])))
			if(model == "ouf") p.s.null["kappa"] <-  tau.hat.null["tau.v"] / tau.hat.null["tau.z"]
		}
		
		ll.null <- getLikelihood.res(p.s = p.s.null, T=FIT$T, Z.res=Z.res.null, model=model)
		k.null <- k-4
		aic.null <-  -2*ll.null + 2*k.null
		test.table <- data.frame(ll = c(ll.null, ll), k = c(k.null, k), aic = c(aic.null, aic))
		row.names(test.table) <- c("with shift", "without shift")
		
	# construct test table
		
		lrt <- with(test.table, 2*(ll[1] - ll[2]))
		p.value <- 1-pchisq(lrt, 2)

		if(verbose){
			print(test.table)
			cat(paste0("
l.r.t.: ", signif(lrt,3), " with ", 4, " degrees of freedom, ",
								 "p-value: ", signif(p.value,4), "\n"))
								 
			if(p.value < 0.05 & with(test.table, aic[2] < aic[1])) cat("There is almost certainly a significant range shift.\n") else 
			if(p.value > 0.05 & with(test.table, aic[2] > aic[1])) cat("There probably wasn't a range shift.\n") else cat("Tough call - there's disagreement in the criteria.\n")
		}
		invisible(list(aic.table = test.table,  lrt = lrt, df = 4, p.value = p.value))
 }

#' @describeIn test_rangeshift Compares a three range fitted model in which the first and third ranges have the same centroid against a model where the first and third centroid are different.
#' @export
test_return <- function(FIT, verbose = TRUE){

	method <- FIT$method
	X <- FIT$X
	Y <- FIT$Y
	T <- FIT$T
	model <- FIT$model
	p.mu.hat <- FIT$p.hat
  n <- length(T)
	Z.res <- FIT$Z.res
	ll <- FIT$ll
	aic <- FIT$aic
	df <- FIT$df
	
	# now ... fit the model with only ONE return
	
	p0.mu.return <- p.mu.hat[- (which(names(p.mu.hat) %in% c('A', 'x3', 'y3', 'tau.z', 'tau.v')))]
	p.mu.return.fit <- getP.mu(T, X, Y, p.m0 = p0.mu.return, FUN = getMu.return, bounds = TRUE)
	
	XY.return.hat <- getMu.return(T, p.mu.return.fit$par)
	Z.return.hat <- XY.return.hat[,1] + 1i* XY.return.hat[,2]
	Z.res.return <- X + 1i*Y - Z.return.hat
	tau.fit.return <- try(getTau(Z.res.return, T=T, model=model, method = method))
	ll.return <- tau.fit.return$ll
	aic.return <- -2*ll.return - 2*(FIT$df - 2)
	
	# construct test table
		
	test.table <- data.frame(ll = c(ll.return, ll), 
												   df = c(df - 2, df), 
													 aic = c(aic.return, aic))
													 
	row.names(test.table) <- c("two range with return", "three ranges")
	  
	lrt <- with(test.table, 2*(ll[2] - ll[1]))
	p.value <- 1-pchisq(lrt, 2)

	if(verbose){
		print(test.table)
		cat(paste0("
l.r.t.: ", signif(lrt,3), " with ", 4, " degrees of freedom, ",
							 "p-value: ", signif(p.value,4), "\n"))
							 
		if(p.value < 0.05 & with(test.table, aic[2] < aic[1])) cat("A model with three unique ranges is better.\n") else 
		if(p.value > 0.05 & with(test.table, aic[2] > aic[1])) cat("The second migration was most likely a return migration.\n") else cat("Tough call - there's disagreement in the criteria.\n")
	}
	
	invisible(list(aic.table = test.table, lrt = lrt, df = 2, p.value = p.value))
}


#' @describeIn test_rangeshift Compare a three range model with an apparent stopover (shorter intermediate range), and see if a more parsimonious model excludes the stopover.
#' @export
 test_stopover <- function(FIT, verbose = TRUE){
	
	method <- FIT$method
	p.hat <- FIT$p.hat
	ll <- FIT$ll
	model <- FIT$model 
	Z <- with(FIT, X + 1i*Y)
	
	# fit two way
  p.m0.null <- c(p.hat[c("t1", "dt1", "x1", "y1")], x2 = t(p.hat["x3"]), y2 = t(p.hat["y3"]))
	names(p.m0.null)[names(p.m0.null) == "dt1"] <- "dt"
		
	FIT.nostopover	<- with(FIT, estimate_shift(T, X, Y, n.clust = 2, p.m0 = p.m0.null,
													CI = FALSE, model = model, method = method))
	
	test.table <- data.frame(ll = c(FIT.nostopover$ll, ll), 
	           aic = c(FIT.nostopover$aic, FIT$aic))
	row.names(test.table) <- c("with stopover", "without stopover")
	  
	lrt <- with(test.table, 2*(ll[2] - ll[1]))
	p.value <- 1-pchisq(lrt, 4)

	if(verbose){
		print(test.table)
		cat(paste0("
l.r.t.: ", signif(lrt,3), " with ", 4, " degrees of freedom, ",
							 "p-value: ", signif(p.value,4), "\n"))
		if(p.value < 0.05 & with(test.table, aic[2] < aic[1])) cat("A three range model (with stopover) is a better model.\n") else 
		if(p.value > 0.05 & with(test.table, aic[2] > aic[1])) cat("The stopover is likely illusory.\n") else cat("Tough call - there's disagreement in the criteria.\n")
	}
	
	invisible(list(aic.table = test.table, lrt = lrt, df = 4, p.value = p.value))
 }
 