#' Compute Range Shift Index
#' 
#' 
#' The range shift index is a dimensionless measure of the distance of the centroids of two ranges divided by the diameter of the 95\% area. This function uses the 95\% confidence intervals from a range shift fit to calculate a point estimate and 95\% confidence intervals of the RSI. 
#'
#' @param FIT  a rnage shift object, outputted by \code{\link{estimate.shift}}
#' @param {n1,n2} the indices of the ranges to estimate from and to, i.e., for single shift, 1 and 2.  For three ranges (two shifts) it can be 1 and 2, 2 and 3, or 1 and 3 - if the ultimate shift is the one of interest. 
#' @return returns a data frame reporting the distance traveled, the RSI and respective bootstrapped confidence intervals. 


getRSI <- function(FIT, n1=1, n2=2, nboot = 1e3){

	p.hat <- FIT$p.hat
	par.CI <- FIT$p.CI
	
	getMCsample <- function(par, n=1e3)
		{
			p.CI <- (par.CI %>% t %>% as.data.frame)[par][,1]
			p.sd <- mean(p.CI[1]-p.CI[2], p.CI[3]-p.CI[1])/2
			list(hat = p.CI[1], MC = rnorm(n, p.CI[1], p.sd))
		}

		x1 <- paste0("x",n1)
		x2 <- paste0("x",n2)
		y1 <- paste0("y",n1)
		y2 <- paste0("y",n2)
	
		if(is.null(par.CI)){
				z1 <- p.hat[x1] + 1i*p.hat[y1]
				Z2 <- p.hat[x2] + 1i*p.hat[y2]
				D.hat <- Mod(z2 - z1)
				RSI.hat <- D.hat / sqrt(4*p.hat["A"]/pi)
				return(c(D = D.hat, RSI = RSI.hat))
		} else {
		
			A.mc <- getMCsample("A", nboot)$MC
			Z1.mc <- getMCsample(x1, nboot)$MC + 1i *  getMCsample(y1, nboot)$MC
			Z2.mc <- getMCsample(x2, nboot)$MC + 1i *  getMCsample(y2, nboot)$MC

			A.hat <- getMCsample("A", nboot)$hat
			Z1.hat <- getMCsample(x1, nboot)$hat + 1i *  getMCsample(y1, nboot)$hat
			Z2.hat <- getMCsample(x2, nboot)$hat + 1i *  getMCsample(y2, nboot)$hat		
			
			D.mc <- Mod(Z2.mc - Z1.mc)
			RSI.mc <- D.mc / sqrt(4*A.mc/pi)
			RSI.hat <- Mod(Z2.hat - Z1.hat) / sqrt(4*A.hat/pi)
			return(data.frame(
								D = c(p.hat =  Mod(Z2.hat - Z1.hat), CI.low = as.numeric(quantile(D.mc,   0.025)), CI.high = as.numeric(quantile(D.mc,   0.975))),
							RSI = c(p.hat = RSI.hat, 						   CI.low = as.numeric(quantile(RSI.mc, 0.025)), CI.high = as.numeric(quantile(RSI.mc, 0.975)))) %>% t
							)			
  }
	}
