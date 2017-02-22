## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache=FALSE, fig.width = 6)
#  html_document:
#    number_sections: yes
#    toc: yes
#    toc_depth: 2

#output:
#  pdf_document:
#    toc: yes
#    number_sections: yes
#    toc_depth: 2
		
#output: rmarkdown::html_vignette
#vignette: >
#  %\VignetteIndexEntry{marcher}
#  %\VignetteEngine{knitr::rmarkdown}
#  \usepackage[utf8]{inputenc}

## ---- message=FALSE, warning=FALSE---------------------------------------
library(marcher)

## ------------------------------------------------------------------------
time <- 1:100

## ------------------------------------------------------------------------
mean.pars <- c(x1 = 0, y1 = 0, x2 = 10, y2 = 10, t1 = 45, dt = 10)

## ------------------------------------------------------------------------
Mean <- getMu(T = time, p.m = mean.pars)

## ----MeanPlot, echo=-1, fig.height = 3-----------------------------------
par(mar = c(2,4,0,2), oma = c(2,0,2,0), bty="l", xpd=NA)
scan.track(time = time, x = Mean)

## ------------------------------------------------------------------------
taus <- c(tau.z = 1, tau.v = 0.5)
area <- 50
SimTrack <- simulate.shift(T = time, tau = taus, mu = Mean, A = area)

## ------------------------------------------------------------------------
head(SimTrack)

## ----FirstSimTrack, echo=-1, fig.height = 4------------------------------
par(mar = c(2,4,0,2), oma = c(2,0,2,0), bty="l", xpd=NA)
with(SimTrack, scan.track(time= T, x = X, y = Y))

## ----SimulatedTracks, echo=-1, fig.height = 3----------------------------
par(bty="l"); set.seed(2015)
MWN.sim <- simulate.shift(T = time, tau = NULL, mu = Mean, A = area) %T>% scan.track()
title("No Autocorrelation: MWN", outer = TRUE)
MOU.sim <- simulate.shift(T = time, tau = c(tau.z = 5), mu = Mean, A = area)  %T>% scan.track()
title("Position Autocorrelation: MOU", outer = TRUE)
MOUF.sim <- simulate.shift(T = time, tau = c(tau.z = 5, tau.v = 1), mu = Mean, A = area) %T>% scan.track
title("Position and Velocity Autocorrelation: MOUF", outer = TRUE)

## ---- echo=-1, fig.height = 3--------------------------------------------
par(bty="l")

time.random <- sort(runif(100,0,100))
mean.random <- getMu(T = time.random, p.m = mean.pars)
MOUF.sim.random <- simulate.shift(T = time.random, tau = c(tau.z = 5, tau.v = 2), mu = mean.random, A = area) %T>% scan.track;
title("MOUF: random times", outer = TRUE)

## ---- echo=-1, fig.height = 3--------------------------------------------
	par(bty="l"); set.seed(2015)
  p.m <- c(x1 = 0, x2 = 5, x3 = 10, y1 = 0, y2 = 5, y3 = 0, t1 = 40, t2 = 130, dt1 = 6, dt2 = 4)
  A <- 20; T <- 1:200
  Mu <- getMu.multi(T, p.m)
  MOU.3range <- simulate.shift(T, tau=c(tau.z = 5), mu=Mu, A=A) %T>% scan.track
	title("MOU: Three Ranges", outer = TRUE)

## ----examplefits---------------------------------------------------------
MWN.fit <- with(MWN.sim, estimate.shift(T=T, X=X, Y=Y))
summary(MWN.fit)

## ----quickfit1, echo=-1, cache=FALSE-------------------------------------
par(bty="l")
(p.m0 <- with(MWN.sim, quickfit(T,X,Y, dt = 1)))

## ----quickfit2, echo=-1, cache=FALSE-------------------------------------
par(bty="l"); set.seed(1979)
(p.m0 <- with(MOU.3range, quickfit(T,X,Y, dt = 1, n.clust = 3)))

## ---- eval=FALSE---------------------------------------------------------
#  (p.m0 <- with(MWN.fit, locate.shift(T,X,Y, dt = 1)))

## ---- eval=FALSE---------------------------------------------------------
#  estimate.shift(T, X, Y, p.m0 = p.m0)

## ----ThreeResidualFits, eval=TRUE, cache=FALSE, echo = -1, fig.height = 3----
par(mfrow= c(1,3), bty="l", echo=-1)
data(SimulatedTracks)
MWN.res <- with(MWN.sim, estimate.shift(T, X, Y, model = "wn"))$Z.res
MOU.res <- with(MOU.sim, estimate.shift(T, X, Y, model = "wn"))$Z.res
MOUF.res <- with(MOUF.sim, estimate.shift(T, X, Y, model = "wn"))$Z.res

plot(MWN.res, asp=1, type="o")
plot(MOU.res, asp=1, type="o")
plot(MOUF.res, asp=1, type="o")

## ----Michela, cache = FALSE, echo=-1, fig.height = 3---------------------
par(bty="l")
data(Michela)
with(Michela, scan.track(time = time, x = x, y = y))

## ------------------------------------------------------------------------
# first day of first year
day1 <- ymd(paste(min(year(Michela$time)), 1, 1))
# find time difference and round (and make numeric)
days <- as.numeric(round(difftime(Michela$time, day1, unit = "day")))

## ----Michela.fit, cache = FALSE, fig.height = 3, echo = -1---------------
par(bty="l")
# p.m0 <- with(Michela, locate.shift(time = day, x = x, y = y, n.clust = 3))
p.m0 <- c(x1 = 653.6,  x2 = 658.85, x3 = 653.8, 
          y1 = 5094.8, y2 = 5096,   y3 = 5094.8, 
					t1 = 118, t2 = 318, dt1 = 10.5, dt2 = 15.8)
Michela.fit <- with(Michela, 
							 estimate.shift(day, x, y, n.clust = 3, model = "ou", method = "like", p.m0 = p.m0))
plot(Michela.fit)

## ----message = FALSE, warning = FALSE------------------------------------
ymd("2005 1 1") + ddays(Michela.fit$p.hat[c('t1','t2')])

## ------------------------------------------------------------------------
getRSI(Michela.fit, 1, 2)

## ------------------------------------------------------------------------
data.frame(FirstToSecond = getRSI(Michela.fit, 1,2)[,1],
           SecondToThird= getRSI(Michela.fit, 2,3)[,1], 
           FirstToThird = getRSI(Michela.fit, 1,3)[,1])

## ----ReturnTest, cache = FALSE-------------------------------------------
test.return(Michela.fit)

## ----MichelaAIC----------------------------------------------------------
c(logLik(Michela.fit.early), AIC(Michela.fit.early))

## ----Michela2b.aic-------------------------------------------------------
Michela.fit.early2 <- with(subset(Michela, day < 200), 
					estimate.shift(day, x, y, n.clust = 2, model = "ou", method = "like"))
c(logLik(Michela.fit.early2), AIC(Michela.fit.early2))

## ----Stopover.test-------------------------------------------------------
test.stopover(Michela.fit.early)

## ---- fig.height = 3, echo = -1------------------------------------------
par(bty="l")
# set initial parameters
A <- 20; T <- 1:100; tau <- c(tau.z = 2, tau.v = 0)

# simulate, fit and test clear disperal
Mu <- getMu(T, c(x1 = 0, y1 = 0, x2 = 4, y2 = 4, t1 = 40, dt = 20))
XY.sim <- simulate.shift(T, tau = tau, Mu, A=A)
with(XY.sim, scan.track(time = T, x = X, y = Y))
with(XY.sim, fitNSD(T, X, Y, plotme=TRUE))


# simulate, fit and test no disperal
Mu <- getMu(T, c(x1 = 0, y1 = 0, x2 = 0, y2 = 0, t1 = 40, dt = 20))
XY.sim <- simulate.shift(T, tau = tau, Mu, A=A)
with(XY.sim, scan.track(time = T, x = X, y = Y))
with(XY.sim, fitNSD(T,X,Y, plotme=TRUE))

## ---- echo=FALSE, eval=FALSE---------------------------------------------
#  require(rmarkdown)
#  setwd("C:/Users/Elie/Box Sync/Rpackages/ecomove/marcher/vignettes")
#  render("marcher.rmd", encoding = "UTF-8", html_document(toc = TRUE)); setwd("./")
#  shell("marcher.html")
#  require(knitr)
#  purl("marcher.rmd")

