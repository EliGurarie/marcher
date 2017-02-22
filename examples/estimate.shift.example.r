# load simulated tracks

data(SimulatedTracks)

# white noise fit
MWN.fit <- with(MWN.sim, estimate.shift(T=T, X=X, Y=Y))
summary(MWN.fit)
plot(MWN.fit)

# OUF fit
MOUF.fit <- with(MOUF.sim.random, estimate.shift(T=T, X=X, Y=Y, model = "ouf", method = "like"))
summary(MOUF.fit)
plot(MOUF.fit)

# Three range fit
## it is helpful to have some initital values for these parameters because the automated quickfit() method is unreliable for three ranges
## in the example, we set a seed that seems to work
set.seed(1976)
MOU.3range.fit <- with(MOU.3range, estimate.shift(T=T, X=X, Y=Y, model = "ou", method = "ar", n.clust = 3))
summary(MOU.3range.fit)
plot(MOU.3range.fit)

