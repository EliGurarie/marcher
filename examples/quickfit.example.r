require(marcher)

## Load simulated data
	data(SimulatedTracks)
  
	# plot the MOU simulation
	scan.track(MOU.sim)
	
  # quick fit - setting dt = 10
  (pm.0 <- with(MOU.sim, quickfit(T, X, Y, dt = 10)))
	
	# interactive locator process
	**## Not run:** 
	(with(MOU.sim, locate.shift(T, X, Y))
	## End(**Not run**)
	
	# fit the model
  fit <- with(myMouf, estimate.shift(T, X, Y))
  

## Three cluster example
	
	# plot the three range shift simulation
	scan.track(MOU.3range)
 
	# quick fit 
	## (note - this may not always work!)
  with(MOU.3range, quickfit(T, X, Y, dt = 10, n.clust = 3))
	
	**## Not run:** 
	with(MOU.3range, locate.shift(T, X, Y, n.clust = 3))
	## End(**Not run**)
	
  
  
