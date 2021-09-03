simSCR <- function(N = 50, sigma = 0.5, lambda = 0.5, StudyPeriod = 25, traps, limits, size = c(22, 5, 0.75), psex = 0.5)
{
    locs <- cbind(x = runif(N, limits[['xlim']][1], limits[['xlim']][2]), 
                  y = runif(N, limits[['ylim']][1], limits[['ylim']][2]))
    J <- nrow(traps)
    capt.hist <- data.frame()
    ID = 1
    for(i in 1:N)
    {
        d2 <- (locs[i,1] - traps[,1])^2 + (locs[i,2] - traps[,2])^2
        hkj <- lambda*exp(-d2/(2*sigma^2))
        hk <- sum(hkj)
        ti <- cumsum(rexp(1000, hk))   #Simulate detection times.
        ti <- ti[ti < StudyPeriod]
        nk <- length(ti)
        if(nk == 0) next;
        # Now assign those detection times to a trap.
        obs <- sample(J, nk, prob = hkj/hk, replace = TRUE)
        size.i <- rnorm(1, size[1], size[2])
        capt.hist <- rbind(capt.hist, data.frame('t_obs' = ti, 'trap_obs' = obs, 'ID' = ID, 
                                                 'size' = rnorm(nk, size.i, size[3]), 
												'sex' = rbinom(1, 1, psex),
												'tru_x' = as.numeric(locs[i,1]), 'tru_y' = as.numeric(locs[i,2])))
        ID <- ID + 1
    }
    capt.hist
}

sim_scenario <- function(scen = "A")
{
	if(scen == "A"){
		side <- 15
		coords <- seq(1, 15, length=side)
		traps <- cbind(x=rep(coords, each=side), y=rep(coords, times=side))
		buffer <- 3
		limits <- list('xlim' = c(min(traps[,1]-buffer), max(traps[,1]+buffer)), 
					   'ylim' = c(min(traps[,2]-buffer), max(traps[,2]+buffer)))
		area <- diff(limits[['xlim']])*diff(limits[['ylim']])
		lambda <- lambdaTrue <- 0.5
		sigma <- sigmaTrue <- 0.5
		N <- NTrue <- 45
		StudyPeriod <- 5
	}else{
		traps <- expand.grid(x = 1:6, y = 1:6, KEEP.OUT.ATTRS = FALSE)
		buffer <- 1.25
		limits <- list('xlim' = c(min(traps[,1]-buffer), max(traps[,1]+buffer)), 
					   'ylim' = c(min(traps[,2]-buffer), max(traps[,2]+buffer)))

		## Parameters:
		lambda <- lambdaTrue <- 0.65
		sigma <- sigmaTrue <- 0.5
		N <- NTrue <- 20
		StudyPeriod <- 10
		area <- diff(limits[['xlim']])*diff(limits[['ylim']])
	}
	sim.dat <- simSCR(N = N, sigma = sigma, lambda = lambda, StudyPeriod = StudyPeriod, traps, limits)
	return(list(lambdaTrue = lambdaTrue, sigmaTrue = sigmaTrue, NTrue = NTrue, StudyPeriod = StudyPeriod,
		area = area, sim.dat = sim.dat, traps = traps, limits = limits))
}

simASCR <- function(N = 50, sigma = 0.5, sigma_toa = 0.01, g0 = 1, lambda = 0.5, StudyPeriod = 25, traps, limits, psex = 0.5, nu = 330)
{
    locs <- cbind(x = runif(N, limits[['xlim']][1], limits[['xlim']][2]), 
                  y = runif(N, limits[['ylim']][1], limits[['ylim']][2]))
    J <- nrow(traps)
    obs <- data.frame()
	pop <- data.frame()
	capthist <- toa <- NULL
    ID = 1
    for(i in 1:N)
    {
        d2 <- (locs[i,1] - traps[,1])^2 + (locs[i,2] - traps[,2])^2
        pkj <- 1-exp(-g0*exp(-d2/(2*sigma^2)))
		p. <- 1-prod(1-pkj)
        ti <- cumsum(rexp(1000, lambda))   #Simulate detection times.
        ti <- ti[ti < StudyPeriod]
        nk <- length(ti)
		keep <- 0
		if(nk != 0) {
			# Now assign those detection times to a trap.
			capt <- do.call("rbind", lapply(1:nk, FUN = function(x){rbinom(J,1,pkj)}))
			keep <- rowSums(capt) != 0
			capt <- capt[keep, ]
			if(sum(keep) > 0) { 
				obs <- rbind(obs, data.frame('t_obs' = ti[keep],'ID' = ID,
														'sex' = rbinom(1, 1, psex),
														'tru_x' = as.numeric(locs[i,1]), 'tru_y' = as.numeric(locs[i,2])))
				capthist <- rbind(capthist, capt)
				
				toa.i <- do.call("rbind", lapply(ti[keep], FUN = function(x){x + sqrt(d2)/nu + rnorm(J,0,sigma_toa)}))
				toa <- rbind(toa, toa.i*capt)
			}
		}
		nk <- sum(keep)
		pop <- rbind(pop, data.frame('n_obs' = nk,'ID' = ID*(nk > 0),
									'tru_x' = locs[i,1], 'tru_y' = locs[i,2]))
		if(nk > 0) 	ID <- ID + 1 							
    }
    list(capt = capthist, toa = toa, obs = obs, pop = pop)
}

# Stolen code for making habitat mask but without loading package:
# Convert secr mask and traps objects for use with Nimble
convertMask <- function(secrmask, secrtraps, plot = TRUE) {

  pixWidth <- min(abs(diff(secrmask$x)))
  bbox <- matrix(c(min(secrmask$x), max(secrmask$x), min(secrmask$y), max(secrmask$y)), nrow = 2, ncol = 2)
  # Create 'false origin' so that SW corner of matrix is at [1, 1]
  origin <- bbox[1, ]

  # Get dimensions of the matrix
  nrows <- ceiling((bbox[2, 1] - bbox[1, 1]) / pixWidth) + 1
  ncols <- ceiling((bbox[2, 2] - bbox[1, 2]) / pixWidth) + 1
  habMat <- matrix(0L, nrow=nrows, ncol=ncols)
  # Convert mask x and y to col/row numbers
  dex <- as.matrix(ceiling(sweep(secrmask, 2, origin) / pixWidth)) + 1
  for(i in 1:nrow(dex)) habMat[dex[i,1], dex[i,2]] <- 1L

  # Convert trap coordinates to the new units:
  newtraps <- sweep(secrtraps, 2, origin) / pixWidth + 1

  out <- list(habMat = habMat, 
              trapMat = as.matrix(newtraps), 
              upperLimit = c(x=nrows+1, y=ncols+1),
              pixelWidth = pixWidth,
              area = sum(habMat) * pixWidth^2)
  attr(out, "boundingbox") <- bbox
  attr(out, "origin") <- origin
  attr(out, "pixelWidth") <- pixWidth
  
  if(plot) {
	image(habMat)
  }
  return(out)
}



rowProd <- function(x)
{
	if(is.null(dim(x))) return(x)
	apply(x, 1, prod)
}
