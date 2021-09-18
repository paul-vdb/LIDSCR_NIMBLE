#################################
# SCR as a Marked Poisson Process
# Fisher Model - Spatial Counts from Burgar et al.
#################################

library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)

source("NimbleFunctions.R")

## Load the Fisher data
###--- set directory and load files
FisherMatrix <- read.csv("../data/Fisher_60min2016_matrix.csv", header=T, row.names=1) 
###--- Fisher sampling overlapping with genetic data
n <- as.matrix(FisherMatrix)
dim(n)
capt.count <- n[-c(10,34),]

traps.df <- read.csv("../data/TDF.csv")
traps <- cbind(x = traps.df$Easting, y = traps.df$Northing)

traps <- traps[-c(10,34),]

coord.scale <- 1000
traps <- traps/coord.scale

# Create a buffer around the traps:
pts <- SpatialPoints(traps)
b.r <- buffer(pts, width = 15)
plot(b.r)
points(traps, col = "red", pch = 4)

cellsize <- c(1,1)

e <- as(extent(b.r), "SpatialPolygons")
grd <- as.data.frame(makegrid(e, "regular", cellsize = cellsize))
points(grd, col = "blue")
mask.pts <- SpatialPoints(grd)
mask.pts <- mask.pts[b.r]
mask <- coordinates(mask.pts)
attr(mask, "area") <- prod(cellsize)
area <- prod(cellsize)

plot(b.r)
points(traps, pch = 4, col = "red")

plot(mask)
points(traps, col = "red", pch = 4)

capt <- matrix(0, nrow = sum(capt.count), ncol = nrow(traps))
counts <- rowSums(capt.count)
indices <- NULL
j <- 1
for(i in 1:nrow(traps))
{
	tmp <- counts[i]
	if(tmp !=0){
		indices <- rbind(indices, cbind(j:(j+tmp-1), rep(i, tmp)) )
		j <- j + tmp
	}
}
capt[indices] <- 1
capture.data <- list()
capture.data$capt <- capt
StudyPeriod <- 64

omega <- apply(capt, 1, FUN = function(x){which(x == 1)})

studyArea <- nrow(mask)*area

M <- 200
J <- nrow(traps)
xlim = range(mask[,1])
ylim = range(mask[,2])


# Chandler and Royle Spatial Count model:
SpatialCountAlg2 <- nimbleCode( {
	# Priors:
	sigma ~ dunif(0,50)
	lambda ~ dunif(0,10)
	psi~dunif(0,1)
	tau <- 1/(sigma*sigma)
	
	for(k in 1:M){
		z[k] ~ dbern(psi)
		X[k,1] ~ dunif(xlim[1], xlim[2])
		X[k,2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- pow(X[k,1] - traps[1:J,1], 2) + pow(X[k,2] - traps[1:J,2], 2)
		hkj[k,1:J] <- lambda*exp(-d2[k,1:J] * tau/2 )*z[k]
	}

	# Model for Observation Process:
	for(j in 1:J){
		Hj[j] <- sum(hkj[1:M,j])*Time	
		counts[j] ~ dpois(Hj[j])
	}
	
	N <- sum(z[1:M])
	D <- N/area
} )

# van Dam-Bates Spatial Count Model using the Marked Poisson process formulation.
SC_MPP <- nimbleCode({
    lambda ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
    for(k in 1:M) {
        z[k] ~ dbern(psi)
        X[k, 1] ~ dunif(xlim[1], xlim[2])
        X[k, 2] ~ dunif(ylim[1], ylim[2])
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
        hkj[k,1:J] <- exp(-d2[k,1:J]*tau2)*lambda
        # Hazard rate for animal across all traps.
        Hk[k] <- sum(hkj[k,1:J])*Time
		Hkz[k] <- Hk[k]*z[k]	# Only put z here for purposes of node dependence and speed.
    }

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # trap probability given ID:
        # This one can certainly be just the ones trick for trap y[i].
		pobs[i] <- hkj[ID[i], omega[i]]
        ones[i] ~ dbern(pobs[i])
		ID[i] ~ dID()	# Dummy distribution to declare this as stochastic.
    }
	p <- exp(-sum(Hkz[1:M]))
    one ~ dbern(p)
	
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area
})


# Run the same model from Burgar et al. Spatial Count on fisher.
#----------------------------------------------------------------------
constants.sc <- list(
    J = nrow(traps),
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
	area = diff(xlim)*diff(ylim)/100
	)

data.sc <- list(
	counts = counts,
	z =  rep(NA, M)
)

SCModel <- nimbleModel(SpatialCountAlg2, constants.sc, data.sc)
SCconf <- configureMCMC(SCModel)
SCconf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))
# Use a block update on locations. Saves time.
SCconf$removeSamplers('X')
for(i in 1:M) SCconf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE)
SCRmcmc <- buildMCMC(SCconf)
SCCmodel <- compileNimble(SCModel)
SCCmcmc <- compileNimble(SCRmcmc, project = SCModel)
samples.sc <- runMCMC(SCCmcmc, 30000, nburnin = 10000, nchains = 3, thin = 1)

out.sc <- mcmc.list(list(as.mcmc(samples.sc[[1]]), as.mcmc(samples.sc[[2]]), as.mcmc(samples.sc[[3]])))
plot(out.sc[,c("N", "sigma", "lambda", "psi")])
# save(out.sc, file = "../output/fisher_sc.Rda")
# load("../output/fisher_sc.Rda")

# Run the same model from van Dam-Bates et al. Marked Poisson Process on fisher.
#----------------------------------------------------------------------
constants.mpp <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
    n_obs = length(omega),
	omega = omega,
	area = diff(xlim)*diff(ylim)/100
	)

data.mpp <- list(
    one = 1,
    ones = rep(1, length(omega)),
	z =  rep(NA, M),
	ID = rep(NA, length(omega))
)

# Need to initialize this model as the stochastic node for ID is kind of wrong...
inits <- function(){
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 0.5, 2)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- lambda*exp(-d2/(2*sigma^2))
	ID <- do.call('c', lapply(omega, FUN = function(x) {sample(1:M, 1, prob = hkj[,x])}))
	z <- rep(0, M)
	z[ID] <- 1
	psi <- length(unique(ID))/M
	list(
		lambda = lambda,
		sigma = sigma,
		psi = psi,
		X = X,
		z = z,
		ID = ID
    )
}

MPPModel <- nimbleModel(SC_MPP, constants.mpp, data.mpp, init = inits())
MPPconf <- configureMCMC(MPPModel)
MPPconf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))
# Use a block update on locations. Saves time.
# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
MPPconf$removeSamplers('X')
# for(i in 1:M){ 
	# SCconf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', 
		# silent = TRUE, control = list(adaptive = FALSE, scale = 0.5))
	# }
for(i in 1:M){ MPPconf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'sampler_myX2', silent = TRUE, 
	control = list(xlim = xlim, ylim = ylim, scale = 0.5, J = nrow(traps)))}	
# Optimized z sampler
MPPconf$removeSamplers('z')
MPPconf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)
# van Dam-Bates categorical sampler
MPPconf$removeSamplers('ID')
# Chandler and Royle Alg. 1 sampler.
# MPPconf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
# van Dam-Bates Alg.
MPPconf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
MPPRmcmc <- buildMCMC(MPPconf)
MPPCmodel <- compileNimble(MPPModel)
MPPCmcmc <- compileNimble(MPPRmcmc, project = MPPModel)

# MPPCmcmc$run(1000)
# mvSamples <- MPPCmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples)
# plot(out)

samples.mpp <- runMCMC(MPPCmcmc, 30000, nburnin = 10000, nchains = 3, 
	thin = 1, inits = list(inits(), inits(), inits()))
out.mpp <- mcmc.list(list(as.mcmc(samples.mpp[[1]]), as.mcmc(samples.mpp[[2]]), as.mcmc(samples.mpp[[3]])))
# save(out.mpp, file = "../output/fisher_mpp.Rda")
# load("../output/fisher_mpp.Rda")
plot(out.mpp[,c("D", "sigma", "lambda")])
dev.new()
plot(out.sc[,c("D", "sigma", "lambda")])