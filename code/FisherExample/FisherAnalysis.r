#################################
# SCR as a Marked Poisson Process
# Fisher Model - Spatial Counts from Burgar et al.
#################################
setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code")
library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)

source("NimbleFunctions.R")

## Load the Fisher data
###--- set directory and load files
load("../data/FisherData.Rda")

xlim <- fisher.data[['xlim']] 
ylim <- fisher.data[['ylim']]
traps <- fisher.data[["traps"]]
obs <- fisher.data[["observations"]]
mustlink <- fisher.data[["mustlink"]] 
cannotlink <- fisher.data[["cannotlink"]] 

omega <- obs$TrapNumber
StudyPeriod <- 64
studyArea <- diff(xlim)*diff(ylim)

M <- 400
J <- nrow(traps)

counts <- NULL
for(j in 1:J)
{
	counts <- c(counts, sum(obs$TrapNumber == j))
}

# Chandler and Royle Spatial Count model:
SpatialCountAlg2 <- nimbleCode( {
	# Priors:
	sigma ~ dunif(0,50)
	lambda ~ dunif(0,10)
	psi ~ dbeta(1,1)
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
    tau2 <- 1/(2*sigma^2)	# Just avoid that extra computation for each animal...
    for(k in 1:M) {
        z[k] ~ dbern(psi)
        X[k, 1] ~ dunif(xlim[1], xlim[2])
        X[k, 2] ~ dunif(ylim[1], ylim[2])
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
        pkj[k,1:J] <- exp(-d2[k,1:J]*tau2)
        # Hazard rate for animal across all traps.
        Hk[k] <- sum(pkj[k,1:J])*Time*lambda
		pkz[k] <- exp(-Hk[k]*z[k])	# Only put z here for purposes of node dependence and speed.
		zones[k] ~ dbern(pkz[k])
   }

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # trap probability given ID:
        # This one can certainly be just the ones trick for trap y[i].
		pobs[i] <- pkj[ID[i], omega[i]]
        ones[i] ~ dbern(pobs[i])
		ID[i] ~ dID(lam = lambda)	# Dummy distribution to declare this as stochastic and mulitply by lambda.
    }
	
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area
})

# van Dam-Bates Spatial Count Model using the Marked Poisson process formulation.
SC_MPP_sex <- nimbleCode({
    lambda ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma[1] ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    sigma[2] ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    tau2[1] <- 1/(2*sigma[1]^2)	# Just avoid that extra computation for each animal...
    tau2[2] <- 1/(2*sigma[2]^2)	# Just avoid that extra computation for each animal...
	psex ~ dbeta(1,1)
    for(k in 1:M) {
        z[k] ~ dbern(psi)
		sex[k] ~ dbern(psex)
        X[k, 1] ~ dunif(xlim[1], xlim[2])
        X[k, 2] ~ dunif(ylim[1], ylim[2])
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
        pkj[k,1:J] <- exp(-d2[k,1:J]*tau2[sex[k] + 1])
        # Hazard rate for animal across all traps.
        Hk[k] <- sum(pkj[k,1:J])*Time*lambda
		pkz[k] <- exp(-Hk[k]*z[k])	# Only put z here for purposes of node dependence and speed.
		zones[k] ~ dbern(pkz[k])
   }

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # trap probability given ID:
        # This one can certainly be just the ones trick for trap y[i].
		pobs[i] <- pkj[ID[i], omega[i]]
        ones[i] ~ dbern(pobs[i])
		ID[i] ~ dID(lam = lambda)	# Dummy distribution to declare this as stochastic and mulitply by lambda.
    }
	
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

SCconf$removeSamplers('sigma')
SCconf$addSampler(target = 'sigma', 
	type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))

# Use a block update on locations. Saves time.
SCconf$removeSamplers('X')
for(i in 1:M) SCconf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE, control = list(scale = 1.7, adaptive = FALSE))
SCRmcmc <- buildMCMC(SCconf)
SCCmodel <- compileNimble(SCModel)
SCCmcmc <- compileNimble(SCRmcmc, project = SCModel)

# SCCmcmc$run(10000)
# SCconf$printSamplers()
# scales <- NULL
# for(i in 1:M){
	# scales <- c(scales,valueInCompiledNimbleFunction(SCCmcmc$samplerFunctions[[403+i]], "scale"))
# }
# Note that when doing adaptive sampling it was around scale ~= 1.7. We will set it at that for all chains.
samples.sc <- runMCMC(SCCmcmc, 30000, nburnin = 10000, nchains = 3, thin = 1)

out.sc <- mcmc.list(list(as.mcmc(samples.sc[[1]]), as.mcmc(samples.sc[[2]]), as.mcmc(samples.sc[[3]])))
plot(out.sc[,c("N", "sigma", "lambda", "psi")])
plot(out.sc[,c("N", "D", "sigma")])
summary(out.sc)
effectiveSize(out.sc)
stat <- NULL
for( i in c('sigma', 'lambda', 'psi', 'N', 'D')) stat[[i]] = gelman.diag(mcmc.list(list(out.sc[[1]][,i], out.sc[[2]][,i], out.sc[[3]][,i] )))


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
    zones = rep(1, M),
    ones = rep(1, length(omega)),
	z =  rep(NA, M),
	ID = rep(NA, length(omega))
)

# Need to initialize this model as the stochastic node for ID is kind of wrong...
inits <- function(){
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 1, 5)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
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

######################################
# Chandler and Royle 2013 Algorithm 1:
######################################
MPPModel <- nimbleModel(SC_MPP, constants.mpp, data.mpp, inits = inits())
MPPconf <- configureMCMC(MPPModel)
MPPconf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID'))
# Use a block update on locations. Saves time.
# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
MPPconf$removeSamplers('X')
for(i in 1:M){
	MPPconf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'RW_block', silent = TRUE, control = list(scale = 1.7, adaptive = FALSE))
	}

MPPconf$removeSamplers('sigma')
MPPconf$addSampler(target = 'sigma', 
	type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))

# Optimized z sampler
MPPconf$removeSamplers('z')
MPPconf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)
# van Dam-Bates categorical sampler
MPPconf$removeSamplers('ID')
# Chandler and Royle Alg. 1 sampler.
MPPconf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
# van Dam-Bates Alg.
MPPRmcmc <- buildMCMC(MPPconf)
MPPCmodel <- compileNimble(MPPModel)
MPPCmcmc <- compileNimble(MPPRmcmc, project = MPPModel)

# MPPCmcmc$run(10000)
# mvSamples <- MPPCmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples)
# plot(out[,c('sigma', 'lambda', 'N', 'D')])
# post.id <- samples[,grep("ID", colnames(samples))]
# NActive <- apply(post.id, 1, FUN = function(x){ length(unique(x))})
# which(NActive > samples[,'N'])



samples.mpp.alg1 <- runMCMC(MPPCmcmc, 30000, nburnin = 10000, nchains = 3, 
	thin = 1, inits = list(inits(), inits(), inits()))

post.id.alg1.1 <- samples.mpp.alg1[[1]][,grep("ID", colnames(samples.mpp.alg1[[1]]))]
post.id.alg1.2 <- samples.mpp.alg1[[2]][,grep("ID", colnames(samples.mpp.alg1[[2]]))]
post.id.alg1.3 <- samples.mpp.alg1[[3]][,grep("ID", colnames(samples.mpp.alg1[[3]]))]
NActive1 <- apply(post.id.alg1.1, 1, FUN = function(x){ length(unique(x))})
NActive2 <- apply(post.id.alg1.2, 1, FUN = function(x){ length(unique(x))})
NActive3 <- apply(post.id.alg1.3, 1, FUN = function(x){ length(unique(x))})

out.mpp.alg1 <- mcmc.list(list(as.mcmc(cbind(samples.mpp.alg1[[1]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "NActive" = NActive1)), 
	as.mcmc(cbind(samples.mpp.alg1[[2]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "NActive" =  NActive2)),
	as.mcmc(cbind(samples.mpp.alg1[[3]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "NActive" = NActive3))))
plot(out.mpp.alg1, ask = TRUE)
summary(out.mpp.alg1)
# save(out.mpp.alg1, file = "../output/fisher_MPP_Alg1.Rda")	
# load("../output/fisher_MPP_Alg1.Rda")	

#################################
# New algorithm:
#################################
MPPModel <- nimbleModel(SC_MPP, constants.mpp, data.mpp, inits = inits())
MPPconf <- configureMCMC(MPPModel)
MPPconf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID'))

MPPconf$removeSamplers('sigma')
MPPconf$addSampler(target = 'sigma', 
	type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))

# Use a block update on locations. Saves time.
# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
MPPconf$removeSamplers('X')
for(i in 1:M){
	MPPconf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'RW_block', silent = TRUE, control = list(scale = 1.7, adaptive = FALSE))
	}
# Optimized z sampler
MPPconf$removeSamplers('z')
MPPconf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)
# van Dam-Bates categorical sampler
MPPconf$removeSamplers('ID')
MPPconf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
MPPRmcmc <- buildMCMC(MPPconf)
MPPCmodel <- compileNimble(MPPModel)
MPPCmcmc <- compileNimble(MPPRmcmc, project = MPPModel)

# MPPCmcmc$run(1000)
# MPPconf$printSamplers()
# scales <- NULL
# for(i in 1:M){
	# scales <- c(scales,valueInCompiledNimbleFunction(MPPCmcmc$samplerFunctions[[3+i]], "scale"))
# }
# mvSamples <- MPPCmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples)
# plot(out[,c('sigma', 'lambda', 'D', 'N')])
# post.id <- samples[,grep("ID", colnames(samples))]
# NActive <- apply(post.id, 1, FUN = function(x){ length(unique(x))})
# hist(NActive - samples[,'N'])

samples.mpp <- runMCMC(MPPCmcmc, 30000, nburnin = 10000, nchains = 3, 
	thin = 1, inits = list(inits(), inits(), inits()))

post.id.1 <- samples.mpp[[1]][,grep("ID", colnames(samples.mpp[[1]]))]
post.id.2 <- samples.mpp[[2]][,grep("ID", colnames(samples.mpp[[2]]))]
post.id.3 <- samples.mpp[[3]][,grep("ID", colnames(samples.mpp[[3]]))]
NActive1 <- apply(post.id.1, 1, FUN = function(x){ length(unique(x))})
NActive2 <- apply(post.id.2, 1, FUN = function(x){ length(unique(x))})
NActive3 <- apply(post.id.3, 1, FUN = function(x){ length(unique(x))})

out.mpp <- mcmc.list(list(as.mcmc(cbind(samples.mpp[[1]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "NActive" = NActive1)), 
	as.mcmc(cbind(samples.mpp[[2]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "NActive" =  NActive2)),
	as.mcmc(cbind(samples.mpp[[3]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "NActive" = NActive3))))

effectiveSize(out.mpp)
stat <- NULL
for( i in c('sigma', 'lambda', 'psi', 'N', 'D', 'NActive')) stat[[i]] = gelman.diag(mcmc.list(list(out.mpp[[1]][,i], out.mpp[[2]][,i], out.mpp[[3]][,i] )))
summary(out.mpp)
plot(out.mpp, ask = TRUE)

# save(out.mpp, file = "../output/fisher_mpp.Rda")
load("../output/fisher_mpp.Rda")
plot(out.mpp[,c("D", "sigma")])
dev.new()
plot(out.mpp.alg1[,c("D", "sigma")])
dev.new()
plot(out.sc[,c("D", "sigma")])


# Now add the SPIM sampler
# and run it again!
###############################
# Need to initialize this model as the stochastic node for ID is kind of wrong...
initsSPIM <- function(){
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 1, 5)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	ID <- numeric(length(omega))
	for(i in 1:length(omega))
	{
		cantMatch <- 1-(1:M %in% ID[cannotlink[i,] == 1])*1
		ID[i] <- sample(1:M, 1, prob = hkj[,omega[i]]*cantMatch)
	}
	# Now must link ID:
	indx <- which(rowSums(mustlink) > 1)
	ID[indx] <- ID[indx[1]]
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

MPPModel <- nimbleModel(SC_MPP, constants.mpp, data.mpp, inits = initsSPIM())
MPPconf <- configureMCMC(MPPModel)
MPPconf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID'))

MPPconf$removeSamplers('sigma')
MPPconf$addSampler(target = 'sigma', 
	type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))

# Use a block update on locations. Saves time.
MPPconf$removeSamplers('X')
for(i in 1:M) MPPconf$addSampler(target = paste0('X[', i, ', 1:2]'), 
	type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 1.7))
# Optimized z sampler
MPPconf$removeSamplers('z')
MPPconf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

MPPconf$removeSamplers('ID')
# Only 1 must link:
indx <- which(rowSums(mustlink) > 1)
MPPconf$addSampler(target = 'ID', type = 'mySPIM', scalarComponents = TRUE, control = list(M = M, cannotlink = cannotlink))
MPPconf$removeSamplers(paste0('ID[', indx, ']'))
MPPconf$addSampler(target = paste0('ID[', indx, ']'), type = 'mySPIM', scalarComponents = FALSE, control = list(M = M, cannotlink = cannotlink))
# MPPconf$printSamplers('ID')
MPPRmcmc <- buildMCMC(MPPconf)
MPPCmodel <- compileNimble(MPPModel)
MPPCmcmc <- compileNimble(MPPRmcmc, project = MPPModel)

MPPCmcmc$run(10000)
mvSamples <- MPPCmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:5000),])
plot(out[,c('N', 'D', 'sigma', 'lambda')])
post.id <- samples[-(1:5000),grep("ID", colnames(samples))]
NActiveCollaredMale <- apply(post.id[,obs$collar == 1 & obs$sex == 1], 1, FUN = function(x){ length(unique(x))})
NActiveCollaredFemale <- apply(post.id[,obs$collar == 1 & obs$sex == 2], 1, FUN = function(x){ length(unique(x))})
NActiveMales <- apply(post.id[,obs$sex == 1], 1, FUN = function(x){ length(unique(x))})
NActiveFemales <- apply(post.id[,obs$sex == 2], 1, FUN = function(x){ length(unique(x))})
par(mfrow = c(2,1))
hist(NActiveCollaredMale)
hist(NActiveCollaredFemale)
sum(obs$collar == 1 & obs$sex == 1)
sum(obs$collar == 1 & obs$sex == 2)

hist(NActiveMales)
hist(NActiveFemales)


# Doesn't look great with adaptive sampling of 0.028...
# Let's tune with 0.1
# valueInCompiledNimbleFunction(MPPCmcmc$samplerFunctions[[3]], "scale")

# post.id <- samples[-(1:5000),grep("ID", colnames(samples))]
# NActive <- apply(post.id, 1, FUN = function(x){ length(unique(x))})
# plot(NActive)
# hist(NActive)

samples.spim <- runMCMC(MPPCmcmc, 30000, nburnin = 10000, nchains = 3, 
	thin = 1, inits = list(initsSPIM(), initsSPIM(), initsSPIM()))

post.id.1 <- samples.spim[[1]][,grep("ID", colnames(samples.spim[[1]]))]
post.id.2 <- samples.spim[[2]][,grep("ID", colnames(samples.spim[[2]]))]
post.id.3 <- samples.spim[[3]][,grep("ID", colnames(samples.spim[[3]]))]
NActive1 <- apply(post.id.1, 1, FUN = function(x){ length(unique(x))})
NActive2 <- apply(post.id.2, 1, FUN = function(x){ length(unique(x))})
NActive3 <- apply(post.id.3, 1, FUN = function(x){ length(unique(x))})

out.spim <- mcmc.list(list(as.mcmc(cbind(samples.spim[[1]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "NActive" = NActive1)), 
	as.mcmc(cbind(samples.spim[[2]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "NActive" =  NActive2)),
	as.mcmc(cbind(samples.spim[[3]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "NActive" = NActive3))))
plot(out.spim[, c('D', 'NActive', 'sigma')])
summary(out.spim)
effectiveSize(out.spim)
stat <- NULL
for( i in c('sigma', 'lambda', 'psi', 'N', 'D', 'NActive')) stat[[i]] = gelman.diag(mcmc.list(list(out.spim[[1]][,i], out.spim[[2]][,i], out.spim[[3]][,i] )))
save(out.spim, file =  "../output/fisher_spim.Rda")
load("../output/fisher_spim.Rda")
plot(out.spim[,c('NActive')])
plot(out.spim[,c('D')])
dev.new()
plot(out.mpp[,c('NActive')])
plot(out.spim[,c('D', 'sigma')])


# MPPCmcmc$run(10000)
# mvSamples <- MPPCmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples[-(1:5000),])
plot(out[,c('N', 'D', 'sigma')])
dev.new()
plot(out.mpp[,c('N', 'D', 'sigma')])
dev.new()
plot(out.spim[,c('N', 'D')])
plot(out.spim[,c('sigma', 'lambda')])


# Run the same model from van Dam-Bates et al. Marked Poisson Process on fisher.
#----------------------------------------------------------------------
data.mpp.sex <- list(
    zones = rep(1, M),
    ones = rep(1, length(omega)),
	z =  rep(NA, M),
	sex = rep(NA, M),	
	ID = rep(NA, length(omega))
)

initsSPIMSex <- function(){
	lambda <- runif(1, 0.1, 1)
	psex <- rbeta(1,1,1)
	sigma <- runif(1, 1, 5)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	ID <- numeric(length(omega))
	for(i in 1:length(omega))
	{
		cantMatch <- 1-(1:M %in% ID[cannotlink[i,] == 1])*1
		ID[i] <- sample(1:M, 1, prob = hkj[,omega[i]]*cantMatch)
	}
	# Now must link ID:
	indx <- which(rowSums(mustlink) > 1)
	ID[indx] <- ID[indx[1]]
	z <- rep(0, M)
	z[ID] <- 1
	sex <- do.call('c', lapply(1:M, FUN = function(x){obs$sex[ID == x][1]})) - 1
	sex[is.na(sex)] <- -1
	sex[sex == -1] <- rbinom(n = sum(sex == -1), size = 1, prob = psex)
	psi <- length(unique(ID))/M
	list(
		lambda = lambda,
		sigma = c(sigma, sigma),
		psi = psi,
		X = X,
		z = z,
		ID = ID,
		sex = sex,
		psex = psex
    )
}

MPPModel <- nimbleModel(SC_MPP_sex, constants.mpp, data.mpp.sex, inits = initsSPIMSex())
MPPconf <- configureMCMC(MPPModel)
MPPconf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'psex', 'sex'))

# MPPconf$removeSamplers('sigma')
# MPPconf$addSampler(target = 'sigma', scalarComponents = TRUE,
	# type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))

MPPconf$removeSamplers('sex')
MPPconf$addSampler(target = 'sex',  scalarComponents = TRUE,
	type = 'myBinarySex', silent = TRUE, control = list(obsSex = obs$sex))

# Use a block update on locations. Saves time.
MPPconf$removeSamplers('X')
for(i in 1:M) MPPconf$addSampler(target = paste0('X[', i, ', 1:2]'), 
	type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 1.7))
# Optimized z sampler
MPPconf$removeSamplers('z')
MPPconf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

MPPconf$removeSamplers('ID')
# Only 1 must link:
indx <- which(rowSums(mustlink) > 1)
MPPconf$addSampler(target = 'ID', type = 'mySPIMSex', scalarComponents = TRUE, control = list(M = M, cannotlink = cannotlink, obsSex = obs$sex))
MPPconf$removeSamplers(paste0('ID[', indx, ']'))
MPPconf$addSampler(target = paste0('ID[', indx, ']'), type = 'mySPIMSex', scalarComponents = FALSE, control = list(M = M, cannotlink = cannotlink,  obsSex = obs$sex))
# MPPconf$printSamplers('ID')
MPPRmcmc <- buildMCMC(MPPconf)
MPPCmodel <- compileNimble(MPPModel)
MPPCmcmc <- compileNimble(MPPRmcmc, project = MPPModel)

MPPCmcmc$run(10000)
mvSamples <- MPPCmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:5000),])
plot(out[,c('N', 'D', 'psex')])
plot(out[,c('sigma[1]', 'sigma[2]', 'lambda')])

post.id <- samples[-(1:5000),grep("ID", colnames(samples))]
NActiveCollaredMale <- apply(post.id[,obs$collar == 1 & obs$sex == 1], 1, FUN = function(x){ length(unique(x))})
NActiveCollaredFemale <- apply(post.id[,obs$collar == 1 & obs$sex == 2], 1, FUN = function(x){ length(unique(x))})
NActiveMales <- apply(post.id[,obs$sex == 1], 1, FUN = function(x){ length(unique(x))})
NActiveFemales <- apply(post.id[,obs$sex == 2], 1, FUN = function(x){ length(unique(x))})
par(mfrow = c(2,1))
hist(NActiveCollaredMale)
hist(NActiveCollaredFemale)
sum(obs$collar == 1 & obs$sex == 1)
sum(obs$collar == 1 & obs$sex == 2)
NActive<- apply(post.id, 1, FUN = function(x){ length(unique(x))})
par(mfrow = c(3,1))
hist(NActive)
hist(NActiveMales)
hist(NActiveFemales)


####################################
# Sex and collar as fixed covariates
####################################
# van Dam-Bates Spatial Count Model using the Marked Poisson process formulation.
SC_MPP_sex_collar <- nimbleCode({
    lambda[1] ~ dunif(0, 20) # Detection rate at distance 0
    lambda[2] ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma[1] ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    sigma[2] ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    tau2[1] <- 1/(2*sigma[1]^2)	# Just avoid that extra computation for each animal...
    tau2[2] <- 1/(2*sigma[2]^2)	# Just avoid that extra computation for each animal...
	psex ~ dbeta(1,1)
	# For the collared individuals, we know the sex too! It's actually observed :)
	for(k in 1:n_collar) {
		z[k] ~ dbern(psi)
		sex[k] ~ dbern(psex)
		X[k, 1] ~ dunif(xlim[1], xlim[2])
		X[k, 2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
		pkj[k,1:J] <- exp(-d2[k,1:J]*tau2[sex[k]+1])*lambda[sex[k]+1]
		# Hazard rate for animal across all traps.
		Hk[k] <- sum(pkj[k,1:J])*Time
		pkz[k] <- exp(-Hk[k]*z[k])	# Only put z here for purposes of node dependence and speed.
		zones[k] ~ dbern(pkz[k])
	}
	for(k in (n_collar + 1):M) {
		z[k] ~ dbern(psi)
		sex[k] ~ dbern(psex)
		X[k, 1] ~ dunif(xlim[1], xlim[2])
		X[k, 2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
		pkj[k,1:J] <- exp(-d2[k,1:J]*tau2[sex[k]+1])*lambda[sex[k]+1]
		# Hazard rate for animal across all traps.
		Hk[k] <- sum(pkj[k,1:J])*Time
		pkz[k] <- exp(-Hk[k]*z[k])	# Only put z here for purposes of node dependence and speed.
		zones[k] ~ dbern(pkz[k])
	}

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # trap probability given ID:
        # This one can certainly be just the ones trick for trap y[i].
		pSex[i] <- (sex[ID[i]]+1) == obsSex[i] | obsSex[i] == 0
		pCollar[i] <- (noCollar[ID[i]]) == obsCollar[i] | obsCollar[i] == 0
		pobs[i] <- pkj[ID[i], omega[i]]*pSex[i]*pCollar[i]
        ones[i] ~ dbern(pobs[i])
		ID[i] ~ dID(lam = 1)	# Dummy distribution to declare this as stochastic and mulitply by lambda.
    }
	
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area
})


constants.mpp.cs <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
    n_obs = length(omega),
	omega = omega,
	area = diff(xlim)*diff(ylim)/100,
	n_collar = 14,
	obsCollar = obs$collar, # Unknown = 0, Collared = 1, Uncollared = 2
	obsSex = obs$sex
)

data.mpp.cs <- list(
    zones = rep(1, M),
    ones = rep(1, length(omega)),
	z =  c(rep(1, 14), rep(NA, M-14)),
	ID = rep(NA, length(omega)),
	sex = c(rep(0, 5), rep(1, 9), rep(NA, M-14)),	# Note 0 is male and 1 is female.
	noCollar = c(rep(1, 14), rep(2, M-14))	# No collar = 1. We observed who is collared for animals.
)

initsSexCollar <- function(){
	lambda <- runif(1, 0.1, 1)
	psex <- rbeta(1, 14, 9)	# Based on collared...
	sigma <- runif(1, 1, 5)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	sexCollar <- c(rep(0, 5), rep(1, 9))
	sex <- c(rep(NA, 14), rbinom(M-14, size = 1, prob = psex))
	ID <- numeric(length(omega))
	ID[obs$collar == 1] <- do.call('c',lapply(which(obs$collar == 1), FUN = function(x){sample(1:14, 1, prob = hkj[1:14,omega[x]]*(sexCollar+1 == obs$sex[x] |  obs$sex[x] == 0))}))
	ID[obs$collar == 2] <- do.call('c',lapply(which(obs$collar == 2), FUN = function(x){sample(15:M, 1, prob = hkj[15:M,omega[x]]*(sex[15:M]+1 == obs$sex[x] |  obs$sex[x] == 0))}))
	ID[obs$collar == 0] <- do.call('c',lapply(which(obs$collar == 0), FUN = function(x){sample(1:M, 1,  prob = hkj[1:M,omega[x]]*(c(sexCollar, sex[15:M])+1 == obs$sex[x] |  obs$sex[x] == 0))}))	
	z <- c(rep(NA, 14), rep(0, M-14))
	z[ID[ID > 14]] <- 1
	psi <- rbeta(1, sum(z,na.rm = TRUE) + 14, M - sum(1-z, na.rm = TRUE))	# NA inits...
	list(
		lambda = c(lambda, lambda),
		sigma = c(sigma, sigma),
		psi = psi,
		X = X,
		z = z,
		ID = ID,
		sex = sex,
		psex = psex
    )
}

MPPModel <- nimbleModel(SC_MPP_sex_collar, constants.mpp.cs, data.mpp.cs, inits = initsSexCollar())
MPPconf <- configureMCMC(MPPModel)
MPPconf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'psex', 'sex'))

MPPconf$removeSamplers('sigma')
MPPconf$addSampler(target = 'sigma', scalarComponents = TRUE,
	type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.25))

# Use a block update on locations. Saves time.
MPPconf$removeSamplers('X')
for(i in 1:M) MPPconf$addSampler(target = paste0('X[', i, ', 1:2]'), 
	type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 1.7))
# Optimized z sampler
MPPconf$removeSamplers('z')
MPPconf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

MPPconf$removeSamplers('ID')
MPPconf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
MPPRmcmc <- buildMCMC(MPPconf)
MPPCmodel <- compileNimble(MPPModel)
MPPCmcmc <- compileNimble(MPPRmcmc, project = MPPModel)

MPPCmcmc$run(10000)
mvSamples <- MPPCmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:5000),])
plot(out[,c('N', 'D', 'psex')])
dev.new()
plot(out[,c('sigma[1]', 'sigma[2]', 'lambda[1]', 'lambda[2]')])

post.id <- samples[-(1:5000),grep("ID", colnames(samples))]
post.sex <- samples[-(1:5000),grepl("sex", colnames(samples)) & !grepl("psex", colnames(samples))]
collaredMale <- apply(post.id, 1, FUN = function(x){ sum(1:5 %in% x)})
collaredFemale <- apply(post.id, 1, FUN = function(x){ sum(6:14 %in% x)})

table(post.sex[cbind(1:5000, post.id[,33])])

par(mfrow = c(2,1))
hist(collaredMale)
abline(v = sum(obs$collar == 1 & obs$sex == 1), col = 'red')
hist(collaredFemale)
abline(v = sum(obs$collar == 1 & obs$sex == 2), col = 'red')

NActiveCollaredMale <- apply(post.id[,obs$collar == 2 & obs$sex == 2], 1, FUN = function(x){ unique(sex())})
NActiveCollaredFemale <- apply(post.id[,obs$collar == 1 & obs$sex == 2], 1, FUN = function(x){ length(unique(x))})
NActiveMales <- apply(post.id[,obs$sex == 1], 1, FUN = function(x){ length(unique(x))})
NActiveFemales <- apply(post.id[,obs$sex == 2], 1, FUN = function(x){ length(unique(x))})
par(mfrow = c(2,1))
hist(NActiveCollaredMale)
hist(NActiveCollaredFemale)
sum(obs$collar == 1 & obs$sex == 1)
sum(obs$collar == 1 & obs$sex == 2)
NActive<- apply(post.id, 1, FUN = function(x){ length(unique(x))})
par(mfrow = c(3,1))
hist(NActive)
abline(v = 24, col = 'red')
hist(NActiveMales)
hist(NActiveFemales)


samples.sexcol <- runMCMC(MPPCmcmc, 30000, nburnin = 10000, nchains = 3, 
	thin = 1, inits = list(initsSexCollar(), initsSexCollar(), initsSexCollar()))

post.id.1 <- samples.sexcol[[1]][,grep("ID", colnames(samples.sexcol[[1]]))]
post.id.2 <- samples.sexcol[[2]][,grep("ID", colnames(samples.sexcol[[2]]))]
post.id.3 <- samples.sexcol[[3]][,grep("ID", colnames(samples.sexcol[[3]]))]
NActive1 <- apply(post.id.1, 1, FUN = function(x){ length(unique(x))})
NActive2 <- apply(post.id.2, 1, FUN = function(x){ length(unique(x))})
NActive3 <- apply(post.id.3, 1, FUN = function(x){ length(unique(x))})

out.sexcol <- mcmc.list(list(as.mcmc(cbind(samples.sexcol[[1]][,c('sigma[1]', 'lambda[1]', 'sigma[2]', 'lambda[2]', 'psi', 'N', 'D', 'psex')], "NActive" = NActive1)), 
	as.mcmc(cbind(samples.sexcol[[2]][,c('sigma[1]', 'lambda[1]', 'sigma[2]', 'lambda[2]', 'psi', 'N', 'D', 'psex')], "NActive" =  NActive2)),
	as.mcmc(cbind(samples.sexcol[[3]][,c('sigma[1]', 'lambda[1]', 'sigma[2]', 'lambda[2]', 'psi', 'N', 'D', 'psex')], "NActive" = NActive3))))

plot(out.sexcol[,1:4])
plot(out.sexcol[,5:8])















#########################################
# Incorporating the Marks sex and collar
#########################################
SC_MPP2 <- nimbleCode({
    lambda ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)	# Just avoid that extra computation for each animal...
	psex <- dbeta(1, 1)
    for(k in 1:M) {
        z[k] ~ dbern(psi)
		sex[k] ~ dbern(psex)
        X[k, 1] ~ dunif(xlim[1], xlim[2])
        X[k, 2] ~ dunif(ylim[1], ylim[2])
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
        pkj[k,1:J] <- exp(-d2[k,1:J]*tau2)
        # Hazard rate for animal across all traps.
        Hk[k] <- sum(pkj[k,1:J])*Time*lambda
		pkz[k] <- exp(-Hk[k]*z[k])	# Only put z here for purposes of node dependence and speed.
		zones[k] ~ dbern(pkz[k])
   }

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # trap probability given ID:
        # This one can certainly be just the ones trick for trap y[i].
		obsex ~ dIndicator(indicator = sex[ID[i]])
		pobs[i] <- pkj[ID[i], omega[i]]
        ones[i] ~ dbern(pobs[i])
		ID[i] ~ dID(lam = lambda)	# Dummy distribution to declare this as stochastic and mulitply by lambda.
    }
	
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area
})


