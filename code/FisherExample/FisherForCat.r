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
# edf <- read.csv("C:/Users/Paul/Documents/Advising/FisherDataProcessing/edf.csv", header = TRUE, sep = ",", check.names = FALSE) ## NEED CHECK.NAMES = FALSE otherwise 3d will not work
# tmp <-edf[!duplicated(edf$individualID_2016),]
# table(tmp$sex)
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

####################################
# Sex and collar as fixed covariates
####################################
# van Dam-Bates Spatial Count Model using the Marked Poisson process formulation.
SC_MPP_sex_collar <- nimbleCode({
    lambda ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    tau2<- 1/(2*sigma^2)	# Just avoid that extra computation for each animal...
	psex ~ dbeta(1,1)
	# For the collared individuals, we know the sex too! It's actually observed :)
	for(k in 1:n_collar) {
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
	for(k in (n_collar + 1):M) {
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
		pSex[i] <- (sex[ID[i]]+1) == obsSex[i] | obsSex[i] == 0
		pCollar[i] <- (noCollar[ID[i]]) == obsCollar[i] | obsCollar[i] == 0
		pobs[i] <- pkj[ID[i], omega[i]]*pSex[i]*pCollar[i]
        ones[i] ~ dbern(pobs[i])
		ID[i] ~ dID(lam = lambda)	# Dummy distribution to declare this as stochastic and mulitply by lambda.
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

# Need to fully initialize because I haven't implemented a random generator distribution for dID()...
# Bit of a pain to make sure we match sex and collar correctly.
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
		lambda = lambda,
		sigma = sigma,
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
	type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 2))
# Optimized z sampler
MPPconf$removeSamplers('z')
MPPconf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

MPPconf$removeSamplers('ID')
MPPconf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
MPPRmcmc <- buildMCMC(MPPconf)
MPPCmodel <- compileNimble(MPPModel)
MPPCmcmc <- compileNimble(MPPRmcmc, project = MPPModel)

# Do a test run on the MCMC object.
# MPPCmcmc$run(1000)
# mvSamples <- MPPCmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples[-(1:5000),])
# plot(out[,c('N', 'D', 'psex')])	# proportion of males/females matches SCR data. 0.6
# dev.new()
# plot(out[,c('sigma', 'lambda')])

# The full run...
samples.sexcol <- runMCMC(MPPCmcmc, 60000, nburnin = 20000, nchains = 3, 
	thin = 1, inits = list(initsSexCollar(), initsSexCollar(), initsSexCollar()))

# Compute the number of actually observed fisher to compare with SCR model in paper.
post.id.1 <- samples.sexcol[[1]][,grep("ID", colnames(samples.sexcol[[1]]))]
post.id.2 <- samples.sexcol[[2]][,grep("ID", colnames(samples.sexcol[[2]]))]
post.id.3 <- samples.sexcol[[3]][,grep("ID", colnames(samples.sexcol[[3]]))]
NActive1 <- apply(post.id.1, 1, FUN = function(x){ length(unique(x))})
NActive2 <- apply(post.id.2, 1, FUN = function(x){ length(unique(x))})
NActive3 <- apply(post.id.3, 1, FUN = function(x){ length(unique(x))})

# Do it by sex as well...
cols <- grepl("sex", colnames(samples.sexcol[[1]])) & !grepl("psex", colnames(samples.sexcol[[1]]))
post.sex1 <- samples.sexcol[[1]][,cols]
post.sex2 <- samples.sexcol[[2]][,cols]
post.sex3 <- samples.sexcol[[3]][,cols]
NActiveFemales1 <- do.call('c', lapply(1:nrow(post.id.1),  FUN = function(x){ sum(post.sex1[x,unique(post.id.1[x,])])}))
NActiveFemales2 <- do.call('c', lapply(1:nrow(post.id.1),  FUN = function(x){ sum(post.sex2[x,unique(post.id.2[x,])])}))
NActiveFemales3 <- do.call('c', lapply(1:nrow(post.id.1),  FUN = function(x){ sum(post.sex3[x,unique(post.id.3[x,])])}))


# Complicated way of adding those observed animals to the MCMC object... Apologies.
out.sexcol <- mcmc.list(list(as.mcmc(cbind(samples.sexcol[[1]][,c('sigma', 'lambda', 'psi', 'N', 'D', 'psex')], "NActive" = NActive1)), 
	as.mcmc(cbind(samples.sexcol[[2]][,c('sigma', 'lambda', 'psi', 'N', 'D', 'psex')], "NActive" =  NActive2)),
	as.mcmc(cbind(samples.sexcol[[3]][,c('sigma', 'lambda', 'psi', 'N', 'D', 'psex')], "NActive" = NActive3))))

save(out.sexcol, file = "../output/FisherSexCollar.Rda")

plot(out.sexcol[,1:2])
dev.new()
plot(out.sexcol[,4:5])
plot(out.sexcol[,"NActive"], main = "Number of active clusters")

par(mfrow = c(3,1))
hist(c(NActive1, NActive2, NActive3), xlab = "# Detected Fisher", main = "")
abline(v = 24, col = 'red')
hist(c(NActiveFemales1, NActiveFemales2, NActiveFemales3), xlab = "# Detected Females", main = "")
abline(v = 15, col = 'red')
hist(c(NActive1 - NActiveFemales1, NActive2 -  NActiveFemales2, NActive3 -  NActiveFemales3), xlab = "# Detected Males", main = "")
abline(v = 9, col = 'red')

####################################
# SPIM version but constrain to max 5 males, 9 females.
# collared animals.
####################################
# van Dam-Bates Spatial Count Model using the Marked Poisson process formulation.
SC_Spim <- nimbleCode({
    lambda ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    tau2<- 1/(2*sigma^2)	# Just avoid that extra computation for each animal...
	# For the collared individuals, we know the sex too! It's actually observed :)
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


constants.spim <- list(
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

data.spim <- list(
    zones = rep(1, M),
    ones = rep(1, length(omega)),
	z =  c(rep(1, 14), rep(NA, M-14)),
	ID = rep(NA, length(omega))
)

# Need to fully initialize because I haven't implemented a random generator distribution for dID()...
# Bit of a pain to make sure we match sex and collar correctly.
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
		lambda = lambda,
		sigma = sigma,
		psi = psi,
		X = X,
		z = z,
		ID = ID
    )
}

MPPModel <- nimbleModel(SC_Spim, constants.spim, data.spim, inits = initsSexCollar())
MPPconf <- configureMCMC(MPPModel)
MPPconf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID'))

MPPconf$removeSamplers('sigma')
MPPconf$addSampler(target = 'sigma', scalarComponents = TRUE,
	type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.25))

# Use a block update on locations. Saves time.
MPPconf$removeSamplers('X')
for(i in 1:M) MPPconf$addSampler(target = paste0('X[', i, ', 1:2]'), 
	type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 2))
# Optimized z sampler
MPPconf$removeSamplers('z')
MPPconf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

MPPconf$removeSamplers('ID')
# Only 1 must link:
indx <- which(rowSums(mustlink) > 1)
# This is a stupid thing to do...
animalLinks <- cbind(rep(1, M),  			   # No info
	c(rep(1, 5), rep(0, M-5)),   			   # Male collared
	c(rep(0, 5), rep(1, 14-5), rep(0, M-14)),  # Female collared
	c(rep(1, 14),rep(0, M-14)),                # No sex collared
	c(rep(0, 14), rep(1, M-14))) # UnCollared
	

for(i in 1:length(omega)){
	if(obs[i, 'collar'] == 1 & obs[i, 'sex'] == 1) linki = animalLinks[,2]
	if(obs[i, 'collar'] == 1 & obs[i, 'sex'] == 2) linki = animalLinks[,3]
	if(obs[i, 'collar'] == 2) linki = animalLinks[,5]
	if(obs[i, 'collar'] == 0 & obs[i, 'sex'] == 0) linki = animalLinks[,1] + animalLinks[,5]
	if(obs[i, 'collar'] == 0 & obs[i, 'sex'] == 1) linki = animalLinks[,2] + animalLinks[,5]
	if(obs[i, 'collar'] == 0 & obs[i, 'sex'] == 2) linki = animalLinks[,3] + animalLinks[,5]
	MPPconf$addSampler(target = paste0('ID[', i, ']'), type = 'mySPIM', control = list(M = M, cannotlink = cannotlink, AnimalLink = linki))
}
MPPconf$removeSamplers(paste0('ID[', indx, ']'))
MPPconf$addSampler(target = paste0('ID[', indx, ']'), type = 'mySPIM', scalarComponents = FALSE, control = list(M = M, cannotlink = cannotlink, AnimalLink = animalLinks[,2]))

MPPRmcmc <- buildMCMC(MPPconf)
MPPCmodel <- compileNimble(MPPModel)
MPPCmcmc <- compileNimble(MPPRmcmc, project = MPPModel)

# Do a test run on the MCMC object.
MPPCmcmc$run(10000)
mvSamples <- MPPCmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:5000),])
plot(out[,c('N', 'D')])	# proportion of males/females matches SCR data. 0.6
dev.new()
plot(out[,c('sigma', 'lambda')])

post.id <- samples[-(1:5000),grep("ID", colnames(samples))]
sum(post.id[,2] == post.id[,106])
sqrt(sum((traps[35,] - traps[34,])^2))
cannotlink[34,106]
cannotlink[106,34]

NActive <- apply(post.id, 1, FUN = function(x){ length(unique(x))})
hist(NActive)
