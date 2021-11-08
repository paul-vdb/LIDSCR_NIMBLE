#################################
# SCR as a Marked Poisson Process
# Fisher Model - Spatial Counts from Burgar et al.
#################################
setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/FisherExample")
library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)

source("../Functions/NimbleFunctions.R")

## Load the Fisher data
###--- set directory and load files
load("../../data/FisherData.Rda")

buffer <- 6
traps <- fisher.data[["traps"]]
xlim <- range(traps[,1])+c(-buffer,buffer)
ylim <- range(traps[,2])+c(-buffer,buffer)
area <- diff(xlim)*diff(ylim)/100
obs <- fisher.data[["observations"]]
cannotlink <- fisher.data$cannotlink
mustlink <- fisher.data$mustlink

omega <- obs$TrapNumber
StudyPeriod <- 64
studyArea <- diff(xlim)*diff(ylim)

M <- 400
J <- nrow(traps)

# Chandler and Royle Spatial Count model, Algorithm 2
# Repeat of Burgar et al. 2018 analysis for papper.
counts <- NULL
for(j in 1:J)
{
	counts <- c(counts, sum(obs$TrapNumber == j))
}

Model1 <- nimbleCode( {
	# Priors:
	sigma ~ dunif(0,50)
	lambda ~ dunif(0,10)
	psi ~ dbeta(1,1)
	tau2 <- 1/(2*sigma^2)
	
	for(k in 1:M){
		z[k] ~ dbern(psi)
		X[k,1] ~ dunif(xlim[1], xlim[2])
		X[k,2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- pow(X[k,1] - traps[1:J,1], 2) + pow(X[k,2] - traps[1:J,2], 2)
		hkj[k,1:J] <- lambda*exp(-d2[k,1:J] * tau2 )*z[k]
	}

	# Model for Observation Process:
	for(j in 1:J){
		Hj[j] <- sum(hkj[1:M,j])*Time
		counts[j] ~ dpois(Hj[j])
	}
	
	N <- sum(z[1:M])
	D <- N/area
} )

# Initialize model:
#-------------------------
init1 <- function(){
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 1, 5)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	psi <- rbeta(1,1,1)
	z <- rbinom(M, prob=psi, size = 1)
	z[1:nrow(traps)] <- 1
	X[1:nrow(traps),] <- as.matrix(traps) + cbind(rnorm(nrow(traps), 0, 0.1), rnorm(nrow(traps), 0, 0.1))
	list(
		lambda = lambda,
		sigma = sigma,
		psi = psi,
		X = X,
		z = z
    )
}


# Run the same model from Burgar et al. Spatial Count on fisher.
#----------------------------------------------------------------------
constants.1 <- list(
    J = nrow(traps),
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
	area = diff(xlim)*diff(ylim)/100
	)

data.1 <- list(
	z =  rep(NA, M),
	counts = counts
	)

Rmodel <- nimbleModel(Model1, constants.1, data.1, inits = init1())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))

# Use a block update on locations. Saves time.
conf$removeSamplers('X')
for(i in 1:M)	{
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'RW_block', silent = TRUE)	# Adaptive is okay for this model.
}
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Cmcmc$run(1000)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples)
# plot(out[,c('sigma', 'lambda', 'N', 'D')])

# Note that when doing adaptive sampling it was around scale ~= 1.7. We will set it at that for all chains.
samples1 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(init1(), init1(), init1()))

out1 <- mcmc.list(list(as.mcmc(samples1[[1]]), as.mcmc(samples1[[2]]), as.mcmc(samples1[[3]])))
save(out1, file = "../../output/fisherModel1.Rda")
# load("../../output/fisherModel1.Rda")
# plot(out1[,c("sigma", "lambda", "psi")])
# plot(out1[,c("sigma", "N", "D")])
# mean(diff(samples1[[1]][, 'sigma']) > 0)
# summary(out1)
# effectiveSize(out1)
stat <- NULL
for( i in c('sigma', 'lambda', 'psi', 'N', 'D')) stat[[i]] = gelman.diag(mcmc.list(list(out1[[1]][,i], out1[[2]][,i], out1[[3]][,i] )))
# load("../../output/fisherModel1.Rda")

# Marked Poisson process model incorporating the 
# latent ID variable as unknown. 
#------------------------------------------------
ModelMPP <- nimbleCode({
    lambda ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)	# Just avoid that extra computation for each animal...
    for(k in 1:M) {
        # Standard SCR data augmentation, animal in pop, loc and dist to traps.
		z[k] ~ dbern(psi)
        X[k, 1] ~ dunif(xlim[1], xlim[2])
        X[k, 2] ~ dunif(ylim[1], ylim[2])
        d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2

		#Now the rate for the Poisson process:
		# Don't put lambda here for efficiency of updates.
        pkj[k,1:J] <- exp(-d2[k,1:J]*tau2)
        # Hazard rate for animal across all traps.
        Hk[k] <- sum(pkj[k,1:J])*Time*lambda
		# Only put z here for purposes of node dependence. Bounded 0-1 for ones trick.
		pkz[k] <- exp(-Hk[k]*z[k])	
		zones[k] ~ dbern(pkz[k])
   }

	pID[1:M] <- z[1:M]*lambda	# Puts a prior on ID of z==1, adds lambda to the detections likelihood.

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # trap probability given ID:
        # This one can certainly be just the ones trick for trap y[i].
        omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i])
		ID[i] ~ dID(pID = pID[1:M])	# Dummy distribution to declare this as stochastic and mulitply by lambda.
    }
	
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area
})

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
	area = diff(xlim)*diff(ylim)/100
	)

data.mpp <- list(
    zones = rep(1, M),
	z =  rep(NA, M),
	ID = rep(NA, length(omega)),
	omega = omega	
)

# Need to initialize this model as the stochastic node for ID is kind of wrong...
inits.mpp <- function(){
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 1, 2)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	psi <- rbeta(1,1,1)
	z <- rbinom(M, size = 1, prob = psi)
	ID <- do.call('c', lapply(omega, FUN = function(x) {sample(1:M, 1, prob = z*hkj[,x])}))
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
Rmodel <- nimbleModel(ModelMPP, constants.mpp, data.mpp, inits = inits.mpp())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID'))
# Use a block update on locations. Saves time.
# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'myJAM', silent = TRUE, control = list(scale = 1.5, xlim = xlim, ylim = ylim, temp = 0.2))
}

conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))		

# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)
# van Dam-Bates categorical sampler
conf$removeSamplers('ID')
# Chandler and Royle Alg. 1 sampler.
conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Cmcmc$run(15000, time = TRUE)
# Cmcmc$getTimes()
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples[-(1:5000),])
# plot(out[,c('sigma', 'lambda', 'N', 'D')])
# mean(diff(out[,c('sigma')]) > 0)
# valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[403]], "width")

samples2 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(inits.mpp(), inits.mpp(), inits.mpp()))

post21 <- samples2[[1]][,grep("ID", colnames(samples2[[1]]))]
post22 <- samples2[[2]][,grep("ID", colnames(samples2[[2]]))]
post23 <- samples2[[3]][,grep("ID", colnames(samples2[[3]]))]
NObs1 <- apply(post21, 1, FUN = function(x){ length(unique(x))})
NObs2 <- apply(post22, 1, FUN = function(x){ length(unique(x))})
NObs3 <- apply(post23, 1, FUN = function(x){ length(unique(x))})

out2 <- mcmc.list(list(as.mcmc(cbind(samples2[[1]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "K" = NObs1)), 
	as.mcmc(cbind(samples2[[2]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "K" =  NObs2)),
	as.mcmc(cbind(samples2[[3]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "K" = NObs3))))
save(out2, file = "../../output/fisherModel2.Rda")	
# load("../../output/fisherModel2.Rda")	
# plot(out2[, c('sigma', 'D')])
# MCMCglmm::posterior.mode(out2)

# summary(out2)
# effectiveSize(out2)
# effectiveSize(out3)
# load("../../output/fisherModel3.Rda")	

# Plot 1 and 2 methods together:
out1D <- do.call("c", out1[,"D"])
out2D <- do.call("c", out2[,"D"])
out3D <- do.call("c", out4[,"D"])
plot(density(out2D), main = "Animal Density", ylim = c(0,1.5), xlab = "Density (Ind/100km2)")
lines(density(out1D), col = "blue")
lines(density(out3D), col = 'red')

####################################
# Model 3: Sex as covariates
####################################
Model3 <- nimbleCode({
    lambda ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    tau2<- 1/(2*sigma^2)	# Just avoid that extra computation for each animal...
	psex ~ dbeta(1,1)
	# For the collared individuals, we know the sex too! It's actually observed.
	# So we have observed sex 1:14 and an additional covariate of collar.
	# As a result, we also have z = 1 observed for 1:14
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

	pID[1:M] <- z[1:M]*lambda

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
		# Hard match/no match info.
		pSex[i] <- (sex[ID[i]]+1) == obsSex[i] | obsSex[i] == 0
		keep[i] ~ dbern(pSex[i])
		
		# Trap Prob
        omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i])

		# ID prior based on z, with lambda mixed in.
		ID[i] ~ dID(pID = pID[1:M])
    }
	
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area
})

constants3 <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
    n_obs = length(omega),
	area = diff(xlim)*diff(ylim)/100,
	obsSex = obs$sex
)

data3 <- list(
    zones = rep(1, M),
    omega = omega,
	keep = rep(1, length(omega)),
	z =  rep(NA,M),
	ID = rep(NA, length(omega)),
	sex = rep(NA, M)	# Note 0 is male and 1 is female.
)

# Need to fully initialize because I haven't implemented a random generator distribution for dID()...
# Bit of a pain to make sure we match sex and collar correctly.
init3 <- function(){
	lambda <- runif(1, 0.1, 1)
	psex <- rbeta(1, Nc, Ncf)	# Based on collared...
	sigma <- runif(1, 2, 4)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	sex <- rbinom(M, size = 1, prob = psex)
	ID <- numeric(length(omega))
	ID <- do.call('c',lapply(1:length(omega), FUN = function(x){sample(1:M, 1, prob = hkj[1:M,omega[x]]*(sex+1 == obs$sex[x] |  obs$sex[x] == 0))}))
	z <- rep(0, M)
	z[ID] <- 1
	psi <- rbeta(1, sum(z, na.rm = TRUE), M - sum(1-z, na.rm = TRUE))	# NA inits...
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

###################################
# Chandler and Royle Sampler:
###################################
Rmodel <- nimbleModel(Model3, constants3, data3, inits = init3())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'psex', 'sex', 'z'))

conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'myJAM', silent = TRUE, control = list(scale = 1.5, xlim = xlim, ylim = ylim, temp = 0.2))
	# conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		# type = 'RW_block', silent = TRUE, control = list(scale = 2, adaptive = FALSE))		
}
conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 1))		
# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

conf$removeSamplers('ID')
conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(10000, time = TRUE)
# Cmcmc$getTimes()
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
# valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[3]], "scale")
out <- mcmc(samples[-(1:5000),])
plot(out[,c('sigma', 'lambda', 'N', 'D', 'psex')])
post <- samples[-(1:5000),grep("ID", colnames(samples))]
K <- apply(post, 1, FUN = function(x){length(unique(x))})
plot(samples[-(1:5000),"N"], K)

# The full run...
samples3 <- runMCMC(Cmcmc, niter = 40000, nburnin = 20000, nchains = 3, 
	thin = 1, inits = list(init3(), init3(), init3()))

# Compute the number of actually observed fisher to compare with SCR model in paper.
post1 <- samples3[[1]][,grep("ID", colnames(samples3[[1]]))]
post2 <- samples3[[2]][,grep("ID", colnames(samples3[[2]]))]
post3 <- samples3[[3]][,grep("ID", colnames(samples3[[3]]))]
NObs1 <- apply(post1, 1, FUN = function(x){ length(unique(x))})
NObs2 <- apply(post2, 1, FUN = function(x){ length(unique(x))})
NObs3 <- apply(post3, 1, FUN = function(x){ length(unique(x))})

# Do it by sex as well...
cols <- grepl("sex", colnames(samples3[[1]])) & !grepl("psex", colnames(samples3[[1]]))
post.sex1 <- samples3[[1]][,cols]
post.sex2 <- samples3[[2]][,cols]
post.sex3 <- samples3[[3]][,cols]
NObsF1 <- do.call('c', lapply(1:nrow(post1),  FUN = function(x){ sum(post.sex1[x,unique(post1[x,])])}))
NObsF2 <- do.call('c', lapply(1:nrow(post2),  FUN = function(x){ sum(post.sex2[x,unique(post2[x,])])}))
NObsF3 <- do.call('c', lapply(1:nrow(post3),  FUN = function(x){ sum(post.sex3[x,unique(post3[x,])])}))

# Complicated way of adding those observed animals to the MCMC object... Apologies.
out3 <- mcmc.list(list(as.mcmc(cbind(samples3[[1]][,c('sigma', 'lambda', 'psi', 'N', 'D', 'psex')], "K" = NObs1, "KF" = NObsF1, "KM" = NObs1 - NObsF1)), 
	as.mcmc(cbind(samples3[[2]][,c('sigma', 'lambda', 'psi', 'N', 'D', 'psex')], "K" =  NObs2, "KF" = NObsF2, "KM" = NObs2 - NObsF2)),
	as.mcmc(cbind(samples3[[3]][,c('sigma', 'lambda', 'psi', 'N', 'D', 'psex')], "K" = NObs3, "KF" = NObsF3, "KM" = NObs3 - NObsF3))))
plot(out3[, c("sigma", "D", "K")])
MCMCglmm::posterior.mode(out3)




####################################
# Model 4: Sex and collar as covariates
####################################
Model4 <- nimbleCode({
    lambda ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
    tau2<- 1/(2*sigma^2)	# Just avoid that extra computation for each animal...
	psex ~ dbeta(1,1)
	# For the collared individuals, we know the sex too! It's actually observed.
	# So we have observed sex 1:14 and an additional covariate of collar.
	# As a result, we also have z = 1 observed for 1:14
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

	pID[1:M] <- z[1:M]*lambda

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
		# Hard match/no match info.
		pSex[i] <- (sex[ID[i]]+1) == obsSex[i] | obsSex[i] == 0
		pCollar[i] <- (noCollar[ID[i]]) == obsCollar[i] | obsCollar[i] == 0
		pMatch[i] <- pCollar[i]*pSex[i]
		keep[i] ~ dbern(pMatch[i])
		
		# Trap Prob
        omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i])

		# ID prior based on z, with lambda mixed in.
		ID[i] ~ dID(pID = pID[1:M])
    }
	
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area
})

Nc <- 14
Ncm <- 5
Ncf <- 9

constants4 <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
    n_obs = length(omega),
	area = diff(xlim)*diff(ylim)/100,
	n_collar = Nc,
	obsCollar = obs$collar, # Unknown = 0, Collared = 1, Uncollared = 2
	obsSex = obs$sex
)

data4 <- list(
    zones = rep(1, M),
    omega = omega,
	keep = rep(1, length(omega)),
	z =  c(rep(1, Nc), rep(NA, M-Nc)),
	ID = rep(NA, length(omega)),
	sex = c(rep(0, Ncm), rep(1, Ncf), rep(NA, M-Nc)),	# Note 0 is male and 1 is female.
	noCollar = c(rep(1, Nc), rep(2, M-Nc))	# No collar = 1. We observed who is collared for animals.
)

# Need to fully initialize because I haven't implemented a random generator distribution for dID()...
# Bit of a pain to make sure we match sex and collar correctly.
init4 <- function(){
	lambda <- runif(1, 0.1, 1)
	psex <- rbeta(1, Nc, Ncf)	# Based on collared...
	sigma <- runif(1, 2, 4)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	sexCollar <- c(rep(0, Ncm), rep(1, Ncf))
	sex <- c(rep(NA, Nc), rbinom(M-Nc, size = 1, prob = psex))
	ID <- numeric(length(omega))
	ID[obs$collar == 1] <- do.call('c',lapply(which(obs$collar == 1), FUN = function(x){sample(1:Nc, 1, prob = hkj[1:Nc,omega[x]]*(sexCollar+1 == obs$sex[x] |  obs$sex[x] == 0))}))
	ID[obs$collar == 2] <- do.call('c',lapply(which(obs$collar == 2), FUN = function(x){sample((Nc+1):M, 1, prob = hkj[(Nc+1):M,omega[x]]*(sex[(Nc+1):M]+1 == obs$sex[x] |  obs$sex[x] == 0))}))
	ID[obs$collar == 0] <- do.call('c',lapply(which(obs$collar == 0), FUN = function(x){sample(1:M, 1,  prob = hkj[1:M,omega[x]]*(c(sexCollar, sex[(Nc+1):M])+1 == obs$sex[x] |  obs$sex[x] == 0))}))	
	z <- rep(0, M)
	z[ID] <- 1
	z[1:Nc] <- NA
	psi <- rbeta(1, sum(z, na.rm = TRUE) + Nc, M - sum(1-z, na.rm = TRUE))	# NA inits...
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

###################################
# Chandler and Royle Sampler:
###################################
Rmodel <- nimbleModel(Model4, constants4, data4, inits = init4())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'psex', 'sex', 'z'))

conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'myJAM', silent = TRUE, control = list(scale = 1.5, xlim = xlim, ylim = ylim, temp = 0.2))
}

conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 1))		

# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

conf$removeSamplers('ID')
conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(20000, time = TRUE)
# Cmcmc$getTimes()
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:5000),])
plot(out[,c('sigma', 'lambda', 'N', 'D', 'psex')])
plot(out[,c('sigma', 'D')])

# The full run...
samples4 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(init4(), init4(), init4()))

# Compute the number of actually observed fisher to compare with SCR model in paper.
post1 <- samples4[[1]][,grep("ID", colnames(samples4[[1]]))]
post2 <- samples4[[2]][,grep("ID", colnames(samples4[[2]]))]
post3 <- samples4[[3]][,grep("ID", colnames(samples4[[3]]))]
NObs1 <- apply(post1, 1, FUN = function(x){ length(unique(x))})
NObs2 <- apply(post2, 1, FUN = function(x){ length(unique(x))})
NObs3 <- apply(post3, 1, FUN = function(x){ length(unique(x))})

# Do it by sex as well...
cols <- grepl("sex", colnames(samples4[[1]])) & !grepl("psex", colnames(samples4[[1]]))
post.sex1 <- samples4[[1]][,cols]
post.sex2 <- samples4[[2]][,cols]
post.sex3 <- samples4[[3]][,cols]
NObsF1 <- do.call('c', lapply(1:nrow(post1),  FUN = function(x){ sum(post.sex1[x,unique(post1[x,])])}))
NObsF2 <- do.call('c', lapply(1:nrow(post2),  FUN = function(x){ sum(post.sex2[x,unique(post2[x,])])}))
NObsF3 <- do.call('c', lapply(1:nrow(post3),  FUN = function(x){ sum(post.sex3[x,unique(post3[x,])])}))

# Complicated way of adding those observed animals to the MCMC object... Apologies.
out4 <- mcmc.list(list(as.mcmc(cbind(samples4[[1]][,c('sigma', 'lambda', 'psi', 'N', 'D', 'psex')], "K" = NObs1, "KF" = NObsF1, "KM" = NObs1 - NObsF1)), 
	as.mcmc(cbind(samples4[[2]][,c('sigma', 'lambda', 'psi', 'N', 'D', 'psex')], "K" =  NObs2, "KF" = NObsF2, "KM" = NObs2 - NObsF2)),
	as.mcmc(cbind(samples4[[3]][,c('sigma', 'lambda', 'psi', 'N', 'D', 'psex')], "K" = NObs3, "KF" = NObsF3, "KM" = NObs3 - NObsF3))))
save(out4, file = "../../output/fisherModel4.Rda")
# load("../../output/fisherModel4.Rda")
# effectiveSize(out4)
# load("../../output/fisherModel4cr.Rda")

# plot(out4[, c("N", "K", "KF", "KM")])
# plot(out4[, c("sigma", "lambda", "D")])

# par(mfrow = c(3,1))
# hist(c(NObs1, NObs2, NObs3), xlab = "# Observed Fisher", main =  "")
# hist(c(NObsF1, NObsF2, NObsF3), xlab = "# Observed Female Fisher", main = "")
# hist(c(NObs1 - NObsF1, NObs2 - NObsF2, NObs3 -NObsF3), xlab = "# Observed Male Fisher", main = "")


library(ggplot2)
library(reshape2)
# Visualize the data:
load("../../output/fisherModel1.Rda")
load("../../output/fisherModel2.Rda")
load("../../output/fisherModel4.Rda")

out1.df <- do.call("rbind", lapply(out1, as.data.frame))
out1.df$K <- NA
out1.df$Method = "SC"
out2.df <- do.call("rbind", lapply(out2, as.data.frame))
out2.df$Method = "MPP"
out4.df <- do.call("rbind", lapply(out4, as.data.frame))
out4.df$Method = "MPP w/ sex and collar"

facet_names <- c('sigma~(km)', 'lambda~(day^-1)', 'D~(Ind~100~km^-2)', 'Number~of~detected~animals')
names(facet_names) <- c("sigma", "lambda", "D", "K")
  
all.dat1 <- rbind(out1.df[, c("sigma", "lambda", "D", "K", "Method")], 
	out2.df[, c("sigma", "lambda", "D", "K", "Method")],
	out4.df[, c("sigma", "lambda", "D", "K", "Method")])
all.dat.l <- melt(all.dat1, id.vars = "Method")
all.dat.l <- all.dat.l[!(all.dat.l$variable == "K" & all.dat.l$Method == "SC"),]
all.dat.l$Method <- factor(all.dat.l$Method, levels = c("SC", "MPP", "MPP w/ sex and collar"))

ggplot(data = all.dat.l, aes(y = value, x = Method)) + 
	facet_wrap(~variable, scale = "free", labeller = as_labeller(x = facet_names, label_parsed)) + 
	theme_classic() +
	geom_boxplot(width = 0.1, outlier.alpha = 0, fill = "black", colour = "white") + 
	geom_violin(alpha = 0, adjust = 2) + 
	ylab("") + xlab("")