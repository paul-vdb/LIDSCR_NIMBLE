#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################
#################################
library(coda)
library(nimble)
library(nimbleSCR)

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

source("../Functions/SimData.R")
load("../../data/FisherData.Rda")
traps <- fisher.data$traps
xlim <- range(traps[,1]) + c(-6,6)
ylim <- range(traps[,2]) + c(-6,6)
J <- nrow(traps)
M <- 300

results <- foreach(h=1:250,
			.packages = c("coda", "nimbleSCR", "nimble"))%dopar%{
	source("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/Functions/NimbleFunctions.R")

	### Scenario 1:
	dat <- simSCR(N = 60, NCollar = c(5, 9), sigma = 1.5, lambda = 0.15,
		StudyPeriod = 64, traps = traps, 
		limits = list(xlim = xlim, ylim = ylim), psex = 0.6)
	omega <- dat$trap
	n <- nrow(dat)

	Model23 <- nimbleCode({
		lambda ~ dunif(0, 20) # Detection rate at distance 0
		psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
		sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
		tau2 <- 1/(2*sigma^2)	# Just avoid that extra computation for each animal...
		for(k in 1:M) {
			z[k] ~ dbern(psi)
			X[k, 1] ~ dunif(xlim[1], xlim[2])
			X[k, 2] ~ dunif(ylim[1], ylim[2])
			d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
			pkj[k,1:J] <- exp(-d2[k,1:J]*tau2)*lambda
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
			omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i], pID = 1)
			ID[i] ~ dID(lam = 1)	# Dummy distribution to declare this as stochastic and mulitply by lambda.
		}
		
		# Predicted population size
		N <- sum(z[1:M])
		D <- N/area
	})

	# Run the same model from van Dam-Bates et al. Marked Poisson Process on fisher.
	#----------------------------------------------------------------------
	constants23 <- list(
		J = J,
		xlim = xlim,
		ylim = ylim,
		traps = traps, 
		Time = 64,
		M = M,
		n_obs = length(omega),
		area = diff(xlim)*diff(ylim)/100
		)

	data23 <- list(
		zones = rep(1, M),
		z =  rep(NA, M),
		ID = rep(NA, length(omega)),
		omega = omega	
	)

	# Need to initialize this model as the stochastic node for ID is kind of wrong...
	init23 <- function(){
		lambda <- runif(1, 0.1, 1)
		sigma <- runif(1, 1, 2)
		X <- cbind(runif(M, xlim[1], xlim[2]), 
				  runif(M, ylim[1], ylim[2]))
		d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
		hkj <- exp(-d2/(2*sigma^2))
		psi <- runif(1,0.15,0.35)
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
	
	Rmodel <- nimbleModel(Model23, constants23, data23, inits = init23())
	conf <- configureMCMC(Rmodel)
	conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))
	# Use a block update on locations. Saves time.
	# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
	conf$removeSamplers('X')
	for(i in 1:M){
		conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
			type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 1.5))
		}

	# Optimized z sampler
	conf$removeSamplers('z')
	conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)
	# van Dam-Bates categorical sampler
	conf$removeSamplers('ID')
	# Chandler and Royle Alg. 1 sampler.
	# conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
	# van Dam-Bates Alg.
	conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))

	Rmcmc <- buildMCMC(conf)
	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

	Cmcmc$run(60000)
	mvSamples <- Cmcmc$mvSamples
	samples <- as.matrix(mvSamples)
	out <- mcmc(samples[-(1:30000),])
	# plot(out[,c("lambda", "sigma","N", "D")])
	save(out, file = paste0("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/FisherSimulations/SCR_Scenario_1IDZ_iter_", h, ".Rda"))
	summary(out)
}

stopCluster(cl)

# Adding sex + collar:
### Scenario 1:
# Make the sex and collar scenario now.
dat$sex <- ifelse(rbinom(n, 1, 0.32) == 1, dat$sex + 1, 0)
dat$collar <-  ifelse(rbinom(n, 1, 0.40) == 1, (1-dat$collar) + 1, 0)

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

	# Trap history model.
	# and unobserved animal ID.
	for(i in 1:n_obs) {
		# trap probability given ID:
		# This one can certainly be just the ones trick for trap y[i].
		pSex[i] <- (sex[ID[i]]+1) == obsSex[i] | obsSex[i] == 0
		pCollar[i] <- (noCollar[ID[i]]) == obsCollar[i] | obsCollar[i] == 0
		pID[i] <- pCollar[i]*pSex[i]
		omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i], pID = pID[i])
		ID[i] ~ dID(lam = lambda)	# Dummy distribution to declare this as stochastic and mulitply by lambda.
	}
	
	# Predicted population size
	N <- sum(z[1:M])
	D <- N/area
})


constants4 <- list(
	J = J,
	xlim = xlim,
	ylim = ylim,
	traps = traps, 
	Time = 64,
	M = M,
	n_obs = length(omega),
	area = diff(xlim)*diff(ylim)/100,
	n_collar = 14,
	obsCollar = dat$collar, # Unknown = 0, Collared = 1, Uncollared = 2
	obsSex = dat$sex
)

data4 <- list(
	zones = rep(1, M),
	omega = omega,
	z =  c(rep(1, 14), rep(NA, M-14)),
	ID = rep(NA, length(omega)),
	sex = c(rep(0, 5), rep(1, 9), rep(NA, M-14)),	# Note 0 is male and 1 is female.
	noCollar = c(rep(1, 14), rep(2, M-14))	# No collar = 1. We observed who is collared for animals.
)

# Need to fully initialize because I haven't implemented a random generator distribution for dID()...
# Bit of a pain to make sure we match sex and collar correctly.
init4 <- function(){
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
	ID[dat$collar == 1] <- do.call('c',lapply(which(dat$collar == 1), FUN = function(x){sample(1:14, 1, prob = hkj[1:14,omega[x]]*(sexCollar+1 == dat$sex[x] |  dat$sex[x] == 0))}))
	ID[dat$collar == 2] <- do.call('c',lapply(which(dat$collar == 2), FUN = function(x){sample(15:M, 1, prob = hkj[15:M,omega[x]]*(sex[15:M]+1 == dat$sex[x] |  dat$sex[x] == 0))}))
	ID[dat$collar == 0] <- do.call('c',lapply(which(dat$collar == 0), FUN = function(x){sample(1:M, 1,  prob = hkj[1:M,omega[x]]*(c(sexCollar, sex[15:M])+1 == dat$sex[x] |  dat$sex[x] == 0))}))	
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

Rmodel <- nimbleModel(Model4, constants4, data4, inits = init4())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D','psex'))
# Use a block update on locations. Saves time.
# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 1.5))
	}

# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)
# van Dam-Bates categorical sampler
conf$removeSamplers('ID')
# Chandler and Royle Alg. 1 sampler.
conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
# van Dam-Bates Alg.
# conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(60000)
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:20000),])
# plot(out[,c("lambda", "sigma","N", "D", "psex")])
save(out, file = paste0("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/FisherSimulations/SCR_Scenario_1SexCollar_iter_", h, ".Rda"))




# Do SCR with known ID to compare:
results <- foreach(h=1:100,
			.packages = c("coda", "nimbleSCR", "nimble"))%dopar%{
	source("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/Functions/NimbleFunctions.R")

	### Scenario 1:
	dat <- simSCR(N = 60, NCollar = c(5, 9), sigma = 1.5, lambda = 0.15,
		StudyPeriod = 64, traps = traps, 
		limits = list(xlim = xlim, ylim = ylim), psex = 0.6)
	omega <- dat$trap
	n <- nrow(dat)

	Model <- nimbleCode({
		lambda ~ dunif(0, 20) # Detection rate at distance 0
		psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
		sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
		tau2 <- 1/(2*sigma^2)	# Just avoid that extra computation for each animal...
		for(k in 1:M) {
			z[k] ~ dbern(psi)
			X[k, 1] ~ dunif(xlim[1], xlim[2])
			X[k, 2] ~ dunif(ylim[1], ylim[2])
			d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
			pkj[k,1:J] <- exp(-d2[k,1:J]*tau2)*lambda
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
			omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i], pID = 1)
			ID[i] ~ dID(lam = 1)
		}
		
		# Predicted population size
		N <- sum(z[1:M])
		D <- N/area
	})

	# Run the same model from van Dam-Bates et al. Marked Poisson Process on fisher.
	#----------------------------------------------------------------------
	ID <- dat$ID
	K <- max(ID)
	
	constants <- list(
		J = J,
		xlim = xlim,
		ylim = ylim,
		traps = traps, 
		Time = 64,
		M = M,
		n_obs = length(omega),
		area = diff(xlim)*diff(ylim)/100
		)

	data <- list(
		zones = rep(1, M),
		z =  c(rep(1, K), rep(NA, M-K)),
		omega = omega,
		ID = ID		
	)

	# Need to initialize this model as the stochastic node for ID is kind of wrong...
	inits <- function(){
		lambda <- runif(1, 0.1, 1)
		sigma <- runif(1, 1, 2)
		X <- cbind(runif(M, xlim[1], xlim[2]), 
				  runif(M, ylim[1], ylim[2]))
		d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
		psi <- rbeta(1,1,1)
		z <- rbinom(M, size = 1, prob = psi)
		z[1:K] <- 1
		list(
			lambda = lambda,
			sigma = sigma,
			psi = psi,
			X = X,
			z = z
		)
	}
	
	Rmodel <- nimbleModel(Model, constants, data, inits = inits())
	conf <- configureMCMC(Rmodel)
	conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))
	# Use a block update on locations. Saves time.
	# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
	conf$removeSamplers('X')
	for(i in 1:M){
		conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
			type = 'RW_block', silent = TRUE)
		}

	# Optimized z sampler
	conf$removeSamplers('z')
	conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

	Rmcmc <- buildMCMC(conf)
	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

	Cmcmc$run(60000)
	mvSamples <- Cmcmc$mvSamples
	samples <- as.matrix(mvSamples)
	out <- mcmc(samples[-(1:20000),])
	# plot(out[,c("lambda", "sigma","N", "D")])
	save(out, file = paste0("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/FisherSimulations/SCR_Scenario_1KnownID_iter_", h, ".Rda"))
	summary(out)
}

stopCluster(cl)

