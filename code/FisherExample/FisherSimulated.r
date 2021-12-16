#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################
setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/FisherExample")
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

# Scenario 1:
results <- foreach(h=20:150,
			.packages = c("coda", "nimbleSCR", "nimble"))%dopar%{
	source("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/Functions/NimbleFunctions.R")

	### Scenario 1:
	dat <- simSCR(N = 40, NCollar = c(5, 9), sigma = 2, lambda = 0.1,
		StudyPeriod = 64, traps = traps, 
		limits = list(xlim = xlim, ylim = ylim), psex = 0.6)
	omega <- dat$trap
	n <- nrow(dat)

	# Model 1: LID
	LID <- nimbleCode({
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
			Hk[k] <- sum(pkj[k,1:J])*Time*lambda*z[k]
			zeros[k] ~ dpois(Hk[k])
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
	constants.lid <- list(
		J = J,
		xlim = xlim,
		ylim = ylim,
		traps = traps, 
		Time = 64,
		M = M,
		n_obs = length(omega),
		area = diff(xlim)*diff(ylim)/100
		)

	data.lid <- list(
		zeros = rep(0, M),
		z =  rep(NA, M),
		ID = rep(NA, length(omega)),
		omega = omega	
	)

	# Need to initialize this model as the stochastic node for ID is kind of wrong...
	inits.lid <- function(){
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
	Rmodel <- nimbleModel(LID, constants.lid, data.lid, inits = inits.lid())
	conf <- configureMCMC(Rmodel)
	conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))
	# Use a block update on locations. Saves time.
	# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
	conf$removeSamplers('X')
	for(i in 1:M){
		conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
			type = 'myJAM', silent = TRUE, 
			control = list(scale = 0.75, xlim = xlim, ylim = ylim, temp = 0.2))
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

	Cmcmc$run(60000)
	mvSamples <- Cmcmc$mvSamples
	samples <- as.matrix(mvSamples)
	out <- mcmc(samples[-(1:20000),])
	# plot(out[,c("lambda", "sigma","N", "D")])
	save(out, file = paste0("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/FisherSimulations/Scenario_2_LID_iter_", h, ".Rda"))
	tmp <- summary(out[,c("lambda", "sigma","N", "D")])
	res <- data.frame(cbind(tmp[[1]], tmp[[2]]))
	res$Method <- "LID"
	res$iter <- h

	##############################
	# Add Sex to the model:
	##############################
	LID_Sex <- nimbleCode({
		lambda ~ dunif(0, 20)

		# Home range parameter m=0/f=1:
		sigma ~ dunif(0, 50)

		# convert now to make it simpler in half-normal computation
		tau2 <- 1/(2*sigma^2)

		psi ~ dbeta(1, 1)      # prior on data augmentation Bernoulli vec.

		psex ~ dbeta(1, 1)
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
			Hk[k] <- sum(pkj[k,1:J])*Time*z[k]*lambda
			zeros[k] ~ dpois(Hk[k])
		}

		pID[1:M] <- z[1:M]*lambda

		# Trap history model.
		# and unobserved animal ID.
		for(i in 1:n_obs) {
			# Hard match/no match info.
			pSex[i] <- (sex[ID[i]]+1) == obsSex[i] | obsSex[i] == 0
			keep[i] ~ dbern(pSex[i])
			
			# Trap Prob.
			omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i])

			# ID prior based on z, with lambda mixed in.
			ID[i] ~ dID(pID = pID[1:M])
		}
		
		# Predicted population size
		N <- sum(z[1:M])
		D <- N/area
	})

	obsSex <- dat$sex
	obsSex <- ifelse(rbinom(n, 1, 0.32) == 1, obsSex + 1, 0)

	constants.sex <- list(
		J = J,
		xlim = xlim,
		ylim = ylim,
		traps = traps, 
		Time = 64,
		M = M,
		n_obs = length(omega),
		area = diff(xlim)*diff(ylim)/100,
		obsSex = obsSex
	)

	data.sex <- list(
		zeros = rep(0, M),
		omega = omega,
		keep = rep(1, length(omega)),
		z =  rep(NA, M),
		ID = rep(NA, length(omega)),
		sex = rep(NA, M)	# Note 0 is male and 1 is female.
	)

	# Need to fully initialize because I haven't implemented a random generator distribution for dID()...
	# Bit of a pain to make sure we match sex and collar correctly.
	init.sex <- function(){
		N <- floor(runif(1, 1, M/2))
		psi <- rbeta(1, N, M-N)	# NA inits...	
		lambda <- runif(1, 0.1, 1)
		psex <- rbeta(1, 1, 1)	# Based on collared...
		sigma <- runif(1, 2, 4)
		X <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
		d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
		hkj <- exp(-d2/(2*sigma^2))
		sex <- rbinom(M, size = 1, prob = psex)
		ID <- do.call('c',lapply(1:length(omega), FUN = function(x){sample(1:M, 1, prob = hkj[1:M,omega[x]]*(sex+1 == obsSex[x] |  obsSex[x] == 0))}))
		Hk <- rowSums(hkj)*64*lambda
		p <- exp(-Hk)*psi/(exp(-Hk)*psi + (1-psi))
		z <- rbinom(M, size = 1, prob=p)
		z[ID] <- 1
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
	Rmodel <- nimbleModel(LID_Sex, constants.sex, data.sex, inits = init.sex())
	conf <- configureMCMC(Rmodel)
	conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))

	conf$removeSamplers('X')
	for(i in 1:M){
		conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
			type = 'myJAM', silent = TRUE, 
			control = list(scale = 0.75, xlim = xlim, ylim = ylim, temp = 0.2))
	}

	conf$removeSamplers('sigma')
	conf$addSampler(target = 'sigma', 
			type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))	

	# Optimized z sampler
	conf$removeSamplers('z')
	conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

	conf$removeSamplers('ID')
	conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
	# conf$addSampler('ID', type = 'mySPIM', scalarComponents = TRUE, control = list(M = M, cannotlink = cannotlink))

	Rmcmc <- buildMCMC(conf)
	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

	Cmcmc$run(60000, time = TRUE)
	mvSamples <- Cmcmc$mvSamples
	samples <- as.matrix(mvSamples)
	out <- mcmc(samples[-(1:20000),])
	plot(out[,c("lambda", "sigma","N", "D")])
	save(out, file = paste0("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/FisherSimulations/Scenario_2_LID_Sex_iter_", h, ".Rda"))
	tmp <- summary(out[,c("lambda", "sigma","N", "D")])
	res2 <- data.frame(cbind(tmp[[1]], tmp[[2]]))
	res2$Method <- "LID_Sex"
	res2$iter <- h
	res <- rbind(res, res2)

	#######################
	# Add Sex and Collar
	#######################
	LID_Sex_collar <- nimbleCode({
		
		lambda ~ dunif(0, 20)

		# Home range parameter m=0/f=1:
		sigma ~ dunif(0, 50)	

		# convert now to make it simpler in half-normal computation
		tau2 <- 1/(2*sigma^2)	

		psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.

		psex ~ dbeta(1, 1)

		# For the collared individuals, we observe the sex too.
		# So we have observed sex 1:14 and z = 1 for those animals.
		for(k in 1:M) {
			z[k] ~ dbern(psi)
			sex[k] ~ dbern(psex)
			
			X[k, 1] ~ dunif(xlim[1], xlim[2])
			X[k, 2] ~ dunif(ylim[1], ylim[2])
			d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2

			pkj[k,1:J] <- exp(-d2[k,1:J]*tau2)*lambda
			# Hazard rate for animal across all traps.
			Hk[k] <- sum(pkj[k,1:J])*Time*z[k]

			# Zero count for all.
			zeros[k] ~ dpois(Hk[k])
		}

		# Prior on animal ID is just 1|z = 1 or 0|z = 0. 
		# The reason is in homogeneous time things cancel nicely and this is easiest.
		pID[1:M] <- z[1:M]

		# Trap history model.
		# and unobserved animal ID.
		for(i in 1:n_obs) {
			# Hard match/no match info.
			pSex[i] <- (sex[ID[i]]+1) == obsSex[i] | obsSex[i] == 0
			# We add matchCollar here to make sure collar matches.
			pCollar[i] <- (collar[ID[i]]+1) == obsCollar[i] | obsCollar[i] == 0
			pMatch[i] <- pSex[i]*pCollar[i]
			# Hard 1/0 for animal match to sex and collar.
			keep[i] ~ dbern(pMatch[i])
			
			# Trap Prob:
			omega[i] ~ dTrap(p = pkj[1:M, 1:J], ID = ID[i])

			# ID prior based on z:
			ID[i] ~ dID(pID = pID[1:M])
		}
		
		# Predicted population size
		N <- sum(z[1:M])
		D <- N/area
	})

	# Process collar information per animal.
	########################################
	obsCollar <-  ifelse(rbinom(n, 1, 0.40) == 1, (1-dat$collar) + 1, 0)

	constants.collar <- list(
		J = J,
		xlim = xlim,
		ylim = ylim,
		traps = traps, 
		Time = 64,
		M = M,
		n_obs = length(omega),
		area = diff(xlim)*diff(ylim)/100,
		obsSex = obsSex,
		obsCollar = obsCollar
		)

	data.collar <- list(
		zeros = rep(0, M),
		omega = omega,
		keep = rep(1, length(omega)),
		z =  rep(NA, M),
		ID = rep(NA, length(omega)),
		sex = c(rep(0, 5), rep(1, 9), rep(NA, M-14)),	# Note 0 is male and 1 is female.
		collar = c(rep(0, 14), rep(1, M-14))
	)


	init.collar <- function(){
		collar.true <- c(rep(0, 14), rep(1, M-14))
		lambda <- runif(1, 0.1, 1)
		psex <- rbeta(1, 9, 14)	# Based on collared...
		sigma <- runif(1, 2, 4)
		X <- cbind(runif(M, xlim[1], xlim[2]), 
				  runif(M, ylim[1], ylim[2]))
		d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
		hkj <- exp(-d2/(2*sigma^2))
		sex <- c(rep(0, 5), rep(1, 9), rbinom(M- 14, size = 1, prob = psex))
		ID <- do.call('c',lapply(1:length(omega), 
			FUN = function(x){sample(1:M, 1, prob = hkj[,omega[x]]*(obsCollar[x] == 0 | (obsCollar[x] - 1) == collar.true)*(obsSex[x] == 0 | (obsSex[x] - 1) == sex))}))
		sex[1:14] <- NA
		z <- rep(0, M)
		z[ID] <- 1
		z[1:14] <- 1
		psi <- rbeta(1, sum(z) + 14, M - sum(1-z))	# NA inits...
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
	Rmodel <- nimbleModel(LID_Sex_collar, constants.collar, data.collar, inits = init.collar())
	conf <- configureMCMC(Rmodel)
	conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))

	conf$removeSamplers('X')
	for(i in 1:M){
		conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
			type = 'myJAM', silent = TRUE, control = list(scale = 1, xlim = xlim, ylim = ylim, temp = 0.2))
	}

	conf$removeSamplers('sigma')
	# sigma 1 and 2 slice sampling.
	conf$addSampler(target = 'sigma', 
			type = 'slice', silent = TRUE, scalarComponents = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))		

	# Optimized z sampler
	conf$removeSamplers('z')
	conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)
	
	conf$removeSamplers('ID')
	conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
	Rmcmc <- buildMCMC(conf)
	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

	Cmcmc$run(60000)
	mvSamples <- Cmcmc$mvSamples
	samples <- as.matrix(mvSamples)
	out <- mcmc(samples[-(1:20000),])
	plot(out[,c("lambda", "sigma","N", "D")])
	save(out, file = paste0("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/FisherSimulations/Scenario_2_LID_Sex_Collar_iter_", h, ".Rda"))
	tmp <- summary(out[,c("lambda", "sigma","N", "D")])
	res3 <- data.frame(cbind(tmp[[1]], tmp[[2]]))
	res3$Method <- "LID_Sex_Collar"
	res3$iter <- h
	res <- rbind(res, res3)

	res
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

