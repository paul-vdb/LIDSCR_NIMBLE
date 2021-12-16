#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################
#################################
setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/FrogExample")

# load("../../output/FrogsLatentID2Session.Rda")
# summary(out)

library(coda)
library(nimble)
library(nimbleSCR)

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

source("../Functions/SimData.R")
load("../../data/stacked-lightfooti.Rdata")

results <- foreach(h=1:110,
			.packages = c("coda", "nimbleSCR", "nimble"))%dopar%{
	source("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/Functions/NimbleFunctions.R")
	
	xlim <- range(mask[,1])
	ylim <- range(mask[,2])
	area <- diff(xlim)*diff(ylim)

	# Constants:
	nu <- 330
	J <- nrow(traps)
	Time <- 30

	M <- 300
	N <- 55
	
	### Scenario 2: 2 sessions
	dat1 <- simASCR(N = N, sigma = 2.3, sigma_toa = 0.00055, g0 = 5.75, lambda = 0.28, 
		StudyPeriod = Time, traps = traps, 
		limits = list(xlim = xlim, ylim = ylim))
	dat2 <- simASCR(N = N, sigma = 2.3, sigma_toa = 0.00055, g0 = 5.75, lambda = 0.28, 
		StudyPeriod = Time, traps = traps, 
		limits = list(xlim = xlim, ylim = ylim))

	toa <- rbind(dat1$toa, dat2$toa)
	capt <- rbind(dat1$capt, dat2$capt)
	occ <- c(rep(1, nrow(dat1$toa)), rep(2, nrow(dat2$toa)))
	n <- nrow(capt)

	# Run the model twice. Once for known ID and once for unknown ID.
	inits <- function(){
		ID <- numeric(length(occ))
		ID[occ == 1] <- 1:sum(occ == 1)
		ID[occ == 2] <- 1:sum(occ == 2)
		psi <- runif(1, 0.1, 0.3)
		z <- matrix(0, ncol = 2, nrow = M)
		z[ID[occ==1],1] <- 1
		z[ID[occ==2],2] <- 1
		z[z==0] <- rbinom(sum(1-z), size = 1, prob = psi)
		lambda = runif(1, 0.15, 0.5)
		sigma = runif(1, 1.5, 4)
		g0 = runif(1,5,10)

		dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
		pkj <- (1-exp(-g0*exp(-dmask2/(2*sigma^2))))
		panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x + .Machine$double.eps)%*%t(capt) + log(1-x + .Machine$double.eps)%*%t(1-capt))})
		X <- array(c(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]), 
			runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2])), c(M,2,2))		
		
		for(k in 1:2){
			for(i in 1:M){
				if(sum(ID[occ == k] == i) == 0) next;
				if(sum(ID[occ == k] == i) == 1) pID <- panimal[which(ID[occ == k] == i), ]
				if(sum(ID[occ == k] == i) > 1) pID <- colSums(panimal[which(ID[occ == k] == i), ])
				mpt <- sample(ncol(panimal), 1, prob = exp(pID))
				X[i,,k] <- mask[mpt,]
			}
		}

		list(
			lambda = lambda,
			psi = psi,
			sigma = sigma,
			sigmatoa = runif(1, 0.01, 1),
			g0 = g0,
			X = X,
			ID = ID,
			z = z
		)
	}

	code <- nimbleCode({
		lambda ~ dunif(0, 10) # Detection rate at distance 0
		psi ~ dbeta(1, 1)      # Prior on data augmentation bernoulli vec.
		sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
		tau2 <- 1/(2*sigma^2)
		sigmatoa ~ dunif(0, 1)
		g0 ~ dunif(0, 20)
		for(v in 1:n_occ) {
			for(k in 1:M) {
				z[k,v] ~ dbern(psi)
				X[k, 1, v] ~ dunif(xlim[1], xlim[2])
				X[k, 2, v] ~ dunif(ylim[1], ylim[2])
				d2[k,1:J, v] <- (X[k,1,v]-traps[1:J,1])^2 + (X[k,2,v]-traps[1:J,2])^2
				expTime[k, 1:J, v] <- sqrt(d2[k,1:J,v])/nu
				pkj[k,1:J,v] <- (1-exp(-g0*exp(-d2[k,1:J,v]*tau2)))
				# Hazard rate for animal across all traps.
				Hk[k,v] <-(1-prod(1-pkj[k,1:J,v]))*lambda*Time*z[k,v]
				# The 'undetected' part of the likelihood.
				zeros[k,v] ~ dpois(Hk[k,v])
			}
			
			# Predicted population size
			N[v] <- sum(z[1:M,v])
		}
		
		pID[1:M, 1:n_occ] <- z[1:M,1:n_occ]*lambda
		
		# Trap history model.
		# and unobserved animal ID.
		for(i in 1:n_obs) {
			# Bernoulli capture history for each call that depends on ID
			y[i,1:J] ~ dbinom_vector(size = trials[1:J], pkj[ID[i],1:J, occ[i]])
			# Time of arrival, depends on which traps actually recorded it.
			toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i], 1:J, occ[i]], sd = sigmatoa, y = y[i,1:J])
			# The likelihood needs to be multiplied by lambda for each detection and
			# I need ID to be a stochastic node. 2 birds...
			ID[i] ~ dID(pID = pID[1:M, occ[i]])
		}

		# Derived Variables.
		EN <- psi*M
		D <- EN/area*10000
	})

	constants <- list(
		J = J,
		xlim = xlim,
		ylim = ylim,
		traps = traps, 
		Time = Time,
		M = M,
		n_obs = nrow(capt),
		trials = rep(1, J),
		nu = nu,
		area = area,
		n_occ = max(occ),
		occ = occ)

	data <- list(
		zeros = matrix(0, nrow = M, ncol = max(occ)),
		y = capt,
		toa = toa,
		z = cbind(rep(NA, M), rep(NA, M)),
		ID = rep(NA, nrow(capt))
	)

	Rmodel <- nimbleModel(code, constants, data, inits = inits())

	conf <- configureMCMC(Rmodel)

	conf$setMonitors(c('psi', 'sigma', 'lambda', 'sigmatoa', 'g0', 'EN', 'N', 'D'))

	conf$removeSamplers('sigmatoa')
	conf$addSampler(target = 'sigmatoa', type = 'RW', control = list(log = TRUE, adaptive = TRUE))

	conf$removeSamplers('g0')
	conf$addSampler(target = 'g0', 
			type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))		

	conf$removeSamplers('X')
	for(v in 1:2){
		for(i in 1:M) {
			conf$addSampler(target = paste0('X[', i, ', 1:2,', v, ']'), 
				type = 'myJAM', silent = TRUE, 
				control = list(scale = 1, xlim = xlim, ylim = ylim, temp = 0.2, occ = 2))
		}
	}

	conf$removeSamplers('z')
	conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE, 
		control = list('Noccasion' = 2, 'IDoccasion' = occ))

	conf$removeSamplers('ID')
	conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))

	Rmcmc <- buildMCMC(conf)
	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

	Cmcmc$run(40000)
	mvSamples <- Cmcmc$mvSamples
	samples <- as.matrix(mvSamples)
	out1 <- mcmc(samples[-(1:20000),])
	out1 <- out1[,c('psi', 'sigma', 'lambda', 'sigmatoa', 'g0', 'EN', 'N[1]','N[2]', 'D')]
	# plot(out1[, c("N[1]", "N[2]", "D", "sigma", "sigmatoa")])
	save(out1, file = paste0("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/ASCRSimulations/ASCR_Scenario1_LatentID_iter_", h, ".Rda"))
	# summary(out)
	output1 <- data.frame(do.call("cbind", summary(out1)))

	###############################################
	# Now do the known ID ASCR Model:
	###############################################
	ID <- c(dat1$obs$ID, dat2$obs$ID)
	initsID <- function(){
		psi <- runif(1, 0.1, 0.3)
		z <- matrix(rbinom(M*2, 1, psi), ncol = 2, nrow = M)
		z[ID[occ==1],1] <- NA
		z[ID[occ==2],2] <- NA
		lambda = runif(1, 0.1, 2)
		sigma = runif(1, 3, 5)
		g0 = runif(1,1,10)

		dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
		pkj <- (1-exp(-g0*exp(-dmask2/(2*sigma^2))))
		panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x + .Machine$double.eps)%*%t(capt) + log(1-x + .Machine$double.eps)%*%t(1-capt))})
		X <- array(c(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]), 
			runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2])), c(M,2,2))		
		
		for(k in 1:2){
			for(i in 1:M){
				if(sum(ID[occ == k] == i) == 0) next;
				if(sum(ID[occ == k] == i) == 1) pID <- panimal[which(ID[occ == k] == i), ]
				if(sum(ID[occ == k] == i) > 1) pID <- colSums(panimal[which(ID[occ == k] == i), ])
				mpt <- sample(ncol(panimal), 1, prob = exp(pID))
				X[i,,k] <- mask[mpt,]
			}
		}

		list(
			lambda = lambda,
			psi = psi,
			sigma = sigma,
			sigmatoa = runif(1, 0.01, 1),
			g0 = g0,
			X = X,
			z = z
		)
	}

	constants.id <- list(
		J = J,
		xlim = xlim,
		ylim = ylim,
		traps = traps, 
		Time = Time,
		M = M,
		n_obs = nrow(capt),
		trials = rep(1, J),
		nu = nu,
		area = area,
		n_occ = 2,
		occ = occ
	)

	data.id <- list(
		zeros = matrix(0, nrow = M, ncol = max(occ)),
		ID = ID,
		y = capt,
		toa = toa,
		z = cbind(c(rep(1, max(ID[occ == 1])), rep(NA, M-max(ID[occ == 1]))),
			c(rep(1, max(ID[occ == 2])), rep(NA, M-max(ID[occ == 2]))))
	)
	Rmodel <- nimbleModel(code, constants.id, data.id, inits = initsID())
	conf <- configureMCMC(Rmodel)
	conf$setMonitors(c('psi', 'sigma', 'lambda', 'sigmatoa', 'g0', 'EN', 'N', 'D'))

	conf$removeSamplers('sigmatoa')
	conf$addSampler(target = 'sigmatoa', type = 'RW', control = list(log = TRUE, adaptive = TRUE))

	conf$removeSamplers('g0')
	conf$addSampler(target = 'g0', 
			type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))		

	conf$removeSamplers('X')
	for(v in 1:2){
		for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2,', v, ']'), type = 'RW_block', silent = TRUE)
	}
	Rmcmc <- buildMCMC(conf)
	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

	Cmcmc$run(40000)
	mvSamples <- Cmcmc$mvSamples
	samples <- as.matrix(mvSamples)
	out2 <- mcmc(samples[-(1:20000),])
	out2 <- out2[,c('psi', 'sigma', 'lambda', 'sigmatoa', 'g0', 'EN', 'N[1]','N[2]', 'D')]
	plot(out2)
	save(out2, file = paste0("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/ASCRSimulations/ASCR_Scenario1_KnownID_iter_", h, ".Rda"))	
	output2 <- data.frame(do.call("cbind", summary(out2)))
	
	output1$Method = "Latent ID"
	output2$Method = "Known ID"

	rbind(output1, output2)
}

stopCluster(cl)

library(ascr)
library(secr)
traps2 <- convert.traps(traps)
mask <- make.mask(traps = traps2, buffer = 15, spacing = 0.2, type = "trapbuffer")
A <- attr(mask, "area")
mask <- as.matrix(mask)
attr(mask, "area") <- A
attr(mask, "buffer") <- 15

multi.capt <- list(list("bincapt" = dat1$capt, "toa" = dat1$toa), list("bincapt" = dat2$capt, "toa" = dat2$toa))
ascr.res <- fit.ascr(capt = multi.capt, list(traps, traps), mask = list(mask, mask), 
			detfn = "hhn", survey.length = c(30, 30))
# citation("ascr")
summary(ascr.res)
confint(ascr.res)





library(Rcpp)
## Paul's function.
source("C:/Users/Paul/Documents/GitHub/scr-cuerate/fit-cuerate-scr.r")
sourceCpp("C:/Users/Paul/Documents/GitHub/scr-cuerate/fit-cuerate-scr.cpp")
ids <- c(dat1$obs$ID, dat2$obs$ID) + (occ-1)*100
start2 <- c(400, 0.3, 2, 9, 10)
fit.cuerate <- cuerate.scr.fit(capt, ids, traps, mask, detfn = "hhn",
                               start = start2, toa = toa, trace = FALSE)
fit.cuerate$res['D','Estimate']/2
summary(out2)[[1]]
55/area*10000




















# Run the model twice. Once for known ID and once for unknown ID.
inits <- function(){
	ID <- 1:length(occ)
	ID[occ == 2] <- (1:sum(occ==2))+M/2
	psi <- runif(1, 0.2, 0.5)
	z <- numeric(M)
	z[ID] <- 1
	z[z==0] <- rbinom(sum(1-z), size = 1, prob = psi)
	lambda = runif(1, 0.15, 0.5)
	sigma = runif(1, 1.5, 4)
	g0 = runif(1,5,10)

	dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
	pkj <- (1-exp(-g0*exp(-dmask2/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x + .Machine$double.eps)%*%t(capt) + log(1-x + .Machine$double.eps)%*%t(1-capt))})
	X <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))		

	for(i in 1:M){
		if(sum(ID == i) == 0) next;
		if(sum(ID == i) == 1) pID <- panimal[which(ID == i), ]
		mpt <- sample(ncol(panimal), 1, prob = exp(pID))
		X[i,] <- mask[mpt,]
	}

	list(
		lambda = lambda,
		psi = psi,
		sigma = sigma,
		sigmatoa = runif(1, 0.01, 1),
		g0 = g0,
		X = X,
		ID = ID,
		z = z
	)
}

code <- nimbleCode({
	lambda ~ dunif(0, 10) # Detection rate at distance 0
	psi ~ dbeta(1, 1)      # Prior on data augmentation bernoulli vec.
	sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
	tau2 <- 1/(2*sigma^2)
	sigmatoa ~ dunif(0, 1)
	g0 ~ dunif(0, 20)
	for(k in 1:M) {
		z[k] ~ dbern(psi)
		X[k, 1] ~ dunif(xlim[1], xlim[2])
		X[k, 2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
		expTime[k, 1:J] <- sqrt(d2[k,1:J])/nu
		pkj[k,1:J] <- (1-exp(-g0*exp(-d2[k,1:J]*tau2)))
		# Hazard rate for animal across all traps.
		Hk[k] <-(1-prod(1-pkj[k,1:J]))*lambda*Time*z[k]
		# The 'undetected' part of the likelihood.
		zeros[k] ~ dpois(Hk[k])
	}
	
	# Predicted population size
	N <- sum(z[1:M])
	
	pID[1:M,1] <- z[1:M]*lambda*occz[1:M,1]
	pID[1:M,2] <- z[1:M]*lambda*occz[1:M,2]
	
	# Trap history model.
	# and unobserved animal ID.
	for(i in 1:n_obs) {
		# Bernoulli capture history for each call that depends on ID
		y[i,1:J] ~ dbinom_vector(size = trials[1:J], pkj[ID[i],1:J])
		# Time of arrival, depends on which traps actually recorded it.
		toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i], 1:J], sd = sigmatoa, y = y[i,1:J])
		# The likelihood needs to be multiplied by lambda for each detection and
		# I need ID to be a stochastic node. 2 birds...
		ID[i] ~ dID(pID = pID[1:M, occ[i]])
	}

	# Derived Variables.
	EN <- psi*M
	D <- N/(2*area)*10000
})

constants <- list(
	J = J,
	xlim = xlim,
	ylim = ylim,
	traps = traps, 
	Time = Time,
	M = M,
	n_obs = nrow(capt),
	trials = rep(1, J),
	nu = nu,
	area = area,
	occz = cbind(c(rep(1, M/2), rep(0, M/2)), c(rep(10, M/2), rep(1, M/2))),
	occ = occ)

data <- list(
	zeros = numeric(M),
	y = capt,
	toa = toa,
	z = rep(NA, M),
	ID = rep(NA, nrow(capt))
)

Rmodel <- nimbleModel(code, constants, data, inits = inits())

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('psi', 'sigma', 'lambda', 'sigmatoa', 'g0', 'EN', 'N', 'D'))

conf$removeSamplers('sigmatoa')
conf$addSampler(target = 'sigmatoa', type = 'RW', control = list(log = TRUE, adaptive = TRUE))

conf$removeSamplers('X')
for(i in 1:M) {
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'myJAM', silent = TRUE, 
		control = list(scale = 0.75, xlim = xlim, ylim = ylim, temp = 0.25))
}


conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

conf$removeSamplers('ID')
conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(20000)














