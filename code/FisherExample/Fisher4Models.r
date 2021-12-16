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
library(ggplot2)

source("../Functions/NimbleFunctions.R")

## Load the Fisher data
###--- set directory and load files
load("../../data/FisherData.Rda")

buffer <- 6
traps <- fisher.data[["traps"]]

avg.dists <- apply(traps, 1, FUN = function(x){sqrt((traps[,1] - x[1])^2 + (traps[,2]-x[2])^2)})
avg.min <- apply(avg.dists, 1, FUN = function(x){min(x[x!=0])})
avg <- mean(avg.min)

# traps <- traps[-(1:13),]
xlim <- range(traps[,1])+c(-buffer,buffer)
ylim <- range(traps[,2])+c(-buffer,buffer)
area <- diff(xlim)*diff(ylim)/100

# Lincoln-Petersen: Full Study
24/(16/24)
# Lincoln-Peterson: 2-month genetic
24/(11/17)

obs <- fisher.data[["observations"]]
collars <- fisher.data$collarInfo
# obs$day <- as.numeric(format(obs$DateTime, "%d"))
# obs$month <- as.numeric(format(obs$DateTime, "%m"))

# ggplot(data = obs, aes(x = day, colour=factor(sex))) + geom_density() + facet_wrap(~month)
# ggplot(data = obs, aes(x = day, colour=factor(collar))) + geom_density() + facet_wrap(~month)

# Start: January 20, 2016
start.date <- as.Date("01-01-2016", format = "%d-%m-%Y")
end.date <- as.Date("05-03-2016", format = "%d-%m-%Y")

# 207 detections.
obs <- fisher.data[["observations"]]
obs <- obs[as.Date(obs$DateTime) >= start.date & end.date > as.Date(obs$DateTime),]
nrow(obs)

table(obs$sex)/nrow(obs)*100
table(obs$collar)/nrow(obs)*100
# obs <- obs[obs$sex != 0,]
# nrow(obs)
# hrs <- as.numeric(format(obs$DateTime, "%H")) 
# plot(hrs, jitter(obs$sex))
omega <- obs$TrapNumber
StudyPeriod <- as.numeric(end.date - start.date)
studyArea <- diff(xlim)*diff(ylim)


M <- 400
J <- nrow(traps)

# Chandler and Royle Spatial Count model, Algorithm 2
# Repeat of Burgar et al. 2018 analysis for paper.
counts <- numeric(J)
for(j in 1:J)
{
	counts[j] <- sum(obs$TrapNumber == j)
}

# unique traps:
avg.traps <- function(x)
{
	n <- NULL
	for(i in unique(x))
	{
		n <- c(n, length(unique(omega[x == i])))
	}
	n <- mean(n)
}


Model1 <- nimbleCode( {
	# Priors:
	sigma ~ dunif(0, 50)
	lambda ~ dunif(0, 20)
	psi ~ dbeta(1,1)
	tau2 <- 1/(2*sigma^2)
	
	for(k in 1:M){
		z[k] ~ dbern(psi)
		X[k,1] ~ dunif(xlim[1], xlim[2])
		X[k,2] ~ dunif(ylim[1], ylim[2])
		d2[k,1:J] <- pow(X[k,1] - traps[1:J,1], 2) + pow(X[k,2] - traps[1:J,2], 2)
		hkj[k,1:J] <- lambda*exp(-d2[k,1:J] * tau2 )*z[k]
	}

	### Model for Observation Process:
	for(j in 1:J){
		Hj[j] <- sum(hkj[1:M,j])*Time
		counts[j] ~ dpois(Hj[j])
	}
	# counts[1:J] ~ dPoisSC(lambda = hkj[1:M,1:J], J = J)

	
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
	# conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		# type = 'myJAM', silent = TRUE, 
		# control = list(scale = 1.5, xlim = xlim, ylim = ylim, temp = 0.1))
}
conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE) #, control = list(adaptive = FALSE, scaleWidth = 0.25))	

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Cmcmc$run(10000, time = TRUE)
# mean(Cmcmc$getTimes())
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples[-(1:5000),])
# plot(out[,c('sigma', 'lambda', 'N', 'D')])
# MCMCglmm::posterior.mode(out)
# valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[803]], "width")
# valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[600]], "scale")


# Note that when doing adaptive sampling it was around scale ~= 1.7. We will set it at that for all chains.
samples1 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(init1(), init1(), init1()))

out1 <- mcmc.list(list(as.mcmc(samples1[[1]]), as.mcmc(samples1[[2]]), as.mcmc(samples1[[3]])))
save(out1, file = "../../output/fisherModel1.Rda")
# load("../../output/fisherModel1.Rda")
# plot(out1[,c("sigma", "lambda", "psi")])
# plot(out1[,c("sigma", "N", "D")])
# mean(diff(samples1[[1]][, 'sigma']) > 0)
# MCMCglmm::posterior.mode(out1)
# summary(out1)
# effectiveSize(out1)
# stat <- NULL
# for( i in c('sigma', 'lambda', 'psi', 'N', 'D')) stat[[i]] = gelman.diag(mcmc.list(list(out1[[1]][,i], out1[[2]][,i], out1[[3]][,i] )))
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
    zeros = rep(0, M),
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
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'X'))
# Use a block update on locations. Saves time.
# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'myJAM', silent = TRUE, 
		control = list(scale = 0.5, xlim = xlim, ylim = ylim, temp = 0.2))
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

# Cmcmc$run(25000, time = TRUE)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out.mpp <- mcmc(samples[-(1:10000),])
# plot(out.mpp[,c('sigma', 'lambda', 'N', 'D')])
# mean(diff(out[,c('sigma')]) > 0)
# valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[403]], "width")
# MCMCglmm::posterior.mode(out[,c('sigma', 'lambda', 'N', 'D')])

# post.id <- out[,grep("ID", colnames(out))]
# post.x <- out[, grep("X", colnames(out))]
# x.vals.1 <- do.call('rbind', lapply(1:nrow(out), FUN = function(x){cbind(post.x[x, 1:400][post.id[x,]], post.x[x, 401:800][post.id[x,]])}))

# plot(traps, xlim = xlim, ylim = ylim, col = 'red', pch = 4)
# points(x.vals.1, pch = 20)
# points(traps[omega[1],1], traps[omega[1], 2], col = 'red', pch = 16, cex = 2)

# plot(density(out1[[1]][,"sigma"]))
# lines(density(out2[[1]][,"sigma"]), col = 'red')
# lines(density(out[,"sigma"]), col = 'blue')

samples2 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(inits.mpp(), inits.mpp(), inits.mpp()) )

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
# plot(out2[, c('lambda', 'N')])
# MCMCglmm::posterior.mode(out2)

# summary(out2)
# effectiveSize(out2)
# effectiveSize(out3)
# load("../../output/fisherModel2.Rda")	

# Plot 1 and 2 methods together:
out1D <- do.call("c", out1[,"D"])
out2D <- do.call("c", out2[,"D"])
out3D <- do.call("c", out4[,"D"])
plot(density(out2D), main = "Animal Density", ylim = c(0,1.5), xlab = "Density (Ind/100km2)")
lines(density(out1D), col = "blue")
abline(v = 2, col = 'red', lty = 'dashed')
lines(density(out3D), col = 'red')

####################################
# Model 3: Sex as covariates
####################################
Model3 <- nimbleCode({
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

obsSex <- obs$sex
constants3 <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
    n_obs = length(omega),
	area = diff(xlim)*diff(ylim)/100,
	obsSex = obsSex
)

data3 <- list(
    zeros = rep(0, M),
    omega = omega,
	keep = rep(1, length(omega)),
	z =  rep(NA, M),
	ID = rep(NA, length(omega)),
	sex = rep(NA, M)	# Note 0 is male and 1 is female.
)

# Need to fully initialize because I haven't implemented a random generator distribution for dID()...
# Bit of a pain to make sure we match sex and collar correctly.
init3 <- function(){
	N <- floor(runif(1, 1, M/2))
	psi <- rbeta(1, N, M-N)	# NA inits...	
	lambda <- runif(1, 0.1, 1)
	psex <- rbeta(1, 1, 1)	# Based on collared...
	sigma <- runif(1, 2, 4)
	X <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	sex <- rbinom(M, size = 1, prob = psex)
	ID <- do.call('c',lapply(1:length(omega), FUN = function(x){sample(1:M, 1, prob = hkj[1:M,omega[x]]*(sex+1 == obs$sex[x] |  obs$sex[x] == 0))}))
	Hk <- rowSums(hkj)*StudyPeriod*lambda
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
Rmodel <- nimbleModel(Model3, constants3, data3, inits = init3())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'z', 'sex', 'psex')) # 'psex', 'sex',

conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'myJAM', silent = TRUE, 
		control = list(scale = 0.5, xlim = xlim, ylim = ylim, temp = 0.2))
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

# Cmcmc$run(25000, time = TRUE)
# Cmcmc$getTimes()
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples[-(1:10000),])
# plot(out[,c('sigma', 'lambda')])
# plot(out[,c('N', 'D')])
# MCMCglmm::posterior.mode(out)
# post <- samples[-(1:10000),grep("ID", colnames(samples))]
# K <- apply(post, 1, FUN = function(x){length(unique(x))})
# cols <- grepl("sex", colnames(samples)) & !grepl("psex", colnames(samples))
# post.sex <- samples[-(1:10000), cols]
# KF1 <- do.call('c', lapply(1:nrow(post.sex),  FUN = function(x){ sum(post.sex[x,unique(post[x,])])}))
# hist(K)
# hist(KF1)

# The full run...
samples3 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(init3(), init3(), init3()))

# Compute the number of actually observed fisher to compare with SCR model in paper.
post1 <- samples3[[1]][,grep("ID", colnames(samples3[[1]]))]
post2 <- samples3[[2]][,grep("ID", colnames(samples3[[2]]))]
post3 <- samples3[[3]][,grep("ID", colnames(samples3[[3]]))]
NObs1 <- apply(post1, 1, FUN = function(x){ length(unique(x))})
NObs2 <- apply(post2, 1, FUN = function(x){ length(unique(x))})
NObs3 <- apply(post3, 1, FUN = function(x){ length(unique(x))})

# trapK1 <- apply(post1, 1, avg.traps)
# trapK2 <- apply(post2, 1, avg.traps)
# trapK3 <- apply(post3, 1, avg.traps)

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

save(out3, file = "../../output/fisherModel3.Rda")	

plot(out3[, c("sigma", "lambda", "D")])
plot(out3[, "trapK"])
plot(out3[, c("K", "KM", "KF")])

MCMCglmm::posterior.mode(out3)


####################################
# Model 4: Sex and collar as covariates
####################################
Model4 <- nimbleCode({
	
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
		# pCollar[i] <- (collar[ID[i]]+1) == obsCollar[i] | obsCollar[i] == 0
		pMatch[i] <- matchCollar[ID[i],i]*pSex[i]
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
collars <- fisher.data$collarInfo
matchCollar <- matrix(1, nrow = M, ncol = length(omega))
collars$Sex <- 2 - (collars$Sex == "M")
dates <- obs$DateTime
dates <- as.Date(dates)
obsCollar <- obs$collar

# Process collar matrix:
for(i in 1:length(omega))
{
	indx.match <- collars$Start < dates[i] & collars$Finish >= dates[i]
	if(obsCollar[i] == 1)
	{
		matchCollar[(nrow(collars)+1):M,i] <- 0
		matchCollar[1:nrow(collars), i] <- indx.match*1
	}else{
		if(obsCollar[i] == 2) matchCollar[1:nrow(collars),i] <- (1-indx.match)		
	}
	if(dates[i] > as.Date("2016-01-25", format = "%Y-%m-%d")){
		matchCollar[which(collars$FisherID == "M02"),i] <- 0
	}
}

constants4 <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
    n_obs = length(omega),
	area = diff(xlim)*diff(ylim)/100,
	obsSex = obs$sex,
	obsMark = obs$collar
	)

data4 <- list(
    zeros = rep(0, M),
    omega = omega,
	keep = rep(1, length(omega)),
	z =  c(rep(1, nrow(collars)), rep(NA, M-nrow(collars))),
	ID = rep(NA, length(omega)),
	sex = c(collars$Sex - 1, rep(NA, M - nrow(collars))),	# Note 0 is male and 1 is female.
	matchCollar = matchCollar	# needs to be data as we can't legally dynamically index constants...
)

init4 <- function(){
	lambda <- runif(1, 0.1, 1)
	psex <- rbeta(1, 9, 14)	# Based on collared...
	sigma <- runif(1, 2, 4)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	sexCollar <- collars$Sex-1
	sex <- c(sexCollar, rbinom(M- nrow(collars), size = 1, prob = psex))
	ID <- do.call('c',lapply(1:length(omega), 
		FUN = function(x){sample(1:M, 1, prob = hkj[,omega[x]]*matchCollar[,x]*(obs$sex[x] == 0 | (obs$sex[x] - 1) == sex))}))
	sex[1:nrow(collars)] <- NA
	z <- rep(0, M)
	z[ID] <- 1
	z[1:nrow(collars)] <- NA
	psi <- rbeta(1, sum(z, na.rm = TRUE) + nrow(collars), M - sum(1-z, na.rm = TRUE))	# NA inits...
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

Cmcmc$run(20000, time = TRUE)
Cmcmc$getTimes()
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:10000),])
plot(out[,c('sigma','lambda')])
plot(out[,c('psex')])
plot(out[,c('N', 'D')])

summary(out[,c('sigma', 'lambda', 'psex','N', 'D')])
post1 <- out[,grep("ID", colnames(out))]
NObs1 <- apply(post1, 1, FUN = function(x){ length(unique(x))})
NObs2 <- apply(post1, 1, FUN = function(x){ length(unique(x[x <= 14]))})
NObs3 <- apply(post1, 1, FUN = function(x){ length(unique(x[x <= 5]))})

# out4 <- mcmc(cbind(out[,c('sigma', 'lambda', 'psex','N', 'D')], "K" = NObs1))
# summary(out4)
# load("../../output/geneticFisherMCMC.Rda")

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
save(out4, file = "../../output/fisherModel4_z15.Rda")
# load(file = "../../output/fisherModel4.Rda")

# out4 <- mcmc.list(list(as.mcmc(cbind(samples4[[1]][,c('sigma[1]', 'sigma[2]', 'lambda[1]', 'lambda[2]', 'psi', 'N', 'D', 'psex')], "K" = NObs1, "KF" = NObsF1, "KM" = NObs1 - NObsF1)), 
	# as.mcmc(cbind(samples4[[2]][,c('sigma[1]', 'sigma[2]', 'lambda[1]', 'lambda[2]', 'psi', 'N', 'D', 'psex')], "K" =  NObs2, "KF" = NObsF2, "KM" = NObs2 - NObsF2)),
	# as.mcmc(cbind(samples4[[3]][,c('sigma[1]', 'sigma[2]', 'lambda[1]', 'lambda[2]', 'psi', 'N', 'D', 'psex')], "K" = NObs3, "KF" = NObsF3, "KM" = NObs3 - NObsF3))))

# MCMCglmm::posterior.mode(out4[[1]][,"D"])
# plot(out4[[1]][,c("sigma[1]", "sigma[2]")])
# plot(as.numeric(out4[[1]][,c("sigma[1]")]), as.numeric(out4[[1]][,c("sigma[2]")]))

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
load("../../output/fisherModel3.Rda")
load("../../output/fisherModel4.Rda")
# load("../../output/geneticFisherMCMC.Rda")

out1.df <- do.call("rbind", lapply(out1, as.data.frame))
out1.df$K <- NA
out1.df$Method = "SC"
MCMCglmm::posterior.mode(out1)

out2.df <- do.call("rbind", lapply(out2, as.data.frame))
out2.df$Method = "LID"
MCMCglmm::posterior.mode(out2)

out3.df <- do.call("rbind", lapply(out3, as.data.frame))
out3.df$Method = "LID+Sex"
MCMCglmm::posterior.mode(out3)

out4.df <- do.call("rbind", lapply(out4, as.data.frame))
out4.df$Method = "LID+Sex/Collar"
MCMCglmm::posterior.mode(out4)

# out.df <- data.frame(out)
# out.df$Method = "Genetic SCR"
# out.df$K <- 24
# out.df$lambda <- out.df$lambda/0.42

facet_names <- c('D~(Ind~100~km^-2)', 'Number~of~detected~animals' , 'sigma~(km)', 'lambda~(day^-1)')
names(facet_names) <- c("D", "K", "sigma", "lambda")

all.dat1 <- rbind(out1.df[, c("sigma", "lambda", "D", "K", "Method")], 
	out2.df[, c("sigma", "lambda", "D", "K", "Method")],
	out3.df[, c("sigma", "lambda", "D", "K", "Method")],
	out4.df[, c("sigma", "lambda", "D", "K", "Method")])
	
all.dat.l <- melt(all.dat1, id.vars = "Method")
all.dat.l <- all.dat.l[!(all.dat.l$variable == "K" & all.dat.l$Method == "SC"),]
all.dat.l$Method <- factor(all.dat.l$Method, levels = c("SC", "LID", "LID+Sex", "LID+Sex/Collar"))

all.dat.l$variable <- factor(all.dat.l$variable, names(facet_names))

ggplot(data = all.dat.l, aes(y = value, x = Method)) + 
	facet_wrap(~variable, scale = "free", labeller = as_labeller(x = facet_names, label_parsed)) + 
	theme_classic() +
	geom_boxplot(width = 0.075, outlier.alpha = 0, fill = "grey", colour = "black") + 
	geom_violin(alpha = 0, adjust = 3) + 
	ylab("") + xlab("")
ggsave("../../output/FisherResults/CombinedFisherMCMCOutput.png", dpi = 'print', 
	width = 9.75, height = 7.35, units = 'in')




# Plotting spatial information from survey:
library(secr)
library(ggpubr)

mask <- make.mask(traps = traps, buffer = 6, spacing = 0.5, type = "trapbuffer")
A <- attr(mask, "area")
mask <- as.matrix(mask)
attr(mask, "area") <- A
attr(mask, "buffer") <- 15

d2 <- t(apply(mask, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))

g1 <- 1 - exp(-0.15*exp(-d2/(2*1.5^2))*64)
g2 <- 1 - exp(-0.1*exp(-d2/(2*2^2))*64)

plot.dat.fish <- data.frame(rbind(mask, mask), val = c(g1[,45], g2[,45]), "Scenario" = rep(c("A", "B"), each = nrow(mask)))

p.fish <- ggplot(data = plot.dat.fish, aes(x = x, y = y, z = val, colour = Scenario)) + 
	geom_contour(aes(linetype = Scenario), breaks = c(0.90), size = 1) +
	theme_classic() + xlim(range(mask[,1])) + ylim(range(mask[,2])) +
	geom_point(data = data.frame(traps), aes(x=x, y=y, z= NULL), col = 'red', shape = 3, size = 2) +
	geom_point(data = data.frame(traps)[45,], aes(x=x, y=y, z= NULL), col = 'red', shape = 1, size = 4) +
	xlab("x (km)") + ylab("y (km)") + 
	scale_color_manual(name = "Detection Function", 
		labels = c("A" = expression(sigma==1.5~lambda==0.15), "B" = expression(sigma==2~lambda==0.10)),
		values = c("#99E1D9", "#548687")) + 
	scale_linetype_discrete(name = "Detection Function", labels = c("A" = expression(sigma==1.5~lambda==0.15), "B" = expression(sigma==2~lambda==0.10))) + 
	coord_fixed() + 
	theme(legend.position = "bottom")
p.fish = p.fish + theme(legend.text = element_text(size=12))
ggsave("../../output/FisherResults/FisherLocationExample.png", dpi = 'print', 
	width = 5.5, height = 5.5, units = 'in')

load(file = "../../output/FrogResults/FrogLocationsExample.Rda")
load("../../data/stacked-lightfooti.Rdata")
test.capt <- capt.all$bincapt[capt.all$bincapt[,7] == 104, 1:6]

p.frog <- ggplot(data = plot.dat, aes(x = x, y = y, z = posterior, colour = sigmatval)) + 
	geom_contour(aes(linetype = sigmatval), breaks = c(0.9), size = 1) +
	theme_classic() + xlim(range(mask[,1])) + ylim(range(mask[,2])) +
	geom_point(data = data.frame(traps), aes(x=x, y=y, z= NULL), col = 'red', shape = 3, size = 2) +
	geom_point(data = data.frame(traps)[test.capt[3,] == 1,], aes(x=x, y=y, z= NULL), col = 'red', shape = 1, size = 4) +
	xlab("x (m)") + ylab("y (m)") + 
	scale_color_manual(name = expression(sigma[t]~(sec)),
		values = c("#99E1D9", "#548687", "#414770")) + 
	scale_linetype_discrete(name = expression(sigma[t]~(sec))) + 
	coord_fixed() + theme(legend.position = "bottom")
p.frog = p.frog + theme(legend.text = element_text(size=12))

figure <- ggarrange(p.fish, p.frog,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
figure					
ggexport(figure, filename = "../../output/FisherResults/CombinedFisherFrogLocations.png", dpi = 'print',
	  width = 850, height = 480)
