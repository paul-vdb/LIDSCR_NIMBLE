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
library(secr)

source("NimbleFunctions.R")
source("SimData.R")

## Load the Fisher data
###--- set directory and load files
load("../data/FisherData.Rda")

buffer <- 4
traps <- fisher.data[["traps"]]
xlim <- range(traps[,1])+c(-buffer,buffer)
ylim <- range(traps[,2])+c(-buffer,buffer)
obs <- fisher.data[["observations"]]

omega <- obs$TrapNumber
StudyPeriod <- 64
studyArea <- diff(xlim)*diff(ylim)

traps.df <- data.frame(names = 1:n.traps, traps)
detector <- ifelse(ss, "signal", "proximity")
read.traps(data = traps.df, detector = detector)

M <- 60
J <- nrow(traps)
mask <- make.mask(traps, buffer = buffer)
habitatMask <- convertMask(mask, traps, plot = TRUE)
habMat <- habitatMask$habMat
traps <- habitatMask$trapMat
upperlimit <- habitatMask$upperLimit
scaledMask <- habitatMask$scaledMask


# Marked Poisson process model incorporating the 
# latent ID variable as unknown. 
#------------------------------------------------
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
    Time = StudyPeriod,
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
Rmodel <- nimbleModel(Model23, constants23, data23, inits = init23())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'X', 'z'))
# Use a block update on locations. Saves time.
# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'RW_block', silent = TRUE, control = list(scale = 1.5, adaptive = FALSE))
	}

conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
	type = 'RW', silent = TRUE, control = list(scale = 0.1, adaptive = FALSE))
conf$removeSamplers('lambda')
conf$addSampler(target = 'lambda', 
	type = 'RW', silent = TRUE, control = list(scale = 0.1, adaptive = FALSE))

# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)
# van Dam-Bates categorical sampler
conf$removeSamplers('ID')
# Chandler and Royle Alg. 1 sampler.
conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
# conf$addSampler('ID', type = 'myCat', scalarComponents = TRUE, control = list(M = M, omega = omega))
# van Dam-Bates Alg.
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(15000)
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:10000),])
plot(out[,c('sigma', 'lambda', 'N', 'D')])
# summary(out[,c('sigma', 'lambda', 'N', 'psi', 'D')])
# tmp <- samples[-(1:5000),]
# post.id <- tmp[,grep("ID", colnames(samples))]

# obs <- which(apply(post.id, 1, function(x){ any(x == 41)}))
# plot(traps, xlim = xlim, ylim = ylim, pch = 4, col = 'blue')
# text(traps, label = counts)
# points(tmp[tmp[, 'z[41]']==1, c("X[41, 1]", "X[41, 2]")], col = 'grey')
# points(tmp[obs, c("X[41, 1]", "X[41, 2]")], col = 'red', pch = 16)

# post.id <- samples[,grep("ID", colnames(samples))]
# NActive <- apply(post.id[-(1:5000),], 1, FUN = function(x){ length(unique(x))})
# plot(samples[-(1:5000),'N'] - NActive)

samples2 <- runMCMC(Cmcmc, 60000, nburnin = 20000, nchains = 3, 
	thin = 1, inits = list(init23(), init23(), init23()))

post21 <- samples2[[1]][,grep("ID", colnames(samples2[[1]]))]
post22 <- samples2[[2]][,grep("ID", colnames(samples2[[2]]))]
post23 <- samples2[[3]][,grep("ID", colnames(samples2[[3]]))]
NObs1 <- apply(post21, 1, FUN = function(x){ length(unique(x))})
NObs2 <- apply(post22, 1, FUN = function(x){ length(unique(x))})
NObs3 <- apply(post23, 1, FUN = function(x){ length(unique(x))})

out2 <- mcmc.list(list(as.mcmc(cbind(samples2[[1]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "K" = NObs1)), 
	as.mcmc(cbind(samples2[[2]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "K" =  NObs2)),
	as.mcmc(cbind(samples2[[3]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "K" = NObs3))))
save(out2, file = "../output/fisherModel2.Rda")	
# summary(out2)
# plot(out2, ask = TRUE)
# plot(out2[,c("lambda", "sigma","D")])
# load("../output/fisherModel2.Rda")	
# plot(out2[,c("sigma","N")])
#################################
# New algorithm:
#################################
Rmodel <- nimbleModel(Model23, constants23, data23, inits = init23())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID'))
# Use a block update on locations. Saves time.
# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
conf$removeSamplers('X')
for(i in 1:M){
		conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
			type = 'RW_block', silent = TRUE, control = list(scale = 1.5, adaptive = FALSE))
	}

# conf$removeSamplers('sigma')
# conf$addSampler(target = 'sigma', 
	# type = 'mySigma', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))
# conf$removeSamplers('lambda')
# conf$addSampler(target = 'lambda', 
	# type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))

# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)
# van Dam-Bates categorical sampler
conf$removeSamplers('ID')
conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
# van Dam-Bates Alg.
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Cmcmc$run(20000)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples[-(1:10000),])
# plot(out[,c('sigma', 'lambda', 'N', 'D')])
# post.id <- samples[,grep("ID", colnames(samples))]
# post.z <- samples[,grep("z", colnames(samples))]
# NActive <- apply(post.id, 1, FUN = function(x){ length(unique(x))})
# which(NActive > samples[,'N'])

# samplers <- conf$printSamplers()
# scales <- NULL
# for(i in 1:nrow(){
	# scales <- c(scales,valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[i]], "scale"))
# }

samples3 <- runMCMC(Cmcmc, 60000, nburnin = 20000, nchains = 3, 
	thin = 1, inits = list(init23(), init23(), init23()))

post21 <- samples3[[1]][,grep("ID", colnames(samples3[[1]]))]
post22 <- samples3[[2]][,grep("ID", colnames(samples3[[2]]))]
post23 <- samples3[[3]][,grep("ID", colnames(samples3[[3]]))]
NObs1 <- apply(post21, 1, FUN = function(x){ length(unique(x))})
NObs2 <- apply(post22, 1, FUN = function(x){ length(unique(x))})
NObs3 <- apply(post23, 1, FUN = function(x){ length(unique(x))})

out3 <- mcmc.list(list(as.mcmc(cbind(samples3[[1]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "K" = NObs1)), 
	as.mcmc(cbind(samples3[[2]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "K" =  NObs2)),
	as.mcmc(cbind(samples3[[3]][,c('sigma', 'lambda', 'psi', 'N', 'D')], "K" = NObs3))))
save(out3, file = "../output/fisherModel3.Rda")	
# plot(out3, ask = TRUE)
# plot(out3[,c("lambda", "sigma","N")])
# load("../output/fisherModel3.Rda")	

# effectiveSize(out3)
# stat <- NULL
# for( i in c('sigma', 'lambda', 'psi', 'N', 'D', 'NActive')) stat[[i]] = gelman.diag(mcmc.list(list(out3[[1]][,i], out3[[2]][,i], out3[[3]][,i] )))
# summary(out3)
# plot(out3, ask = TRUE)

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
    Time = StudyPeriod,
    M = M,
    n_obs = length(omega),
	area = diff(xlim)*diff(ylim)/100,
	n_collar = 14,
	obsCollar = obs$collar, # Unknown = 0, Collared = 1, Uncollared = 2
	obsSex = obs$sex
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

Rmodel <- nimbleModel(Model4, constants4, data4, inits = init4())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'psex', 'sex', 'z'))

# conf$removeSamplers('sigma')
# conf$addSampler(target = 'sigma', scalarComponents = TRUE,
	# type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))
# conf$removeSamplers('lambda')
# conf$addSampler(target = 'lambda', scalarComponents = TRUE,
	# type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))
# Use a block update on locations. Saves time.
conf$removeSamplers('X')
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
	type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 1.5))
# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

conf$removeSamplers('ID')
# conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Do a test run on the MCMC object.
# Cmcmc$run(10000)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples[-(1:5000),])
# plot(out[,c('sigma', 'D', 'psex')])	# proportion of males/females matches SCR data. 0.6
# dev.new()
# plot(out[,c('sigma', 'lambda')])
# post1 <- samples[-(1:5000),grep("ID", colnames(samples))]
# Nobs1 <- apply(post1, 1, FUN = function(x){ length(unique(x))})
# cols <- grepl("sex", colnames(samples)) & !grepl("psex", colnames(samples))
# post.sex1 <- samples[-(1:5000),cols]
# post.z <- samples[-(1:5000), grepl("z\\[", colnames(samples))]
# NObsF1 <- do.call('c', lapply(1:nrow(post1),  FUN = function(x){ sum(post.sex1[x,unique(post1[x,])])}))
# NObs <- do.call('c', lapply(1:nrow(post1),  FUN = function(x){ sum(post.z[x,unique(post1[x,])])}))

# The full run...
samples4 <- runMCMC(Cmcmc, 60000, nburnin = 20000, nchains = 3, 
	thin = 1, inits = list(init4(), init4(), init4()))

# Compute the number of actually observed fisher to compare with SCR model in paper.
post1 <- samples4[[1]][,grep("ID", colnames(samples4[[1]]))]
post2 <- samples4[[2]][,grep("ID", colnames(samples4[[2]]))]
post3 <- samples4[[3]][,grep("ID", colnames(samples4[[3]]))]
Nobs1 <- apply(post1, 1, FUN = function(x){ length(unique(x))})
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
save(out4, file = "../output/fisherModel4.Rda")
load("../output/fisherModel4.Rda")
plot(out4[, c("K", "KF", "KM")])

###################################
# Chandler and Royle Santiy Check:
###################################
Rmodel <- nimbleModel(Model4, constants4, data4, inits = init4())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'psex', 'sex', 'z'))

# conf$removeSamplers('sigma')
# conf$addSampler(target = 'sigma', scalarComponents = TRUE,
	# type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))
# conf$removeSamplers('lambda')
# conf$addSampler(target = 'lambda', scalarComponents = TRUE,
	# type = 'RW', silent = TRUE, control = list(adaptive = FALSE, scale = 0.1))
# Use a block update on locations. Saves time.
conf$removeSamplers('X')
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
	type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 1.5))
# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

conf$removeSamplers('ID')
conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# The full run...
samples4 <- runMCMC(Cmcmc, 60000, nburnin = 20000, nchains = 3, 
	thin = 1, inits = list(init4(), init4(), init4()))

# Compute the number of actually observed fisher to compare with SCR model in paper.
post1 <- samples4[[1]][,grep("ID", colnames(samples4[[1]]))]
post2 <- samples4[[2]][,grep("ID", colnames(samples4[[2]]))]
post3 <- samples4[[3]][,grep("ID", colnames(samples4[[3]]))]
Nobs1 <- apply(post1, 1, FUN = function(x){ length(unique(x))})
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
save(out4, file = "../output/fisherModel4CR.Rda")







































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
	type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 1.5))
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
# plot(out[,c('N', 'D', 'sigma', 'lambda')])
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
	type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 1.4))
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


