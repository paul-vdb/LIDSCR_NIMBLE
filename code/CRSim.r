#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################

library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)
library(coda)

source("NimbleFunctions.R")
source("SimData.R")

inits <- function(){
    p <- runif(1, 0.1, 0.7)
    K <- rbinom(1, M, p)
	list(
        lambda = runif(1, 0.1, 2),
        psi = p,
        sigma = runif(1, 0.1, 1.5),
        X = cbind(runif(M, limits[['xlim']][1], limits[['xlim']][2]), 
                  runif(M, limits[['ylim']][1], limits[['ylim']][2])),
		ID = sample(K, nrow(sim.dat), replace = TRUE),
        z = c(rep(1,K), rep(0, M-K))
    )
}

## Load the from Chandler and Royle:
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
sim.dat <- simSCR(N = N, sigma = sigma, lambda = lambda, StudyPeriod = StudyPeriod, traps, limits)
# Constants:
M <- 200
J <- nrow(traps)
Time <- StudyPeriod
xlim <- limits[['xlim']]
ylim <- limits[['ylim']]
traps <- traps

# Data:
y <- sim.dat$trap_obs
n <- nrow(sim.dat)
K <- length(table(sim.dat$ID))

code <- nimbleCode({
    lambda ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
    for(i in 1:M) {
        z[i] ~ dbern(psi)
		# I'm not convinced on the independence sampler for X as you implemented Daniel...
		# I think it could be a really good idea for weird state spaces and we should keep it in the back
		# pocket but for the simulations I don't think it's necessary.
        X[i, 1] ~ dunif(xlim[1], xlim[2])
        X[i, 2] ~ dunif(ylim[1], ylim[2])
        d2[i,1:J] <- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
        hkj[i,1:J] <- exp(-d2[i,1:J]*tau2)*lambda
        # Hazard rate for animal across all traps.
        Hk[i] <- sum(hkj[i,1:J])*Time
		Hkz[i] <- Hk[i]*z[i]
    }
    # Total thinning for all animals and traps.
    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # trap probability given ID:
        # This one can certainly be just the ones trick for trap y[i].
		pobs[i] <- hkj[ID[i], y[i]] + 0.0000000001
        ones[i] ~ dbern(pobs[i])
		ID[i] ~ dID()
    }
	p <- exp(-sum(Hkz[1:M]))
    one ~ dbern(p)
    # Predicted population size
    Nhat <- sum(z[1:M])
})

constants <- list(
    J = J,
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    Time = StudyPeriod,
    M = M,
    n_obs = nrow(sim.dat),
	y = y
	)

data <- list(
    one = 1,
    ones = rep(1, nrow(sim.dat)),
	# z = c(rep(1, K), rep(NA, M-K)),
	z =  rep(NA, M),
	ID = rep(NA, nrow(sim.dat))
)

Rmodel <- nimbleModel(code, constants, data, inits = inits())

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('sigma', 'lambda', 'psi', 'Nhat'))

conf$removeSamplers('X')
# for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'myX', control = list(xlim = limits$xlim, ylim = limits$ylim, J = nrow(traps)))
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE))

conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

conf$removeSamplers('ID')
conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
# conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))

# conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Cmcmc$run(5000)
# mvSamples <- Cmcmc$mvSamples
# samps <- as.matrix(mvSamples)
# samps.mcmc <- mcmc(samps)
# plot(samps.mcmc)

samples <- runMCMC(Cmcmc, 30000, nburnin = 10000, 
	nchains = 3, thin = 1, inits = list(inits(), inits(), inits()))
	
# samples <- runMCMC(Cmcmc, 30000, nburnin = 10000, 
	# nchains = 1, thin = 1, inits = inits())

out <- mcmc.list(list(as.mcmc(samples[[1]]), as.mcmc(samples[[2]]), as.mcmc(samples[[3]])))
plot(out[,c("Nhat", "sigma", "lambda", "psi")])
plot(out.cr[,c("Nhat", "sigma", "lambda", "psi")])


out <- as.mcmc(samples)
plot(out[,c("Nhat", "sigma", "lambda", "psi")])	
