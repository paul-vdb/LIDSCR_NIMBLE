#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################

library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)

source("NimbleFunctions.R")

## Load the from Chandler and Royle:

nopaDat <-
structure(list(y = structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0,
1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1,
1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0,
0, 1, 2, 1, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1,
1, 0, 2, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 2, 1, 1, 1, 0, 0, 0,
0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1,
1, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0), .Dim = c(105L, 3L), .Dimnames = list(
    c("4005", "3055", "3005", "2055", "2005", "1055", "1005",
    "4100", "3150", "3100", "2150", "2100", "1150", "1100", "4105",
    "3155", "3105", "2155", "2105", "1155", "1105", "4200", "3250",
    "3200", "2250", "2200", "1250", "1200", "4205", "3255", "3205",
    "2255", "2205", "1255", "1205", "4300", "3350", "3300", "2350",
    "2300", "1350", "1300", "4305", "3355", "3305", "2355", "2305",
    "1355", "1305", "4400", "3450", "3400", "2450", "2400", "1450",
    "1400", "4405", "3455", "3405", "2455", "2405", "1455", "1405",
    "4500", "3550", "3500", "2550", "2500", "1550", "1500", "4505",
    "3555", "3505", "2555", "2505", "1555", "1505", "4600", "3650",
    "3600", "2650", "2600", "1650", "1600", "4605", "3655", "3605",
    "2655", "2605", "1655", "1605", "4700", "3750", "3700", "2750",
    "2700", "1750", "1700", "4705", "3755", "3705", "2755", "2705",
    "1755", "1705"), c("v1", "v2", "v3"))), X = structure(c(5,
5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8,
8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10,
10, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 13,
13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15,
15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17,
17, 17, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19,
5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9,
10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7,
8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5,
6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10,
11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8,
9, 10, 11, 5, 6, 7, 8, 9, 10, 11), .Dim = c(105L, 2L), .Dimnames = list(
    NULL, c("x", "y"))), xSide = 24, ySide = 16, M = 100, nTraps = 105L,
    nReps = 3L), .Names = c("y", "X", "xSide", "ySide", "M",
"nTraps", "nReps"))


######################
# Process data:
######################
Time <- StudyPeriod <- nopaDat$nReps
J <- nopaDat$nTraps
traps <- nopaDat$X
nj <- nopaDat$y
xlim <- c(1, nopaDat$xSide)+1
ylim <- c(1, nopaDat$ySide)+1
y <- NULL
for(i in 1:J){
	for(j in 1:StudyPeriod)
	y <- c(y, rep(i, nj[i,j]))
}
n <- length(y)
table(y)

M <- 300

code <- nimbleCode({
    lambda ~ dunif(0, 20) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 50)	# Now the prior is directly on sigma to be consistent with literature.
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
    n_obs = length(y),
	y = y
	)

data <- list(
    one = 1,
    ones = rep(1, length(y)),
	# z = c(rep(1, K), rep(NA, M-K)),
	z =  rep(NA, M),
	ID = rep(NA, length(y))
)

inits <- function(){
    p <- runif(1, 0.1, 0.7)
    K <- rbinom(1, M, p)
	list(
        lambda = runif(1, 0.1, 2),
        psi = p,
        sigma = runif(1, 0.1, 1.5),
        X = cbind(runif(M, xlim[1], xlim[2]), 
                  runif(M, ylim[1], ylim[2])),
		ID = sample(K, length(y), replace = TRUE),
        z = c(rep(1,K), rep(0, M-K))
    )
}

Rmodel <- nimbleModel(code, constants, data, inits())

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

samples <- runMCMC(Cmcmc, 30000, nburnin = 10000, 
	nchains = 3, thin = 1, inits = list(inits(), inits(), inits()))
out <- mcmc.list(list(as.mcmc(samples[[1]]), as.mcmc(samples[[2]]), as.mcmc(samples[[3]])))
plot(out[,c("Nhat", "sigma", "lambda", "psi")])

