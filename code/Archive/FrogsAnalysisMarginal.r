#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################
#################################

library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)
library(coda)
library(ggplot2)
library(combinat)

source("NimbleFunctions.R")
source("SimData.R")


load("../data/stacked-lightfooti.Rdata")

inits <- function(){
	ID = 1:n
    p <- runif(1, 0.1, 0.5)
	z = c(rep(1, n), rbinom(M-n, 1, p))
	lambda = runif(1, 0.1, 2)
	sigma = runif(1, 5, 10)
	lam0 = runif(1,1,10)

	dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
	pkj <- (1-exp(-lam0*exp(-dmask2/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x)%*%t(capt) + log(1-x)%*%t(1-capt))})
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
			  
	for(i in 1:M){
		if(sum(ID == i) == 0) next;
		if(sum(ID == i) == 1) pID <- panimal[ID == i, ]
		if(sum(ID == i) > 1) pID <- colSums(panimal[ID == i, ])
		mpt <- sample(ncol(panimal), 1, prob = exp(pID))
		X[i,] <- mask[mpt,]
	}
	sigmatoa = runif(1, 0.8, 1)
	list(
        lambda = lambda,
        psi = p,
        sigma = sigma,
		sigmatoa = sigmatoa,
		lam0 = lam0,
		X=X,
		z=z
    )
}

initsTRUE <- function(){
	ID <- capt.all$bincapt[, 7]
	ID <- ID[keep]
	ID <- as.integer(as.factor(ID))
    p <- runif(1, 0.1, 0.5)
	z = c(rep(1, max(ID)), rbinom(M-max(ID), 1, p))
	lambda = 0.3
	sigma = 2.3
	lam0 = 5.8

	dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
	pkj <- (1-exp(-lam0*exp(-dmask2/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x + .Machine$double.eps)%*%t(capt) + log(1-x + .Machine$double.eps)%*%t(1-capt))})
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
			  
	for(i in 1:M){
		if(sum(ID == i) == 0) next;
		if(sum(ID == i) == 1) pID <- panimal[ID == i, ]
		if(sum(ID == i) > 1) pID <- colSums(panimal[ID == i, ])
		mpt <- sample(ncol(panimal), 1, prob = exp(pID))
		X[i,] <- mask[mpt,]
	}
	
	sigmatoa = 0.002

	list(
        lambda = lambda,
        psi = p,
        sigma = sigma,
		sigmatoa = sigmatoa,
		lam0 = lam0,
		X=X,
		z=z
    )
}


xlim <- range(mask[,1])
ylim <- range(mask[,2])
area <- diff(xlim)*diff(ylim)
toa <- capt.all$toa
capt <- capt.all$bincapt[,1:6]
tmin <- apply(toa, 1, max)
keep <- which(tmin < 1200)
toa <- toa[keep,]
ID <- capt.all$bincapt[, 7]
ID <- ID[keep]
ID <- as.integer(as.factor(ID))
IDBen <- ID
capt <- capt[keep,]

# Constants:
M <- 200
nu <- 330
J <- nrow(traps)
n <- nrow(capt)
Time <- 30
mint <- min(toa[toa != 0])
toa<- toa - mint + 1	# That add one is to make sure they can't go negative for time of calling.
toa <- toa*capt
tmink <- tmin[keep] - mint

# tmp <- do.call('rbind', permn(c(0,1,0,1,0,0)))
# rowSums(tmp)

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
	y = capt,
	toa = toa,
	area = area)

data <- list(
	one = 1,
	ones = rep(1, nrow(capt)),
	z = rep(NA, M)
)

codeMarg <- nimbleCode({
    lambda ~ dunif(0, 10) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # Prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
	#tautoa ~ dunif(10,10000000)
	sigmatoa ~ dunif(0,1) # 1/sqrt(tautoa)
	lam0 ~ dunif(0, 20)
    for(i in 1:M) {
        z[i] ~ dbern(psi)
        X[i, 1] ~ dunif(xlim[1], xlim[2])
        X[i, 2] ~ dunif(ylim[1], ylim[2])
		# X[i, 1:2] ~ dX()
        d2[i,1:J] <- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
		expTime[i,1:J] <- sqrt(d2[i,1:J])/nu
        # pkj[i,1:J] <- lam0*exp(-d2[i,1:J]*tau2)
        pkj[i,1:J] <- (1-exp(-lam0*exp(-d2[i,1:J]*tau2)))
        # Hazard rate for animal across all traps.
        Hk[i] <- (1-prod(1-pkj[i,1:J]))*lambda*Time
		Hkz[i] <- Hk[i]*z[i]
    }	
    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
		ones[i] ~ dMarg(toa = toa[i,1:J], y = y[i, 1:J], J=J, size = trials[1:J], pkj = pkj, sd = sigmatoa, expTime = expTime, M = M)
    }
	p <- exp(-sum(Hkz[1:M]))*lambda^n_obs
    one ~ dbern(p)
    # Predicted population size
    Nhat <- sum(z[1:M])
	D <- Nhat/area*10000
})

Rmodel <- nimbleModel(codeMarg, constants, data, inits = initsTRUE())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'sigmatoa', 'lam0', 'Nhat', 'D', 'z', 'X'))
conf$removeSamplers('X')
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE)

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(10000)
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples)
plot(out[, c("lambda", "sigmatoa")])
plot(out[, c("sigma", "lam0")])
plot(out[, c("D", "Nhat")])
summary(out[, c("D", "Nhat")])
summary(out[, c("lambda", "sigmatoa")])
summary(out[, c("sigma", "lam0")])

