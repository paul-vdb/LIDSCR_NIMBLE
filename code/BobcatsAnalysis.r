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

## Load the from Augustine:
load("../data/bobcats.Rdata")
secrmask <- data.frame(bobcats$mask)
names(secrmask) <- c("x", "y")
secrtraps <- bobcats$traps
habitatMask <- convertMask(secrmask, secrtraps, plot = TRUE)
habMat <- habitatMask$habMat
traps <- habitatMask$trapMat
upperlimit <- habitatMask$upperLimit
Time <- 187
M <- 150
ID <- bobcats$mark
ID[grep("L|R", bobcats$mark)] <- NA
ID <- as.numeric(ID)
J <- nrow(traps)
y <- apply(bobcats$capt, 1, FUN = function(x){which(x == 1)})
Time <- 187
mustlink <- matrix(0, nrow = length(y), ncol = length(y))+diag(1, length(y))
cannotlink <- matrix(0, nrow = length(y), ncol = length(y))

for(i in 1:length(y))
{
	marki <- bobcats$mark[i]
	LR <- gsub("[0-9]", "", marki)
	num <-  gsub("[[:alpha:]]", "", marki)
	ml <- bobcats$mark == marki
	mustlink[ml,i] <- 1
	if(LR == ""){
		cannotlink[!ml,i] <- 1
	}else{
		cannotlink[!ml & grepl(LR, bobcats$mark)] <- 1
	}
}


inits <- function(){
    p <- runif(1, 0.1, 0.7)
    kobs <- max(ID, na.rm = TRUE)
	id.known <- as.integer(factor(bobcats$mark))
	id.lr <- as.integer(factor(id.known[is.na(ID)]))
	K <- max(id.known)
	ID[is.na(ID)] <- id.lr + kobs
	list(
        lambda = runif(1, 0.1, 2),
        psi = p,
        sigma = runif(1, 0.1, 1.5),
        X = cbind(runif(M, 1, upperlimit[1]), 
                  runif(M, 1, upperlimit[2])),
        ID = ID,
        z = c(rep(1,K), rep(0, M-K))
    )
}

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
        X[i, 1] ~ dunif(1, upperlimit[1])
        X[i, 2] ~ dunif(1, upperlimit[2])
		pOK[i] <- habMat[trunc(X[i, 1]), trunc(X[i, 2])] # habitat check
		OK[i] ~ dbern(pOK[i]) # OK[i] = 1, the ones trick		
		
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
    upperlimit = upperlimit,
    traps = traps, 
    Time = Time,
    M = M,
    n_obs = length(y),
	y = y
	)

data <- list(
    one = 1,
    ones = rep(1, length(y)),
	OK = rep(1, M),
	z =  c(rep(1, max(ID, na.rm = TRUE)), rep(NA, M - max(ID, na.rm = TRUE))),
	ID = ID,
	habMat = habMat)

Rmodel <- nimbleModel(code, constants, data, inits = inits())

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('sigma', 'lambda', 'psi', 'Nhat'))

conf$removeSamplers('X')
# for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'myX', control = list(xlim = limits$xlim, ylim = limits$ylim, J = nrow(traps)))
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE) #, control = list(adaptive = FALSE)

conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

conf$removeSamplers('ID')
# conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
# conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
marks <- unique(bobcats$mark)
marks <- marks[grep("L|R", marks)]
for(i in 1:length(marks)){
	add <- which(bobcats$mark == marks[i])
	add.comma <- paste(add, collapse = ",")
	conf$addSampler(target = paste0('ID[c(', add.comma, ')]'), type = 'mySPIM', scalarComponents = TRUE, control = list(M = M, cannotlink = cannotlink))
}
# conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(5000)
mvSamples <- Cmcmc$mvSamples
samps <- as.matrix(mvSamples)
samps.mcmc <- mcmc(samps)
plot(samps.mcmc)
