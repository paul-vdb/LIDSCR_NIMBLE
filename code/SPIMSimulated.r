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
library(ggplot2)


source("NimbleFunctions.R")
source("SimData.R")

## Load the from Ben Augustine:
side <- 8
coords <- seq(1, side, length=side)
traps <- cbind(x=rep(coords, each=side), y=rep(coords, times=side))
buffer <- 2
limits <- list('xlim' = c(min(traps[,1]-buffer), max(traps[,1]+buffer)), 
			   'ylim' = c(min(traps[,2]-buffer), max(traps[,2]+buffer)))
area <- diff(limits[['xlim']])*diff(limits[['ylim']])
lambda <- lambdaTrue <- 0.5
sigma <- sigmaTrue <- 0.5
N <- NTrue <- 48
StudyPeriod <- 5
sim.dat <- simSCR(N = N, sigma = sigma, lambda = lambda, StudyPeriod = StudyPeriod, traps, limits)

# Now add right or left flank information:
y <- sim.dat$trap_obs
n <- length(y)
mark <- paste0(sim.dat$ID, c("L", "R")[rbinom(n, 1, 0.5) + 1])

# Constants:
M <- 200
J <- nrow(traps)
Time <- StudyPeriod
xlim <- limits[['xlim']]
ylim <- limits[['ylim']]
traps <- traps

# Process marks
mustlink <- diag(1, n)
cannotlink <- matrix(0, nrow = n, ncol = n)
for(i in 1:n)
{
	marki <- mark[i]
	LR <- gsub("[0-9]", "", marki)
	num <-  gsub("[[:alpha:]]", "", marki)
	ml <- mark == marki
	mustlink[ml,i] <- 1
	mustlink[i,ml] <- 1	
	if(LR == ""){
		cannotlink[!ml,i] <- 1
		cannotlink[i,!ml] <- 1		
	}else{
		cannotlink[!ml & grepl(LR, mark), i] <- 1
		cannotlink[i,!ml & grepl(LR, mark)] <- 1
		cannotlink[grep("L|R", mark, invert = TRUE), i] <- 1
		cannotlink[i, grep("L|R", mark, invert = TRUE)] <- 1		
	}
}


inits <- function(){
    p <- runif(1, 0.1, 0.7)
	ID <- as.integer(as.factor(mark))
	list(
        lambda = runif(1, 0.1, 2),
        psi = p,
        sigma = runif(1, 0.1, 1.5),
        X = cbind(runif(M, xlim[1],xlim[2]),runif(M, ylim[1], ylim[2])),
        ID = ID,
        z = c(rep(1, max(ID)), rep(0, M-max(ID)))
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
    Time = Time,
    M = M,
    n_obs = length(y),
	y = y)

data <- list(
    one = 1,
    ones = rep(1, length(y)),
	z =  rep(NA, M),
	ID = rep(NA, n))

Rmodel <- nimbleModel(code, constants, data, inits = inits())

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('sigma', 'lambda',  'psi', 'Nhat', 'ID'))

conf$removeSamplers('X')
# for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'myX', control = list(xlim = limits$xlim, ylim = limits$ylim, J = nrow(traps)))
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE) #, control = list(adaptive = FALSE)

conf$removeSamplers('z')
# Careful how to add sampler back!!
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)
# conf$printSamplers("z")

conf$removeSamplers('ID')
# conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
# conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
marktypes <- unique(mark)
for(i in 1:length(marktypes)){
	add <- which(mark == marktypes[i])
	add.names <- paste0("ID[",add, "]")
	conf$addSampler(target = add.names, type = 'mySPIM', scalarComponents = FALSE, control = list(M = M, cannotlink = cannotlink))
}
# conf$printSamplers("ID")
# conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(10000)
mvSamples <- Cmcmc$mvSamples
samps <- as.matrix(mvSamples)
post.id <- samps[-(1:5000),grep("ID", colnames(samps))]
NActive <- apply(post.id, 1, FUN = function(x){ length(unique(x))})
hist(NActive)
samps.mcmc <- mcmc(samps[-(1:5000),c("sigma", "lambda", "Nhat", 'psi')])
plot(samps.mcmc)

id.tmp <- post.id[2500,]
indx <- which(id.tmp == id.tmp[86])
bobcats$mark[indx]
y[indx]

# samps <- cbind(samps, samps[,"sigma"]*attr(habitatMask, "pixelWidth"))
# colnames(samps)[ncol(samps)] <- "sigma_scaled"
# samps.mcmc <- mcmc(samps[-(1:5000),c("sigma", "lambda", "Nhat", 'psi', 'sigma_scaled')])
# plot(samps.mcmc)
# summary(samps.mcmc)
# plot(samps.mcmc[,c("lambda", "sigma_scaled", "Nhat")])


samples <- runMCMC(Cmcmc, niter = 10000, nburnin = 5000, nchains = 3, thin = 1, inits = list(inits(), inits(), inits()))
out <- mcmc.list(list(as.mcmc(samples[[1]]), as.mcmc(samples[[2]]), as.mcmc(samples[[3]])))
plot(out[, c("Nhat", "sigma_scaled", "lambda")])
summary(out)

indx <- grep('ID', colnames(samps))
id.out <- samps[, indx]
z.out <- samps[,grep('z', colnames(samps))]
apply(id.out, 2, FUN = function(x){(sweep(id.out, 2, x) == 0)}

tmp <- sweep(id.out, 1, id.out[,30]) == 0
colSums(tmp)
bobcats$mark[30:34]
which(id.out[,c("ID[30]")] ==  id.out[,c("ID[31]")])
cannotlink[30,31:34]

no_link <- sum(cannotlink[model[['ID']] == i, nodeIndex])
table(z.out[cbind(1:nrow(samps), id.out[,72])])


# Single Side SCR0:
indxL <- !grepl("R", bobcats$mark)
yL <- y[indxL]
IDL <- as.integer(as.factor(bobcats$mark[indxL]))
constants <- list(
    J = J,
    upperlimit = upperlimit,
    traps = traps, 
    Time = Time,
    M = M,
    n_obs = length(yL),
	y = yL,
	pixelWidth = attr(habitatMask, "pixelWidth"),
	captureType = captureType[indxL]
	)

data <- list(
    one = 1,
    ones = rep(1, length(yL)),
	OK = rep(1, M),
	z =  c(rep(1, max(IDL, na.rm = TRUE)), rep(NA, M - max(IDL, na.rm = TRUE))),
	ID = IDL,
	habMat = habMat)



Rmodel <- nimbleModel(code, constants, data, inits = inits())

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('sigma', 'lambda',  'psi', 'Nhat', 'sigma_scaled', 'ID'))

conf$removeSamplers('X')
# for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'myX', control = list(xlim = limits$xlim, ylim = limits$ylim, J = nrow(traps)))
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE) #, control = list(adaptive = FALSE)

conf$removeSamplers('z')
# Careful how to add sampler back!!
conf$addSampler(paste0('z[', max(IDL) + 1, ':', M, ']'), type = 'myBinary', scalarComponents = TRUE)
# conf$printSamplers("z")

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(30000)
mvSamples <- Cmcmc$mvSamples
samps <- as.matrix(mvSamples)
samps.mcmc <- mcmc(samps[-(1:5000),c("sigma", "lambda[1]", "lambda[2]", "Nhat", 'psi', 'sigma_scaled')])
plot(samps.mcmc)
