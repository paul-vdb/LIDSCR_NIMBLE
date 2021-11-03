#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################
# setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/FrogExample")
library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)
library(coda)
library(ggplot2)
library(secr)
library(ascr)

source("../Functions/NimbleFunctions.R")
source("../Functions/SimData.R")

load("../../data/stacked-lightfooti.Rdata")

# Change the mask to finer resolution:
# traps2 <-convert.traps(lightfooti$traps)
traps2 <- convert.traps(traps)
mask <- make.mask(traps = traps2, buffer = 15, spacing = 0.2, type = "trapbuffer")
A <- attr(mask, "area")
mask <- as.matrix(mask)
attr(mask, "area") <- A
attr(mask, "buffer") <- 15


inits <- function(){
	ID <- numeric(length(occ))
	ID[occ == 1] <- 1:sum(occ == 1)
	ID[occ == 2] <- 1:sum(occ == 2)
	z <- matrix(0, ncol = 2, nrow = M)
	z[ID[occ==1],1] <- 1
	z[ID[occ==2],2] <- 1
    p <- runif(1, 0.1, 0.3)
	lambda = runif(1, 0.1, 2)
	sigma = runif(1, 3, 5)
	g0 = runif(1,1,10)

	dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
	pkj <- (1-exp(-g0*exp(-dmask2/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x)%*%t(capt) + log(1-x)%*%t(1-capt))})
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
        psi = p,
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
			Hk[k,v] <-(1-prod(1-pkj[k,1:J,v]))*lambda*Time
			pkz[k,v] <- exp(-Hk[k,v]*z[k,v])
			zones[k,v] ~ dbern(pkz[k,v])
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

xlim <- range(mask[,1])
ylim <- range(mask[,2])
area <- diff(xlim)*diff(ylim)
toa <- capt.all$toa
capt <- capt.all$bincapt[,1:6]
tmin <- apply(toa, 1, max)
occ <- 1+(tmin > 1200)
ID <- capt.all$bincapt[, 7]
ID <- as.integer(as.factor(ID))
IDBen <- ID

# indx2 <- rowSums(capt) > 1
# capt <- capt[indx2,] 
# toa <- toa[indx2,]
# occ <- occ[indx2]

# Constants:
M <- 200
nu <- 330
J <- nrow(traps)
n <- nrow(capt)
Time <- 30

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
	n_occ = 2,
	occ = occ,
	mi = rowSums(capt))

data <- list(
	zones = cbind(rep(1,M), rep(1,M)),
    y = capt,
	toa = toa,
	z = cbind(rep(NA, M), rep(NA, M)),
	ID = rep(NA, nrow(capt))
)

Rmodel <- nimbleModel(code, constants, data, inits = inits())

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('psi', 'sigma', 'lambda', 'sigmatoa', 'g0', 'N', 'D', 'ID'))

conf$removeSamplers('X')
for(v in 1:2){
	for(i in 1:M) {
		conf$addSampler(target = paste0('X[', i, ', 1:2,', v, ']'), type = 'RW_block', silent = TRUE, 
			control = list(scale = 2, adaptive = FALSE))
	}	
}

conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE, 
	control = list('Noccasion' = 2, 'IDoccasion' = occ))

conf$removeSamplers('ID')
# Sampler from Chandler and Royle 2013
conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Make sure it runs...
# Cmcmc$run(25000)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# out <- mcmc(samples[-(1:5000), c("psi", "sigma", "sigmatoa", "lambda", "g0","N[1]", "N[2]", "D")])
# plot(out, ask = TRUE)
# summary(out)[[1]]["D", "Mean"]*summary(out)[[1]]["lambda", "Mean"]*60
# quantile(out[, "D"]*out[, "lambda"]*60, c(0.0275, 0.975))

# post.id <- samples[-(1:5000),grep("ID", colnames(samples))]
# post.id1 <- post.id[,occ == 1]
# post.id2 <- post.id[,occ == 2]
# NActive <- apply(post.id, 1, FUN = function(x){ length(unique(x))})
# NActive1 <- apply(post.id1, 1, FUN = function(x){ length(unique(x))})
# NActive2 <- apply(post.id2, 1, FUN = function(x){ length(unique(x))})
# par(mfrow = c(2,1))
# hist(NActive1)
# abline(v = 14, col = 'red')
# hist(NActive2)
# abline(v = 11, col = 'red')

# par(mfrow = c(2,1))
# hist(samples[-(1:5000),"N[1]"])
# hist(samples[-(1:5000),"N[2]"])


# plot(density(as.numeric(out[,"D"])), main = "Frogs", xlab = "D (Ind/Ha)")
# abline(v = c(358), col = 'red')
# abline(v = c(240,534), col = 'grey')

mcmc.out <- runMCMC(Cmcmc, nburnin = 20000, niter = 50000, nchains = 3, 
	inits = list(inits(), inits(), inits()))	

out.list <- list()
for( i in 1:3 )
{
	post.id <- mcmc.out[[i]][,grep("ID", colnames(mcmc.out[[i]]))]
	post.id1 <- post.id[,occ == 1]
	post.id2 <- post.id[,occ == 2]
	NActive  <- apply(post.id, 1, FUN = function(x){ length(unique(x))})	
	NActive1 <- apply(post.id1, 1, FUN = function(x){ length(unique(x))})
	NActive2 <- apply(post.id2, 1, FUN = function(x){ length(unique(x))})
	out.list[[i]] <- mcmc(cbind(mcmc.out[[i]][, c("psi", "sigma", "sigmatoa", "lambda", "g0", "N[1]", "N[2]", "D")], K = NActive, K1 = NActive1, K2 = NActive2))
}
out <- as.mcmc.list(out.list)
summary(out)
plot(out, ask = TRUE)
# save(out, file = "../../output/FrogsLatentID2Session_HOLYSHITITWORKED.Rda")
save(out, file = "../../output/FrogsLatentID2Session.Rda")
load( "../../output/FrogsLatentID2Session.Rda")

###########
# Known ID ASCR from Stevenson 2020
###########

initsID <- function(){
    p <- runif(1, 0.1, 0.3)
	z <- matrix(rbinom(M*2, 1, p), ncol = 2, nrow = M)
	z[IDBen[occ==1],1] <- NA
	z[IDBen[occ==2],2] <- NA
	lambda = runif(1, 0.1, 2)
	sigma = runif(1, 3, 5)
	g0 = runif(1,1,10)

	dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
	pkj <- (1-exp(-g0*exp(-dmask2/(2*sigma^2))))
	panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x)%*%t(capt) + log(1-x)%*%t(1-capt))})
	X <- array(c(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]), 
		runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2])), c(M,2,2))		
	
	for(k in 1:2){
		for(i in 1:M){
			if(sum(IDBen[occ == k] == i) == 0) next;
			if(sum(IDBen[occ == k] == i) == 1) pID <- panimal[which(IDBen[occ == k] == i), ]
			if(sum(IDBen[occ == k] == i) > 1) pID <- colSums(panimal[which(IDBen[occ == k] == i), ])
			mpt <- sample(ncol(panimal), 1, prob = exp(pID))
			X[i,,k] <- mask[mpt,]
		}
	}

	list(
        lambda = lambda,
        psi = p,
        sigma = sigma,
		sigmatoa = runif(1, 0.01, 1),
		g0 = g0,
		X = X,
		z = z
    )
}

code.id <- nimbleCode({
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
		}
		# Total thinning for all animals and traps.
		H[v] <- sum(Hk[1:M,v])		
		# Predicted population size
		N[v] <- sum(z[1:M,v])
	}
	
    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # Bernoulli capture history for each call that depends on ID
		y[i,1:J] ~ dbinom_vector(size = trials[1:J], pkj[ID[i],1:J, occ[i]])
		# Time of arrival, depends on which traps actually recorded it.
		toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i], 1:J, occ[i]], sd = sigmatoa, y = y[i,1:J])
    }

	# Animal Process model:
	D <- psi*M/area*10000
	p1 <- exp(-sum(H[1:n_occ]))*lambda^n_obs
	one ~ dbern(p1)
})

IDBen[occ == 1] <- as.integer(as.factor(IDBen[occ == 1]))
IDBen[occ == 2] <- as.integer(as.factor(IDBen[occ == 2]))

IDBen[occ == 1] <- 1:sum(occ == 1)
IDBen[occ == 2] <- 1:sum(occ == 2)

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
	occ = occ,
	ID = IDBen
)

data.id <- list(
	one = 1,
    y = capt,
	toa = toa,
	z = cbind(c(rep(1, max(IDBen[occ == 1])), rep(NA, M-max(IDBen[occ == 1]))),
		c(rep(1, max(IDBen[occ == 2])), rep(NA, M-max(IDBen[occ == 2]))))
)

Rmodel <- nimbleModel(code.id, constants.id, data.id, inits = initsID())

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('sigma', 'lambda', 'sigmatoa', 'g0', 'N', 'D'))

conf$removeSamplers('X')
for(v in 1:2){
	for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2,', v, ']'), type = 'RW_block', silent = TRUE)
}

conf$removeSamplers('sigmatoa')
conf$addSampler(target = 'sigmatoa', type = 'RW', control = list(log = TRUE))

# conf$removeSamplers('z')
# conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE, 
	# control = list('Noccasion' = 2, 'IDoccasion' = occ))

conf$removeSamplers('sigmatoa')
conf$addSampler(target = 'sigmatoa', type = 'RW', control = list(log = TRUE))

# conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

mcmc.out.id <- runMCMC(Cmcmc, nburnin = 10000, niter = 30000, nchains = 3, 
	inits = list(initsID(), initsID(), initsID()))	

out.id <- as.mcmc.list(list(mcmc(mcmc.out.id[[1]])[, c("sigma", "sigmatoa", "lambda", "g0", "N[1]", "N[2]", "D")], 
	mcmc(mcmc.out.id[[2]])[, c("sigma", "sigmatoa", "lambda", "g0", "N[1]", "N[2]", "D")], 
	mcmc(mcmc.out.id[[3]])[, c("sigma", "sigmatoa", "lambda", "g0", "N[1]", "N[2]", "D")]))
summary(out.id)
plot(out.id, ask = TRUE)

save(out.id, file = "../../output/FrogsKnownID2Session.Rda")
load("../../output/FrogsKnownID2Session.Rda")

###################
library(lattice)
lapply(posterior, function(x) varnames(x) <- names)
out.D <- mcmc.list(list(mcmc(cbind(mcmc.out[[1]][,"D"], mcmc.out.id[[1]][,"D"])),
	mcmc(cbind(mcmc.out[[2]][,"D"], mcmc.out.id[[2]][,"D"])),
	mcmc(cbind(mcmc.out[[3]][,"D"], mcmc.out.id[[3]][,"D"]))))
varnames(out.D) <- c("LatentID","KnownID"))
xyplot(out.D)
	default.scales = list(y = list(relation = "free", limits = list(c(0,1.5), 'sigma' = c(0, 10), c(0, 300)))))
plot(out.D)
###################


############################
# Now run it with ASCR, 
# but we need to increase the mask:
############################
multi.capt <- list(list("bincapt" = capt[occ==1,], "toa" = toa[occ==1,]), list("bincapt" = capt[occ ==2,], "toa" = toa[occ==2,]))
ascr.res <- fit.ascr(capt = multi.capt, list(traps, traps), mask = list(mask, mask), 
			detfn = "hhn", survey.length = c(30, 30))
# citation("ascr")
summary(ascr.res)
confint(ascr.res)
#Output:
#           D      lambda0        sigma    sigma.toa 
# 9.200900e+02 4.966609e+00 2.541285e+00 3.486909e-04 



Cmcmc$run(10000)
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:5000),])
plot(out[, c("lambda", "sigmatoa")])
plot(out[, c("sigma", "g0")])
plot(out[, c("D", "Nhat")])
summary(out[, c("D", "Nhat")])
summary(out[, c("lambda", "sigmatoa")])
summary(out[, c("sigma", "g0")])

# Demonstrate the problem:
post.x <- samples[-(1:5000),grep("X", colnames(samples))]
post.x1 <- post.x[,grep("1]", colnames(post.x))]
post.x2 <- post.x[,grep("2]", colnames(post.x))]

post.id <- samples[-(1:5000),grep("ID", colnames(samples))]
NActive <- apply(post.id, 1, FUN = function(x){ length(unique(x))})
hist(NActive)

ID <- capt.all$bincapt[, 7]
ID <- ID[keep]
ID <- as.integer(as.factor(ID))
i1 <- 22
i2 <- 60
x1 <- data.frame(x = post.x1[cbind(1:nrow(post.id), post.id[,i1])], y= post.x2[cbind(1:nrow(post.id), post.id[,i1])])
x14 <-  data.frame(x = post.x1[cbind(1:nrow(post.id), post.id[,i2])], y= post.x2[cbind(1:nrow(post.id), post.id[,i2])])

# This is two obvious detections that should be matched and it is working great.
ggplot(data = data.frame(traps), aes(x=x,y=y)) + geom_point(shape = 4) + 
	theme_classic() + geom_point(data = x1, aes(x=x, y=y), col = "red", alpha = 0.1) + 
	geom_point(data = x14, aes(x=x, y=y), col = "blue", alpha = 0.1) + 
	geom_point(data = data.frame(traps)[capt[i1,] == 1, ], aes(x=x,y=y), shape = 2, col = "red", size= 3) +
	geom_point(data = data.frame(traps)[capt[i2,] == 1, ], aes(x=x,y=y), shape = 3, col = "blue", size= 3)	
sum(post.id[,1] == post.id[,61])/nrow(post.id)

# These are two not so obvious mathces that may actually not match 
# but they should certainly not look like they do below.
# This doesn't match with the posterior call location of an animal heard at a single trap at all.
# I can't make any sense of why this matching would be sooo bad.
x78 <- data.frame(x = post.x1[cbind(1:nrow(post.id), post.id[,78])], y= post.x2[cbind(1:nrow(post.id), post.id[,78])])
x85 <-  data.frame(x = post.x1[cbind(1:nrow(post.id), post.id[,85])], y= post.x2[cbind(1:nrow(post.id), post.id[,85])])

ggplot(data = data.frame(traps), aes(x=x,y=y)) + geom_point(shape = 4) + 
	theme_classic() + geom_point(data = x78, aes(x=x, y=y), col = "red", alpha = 0.1) + 
	geom_point(data = x85, aes(x=x, y=y), col = "blue", alpha = 0.1) +
	geom_point(data = data.frame(traps)[capt[78,] == 1, ], aes(x=x,y=y), shape = 2, col = "red", size= 3) +
	geom_point(data = data.frame(traps)[capt[85,] == 1, ], aes(x=x,y=y), shape = 3, col = "blue", size= 3)
sum(post.id[,78] == post.id[,85])/nrow(post.id)











# SINGLE SESSION FROG ASCR:
code <- nimbleCode({
    lambda ~ dunif(0, 10) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # Prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
	sigmatoa ~ dunif(0,1) #1/sqrt(tautoa)
	g0 ~ dunif(0, 50)
    for(i in 1:M) {
        z[i] ~ dbern(psi)
        X[i, 1] ~ dunif(xlim[1], xlim[2])
        X[i, 2] ~ dunif(ylim[1], ylim[2])
        d2[i,1:J] <- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
		expTime[i,1:J] <- sqrt(d2[i,1:J])/nu
        pkj[i,1:J] <- (1-exp(-g0*exp(-d2[i,1:J]*tau2)))
        # Hazard rate for animal across all traps.
        Hk[i] <-(1-prod(1-pkj[i,1:J]))*lambda*Time
		Hkz[i] <- Hk[i]*z[i]
    }	
    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # Bernoulli capture history for each call that depends on ID
		y[i,1:J] ~ dbinom_vector(size = trials[1:J], pkj[ID[i],1:J])
		# Time of arrival, depends on which traps actually recorded it.
		toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i],1:J], sd = sigmatoa, y = y[i,1:J])
		ID[i] ~ dID()
    }
	p <- exp(-sum(Hkz[1:M]))*lambda^n_obs
    one ~ dbern(p)
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area*10000
})
