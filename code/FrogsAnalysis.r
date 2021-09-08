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
		ID = ID,
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
		ID = ID,
		z=z
    )
}

code <- nimbleCode({
    lambda ~ dunif(0, 10) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # Prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
	sigmatoa ~ dunif(0,1)
	# sigmatoa <- 0.001
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
        # Bernoulli capture history for each call that depends on ID
		y[i,1:J] ~ dbinom_vector(size = trials[1:J], pkj[ID[i],1:J])
		# Time of arrival, depends on which traps actually recorded it.
		toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i],1:J], sd = sigmatoa, y = y[i,1:J])
		ID[i] ~ dID()
    }
	p <- exp(-sum(Hkz[1:M]))*lambda^n_obs
    one ~ dbern(p)
    # Predicted population size
    Nhat <- sum(z[1:M])
	D <- Nhat/area*10000
})

xlim <- range(mask[,1])
ylim <- range(mask[,2])
area <- diff(xlim)*diff(ylim)
toa <- capt.all$toa
tmin <- apply(toa, 1, max)
keep <- which(tmin < 1200)
toa <- toa[keep,]
ID <- capt.all$bincapt[, 7]
ID <- ID[keep]
ID <- as.integer(as.factor(ID))
IDBen <- ID
capt <- capt.all$bincapt[,1:6]
capt <- capt[keep,]

# Constants:
M <- 200
nu <- 330
J <- nrow(traps)
n <- nrow(capt)
Time <- 30
mint <- min(toa[toa != 0])
toa[toa != 0]<- toa[toa != 0] - mint + 1	# That add one is to make sure they can't go negative for time of calling.
tmink <- tmin[keep] - mint

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
	area = area)

data <- list(
	one = 1,
    y = capt,
	toa = toa,
	z = rep(NA, M),
	ID = rep(NA, nrow(capt))
	# ID = IDBen
)

Rmodel <- nimbleModel(code, constants, data, inits = initsTRUE())

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('sigma', 'lambda', 'sigmatoa', 'lam0', 'Nhat', 'D', 'ID', 'z', 'X'))

conf$removeSamplers('X')
# for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'myX', control = list(xlim = xlim, ylim = ylim, J = nrow(traps)))
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE)

# conf$removeSamplers(c('sigma', 'lam0'))
# conf$addSampler(target = c('sigma', 'lam0'), type = 'RW_block', silent = TRUE)

conf$removeSamplers('sigmatoa')
conf$addSampler(target = 'sigmatoa', type = 'RW', control = list(log = TRUE))

# conf$printSamplers()

conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

conf$removeSamplers('ID')
# conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
# conf$addSampler('ID', type = 'myIDZASCR2', scalarComponents = TRUE, control = list(M = M, J = J))

# conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(10000)
mvSamples <- Cmcmc$mvSamples
samples <- as.matrix(mvSamples)
out <- mcmc(samples[-(1:2500),])
plot(out[, c("lambda", "sigmatoa", "Nhat")])
plot(out[, c("sigma", "lam0")])

# Demonstrate the problem:
post.x <- samples[-(1:5000),grep("X", colnames(samples))]
post.x1 <- post.x[,grep("1]", colnames(post.x))]
post.x2 <- post.x[,grep("2]", colnames(post.x))]

post.id <- samples[-(1:5000),grep("ID", colnames(samples))]

ID <- capt.all$bincapt[, 7]
ID <- ID[keep]
ID <- as.integer(as.factor(ID))
x1 <- data.frame(x = post.x1[cbind(1:nrow(post.id), post.id[,1])], y= post.x2[cbind(1:nrow(post.id), post.id[,1])])
x14 <-  data.frame(x = post.x1[cbind(1:nrow(post.id), post.id[,14])], y= post.x2[cbind(1:nrow(post.id), post.id[,14])])

# This is two obvious detections that should be matched and it is working great.
ggplot(data = data.frame(traps), aes(x=x,y=y)) + geom_point(shape = 4) + 
	theme_classic() + geom_point(data = x1, aes(x=x, y=y), col = "red", alpha = 0.1) + 
	geom_point(data = x14, aes(x=x, y=y), col = "blue", alpha = 0.1) + 
	geom_point(data = data.frame(traps)[capt[1,] == 1, ], aes(x=x,y=y), shape = 2, col = "red", size= 3) +
	geom_point(data = data.frame(traps)[capt[14,] == 1, ], aes(x=x,y=y), shape = 3, col = "blue", size= 3)	
sum(post.id[,1] == post.id[,14])/nrow(post.id)


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
