#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################
setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/FrogExample")
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
			Hk[k,v] <-(1-prod(1-pkj[k,1:J,v]))*lambda*Time*z[k,v]
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
	zeros = cbind(rep(0,M), rep(0,M)),
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
		conf$addSampler(target = paste0('X[', i, ', 1:2,', v, ']'), 
			type = 'myJAM', silent = TRUE, 
			control = list(scale = 1, xlim = xlim, ylim = ylim, temp = 0.2, occ = 2))
	}
}

conf$removeSamplers('g0')
conf$addSampler(target = 'g0', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))		

conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE, 
	control = list('Noccasion' = 2, 'IDoccasion' = occ))

# conf$printSamplers()

conf$removeSamplers('ID')
# Sampler from Chandler and Royle 2013
conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Make sure it runs...
# Cmcmc$run(20000)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# samples <- cbind(samples, CD = samples[,"D"]*samples[, "lambda"]*60)
# out <- mcmc(samples[-(1:10000), c("psi", "sigma", "sigmatoa", "lambda", "g0","N[1]", "N[2]", "D", "CD")])
# summary(out)

# Plot the deteciton function:

# post.id <- samples[-(1:10000),grep("ID", colnames(samples))]
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
# summary(out)
save(out, file = "../../output/FrogResults/LIDASCR.Rda")

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

ID <- capt.all$bincapt[, 7]
IDBen[occ == 1] <- as.integer(as.factor(ID[occ == 1]))
IDBen[occ == 2] <- as.integer(as.factor(ID[occ == 2]))
K <- c(max(IDBen[occ == 1]), max(IDBen[occ == 2]))

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
	zeros = matrix(0, nrow = M, ncol = 2),
    y = capt,
	toa = toa,
	z = cbind(c(rep(1, K[1]), rep(NA, M-K[1])),
		c(rep(1, K[2]), rep(NA, M-K[2]))),
	ID = IDBen
)

Rmodel <- nimbleModel(code, constants.id, data.id, inits = initsID())

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('sigma', 'lambda', 'sigmatoa', 'g0', 'N', 'D'))

conf$removeSamplers('X')
for(v in 1:2){
	for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2,', v, ']'), type = 'RW_block', silent = TRUE)
}

conf$removeSamplers('g0')
conf$addSampler(target = 'g0', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.5))		

# conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Cmcmc$run(5000)
# mvSamples <- Cmcmc$mvSamples
# samples <- as.matrix(mvSamples)
# samples <- cbind(samples, CD = samples[,"D"]*samples[, "lambda"]*60)
# out.id <- mcmc(samples[-(1:10000),])
# plot(out.id[,c('D', 'sigma', 'lambda')])
# plot(out.id[,c('g0', 'sigmatoa')])
valueInCompiledNimbleFunction(Cmcmc$samplerFunctions[[381]], "scale")

mcmc.out.id <- runMCMC(Cmcmc, nburnin = 20000, niter = 50000, nchains = 3, 
	inits = list(initsID(), initsID(), initsID()))	

out.id <- as.mcmc.list(list(mcmc(mcmc.out.id[[1]])[, c("sigma", "sigmatoa", "lambda", "g0", "N[1]", "N[2]", "D")], 
	mcmc(mcmc.out.id[[2]])[, c("sigma", "sigmatoa", "lambda", "g0", "N[1]", "N[2]", "D")], 
	mcmc(mcmc.out.id[[3]])[, c("sigma", "sigmatoa", "lambda", "g0", "N[1]", "N[2]", "D")]))
# summary(out.id)
# plot(out.id, ask = TRUE)
# save(out.id, file = "../../output/FrogResults/IDASCR.Rda")
load("../../output/FrogResults/IDASCR.Rda")



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
# Detection function: Hazard halfnormal 
# Information types: Times of arrival
 
 # Parameters: 
            # Estimate Std. Error
# D         1.2489e+02    10.1940
# lambda0   6.5932e+00     1.0356
# sigma     2.2916e+00     0.0943
# sigma.toa 9.6946e-04     0.0001
# ---                            
                               
# esa.1     2.4556e-02     0.0009
# esa.2     2.4556e-02     0.0009
# > confint(ascr.res)
                 # 2.5 %       97.5 %
# D         1.049065e+02 1.448663e+02
# lambda0   4.563476e+00 8.622954e+00
# sigma     2.106834e+00 2.476330e+00
# sigma.toa 7.833584e-04 1.155556e-03

cdascr <- data.frame(cbind(Mean = summary(ascr.res)[[1]], confint(ascr.res)))
cdascr$variable <- c("CD", "g0", "sigma", "sigmatoa")
cdascr$Method <- "CD"
# cdascr["lambda0", "Mean"] <- cdascr["lambda0", "Mean"]
# save(cdascr, file = "../../output/FrogResults/ascr_mle.Rda")
load("../../output/FrogResults/ascr_mle.Rda")

library(reshape2)
library(ggplot2)
# Now plot results:

out.lid <- do.call("rbind", lapply(out, as.data.frame))
out.lid$Method = "LID"
out.lid$CD <- out.lid$D*out.lid$lambda
out.kid <- do.call("rbind", lapply(out.id, as.data.frame))
out.kid$Method = "ID"
out.kid$CD <- out.kid$D*out.kid$lambda
pars <- c("g0", "sigma", "lambda", "D", "CD", "sigmatoa")

output <- rbind(out.lid[,c(pars,"Method")], out.kid[,c(pars,"Method")])

all.out <- melt(output, id.vars = "Method")
all.out$value <- as.numeric(all.out$value)

save(all.out, file = "../../output/FrogResults/CombinedFrogMCMCOutput.Rda")

facet_names <- c('g[0]', 'sigma~(m)', 'lambda~(sec^{-1})', 'D~(Ind~Ha^{-1})', 'sigma[t]~(s)', "Call~Density~(Calls~Ha^-1~Second^-1)" )
names(facet_names) <- c("g0", "sigma", "lambda", "D", "sigmatoa", "CD")


ggplot(data = all.out, aes(y = value, x = Method)) + 
	facet_wrap(~variable, scale = "free", labeller = as_labeller(x = facet_names, label_parsed)) + 
	geom_point(data = cdascr, aes(x = Method, y = Mean)) +
	geom_errorbar(data = cdascr, aes(y = Mean, ymin = `X2.5..`, ymax = `X97.5..`), width = 0.2) +
	theme_classic() +
	geom_boxplot(width = 0.1, outlier.alpha = 0, fill = "black", colour = "white") + 
	geom_violin(alpha = 0, adjust = 2) + 
	ylab("") + xlab("")
ggsave("../../output/FrogResults/CombinedFrogMCMCOutput.png", dpi = 'print')

out.lid.l <- melt(out.lid)

ggplot(data = out.lid.l[out.lid.l$variable %in% c("K1", "K2"),], aes(x = value)) + 
	geom_histogram(fill = "black", colour = "black") +
	facet_wrap(~variable) +
	theme_classic() + 
	geom_vline(data = data.frame(variable = c("K1", "K2"), 
		value = c(14,11)), aes(xintercept = value), col = 'red')









###
# Posterior plots of location given parameter values:
###
library(ggplot2)
post.mask <- function(capt, toa = NULL, pars, Time = 30, mask, traps, detfn = 'hn', speed.sound = 330){
	d2 <- t(apply(mask, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	
	sigma <- pars['sigma']
	lambda <- pars['lambda']
	g0 <- pars['g0']
	sigmatoa <- pars['sigmatoa']	

	if(g0>0) detfn <- "hhn"

	if(detfn == "hhn")
	{
		pkj <- 1-exp(-g0*exp(-d2/(2*sigma^2)))
	}else{
		print("half-normal detection function")
		pkj <- g0*exp(-d2/(2*sigma^2))
	}
	
	l.tmp <- rowSums(log(1-pkj))
	p. <- 1-exp(l.tmp)

	# Time of arrival function:
	getTOA <- function(x, y, d2, nu = 330, sd = 0.001)
	{
		m <- sum(y)
		if(m == 1){
			logProb <- 0
		}else{
			logProb <- apply(sqrt(d2[,y==1])/nu, 1, FUN = function(z){
				td <- x[y==1] - z;
				etd <- sum(td)/m;
				sum((td - etd)^2)/(2*sd^2)
			})
			logProb <- (1-m)*log(sd) - logProb
		}
		return(logProb)
	}

	ll <- numeric(nrow(mask))
	k <- nrow(capt)
	if(is.null(dim(capt))) {
		k <- 1
		ll <- rowSums(t(apply(pkj, 1, FUN = function(x){dbinom(x = capt, size = 1, prob = x, log = TRUE)})))
		# ll <- ll + getTOA(toa, capt, d2, nu = speed.sound, sd = sigmatoa)
	}else{
		for(i in 1:k)
		{
			ll <- ll + rowSums(t(apply(pkj, 1, FUN = function(x){dbinom(x = capt[i,], size = 1, prob = x, log = TRUE)})))
			ll <- ll + getTOA(toa[i,], capt[i,], d2, nu = speed.sound, sd = sigmatoa)
		}
	}
	ll <- ll - lambda*Time*p. + k*log(lambda)
	return(exp(ll))
}

test.capt <- capt.all$bincapt[capt.all$bincapt[,7] == 1, 1:6]
test.toa <- capt.all$toa[capt.all$bincapt[,7] == 1,]
pars <- c("sigma" = 2.3, "lambda" = 0.3, "g0" = 5.75, "sigmatoa" = 0.001)
post.locs <- post.mask(capt = test.capt, toa = test.toa, pars, Time = 30, mask, traps, detfn = 'hn')

ggplot(data = data.frame(mask), aes(x = x, y = y, fill = post.locs)) + geom_tile() +
	geom_point(data = data.frame(traps), aes(x=x, y=y, fill = NULL), col = 'red') + 
	theme_classic()
ggplot(data = data.frame(mask), aes(x = x, y = y, z = post.locs)) + geom_contour() +
	geom_point(data = data.frame(traps), aes(x=x, y=y, z= NULL), col = 'red') + 
	theme_classic()
	