#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################
#################################
library(coda)
library(nimble)
library(nimbleSCR)

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

source("SimData.R")
load("../data/stacked-lightfooti.Rdata")

results <- foreach(h=1:100,
			.packages = c("coda", "nimbleSCR", "nimble"))%dopar%{
	source("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/NimbleFunctions.R")

	### Scenario 1:
	dat <- simASCR(N = 60, sigma = 2.3, sigma_toa = 0.00055, g0 = 5.75, lambda = 0.27, 
		StudyPeriod = 30, traps = traps, 
		limits = list(xlim = range(mask[,1]), ylim = range(mask[,2])))

	inits <- function(){
		ID = 1:n
		p <- runif(1, 0.1, 0.5)
		z = c(rep(1, n), rbinom(M-n, 1, p))
		lambda = runif(1, 0.1, 2)
		sigma = runif(1, 5, 10)
		g0 = runif(1,1,10)

		dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
		pkj <- (1-exp(-g0*exp(-dmask2/(2*sigma^2))))
		panimal <- apply(pkj, 1, FUN = function(x){colSums(log(x)%*%t(capt) + log(1-x)%*%t(1-capt))})
		X <- cbind(runif(M, xlim[1], xlim[2]), 
				  runif(M, ylim[1], ylim[2]))
				  
		for(i in 1:M){
			if(sum(ID == i) == 0) next;
			pID <- panimal[ID == i, ]
			mpt <- sample(ncol(panimal), 1, prob = exp(panimal[ID == i, ]))
			X[i,] <- mask[mpt,]
		}
		sigmatoa = runif(1, 0.8, 1)
		list(
			lambda = lambda,
			psi = p,
			sigma = sigma,
			sigmatoa = sigmatoa,
			g0 = g0,
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
		sigmatoa ~ dunif(0,1) #1/sqrt(tautoa)
		g0 ~ dunif(0, 50)
		for(i in 1:M) {
			z[i] ~ dbern(psi)
			X[i, 1] ~ dunif(xlim[1], xlim[2])
			X[i, 2] ~ dunif(ylim[1], ylim[2])
			# X[i, 1:2] ~ dX()
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

	xlim <- range(mask[,1])
	ylim <- range(mask[,2])
	area <- diff(xlim)*diff(ylim)
	toa <- dat$toa
	capt <- dat$capt
	ID <- dat$obs$ID

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
		area = area)

	data <- list(
		one = 1,
		y = capt,
		toa = toa,
		z = rep(NA, M),
		ID = rep(NA, nrow(capt))
	)

	Rmodel <- nimbleModel(code, constants, data, inits = inits())

	conf <- configureMCMC(Rmodel)

	conf$setMonitors(c('sigma', 'lambda', 'sigmatoa', 'g0', 'N', 'D', 'ID'))

	# conf$removeSamplers('X')
	# for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE, control = list(scale = 0.05, adaptive = FALSE))
	# for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'sampler_myX2', silent = TRUE, 
		# control = list(xlim = xlim, ylim = ylim, scale = 0.25, J = nrow(traps)))

	conf$removeSamplers('sigmatoa')
	conf$addSampler(target = 'sigmatoa', type = 'RW', control = list(log = TRUE))

	conf$removeSamplers('z')
	conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

	conf$removeSamplers('ID')
	conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))

	Rmcmc <- buildMCMC(conf)

	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

	Cmcmc$run(40000)
	mvSamples <- Cmcmc$mvSamples
	samples <- as.matrix(mvSamples)
	out <- mcmc(samples[-(1:10000),])
	save(out, file = paste0("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/ASCRSimulations/ASCR_Scenario_1_iter_", h, ".Rda"))
	summary(out)
}

# Demonstrate the problem:
post.x <- samples[-(1:5000),grep("X", colnames(samples))]
post.x1 <- post.x[,grep("1]", colnames(post.x))]
post.x2 <- post.x[,grep("2]", colnames(post.x))]

post.id <- samples[-(1:5000),grep("ID", colnames(samples))]
NActive <- apply(post.id, 1, FUN = function(x){ length(unique(x))})
hist(NActive)


x1 <- data.frame(x = post.x1[cbind(1:nrow(post.id), post.id[,1])], y= post.x2[cbind(1:nrow(post.id), post.id[,1])])
x14 <-  data.frame(x = post.x1[cbind(1:nrow(post.id), post.id[,2])], y= post.x2[cbind(1:nrow(post.id), post.id[,2])])

# This is two obvious detections that should be matched and it is working great.
ggplot(data = data.frame(traps), aes(x=x,y=y)) + geom_point(shape = 4) + 
	theme_classic() + geom_point(data = x1, aes(x=x, y=y), col = "red", alpha = 0.1) + 
	geom_point(data = x14, aes(x=x, y=y), col = "blue", alpha = 0.1) + 
	geom_point(data = data.frame(traps)[capt[1,] == 1, ], aes(x=x,y=y), shape = 2, col = "red", size= 3) +
	geom_point(data = data.frame(traps)[capt[2,] == 1, ], aes(x=x,y=y), shape = 3, col = "blue", size= 3)	
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


Rmodel$sigmatoa <- 0.0003
Rmodel$calculate()

target <- 'ID[78]'
logprobs <- numeric(M)
calcNodes <- Rmodel$getDependencies(target)
for(i in 1:M)
{
	Rmodel[[target]] <- i
	logprobs[i] <- Rmodel$calculate(calcNodes)
}
p <- exp(logprobs - max(logprobs))
ggplot(data = data.frame(traps), aes(x=x,y=y)) + geom_point(shape = 4) + 
	theme_classic() + geom_point(data = data.frame(Rmodel$X), 
		aes(x=X1, y=X2, colour = p))


# Is there ever a z = 0 with an ID.
ids <- samples[,grep("ID", colnames(samples))]
zs <- samples[,grep("z", colnames(samples))]
ones <- do.call('rbind', lapply(1:nrow(ids), FUN = function(x){zs[x, ids[x,]]}))

# Show Mask info:
ind <- 3
dmask2 <- t(apply(mask, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))
pkj <- (1-exp(-6.5*exp(-dmask2/(2*2.2^2))))
pcapt <- apply(pkj, 1, FUN = function(x){dbinom_vector(capt[ind,], size = rep(1,J),x)})
ptoa <- apply(sqrt(dmask2)/330, 1, FUN = function(x){dnorm_vector_marg(x = toa[ind, ], mean = x, sd = 0.001, y = capt[ind,1:J])})
p <- 1-rowProd(1-pkj)
ggplot(data = data.frame(mask), aes(x=x,y=y)) + geom_tile(aes(fill=ptoa*pcapt)) + 
	geom_point(data = data.frame(traps), aes(x=x,y=y), shape = 16, col = "red", size= 3) + 
	geom_point(data = data.frame(traps)[capt[ind,]==1,], aes(x=x,y=y), shape = 16, col = "green", size= 3) 
	# geom_point(data = x85, aes(x=x, y=y), col = "purple", alpha = 0.1)



d <- seq(0, 20, by = 0.1)
p <- (1-exp(-3.93*exp(-d^2/(2*2.76^2))))
p2 <- (1-exp(-7.5*exp(-d^2/(2*2.76^2))))
p3 <- (1-exp(-7.5*exp(-d^2/(2*2.2^2))))
plot(d, p, type = 'l')
lines(d, p2, type = 'l', col = 'red')
lines(d, p3, type = 'l', col = 'blue')



code_2 <- nimbleCode({
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
		p0[i] <- (1-prod(1-pkj[i,1:J]))
		for(j in 1:J) p1[j,i] <- dbinom_vector(x = mat1[j,], size = trials[1:J], prob = pkj[i,1:J], log = 0)
        Hk[i] <- (1-p0[i]-sum(p1[1:J,i]))*lambda*Time
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
