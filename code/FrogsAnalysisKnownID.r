#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################

library(coda)
library(foreach)
library(doParallel)
library(nimble)
library(nimbleSCR)


load("stacked-lightfooti.Rdata")

#############################################
# Part 1:
# Simulate some ASCR data in continous time:
#############################################
# Simulate some data:
simASCR <- function(N = 50, sigma = 0.5, sigma_toa = 0.01, g0 = 1, lambda = 0.5, StudyPeriod = 25, traps, limits, psex = 0.5, nu = 330)
{
    locs <- cbind(x = runif(N, limits[['xlim']][1], limits[['xlim']][2]), 
                  y = runif(N, limits[['ylim']][1], limits[['ylim']][2]))
    J <- nrow(traps)
    obs <- data.frame()
	capthist <- toa <- data.frame()
    ID = 1
    for(i in 1:N)
    {
        d2 <- (locs[i,1] - traps[,1])^2 + (locs[i,2] - traps[,2])^2
        pkj <- g0*exp(-d2/(2*sigma^2))
		p. <- 1-prod(1-pkj)
        ti <- cumsum(rexp(1000, lambda))   #Simulate detection times.
        ti <- ti[ti < StudyPeriod]
        nk <- length(ti)
		if(nk == 0) next; 
        # Now assign those detection times to a trap.
		capt <- do.call("rbind", lapply(1:nk, FUN = function(x){rbinom(J,1,pkj)}))
		keep <- rowSums(capt) != 0
        capt <- capt[keep, ]
		if(sum(keep) == 0) next; 
		if(any(is.na(ti[keep]))) stop;
        obs <- rbind(obs, data.frame('t_obs' = ti[keep],'ID' = ID,
												'sex' = rbinom(1, 1, psex),
												'tru_x' = as.numeric(locs[i,1]), 'tru_y' = as.numeric(locs[i,2])))
		capthist <- rbind(capthist, capt)
		
		toa.i <- do.call("rbind", lapply(ti[keep], FUN = function(x){x + sqrt(d2)/nu + rnorm(J,0,sigma_toa)}))
		toa <- rbind(toa, toa.i)
		
        ID <- ID + 1
    }
    list(capt = capthist, toa = toa, obs = obs)
}


###
# Samplers to speed things up.
##

sampler_myBinary <- nimbleFunction(
    name = 'sampler_myBinary',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
        if((model[[target]] == 1) & (any(model[['ID']] == nodeIndex))) return()
        currentLogProb <- model$getLogProb(calcNodes)
        model[[target]] <<- 1 - model[[target]]
        otherLogProb <- model$calculate(calcNodes)
        acceptanceProb <- 1/(exp(currentLogProb - otherLogProb) + 1)
        jump <- (!is.nan(acceptanceProb)) & (runif(1,0,1) < acceptanceProb)
        if(jump) {
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)


sampler_myCategorical <- nimbleFunction(
    name = 'sampler_myCategorical',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        k <- length(model$getParam(target, 'prob'))
        probs <- numeric(k)
        logProbs <- numeric(k)
    },
    run = function() {
        currentValue <- model[[target]]
        logProbs[currentValue] <<- model$getLogProb(calcNodes)
        for(i in 1:k) {
            if(i != currentValue) {
                if(model[['z']][i] == 0) {
                    logProbs[i] <<- -Inf
                } else {
                    model[[target]] <<- i
                    logProbPrior <- model$calculate(target)
                    if(logProbPrior == -Inf) {
                        logProbs[i] <<- -Inf
                    } else {
                        if(is.nan(logProbPrior)) {
                            logProbs[i] <<- -Inf
                        } else {
                            logProbs[i] <<- logProbPrior + model$calculate(calcNodesNoSelf)
                            if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
                        }
                    }
                }
            }
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
            model[[target]] <<- newValue
            model$calculate(calcNodes)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

sampler_myTOA <- nimbleFunction(
    name = 'sampler_myTOA',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    },
    run = function() {
		n <- length(model[['ID']])
		J <- ncol(model[['y']])
		ssq <- 0
		ntot <- 0
		for(i in 1:n)
		{
			m <- sum(y)
			if(m != 1)
			{
				y <- model[['y']][i,1:J]
				x <- model[['toa']][i,1:J]
				tdiff <- y*(x - model[['expTime']][i,1:J])
				etdiff <- sum(tdiff)/m
				ssq <- ssq + sum((tdiff - etdiff)^2)
				ntot <- ntot + (1-m)
			}
		}
			
		model[[target]] <<- rgamma(1, shape = 0.01 - ntot/2, scale = 0.01 + 2/ssq))
		model$calculate(calcNodes)
		nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
		nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
		nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        },
    methods = list( reset = function() { } )
)


# Sampler to condition on my capture history for time of arrival and to speed things up.
dnorm_vector <- nimbleFunction(
  run = function( x = double(1),
                  mean = double(1),
                  sd = double(0),
				  y = double(1),
                  log = integer(0, default = 0)
                  ) {
    returnType(double(0))
    logProb <- sum(y*dnorm(x, mean = mean, sd = sd, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  })

rnorm_vector <- nimbleFunction(
  run = function( n = integer(0, default = 1),
                  mean = double(1),
                  sd = double(0),
				  y = double(1)
  ) {
    returnType(double(1))
    return(rnorm(length(y), mean = mean, sd = sd))
  })

  # Sampler without known detection times:
dnorm_vector_marg <- nimbleFunction(
  run = function( x = double(1),
                  mean = double(1),
                  sd = double(0),
				  y = double(1),
                  log = integer(0, default = 0)
                  ) {
    returnType(double(0))
	# This is not really a normal distribution but a small correction and it's proportional.
	# Fingers crossed...
	m <- sum(y)
	if(m == 1){
		if(log) return(0) else return(1)
	}
	tdiff <- y*(x - mean)
	etdiff <- sum(tdiff)/m
    logProb <- sum(y*dnorm(tdiff, mean = etdiff, sd = sd, log = TRUE)) + m*log(sd) + (1-m)*log(sd) # Normal correction factor.
    if(log) return(logProb) else return(exp(logProb))
  })

rnorm_vector_marg <- nimbleFunction(
  run = function( n = integer(0, default = 1),
                  mean = double(1),
                  sd = double(0),
				  y = double(1)
  ) {
    returnType(double(1))
    return(rnorm(length(y), mean = mean, sd = sd))
  })


# registerDistributions(
    # list(dnorm_vector = list(BUGSdist = 'dnorm_vector()',
                   # types = c('value = double(1)'))))


inits <- function(){
    p <- runif(1, 0.1, 0.2)
    # K <- rbinom(1, M, p)
	K <- max(ID)
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
	d2 <- t(apply(X, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))

	list(
        lambda = lambda,
        psi = p,
        sigma = sigma,
		sigmatoa = runif(1, 0.0001, 0.01),
		lam0 = lam0,
		X=X#,
		# tcall = abs(tmink - sqrt(d2[ID,1])/nu)
        # ID = sample(K, n, replace = TRUE)
    )
}

# inits <- function(){
    # p <- runif(1, 0.1, 0.2)
    ## K <- rbinom(1, M, p)
	# K <- max(ID)
	# X <- cbind(runif(M, xlim[1], xlim[2]), 
			  # runif(M, ylim[1], ylim[2]))		  
	# d2 <- t(apply(X, 1, FUN = function(x){(traps[,1]-x[1])^2 + (traps[,2] - x[2])^2}))

	# list(
        # lambda = runif(1,1,10),
        # psi = p,
        # sigma = runif(1, 5, 10),
		# sigmatoa = runif(1, 1, 2),
		# lam0 = runif(1, 5, 20),
		# X=X,
		# tcall = abs(tmink - sqrt(d2[ID,1])/nu)
        ## ID = sample(K, n, replace = TRUE)
    # )
# }


code <- nimbleCode({
    lambda ~ dunif(0,10) # Detection rate at distance 0
    psi ~ dbeta(1, 1)      # Prior on data augmentation bernoulli vec.
    sigma ~ dunif(0, 10)	# Now the prior is directly on sigma to be consistent with literature.
    tau2 <- 1/(2*sigma^2)
	# tautoa ~ dgamma(shape = 0.01, rate = 0.01)
	sigmatoa ~ dunif(0, 1)# 1/sqrt(tautoa)
	lam0 ~ dunif(0, 20)
    for(i in 1:M) {
        z[i] ~ dbern(psi)
        X[i, 1] ~ dunif(xlim[1], xlim[2])
        X[i, 2] ~ dunif(ylim[1], ylim[2])
        d2[i,1:J] <- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
		expTime[i,1:J] <- sqrt(d2[i,1:J])/nu
        # pkj[i,1:J] <- exp(-d2[i,1:J]*tau2)*z[i]
        pkj[i,1:J] <- (1-exp(-lam0*exp(-d2[i,1:J]*tau2)))*z[i]		
        # Hazard rate for animal across all traps.
        pk[i] <- (1-prod(1-pkj[i,1:J]))
    }
    # Total thinning for all animals and traps.
    p <- sum(pk[1:M])
	
    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
        # Bernoulli capture history for each call that depends on ID
		y[i,1:J] ~ dbinom_vector(size = trials[1:J], pkj[ID[i],1:J])
		# Time of arrival, depends on which traps actually recorded it.
		toa[i, 1:J] ~ dnorm_vector_marg(mean = expTime[ID[i],1:J], sd = sigmatoa, y = y[i,1:J])
		# tcall[i] ~ dunif(0, Time)
    }
    p1 <- exp(-lambda*p*Time)*lambda^n
	one ~ dbern(p1)
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
capt <- capt.all$bincapt[,1:6]
capt <- capt[keep,]

# Constants:
M <- 100
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
	ID = ID,
	area = area)

data <- list(
    n = nrow(capt),
    y = capt,
	toa = toa,
	# tcall = rep(NA, nrow(capt)),
	z = c(rep(1, max(ID)), rep(NA, M - max(ID))),
	one = 1
)

Rmodel <- nimbleModel(code, constants, data, inits = inits())
Rmodel$calculate()
Rmodel$simulate()

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('sigma', 'lambda', 'sigmatoa', 'lam0', 'Nhat', 'D', 'X', 'expTime'))

conf$removeSamplers('X')
# for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'myX', control = list(xlim = limits$xlim, ylim = limits$ylim, J = nrow(traps)))
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE)

# conf$printSamplers()

# conf$removeSamplers('z')
# for(i in (max(ID)+1):M) conf$addSampler(paste0('z[', i, ']'), type = 'myBinary', scalarComponents = TRUE)

# conf$removeSamplers('ID')
# conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000, nchains = 3, thin = 3, inits = list(inits(), inits(), inits()))
# samples <- runMCMC(Cmcmc, niter = 5000, nburnin = 1000, nchains = 1, thin = 10, inits =inits())

out <- mcmc.list(list(as.mcmc(samples[[1]]), as.mcmc(samples[[2]]), as.mcmc(samples[[3]])))
plot(out[,1:3])
dev.new()
plot(out[,4:6])

plot(out[,'sigmatoa'])
summary(out[,'sigmatoa'])

save(out, file = "frogsMarked.Rda")

# out <- mcmc(samples)
plot(out[,which(colnames(out[[1]]) %in% c("Nhat", "sigma", "sigmatoa", "psi"))])

library(ggplot2)
library(viridis)
post.x <- samples[[1]][,grep("X", colnames(samples[[1]]))]
post.x1 <- post.x[,grep("1]", colnames(post.x))]
post.x2 <- post.x[,grep("2]", colnames(post.x))]

ID <- capt.all$bincapt[, 7]
ID <- ID[keep]
ID <- as.integer(as.factor(ID))
x1 <- data.frame(x = post.x1[,1], y= post.x2[,1])
x2 <-  data.frame(x = post.x1[,2], y= post.x2[,2])
x3 <-  data.frame(x = post.x1[,3], y= post.x2[,3])
x4 <-  data.frame(x = post.x1[,4], y= post.x2[,4])
x5 <-  data.frame(x = post.x1[,5], y= post.x2[,5])

cols <- viridis(max(ID))

p <- ggplot(data = data.frame(traps), aes(x=x,y=y)) + geom_point(shape = 4) + 
	theme_classic() 
	
for( i in 1:max(ID))
{
	p <- p + geom_point(data = data.frame(x = post.x1[,i], y= post.x2[,i]), aes(x=x, y=y), alpha = 0.1, col = cols[i])
	meanx <- mean(post.x1[,i])
	meany <-  mean(post.x2[,i])
	p <- p + geom_text(data = data.frame(x = meanx, y=meany, lab = paste(i)), aes(x = x, y=y, label = lab), col = 'red', size = 5)
}
plot(p)	




#########
# Plotting the results
#########

library(ggplot2)
library(coda)
library(lattice)
library(bayesplot)
library(gridExtra)
load("frogsMarked.Rda")
out.marked <- out
load("frogsUnmarked.Rda")
out.unmarked <- out
mcmc_dens(out.marked, pars = c("D"))
mcmc_dens(out.marked, pars = c("D"))
marked <-  rbind(as.data.frame(out.marked[[1]]), as.data.frame(out.marked[[2]]), as.data.frame(out.marked[[3]]))
unmarked <-  rbind(as.data.frame(out.unmarked[[1]]), as.data.frame(out.unmarked[[2]]), as.data.frame(out.unmarked[[3]]))



p1 <- ggplot(data = marked, aes(x = D)) + geom_density(alpha = 0.5, colour = "blue") + 
	geom_density(data = unmarked, aes(x = D), colour = "red", alpha = 0.5) + theme_bw() + 
	xlab("Individuals per Hectare")
	
p2 <- ggplot(data = marked, aes(x = sigma)) + geom_density(alpha = 0.5, colour = "blue") + 
	geom_density(data = unmarked, aes(x = sigma), colour = "red", alpha = 0.5) + theme_bw() + 
	xlab(expression(sigma))

p3 <- ggplot(data = marked, aes(x = lambda)) + geom_density(alpha = 0.5, colour = "blue") + 
	geom_density(data = unmarked, aes(x = lambda), colour = "red", alpha = 0.5) + theme_bw()+
	xlab(expression(mu))

p4 <- ggplot(data = marked, aes(x = lam0)) + geom_density(alpha = 0.5, colour = "blue") + 
	geom_density(data = unmarked, aes(x = lam0), colour = "red", alpha = 0.5) + theme_bw()	+ 
	xlab(expression(lambda[0]))

grid.arrange(p1, p2, p3, p4, ncol=2)	