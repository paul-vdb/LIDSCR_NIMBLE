#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################
setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/Bernoulli")
#################################
library(coda)
library(nimble)
library(nimbleSCR)
library(poisbinom)

source("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/Functions/NimbleFunctions.R")

load("../../data/MartenData.Rda")
xlim = marten.data$xlim
ylim = marten.data$ylim
area <- diff(xlim)*diff(ylim)/100	# Density reported per 100 sq km
traps <- marten.data$traps
obs <- marten.data$observation

M <- 200

y <- obs$TrapNumber

Model.a <- nimbleCode({
    lambda ~ dunif(0, 20)

	# Home range parameter m=0/f=1:
    sigma ~ dunif(0, 50)

	# convert now to make it simpler in half-normal computation
    tau2 <- 1/(2*sigma^2)

    psi ~ dbeta(1, 1)      # prior on data augmentation Bernoulli vec.

	# For the collared individuals, we know the sex too! It's actually observed.
	# So we have observed sex 1:14 and an additional covariate of collar.
	# As a result, we also have z = 1 observed for 1:14
	for(k in 1:M) {
		z[k] ~ dbern(psi)

		X[k, 1] ~ dunif(xlim[1], xlim[2])
		X[k, 2] ~ dunif(ylim[1], ylim[2])
		
		d2[k,1:J] <- (X[k,1]-traps[1:J,1])^2 + (X[k,2]-traps[1:J,2])^2
		gkj[k,1:J] <- exp(-d2[k,1:J]*tau2)
		
		# Need probability of switching off/on an animal to a trap.
		pkj[k, 1:J] <- exp(-gkj[k,1:J]*lambda*z[k])
		p[k, 1:J] <- (1-pkj[k, 1:J])/pkj[k, 1:J]

		# Hazard rate for animal across all traps.
		Hk[k] <- sum(gkj[k,1:J])*lambda*z[k]		
		# Total zero process
		zeros[k] ~ dpois(Hk[k])
	}

    # Trap history model.
    # and unobserved animal ID.
    for(i in 1:n_obs) {
		# Trap Prob.
        omega[i] ~ dTrap(p = p[1:M, 1:J], ID = ID[i])
		ID[i] ~ dID(pID = z[1:M])	# Dummy distribution to declare this as stochastic and mulitply by lambda.		
    }
	
    # Predicted population size
    N <- sum(z[1:M])
	D <- N/area
})

constants <- list(
    J = nrow(traps),
    xlim = xlim,
    ylim = ylim,
    traps = traps, 
    M = M,
    n_obs = length(y),
	area = diff(xlim)*diff(ylim)/100
)

data <- list(
	omega = y,
    zeros = rep(0, M),
	z =  rep(NA, M),
	ID = rep(NA, length(y))
)

# Need to fully initialize because I haven't implemented a random generator distribution for dID()...
# Bit of a pain to make sure we match sex and collar correctly.
init <- function(){
	N <- floor(runif(1, 1, M/2))
	psi <- rbeta(1, N, M-N)	# NA inits...	
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 2, 4)
	X <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - traps[,1])^2 + (x[2] - traps[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	ID <- do.call('c',lapply(1:length(y), FUN = function(x){sample(1:M, 1, prob = hkj[1:M,y[x]])}))
	for(j in 1:nrow(traps))
	{
		bigI <- which(y==j)
		ids <- ID[bigI]
		indx <- duplicated(ids)
		if(sum(indx) == 0) next;
		for(k in 1:sum(indx))
		{
			prob <- hkj[1:M,j]
			prob[ids[indx]] <- 0
			ID[bigI[indx]] <- sample(1:M, 1, prob = prob)
		}
	}
	Hk <- rowSums(hkj)*lambda
	p <- exp(-Hk)*psi/(exp(-Hk)*psi + (1-psi))
	z <- rbinom(M, size = 1, prob=p)
	z[ID] <- 1
	list(
		lambda = lambda,
		sigma = sigma,
		psi = psi,
		X = X,
		z = z,
		ID = ID
		)
}

Rmodel <- nimbleModel(Model.a, constants, data, inits = init())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID'))
# Use a block update on locations. Saves time.
# Turn off adaptive samping and fix the scale of the sampler to something reasonable.
conf$removeSamplers('X')
for(i in 1:M){
	conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'RW_block', silent = TRUE, 
		control = list(scale = 0.5, adaptive = FALSE))
}

conf$removeSamplers('sigma')
conf$addSampler(target = 'sigma', 
		type = 'slice', silent = TRUE, control = list(adaptive = FALSE, scaleWidth = 0.2))	

# van Dam-Bates categorical sampler
conf$removeSamplers('ID')
# Chandler and Royle Alg. 1 sampler.
conf$addSampler('ID', type = 'myCategoricalBernoulli', scalarComponents = TRUE, 
	control = list(M = M, omega = y))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(10000)
mvSamplesa <- Cmcmc$mvSamples
samplesa <- as.matrix(mvSamplesa)
out.a <- mcmc(samplesa[-(1:5000),])
plot(out.a[,c('sigma', 'lambda', 'N', 'D')])
D <- out.a[, "N"]/diff(xlim)*diff(ylim)/100
plot(D)

plot(density(out.mpp[, 'N']), col = 'blue', main = "Comparing Methods", xlab="N", ylab = "Density")
lines(density(out.p[, 'N']))
lines(density(out.a[, 'N']), col = "red")
legend("topright",legend = c("Poisson-Binomial", "Poisson Approx", "Allocation Model"), col = c("blue", "black", "red"), lty = c(1,1,1))
abline(v = 15, col = 'black', lty = 'dashed')

plot(out.mpp[,c('sigma', 'lambda', 'N')])
plot(out.1[,c('sigma', 'lambda', 'N')])

par(mfrow = c(3,1))
plot(density(out.mpp[, 'sigma']), main = 'sigma')
lines(density(out.a[, 'sigma']), col = 'blue')
abline(v = 0.75, col = 'red')
legend('topright', legend = c("Poisson Binomial", "Allocation"), col = c("black", "blue"), lty = 1)

plot(density(out.mpp[, 'lambda']), main = 'lambda')
lines(density(out.a[, 'lambda']), col = 'blue')
abline(v = 2, col = 'red')

plot(density(out.mpp[, 'N']), main = 'N')
lines(density(out.a[, 'N']), col = 'blue')
abline(v = 30, col = 'red')


p <- rbeta(100, 1, 5)
nj <- 0:100
plot(dpoisbinom(nj, p))

approx.func <- function(p, iters = 1000)
{
	n <- NULL
	for(i in 1:iters)
	{
		n <- c(n, sum(rbinom(length(p), 1, prob = p)))
	}
	return(n)
}

tmp <- approx.func(p, 10000)
plot(dpoisbinom(nj, p), type = 'l')
lines(density(n), col = 'red')
