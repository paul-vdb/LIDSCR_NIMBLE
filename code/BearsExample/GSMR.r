#################################
# gSMR as a Marked Poisson Process
#################################
# setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/BearsExample")
library(sp)
library(coda)
library(raster)
library(nimble)
library(nimbleSCR)

source("../Functions/NimbleFunctions.R")
source("../Functions/GSMRSIM.R")



df <- fnc.create.SMR.data(N = 100, K.trap = 5, K.camera = 3, sigma = 0.04, g0trap = 0.5, lam0 = 0.7, 
                          n.trap = 16, n.camera = 100, trap.design = 'linear',
                          xlims = 0:1, ylims = 0:1, xlim.R = c(0.2, 0.8), ylim.R = c(0.2, 0.8), obsmod = "pois", 
                          n.collar = 1000, nlocs = 0) 
M <- 300
y.trap <- df$y.trap
y.camera.marked <- df$y.camera.marked 
n.camera.unmarked <- df$n.camera.unmarked

# Process the camera known ID detections.
omega <- NULL
ID <- NULL
for(i in 1:nrow(y.camera.marked))
{
	capti <- apply(y.camera.marked[i,,], 2, FUN = function(x){which(x > 0)})
	if(is.list(capti)) {
		capti <- do.call("c",capti)
	}else{
		capti <- as.numeric(capti)
	}
	omega <- c(omega, capti)
	ID <- c(ID, rep(i, length(capti)))
}

# Now process the unknown ID detections.
for(i in 1:ncol(n.camera.unmarked))
{
	capti <- which(n.camera.unmarked[,i] > 0)
	omega <- c(omega, capti)
	ID <- c(ID, rep(NA, length(capti)))
}


X.trap <- df$X.trap
X.camera <- df$X.camera

M <- 300
y <- array(0, dim = c(M, 16, 5))
y[1:nrow(y.trap), ,] <- y[1:nrow(y.trap), ,]  + y.trap

n_obs <- length(omega)
kmark <- nrow(y.camera.marked)

gSMR <- nimbleCode({
	sigma ~ dunif(0,1000) # uninformative prior
	psi ~ dbeta(1,1)
	lambda ~ dunif(0,20)
	g0 ~ dbeta(1,1)
	
	for(i in 1:M){
		z[i] ~ dbern(psi)
		X[i,1]~dunif(xlim[1],xlim[2])
		X[i,2]~dunif(ylim[1],ylim[2])
		
		d2.trap[i,1:J]<- (X[i,1]-X.trap[1:Jtrap,1])^2 + (X[i,2]-X.trap[1:Jtrap,2])^2
		d2.camera[i,1:J]<- (X[i,1]-X.camera[1:Jcamera,1])^2 + (X[i,2]-X.camera[1:Jcamera,2])^2
		p[i,1:J]<- g0*exp(-d2.trap[i, 1:J]/(2*sigma*sigma))
		lam[i,1:J] <- lambda*exp(-d2.camera[i, 1:Jcamera]/(2*sigma*sigma))
		y[i,1:J] ~ dbinom_vector(size = trials[1:Jtrap], prob = p[i,1:J])

		# Counting process for camera traps
		Hk[k] <- sum(pkj[k,1:Jcamera])*Time*lambda
		pkz[k] <- exp(-Hk[k]*z[k])	# Only put z here for purposes of node dependence and speed.
		zones[k] ~ dbern(pkz[k])
	}
	for(i in 1:n_obs)
	{
		pID[i] <- (obsMark[i] == animalMark[i])
        omega[i] ~ dTrap(p = pkj[1:M, 1:Jcamera], ID = ID[i], pID = pID[i])
		ID[i] ~ dID(lam = lambda)	# Dummy distribution to declare this as stochastic and mulitply by lambda.		
	}

    # Predicted population size
    N <- sum(z[1:M])
})

constants4 <- list(
    Jtrap = nrow(X.trap),
    Jcamera = nrow(X.camera),
    xlim = c(0,1),
    ylim = c(0,1),
    X.trap = X.trap,
	X.camera = X.camera,
    Time = 3,
    M = M,
    n_obs = length(omega),
	n_collar = 14,
	obsMark = !is.na(ID)*1,
	AnimalMark = c(rep(1, kmark), rep(0, M-kmark))
)

data4 <- list(
    zones = rep(1, M),
    omega = omega,
	z =  c(rep(1, kmark), rep(NA, M-kmark)),
	ID = ID,
	y = y
)

# Need to fully initialize because I haven't implemented a random generator distribution for dID()...
# Bit of a pain to make sure we match sex and collar correctly.
inits <- function(){
	g0 <- rbeta(1,1,1)
	lambda <- runif(1, 0.1, 1)
	sigma <- runif(1, 0.1, 1)
	X <- cbind(runif(M, xlim[1], xlim[2]), 
			  runif(M, ylim[1], ylim[2]))
	d2 <- t(apply(X, 1, FUN = function(x){(x[1] - X.camera[,1])^2 + (x[2] - X.camera[,2])^2}))
	hkj <- exp(-d2/(2*sigma^2))
	sexCollar <- c(rep(0, 5), rep(1, 9))
	sex <- c(rep(NA, 14), rbinom(M-14, size = 1, prob = psex))
	ID <- numeric(length(omega))
	z <- c(rep(NA, 14), rep(0, M-14))
	z[ID[ID > 14]] <- 1
	psi <- rbeta(1, sum(z,na.rm = TRUE) + 14, M - sum(1-z, na.rm = TRUE))	# NA inits...
	list(
		lambda = lambda,
		sigma = sigma,
		psi = psi,
		X = X,
		z = z,
		ID = ID,
		sex = sex,
		psex = psex
    )
}

###################################
# Chandler and Royle Sampler:
###################################
Rmodel <- nimbleModel(Model4, constants4, data4, inits = init4())
conf <- configureMCMC(Rmodel)
conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D', 'ID', 'psex', 'sex', 'z'))

conf$removeSamplers('X')
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
	type = 'RW_block', silent = TRUE, control = list(adaptive = FALSE, scale = 1.5))
# Optimized z sampler
conf$removeSamplers('z')
conf$addSampler('z', type = 'myBinary', scalarComponents = TRUE)

conf$removeSamplers('ID')
conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# The full run...
samples4 <- runMCMC(Cmcmc, niter = 100000, nburnin = 40000, nchains = 3, 
	thin = 1, inits = list(init4(), init4(), init4()))
