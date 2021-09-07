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

## Load the from Augustine:
load("../data/bobcats.Rdata")
secrmask <- data.frame(bobcats$mask)
names(secrmask) <- c("x", "y")
secrtraps <- bobcats$traps
habitatMask <- convertMask(secrmask, secrtraps, plot = TRUE)
habMat <- habitatMask$habMat
traps <- habitatMask$trapMat
upperlimit <- habitatMask$upperLimit
scaledMask <- habitatMask$scaledMask
Time <- 187
M <- 150
marks <- bobcats$mark
IDbobcat <- marks
IDbobcat[grep("L|R", marks)] <- NA
IDbobcat <- as.numeric(IDbobcat)
kobs <- max(IDbobcat, na.rm = TRUE)
J <- nrow(traps)
y <- apply(bobcats$capt, 1, FUN = function(x){which(x == 1)})
Time <- 187
mustlink <- matrix(0, nrow = length(y), ncol = length(y))+diag(1, length(y))
cannotlink <- matrix(0, nrow = length(y), ncol = length(y))
n <- length(y)
for(i in 1:n)
{
	marki <- bobcats$mark[i]
	LR <- gsub("[0-9]", "", marki)
	num <-  gsub("[[:alpha:]]", "", marki)
	ml <- bobcats$mark == marki
	mustlink[ml,i] <- 1
	mustlink[i,ml] <- 1	
	if(LR == ""){
		cannotlink[!ml,i] <- 1
		cannotlink[i,!ml] <- 1		
	}else{
		cannotlink[!ml & grepl(LR, bobcats$mark), i] <- 1
		cannotlink[i,!ml & grepl(LR, bobcats$mark)] <- 1
		cannotlink[grep("L|R", bobcats$mark, invert = TRUE), i] <- 1
		cannotlink[i, grep("L|R", bobcats$mark, invert = TRUE)] <- 1		
	}
}

inits <- function(){
    p <- runif(1, 0.1, 0.7)
	id.known <- factor(bobcats$mark)
	id.lr <- as.integer(factor(grep("L|R", id.known, value = TRUE)))
	ID <- IDbobcat
	ID[is.na(IDbobcat)] <- id.lr + kobs
	ID[!is.na(IDbobcat)] <- NA	# True Values should be NA...
	list(
        lambda = runif(1, 0.1, 2),
        psi = p,
        sigma = runif(1, 0.1, 1.5),
        X = scaledMask[sample(nrow(scaledMask), M),],
        ID = ID,
        z = c(rep(NA, kobs), rep(1,max(id.lr)), rep(0, M-kobs-max(id.lr)))
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
	z =  c(rep(1, max(IDbobcat, na.rm = TRUE)), rep(NA, M - max(IDbobcat, na.rm = TRUE))),
	ID = IDbobcat,
	habMat = habMat)

Rmodel <- nimbleModel(code, constants, data, inits = inits())

conf <- configureMCMC(Rmodel)

conf$setMonitors(c('sigma', 'lambda', 'psi', 'Nhat'))

conf$removeSamplers('X')
# for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'myX', control = list(xlim = limits$xlim, ylim = limits$ylim, J = nrow(traps)))
for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), type = 'RW_block', silent = TRUE) #, control = list(adaptive = FALSE)

conf$removeSamplers('z')
# Careful how to add sampler back!!
conf$addSampler('z[16:150]', type = 'myBinary', scalarComponents = TRUE)
# conf$printSamplers("z")

conf$removeSamplers('ID')
# conf$addSampler('ID', type = 'myIDZ', scalarComponents = TRUE, control = list(M = M))
# conf$addSampler('ID', type = 'myCategorical', scalarComponents = TRUE, control = list(M = M))
mark <- unique(bobcats$mark)
mark <- mark[grep("L|R", mark)]
for(i in 1:length(marks)){
	add <- which(bobcats$mark == mark[i])
	add.names <- paste0("ID[",add, "]")
	conf$addSampler(target = add.names, type = 'mySPIM', scalarComponents = FALSE, control = list(M = M, cannotlink = cannotlink))
}
# conf$printSamplers("ID")

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Cmcmc$run(10000)
mvSamples <- Cmcmc$mvSamples
samps <- as.matrix(mvSamples)
samps <- cbind(samps, samps[,"sigma"]*attr(habitatMask, "pixelWidth"))
colnames(samps)[ncol(samps)] <- "sigma_scaled"
samps.mcmc <- mcmc(samps[-(1:5000),c("sigma", "lambda", "Nhat", 'psi', 'sigma_scaled')])
plot(samps.mcmc)
summary(samps.mcmc)


post.x <- samps[-(1:5000),grep("X", colnames(samps))]
post.x1 <- post.x[,grep("1]", colnames(post.x))]
post.x2 <- post.x[,grep("2]", colnames(post.x))]
post.id <- samps[-(1:5000),grep("ID", colnames(samps))]
x1 <- data.frame(x = post.x1[cbind(1:nrow(post.id), post.id[,85])], y= post.x2[cbind(1:nrow(post.id), post.id[,85])])
allx <- cbind(post.x1[,20], post.x2[,20])
ggplot(data = data.frame(traps), aes(x=X,y=Y)) + geom_point(shape = 4) + 
	theme_classic() + geom_line(data = x1, aes(x=x, y=y), col = "red", alpha = 0.1)
ggplot(data = data.frame(traps), aes(x=X,y=Y)) + geom_point(shape = 4) + 
	theme_classic() + geom_line(data = data.frame(allx), aes(x=X1, y=X2), col = "red", alpha = 0.1)

sum(post.id[,1] == post.id[,14])/nrow(post.id)
