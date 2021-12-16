#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################
#################################
setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/FisherExample")

library(coda)
library(nimble)
library(nimbleSCR)
library(secr)

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

source("../Functions/SimData.R")
load("../../data/FisherData.Rda")
traps <- fisher.data$traps
xlim <- range(traps[,1]) + c(-6,6)
ylim <- range(traps[,2]) + c(-6,6)
J <- nrow(traps)
M <- 300
area <- diff(xlim)*diff(ylim)/100

results <- foreach(h=1:250,
			.packages = c("coda", "nimbleSCR", "nimble"))%dopar%{
	source("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/code/Functions/NimbleFunctions.R")

	### Scenario 1:
	dat <- simSCR(N = 60, NCollar = c(5, 9), sigma = 2.5, lambda = 0.15/30/0.4,
		StudyPeriod = 120, traps = traps, 
		limits = list(xlim = xlim, ylim = ylim), psex = 0.6)
	omega <- dat$trap
	nj <- numeric(nrow(traps))
	for(j in 1:nrow(traps)) nj[j] <- nrow(dat[dat$trap_obs == j,])

	# Now thin 40% of detections.
	dat <- dat[rbinom(nrow(dat), size = 1, prob = 0.4) == 1,]
	dat$occ <- dat$t_obs %/% 30 + 1
	M <- 300
	y <- array(0, dim = c(M, nrow(traps), 4))

	K <- length(unique(dat$ID))
	for(i in 1:K)
	{
		y[cbind(i,dat[dat$ID == i,]$trap_obs, dat[dat$ID == i,]$occ)] <- 1
	}
	y_all <- apply(y, 1:2, sum)

	SCR_bern <- nimbleCode({
		sigma ~ dunif(0,1000) # uninformative prior
		psi ~ dbeta(1,1)
		psex ~ dbeta(1,1)
		lambda ~ dunif(0,20)
		
		for(i in 1:M){
			z[i] ~ dbern(psi)
			X[i,1]~dunif(xlim[1],xlim[2])
			X[i,2]~dunif(ylim[1],ylim[2])
			d2[i,1:J]<- (X[i,1]-traps[1:J,1])^2 + (X[i,2]-traps[1:J,2])^2
			# Because detectability didn't vary per trap night we aggregated!
			# For comparision to Poisson we will use a hazard half normal detection function.
			# That way we are using the same lambda, encounter rate.
			# See Augustine paper as he does this as well.
			Hkj[i, 1:J] <- 30*lambda*exp(-d2[i, 1:J]/(2*sigma*sigma))*z[i]
			p[i,1:J]<- (1-exp(-Hkj[i,1:J]*0.4))
			#From Daniel Turek in nimbleSCR package. Fast binomial! Avoids loopin.
			y[i,1:J] ~ dbinom_vector(size = trials[1:J], prob = p[i,1:J])
		}
		# for( j in 1:J )
		# {
			# Hj[j] <- sum(Hkj[1:M, j]*4)
			# nj[j] ~ dpois(Hj[j])
		# }
		N <- sum(z[1:M])
		D <- N/area
	})

	constants<- list(
		J = nrow(traps),
		area = area,
		xlim = xlim,
		ylim = ylim,
		traps = traps, 
		M = M,
		trials = rep(4, nrow(traps))
		)
	data <- list(
		z =  c(rep(1, K), rep(NA, M-K)),
		y = y_all#,
		#nj = nj
		)

	inits <- function()
	{
		psi <- rbeta(1, 1+K, M-K)
		z <- numeric(M)
		z[1:K] <- NA
		X <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
		for( i in 1:K ){
			X[i,] <- as.numeric(traps[which.max(y_all[i,]),]) + rnorm(2,0,0.1)
		}
		list(X = X, z = z, sigma = runif(1,2, 3), lambda = runif(1,0.0005,0.001), psi = psi)
	}

	Rmodel <- nimbleModel(SCR_bern, constants, data, inits = inits())
	Rmodel$calculate()
	Rmodel$calculate("y[2,1:64]")
	conf <- configureMCMC(Rmodel)
	conf$setMonitors(c('sigma', 'lambda', 'psi', 'N', 'D'))

	# Use a block update on locations. Saves time.
	conf$removeSamplers('X')
	for(i in 1:M) conf$addSampler(target = paste0('X[', i, ', 1:2]'), 
		type = 'RW_block', silent = TRUE)

	Rmcmc <- buildMCMC(conf)
	Cmodel <- compileNimble(Rmodel)
	Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

	Cmcmc$run(20000)
	mvSamples <- Cmcmc$mvSamples
	samples <- as.matrix(mvSamples)
	out <- mcmc(samples[-(1:10000),])
	# plot(out[,c("lambda", "sigma","N", "D")])
	save(out, file = paste0("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/FisherSimulations/SCR_Scenario_2Genetic_iter_", h, ".Rda"))
	summary(out)
}

stopCluster(cl)



cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)

results <- foreach(h=1:250,
			.packages = c("secr"))%dopar%{
	### Scenario 1:
	dat <- simSCR(N = 60, NCollar = c(5, 9), sigma = 2.5, lambda = 0.15/30/0.4,
		StudyPeriod = 120, traps = traps, 
		limits = list(xlim = xlim, ylim = ylim), psex = 0.6)
	omega <- dat$trap
	dat$occ <- dat$t_obs %/% 30 + 1
		
	names(traps) <- c('x','y')
	traps <- cbind("trapID" = paste0("T",1:nrow(traps)), traps)
	traps <- read.traps(data = traps, detector = "proximity", trapID = "trapID")
	mesh <- make.mask(traps, buffer = 6)
	capt.dat <- data.frame("session" = 1, "ID" = dat$ID, "occasion" = dat$t_obs %/% 30 + 1, "trapID" = paste0("T",dat$trap_obs))
	capt.dat <- capt.dat[!duplicated(paste(capt.dat$ID, capt.dat$occasion, capt.dat$trapID)),]
	capt <- make.capthist(capt.dat, traps, fmt = "trapID", noccasions = 4)
	fit <- secr.fit(capthist = capt, mask = mesh, binomN = 1, detectfn="hhn", trace = FALSE)
	res1 <- summary(fit)$predicted
	res1$Scenario = "Unthinned"
	res1$parameter = rownames(res1)

	dat <- dat[rbinom(nrow(dat), size = 1, prob = 0.4) == 1,]
	capt.dat <- data.frame("session" = 1, "ID" = dat$ID, "occasion" = dat$t_obs %/% 30 + 1, "trapID" = paste0("T",dat$trap_obs))
	capt.dat <- capt.dat[!duplicated(paste(capt.dat$ID, capt.dat$occasion, capt.dat$trapID)),]
	capt <- make.capthist(capt.dat, traps, fmt = "trapID", noccasions = 4)
	fit <- secr.fit(capthist = capt, mask = mesh, binomN = 1, detectfn="hhn",trace = FALSE)
	res2 <- summary(fit)$predicted
	res2$Scenario = "Thinned"
	res2$parameter = rownames(res2)
	rbind(res1,res2)
}

stopCluster(cl)

res <- do.call('rbind', results)
true.vals <- data.frame(parameter = c("D", "sigma", "lambda0"), 
		true.val =c(300,2.5, 0.0125*30))
library(ggplot2)
ggplot(data = res, aes(x = Scenario, y = estimate)) + 
	geom_boxplot() +
	geom_hline(data = true.vals, aes(yintercept = true.val), col = 'red') +
	facet_wrap(~ parameter, scales = "free_y")+
	theme_bw()
hist(res[res$parameter == "D","estimate"])
abline(v = 300, col ='red')