#################################
# SCR as a Marked Poisson Process
# Completely latent ID model
#################################

library(sp)
library(coda)
library(raster)

source("FunctionsMCMC.R")

## Load the from Augustine:
load("../data/bobcats.Rdata")

y <- apply(bobcats$capt, 1, FUN = function(x){which(x == 1)})
Time <- 187

traps <- bobcats$traps
xlim <- range(bobcats$mask[,1])
ylim <- range(bobcats$mask[,2])
A <- diff(xlim)*diff(ylim)

mustlink <- matrix(0, nrow = length(y), ncol = length(y))+diag(1, length(y))
cannotlink <- matrix(0, nrow = length(y), ncol = length(y))

for(i in 1:length(y))
{
	marki <- bobcats$mark[i]
	LR <- gsub("[0-9]", "", marki)
	num <-  gsub("[[:alpha:]]", "", marki)
	ml <- bobcats$mark == marki
	mustlink[ml,i] <- 1
	if(LR == ""){
		cannotlink[!ml,i] <- 1
	}else{
		cannotlink[!ml & grepl(LR, bobcats$mark)] <- 1
	}
}

out <- sample_mcmc(y, traps, ylim, xlim, M = 150, Time, tune = c(0.1, 0.01, 2), 
	iterations = 30000, nchains = 3, CR = FALSE, nburnin = 10000,
	priorSigma = NULL, priorLambda = NULL, priorPsi = c(1,1), 
	link = list(mustlink=mustlink, cannotlink=cannotlink), mask = mask)
plot(out[, c("N","NActive")])
plot(out[, c("sigma", "lambda")])
save(out, file = "../ExampleOutput/BobcatMCMCOutput.Rda")
out2 <- sample_mcmc(y, traps, ylim, xlim, M = 150, Time, tune = c(0.1, 0.01, 2), 
	iterations = 20000, nchains = 1, CR = FALSE, nburnin = 5000,
	priorSigma = NULL, priorLambda = NULL, priorPsi = c(1,1), link = list(mustlink=mustlink, cannotlink=cannotlink))
plot(out2[, c("sigma", "N","NActive")])
plot(out2[, c("lambda", "psi")])
summary(out2)

library(data.table)
library(ggplot2)
dat <- data.table(out[[1]])
dat <- dat[, .(iter = .N, maxll = max(ll, na.rm = TRUE), chain = "A"), by = "N"]
setorder(dat, N)
ggplot(data = dat, aes(x=N, y=maxll)) + geom_line() + geom_point(aes(colour = iter))
ggplot(data = dat, aes(x=N, y=ll, group = N)) + geom_boxplot()
 + geom_point(aes(colour = iter))

plot(as.numeric(out[[1]][,"sigma"]), as.numeric(out[[1]][, "N"]))



# Simulate SPIM data:
# N=50
# p01=0.52
# p02=0.26
# lam01=-log(1-p01)
# lam02=-log(1-p02)
# 1-exp(-c(lam01,lam02))
# sigma=0.50
# K=6
# buff=2
# X<- expand.grid(3:8,3:8)
# X=cbind(X,1)
# data <- sim2side(N = 50, lam01 = 0.1, lam02 = 0.1, sigma = 0.5, K = 10,
  # X = X, buff = 3, obstype = "poisson")
# sim2side.plot(data)
