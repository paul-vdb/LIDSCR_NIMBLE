library(coda)
library(ggplot2)
setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/ASCRSimulations")

files <- dir()
# files <- grep("Scenario2", files, value = TRUE)
# files <- files[!grepl("Scenario2", files)]
# files <- grep("Scenario3", files, value = TRUE)

res <- NULL
modes <- data.frame()
pars <- c('sigmatoa', 'sigma', 'lambda', 'g0', 'N[1]', 'N[2]', 'D')
for( i in files) 
{
	load(i)
	if(grepl("Known", i)){
		out <- out2
		Method = "Known ID"
	}else{
		out <- out1
		Method = "Latent ID"
	}
		Dmode <- MCMCglmm::posterior.mode(out)[pars]


	tmp <- cbind(data.frame(summary(out[, pars])[[1]]), 
		data.frame(summary(out[, pars])[[2]]),
		mode = as.numeric(Dmode))
	tmp$parameter <- rownames(tmp)
	tmp$Method = Method
	tmp$Scenario <- 1
	if(grepl("Scenario2", i)) tmp$Scenario <- 2
	if(grepl("Scenario3", i)) tmp$Scenario <- 3
	res <- rbind(res, tmp)
}
# plot(out[,c("N", "sigma", "lambda", "g0")])

# scenario.vals <- data.frame(parameter = c('sigmatoa', 'sigma', 'lambda', 'g0', 'N[1]', 'N[2]', 'D'),
	# val = c(0.00055, 2.3, 0.28, 5.75, 55, 55, 408.1633))
scenario1 <- data.frame(parameter = c('sigmatoa', 'sigma', 'lambda', 'g0', 'N[1]', 'N[2]', 'D'),
	val = c(0.00055, 2.3, 0.28, 5.75, 55, 55, 408.1633), Scenario = 1)
scenario2 <- data.frame(parameter = c('sigmatoa', 'sigma', 'lambda', 'g0', 'N[1]', 'N[2]', 'D'),
	val = c(0.001, 2.3, 0.28, 5.75, 55, 55, 408.1633), Scenario = 2)
scenario3 <- data.frame(parameter = c('sigmatoa', 'sigma', 'lambda', 'g0', 'N[1]', 'N[2]', 'D'),
	val = c(0.05, 2.3, 0.28, 5.75, 55, 55, 408.1633), Scenario = 3)
scenarios <- rbind(scenario1, scenario2, scenario3)

ggplot(data=res, aes(y = Mean, x = Method)) + facet_wrap(~parameter, scales = "free") + 
	geom_violin() + 
	geom_boxplot(width = 0.25) +
	geom_hline(data = scenario.vals, aes(yintercept=val), col = 'red') + 
	theme_classic() + ggtitle("Posterior Means")

ggplot(data=res, aes(y = X50., x = Method)) + facet_wrap(~parameter, scales = "free") + 
	geom_violin() + 
	geom_boxplot(width = 0.25) +
	geom_hline(data = scenario.vals, aes(yintercept=val), col = 'red') + 
	theme_classic() + ggtitle("Posterior Medians")

ggplot(data=res, aes(y = mode, x = Method)) + facet_wrap(~parameter+Scenario, scales = "free") + 
	geom_violin() + 
	geom_boxplot(width = 0.25) +
	geom_hline(data = scenarios, aes(yintercept=val), col = 'red') + 
	theme_classic() + ggtitle("Posterior Modes")

x_names <- c("sigma_t = 0.00055","sigma_t = 0.001","sigma_t = 0.05")
res$tmp <- res$Scenario
res$Scenario <- factor(res$Scenario, label = x_names)
ggplot(data=res[res$parameter == "D",], aes(y = mode, x = Method)) + facet_wrap(~Scenario) + 
	geom_violin() + 
	geom_boxplot(width = 0.15) +
	geom_hline(aes(yintercept=408.16330), col = 'red') + 
	theme_classic() + ggtitle("Posterior Modes") 

ggplot(data=res[res$parameter == "D",], aes(y = X50., x = Method)) + facet_wrap(~Scenario) + 
	geom_violin() + 
	geom_boxplot(width = 0.15) +
	geom_hline(aes(yintercept=408.16330), col = 'red') + 
	theme_classic() + ggtitle("Posterior Medians") 

ggplot(data=res[res$parameter == "D",], aes(y = Mean, x = Method)) + facet_wrap(~Scenario) + 
	geom_violin() + 
	geom_boxplot(width = 0.15) +
	geom_hline(aes(yintercept=408.16330), col = 'red') + 
	theme_classic() + ggtitle("Posterior Means") 


facet_names <- c("Scenario 1", "Scenario 2")
names(facet_names) <- c("1", "2")
ggplot(data=res[res$parameter == "D",], aes(y = mode, x = Method)) + 
	facet_wrap(~Scenario, scales = "fixed", labeller = as_labeller(x = facet_names)) + 
	geom_violin() + 
	geom_boxplot(width = 0.25) +
	geom_hline(aes(yintercept=408.1633), col = 'red') + 
	theme_classic() + ggtitle("Posterior Modes")


modes <- data.frame(parameter = names(Dmode1), Mean = Dmode1, Method = "Known ID")
modes <- rbind(modes, data.frame(parameter = names(Dmode2), Mean = Dmode2, Method = "Latent ID"))
pars <- c('sigmatoa', 'sigma', 'lambda', 'g0', 'D')
ggplot(data=modes[modes$parameter %in% pars,], aes(y = Mean, x = Method)) + facet_grid(Scenario~parameter, scales = "free") + 
	geom_violin() + 
	geom_boxplot(width = 0.25) +
	geom_hline(data = scenarios, aes(yintercept=val), col = 'red') + 
	theme_classic() + ggtitle("Posterior Modes")

ggplot(data=modes[modes$parameter %in% "D",], aes(y = Mean, x = Method, col = Scenarios)) +
	geom_violin() + 
	geom_boxplot(width = 0.25) +
	geom_hline(data = scenarios, aes(yintercept=val), col = 'red') + 
	theme_classic() + ggtitle("Posterior Modes")


par(mfrow = c(2,2))
boxplot(Dmode1[names(Dmode1) == "D"], ylab = "Posterior Mode D", main = "Known ID")
abline(h = 408.1633, col = 'red')
boxplot(res[res$parameter == "D" & res$Method == "Known ID", "Mean"], ylab = "Posterior Mean D", main = "Known ID")
abline(h = 408.1633, col = 'red')
boxplot(Dmode2[names(Dmode2) == "D"], ylab = "Posterior Mode D", main = "Latent ID")
abline(h = 408.1633, col = 'red')
boxplot(res[res$parameter == "D" & res$Method == "Latent ID", "Mean"], ylab = "Posterior Mean D", main = "Latent ID")
abline(h = 408.1633, col = 'red')

par(mfrow = c(2,1))
boxplot(Dmode1[names(Dmode1) == "N[1]"])
abline(h = 55, col = 'red')
boxplot(Dmode2[names(Dmode2) == "N[1]"])
abline(h = 55, col = 'red')


library(data.table)
res.dt <- data.table(res)
res.dt[, .(avg = mean(X50.)), by = c("Method", "parameter")]

N <- res[res$parameter == "N",]
hist(N[, "Mean"])
abline(v = 60, col = 'red')
mean(N[, "Mean"])
median(N[, "Mean"])

boxplot((N[, "Mean"] - 60)/60)
abline(h = 0, col = 'red')

g0 <- res[res$parameter == "g0",]
hist(g0[, "Mean"])
abline(v = 5.75, col = 'red')
mean(g0[, "Mean"])
median(g0[, "Mean"])
boxplot((g0[, "Mean"] - 5.75)/5.75)
abline(h = 0, col = 'red')

sigma <- res[rownames(res) == "sigma",]
hist(sigma[, "Mean"])
abline(v = 2.3, col = 'red')
mean(sigma[, "Mean"])
median(sigma[, "Mean"])
boxplot((sigma[, "Mean"] - 2.3)/2.3)
abline(h = 0, col = 'red')

lambda <- res[rownames(res) == "lambda",]
hist(lambda[, "Mean"])
abline(v = 0.27, col = 'red')
mean(lambda[, "Mean"])
median(lambda[, "Mean"])
boxplot((lambda[, "Mean"] - 0.27)/0.27)
abline(h = 0, col = 'red')

sigmatoa <- res[rownames(res) == "sigmatoa",]
hist(sigmatoa[, "Mean"])
abline(v = 0.00055, col = 'red')
mean(sigmatoa[, "Mean"])
median(sigmatoa[, "Mean"])

true.vals <- data.frame(parameter = c("N", "g0", "sigma", "lambda", "sigmatoa"), values = c(60, 5.75, 2.3, 0.27, 0.00055))
res <- merge(res, true.vals, by = "parameter")
res$bias <- (res$Mean - res$values)/res$values*100
res$biasmed <- (res$X50. - res$values)/res$values*100
res$parameter <- as.factor(res$parameter)
ggplot(data = res, aes(x = parameter, y = bias)) + geom_boxplot() +
	geom_hline(yintercept = 0, col = 'red') + theme_classic() + ylab("% Bias")
ggplot(data = res, aes(x = parameter, y = biasmed)) + geom_boxplot() +
	geom_hline(yintercept = 0, col = 'red') + theme_classic() + ylab("% Bias")



# Scenario 2:
files <- dir()
files <- files[grepl("Scenario_2", files)]
res <- NULL
for( i in files) 
{
	load(i)
	tmp <- cbind(data.frame(summary(out[, c('sigmatoa', 'sigma', 'lambda', 'g0', 'D', 'EN', 'psi')])[[1]]), 
		data.frame(summary(out[, c('sigmatoa', 'sigma', 'lambda', 'g0', 'D', 'EN', 'psi')])[[2]]))
	tmp$parameter <- rownames(tmp)
	res <- rbind(res, tmp)
}
# plot(out)
hist(res[res$parameter == "EN", "Mean"])
abline(v = 60, col = "red")
hist(res[res$parameter == "EN", "X50."])
abline(v = 60, col = "red")
boxplot(res[res$parameter == "EN", "X50."])
abline(h = 60, col = "red")


hist(res[res$parameter == "psi", "Mean"])
abline(v = 0.30, col = "red")
hist(res[res$parameter == "sigma", "Mean"])
abline(v = 2.3, col = 'red')
hist(res[res$parameter == "g0", "Mean"])
abline(v = 5.75, col = 'red')
hist(res[res$parameter == "sigmatoa", "Mean"])
abline(v = 0.00055, col = 'red')
hist(res[res$parameter == "lambda", "Mean"])
abline(v = 0.27, col = 'red')










d2 <- t(apply(mask, 1, FUN = function(x){(traps[,1] - x[1])^2 + (traps[,2] - x[2])^2}))
pkj <- 1-exp(-5.5*exp(-d2/(2*2.2^2)))

post.x <- list()
x.mle <- NULL
for(i in 1:nrow(capt)){
	post.x[[i]] <- apply(d2, 1, FUN = function(x){dnorm_vector_marg(toa[i,], sqrt(x)/330, 0.00045, capt[i,], 0)})
	post.x[[i]] <- post.x[[i]]*apply(d2, 1, FUN = function(x){dbinom_vector(capt[i,], rep(1,6),  1-exp(-5.5*exp(-x/(2*2.2^2))), 0)})
	x.mle <- c(x.mle, which.max(post.x[[i]]))
}

ll <- function(pars)
{
	ll <- 0
	for( i in 1:nrow(capt)) ll <- ll + dnorm_vector_marg(toa[i,], sqrt(d2[x.mle[i],])/330, pars, capt[i,], 1)
	-ll
}
optim(0.001, ll)

library(ggplot2)
ggplot(data = data.frame(mask), aes(x = x, y=y)) + geom_tile(aes(fill = post.x[[indx[1]]])) + 
	ggtitle("indx1") + geom_point(data = data.frame(traps), aes(x=x,y=y), col = 'red') 
dev.new() 
ggplot(data = data.frame(mask), aes(x = x, y=y)) + geom_tile(aes(fill = post.x[[indx[2]]])) + 
	ggtitle("indx2") + geom_point(data = data.frame(traps), aes(x=x,y=y), col = 'red') 
dev.new() 	
ggplot(data = data.frame(mask), aes(x = x, y=y)) + geom_tile(aes(fill = post.x[[indx[3]]]))+ 
	ggtitle("indx3") + geom_point(data = data.frame(traps), aes(x=x,y=y), col = 'red') 
dev.new() 	
ggplot(data = data.frame(mask), aes(x = x, y=y)) + geom_tile(aes(fill = post.x[[indx[9]]])) + 
	ggtitle("indx4") + geom_point(data = data.frame(traps), aes(x=x,y=y), col = 'red') 

# All obs:
all.x <- 1+numeric(nrow(mask))
for(i in indx) {all.x <- all.x * post.x[[i]]}
tmp <- all.x[which(all.x > 0)]
max(tmp)
quantile(tmp, c(0.0275, 0.5, 0.975))
ggplot(data = data.frame(mask), aes(x = x, y=y)) + geom_tile(aes(fill = log(all.x))) + 
	ggtitle("indx4") + geom_point(data = data.frame(traps), aes(x=x,y=y), col = 'red') 
