library(coda)
library(ggplot2)
library(tidyverse)
setwd("C:/Users/Paul/Documents/GitHub/LIDSCR_NIMBLE/output/FisherSimulations")

files <- dir()
files <- grep("Scenario_1_LID|Scenario_2_LID", files, value = TRUE)
# files <- files[!grepl("Scenario2", files)]
# files <- grep("Scenario3", files, value = TRUE)

res <- NULL
modes <- data.frame()
pars <- c('D', 'N', 'sigma', 'lambda')
for( i in files) 
{
	load(i)
	Method <- "LID"
	if(grepl("LID_Sex", i))	Method <- "LID+Sex"
	if(grepl("LID_Sex_Collar", i)) Method <- "LID+Sex/Collar"
	Scenario <- 1
	if(grepl("Scenario_2", i)) Scenario <- 2
	Dmode <- MCMCglmm::posterior.mode(out)[pars]
	tmp <- cbind(data.frame(summary(out[, pars])[[1]]), 
		data.frame(summary(out[, pars])[[2]]),
		mode = as.numeric(Dmode))
	tmp$sd.chain <- apply(out[, pars], 2, sd)
	tmp$parameter <- rownames(tmp)
	tmp$Method <- Method
	tmp$Scenario <- Scenario
	tmp$iter <- i
	res <- rbind(res, tmp)
}

failed.conv <- which(res$sd.chain == 0)
res <- res[!res$iter %in% res$iter[failed.conv], ]

area <- 20.46328

scenario1 <- data.frame(parameter = c('D', 'N', 'sigma', 'lambda'),
	val = c(50/area, 50, 1.5, 0.15), Scenario = 1)
scenario2 <- data.frame(parameter = c('D', 'N', 'sigma', 'lambda'),
	val = c(40/area, 40, 2, 0.10), Scenario = 2)	
scenarios <- rbind(scenario1, scenario2)

facet_names <- c('N', 'sigma~(km)', 'lambda~(day^-1)')
names(facet_names) <- c("N", "sigma", "lambda")
# , labeller = as_labeller(x = facet_names, label_parsed)

# Just keep the first 100 of each simulation:
files <- res[!duplicated(res$iter),]
iter.keep <- files %>% group_by(Method, Scenario) %>% slice(1:100) %>% ungroup() %>% select(iter)
res.true <- res
res <- res[res$iter %in% c(iter.keep[[1]]),]
res$parameter <- factor(res$parameter, levels = c("D", "N", "sigma", "lambda"))

res <- res %>% left_join(scenarios) %>% mutate(MeanBias = (Mean - val)/val*100, 
		MedianBias = (X50. - val)/val*100, ModeBias = (mode - val)/val*100)

res %>% filter(MeanBias > 200)
res %>% filter(iter == 'Scenario_1_LID_Sex_iter_136.Rda')

ggplot(data=res %>% filter(parameter != "D"), aes(y = MeanBias, fill = Method, x = factor(Scenario))) + 
	facet_wrap(~parameter, labeller = as_labeller(x = facet_names, label_parsed)) + 
	geom_boxplot(width = 0.75) +
	geom_hline(yintercept = 0, col = 'red') +
	theme_classic() + ylab("Posterior Mean Bias (%)") + xlab("Scenario") + 
	# scale_fill_grey(start = 0.9, end = 0.25)
	scale_fill_brewer(palette = 1, direction = 1)
ggsave("../../output/FisherResults/FisherSimulationsMean.png", dpi = 'print', 
	width = 9.75, height = 7.35, units = 'in')

ggplot(data=res %>% filter(parameter != "D"), aes(y = MedianBias, fill = Method, x = factor(Scenario))) + 
	facet_wrap(~parameter, labeller = as_labeller(x = facet_names, label_parsed)) + 
	geom_boxplot(width = 0.75) +
	geom_hline(yintercept = 0, col = 'red') +
	theme_classic() + ylab("Posterior Median Bias (%)") + xlab("Scenario") + 
	# scale_fill_grey(start = 0.9, end = 0.25)
	scale_fill_brewer(palette = 1, direction = 1)
ggsave("../../output/FisherResults/FisherSimulationsMedian.png", dpi = 'print', 
	width = 9.75, height = 7.35, units = 'in')

ggplot(data=res %>% filter(parameter != "D"), aes(y = ModeBias, fill = Method, x = factor(Scenario))) + 
	facet_wrap(~parameter, labeller = as_labeller(x = facet_names, label_parsed)) + 
	geom_boxplot(width = 0.75) +
	geom_hline(yintercept = 0, col = 'black') +
	# stat_summary(fun.y="mean", color="red", shape=4)	
	theme_classic() + ylab("Posterior Mode Bias (%)") + xlab("Scenario") + 
	# scale_fill_grey(start = 0.9, end = 0.25)
	scale_fill_brewer(palette = 1, direction = 1)
ggsave("../../output/FisherResults/FisherSimulationsMode.png", dpi = 'print', 
	width = 9.75, height = 7.35, units = 'in')

res.long <- res %>% dplyr::select(Method, Scenario, iter, parameter, MeanBias, ModeBias, MedianBias) %>%
	gather("variable", "value", MeanBias, ModeBias, MedianBias)

facet_names <- c('N', 'sigma~(km)', 'lambda~(day^-1)', 'Mean', 'Median', 'Mode')
names(facet_names) <- c("N", "sigma", "lambda", 'MeanBias', 'MedianBias',  'ModeBias')

ggplot(data=res.long %>% filter(parameter != "D"), aes(y = value, fill = Method, x = factor(Scenario))) + 
	facet_grid(variable~parameter, labeller = as_labeller(x = facet_names, label_parsed)) + 
	geom_boxplot(width = 0.75) +
	geom_hline(yintercept = 0, col = 'black') +
	# stat_summary(fun.y="mean", color="red", shape=4)	
	theme_classic() + ylab("Relative Posterior Bias (%)") + xlab("Scenario") + 
	ylim(-75,150)  +
	# scale_fill_grey(start = 0.9, end = 0.25)
	scale_fill_brewer(palette = 1, direction = 1)
ggsave("../../output/FisherResults/FisherSimulationsAllPosterior.png", dpi = 'print', 
	width = 9.75, height = 7.35, units = 'in')
