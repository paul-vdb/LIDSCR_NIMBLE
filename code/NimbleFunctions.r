dX <- nimbleFunction(
    run = function(x = double(1), log = integer(0, default = 0)) {
        returnType(double(0))
        return(0)
    }
)

rX <- nimbleFunction(
    run = function(n = integer(0)) {
        print('Error: rX() function not implemented')
        returnType(double(1))
        return(c(0, 0))
    }
)

# Stupid holder function for sampling ID.
dID <- nimbleFunction(
    run = function(x = integer(0), log = integer(0, default = 0)) {
        returnType(double(0))
		if(log) return(0) else return(1)
    }
)

rID <- nimbleFunction(
    run = function(n = integer(0)) {
        print('Error: rX() function not implemented')
        returnType(integer(0))
        return(1)
    }
)

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

dbinom_vector_ascr <- nimbleFunction(
  run = function( x = double(1),
                  size = double(1),
                  prob = double(1), 
				  pcapt = double(0),
                  log = integer(0, default = 0)
                  ) {
    returnType(double(0))
    logProb <- sum(dbinom(x, prob = prob, size = size, log = TRUE)) - log(pcapt)
    if(log) return(logProb) else return(exp(logProb))
  })

#' @rdname dbinom_vector
#' @export
rbinom_vector_ascr <- nimbleFunction(
  run = function( n = integer(0, default = 1),
                  size = double(1),
                  prob = double(1),
				  pcapt = double(0)
  ) {
    returnType(double(1))
    return(rbinom(length(size), prob = prob, size = size))
  })

registerDistributions(
    list(dX = list(BUGSdist = 'dX()',
                   types = c('value = double(1)')),
		 dID = list(BUGSdist = 'dID()',
                   types = c('value = integer(0)')),
		 dbinom_vector_ascr = list(BUGSdist = 'dbinom_vector_ascr()',
				   types =c('value = double(1)')),
		 dnorm_vector_marg = list(BUGSdist = 'dnorm_vector_marg(mean, sd, y)',
				   types = c('value = double(1)', 'mean = double(1)', 'sd = double(0)', 'y = double(1)'))		   
				   ) )

sampler_myX <- nimbleFunction(
    name = 'sampler_myX',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        xlim <- control$xlim
        ylim <- control$ylim
        calcNodesAll <- model$getDependencies(target)
        calcNodesNoIDs <- model$getDependencies(target)
        if(model$getDistribution(target) != 'dX') stop('myX sampler error')
        nodeIndex <- as.numeric(gsub('^X\\[([[:digit:]]+), .*$', '\\1', target))
        zNode <- paste0('z[', nodeIndex, ']')
        dNode <- paste0('d2[', nodeIndex, ', 1:', control$J, ']')
        copyNodes <- c(target, dNode)
    },
    run = function() {
        newx <- runif(1, xlim[1], xlim[2])
        newy <- runif(1, ylim[1], ylim[2])
        if(model[[zNode]] == 0) {
            logMHR <- 0
            jump <- decide(logMHR)
            if(jump) {
                model[[target]] <<- c(newx, newy)
                model$calculate(copyNodes)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = copyNodes, logProb = FALSE)
            }
            return()
        }
        anyID <- any(model[['ID']] == nodeIndex)
        if(anyID) { lpcurrent <- model$getLogProb(calcNodesAll)
        } else { lpcurrent <- model$getLogProb(calcNodesNoIDs) }
        model[[target]] <<- c(newx, newy)
        if(anyID) { lpprop <- model$calculate(calcNodesAll)
        } else { lpprop <- model$calculate(calcNodesNoIDs) }
        logMHR <- lpprop - lpcurrent
        jump <- decide(logMHR)
        if(jump) {
            if(anyID) {
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesAll, logProb = TRUE)
            } else {
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoIDs, logProb = TRUE)
            }
        } else {
            if(anyID) {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesAll, logProb = TRUE)
            } else {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoIDs, logProb = TRUE)
            }
        }
    },
    methods = list( reset = function() {} )
)

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
        k <- control$M
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
					logProbs[i] <<- model$calculate(calcNodes)
                    if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
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


sampler_myIDZ <- nimbleFunction(
    name = 'sampler_myIDZ',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		calcNodesZ <- model$getDependencies('z')
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        k <- control$M
        probs <- numeric(k)
        logProbs <- numeric(k)
    },
    run = function() {
		psi <- model[['psi']]
        currentValue <- model[[target]]
		logProbs[currentValue] <<- model$getLogProb(calcNodes)
		n_currentValue <- sum(model[['ID']] == currentValue)
		if(n_currentValue == 1)
		{
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi) - log(1-psi) - model[['Hk']][currentValue]
			model[['z']][currentValue] <<- 0
		}
        for(i in 1:k) {
			if(i != currentValue){
				if(model[['z']][i] == 1) {
					nk <- sum(model[['ID']] == i)
					if(nk == 0)
					{
						logProbs[i] <<- -Inf
					}else { 
						model[[target]] <<- i
						logProbs[i] <<- model$calculate(calcNodes)
						if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
					}
				}else {
					model[[target]] <<- i
					logProbs[i] <<- model$calculate(calcNodes) + log(psi) - log(1-psi) - model[['Hk']][i]
					if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
				}
            }
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
            model[[target]] <<- newValue
			if(model[['z']][newValue] == 0){
				model[['z']][newValue] <<- 1
				model$calculate(calcNodesZ)
			}
            model$calculate(calcNodes)	# I've made ID independent of z so this shouldn't double.
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
			model[['z']][newValue] <<- 1
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

sampler_mySPIM <- nimbleFunction(
    name = 'sampler_mySPIM',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		calcNodesZ <- model$getDependencies('z')
        targetNodesAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)		
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        k <- control$M
        probs <- numeric(k)
        logProbs <- numeric(k)
		cannotlink <- extractControlElement(control, 'cannotlink', 'identity')	# n x n matrix where a 1 indicates i cannot link to detection j.
        nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', targetNodesAsScalar[1]))
		n_grp <- length(targetNodesAsScalar)
    },
    run = function() {
		psi <- model[['psi']]
        currentValue <- model[["ID"]][nodeIndex]
		logProbs[currentValue] <<- model$getLogProb(calcNodes)
		n_currentValue <- sum(model[['ID']] == currentValue)
		no_link <- 0
		if(n_currentValue == n_grp)
		{
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi) - log(1-psi) - model[['Hk']][currentValue]
			model[['z']][currentValue] <<- 0
		}
        for(i in 1:k) {
			if(i != currentValue){
				if(model[['z']][i] == 1) {
					nk <- sum(model[['ID']] == i)
					if(nk == 0)
					{
						logProbs[i] <<- -Inf
					}else {
						no_link <- sum(cannotlink[model[['ID']] == i, nodeIndex])	# This is the check to see if there are cannot links.
						if(no_link > 0)
						{
							logProbs[i] <<- -Inf
						}else {
							# values(model, targetNodesAsScalar) <<- i
							values(model, targetNodesAsScalar) <<- rep(i, n_grp) 					
							logProbs[i] <<- model$calculate(calcNodes)
							if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
						}
					}
				}else {
					# values(model, targetNodesAsScalar) <<- i
					values(model, targetNodesAsScalar) <<- rep(i, n_grp) ##  replace with this					
					logProbs[i] <<- model$calculate(calcNodes) + log(psi) - log(1-psi) - model[['Hk']][i]
					if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
				}
            }
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
            # values(model, targetNodesAsScalar) <<- newValue
			values(model, targetNodesAsScalar) <<- rep(newValue, n_grp) ##  replace with this			
			if(model[['z']][newValue] == 0){
				model[['z']][newValue] <<- 1
				model$calculate(calcNodesZ)
			}
            model$calculate(calcNodes)	# I've made ID independent of z so this shouldn't double effort.
            nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
			model[['z']][newValue] <<- 1
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

