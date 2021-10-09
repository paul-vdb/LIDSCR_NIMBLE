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
    run = function(x = integer(0), lam = double(0), log = integer(0, default = 0)) {
        returnType(double(0))
		if(log) return(log(lam)) else return(lam)
    }
)

rID <- nimbleFunction(
    run = function(n = integer(0), lam = double(0)) {
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
	m <- sum(y)
	if(m == 1){
		logProb <- 0
	}else{
		td <- (x - mean)*y
		etd <- sum(td)/m
		logProb <- (1-m)*log(sd) - sum(y*(td - etd)^2)/(2*sd^2)
	}	
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


dMarg <- nimbleFunction(
  run = function( x = double(0),
				  toa = double(1),
				  y = double(1),
				  J = double(0),
                  size = double(1),
                  pkj = double(2),
				  sd = double(0),
				  expTime = double(2),
				  M = double(0),
                  log = integer(0, default = 0)
                  ) {
    returnType(double(0))
		likelihood <- 0
		for(k in 1:M) {
			lpcapt <- dbinom_vector(x = y, size = size, prob = pkj[k,1:J], log = 1)
			lptoa <- dnorm_vector_marg(x = toa, mean = expTime[k,1:J], sd = sd, y = y, log = 1)
			likelihood <- likelihood + exp(lptoa + lpcapt)
		}
    if(log) return(log(likelihood)) else return(likelihood)
  })


registerDistributions(
    list(dX = list(BUGSdist = 'dX()',
                   types = c('value = double(1)')),
		 dID = list(BUGSdist = 'dID(lam)',
                   types = c('value = integer(0)', 'lam = double(0)')),
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

sampler_myX2 <- nimbleFunction(
    name = 'sampler_myX2',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        xlim <- control$xlim
        ylim <- control$ylim
		nv <- extractControlElement(control, 'Noccasion', 1)		
		v <- extractControlElement(control, 'occasion', 1)		
        calcNodesAll <- model$getDependencies(target)
        nodeIndex <- as.numeric(gsub("\\D*(\\d+).*", '\\1', target))
        if(nv == 1){
			zNode <- paste0('z[', nodeIndex, ']')
		}else{
			zNode <- paste0('z[', nodeIndex, ',', v,']')
		}		
		scale <- extractControlElement(control, 'scale', 1)		
    },
    run = function() {
		if((model[[zNode]] == 0)) {
				model[[target]] <<- c(runif(1, xlim[1],xlim[2]), runif(1, ylim[1], ylim[2]))
				model$calculate(calcNodesAll)
				jump <- TRUE
        }else {
			lpcurrent <- model$getLogProb(calcNodesAll)	
			XCurrent <- model[[target]]	
			model[[target]] <<- XCurrent + rnorm(2, 0, sd = scale)
			prior <- model$calculateDiff(target)
			if(prior == -Inf)
			{
				jump <- FALSE
			}else{
				
				lpprop <- model$calculate(calcNodesAll)
				logMHR <- lpprop - lpcurrent + prior
				jump <- decide(logMHR)
			}
		}
		if(jump) {
				nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesAll, logProb = TRUE)
		} else {
				nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesAll, logProb = TRUE)
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
        nodeIndex <- as.numeric(sub("\\D*(\\d+).*", "\\1", target))
		nv <- extractControlElement(control, 'Noccasion', 1)		
		occ <- extractControlElement(control, 'IDoccasion', 1)
		n_obs <- length(model[['ID']])
		if(nv == 1){ 
			IDMatch <- 1:n_obs
		}else{
			v <- as.numeric(sub(".*,\\D*(\\d+).*", "\\1", target))
			IDMatch <- which(occ == v)
		}
    },
    run = function() {
        if((model[[target]] == 1) & (any(model[['ID']][IDMatch] == nodeIndex))) return()
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
        M <- control$M
        probs <- numeric(M)
        logProbs <- numeric(M)
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
			model[['z']][currentValue] <<- 0
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi)-model[['Hk']][currentValue]
		}
        for(k in 1:M) {
			if(k != currentValue){
				if(model[['z']][k] == 1) {
					nk <- sum(model[['ID']] == k)
					if(nk == 0)
					{
						logProbs[k] <<- -Inf
					}else {
						no_link <- sum(cannotlink[model[['ID']] == k, nodeIndex])	# This is the check to see if there are cannot links.
						if(no_link > 0)
						{
							logProbs[k] <<- -Inf
						}else {
							values(model, targetNodesAsScalar) <<- rep(k, n_grp) 					
							logProbs[k] <<- model$calculate(calcNodes) + log(1-psi)
							if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
						}
					}
				}else {
					values(model, targetNodesAsScalar) <<- rep(k, n_grp)
					logProbs[k] <<- model$calculate(calcNodes) + log(psi)-model[['Hk']][k]
					if(is.nan(logProbs[k])) logProbs[k] <<- -Inf
				}
            }
        }
		# Note that logProbs of z=1 and nk=0 is -Inf, or it had better be!y
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
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
			model[['z']][currentValue] <<- 1
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

sampler_mySigmaToa <- nimbleFunction(
    name = 'sampler_mySigmaToa',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodesAll <- model$getDependencies(target)
		mi <- control$mi
		n <- length(mi)
		J <- control$J
    },
    run = function() {
		alpha <- 1
		ssq <- 0
		mdiff <- 0
		tdiff <- numeric(J)
		for(i in 1:n){
			if(mi[i] > 1){
				tdiff <- (model[['toa']][i,1:J] - model[['expTime']][model[['ID']][i],1:J])*model[['y']][i,1:J]
				mdiff <- sum(tdiff[1:J])/mi[i]
				ssq <- ssq + sum(model[['y']][i,1:J]*(tdiff[1:J] - mdiff)^2)
				alpha = alpha + (mi[i]-1)/2
			}
		}
	model[[target]]	<<- 1/sqrt(rgamma(1, shape = alpha, rate = ssq/2))
	model$calculate(calcNodesAll)
	nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesAll, logProb = TRUE)
    },
    methods = list( reset = function() {} )
)

##################
# Generalized for Multi Session:
##################
sampler_myIDZ <- nimbleFunction(
    name = 'sampler_myIDZ',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        M <- control$M
		nv <- extractControlElement(control, 'Noccasion', 1)		
		v <- extractControlElement(control, 'occasion', 1)
		if(nv == 1){ 
			zNodes <- paste0('z[1:', M, ']')
			Hkv <- 'Hk'
		}else{
			zNodes <- paste0('z[1:', M,',', v,']')
			Hkv <- paste0('Hk[1:', M, ',', v,']')
		}
		calcNodesZ <- model$getDependencies(zNodes)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        probs <- numeric(M)
        logProbs <- numeric(M)
	    nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
    },
    run = function() {
		psi <- model[['psi']]
        currentValue <- model[[target]]
		logProbs[currentValue] <<- model$getLogProb(calcNodes)
		n_currentValue <- sum(model[['ID']] == currentValue)
		if(n_currentValue == 1)
		{
			logProbs[currentValue] <<- logProbs[currentValue] + log(psi) - model[[Hkv]][currentValue]
			model[[zNodes]][currentValue] <<- 0
		}
        for(k in 1:M) {
			if(k != currentValue){
				# Start with check if it can match and assign 0 to 1 prob or 1 already prob.
				if(model[[zNodes]][k] == 0){
					logProbs[k] <<- log(psi)-model[[Hkv]][k]				
				}else{
					if(sum(model[['ID']] == k) == 0){
						logProbs[k] <<- -Inf
					}else{
						logProbs[k] <<- log(1-psi)
					}
				}
				# If it is not a zero-inflated value then find the full conditional.
				if(logProbs[k] != -Inf)
				{
					model[[target]] <<- k
					logProbs[k] <<- logProbs[k] + model$calculate(calcNodes)		
				}
				
				# If it's a bad number make it zero.
				if(is.nan(logProbs[k])){
						logProbs[k] <<- -Inf
				}
			}
        }
        logProbs <<- logProbs - max(logProbs)
        probs <<- exp(logProbs)
        newValue <- rcat(1, probs)
        if(newValue != currentValue) {
            model[[target]] <<- newValue
            model$calculate(calcNodes)			
			if(model[[zNodes]][newValue] == 0){
				model[[zNodes]][newValue] <<- 1
				model$calculate(calcNodesZ)
			}
            nimCopy(from = model, to = mvSaved, row = 1, nodes = c(target, zNodes), logProb = TRUE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        } else {
			model[[zNodes]][currentValue] <<- 1
            nimCopy(from = mvSaved, to = model, row = 1, nodes = c(target, zNodes), logProb = TRUE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        }
    },
    methods = list( reset = function() { } )
)

## Depricated single session version of the IDZ sampler.
# sampler_myIDZ <- nimbleFunction(
    # name = 'sampler_myIDZ',
    # contains = sampler_BASE,
    # setup = function(model, mvSaved, target, control) {
        # calcNodes <- model$getDependencies(target)
        # calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
		# calcNodesZ <- model$getDependencies('z')
        # isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)
        # calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        # calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        # k <- control$M
        # probs <- numeric(k)
        # logProbs <- numeric(k)
	    # nodeIndex <- as.numeric(gsub('[^[:digit:]]', '', target))
		##check <- paste0("y[", nodeIndex,",1:6]")
		##count <- 0
		# zAdjust <- numeric(k)
    # },
    # run = function() {
		##count <<- count+1
		# psi <- model[['psi']]
        # currentValue <- model[[target]]
		# logProbs[currentValue] <<- model$getLogProb(calcNodes)
		# n_currentValue <- sum(model[['ID']] == currentValue)
		# if(n_currentValue == 1)
		# {
			# zAdjust[currentValue] <<- 1
			# model[['z']][currentValue] <<- 0
		# }
        # for(i in 1:k) {
			# if(i != currentValue){
				# if((sum(model[['ID']] == i) == 0 & model[['z']][i] == 1))
				# {
					# logProbs[i] <<- -Inf
					# zAdjust[i] <<- 0
				# }else {
					# model[[target]] <<- i
					# logProbs[i] <<- model$calculate(calcNodes)
					# if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
					
					## Keep track of which ones are existing or not.
					# zAdjust[i] <<- 1-model[['z']][i]
				# }
			# }
        # }
		# logProbs <<- logProbs + (1-zAdjust)*log(1-psi) + zAdjust*(log(psi)-model[['Hk']])
        # logProbs <<- logProbs - max(logProbs)
        # probs <<- exp(logProbs)
        # newValue <- rcat(1, probs)
        # if(newValue != currentValue) {
            # model[[target]] <<- newValue
            # model$calculate(calcNodes)	# I've made ID independent of z so this shouldn't double.			
			##if(count %% 100 == 0 & nodeIndex == 85) {
			##	print(model$getLogProb(check))	# REMOVE LATER
			##}
			# if(model[['z']][newValue] == 0){
				# model[['z']][newValue] <<- 1
				# model$calculate(calcNodesZ)
			# }
            # nimCopy(from = model, to = mvSaved, row = 1, nodes = c(target, 'z'), logProb = TRUE)
            # nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            # nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        # } else {
			# model[['z']][newValue] <<- 1
            # nimCopy(from = mvSaved, to = model, row = 1, nodes = c(target, 'z'), logProb = TRUE)
            # nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
            # nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
        # }
    # },
    # methods = list( reset = function() { } )
# )


