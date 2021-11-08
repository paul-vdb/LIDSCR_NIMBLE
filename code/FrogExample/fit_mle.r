## Main fitting function.
## - Arguments:
##   - capthist:     Capture history matrix.
##   - ids:          A vector of IDs, where each element indicates which individual is associated with each capture history.
##   - traps:        Matrix of detector locations.
##   - mask:         Mask object with points for numerical integration over activity centre locations. Must include an attribute called "area", providing the area covered by a single point.
##   - detfn:        The detection function to use; either "hn" (halfnormal) or "hhn" (hazard halfnormal).
##   - start:        Start values for numerical maximisation. Order of parameters is D, lambda0/g0, sigma, lambda_c, sigma_t. 
##   - toa:          Time of arrival matrix in the same structure as the capthist object (optional).
##   - speed_sound:  The speed of sound in metres per second.
stacked.scr.fit <- function(capthist, ids, traps, mask, detfn, start, toa = NULL, speed_sound = 330, trace = FALSE){
	## ids must be 1:N animals
	ids = as.numeric(factor(ids))
	
    ## Number of traps.
    nTraps <- nrow(traps)
    ## Number of mask points.
    nMask <- nrow(mask)
    ## Area of a single mask point.
    aMask <- attr(mask, "area")
    ## Indicator for whether or not times of arrival are used.
    use_toa <- !is.null(toa)
    ## Calculating distances between mask points and detectors.
    maskDists <- eucdist(mask, traps)
    if (use_toa){
        ## Creating TOA sum of squares matrix.
        toa_ssq  <- make_toa_ssq(toa, eucdist(traps, mask), speed_sound)
    } else {
        ## Dummy object if not used.
        toa_ssq <- matrix(0, nrow = 1, ncol = 1)
    }
    ## Sorting capture histories in numerical order by ID.
    uniqueIDs <- sort(unique(ids))
    nAnimals <- length(uniqueIDs)

    ## Indicator for detection function.
    if (detfn == "hn"){
        hn <- TRUE
    } else if (detfn == "hhn"){
        hn <- FALSE
    } else {
        stop("The argument detfn must either be 'hn' or 'hhn'.")
    }
    ## Converting parameters to link scale.
    start.link <- numeric(length(start))
    start.link[c(1, 3, 4)] <- log(start[c(1, 3, 4)])
    if (hn){
        start.link[2] <- qlogis(start[2])
    } else {
        start.link[2] <- log(start[2])
    }
    if (use_toa){
        start.link[5] <- log(start[5])
    }
							 			  
    ## Fitting model.
    fit <- nlminb(start.link, scr_nll_stacked,
                  caps = capthist,
                  aMask = aMask,
                  maskDists = maskDists,
				  ID = ids,
                  toa = toa,
                  toa_ssq = toa_ssq,
                  use_toa = use_toa,
                  hn = hn,
                  trace = trace)
    ## Approximating Hessian.
    hess <- optimHess(fit$par, scr_nll_stacked,
                      caps = capthist,
                      aMask = aMask,
                      maskDists = maskDists,
                      ID = ids,
                      toa = toa,
                      toa_ssq = toa_ssq,
                      use_toa = use_toa,
                      hn = hn,
                      trace = trace)
	
	## Calculating confidence intervals
    ## - Using the (sqrt of) diagonals of (-ve) Hessian obtained from optim
    ##    - i.e. Information matrix
    ## - Wald CIs calculated by sapply() loop
    ##    - Loops through each of fitted parameters and calculates lower/upper bounds
    ## Note: fitted pars must be on LINK scale
    ##     : if matrix is singular, none of the SEs or CIs are calculated (inherits/try statement)
    fittedPars = fit$par
    if(inherits(try(solve(hess), silent = TRUE), "try-error")) {
        ## Hessian is singular
        warning("Warning: singular hessian")
        ## SE and Wald CIs not calculated
        se = NA
        waldCI = matrix(NA, nrow = length(fittedPars), ncol = 2)
        ## But columns still need to be returned (if/when simulations are run)
        cnames = c("Estimate", "SE", "Lower", "Upper")
    } else {
        ## Calculating basic Wald CIs
        se = sqrt(diag(solve(hess)))
        waldCI = t(sapply(1:length(fittedPars),
                          function(i) fittedPars[i] + (c(-1, 1) * (qnorm(0.975) * se[i]))))
        ## Back-transforming the confidence limits, depending on whether we're using lambda0 or g0
        if(hn) {
            waldCI = rbind(exp(waldCI[1, ]),
                           plogis(waldCI[2, ]),
                           exp(waldCI[3:length(fittedPars), ]))
        } else {
            waldCI = exp(waldCI)
        }
        
        ## Using the delta method to get the standard errors
        ## - G = jacobian matrix of partial derivatives of back-transformed
        ##    - i.e. log(D) -> exp(D) -- deriv. --> exp(D)
        ##    - Note: 1st deriv of plogis (CDF) = dlogis (PDF)
        if (hn){
            G.mult <- c(exp(fittedPars[1]),
                           dlogis(fittedPars[2]),
                           exp(fittedPars[3:length(fittedPars)]))
        } else {
            G.mult <- c(exp(fittedPars[1]),
                           exp(fittedPars[2]),
                           exp(fittedPars[3:length(fittedPars)]))
        }
        G = diag(length(fittedPars)) * G.mult
        se = sqrt(diag(G %*% solve(hess) %*% t(G)))
    }
    
    
    ## Returning the fitted parameters in a named vector
    ## - First checks to see if TOA is being used,
    ##    then inserts par names in front of "sigma_toa"
    parNames = NULL
    if(!(is.null(toa) || is.na(toa))) {
        parNames = "sigma_toa"
    }
    if (hn) {
            parNames = c("D", "g0", "sigma", "lambda", parNames)
            fittedPars = c(exp(fittedPars[1]),
                           plogis(fittedPars[2]),
                           exp(fittedPars[3:length(fittedPars)]))
        }else {
        parNames = c("D", "lambda0", "sigma", "lambda", parNames)
        
        fittedPars = exp(fittedPars)
    }
    results = cbind(fittedPars, se, waldCI)
    dimnames(results) = list(parNames, c("Estimate", "SE", "Lower", "Upper"))
    
	out <- list(results = results, capthist = capthist, 
		mask = mask, aMask = aMask, maskDists = maskDists, 
		speed_sound = speed_sound, ids = ids,
		traps = traps, detfn = detfn, toa = toa, toa_ssq = toa_ssq)
	
	return(out)
}


## Main fitting function.
## - Arguments:
##   - capthist:     Capture history matrix.
##   - ids:          A vector of IDs, where each element indicates which individual is associated with each capture history.
##   - traps:        Matrix of detector locations.
##   - mask:         Mask object with points for numerical integration over activity centre locations. Must include an attribute called "area", providing the area covered by a single point.
##   - detfn:        The detection function to use; either "hn" (halfnormal) or "hhn" (hazard halfnormal).
##   - start:        Start values for numerical maximisation. Order of parameters is D, lambda0/g0, sigma, lambda_c, sigma_t. 
##   - toa:          Time of arrival matrix in the same structure as the capthist object (optional).
##   - speed_sound:  The speed of sound in metres per second.
ascr.fit <- function(capthist, traps, mask, detfn, start, toa = NULL, speed_sound = 330, trace = FALSE){
	
    ## Number of traps.
    nTraps <- nrow(traps)
    ## Number of mask points.
    nMask <- nrow(mask)
    ## Area of a single mask point.
    aMask <- attr(mask, "area")
    ## Indicator for whether or not times of arrival are used.
    use_toa <- !is.null(toa)
    ## Calculating distances between mask points and detectors.
    maskDists <- eucdist(mask, traps)
    if (use_toa){
        ## Creating TOA sum of squares matrix.
        toa_ssq  <- make_toa_ssq(toa, eucdist(traps, mask), speed_sound)
    } else {
        ## Dummy object if not used.
        toa_ssq <- matrix(0, nrow = 1, ncol = 1)
    }

    ## Indicator for detection function.
    if (detfn == "hn"){
        hn <- TRUE
    } else if (detfn == "hhn"){
        hn <- FALSE
    } else {
        stop("The argument detfn must either be 'hn' or 'hhn'.")
    }
    ## Converting parameters to link scale.
    start.link <- numeric(length(start))
    start.link <- log(start)
    if (hn){
        start.link[1] <- qlogis(start[1])
    } else {
        start.link[1] <- log(start[1])
    }
    if (use_toa){
        start.link[3] <- log(start[3])
    }
	
	neg_call_ll <- function(pars)
	{
		call.mat <- ascr_calls(pars, caps = capthist, 
			aMask = aMask,
			maskDists = maskDists,
			toa = toa,
			toa_ssq = toa_ssq,
			use_toa = use_toa,
			hn = hn)
		-sum(log(colSums(exp(call.mat))*aMask))		
	}

    ## Fitting model.
    fit <- nlminb(start.link, neg_call_ll)
    ## Approximating Hessian.
    hess <- optimHess(fit$par, neg_call_ll)
	
	## Calculating confidence intervals
    ## - Using the (sqrt of) diagonals of (-ve) Hessian obtained from optim
    ##    - i.e. Information matrix
    ## - Wald CIs calculated by sapply() loop
    ##    - Loops through each of fitted parameters and calculates lower/upper bounds
    ## Note: fitted pars must be on LINK scale
    ##     : if matrix is singular, none of the SEs or CIs are calculated (inherits/try statement)
    fittedPars = fit$par
    if(inherits(try(solve(hess), silent = TRUE), "try-error")) {
        ## Hessian is singular
        warning("Warning: singular hessian")
        ## SE and Wald CIs not calculated
        se = NA
        waldCI = matrix(NA, nrow = length(fittedPars), ncol = 2)
        ## But columns still need to be returned (if/when simulations are run)
        cnames = c("Estimate", "SE", "Lower", "Upper")
    } else {
        ## Calculating basic Wald CIs
        se = sqrt(diag(solve(hess)))
        waldCI = t(sapply(1:length(fittedPars),
                          function(i) fittedPars[i] + (c(-1, 1) * (qnorm(0.975) * se[i]))))
        ## Back-transforming the confidence limits, depending on whether we're using lambda0 or g0
        if(hn) {
            waldCI = rbind(plogis(waldCI[1, ]),
							exp(waldCI[2, ]),
                           exp(waldCI[3, ]))
        } else {
            waldCI = exp(waldCI)
        }
        
        ## Using the delta method to get the standard errors
        ## - G = jacobian matrix of partial derivatives of back-transformed
        ##    - i.e. log(D) -> exp(D) -- deriv. --> exp(D)
        ##    - Note: 1st deriv of plogis (CDF) = dlogis (PDF)
        if (hn){
            G.mult <- c(dlogis(fittedPars[1]),
						exp(fittedPars[2:3]))
        } else {
            G.mult <- exp(fittedPars)
        }
        G = diag(length(fittedPars)) * G.mult
        se = sqrt(diag(G %*% solve(hess) %*% t(G)))
    }
    
    
    ## Returning the fitted parameters in a named vector
    ## - First checks to see if TOA is being used,
    ##    then inserts par names in front of "sigma_toa"
    parNames = NULL
    if(!(is.null(toa) || is.na(toa))) {
        parNames = "sigma_toa"
    }
    if (hn) {
            parNames = c("g0", "sigma", parNames)
            fittedPars = c(plogis(fittedPars[1]),
                           exp(fittedPars[2]),
                           exp(fittedPars[3:length(fittedPars)]))
        }else {
        parNames = c("lambda0", "sigma", parNames)
        
        fittedPars = exp(fittedPars)
    }
    results = cbind(fittedPars, se, waldCI)
    dimnames(results) = list(parNames, c("Estimate", "SE", "Lower", "Upper"))
    results
}

