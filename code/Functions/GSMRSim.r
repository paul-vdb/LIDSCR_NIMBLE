##################################################
# Add function from gSMR paper:
##################################################
e2dist <- function (x, y) {  # Function from scrbook package to calculate the distance between locations in 2 matrices.
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

e2dist.2 <- function (x, y) {  # Function from scrbook package to calculate the squared distance between locations in 2 matrices.
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- (x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

###############################################################-
### 1.  Function to generate generalized SMR data with linear, gridded, or random trap (within the state-space), trap designs ####
###############################################################-

fnc.create.SMR.data <- function (N = N, K.trap = K.trap, K.camera = K.camera, sigma = sigma, g0trap = g0trap, lam0 = lam0, 
                                 n.trap = n.trap, n.camera = n.camera, trap.design = trap.design,
                                 xlims = xlims, ylims = ylims, xlim.R = xlim.R, ylim.R = ylim.R, obsmod = c("pois", "bern"),
                                 n.collar = 10000, nlocs = 100, plot.map = TRUE) {
  # N = True number of animals in the state-space
  # K.trap = number of marking (trapping) occasions
  # K = number of resighting (eg camera) occasions
  # sigma = scale parameter for home range size
  # g0trap = probability of detection at the home range centre for marking
  # lam0 = encounter rate at home range centre for resighting
  # n.trap = total number of traps for marking.  Used to create a square grid of traps
  # n.camera = total number of cameras for resighting.  Used to create a square grid of detectors
  # trap.design = one of 'linear', 'grid', 'random'.  While we ran simulations with 'random', we did not include results in manuscript
  # xlims = min and max coordinates of easting for the state-space
  # ylims = min and max coordinates of northing for the state-space
  # xlim.R = min and max coorindates of easting for the study area
  # ylim.R = min and max coordinates of northing for the study area.
  # obsmod = type of observation model for the resighting process.  We used 'pois' for the manuscript.
  # n.collar = maximum number of animals to fit with telemetry tags or collars.  Telemetry data is created for the lesser of n.collar or n.marked
  # nlocs = number of telemetry locations for each marked animal.  We used 100 for each marked animal in the manuscript

  # VALUES RETURNED
  # y.trap.all = all capture observations including known animals that were undetected.  Array of individuals x trap x trap.occasion
  # y.trap.marked = SCR observations of captured animals.  Array of individuals x trap x trap.occasion
  # y.camera.all = all resight observations including animals not detected and/or not trapped.  Array of individuals x detector x occasion
  # y.camera.marked = resight encounter history of marked animals.  Array of individuals x detector x occasion
  # n.camera.unmarked = number of detections of unmarked animals at location j (row) and occasion k (column)
  # i.marked = row numbers of marked animals in y.trap.all and y.camera.all
  # i.unmarked = row numbers of unmarked animals in y.camera.all
  # telemetry.array = Telemetry locations of marked individuals.  Array of individual x location number x coordinates(x,y).
  #                   # When telemetry data is collected from a subset of marked individuals.  The first n.collar marked individuals receive telemetry data.
  # X.trap = locations of traps for marking.  Matrix of trap x coordinates
  # X.camera = locations of detectors for resighting.  Matrix of detector x coordinates

  if ( ! trap.design %in% c('grid', 'linear', 'random')){
    stop('trap.design should be one of:  grid, linear, or random')
  }
  
  obsmod <- match.arg(obsmod)
  n.trap.row <- sqrt(n.trap)
  n.camera.row <- sqrt(n.camera)
  if (trap.design == 'grid'){
    coor0 <- seq(xlim.R[1], xlim.R[2], length = sqrt(n.trap))
    X.trap <- cbind(rep(coor0, each=length(coor0)),rep(coor0, times=length(coor0))) # Nets for a square grid of traps
  } 
  
  if (trap.design == 'linear'){
    coor0 <- seq(xlim.R[1], xlim.R[2], length = n.trap)
    X.trap <- cbind(coor0, 0.5)
  }
  
  if (trap.design == 'random') {
    n.rand <- n.trap
    X.trap <- cbind( runif(n.rand, xlims[1], xlims[2]), runif(n.rand, xlims[1], xlims[2])    )
  }
  J.trap <- nrow(X.trap)  # nN
  
  ## Camera coordinates
  coor0 <- seq(xlim.R[1], xlim.R[2], length = n.camera.row)
  X.camera <- cbind(rep(coor0, each=length(coor0)),rep(coor0, times=length(coor0))) # Nets for a square grid of traps
  t.jitter <- (coor0[2] - coor0[1])/3
  J.camera <- nrow(X.camera)  # nN

    # Activity Centers
    sx <- runif(N, xlims[1], xlims[2])
    sy <- runif(N, ylims[1], ylims[2])
    S <- cbind(sx, sy)
    #### MARK DATA
    D.trap <- e2dist(S, X.trap)
    ptrap <- g0trap * exp(-(D.trap * D.trap)/(2 * sigma * sigma))
    y.trap.all <-array(0, c(N, J.trap, K.trap))
    for (i in 1:N){
      for (j in 1:J.trap){
        y.trap.all[i, j, ] <- rbinom(K.trap, 1, ptrap[i, j])
      }
    }
    n.trap.ind <- apply(y.trap.all, 1, sum)    # number of captures per individual
    marked <- ifelse(n.trap.ind > 0, 1, 0) # is each animal marked (0 or 1)
    i.marked = (1:N)[marked == 1]          # ID for marked individuals
    n.marked <- sum(marked)                # number captured and marked
    y.trap <- y.trap.all[marked == 1, , ]  # capture-recapture data for marked animals

    
    #### RESIGHT DATA
    D <- e2dist(S, X.camera)               # Distance between each home range center (row) and camera (column)
    lam <- lam0 * exp(-(D * D)/(2 * sigma * sigma))   # Encounter rates
    y.camera.all <- array(NA, c(N, J.camera, K.camera))          # Array for resighting data
    for (i in 1:N) {
      for (j in 1:J.camera) {
        if (identical(obsmod, "bern")) {
          y.camera.all[i, j, ] <- rbinom(K.camera, 1, lam[i, j])
        }
        else if (identical(obsmod, "pois")) {
          y.camera.all[i, j, ] <- rpois(K.camera, lam[i, j])
        }
      }
    }
    y.camera.unmarked <- y.camera.all * (1 - marked)
    i.unmarked <- (1:N)[rowSums(y.camera.unmarked) > 0] 
    n.camera.unmarked <- apply(y.camera.unmarked, c(2, 3), sum) # Sum of detections by Camera (row) and Occasion (column)
    y.camera.marked <- y.camera.all[marked == 1, , ]       # Resight data of marked animals
    n.ind <- apply(y.camera.all, c(1), sum)               # Number of detections per individual
    det.camera <- ifelse(n.ind > 0, 1, 0)                  # was each individual detected by camera.  Used in plot below.

    # Telemetry data
    n.collar <- min(c(n.marked, n.collar))
    telemetry.array <- array(NA, dim=c(n.collar, nlocs, 2))
    if (nlocs > 0 & n.collar > 0) {
      for (i in 1:n.collar) {
        telemetry.array[i, , 1] <- rnorm(nlocs, S[i.marked[i], 1], sigma)
        telemetry.array[i, , 2] <- rnorm(nlocs, S[i.marked[i], 2], sigma)
      }        
    }
    
    # Plot marked and unmarked animals
    if (plot.map == TRUE){
      par(mfrow = c(1,1))
      plot(S, col = 'red', pch = 19, xlim = xlims, ylim = ylims, asp = 1, xlab = 'X', ylab = 'Y',
           main = paste('Spatial Mark-Resight Simulated Data \n g0trap =', g0trap,
                        ' lam0camera =', lam0, ' sigma =', sigma, ' N =', N, ' \n Trap Design:', trap.design,  sep = ' '))
      rect(xleft = xlims[1], xright = xlims[2], ybottom = ylims[1], ytop = ylims[2], col = 'gray80' )  # xlim.S = c(0, 1), ylim.S = c(0, 1), xlim.R = c(0.2, 0.8), ylim.R = c(0.2, 0.8)
      rect(xleft = xlim.R[1], xright = xlim.R[2], ybottom = ylim.R[1], ytop = ylim.R[2], col = 'gray70' )  # xlim.S = c(0, 1), ylim.S = c(0, 1), xlim.R = c(0.2, 0.8), ylim.R = c(0.2, 0.8)
      lapply(seq(0, 1, by = 0.1), function(x) lines(c(x, x), c(0, 1), col = 'white'))
      lapply(seq(0, 1, by = 0.1), function(x) lines(c(0, 1), c(x, x), col = 'white'))
      points(S, col = 'red', pch = 19)
      points(X.trap, col='blue', pch = 1, cex = 1.2)
      points(X.camera, col = 'forestgreen', pch = '+', cex = 0.8)
      points(S[marked == 1, 1], S[marked == 1, 2], pch=19, col='purple', cex = 1)
      points(S[det.camera == 1, 1], S[det.camera == 1, 2], pch=1, col='yellow', cex = 0.5)
      legend(x = 0.7, y = 1, cex = 0.7, legend = c( 'Traps for Marking', 'Cameras', 'Home Range Centre', '  Marked', '  Detected Camera'), pch = c(1, 3, 19, 19, 1), col = c('blue','forestgreen', 'red','purple', 'yellow' ), bg = 'gray80')
    } # End of Plot
    # Collect Results
    list(y.trap.all = y.trap.all, y.trap = y.trap, i.marked = i.marked, y.camera.all = y.camera.all, y.camera.marked = y.camera.marked, i.unmarked = i.unmarked,
         n.camera.unmarked = n.camera.unmarked,
         telemetry.array = telemetry.array, X.trap = X.trap, X.camera = X.camera)

}   # END OF FUNCTION  fnc.create.SMR.data
