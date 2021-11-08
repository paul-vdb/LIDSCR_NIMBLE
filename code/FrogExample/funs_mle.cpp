#include <math.h>
#include <iostream>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;



// ============================== //
//          Stacked-AC NLL          //
// ============================== //
/*
 * Calculates the log-likelihood for cue-based (acoustic) SCR data
 * - Arguments:
 *    - pars:       Vector of (density, g0/lambda0, sigma, lambda_c [, sigma_toa]).
 *    - caps:       Matrix of captures; last column is animal ID
 *    - nTraps:     Number of detectors.
 *    - toa_ssq:    Matrix of times of arrival for each mask point.
 *    - nMask:      Number of mask points.
 *    - aMask:      Area covered by a single mask point.
 *    - maskDists:  Matrix of distances between each trap and mask point.
 *    - nCalls:     Vector of counts; number of calls by each animal ID.
 *    - use_toa:    TRUE/FALSE; whether TOA matrix is used.
 */
// [[Rcpp::export]]
double scr_nll_stacked(const NumericVector& pars,
		       const arma::mat& caps,
		       const double& aMask,
		       const arma::mat& maskDists,
		       const arma::vec& ID,
		       const arma::mat& toa,
		       const arma::mat& toa_ssq,
		       const bool& use_toa,
		       const bool& hn,
		       const bool& trace) {
  /*
   *  Storing/initialising (starting) parameter values.
   *  - Note that parameters are back-transformed
   *  - Also note that if use_toa = FALSE, there will be no sigma_toa to estimate
   */
  double D = exp(pars[0]);
  double g0;
  /*
   * Detection function intercept is on [0, 1] using halfnormal, [0,
   * Inf) if hazard halfnormal.
   */
  if (hn){
    g0 = R::plogis(pars[1], 0, 1, 1, 0);
  } else {
    g0 = exp(pars[1]);
  }
  double sigma = exp(pars[2]);
  double lambda_c = exp(pars[3]);
  double sigma_toa = 0;
  if (use_toa){
    sigma_toa = exp(pars[4]);
  }
  // Number of animals, calls, mask points
  int nCalls = caps.n_rows;
  arma::vec animal_ID = unique(ID);
  int nAnimals = animal_ID.n_rows;
  
  int nMask = maskDists.n_rows;
  int nTraps = caps.n_cols;
  
  /*
  * Constructing a detection probability matrix.
  * - Element (i, j) gives prob. of animal @ ith mask pt. being detected @ jth trap.
  * - Line that fills in maskProbs(i, j) is the Hazard Half-Normal function (HHN)
  */
  arma::mat maskProbs(nMask, nTraps);
  for(int i = 0; i < nMask; i++) {
    for(int j = 0; j < nTraps; j++) {
      if (hn){
	maskProbs(i, j) = g0 * exp(-pow(maskDists(i, j), 2.0) / (2 * pow(sigma, 2.0))) + DBL_MIN;
      } else {
	maskProbs(i, j) = 1 - exp(-(g0 * exp(-pow(maskDists(i, j), 2.0) / (2 * pow(sigma, 2.0))))) + DBL_MIN;
      }
    }
  }

  /*
  * Constructing a detection probability vector
  * - ith element = P(animal @ ith mask pt. is detected by >= 1 trap)
  */
  arma::vec pAvoid(nMask);
  for(int i = 0; i < nMask; i++) {
    pAvoid[i] = 1 - maskProbs(i, 0);
    for(int j = 1; j < nTraps; j++) {
      pAvoid[i] *= 1 - maskProbs(i, j);
    }
  }

  /*
  *  Probability of detecting one specific call emitted from s
  *
  */
  arma::vec pDetected = 1 - pAvoid;

  /*
  * Probability of detecting at least one call on at least one microphone
  *  from an individual located at s
  */
  arma::vec pAnimal = 1 - exp(-lambda_c * pDetected);

  // ========================================= //
  // ========================================= //
  arma::mat fCapt(nMask, nAnimals);
  arma::mat lpTrap = log(maskProbs + DBL_MIN);
  arma::mat lpMiss = log(1-maskProbs + DBL_MIN);
  arma::vec lpDetected = log(pDetected + DBL_MIN);

  // Looping through all calls	
  for (int i = 0; i < nCalls; i++){

     // Bernoulli capture history
     //----------------------------
     fCapt.col(ID(i) - 1) += lpTrap*caps.row(i).t() + lpMiss*(1-caps.row(i)).t(); 
     fCapt.col(ID(i) - 1) -= lpDetected;

	 // Time of Arrival:
     //----------------------------
	 arma::uvec mic_obs = find(caps.row(i) > 0);
	 int mi = mic_obs.n_elem;
     
	 if(mi != 1){ 
	   fCapt.col(ID(i) - 1) += (1 - mi) * log(sigma_toa/1000) - (toa_ssq.row(i).t() / (2 * pow(sigma_toa/1000, 2)));
	 }
    
	 // Add the count process for number of calls:
	 fCapt.col(ID(i) - 1) += log(lambda_c) + lpDetected; // For each poisson this is the additional part per call.
  }

    // Now we add the full poisson bit to the whole thing...
	fCapt -= arma::repmat(pDetected, 1, nAnimals)*lambda_c;
	
	// Now approximate the integral:
	arma::mat efCapt = exp(fCapt);
	arma::rowvec sum_fCapt = sum(efCapt, 0);	// Don't need to add the area of the mask as it cancels with in the f_x term.
	
    /*
    * Log-likelihood contribution from all capture histories
    * - Calculated by log of sum of individual likelihood contributions.
    */
    double log_lik = sum(log(sum_fCapt + DBL_MIN));

    // Calculating effective survey area.
    double esa = aMask * sum(pAnimal);

    // Log-likelihood contribution from number of animals detected.
    log_lik += R::dpois(nAnimals, D * esa, 1);

    // Add animal locations probabilties:
	log_lik -= nAnimals * log(sum(pAnimal));
	
   // Printing out parameter values and log-likelihood, if required.
   if (trace){
     Rcout << "D: " << D << ", g0: " << g0 << ", sigma: " << sigma << ", lambda_c: " << lambda_c << ", sigma_toa: " << sigma_toa << ", LL: " << log_lik << std::endl;
   }
   // Returning log-likelihood
   return -log_lik;
}



// ============================== //
//          Traditional ASCR      //
// ============================== //
/*
 * Calculates the log-likelihood for cue-based (acoustic) SCR data
 * - Arguments:
 *    - pars:       Vector of (density, g0/lambda0, sigma, lambda_c [, sigma_toa]).
 *    - caps:       Matrix of captures; last column is animal ID
 *    - nTraps:     Number of detectors.
 *    - toa_ssq:    Matrix of times of arrival for each mask point.
 *    - nMask:      Number of mask points.
 *    - aMask:      Area covered by a single mask point.
 *    - maskDists:  Matrix of distances between each trap and mask point.
 *    - nCalls:     Vector of counts; number of calls by each animal ID.
 *    - use_toa:    TRUE/FALSE; whether TOA matrix is used.
 */
// [[Rcpp::export]]
arma::mat ascr_calls(NumericVector pars,
		       const arma::mat& caps,
		       const double& aMask,
		       const arma::mat& maskDists,
		       const arma::mat& toa,
		       const arma::mat& toa_ssq,
		       const bool& use_toa,
		       const bool& hn) {
  /*
   *  Storing/initialising (starting) parameter values.
   *  - Note that parameters are back-transformed
   *  - Also note that if use_toa = FALSE, there will be no sigma_toa to estimate
   */
  double g0;
  /*
   * Detection function intercept is on [0, 1] using halfnormal, [0,
   * Inf) if hazard halfnormal.
   */
  if (hn){
    g0 = R::plogis(pars[0], 0, 1, 1, 0);
  } else {
    g0 = exp(pars[0]);
  }
  double sigma = exp(pars[1]);
  double sigma_toa = 0;
  if (use_toa){
    sigma_toa = exp(pars[2]);
  }
  // Number of animals, calls, mask points
  int nCalls = caps.n_rows;
  
  int nMask = maskDists.n_rows;
  int nTraps = caps.n_cols;
  
  /*
  * Constructing a detection probability matrix.
  * - Element (i, j) gives prob. of animal @ ith mask pt. being detected @ jth trap.
  * - Line that fills in maskProbs(i, j) is the Hazard Half-Normal function (HHN)
  */
  arma::mat maskProbs(nMask, nTraps);
  for(int i = 0; i < nMask; i++) {
    for(int j = 0; j < nTraps; j++) {
      if (hn){
	maskProbs(i, j) = g0 * exp(-pow(maskDists(i, j), 2.0) / (2 * pow(sigma, 2.0))) + DBL_MIN;
      } else {
	maskProbs(i, j) = 1 - exp(-(g0 * exp(-pow(maskDists(i, j), 2.0) / (2 * pow(sigma, 2.0))))) + DBL_MIN;
      }
    }
  }

  /*
  * Constructing a detection probability vector
  * - ith element = P(animal @ ith mask pt. is detected by >= 1 trap)
  */
  arma::vec pAvoid(nMask);
  for(int i = 0; i < nMask; i++) {
    pAvoid[i] = 1 - maskProbs(i, 0);
    for(int j = 1; j < nTraps; j++) {
      pAvoid[i] *= 1 - maskProbs(i, j);
    }
  }

  /*
  *  Probability of detecting one specific call emitted from s
  *
  */
  arma::vec pDetected = 1 - pAvoid;

  // ========================================= //
  // ========================================= //
  arma::mat fCapt(nMask, nCalls);
  arma::mat lpTrap = log(maskProbs + DBL_MIN);
  arma::mat lpMiss = log(1-maskProbs + DBL_MIN);
  arma::vec lpDetected = log(pDetected + DBL_MIN);

  // Looping through all calls	
  for (int i = 0; i < nCalls; i++){

     // Bernoulli capture history
     //----------------------------
     fCapt.col(i) += lpTrap*caps.row(i).t() + lpMiss*(1-caps.row(i)).t(); 
     // fCapt.col(i) -= lpDetected; // This term cancels with f_x numerator and is not included for now.

	 // Time of Arrival:
     //----------------------------
	 arma::uvec mic_obs = find(caps.row(i) > 0);
	 int mi = mic_obs.n_elem;
     
	 if(mi != 1){ 
	   fCapt.col(i) += (1 - mi) * log(sigma_toa/1000) - (toa_ssq.row(i).t() / (2 * pow(sigma_toa/1000, 2)));
	 }
   }
    // Calculating effective survey area.
    double esa = aMask * sum(pDetected);

    // Now we add f_x (Cancels with truncation term.
	// fCapt += arma::repmat(pDetected, 1, nCalls)*lambda_c;
	
	// Denominator for f_x
	fCapt -= log(esa);

   // Returning log-likelihood
   return fCapt;
}
// =================================================================================== //
// =================================================================================== //


// A function to create the TOA SSQ matrix.
// [[Rcpp::export]]
NumericMatrix make_toa_ssq(const NumericMatrix& capt, const NumericMatrix& dists, const double& sound_speed){
  int n = capt.nrow();
  int n_traps = capt.ncol();
  int n_mask = dists.ncol();
  NumericMatrix out(n, n_mask);
  int n_dets;
  int index;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n_mask; j++){
      n_dets = 0;
      for (int k = 0; k < n_traps; k++){
	if (capt(i, k) > 0) n_dets++;
      }
      NumericVector delts(n_dets);
      index = 0;
      for (int k = 0; k < n_traps; k++){
	if (capt(i, k) > 0){
	  delts(index) = capt(i, k) - dists(k, j)/sound_speed;
	  index++;
	}
	out(i, j) = sum(pow(delts - mean(delts), 2));
      }
    }
  }  
  return out;
}


/*
 * Calculating the Euclidean distance between a point and each trap.
 * Returns a vector of distances.
 */
// [[Rcpp::export]]
NumericMatrix eucdist(NumericMatrix points,
                          NumericMatrix traplocations) {
  NumericMatrix dists(points.nrow(), traplocations.nrow());
  for(int i = 0;  i < points.nrow(); i++) {
    for(int j = 0; j < traplocations.nrow(); j++) {
      dists(i, j) = sqrt(pow(points(i, 0) - traplocations(j, 0), 2.0)
                         + pow(points(i, 1) - traplocations(j, 1), 2.0));
    }
  }
  return dists;
}
