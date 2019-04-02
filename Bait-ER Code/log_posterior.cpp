// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <log_likelihood_trajectories.cpp>

using namespace arma;
using namespace Rcpp;


/*
 Calculates de log posterior
 used the log_likelihood and adds the log_prior
 prior of sigma is exponential with rate prior_rate

 prior_rate: the rate of the exponential distribution, prior of sigma
 
*/

// [[Rcpp::export]]
double log_posterior(double N, double sigma, mat trajectories_matrix,int number_replicates, int number_time_points, vec time, vec prior_parameters) {
  
  // adds the log likelihood to the log prior
  double log_total = R::dgamma(1.0+sigma, prior_parameters(0), 1/prior_parameters(1), 1) + log_likelihood_trajectories(N,sigma,number_replicates,number_time_points,time,trajectories_matrix);
  
  return log_total;
}
