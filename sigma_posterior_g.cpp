// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <log_posterior.cpp>
#include <quantile_interpolation.cpp>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
List sigma_posterior(double N, vec range_sigma, int number_replicates, int number_time_points,  int delta, mat trajectories_matrix, vec times, double prior_rate) {
  
  vec sigma = linspace(range_sigma(0),range_sigma(1),delta);
  vec posterior_sigma(delta);
  
  for(int i=0; i<delta; i++) {
    posterior_sigma(i) = log_posterior(N,sigma(i),trajectories_matrix,number_replicates,number_time_points,time,prior_rate);
    std::cout << i << " of " << delta << " \n";
  }
  
  vec prob_posterior_sigma = exp(posterior_sigma);
  prob_posterior_sigma = prob_posterior_sigma/sum(prob_posterior_sigma);
  
  double sigma_mean   =  sum(sigma%prob_posterior_sigma);
  double sigma_median =  quantile_interpolation(sigma,prob_posterior_sigma,0.5);
  double sigma_05     =  quantile_interpolation(sigma,prob_posterior_sigma,0.05);
  double sigma_95     =  quantile_interpolation(sigma,prob_posterior_sigma,0.95);
  
  vec credibility_interval(2);
  credibility_interval(0) = sigma_05;
  credibility_interval(0) = sigma_95;
  
  return List::create(Named("sigma") = sigma,
                      Named("sigma_log_posterior") = posterior_sigma,
                      Named("sigma_mean") = sigma_mean,
                      Named("sigma_median")=sigma_median,
                      Named("sigma_credibility_interval")= credibility_interval);
}

