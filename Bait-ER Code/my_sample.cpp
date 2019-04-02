// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

/*
Simple verison of function sample
 */

// [[Rcpp::export]]
int my_sample(vec pomo_states, vec prob) {
  
  // normalizes the prob vector
  vec prob1 = prob/sum(prob);
  
  // sample a value using the uniform distribution
  NumericVector u = runif(1);
  
  // cumulative probability
  vec prob2 = cumsum(prob1);
  
  // returns the value sampled
  uvec index = find(prob2>u(0));
  return pomo_states(index(0));
  
}

