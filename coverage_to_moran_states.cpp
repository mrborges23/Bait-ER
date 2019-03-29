// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

/*
 Calculates quantiles using interpolation
 */

// [[Rcpp::export]]
vec coverage_to_moran_states(double N, double coverage_allele_A, double total_coverage) {
  
  int number_pomo_states = N+1;
  vec conditional_probability(N+1);
  conditional_probability(N+1).zeros();
  
  for (int i=0; i< number_pomo_states; i++ ) {
    
    for (int j=0; j<C; j++) {
      
      
    }
    
    
  }
  
  
}

