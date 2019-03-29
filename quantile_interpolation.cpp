// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

/*
Calculates quantiles using interpolation
 */
// [[Rcpp::export]]
double quantile_interpolation(vec sigma, vec prob_post,double quantile) {
  
  vec c_prob_post = cumsum(prob_post);
  vec a_prob_post = abs(c_prob_post-quantile);
  uvec index      = sort_index(a_prob_post);

  double sigma1 = sigma(index(0));
  double sigma2 = sigma(index(1));
  
  double c_prob_post1 = c_prob_post(index(0));
  double c_prob_post2 = c_prob_post(index(1));
  
  double sigma_asterisk = (quantile-c_prob_post1)*(sigma2-sigma1)/(c_prob_post2-c_prob_post1)+sigma1;
  return sigma_asterisk;
}
