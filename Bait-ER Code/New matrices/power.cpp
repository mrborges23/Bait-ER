// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <log_posterior.cpp>
#include <quantile_interpolation.cpp>

using namespace arma;
using namespace Rcpp;

/*
Calculates the power matrix
*/

// [[Rcpp::export]]
mat power_m(double N, mat matrix) {
  
  mat p_matrix = pow(matrix,N);
  return p_matrix;
  
}

