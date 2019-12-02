// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

/*
 
 Conditional probability
 calculates the coditional probability p(vector2|vector1), where vector 1 and 2 are two probability vectors
 of observing the pomo states in two diferent time points t1 and t2, respectively
 
 observed vectors: vector 1 and vector 2
 time elapsed between observing the vector 1 and 2: time = t2-t1
 rate_matrix: the pomo rate matrix, depending on N and sigma
 
 */

// [[Rcpp::export]]
mat exp_matrix(mat rate_matrix, double time) {
  mat pomo_matrix = expmat(rate_matrix*time); 
  return pomo_matrix;
}