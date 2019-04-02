// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

/*
 
 Conditional probability
 calculates the coditional probability p(vector2|vector1), where vector 1 and 2 are two probability vectors

 observed vectors: vector 1 and vector 2
 prob_matrix: the transpose of the pomo rate matrix, in the increment of time between t1 and t2
 
 */

// [[Rcpp::export]]
vec trajectory_probability(mat prob_matrix, vec vector1, vec vector2) {
  // calculates the conditional probabilities of observing the pomo state vector vector2 in time t2 
  // given that the pomo state vector in t1 was vector1
 
  vec pomo_vector = prob_matrix*vector1; 
  vec probability_vector = pomo_vector%vector2;
  
  return probability_vector;
}
