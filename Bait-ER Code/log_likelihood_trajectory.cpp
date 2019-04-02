// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <trajectory_probability.cpp>

using namespace arma;
using namespace Rcpp;

/*
 Trajectory likelihood
 calculates the likelihood of observing several allele trajectory within the same replicate
 it accounts for noise coming from the coverage readings

 data consists of state vector containing probabilities, each column representing one time point
 and each row represents a pomo state
 
 Parameters:
 N: population size
 sigma: selection coefficient
 trajectory_matrix: congregates the vector of observed states
 times:

 */

// [[Rcpp::export]]
double log_likelihood_trajectory(double N, mat trajectory_matrix, vec indexes, mat prob_matrices) {
  
  // useful quantities
  vec vector1 =trajectory_matrix.col(0);
  int number_time_points = trajectory_matrix.n_cols;
  
  for (int i=0; i<(number_time_points-1); i++){
    
    // nth observed vector
    vec vector2 = trajectory_matrix.col(i+1);
    
    // gets the transpose prob_matrix
    int index = indexes(i);
    mat prob_matrix = prob_matrices(span(index*N+index,(index+1)*N+index),span(0,N));
    
    // calculates conditional probabilities
    vector1 = trajectory_probability(prob_matrix,vector1,vector2);
  }
  
  // returns the log conditional probability of observing the trajectory given 
  // a particular selection coefficient sigma 
  double log_probability = log(sum(vector1));
  return log_probability;
}
