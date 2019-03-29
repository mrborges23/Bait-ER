// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <log_likelihood_trajectory.cpp>
#include <prob_matrix.cpp>

using namespace arma;
using namespace Rcpp;


/*

 Calculates log likelihood over several replicates
 
 trajectories_matrix: columns represent each time point
 rows represent the pomo states
 the replicates are row binded

*/

// [[Rcpp::export]]
double log_likelihood_trajectories(double N, double sigma, int number_replicates,int number_time_points, vec time, mat trajectories_matrix) {
  
  List list_prob_matrix = prob_matrix(N,sigma,time,number_time_points);
  
  vec indexes = list_prob_matrix(0);
  mat prob_matrices = list_prob_matrix(1);
  
  vec log_likelihood_replicates(number_replicates);
  log_likelihood_replicates.zeros();
  
  for (int i=0; i < number_replicates; i++){
    
    // preparing trajectories_matrix, for replicate number i
    mat trajectory_matrix = trajectories_matrix.rows(i*N+i,(i+1)*N+i);
    
    // calculates the likelihood of each trajectory
    log_likelihood_replicates(i) = log_likelihood_trajectory(N, trajectory_matrix, indexes, prob_matrices);
  }
  
  // total likelihood
  return sum(log_likelihood_replicates);
}


      