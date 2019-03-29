// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <pomo_rate_matrix.cpp>

using namespace arma;
using namespace Rcpp;


/*
 Prepares the prob_matrix file for the analysis, 
 avoids calculating the exp matrix an unecessary number of times if the incremente between any two time points is equal
*/

// [[Rcpp::export]]
List prob_matrix(double N, double sigma, vec time, int number_time_points) {
  
  // rate matrix
  mat rate_matrix = pomo_rate_matrix(N,sigma);
  
  // calculates increments
  vec increments(number_time_points-1);
  vec indexes(number_time_points-1);
  indexes.zeros();
  
  for (int i=0; i<(number_time_points-1); i++){
    increments(i) = time(i+1) - time(i);
  }
  
  // evaluates prob_matrix for all the unique increments
  vec unique_increments = unique(increments);
  int number_unique_increments = unique_increments.n_elem;
  
  mat t_prob_matrix_keeper((N+1)*number_unique_increments,N+1);
  for (int i=0; i < number_unique_increments; i++){
      double delta = unique_increments(i);
      t_prob_matrix_keeper(span(i*N+i,(i+1)*N+i),span(0,N)) = expmat(rate_matrix*delta).t();
      indexes.elem(find(increments == delta)) += i;
  }
  
  return List::create(Named("index") = indexes,
                      Named("prob_matrix") = t_prob_matrix_keeper);
}

