// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

/*
PoMo rate matrix
continuous time: may be usefull in some cases
without mutation: can be ignored for few generations
with selection 
stationarity is not assumed

population size: N
selection coefficient: sigma

*/

// [[Rcpp::export]]
mat pomo_rate_matrix(double N, double sigma) {
  
  // useful quantities
  int number_pomo_states = N+1;
  mat pomo_matrix(number_pomo_states,number_pomo_states);
  pomo_matrix.zeros();
  
  // populates the matrix
  // pomo states: {Na}, {1A,(N-1)a}, ..., {(N-1)A,1a}, {NA}
  for (int i=1; i<(number_pomo_states-1); i++ ){
    
    //allele A increases in frequency (a decreases)
    pomo_matrix(i,i+1) = i*(N-i)*(1+sigma)/N;
    
    //allele A decreases in frequency (a increases)
    pomo_matrix(i,i-1) = i*(N-i)/N;
    
  }
  
  // populates the diagonal elements
  pomo_matrix.diag() = -sum(pomo_matrix, 1); 
  
  return pomo_matrix;
}

