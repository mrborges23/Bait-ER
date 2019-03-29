// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <moran_states_distribution.cpp>

using namespace arma;
using namespace Rcpp;
/*
 
 
 
 */
// [[Rcpp::export]]
mat reads_to_moran_states(double N, List reads) {
  
  // reads the list
  mat allele_coverage = reads["allele_coverage"];
  mat total_coverage  = reads["total_coverage"];

  // useful quantitites
  int number_pomo_states = N+1;
  int number_replicates  = allele_coverage.n_rows;
  int number_time_points = allele_coverage.n_cols;

  mat paths(number_pomo_states*number_replicates,number_time_points);
  
  // populates paths
  for (int i=0; i<number_time_points; i++){
    
    for (int j=0; j<number_replicates; j++){
      
      vec distribution = moran_states_distribution(N,allele_coverage(j,i),total_coverage(j,i));
      paths(span(j*N+j,(j+1)*N+j),i) = distribution;
        
    }
    
  }

  return paths;

}
