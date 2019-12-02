// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <moran_states_distribution.cpp>

using namespace arma;
using namespace Rcpp;
using namespace std;

/*

 Calculates the probability of the moran states given an obaerved count 
 p({nA,(N-N)a}!c)=Binomial(n/N,C) for given set of time points, replicates 
 C
 the vectors allele_counts and total_counts should be organized first per generation 
 and then per replicate. Example for 2 generations and 3 replicates:
 F0R1 F0R2 F0R3 F1R1 F1R2 F1R3
 
*/

// [[Rcpp::export]]
mat counts_to_moran_states(Row<int> allele_counts, Row<int> total_counts,int number_time_points, int number_replicates, double N) {

  // usefull quantitie
  mat allele_trajectories((N+1)*number_replicates,number_time_points);
  allele_trajectories.zeros();
    
  vec moran_distribution(N+1); 
  moran_distribution.zeros();
    
  // calculates the probability of each moran state given the observed counts
  int index = 0; 
    
  for (int t=0; t<number_time_points; t++){
    for (int r=0; r<number_replicates; r++){
      
      double allele_coverage = allele_counts(index);
      double total_coverage  = total_counts(index);
      
      moran_distribution = moran_states_distribution(N, allele_coverage, total_coverage);
      allele_trajectories(span(r*N+r,(r+1)*N+r),t) = moran_distribution;
      index = index+1;
    }
  }
  
  return allele_trajectories;
  
}
