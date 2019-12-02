// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

/*
 Calculates the probability of the moran states given a value of coverage
 p({nA,(N-N)a}!c)=Binomial(n/N,C)
 
 N: population size
 allele_coverage: coverage of allele A 
 total_coverage: total coverage
 
 Example:
 N               <-100
 allele_coverage <- 0
 total_coverage  <- 20
 moran_states_distribution(N,allele_coverage,total_coverage)
 */

// [[Rcpp::export]]
NumericVector moran_states_distribution(double N, int allele_coverage, int total_coverage) {
  
  //usefull quantities
  int number_pomo_states =  N+1;
  NumericVector moran_distribution(number_pomo_states);
  
  // populates the moran distribution using the inverse binomial sampling
  for (int i=0; i<number_pomo_states;i++){
    moran_distribution(i) = R::dbinom(allele_coverage,total_coverage,i/N,false);
  }
  
  moran_distribution = moran_distribution/sum(moran_distribution);
  
  return moran_distribution;
}
