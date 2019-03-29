// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <counts_to_moran_states.cpp>
#include <sigma_posterior2.cpp>

using namespace arma;
using namespace Rcpp;

/*

 Calculates the posterior of sigma for several sites
 Future plans: add a parallel for from Rcpp parallel
  
*/

// [[Rcpp::export]]
mat sites_sigma_posterior(List sync_file, int number_sites, double N, vec time, int number_time_points, int number_replicates, double prior_rate) {
  
  // usefull quantities
  mat output(number_sites,15);
  output.zeros();
  
  mat trajectories_matrix((N+1)*number_replicates,number_time_points);
  
  Mat<int> allele_counts = sync_file["allele1_counts"];
  Mat<int> total_counts  = sync_file["total_coverage"];
  
  Row<int> acounts(number_time_points*number_replicates);
  Row<int> tcounts(number_time_points*number_replicates); 
  
  for (int i = 0; i < number_sites; i++){
    cout << "Site: " << i+1 << "\n"; 
    // gets information from the allele and total counts for each replicates/time point
    acounts = allele_counts.row(i);
    tcounts = total_counts.row(i);
    
    // calculates the vitual trajectories
    trajectories_matrix = counts_to_moran_states(acounts,tcounts,number_time_points,number_replicates,N);
    
    // calculates summary statistics from the posterior of sigma  
    output.row(i) = sigma_posterior2(N,time,number_time_points,number_replicates,trajectories_matrix,prior_rate);
        
  }
  
  return output;
  
}

