// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <fstream>

#include <RcppArmadillo.h>

#include <counts_to_moran_states.cpp>
#include <sigma_posterior2.cpp>

using namespace arma;
using namespace Rcpp;
using namespace std;

/*

Calculates the posterior of sigma for several sites
Future plans: add a parallel for from Rcpp parallel

*/

// [[Rcpp::export]]
mat sites_sigma_posterior2(List sync_file, int initial_site, int final_site, string output_file,  double N, vec time, int number_time_points, int number_replicates, double prior_rate) {
  
  ofstream file_output;
  file_output.open(output_file);
  
  // usefull quantities
  mat output(final_site-initial_site+1,15);
  output.zeros();
  
  mat trajectories_matrix((N+1)*number_replicates,number_time_points);
  
  StringVector chr  = sync_file["chromosome"];
  Col<int> pos      = sync_file["position"];
  StringVector ref  = sync_file["reference_allele"];
  StringVector a1   = sync_file["allele1"];
  StringVector a2   = sync_file["allele2"];
  
  Mat<int> allele_counts = sync_file["allele1_counts"];
  Mat<int> total_counts  = sync_file["total_coverage"];
  

  Row<int> acounts(number_time_points*number_replicates);
  Row<int> tcounts(number_time_points*number_replicates); 
  
  for (int i = (initial_site-1); i < final_site; i++){
    cout << "Site: " << i+1 << "\n"; 
    // gets information from the allele and total counts for each replicates/time point
    acounts = allele_counts.row(i);
    tcounts = total_counts.row(i);
    
    // calculates the vitual trajectories
    trajectories_matrix = counts_to_moran_states(acounts,tcounts,number_time_points,number_replicates,N);

    // calculates summary statistics from the posterior of sigma  
    int j = i-initial_site+1;
    output.row(j) = sigma_posterior2(N,time,number_time_points,number_replicates,trajectories_matrix,prior_rate);
    
    file_output << chr(i)       << "\t" << pos(i)       << "\t" << ref(i)       << "\t" << a1(i)        << "\t"<< a2(i)         << "\t"
                << output(j,0)  << "\t" << output(j,1)  << "\t" << output(j,2)  << "\t" << output(j,3)  << "\t" << output(j,4)  << "\t"
                << output(j,5)  << "\t" << output(j,6)  << "\t" << output(j,7)  << "\t" << output(j,8)  << "\t" << output(j,9)  << "\t" 
                << output(j,10) << "\t" << output(j,11) << "\t" << output(j,12) << "\t" << output(j,13) << "\t" << output(j,14) << "\n";
  }
  
  return output;
  
}

