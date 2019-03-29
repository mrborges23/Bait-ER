// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <split.cpp>
#include <moran_states_distribution.cpp>
#include <quantile_interpolation.cpp>
#include <log_posterior.cpp>

using namespace arma;
using namespace Rcpp;
using namespace std;

/*
 Reads sync file
 
 Sync files contain 3 + n columns with: 
 col 1: chromosome (reference contig)
 col 2: position (in the reference contig)
 col 3: reference allele 
 col >3: sync entries for allele frequencies for all populations in the form 
 A-count:T-count:C-count:G-count:N-count:deletion-count. 
 
 Sync files originally don’t have a header but headers are accepted when specified with header=T.
 Sync files are  specified in Kofler et al. (2011). 
 
 
 This code was take from https://gist.github.com/hadley/6353939.
 */


// [[Rcpp::export]]
DataFrame read_sync_file_g(string path,int number_sites,int number_time_points, int number_replicates, double N, vec time, double prior_rate) {
  
  string line;
  StringVector lines_sync_file;
  ifstream myfile(path.c_str());

  // Tests file open
  if(!myfile) {
    throw range_error("Error opening output file");
  }
  
  // Reads the file
  while (std::getline(myfile, line)) {
    lines_sync_file.push_back(line);
  }
  
  // quantities to keep
  CharacterVector chromossome(number_sites);
  IntegerVector position(number_sites); 
  StringVector reference_allele(number_sites);
  StringVector allele_1(number_sites);
  StringVector allele_2(number_sites);
  NumericVector posterior_mean(number_sites);
  NumericVector posterior_median(number_sites);
  NumericVector posterior_q05(number_sites);
  NumericVector posterior_q95(number_sites);
  NumericVector bayes_factor(number_sites);
  
  for (int i=1; i<(number_sites+1); i++){
    
    cout << "Site number " << i << "\n";
    
    StringVector vector_sync_file = my_split(string(lines_sync_file(i)), '\t');
    
    // Extracts information from chromossome, position and counts
    chromossome(i-1) = string(vector_sync_file(0));
    position(i-1) = stoi(string(vector_sync_file(1))) ; 
    reference_allele(i-1) = string(vector_sync_file(2));
    
    // creates counts matrix
    int total_count_entries = number_time_points*number_replicates;
    mat counts_matrix(4,total_count_entries);
    counts_matrix.zeros();
    
    for (int j = 0; j<total_count_entries; j++){
      StringVector counts = my_split(string(vector_sync_file(j+3)), ':');
      
      counts_matrix(0,j) = stod(string(counts(0)));
      counts_matrix(1,j) = stod(string(counts(1)));
      counts_matrix(2,j) = stod(string(counts(2)));
      counts_matrix(3,j) = stod(string(counts(3)));
      
    }
    
    vec sum_row = sum(counts_matrix,1);
    
    // tests if the site is biallelic
    uvec test_biallelic = find(sum_row == 0);
    if (test_biallelic.n_elem < 2) {
      cout << "Site is triallelic, but the algorithm will procceed by consedering the two alleles with the highest counts." << "\n";
    }
    
    uvec allele1 = find(sort_index(sum_row)==3);
    uvec allele2 = find(sort_index(sum_row)==2);
    
    mat allele_trajectories((N+1)*number_replicates,number_time_points);
    allele_trajectories.zeros();
    
    vec moran_distribution(N+1); 
    moran_distribution.zeros();
    
    int index = 0; 
  
    for (int t=0; t<number_time_points; t++){
      for (int r=0; r<number_replicates; r++){
        double allele_coverage = counts_matrix(allele1(0),index);
        double total_coverage  = counts_matrix(allele1(0),index)+counts_matrix(allele2(0),index);
        moran_distribution = moran_states_distribution(N, allele_coverage, total_coverage);
        allele_trajectories(span(r*N+r,(r+1)*N+r),t) = moran_distribution/sum(moran_distribution);
        index = index+1;
      }
    }

    // usefull quantities
    vec pomo_states = regspace(0,1,N);
    StringVector bases(4);
    bases(0) = "A";
    bases(1) = "T";
    bases(2) = "C";
    bases(3) = "G";
      
    // calculates increments
    mat increments; 
    increments.zeros(number_replicates,number_time_points-1);
    
    double n1,n2,t1,t2;
    
    for (int i=0;i<number_replicates; i++){
      for (int j=1;j<number_time_points;j++){
        
        n1 = sum(pomo_states%allele_trajectories(span(i*N+i,(i+1)*N+i),j-1));
        n2 = sum(pomo_states%allele_trajectories(span(i*N+i,(i+1)*N+i),j));
        t1 = time(j-1);
        t2 = time(j);
        
        increments(i,j-1) = (n2/n1-1)/(t2-t1);
        
      }
    }
    
    // calculates the empirical average and sd of sigma
    double m_sigma = mean(mean(increments));
    double sd_sigma = sqrt(accu(square(m_sigma-increments))/(N-1));
    
    if (m_sigma < 0) {
      
      allele_trajectories = 1 - allele_trajectories;
      
      increments.zeros(number_replicates,number_time_points-1);
      
      double n1,n2,t1,t2;
      
      for (int i=0;i<number_replicates; i++){
        for (int j=1;j<number_time_points;j++){
          
          n1 = sum(pomo_states%allele_trajectories(span(i*N+i,(i+1)*N+i),j-1));
          n2 = sum(pomo_states%allele_trajectories(span(i*N+i,(i+1)*N+i),j));
          t1 = time(j-1);
          t2 = time(j);
          
          increments(i,j-1) = (n2/n1-1)/(t2-t1);
          
        }
      }
      
      // calculates the empirical average and sd of sigma
      double m_sigma = mean(mean(increments));
      double sd_sigma = sqrt(accu(square(m_sigma-increments))/(N-1));
      
    }
    
    // creates the grid
    // eventually put as input
    vec grid;
    grid << -7 << -4 << -2 << -1 << -0.5 <<  0 << 0.5 << 1 << 2 << 4 << 7;
    grid = m_sigma + grid*sd_sigma;
    
    uvec m_one = find(grid < -1);
    grid(m_one) = randu(m_one.n_elem)*(m_sigma+1)-1;
    
    vec log_likelihood(grid.n_elem);
    for (int i =0; i < grid.n_elem; i++){
      log_likelihood(i) = log_posterior(N,grid(i),allele_trajectories,number_replicates,number_time_points,time,prior_rate);
    }
    
    // polynomial regression, degree 9
    vec coefficients = polyfit(grid,log_likelihood,9);
  
    vec x = linspace(min(grid),max(grid),100);
    vec y = coefficients(9)+ coefficients(8)*x + coefficients(7)*pow(x,2)+coefficients(6)*pow(x,3)+coefficients(5)*pow(x,4)+coefficients(4)*pow(x,5)+coefficients(3)*pow(x,6)+coefficients(2)*pow(x,7)+coefficients(1)*pow(x,8)+coefficients(0)*pow(x,9);
    
    // mean
    double mean = sum(x%exp(y)/sum(exp(y)));
    
    // calculates the bazes factor against neutrality
    vec test_sigma;
    test_sigma << 0 << mean;
    vec test_posterior = coefficients(0)+ coefficients(1)*test_sigma  + coefficients(2)*pow(test_sigma ,1)+coefficients(3)*pow(test_sigma ,2)+coefficients(4)*pow(test_sigma ,3)+coefficients(5)*pow(test_sigma ,4)+coefficients(6)*pow(test_sigma ,5)+coefficients(7)*pow(test_sigma ,6)+coefficients(8)*pow(test_sigma ,7)+coefficients(1)*pow(test_sigma ,8)+coefficients(9)*pow(test_sigma ,9);
    
    double c_bayes_factor = 2*(test_posterior(1) - test_posterior(0));
    
    // quantile computation
    vec exp_y = exp(y)/sum(exp(y));
    double median =  quantile_interpolation(x,exp_y,0.5);
    double q_05   =  quantile_interpolation(x,exp_y,0.05);
    double q_95   =  quantile_interpolation(x,exp_y,0.95);
    
    posterior_mean(i-1) = mean;
    posterior_median(i-1) = median;
    posterior_q05(i-1) = q_05;
    posterior_q95(i-1) = q_95;
    bayes_factor(i-1) = c_bayes_factor;

  }
  
  return DataFrame::create(_["chromosome"]= chromossome, _["position"]= position,
                           _["reference_allele"]= reference_allele, _["posterior_mean"]= posterior_mean,
                           _["posterior_median"]= posterior_median, _["posterior_q05"]= posterior_q05,
                           _["posterior_q95"]= posterior_q95, _["bayes_factor"]= bayes_factor);

}


