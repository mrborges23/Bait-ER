// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <split.cpp>

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

*/


// [[Rcpp::export]]
List read_sync_file(string path,int number_sites,int number_time_points, int number_replicates) {
  
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
  
  // useful quantities
  int total_count_entries = number_time_points*number_replicates;
  
  StringVector vector_sync_file;
  StringVector counts(total_count_entries+3);

  Mat<int> counts_matrix(4,total_count_entries);

  Col<int> sum_row;
  uvec test_biallelic;
  uvec sort_sum_row;
  int a1;
  int a2;
  
  // nucleotide bases vector
  StringVector nuc_bases(4);
  nuc_bases(0) = "A"; nuc_bases(1) = "T"; nuc_bases(2) = "C"; nuc_bases(3) = "G";
  
  // quantities to keep
  CharacterVector chromossome(number_sites);
  IntegerVector position(number_sites); 
  StringVector reference_allele(number_sites);
  
  StringVector allele1(number_sites);
  StringVector allele2(number_sites);
  
  Mat<int> counts_allele1(number_sites,total_count_entries);
  Mat<int> total_coverage(number_sites,total_count_entries);
  
  IntegerVector info(number_sites);
  
  // Goes line by line and extract iformation that fills the list of 
  // vectors and matrices in // quantitites to keep
  for (int i=1; i<(number_sites+1); i++){
    
    vector_sync_file = my_split(string(lines_sync_file(i)), '\t');
    
    // Extracts information from chromossome, position and reference allele
    chromossome(i-1) = string(vector_sync_file(0));
    position(i-1) = stoi(string(vector_sync_file(1))) ; 
    reference_allele(i-1) = string(vector_sync_file(2));
    
    // creates counts matrix
    // tranforms the strings countsA:countsT:countsC:countsG in measures of coverage
    counts_matrix.zeros();
    
    for (int j = 0; j<total_count_entries; j++){
      
      counts = my_split(string(vector_sync_file(j+3)), ':');
      
      counts_matrix(0,j) = stoi(string(counts(0)));
      counts_matrix(1,j) = stoi(string(counts(1)));
      counts_matrix(2,j) = stoi(string(counts(2)));
      counts_matrix(3,j) = stoi(string(counts(3)));
      
    }

    sum_row = sum(counts_matrix,1);
    test_biallelic = find(sum_row > 0);

    // tests if the site
    // 0: has no counts for any of the nucleotide bases
    // 1: is monomorphic
    // 2: is biallelic
    // 3: is triallelic
    // this information is returned in vector info
    info(i-1) = test_biallelic.n_elem;
    
    sort_sum_row = sort_index(sum_row);
    
    // Finds the two alleles that have the 
    a1 = sort_sum_row(3);
    a2 = sort_sum_row(2);
    
    allele1(i-1) = nuc_bases(a1);
    allele2(i-1) = nuc_bases(a2);
    
    // fills the total coverage and the allele coverage matrices
    counts_allele1(i-1,span(0,total_count_entries-1)) = counts_matrix(a1,span(0,total_count_entries-1));
    total_coverage(i-1,span(0,total_count_entries-1)) = counts_matrix(a1,span(0,total_count_entries-1))+counts_matrix(a2,span(0,total_count_entries-1));
  }  
  
  // returns a list
  return List::create(_["chromosome"]= chromossome, 
                      _["position"]= position,
                      _["reference_allele"]= reference_allele, 
                      _["allele1"]= allele1,
                      _["allele2"]= allele2, 
                      _["allele1_counts"]= counts_allele1,
                      _["total_coverage"]= total_coverage, 
                      _["info"]= info);
}
    
