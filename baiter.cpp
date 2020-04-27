#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <armadillo>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/binomial.hpp>

using namespace std;
using namespace arma;


// global variables

char nuc_bases[] = "ATCGND";

int n_replicates;
int n_time_points;
int n_loci;
int header;  
int order;

double N;

vec time_vector;
vec prior_parameters;

//functions declaration
void   read_sync_file_line          (ifstream &inFile, string &info, Col<int> &allele_counts, Col<int> &total_counts ) ;
vec    trajectory_probability       (mat &prob_matrix, vec &vector1, vec &vector2);
void   moran_states_distribution    (vec &moran_distribution, double &allele_coverage, double &total_coverage) ;
void   counts_to_moran_states       (mat &allele_trajectories, Col<int> &allele_counts, Col<int> &total_counts);

mat    pomo_rate_matrix             (double &sigma);
double log_likelihood_trajectory    (mat &trajectory_matrix, vec &indexes, mat &prob_matrices);
double log_likelihood_trajectories  (double &sigma, mat &trajectories_matrix);
double log_posterior                (double &sigma, mat &trajectories_matrix);
void   sigma_posterior2             (mat &trajectories_matrix, string &info);

// functions for gamma and binomial distributions
double dbinom(double k,double C, double p);



int main (int argc, char *argv[]){

  // some output
  cout << "\n  You are running Bait-ER!\n \n" << "  Please do not forget to cite us: \n \n" << "  If you find bugs, please report them on our GitHub: \n   - github.com/mrborges23/Bait-ER\n \n";   


  ifstream control(argv[1]);

  if (!control.is_open()) {
    cerr << "  Unable to open " << argv[1] << ". \n \n" << endl;
    return 1;
  }

  string name;

  //these are to be inputed by the user
  // name of sync_file
  // order of counts in important: check again
  // n_replicates and loci

  string sync_file;
  string output_file;

  string time_points;
  string prior_vector;

  control >> name >> sync_file;
  control >> name >> header;  
  control >> name >> order;
  control >> name >> n_replicates;
  control >> name >> time_points;
  control >> name >> n_loci;
  control >> name >> N;
  control >> name >> prior_vector;
  control >> name >> output_file;

  control.close();

  // vector time and n_time_points
  time_vector = vec(time_points);
  n_time_points = time_vector.size();

  // vector of prior parameters
  prior_parameters = vec(prior_vector);

  // some info is printed in the terminal
  cout << "  Information received: \n   - sync file:                 " << sync_file <<"\n   - number loci:               " << n_loci << "\n   - number replicates:         " << n_replicates << "\n   - time points:               " << time_points << "\n   - number time points:        " << n_time_points << "\n   - effective population size: " << N << "\n   - prior parameters:          " << prior_vector << "\n\n";

  // output file
  ofstream outFile;
  outFile.open (output_file);
  outFile << "chromosome\tposition\treference\tsigma\tlogBF\talpha\tbeta\n";

  //opens sync file
  //checks whether sync file is open
  ifstream inFile;
  inFile.open(sync_file, ios::in);
  if (! inFile) {
    cerr << "  Bait-ER could not find " << sync_file << ". Make sure you have " << sync_file << " in the working directory." << endl;
    return 1;
  }

  // checks whether the number of replicates and time points respect 
  // the number of columns in the sync file
  string line,value;
  getline(inFile, line);         
  stringstream s;
  s << line;                   
  int n_cols = 0;    
  while(s >> value) { 
    n_cols++; 
  }
  inFile.close();
  if (n_replicates*n_time_points != (n_cols-3) ) {
    cerr << "  Unexpected number of columns in the sync file.\n  Make sure that your sync file has 3+n_replicates*n_time_points columns. \n" << endl;
    return 1;
  }


  //skip the first row if there is an header
  inFile.open(sync_file, ios::in);
  string H[n_replicates*n_time_points+3];
  if (header == 1){
    for(int j = 0; j < (n_replicates*n_time_points+3); j++){
      inFile >> H[j];
    }
  }

  cout << "  Bait-ER has started!\n";

  // some useful variables
  string info1; 
  string info2; 

  Col<int> allele_counts(n_replicates*n_time_points);
  Col<int> total_counts(n_replicates*n_time_points);
  mat allele_trajectories((N+1)*n_replicates,n_time_points);
  //reads sync file line by line
  for(int i = 0; i < n_loci; i++){
    
    // reset the objects
    allele_counts.zeros();
    total_counts.zeros();
    allele_trajectories.zeros();

    // update info1 and counts
    read_sync_file_line(inFile, info1, allele_counts, total_counts ) ;
    
    // update trajectories matrix
    counts_to_moran_states(allele_trajectories,allele_counts,total_counts);
    
    // posterior estimates and statistics
    sigma_posterior2(allele_trajectories,info2);

    outFile << info1 << info2 << "\n"; 
  }

  // close in and output files
  inFile.close();
  outFile.close();




  // final message
  cout << "  Bait-ER has finished! You can now analyse your output file: " << output_file << ".\n\n";

  return 0;
}


/*
PoMo rate matrix

we assume a time continuous process as may be usefull in E&R designs
we model the rates based on a moran dynamics without mutation, 
which can be ignored for few generations, but with selection

N  population size
sigma selection coefficient
*/

mat pomo_rate_matrix( double &sigma) {
  
  // useful quantities
  int number_pomo_states = N+1;
  mat pomo_matrix(number_pomo_states,number_pomo_states);
  pomo_matrix.zeros();
  
  // populates the matrix
  // pomo states: {Na}, {1A,(N-1)a}, ..., {(N-1)A,1a}, {NA}
  for (int i=1; i<(number_pomo_states-1); i++ ){
    
    //allele A increases in frequency (a decreases)
    pomo_matrix(i,i+1) = i*(N-i)*(1.0+sigma)/N;     
    //allele A decreases in frequency (a increases)
    pomo_matrix(i,i-1) = i*(N-i)/N;

    //diagonal
    pomo_matrix(i,i) = -i*(N-i)*(2.0+sigma)/N;
    
  }
  
  return pomo_matrix;
}




/*
Trajectory probability
 
calculates the conditional probability p(vector2|vector1), where vector 1 and 2 are two probability vectors
observed vectors: vector 1 and vector 2
prob_matrix: the transpose of the pomo rate matrix, in the increment of time between t1 and t2
 
 */

vec trajectory_probability(mat &prob_matrix, vec &vector1, vec &vector2) {

  // calculates the conditional probabilities of observing the pomo state vector vector2 in time t2 
  // given that the pomo state vector in t1 was vector1
  vec pomo_vector = prob_matrix*vector1; 
  vec probability_vector = pomo_vector%vector2;
  
  return probability_vector;
}


/*
 Trajectory likelihood

 calculates the likelihood of observing several allele trajectory within the same replicate
 it accounts for noise coming from the coverage readings
 data consists of state vector containing probabilities, each column representing one time point
 and each row represents a pomo state
 
 Parameters:
 N: population size
 sigma: selection coefficient
 trajectory_matrix: congregates the vector of observed states
 times:
 */

double log_likelihood_trajectory(  mat &trajectory_matrix, vec &indexes, mat &prob_matrices) {
  
  // useful quantities
  vec vector1 =trajectory_matrix.col(0);
  
  for (int i=0; i<(n_time_points-1); i++){
    
    // nth observed vector
    vec vector2 = trajectory_matrix.col(i+1);
    
    // gets the transpose prob_matrix
    int index = indexes(i);
    mat prob_matrix = prob_matrices(span(index*N+index,(index+1)*N+index),span(0,N));
    
    // calculates conditional probabilities
    vector1 = trajectory_probability(prob_matrix,vector1,vector2);
  }
  
  // returns the log conditional probability of observing the trajectory given 
  // a particular selection coefficient sigma 
  double log_probability = log(sum(vector1));
  return log_probability;
}

/*
 Calculates log likelihood over several replicates
 
 trajectories_matrix: columns represent each time point
 rows represent the pomo states
 the replicates are row binded
*/

double log_likelihood_trajectories( double &sigma,  mat &trajectories_matrix) {

  //avoids calculating the exp matrix unecessary number of times as it 
  //checks whether the increments between any two time points is equal
  vec indexes(n_time_points-1);
  indexes.zeros();

  mat rate_matrix = pomo_rate_matrix(sigma);

  // calculates increments
  vec increments(n_time_points-1);
  
  for (int i=0; i<(n_time_points-1); i++){
    increments(i) = time_vector(i+1) - time_vector(i);
  }
  
  // evaluates prob_matrix for all the unique increments
  vec unique_increments = unique(increments);
  int number_unique_increments = unique_increments.n_elem;

  //calculates the probability matrix 
  mat prob_matrix_keeper((N+1)*number_unique_increments,N+1);
  for (int i=0; i < number_unique_increments; i++){
      double delta = unique_increments(i);
      prob_matrix_keeper(span(i*N+i,(i+1)*N+i),span(0,N)) = expmat(rate_matrix*delta).t();
      indexes.elem(find(increments == delta)) += i;
  }

  vec log_likelihood_replicates(n_replicates);
  log_likelihood_replicates.zeros();
  
  for (int i=0; i < n_replicates; i++){
    
    // preparing trajectories_matrix, for replicate number i
    mat trajectory_matrix = trajectories_matrix.rows(i*N+i,(i+1)*N+i);

    // calculates the likelihood of each trajectory
    log_likelihood_replicates(i) = log_likelihood_trajectory( trajectory_matrix, indexes, prob_matrix_keeper);
  }

  // total likelihood
  return sum(log_likelihood_replicates);
}




/*
 Calculates de log posterior
 used the log_likelihood and adds the log_prior
 prior of sigma is exponential with rate prior_rate
 prior_rate: the rate of the exponential distribution, prior of sigma
 
*/

double log_posterior(  double &sigma, mat &trajectories_matrix ) {

   // adds the log likelihood to the log prior
   double log_total = log_likelihood_trajectories(sigma , trajectories_matrix)
                      + prior_parameters(0)*log(prior_parameters(1))+(prior_parameters(0)-1.0)*log(1.0+sigma)-prior_parameters(1)*(1.0+sigma)-lgamma(prior_parameters(0));
                      

  return log_total;
}


/*
Calculates a bayes factor against sigma=0
Constitutes a test of neutrality
*/

void sigma_posterior2( mat &trajectories_matrix, string &info) {
  
  // calculates increments
  mat increments; 
  increments.zeros(n_replicates,n_time_points-1);
  vec pomo_states = regspace(0,1,N);
  
  double n1,n2,t1,t2;
  
  for (int i=0;i<n_replicates; i++){
    for (int j=1;j<n_time_points;j++){
      
      n1 = sum(pomo_states%trajectories_matrix(span(i*N+i,(i+1)*N+i),j-1));
      n2 = sum(pomo_states%trajectories_matrix(span(i*N+i,(i+1)*N+i),j));
      t1 = time_vector(j-1);
      t2 = time_vector(j);
      
      increments(i,j-1) = (n2/n1-1)/(t2-t1);
      
    }
  }
  
  // calculates the empirical average and sd of sigma
  double m_sigma  = mean(mean(increments));
  double sd_sigma = mean(mean(abs(increments-m_sigma)));
  
  //cout << m_sigma << " " << sd_sigma << " "; 

  // sometime_vectors de emprirical variance gets too small and numerical problems 
  // on fiting alpha and beta arise (e.g. alpha is negative or beta is null) 
  /*
  if (sd_sigma < 0.001) {
    sd_sigma = 0.001;
  }
  // code ready to exclud
  if (sd_sigma < 0.001) {
  	info = "\tNA\tNA\tNA\tNA"; 
    return;
  }
  */




  // fit the gamma distribution
  vec grid;
  grid << m_sigma-2*sd_sigma << m_sigma-sd_sigma << m_sigma << m_sigma+sd_sigma << m_sigma+2*sd_sigma ;

  // calculate the log_likelihood

  vec log_likelihood(grid.n_elem);
  for (int i =0; i < grid.n_elem; i++){
    log_likelihood(i) = log_posterior(grid(i),trajectories_matrix);
  }
  
  /*
  for (int i =0; i < grid.n_elem; i++){
    cout << grid(i) << " "; 
  }
  for (int i =0; i < grid.n_elem; i++){
    cout << log_likelihood(i) << " "; 
  }
  */

  vec xi = grid + 1.0;
  vec yi = log_likelihood;
  
  double si1 = arma::mean(xi);
  double si2 = arma::mean(yi);
  double si3 = arma::mean(xi%yi);
  double si4 = arma::mean(log(xi));
  double si5 = arma::mean(xi%log(xi));
  double si6 = arma::mean(yi%log(xi));
  double si7 = arma::mean(log(xi)%log(xi));
  double si8 = arma::mean(xi%xi);

  // least square estimates of alpha and beta
  double d      =  (si7*si1*si1-2*si4*si5*si1+si5*si5+si4*si4*si8-si7*si8);
  double alpha  = -((-(-si2-si4)*si4-si6-si7)*(si1*si1-si8)-(-si3-si1*(-si2-si4)-si5)*(si1*si4-si5))/d;
  double beta   = -(si3*si4*si4-si2*si5*si4-si1*si6*si4+si5*si6+si1*si2*si7-si3*si7)/d;

  // another step to avoid numerical problems
  /*
  if (alpha < 0) {

    sd_sigma = 0.01;

    grid << m_sigma-2*sd_sigma << m_sigma-sd_sigma << m_sigma << m_sigma+sd_sigma << m_sigma+2*sd_sigma ;

    // calculate the log_likelihood
    for (int i =0; i < grid.n_elem; i++){
      log_likelihood(i) = log_posterior(grid(i),trajectories_matrix);
    }

    xi = grid + 1.0;
    yi = log_likelihood;
  
    si1 = arma::mean(xi);
    si2 = arma::mean(yi);
    si3 = arma::mean(xi%yi);
    si4 = arma::mean(log(xi));
    si5 = arma::mean(xi%log(xi));
    si6 = arma::mean(yi%log(xi));
    si7 = arma::mean(log(xi)%log(xi));
    si8 = arma::mean(xi%xi);

    // least square estimates of alpha and beta
    d      =  (si7*si1*si1-2*si4*si5*si1+si5*si5+si4*si4*si8-si7*si8);
    alpha  = -((-(-si2-si4)*si4-si6-si7)*(si1*si1-si8)-(-si3-si1*(-si2-si4)-si5)*(si1*si4-si5))/d;
    beta   = -(si3*si4*si4-si2*si5*si4-si1*si6*si4+si5*si6+si1*si2*si7-si3*si7)/d;

  }
  */

  // posterior statistics
  double average           = alpha/beta-1;
  double log_bayes_factor  = log(boost::math::gamma_q(alpha,beta)) - log(boost::math::gamma_p(alpha,beta)); 
  //+                          log(boost::math::gamma_p(prior_parameters(0),prior_parameters(1))) - log(boost::math::gamma_q(prior_parameters(0),prior_parameters(1)));

  //cout << alpha << " " << beta << " " << average << " " << log_bayes_factor << " " << endl;
  // exports posterior as a string
  info = "\t" + to_string(average) + "\t" + to_string(log_bayes_factor) + "\t" + to_string(alpha) + "\t" + to_string(beta); 

}

/*

 Calculates the probability of the moran states given an observed count 
 p({nA,(N-N)a}!c)=Binomial(n/N,C) for given set of time points, replicates 
 C
 the vectors allele_counts and total_counts should be organized first per generation 
 and then per replicate. Example for 2 generations and 3 replicates:
 F0R1 F0R2 F0R3 F1R1 F1R2 F1R3
 
*/


void counts_to_moran_states(mat &allele_trajectories, Col<int> &allele_counts, Col<int> &total_counts) {

  // usefull quantitie
  vec moran_distribution(N+1); 
  
  // calculates the probability of each moran state given the observed counts  
  if (order == 0) {
  int index = 0; 
    for (int t=0; t<n_time_points; t++){
      for (int r=0; r<n_replicates; r++){
      
        double allele_coverage = allele_counts(index);
        double total_coverage  = total_counts(index);

        moran_states_distribution(moran_distribution, allele_coverage, total_coverage);
        allele_trajectories(span(r*N+r,(r+1)*N+r),t) = moran_distribution;
        index = index+1;

      }
    }
  }

  if (order == 1) {
  int index = 0; 
    for (int r=0; r<n_replicates; r++){
      for (int t=0; t<n_time_points; t++){
      
        double allele_coverage = allele_counts(index);
        double total_coverage  = total_counts(index);
      
        moran_states_distribution(moran_distribution,allele_coverage, total_coverage);
        allele_trajectories(span(r*N+r,(r+1)*N+r),t) = moran_distribution;
        index = index+1;

      }
    }
  }
  
}

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

void moran_states_distribution(vec &moran_distribution,  double &allele_coverage, double &total_coverage) {
  
  // populates the moran distribution using the inverse binomial sampling
  for (int i=0; i<(N+1) ;i++){
    moran_distribution(i) = dbinom(allele_coverage,total_coverage,i/N);
  }
  
  moran_distribution = moran_distribution/sum(moran_distribution);

}

/*
  Reads a line of the sync file
*/

void read_sync_file_line(ifstream &inFile, string &info, Col<int> &allele_counts, Col<int> &total_counts ) {

  string chromosome;
  string position;
  string counts_str;

  inFile >> chromosome;
  inFile >> position;
  inFile >> counts_str;

  // goes through all the cols and saves content to C[j]
  Mat<int> counts_matrix(n_replicates*n_time_points,4);

  for(int j = 0; j < (n_replicates*n_time_points); j++){

    // counts_str = "c1:c2:c3:C4:c5:c6"
    // count_str  = "c1"
    inFile >> counts_str;
    stringstream reader(counts_str);
    string count_str;

    for (int k = 0; k < 4; k++){
      getline(reader, count_str, ':');
      counts_matrix(j,k) = stoi(count_str);
    }

  }

  // sums the counts per allele
  // sorts the counts
  Row<int> sum_row = sum(counts_matrix,0);

  uvec sort_sum_row = sort_index(sum_row);
    
  // if the site is biallelic (most of them)
  string s1,s2;
  s1 = nuc_bases[sort_sum_row(3)];     
  s2 = nuc_bases[sort_sum_row(2)];

  // if the site is triallelic
  if (sum_row(sort_sum_row(1)) >1){
    s2 = s2 + nuc_bases[sort_sum_row(1)];
  }

  // if the site is triallelic
  if (sum_row(sort_sum_row(0)) >1){
    s2 = s2 + nuc_bases[sort_sum_row(0)];
  }


  //cout << chromosome << " " << position << " " << nuc_bases[a1] << nuc_bases[a2] << " ";

  info = chromosome + "\t" + position + "\t" + s1 + ":" + s2;   
  // fills the total coverage and the allele coverage matrices
  allele_counts = counts_matrix.col(sort_sum_row(3));
  total_counts  = sum(counts_matrix,1);

}





/*
 Some usefull distributions
*/

double dbinom(double k, double C,  double p){
  double probability = boost::math::binomial_coefficient<double>(C,k)*pow(p,k)*pow(1-p,C-k);
  return probability;
}

