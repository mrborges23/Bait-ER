// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <log_posterior.cpp>

using namespace arma;
using namespace Rcpp;

/*

Calculates a bayes factor of sigma<0 and sigma>0
Constitutes a test of neutrality

*/

// [[Rcpp::export]]
List sigma_posterior(double N, vec time, int number_time_points, int number_replicates,mat trajectories_matrix, vec prior_parameters) {
  
  // calculates increments
  mat increments; 
  increments.zeros(number_replicates,number_time_points-1);
  vec pomo_states = regspace(0,1,N);
  
  double n1,n2,t1,t2;
  
  for (int i=0;i<number_replicates; i++){
    for (int j=1;j<number_time_points;j++){
      
      n1 = sum(pomo_states%trajectories_matrix(span(i*N+i,(i+1)*N+i),j-1));
      n2 = sum(pomo_states%trajectories_matrix(span(i*N+i,(i+1)*N+i),j));
      t1 = time(j-1);
      t2 = time(j);
      
      increments(i,j-1) = (n2/n1-1)/(t2-t1);
      
    }
  }
  
  // calculates the empirical average and sd of sigma
  double m_sigma  = mean(mean(increments));
  double sd_sigma = mean(mean(abs(increments-m_sigma)));
  
  // fit tht gamma distribution
  vec grid;
  grid << m_sigma-2*sd_sigma << m_sigma-sd_sigma << m_sigma << m_sigma+sd_sigma << m_sigma+2*sd_sigma ;
  
  // calculate the log_likelihood
  vec log_likelihood(grid.n_elem);
  for (int i =0; i < grid.n_elem; i++){
    log_likelihood(i) = log_posterior(N,grid(i),trajectories_matrix,number_replicates,number_time_points,time,prior_parameters);
  }
   
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
  double d     =  (si7*si1*si1-2*si4*si5*si1+si5*si5+si4*si4*si8-si7*si8);
  double alpha = -((-(-si2-si4)*si4-si6-si7)*(si1*si1-si8)-(-si3-si1*(-si2-si4)-si5)*(si1*si4-si5))/d;
  double beta  = -(si3*si4*si4-si2*si5*si4-si1*si6*si4+si5*si6+si1*si2*si7-si3*si7)/d;
  
  // posterior statistics
  double mean              = alpha/beta-1;
  double var               = alpha/(beta*beta);
  double median            = R::qgamma(0.50,alpha,1/beta,1,0)-1;
  double q05               = R::qgamma(0.05,alpha,1/beta,1,0)-1;
  double q95               = R::qgamma(0.95,alpha,1/beta,1,0)-1;
  double log_bayes_factor  = R::pgamma(1,alpha,1/beta,0,1) - R::pgamma(1,alpha,1/beta,1,1)- 
                             R::pgamma(1,prior_parameters(0),1/prior_parameters(1),0,1) + R::pgamma(1,prior_parameters(0),1/prior_parameters(1),1,1);
 
  // exports posterior statistics as a list
  return List::create(Named("mean")              = mean,
                      Named("variance")          = var,
                      Named("median")            = median,
                      Named("q05")               = q05,
                      Named("q95")               = q95,
                      Named("log_bayes_factor")  = log_bayes_factor,
                      Named("alpha")             = alpha,
                      Named("beta")              = beta );
  
  
}

