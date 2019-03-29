// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <log_posterior.cpp>
#include <quantile_interpolation.cpp>

using namespace arma;
using namespace Rcpp;

/*

Calculates a bayes factor against sigma=0
Constitutes a test of neutrality

*/

// [[Rcpp::export]]
rowvec sigma_posterior2(double N, vec time, int number_time_points, int number_replicates,mat trajectories_matrix, double prior_rate) {
  
  // usefull quantities
  vec pomo_states = regspace(0,1,N);
  rowvec output(15);
  output.zeros();
  
  // calculates increments
  mat increments; 
  increments.zeros(number_replicates,number_time_points-1);
  
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
  double m_sigma = mean(mean(increments));
  double sd_sigma = sqrt(accu(square(m_sigma-increments))/(N-1));

  // creates the grid
  // eventually put as input
  vec grid;
  grid << -5 << -3 << -2 << -1 << -0.2 << 0 << 0.2 << 1 << 2 << 3 << 5;
  grid = m_sigma + grid*sd_sigma;

  vec log_likelihood(grid.n_elem);
  for (int i =0; i < grid.n_elem; i++){
    log_likelihood(i) = log_posterior(N,grid(i),trajectories_matrix,number_replicates,number_time_points,time,prior_rate);
  }

  // if loglik returns nonfinite values
  // find them and recalculate them in a suitable range
  uvec non_finite_ll = find_nonfinite(log_likelihood);
  int number_non_finite_ll = non_finite_ll.n_elem;
  
  if (number_non_finite_ll > 0){
    
    // finds finite elements
    // establishes a range for those values that we were able to assess the log-likelihood
    uvec finite_ll = find_finite(log_likelihood);
    double sigma_min = grid(finite_ll(0));
    double sigma_max = grid(finite_ll(finite_ll.n_elem-1));
    
    // recalculate grid
    grid(non_finite_ll) = randu(number_non_finite_ll)*(sigma_max-sigma_min)+sigma_min;
    
    // recalculate log_lik
    for (int i=0; i<number_non_finite_ll; i++){
      log_likelihood(non_finite_ll(i)) = log_posterior(N,grid(non_finite_ll(i)),trajectories_matrix,number_replicates,number_time_points,time,prior_rate);
    }
  }

  // polynomial regression, degree 9
  // several analysis showes this degree is approapriated
  vec coefficients = polyfit(grid,log_likelihood,9);

  vec x  = linspace(min(grid),max(grid),100);
  vec y  = coefficients(9)+ coefficients(8)*x + coefficients(7)*pow(x,2)+coefficients(6)*pow(x,3)+coefficients(5)*pow(x,4)+coefficients(4)*pow(x,5)+coefficients(3)*pow(x,6)+coefficients(2)*pow(x,7)+coefficients(1)*pow(x,8)+coefficients(0)*pow(x,9);
  vec ny = exp(y)/sum(exp(y));
  
  // mean
  double mean     = sum(x%ny);
  double variance = sum(x%x%ny)-mean*mean;
  
  vec x2           = x+1;
  double mean2     = sum(x2%ny);
  double variance2 = sum(x2%x2%ny)-mean2*mean2;
  
  // calculates the bayes factor against neutrality
  double log_bayes_factor       = R::pgamma(1,mean2*mean2/variance2,variance2/mean2,0,1) - R::pgamma(1,mean2*mean2/variance2,variance2/mean2,1,1);
  
  // quantile computation
  double median =  quantile_interpolation(x,ny,0.5);
  double q_05   =  quantile_interpolation(x,ny,0.05);
  double q_95   =  quantile_interpolation(x,ny,0.95);
  
  // exports posterior statistics as a vector
  output(0) = mean;
  output(1) = median;
  output(2) = q_05;
  output(3) = q_95;
  output(4) = log_bayes_factor;
  output(span(5,14)) = coefficients.t();

  return output;
  
}

