rm(list=ls())

library("Rcpp")
library("RcppArmadillo")
library("RcppParallel")

sourceCpp("Bait-ER Code/read_sync_file.cpp")
sourceCpp("Bait-ER Code/allele_trajectory_simulator.cpp")
sourceCpp("Bait-ER Code/reads_to_moran_states.cpp")
sourceCpp("Bait-ER Code/sigma_posterior.cpp")
sourceCpp("Bait-ER Code/sites_sigma_posterior.cpp")

source("Bait-ER Code/plots.R")


# SIMULATED DATA

# populatio parameters
Ne                  <- 300
sigma               <- 10/300
initial_frequency   <- 0.5

# experimental design parameters
times               <- c(0,10,20,30,40)
number_replicates   <- 10
coverage            <- 200

# initial state
prob_vector_initial     <- rep(0,Ne+1)
prob_vector_initial[floor(initial_frequency*Ne)+1] <- 1

reads               <- allele_trajectory_simulator(Ne,sigma,times,prob_vector_initial,number_replicates,coverage)

# Estimates sigma from simulated data
N                   <- 300
trajectories_matrix <- reads_to_moran_states(N,reads) 
number_time_points  <- length(times)
prior_parameters    <- c(1,rate[i])
post_sigma          <- sigma_posterior(N,times,number_time_points,number_replicates,trajectories_matrix,prior_parameters)

# plots virtual trajectories and the posterior of sigma (along with some summary statistics)
par(mfrow=c(1,2))
plot_paths(N,trajectories_matrix,times)
plot_posterior_sigma(post_sigma)


# SYNC_FILE: YEAST DATA

# information regarding the data set experimental design
file               <- "Data/data_sync.txt"
time               <- c(0,20,40,60)
number_sites       <- 100
number_replicates  <- 12
number_time_points <- 4

# estimated sigma for several sites
N          <- 300
prior_parameters    <- c(1,1)
sync_file  <- read_sync_file(file,number_sites,number_time_points,number_replicates)
output     <- sites_sigma_posterior(sync_file,1,10,"Data/output.txt",N,time,number_time_points,number_replicates,pprior_parameters)

