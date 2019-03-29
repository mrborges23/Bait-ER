rm(list=ls())

library("Rcpp")
library("RcppArmadillo")
library("RcppParallel")

sourceCpp("read_sync_file.cpp")
sourceCpp("allele_trajectory_simulator.cpp")
sourceCpp("reads_to_moran_states.cpp")
sourceCpp("sigma_posterior.cpp")
sourceCpp("sites_sigma_posterior.cpp")
sourceCpp("sites_sigma_posterior2.cpp")

source("plots.R")


# SIMULATED DATA

# populatio parameters
Ne                  <- 300
sigma               <- 10/300
initial_frequency   <- 0.5

# experimental design parameters
times               <- c(0,10,20,30,40)
number_replicates   <- 5
coverage            <- 40

# initial state
prob_vector_initial     <- rep(0,Ne+1)
prob_vector_initial[floor(initial_frequency*Ne)+1] <- 1

reads               <- allele_trajectory_simulator(Ne,sigma,times,prob_vector_initial,number_replicates,coverage)

# Estimates sigma from simulated data
N                   <- 300
trajectories_matrix <- reads_to_moran_states(N,reads) 
number_time_points  <- length(times)
prior_rate          <- 0.0001
post_sigma          <- sigma_posterior(N,times,number_time_points,number_replicates,trajectories_matrix,prior_rate)

# plots virtual trajectories and the posterior of sigma (along with some summary statistics)
par(mfrow=c(1,2))
plot_paths(N,trajectories_matrix,times)
plot_posterior_sigma(post_sigma)


library(poolSeq)
simTraj <- reads$allele_coverage/reads$total_coverage
out     <- estimateSH(simTraj, Ne, t=times, haploid=TRUE,h=0.5, simulate.p.value=TRUE)
out$p.value
out$s


# SYNC_FILE: YEAST DATA

# iformation regarding the data set experimental design
file               <- "sync_file_chr1.txt"
time               <- c(0,180,360,540)
number_sites       <- 100
number_replicates  <- 12
number_time_points <- 4

# estimated sigma for several sites
N          <- 700
prior_rate <- 0.00001
sync_file  <- read_sync_file(file,number_sites,number_time_points,number_replicates)
output     <- sites_sigma_posterior2(sync_file,1,10,"output_chr1.txt",N,time,number_time_points,number_replicates,prior_rate)
