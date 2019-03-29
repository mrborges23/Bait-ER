# Plots allele trajectories
# delta is a graphical parameter


plot_paths <- function(N,paths,times){
  
  number_replicates <- nrow(paths)/(N+1)
  time_points       <- ncol(paths)
  pomo_states       <- 0:N
  
  plot(NA,xlim=c(times[1],times[time_points]),ylim=c(0,N),ylab="Moran states",xlab="Time points",main="Virtual trajectories")
  for (i in 1:number_replicates) {
    
    # Path for a replicate
    path <- paths[((i-1)*N+i):(i*N+i),]
    
    virtual_trajectory <- rep(NA,time_points)
    for (j in 1:time_points){
      
      # Normalizes the conditional probabilities and calculates the respective colour
      vector <- path[,j]
      virtual_trajectory[j] = sum(pomo_states*vector)
      
    }
    #Plots the most probable trajectory for each replicate
    lines(times,virtual_trajectory,col=i)

  }
  legend("topleft", legend=paste0("replicate ",1:number_replicates),
         col=1:number_replicates, lty=1, cex=0.8, bty="n")
}

# Plot the posterior distribution of sigma

plot_posterior_sigma <- function(post_sigma){
  
  poly_coef <- post_sigma$polynomial_coefficients
  mean      <- post_sigma$mean
  sigma_sd  <- (post_sigma$q95-mean)/2
  sigma     <- seq(mean-7*sigma_sd,mean+7*sigma_sd,length.out=200)
  prob_post <- exp(poly_coef[10]+
                     poly_coef[9]*sigma+
                     poly_coef[8]*sigma^2+
                     poly_coef[7]*sigma^3+
                     poly_coef[6]*sigma^4+
                     poly_coef[5]*sigma^5+
                     poly_coef[4]*sigma^6+
                     poly_coef[3]*sigma^7+
                     poly_coef[2]*sigma^8+
                     poly_coef[1]*sigma^9)
  
  # Plotting the posterior disitrbution of sigma
  plot(sigma,prob_post,ylim=c(min(prob_post),max(prob_post)*1.1),type="l",ylab="posterior probabilities",col="gray",main="Posterior distribution of sigma")
  abline(v=0,col="blue",lty=2)
  
  # Calculating some summary statistics
  mean   <- post_sigma$mean
  median <- post_sigma$median
  q05    <- post_sigma$q05
  q95    <- post_sigma$q95
  
  # Plotting the summary statistics
  y_coordinate <- min(prob_post) + (max(prob_post)-min(prob_post))/2
  segments(q05 ,y_coordinate,q95,y_coordinate,col="red")
  
  segments(q05 ,y_coordinate,q95,y_coordinate,col="red")
  segments(q05 ,y_coordinate,q95,y_coordinate,col="red")
  
  
  points(median,y_coordinate,col="red",pch=4)
  abline(v=mean ,lty=2)
  
  legend(x=min(sigma),y=max(prob_post)*1.1, legend="credibility interval 95%",
         col="red", lty=1, cex=0.8,bty="n")
  
  legend(x=min(sigma),y=max(prob_post)*1.05, legend="mean sigma",
         col="black", lty=2, cex=0.8,bty="n")
  
  legend(x=min(sigma),y=max(prob_post)*1, legend="sigma = 0",
         col="blue", lty=2, cex=0.8,bty="n")
  
}




