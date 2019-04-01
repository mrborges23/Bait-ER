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
  
  # Calculating some summary statistics
  mean   <- post_sigma$mean
  median <- post_sigma$median
  q05    <- post_sigma$q05
  q95    <- post_sigma$q95
  
  alpha  <- post_sigma$alpha
  beta   <- post_sigma$beta
  
  x <- seq(qgamma(0.00001,shape=alpha,rate=beta),qgamma(0.99999,shape=alpha,rate=beta),length.out=500)
  y <- dgamma(x,shape=alpha,rate=beta)
  
  # Plotting the posterior disitrbution of sigma
  plot(x-1,y,type="l",ylab="posterior probabilities",col="gray",main="Posterior distribution of sigma",xlab="sigma")
  abline(v=0,col="blue",lty=2)

  # Plotting the summary statistics
  y_coordinate <- min(y) + (max(y)-min(y))/2
  segments(q05 ,y_coordinate,q95,y_coordinate,col="red")
  
  segments(q05 ,y_coordinate,q95,y_coordinate,col="red")
  segments(q05 ,y_coordinate,q95,y_coordinate,col="red")
  
  
  points(median,y_coordinate,col="red",pch=4)
  abline(v=mean ,lty=2)
  
  legend(x=min(sigma),y=max(y)*1.0, legend="credibility interval 95%",
         col="red", lty=1, cex=0.8,bty="n")
  
  legend(x=min(sigma),y=max(y)*0.95, legend="mean sigma",
         col="black", lty=2, cex=0.8,bty="n")
  

  
}




