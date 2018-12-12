############################
# Summarizes the preferential
# sampling x prevalence
# simulation results
############################

library(plyr)
library(grid)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
outputs <- load_sim_outputs(tag='size')
sizes <- c(75, 123, 619, 1689, 5495)


# log odds scatterplots
par(mfrow=c(3,2))
xl <- c(-15, 5)
yl <- c(-15, 5)
rmses <- c()
for (s in sizes){
  
  o <- get_output_general(outputs, tag=paste('size_', s, sep=""))
  true_params <- load_params(paste('true_params_size_', s, '.json', sep=''))
  lodds <- calc_log_odds_output(o, true_params)
  lodds_true <- calc_log_odds_true_general(true_params)
  rmse <- sqrt(mean((lodds-lodds_true)^2))
  rmses <- c(rmses, rmse)
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl); abline(0, 1, col=2)

}


# log odds rmses plot
par(mfrow=c(1,1))
plot(x=sizes, y=rmses, type='l', ylim=c(0, 1.5), col='2', xlab='Number of collected specimen')
points(x=sizes, y=rmses, col='2')


# log odds rmses table
rows <- list()
counter <- 1
for ( i in 1:length(sizes)){
  s <- sizes[i]
  rmse <- round(rmses[i], 3)
  new_row <- list(
    size=s,
    rmse=rmse
  )
  rows[[counter]] <- new_row
  counter <- counter + 1
}
df_rmse <- ldply(rows, 'data.frame')
write_latex_table(df_rmse, "latex_size_rmse.txt", 
                  path="/Users/brianconroy/Documents/research/project1/simulations_size/")


# traceplots
s <- 75
for (s in sizes){
  
  true_params <- load_params(paste('true_params_size_', s, '.json', sep=''))
  fname <- paste("size_traces_", s, ".png", sep="")
  png(paste("/Users/brianconroy/Documents/research/project1/simulations_size/", fname, sep=""),
      width=900, height=700, res=100)
  par(mfrow=c(2,3))
  o <- get_output_general(outputs, tag=paste('size_', s, sep=""))
  plot_traces_general(o, true_params)
  dev.off()
  
}


# w scatterplots
par(mfrow=c(3,2))
counter <- 1
for (s in sizes){
  
  output <- get_output_general(outputs, tag=paste('size_', s, sep=""))
  truevals <- load_params(paste('true_params_size_', s, '.json', sep=''))
  w.hat <- colMeans(output$samples.w)
  plot(x=truevals$W, y=w.hat, xlab="True W", ylab="Estimated W"); abline(0, 1, col='2')
  
}
