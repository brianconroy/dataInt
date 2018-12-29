###############################
# Summarizes simulation results
# of the normal vs. truncated
# normal prior comparison
###############################

library(plyr)
library(grid)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
outputs <- load_sim_outputs(tag='priorcompare')
priors <- c("normal_uncal", "normal_cal", "tnormal_uncal", "tnormal_cal")
dst <- "/Users/brianconroy/Documents/research/project1/simulations_priorcompare/"


# log odds scatterplots
par(mfrow=c(2,2))
xl <- c(-22, 10)
yl <- c(-22, 10)
rmses <- c()
true_params <- load_params('true_params_priorcompare.json')
for (p in priors){
  
  o <- get_output_general(outputs, tag=paste('priorcompare_', p, sep=""))
  lodds <- calc_log_odds_output(o, true_params)
  lodds_true <- calc_log_odds_true_general(true_params)
  rmse <- sqrt(mean((lodds-lodds_true)^2))
  rmses <- c(rmses, rmse)
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl); abline(0, 1, col=2)

}


# log odds rmses table
rows <- list()
counter <- 1
for ( i in 1:length(priors)){
  p <- priors[i]
  rmse <- round(rmses[i], 3)
  new_row <- list(
    configuration=p,
    rmse=rmse
  )
  rows[[counter]] <- new_row
  counter <- counter + 1
}
df_rmse <- ldply(rows, 'data.frame')
write_latex_table(df_rmse, "latex_priorcompare_rmse.txt", path=dst)


# traceplots
true_params <- load_params('true_params_priorcompare.json')
for (p in priors){
  
  fname <- paste("priorcompare_traces_", p, ".png", sep="")
  png(paste(dst, fname, sep=""),
      width=900, height=700, res=100)
  par(mfrow=c(2,3))
  o <- get_output_general(outputs, tag=paste('priorcompare_', p, sep=""))
  plot_traces_general(o, true_params)
  dev.off()
  
}


# w scatterplots
par(mfrow=c(2,2))
counter <- 1
truevals <- load_params('true_params_priorcompare.json')
for (p in priors){
  
  o <- get_output_general(outputs, tag=paste('priorcompare_', p, sep=""))
  w.hat <- colMeans(o$samples.w)
  plot(x=truevals$W, y=w.hat, xlab="True W", ylab="Estimated W"); abline(0, 1, col='2')
  
}


# export tuning parameters to LaTeX
rows_all <- c()
for (p in priors){
  
  o <- get_output_general(outputs, tag=paste('priorcompare_', p, sep=""))
  rows_o <- summarize_mcmc_pscc(o, p)
  rows_all <- c(rows_o, rows_all)
  
}
df_mcmc <- ldply(rows_all, 'data.frame')
write_latex_table(df_mcmc, "latex_priorcompare_mcmc.txt", path=dst)
