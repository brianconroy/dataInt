#########################
# Summarizes the sampling
# size simulations

# key finding: need to set
# good initial values as 
# size > 619
#########################

library(plyr)
library(grid)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
outputs <- load_sim_outputs(tag='size')
sizes <- c(75, 123, '619_fail', '619_success', 1689, 5495)
dst <- "/Users/brianconroy/Documents/research/project1/simulations_size/"


# log odds scatterplots
par(mfrow=c(3,2))
xl <- c(-15, 5)
yl <- c(-15, 5)
rmses <- c()
labels <- c("A)", "B)", "C)", "D)", "E)", "F)")
counter <- 1
for (s in sizes){
  
  o <- get_output_general(outputs, tag=paste('size_', s, sep=""))
  true_params <- load_params(paste('true_params_size_', s, '.json', sep=''))
  lodds <- calc_log_odds_output(o, true_params)
  lodds_true <- calc_log_odds_true_general(true_params)
  rmse <- sqrt(mean((lodds-lodds_true)^2))
  rmses <- c(rmses, rmse)
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl, main=labels[counter]); abline(0, 1, col=2)
  counter <- counter + 1
  
}


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
write_latex_table(df_rmse, "latex_size_rmse.txt", path=dst)


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
labels <- c("A)", "B)", "C)", "D)", "E)", "F)")
for (s in sizes){
  
  output <- get_output_general(outputs, tag=paste('size_', s, sep=""))
  truevals <- load_params(paste('true_params_size_', s, '.json', sep=''))
  w.hat <- colMeans(output$samples.w)
  plot(x=truevals$W, y=w.hat, xlab="True W", ylab="Estimated W", main=labels[counter]); abline(0, 1, col='2')
  counter <- counter + 1
  
}


# export tuning parameters to LaTeX
rows_all <- c()
for (s in sizes){
  
  o <- get_output_general(outputs, tag=paste('size_', s, sep=""))
  rows_o <- summarize_mcmc_pscc(o, s)
  rows_all <- c(rows_o, rows_all)
  
}
df_mcmc <- ldply(rows_all, 'data.frame')
write_latex_table(df_mcmc, "latex_size_mcmc.txt", path=dst)

# save params to table
sizes <- c(75, 123, '619_success', 1689, 5495)
rows <- list()
counter <- 1
for (s in sizes){
  params <- load_sim_params_size(s)
  rows[[counter]] <- list(
    size=params$size,
    prevalence=params$prev,
    beta.case=params$beta.case,
    beta.ctrl=params$beta.ctrl,
    alpha.case=params$alpha.case,
    alpha.ctrl=params$alpha.ctrl
  )
  counter <- counter + 1
}
write_latex_table(ldply(rows, 'data.frame'), "latex_size_params.txt", path=dst)


# table of bias for alpha params
# summary table of bias
rows_alpha <- list()
counter <- 1
for (s in sizes){
  true_params <- load_params(paste('true_params_size_', s, '.json', sep=''))
  o <- get_output_general(outputs, tag=paste('size_', s, sep=""))
  rows_alpha[[counter]] <- list(
    size=s,
    parameter="alpha (case)",
    estimate=round(mean(o$samples.alpha.ca), 3),
    bias=round(mean(o$samples.alpha.ca), 3) - true_params$Alpha.case
  )
  counter <- counter + 1
}
for (s in sizes){
  true_params <- load_params(paste('true_params_size_', s, '.json', sep=''))
  o <- get_output_general(outputs, tag=paste('size_', s, sep=""))
  rows_alpha[[counter]] <- list(
    size=s,
    parameter="alpha (control)",
    estimate=round(mean(o$samples.alpha.co), 3),
    bias=round(mean(o$samples.alpha.co), 3) - true_params$Alpha.ctrl
  )
  counter <- counter + 1
}
write_latex_table(ldply(rows_alpha, 'data.frame'), "latex_alpha_bias.txt", path=dst)


rows_beta <- list()
counter <- 1
for (s in sizes){
  true_params <- load_params(paste('true_params_size_', s, '.json', sep=''))
  o <- get_output_general(outputs, tag=paste('size_', s, sep=""))
  for (i in 1:3){
    rows_beta[[counter]] <- list(
      size=s,
      parameter=paste("beta (case)", i-1),
      estimate=round(mean(o$samples.beta.ca[,i]), 3),
      bias=round(mean(o$samples.beta.ca[,i]), 3) - true_params$beta.case[i]
    )
    counter <- counter + 1
  }
  for (i in 1:3){
    rows_beta[[counter]] <- list(
      size=s,
      parameter=paste("beta (control)", i-1),
      estimate=round(mean(o$samples.beta.co[,i]), 3),
      bias=round(mean(o$samples.beta.co[,i]), 3) - true_params$beta.ctrl[i]
    )
    counter <- counter + 1
  }
  
}
rows_beta <- ldply(rows_beta, 'data.frame')
boxplot(bias ~ size, data=rows_beta, names=c(75, 123, '619 (uncal)', 
                                             '619 (cal)', 1689, 5495), ylab='bias', xlab='size')
