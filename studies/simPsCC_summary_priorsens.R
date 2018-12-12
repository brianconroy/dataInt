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
outputs <- load_sim_outputs_priorsens()


##################
# Alpha iteration
##################

# alpha case traceplots
param <- "alpha_case"
sampling <- "medium"
prevalence <- "medium"
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
fname <- paste("priorsens_traces_", param, ".png", sep="")
png(paste("/Users/brianconroy/Documents/research/project1/simulations_prior/", fname, sep=""),
    width=900, height=700, res=100)
par(mfrow=c(2,3))
for (i in 1:5){
 
  o <- get_output_priorsens(outputs, 'alpha', i)
  padded_plot(o$samples.alpha.ca, true_params$Alpha.case, '')
  
}
dev.off()


# alpha control traceplots
param <- "alpha_ctrl"
sampling <- "medium"
prevalence <- "medium"
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
fname <- paste("priorsens_traces_", param, ".png", sep="")
png(paste("/Users/brianconroy/Documents/research/project1/simulations_prior/", fname, sep=""),
    width=900, height=700, res=100)
par(mfrow=c(2,3))
for (i in 1:5){
  
  o <- get_output_priorsens(outputs, 'alpha', i)
  padded_plot(o$samples.alpha.co, true_params$Alpha.ctrl, '')
  
}
dev.off()


# log odds scatter plots
s <- "medium"
p <- "medium"
par(mfrow=c(2,3))
xl <- c(-22, 7)
yl <- c(-22, 7)
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
lodds_true_general <- calc_log_odds_true(true_params)
for (i in 1:5){
  o <- get_output_priorsens(outputs, 'alpha', i)
  lodds <- calc_log_odds_output(o, true_params)
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl); abline(0, 1, col=2)
}


##################
# Phi iteration
##################

# Phi traceplots
param <- "phi"
sampling <- "medium"
prevalence <- "medium"
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
fname <- paste("priorsens_traces_", param, ".png", sep="")
png(paste("/Users/brianconroy/Documents/research/project1/simulations_prior/", fname, sep=""),
    width=900, height=700, res=100)
par(mfrow=c(2,3))
for (i in 1:5){
  
  o <- get_output_priorsens(outputs, 'phi', i)
  padded_plot(o$samples.phi, true_params$Phi, '')
  
}
dev.off()


# log odds scatter plots
s <- "medium"
p <- "medium"
par(mfrow=c(2,3))
xl <- c(-22, 7)
yl <- c(-22, 7)
lodds_true <- calc_log_odds_true(s, p)
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
for (i in 1:5){
  o <- get_output_priorsens(outputs, 'phi', i)
  lodds <- calc_log_odds_output(o, true_params)
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl); abline(0, 1, col=2)
}


##################
# Theta iteration
##################

# Theta traceplots
param <- "theta"
sampling <- "medium"
prevalence <- "medium"
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
fname <- paste("priorsens_traces_", param, ".png", sep="")
png(paste("/Users/brianconroy/Documents/research/project1/simulations_prior/", fname, sep=""),
    width=900, height=700, res=100)
par(mfrow=c(2,3))
yl <- c(0, 22)
for (i in 1:5){
  
  o <- get_output_priorsens(outputs, 'phi', i)
  plot(o$samples.theta, ylim=yl, type='l'); abline(h=true_params$Theta, col='2')
  
}
dev.off()


# log odds scatter plots
s <- "medium"
p <- "medium"
par(mfrow=c(2,3))
xl <- c(-22, 7)
yl <- c(-22, 7)
lodds_true <- calc_log_odds_true(s, p)
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
for (i in 1:5){
  o <- get_output_priorsens(outputs, 'theta', i)
  lodds <- calc_log_odds_output(o, true_params)
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl); abline(0, 1, col=2)
}