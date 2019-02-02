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
dst <- "/Users/brianconroy/Documents/research/project1/simulations_prior/"

##################
# Alpha iteration
##################

# alpha case traceplots
param <- "alpha_case"
sampling <- "medium"
prevalence <- "medium"
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
fname <- paste("priorsens_traces_", param, ".png", sep="")
png(paste(dst, fname, sep=""), width=900, height=700, res=100)
par(mfrow=c(2,3))
for (i in 1:5){
 
  o <- get_output_priorsens(outputs, 'alpha', i)
  plot(o$samples.alpha.ca, type='l', ylim=c(0, 1.5)); abline(h=true_params$Alpha.case, col=2)
  
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
  plot(o$samples.alpha.co, type='l', ylim=c(-1.5, 0)); abline(h=true_params$Alpha.ctrl, col=2)
  
}
dev.off()


# export tuning parameters to LaTeX
rows_all <- c()
for (i in 1:5){
  
  o <- get_output_priorsens(outputs, 'alpha', i)
  rows_o <- summarize_mcmc_pscc(o, paste('alpha prior', i))
  rows_all <- c(rows_o, rows_all)
  
}
df_mcmc <- ldply(rows_all, 'data.frame')
write_latex_table(df_mcmc, "latex_priorsens_mcmc_alpha.txt", path=dst)


# export summary table of parameter biases
sampling <- "medium"
prevalence <- "medium"
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
param <- "alpha"
rows_alpha <- list()
counter <- 1
for (i in 1:5){
  
  priors <- load_priors(param, i)
  o <- get_output_priorsens(outputs, param, i)
  rows_alpha[[counter]] <- list(
    parameter='alpha (case)',
    prior_variance=priors$prior_alpha,
    estimate=round(mean(o$samples.alpha.ca), 3),
    posterior_var=round(var(o$samples.alpha.ca), 3),
    bias=round(mean(o$samples.alpha.ca), 3)-true_params$Alpha.case
  )
  counter <- counter + 1

}
for (i in 1:5){
  
  priors <- load_priors(param, i)
  o <- get_output_priorsens(outputs, param, i)
  rows_alpha[[counter]] <- list(
    parameter='alpha (control)',
    prior_variance=priors$prior_alpha,
    estimate=round(mean(o$samples.alpha.co), 3),
    posterior_var=round(var(o$samples.alpha.co), 3),
    bias=round(mean(o$samples.alpha.co), 3)-true_params$Alpha.ctrl
  )
  counter <- counter + 1
  
}
alpha_df <- ldply(rows_alpha, 'data.frame')
write_latex_table(alpha_df, "latex_priorsens_alpha_bias.txt", path=dst)

# log odds scatter plots
s <- "medium"
p <- "medium"
par(mfrow=c(2,3))
xl <- c(-22, 7)
yl <- c(-22, 7)
param <- "alpha"
true_params <- load_params(paste('true_params_', s, '_', p, '.json', sep=''))
lodds_true_general <- calc_log_odds_true_general(true_params)
rmses <- list()
for (i in 1:5){
  o <- get_output_priorsens(outputs, 'alpha', i)
  lodds <- calc_log_odds_output(o, true_params)
  plot(x=lodds_true_general, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl); abline(0, 1, col=2)

  priors <- load_priors(param, i)
  rmses[[i]] <- list(
    prior_variance=priors$prior_alpha,
    rmse=round(sqrt(mean((lodds-lodds_true_general)^2)), 3)
  )
  
}
write_latex_table(ldply(rmses, 'data.frame'), "latex_priorsens_alpha_rmse.txt", path=dst)

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
  plot(o$samples.phi, type='l', ylim=c(0,30)); abline(h=true_params$Phi, col=2)
  
}
dev.off()


# export tuning parameters to LaTeX
rows_all <- c()
for (i in 1:5){
  
  o <- get_output_priorsens(outputs, 'phi', i)
  rows_o <- summarize_mcmc_pscc(o, paste('phi prior', i))
  rows_all <- c(rows_o, rows_all)
  
}
df_mcmc <- ldply(rows_all, 'data.frame')
write_latex_table(df_mcmc, "latex_priorsens_mcmc_phi.txt", path=dst)


# export summary table of parameter biases
sampling <- "medium"
prevalence <- "medium"
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
param <- "phi"
rows_phi <- list()
counter <- 1
for (i in 1:5){
  
  priors <- load_priors(param, i)
  o <- get_output_priorsens(outputs, param, i)
  rows_phi[[counter]] <- list(
    Prior_Variance=ig_var(priors$prior_phi[1], priors$prior_phi[2]),
    Estimate=round(mean(o$samples.phi), 3),
    Posterior_Variance=round(var(o$samples.phi), 3),
    Bias=round(mean(o$samples.phi), 3)-true_params$Phi
  )
  counter <- counter + 1
  
}
write_latex_table(ldply(rows_phi, 'data.frame'), "latex_priorsens_phi_bias.txt", path=dst)


# log odds scatter plots
s <- "medium"
p <- "medium"
par(mfrow=c(2,3))
xl <- c(-22, 7)
yl <- c(-22, 7)
rmses <- list()
lodds_true <- calc_log_odds_true(s, p)
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
for (i in 1:5){
  o <- get_output_priorsens(outputs, 'phi', i)
  lodds <- calc_log_odds_output(o, true_params)
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl); abline(0, 1, col=2)

  priors <- load_priors(param, i)
  rmses[[i]] <- list(
    prior_variance=ig_var(priors$prior_phi[1], priors$prior_phi[2]),
    rmse=round(sqrt(mean((lodds-lodds_true_general)^2)), 3)
  )
}
write_latex_table(ldply(rmses, 'data.frame'), "latex_priorsens_phi_rmse.txt", path=dst)


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


# export tuning parameters to LaTeX
rows_all <- c()
for (i in 1:5){
  
  o <- get_output_priorsens(outputs, 'theta', i)
  rows_o <- summarize_mcmc_pscc(o, paste('theta prior', i))
  rows_all <- c(rows_o, rows_all)
  
}
df_mcmc <- ldply(rows_all, 'data.frame')
write_latex_table(df_mcmc, "latex_priorsens_mcmc_theta.txt", path=dst)


# export summary table of parameter biases
sampling <- "medium"
prevalence <- "medium"
true_params <- load_params(paste('true_params_', sampling, '_', prevalence, '.json', sep=''))
param <- "theta"
rows_theta <- list()
counter <- 1
for (i in 1:5){
  
  priors <- load_priors(param, i)
  o <- get_output_priorsens(outputs, param, i)
  rows_theta[[counter]] <- list(
    Prior_Variance=g_var(priors$prior_theta[1], priors$prior_theta[2]),
    Estimate=round(mean(o$samples.theta), 3),
    Posterior_Variance=round(var(o$samples.theta), 3),
    Bias=round(mean(o$samples.theta), 3)-true_params$Theta
  )
  counter <- counter + 1
  
}
write_latex_table(ldply(rows_theta, 'data.frame'), "latex_priorsens_theta_bias.txt", path=dst)



# log odds scatter plots
s <- "medium"
p <- "medium"
par(mfrow=c(2,3))
xl <- c(-22, 7)
yl <- c(-22, 7)
param <- "theta"
rmses <- list()
lodds_true <- calc_log_odds_true(s, p)
true_params <- load_params(paste('true_params_', s, '_', p, '.json', sep=''))
for (i in 1:5){
  o <- get_output_priorsens(outputs, 'theta', i)
  lodds <- calc_log_odds_output(o, true_params)
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl); abline(0, 1, col=2)

  priors <- load_priors(param, i)
  rmses[[i]] <- list(
    prior_variance=g_var(priors$prior_theta[1], priors$prior_theta[2]),
    rmse=round(sqrt(mean((lodds-lodds_true_general)^2)), 3)
  )
}
write_latex_table(ldply(rmses, 'data.frame'), "latex_priorsens_theta_rmse.txt", path=dst)

