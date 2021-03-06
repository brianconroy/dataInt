###############################
# Summarizes simulation results
# of the multispecies iteration
# study (modeling option 1)
###############################

library(plyr)
library(grid)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')

caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
outputs <- load_sim_outputs(tag='simMulti_opt1_iteration')
dst <- "/Users/brianconroy/Documents/research/project2/simulation_1_iteration/"
species_ns <- c(2, 4, 6, 8)

###########
# log odds 
#   scatter
#   box
###########

data_store <- load_output('simMulti_data_iteration.json')
true_params <- data_store
true_params$data <- list(
  locs=data_store$locs
)

rmses <- list()
counter <- 1
for (n_s in species_ns){
  o <- get_output_general(outputs, tag=paste('simMulti_opt1_iteration_', n_s, sep=""))
  for (species in 1:n_s){
    lodds <- calc_log_odds_multi(o, true_params, species, output_type='_multi')
    lodds_true <- calc_log_odds_true_multi_general(true_params, species)
    rmse <- round(sqrt(mean((lodds-lodds_true)^2)), 3)
    plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col=2)
    rmses[[counter]] <- list(species=species, total_species=n_s, rmse=rmse)
    counter <- counter + 1
  }
}
rmse_df <- ldply(rmses, 'data.frame')
boxplot(rmse ~ total_species, data=rmse_df, xlab='Total Species', ylab='RMSE')

# summary error table
rows <- list()
counter <- 1
for (n_s in species_ns){
  rmse_s <- rmse_df[rmse_df$total_species == n_s,]$rmse
  rows[[counter]] <- list(
    min=min(rmse_s),
    q1=quantile(rmse_s, p=0.25),
    mean=round(mean(rmse_s), 3),
    median=round(mean(rmse_s), 3),
    q3=quantile(rmse_s, p=0.75),
    max=max(rmse_s)
  )
  counter <- counter + 1
}
write_latex_table(ldply(rows, 'data.frame'), "latex_simMulti_iteration_rmses.txt", path=dst)

# scatterplots for just species 1
#   illustrative purposes
labs <- c("A)", "B)", "C)", "D)")
par(mfrow=c(2,2))
species <- 1
counter <- 1
for (n_s in species_ns){
  lab <- labs[counter]
  o <- get_output_general(outputs, tag=paste('simMulti_opt1_iteration_', n_s, sep=""))
  lodds <- calc_log_odds_multi(o, true_params, species, output_type='_multi')
  lodds_true <- calc_log_odds_true_multi_general(true_params, species)
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', main=lab); abline(0, 1, col=2)
  counter <- counter + 1
}

###
# W
###
labs <- c("A)", "B)", "C)", "D)")
counter <- 1
W <- true_params$W
for (n_s in species_ns){
  lab <- labs[counter]
  o <- get_output_general(outputs, tag=paste('simMulti_opt1_iteration_', n_s, sep=""))
  w.hat <- colMeans(o$samples.w)
  plot(x=W, y=w.hat, xlab='True W', ylab='Estimated W', main=lab); abline(0, 1, col=2)
}

############
# alpha bias 
# estimates
############

alpha_est <- list()
counter <- 1
for (n_s in species_ns){
  o <- get_output_general(outputs, tag=paste('simMulti_opt1_iteration_', n_s, sep=""))
  for (species in 1:n_s){
    alpha_est[[counter]] <- list(
      species=species,
      total_species=n_s,
      parameter="alpha (case)",
      estimate=round(mean(o$samples.alpha.ca[species,,]), 3),
      true_value=true_params$alpha.cases[species,],
      bias=round(mean(o$samples.alpha.ca[species,,]), 3)-true_params$alpha.cases[species,]
    )
    counter <- counter + 1
  }
}
for (n_s in species_ns){
  o <- get_output_general(outputs, tag=paste('simMulti_opt1_iteration_', n_s, sep=""))
  for (species in 1:n_s){
    alpha_est[[counter]] <- list(
      species=species,
      total_species=n_s,
      parameter="alpha (control)",
      estimate=round(mean(o$samples.alpha.co[species,,]), 3),
      true_value=true_params$alpha.ctrls[species,],
      bias=round(mean(o$samples.alpha.co[species,,]), 3)-true_params$alpha.ctrls[species,]
    )
    counter <- counter + 1
  }
}
alpha_df <- ldply(alpha_est, 'data.frame')
boxplot(bias ~ total_species, data=alpha_df[alpha_df$parameter=='alpha (case)',],
        xlab='Total Species', ylab='Bias')
boxplot(bias ~ total_species, data=alpha_df[alpha_df$parameter=='alpha (control)',],
        xlab='Total Species', ylab='Bias')
