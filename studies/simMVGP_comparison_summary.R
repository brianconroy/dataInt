###############################
# Summarizes simulation results
# of the multispecies comparison
# study
###############################


library(plyr)
library(grid)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
outputs <- load_sim_outputs(tag='simMVGP')
dst <- "/Users/brianconroy/Documents/research/project2/simulations_comparison/"
species <- c(1, 2)
levels <- c("none", "medium", "high")


#### summarize parameters used to generate data
param_rows <- list()
counter <- 1
for (l in levels){
  params <- load_output(paste('simMVGP_comparison_params_', l, '.json', sep=''))
  param_rows[[counter]] <- list(Correlation=l, 
                          Parameter='$\\alpha_+$', 
                          Value_Species1=as.character(params$Alpha.cases[1]), 
                          Value_Species2=as.character(params$Alpha.cases[2]))
  param_rows[[counter+1]] <- list(Correlation=l, 
                          Parameter='$\\alpha_-$', 
                          Value_Species1=as.character(params$Alpha.ctrls[1]), 
                          Value_Species2=as.character(params$Alpha.ctrls[2]))
  param_rows[[counter+2]] <- list(Correlation=l, 
                          Parameter='$\\beta_+$', 
                          Value_Species1=paste('(', paste(params$beta.cases[1,], collapse=', '), ')', sep=''), 
                          Value_Species2=paste('(', paste(params$beta.cases[2,], collapse=', '), ')', sep=''))
  param_rows[[counter+3]] <- list(Correlation=l, 
                          Parameter='$\\beta_-$', 
                          Value_Species1=paste('(', paste(params$beta.ctrls[1,], collapse=', '), ')', sep=''), 
                          Value_Species2=paste('(', paste(params$beta.ctrls[2,], collapse=', '), ')', sep=''))
  counter <- counter + 4
  
}
write_latex_table(ldply(param_rows, 'data.frame'), 'latex_sim_params.txt', dst)


#### summarize simulated case/control counts
rows_count <- list()
for (l in levels){
  
  data <- load_output(paste('simMVGP_comparison_data_', l, '.json', sep=''))
  locs1 <- list(
    status=data$locs$status[[1]],
    cells=data$locs$cells[[1]],
    coords=data$locs$coords[[1]],
    ids=data$locs$ids[[1]]
  )
  locs2 <- list(
    status=data$locs$status[[2]],
    cells=data$locs$cells[[2]],
    coords=data$locs$coords[[2]],
    ids=data$locs$ids[[2]]
  )
  case.data1 <- list(
    y=data$case.data$y[[1]],
    x.standardised=data$case.data$x.standardised[[1]],
    x=data$case.data$x[[1]],
    p=data$case.data$p[[1]]
  )
  case.data2 <- list(
    y=data$case.data$y[[2]],
    x.standardised=data$case.data$x.standardised[[2]],
    x=data$case.data$x[[2]],
    p=data$case.data$p[[2]]
  )
  ctrl.data1 <- list(
    y=data$ctrl.data$y[[1]],
    x.standardised=data$ctrl.data$x.standardised[[1]],
    x=data$ctrl.data$x[[1]],
    p=data$ctrl.data$p[[1]]
  )
  ctrl.data2 <- list(
    y=data$ctrl.data$y[[2]],
    x.standardised=data$ctrl.data$x.standardised[[2]],
    x=data$ctrl.data$x[[2]],
    p=data$ctrl.data$p[[2]]
  )
  data <- list(
    locs=list(locs1, locs2),
    case.data=list(case.data1, case.data2),
    ctrl.data=list(ctrl.data1, ctrl.data2)
  )
  
  rows_count[[counter]] <- list(
    Correlation=l,
    Species=1,
    Cases=sum(case.data1$y),
    Controls=sum(ctrl.data1$y),
    Prevalence=round(sum(case.data1$y)/sum(case.data1$y + ctrl.data1$y), 3)
  )
  rows_count[[counter+1]] <- list(
    Correlation=l,
    Species=2,
    Cases=sum(case.data2$y),
    Controls=sum(ctrl.data2$y),
    Prevalence=round(sum(case.data2$y)/sum(case.data2$y + ctrl.data2$y), 3)
  )
  counter <- counter + 2
}
write_latex_table(ldply(rows_count, 'data.frame'),'latex_sim_counts.txt', dst)


#### log odds scatterplots
par(mfrow=c(2,3))
xl <- c(-15, 7)
yl <- c(-15, 7)
rmses <- list()
counter <- 1
for (l in levels){
  params <- load_output(paste('simMVGP_comparison_params_', l, '.json', sep=''))
  data <- load_output(paste('simMVGP_comparison_data_', l, '.json', sep=''))
  for (s in species){
    # MVGP
    o <- get_output_general(outputs, tag=paste(l, 'mvgp', sep="_"))
    lodds <- calc_lodds_mvgp(o, data, s)
    lodds_true <- calc_lodds_true_multi(params, data, s)
    rmse <- round(sqrt(mean((lodds-lodds_true)^2)), 3)
    plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col=2)
    rmses[[counter]] <- list(Correlation=l, Species=s, Model='MVGP', rmse=rmse)
    counter <- counter + 1
    
    # Separaete species Models
    o <- get_output_general(outputs, tag=paste(l, '_species', s, sep=""))
    lodds_separate <- calc_log_odds_species(o, data, s)
    rmse <- round(sqrt(mean((lodds_separate-lodds_true)^2)), 3)
    plot(x=lodds_true, y=lodds_separate, xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col=2)
    rmses[[counter]] <- list(Correlation=l, Species=s, Model='Separate', rmse=rmse)
    counter <- counter + 1
    
    # Pooled model
    o <- get_output_general(outputs, tag=paste(l, '_pooled', sep=""))
    lodds_pooled <- calc_log_odds_species(o, data, s)
    rmse <- round(sqrt(mean((lodds_pooled-lodds_true)^2)), 3)
    plot(x=lodds_true, y=lodds_pooled, xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col=2)
    rmses[[counter]] <- list(Correlation=l, Species=s, Model='Pooled', rmse=rmse)
    counter <- counter + 1
  }
}
write_latex_table(ldply(rmses, 'data.frame'), "latex_simMVGP_rmses.txt", path=dst)


##### Posterior variance comparison


summarize_post_var <- function(post_var, species, model, level){
  
  s1 <- as.numeric(summary(post_var))
  return(
    list(
      Correlation=level,
      Model=model,
      Species=species,
      Mean=s1[4],
      Median=s1[3]
    )
  )
  
}


plot_post_var <- function(post_var, r, main=''){
  
  r[][!is.na(r[])] <- post_var
  plot(r, main=main)
  
}


#### Summarize and plot posterior variances
postvar_summaries <- list()
lpostvar_summaries <- list()
counter <- 1
for (l in levels){
  
  data <- load_output(paste('simMVGP_comparison_data_', l, '.json', sep=''))
  
  o_mvgp <- get_output_general(outputs, tag=paste(l, 'mvgp', sep="_"))
  o_rodent_sep <- get_output_general(outputs, tag=paste(l, '_species', 1, sep=""))
  o_coyote_sep <- get_output_general(outputs, tag=paste(l, '_species', 2, sep=""))
  o_pooled <- get_output_general(outputs, tag=paste(l, '_pooled', sep=""))
  
  X_rodent <- load_x_standard(as.logical(data$locs$status[[1]]), agg_factor=8)
  X_coyote <- load_x_standard(as.logical(data$locs$status[[2]]), agg_factor=8)
  risk_rodent_sep <- calc_posterior_risk(o_rodent_sep, X_rodent)
  risk_coyote_sep <- calc_posterior_risk(o_coyote_sep, X_coyote)
  risk_rodent_mvgp <- calc_posterior_risk_multi(o_mvgp, X_rodent, species=1)
  risk_coyote_mvgp <- calc_posterior_risk_multi(o_mvgp, X_coyote, species=2)
  risk_rodent_pooled <- calc_posterior_risk(o_pooled, X_rodent)
  risk_coyote_pooled <- calc_posterior_risk(o_pooled, X_coyote)
  
  lodds_rodent_sep <- calc_posterior_lodds(o_rodent_sep, X_rodent)
  lodds_coyote_sep <- calc_posterior_lodds(o_coyote_sep, X_coyote)
  lodds_rodent_pooled <- calc_posterior_lodds(o_pooled, X_rodent)
  lodds_coyote_pooled <- calc_posterior_lodds(o_pooled, X_coyote)
  lodds_rodent_mvgp <- calc_posterior_lodds_multi(o_mvgp, X_rodent, species=1)
  lodds_coyote_mvgp <- calc_posterior_lodds_multi(o_mvgp, X_coyote, species=2)
  
  lpostvar_rodent_sep <- apply(lodds_rodent_sep, 2, var)
  lpostvar_coyote_sep <- apply(lodds_coyote_sep, 2, var)
  lpostvar_rodent_pooled <- apply(lodds_rodent_pooled, 2, var)
  lpostvar_coyote_pooled <- apply(lodds_coyote_pooled, 2, var)
  lpostvar_rodent_mvgp <- apply(lodds_rodent_mvgp, 2, var)
  lpostvar_coyote_mvgp <- apply(lodds_coyote_mvgp, 2, var)
  
  postvar_rodent_sep <- apply(risk_rodent_sep, 2, var)
  postvar_coyote_sep <- apply(risk_coyote_sep, 2, var)
  postvar_rodent_pooled <- apply(risk_rodent_pooled, 2, var)
  postvar_coyote_pooled <- apply(risk_coyote_pooled, 2, var)
  postvar_rodent_mvgp <- apply(risk_rodent_mvgp, 2, var)
  postvar_coyote_mvgp <- apply(risk_coyote_mvgp, 2, var)
  
  postvar_summaries[[counter]] <- summarize_post_var(postvar_rodent_sep, 'rodent', 'separate', l)
  postvar_summaries[[counter+1]] <- summarize_post_var(postvar_rodent_mvgp, 'rodent', 'mvgp', l)
  postvar_summaries[[counter+2]] <- summarize_post_var(postvar_coyote_sep, 'coyote', 'separate', l)
  postvar_summaries[[counter+3]] <- summarize_post_var(postvar_coyote_mvgp, 'coyote', 'mvgp', l)
  
  lpostvar_summaries[[counter]] <- summarize_post_var(lpostvar_rodent_sep, 'rodent', 'separate', l)
  lpostvar_summaries[[counter+1]] <- summarize_post_var(lpostvar_rodent_mvgp, 'rodent', 'mvgp', l)
  lpostvar_summaries[[counter+2]] <- summarize_post_var(lpostvar_coyote_sep, 'coyote', 'separate', l)
  lpostvar_summaries[[counter+3]] <- summarize_post_var(lpostvar_coyote_mvgp, 'coyote', 'mvgp', l)
  
  counter <- counter + 4
  
  par(mfrow=c(2,3))
  plot_post_var(postvar_rodent_mvgp, caPr.disc[[1]], main='A)')
  plot_post_var(postvar_rodent_sep, caPr.disc[[1]], main='B)')
  plot_post_var(postvar_rodent_pooled, caPr.disc[[1]], main='C)')
  plot_post_var(postvar_coyote_mvgp, caPr.disc[[1]], main='D)')
  plot_post_var(postvar_coyote_sep, caPr.disc[[1]], main='E)')
  plot_post_var(postvar_coyote_pooled, caPr.disc[[1]], main='F)')

}
ldply(postvar_summaries, 'data.frame')
ldply(lpostvar_summaries, 'data.frame')


#### Calculate parameter estimates/bias/posterior variance
for (l in c('none', 'medium', 'high')){
  
  params <- load_output(paste('simMVGP_comparison_params_', l, '.json', sep=''))
  data <- load_output(paste('simMVGP_comparison_data_', l, '.json', sep=''))
  
  
  for (s in species){
    
    estimates <- list()
    
    # MVGP
    o <- get_output_general(outputs, tag=paste(l, 'mvgp', sep="_"))
    estimates <- c(estimates, summarize_multi_params(o, params, s))
    
    # Separaete species models
    o_sep <- get_output_general(outputs, tag=paste(l, '_species', s, sep=""))
    estimates <- c(estimates, summarize_params(o_sep, params, s, model='Separate'))
    
    # Pooled model
    o_pooled <- get_output_general(outputs, tag=paste(l, '_pooled', sep=""))
    estimates <- c(estimates, summarize_params(o_pooled, params, s, model='Pooled'))
    
    estimates_df <- ldply(estimates, data.frame)
    estimates_df <- estimates_df[with(estimates_df, order(Species, Parameter, Model)),]
    estimates_df$Parameter <- as.character(estimates_df$Parameter)
    estimates_df <- replace_vals(estimates_df, column='Parameter', val='Beta 0 (case)', '$\\beta_{0, +}$')
    estimates_df <- replace_vals(estimates_df, column='Parameter', val='Beta 1 (case)', '$\\beta_{1, +}$')
    estimates_df <- replace_vals(estimates_df, column='Parameter', val='Beta 2 (case)', '$\\beta_{2, +}$')
    
    estimates_df <- replace_vals(estimates_df, column='Parameter', val='Beta 0 (control)', '$\\beta_{0, -}$')
    estimates_df <- replace_vals(estimates_df, column='Parameter', val='Beta 1 (control)', '$\\beta_{1, -}$')
    estimates_df <- replace_vals(estimates_df, column='Parameter', val='Beta 2 (control)', '$\\beta_{2, -}$')
    
    estimates_df <- replace_vals(estimates_df, column='Parameter', val='Alpha (case)', '$\\alpha_+$')
    estimates_df <- replace_vals(estimates_df, column='Parameter', val='Alpha (control)', '$\\alpha_-$')
    
    write_latex_table(estimates_df, paste("latex_simMVGP_estimates_", l, '_', s, ".txt", sep=""), path=dst)
    
  }
  
}
