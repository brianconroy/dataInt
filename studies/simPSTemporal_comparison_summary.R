############################
# Fit the spatiotemporal
# shared latent process model
############################


library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


dst <- "/Users/brianconroy/Documents/research/project3/simulation_comparison/"
outputs <- load_sim_outputs(tag="simPsTemporal")


#### Prism Principal Components
years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
caPr_all <- list()
caPr.disc_all <- list()
for (h in years){
  caPr_y <- load_prism_pcs_time(h)
  caPr.disc_y <- aggregate(caPr_y, fact=7) 
  caPr_all <- c(caPr_all, caPr_y)
  caPr.disc_all <- c(caPr.disc_all, caPr.disc_y)
}
plot(caPr_all[[1]])
plot(caPr.disc_all[[1]])


#### Plot pcs 1 and 2, separately, rescaled
pcs1 <- list()
pcs2 <- list()
for (i in 1:length(years)){
  pcs1[[i]] <- caPr_all[[i]][[1]]
  pcs2[[i]] <- caPr_all[[i]][[2]]
}
pcs1_re <- equalize_scales4(pcs1)
pcs2_re <- equalize_scales4(pcs2)

#### Figure: pcs 1
par(mfrow=c(3,3))
for (i in 1:length(years)){
  plot(pcs1_re[[i]], main=years[i])
}

#### Figure: pcs 1
par(mfrow=c(3,3))
for (i in 1:length(years)){
  plot(pcs2_re[[i]], main=years[i])
}


#### Table: summarize levels of U
levels <- c("increasing", "decreasing", "alternating")
params_u <- list()
counter <- 1
for (l in levels){
  sim_name <- paste("simPsTemporal", l, sep="_")
  params <- load_output(paste(sim_name, '_params.json', sep=''))
  params_u[[counter]] <- list(Pattern=l, U=paste("(", paste(params$U, collapse=", "), ")", sep=""))
  counter <- counter + 1
}
params_u <- ldply(params_u, 'data.frame')
params_u$U <- as.character(params_u$U)
params_u$Pattern <- as.character(params_u$Pattern)
write_latex_table(params_u, "latex_params_u.txt", path=dst)


#### Table: summarize other simulation parameters used
levels <- c("increasing", "decreasing", "alternating")
params_other <- list()
counter <- 1
for (l in levels){
  sim_name <- paste("simPsTemporal", l, sep="_")
  params <- load_output(paste(sim_name, '_params.json', sep=''))
  params_other[[counter]] <- list(
    Pattern=l, 
    Beta.case=paste("(", paste(params$beta.case, collapse=", "), ")", sep=""),
    Beta.ctrl=paste("(", paste(params$beta.ctrl, collapse=", "), ")", sep=""),
    alpha.case=params$Alpha.case,
    alpha.ctrl=params$Alpha.ctrl
  )
  counter <- counter + 1
}
write_latex_table(ldply(params_other, 'data.frame'), "latex_params_other.txt", path=dst)


#### Table: summarize number of observation sites per year
levels <- c("increasing", "decreasing", "alternating")
counts_loc <- list()
counter <- 1
for (l in levels){
  sim_name <- paste("simPsTemporal", l, sep="_")
  data_ <- load_output(paste(sim_name, '_data.json', sep=''))
  data <- reformat_saved_data(data_)
  row_l <- list(Pattern=l)
  for (i in 1:length(data$locs)){
    row_l[as.character(years[i])] <- sum(data$locs[[i]]$status)
  }
  counts_loc[[counter]] <- row_l
  counter <- counter + 1
}
write_latex_table(ldply(counts_loc, 'data.frame'), "latex_counts_loc.txt", path=dst)


#### Table: summarize counts of observed specimen and prevalences per year
levels <- c("increasing", "decreasing", "alternating")
counts_prev <- list()
counter <- 1
for (l in levels){
  sim_name <- paste("simPsTemporal", l, sep="_")
  data_ <- load_output(paste(sim_name, '_data.json', sep=''))
  data <- reformat_saved_data(data_)
  row_l <- list(Pattern=l)
  for (i in 1:length(data$locs)){
    ncases <- sum(data$case.data[[i]]$y)
    nctrls <- sum(data$ctrl.data[[i]]$y)
    prev <- round(ncases/(ncases+nctrls), 2)
    # row_l[as.character(years[i])] <- paste(prev, ' (', ncases + nctrls, ')', sep='')
    row_l[as.character(years[i])] <- prev
  }
  counts_prev[[counter]] <- row_l
  counter <- counter + 1
}
write_latex_table(ldply(counts_prev, 'data.frame'), "latex_counts_prev.txt", path=dst)


##################
# RMSEs in
# log disease odds
##################


rmses_all <- list()
counter <- 1
levels <- c("increasing", "decreasing", "alternating")
for (level in levels){
  
  sim_name <- paste("simPsTemporal", level, sep="_")
  for (o in outputs){
    if (o$description == sim_name){
      o_tmp <- o
    }
  }
  o_agg <- get_output_general(outputs, tag=paste(level, "_aggregated_ps", sep=""))
  
  
  #### Load simulation parameters
  params <- load_output(paste(sim_name, '_params.json', sep=''))
  alpha.case <- params$Alpha.case
  alpha.ctrl <- params$Alpha.ctrl
  beta.case <- params$beta.case
  beta.ctrl <- params$beta.ctrl
  W <- params$W
  U <- params$U
  Theta <- params$Theta
  Phi <- params$Phi
  
  cells.all <- c(1:ncell(caPr.disc_all[[1]]))[!is.na(values(caPr.disc_all[[1]][[1]]))]
  coords <- xyFromCell(caPr.disc_all[[1]], cell=cells.all)
  D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
  
  
  #### Load Data
  data_ <- load_output(paste(sim_name, '_data.json', sep=''))
  data <- reformat_saved_data(data_)
  data_pooled <- pool_temporal_data(data)
  
  
  #### Calculate errors in each time window for each model
  lodds_est_pooled <- calc_est_lodds_pool_ps(o_agg, data_pooled)
  rmse_temporal <- c()
  rmse_aggregate <- c()
  rmse_pooled <- c()
  rmse_poisson <- c()
  for (y_idx in 1:length(years)){
    year <- years[y_idx]
    
    # calc true log odds
    lodds_true <- calc_true_lodds_time(params, data, year, y_idx)
    
    # calc estimated log odds
    lodds_est_tmp <- calc_est_lodds_time(o_tmp, data, year, y_idx)
    lodds_est_agg <- calc_est_lodds_time_ps(o_agg, data, year, y_idx)
    lodds_est_pois <- calc_est_lodds_time_poisson(data, year, y_idx)
    
    rmse_tmp <- sqrt(mean((lodds_true - lodds_est_tmp)^2))
    rmse_agg <- sqrt(mean((lodds_true - lodds_est_agg)^2))
    rmse_p <- sqrt(mean((lodds_true - lodds_est_pooled)^2))
    rmse_pois <- sqrt(mean((lodds_true - lodds_est_pois)^2))
    
    rmse_temporal <- c(rmse_temporal, rmse_tmp)
    rmse_aggregate <- c(rmse_aggregate, rmse_agg)
    rmse_pooled <- c(rmse_pooled, rmse_p)
    rmse_poisson <- c(rmse_poisson, rmse_pois)
  }
  
  #### store rmses in table
  rmses_all[[counter]] <- make_rmse_row(rmse_poisson, level, years, "poisson")
  rmses_all[[counter+1]] <- make_rmse_row(rmse_aggregate, level, years, "ps")
  rmses_all[[counter+2]] <- make_rmse_row(rmse_pooled, level, years, "ps pooled")
  rmses_all[[counter+3]] <- make_rmse_row(rmse_temporal, level, years, "ps temporal")
  counter <- counter + 4
  
  par(mfrow=c(2,2))
  ul <- 2.75
  barplot(rmse_poisson, ylim=c(0, ul), main='A)', ylab='RMSE', xlab='Year Index', names.arg=years)
  barplot(rmse_aggregate, ylim=c(0,ul), main='B)', ylab='RMSE', xlab='Year Index', names.arg=years)
  barplot(rmse_pooled, ylim=c(0,ul), main='C)', ylab='RMSE', xlab='Year Index', names.arg=years)
  barplot(rmse_temporal, ylim=c(0,ul), main='D)', ylab='RMSE', xlab='Year Index', names.arg=years)
}
write_latex_table(ldply(rmses_all, 'data.frame'), "latex_rmses_all.txt", path=dst)


#########################
# RMSEs in
# random effects
# and temporal parameters
#########################


#### Figure: w scatterplots
par(mfrow=c(1,3))
levels <- c("increasing", "decreasing", "alternating")
labels <- c("A)", "B)", "C)")
w_rmses <- c()
counter <- 1
for (level in levels){
  
  sim_name <- paste("simPsTemporal", level, sep="_")
  for (o in outputs){
    if (o$description == sim_name){
      o_tmp <- o
    }
  }
  
  params <- load_output(paste(sim_name, '_params.json', sep=''))
  W <- params$W
  w.hat <- colMeans(o_tmp$samples.w)
  w_rmses <- c(w_rmses, sqrt(mean((W-w.hat)^2)))
  
  plot(x=W, 
    y=w.hat, 
    main=labels[counter], 
    xlab='Estimated W', 
    ylab='True W'); abline(0,1,col=2)
  counter <- counter + 1
  
}


#### Figure: w + u scatterplots
par(mfrow=c(3,7))
levels <- c("increasing", "decreasing", "alternating")
wu_rmses <- list()
counter <- 1
for (level in levels){
  
  sim_name <- paste("simPsTemporal", level, sep="_")
  for (o in outputs){
    if (o$description == sim_name){
      o_tmp <- o
    }
  }
  
  params <- load_output(paste(sim_name, '_params.json', sep=''))
  W <- params$W
  U <- params$U
  w.hat <- colMeans(o_tmp$samples.w)
  u.hat <- colMeans(o_tmp$samples.u)
  
  wu_rmses_y <- list(Pattern=level)
  for (i in 1:length(years)){
    wu_rmses_y[[as.character(years[i])]] <- round(sqrt(mean((W + U[i] - w.hat - u.hat[i])^2)), 3)
    plot(y=W + U[i],
         x=w.hat + u.hat[i],
         ylab='True Value',
         xlab='Estimated Value',
         main=years[i]); abline(0,1,col=2)
  }
  wu_rmses[[counter]] <- wu_rmses_y
  
  counter <- counter + 1
  
}
#### Table: w + u RMSEs
write_latex_table(ldply(wu_rmses, 'data.frame'), "latex_rmses_w_u.txt", path=dst)


#### Table: other parameter estimates
levels <- c("increasing", "decreasing", "alternating")
param_ests <- list()
counter <- 1
for (level in levels){
  
  sim_name <- paste("simPsTemporal", level, sep="_")
  for (o in outputs){
    if (o$description == sim_name){
      o_tmp <- o
    }
  }
  
  params <- load_output(paste(sim_name, '_params.json', sep=''))
  Theta <- params$Theta
  Phi <- params$Phi
  alpha.case <- params$Alpha.case
  alpha.ctrl <- params$Alpha.ctrl
  theta.hat <- mean(o_tmp$samples.theta)
  phi.hat <- mean(o_tmp$samples.phi)
  alpha.ca.hat <- mean(o_tmp$samples.alpha.ca)
  alpha.co.hat <- mean(o_tmp$samples.alpha.co)
    
  param_ests[[counter]] <- list(
    Trend=level,
    theta.hat=round(theta.hat, 3),
    theta.bias=round(theta.hat-Theta, 3),
    phi.hat=round(phi.hat, 3),
    phi.bias=round(phi.hat-Phi, 3),
    alpha.ca.hat=round(alpha.ca.hat, 3),
    alpha.ca.bias=round(alpha.ca.hat-alpha.case, 3),
    alpha.co.hat=round(alpha.co.hat, 3),
    alpha.co.bias=round(alpha.co.hat-alpha.ctrl, 3)
  )
  counter <- counter + 1
  
}
write_latex_table(ldply(param_ests, 'data.frame'), "latex_param_ests.txt", path=dst)
