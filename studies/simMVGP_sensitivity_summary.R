###############################
# Summarizes simulation results
# of the multispecies robustness
# study
###############################


library(plyr)
library(grid)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


dst <- "/Users/brianconroy/Documents/research/project2/simulations_sensitivity/"
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)


#### Figure: prism pcs
par(mfrow=c(1,2))
plot(caPr.disc[[1]], main='A)')
plot(caPr.disc[[2]], main='B)')


#### Table: summary of simulated locations
sim_locs <- list()
for (l in c("0", "01", "02", "1", "2", "3")){
  params <- load_output(paste('simMVGP_sensitivity_params_', l, '.json', sep=''))
  data <- load_output(paste('simMVGP_sensitivity_data_', l, '.json', sep=''))
  if (l %in% c("0", "01", "02")){
    phis <- "(1, 1)"
  } else{
    phis <- "(3, 3)"
  }
  l1 <- sum(data$locs$status[[1]])
  l2 <- sum(data$locs$status[[2]])
  
  sim_locs[[counter]] <- list(
    Phis=phis,
    Thetas=paste("(", params$Theta1, ", ", params$Theta2, ")", sep=""),
    species=1,
    Locations=l1
  )
  sim_locs[[counter+1]] <- list(
    Phis=phis,
    Thetas=paste("(", params$Theta1, ", ", params$Theta2, ")", sep=""),
    species=2,
    Locations=l2
  )
  counter <- counter + 2
}
write_latex_table(ldply(sim_locs, 'data.frame'), "latex_locs.txt", path=dst)


#### Table: summarize simulation params
sim_params <- list()
counter <- 1
for (s in c(1, 2)){
  for (l in c("0", "01", "02", "1", "2", "3")){
    params <- load_output(paste('simMVGP_sensitivity_params_', l, '.json', sep=''))
    if (l %in% c("0", "01", "02")){
      phis <- "(1, 1)"
    } else{
      phis <- "(3, 3)"
    }
    sim_params[[counter]] <- list(
      Phis=phis,
      Thetas=paste("(", params$Theta1, ", ", params$Theta2, ")", sep=""),
      species=s,
      Betacase=paste("(", paste(params$beta.cases[s,], collapse=", "),")", sep=""),
      Alphacase=params$Alpha.cases[s,],
      Betactrl=paste("(", paste(params$beta.ctrls[s,], collapse=", "), ")", sep=""),
      Alphactrl=params$Alpha.ctrls[s,]
    )
    counter <- counter + 1
  }
}
write_latex_table(ldply(sim_params, 'data.frame'), "latex_params.txt", path=dst)


#### Table: summarize simulation counts
sim_counts <- list()
counter <- 1
for (l in c("0", "01", "02", "1", "2", "3")){
  params <- load_output(paste('simMVGP_sensitivity_params_', l, '.json', sep=''))
  data <- load_output(paste('simMVGP_sensitivity_data_', l, '.json', sep=''))
  if (l %in% c("0", "01", "02")){
    phis <- "(1, 1)"
  } else{
    phis <- "(3, 3)"
  }
  y1_case <- sum(data$case.data$y[[1]])
  y2_case <- sum(data$case.data$y[[2]])
  y1_ctrl <- sum(data$ctrl.data$y[[1]])
  y2_ctrl <- sum(data$ctrl.data$y[[2]])
  p1 <- y1_case/(y1_case + y1_ctrl)
  p2 <- y2_case/(y2_case + y2_ctrl)
  
  sim_counts[[counter]] <- list(
    Phis=phis,
    Thetas=paste("(", params$Theta1, ", ", params$Theta2, ")", sep=""),
    species=1,
    Case_Count=y1_case,
    Control_Count=y1_ctrl,
    Prevalence=round(p1, 2)
  )
  sim_counts[[counter+1]] <- list(
    Phis=phis,
    Thetas=paste("(", params$Theta1, ", ", params$Theta2, ")", sep=""),
    species=2,
    Case_Count=y2_case,
    Control_Count=y2_ctrl,
    Prevalence=round(p2, 2)
  )
  counter <- counter + 2
}
write_latex_table(ldply(sim_counts, 'data.frame'), "latex_counts.txt", path=dst)


outputs <- load_sim_outputs(tag='simMVGP_sensitivity')
species <- c(1, 2)


####  Figure & Table: rmse
# levels <- c("0", "01", "02")
levels <- c("1", "2", "3")
par(mfrow=c(2,3))
rmses <- list()
counter <- 1
for (s in species){
  for (l in levels){
    params <- load_output(paste('simMVGP_sensitivity_params_', l, '.json', sep=''))
    data <- load_output(paste('simMVGP_sensitivity_data_', l, '.json', sep=''))
    o <- get_output_general(outputs, tag=paste(l, 'mvgp', sep="_"))
    lodds <- calc_lodds_mvgp(o, data, s, agg_factor=8)
    lodds_true <- calc_lodds_true_multi(params, data, s)
    rmse <- round(sqrt(mean((lodds-lodds_true)^2)), 3)
    plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col=2)
    
    if (l %in% c("0", "01", "02")){
      phis <- "(1, 1)"
    } else{
      phis <- "(3, 3)"
    }
    
    rmses[[counter]] <- list(
      Phis=phis,
      Thetas=paste("(", params$Theta1, ", ", params$Theta2, ")", sep=""),
      Species=s, 
      rmse=rmse)
    counter <- counter + 1
  }
}
write_latex_table(ldply(rmses, 'data.frame'), "latex_simMVGP_rmses_01.txt", path=dst)


#### w
# levels <- c("0", "01", "02")
levels <- c("1", "2", "3")
par(mfrow=c(2,3))
rmses_w <- list()
counter <- 1
for (s in species){
  for (l in levels){
    params <- load_output(paste('simMVGP_sensitivity_params_', l, '.json', sep=''))
    data <- load_output(paste('simMVGP_sensitivity_data_', l, '.json', sep=''))
    o <- get_output_general(outputs, tag=paste(l, 'mvgp', sep="_"))
    w.hat <- colMeans(o$samples.w)
    w.hat <- w.hat[seq(s, length(w.hat), by=2)]
    W <- params$W
    W <- W[seq(s, length(W), by=2)]
    plot(x=W, y=w.hat, xlab='True W', ylab='EstimatedW'); abline(0, 1, col=2)
    
    if (l %in% c("0", "01", "02")){
      phis <- "(1, 1)"
    } else{
      phis <- "(3, 3)"
    }
    
    rmse <- round(sqrt(mean((W-w.hat)^2)), 3)
    rmses_w[[counter]] <- list(
      Phis=phis,
      Thetas=paste("(", params$Theta1, ", ", params$Theta2, ")", sep=""),
      Species=s, 
      rmse=rmse)
    counter <- counter + 1
    
  }
}
write_latex_table(ldply(rmses_w, 'data.frame'), "latex_simMVGP_rmses_w_1.txt", path=dst)

