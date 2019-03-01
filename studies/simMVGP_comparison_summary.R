###############################
# Summarizes simulation results
# of the multispecies comparison
# study (modeling option 1)
###############################

library(plyr)
library(grid)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')

# compare on basis of 
#   estimated log disease odds for each species
#   estimated w
#   bias parameter estimates

caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=8)
outputs <- load_sim_outputs(tag='simMVGP')
dst <- "/Users/brianconroy/Documents/research/project2/simulations_comparison/"
species <- c(1, 2)
levels <- c("none", "medium", "high")

#######################
# risk scatterplots
#######################

# 4 models, 3 correlation levels

# row 1: no correlation
#   mvgp1, species1, mvgp2, species2
# ...

# rmse table
  # correlation | species | model | lodds rmse

par(mfrow=c(3,4))
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
    
    # Multispecies Model
    o <- get_output_general(outputs, tag=paste(l, '_species', s, sep=""))
    lodds_separate <- calc_log_odds_species(o, data, s)
    rmse <- round(sqrt(mean((lodds_separate-lodds_true)^2)), 3)
    plot(x=lodds_true, y=lodds_separate, xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col=2)
    rmses[[counter]] <- list(Correlation=l, Species=s, Model='Separate', rmse=rmse)
    counter <- counter + 1
  }
}
write_latex_table(ldply(rmses, 'data.frame'), "latex_simMVGP_rmses.txt", path=dst)



# fix: same W
# fix: species 2


perturbs <- c("low", "high")
models <- c("_multi", "_separate_species2", "_pooled")
species <- 2

par(mfrow=c(2,3))
xl <- c(-20, 7)
yl <- c(-20, 7)
rmses2 <- c()
labs <- c("A)", "B)", "C)", "D)", "E)", "F)")
counter <- 1
for (p in perturbs){
  for (m in models){
    lab <- labs[counter]
    true_params <- load_params(paste("true_params_simMulti_opt1_comparison_", p, ".json", sep=""))
    o <- get_output_general(outputs, tag=paste('simMulti_opt1_comparison_', p, m, sep=""))
    lodds <- calc_log_odds_multi(o, true_params, species, m)
    lodds_true <- calc_log_odds_true_multi(true_params, species)
    rmse <- round(sqrt(mean((lodds-lodds_true)^2)), 3)
    plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl, main=lab); abline(0, 1, col=2)
    rmses2[[counter]] <- list(perturbation=p, model=m, w="same", rmse=rmse)
    counter <- counter + 1
  }
}


# fix: different W
# fix: species 1


perturbs <- c("low_different")
models <- c("_multi", "_separate1", "_pooled")
species <- 1

par(mfrow=c(2,3))
xl <- c(-20, 7)
yl <- c(-20, 7)
rmses3 <- c()
labs <- c("A)", "B)", "C)")
counter <- 1
for (p in perturbs){
  for (m in models){
    lab <- labs[counter]
    true_params <- load_params(paste("true_params_simMulti_opt1_comparison_", p, ".json", sep=""))
    o <- get_output_general(outputs, tag=paste('simMulti_opt1_comparison_', p, m, sep=""))
    lodds <- calc_log_odds_multi(o, true_params, species, m)
    lodds_true <- calc_log_odds_true_multi(true_params, species)
    rmse <- round(sqrt(mean((lodds-lodds_true)^2)), 3)
    plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl, main=lab); abline(0, 1, col=2)
    rmses3[[counter]] <- list(perturbation=p, model=m, w="different", rmse=rmse)
    counter <- counter + 1
  }
}

# fix: different W
# fix: species 2

models <- c("_multi", "_separate_species2", "_pooled")
species <- 2
xl <- c(-20, 7)
yl <- c(-20, 7)
rmses4 <- c()
labs <- c("D)", "E)", "F)")
counter <- 1
for (p in perturbs){
  for (m in models){
    lab <- labs[counter]
    true_params <- load_params(paste("true_params_simMulti_opt1_comparison_", p, ".json", sep=""))
    o <- get_output_general(outputs, tag=paste('simMulti_opt1_comparison_', p, m, sep=""))
    lodds <- calc_log_odds_multi(o, true_params, species, m)
    lodds_true <- calc_log_odds_true_multi(true_params, species)
    rmse <- round(sqrt(mean((lodds-lodds_true)^2)), 3)
    plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl, main=lab); abline(0, 1, col=2)
    rmses4[[counter]] <- list(perturbation=p, model=m, w="different", rmse=rmse)
    counter <- counter + 1
  }
}


rmses <- rbind(
  ldply(rmses1, 'data.frame'),
  ldply(rmses2, 'data.frame'),
  ldply(rmses3, 'data.frame'),
  ldply(rmses4, 'data.frame')
)
write_latex_table(rmses, "latex_simMulti_opt1_rmse.txt", path=dst)


# traceplots
perturbs <- c("low", "high")
models <- c("_multi", "_separate1", "_pooled")
species <- c(1)

for (p in perturbs){
  for (m in models){
    for (s in species){
      true_params <- load_params(paste("true_params_simMulti_opt1_comparison_", p, ".json", sep=""))
      o <- get_output_general(outputs, tag=paste('simMulti_opt1_comparison_', p, m, sep=""))
      
      fname <- paste("multi_opt1_traces_", p, m, "_", s, ".png", sep="")
      png(paste(dst, fname, sep=""),
          width=900, height=700, res=100)
      plot_traces_multi(o, true_params, s, m)
      dev.off()
    }
  }
}

perturbs <- c("low", "high")
models <- c("_multi", "_separate_species2", "_pooled")
species <- c(2)

for (p in perturbs){
  for (m in models){
    for (s in species){
      true_params <- load_params(paste("true_params_simMulti_opt1_comparison_", p, ".json", sep=""))
      o <- get_output_general(outputs, tag=paste('simMulti_opt1_comparison_', p, m, sep=""))
      
      fname <- paste("multi_opt1_traces_", p, m, "_", s, ".png", sep="")
      png(paste(dst, fname, sep=""),
          width=900, height=700, res=100)
      plot_traces_multi(o, true_params, s, m)
      dev.off()
    }
  }
}


perturbs <- c("low_different")
# models <- c("_multi", "_separate1", "_pooled")
# species <- 1

models <- c("_multi", "_separate_species2", "_pooled")
species <- 2

for (p in perturbs){
  for (m in models){
    for (s in species){
      true_params <- load_params(paste("true_params_simMulti_opt1_comparison_", p, ".json", sep=""))
      o <- get_output_general(outputs, tag=paste('simMulti_opt1_comparison_', p, m, sep=""))
      
      fname <- paste("multi_opt1_traces_", p, m, "_", s, ".png", sep="")
      png(paste(dst, fname, sep=""),
          width=900, height=700, res=100)
      plot_traces_multi(o, true_params, s, m)
      dev.off()
    }
  }
}


#######################
# param estimate tables
#######################


make_row <- function(p, w, species, m, parameter, est, true_val){
  
  return(list(
    perturbation=p,
    w=w,
    species=species,
    model=m,
    parameter=parameter,
    estimate=est,
    bias=est-true_val
  ))
  
}


# perturbation | w | species | model | parameter | estimate | bias
perturbs <- c("low", "high", "low_different")
models <- c("_multi", "_separate1", "_pooled")
species <- 1

rows <- list()
dec <- 3
counter <- 1
for (p in perturbs){
  for (m in models){
    true_params <- load_params(paste("true_params_simMulti_opt1_comparison_", p, ".json", sep=""))
    o <- get_output_general(outputs, tag=paste('simMulti_opt1_comparison_', p, m, sep=""))
    if (m == "_multi"){
      bc0 <- round(mean(o$samples.beta.ca[1,,1]), dec)
      bc1 <- round(mean(o$samples.beta.ca[1,,2]), dec)
      bc2 <- round(mean(o$samples.beta.ca[1,,3]), dec)
      bco0 <- round(mean(o$samples.beta.co[1,,1]), dec)
      bco1 <- round(mean(o$samples.beta.co[1,,2]), dec)
      bco2 <- round(mean(o$samples.beta.co[1,,3]), dec)
      alpha.ca <- round(mean(o$samples.alpha.ca[species,,]), dec)
      alpha.co <- round(mean(o$samples.alpha.co[species,,]), dec)
    } else {
      bc0 <- round(mean(o$samples.beta.ca[,1]), dec)
      bc1 <- round(mean(o$samples.beta.ca[,2]), dec)
      bc2 <- round(mean(o$samples.beta.ca[,3]), dec)
      bco0 <- round(mean(o$samples.beta.co[,1]), dec)
      bco1 <- round(mean(o$samples.beta.co[,2]), dec)
      bco2 <- round(mean(o$samples.beta.co[,3]), dec)
      alpha.ca <- round(mean(o$samples.alpha.ca), dec)
      alpha.co <- round(mean(o$samples.alpha.co), dec)
    }
    phi <- round(mean(o$samples.phi), dec)
    theta <- round(mean(o$samples.theta), dec)
    
    bc0_ <- true_params$beta.case1[1]
    bc1_ <- true_params$beta.case1[2]
    bc2_ <- true_params$beta.case1[3]
    bco0_ <- true_params$beta.ctrl1[1]
    bco1_ <- true_params$beta.ctrl1[2]
    bco2_ <- true_params$beta.ctrl1[3]
    alpha.ca_ <- true_params$Alpha.case1
    alpha.co_ <- true_params$Alpha.ctrl1
    phi_ <- true_params$Phi
    theta_ <- true_params$Theta
    
    if (p == 'low_different'){
      w='different'
      p_='low'
    } else {
      w='same'
      p_=p
    }
    
    rows[[counter]] <- make_row(p_, w=w, species, m, parameter="beta0 (case)", bc0, bc0_)
    rows[[counter+1]] <- make_row(p_, w=w, species, m, parameter="beta1 (case)", bc1, bc1_)
    rows[[counter+2]] <- make_row(p_, w=w, species, m, parameter="beta2 (case)", bc2, bc2_)
    rows[[counter+3]] <- make_row(p_, w=w, species, m, parameter="beta2 (control)", bco0, bco0_)
    rows[[counter+4]] <- make_row(p_, w=w, species, m, parameter="beta2 (control)", bco1, bco1_)
    rows[[counter+5]] <- make_row(p_, w=w, species, m, parameter="beta2 (control)", bco2, bco2_)
    rows[[counter+6]] <- make_row(p_, w=w, species, m, parameter="alpha (case)", alpha.ca, alpha.ca_)
    rows[[counter+7]] <- make_row(p_, w=w, species, m, parameter="alpha (control)", alpha.co, alpha.co_)
    rows[[counter+8]] <- make_row(p_, w=w, species, m, parameter="theta", theta, theta_)
    rows[[counter+9]] <- make_row(p_, w=w, species, m, parameter="phi", phi, phi_)
    counter <- counter + 10
    
  }
 
}


models <- c("_multi", "_separate2", "_pooled")
species <- 2
for (p in perturbs){
  for (m in models){
    true_params <- load_params(paste("true_params_simMulti_opt1_comparison_", p, ".json", sep=""))
    o <- get_output_general(outputs, tag=paste('simMulti_opt1_comparison_', p, m, sep=""))
    if (m == "_multi"){
      bc0 <- round(mean(o$samples.beta.ca[1,,1]), dec)
      bc1 <- round(mean(o$samples.beta.ca[1,,2]), dec)
      bc2 <- round(mean(o$samples.beta.ca[1,,3]), dec)
      bco0 <- round(mean(o$samples.beta.co[1,,1]), dec)
      bco1 <- round(mean(o$samples.beta.co[1,,2]), dec)
      bco2 <- round(mean(o$samples.beta.co[1,,3]), dec)
      alpha.ca <- round(mean(o$samples.alpha.ca[species,,]), dec)
      alpha.co <- round(mean(o$samples.alpha.co[species,,]), dec)
    } else {
      bc0 <- round(mean(o$samples.beta.ca[,1]), dec)
      bc1 <- round(mean(o$samples.beta.ca[,2]), dec)
      bc2 <- round(mean(o$samples.beta.ca[,3]), dec)
      bco0 <- round(mean(o$samples.beta.co[,1]), dec)
      bco1 <- round(mean(o$samples.beta.co[,2]), dec)
      bco2 <- round(mean(o$samples.beta.co[,3]), dec)
      alpha.ca <- round(mean(o$samples.alpha.ca), dec)
      alpha.co <- round(mean(o$samples.alpha.co), dec)
    }
    phi <- round(mean(o$samples.phi), dec)
    theta <- round(mean(o$samples.theta), dec)
    
    bc0_ <- true_params$beta.case2[1]
    bc1_ <- true_params$beta.case2[2]
    bc2_ <- true_params$beta.case2[3]
    bco0_ <- true_params$beta.ctrl2[1]
    bco1_ <- true_params$beta.ctrl2[2]
    bco2_ <- true_params$beta.ctrl2[3]
    alpha.ca_ <- true_params$Alpha.case2
    alpha.co_ <- true_params$Alpha.ctrl2
    phi_ <- true_params$Phi
    theta_ <- true_params$Theta
    
    if (p == 'low_different'){
      w='different'
      p_='low'
    } else {
      w='same'
      p_=p
    }
    
    rows[[counter]] <- make_row(p_, w=w, species, m, parameter="beta0 (case)", bc0, bc0_)
    rows[[counter+1]] <- make_row(p_, w=w, species, m, parameter="beta1 (case)", bc1, bc1_)
    rows[[counter+2]] <- make_row(p_, w=w, species, m, parameter="beta2 (case)", bc2, bc2_)
    rows[[counter+3]] <- make_row(p_, w=w, species, m, parameter="beta2 (control)", bco0, bco0_)
    rows[[counter+4]] <- make_row(p_, w=w, species, m, parameter="beta2 (control)", bco1, bco1_)
    rows[[counter+5]] <- make_row(p_, w=w, species, m, parameter="beta2 (control)", bco2, bco2_)
    rows[[counter+6]] <- make_row(p_, w=w, species, m, parameter="alpha (case)", alpha.ca, alpha.ca_)
    rows[[counter+7]] <- make_row(p_, w=w, species, m, parameter="alpha (control)", alpha.co, alpha.co_)
    rows[[counter+8]] <- make_row(p_, w=w, species, m, parameter="theta", theta, theta_)
    rows[[counter+9]] <- make_row(p_, w=w, species, m, parameter="phi", phi, phi_)
    counter <- counter + 10
    
  }
  
}
df <- ldply(rows, 'data.frame')
write_latex_table(df, "latex_multi_params.txt", path=dst)

boxplot(bias ~ model, 
        data=df[! df$parameter %in% c('theta', 'phi') & df$w == 'same',], 
        names=c('Multispecies', 'Separate (1)', 'Pooled', 'Separate (2)'),
        xlab='Model',
        ylab='Bias')
boxplot(bias ~ model, 
        data=df[! df$parameter %in% c('theta', 'phi') & df$w != 'same',],
        names=c('Multispecies', 'Separate (1)', 'Pooled', 'Separate (2)'),
        xlab='Model',
        ylab='Bias')

sums <- c()
for (m in c('_multi', '_separate1', '_separate2', '_pooled')){
  rowsm <- df[df$model == m & !df$parameter %in% c('theta', 'phi') & df$w =='same',]
  sums <- rbind(sums, c(model=m, w='same', round(summary(rowsm$bias), 3)))
}

for (m in c('_multi', '_separate1', '_separate2', '_pooled')){
  rowsm <- df[df$model == m & !df$parameter %in% c('theta', 'phi') & df$w =='different',]
  sums <- rbind(sums, c(model=m, w='different', round(summary(rowsm$bias), 3)))
}
sums <- data.frame(sums)
write_latex_table(sums, "latex_multi_bias_summary.txt", path=dst)
