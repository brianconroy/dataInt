###############################
# Simulation study to
# compare the proposed
# model against benchmarks
# (poisson and spatial poisson)
###############################


library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


sampling <- "low"
src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration/"
sim_name <- paste("sim_iteration_v2_", sampling, sep="")
n_sims <- 25
agg_factor <- 10


#### Prism Principal Components
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=agg_factor)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


######
# RMSE
######


rmse_ps <- c()
rmse_pr <- c()
n_cells_obs <- c()
for (i in 1:n_sims){
  print(i)
  
  params <- load_output(paste("params_low_", i, ".json", sep=""), src=src)
  Theta <- params$Theta
  Phi <- params$Phi
  Alpha.case <- params$alpha.case
  Alpha.ctrl <- params$alpha.ctrl
  beta.case <- params$beta.case
  beta.ctrl <- params$beta.ctrl
  W <- params$W
  
  output <- load_output(paste("output_", sim_name, "_", i, ".json", sep=""), src=src)
 
  # Check estimated log odds
  w.hat <- colMeans(output$samples.w)
  beta_ca_h <- colMeans(output$samples.beta.ca)
  beta_co_h <- colMeans(output$samples.beta.co)
  alpha_ca_h <- mean(output$samples.alpha.ca)
  alpha_co_h <- mean(output$samples.alpha.co)
  phi_h <- mean(output$samples.phi)
  theta_h <- mean(output$samples.theta)
  
  # RMSE
  X.standard <- load_x_standard2(as.logical(data$locs$status), agg_factor=agg_factor)
  lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
  lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
  plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_ps <- c(rmse_ps, sqrt(mean((lodds.true-lodds.ps)^2)))
  
  # Reference model (Poisson Regression)
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  X.ca <- data$case.data$x.standardised
  Y.ca <- data$case.data$y
  X.co <- data$ctrl.data$x.standardised
  Y.co <- data$ctrl.data$y
  # Cases
  rmodel.ca <- glm(Y.ca ~ X.ca-1, family='poisson')
  beta_ca_r <- coefficients(rmodel.ca)
  # Controls
  rmodel.co <- glm(Y.co ~ X.co-1, family='poisson')
  beta_co_r <- coefficients(rmodel.co)
  # Log odds
  lodds.r <- X.standard %*% beta_ca_r - X.standard %*% beta_co_r
  plot(x=lodds.true, y=lodds.r, main='reference', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_pr <- c(rmse_pr, sqrt(mean((lodds.true-lodds.r)^2)))
  
  n_cells_obs <- c(n_cells_obs, sum(data$locs$status))
  
}

boxplot(rmse_ps) 
summary(rmse_ps)
summary(rmse_pr)

df <- cbind(rmse_ps, rmse_pr)
boxplot(df, col=rainbow(2, s=0.5))
plot(x=n_cells_obs, y=rmse_ps)
