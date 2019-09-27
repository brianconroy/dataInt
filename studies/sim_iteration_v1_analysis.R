###############################
# Simulation study to
# compare the proposed
# model against benchmarks
# (poisson and spatial poisson)
###############################


library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration_v1/"
agg_factor <- 11
n_sims <- 25

#### Prism Principal Components
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=agg_factor)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


######
# RMSE: low 
######


sampling <- "low"
sim_name <- paste("sim_iteration_v1_", sampling, sep="")

rmse_ps_low <- c()
rmse_pr_low <- c()
rmse_sp_low <- c()
n_cells_obs_low <- c()
prevalences_low <- c()
for (i in 1:n_sims){
  print(i)
  
  data <- load_output(paste("data_low_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_", sampling, "_", i, ".json", sep=""), src=src)
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
  rmse_ps_low <- c(rmse_ps_low, sqrt(mean((lodds.true-lodds.ps)^2)))
  
  # Reference model (Poisson Regression)
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
  rmse_pr_low <- c(rmse_pr_low, sqrt(mean((lodds.true-lodds.r)^2)))
  
  # Reference model (Spatial Poisson Regression)
  # output_sp_ca <- load_output(paste("output.sp_ca_", sim_name, "_", i, ".json", sep=""), src=src)
  # output_sp_co <- load_output(paste("output.sp_co_", sim_name, "_", i,  ".json", sep=""), src=src)
  # kriged_w_ca <- load_output(paste("output.krige_ca_", sim_name, "_", i, ".json", sep=""), src=src)
  # kriged_w_co <- load_output(paste("output.krige_co_", sim_name, "_", i, ".json", sep=""), src=src)
  # w.hat_spca <- colMeans(output_sp_ca$samples.w)
  # w.hat_spco <- colMeans(output_sp_co$samples.w)
  # w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(data$locs$status))
  # w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, as.logical(data$locs$status))
  # beta_co_sp <- colMeans(output_sp_co$samples.beta)
  # beta_ca_sp <- colMeans(output_sp_ca$samples.beta)
  # lodds_sp <- X.standard %*% beta_ca_sp + w_ca_est - X.standard %*% beta_co_sp - w_co_est
  # rmse_sp_low <- c(rmse_sp_low, sqrt(mean((lodds.true-lodds_sp)^2)))
  # plot(x=lodds.true, y=lodds_sp, main='reference 2', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  # 
  n_cells_obs_low <- c(n_cells_obs_low, sum(data$locs$status))
  prevalences_low <- c(prevalences_low, sum(data$case.data$y)/(sum(data$case.data$y + data$ctrl.data$y)))
  
}

# Initial inspection
summary(rmse_ps_low)
summary(rmse_sp_low)
summary(rmse_pr_low)

boxplot(cbind(rmse_ps_low, rmse_sp_low, rmse_pr_low), col=rainbow(3, s=0.5))
plot(x=n_cells_obs_low, y=rmse_ps_low)
points(x=n_cells_obs_low, y=rmse_pr_low, col=2)


######
# RMSE: high
######


sampling <- "high"
src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_iteration/"
sim_name <- paste("sim_iteration_v1_", sampling, sep="")

rmse_ps_high <- c()
rmse_pr_high <- c()
rmse_sp_high <- c()
n_cells_obs_high <- c()
prevalences_high <- c()
for (i in 1:n_sims){
  print(i)
  
  data <- load_output(paste("data_high_", i, ".json", sep=""), src=src)
  params <- load_output(paste("params_high_", i, ".json", sep=""), src=src)
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
  
  # # RMSE
  X.standard <- load_x_standard2(as.logical(data$locs$status), agg_factor=agg_factor)
  lodds.true <- X.standard %*% beta.case + Alpha.case * W - X.standard %*% beta.ctrl - Alpha.ctrl * W
  lodds.ps <- X.standard %*% beta_ca_h + alpha_ca_h * w.hat - X.standard %*% beta_co_h - alpha_co_h * w.hat
  plot(x=lodds.true, y=lodds.ps, main='A)', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  rmse_ps_high <- c(rmse_ps_high, sqrt(mean((lodds.true-lodds.ps)^2)))
  
  # Reference model (Poisson Regression)
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
  rmse_pr_high <- c(rmse_pr_high, sqrt(mean((lodds.true-lodds.r)^2)))
  
  # Reference model (Spatial Poisson Regression)
  output_sp_ca <- load_output(paste("output.sp_ca_", sim_name, "_", i, ".json", sep=""), src=src)
  output_sp_co <- load_output(paste("output.sp_co_", sim_name, "_", i,  ".json", sep=""), src=src)
  kriged_w_ca <- load_output(paste("output.krige_ca_", sim_name, "_", i, ".json", sep=""), src=src)
  kriged_w_co <- load_output(paste("output.krige_co_", sim_name, "_", i, ".json", sep=""), src=src)
  w.hat_spca <- colMeans(output_sp_ca$samples.w)
  w.hat_spco <- colMeans(output_sp_co$samples.w)
  w_ca_est <- combine_w(w.hat_spca, kriged_w_ca$mu.new, as.logical(data$locs$status))
  w_co_est <- combine_w(w.hat_spco, kriged_w_co$mu.new, as.logical(data$locs$status))
  beta_co_sp <- colMeans(output_sp_co$samples.beta)
  beta_ca_sp <- colMeans(output_sp_ca$samples.beta)
  lodds_sp <- X.standard %*% beta_ca_sp + w_ca_est - X.standard %*% beta_co_sp - w_co_est
  rmse_sp_high <- c(rmse_sp_high, sqrt(mean((lodds.true-lodds_sp)^2)))
  plot(x=lodds.true, y=lodds_sp, main='reference 2', xlab='True Log Odds', ylab='Estimated Log Odds'); abline(0, 1, col='2')
  
  
  n_cells_obs_high <- c(n_cells_obs_high, sum(data$locs$status))
  prevalences_high <- c(prevalences_high, sum(data$case.data$y)/(sum(data$case.data$y + data$ctrl.data$y)))
  
}

# initial inspection
summary(rmse_ps_high)
summary(rmse_sp_high)
summary(rmse_pr_high)

boxplot(cbind(rmse_ps_high, rmse_sp_high, rmse_pr_high), col=rainbow(3, s=0.5))
plot(x=n_cells_obs_high, y=rmse_ps_high)
points(x=n_cells_obs_high, y=rmse_pr_high, col=2)


#########
# Figure: RMSE boxplots
#########

## All datasets

# High
df2 <- data.frame(cbind(c(rmse_ps_high, rmse_sp_high, rmse_pr_high)))
names(df2) <- "RMSE"
df2$Model <- c(rep("PS", length(rmse_ps_high)), 
               rep("SP", length(rmse_ps_high)),
               rep("PR", length(rmse_pr_high)))
df2$Sampling <- rep("High", 3*length(rmse_ps_high))

# Low
df1 <- data.frame(cbind(c(rmse_ps_low, rmse_sp_low, rmse_pr_low)))
names(df1) <- "RMSE"
df1$Model <- c(rep("PS", length(rmse_ps_high)), 
               rep("SP", length(rmse_ps_high)),
               rep("PR", length(rmse_pr_high)))
df1$Sampling <- rep("Low", 3*length(rmse_ps_high))

# Combined
df <- rbind(df1, df2)
df$Model <- factor(df$Model,levels=c("PS", "SP", "PR"))

# Plot
p1 <- ggplot() + 
  geom_boxplot(data = df, mapping = aes(Sampling, RMSE, fill=Model)) + 
  scale_x_discrete(limits=c("Low", "High")) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("A)") + 
  ylim(0,10.5)

## Datasets with >= 75 observed cells

high_geq75 <- n_cells_obs_high >= 75
low_geq75 <- n_cells_obs_low >= 75

# High
df2_75 <- data.frame(cbind(c(rmse_ps_high[high_geq75], rmse_sp_high[high_geq75], rmse_pr_high[high_geq75])))
names(df2_75) <- "RMSE"
df2_75$Model <- c(rep("PS", sum(high_geq75)), 
                  rep("SP", sum(high_geq75)),
                  rep("PR", sum(high_geq75)))
df2_75$Sampling <- rep("High", 3*sum(high_geq75))

# Low
df1_75 <- data.frame(cbind(c(rmse_ps_low[low_geq75], rmse_sp_low[low_geq75], rmse_pr_low[low_geq75])))
names(df1_75) <- "RMSE"
df1_75$Model <- c(rep("PS", sum(low_geq75)), 
                  rep("SP", sum(low_geq75)),
                  rep("PR", sum(low_geq75)))
df1_75$Sampling <- rep("Low", 3*sum(low_geq75))

# Combined
df_75 <- rbind(df1_75, df2_75)
df_75$Model <- factor(df_75$Model,levels=c("PS", "SP", "PR"))

# Plot
p2 <- ggplot() + 
  geom_boxplot(data = df_75, mapping = aes(Sampling, RMSE, fill=Model)) + 
  scale_x_discrete(limits=c("Low", "High")) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("B)") + 
  ylim(0,10.5)

grid.arrange(p1, p2, ncol=2)

#########
# Figure: RMSE vs observed cells and prevalences
#########

# Observed cells
df1 <- data.frame(cbind(n_cells_obs_high, rmse_ps_high))
names(df1) <- c("N", "RMSE")
df1$Sampling <- "High"

df2 <- data.frame(cbind(n_cells_obs_low, rmse_ps_low))
names(df2) <- c("N", "RMSE")
df2$Sampling <- "Low"

df <- rbind(df1, df2)

p1 <- ggplot(df, aes(x=N, y=RMSE, group=Sampling)) +
  geom_line(aes(color=Sampling)) +
  geom_point(aes(color=Sampling)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("A)")

# Prevalences
df_prev1 <- data.frame(cbind(prevalences_high, rmse_ps_high))
names(df_prev1) <- c("Prevalence", "RMSE")
df_prev1$Sampling <- "High"

df_prev2 <- data.frame(cbind(prevalences_low, rmse_ps_low))
names(df_prev2) <- c("Prevalence", "RMSE")
df_prev2$Sampling <- "Low"

df_prev <- rbind(df_prev1, df_prev2)

p2 <- ggplot(df_prev, aes(x=Prevalence, y=RMSE, group=Sampling)) +
  geom_line(aes(color=Sampling)) +
  geom_point(aes(color=Sampling)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  ggtitle("B)")

grid.arrange(p1, p2, ncol=2)


#########
# Figure: Plot of study region covariates
#########
par(mfrow=c(1,2))
plot(caPr.disc[[1]], main='')
plot(caPr.disc[[2]], main='')


########
# Table: RMSE summary
########
rmse_table <- list()
rmse_table[[1]] <- list(
  Sampling="Low",
  Model="PS",
  Mean_RMSE=round(mean(rmse_ps_low), 3),
  SD_RMSE=round(sd(rmse_ps_low), 3),
  Mean_RMSE_above=round(mean(rmse_ps_low[n_cells_obs_low >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_ps_low[n_cells_obs_low >= 75]), 3)
)
rmse_table[[2]] <- list(
  Sampling="Low",
  Model="SP",
  Mean_RMSE=round(mean(rmse_sp_low), 3),
  SD_RMSE=round(sd(rmse_sp_low), 3),
  Mean_RMSE_above=round(mean(rmse_sp_low[n_cells_obs_low >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_sp_low[n_cells_obs_low >= 75]), 3)
)
rmse_table[[3]] <- list(
  Sampling="Low",
  Model="PR",
  Mean_RMSE=round(mean(rmse_pr_low), 3),
  SD_RMSE=round(sd(rmse_pr_low), 3),
  Mean_RMSE_above=round(mean(rmse_pr_low[n_cells_obs_low >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_pr_low[n_cells_obs_low >= 75]), 3)
)
rmse_table[[4]] <- list(
  Sampling="High",
  Model="PS",
  Mean_RMSE=round(mean(rmse_ps_high), 3),
  SD_RMSE=round(sd(rmse_ps_high), 3),
  Mean_RMSE_above=round(mean(rmse_ps_high[n_cells_obs_high >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_ps_high[n_cells_obs_high >= 75]), 3)
)
rmse_table[[5]] <- list(
  Sampling="High",
  Model="SP",
  Mean_RMSE=round(mean(rmse_sp_high), 3),
  SD_RMSE=round(sd(rmse_sp_high), 3),
  Mean_RMSE_above=round(mean(rmse_sp_high[n_cells_obs_high >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_sp_high[n_cells_obs_high >= 75]), 3)
)
rmse_table[[6]] <- list(
  Sampling="High",
  Model="PR",
  Mean_RMSE=round(mean(rmse_pr_high), 3),
  SD_RMSE=round(sd(rmse_pr_high), 3),
  Mean_RMSE_above=round(mean(rmse_pr_high[n_cells_obs_high >= 75]), 3),
  SD_RMSE_above=round(sd(rmse_pr_high[n_cells_obs_high >= 75]), 3)
)
write_latex_table(ldply(rmse_table, "data.frame"), fname="rmse_table.txt", 
                  path="/Users/brianconroy/Documents/research/paper/")


########
# Table: Dataset summary
########
table_dataset <- list()
table_dataset[[1]] <- list(
  Sampling="Low",
  Quantity="Observed Cells",
  Mean=round(mean(n_cells_obs_low), 3),
  SD=round(sd(n_cells_obs_low), 3),
  Q1=quantile(n_cells_obs_low, 0.25),
  Q2=quantile(n_cells_obs_low, 0.75)
)
table_dataset[[2]] <- list(
  Sampling="Low",
  Quantity="Prevalence",
  Mean=round(mean(prevalences_low), 3),
  SD=round(sd(prevalences_low), 3),
  Q1=round(quantile(prevalences_low, 0.25), 3),
  Q2=round(quantile(prevalences_low, 0.75), 3)
)
table_dataset[[3]] <- list(
  Sampling="High",
  Quantity="Observed Cells",
  Mean=round(mean(n_cells_obs_high), 3),
  SD=round(sd(n_cells_obs_high), 3),
  Q1=quantile(n_cells_obs_high, 0.25),
  Q2=quantile(n_cells_obs_high, 0.75)
)
table_dataset[[4]] <- list(
  Sampling="High",
  Quantity="Prevalence",
  Mean=round(mean(prevalences_high), 3),
  SD=round(sd(prevalences_high), 3),
  Q1=round(quantile(prevalences_high, 0.25), 3),
  Q2=round(quantile(prevalences_high, 0.75), 3)
)
write_latex_table(ldply(table_dataset, "data.frame"), fname="dataset_table.txt", 
                  path="/Users/brianconroy/Documents/research/paper/")
