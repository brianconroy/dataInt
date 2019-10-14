############################
# Fit the spatiotemporal
# shared latent process model
############################


library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_mvgp/"
agg_factor <- 12
n_sims <- 25
sim_name <- "sim_mvgp"


#### Prism Principal Components
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=agg_factor)
n_values(caPr.disc[[1]])
plot(caPr.disc)
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc, cell=cells.all)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


#####################
#### RMSE in log odds
#####################
rmse_mvgp <- array(NA, c(n_sims, 2))
n_cells_obs <- array(NA, c(n_sims, 2))
prevalences <- array(NA, c(n_sims, 2))
bias_alpha_ca <- array(NA, c(n_sims, 2))
bias_alpha_co <- array(NA, c(n_sims, 2))
for (i in 1:14){
  print(paste("dataset", i))
  data <- reformat_saved_mvgp(load_output(paste("data_", i, ".json", sep=""), src=src))
  params <- load_output(paste("params_", i, ".json", sep=""), src=src)
  
  output_mvgp <- load_output(paste("output_", sim_name, "_", i, ".json", sep=""), src=src)
  
  par(mfrow=c(1,2))
  for (k in 1:2){
    # true log odds
    lodds_true <- calc_lodds_true_multi(params, data, species=k, agg_factor=agg_factor)
    
    # proposed model
    lodds <- calc_lodds_mvgp(output_mvgp, data, species=k, agg_factor=agg_factor)
    rmse_mvgp[i,k] <- sqrt(mean((lodds-lodds_true)^2))
    n_cells_obs[i,k] <- sum(data$locs[[k]]$status)
    prevalences[i,k] <- sum(data$case.data[[k]]$y)/sum(data$case.data[[k]]$y + data$ctrl.data[[k]]$y)
    plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', main=paste("species", k)); abline(0, 1, col=2)
    
    # parameter biases
    bias_alpha_ca[i,k] <- mean(output_mvgp$samples.alpha.ca[k,,]) - params$alpha.cases[k]
    bias_alpha_co[i,k] <- mean(output_mvgp$samples.alpha.co[k,,]) - params$alpha.ctrls[k]
  }
}

# Summary of datasets
data_summary <- list()
data_summary[[1]] <- list(
  Species=1,
  Attribute="N Cells",
  Mean=round(mean(n_cells_obs[,1]), 2),
  Median=round(median(n_cells_obs[,1]), 2),
  Q1=round(quantile(n_cells_obs[,1], 0.25), 2),
  Q3=round(quantile(n_cells_obs[,1], 0.75), 2)
)
data_summary[[2]] <- list(
  Species=1,
  Attribute="Prevalence",
  Mean=round(mean(prevalences[,1]), 2),
  Median=round(median(prevalences[,1]), 2),
  Q1=round(quantile(prevalences[,1], 0.25), 2),
  Q3=round(quantile(prevalences[,1], 0.75), 2)
)
data_summary[[3]] <- list(
  Species=2,
  Attribute="N Cells",
  Mean=round(mean(n_cells_obs[,2]), 2),
  Median=round(median(n_cells_obs[,2]), 2),
  Q1=round(quantile(n_cells_obs[,2], 0.25), 2),
  Q3=round(quantile(n_cells_obs[,2], 0.75), 2)
)
data_summary[[4]] <- list(
  Species=2,
  Attribute="Prevalence",
  Mean=round(mean(prevalences[,2]), 2),
  Median=round(median(prevalences[,2]), 2),
  Q1=round(quantile(prevalences[,2], 0.25), 2),
  Q3=round(quantile(prevalences[,2], 0.75), 2)
)

# Boxplot

# RMSE summary  table

# Lineplot: RMSE vs N
df_n <- c()
for (i in 1:n_sims){
  for (k in 1:ncol(rmse_mvgp)){
    df_n <- rbind(df_n, cbind(k, rmse_mvgp[i,k], n_cells_obs[i,k]))
  }
}
df_n <- data.frame(df_n)
names(df_n) <- c("Species", "RMSE", "N")
df_n$Species <- as.character(df_n$Species)

ggplot(df_n, aes(x=N, y=RMSE, group=Species)) +
  geom_line(aes(color=Species)) +
  geom_point(aes(color=Species)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  geom_hline(yintercept=0, color="red") # + ylim(-3.5, 2.5)

# Biases in alphas
df_aca <- c()
for (i in 1:n_sims){
  for (k in 1:ncol(rmse_mvgp)){
    df_aca <- rbind(df_aca, cbind(k, rmse_mvgp[i,k], n_cells_obs[i,k]))
  }
}
df_aca <- data.frame(df_aca)
names(df_aca) <- c("Species", "Bias", "N")
df_aca$Species <- as.character(df_aca$Species)


# Biases in T params

