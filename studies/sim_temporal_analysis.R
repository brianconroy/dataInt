############################
# Fit the spatiotemporal
# shared latent process model
############################

# ToDo:
# line plots (by trend)
# comparative box plots (by trend)


library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(plyr)
library(mvtnorm)
library(R.utils)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


src <- "/Users/brianconroy/Documents/research/dataInt/output/sim_temporal/"
dst <- "/Users/brianconroy/Documents/research/project3/sim_temporal/"
agg_factor <- 12
n_sims <- 25


#### Principal Components
years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
caPr_all <- list()
caPr.disc_all <- list()
for (h in years){
  caPr_y <- load_prism_pcs_time(h)
  caPr.disc_y <- aggregate(caPr_y, fact=agg_factor) 
  caPr_all <- c(caPr_all, caPr_y)
  caPr.disc_all <- c(caPr.disc_all, caPr.disc_y)
}
plot(caPr_all[[1]])
plot(caPr.disc_all[[1]])
print(n_values(caPr.disc_all[[1]][[1]]))
print(mean(area(caPr.disc_all[[1]][[1]])[]))

cells.all <- c(1:ncell(caPr.disc_all[[1]]))[!is.na(values(caPr.disc_all[[1]][[1]]))]
coords <- xyFromCell(caPr.disc_all[[1]], cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))


# plot principal components
pcs1 <- list()
pcs2 <- list()
for (i in 1:length(years)){
  pcs1[[i]] <- caPr.disc_all[[i]][[1]]
  pcs2[[i]] <- caPr.disc_all[[i]][[2]]
}
pcs1_eq <- equalize_scales4(pcs1)
pcs2_eq <- equalize_scales4(pcs2)
par(mfrow=c(2,4))
for (i in 1:7){
  plot(pcs1_eq[[i]], main=years[i])
}
par(mfrow=c(1,1))
par(mfrow=c(2,4))
for (i in 1:7){
  plot(pcs2_eq[[i]], main=years[i])
}
par(mfrow=c(1,1))


####################
# Make main analytics tables
# year | trend | N | rmse
####################
rmses_temporal <- array(NA, c(3*n_sims*length(years), 4))
rmses_pooled <- array(NA, c(3*n_sims*length(years), 4))
rmses_glm <- array(NA, c(3*n_sims*length(years), 4))
prevalences <- array(NA, c(3*n_sims*length(years), 3))
n_obs_cells <- array(NA, c(3*n_sims*length(years), 3))
biases <- array(NA, c(3*n_sims, 4))
rmses_w <- array(NA, c(3*n_sims, 3))
rmses_uw <- array(NA, c(3*n_sims*length(years), 4))
counter <- 1
counter2 <- 1
for (level in c("Increasing", "Decreasing", "Alternating")){
  for (i in 1:n_sims){
    
    print(paste("dataset", i))
    data <- reformat_saved_data(load_output(paste("data_", tolower(level), "_", i, ".json", sep=""), src=src))
    params <- load_output(paste("params_", tolower(level), "_", i, ".json", sep=""), src=src)
    
    sim_name <- paste("simPsTemporal", tolower(level), sep="_")
    output <- load_output(paste("output_", sim_name, "_", i, ".json", sep=""), src=src)
    output_pooled <- load_output(paste("output_pooled_", sim_name, "_", i, ".json", sep=""), src=src)
    
    #### Calculate RMSEs estimated log odds
    for (y_idx in 1:length(years)){
      year <- years[y_idx]
      params$alpha.case=params$Alpha.case
      params$alpha.ctrl=params$Alpha.ctrl
      lodds_true <- calc_true_lodds_time(params, data, year, y_idx, agg_factor=agg_factor)
      lodds_est_tmp <- calc_est_lodds_time(output, data, year, y_idx, agg_factor=agg_factor)
      lodds_est_pooled <- calc_est_lodds_time_ps(output_pooled, data, year, y_idx, agg_factor=agg_factor)
      lodds_est_glm <- calc_est_lodds_time_poisson(data, year, y_idx, agg_factor=agg_factor)
      n <- sum(data$locs[[y_idx]]$status)
      prev <- sum(data$case.data[[y_idx]]$y)/sum(data$case.data[[y_idx]]$y + data$ctrl.data[[y_idx]]$y)
      # store RMSEs
      rmses_temporal[counter,] <- c(year, level, n, sqrt(mean((lodds_true - lodds_est_tmp)^2)))
      rmses_pooled[counter,] <- c(year, level, n, sqrt(mean((lodds_true - lodds_est_pooled)^2)))
      rmses_glm[counter,] <- c(year, level, n, sqrt(mean((lodds_true - lodds_est_glm)^2)))
      u.hat <- colMeans(output$samples.u)[y_idx]
      w.hat <- colMeans(output$samples.w)
      rmses_uw[counter,] <- c(year, level, n, sqrt(mean((params$U[y_idx] + params$W - u.hat - w.hat)^2)))
      # store other details
      n_obs_cells[counter,] <- c(year, level, n)
      prevalences[counter,] <- c(year, level, prev)
      # update counter
      counter <- counter + 1
    }
    rmses_w[counter2,] <- c(level, n, sqrt(mean((params$W - w.hat)^2)))
    biases[counter2,1] <- mean(output$samples.alpha.ca) - params$Alpha.case
    biases[counter2,2] <- mean(output$samples.alpha.co) - params$Alpha.ctrl
    biases[counter2,3] <- mean(output$samples.theta) - params$Theta
    biases[counter2,4] <- mean(output$samples.phi) - params$Phi
    counter2 <- counter2 +  1
  }
}
# Reformat
rmses_temporal <- data.frame(rmses_temporal)
rmses_pooled <- data.frame(rmses_pooled)
rmses_glm <- data.frame(rmses_glm)
n_obs_cells <- data.frame(n_obs_cells)
prevalences <- data.frame(prevalences)
names(rmses_temporal) <- c("Year", "Trend", "N", "RMSE")
names(rmses_pooled) <- c("Year", "Trend", "N", "RMSE")
names(rmses_glm) <- c("Year", "Trend", "N", "RMSE")
names(n_obs_cells) <- c("Year", "Trend", "N")
names(prevalences) <- c("Year", "Trend", "Prevalence")
rmses_temporal$RMSE <- as.numeric(as.character(rmses_temporal$RMSE))
rmses_pooled$RMSE <- as.numeric(as.character(rmses_pooled$RMSE))
rmses_glm$RMSE <- as.numeric(as.character(rmses_glm$RMSE))
rmses_temporal$N <- as.numeric(as.character(rmses_temporal$N))
rmses_pooled$N <- as.numeric(as.character(rmses_pooled$N))
rmses_glm$N <- as.numeric(as.character(rmses_glm$N))
rmse_combined <- rbind(rmses_temporal, rmses_pooled, rmses_glm)
rmse_combined$Model <- c(rep("Spatiotemporal", nrow(rmses_temporal)), rep("Pooled", nrow(rmses_temporal)), rep("GLM", nrow(rmses_temporal)))
n_obs_cells$N <- as.numeric(as.character(n_obs_cells$N))
prevalences$Prevalence <- as.numeric(as.character(prevalences$Prevalence))
rmses_w <- data.frame(rmses_w)
names(rmses_w) <- c("Trend", "N", "RMSE")
rmses_w$RMSE <- as.numeric(as.character(rmses_w$RMSE))
rmses_uw <- data.frame(rmses_uw)
names(rmses_uw) <- c("Year", "Trend", "N", "RMSE")
rmses_uw$RMSE <- as.numeric(as.character(rmses_uw$RMSE))
rmses_uw$N <- as.numeric(as.character(rmses_uw$N))
biases <- data.frame(biases)
biases$Trend <- c(rep("Increasing", n_sims), rep("Decreasing", n_sims), rep("Alternating", n_sims))
names(biases) <- c("alpha_ca", "alpha_co", "theta", "phi", "Trend")

#########
# Figure:  RMSE boxplots
#########
rmse_increasing <- rmse_combined[rmse_combined$Trend == "Increasing",]
rmse_decreasing <- rmse_combined[rmse_combined$Trend == "Decreasing",]
rmse_alternating <- rmse_combined[rmse_combined$Trend == "Alternating",]

p1 <- ggplot() + 
  geom_boxplot(data = rmse_increasing, mapping = aes(Year, RMSE, fill=Model)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylim(0, 10) + 
  ggtitle("A)") +
  geom_hline(yintercept=0, color="red")
p2 <- ggplot() + 
  geom_boxplot(data = rmse_decreasing, mapping = aes(Year, RMSE, fill=Model)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylim(0, 10) + 
  ggtitle("B)") +
  geom_hline(yintercept=0, color="red")
p3 <- ggplot() + 
  geom_boxplot(data = rmse_alternating, mapping = aes(Year, RMSE, fill=Model)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylim(0, 10) + 
  ggtitle("C)") +
  geom_hline(yintercept=0, color="red")
grid.arrange(p1,p2,p3,ncol=2)

########
# Table: RMSE summary
# Trend | Year 1 (Mean Rmse) | ....
########
rmse_summary <- c()
for (df in list(rmse_increasing, rmse_decreasing, rmse_alternating)){
  agg_mean <- aggregate(df$RMSE,
                   by=list(Year=df$Year, Model=df$Model), 
                   FUN=function(x){mean(x, na.rm=T)})
  agg_sd <- aggregate(df$RMSE,
                   by=list(Year=df$Year, Model=df$Model), 
                   FUN=function(x){sd(x, na.rm=T)})
  for (model in c("Spatiotemporal", "Pooled", "GLM")){
    agg_m_mean <- agg_mean[agg_mean$Model == model,]
    agg_m_sd <- agg_sd[agg_sd$Model == model,]
    rmse_summary <- rbind(rmse_summary, round(agg_m_mean$x, 3))
    rmse_summary <- rbind(rmse_summary, round(agg_m_sd$x, 3))
  }
}
rmse_summary <- data.frame(rmse_summary)
rmse_summary <- cbind(c(rep(c("Mean", "SD"), 9)), rmse_summary)
rmse_summary <- cbind(c(rep(c("Spatiotemporal", "Spatiotemporal", "Pooled", "Pooled", "GLM", "GLM"), 3)), rmse_summary)
rmse_summary <- cbind(c(rep("Increasing", 6), rep("Decreasing", 6), rep("Alternating", 6)), rmse_summary)
names(rmse_summary) <- c("Trend", "Model", "Metric", years)
write_latex_table(rmse_summary, "rmse_summary", dst)

# Summarize mean observed cells over time for each trend
for (trend in c("Increasing", "Decreasing", "Alternating")){
  obs_trend <- n_obs_cells[n_obs_cells$Trend == trend,]
  summary_trend <- c()
  for (i in 1:length(years)){
    summary_trend <- c(summary_trend, round(mean(obs_trend[,i]), 2))
  }
  print(summary_trend)
}

# Summarize mean prevalences over time for each trend
for (trend in c("Increasing", "Decreasing", "Alternating")){
  prev_trend <- prevalences[prevalences$Trend == trend,]
  summary_trend <- c()
  for (i in 1:length(years)){
    summary_trend <- c(summary_trend, round(mean(prev_trend[,i]), 3))
  }
  print(summary_trend)
}

# Make lineplots
p1 <- ggplot(rmse_increasing[rmse_increasing$Model=="Spatiotemporal",], aes(x=N, y=RMSE, fill=Year)) +
  geom_line(aes(color=Year)) +
  geom_point() + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("A)") + 
  ylim(c(0,15)) + 
  geom_hline(yintercept=0, color="red")
p2 <- ggplot(rmse_decreasing[rmse_decreasing$Model=="Spatiotemporal",], aes(x=N, y=RMSE, fill=Year)) +
  geom_line(aes(color=Year)) +
  geom_point() + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("B)") + 
  ylim(c(0,15)) + 
  geom_hline(yintercept=0, color="red")
p3 <- ggplot(rmse_alternating[rmse_alternating$Model=="Spatiotemporal",], aes(x=N, y=RMSE, fill=Year)) +
  geom_line(aes(color=Year)) +
  geom_point() + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("C)") + 
  ylim(c(0,15)) + 
  geom_hline(yintercept=0, color="red")
grid.arrange(p1,p2,p3,ncol=2)

# Make W, U+W plots
p1 <- ggplot() + 
  geom_boxplot(data = rmses_w, mapping = aes(Trend, RMSE)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  geom_hline(yintercept=0, color="red") + 
  ylim(0,12) + 
  ggtitle("A)")
p2 <- ggplot() + 
  geom_boxplot(data = rmses_uw, mapping = aes(Year, RMSE, fill=Trend)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  geom_hline(yintercept=0, color="red") + 
  ylim(0,12) + 
  ggtitle("B)")
grid.arrange(p1,p2,ncol=2)

# summary of RMSE in W
summary_rmse_w <- list()
counter <- 1
for (trend in c("Increasing", "Decreasing", "Alternating")){
  rmses_w_t <- rmses_w[rmses_w$Trend == trend,]
  summary_rmse_w[[counter]] <- list(
    Trend=trend,
    Mean=round(mean(rmses_w_t$RMSE), 3),
    SD=round(sd(rmses_w_t$RMSE), 3),
    Median=round(median(rmses_w_t$RMSE), 3),
    Q1=round(quantile(rmses_w_t$RMSE, 0.25), 3),
    Q3=round(quantile(rmses_w_t$RMSE, 0.75), 3)
  )
  counter <- counter + 1
}
write_latex_table(ldply(summary_rmse_w, 'data.frame'), 'summary_rmse_w', dst)

# summary of RMSE in u + W
summary_rmse_uw <- c()
for (trend in c("Increasing", "Decreasing", "Alternating")){
  rmses_uw_t <- rmses_uw[rmses_uw$Trend == trend,]
  agg_mean <- aggregate(rmses_uw_t$RMSE,
                        by=list(Year=rmses_uw_t$Year), 
                        FUN=function(x){mean(x, na.rm=T)})
  agg_sd <- aggregate(rmses_uw_t$RMSE,
                      by=list(Year=rmses_uw_t$Year), 
                      FUN=function(x){sd(x, na.rm=T)})
  summary_rmse_uw <- rbind(summary_rmse_uw, round(agg_mean$x, 3))
  summary_rmse_uw <- rbind(summary_rmse_uw, round(agg_sd$x, 3))
}
summary_rmse_uw <- data.frame(summary_rmse_uw)
summary_rmse_uw <- cbind(rep(c("Mean", "SD"), 3), summary_rmse_uw)
summary_rmse_uw <- cbind(c(rep(c("Increasing"), 2), 
                         rep(c("Decreasing"), 2), 
                         rep(c("Alternating"), 2)), summary_rmse_uw)
names(summary_rmse_uw) <- c("Trend", "Metric", years)
write_latex_table(summary_rmse_uw, 'summary_rmse_uw', dst)

# Biases in other param estimates
bias_summary <- list()
counter <- 1
for (trend in c("Increasing", "Decreasing", "Alternating")){
  bias_t <- biases[biases$Trend == trend,]
  mean_biases <- round(colMeans(bias_t[,1:4]), 3)
  sd_biases <- round(apply(bias_t[,1:4], 2, sd), 3)
  bias_summary[[counter]] <- list(
    Trend=trend,
    Metric="Mean",
    Alpha_ca=mean_biases[1],
    Alpha_co=mean_biases[2],
    Theta=mean_biases[3],
    Phi=mean_biases[4]
  )
  bias_summary[[counter+1]] <- list(
    Trend=trend,
    Metric="SD",
    Alpha_ca=sd_biases[1],
    Alpha_co=sd_biases[2],
    Theta=sd_biases[3],
    Phi=mean_biases[4]
  )
  counter <- counter + 2
}
write_latex_table(ldply(bias_summary, 'data.frame'), 'bias_summary', dst)
