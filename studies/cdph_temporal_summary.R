#############################
# Summarize the discrete time
# spatiotemporal model fit
# to Sciurid data
#############################

library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(mvtnorm)
library(R.utils)
sourceDirectory('Documents/research/dataInt/R/')


analysis_name <- "cdph_temporal_analysis"
agg_factor <- 7
dst <- "/Users/brianconroy/Documents/research/project3/analysis/p3_figs_analysis/"
src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")
output <- load_output(paste("output_", analysis_name, ".json", sep=""))
output_agg <- load_output(paste('output_', paste(analysis_name, "_aggregated_ps", sep=""), ".json", sep=""))


#### Prism Principal Components
years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
caPr_all <- list()
caPr.disc_all <- list()
for (h in years){
  caPr_y <- load_prism_pcs_time(h)
  caPr.disc_y <- aggregate(caPr_y, fact=agg_factor) 
  caPr_all <- c(caPr_all, caPr_y)
  caPr.disc_all <- c(caPr.disc_all, caPr.disc_y)
  plot(caPr_y)
}
plot(caPr_all[[1]])
plot(caPr.disc_all[[1]])
loc.disc_y <- caPr.disc_all[[1]]


#### Assemble data
locs <- list()
years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
case.data <- list()
ctrl.data <- list()
bin_width <- 5

for (i in 1:length(years)){
  y <- years[i]
  y_ub <- y + (bin_width - 1)/2
  y_lb <- y - (bin_width - 1)/2
  
  rodents_y <- rodents[rodents$Year <= y_ub & rodents$Year >= y_lb,]
  loc.disc_y <- caPr.disc_all[[i]][[1]]
  caPr.disc_y <- caPr.disc_all[[i]]
  data_y <- assemble_data(rodents_y, loc.disc_y, caPr.disc_y)
  
  locs[[i]] <- data_y$loc
  case.data[[i]] <- data_y$case.data
  ctrl.data[[i]] <- data_y$ctrl.data
}


cells.all <- c(1:ncell(caPr.disc_all[[1]]))[!is.na(values(caPr.disc_all[[1]][[1]]))]
coords <- xyFromCell(caPr.disc_all[[1]], cell=cells.all)
D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
data <- list(locs=locs, case.data=case.data, ctrl.data=ctrl.data)


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

#### Figure: pcs 2
par(mfrow=c(3,3))
for (i in 1:length(years)){
  plot(pcs2_re[[i]], main=years[i])
}


#### Table: summarize specimen counts/prevalences/observation sites per year
prevs <- list()
counter <- 1
for (i in 1:length(data$locs)){
  ncases <- sum(data$case.data[[i]]$y)
  nctrls <- sum(data$ctrl.data[[i]]$y)
  prev <- round(ncases/(ncases+nctrls), 2)
  row_i <- list(
    Time_Interval=years[i],
    Observation_Sites=sum(data$locs[[i]]$status),
    Total_Specimen=ncases+nctrls,
    Prevalence=prev
  )
  prevs[[i]] <- row_i
  counter <- counter + 1
}
write_latex_table(ldply(prevs, 'data.frame'), "latex_counts_prev.txt", path=dst)


#### Model fitting details
## Temporal model
print(output$n.sample)
print(output$burnin)
print(ncol(D))
print(mean(area(caPr.disc_all[[1]][[1]])[]))
## Reference model
print(output_agg$n.sample)
print(output_agg$burnin)


#### Figure: Plot temporal risk surfaces
temporal_risks <- calc_temporal_risks(output, caPr.disc_all, caPr_all, agg_factor, years)

rescaled <- equalize_scales4(temporal_risks$r_risks_high)
par(mfrow=c(2,4))
for (i in 1:length(years)){
  plot(rescaled[[i]], main=years[i])
}


#### Figure: Plot risk surface from reference method
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=agg_factor)
loc.disc <- caPr.disc[[1]]
all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
risk_agg <- calc_risk_cdph_output(
  output_agg, 
  data, 
  caPr.disc, 
  all_ids, 
  agg_factor=agg_factor)$r_risk_high
plot(risk_agg)


#### Figure: compare prevalences and process
w.hat <- colMeans(output$samples.w)
u.hat <- colMeans(output$samples.u)
prevalences <- c()
n_obs <- c()
for (i in 1:length(years)){
  prevalences <- c(prevalences,
                   sum(data$case.data[[i]]$y)/(sum(data$case.data[[i]]$y) + sum(data$ctrl.data[[i]]$y)))
  n_obs <- c(n_obs, sum(data$locs[[i]]$status))
}

# quantiles
samples <- array(NA, c(nrow(output$samples.w), length(years)))
for (i in 1:nrow(output$samples.w)){
  w.hat_i <- mean(output$samples.w[i,])
  for (t in 1:length(years)){
    samples[i,t] <- w.hat_i + output$samples.u[i,t]
  }
}
q25 <- c()
q75 <- c()
for (t in 1:ncol(samples)){
  q25 <- c(q25, quantile(samples[,t], probs=c(0.25)))
  q75 <- c(q75, quantile(samples[,t], probs=c(0.75)))
}

dat <- cbind(years, mean(w.hat) + u.hat, q25, q75)
dat <- data.frame(dat)
names(dat) <- c('years', 'mean', 'q25', 'q75')

p1 <- ggplot(dat, aes(x=years, y=mean)) + 
  geom_ribbon(aes(ymin=q25, ymax=q75), alpha=0.2) + 
  geom_line() + 
  scale_x_continuous(breaks=years) + 
  xlab("Year") + ylab("Mean") + ggtitle("A)")

p2 <- ggplot(data.frame(cbind(years, prevalences)), aes(x=years, y=prevalences)) + 
  geom_line() + 
  scale_x_continuous(breaks=years) +
  xlab("Year") + ylab("Prevalence") + ggtitle("B)")

p3 <- ggplot(data.frame(cbind(years, n_obs)), aes(x=years, y=n_obs)) + 
  geom_line() + 
  scale_x_continuous(breaks=years) +
  xlab("Year") + ylab("Observed Cells") + ggtitle("C)")

grid.arrange(p1, p2, p3)


#### Figure: scatterplots of values for spatiotemp model vs aggregate model
par(mfrow=c(2,4))
for (t in 1:length(years)){
  r_t <- temporal_risks$risks_high[[t]]
  r_a <- risk_agg[][!is.na(risk_agg[])]
  plot(x=r_a, y=r_t, main=years[t], ylab='Temporal Risk', xlab='Aggregate Risk'); abline(0,1,col=2)
  print(round(mean(r_t - r_a), 3))
}


#### Figure: significance value comparison
rs <- calc_significance_rasters_ds(rodents, output_agg, caPr, threshold=0.05)
plot(rs$r_fracs)
plot(rs$r_inds)
plot(rs$r_inds_95)
plot_risk_overlay(rs$r_inds_95, rodents)
plot_risk_overlay(rs$r_inds, rodents)

rs_temporal <- calc_significance_rasters_ds_temporal(data, output, caPr_all, caPr.disc_all, threshold=0.05)
par(mfrow=c(2,4))
for (t in 1:length(years)){
  plot(rs_temporal[[t]]$r_inds_95, main=years[t])
}

baseline <- sum(rs$r_inds_95[][!is.na(rs$r_inds_95[])])
frac_baseline <- c()
for (t in 1:length(years)){
  frac_baseline <- c(
    frac_baseline,
    round(sum(rs_temporal[[t]]$r_inds_95[][!is.na(rs_temporal[[t]]$r_inds_95[])])/baseline, 3)
  )
}


