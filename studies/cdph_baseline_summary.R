###############################
# Summarizes simulation results
# of the baseline (scurid)
# cdph analysis
###############################

library(plyr)
library(grid)
library(mvtnorm)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


dst <- "/Users/brianconroy/Documents/research/project1/cdph_baseline/"
caPr <- load_prism_pcs()
caPr.disc <- aggregate(caPr, fact=5)
N <- n_values(caPr.disc[[1]])
print(N)
print(mean(area(caPr.disc[[1]])[]))
plot(caPr.disc)

src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")
output <- load_output("output_cdph_baseline.json")

##################
# data description
##################

coords_all <- cbind(matrix(rodents$Lon_Add_Fix), rodents$Lat_Add_Fix)
loc.disc <- caPr.disc[[1]]
cells_all <- cellFromXY(loc.disc, coords_all)
cells_all <- cells_all[!is.na(cells_all[])]
cells_obs <- sort(unique(cells_all))

loc.disc[][!is.na(loc.disc[])] <- 0
loc.disc[][cells_obs] <- 5

plot(loc.disc)
plot(rasterToPolygons(loc.disc), add=T, border='black', lwd=1) 
points(coords_all, col='2')

##################
# data description
##################

# positive counts at each cell
rodents_pos <- rodents[rodents$Res == 'POS',]
coords_pos <- cbind(matrix(rodents_pos$Lon_Add_Fix), rodents_pos$Lat_Add_Fix)
cells_pos <- cellFromXY(loc.disc, coords_pos)
counts_pos <- data.frame(table(cells_pos))
names(counts_pos) <- c('cell', 'count_pos')

# negative counts at each cell
rodents_neg <- rodents[rodents$Res == 'NEG',]
coords_neg <- cbind(matrix(rodents_neg$Lon_Add_Fix), rodents_neg$Lat_Add_Fix)
cells_neg <- cellFromXY(loc.disc, coords_neg)
counts_neg <- data.frame(table(cells_neg))
names(counts_neg) <- c('cell', 'count_neg')

# combine counts
counts_all <- merge(counts_pos, counts_neg, by='cell', all=T)
counts_all$cell <- as.numeric(as.character(counts_all$cell))
counts_all[is.na(counts_all$count_pos),]$count_pos <- 0
counts_all[is.na(counts_all$count_neg),]$count_neg <- 0
counts_all <- counts_all[with(counts_all, order(cell)),]

# location data
all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
locs <- list(
  cells=cells_obs,
  status=1 * c(all_ids %in% cells_obs),  
  coords=xyFromCell(loc.disc, cells_obs)
)
locs$ids <- c(1:length(all_ids))[as.logical(locs$status)]
plot(loc.disc)
points(locs$coords)

# case data
cov.disc <- caPr.disc
x1 <- cov.disc[[1]][][locs$cells]
x2 <- cov.disc[[2]][][locs$cells]
x1.standardised <- (x1 - mean(x1))/sd(x1)
x2.standardised <- (x2 - mean(x2))/sd(x2)
x <- cbind(1, x1, x2)
x.standardised <- cbind(1, x1.standardised, x2.standardised)

case.data <- list(
  y=counts_all$count_pos,
  x.standardised=x.standardised,
  x=x,
  p=3
)
print(sum(case.data$y))

# control data
ctrl.data <- list(
  y=counts_all$count_neg,
  x.standardised=x.standardised,
  x=x,
  p=3
)

coords <- xyFromCell(caPr.disc, cell=all_ids)
d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
data <- list(loc=locs, case.data=case.data, ctrl.data=ctrl.data)

# table of summary metrics
count_sums <- rbind(summary(case.data$y), summary(ctrl.data$y))
count_sums <- cbind(c('Plague Positive', 'Plague Negative'), count_sums)
count_sums <- data.frame(count_sums)
count_sums$Total <- c(sum(case.data$y), sum(ctrl.data$y))
names(count_sums) <- c("Disease Status", "Min", "1st Quartile",
                       "Median", "Mean", "3rd Quartile", "Max", "Total")
write_latex_table(count_sums, "cdph_baseline_count_summary.txt", path=dst)


##############
# patch output
##############
prior_alpha_ca_mean <- 0
prior_alpha_co_mean <- 0
prior_alpha_ca_var <- 6
prior_alpha_co_var <- 6
output <- burnin_after(output, n.burn=200)
output <- continueMCMC(data, output, n.sample=3000)
save_output(output, "output_cdph_baseline.json")


###################
# model convergence
###################


# traceplots
fname <- paste("cdph_baseline_traces", ".png", sep="")
png(paste(dst, fname, sep=""),
     width=900, height=700, res=100)
par(mfrow=c(2,3))
plot_traces_general(output)



# log odds scatterplots
par(mfrow=c(2,2))
xl <- c(-22, 10)
yl <- c(-22, 10)
rmses <- c()
true_params <- load_params('true_params_priorcompare.json')
labs <- c("A)", "B)", "C)", "D)")
for (i in 1:length(priors)){
  
  p <- priors[i]
  lab <- labs[i]
  o <- get_output_general(outputs, tag=paste('priorcompare_', p, sep=""))
  lodds <- calc_log_odds_output(o, true_params)
  lodds_true <- calc_log_odds_true_general(true_params)
  rmse <- sqrt(mean((lodds-lodds_true)^2))
  rmses <- c(rmses, rmse)
  plot(x=lodds_true, y=lodds, xlab='True Log Odds', ylab='Estimated Log Odds', xlim=xl, ylim=yl, main=lab); abline(0, 1, col=2)

}


# log odds rmses table
rows <- list()
counter <- 1
for ( i in 1:length(priors)){
  p <- priors[i]
  rmse <- round(rmses[i], 3)
  new_row <- list(
    configuration=p,
    rmse=rmse
  )
  rows[[counter]] <- new_row
  counter <- counter + 1
}
df_rmse <- ldply(rows, 'data.frame')
write_latex_table(df_rmse, "latex_priorcompare_rmse.txt", path=dst)


# traceplots
true_params <- load_params('true_params_priorcompare.json')
for (p in priors){
  
  fname <- paste("priorcompare_traces_", p, ".png", sep="")
  png(paste(dst, fname, sep=""),
      width=900, height=700, res=100)
  par(mfrow=c(2,3))
  o <- get_output_general(outputs, tag=paste('priorcompare_', p, sep=""))
  plot_traces_general(o, true_params)
  dev.off()
  
}


# w scatterplots
par(mfrow=c(2,2))
counter <- 1
truevals <- load_params('true_params_priorcompare.json')
for (p in priors){
  
  o <- get_output_general(outputs, tag=paste('priorcompare_', p, sep=""))
  w.hat <- colMeans(o$samples.w)
  plot(x=truevals$W, y=w.hat, xlab="True W", ylab="Estimated W"); abline(0, 1, col='2')
  
}


# export tuning parameters to LaTeX
rows_all <- c()
for (p in priors){
  
  o <- get_output_general(outputs, tag=paste('priorcompare_', p, sep=""))
  rows_o <- summarize_mcmc_pscc(o, p)
  rows_all <- c(rows_o, rows_all)
  
}
df_mcmc <- ldply(rows_all, 'data.frame')
write_latex_table(df_mcmc, "latex_priorcompare_mcmc.txt", path=dst)


# summary table of alpha estimates
# model # parameter # estimate # bias
# (start with alpha case)
true_params <- load_params('true_params_priorcompare.json')
Alpha.case <- true_params$Alpha.case
Alpha.ctrl <- true_params$Alpha.ctrl
rows <- list()
counter <- 1
for (p in priors){
  o <- get_output_general(outputs, tag=paste('priorcompare_', p, sep=""))
  rows[[counter]] <- list(
    model=p,
    parameter="alpha (case)",
    estimate=round(mean(o$samples.alpha.ca), 3),
    true_value=Alpha.case,
    bias=round(mean(o$samples.alpha.ca), 3)-Alpha.case
  )
  counter <- counter + 1
}
for (p in priors){
  o <- get_output_general(outputs, tag=paste('priorcompare_', p, sep=""))
  rows[[counter]] <- list(
    model=p,
    parameter="alpha (control)",
    estimate=round(mean(o$samples.alpha.co), 3),
    true_value=Alpha.ctrl,
    bias=round(mean(o$samples.alpha.co), 3)-Alpha.ctrl
  )
  counter <- counter + 1
}
write_latex_table(ldply(rows, 'data.frame'), "latex_priorcompare_alpha_estimates.txt", path=dst)
