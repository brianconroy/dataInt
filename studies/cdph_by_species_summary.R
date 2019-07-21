###########################
# Summarizes results of the 
# by species (scurid)
# cdph analysis
###########################

library(plyr)
library(grid)
library(mvtnorm)
library(ggplot2)
library(R.utils)
library(gridExtra)
sourceDirectory('Documents/research/dataInt/R/')


dst <- "/Users/brianconroy/Documents/research/project2/cdph_by_species/"
src <- "/Users/brianconroy/Documents/research/cdph/data/"
caPr <- load_prism_pcs2()
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")


#### Significance maps
output_base <- load_output("output_cdph_baseline.json")
rs <- calc_significance_rasters(rodents, output_base, caPr)
plot(rs$r_fracs)
plot(rs$r_inds)


species <- 'Chipmunk, M'
rodents_species <- rodents[rodents$Short_Name == species,]
analysis_name <- gsub(',', '', gsub(' ', '_', paste('analysis', paste(species, collapse="_"), sep='_'), fixed=T))
output_species <- load_output(paste("cdph_", analysis_name, ".json", sep=""))
rs <- calc_significance_rasters(rodents_species, output_species, caPr)
plot(rs$r_fracs)
plot(rs$r_inds)


groupings <- list(
  c('Pine Squirrel'),
  c('Chipmunk, YP'),
  c('Chipmunk, LP'),
  c('Chipmunk, S'),
  c('Chipmunk, M'),
  c('CA G Sq'),
  c('GM G Sq')
)


null_alphas <- c('GM G Sq', 'Pine Squirrel', 'CA G Sq')


for (species in groupings){
  
  rodents_species <- rodents[rodents$Short_Name == species,]
  analysis_name <- gsub(',', '', gsub(' ', '_', paste('analysis', paste(species, collapse="_"), sep='_'), fixed=T))
  output_species <- load_output(paste("cdph_", analysis_name, ".json", sep=""))
  na <- F
  if (species %in% null_alphas){ na <- T} 
  rs <- calc_significance_rasters(rodents_species, output_species, caPr, null_alphas=na)
  plot(rs$r_fracs, main=species)
  plot(rs$r_inds, main=species)
  
}

caPr.disc <- aggregate(caPr, fact=agg_factor)
loc.disc <- caPr.disc[[1]]
all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
N <- n_values(caPr.disc[[1]])
print(N)
print(mean(area(caPr.disc[[1]])[]))
plot(caPr.disc)


#### risk maps for single species views
for (species in groupings){
  
  results <- calc_risk_cdph(species, rodents, caPr.disc, all_ids, agg_factor=agg_factor, null_alphas=T)
  plot(results$r_risk_high, main=species)
  
  #### Figure: risk map with cases overlayed
  if (paste(species, collapse="") == 'all_but_ds'){
    all_species <- unique(rodents$Short_Name)
    species_group <- as.character(all_species[all_species != 'Pine Squirrel'])
    rodents_species <- rodents[rodents$Short_Name %in% species_group,]
  } else {
    rodents_species <- rodents[rodents$Short_Name %in% species,]
  }
  plot_risk_overlay(results, rodents_species)
  
  #### Figure: covariate contribution to log odds, random field contribution to log odds
  plot_cov_vs_w(results, caPr)
  
  #### posterior variance
  rodents_species <- rodents[rodents$Short_Name %in% species,]
  data <- assemble_data(rodents_species, loc.disc, caPr.disc)
  X_rodent <- load_x_standard2(as.logical(data$loc$status), agg_factor=agg_factor)
  analysis_name <- gsub(',', '', gsub(' ', '_', paste('analysis', paste(species, collapse="_"), sep='_'), fixed=T))
  output <- load_output(paste("cdph_", analysis_name, ".json", sep=""))
  risk_rodent_sep <- calc_posterior_risk(output, X_rodent)
  
  postvar_rodent_sep <- apply(risk_rodent_sep, 2, var)
  r <- caPr.disc[[1]]
  r[][!is.na(r[])] <- postvar_rodent_sep
  plot(r)
  
  xy <- data.frame(xyFromCell(r, 1:ncell(r)))
  v <- getValues(r)
  
  tps <- Tps(xy, v)
  p <- raster(caPr[[2]])
  p <- interpolate(p, tps)
  p <- mask(p, caPr[[1]])
  plot(p)
  
}


#### save rasters
for (s in 1:5){
  species <- groupings[[s]]
  r_risk_high <- calc_risk_cdph(species, rodents, caPr.disc, all_ids)
  writeRaster(
    r_risk_high, 
    filename=paste("/Users/brianconroy/Desktop/rasters/", gsub(" ", "_", species), sep=""), 
    format="GTiff")
}


#### risk maps for species comparisons, on the same scale
plot(calc_risk_cdph('CA G Sq', rodents, caPr.disc, all_ids))
plot(calc_risk_cdph('GM G Sq', rodents, caPr.disc, all_ids))


#### Chipmunks
r_risk_yp <- calc_risk_cdph('Chipmunk, YP', rodents, caPr.disc, all_ids)
r_risk_lp <- calc_risk_cdph('Chipmunk, LP', rodents, caPr.disc, all_ids)
r_risk_s <- calc_risk_cdph('Chipmunk, S', rodents, caPr.disc, all_ids)
# r_risk_m <- calc_risk_cdph('Chipmunk, M', rodents, caPr.disc, all_ids)

plot(r_risk_yp$r_risk_high)
plot(r_risk_s$r_risk_high)
plot(r_risk_lp$r_risk_high)

r_list <- list(r_risk_yp, r_risk_lp, r_risk_s)
r_list_new <- equalize_scales3(r_list)

plot(r_list_new[[1]]$r_risk_high)
rodents_species <- rodents[rodents$Short_Name == 'Chipmunk, YP',]
plot_risk_overlay(r_list_new[[1]], rodents_species)
plot_cov_vs_w(r_list_new[[1]], caPr)

plot(r_list_new[[2]]$r_risk_high)
rodents_species <- rodents[rodents$Short_Name == 'Chipmunk, LP',]
plot_risk_overlay(r_list_new[[2]], rodents_species)
plot_cov_vs_w(r_list_new[[2]], caPr)

plot(r_list_new[[3]]$r_risk_high)
rodents_species <- rodents[rodents$Short_Name == 'Chipmunk, S',]
plot_risk_overlay(r_list_new[[3]], rodents_species)
plot_cov_vs_w(r_list_new[[3]], caPr)

#### CA GS and T Amoenus together and separately
r_risk_yp <- calc_risk_cdph('Chipmunk, YP', rodents, caPr.disc, all_ids, agg_factor=agg_factor)
r_risk_cgs <- calc_risk_cdph('CA G Sq', rodents, caPr.disc, all_ids, agg_factor=agg_factor, null_alphas = T)
r_risk_yp_cgs <- calc_risk_cdph(c('CA G Sq', 'Chipmunk, YP'), rodents, caPr.disc, all_ids, agg_factor=agg_factor)
r_list2 <- list(r_risk_yp$r_risk_high, r_risk_cgs$r_risk_high, r_risk_yp_cgs$r_risk_high)
r_list_new2 <- equalize_scales2(r_list2)
plot(r_list_new2[[1]])
plot(r_list_new2[[2]])
plot(r_list_new2[[3]])

writeRaster(
  r_list_new2[[3]], 
  filename="/Users/brianconroy/Desktop/rasters/CGS_and_Chipmunk_YP", 
  format="GTiff")




#### Douglas squirrels together and separately
r_risk_ds <- calc_risk_cdph('Pine Squirrel', rodents, caPr.disc, all_ids)
r_risk_allbutds <- calc_risk_cdph('all_but_ds', rodents, caPr.disc, all_ids)
r_list3 <- list(r_risk_ds, r_risk_allbutds)
r_list_new3 <- equalize_scales3(r_list3)
plot(r_list_new3[[1]]$r_risk_high)
plot(r_list_new3[[2]]$r_risk_high)




plot(r_list_new3[[1]])
rodents_species <- rodents[rodents$Short_Name == 'Pine Squirrel',]
plot_risk_overlay(r_risk_ds, rodents_species)
plot_cov_vs_w(r_risk_ds, caPr)

writeRaster(
  r_list_new3[[2]], 
  filename="/Users/brianconroy/Desktop/rasters/All_but_Douglas_Squirrel", 
  format="GTiff")


#### summarize counts by species/group
rows <- list()
counter <- 1
for (species in groupings){
  print(species)
  if (paste(species, collapse="") == 'all_but_ds'){
    all_species <- unique(rodents$Short_Name)
    species_group <- as.character(all_species[all_species != 'Pine Squirrel'])
    rodents_species <- rodents[rodents$Short_Name %in% species_group,]
  } else {
    rodents_species <- rodents[rodents$Short_Name %in% species,]
  }
  ca <- table(rodents_species$Res)['POS']
  co <- table(rodents_species$Res)['NEG']
  rows[[counter]] <- list(
    Species=paste(species, collapse="|"),
    Cases=ca,
    Controls=co,
    Total=ca+co,
    Prevalence=round(ca/(ca+co), 3)
  )
  counter <- counter + 1
}
write_latex_table(ldply(rows, 'data.frame'), 'latex_species_counts.txt', dst)


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

# parameter estimates
params <- summarize_ps_params(output)
ite_latex_table(params, 'cdph_baseline_params.txt', dst)

# random field (downscaled)
w.hat <- colMeans(output$samples.w)
rw <- caPr.disc[[1]]
rw[][!is.na(rw[])] <- w.hat

xy <- data.frame(xyFromCell(rw, 1:ncell(rw)))
v <- getValues(rw)

tps <- Tps(xy, v)
p <- raster(caPr[[2]])
p <- interpolate(p, tps)
p <- mask(p, caPr[[1]])
w.hat_ds <- p[][!is.na(p[])]

par(mfrow=c(1,2))
plot(rw)
plot(p)

par(mfrow=c(1,2))
plot(rw)
plot(rw)
points(locs$coords, col=1, pch=16)

alpha.ca.hat <- mean(output$samples.alpha.ca)
alpha.co.hat <- mean(output$samples.alpha.co)
beta.ca.hat <- colMeans(output$samples.beta.ca)
beta.co.hat <- colMeans(output$samples.beta.co)

# risk map (low resolution)
X_low <- load_x_ca(factor=5)
lodds_low <- X_low %*% beta.ca.hat + alpha.ca.hat * w.hat - X_low %*% beta.co.hat - alpha.co.hat * w.hat
risk_low <- calc_risk(lodds_low)

r_lodds_low <- caPr.disc[[1]]
r_lodds_low[][!is.na(r_lodds_low[])] <- lodds_low
plot(r_lodds_low)

r_risk_low <- caPr.disc[[1]]
r_risk_low[][!is.na(r_risk_low[])] <- risk_low
plot(r_risk_low)

# risk map (downscaled)
X_high <- load_x_ca()
lodds_high <- X_high %*% beta.ca.hat + alpha.ca.hat * w.hat_ds - X_high %*% beta.co.hat - alpha.co.hat * w.hat_ds
risk_high <- calc_risk(lodds_high)

r_lodds_high <- caPr[[2]]
r_lodds_high[][!is.na(r_lodds_high[])] <- lodds_high
plot(r_lodds_high)

r_risk_high <- caPr[[2]]
r_risk_high[][!is.na(r_risk_high[])] <- risk_high

#pal <- colorRampPalette(c("blue","red"))
plot(r_risk_high)#, col=pal(15))

#################
# Compare to 
# reference model
#################

mod.ca <- glm(case.data$y ~ case.data$x.standardised - 1, family='poisson')
mod.co <- glm(ctrl.data$y ~ ctrl.data$x.standardised - 1, family='poisson')

beta.ca.hat_p <- unname(coefficients(mod.ca))
beta.co.hat_p <- unname(coefficients(mod.co))

lodds_low_p <- X_low %*% beta.ca.hat_p - X_low %*% beta.co.hat_p
r_lodds_low_p <- caPr.disc[[1]]
r_lodds_low_p[][!is.na(r_lodds_low_p[])] <- lodds_low_p
plot(r_lodds_low_p)

lodds_high_p <- X_high %*% beta.ca.hat_p - X_high %*% beta.co.hat_p
risk_high_p <- calc_risk(lodds_high_p)

r_lodds_high_p <- caPr[[2]]
r_lodds_high_p[][!is.na(r_lodds_high_p[])] <- lodds_high_p
plot(r_lodds_high_p)

r_risk_high_p <- caPr[[2]]
r_risk_high_p[][!is.na(r_risk_high_p[])] <- risk_high_p

par(mfrow=c(1,2))
pal <- colorRampPalette(c("blue","red"))
plot(r_risk_high_p, main='A)', 
     breaks=c(0.01, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21), col=pal(7))
plot(r_risk_high, main='B)',
     breaks=c(0.01, 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21), col=pal(7))

breaks=c(0, 5, 10, 15, 21, 25)

par(mfrow=c(1,1))
plot(y=r_risk_high_p[], 
     x=r_risk_high[], 
     xlab='Risk (Preferential Sampling)', 
     ylab='Risk (Poisson)'
     ); abline(0, 1, col=2)

# summary table of risk differences
r_diff <- abs(r_risk_high_p[][!is.na(r_risk_high_p[])] - r_risk_high[][!is.na(r_risk_high[])])
r_diff <- 100*r_diff/r_risk_high_p[][!is.na(r_risk_high_p[])]
r_diff_summary <- round(matrix(summary(r_diff), ncol=6), 4)
r_diff_summary <- data.frame(r_diff_summary)
names(r_diff_summary) <- c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Max")
print(round(r_diff_summary), 3)

# comparison table of parameter estimates
param_comp <- compare_params(beta.ca.hat_p, beta.co.hat_p, beta.ca.hat, beta.co.hat)
write_latex_table(param_comp, 'cdph_baseline_param_comparison.txt', dst)
