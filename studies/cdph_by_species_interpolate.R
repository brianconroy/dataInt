###########################
# Summarizes results of the 
# by species (Scurid)
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

##########
# baseline
##########

output_base <- load_output("output_cdph_baseline.json")
samples <- output_base$samples.w
r_pred <- caPr[[1]]
r_train <- aggregate(r_pred, fact=5)
bws <- c(0.08, 0.09)
samples_int <- interpolate_w_batched(samples, bws, r_train, r_pred, batch_size=500)
# save_output(samples_int, "cdph_baseline_interpolated_w.json")







# w.hat <- colMeans(samples_int)
# var.hat <- apply(samples_int, 2, var)
# plot(overlay(w.hat, r_pred))
# plot(overlay(var.hat, r_pred))
# plot(overlay(apply(output_base$samples.w, 2, var), r_train))



rs <- calc_significance_rasters_ds(rodents, output_base, caPr, threshold=0.05)
plot(rs$r_fracs)
plot(rs$r_inds)
plot(rs$r_inds_95)
plot_risk_overlay(rs$r_inds_95, rodents)
plot_risk_overlay(rs$r_inds, rodents)


groupings <- list(
  c('Pine Squirrel'),
  c('Chipmunk, YP'),
  c('Chipmunk, LP'),
  c('Chipmunk, S'),
  c('Chipmunk, M'),
  c('CA G Sq'),
  c('GM G Sq')
)

null_alphas <- c('GM G Sq', 'CA G Sq')

for (species in groupings){
  
  rodents_species <- rodents[rodents$Short_Name == species,]
  analysis_name <- gsub(',', '', gsub(' ', '_', paste('analysis', paste(species, collapse="_"), sep='_'), fixed=T))
  output_species <- load_output(paste("cdph_", analysis_name, ".json", sep=""))
  na <- F
  if (species %in% null_alphas){ na <- T} 
  rs <- calc_significance_rasters_ds(rodents_species, output_species, caPr, threshold=0.05, null_alphas=na)
  plot(rs$r_inds, main=species)
  plot_risk_overlay(rs$r_inds, rodents_species, main=species)
  plot(rs$r_inds_95, main=species)
  
}

caPr.disc <- aggregate(caPr, fact=agg_factor)
loc.disc <- caPr.disc[[1]]
all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
N <- n_values(caPr.disc[[1]])
print(N)
print(mean(area(caPr.disc[[1]])[]))
plot(caPr.disc)


##############
#### Risk maps
##############


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


#################
#### Save rasters
#################


for (s in 1:5){
  species <- groupings[[s]]
  r_risk_high <- calc_risk_cdph(species, rodents, caPr.disc, all_ids)
  writeRaster(
    r_risk_high, 
    filename=paste("/Users/brianconroy/Desktop/rasters/", gsub(" ", "_", species), sep=""), 
    format="GTiff")
}


#################
#### Rescale maps
#################


plot(calc_risk_cdph('CA G Sq', rodents, caPr.disc, all_ids))
plot(calc_risk_cdph('GM G Sq', rodents, caPr.disc, all_ids))

r_risk_yp <- calc_risk_cdph('Chipmunk, YP', rodents, caPr.disc, all_ids)
r_risk_lp <- calc_risk_cdph('Chipmunk, LP', rodents, caPr.disc, all_ids)
r_risk_s <- calc_risk_cdph('Chipmunk, S', rodents, caPr.disc, all_ids)

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


#####################
#### Summarize counts
#####################


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
