##############################
# New method to calculate high
# resolution significance maps
##############################


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


######################
#### Significance maps
######################


output <- load_output("output_cdph_baseline.json")
sig_maps <- calc_significance_rasters_ds2(output, caPr, threshold=0.05)
plot(sig_maps$r_fracs)
plot(sig_maps$r_inds)
plot(sig_maps$r_inds_95)


# compare
sig_maps2 <- calc_significance_rasters_ds(rodents, output, caPr, threshold=0.05)
plot(sig_maps2$r_fracs)
plot(sig_maps2$r_inds)
plot(sig_maps2$r_inds_95)
