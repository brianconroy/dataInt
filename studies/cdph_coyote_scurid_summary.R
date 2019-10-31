###############################
# Summarizes simulation results
# of the joint (sciurid + coyote)
# cdph analysis
###############################

library(plyr)
library(grid)
library(mvtnorm)
library(ggplot2)
library(R.utils)
library(gridExtra)
library(MCMCpack)
sourceDirectory('Documents/research/dataInt/R/')


dst <- "/Users/brianconroy/Documents/research/project2/analysis/"
caPr <- load_prism_pcs2()
caPr.disc <- aggregate(caPr, fact=6)
N <- n_values(caPr.disc[[1]])
print(N)
print(mean(area(caPr.disc[[1]])[]))
plot(caPr.disc)
loc.disc <- caPr.disc[[1]]

analysis_name <- "cdph_coyote_scurid"
outputs <- load_sim_outputs(tag=analysis_name)
data <- load_output(paste(analysis_name, '_data.json', sep=''))

src <- "/Users/brianconroy/Documents/research/cdph/data/"

########
# Figure
########

# gridded summary of sample sites
loc.disc <- caPr.disc[[1]]
all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
rodent_coords_all <- cbind(matrix(rodents$Lon_Add_Fix), rodents$Lat_Add_Fix)
rodent_coords_distinct <- unique(rodent_coords_all)
rodent_cells_distinct <- cellFromXY(loc.disc, rodent_coords_distinct)
rodent_coord_to_cell <- cbind(rodent_coords_distinct, rodent_cells_distinct)
rodent_coord_to_cell <- data.frame(rodent_coord_to_cell)
names(rodent_coord_to_cell) <- c("x", "y", "cell")
rodent_cell_counts <- table(rodent_coord_to_cell$cell)
rodent_cell_counts <- data.frame(cbind(names(rodent_cell_counts), unname(rodent_cell_counts)))
names(rodent_cell_counts) <- c("cell", "count")
rodent_cells_all <- data.frame(all_ids)
names(rodent_cells_all) <- c("cell")
rodent_cells_all_counts <- merge(rodent_cells_all, rodent_cell_counts, all.x = T)
rodent_cells_all_counts$count <- as.numeric(as.character(rodent_cells_all_counts$count))
rodent_cells_all_counts[is.na(rodent_cells_all_counts$count),]$count <- 0
loc.disc[][!is.na(loc.disc[])] <- rodent_cells_all_counts$count
plot(loc.disc)

coyote_loc.disc <- caPr.disc[[1]]
all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
coyote_coords_all <- cbind(matrix(coyotes$Long_QC), coyotes$Lat_QC)
coyote_coords_distinct <- unique(coyote_coords_all)
coyote_cells_distinct <- cellFromXY(loc.disc, coyote_coords_distinct)
coyote_coord_to_cell <- cbind(coyote_coords_distinct, coyote_cells_distinct)
coyote_coord_to_cell <- data.frame(coyote_coord_to_cell)
names(coyote_coord_to_cell) <- c("x", "y", "cell")
coyote_cell_counts <- table(coyote_coord_to_cell$cell)
coyote_cell_counts <- data.frame(cbind(names(coyote_cell_counts), unname(coyote_cell_counts)))
names(coyote_cell_counts) <- c("cell", "count")
coyote_cells_all <- data.frame(all_ids)
names(coyote_cells_all) <- c("cell")
coyote_cells_all_counts <- merge(coyote_cells_all, coyote_cell_counts, all.x = T)
coyote_cells_all_counts$count <- as.numeric(as.character(coyote_cells_all_counts$count))
coyote_cells_all_counts[is.na(coyote_cells_all_counts$count),]$count <- 0
coyote_loc.disc[][!is.na(coyote_loc.disc[])] <- coyote_cells_all_counts$count
plot(coyote_loc.disc)

par(mfrow=c(1,2))
coyote_loc.disc[][1] <- max(loc.disc[], na.rm=T)
plot(loc.disc, main='A)')
plot(coyote_loc.disc, main='B)')

#########
# Figure: Coyote and Sciurid overlay
#########


us <- getData("GADM", country="USA", level=2)
ca <- us[us$NAME_1 == 'California',]
plot(ca)
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")
coyotes <- read.csv(paste(src, "CDPH_coyote_recoded_full.csv", sep=""), header=T, sep=",")
c_coords <- cbind(coyotes$Long_QC, coyotes$Lat_QC)
r_coords <- cbind(rodents$Lon_Add_Fix, rodents$Lat_Add_Fix)
points(c_coords, col=rgb(0,0,1,0.25), pch=16, cex=0.75)
points(r_coords, col=rgb(1,0,0,0.25), pch=16, cex=0.75)


#########
# Figure: PRISM pcs
#########


par(mfrow=c(1,2))
plot(caPr[[1]], main='A)')
plot(caPr[[2]], main='B)')


###############
# patch output
###############


# cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
# coords <- xyFromCell(caPr.disc, cell=cells.all)
# D <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
# data_input <- load_mvgp_data(paste(analysis_name, '_data.json', sep=''))
# 
# prior_theta <- c(3, 2)
# prior_alpha_ca_mean <- c(0, 0)
# prior_alpha_co_mean <- c(0, 0)
# prior_alpha_ca_var <- c(4, 4)
# prior_alpha_co_var <- c(4, 4)
# 
# o_mvgp <- get_output_general(outputs, tag=paste(analysis_name, 'mvgp', sep="_"))
# 
# o_mvgp <- continueMCMC_mvgp(data_input, D, o_mvgp, n.sample=2000)
# 
# o_mvgp <- burnin_mvgp(o_mvgp, n.burn=400)
# 
# plot(o_mvgp$samples.theta, type='l')
# plot(o_mvgp$samples.t[,1], type='l')
# plot(o_mvgp$samples.t[,2], type='l')
# plot(o_mvgp$samples.t[,4], type='l')
# 
# plot(o_mvgp$samples.alpha.ca[1,,], type='l')
# plot(o_mvgp$samples.alpha.co[1,,], type='l')
# plot(o_mvgp$samples.alpha.ca[2,,], type='l')
# plot(o_mvgp$samples.alpha.co[2,,], type='l')
# 
# save_output(o_mvgp, paste(o_mvgp$description, ".json", sep=""))


#########
# Figure: risk maps
#########


o_mvgp <- get_output_general(outputs, tag=paste(analysis_name, 'mvgp', sep="_"))
lodds_rodent_mvgp <- calc_lodds_mvgp2(o_mvgp, data, 1, agg_factor=6)
lodds_coyote_mvgp <- calc_lodds_mvgp2(o_mvgp, data, 2, agg_factor=6)
risk_rodent_mvgp <- calc_risk(lodds_rodent_mvgp)
risk_coyote_mvgp <- calc_risk(lodds_coyote_mvgp)

#### downscale
w.hat <- colMeans(o_mvgp$samples.w)
w.hat1 <- w.hat[seq(1, length(w.hat), by=2)]
w.hat2 <- w.hat[seq(2, length(w.hat), by=2)]
ds1 <- downscale(w.hat1, caPr.disc, caPr)
ds2 <- downscale(w.hat2, caPr.disc, caPr)
w.hat_ds1 <- ds1$w.hat_ds
w.hat_ds2 <- ds2$w.hat_ds

alpha.ca.hat1 <- mean(o_mvgp$samples.alpha.ca[1,,])
alpha.ca.hat2 <- mean(o_mvgp$samples.alpha.ca[2,,])
alpha.co.hat1 <- mean(o_mvgp$samples.alpha.co[1,,])
alpha.co.hat2 <- mean(o_mvgp$samples.alpha.co[2,,])
beta.ca.hat1 <- colMeans(o_mvgp$samples.beta.ca[1,,])
beta.ca.hat2 <- colMeans(o_mvgp$samples.beta.ca[2,,])
beta.co.hat1 <- colMeans(o_mvgp$samples.beta.co[1,,])
beta.co.hat2 <- colMeans(o_mvgp$samples.beta.co[2,,])

lodds_rodent_ds <- calc_lodds_ds(alpha.ca.hat1, alpha.co.hat1, beta.ca.hat1, beta.co.hat1, w.hat_ds1)
lodds_coyote_ds <- calc_lodds_ds(alpha.ca.hat2, alpha.co.hat2, beta.ca.hat2, beta.co.hat2, w.hat_ds2)

risk_rodent_ds <- calc_risk_ds(alpha.ca.hat1, alpha.co.hat1, beta.ca.hat1, beta.co.hat1, w.hat_ds1)
risk_coyote_ds <- calc_risk_ds(alpha.ca.hat2, alpha.co.hat2, beta.ca.hat2, beta.co.hat2, w.hat_ds2)

r_r_ds <- caPr[[1]]
r_r_ds[][!is.na(r_r_ds[])] <- risk_rodent_ds
r_c_ds <- caPr[[1]]
r_c_ds[][!is.na(r_c_ds[])] <- risk_coyote_ds


#########
# Figure:  coyote MVGP risk map
#########


plot(r_c_ds)


#########
# Figure: sciurid MVGP risk map
#########


plot(r_r_ds)


# risk maps with cases overlayed
src <- "/Users/brianconroy/Documents/research/cdph/data/"
rodents <- read.csv(paste(src, "CDPH_scurid_updated_full.csv", sep=""), header=T, sep=",")
coyotes <- read.csv(paste(src, "CDPH_coyote_recoded_full.csv", sep=""), header=T, sep=",")


#########
# Figure: coyote overlay
#########


c_cases <- coyotes[coyotes$Res == 'POS',]
c_coords <- cbind(c_cases$Long_QC, c_cases$Lat_QC)
c_ctrls <- coyotes[coyotes$Res == 'NEG',]
c_coords_ctrl <- cbind(c_ctrls$Long_QC, c_ctrls$Lat_QC)
plot(r_c_ds)
points(c_coords_ctrl, col=rgb(0,0,1,0.05), pch=16, cex=0.65)
points(c_coords, col=rgb(1,0,0,0.2), pch=16, cex=0.65)


#########
# Figure: rodent overlay
#########


r_cases <- rodents[rodents$Res == 'POS',]
r_coords <- cbind(r_cases$Lon_Add_Fix, r_cases$Lat_Add_Fix)
r_ctrls <- rodents[rodents$Res == 'NEG',]
r_coords_ctrl <- cbind(r_ctrls$Lon_Add_Fix, r_ctrls$Lat_Add_Fix)
plot(r_r_ds)
points(r_coords_ctrl, col=rgb(0,0,1,0.05), pch=16, cex=0.65)
points(r_coords, col=rgb(1,0,0,0.2), pch=16, cex=0.65)


#########
# Figure: side by side of risk map
#########


#### with case and control gridded counts
rodent_case.disc <- caPr.disc[[1]]
all_ids <- c(1:length(rodent_case.disc[]))[!is.na(rodent_case.disc[])]
rodent_coords_case <- cbind(rodents[rodents$Res == 'POS',]$Lon_Add_Fix, rodents[rodents$Res == 'POS',]$Lat_Add_Fix)
rodent_cells_case <- cellFromXY(rodent_case.disc, rodent_coords_case)
rodent_coord_to_cell_case <- cbind(rodent_coords_case, rodent_cells_case)
rodent_coord_to_cell_case <- data.frame(rodent_coord_to_cell_case)
names(rodent_coord_to_cell_case) <- c("x", "y", "cell")

rodent_cell_counts_case <- table(rodent_coord_to_cell_case$cell)
rodent_cell_counts_case <- data.frame(cbind(names(rodent_cell_counts_case), unname(rodent_cell_counts_case)))
names(rodent_cell_counts_case) <- c("cell", "count")

cells_all <- data.frame(all_ids)
names(cells_all) <- c("cell")
rodent_cells_all_counts_case <- merge(cells_all, rodent_cell_counts_case, all.x = T)
rodent_cells_all_counts_case$count <- as.numeric(as.character(rodent_cells_all_counts_case$count))
rodent_cells_all_counts_case[is.na(rodent_cells_all_counts_case$count),]$count <- 0
rodent_case.disc[][!is.na(rodent_case.disc[])] <- rodent_cells_all_counts_case$count
plot(rodent_case.disc)

rodent_ctrl.disc <- caPr.disc[[1]]
all_ids <- c(1:length(rodent_ctrl.disc[]))[!is.na(rodent_ctrl.disc[])]
rodent_coords_ctrl <- cbind(rodents[rodents$Res == 'NEG',]$Lon_Add_Fix, rodents[rodents$Res == 'NEG',]$Lat_Add_Fix)
rodent_cells_ctrl <- cellFromXY(rodent_ctrl.disc, rodent_coords_ctrl)
rodent_coord_to_cell_ctrl <- cbind(rodent_coords_ctrl, rodent_cells_ctrl)
rodent_coord_to_cell_ctrl <- data.frame(rodent_coord_to_cell_ctrl)
names(rodent_coord_to_cell_ctrl) <- c("x", "y", "cell")

rodent_cell_counts_ctrl <- table(rodent_coord_to_cell_ctrl$cell)
rodent_cell_counts_ctrl <- data.frame(cbind(names(rodent_cell_counts_ctrl), unname(rodent_cell_counts_ctrl)))
names(rodent_cell_counts_ctrl) <- c("cell", "count")

cells_all <- data.frame(all_ids)
names(cells_all) <- c("cell")
rodent_cells_all_counts_ctrl <- merge(cells_all, rodent_cell_counts_ctrl, all.x = T)
rodent_cells_all_counts_ctrl$count <- as.numeric(as.character(rodent_cells_all_counts_ctrl$count))
rodent_cells_all_counts_ctrl[is.na(rodent_cells_all_counts_ctrl$count),]$count <- 0
rodent_ctrl.disc[][!is.na(rodent_ctrl.disc[])] <- rodent_cells_all_counts_ctrl$count
plot(rodent_ctrl.disc)

par(mfrow=c(2,2))
plot(rodent_case.disc, main='A)')
plot(rodent_ctrl.disc, main='B)')
plot(r_r_ds, main='C)')

#########
# Coyotes


coyote_case.disc <- caPr.disc[[1]]
all_ids <- c(1:length(coyote_case.disc[]))[!is.na(coyote_case.disc[])]
coyote_coords_case <- cbind(coyotes[coyotes$Res == 'POS',]$Long_QC, coyotes[coyotes$Res == 'POS',]$Lat_QC)
coyote_cells_case <- cellFromXY(coyote_case.disc, coyote_coords_case)
coyote_coord_to_cell_case <- cbind(coyote_coords_case, coyote_cells_case)
coyote_coord_to_cell_case <- data.frame(coyote_coord_to_cell_case)
names(coyote_coord_to_cell_case) <- c("x", "y", "cell")

coyote_cell_counts_case <- table(coyote_coord_to_cell_case$cell)
coyote_cell_counts_case <- data.frame(cbind(names(coyote_cell_counts_case), unname(coyote_cell_counts_case)))
names(coyote_cell_counts_case) <- c("cell", "count")

cells_all <- data.frame(all_ids)
names(cells_all) <- c("cell")
coyote_cells_all_counts_case <- merge(cells_all, coyote_cell_counts_case, all.x = T)
coyote_cells_all_counts_case$count <- as.numeric(as.character(coyote_cells_all_counts_case$count))
coyote_cells_all_counts_case[is.na(coyote_cells_all_counts_case$count),]$count <- 0
coyote_case.disc[][!is.na(coyote_case.disc[])] <- coyote_cells_all_counts_case$count
plot(coyote_case.disc)

coyote_ctrl.disc <- caPr.disc[[1]]
all_ids <- c(1:length(coyote_ctrl.disc[]))[!is.na(coyote_ctrl.disc[])]
coyote_coords_ctrl <- cbind(coyotes[coyotes$Res == 'NEG',]$Long_QC, coyotes[coyotes$Res == 'NEG',]$Lat_QC)
coyote_cells_ctrl <- cellFromXY(coyote_ctrl.disc, coyote_coords_ctrl)
coyote_coord_to_cell_ctrl <- cbind(coyote_coords_ctrl, coyote_cells_ctrl)
coyote_coord_to_cell_ctrl <- data.frame(coyote_coord_to_cell_ctrl)
names(coyote_coord_to_cell_ctrl) <- c("x", "y", "cell")

coyote_cell_counts_ctrl <- table(coyote_coord_to_cell_ctrl$cell)
coyote_cell_counts_ctrl <- data.frame(cbind(names(coyote_cell_counts_ctrl), unname(coyote_cell_counts_ctrl)))
names(coyote_cell_counts_ctrl) <- c("cell", "count")

cells_all <- data.frame(all_ids)
names(cells_all) <- c("cell")
coyote_cells_all_counts_ctrl <- merge(cells_all, coyote_cell_counts_ctrl, all.x = T)
coyote_cells_all_counts_ctrl$count <- as.numeric(as.character(coyote_cells_all_counts_ctrl$count))
coyote_cells_all_counts_ctrl[is.na(coyote_cells_all_counts_ctrl$count),]$count <- 0
coyote_ctrl.disc[][!is.na(coyote_ctrl.disc[])] <- coyote_cells_all_counts_ctrl$count
plot(coyote_ctrl.disc)

par(mfrow=c(2,2))
plot(coyote_case.disc, main='A)')
plot(coyote_ctrl.disc, main='B)')
plot(r_c_ds, main='C)')




# covariate contribution to log odds, random field contribution to log odds
#########
# Figure:  rodents
#########


X_high <- load_x_ca2()
cov_rodent <- X_high %*% beta.ca.hat1 - X_high %*% beta.co.hat1
w_rodent <- alpha.ca.hat1 * w.hat_ds1 - alpha.co.hat1 * w.hat_ds1
r_cov_rodent <- caPr[[1]]
r_cov_rodent[][!is.na(r_cov_rodent[])] <- cov_rodent
r_w_rodent <- caPr[[1]]
r_w_rodent[][!is.na(r_w_rodent[])] <- w_rodent
par(mfrow=c(1,2))
plot(r_cov_rodent, main='A)')
pal <- colorRampPalette(c("blue","red"))
plot(r_w_rodent, main='B)', col=pal(20))
summary(w_rodent[])


#########
# Figure:  coyotes
#########


X_high <- load_x_ca2()
cov_coyote <- X_high %*% beta.ca.hat2 - X_high %*% beta.co.hat2
w_coyote <- alpha.ca.hat2 * w.hat_ds2 - alpha.co.hat2 * w.hat_ds2
r_cov_coyote <- caPr[[1]]
r_cov_coyote[][!is.na(r_cov_coyote[])] <- cov_coyote
r_w_coyote <- caPr[[1]]
r_w_coyote[][!is.na(r_w_coyote[])] <- w_coyote
par(mfrow=c(1,2))
plot(r_cov_coyote, main='A)')
pal <- colorRampPalette(c("blue","red"))
plot(r_w_coyote, main='B)', col=pal(20))


#########################
# Table: Alpha estimates
#########################


a_params <- list()
a_params[[1]] <- list(
  Parameter="Alpha (Case, Rodents)",
  Estimate=round(mean(o_mvgp$samples.alpha.ca[1,,1]), 3),
  Variance=round(var(o_mvgp$samples.alpha.ca[1,,1]), 3)
)
a_params[[2]] <- list(
  Parameter="Alpha (Control, Rodents)",
  Estimate=round(mean(o_mvgp$samples.alpha.co[1,,1]), 3),
  Variance=round(var(o_mvgp$samples.alpha.co[1,,1]), 4)
)
a_params[[3]] <- list(
  Parameter="Alpha (Case, Coyotes)",
  Estimate=round(mean(o_mvgp$samples.alpha.ca[2,,1]), 3),
  Variance=round(var(o_mvgp$samples.alpha.ca[2,,1]), 3)
)
a_params[[4]] <- list(
  Parameter="Alpha (Control, Coyotes)",
  Estimate=round(mean(o_mvgp$samples.alpha.co[2,,1]), 3),
  Variance=round(var(o_mvgp$samples.alpha.co[2,,1]), 4)
)
write_latex_table(ldply(a_params, 'data.frame'), "cdph_mvgp_alpha_params.txt", path=dst)


##########################
# Table: T matrix estimate
##########################


t_params <- list()
t_params[[1]] <- list(
  Parameter="T (1,1)",
  Estimate=round(mean(o_mvgp$samples.t[,1]), 3),
  Variance=round(var(o_mvgp$samples.t[,1]), 3)
)
t_params[[2]] <- list(
  Parameter="T (1,2)",
  Estimate=round(mean(o_mvgp$samples.t[,2]), 3),
  Variance=round(var(o_mvgp$samples.t[,2]), 3)
)
t_params[[3]] <- list(
  Parameter="T (2,2)",
  Estimate=round(mean(o_mvgp$samples.t[,4]), 3),
  Variance=round(var(o_mvgp$samples.t[,4]), 3)
)
write_latex_table(ldply(t_params, 'data.frame'), "cdph_mvgp_t_params.txt", path=dst)


#########
# Figure: compare log odds (not downscaled)
#########


lodds_rodent_mvgp <- calc_lodds_mvgp2(o_mvgp, data, 1, agg_factor=6)
lodds_coyote_mvgp <- calc_lodds_mvgp2(o_mvgp, data, 2, agg_factor=6)

o_coyote_sep <- get_output_general(outputs, tag=paste(analysis_name, 'coyote', sep="_"))
lodds_coyote_sep <- calc_log_odds_species2(o_coyote_sep, data, species=2, agg_factor=6)

o_rodent_sep <- get_output_general(outputs, tag=paste(analysis_name, 'rodent', sep="_"))
lodds_rodent_sep <- calc_log_odds_species2(o_rodent_sep, data, species=1, agg_factor=6)

par(mfrow=c(1,2))
plot(y=lodds_rodent_mvgp, x=lodds_rodent_sep, ylab='Log Odds (MVGP)', xlab='Log Odds (separate)', main='A)'); abline(0,1,col=2)
plot(y=lodds_coyote_mvgp, x=lodds_coyote_sep, ylab='Log Odds (MVGP)', xlab='Log Odds (separate)', main='B)'); abline(0,1,col=2)

summary(lodds_coyote_mvgp-lodds_coyote_sep)
summary(lodds_rodent_mvgp-lodds_rodent_sep)


#########
# Figure: Compare Risk (not downscaled)
#########


lodds_rodent_mvgp <- calc_lodds_mvgp2(o_mvgp, data, 1, agg_factor=6)
lodds_coyote_mvgp <- calc_lodds_mvgp2(o_mvgp, data, 2, agg_factor=6)
risk_rodent_mvgp <- calc_risk(lodds_rodent_mvgp)
risk_coyote_mvgp <- calc_risk(lodds_coyote_mvgp)

lodds_coyote_sep <- calc_log_odds_species2(o_coyote_sep, data, species=2, agg_factor=6)
lodds_rodent_sep <- calc_log_odds_species2(o_rodent_sep, data, species=1, agg_factor=6)
risk_coyote_sep <- calc_risk(lodds_coyote_sep)
risk_rodent_sep <- calc_risk(lodds_rodent_sep)

par(mfrow=c(1,2))
plot(y=risk_rodent_mvgp, x=risk_rodent_sep, ylab='Risk (MVGP)', xlab='Risk (Separate)', main='A)'); abline(0,1,col=2)
plot(y=risk_coyote_mvgp, x=risk_coyote_sep, ylab='Risk (MVGP)', xlab='Risk (Separate)', main='B)'); abline(0,1,col=2)

summary(risk_coyote_mvgp-risk_coyote_sep)
summary(risk_rodent_mvgp-risk_rodent_sep)


#####################################
# Figure: Posterior log odds variance
#####################################


X_rodent <- load_x_standard2(as.logical(data$locs$status[[1]]), agg_factor=6)
X_coyote <- load_x_standard2(as.logical(data$locs$status[[2]]), agg_factor=6)
lodds_rodent_sep <- calc_posterior_lodds(o_rodent_sep, X_rodent)
lodds_coyote_sep <- calc_posterior_lodds(o_coyote_sep, X_coyote)
lodds_rodent_mvgp <- calc_posterior_lodds_multi(o_mvgp, X_rodent, species=1)
lodds_coyote_mvgp <- calc_posterior_lodds_multi(o_mvgp, X_coyote, species=2)

postvar_lodds_rodent_sep <- apply(lodds_rodent_sep, 2, var)
postvar_lodds_coyote_sep <- apply(lodds_coyote_sep, 2, var)
postvar_lodds_rodent_mvgp <- apply(lodds_rodent_mvgp, 2, var)
postvar_lodds_coyote_mvgp <- apply(lodds_coyote_mvgp, 2, var)

par(mfrow=c(1,2))
plot(x=postvar_lodds_coyote_sep, y=postvar_lodds_coyote_mvgp,
     xlab='Posterior Log Odds Variance (Separate)',
     ylab='Posterior Log Odds Variance (MVGP)',
     main='A)'
); abline(0, 1, col=2)
plot(x=postvar_lodds_rodent_sep, y=postvar_lodds_rodent_mvgp, 
     xlab='Posterior Log Odds Variance (Separate)',
     ylab='Posterior Log Odds Variance (MVGP)',
     main='B)'
); abline(0, 1, col=2)

# number of cells where mvgp beats separate
sum(postvar_lodds_rodent_sep > postvar_lodds_rodent_mvgp)/length(postvar_lodds_rodent_mvgp)
sum(postvar_lodds_coyote_sep > postvar_lodds_coyote_mvgp)/length(postvar_lodds_coyote_mvgp)

# mean, median, min, max
summary(postvar_lodds_rodent_sep)
summary(postvar_lodds_rodent_mvgp)

summary(postvar_lodds_coyote_sep)
summary(postvar_lodds_coyote_mvgp)


#########################
# Posterior risk variance
#########################


X_rodent <- load_x_standard2(as.logical(data$locs$status[[1]]), agg_factor=6)
X_coyote <- load_x_standard2(as.logical(data$locs$status[[2]]), agg_factor=6)
risk_rodent_sep <- calc_posterior_risk(o_rodent_sep, X_rodent)
risk_coyote_sep <- calc_posterior_risk(o_coyote_sep, X_coyote)
risk_rodent_mvgp <- calc_posterior_risk_multi(o_mvgp, X_rodent, species=1)
risk_coyote_mvgp <- calc_posterior_risk_multi(o_mvgp, X_coyote, species=2)

postvar_rodent_sep <- apply(risk_rodent_sep, 2, var)
postvar_coyote_sep <- apply(risk_coyote_sep, 2, var)
postvar_rodent_mvgp <- apply(risk_rodent_mvgp, 2, var)
postvar_coyote_mvgp <- apply(risk_coyote_mvgp, 2, var)

par(mfrow=c(1,2))
plot(x=postvar_coyote_sep, y=postvar_coyote_mvgp,
     xlab='Posterior Risk Variance (Separate)',
     ylab='Posterior Risk Variance (MVGP)',
     main='A)'
); abline(0, 1, col=2)
plot(x=postvar_rodent_sep, y=postvar_rodent_mvgp, 
     xlab='Posterior Risk Variance (Separate)',
     ylab='Posterior Risk Variance (MVGP)',
     main='B)'
); abline(0, 1, col=2)

summary(postvar_rodent_sep)
summary(postvar_rodent_mvgp)

summary(postvar_coyote_sep)
summary(postvar_coyote_mvgp)

r1 <- caPr.disc[[1]]
r1[][!is.na(r1[])] <- postvar_rodent_mvgp
plot(r1)

r2 <- caPr.disc[[1]]
r2[][!is.na(r2[])] <- postvar_coyote_mvgp
plot(r2)


#### compare to separate model downscaled risk maps
risk_rodent_sep <- calc_risk(lodds_rodent_sep)
risk_coyote_sep <- calc_risk(lodds_coyote_sep)
w.hat_r_sep <- colMeans(o_rodent_sep$samples.w)
w.hat_c_sep <- colMeans(o_coyote_sep$samples.w)

ds_r_sep <- downscale(w.hat_r_sep, caPr.disc, caPr)
ds_c_sep <- downscale(w.hat_c_sep, caPr.disc, caPr)
w.hat_r_sep_ds <- ds_r_sep$w.hat_ds
w.hat_c_sep_ds <- ds_c_sep$w.hat_ds

alpha.ca.hat_r <- mean(o_rodent_sep$samples.alpha.ca)
alpha.co.hat_r <- mean(o_rodent_sep$samples.alpha.co)
beta.ca.hat_r <- colMeans(o_rodent_sep$samples.beta.ca)
beta.co.hat_r <- colMeans(o_rodent_sep$samples.beta.co)

alpha.ca.hat_c <- mean(o_coyote_sep$samples.alpha.ca)
alpha.co.hat_c <- mean(o_coyote_sep$samples.alpha.co)
beta.ca.hat_c <- colMeans(o_coyote_sep$samples.beta.ca)
beta.co.hat_c <- colMeans(o_coyote_sep$samples.beta.co)

lodds_rodent_ds_sep <- calc_lodds_ds(alpha.ca.hat_r, alpha.co.hat_r, beta.ca.hat_r, beta.co.hat_r, w.hat_r_sep_ds)
lodds_coyote_ds_sep <- calc_lodds_ds(alpha.ca.hat_c, alpha.co.hat_c, beta.ca.hat_c, beta.co.hat_c, w.hat_c_sep_ds)

risk_rodent_ds_sep <- calc_risk_ds(alpha.ca.hat_r, alpha.co.hat_r, beta.ca.hat_r, beta.co.hat_r, w.hat_r_sep_ds)
risk_coyote_ds_sep <- calc_risk_ds(alpha.ca.hat_c, alpha.co.hat_c, beta.ca.hat_c, beta.co.hat_c, w.hat_c_sep_ds)

r_r_ds_sep <- caPr[[1]]
r_r_ds_sep[][!is.na(r_r_ds_sep[])] <- risk_rodent_ds_sep
r_c_ds_sep <- caPr[[1]]
r_c_ds_sep[][!is.na(r_c_ds_sep[])] <- risk_coyote_ds_sep




#########
# Figure: ds lodds comparisons
#########


par(mfrow=c(1,2))
plot(y=lodds_rodent_ds, x=lodds_rodent_ds_sep, ylab='Log Odds (MVGP)', xlab='Log Odds (separate)', main='A)'); abline(0,1,col=2)
plot(y=lodds_coyote_ds, x=lodds_coyote_ds_sep, ylab='Log Odds (MVGP)', xlab='Log Odds (separate)', main='B)'); abline(0,1,col=2)


summary(lodds_coyote_ds - lodds_coyote_ds_sep)
summary(lodds_rodent_ds - lodds_rodent_ds_sep)


#########
# Figure: ds risk comparisons
#########


par(mfrow=c(1,2))
plot(y=risk_rodent_ds, x=risk_rodent_ds_sep, ylab='Risk (MVGP)', xlab='Risk (separate)', main='A)'); abline(0,1,col=2)
plot(y=risk_coyote_ds, x=risk_coyote_ds_sep, ylab='Risk (MVGP)', xlab='Risk (separate)', main='B)'); abline(0,1,col=2)


summary(risk_coyote_ds - risk_coyote_ds_sep)
summary(risk_rodent_ds - risk_rodent_ds_sep)

summary(100*(risk_rodent_ds - risk_rodent_ds_sep)/risk_rodent_ds)


#### Coyote risk map comparison: mvgp vs ps
c_eq <- equalize_scales(r_c_ds, r_c_ds_sep)
par(mfrow=c(1,2))
plot(c_eq[[1]], main='A)')
plot(c_eq[[2]], main='B)')


#### Rodent risk map comparison: mvgp vs ps
r_eq <- equalize_scales(r_r_ds, r_r_ds_sep)
par(mfrow=c(1,2))
plot(r_eq[[1]], main='A)')
plot(r_eq[[2]], main='B)')


#### T Matrix Estimate
#### parameter | estimate | posterior variance
par(mfrow=c(2,2))
plot(o_mvgp$samples.t[,1], type='l')
plot(o_mvgp$samples.t[,2], type='l')
plot(o_mvgp$samples.t[,4], type='l')

t_params <- list()
t_params[[1]] <- list(
  Parameter="T (1,1)",
  Estimate=round(mean(o_mvgp$samples.t[,1]), 3),
  Variance=round(var(o_mvgp$samples.t[,1]), 3)
)
t_params[[2]] <- list(
  Parameter="T (1,2)",
  Estimate=round(mean(o_mvgp$samples.t[,2]), 3),
  Variance=round(var(o_mvgp$samples.t[,2]), 3)
)
t_params[[3]] <- list(
  Parameter="T (2,2)",
  Estimate=round(mean(o_mvgp$samples.t[,4]), 3),
  Variance=round(var(o_mvgp$samples.t[,4]), 3)
)
write_latex_table(ldply(t_params, 'data.frame'), "cdph_mvgp_t_params.txt", path=dst)
