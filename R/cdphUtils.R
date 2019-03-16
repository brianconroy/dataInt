

downscale <- function(w.est, caPr.disc, caPr){
  
  rw <- caPr.disc[[1]]
  rw[][!is.na(rw[])] <- w.est
  xy <- data.frame(xyFromCell(rw, 1:ncell(rw)))
  v <- getValues(rw)
  tps <- Tps(xy, v)
  p <- raster(caPr[[2]])
  p <- interpolate(p, tps)
  p <- mask(p, caPr[[1]])
  w.hat_ds <- p[][!is.na(p[])]
  return(list(p=p, w.hat_ds=w.hat_ds))
  
}


assemble_data <- function(rodents, loc.disc, caPr.disc){
  
  coords_all <- cbind(matrix(rodents$Lon_Add_Fix), rodents$Lat_Add_Fix)
  loc.disc <- caPr.disc[[1]]
  cells_all <- cellFromXY(loc.disc, coords_all)
  cells_all <- cells_all[!is.na(cells_all[])]
  cells_obs <- sort(unique(cells_all))
  
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
  if (sum(is.na(counts_all$count_pos))){
    counts_all[is.na(counts_all$count_pos),]$count_pos <- 0
  } 
  if (sum(is.na(counts_all$count_neg))){
    counts_all[is.na(counts_all$count_neg),]$count_neg <- 0
  }
  counts_all <- counts_all[with(counts_all, order(cell)),]
  
  # location data
  all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
  locs <- list(
    cells=cells_obs,
    status=1 * c(all_ids %in% cells_obs),  
    coords=xyFromCell(loc.disc, cells_obs)
  )
  locs$ids <- c(1:length(all_ids))[as.logical(locs$status)]
  
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
  
  data <- list(loc=locs, case.data=case.data, ctrl.data=ctrl.data)
  
  return(data)
  
}


assemble_data_coyotes <- function(coyotes, caPr.disc){
  
  coords_all <- cbind(matrix(coyotes$Long_QC), coyotes$Lat_QC)
  loc.disc <- caPr.disc[[1]]
  cells_all <- cellFromXY(loc.disc, coords_all)
  cells_all <- cells_all[!is.na(cells_all[])]
  cells_obs <- sort(unique(cells_all))
  
  # positive counts at each cell
  coyotes_pos <- coyotes[coyotes$Res == 'POS',]
  coords_pos <- cbind(matrix(coyotes_pos$Long_QC), coyotes_pos$Lat_QC)
  cells_pos <- cellFromXY(loc.disc, coords_pos)
  counts_pos <- data.frame(table(cells_pos))
  names(counts_pos) <- c('cell', 'count_pos')
  
  # negative counts at each cell
  coyotes_neg <- coyotes[coyotes$Res == 'NEG',]
  coords_neg <- cbind(matrix(coyotes_neg$Long_QC), coyotes_neg$Lat_QC)
  cells_neg <- cellFromXY(loc.disc, coords_neg)
  counts_neg <- data.frame(table(cells_neg))
  names(counts_neg) <- c('cell', 'count_neg')
  
  # combine counts
  counts_all <- merge(counts_pos, counts_neg, by='cell', all=T)
  counts_all$cell <- as.numeric(as.character(counts_all$cell))
  if (sum(is.na(counts_all$count_pos))){
    counts_all[is.na(counts_all$count_pos),]$count_pos <- 0
  } 
  if (sum(is.na(counts_all$count_neg))){
    counts_all[is.na(counts_all$count_neg),]$count_neg <- 0
  }
  counts_all <- counts_all[with(counts_all, order(cell)),]
  
  # location data
  all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
  locs <- list(
    cells=cells_obs,
    status=1 * c(all_ids %in% cells_obs),  
    coords=xyFromCell(loc.disc, cells_obs)
  )
  locs$ids <- c(1:length(all_ids))[as.logical(locs$status)]
  
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
  
  data <- list(loc=locs, case.data=case.data, ctrl.data=ctrl.data)
  
  return(data)
  
}


summarize_ps_params <- function(output){
  
  params <- list()
  params[[1]] <- list(
    Name='Beta 0 (case)',
    Posterior_Mean=round(mean(output$samples.beta.ca[,1]), 3),
    Posterior_Median=round(median(output$samples.beta.ca[,1]), 3),
    Posterior_Variance=round(var(output$samples.beta.ca[,1]), 3)
  )
  params[[2]] <- list(
    Name='Beta 1 (case)',
    Posterior_Mean=round(mean(output$samples.beta.ca[,2]), 3),
    Posterior_Median=round(median(output$samples.beta.ca[,2]), 3),
    Posterior_Variance=round(var(output$samples.beta.ca[,2]), 3)
  )
  params[[3]] <- list(
    Name='Beta 2 (case)',
    Posterior_Mean=round(mean(output$samples.beta.ca[,3]), 3),
    Posterior_Median=round(median(output$samples.beta.ca[,3]), 3),
    Posterior_Variance=round(var(output$samples.beta.ca[,3]), 3)
  )
  params[[4]] <- list(
    Name='Beta 0 (control)',
    Posterior_Mean=round(mean(output$samples.beta.co[,1]), 3),
    Posterior_Median=round(median(output$samples.beta.co[,1]), 3),
    Posterior_Variance=round(var(output$samples.beta.co[,1]), 3)
  )
  params[[5]] <- list(
    Name='Beta 1 (control)',
    Posterior_Mean=round(mean(output$samples.beta.co[,2]), 3),
    Posterior_Median=round(median(output$samples.beta.co[,2]), 3),
    Posterior_Variance=round(var(output$samples.beta.co[,2]), 3)
  )
  params[[6]] <- list(
    Name='Beta 2 (control)',
    Posterior_Mean=round(mean(output$samples.beta.co[,3]), 3),
    Posterior_Median=round(median(output$samples.beta.co[,3]), 3),
    Posterior_Variance=round(var(output$samples.beta.co[,3]), 3)
  )
  params[[7]] <- list(
    Name='Alpha (case)',
    Posterior_Mean=round(mean(output$samples.alpha.ca), 3),
    Posterior_Median=round(median(output$samples.alpha.ca), 3),
    Posterior_Variance=round(var(output$samples.alpha.ca), 3)
  )
  params[[8]] <- list(
    Name='Alpha (control)',
    Posterior_Mean=round(mean(output$samples.alpha.co), 3),
    Posterior_Median=round(median(output$samples.alpha.co), 3),
    Posterior_Variance=round(var(output$samples.alpha.co), 3)
  )
  params[[9]] <- list(
    Name='Range',
    Posterior_Mean=round(mean(output$samples.theta), 3),
    Posterior_Median=round(median(output$samples.theta), 3),
    Posterior_Variance=round(var(output$samples.theta), 3)
  )
  params[[10]] <- list(
    Name='Marginal Variance',
    Posterior_Mean=round(mean(output$samples.phi), 3),
    Posterior_Median=round(median(output$samples.phi), 3),
    Posterior_Variance=round(var(output$samples.phi), 3)
  )
  return(ldply(params, 'data.frame'))
  
}


#' compare_params
#'
#' @param beta.ca.hat_p (numeric) vector of parameter estimates from case poisson model
#' @param beta.co.hat_p (numeric) vector of parameter estimates from control poisson model
#' @param beta.ca.hat (numeric) vector of parameter estimates from case preferential sampling model
#' @param beta.co.hat (numeric) vector of parameter estimates from control preferential sampling model
#'
#' @return
#' @export
#'
#' @examples
compare_params <- function(beta.ca.hat_p, beta.co.hat_p, beta.ca.hat, beta.co.hat, beta_ca_sp, beta_co_sp){
  
  param_comp <- list()
  param_comp[[1]] <- list(
    Parameter='Beta 0 (case)',
    Model='Poisson',
    Estimate=round(beta.ca.hat_p[1], 3)
  )
  param_comp[[2]] <- list(
    Parameter='Beta 0 (case)',
    Model='Spatial Poisson',
    Estimate=round(beta_ca_sp[1], 3)
  )
  param_comp[[3]] <- list(
    Parameter='Beta 0 (case)',
    Model='Preferential Sampling',
    Estimate=round(beta.ca.hat[1], 3)
  )
  param_comp[[4]] <- list(
    Parameter='Beta 1 (case)',
    Model='Poisson',
    Estimate=round(beta.ca.hat_p[2], 3)
  )
  param_comp[[5]] <- list(
    Parameter='Beta 1 (case)',
    Model='Spatial Poisson',
    Estimate=round(beta_ca_sp[2], 3)
  )
  param_comp[[6]] <- list(
    Parameter='Beta 1 (case)',
    Model='Preferential Sampling',
    Estimate=round(beta.ca.hat[2], 3)
  )
  param_comp[[7]] <- list(
    Parameter='Beta 2 (case)',
    Model='Poisson',
    Estimate=round(beta.ca.hat_p[3], 3)
  )
  param_comp[[8]] <- list(
    Parameter='Beta 2 (case)',
    Model='Spatial Poisson',
    Estimate=round(beta_ca_sp[3], 3)
  )
  param_comp[[9]] <- list(
    Parameter='Beta 2 (case)',
    Model='Preferential Sampling',
    Estimate=round(beta.ca.hat[3], 3)
  )
  
  param_comp[[10]] <- list(
    Parameter='Beta 0 (control)',
    Model='Poisson',
    Estimate=round(beta.co.hat_p[1], 3)
  )
  param_comp[[11]] <- list(
    Parameter='Beta 0 (control)',
    Model='Spatial Poisson',
    Estimate=round(beta_co_sp[1], 3)
  )
  param_comp[[12]] <- list(
    Parameter='Beta 0 (control)',
    Model='Preferential Sampling',
    Estimate=round(beta.co.hat[1], 3)
  )
  param_comp[[13]] <- list(
    Parameter='Beta 1 (control)',
    Model='Poisson',
    Estimate=round(beta.co.hat_p[2], 3)
  )
  param_comp[[14]] <- list(
    Parameter='Beta 1 (control)',
    Model='Spatial Poisson',
    Estimate=round(beta_co_sp[2], 3)
  )
  param_comp[[15]] <- list(
    Parameter='Beta 1 (control)',
    Model='Preferential Sampling',
    Estimate=round(beta.co.hat[2], 3)
  )
  param_comp[[16]] <- list(
    Parameter='Beta 2 (control)',
    Model='Poisson',
    Estimate=round(beta.co.hat_p[3], 3)
  )
  param_comp[[17]] <- list(
    Parameter='Beta 2 (control)',
    Model='Spatial Poisson',
    Estimate=round(beta_co_sp[3], 3)
  )
  param_comp[[18]] <- list(
    Parameter='Beta 2 (control)',
    Model='Preferential Sampling',
    Estimate=round(beta.co.hat[3], 3)
  )
  return(ldply(param_comp, 'data.frame'))
  
}
