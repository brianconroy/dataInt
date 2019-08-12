

calc_significance_rasters <- function(rodents_data, output, caPr, null_alphas=F){
  
  if (ncol(output$samples.w) == 788){
    agg_factor <- 6
  } else{
    agg_factor <- 5
  }
  
  caPr.disc <- aggregate(caPr, fact=agg_factor)
  loc.disc <- caPr.disc[[1]]
  all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
  N <- n_values(caPr.disc[[1]])
  
  data <- assemble_data(rodents_data, loc.disc, caPr.disc)
  X_rodent <- load_x_standard2(as.logical(data$loc$status), agg_factor=agg_factor)
  risk_rodent <- calc_posterior_risk(output, X_rodent, null_alphas=null_alphas)
  
  threshold <- 0.05
  fracs <- apply(risk_rodent, 2, function(x){sum(x > threshold)/nrow(risk_rodent)})
  r <- caPr.disc[[1]]
  r[][!is.na(r[])] <- fracs
  
  indicators <- fracs > 0.95
  r2 <- caPr.disc[[1]]
  r2[][!is.na(r2[])] <- indicators
  
  return(list(
    r_fracs=r,
    r_inds=r2
  ))
  
}


calc_posterior_risk <- function(output, samples_int){
  
  samples.beta.ca <- output$samples.beta.ca
  samples.beta.co <- output$samples.beta.co
  samples.alpha.ca <- output$samples.alpha.ca
  samples.alpha.co <- output$samples.alpha.co
  X_high <- load_x_ca2()
  n.samples <- nrow(samples_int)
  samples.risk <- c()
  
  samples.l <- lapply(1:n.samples, function(i){
    w_i <- samples_int[i,]
    beta.ca_i <- samples.beta.ca[i,]
    beta.co_i <- samples.beta.co[i,]
    alpha.ca_i <- samples.alpha.ca[i]
    alpha.co_i <- samples.alpha.co[i]
    lodds_i <- X_high %*% beta.ca_i + alpha.ca_i * w_i - X_high %*% beta.co_i - alpha.co_i * w_i
    t(calc_risk(lodds_i))
  })
  samples.risk <- do.call(rbind, samples.l)
  
  return(samples.risk)
  
}


calc_posterior_risk_temporal <- function(output, samples_int){
  
  samples.beta.ca <- output$samples.beta.ca
  samples.beta.co <- output$samples.beta.co
  samples.alpha.ca <- output$samples.alpha.ca
  samples.alpha.co <- output$samples.alpha.co
  samples.u <- output$samples.u
  n.samples <- nrow(samples_int)
  nyears <- ncol(samples.u)
  years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
  
  X_highs <- lapply(years, function(y){
    load_x_time(year=y, agg_factor=1)
  })
  
  samples.temporal <- lapply(1:nyears, function(t){
    print(t)
    X_high <- X_highs[[t]]
    samples.l <- lapply(1:n.samples, function(i){
      w_i <- samples_int[i,]
      beta.ca_i <- samples.beta.ca[i,]
      beta.co_i <- samples.beta.co[i,]
      alpha.ca_i <- samples.alpha.ca[i]
      alpha.co_i <- samples.alpha.co[i]
      u.i <- samples.u[i,t]
      lodds_i <- X_high %*% beta.ca_i + alpha.ca_i * (u.i + w_i) - X_high %*% beta.co_i - alpha.co_i * (u.i + w_i)
      t(calc_risk(lodds_i))
    })
    do.call(rbind, samples.l)
  })
  
  return(samples.temporal)
  
}


calc_significance <- function(samples.risk, r, threshold=0.05){
  
  pvals <- sapply(1:ncol(samples.risk), function(i){
    r_i <- samples.risk[,i]
    sum(r_i > threshold)/length(r_i)
  })
  
  rp <- overlay(pvals, r)
  
  inds_95 <- 1 * (pvals > 0.95)
  inds_50 <- 1 * (pvals > 0.5)
  inds_25 <- 1 * (pvals > 0.25)
  # normalize color scales
  if (sum(inds_95) == 0){inds_95[1]=3}
  if (sum(inds_50) == 0){inds_95[2]=2}
  if (sum(inds_25) == 0){inds_95[3]=1}
  
  inds <- inds_95 + inds_50 + inds_25
  r_inds <- overlay(inds, r)
  r_inds_95 <- overlay(inds_95, r)
  
  return(list(
    r_fracs=rp,
    r_inds=r_inds,
    r_inds_95=r_inds_95
  ))
  
}


calc_significance_temporal <- function(samples.temporal, r, threshold=0.05){
  
  pvals <- lapply(1:length(samples.temporal), function(t){
      sapply(1:ncol(samples.temporal[[t]]), function(i){
        r_i <- samples.temporal[[t]][,i]
        sum(r_i > threshold)/length(r_i)
      })
  })
  
  rp <- lapply(1:length(samples.temporal), function(t){
    overlay(pvals[[t]], r)
  })
  
  inds_95 <- lapply(1:length(samples.temporal), function(t){ 1 * (pvals[[t]] > 0.95)})
  inds_50 <- lapply(1:length(samples.temporal), function(t){ 1 * (pvals[[t]] > 0.5)})
  inds_25 <- lapply(1:length(samples.temporal), function(t){ 1 * (pvals[[t]] > 0.25)})
  
  inds <- lapply(1:length(samples.temporal), function(t){
    inds_95[[t]] + inds_50[[t]] + inds_25[[t]]
  })
  r_inds <- lapply(1:length(samples.temporal), function(t){
    overlay(inds[[t]], r)
  })
  r_inds_95 <- lapply(1:length(samples.temporal), function(t){
    overlay(inds_95[[t]], r)
  })
    
  return(list(
    r_fracs=rp,
    r_inds=r_inds,
    r_inds_95=r_inds_95
  ))
  
}


calc_significance_rasters_ds <- function(rodents_data, output, caPr, threshold, null_alphas=F){
  
  if (ncol(output$samples.w) == 788){
    agg_factor <- 6
  } else if (ncol(output$samples.w) == 584){
    agg_factor <- 7
  }else{
    agg_factor <- 5
  }
  
  caPr.disc <- aggregate(caPr, fact=agg_factor)
  loc.disc <- caPr.disc[[1]]
  all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
  N <- n_values(caPr.disc[[1]])
  
  # sample from posterior distributions at high resolution
  data <- assemble_data(rodents_data, loc.disc, caPr.disc)
  X_rodent <- load_x_standard2(as.logical(data$loc$status), agg_factor=agg_factor)
  risk_rodent <- calc_posterior_risk(output, X_rodent, null_alphas=null_alphas)
  
  # downscale posterior variances
  postvar_rodent <- apply(risk_rodent, 2, var)
  r <- overlay(postvar_rodent, caPr.disc[[1]])
  xy <- data.frame(xyFromCell(r, 1:ncell(r)))
  v <- getValues(r)
  tps <- Tps(xy, v)
  p <- raster(caPr[[2]])
  p <- interpolate(p, tps)
  p <- mask(p, caPr[[1]])
  postvar_ds <- p[][!is.na(p)[]]
  postvar_ds[postvar_ds <= 0] <- 1e-5
  
  # calculate risk estimates at high resolution
  results_ds <- calc_risk_cdph_output(
    output, 
    data, 
    caPr.disc, 
    all_ids, 
    agg_factor=agg_factor, 
    null_alphas=null_alphas)$r_risk_high
  risk_est <- results_ds[][!is.na(results_ds[])]

  # calculate pvalues
  pvals <- c()
  for (i in 1:length(risk_est)){
    r_i <- risk_est[i]
    v_i <- postvar_ds[i]
    pvals <- c(pvals, 1 - pnorm(threshold, mean=r_i, sd=sqrt(v_i)))
  }
  rp <- overlay(pvals, caPr[[1]])

  inds_95 <- 1 * (pvals > 0.95)
  inds_50 <- 1 * (pvals > 0.5)
  inds_25 <- 1 * (pvals > 0.25)
  # normalize color scales
  if (sum(inds_95) == 0){inds_95[1]=3}
  if (sum(inds_50) == 0){inds_95[2]=2}
  if (sum(inds_25) == 0){inds_95[3]=1}
  
  inds <- inds_95 + inds_50 + inds_25
  r_inds <- overlay(inds, caPr[[1]])
  r_inds_95 <- overlay(inds_95, caPr[[1]])
  
  return(list(
    r_fracs=rp,
    r_inds=r_inds,
    r_inds_95=r_inds_95
  ))
  
}


calc_significance_rasters_ds_temporal <- function(data, output, caPr_all, caPr.disc_all, threshold){
  
  if (ncol(output$samples.w) == 788){
    agg_factor <- 6
  } else if (ncol(output$samples.w) == 584){
    agg_factor <- 7
  }else{
    agg_factor <- 5
  }
  
  loc.disc <- caPr.disc_all[[1]][[1]]
  all_ids <- c(1:length(loc.disc[]))[!is.na(loc.disc[])]
  N <- n_values(caPr.disc[[1]])
  years <- c(1983, 1988, 1993, 1998, 2003, 2008, 2013)
  
  risks_all_t <- calc_temporal_risks(output, caPr.disc_all, caPr_all, agg_factor, years)
  
  results <- list()
  counter <- 1
  for (t in 1:ncol(output$samples.u)){
    
    # sample from posterior distributions at high resolution
    data_y <- list(locs=data$locs[[t]], case.data=data$case.data[[t]], ctrl.data=data$ctrl.data[[t]])
    X_rodent <- load_x_time(year=years[t], agg_factor=agg_factor)
    risk_rodent_y <- calc_posterior_risk(output, X_rodent, null_alphas=F)
    
    # downscale posterior variances
    postvar_rodent <- apply(risk_rodent_y, 2, var)
    r <- overlay(postvar_rodent, caPr.disc[[1]])
    xy <- data.frame(xyFromCell(r, 1:ncell(r)))
    v <- getValues(r)
    tps <- Tps(xy, v)
    p <- raster(caPr[[2]])
    p <- interpolate(p, tps)
    p <- mask(p, caPr[[1]])
    postvar_ds <- p[][!is.na(p)[]]
    postvar_ds[postvar_ds <= 0] <- 1e-5
    
    # get risk estimates at high resolution
    risk_est <- risks_all_t$risks_high[[t]]
    
    # calculate pvalues
    pvals <- c()
    for (i in 1:length(risk_est)){
      r_i <- risk_est[i]
      v_i <- postvar_ds[i]
      pvals <- c(pvals, 1 - pnorm(threshold, mean=r_i, sd=sqrt(v_i)))
    }
    rp <- overlay(pvals, caPr_all[[1]][[1]])
    
    inds_95 <- 1 * (pvals > 0.95)
    inds_50 <- 1 * (pvals > 0.5)
    inds_25 <- 1 * (pvals > 0.25)
    
    inds <- inds_95 + inds_50 + inds_25
    r_inds <- overlay(inds, caPr[[1]])
    r_inds_95 <- overlay(inds_95, caPr[[1]])
    
    results[[counter]] <- list(
      r_fracs=rp,
      r_inds=r_inds,
      r_inds_95=r_inds_95
    )
    counter <- counter + 1
    
  }
  
  return(results)
  
}


calc_temporal_risks <- function(output, caPr.disc_all, caPr_all, agg_factor, years){
  
  # interpolate w
  w.hat <- colMeans(output$samples.w)
  rw <- caPr.disc_all[[1]][[1]]
  rw[][!is.na(rw[])] <- w.hat
  
  xy <- data.frame(xyFromCell(rw, 1:ncell(rw)))
  v <- getValues(rw)
  
  tps <- Tps(xy, v)
  p <- raster(caPr_all[[1]][[2]])
  p <- interpolate(p, tps)
  p <- mask(p, caPr_all[[1]][[1]])
  w.hat_ds <- p[][!is.na(p[])]
  
  # calculate risks
  alpha.ca.hat <- mean(output$samples.alpha.ca)
  alpha.co.hat <- mean(output$samples.alpha.co)
  beta.ca.hat <- colMeans(output$samples.beta.ca)
  beta.co.hat <- colMeans(output$samples.beta.co)
  u.hat <- colMeans(output$samples.u)
  
  risks_low <- list()
  risks_high <- list()
  r_risks_low <- list()
  r_risks_high <- list()
  for (i in 1:length(years)){
    
    y <- years[i]
    u.hat_y <- u.hat[i]
    
    # low resolution risks
    X_low <- load_x_time(year=y, agg_factor=agg_factor)
    lodds_low <- X_low %*% beta.ca.hat + alpha.ca.hat * (u.hat_y + w.hat) - X_low %*% beta.co.hat - alpha.co.hat * (u.hat_y + w.hat)
    risk_low <- calc_risk(lodds_low)
    
    r_lodds_low <- caPr.disc_all[[i]][[1]]
    r_lodds_low[][!is.na(r_lodds_low[])] <- lodds_low
    
    r_risk_low <- r_lodds_low
    r_risk_low[][!is.na(r_risk_low[])] <- risk_low
    
    risks_low[[i]] <- risk_low
    r_risks_low[[i]] <- r_risk_low
    
    # high resolution
    X_high <- load_x_time(year=y, agg_factor=1)
    lodds_high <- X_high %*% beta.ca.hat + alpha.ca.hat * (u.hat_y + w.hat_ds) - X_high %*% beta.co.hat - alpha.co.hat * (u.hat_y + w.hat_ds)
    risk_high <- calc_risk(lodds_high)
    
    r_lodds_high <- caPr_all[[i]][[1]]
    r_lodds_high[][!is.na(r_lodds_high[])] <- lodds_high
    
    r_risk_high <- r_lodds_high
    r_risk_high[][!is.na(r_risk_high[])] <- risk_high
    
    risks_high[[i]] <- risk_high
    r_risks_high[[i]] <- r_risk_high
  }
  
  return(
    list(
      risks_low=risks_low,
      risks_high=risks_high,
      r_risks_low=r_risks_low,
      r_risks_high=r_risks_high
    )
  )
  
}


calc_lodds_ds <- function(alpha.ca, alpha.co, beta.ca, beta.co, w.hat_ds){
  
  X_high <- load_x_ca2()
  lodds_high <- X_high %*% beta.ca + alpha.ca * w.hat_ds - X_high %*% beta.co - alpha.co * w.hat_ds
  return(lodds_high)
  
}


calc_risk_ds <- function(alpha.ca, alpha.co, beta.ca, beta.co, w.hat_ds){
  
  X_high <- load_x_ca2()
  lodds_high <- X_high %*% beta.ca + alpha.ca * w.hat_ds - X_high %*% beta.co - alpha.co * w.hat_ds
  risk_high <- calc_risk(lodds_high)
  return(risk_high)
  
}


plot_risk_overlay <- function(results, rodents_species, main=''){
  
  if ("r_risk_high" %in% names(results)){
    plot(results$r_risk_high, main=main)
  } else{
    plot(results, main=main)
  }
  r_cases <- rodents_species[rodents_species$Res == 'POS',]
  r_coords <- cbind(r_cases$Lon_Add_Fix, r_cases$Lat_Add_Fix)
  r_ctrls <- rodents_species[rodents_species$Res == 'NEG',]
  r_coords_ctrl <- cbind(r_ctrls$Lon_Add_Fix, r_ctrls$Lat_Add_Fix)
  points(r_coords_ctrl, col=rgb(0,0,1,0.05), pch=16, cex=0.65)
  points(r_coords, col=rgb(1,0,0,0.2), pch=16, cex=0.65)
  
}


plot_cov_vs_w <- function(results, caPr){
  
  X_high <- load_x_ca2()
  cov_rodent <- X_high %*% results$beta.ca.hat - X_high %*% results$beta.co.hat
  w_rodent <- (results$alpha.ca.hat -  results$alpha.co.hat) * results$w.hat_ds
  r_cov_rodent <- caPr[[1]]
  r_cov_rodent[][!is.na(r_cov_rodent[])] <- cov_rodent
  r_w_rodent <- caPr[[1]]
  r_w_rodent[][!is.na(r_w_rodent[])] <- w_rodent
  par(mfrow=c(1,2))
  plot(r_cov_rodent, main='A)')
  pal <- colorRampPalette(c("blue","red"))
  plot(r_w_rodent, main='B)', col=pal(20))
  par(mfrow=c(1,1))
  
}


calc_risk_cdph <- function(species, rodents, caPr.disc, all_ids, agg_factor=5, null_alphas=FALSE){
  
  print(species)
  if (paste(species, collapse="") == 'all_but_ds'){
    all_species <- unique(rodents$Short_Name)
    species_group <- as.character(all_species[all_species != 'Pine Squirrel'])
    rodents_species <- rodents[rodents$Short_Name %in% species_group,]
  } else {
    rodents_species <- rodents[rodents$Short_Name %in% species,]
  }
  analysis_name <- gsub(',', '', gsub(' ', '_', paste('analysis', paste(species, collapse="_"), sep='_'), fixed=T))
  output <- load_output(paste("cdph_", analysis_name, ".json", sep=""))
  data <- assemble_data(rodents_species, loc.disc, caPr.disc)
  
  coords <- xyFromCell(caPr.disc, cell=all_ids)
  d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
  
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
  
  if (null_alphas){
    alpha.ca.hat <- 0
    alpha.co.hat <- 0
  } else{
    alpha.ca.hat <- mean(output$samples.alpha.ca)
    alpha.co.hat <- mean(output$samples.alpha.co)
  }
  beta.ca.hat <- colMeans(output$samples.beta.ca)
  beta.co.hat <- colMeans(output$samples.beta.co)
  
  # risk map (low resolution)
  X_low <- load_x_ca2(factor=agg_factor)
  lodds_low <- X_low %*% beta.ca.hat + alpha.ca.hat * w.hat - X_low %*% beta.co.hat - alpha.co.hat * w.hat
  risk_low <- calc_risk(lodds_low)
  
  r_lodds_low <- caPr.disc[[1]]
  r_lodds_low[][!is.na(r_lodds_low[])] <- lodds_low
  
  r_risk_low <- caPr.disc[[1]]
  r_risk_low[][!is.na(r_risk_low[])] <- risk_low
  
  # risk map (downscaled)
  X_high <- load_x_ca2()
  lodds_high <- X_high %*% beta.ca.hat + alpha.ca.hat * w.hat_ds - X_high %*% beta.co.hat - alpha.co.hat * w.hat_ds
  risk_high <- calc_risk(lodds_high)
  
  r_lodds_high <- caPr[[2]]
  r_lodds_high[][!is.na(r_lodds_high[])] <- lodds_high
  
  r_risk_high <- caPr[[2]]
  r_risk_high[][!is.na(r_risk_high[])] <- risk_high
  
  response <- list(
    r_risk_high=r_risk_high,
    alpha.ca.hat=alpha.ca.hat,
    alpha.co.hat=alpha.co.hat,
    beta.ca.hat=beta.ca.hat,
    beta.co.hat=beta.co.hat,
    w.hat_ds=w.hat_ds
  )
  return(response)
  
}


calc_risk_cdph_output <- function(output, data, caPr.disc, all_ids, agg_factor=5, null_alphas=FALSE){
  
  coords <- xyFromCell(caPr.disc, cell=all_ids)
  d <- as.matrix(dist(coords, diag=TRUE, upper=TRUE))
  
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
  
  if (null_alphas){
    alpha.ca.hat <- 0
    alpha.co.hat <- 0
  } else{
    alpha.ca.hat <- mean(output$samples.alpha.ca)
    alpha.co.hat <- mean(output$samples.alpha.co)
  }
  beta.ca.hat <- colMeans(output$samples.beta.ca)
  beta.co.hat <- colMeans(output$samples.beta.co)
  
  # risk map (low resolution)
  X_low <- load_x_ca2(factor=agg_factor)
  lodds_low <- X_low %*% beta.ca.hat + alpha.ca.hat * w.hat - X_low %*% beta.co.hat - alpha.co.hat * w.hat
  risk_low <- calc_risk(lodds_low)
  
  r_lodds_low <- caPr.disc[[1]]
  r_lodds_low[][!is.na(r_lodds_low[])] <- lodds_low
  
  r_risk_low <- caPr.disc[[1]]
  r_risk_low[][!is.na(r_risk_low[])] <- risk_low
  
  # risk map (downscaled)
  X_high <- load_x_ca2()
  lodds_high <- X_high %*% beta.ca.hat + alpha.ca.hat * w.hat_ds - X_high %*% beta.co.hat - alpha.co.hat * w.hat_ds
  risk_high <- calc_risk(lodds_high)
  
  r_lodds_high <- caPr[[2]]
  r_lodds_high[][!is.na(r_lodds_high[])] <- lodds_high
  
  r_risk_high <- caPr[[2]]
  r_risk_high[][!is.na(r_risk_high[])] <- risk_high
  
  response <- list(
    r_risk_high=r_risk_high,
    alpha.ca.hat=alpha.ca.hat,
    alpha.co.hat=alpha.co.hat,
    beta.ca.hat=beta.ca.hat,
    beta.co.hat=beta.co.hat,
    w.hat_ds=w.hat_ds
  )
  return(response)
  
}


equalize_scales <- function(r1, r2){
  
  v1 <- r1[][!is.na(r1[])]
  v2 <- r2[][!is.na(r2[])]
  r_max <- max(v1, v2)
  r_min <- min(v1, v2)
  if (sum(v1 == r_max) == 0){
    v1[length(v1)] <- r_max
  } else{
    v2[length(v1)] <- r_max
  }
  if (sum(v1 == r_min) == 0){
    v1[1] <- r_min
  } else{
    v2[1] <- r_min
  }
  r1[][!is.na(r1[])] <- v1
  r2[][!is.na(r2[])] <- v2
  return(list(r1, r2))
  
}


equalize_scales2 <- function(r_list){
  
  vs <- list()
  r_max <- 0
  r_min <- 1
  r_list_new <- list()
  for (i in 1:length(r_list)){
    r <- r_list[[i]]
    r_vals <- r[][!is.na(r[])]
    vs[[i]] <- r_vals
    r_max <- max(r_max, r_vals)
    r_min <- min(r_min, r_vals)
  }
  for (i in 1:length(r_list)){
    r <- r_list[[i]]
    r_vals <- r[][!is.na(r[])]
    if (sum(r_vals == r_max) == 0){
      r_vals[length(r_vals)] <- r_max
    }
    if (sum(r_vals == r_min) == 0){
      r_vals[1] <- r_min
    }
    r_new <- r
    r_new[][!is.na(r_new[])] <- r_vals 
    r_list_new[[i]] <- r_new
  }
  return(r_list_new)
  
}


equalize_scales3 <- function(r_list){
  
  vs <- list()
  r_max <- 0
  r_min <- 1
  r_list_new <- list()
  for (i in 1:length(r_list)){
    result_i <- r_list[[i]]
    r <- result_i$r_risk_high
    r_vals <- r[][!is.na(r[])]
    vs[[i]] <- r_vals
    r_max <- max(r_max, r_vals)
    r_min <- min(r_min, r_vals)
  }
  for (i in 1:length(r_list)){
    result_i <- r_list[[i]]
    r <- result_i$r_risk_high
    r_vals <- r[][!is.na(r[])]
    if (sum(r_vals == r_max) == 0){
      r_vals[length(r_vals)] <- r_max
    }
    if (sum(r_vals == r_min) == 0){
      r_vals[1] <- r_min
    }
    r_new <- r
    r_new[][!is.na(r_new[])] <- r_vals 
    result_i$r_risk_high <- r_new
    r_list_new[[i]] <- result_i
  }
  return(r_list_new)
  
}


equalize_scales4 <- function(r_list){
  
  vs <- list()
  r_max <- 0
  r_min <- 1
  r_list_new <- list()
  for (i in 1:length(r_list)){
    r <- r_list[[i]]
    r_vals <- r[][!is.na(r[])]
    vs[[i]] <- r_vals
    r_max <- max(r_max, r_vals)
    r_min <- min(r_min, r_vals)
  }
  for (i in 1:length(r_list)){
    r <- r_list[[i]]
    r_vals <- r[][!is.na(r[])]
    if (sum(r_vals == r_max) == 0){
      r_vals[length(r_vals)] <- r_max
    }
    if (sum(r_vals == r_min) == 0){
      r_vals[1] <- r_min
    }
    r_new <- r
    r_new[][!is.na(r_new[])] <- r_vals 
    r_list_new[[i]] <- r_new
  }
  return(r_list_new)
  
}


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


compare_params <- function(output, output.sp_ca, output.sp_co, mod.ca, mod.co){
  
  beta.ca.hat <- colMeans(output$samples.beta.ca)
  beta.co.hat <- colMeans(output$samples.beta.co)
  beta.ca.var <- apply(output$samples.beta.ca, 2, var)
  beta.co.var <- apply(output$samples.beta.co, 2, var)
  
  beta_ca_sp <- colMeans(output.sp_ca$samples.beta)
  beta_co_sp <- colMeans(output.sp_co$samples.beta)
  beta_ca_sp_var <- apply(output.sp_ca$samples.beta, 2, var)
  beta_co_sp_var <- apply(output.sp_co$samples.beta, 2, var)
  
  beta.ca.hat_p <- unname(coefficients(mod.ca))
  beta.co.hat_p <- unname(coefficients(mod.co))
  beta.ca.var_p <- diag(vcov(mod.ca))
  beta.co.var_p <- diag(vcov(mod.co))
  
  param_comp <- list()
  param_comp[[1]] <- list(
    Parameter='Beta 0 (case)',
    Model='Poisson',
    Estimate=round(beta.ca.hat_p[1], 3),
    Variance=round(beta.ca.var_p[1], 5)
  )
  param_comp[[2]] <- list(
    Parameter='Beta 0 (case)',
    Model='Spatial Poisson',
    Estimate=round(beta_ca_sp[1], 3),
    Variance=round(beta_ca_sp_var[1], 5)
  )
  param_comp[[3]] <- list(
    Parameter='Beta 0 (case)',
    Model='Preferential Sampling',
    Estimate=round(beta.ca.hat[1], 3),
    Variance=round(beta.ca.var[1], 3)
  )
  param_comp[[4]] <- list(
    Parameter='Beta 1 (case)',
    Model='Poisson',
    Estimate=round(beta.ca.hat_p[2], 3),
    Variance=round(beta.ca.var_p[2], 5)
  )
  param_comp[[5]] <- list(
    Parameter='Beta 1 (case)',
    Model='Spatial Poisson',
    Estimate=round(beta_ca_sp[2], 3),
    Variance=round(beta_ca_sp_var[2], 5)
  )
  param_comp[[6]] <- list(
    Parameter='Beta 1 (case)',
    Model='Preferential Sampling',
    Estimate=round(beta.ca.hat[2], 3),
    Variance=round(beta.ca.var[2], 3)
  )
  param_comp[[7]] <- list(
    Parameter='Beta 2 (case)',
    Model='Poisson',
    Estimate=round(beta.ca.hat_p[3], 3),
    Variance=round(beta.ca.var_p[3], 5)
  )
  param_comp[[8]] <- list(
    Parameter='Beta 2 (case)',
    Model='Spatial Poisson',
    Estimate=round(beta_ca_sp[3], 3),
    Variance=round(beta_ca_sp_var[3], 5)
  )
  param_comp[[9]] <- list(
    Parameter='Beta 2 (case)',
    Model='Preferential Sampling',
    Estimate=round(beta.ca.hat[3], 3),
    Variance=round(beta.ca.var[3], 3)
  )
  
  param_comp[[10]] <- list(
    Parameter='Beta 0 (control)',
    Model='Poisson',
    Estimate=round(beta.co.hat_p[1], 3),
    Variance=round(beta.co.var_p[1], 5)
  )
  param_comp[[11]] <- list(
    Parameter='Beta 0 (control)',
    Model='Spatial Poisson',
    Estimate=round(beta_co_sp[1], 3),
    Variance=round(beta_co_sp_var[1], 5)
  )
  param_comp[[12]] <- list(
    Parameter='Beta 0 (control)',
    Model='Preferential Sampling',
    Estimate=round(beta.co.hat[1], 3),
    Variance=round(beta.co.var[1], 3)
  )
  param_comp[[13]] <- list(
    Parameter='Beta 1 (control)',
    Model='Poisson',
    Estimate=round(beta.co.hat_p[2], 3),
    Variance=round(beta.co.var_p[2], 5)
  )
  param_comp[[14]] <- list(
    Parameter='Beta 1 (control)',
    Model='Spatial Poisson',
    Estimate=round(beta_co_sp[2], 3),
    Variance=round(beta_co_sp_var[2], 5)
  )
  param_comp[[15]] <- list(
    Parameter='Beta 1 (control)',
    Model='Preferential Sampling',
    Estimate=round(beta.co.hat[2], 3),
    Variance=round(beta.co.var[2], 3)
  )
  param_comp[[16]] <- list(
    Parameter='Beta 2 (control)',
    Model='Poisson',
    Estimate=round(beta.co.hat_p[3], 3),
    Variance=round(beta.co.var_p[3], 5)
  )
  param_comp[[17]] <- list(
    Parameter='Beta 2 (control)',
    Model='Spatial Poisson',
    Estimate=round(beta_co_sp[3], 3),
    Variance=round(beta_co_sp_var[3], 5)
  )
  param_comp[[18]] <- list(
    Parameter='Beta 2 (control)',
    Model='Preferential Sampling',
    Estimate=round(beta.co.hat[3], 3),
    Variance=round(beta.co.var[3], 3)
  )
  return(ldply(param_comp, 'data.frame'))
  
}
