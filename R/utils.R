

gen_sim_name <- function(sampling, prevalence){
  
  return(paste("prefSampleGpCC", sampling, prevalence, sep="_"))
  
}


load_sim_params <- function(sampling, prevalence){
  
  sim_params <- read.table("/Users/brianconroy/Documents/research/dataInt/output/simParams.txt",
                           stringsAsFactors=F, header=T)
  params <- sim_params[sim_params$sampling == sampling & sim_params$prevalence == prevalence,]
  return(params)
  
}


load_sim_params_size <- function(size){
  
  sim_params <- read.table("/Users/brianconroy/Documents/research/dataInt/output/simParams_size.txt",
                           stringsAsFactors=F, header=T)
  params <- sim_params[sim_params$size == size,]
  return(params)
  
}


load_sim_params_priorcompare <- function(){
  
  sim_params <- read.table("/Users/brianconroy/Documents/research/dataInt/output/simParams_priorcompare.txt",
                           stringsAsFactors=F, header=T)
  return(sim_params)
  
}


ig_var <- function(a, b){
  
  return(b^2/((a-1)^2 * (a-2)))
  
}


iterate_ig_variance <- function(phi){
  
  vars <- list()
  c <- 1
  for (b in 40:400){
    a <- b/phi + 1
    vars[[c]] <- list(a=a, b=b, v=ig_var(a, b))
    c <- c + 1
  }
  return(ldply(vars, 'data.frame'))

}


g_var <- function(a, b){
  
  return(a*(b^2))

}


iterate_g_variance <- function(theta, lower, upper){
  
  if (lower < 1){
    rng <- c(seq(lower, 1, by=0.25), 1:upper)
  } else{
    rng <- lower:upper
  }
  
  vars <- list()
  c <- 1
  for (b in rng){
    a <- round(theta/b, 3)
    vars[[c]] <- list(a=a, b=b, v=g_var(a, b))
    c <- c + 1
  }
  return(ldply(vars, 'data.frame'))
  
}


n_values <- function(r){
  
  return(length(r[][!is.na(r[])]))
  
}


#' load_prism_pcs
#'
#' @return raster stack of first 2 PRISM principal components
#' @export
#'
#' @examples
load_prism_pcs <- function(){
  
  file <- "Documents/research/dataInt/data/prism_pcas_ca.grd"
  return(stack(file))
  
}


#' summarize_param
#'
#' @param name parameter name
#' @param true true parameter value
#' @param estimated estimated parameter value
#'
#' @return list summarizing estimate and bias
#' @export
#'
#' @examples
summarize_param <- function(name, true, estimated){
  
  summ <- list()
  summ$name <- name
  summ$true <- true
  summ$estimated <- round(estimated, 2)
  summ$bias <- round(estimated - true, 2)
  return(summ)
  
}


#' overlay
#'
#' @param vals vector of ordered cell values
#' @param r raster
#'
#' @return new raster with replaced values
#' @export
#'
#' @examples
overlay <- function(vals, r){
  
  r[][!is.na(r[])] <- vals
  return(r)
  
}


#' combine_w
#' 
#' combines estimated and kriged random effects 
#' according to cell observation status 
#' 
#' @param w_est estimated random effects of observed sites
#' @param w_krige kriged random effects
#' @param status vector of TRUE/FALSE values indicating observation 
#'
#' @return vector of combined estimated and kriged random effects
#' @export
#'
#' @examples
combine_w <- function(w_est, w_krige, status){
  
  w_comb <- c()
  i1 <- 1
  i2 <- 1
  for (i in 1:length(status)){
    if (status[i]){
      w_comb <- c(w_comb, w_est[i1])
      i1 <- i1 + 1
    } else {
      w_comb <- c(w_comb, w_krige[i2])
      i2 <- i2 + 1
    }
  }
  return(w_comb)
  
}


#' view_tr
#'
#' @param samples vector of mcmc samples 
#' @param trueval optional true parameter value displayed as horizontal line
#' @param title optional plot title
#' @param ylim optional y axis limits
#'
#' @return
#' @export
#'
#' @examples
view_tr <- function(samples, trueval=NULL, title='', ylim=NULL){
  
  if (is.null(ylim)){
    plot(samples, type='l', main=title)
  } else {
    plot(samples, type='l', main=title, ylim=ylim)
  }
  
  if (!is.null(trueval)){
    abline(h=trueval, col='2')
  }
  
}


#' view_tr_w
#' 
#' multipanel (4x4) view of random effect traceplots
#'
#' @param samples matrix of mcmc samples
#' @param w_true optional vector of true random effect values
#' @param page view offset
#'
#' @return
#' @export
#'
#' @examples
view_tr_w <- function(samples, w_true=NULL, page=1){
  
  par(mfrow=c(4, 4))
  for(i in 1:16){
    plot(samples[,page*i], type='l', ylab='w')
    if (!is.null(w_true)){
      abline(h=w_true[page*i], col='2')
    }
  }
  par(mfrow=c(1, 1))
  
}


summarizeSamps <- function(samples, beta){
  
  resultsList <- list()
  
  for (i in 1:length(beta)){
    
    mu <- round(mean(samples[,i]), 3)
    resultsList[paste('mu.b', i-1, sep='')] <- mu
    resultsList[paste('sd.b', i-1, sep='')] <- round(sd(samples[,i]), 3)
    resultsList[paste('pbias.b', i-1, sep='')] <- round(100*(mu - beta[i])/beta[i], 3)
    
  }
  
  return(resultsList)
  
}


riskMSE <- function(surf.true, surf.est){
  
  inds <- !is.na(values(surf.true))
  diffs <- values(surf.true)[inds] - values(surf.est)[inds]
  return(round(mean(diffs^2), 3))
  
}


riskMAE <- function(surf.true, surf.est){
  
  inds <- !is.na(values(surf.true))
  diffs <- values(surf.true)[inds] - values(surf.est)[inds]
  return(round(mean(abs(diffs)),3))
  
}


summarizeSurf <- function(surf.true, surf.est, mod.name){
  
  summ <- list()
  summ$model <- mod.name
  summ$mse <- riskMSE(surf.true, surf.est)
  summ$mae <- riskMAE(surf.true, surf.est)
  return(summ)
  
}


calcIntens <- function(cov, params){
  
  inds <- !is.na(cov[[1]][])
  intens <- rep(params[1], sum(inds))
  for (i in 1:length(names(cov))){
    cov.i <- cov[[i]][][inds]
    mu.i <- mean(cov.i)
    sd.i <- sd(cov.i)
    intens <- intens + params[i+1] * (cov.i - mu.i)/sd.i
  }
  intens.surf <- cov[[1]]
  intens.surf[][inds] <- intens
  return(intens.surf)
  
}


rastDiff <- function(r1, r2){
  
  r.new <- r1
  inds <- !is.na(values(r.new))
  values(r.new)[inds] <- values(r.new)[inds] - values(r2)[inds]
  return(r.new)
  
}


# discards the first n samples of an mcmc output
trimSamps <- function(output, n){
  
  
  new.out <- output
  for (x in names(output)){
    if (grepl("samples", x)){
      
      if (grepl("rho", x)){
        
        samps <- new.out$samples.rho
        D <- dim(new.out$samples.rho)[1]
        n.keep <- dim(samps)[2] - n
        trimmed.rho <- samps[,(n+1):dim(samps)[2],]
        trimmed <- array(NA, c(D, n.keep, 1))
        trimmed[1,,] <- trimmed.rho
        
      } else if (grepl("phi", x)) {
        
        samps <- new.out$samples.phi
        D <- dim(new.out$samples.phi)[1]
        n.keep <- dim(samps)[2] - n
        trimmed.phi <- samps[,(n+1):dim(samps)[2],]
        trimmed <- array(NA, c(D, n.keep, dim(samps)[3]))
        trimmed[1,,] <- trimmed.phi
        
      } else {
        
        samps <- new.out[[x]]
        trimmed <- samps[(n+1):nrow(samps),]
        if (ncol(samps) == 1){
          trimmed <- matrix(trimmed)
        }
        
      }
      new.out[[x]] <- trimmed
    }
  }
  
  
  return(new.out)
  
  
}


# pools output of the carLeroux models
pool.results <- function(results){
  results.new <- list()
  rho.present <- 'samples.rho' %in% names(results[[1]])
  new.beta <- results[[1]]$samples.beta
  new.tau2 <- results[[1]]$samples.tau2
  if (rho.present){
    D <- dim(results[[1]]$samples.rho)[1]
    n.keep <- dim(results[[1]]$samples.rho)[2]
    new.rho <- array(NA, c(D, length(results)*n.keep, 1))
    for (d in 1:D){
      new.rho[d,1:n.keep,] <- results[[1]]$samples.rho[d,,]
    }
  }
  
  if (length(results) > 2){
    for (i in 2:length(results)){
      new.beta <- rbind(new.beta, results[[i]]$samples.beta)  
      new.tau2 <- rbind(new.tau2, results[[i]]$samples.tau2)
      if (rho.present){
        for (d in 1:D){
          new.rho[d,(n.keep*(i-1) + 1):(n.keep*i),] <- results[[i]]$samples.rho[d,,]
        }
      }
    }
  }
  
  results.new$samples.beta <- new.beta
  results.new$samples.tau2 <- new.tau2
  if (rho.present){
    results.new$samples.rho <- new.rho
  }
  results.new$pooled.accept <- rowMeans(sapply(results, function(x){x$accept}))
  
  
  return(results.new)
  
}


# pools output of the preferential sampling case control model
pool.results.pscc <- function(results){
  
  results.new <- list()
  new.beta.l <- results[[1]]$samples.beta.l
  new.alpha <- results[[1]]$samples.alpha
  new.beta.case <- results[[1]]$samples.beta.case
  new.beta.control <- results[[1]]$samples.beta.control
  new.accept <- results[[1]]$accept
  
  for (i in 2:length(results)){
    new.beta.l <- rbind(new.beta.l, results[[i]]$samples.beta.l)  
    new.alpha <- rbind(new.alpha, results[[i]]$samples.alpha)
    new.beta.case <- rbind(new.beta.case, results[[i]]$samples.beta.case)
    new.beta.control <- rbind(new.beta.control, results[[i]]$samples.beta.control)
    new.accept <- new.accept + results[[i]]$accept
  }
  
  new.accept <- new.accept/length(results)
  results.new$samples.beta.l <- new.beta.l
  results.new$samples.beta.case <- new.beta.case
  results.new$samples.beta.control <- new.beta.control
  results.new$samples.alpha <- new.alpha
  results.new$accept <- new.accept
  
  return(results.new)
  
}


# pools output of the case control differential 
# preferential sampling model
pool.results.psccd <- function(results){
  
  results.new <- list()
  new.beta.l <- results[[1]]$samples.beta.l
  new.alpha.case <- results[[1]]$samples.alpha.case
  new.alpha.control <- results[[1]]$samples.alpha.control
  new.beta.case <- results[[1]]$samples.beta.case
  new.beta.control <- results[[1]]$samples.beta.control
  new.accept <- results[[1]]$accept
  
  for (i in 2:length(results)){
    new.beta.l <- rbind(new.beta.l, results[[i]]$samples.beta.l)  
    new.alpha.case <- rbind(new.alpha.case, results[[i]]$samples.alpha.case)
    new.alpha.control <- rbind(new.alpha.control, results[[i]]$samples.alpha.control)
    new.beta.case <- rbind(new.beta.case, results[[i]]$samples.beta.case)
    new.beta.control <- rbind(new.beta.control, results[[i]]$samples.beta.control)
    new.accept <- new.accept + results[[i]]$accept
  }
  
  new.accept <- new.accept/length(results)
  results.new$samples.beta.l <- new.beta.l
  results.new$samples.beta.case <- new.beta.case
  results.new$samples.beta.control <- new.beta.control
  results.new$samples.alpha.case <- new.alpha.case
  results.new$samples.alpha.control <- new.alpha.control
  results.new$accept <- new.accept
  
  return(results.new)
  
}


# pools output of the preferential sampling model
pool.results.ps <- function(results){
  
  results.new <- list()
  new.beta.l <- results[[1]]$samples.beta.l
  new.beta.c <- results[[1]]$samples.beta.c
  new.alpha <- results[[1]]$samples.alpha
  new.accept <- results[[1]]$accept
  
  for (i in 2:length(results)){
    new.beta.l <- rbind(new.beta.l, results[[i]]$samples.beta.l)  
    new.beta.c <- rbind(new.beta.c, results[[i]]$samples.beta.c)
    new.alpha <- rbind(new.alpha, results[[i]]$samples.alpha)
    new.accept <- new.accept + results[[i]]$accept
  }
  
  new.accept <- new.accept/length(results)
  results.new$samples.beta.l <- new.beta.l
  results.new$samples.beta.c <- new.beta.c
  results.new$samples.alpha <- new.alpha
  results.new$accept <- new.accept
  
  return(results.new)
  
}


pool.results.car <- function(results){
  results.new <- list()
  new.beta <- results[[1]]$samples.beta
  new.tau2 <- results[[1]]$samples.tau2
  
  for (i in 2:length(results)){
    new.beta <- rbind(new.beta, results[[i]]$samples.beta)  
    new.tau2 <- rbind(new.tau2, results[[i]]$samples.tau2)
  }
  
  results.new$samples.beta <- new.beta
  results.new$samples.tau2 <- new.tau2
  
  return(results.new)
  
}


## arguments:
#   resultList: unnamed list of outputs of a leroux model
viewResults <- function(resultList, params, overall.main=''){
  if ('beta0' %in% params){
    viewTracesInt(resultList, 'beta0', line=beta0, title.type='pooled', overall.title=overall.main)
    viewHistsInt(resultList, 'beta0', line=beta0, title.type='pooled', overall.title=overall.main)
  }
  if ('beta1' %in% params){
    viewTracesInt(resultList, 'beta1', line=beta1, title.type='pooled', overall.title=overall.main)
    viewHistsInt(resultList, 'beta1', line=beta1, title.type='pooled', overall.title=overall.main)
  }
  if ('beta2' %in% params){
    viewTracesInt(resultList, 'beta2', line=beta2, title.type='pooled', overall.title=overall.main)
    viewHistsInt(resultList, 'beta2', line=beta2, title.type='pooled', overall.title=overall.main)
  }
  if ('tau2.1' %in% params){
    viewTracesInt(resultList, 'tau2', line=tau.sq.1, d=1, title.type='pooled', overall.title=overall.main)
    viewHistsInt(resultList, 'tau2', line=tau.sq.1, d=1, title.type='pooled', overall.title=overall.main)
  }
  if ('rho.1' %in% params){
    viewTracesInt(resultList, 'rho.1', line=rho.1, d=1, title.type='pooled', overall.title=overall.main)
    viewHistsInt(resultList, 'rho.1', line=rho.1, d=1, title.type='pooled', overall.title=overall.main)
  }
}


list2csv <- function(nestedList, outfile){
  
  outdir <- '/Users/brianconroy/Documents/research/dataInt/'
  df <- c()
  for (i in 1:length(nestedList)){
    df <- rbind(df, nestedList[[i]])
  }
  write.table(df, file=paste(outdir, outfile, sep=''), row.names=F, sep=',')
  
}


## arguments
#   output: list, output of a carLeroux function
#   params: list of parameter names to values
summarize <- function(output, params, dic=TRUE){
  
  summary <- list()
  
  counter <- 1
  for (n in names(params)){
    n.summary <- list()
    n.summary[['parameter']] <- n
    n.summary[['true']] <- params[[n]]
    n.summary[['posterior.mean']] <- postMean(n, output)
    n.summary[['posterior.median']] <- postMedian(n, output)
    n.summary[['posterior.sd']] <- postSd(n, output)
    bias.p <- bias(output, n, params[[n]])
    bias.perc.p <- bias(output, n, params[[n]], type='perc')
    n.summary[['bias']] <- bias.p
    n.summary[['percbias']] <- bias.perc.p
    if (dic){
      n.summary[['DIC']] <- output$fit$DIC
      n.summary[['percent deviance']] <- output$fit$percent_dev
    }
    summary[[counter]] <- n.summary
    counter <- counter + 1
  }
  
  if ('rho' %in% names(params)){
    summary[['accept.beta']] <- output$accept[1]
    summary[['accept.phi']] <- output$accept[2]
    summary[['accept.rho']] <- output$accept[3]
  }
  
  return(summary)

}


getSamples <- function(param, output){
  
  if (param == 'beta0'){
    samples <- output$samples.beta[,1]
  } else if (param == 'beta1'){
    samples <- output$samples.beta[,2]
  } else if (param == 'beta2'){
    samples <- output$samples.beta[,3]
  } else if (param == 'beta3'){
    samples <- output$samples.beta[,4]
  } else if (param == 'beta4'){
    samples <- output$samples.beta[,5]
  } else if (param == 'tau2'){
    samples <- output$samples.tau2
  } else if (param == 'tau.sq.1'){
    samples <- output$samples.tau2[,1]
  } else if (param == 'tau.sq.2'){
    samples <- output$samples.tau2[,2]
  } else if (param == 'sigma2'){
    samples <- output$samples.sigma2
  } else if (param == 'rho'){
    samples <- output$samples.rho
  } else if (param == 'rho.1'){
    samples <- output$samples.rho[1,,]
  } else if (param == 'rho.2'){
    samples <- output$samples.rho[2,,]
  } else if (param == 'beta0.loc'){
    samples <- output$samples.beta.l[,1]
  } else if (param == 'beta1.loc'){
    samples <- output$samples.beta.l[,2]
  } else if (param == 'beta2.loc'){
    samples <- output$samples.beta.l[,3]
  } else if (param == 'beta0.cond'){
    samples <- output$samples.beta.c[,1]
  } else if (param == 'beta1.cond'){
    samples <- output$samples.beta.c[,2]
  } else if (param == 'beta2.cond'){
    samples <- output$samples.beta.c[,3]
  } else if (param == 'alpha'){
    samples <- output$samples.alpha
  } else if (param == 'alpha.case'){
    samples <- output$samples.alpha.case
  } else if (param == 'alpha.control'){
    samples <- output$samples.alpha.control
  } else if (param == 'beta0.case'){
    samples <- output$samples.beta.case[,1]
  } else if (param == 'beta1.case'){
    samples <- output$samples.beta.case[,2]
  } else if (param == 'beta2.case'){
    samples <- output$samples.beta.case[,3]
  } else if (param == 'beta0.control'){
    samples <- output$samples.beta.control[,1]
  } else if (param == 'beta1.control'){
    samples <- output$samples.beta.control[,2]
  } else if (param == 'beta2.control'){
    samples <- output$samples.beta.case[,2]
  } 
  
  return(samples)
  
}


postSd <- function(param, output){
  
  samples <- getSamples(param, output)
  return(round(sd(samples), 3))
  
}


postMean <- function(param, output){
  
  samples <- getSamples(param, output)
  return(round(mean(samples), 3))
  
}


postMedian <- function(param, output){
  
  samples <- getSamples(param, output)
  return(round(median(samples), 3)) 
  
}


## arguments
#   output: output of a CAR model
#   param: name of the parameter
#   trueval: true parameter value
#   type: 'abs' or 'perc' (percent)
#   round: rounds bias to 3 digits
bias <- function(output, param, trueval, type='abs', round=TRUE){
  
  samples <- getSamples(param, output)
  bs <- mean(samples) - trueval
  
  if ((type) == 'perc'){
    bs <- 100*bs/trueval
  }
  
  if (round){
    bs <- round(bs, 3)
  }
  
  return(bs)
  
}


# creates traceplots from a list of results,
# i.e. outputs of either carBYM / carLeroux. 
# param specifies which parameter to plot
viewTraces <- function(results, param, line=NULL, titletype='ini'){
  
  for (r in results){
    
    samples <- getSamples(param, r)
    
    if (titletype=='ini'){
      title <- paste('beta0: ', round(r$beta.initial[1], 2),
                    'beta1:', round(r$beta.initial[2], 2),
                    'tau2:', round(r$tau2.initial, 2))
    } else {
      title <- param
    }
    
    plot(samples, type='l', main=title)
    if (!is.null(line)){
      abline(h=line, col='2')
    }
    
  }
  
}


## viewTracesInt
  # create traceplots from the output of the 
  # data integration CAR models
## arguments
  # results: list (output of data integration)
  # d: index for which dataset to visualize
  # title.type: either 'tune' or 'initial'
  # overall.title: overall title
viewTracesInt <- function(results, param, d=NULL, line=NULL, title.type='tune', overall.title=NULL){
  
  for (r in results){
    
    samples <- getSamples(param, output)
    
    if (title.type == 'tune'){
      title <- paste('tune beta:', r$tune.beta, '| tune phi:', r$tune.phi)
    } else if (title.type == 'pooled'){
      title <- param
    } else{
      taus <- sapply(r$tau2.initial, function(x){round(x,2)})
      title <- paste(c('beta0:', round(r$beta.initial[1], 2),
                     '| beta1:', round(r$beta.initial[2], 2),
                     '| tau2:', taus), collapse=" ")
    }
    
    plot(samples, type='l', main=title)
    if (!is.null(line)){
      abline(h=line, col='2')
    } 
    
  }
  
  if (!is.null(overall.title)){
    title(overall.title, outer=TRUE)
  }

}


# creates histograms from a list of results,
# i.e. outputs of either carBYM / carLeroux. 
# param specifies which parameter to plot
viewHists <- function(results, param, line=NULL){
  
  
  for (r in results){
    
    samples <- getSamples(param, r) 
    hist(samples, main=param)
    if (!is.null(line)){
      abline(v=line, col='2')
    }
    
  }
}


## viewTracesInt
  # create sample histograms from the output of the 
  # data integration CAR models
## arguments
  # results: list (output of data integration)
  # d: index for which dataset to visualize
  # title.type: either 'tune' or 'initial'
  # overall.title: overall title
viewHistsInt <- function(results, param, d=NULL, line=NULL, title.type='tune', overall.title=NULL){
  
  for (r in results){
    
    if (param == 'beta0'){
      samples <- r$samples.beta[,1]
    } else if (param == 'beta1'){
      samples <- r$samples.beta[,2]
    } else if (param == 'beta2'){
      samples <- r$samples.beta[,3]
    } else if (param == 'beta3'){
      samples <- r$samples.beta[,4]
    } else if (param == 'beta4'){
      samples <- r$samples.beta[,5]
    } else if (param == 'tau2'){
      samples <- r$samples.tau2[,d]
    } else if (param == 'sigma2'){
      samples <- r$samples.sigma2
    } else if (param == 'rho.1'){
      samples <- r$samples.rho[1,,]
    } else if (param == 'rho.2'){
      samples <- r$samples.rho[2,,]
    }
    
    if (title.type == 'tune'){
      title <- paste('tune beta:', r$tune.beta, '| tune phi:', r$tune.phi)
    } else if (title.type == 'pooled'){
      title <- param
    } else{
      taus <- sapply(r$tau2.initial, function(x){round(x,2)})
      title <- paste(c('beta0:', round(r$beta.initial[1], 2),
                       '| beta1:', round(r$beta.initial[2], 2),
                       '| tau2:', taus), collapse=" ")
    }
    
    hist(samples, main=title)
    if (!is.null(line)){
      abline(v=line, col='2')
    }
    
  }
  
  if (!is.null(overall.title)){
    title(overall.title, outer=TRUE)
  }
  
}


# exploratory function: plot any variable of the bioclim dataset (CA only)
# var: character name of worldclim variable to be plotted
# wc: worldclim RasterStack
plotwc <- function(varIndex, wc){
  
  datVar <- wc[[varIndex]]
  plot(datVar)
  bounds <- us_boundaries(type='state', resolution='low')
  plot(st_geometry(bounds), add=T)
  
}


# x: RasterLayer
plotstates <- function(x){
  
  plot(x)
  bounds <- us_boundaries(type='state', resolution='low')
  plot(st_geometry(bounds), add=T)
  
}


# arguments
  # x: raster stack to crop and mask
  # ext: extent
  # spdf: spatial polygons data frame defining regions to keep
cropmask <- function(x, ext, spdf){
  
  dat_sub <- crop(x, ext)
  return(mask(dat_sub, spdf))
  
}


# create a raster from scratch with the supplied values
inirast <- function(vals, rname='samp1'){
  
  r <- raster(nrow=nrow(dat_ca),
              ncol=ncol(dat_ca),
              xmn=-130,
              xmx=-100,
              ymn=25,
              ymx=50)
  r[] <- vals
  rm <- mask(r, ca)
  names(rm) <- c(rname)
  return(rm)
  
}


get_ca_data <- function(vars=NULL){
  
  dat_bio <- getData('worldclim', var='bio', res=10)
  if (!is.null(var)){
    dat_bio <- dat_bio[[vars]]
  }
  dat_sub <- crop(dat_bio, extent(-125, -110, 25, 50))
  us <- readRDS('Documents/research/dataInt/data/GADM_2.8_USA_adm1.rds')
  ca <- us[match("CALIFORNIA", toupper(us$NAME_1)), ]
  dat_ca <- mask(dat_sub, ca)
  return(dat_ca)
  
}

