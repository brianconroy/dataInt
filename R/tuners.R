
## update.proposal1
#   updates the proposal standard deviation 
#   based on acceptance rate. update falls
#   within the min and max arguments
## arguments
#   accept: number of acceptances
#   total: number of total MH steps
#   sd: current proposal standard deviation
#   min: minimum desired acceptance rate
#   max: maximum desired acceptance rate
tuners.sd.1 <- function(accept, total, sd, min, max){
  
  
  rate <- 100*accept/total
  if(rate > max){
    sd <- sd + 0.1*sd
  }else if(rate < min)              
  {
    sd <- sd - 0.1*sd
  }else{}
  
  return(sd)
  
  
}


tuners.sd.2 <- function(accept, total, sd, min, max, sd.max){
  
  rate <- 100*accept/total
  if(rate > max){
    sd <- sd + 0.1*sd
    sd[which(sd>sd.max)] <- sd.max
  }else if(rate < min){
    sd <- sd - 0.1*sd
  }
  
  return(sd)
  
  
}


tuneMCMC <- function(data, mcmcFunction, tune.beta, tune.phi, tune.rho=list(NULL), 
                     tune.tau.prior, tune.beta.prior, n.sample, 
                     plot=FALSE, self.tune=FALSE, fix.rho=TRUE, rho.list=NULL){

  
  results <- list()
  counter <- 1
  for (i in 1:length(tune.beta)){
    for(j in 1:length(tune.phi)){
      for (k in 1:length(tune.tau.prior)){
        for (l in 1:length(tune.beta.prior)){
          for (m in tune.rho){
            print("***************")
            t.beta <- tune.beta[i]
            t.phi <- tune.phi[j]
            t.rho <- m
            prior.tau <- tune.tau.prior[[k]]
            prior.beta <- tune.beta.prior[[l]]
            print(paste("tuning beta prop:", t.beta, "phi prop:", t.phi))
            print(paste("phi prior params:", prior.tau[1], prior.tau[2]))
            print(paste("beta prior params:", prior.beta[1], prior.beta[2]))
            
            output <- mcmcFunction(data, t.beta, t.phi, t.rho, prior.var.beta=prior.beta, 
                                   prior.tau2=prior.tau, n.sample=n.sample, self.tune=self.tune, 
                                   fix.rho=fix.rho, rho.list=rho.list)
            output$tune.beta <- t.beta
            output$tune.phi <- t.phi
            output$tune.tau.prior <- prior.tau
            output$tune.beta.prior <- prior.beta
            
            results[[counter]] <- output
            counter <- counter + 1
            
            print("acceptance rates (beta, phi)")
            print(output$accept)
            print(paste("beta0", mean(output$samples.beta[,1])))
            print(paste("beta1", mean(output$samples.beta[,2])))
            print(paste("tau2", mean(output$samples.tau2)))
            
            if (plot){
              plot(output$samples.beta[,1], type='l', 
                   main=paste('beta0 | ', 'phi tune:', t.phi, '| beta tune:', t.beta,
                              '| beta prior:', prior.beta[1], '| tau2 prior:', prior.tau[1]))
              abline(h=beta0, col='2')
            
            }
          
          }
        }
      }
    }
  }
  
  return(results)
  
  
}


bestTune <- function(tuneResults, params, truevals){
  
  
  best.DIC <- Inf
  best.configs <- list()
  best.output <- list()
  for (r in tuneResults){
    r.summary <- summarize(r, params)
    if (r.summary$DIC < best.DIC){
      best.DIC <- r.summary$DIC
      best.configs <- list(tune.beta=r$tune.beta,
                           tune.phi=r$tune.phi)
      best.output <- r
    }
  }
  
  return(list(best.configs=best.configs, 
              best.DIC=best.DIC,
              output=best.output))
  
  
}


viewOutput <- function(output, type='thinned', truevals, overall.title=NULL, view.hist=TRUE, view.trace=TRUE){

  
  if (type=='preferential.alpha'){
    if (view.trace){
      for (n in names(truevals)){
        viewTraces(list(output), n, line=unlist(truevals[n]), titletype='param')
      }
    }
    if (view.hist){
      for (n in names(truevals)){
        viewHists(list(output), n, line=unlist(truevals[n]))
      }
    }
  }
  
  if (type=='thinned'){
    overall.main <- "integrated thinned model"
    par(mfrow=c(3,2))
    viewTraces(list(output), 'beta0', line=truevals$beta0, titletype='param')
    viewTraces(list(output), 'beta1', line=truevals$beta1, titletype='param')
    viewTraces(list(output), 'beta2', line=truevals$gamma.1, titletype='param')
    viewTraces(list(output), 'beta3', line=truevals$gamma.2, titletype='param')
  } else if (type=='preferential') {
    if (view.trace){
      viewTraces(list(output), 'beta1.loc', line=truevals$beta1.loc, titletype='param')
      viewTraces(list(output), 'beta0.cond', line=truevals$beta0.cond, titletype='param')
      viewTraces(list(output), 'beta1.cond', line=truevals$beta1.cond, titletype='param')  
    }
    if (view.hist){
      viewHists(list(output), 'beta1.loc', line=truevals$beta1.loc)
      viewHists(list(output), 'beta0.cond', line=truevals$beta0.cond)
      viewHists(list(output), 'beta1.cond', line=truevals$beta1.cond)  
    }
  } else if (type=='preferential.cc') {
    if (view.trace){
      viewTraces(list(output), 'beta1.loc', line=truevals$beta1.loc, titletype='param')
      viewTraces(list(output), 'alpha', line=truevals$alpha, titletype='param')
      viewTraces(list(output), 'beta0.case', line=truevals$beta0.case, titletype='param')
      viewTraces(list(output), 'beta1.case', line=truevals$beta1.case, titletype='param')
      viewTraces(list(output), 'beta0.control', line=truevals$beta0.control, titletype='param')
      viewTraces(list(output), 'beta1.control', line=truevals$beta1.control, titletype='param')
    }
    if (view.hist){
      viewHists(list(output), 'beta1.loc', line=truevals$beta1.loc)
      viewHists(list(output), 'alpha', line=truevals$alpha)
      viewHists(list(output), 'beta0.case', line=truevals$beta0.case)
      viewHists(list(output), 'beta1.case', line=truevals$beta1.case)
      viewHists(list(output), 'beta0.control', line=truevals$beta0.control)
      viewHists(list(output), 'beta1.control', line=truevals$beta1.control)
    }
  } else if (type=='preferential.cc.diff') {
    if (view.trace){
      viewTraces(list(output), 'beta1.loc', line=truevals$beta1.loc, titletype='param')
      viewTraces(list(output), 'alpha.case', line=truevals$alpha.case, titletype='param')
      viewTraces(list(output), 'alpha.control', line=truevals$alpha.control, titletype='param')
      viewTraces(list(output), 'beta0.case', line=truevals$beta0.case, titletype='param')
      viewTraces(list(output), 'beta1.case', line=truevals$beta1.case, titletype='param')
      viewTraces(list(output), 'beta0.control', line=truevals$beta0.control, titletype='param')
      viewTraces(list(output), 'beta1.control', line=truevals$beta1.control, titletype='param')
    }
    if (view.hist){
      viewHists(list(output), 'beta1.loc', line=truevals$beta1.loc)
      viewHists(list(output), 'alpha.case', line=truevals$alpha.case)
      viewHists(list(output), 'alpha.control', line=truevals$alpha.control)
      viewHists(list(output), 'beta0.case', line=truevals$beta0.case)
      viewHists(list(output), 'beta1.case', line=truevals$beta1.case)
      viewHists(list(output), 'beta0.control', line=truevals$beta0.control)
      viewHists(list(output), 'beta1.control', line=truevals$beta1.control)
    }
  }# else {
  #   overall.main <- "integrated model"
  #   par(mfrow=c(2,2))
  #   viewTraces(list(output), 'beta0', line=truevals$beta0, titletype='param')
  #   viewTraces(list(output), 'beta1', line=truevals$beta1, titletype='param')
  # }

  if (!(type %in% c('preferential', 'preferential.alpha', 
                    'preferential.cc', 'preferential.cc.diff'))){
    viewTraces(list(output), 'tau.sq.1', line=truevals$tau.sq.1, titletype='param')
    viewTraces(list(output), 'tau.sq.2', line=truevals$tau.sq.2, titletype='param')
  }
  
  if (!is.null(overall.title)){
    title(overall.title, outer=TRUE)
  }
  
  
}
