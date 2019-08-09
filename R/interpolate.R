library(np)
library(tictoc)
library(foreach)
library(doParallel)


interpolate_w_batched  <- function(samples, bws, r_train, r_pred, batch_size=500){
  
  samples_int <- c()
  
  for (batch in 1:(nrow(output_base$samples.w)/batch_size)){
    
    print(paste(batch))
    samples_batch <- interpolate_w(samples[(batch_size * (batch-1) + 1):(batch_size * batch),], bws, r_train, r_pred)
    samples_int <- rbind(samples_int, samples_batch)

  }
  
  return(samples_int)
  
}


interpolate_w <- function(samples, bws, r_train, r_pred, out_file=NULL){
  
  txdat <- data.frame(xyFromCell(r_train, (1:ncell(r_train))[!is.na(r_train[])]))
  x <- txdat[,1]
  y <- txdat[,2]
  df_new <- data.frame(xyFromCell(r_pred, (1:ncell(r_pred))[!is.na(r_pred[])]))
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  print(Sys.time())
  print(paste("interpolating", nrow(samples), "samples"))
  tic("smoothing")
  samples_pred <- foreach(i=1:nrow(samples), .combine=rbind) %dopar% {
    library(np)
    z <- samples[i,]
    model.np = npreg(bws=c(0.08, 0.09),
                     formula=z~x+y,
                     regtype="lc",
                     ckertype="gaussian")
    pred <- predict(model.np,
                    newdata=df_new)
    pred
  }
  toc()
  
  stopCluster(cl)
  
  if (!is.null(out_file)){
    save_output(samples_pred, out_file)
  }
  
  return(samples_pred)
  
}
