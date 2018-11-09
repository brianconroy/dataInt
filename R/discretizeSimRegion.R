
## discretize simulation region
# return a list containing:
# W: adjacency matrix

# ToDO: remove this: vals: values of the raster

# raster: RasterBrick of covariates over the study region,
# nNeighbors: number of neighbors for each grid in the study region

## arguments
# param caWin: window polygon boundary for california
# param caWc: RasterBrick cropped and masked around california
# param factor: aggregation factor to discretize the study region
discretizeSimRegion <- function(caWin, caWc, factor=4){
  
  
  # reduce raster resolution
  # i.e. discretize the study region
  caWc.disc <- aggregate(caWc, fact=factor)
  vals <- values(caWc.disc[[1]])
  
  
  # how to count number of neighbors
  # iterate over matrix of raster vals
  valMatrix <- matrix(vals, 
                      nrow=nrow(caWc.disc), 
                      ncol=ncol(caWc.disc),
                      byrow=TRUE)
  nNeighbors <- c()
  for(i in 1:nrow(valMatrix)){
    for(j in 1:ncol(valMatrix)){
      pos <- c('left', 'right', 'bottom', 'top')
      nbhInds <- sapply(pos, function(x) {getValIndicator(valMatrix, i, j, x)})
      numNbhs <- sum(nbhInds)
      nNeighbors <- c(nNeighbors, numNbhs)
    }
  }
  
  
  # discard values of islands
  inds <- !is.na(vals)
  inds.0 <- nNeighbors == 0
  inds.discard <- inds & inds.0
  vals[inds.discard] <- NA
  values(caWc.disc)[inds.discard] <- NA

  
  nbhRast <- caWc.disc
  nbhRast[] <- nNeighbors
  
  
  # Create the adjacency matrix
  # same concept as above
  # get neighbor indices
  # only non NAs
  getNbhIndices <- function(m, i, j){
    pos <- c('left', 'right', 'bottom', 'top')
    nbhInds <- sapply(pos, function(x) {getValIndicator(m, i, j, x)})
    posNbh <- names(nbhInds[nbhInds])
    inds <- lapply(posNbh, function(x){position2Index(x, i, j)})
    return(inds)
  }
  
  
  index2CellVal <- function(m, ij){
    i <- ij[1]
    j <- ij[2]
    preceding <- ncol(m)*(i-1)
    return(preceding + j)
  }
  
  
  Wna <- matrix(rep(0, ncell(caWc.disc)^2), ncol=ncell(caWc.disc))
  for (i in 1:nrow(caWc.disc)){
    for (j in 1:ncol(caWc.disc)){
      currCell <- index2CellVal(caWc.disc, c(i, j))
      nbs <- getNbhIndices(valMatrix, i, j)
      cellVals <- sapply(nbs, function(x){index2CellVal(caWc.disc, x)})
      for (v in cellVals){
        Wna[currCell, v] <- 1
      }
    }
  }
  
  
  # which cell indices aren't NA?
  indRetain <- c(1:length(vals))[!is.na(vals)]
  W <- matrix(rep(0, length(indRetain)^2), ncol=length(indRetain))
  for(i in 1:length(indRetain)){
    newRow <- Wna[indRetain[i], c(indRetain)]
    W[i,] <- newRow
  }
  
  
  sim <- list()
  sim$raster <- caWc.disc
  sim$W <- W
  sim$nNeighbors <- nNeighbors[!is.na(values(caWc.disc))[,1]]
  return(sim)
  
  
}


getVal <- function(m, i, j){
  if (j > ncol(m) || j < 1){
    return(NaN)
  } else if (i > nrow(m) || i < 1){
    return(NaN)
  } else {
    return(m[i, j])
  }
}


# get the row, column indices at some position
# (right, left, etc) with respect to a given
# position (i, j)
position2Index <- function(position, i, j){
  if (position == 'right'){
    r <- i
    c <- j + 1
  } else if (position == 'left'){
    r <- i
    c <- j - 1
  } else if (position == 'bottom right'){
    r <- i + 1
    c <- j + 1
  } else if (position == 'bottom left'){
    r <- i + 1
    c <- j - 1
  } else if (position == 'bottom'){
    r <- i + 1
    c <- j
  } else if (position == 'top'){
    r <- i - 1
    c <- j
  } else if (position == 'top right'){
    r <- i - 1
    c <- j + 1
  } else if (position == 'top left'){
    r <- i - 1
    c <- j - 1
  }
  return(c(r, c))
}


getValIndicator <- function(m, i, j, position){
  inds <- position2Index(position, i, j)
  r <- inds[1]
  c <- inds[2]
  targetVal <- getVal(m, r, c)
  return(!is.na(targetVal))
}
