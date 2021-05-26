#############################################
###  helper function for the algorithm   ####
#############################################

### get the index of each node in the mixed DAG, including two distributions

CheckIndex <- function(CI.X, CI.N){
  CI.p <- dim(CI.X)[2]
  getNumber <- apply(CI.X, 2, unique)
  CI.IndexPois <- which(sapply(1:CI.p, function(i) any(unique(getNumber[[i]]) > CI.N)) == T)
  
  if(length(CI.IndexPois) == 0){
    CI.IndexBino <- seq(CI.p)
  }else{
    CI.IndexBino <- seq(CI.p)[-CI.IndexPois]
  }
  
  return(list(IndexBino = CI.IndexBino, IndexPois = CI.IndexPois))
}





## transform each binomial variable to binary ones
ConvBinomial <- function(conv.X, conv.N){
  conv.n <- length(conv.X)
  conv.vX <- NULL
  conv.list <- lapply(c(1:conv.n), function(i) c(rep(0, conv.N-conv.X[i]), rep(1, conv.X[i])))
  conv.vec <- unlist(conv.list)
  
  return(conv.vec)
}


## transfrom the square of each binomial variable to binary ones
Conv2one <- function(conv.X, conv.N){
  conv.n <- length(conv.X)
  conv.vX <- NULL
  conv.list <- lapply(c(1:conv.n), function(i) c(rep(0, conv.N^2-conv.X[i]), rep(1, conv.X[i])))
  conv.vec <- unlist(conv.list)
  
  return(conv.vX)
}



## transform the poisson variable as the binomial is changed into binary
PoissonNew <- function(pois.X, pois.N){
  pois.n <- length(pois.X)
  poisNew <- lapply(c(1:pois.N), function(i) rep(pois.X[i], pois.N))
  poisVec <- unlist(poisNew)
  
  return(poisVec)
}




PoissonNew2 <- function(pois.X, pois.N){
  pois.n <- length(pois.X)
  poisNew <- lapply(c(1:pois.n), function(i) rep(pois.X[i], pois.N))
  poisVec <- unlist(poisNew)
  
  return(poisVec)
}




## transform the poisson varianble with the square for the binomial ones
Poisson2rep <- function(pois.X, pois.N){
  pois.n <- length(pois.X)
  poisNew <- lapply(c(1:pois.N), function(i) rep(pois.X[i], pois.N^2))
  poisVec <- unlist(poisNew)
  
  return(poisNew)
}


### find the parents of the node 
potentialParents <- function(x, index){
  parents <- NULL
  T <- length(x)
  for (pp in (index+1):T) {
    NewParents <- x[[pp]]
    parents <- c(parents, NewParents)
  }
  return(parents)
} 





### find the direction in the binomial case, we need to deal with the binary variables
TransFactor <- function(TrF.X){
  TrF.p <- dim(TrF.X)[2]
  X.final <- sapply(1:TrF.p, function(i) as.factor(TrF.X[, i]))
  return(X.final)
}




##############################################
### some compariasion functions for later use
##############################################


hammingDistance <- function(G1, G2){
  allMistakeOne <- FALSE
  if(allMistakeOne){
    Gtmp <- (G1 + G2) %% 2
    Gtmp <- Gtmp + t(Gtmp)
    nrRevesals <- sum(Gtmp == 2)/2
    nrIncDel <- sum(Gtmp == 1)/2
    hammingDis <- nrRevesals + nrIncDel
    
  }else{
    hammingDis <- sum(abs(G1 - G2))
    hammingDis <- hammingDis - 0.5 *sum( G1 * t(G1) * (1-G2) * t(1-G2) + G2 * t(G2) * (1-G1) * t(1-G1))
  }
  return(hammingDis)
}


charaToNume <- function(pre.adj, p){
  # pre.adj is the the character
  
  nume.adj <- as.numeric(as.factor(pre.adj))
  
  parentsV <- nume.adj[1:(length(nume.adj)/2)] 
  childsV <- nume.adj[- (1:(length(nume.adj)/2))] 
  edge <- length(parentsV) 
  nume.final <- matrix(0, p, p)
  for (i in 1:edge) {
    parV <- parentsV[i]
    chiV <- childsV[i]
    nume.final[parV, chiV] <- 1
  }
  
  return(nume.final)
}





##  check distribution for one node in  layer_At function
To_dist <- function(node.index, IndexBino, IndexPois){
  
  if(node.index %in% IndexBino){
    left.index <- "bino"
  }else{
    left.index <- "pois"
  }
  
  return(left.index)
}

  
 

  