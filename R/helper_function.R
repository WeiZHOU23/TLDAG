
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
  #Final.X <- as.matrix(data.frame(TrF.X))
  #return(Final.X)
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





## transform each binomial variable to binary ones
ConvBinomial <- function(conv.X,conv.N){
  conv.n <- length(conv.X)
  conv.vX <- NULL
  #conv.vX <- rep(0,conv.n*conv.N)
  #  conv.vX <- rep(0, 10)
  #con.vX[sample(1:(conv.N*conv.n), sum(conv.X), replace = F)]<- 1
  for (i in 1:conv.n) {
    conv.vX <- c(conv.vX, c(rep(0, conv.N-conv.X[i]), rep(1, conv.X[i])))
  }
  #conv.X <- conv.vX
  return(conv.vX)
}



## transform the poisson variable as the binomial is changed into binary
PoissonNew <- function(pois.X, pois.N){
  pois.n <- length(pois.X)
  poisNew <- NULL
  for (i in 1:pois.n) {
    poisNew <- c(poisNew, rep(pois.X[i], pois.N))
  }
  
  return(poisNew)
}