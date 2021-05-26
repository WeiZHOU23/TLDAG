################################################################################
###### To calculate kappa value to determine the tuning parameter \epsilon   ###
################################################################################


agree.twosets <- function(AG.set1, AG.set2, AG.left){
  p.AGtot<-length(AG.left)
  if(length(AG.set1)+length(AG.set2) ==0 || length(AG.set1)+length(AG.set2)==2*p.AGtot )
    AG.kap <- -1 
  else{
    n11 <- length(intersect(AG.set1, AG.set2)) 
    n12 <- length(setdiff(AG.set1,  intersect(AG.set1, AG.set2))) 
    n21 <- length(setdiff(AG.set2,  intersect(AG.set1, AG.set2)))
    n22 <- length(intersect(setdiff(AG.left, AG.set1), setdiff(AG.left, AG.set2)))
    
    
    AG.kap <- ( (n11+n22)/p.AGtot - ((n11+n12)*(n11+n21)+(n12+n22)*(n21+n22)) / (p.AGtot*p.AGtot) ) / 
      (1 -  ((n11+n12)*(n11+n21)+(n12+n22)*(n21+n22)) / (p.AGtot*p.AGtot) )
    
  }
  return(AG.kap)
}


#################################################
Kappa.cal <- function(KC.layer, KC.X, KC.n, KC.p, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.bino.N, KC.Distri = c("pois","mix"), KC.IndexPois, KC.IndexBino, KC.seed){
  
  KC.n1 <- as.integer(KC.n/2)
  KC.n2 <- KC.n - KC.n1
  KC.done <- unlist(KC.layer) ## The set S
  set.seed(KC.seed) #### Add seed here
  KC.index <- sample(1:KC.n, KC.n, replace = FALSE)
  data.perm <- KC.X[KC.index, ] 
  KC.data1 <- data.perm[1:KC.n1, ]
  KC.data2 <- data.perm[-(1:KC.n1), ]
  
  
  if(KC.Distri == "pois"){
    ## only for poission distribution
    
    if(length(KC.done) == 0){ 
      ## A_0 for roots
      
      KC.set1 <- layer_A0(KC.data1, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.bino.N, KC.Distri)
      KC.set2 <- layer_A0(KC.data2, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.bino.N, KC.Distri)
      
    }else if( length(KC.done) >=1 & length(KC.done) < (KC.p-1)){  
      ## layer A_1 to A_{T-1}
      
      KC.set1 <- layer_At(KC.layer, KC.data1, KC.data1, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.bino.N, KC.Distri, KC.IndexPois, KC.IndexBino)
      KC.set2 <- layer_At(KC.layer, KC.data2, KC.data1, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.bino.N, KC.Distri, KC.IndexPois, KC.IndexBino)
    }
  }else{
    ## for mixed distributions
    
    if(length(KC.IndexBino) != 0){
      KC.X.Binary <- matrix(0, KC.n*KC.bino.N, KC.p)
      KC.X.Binary[, KC.IndexBino] <- sapply(1:length(KC.IndexBino), function(i) ConvBinomial(KC.X[, KC.IndexBino[i]], KC.bino.N))
      
      if(length(KC.IndexPois) != 0){
        
        KC.X.Binary[, KC.IndexPois] <- sapply(1:length(KC.IndexPois), function(i) PoissonNew2(KC.X[, KC.IndexPois[i]], KC.bino.N))
      }
      
      KC.XBino1 <- KC.X.Binary[1:(KC.n1*KC.bino.N), ]
      KC.XBino2 <- KC.X.Binary[-(1:(KC.n1*KC.bino.N)), ]
    }else{

      KC.X.Binary <- KC.X
      KC.XBino1 <- KC.data1
      KC.XBino2 <- KC.data2
    }
    
    if(length(KC.done) == 0){ 
      ## A_0 for roots
      KC.set1 <- layer_A0(KC.data1, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.bino.N, KC.Distri)
      KC.set2 <- layer_A0(KC.data2, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.bino.N, KC.Distri)
      
    }else if( length(KC.done) >=1 & length(KC.done) < (KC.p-1)){  
      ## layer A_1 to A_{T-1}
      KC.set1 <- layer_At(KC.layer, KC.data1, KC.XBino1, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.bino.N, KC.Distri, KC.IndexPois, KC.IndexBino)
      KC.set2 <- layer_At(KC.layer, KC.data2, KC.XBino2, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.bino.N, KC.Distri, KC.IndexPois, KC.IndexBino)
      
    }
  }
  
  return(list(KC.set1, KC.set2))
}



#################################################
Kapp.compute<-function(kap.re1, kap.re2, kap.lam, kap.left){
  kap.value<-0
  
  kap.re1[ which(abs(kap.re1 - 1) <= kap.lam)] <- 1
  kap.re1[ which(abs(kap.re1 - 1) > kap.lam)] <- 0
  kap.re2[ which(abs(kap.re2 - 1) <= kap.lam)] <- 1
  kap.re2[ which(abs(kap.re2 - 1) > kap.lam)] <- 0
  kap.value<-sum(sapply(1:nrow(kap.re1), function(i) agree.twosets(which(kap.re1[i,]==1),
                                                                   which(kap.re2[i,]==1), kap.left)))
  
  return(kap.value)
}



#################################################
Kappa <- function(Ka.layer, Ka.X, Ka.grid.t, Ka.beta1.pois, Ka.beta2.pois, Ka.beta1.bino, Ka.beta2.bino, Ka.N, Ka.Dist = c("pois","mix"), Ka.IndexPois, Ka.IndexBino){
  
  Ka.n <- dim(Ka.X)[1]
  Ka.p <- dim(Ka.X)[2]
  Ka.done <- unlist(Ka.layer)
  # print(Ka.done)
  if(is.null(Ka.done) == T){
    Ka.left <- seq(Ka.p)
  }else{
    Ka.left <- seq(Ka.p)[-Ka.done]
  }
  
  Kappa.temp <- sapply(c(1,2,3,4,5), Kappa.cal, KC.layer=Ka.layer, KC.X=Ka.X, KC.n=Ka.n, KC.p=Ka.p, 
                     KC.beta1.pois = Ka.beta1.pois, KC.beta2.pois = Ka.beta2.pois, KC.beta1.bino = Ka.beta1.bino, KC.beta2.bino = Ka.beta2.bino, KC.bino.N = Ka.N,
                     KC.Distri = Ka.Dist, KC.IndexPois = Ka.IndexPois, KC.IndexBino = Ka.IndexBino)
  kappa.temp1 <- matrix(c(Kappa.temp[[1]],Kappa.temp[[3]],Kappa.temp[[5]],Kappa.temp[[7]],Kappa.temp[[9]]), 5, length(Kappa.temp[[1]]), byrow = T)
  kappa.temp2 <- matrix(c(Kappa.temp[[2]],Kappa.temp[[4]],Kappa.temp[[6]],Kappa.temp[[8]],Kappa.temp[[10]]), 5, length(Kappa.temp[[2]]), byrow = T)
  Ka.result <- sapply(Ka.grid.t, Kapp.compute, kap.re1=kappa.temp1, kap.re2=kappa.temp2, kap.left=Ka.left)
  return(Ka.result/5)
}


#################################################
TuningPara <- function(Tune.layer, Tune.X, Tune.grid.t, Tune.alpha, Tune.beta1.pois, Tune.beta2.pois, Tune.beta1.bino, Tune.beta2.bino, Tune.N, Tune.Dist = c("pois","mix"), Tune.IndexPois, Tune.IndexBino){

  kappa.vec <- Kappa(Tune.layer, Tune.X, Tune.grid.t, Tune.beta1.pois, Tune.beta2.pois,Tune.beta1.bino, Tune.beta2.bino, Tune.N, Tune.Dist, Tune.IndexPois, Tune.IndexBino)
  
  if(any(is.na(kappa.vec))){
    kappa.vec <- kappa.vec[-which(is.nan(kappa.vec))]
  }
  if(any(is.infinite(kappa.vec))){
    kappa.vec <- kappa.vec[-which(is.infinite(kappa.vec))]
  }
  
  if(max(kappa.vec) == 0){
    # the maximum of kappa.vec is 0
    i.opt.kappa <- Tune.grid.t[which( kappa.vec >=  max(kappa.vec)*(1-Tune.alpha) )[1]]
  }else{
    i.opt.kappa <- Tune.grid.t[which( kappa.vec / max(kappa.vec) >= (1-Tune.alpha) )[1]]
  }
  result <- list(alpha = Tune.alpha, kappa.values = kappa.vec, thres = i.opt.kappa)
  
  return(result)
}



