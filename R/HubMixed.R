###################################
###  TLDAG for hub mixed graph  ###
###################################




source("..GDS/util_DAGs/computeScoreSEMGauss.R")
source("../GDS/inferDAG/gesWrap.R")
source("../R/DLiNGAM_p.R")
sourceCpp('direct_lingam_funcs.cpp')
library(bnlearn)
library(pcalg)
library(Rcpp)
library(RcppArmadillo)





#################################################
## Calculate the Kappa coefficient of Two Sets
#################################################



agree.twosets <- function(AG.set1, AG.set2, AG.left){
  p.AGtot<-length(AG.left)
  if(length(AG.set1)+length(AG.set2) ==0 || length(AG.set1)+length(AG.set2)==2*p.AGtot )
    AG.kap <- -1 
  else{
    n11 <- length(intersect(AG.set1, AG.set2)) #### 
    n12 <- length(setdiff(AG.set1,  intersect(AG.set1, AG.set2))) 
    n21 <- length(setdiff(AG.set2,  intersect(AG.set1, AG.set2)))
    n22 <- length(intersect(setdiff(AG.left, AG.set1), setdiff(AG.left, AG.set2)))
    
    
    AG.kap <- ( (n11+n22)/p.AGtot - ((n11+n12)*(n11+n21)+(n12+n22)*(n21+n22)) / (p.AGtot*p.AGtot) ) / 
      (1 -  ((n11+n12)*(n11+n21)+(n12+n22)*(n21+n22)) / (p.AGtot*p.AGtot) )
    
  }
  return(AG.kap)
}


#################################################
Kappa.cal <- function(KC.layer, KC.X, KC.n, KC.p, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.N, KC.index){
  KC.n1 <- as.integer(KC.n/2)
  KC.n2 <- KC.n - KC.n1
  KC.done <- unlist(KC.layer) ## The set S
  set.seed(KC.index) #### Add seed here
  KC.index <- sample(1:KC.n, KC.n, replace = FALSE)
  data.perm <- KC.X[KC.index, ] 
  KC.data1 <- data.perm[1:KC.n1, ]
  KC.data2 <- data.perm[-(1:KC.n1), ]
  
  KC.IndexPois <- CheckIndex(KC.X, KC.N)$IndexPois
  KC.IndexBino <- CheckIndex(KC.X, KC.N)$IndexBino
  
  if(length(KC.IndexBino) != 0){
    KC.X.Binary <- matrix(0, KC.n*KC.N, KC.p)
    KC.X.Binary[, KC.IndexBino] <- sapply(1:length(KC.IndexBino), function(i) ConvBinomial(KC.X[, KC.IndexBino[i]], KC.N))
    
    if(length(KC.IndexPois) != 0){
      KC.X.Binary[, KC.IndexPois] <- sapply(1:length(KC.IndexPois), function(i) PoissonNew(KC.X[, KC.IndexPois[i]], KC.N))
    }
  }else{
    
    # all nodes are poisson nodes
    KC.X.Binary <- matrix(0, KC.n*KC.N, KC.p)
  }
  
  
  KC.XBino1 <- KC.X.Binary[1:(KC.n1*KC.N), ]
  KC.XBino2 <- KC.X.Binary[-(1:(KC.n1*KC.N)), ]
  
  if(length(KC.done) == 0){ 
    ## A_0 for roots
    KC.set1 <- layer_A0(KC.data1, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.N)
    KC.set2 <- layer_A0(KC.data2, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.N)
    
  }else if( length(KC.done) >=1 & length(KC.done) < (KC.p-1)){ 
    ## layer A_1 to A_{T-1}
    KC.set1 <- layer_At(KC.layer, KC.data1, KC.XBino1, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.N)
    KC.set2 <- layer_At(KC.layer, KC.data2, KC.XBino2, KC.beta1.pois, KC.beta2.pois, KC.beta1.bino, KC.beta2.bino, KC.N)
    
  }
  return(list(KC.set1, KC.set2))
}


#################################################
Kapp.compute<-function(kap.re1, kap.re2, kap.lam, kap.left){
  kap.value<-0
  
  kap.re1[ which(abs(kap.re1 - 1) <= kap.lam)]<-1
  kap.re1[ which(abs(kap.re1 - 1) > kap.lam)]<-0
  kap.re2[ which(abs(kap.re2 - 1) <= kap.lam)]<-1
  kap.re2[ which(abs(kap.re2 - 1) > kap.lam)]<-0
  kap.value<-sum(sapply(1:nrow(kap.re1), function(i) agree.twosets(which(kap.re1[i,]==1),
                                                                   which(kap.re2[i,]==1), kap.left)))
  
  return(kap.value)
}



#################################################
Kappa <- function(Ka.layer, Ka.X, Ka.grid.t, Ka.beta1.pois, Ka.beta2.pois, Ka.beta1.bino, Ka.beta2.bino, Ka.N){
  Ka.n <- dim(Ka.X)[1]
  Ka.p <- dim(Ka.X)[2]
  Ka.done<-unlist(Ka.layer)
  # print(Ka.done)
  if(is.null(Ka.done)==T){
    Ka.left<-seq(Ka.p)
  }else{
    Ka.left<-seq(Ka.p)[-Ka.done]
  }
  Kappa.temp<-sapply(c(1,2,3,4,5), Kappa.cal, KC.layer=Ka.layer, KC.X=Ka.X, KC.n=Ka.n, KC.p=Ka.p, 
                     KC.beta1.pois=Ka.beta1.pois, KC.beta2.pois=Ka.beta2.pois, KC.beta1.bino=Ka.beta1.bino, KC.beta2.bino=Ka.beta2.bino, 
                     KC.N=Ka.N)
  kappa.temp1<-matrix(c(Kappa.temp[[1]],Kappa.temp[[3]],Kappa.temp[[5]],Kappa.temp[[7]],Kappa.temp[[9]]),5,length(Kappa.temp[[1]]),byrow = T)
  kappa.temp2<-matrix(c(Kappa.temp[[2]],Kappa.temp[[4]],Kappa.temp[[6]],Kappa.temp[[8]],Kappa.temp[[10]]),5,length(Kappa.temp[[2]]),byrow = T)
  Ka.result<-sapply(Ka.grid.t, Kapp.compute, kap.re1=kappa.temp1, kap.re2=kappa.temp2, kap.left=Ka.left)
  return(Ka.result/5)
}


#################################################
TuningPara <- function(Tune.layer, Tune.X, Tune.grid.t, Tune.alpha, Tune.beta1.pois, Tune.beta2.pois, Tune.beta1.bino, Tune.beta2.bino, Tune.N){
  kappa.vec <- Kappa(Tune.layer, Tune.X, Tune.grid.t, Tune.beta1.pois, Tune.beta2.pois,Tune.beta1.bino, Tune.beta2.bino, Tune.N)
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



#############################################################################
## get the index of each node in the mixed DAG, including two distributions

CheckIndex <- function(CI.X, CI.N){
  CI.p <- dim(CI.X)[2]
  getNumber <- apply(CI.X, 2, unique)
  #getNumberTF <- sapply(1:p, function(i) unique(getNumber[[i]]) <= N)   ## binomial node only takes values 0/1 
  #IndexBino <- which(getNumberTF==TRUE)
  CI.IndexPois <- which(sapply(1:CI.p, function(i) any(unique(getNumber[[i]]) > CI.N)) == T)
  
  if(length(CI.IndexPois) == 0){
    CI.IndexBino <- seq(CI.p)
  }else{
    CI.IndexBino <- seq(CI.p)[-CI.IndexPois]
  }
  
  return(list(IndexBino = CI.IndexBino, IndexPois = CI.IndexPois))
}

## calculate the unconditional ratio in the layer A_0
unconRatio <- function(uR.X, uR.beta_1, uR.beta_2){
  uR.n <- dim(uR.X)[1]
  uR.p <- dim(uR.X)[2]
  uR.T <- matrix(0, uR.n, uR.p)
  
  uR.w <- (uR.beta_1 + uR.beta_2 * colMeans(uR.X))^(-1)
  uR.w2 <- diag(uR.w)
  if(uR.p == 1){
    uR.T <- uR.X * uR.w
  }else{
    uR.T <- uR.X %*% uR.w2
  }
  
  uR.ratio <- c()
  uR.ratio <- apply(uR.T, 2, var) / colMeans(uR.T)
  
  #uR.adjustratio <- uR.ratio / min(uR.ratio)
  return(uR.ratio)
}

##################################
## roots for the layer A_0     ###
##################################

layer_A0 <- function(A0.X, A0.beta1.pois, A0.beta2.pois, A0.beta1.bino, A0.beta2.bino, A0.N){
  A0.n <- dim(A0.X)[1]
  A0.p <- dim(A0.X)[2]
  
  A0.IndexPois <- CheckIndex(A0.X, A0.N)$IndexPois
  A0.IndexBino <- CheckIndex(A0.X, A0.N)$IndexBino
  
  
  A0.XPois <- as.matrix(A0.X[, A0.IndexPois])
  A0.XBino <- as.matrix(A0.X[, A0.IndexBino])
  
  
  A0.ratio <- rep(0, A0.p)
  A0.ratio[A0.IndexPois] <- unconRatio(A0.XPois, uR.beta_1=A0.beta1.pois, uR.beta_2=A0.beta2.pois)
  A0.ratio[A0.IndexBino] <- unconRatio(A0.XBino, uR.beta_1=A0.beta1.bino, uR.beta_2=A0.beta2.bino)
  
  A0.final.ratio <- A0.ratio / min(A0.ratio)
  return(A0.final.ratio)
}



#####################################################
## calculate ratio for the nodes in the layer A_t ###
#####################################################

layer_At <-function(At.layer, At.X, At.XBino, At.beta1.pois, At.beta2.pois, At.beta1.bino, At.beta2.bino=-1/At.N, At.N){
  # AT.R is the the number of failness in the negative binomial distribution
  At.n1 <- dim(At.X)[1]
  At.p <- dim(At.X)[2]
  At.Sdone <- unlist(At.layer)
  At.left <- seq(At.p)[-At.Sdone]
  At.n.left <- length(At.left)
  #At.con.mean = At.con.var = c();
  At.IndexBino  <- CheckIndex(At.X, At.N)$IndexBino
  At.IndexPois  <- CheckIndex(At.X, At.N)$IndexPois
  
  At.con.var = At.con.mean = c();
  for (k in 1:At.n.left) {
    if(At.left[k] %in% At.IndexBino){
      
      glmfit <- glm(At.XBino[, At.left][, k] ~ At.XBino[, At.Sdone], family = "binomial", control=list(maxit=100))
      glmfit.coef <- as.matrix(glmfit$coefficients)
      
      if(any(is.na(glmfit.coef))){
        na.index <- which(is.na(glmfit.coef) == TRUE)
        glmfit.coef[na.index] <- 0
      }
      
      glm_St <- cbind(At.XBino[, At.left][, k], At.XBino[, At.Sdone])
      exp.term <- exp(glm_St %*% glmfit.coef)
      con.prob <- mean(exp.term^(At.XBino[, At.left][, k]) * (factorial(1)/(factorial(At.XBino[, At.left][, k])*factorial(1-At.XBino[, At.left][, k])))
                       / (1 + exp.term)^(At.XBino[, At.left][, k]))
      con.mean.X <-  con.prob * sum(At.X[, At.left][, k])
      w.Tt <- (At.beta1.bino + At.beta2.bino * con.mean.X)^(-1)
      At.Tt <- At.X[, At.left][, k] * w.Tt
      At.con.mean[k] <- con.prob * w.Tt * sum(At.Tt)
      
      sec.moment <- con.prob * w.Tt * sum(At.Tt^2)
      At.con.var[k] <- sec.moment - At.con.mean[k]^2
      
    }else{
      glmfit <- glm(At.X[, At.left][, k] ~ At.X[, At.Sdone], family = "poisson", control=list(maxit=100))
      con.mean1Pois <- mean(exp(predict(glmfit)))
      wPois <- (At.beta1.pois + At.beta2.pois * con.mean1Pois)^(-1)
      TtPois <- At.X[, At.left][, k] * wPois
      Tt2Pois <- (At.X[, At.left][, k] * wPois)^2
      
      fitTPois <- glm(TtPois ~ At.X[, At.Sdone], family = "poisson", control=list(maxit=100)) 
      At.con.meanT <- mean(exp(predict(fitTPois)))
      
      At.con.mean[k] <- mean(TtPois)
      
      fitT2Pois <- glm(Tt2Pois ~ At.X[, At.Sdone], family = "poisson", control=list(maxit=100)) 
      sec.moment <- mean(exp(predict(fitT2Pois)))
      At.con.var[k] <- sec.moment - At.con.meanT^2
      
    }
  }
  
  At.ratio <- At.con.var / At.con.mean
  At.final.ratio <- At.ratio / min(At.ratio)
  return(At.final.ratio) 
}





#################################################
###   The proposed method, denoted by TLDAG
#################################################

EDAG <- function(ED.X, ED.grid.t, ED.beta1.pois, ED.beta2.pois, ED.beta1.bino, ED.beta2.bino, ED.N){
  ptm <- proc.time() 
  ED.n <- dim(ED.X)[1]
  ED.p <- dim(ED.X)[2]
  ED.TF <- matrix(0, ED.n, ED.p)
  ED.done <- NULL
  ED.layer <- list()
  ED.ratio.t <- list()
  ED.thres <- list()
  ED.final.layer<-list()
  
  
  
  ED.IndexPois <- CheckIndex(ED.X, ED.N)$IndexPois
  ED.IndexBino <- CheckIndex(ED.X, ED.N)$IndexBino
  
  if(length(ED.IndexBino) != 0){
    ED.X.Binary <- matrix(0, ED.n*ED.N, ED.p)
    ED.X.Binary[, ED.IndexBino] <- sapply(1:length(ED.IndexBino), function(i) ConvBinomial(ED.X[, ED.IndexBino[i]], ED.N))
    
    if(length(ED.IndexPois) != 0){
      ED.X.Binary[, ED.IndexPois] <- sapply(1:length(ED.IndexPois), function(i) PoissonNew(ED.X[, ED.IndexPois[i]], ED.N))
    }
  }else{
    
    ED.X.Binary <- ED.X
  }
  
  
  
  ##################################################################################  
  # Step1: layer and ordering
  # roots layer A_0
  ############## Node identification ##################
  
  ED.A0.thresh <- TuningPara(ED.layer, ED.X, ED.grid.t, 0.1, ED.beta1.pois, ED.beta2.pois, ED.beta1.bino, ED.beta2.bino, ED.N)$thres
  ED.ratio.final<-layer_A0(ED.X, ED.beta1.pois, ED.beta2.pois, ED.beta1.bino, ED.beta2.bino, ED.N)
  ED.layer[[1]] <- which(abs(ED.ratio.final - 1) <= ED.A0.thresh)
  ED.thres[[1]] <- ED.A0.thresh
  ED.ratio.t[[1]]<-ED.ratio.final
  
  #######################################################################################
  # identify nodes in layer A_t
  if( length(unlist(ED.layer)) == ED.p){
    EDS.result <- matrix(0, ED.p, ED.p)
    EDS.layernodes <- unlist(ED.layer)
  }else{
    
    cont <- TRUE
    while( cont ){
      for (lay in 2:ED.p) {
        if( length(unlist(ED.layer)) < ED.p-1){   
          
          ED.left <- seq(ED.p)[-unlist(ED.layer)]
          print(lay)
          ####################################################################
          ED.AT.thresh <- TuningPara(ED.layer, ED.X, ED.grid.t, 0.1, ED.beta1.pois, ED.beta2.pois, ED.beta1.bino, ED.beta2.bino, ED.N)$thres
          ED.thres[[lay]] <- ED.AT.thresh
          ED.AT.ratio.final<-layer_At(ED.layer, ED.X, ED.X.Binary, ED.beta1.pois, ED.beta2.pois, ED.beta1.bino, ED.beta2.bino, ED.N)
          ED.layer[[lay]] <- ED.left[which(abs(ED.AT.ratio.final - 1) <= ED.AT.thresh)]
          ED.ratio.t[[lay]]<-ED.AT.ratio.final
          #########################################################################
          
        }else if(length(unlist(ED.layer)) == ED.p-1){
          ## sink node is in A_T 
          ED.left <- seq(ED.p)[-unlist(ED.layer)]
          ED.n.left <- length(ED.left)
          ED.layer[[length(ED.layer) + 1]] <- ED.left
          
        }else (length(unlist(ED.layer)) >= ED.p)
        cont <- FALSE
        #break
      }
    }
    
    ED.final.layer <- list(ED.layer, ED.ratio.t, ED.thres)  # layer[[1]] indicates the topological layer
    
    
    ###############################################
    ##Step2:Sparse Estimation for Direction Edges
    ###############################################
    
    
    EDS.layernodes <- rev(ED.final.layer[[1]]) ### Selected layer structure with nodes in each layer
    EDS.layernum  <- length(ED.final.layer[[1]]) ### number of layer in the DAG
    EDS.result <- matrix(0, ED.p, ED.p)
    
    
    for (ii in 1:(EDS.layernum-1)) {   
      for (k in 1:length(EDS.layernodes[[ii]])) {
        EDS.now <- EDS.layernodes[[ii]][k]  ## each element in each layer ### the node as a child
        ## now turn to reconstruct the directed structure among SE.now and all the nodes in the upper layers
        
        #for (j in (ii+1):EDS.layernum) {
        
        #EDS.this <- EDS.layernodes[[j]]  ## each layer ### should be the nodes in the upper layer
        EDS.this <- potentialParents(EDS.layernodes, ii)
        
        if(length(EDS.this) > 1){
          
          if(EDS.now %in% ED.IndexBino){
            EDS.nfold <- 5
            EDS.fold0 <- sample.int(sum(ED.X.Binary[, EDS.now]==0)) %% EDS.nfold
            EDS.fold1 <- sample.int(sum(ED.X.Binary[, EDS.now]==1)) %% EDS.nfold
            EDS.foldid <- numeric(length(ED.X.Binary[, EDS.now]))
            EDS.foldid[ED.X.Binary[, EDS.now]==0] <- EDS.fold0
            EDS.foldid[ED.X.Binary[, EDS.now]==1] <- EDS.fold1
            EDS.foldid <- EDS.foldid + 1
            EDS.fraction <- table(ED.X.Binary[, EDS.now])/length(ED.X.Binary[, EDS.now])
            EDS.weights <- 1 - EDS.fraction[as.character(ED.X.Binary[, EDS.now])]
            
            EDS.lassom <- glmnet::cv.glmnet(ED.X.Binary[, EDS.this], ED.X.Binary[, EDS.now], family = "binomial", foldid = EDS.foldid, weights = EDS.weights )
            EDS.bfit <- coefficients(EDS.lassom)[-1]
            
          }else{
            EDS.lassomPois <- glmnet::cv.glmnet(ED.X[, EDS.this], ED.X[, EDS.now], family = "poisson", nfolds = 5 )
            EDS.bfit <- coefficients(EDS.lassomPois)[-1]
            
          }
          
          for (mm in 1:length(EDS.this)) {
            if(EDS.bfit[mm] != 0)
              EDS.result[EDS.this[mm], EDS.now] <- 1
          }
          
        }else{
          if(EDS.now %in% ED.IndexBino){
            EDS.lmod <- glm(ED.X.Binary[, EDS.now] ~ ED.X.Binary[, EDS.this], family = "binomial", control=list(maxit=100))
            EDS.bfit <- as.vector(EDS.lmod$coefficients[2])
          }else{
            lmodPois <- glm(ED.X[, EDS.now] ~ ED.X[, EDS.this], family = "poisson", control=list(maxit=100))
            EDS.bfit <- as.vector(lmodPois$coefficients[2])
          }
          
          if(EDS.bfit != 0)
            EDS.result[EDS.this, EDS.now] <- 1
        }
      }
    }   
  }
  #elapse <- proc.time() - ptm
  return(list(adj = EDS.result, OrderLayer = rev(EDS.layernodes)))
}



######################################################
### the ODS algorithm in Park and Raskutti (2018)  ###
######################################################


ODSGLMLasso <- function(X, beta1.pois, beta2.pois, beta1.bino, beta2.bino, N){
  ptm<-proc.time()
  n <- dim(X)[1]
  p <- dim(X)[2]
  IndexBino <- CheckIndex(X, N)$IndexBino
  IndexPois <- CheckIndex(X, N)$IndexPois
  
  XBino <- as.matrix(X[, IndexBino])
  XPois <- as.matrix(X[, IndexPois])
  
  
  if(length(IndexBino) != 0){
    X.all <- matrix(0, n*N, p)
    X.all[, IndexBino] <- sapply(1:length(IndexBino), function(i) ConvBinomial(X[, IndexBino[i]], N))
    
    if(length(IndexPois) != 0){
      X.all[, IndexPois] <- sapply(1:length(IndexPois), function(i) PoissonNew(X[, IndexPois[i]], N))
    }
  }
  
  #########################
  #Step1: NormalizedGraph
  #########################
  
  NormalizedGraph <- matrix(0, p, p)
  for (i in 1:p) {
    left <- seq(p)[-i]
    
    if(i %in% IndexBino){
      nfold <- 5
      fold0 <- sample.int(sum(X.all[, i]==0)) %% nfold
      fold1 <- sample.int(sum(X.all[, i]==1)) %% nfold
      foldid <- numeric(length(X.all[, i]))
      foldid[X.all[, i]==0] <- fold0
      foldid[X.all[, i]==1] <- fold1
      foldid <- foldid + 1
      fraction <- table(X.all[, i])/length(X.all[, i])
      weights <- 1 - fraction[as.character(X.all[, i])]
      
      lassomBino <- glmnet::cv.glmnet(X.all[, left], X.all[, i], family = "binomial", foldid = foldid, weights = weights)
      bfit <- coefficients(lassomBino)[-1]
      
    }else{
      lassomPois <- glmnet::cv.glmnet(X[, left], X[, i], family = "poisson", alpha = 1, nfolds = 5)
      bfit <- coefficients(lassomPois)[-1]
      
    }
    
    for (j in 1:length(bfit)) {
      if(bfit[j] != 0)
        NormalizedGraph[i, left[j]] <- 1
    }
  }
  ## return(NormalizedGraph)
  
  ###################
  #Step2: ordering
  ###################
  
  ODS.order <- c()
  ratioBino <- RatioODS(XBino, beta_1 = beta1.bino, beta_2 = beta2.bino)
  ratioPois <- RatioODS(XPois, beta_1 = beta1.pois, beta_2 = beta2.pois)
  
  S1 <- rep(0,p)
  S1[IndexBino] <- ratioBino
  S1[IndexPois] <- ratioPois
  ODS.order[1] <- seq(p)[which.min(S1)]
  #ODS.left[[1]] <- ODS.nodes[-which.min(S1)]
  
  # m-th element of the ordering m=2 to p-1
  # find condidate parents
  left.order <- seq(p)[-ODS.order]
  norgraph <- NormalizedGraph
  #dat <- X[, left.order]
  #kn <- c(left.order[1]:left.order[length(left.order)])
  C <- list()
  #S <- matrix(0, p-2, length(left.order))
  c0 <- 0.005
  
  # calculate scores
  for (m in 2:(p-1)) {
    ODS.order <- ODS.order  #### renew the ordering
    left.order <- seq(p)[-ODS.order]
    S <- c()
    
    for (k1 in 1:length(left.order) ) { # k for the original variable index
      k <- left.order[k1]
      neigh <- union(which( norgraph[k, ] != 0), which( norgraph[, k] != 0))   #k'neighborhood 
      C[[k]] <- intersect(neigh, ODS.order)   # given parents
      kk <- which(left.order == k)
      
      # calculate conditinal mean/variance  
      if(length(C[[k]]) ==0 ){  
        
        if(k %in% IndexBino){
          S[kk] <- RatioODS(as.matrix(X[, k]), beta_1 = beta1.bino, beta_2 = beta2.bino)
        }else{
          S[kk] <- RatioODS(as.matrix(X[, k]), beta_1 = beta1.pois, beta_2 = beta2.pois)
        }
        
      }else{
        if(k %in% IndexBino){
          
          glmfit <- glm(X.all[, k] ~ X.all[, C[[k]]], family = "binomial", control=list(maxit=100))
          glmfit.coef <- as.matrix(glmfit$coefficients)
          
          if(any(is.na(glmfit.coef))){
            na.index <- which(is.na(glmfit.coef) == TRUE)
            glmfit.coef[na.index] <- 0
          }
          
          glm_St <- cbind(X.all[, k], X.all[, C[[k]]])
          exp.term <- exp(glm_St %*% glmfit.coef)
          con.prob <- mean(exp.term^(X.all[, k]) * (factorial(1)/(factorial(X.all[, k])*factorial(1-X.all[, k])))
                           / (1 + exp.term)^(X.all[, k]))
          con.mean <-  con.prob * sum(X[, k])
          sec.moment <- con.prob * sum(X[, k]^2)
          con.var <- sec.moment - con.mean^2
          
          w <- (beta1.bino + beta2.bino * con.mean)^(-1)
          nx <- n*length(C[[k]])
          nC <- sum( nx * ifelse(nx >= c0*n, 1, 0) )   
          kk <- which(left.order == k)
          S[kk] <- sum( (w^2 * con.var - w * con.mean) * nx / nC ) 
          
        }else{
          glmPois <- glm( X[, k] ~ X[, C[[k]]], family = "poisson", control = list(maxit=100)) 
          glm2Pois <- glm( X[, k]^2 ~ X[, C[[k]]], family = "poisson", control = list(maxit=100))
          con.mean <- mean(exp(predict(glmPois)))
          sec.moment <- mean(exp(predict(glm2Pois)))
          con.var <- sec.moment - con.mean^2
          
          w <- (beta1.pois + beta2.pois * con.mean)^(-1)
          nx <- n*length(C[[k]])
          nC <- sum( nx * ifelse(nx >= c0*n, 1, 0) )   
          S[kk] <- sum( (w^2 * con.var - w * con.mean) * nx / nC ) 
          
        }
      } 
    }
    
    ODS.order[m] <- left.order[which.min(S)]
    #ODS.left[[m]] <- ODS.left[[m-1]][-which.min(S)]
  }
  
  ODS.order[p] <- seq(p)[-ODS.order]
  
  ########################
  ##Step3: directed edges
  ########################
  
  ordering <- ODS.order
  rr <- rev(ordering)
  result <- matrix(0, p, p)
  
  for (ii in 1:(p-1)){
    now <- rr[ii]
    this <- sort(rr[(ii+1):p])
    if (length(this) > 1){
      
      if(now %in% IndexBino){
        nfold <- 5
        fold0 <- sample.int(sum(X.all[, now]==0)) %% nfold
        fold1 <- sample.int(sum(X.all[, now]==1)) %% nfold
        foldid <- numeric(length(X.all[, now]))
        foldid[X.all[, now]==0] <- fold0
        foldid[X.all[, now]==1] <- fold1
        foldid <- foldid + 1
        fraction <- table(X.all[, now])/length(X.all[, now])
        weights <- 1 - fraction[as.character(X.all[, now])]
        
        lassomBino <- glmnet::cv.glmnet(X.all[, left], X.all[, now], family = "binomial", foldid = foldid, weights = weights)
        bfit <- coefficients(lassomBino)[-1]
        
      }else{
        lassomPois <- glmnet::cv.glmnet(X[, left], X[, now], family = "poisson", alpha = 1, nfolds = 5)
        bfit <- coefficients(lassomPois)[-1]
        
      }
      
      for (jj in 1:length(this)){
        if( bfit[jj] != 0 )
          result[this[jj],now] <- 1
      }
    }else {
      # deal with the last two nodes
      if(now %in% IndexBino){
        lmodBino <- glm(X.all[, now] ~ X.all[, this], family = "binomial", control = list(maxit=100))
        coeff <- as.vector(lmodBino$coefficients[2])
        
      }else{
        lmodPois <- glm(X[, now] ~ X[, this], family = "poisson", control = list(maxit=100))
        coeff <- as.vector(lmodPois$coefficients[2])
        
      }
      
      if (coeff != 0) 
        result[this,now] <- 1
    }
  }
  return(list(adj=result,ordering=rev(rr)))
}



## a helper function for ODS algorithm

RatioODS <- function(X, beta_1, beta_2){
  con.mean1 <- apply(X, 2, mean)
  con.var1  <- apply(X, 2, var)
  w1 <- (beta_1 + beta_2 * con.mean1)^(-1)
  S1 <- w1^2 * con.var1 - w1 * con.mean1
  
  return(S1)
}












#############################
########## Experiment
#############################

## N: Binomial parameter
do_one <- function(do.n, do.p, do.N, do.sparsity=c("dense","sparse"), do.grid.t, do.theta1.pois, do.theta2.pois, do.theta1.pois.beta, do.theta2.pois.beta, do.theta1.bino, do.theta2.bino, do.seed){
  do.result <- matrix(0, 6, 5) ##
  # Data generation
  #dat <- hubgraph(n, p)
  do.dat <- hubMixedgraph(do.n, do.p, do.N, do.sparsity, do.theta1.pois, do.theta2.pois, do.theta1.pois.beta, do.theta2.pois.beta, do.theta1.bino, do.theta2.bino, do.seed)
  do.truth <- do.dat$truth
  do.X <- do.dat$X
  
  # EDAG
  ptm <- proc.time()
  do.sresult <- EDAG(do.X, do.grid.t, 1, 0, 1, -1/do.N, do.N)
  do.sx <- do.sresult$adj
  do.result[1, 1] <- (proc.time() - ptm)[1]
  do.result[2, 1] <- hammingDistance(do.sx, do.truth)/(do.p*(do.p-1))
  do.result[3, 1] <- ifelse(sum(do.truth), sum(do.truth*do.sx)/sum(do.truth), 0) #recall
  do.result[4, 1] <- ifelse(sum(do.truth), sum(do.truth*do.sx)/(sum(do.sx)+1e-6), 0) # precision
  do.result[5, 1] <- 2*do.result[3, 1]*do.result[4,1]/(do.result[3, 1]+do.result[4,1]+1e-6)  # F1 score
  do.result[6, 1] <- ifelse(sum(do.sx), 1 - sum(do.truth*do.sx)/(sum(do.sx)+1e-6), 0) #FDR
  #sxtd <- sx
  
  
  # ODSGLMLasso
  ptm <- proc.time()
  sresult <- ODSGLMLasso(do.X, 1, 0, 1, -1/do.N, do.N)
  sx <- sresult$adj 
  do.result[1, 2] <- (proc.time() - ptm)[1]
  do.result[2, 2] <- hammingDistance(sx, do.truth)/(do.p*(do.p-1))     ##normalized hammdis
  do.result[3, 2] <- ifelse(sum(do.truth), sum(do.truth*sx)/sum(do.truth), 0) #recall
  do.result[4, 2] <- ifelse(sum(do.truth), sum(do.truth*sx)/(sum(sx)+1e-6), 0) # precision
  do.result[5, 2] <- 2*do.result[3, 2]*do.result[4,2]/(do.result[3, 2]+do.result[4,2]+1e-6)  # F1-score 
  do.result[6, 2] <- ifelse(sum(sx), 1 - sum(do.truth*sx)/(sum(sx)+1e-6), 0) #FDR
  
  
  # MRS
  ptm <- proc.time()
  Mresult <- MRS(do.X)
  Mx <- Mresult$adj 
  do.result[1, 3] <- (proc.time() - ptm)[1]
  do.result[2, 3] <- hammingDistance(Mx, do.truth)/(do.p*(do.p-1))     ##normalized hammdis
  do.result[3, 3] <- ifelse(sum(do.truth), sum(do.truth*Mx)/sum(do.truth), 0) #recall
  do.result[4, 3] <- ifelse(sum(do.truth), sum(do.truth*Mx)/(sum(Mx)+1e-6), 0) # precision
  do.result[5, 3] <- 2*do.result[3, 3]*do.result[4,3]/(do.result[3, 3]+do.result[4,3]+1e-6)  # F1-score 
  do.result[6, 3] <- ifelse(sum(Mx), 1 - sum(do.truth*Mx)/(sum(Mx)+1e-6), 0) #FDR
  
  
  # DirectLINGAM
  ptm <- proc.time()
  lx <- DirectLINGAM(do.X)
  do.result[1, 4] <- (proc.time() - ptm)[1]
  do.result[2, 4] <- hammingDistance(lx, do.truth)/(do.p*(do.p-1))     ##normalized hammdis
  do.result[3, 4] <- ifelse(sum(do.truth), sum(do.truth*lx)/sum(do.truth), 0) #recall
  do.result[4, 4] <- ifelse(sum(do.truth), sum(do.truth*lx)/(sum(lx)+1e-6), 0) # precision
  do.result[5, 4] <- 2*do.result[3, 3]*do.result[4,3]/(do.result[3, 3]+do.result[4,3]+1e-6)  # F1-score 
  do.result[6, 4] <- ifelse(sum(lx), 1 - sum(do.truth*lx)/(sum(lx)+1e-6), 0) #FDR
  
  
  # GES
  ptm <- proc.time()
  gx.log <- gesWrap(do.X)$Adj
  gx <- matrix(as.numeric(gx.log), do.p, do.p)
  do.result[1, 5] <- (proc.time() - ptm)[1]
  do.result[2, 5] <- hammingDistance(gx, do.truth)/(do.p*(do.p-1))     ##normalized hammdis
  do.result[3, 5] <- ifelse(sum(do.truth), sum(do.truth*gx)/sum(do.truth), 0) #recall
  do.result[4, 5] <- ifelse(sum(do.truth), sum(do.truth*gx)/(sum(gx)+1e-6), 0) # precision
  do.result[5, 5] <- 2*do.result[3, 4]*do.result[4,4]/(do.result[3, 4]+do.result[4,4]+1e-6)  # F1-score 
  do.result[6, 5] <- ifelse(sum(gx), 1 - sum(do.truth*gx)/(sum(gx)+1e-6), 0) #FDR
  
  
  # MMHC
  ptm <- proc.time()
  mmhc.cha <- mmhc(data.frame(do.X))$arcs
  mx <- charaToNume(mmhc.cha, do.p)
  do.result[1, 6] <- (proc.time() - ptm)[1]
  do.result[2, 6] <- hammingDistance(mx, do.truth)/(do.p*(do.p-1))     ##normalized hammdis
  do.result[3, 6] <- ifelse(sum(do.truth), sum(do.truth*mx)/sum(do.truth), 0) #recall
  do.result[4, 6] <- ifelse(sum(do.truth), sum(do.truth*mx)/(sum(mx)+1e-6), 0) # precision
  do.result[5, 6] <- 2*do.result[3, 5]*do.result[4, 5]/(do.result[3, 5]+do.result[4, 5]+1e-6)  # F1-score 
  do.result[6, 6] <- ifelse(sum(mx), 1 - sum(do.truth*mx)/(sum(mx)+1e-6), 0) #FDR
  
  
  return(list(DAG=do.truth, EDAGs=do.sresult$adj, EDAGl=do.sresult$OrderLayer, ODSs=sresult$adj, ODSl=sresult$ordering, Perf=do.result))
}



simu_replicate <- function(simr.n, simr.p, simr.N, simr.sparsity=c("dense","sparse"), simr.times, simr.grid.t, simr.theta1.pois, simr.theta2.pois, simr.theta1.pois.beta, simr.theta2.pois.beta, simr.theta1.bino, simr.theta2.bino){
  simr.results <- list()
  simr.sd <- list()  ## standard deviation
  
  for (i in 1:simr.times) {
    print(i)
    #results[[i]] <- do_one(n = n, p = p, grid.t)
    simr.results[[i]] <- do_one(simr.n, simr.p, simr.N, simr.sparsity, simr.grid.t, simr.theta1.pois, simr.theta2.pois, simr.theta1.pois.beta, simr.theta2.pois.beta, simr.theta1.bino, simr.theta2.bino, i)$Perf
    
  } 
  simr.result  <- Reduce('+', simr.results)
  Metric = simr.result/simr.times   ## mean of all the metrics 
  simr.sd <- apply(array(unlist(simr.results), c(6, 5, simr.times)), c(1,2), sd)  ## standard deviation using 3Darray
  
  return(list(Metrics = round(Metric, 2), SD = round(simr.sd, 2)))
  
}





#####################################
## generate data for hub graph (d=1)  
#####################################




hubMixedgraph <- function(n, p, N, sparsity=c("dense","sparse"), theta1.pois, theta2.pois, theta1.pois.beta, theta2.pois.beta, theta1.bino, theta2.bino, hub.seed){
  ordering <- seq(p)
  truth <- matrix(0, p, p)
  hub.X <- matrix(0, n, p)
  
  
  ## generate the parameters
  set.seed(hub.seed)
  theta.pois <- runif(p, theta1.pois, theta2.pois)
  set.seed(hub.seed)
  beta.pois <- matrix(runif(p*p, theta1.pois.beta, theta2.pois.beta), p ,p)
  
  set.seed(hub.seed)
  theta.bino <- runif(p, theta1.bino, theta2.bino)
  set.seed(hub.seed)
  beta.bino  <- matrix(runif(p*p, theta1.bino, theta2.bino), p ,p)
  
  ## generate the root node as poisson 
  set.seed(hub.seed)
  hub.index <- sample(2, 1, replace = FALSE)    
  
  # randomly root node as the only parent node
  if(hub.index == 1){
    #poisson root node
    A0.lambda <- exp(theta.pois[1])
    set.seed(hub.seed)
    hub.X[, 1] <- rpois(n, lambda = A0.lambda)
    
  }else{
    #binomial root node
    hub.A0.prob <- exp(theta.bino[1]) / (1 + exp(theta.bino[1]))
    set.seed(hub.seed)
    hub.X[, 1]  <- rbinom(n, size = N, prob = hub.A0.prob)
    
  }
  
  
  
  if(sparsity == "dense"){    # for Example 2 the hub graph for mixed DAG
    truth[1, 2:p] <- 1
    child <- c(2:p)
    
    
    #generate the sink nodes
    set.seed(hub.seed)
    hub.sink.index <- sample(2, p - 1, replace = TRUE)
    for (j in 2:p) {
      if(hub.sink.index[j-1] == 1){
        #poisson sink nodes
        A1.lambda <- exp(theta.pois[j] + hub.X[, 1] * beta.pois[j, 1])
        set.seed(hub.seed)
        hub.X[, j] <- rpois(n, lambda = A1.lambda)
        
      }else{
        #binomial sink nodes
        g.prob <- exp(theta.bino[j] + beta.bino[j, 1] * hub.X[, 1]) / (1+exp(theta.bino[j] + beta.bino[j, 1] * hub.X[, 1]))
        set.seed(hub.seed)
        hub.X[, j] <- rbinom(n, size = N, prob = g.prob)
        
      }
    }
    
    
    
  }else{
    set.seed(hub.seed)
    direct.edge <- sample(p-1, ceiling(p/2)-1, replace = FALSE) + 1  ## randomly choose p/2 edges for hub graph 
    truth[1, direct.edge] <- 1  
    child <- direct.edge
    root <- setdiff(c(2:p), intersect(c(2:p), direct.edge))   # other half root nodes except X1
    direct_length <- length(direct.edge)
    root_length <- length(root)
    
    #generate the half root nodes
    set.seed(hub.seed)
    hub.root.index <- sample(2, root_length, replace = TRUE)
    for (k in 1:root_length) {
      root.index <- root[k]
      if(hub.root.index[k] == 1){
        #poisson root node
        Ak.lambda <- exp(theta.pois[root.index])
        set.seed(hub.seed)
        hub.X[, root.index] <- rpois(n, lambda = Ak.lambda)
        
      }else{
        #binomial root node
        hub.Ak.prob <- exp(theta.bino[root.index]) / (1 + exp(theta.bino[root.index]))
        set.seed(hub.seed)
        hub.X[, root.index]  <- rbinom(n, size = N, prob = hub.Ak.prob)
        
      }
    }
    
    #genrate the sink nodes
    set.seed(hub.seed)
    hub.sink.index <- sample(2, direct_length, replace = TRUE)
    for (j in 1:direct_length) {
      direct.index <- direct.edge[j]
      if(hub.sink.index[j] == 1){
        #poisson sink nodes
        A1.lambda <- exp(theta.pois[direct.index] + hub.X[, 1] * beta.pois[direct.index, 1])
        set.seed(hub.seed)
        hub.X[, direct.index] <- rpois(n, lambda = A1.lambda)
        
      }else{
        #binomial sink nodes
        g.prob <- exp(theta.bino[direct.index] + beta.bino[direct.index, 1] * hub.X[, 1]) / (1+exp(theta.bino[direct.index] + beta.bino[direct.index, 1] * hub.X[, 1]))
        set.seed(hub.seed)
        hub.X[, direct.index] <- rbinom(n, size = N, prob = g.prob)
        
      }
    }
    
  }
  
  #return(list(X = hub.X, truth = truth, TO = ordering, Index = c(hub.index, hub.sink.index)))
  return(list(X = hub.X, truth = truth, TO = ordering, Index = hub.sink.index))
}








##### Examples
grid.t <- 10^(-2 + 0.15*seq(0, 60, by=1))

a1<-simu_replicate(simr.n=200, simr.p=5, simr.N=4, simr.sparsity="sparse", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                   simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)
a3<-simu_replicate(simr.n=500, simr.p=5, simr.N=4, simr.sparsity="sparse", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                   simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)


b1 <- simu_replicate(simr.n=200, simr.p=20, simr.N=4, simr.sparsity="sparse", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                     simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)
b3 <- simu_replicate(simr.n=500, simr.p=20, simr.N=4, simr.sparsity="sparse", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                     simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)


c1 <- simu_replicate(simr.n=200, simr.p=100, simr.N=4, simr.sparsity="sparse", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                     simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)
c3 <- simu_replicate(simr.n=500, simr.p=100, simr.N=4, simr.sparsity="sparse", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                     simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)



