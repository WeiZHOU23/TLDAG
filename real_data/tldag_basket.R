
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
Kappa.cal <- function(KC.layer, KC.X, KC.n, KC.p, KC.beta_1, KC.beta_2, KC.index){
  KC.n1 <- as.integer(KC.n/2)
  KC.n2 <- KC.n - KC.n1
  KC.done <- unlist(KC.layer) ## The set S
  set.seed(KC.index) #### Add seed here
  KC.index <- sample(1:KC.n, KC.n, replace = FALSE)
  data.perm <- KC.X[KC.index, ] 
  KC.data1 <- data.perm[1:KC.n1, ]
  KC.data2 <- data.perm[-(1:KC.n1), ]
  if(length(KC.done) == 0){ ## A_0 for roots
    KC.set1 <- layer_A0(KC.data1, KC.beta_1, KC.beta_2 )
    KC.set2 <- layer_A0(KC.data2, KC.beta_1, KC.beta_2 )
  }else if( length(KC.done) >=1 & length(KC.done) < (KC.p-1)){  ## layer A_1 to A_{T-1}
    KC.set1 <- layer_At(KC.layer, KC.data1, KC.beta_1, KC.beta_2)
    KC.set2 <- layer_At(KC.layer, KC.data2, KC.beta_1, KC.beta_2)
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
Kappa <- function(Ka.layer, Ka.X, Ka.grid.t, Ka.beta_1, Ka.beta_2){
  Ka.n <- dim(Ka.X)[1]
  Ka.p <- dim(Ka.X)[2]
  Ka.done<-unlist(Ka.layer)
  # print(Ka.done)
  if(is.null(Ka.done)==T){
    Ka.left<-seq(Ka.p)
  }else{
    Ka.left<-seq(Ka.p)[-Ka.done]
  }
  Kappa.temp<-sapply(c(1,2,3,4,5), Kappa.cal, KC.layer=Ka.layer, KC.X=Ka.X, KC.n=Ka.n, KC.p=Ka.p, KC.beta_1=Ka.beta_1, KC.beta_2=Ka.beta_2)
  #kappa.temp1<-matrix(c(Kappa.temp[[1]],Kappa.temp[[3]],Kappa.temp[[5]],Kappa.temp[[7]],Kappa.temp[[9]]),5, length(Kappa.temp[[1]]),byrow = F)
  #kappa.temp2<-matrix(c(Kappa.temp[[2]],Kappa.temp[[4]],Kappa.temp[[6]],Kappa.temp[[8]],Kappa.temp[[10]]), 5, length(Kappa.temp[[2]]),byrow = F)
  kappa.temp1<-matrix(c(Kappa.temp[[1]],Kappa.temp[[3]],Kappa.temp[[5]],Kappa.temp[[7]],Kappa.temp[[9]]),5, length(Kappa.temp[[1]]),byrow = T)
  kappa.temp2<-matrix(c(Kappa.temp[[2]],Kappa.temp[[4]],Kappa.temp[[6]],Kappa.temp[[8]],Kappa.temp[[10]]), 5, length(Kappa.temp[[2]]),byrow = T)
  Ka.result<-sapply(Ka.grid.t, Kapp.compute, kap.re1=kappa.temp1, kap.re2=kappa.temp2, kap.left=Ka.left)
  return(Ka.result/5)
}


#################################################
TuningPara <- function(Tune.layer, Tune.X, Tune.grid.t, Tune.alpha, Tune.beta_1, Tune.beta_2){
  kappa.vec <- Kappa(Tune.layer, Tune.X, Tune.grid.t, Tune.beta_1, Tune.beta_2)
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
#########################################################



####################################################
## roots for the layer A_0  ##
layer_A0 <- function(A0.X, A0.beta_1, A0.beta_2){
  A0.n <- dim(A0.X)[1]
  A0.p <- dim(A0.X)[2]
  A0.TF <- matrix(0, A0.n, A0.p)
  A0.TF<-t( (A0.beta_1 + A0.beta_2 *colMeans(A0.X))^(-1) ) %x% rep(1,A0.n)*A0.X
  A0.ratio <- c()
  A0.ratio <-  apply(A0.TF,2,var) /colMeans(A0.TF)
  A0.Adjustratio <- A0.ratio / min(A0.ratio)
  return(A0.Adjustratio)
}

########################################################################################


layer_At <-function(At.layer, At.X, At.beta_1, At.beta_2){
  At.n1 <- dim(At.X)[1]
  At.p <- dim(At.X)[2]
  At.Sdone <- unlist(At.layer)
  At.left <- seq(At.p)[-At.Sdone]
  At.n.left <- length(At.left)
  At.con.mean = At.con.var = c();
  for (k in 1:At.n.left) {
    glmfit <- glm(At.X[, At.left][, k] ~ At.X[, At.Sdone], family = "poisson", control=list(maxit=100))
    con.mean1Pois <- mean(exp(predict(glmfit)))
    Tt2Pois <- (At.X[, At.left][, k])^2
    At.con.mean[k] <- mean(At.X[, At.left][, k])
    
    
    fitT2Pois <- glm(Tt2Pois ~ At.X[, At.Sdone], family = "poisson", control=list(maxit=100)) 
    sec.moment <- mean(exp(predict(fitT2Pois)))
    At.con.var[k] <- sec.moment - con.mean1Pois^2
    
  }
  At.ratio <- At.con.var / At.con.mean
  final.result <- At.ratio / min(At.ratio)
  return(final.result) 
}



### helper functions
# lasso regression to find the parents of the node 
potentialParents <- function(x, index){
  parents <- NULL
  T <- length(x)
  for (kk in (index+1):T) {
    NewParents <- x[[kk]]
    parents <- c(parents, NewParents)
  }
  return(parents)
} 


#######################################################

EDAG <- function(ED.X, ED.grid.t, ED.beta_1, ED.beta_2, tuning = c("small","large","adaptive"), c0, lambda.index){
  ptm <- proc.time() 
  ED.n <- dim(ED.X)[1]
  ED.p <- dim(ED.X)[2]
  ED.TF <- matrix(0, ED.n, ED.p)
  ED.done <- NULL
  ED.layer <- list()
  ED.ratio.t <- list()
  ED.thres <- list()
  ED.final.layer<-list()
  ##################################################################################  
  # Step1: layer and ordering
  # roots layer A_0
  ############## Node identification ##################
  EDA0.thresh <- TuningPara(ED.layer, ED.X, ED.grid.t,  0.1, ED.beta_1, ED.beta_2)$thres
  ED.ratio.final<-layer_A0(ED.X, ED.beta_1, ED.beta_2)
  ED.layer[[1]] <- which(abs(ED.ratio.final - 1) <= EDA0.thresh)
  ED.thres[[1]] <- EDA0.thresh
  ED.ratio.t[[1]]<-ED.ratio.final
  #######################################################################################
  # identify nodes in layer A_t
  cont <- TRUE
  while( cont ){
    for (lay in 2:ED.p) {
      if( length(unlist(ED.layer)) < ED.p-1){   
        
        ED.left <- seq(ED.p)[-unlist(ED.layer)]
        print(lay)
        ####################################################################
        EDAT.thresh <- TuningPara(ED.layer, ED.X, ED.grid.t, 0.1, ED.beta_1, ED.beta_2)$thres
        ED.thres[[lay]] <- EDAT.thresh
        ED.ATratio.final<-layer_At(ED.layer, ED.X, ED.beta_1, ED.beta_2)
        ED.layer[[lay]] <- ED.left[which(abs(ED.ATratio.final - 1) <= EDAT.thresh)]
        ED.ratio.t[[lay]]<-ED.ATratio.final
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
  
  ED.final.layer <- list(ED.layer, ED.ratio.t, ED.thres)  ##list 其中layer[[1]]表示topological layer
  
  
  ##Step2:Sparse Estimation for Direction Edges
  EDS.layernodes <- rev(ED.final.layer[[1]]) ### Selected layer structure with nodes in each layer
  EDS.layernum  <- length(ED.final.layer[[1]]) ### number of layer in the DAG
  EDS.result <- matrix(0, ED.p, ED.p)
  
  
  
  for (ii in 1:(EDS.layernum-1)) {   ###
    for (k in 1:length(EDS.layernodes[[ii]])) {
      EDS.now <- EDS.layernodes[[ii]][k]  ## each element in each layer ### the node as a child
      ## now turn to reconstruct the directed structure among SE.now and all the nodes in the upper layers
      
      #for (j in (ii+1):EDS.layernum) {
      
      #EDS.this <- EDS.layernodes[[j]]  ## each layer ### should be the nodes in the upper layer
      EDS.this <- potentialParents(EDS.layernodes, ii)
      
      if(length(EDS.this) > 1){
        
        if(tuning == "small"){
          EDS.lassom <- glmnet::cv.glmnet(ED.X[, EDS.this], ED.X[, EDS.now], family = "poisson", nfolds = 5)
          EDS.bfit <- coefficients(EDS.lassom)[-1]
        }
        if(tuning == "large"){
          large.lambda <- largeTuningPara(ED.X[, EDS.this], ED.X[, EDS.now], c0)[2]
          EDS.lassom <- glmnet::glmnet(ED.X[, EDS.this], ED.X[, EDS.now], lambda = large.lambda)
          EDS.bfit <- coefficients(EDS.lassom)[-1]
        }
        if(tuning == "adaptive"){
          EDS.lassoma <- glmnet::glmnet(ED.X[, EDS.this], ED.X[, EDS.now], lambda = grid.tc[lambda.index])
          EDS.bfit <- coefficients(EDS.lassoma)[-1]
          
        }
        
        
        
        for (mm in 1:length(EDS.this)) {
          if(EDS.bfit[mm] != 0)
            EDS.result[EDS.this[mm], EDS.now] <- 1
        }
      }else{
        EDS.lmod <- glm(ED.X[, EDS.now] ~ ED.X[, EDS.this], family = "poisson", control=list(maxit=100))
        EDS.bfit <- as.vector(EDS.lmod$coefficients[2])
        
        if(EDS.bfit != 0)
          EDS.result[EDS.this, EDS.now] <- 1
      }
    }
  }   
  #}
  #elapse <- proc.time() - ptm
  return(list(adj = EDS.result, OrderLayer = rev(EDS.layernodes)))
}
#########################################################################################




largeTuningPara <- function(x, y, c = 2.5){
  
  cv.lasso <- glmnet::cv.glmnet(x, y, family = "poisson", nfolds = 5)
  lambda_range <- cv.lasso$lambda
  res <- cbind(lambda_range, NaN)
  for (i in 1:length(lambda_range)) {
    res[i, 2] <- mean(predict(cv.lasso, newx = x, s = lambda_range[i]))
  }
  bestlam <- cv.lasso$lambda.min
  mini.mse <- mean(predict(cv.lasso, newx = x, s = bestlam))
  large.mse <- c * mini.mse
  large.mse.index <- which.min(abs(res[, 2] - large.mse))
  #large.tunpara <- as.numeric(res[1,][large.mse.index])
  large.tunpara <- as.numeric(res[large.mse.index, ][1])
  
  return(c(bestlam,large.tunpara))
} 


