########################################################
##### Main Functions for Random Graph for Binomial DAG
########################################################
#install.packages("pscl")
#install.packages("doParallel")
#install.packages("foreach")
#install.packages("doSNOW")
#install.packages("glmnet")
#install.packages("RcppEigen")
#install.packages("iterators")
#install.packages("matrixcalc")
library(pscl)
library(glmnet)
library(doParallel)
library(matrixcalc)
library(foreach)
library(doSNOW)
library(iterators)
library(RcppEigen)



##############################
### main functions
#############################


#################################################
### Tuning Parameter by Sun et al (2013 JMLR) ###
#################################################
## Helper Function
## Calculate the Kappa coefficient of Two Sets

##########################
### 并行方法计算 kappa  ##
##########################

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
Kappa.cal <-function(KC.layer, KC.X, KC.n, KC.p, KC.beta_1, KC.beta_2, KC.index){
  KC.n1 <- as.integer(KC.n/2)
  KC.n2 <- KC.n - KC.n1
  KC.done <- unlist(KC.layer) ## The set S
  set.seed(KC.index) #### Add seed here
  KC.index<-sample(1:KC.n, KC.n, replace = FALSE)
  data.perm <- KC.X[KC.index, ] 
  KC.data1 <- data.perm[1:KC.n1, ]
  KC.data2 <- data.perm[-(1:KC.n1), ]
   if(length(KC.done) == 0){ ## A_0 for roots
    KC.set1 <- layer_A0(KC.data1, KC.beta_1, KC.beta_2 )
    KC.set2 <- layer_A0(KC.data2, KC.beta_1, KC.beta_2 )
  }else if( length(KC.done) >=1 & length(KC.done) < (KC.p-1)){  ## layer A_1 to A_{T-1}
    KC.set1 <- layer_At(KC.layer, KC.data1, KC.beta_1, KC.beta_2 )
    KC.set2 <- layer_At(KC.layer, KC.data2, KC.beta_1, KC.beta_2 )
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
    kappa.temp1<-matrix(c(Kappa.temp[[1]],Kappa.temp[[3]],Kappa.temp[[5]],Kappa.temp[[7]],Kappa.temp[[9]]),5, length(Kappa.temp[[1]]),byrow = F)
    kappa.temp2<-matrix(c(Kappa.temp[[2]],Kappa.temp[[4]],Kappa.temp[[6]],Kappa.temp[[8]],Kappa.temp[[10]]), 5, length(Kappa.temp[[2]]),byrow = F)
    Ka.result<-sapply(Ka.grid.t, Kapp.compute, kap.re1=kappa.temp1, kap.re2=kappa.temp2, kap.left=Ka.left)
   return(Ka.result/5)
}
#################################################

TuningPara <- function(Tune.layer, Tune.X, Tune.grid.t, Tune.alpha, Tune.beta_1, Tune.beta_2){
 kappa.vec <- Kappa(Tune.layer, Tune.X, Tune.grid.t, Tune.beta_1, Tune.beta_2)
  if(is.nan(kappa.vec) || is.infinite(kappa.vec)){
    kappa.vec <- kappa.vec[-which(is.nan(kappa.vec))]
    kappa.vec <- kappa.vec[-which(is.infinite(kappa.vec))]
  }
   if(max(kappa.vec) == 0){
    # the maximum of kappa.vec is 0
    i.opt.kappa <- grid.t[which( kappa.vec >=  max(kappa.vec)*(1-Tune.alpha) )[1]]
  }else{
    i.opt.kappa <- grid.t[which( kappa.vec / max(kappa.vec) >= (1-Tune.alpha) )[1]]
  }
  result <- list(alpha = Tune.alpha, kappa.values = kappa.vec, thres = i.opt.kappa)
  return(result)
}
#########################################################


####################################################
## roots for the layer A_0###
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
########################################################################################



layer_At <-function(At.layer, At.X, At.beta_1, At.beta_2){
  At.n1 <- dim(At.X)[1]
  At.p <- dim(At.X)[2]
  At.Sdone <- unlist(At.layer)
  At.left <- seq(At.p)[-At.Sdone]
  At.n.left <- length(At.left)
  At.con.mean = At.con.var = c();
  for (k in 1:At.n.left) {
    glmfit <- glm(At.X[, At.left][, k] ~ At.X[, At.Sdone], family = "binomial", control=list(maxit=100))
    glm.probs <- predict(glmfit, type = "response")
    
    con.mean1 <- t(glm.probs) %*% At.X[, At.left][, k] 
    At.Tt <- ((At.beta_1 + At.beta_2 * con.mean1)^(-1))* At.X[, At.left][,k] 
    At.con.mean[k] <- glm.probs %*% At.Tt 
    At.con.var[k] <-  glm.probs %*% At.Tt^2 - At.con.mean[k]^2
  }
  At.ratio <- At.con.var / At.con.mean
  final.result <- At.ratio / min(At.ratio)
  return(final.result) 
}



EDAG <- function(ED.X, ED.grid.t, ED.beta_1, ED.beta_2){
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
 # while( cont ){
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
        ED.n.left <- length(left)
        ED.layer[[length(layer) + 1]] <- ED.n.left
        
      }else (length(unlist(ED.layer)) >= ED.p)
      #cont <- FALSE
      break
   }
  #}
  
  ED.final.layer <- list(ED.layer, ED.ratio.t, ED.thres)  ##list 其中layer[[1]]表示topological layer

  
  ##Step2:Sparse Estimation for Direction Edges
  EDS.layernodes <- rev(ED.final.layer[[1]]) ### Selected layer structure with nodes in each layer
  EDS.layernum  <- length(ED.final.layer[[1]]) ### number of layer in the DAG
  EDS.result <- matrix(0, ED.p, ED.p)

  
  
  for (ii in 1:(EDS.layernum-1)) {   ###
    for (k in 1:length(EDS.layernodes[[ii]])) {
      EDS.now <- EDS.layernodes[[ii]][k]  ## each element in each layer ### the node as a child
      ## now turn to reconstruct the directed structure among SE.now and all the nodes in the upper layers
      
      for (j in (ii+1):EDS.layernum) {
        
        EDS.this <- EDS.layernodes[[j]]  ## each layer ### should be the nodes in the upper layer
        if(length(EDS.this) > 1){
          EDS.nfold <- 5
          EDS.fold0 <- sample.int(sum(ED.X[, EDS.now]==0)) %% EDS.nfold
          EDS.fold1 <- sample.int(sum(ED.X[, EDS.now]==1)) %% EDS.nfold
          EDS.foldid <- numeric(length(ED.X[, EDS.now]))
          EDS.foldid[ED.X[, EDS.now]==0] <- EDS.fold0
          EDS.foldid[ED.X[, EDS.now]==1] <- EDS.fold1
          EDS.foldid <- EDS.foldid + 1
          EDS.fraction <- table(ED.X[, EDS.now])/length(ED.X[, EDS.now])
          EDS.weights <- 1 - EDS.fraction[as.character(ED.X[, EDS.now])]
          
          EDS.lassom <- glmnet::cv.glmnet(ED.X[, EDS.this], ED.X[, EDS.now], family = "binomial", foldid = EDS.foldid, weights = EDS.weights )
          EDS.bfit <- coefficients(EDS.lassom)[-1]
          for (mm in 1:length(EDS.this)) {
            if(EDS.bfit[mm] != 0)
              EDS.result[EDS.this[mm], EDS.now] <- 1
          }
        }else{
          EDS.lmod <- glm(ED.X[, EDS.now] ~ ED.X[, EDS.this], family = "binomial", control=list(maxit=100))
          EDS.bfit <- as.vector(EDS.lmod$coefficients[2])
          
          if(EDS.bfit != 0)
            EDS.result[EDS.this, EDS.now] <- 1
        }
      }
      
      
      
      
    }   
  }
  #elapse <- proc.time() - ptm
  return(list(adj = EDS.result, OrderLayer = rev(EDS.layernodes)))
}
#########################################################################################



######################################################################################3
### ODS(GLMLasso)
ODSGLMLasso <- function(X, beta_1, beta_2){
  ptm<-proc.time()
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Step1: NormalizedGraph
  NormalizedGraph <- matrix(0, p, p)
  for (i in 1:p) {
    left <- seq(p)[-i]
    
    nfold <- 5
    fold0 <- sample.int(sum(X[, i]==0)) %% nfold
    fold1 <- sample.int(sum(X[, i]==1)) %% nfold
    foldid <- numeric(length(X[, i]))
    foldid[X[, i]==0] <- fold0
    foldid[X[, i]==1] <- fold1
    foldid <- foldid + 1
    fraction <- table(X[, i])/length(X[, i])
    weights <- 1 - fraction[as.character(X[, i])]
    
    lassom <- glmnet::cv.glmnet(X[, left], X[, i], family = "binomial", foldid = foldid, weights = weights)
    bfit <- coefficients(lassom)[-1]
    for (j in 1:length(bfit)) {
      if(bfit[j] != 0)
        NormalizedGraph[i, left[j]] <- 1
    }
  }
  ## return(NormalizedGraph)
  
  #Step2: ordering
  ODS.order <- c()
  
  ## first element of the ordering
  con.mean1 <- apply(X, 2, mean)
  con.var1  <- apply(X, 2, var)
  w1 <- (beta_1 + beta_2* con.mean1)^(-1)
  S1 <- w1^2 * con.var1 - w1 * con.mean1
  ODS.order[1] <- seq(p)[which.min(S1)]
  
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
    ODS.order <- ODS.order  ####解决S每次相同的问题, renew the ordering
    left.order <- seq(p)[-ODS.order]
    S <- c()
    
    for (k1 in 1:length(left.order) ) { # k 代表原始的变量
      k <- left.order[k1]
      neigh <- union(which( norgraph[k, ] != 0), which( norgraph[, k] != 0))   #k'neighborhood 是否要考虑第k列不为0，两者的union, YES!
      C[[k]] <- intersect(neigh, ODS.order)   # given parents
      
      # calculate conditinal mean/variance  
      if(length(C[[k]]) ==0 ){  
        con.mean <- mean(X[, k])
        con.var  <- var(X[, k])
        
        w <- (beta_1 + beta_2 * con.mean)^(-1)
        kk <- which(left.order == k)
        S[kk] <- w^2 * con.var - w * con.mean
        
      }else{
        #para.mean <- glm( X[, k] ~ X[, C[[k]]], family = "binomial")$coefficients 
        #para.var  <- glm( X[, k]^2 ~ X[, C[[k]]], family = "binomial")$coefficients
        glm.mean <- glm( X[, k] ~ X[, C[[k]]], family = "binomial", control=list(maxit=100))
        glm.probm <- predict(glm.mean, type = "response")
        con.mean <- X[, k] %*% glm.probm
        
        glm.var  <- glm( X[, k]^2 ~ X[, C[[k]]], family = "binomial", control=list(maxit=100))
        glm.probv <- predict(glm.var, type = "response")
        con.mean2 <- X[, k]^2 %*% glm.probv
        con.var <- con.mean2 - con.mean^2
        
        w <- (beta_1 + beta_2 * con.mean)^(-1)
        nx <- n*length(C[[k]])
        nC <- sum( nx * ifelse(nx >= c0*n, 1, 0) )   
        kk <- which(left.order == k)
        S[kk] <- sum( (w^2 * con.var - w * con.mean) * nx / nC ) 
        
      } 
    }
    
    ODS.order[m] <- left.order[which.min(S)]
    
  }
  
  ODS.order[p] <- seq(p)[-ODS.order]
  ## return(ODS.order)
  
  ##Step3: directed edges
  ordering <- ODS.order
  rr <- rev(ordering)
  result <- matrix(0, p, p)
  
  for (ii in 1:(p-1)){
    now <- rr[ii]
    this <- sort(rr[(ii+1):p])
    if (length(this) > 1){
      nfold <- 5
      fold0 <- sample.int(sum(X[, now]==0)) %% nfold
      fold1 <- sample.int(sum(X[, now]==1)) %% nfold
      foldid <- numeric(length(X[, now]))
      foldid[X[, now]==0] <- fold0
      foldid[X[, now]==1] <- fold1
      foldid <- foldid + 1
      fraction <- table(X[, now])/length(X[, now])
      weights <- 1 - fraction[as.character(X[, now])]
      
      lassom <- glmnet::cv.glmnet(X[, this], X[, now], family = "binomial",foldid = foldid, weights = weights)
      bfit <- coefficients(lassom)[-1]
      
      for (jj in 1:length(this)){
        if( bfit[jj] != 0 )
          result[this[jj],now] <- 1
      }
    }else {
      # deal with the last two nodes
      lmod <- glm(X[, now] ~ X[, this], family = "binomial", control = list(maxit=100))
      coeff <- as.vector(lmod$coefficients[2])
      if (coeff != 0) 
        result[this,now] <- 1
    }
  }
  return(list(adj=result,ordering=rev(rr)))
}

####################################


##############################################
### some compariasion functions for later use
##############################################
### hamming distance from EqVarDAG
# hammingDistance(G1,G2)
#
# Computes Hamming Distance between DAGs G1 and G2 with SHD(->,<-) = 1!!!!
#
# INPUT:  G1, G2     adjacency graph containing only zeros and ones: (i,j)=1 means edge from X_i to X_j.
#         
# OUTPUT: hammingDis Hamming Distance between G1 and G2
#
# Copyright (c) 2012-2013  Jonas Peters [peters@stat.math.ethz.ch]

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

##
#library(pcalg)
#shd(g1, g2)



#############################
########## Experiment
#############################

## N: Binomial parameter
do_one <- function(do.n, do.p, do.grid.t, do.at, do.T, do.d, do.fixdegree, do.theta1, do.theta2, do.N,do.seed){
  do.result <- matrix(0, 6, 2) ##
  # Data generation
  #dat <- hubgraph(n, p)
  do.dat <- get_DAGstru(do.n, do.p, do.at, do.T, do.d, do.fixdegree, do.theta1, do.theta2, do.N, do.seed) # N is the size for binomial
  do.truth <- do.dat$DAG
  do.X <- do.dat$X
  
  # EDAG
  ptm <- proc.time()
  do.sresult <- EDAG(do.X, do.grid.t, 1, -1/do.N)
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
  sresult <- ODSGLMLasso(do.X, beta_1 = 1, beta_2 = -1/do.N)
  do.sx <- do.sresult$adj ## ？？？The same？
  do.result[1, 2] <- (proc.time() - ptm)[1]
  do.result[2, 2] <- hammingDistance(do.sx, do.truth)/(do.p*(do.p-1))     ##normalized hammdis
  do.result[3, 2] <- ifelse(sum(do.truth), sum(do.truth*do.sx)/sum(do.truth), 0) #recall
  do.result[4, 2] <- ifelse(sum(do.truth), sum(do.truth*do.sx)/(sum(do.sx)+1e-6), 0) # precision
  do.result[5, 2] <- 2*do.result[3, 2]*do.result[4,2]/(do.result[3, 2]+do.result[4,2]+1e-6)  # F1-score 
  do.result[6, 2] <- ifelse(sum(do.sx), 1 - sum(do.truth*do.sx)/(sum(do.sx)+1e-6), 0) #FDR
  
  return(do.result)
}



####### repeat

#sim_run <- function(SR.n, SR.p, SR.grid.t, SR.simu.times, SR.at, SR.T, SR.d, SR.fixdegree = T, SR.t1, SR.t2, SR.N, SR.pa){
  # p: number of nodes
  # n: sample size
  # N: simulation times
#  if(SR.pa == T){
#    cl <- makeCluster(detectCores()-1)
#    on.exit(stopCluster(cl))
#    registerDoSNOW(cl)
#    print(detectCores())
#    results <- foreach(i = 1:SR.simu.times, .export = ls(globalenv()), .packages = c("glmnet", "matrixcalc","doParallel","doSNOW")) %dopar% {
#      do_one(SR.n, SR.p, SR.grid.t, SR.at, SR.T, SR.d, SR.fixdegree = T, SR.t1, SR.t2, SR.N)
 #   }
#    stopCluster(cl)
#  }else{
#    SR.pb <- txtProgressBar(min = 0, max = SR.simu.times, style = 3)
#    SR.result <- matrix(0, 6, 2)
#    for (i in 1:SR.simu.times) {
#      results <- foreach(i = 1:SR.simu.times, .export = ls(globalenv()),.packages = c("glmnet", "matrixcalc","doParallel","doSNOW")) %do% {
 #       do_one(SR.n, SR.p, SR.grid.t, SR.at, SR.T, SR.d, SR.fixdegree = T, SR.theta1, theta2, SR.N, SR.beta_1, SR.beta_2) ###?? where beta1 and beta2
#      }
#      setTxtProgressBar(SR.pb, i)
#    }
#    close(SR.pb)
#  }
#  SR.result <- Reduce('+', SR.results)
  
#  return(SR.result/SR.simu.times)
#}


#simulation_randomBino <- function(sran.n, sran.p, sran.grid.t, sran.simu.times, sran.at, sran.T, sran.d, sran.fixdegree = T, sran.theta1, sran.theta2, sran.N){
  #set.seed(1)
#  cat('n=',sran.n, "p=",sran.p, "N=",sran.simu.times)
#  sran.result <- sim_run(sran.n, sran.p, sran.grid.t, sran.simu.times, sran.at, sran.T, sran.d, sran.fixdegree = T, sran.theta1, sran.theta2, sran.N, T)
#  sink(paste("p=",sran.p,"n=",sran.n,"_randomBino.txt"))
#  cat('Setting', 'p=',sran.p, "n=",sran.n,"N=",sran.simu.times, "\n")
#  cat('\tEDAG\tODS\n')
#  cat('Avg. time\t', round(sran.result[1,],4),'\n')
#  cat('Avg. Dist\t', round(sran.result[2,],4),'\n')
#  cat('Avg. Reca\t', round(sran.result[3,],4),'\n')
#  cat('Avg. Prec\t', round(sran.result[4,],4),'\n')
#  cat('Avg. fsco\t', round(sran.result[5,],4),'\n')
#  cat('Avg. FDR \t', round(sran.result[6,],4),'\n')
  #cat('SD . Dist\t', round(result[7,],4),'\n')
#  sink()
  
#}

#cat("done.\n")



simu_replicate <- function(simr.n, simr.p, simr.grid.t, simr.at, simr.T, simr.d, simr.fixdegree = T, simr.theta1, simr.theta2, simr.N){
  simr.results <- list()
  for (i in 1:50) {
    print(i)
    #results[[i]] <- do_one(n = n, p = p, grid.t)
    simr.results[[i]] <- do_one(simr.n, simr.p, simr.grid.t, simr.at, simr.T, simr.d, simr.fixdegree, simr.theta1, simr.theta2, simr.N, i)
  }                   
  simr.result  <- Reduce('+', simr.results)
  return(simr.result/50)
}

###################################################
## generate data for random graph for binomial DAG  
###################################################



randomBinomialDag <- function(RBD.p, RBD.at, RBD.T, RBD.d=2, RBD.fixdegree=T,RBD.seed){
  # at is the number of nodes in each layer
  # T is the number of layer minus 1
  # d is the number of parents for each node except for the roots
  RBD.numberLayer <- RBD.T + 1
  RBD.indegree <- RBD.d
  ### generate node in the first layer
  RBD.adjaMatrix <- matrix(0, RBD.p, RBD.p)
  RBD.layer <- list()
  set.seed(RBD.seed)
  RBD.layer[[1]] <- sample(RBD.p, RBD.at[1], replace = FALSE)
  
  #sample does not work properly when choosing from sets with one element, we deal with A1
  RBD.roots <- RBD.layer[[1]]
  RBD.left <- seq(RBD.p)[-RBD.roots]
  set.seed(RBD.seed)
  RBD.layer[[2]] <- sample(RBD.left, RBD.at[2], replace = FALSE)
  ### generate the directed structure between layer A0 and A1
  RBD.numberParents <- RBD.d
  set.seed(RBD.seed)
  RBD.Parents <- sample(x = RBD.roots, size = RBD.d, replace = FALSE)
  RBD.node <- RBD.layer[[2]]
  RBD.adjaMatrix[RBD.Parents, RBD.node] <- rep(1, length(RBD.Parents)*RBD.at[2])
  ### Note that for any nodes in A1, they share the same parent
  
  
  # from A2 layer 
  for (t in 3:RBD.numberLayer) {
    RBD.numberNode <- RBD.at[t]
    RBD.possibleParents <- unlist(RBD.layer)
    RBD.left <- seq(RBD.p)[-RBD.possibleParents]
    set.seed(RBD.seed)
    RBD.layer[[t]] <- sample(RBD.left, RBD.at[t], replace = FALSE)
    
    for (j in 1:RBD.numberNode) {
      if(RBD.fixdegree==T){
        RBD.numberParents <- RBD.indegree   ## 固定indegree
      }else{
        #draw <- c(d : length(possibleParents))  
        RBD.draw <- c(min(length(RBD.possibleParents), RBD.d): max(length(RBD.possibleParents), RBD.d))
        set.seed(RBD.seed)
        RBD.numberParents <- sample(x = RBD.draw, 1 )  ## random indegree 
      }
      set.seed(RBD.seed)
      RBD.Parents <- sample(x  = RBD.possibleParents, size = RBD.numberParents, replace = FALSE)
      RBD.child <- RBD.layer[[t]][j]
      RBD.adjaMatrix[RBD.Parents, RBD.child] <- rep(1, RBD.numberParents)
    }
  }
  return(list(DAG = RBD.adjaMatrix, layer = RBD.layer))
}


##### generate data #########################################################

get_DAGdata <- function(getD.n, getD.p, getD.at, getD.T, getD.d, fgetD.ixdegree , getD.theta1, getD.theta2, getD.N, getD.seed){
  # N is size for binomial parameter
  set.seed(getD.seed)
  getD.dat <- randomBinomialDag(getD.p, getD.at, getD.T, getD.d, fgetD.ixdegree , getD.seed)
  getD.DAG <- getD.dat$DAG
  getD.layer <- getD.dat$layer
  set.seed(getD.seed)
  getD.theta <- runif(getD.p, getD.theta1, getD.theta2)
  getD.beta  <- matrix(runif(getD.p*getD.p, getD.theta1, getD.theta2), getD.p, getD.p)
  getD.X <- matrix(0, getD.n, getD.p)
  
  ## generate data for roots
  for (i in 1:getD.at[1]) {
    getD.A0node <- getD.layer[[1]]
    getD.root <- getD.A0node[i]
    #X[, root] <- rpois(n, theta[root])
    getD.g.root <- exp(getD.theta[i])/(1+exp(getD.theta[i])) ### probablity should be logit^{-1}
    set.seed(getD.seed)
    getD.X[, getD.root] <- rbinom(getD.n, size = getD.N, prob = getD.g.root) 
  }
  ## generate data for A1 layer
  for (k in 2:(getD.T+1)) {
    
    for (j in 1:getD.at[k]) {
      getD.node <- getD.layer[[k]]
      getD.index <- getD.node[j]
      getD.parents <- which(getD.DAG[, getD.node[j]] != 0)
      
      if(sum(getD.DAG[, getD.index] != 0) == 1){
        #g[j] <- exp(theta[index] + X[, parents]*beta[index, parents])/(1+exp(theta[index] + X[, parents]*beta[index, parents]))
        getD.g <- exp(getD.theta[getD.index] + getD.X[, getD.parents]*getD.beta[getD.index, getD.parents])/(1+exp(getD.theta[getD.index] + getD.X[, getD.parents]*getD.beta[getD.index, getD.parents]))
        set.seed(getD.seed)
        getD.X[, getD.index] <- rbinom(getD.n, size = getD.N, prob = getD.g)  
      }else{
        getD.g<- exp(getD.theta[getD.index] + getD.X[, getD.parents] %*% getD.beta[getD.index, getD.parents])/(1+exp(getD.theta[getD.index] + getD.X[, getD.parents] %*% getD.beta[getD.index, getD.parents]))
        set.seed(getD.seed)
         getD.X[, getD.index] <- rbinom(getD.n, size = getD.N, prob = getD.g)
      }
    }
  }
  ##这里需要check生成的data X 是否变量下的每一类只有一个观测值, remove these observations
  #if(length(onefactor(X)) != 0){
  #  X <- X[-onefactor(X), ]
  #}else{
  #  X <- X
  #}
  ## convert binomial to binary case
  getD.X <- sapply(1:getD.p, function(i) ConvBinomial(getD.X[, i], getD.N))
  return(list(X = getD.X, DAG = getD.DAG, layer = getD.layer))
}


get_DAGstru <- function(get.n, get.p, get.at, get.T, get.d, get.fixdegree , get.theta1, get.theta2, get.N, get.seed){
  get.contin <- TRUE
  while(get.contin){
    get.DAGdata <- get_DAGdata(get.n, get.p, get.at, get.T, get.d, get.fixdegree , get.theta1, get.theta2, get.N, get.seed)
    #DAGdata<-get_DAGdata(n=200, p=20, at=c(10,5,3,2), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
    #DAGdata <- get_DAGdata(n=200, p=100, at=c(70, 20, 5, 5), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
    get.X <- get.DAGdata$X
    get.zeroNumber <- apply(get.X, 2, function(i)sum(i==0))
    get.oneNumber <- apply(get.X, 2, function(i)sum(i==1))
    #min(zeroNumber) >= 2 & min(oneNumber) >= 2
    #zeroNumber;oneNumber
    #if(min(zeroNumber) >= 8 & min(oneNumber) >= 8){  ## p=5可以设置min(Number)>=8 以此满足glmnet函数的要求 
    if(min(get.zeroNumber) > 2 & min(get.oneNumber) > 2){
      get.contin <- FALSE
    }
  }
  return(get.DAGdata)
}


## transform each binomial variable to binary ones
ConvBinomial <- function(conv.X,conv.N){
  conv.n <- length(conv.X)
  conv.vX <- rep(0,conv.n*conv.N)
  #con.vX[sample(1:(conv.N*conv.n), sum(conv.X), replace = F)]<- 1
  for (i in 1:conv.n) {
   conv.vX <- c(conv.vX, c(rep(0, conv.N-conv.X[i]), rep(1, conv.X[i])))
  }
  #conv.X <- conv.vX
  return(conv.vX)
}


# To check whether one class has 1 observation
#onefactor <- function(one.X){
#  one.n <- dim(one.X)[1]
#  one.p <- dim(one.X)[2]
#  one.classNumber <- apply(one.X, 2, unique)
#  one.loca <- list()
#  one.subloca <- list()
#  for(pp in 1:one.p){
#    for (m in 1:length(one.classNumber[[pp]])) {
#      onefactor <- which(one.X[, pp] == one.classNumber[[pp]][m])
      
#      if(length(onefactor) == 1){
#        one.subloca[[m]] <- onefactor
#      } 
#    }
#    one.loca[[pp]] <- one.subloca
#  }
#  one.oneObser <- unique(unlist(one.loca))
#  return(one.oneObser)
#}


################################################
############# Experiments  #####################
################################################



# XX<-get_DAGstru(n=200, p=5, at=c(2,1,2), T=2, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# XX<-get_DAGstru(n=200, p=10, at=c(2,3,5), T=2, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# XX<-get_DAGstru(n=200, p=20, at=c(2,3,5,10), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# XX<-get_DAGstru(n=200, p=50, at=c(5,10,15,20), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# XX<-get_DAGstru(n=200, p=50, at=c(30,10,5,5), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# XX<-get_DAGstru(n=200, p=100, at=c(10,20,30,40), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# XX<-get_DAGstru(n=200, p=100, at=c(70,10,10,10), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4) 

# X <- XX$X
# truth <- XX$DAG 
# layer <- XX$layer

# N<-4
# beta_1=1; beta_2=-1/N
# edag <- EDAG(X, grid.t, beta_1=1, beta_2=-1/N)
# ods <- ODSGLMLasso(X, beta_1=1, beta_2=-1/N)


# do_one(n=200, p=5, grid.t, at=c(2,1,2), T=2, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# do_one(n=200, p=10, grid.t, at=c(2,1,2,5), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# do_one(n=200, p=20, grid.t, at=c(10,5,3,2), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# do_one(n=200, p=50, grid.t, at=c(10,10,10,10,10), T=4, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# do_one(n=200, p=100, grid.t, at=c(70,10,10,10), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)

# simulation_randomBino(n = 200, p = 5, grid.t, simu.times = 50, at=c(2,1,2), T=2, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# simulation_randomBino(n = 300, p = 5, grid.t, simu.times = 50, at=c(2,1,2), T=2, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# simulation_randomBino(n = 500, p = 5, grid.t, simu.times = 50, at=c(2,1,2), T=2, d=2, fixdegree = T, t1=0.5, t2=1, N=4)

# simulation_randomBino(n = 200, p = 20, grid.t, simu.times = 50, at=c(10,5,3,2), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# simulation_randomBino(n = 300, p = 20, grid.t, simu.times = 50, at=c(10,5,3,2), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# simulation_randomBino(n = 500, p = 20, grid.t, simu.times = 50, at=c(10,5,3,2), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)

# simulation_randomBino(n = 200, p = 100, grid.t, simu.times = 50, at=c(70,20,5,5), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# simulation_randomBino(n = 300, p = 100, grid.t, simu.times = 50, at=c(70,20,5,5), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
# simulation_randomBino(n = 500, p = 100, grid.t, simu.times = 50, at=c(70,20,5,5), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)

grid.t <- 10^(-2+0.15*seq(0,60,by=1))[-1] #####?? why [-1]

## p =5 replicate 50 times

simu_replicate(200, 5, grid.t, c(2,1,2), 2, 2,  T, 0.5, 1, 4)

simu_replicate(n = 200, p = 5, grid.t, at=c(2,1,2), T=2, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
simu_replicate(n = 300, p = 5, grid.t, at=c(2,1,2), T=2, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
simu_replicate(n = 500, p = 5, grid.t, at=c(2,1,2), T=2, d=2, fixdegree = T, t1=0.5, t2=1, N=4)


## p = 20  running for 50 times
simu_replicate(n = 200, p = 20, grid.t, at=c(10,5,3,2), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
simu_replicate(n = 300, p = 20, grid.t, at=c(10,5,3,2), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
simu_replicate(n = 500, p = 20, grid.t, at=c(10,5,3,2), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)


## p = 100  running for 50 times
simu_replicate(n = 200, p = 100, grid.t, at=c(70,20,5,5), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
simu_replicate(n = 300, p = 100, grid.t, at=c(70,20,5,5), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)
simu_replicate(n = 500, p = 100, grid.t, at=c(70,20,5,5), T=3, d=2, fixdegree = T, t1=0.5, t2=1, N=4)


