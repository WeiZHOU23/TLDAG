##############################################
### TLDAG for dense hub graph for poisson  ###
##############################################



source("../GDS/util_DAGs/computeScoreSEMGauss.R")
source("../GDS/inferDAG/gesWrap.R")
source("../R/DLiNGAM_p.R")
source("../R/helper_function.R")
sourceCpp('direct_lingam_funcs.cpp')
library(bnlearn)
library(pcalg)
library(Rcpp)
library(RcppArmadillo)



################################################################################
###### To calculate kappa value to determine the tuning parameter \alpha
################################################################################

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
## roots for the layer A_0  

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
## the ratio values for the node in the layer A_t

layer_At <-function(At.layer, At.X, At.beta_1, At.beta_2){
  At.n1 <- dim(At.X)[1]
  At.p <- dim(At.X)[2]
  At.Sdone <- unlist(At.layer)
  At.left <- seq(At.p)[-At.Sdone]
  At.n.left <- length(At.left)
  At.con.mean = At.con.var = c();
  for (k in 1:At.n.left) {
    glmfit <- glm(At.X[, At.left][, k] ~ At.X[, At.Sdone], family = "poisson", control=list(maxit=100))
    con.mean1Pois <- mean(predict(glmfit))
    wPois <- (At.beta_1 + At.beta_2 * con.mean1Pois)^(-1)
    TtPois <- At.X[, At.left][, k] * wPois
    Tt2Pois <- (At.X[, At.left][, k] * wPois)^2
    
    At.con.mean[k] <- mean(TtPois)    
    
    fitT2Pois <- glm(Tt2Pois ~ At.X[, At.Sdone], family = "poisson", control=list(maxit=100)) 
    sec.moment <- mean(exp(predict(fitT2Pois)))
    At.con.var[k] <- sec.moment - At.con.mean[k]^2
    
  }
  At.ratio <- At.con.var / At.con.mean
  final.result <- At.ratio / min(At.ratio)
  return(final.result) 
}





#################################################
## the main function for the proposed method 
#################################################

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
  
  ED.final.layer <- list(ED.layer, ED.ratio.t, ED.thres)  
  
  ############################################################
  ##Step2:Sparse Estimation for Direction Edges
  ############################################################
  
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
        
        EDS.lassom <- glmnet::cv.glmnet(ED.X[, EDS.this], ED.X[, EDS.now], family = "poisson", nfolds = 5)
        EDS.bfit <- coefficients(EDS.lassom)[-1]
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



####################################################
### the ODS algorithm in Park and Raskutti (2018)
####################################################

ODSGLMLasso <- function(X, beta_1, beta_2){
  ptm<-proc.time()
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  ####################################
  #Step1: NormalizedGraph
  ####################################
  
  NormalizedGraph <- matrix(0, p, p)
  eachtime <- list()
  for (i in 1:p) {
    left <- seq(p)[-i]
    ptm <- proc.time()
    lassom <- glmnet::cv.glmnet(X[, left], X[, i], family = "poisson", alpha = 1, nfolds = 5)
    bfit <- coefficients(lassom)[-1]
    for (j in 1:length(bfit)) {
      if(bfit[j] != 0)
        NormalizedGraph[i, left[j]] <- 1
    }
    eachtime[[i]] <- (proc.time() - ptm)[1]
  }
  ## return(NormalizedGraph)
  
  
  ####################################
  #Step2: ordering
  ####################################
  
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
    ODS.order <- ODS.order  #### renew the ordering
    left.order <- seq(p)[-ODS.order]
    S <- c()
    
    for (k1 in 1:length(left.order) ) { # k  the original variable index
      k <- left.order[k1]
      neigh <- union(which( norgraph[k, ] != 0), which( norgraph[, k] != 0))   # k'neighborhood 
      C[[k]] <- intersect(neigh, ODS.order)   # given parents
      
      # calculate conditinal mean/variance  
      if(length(C[[k]]) ==0 ){  
        con.mean <- mean(X[, k])
        con.var  <- var(X[, k])
        
        w <- (beta_1 + beta_2 * con.mean)^(-1)
        kk <- which(left.order == k)
        S[kk] <- w^2 * con.var - w * con.mean
        
      }else{
        glmPois <- glm( X[, k] ~ X[, C[[k]]], family = "poisson", control = list(maxit=100)) 
        glm2Pois <- glm( X[, k]^2 ~ X[, C[[k]]], family = "poisson", control = list(maxit=100))
        con.mean <- mean(exp(predict(glmPois)))
        sec.moment <- mean(exp(predict(glm2Pois)))
        con.var <- sec.moment - con.mean^2
        
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
  
  ########################################
  ##Step3: directed edges
  ########################################
  
  ordering <- ODS.order
  rr <- rev(ordering)
  result <- matrix(0, p, p)
  
  for (ii in 1:(p-1)){
    now <- rr[ii]
    this <- sort(rr[(ii+1):p])
    if (length(this) > 1){
      
      lassom <- glmnet::cv.glmnet(X[, this], X[, now], family = "poisson", alpha = 1, nfolds = 5)
      bfit <- coefficients(lassom)[-1]
      
      for (jj in 1:length(this)){
        if( bfit[jj] != 0 )
          result[this[jj],now] <- 1
      }
    }else {
      # deal with the last two nodes
      lmod <- glm(X[, now] ~ X[, this], family = "poisson", control = list(maxit=100))
      coeff <- as.vector(lmod$coefficients[2])
      if (coeff != 0) 
        result[this,now] <- 1
    }
  }
  return(list(adj=result,ordering=rev(rr)))
}











#############################
########## Experiment
#############################


do_one <- function(do.n, do.p, do.grid.t, do.seed){
  do.result <- matrix(0, 6, 6) ## all metrics
  # Data generation
  #dat <- hubgraph(n, p)
  do.dat <- hubgraph(do.n, do.p, do.seed)
  do.truth <- do.dat$truth
  do.X <- do.dat$X
  
  # TLDAG
  ptm <- proc.time()
  do.sresult <- EDAG(do.X, do.grid.t, 1, 0)
  do.sx <- do.sresult$adj
  do.result[1, 1] <- (proc.time() - ptm)[1]
  do.result[2, 1] <- hammingDistance(do.sx, do.truth)/(do.p*(do.p-1))
  do.result[3, 1] <- ifelse(sum(do.truth), sum(do.truth*do.sx)/sum(do.truth), 0) #recall
  do.result[4, 1] <- ifelse(sum(do.truth), sum(do.truth*do.sx)/(sum(do.sx)+1e-6), 0) # precision
  do.result[5, 1] <- 2*do.result[3, 1]*do.result[4,1]/(do.result[3, 1]+do.result[4,1]+1e-6)  # F1 score
  do.result[6, 1] <- ifelse(sum(do.sx), 1 - sum(do.truth*do.sx)/(sum(do.sx)+1e-6), 0) #FDR
  
  
  
  # ODSGLMLasso
  ptm <- proc.time()
  sresult <- ODSGLMLasso(do.X, beta_1 = 1, beta_2 = 0)
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
  
  #do.result1 <- round(do.result, 2)
  return(do.result)
}





############### repeat ################################


simu_replicate <- function(simr.n, simr.p, simu.times, simr.grid.t){
  simr.results <- list()
  simr.sd <- list()  ## standard deviation
  
  for (i in 1:simu.times) {
    print(i)
    #results[[i]] <- do_one(n = n, p = p, grid.t)
    simr.results[[i]] <- do_one(simr.n, simr.p, simr.grid.t, i)
  }                   
  simr.result  <- Reduce('+', simr.results)
  Metric = simr.result/simu.times   ## mean of all the metrics 
  simr.sd <- apply(array(unlist(simr.results), c(6, 5, simu.times)), c(1,2), sd)  ## standard deviation using 3Darray
  
  return(list(Metrics = round(Metric, 2), SD = round(simr.sd, 2)))
}



#####################################
## generate data for hub graph (d=1)  
#####################################

hubgraph <- function(n, p, h.seed){
  ordering <- seq(p)
  truth <- matrix(0,p,p)
  truth[1, 2:p] <- 1
  
  set.seed(h.seed)
  theta <- runif(p, 1, 3)
  set.seed(h.seed)
  beta  <- matrix(runif(p*p, 0.1, 0.5), p ,p)
  X <- matrix(0, n, p)
  set.seed(h.seed)
  X[, 1] <- rpois(n, lambda = exp(theta[1]))
  
  for (j in 2:p) {
    g <- exp(theta[j] + beta[j, 1] * X[, 1])
    X[, j] <- rpois(n, lambda = g)
  }
  
  return(list(X = X, truth = truth, TO = ordering))
}






