###########################################
### TLDAG for the random binomial graph ###
###########################################

source("..GDS/util_DAGs/computeScoreSEMGauss.R")
source("../GDS/inferDAG/gesWrap.R")
source("../R/helper_function.R")
source("../R/DLiNGAM_p.R")
source("helper_function.R")
sourceCpp('direct_lingam_funcs.cpp')
library(bnlearn)
library(pcalg)
library(Rcpp)
library(RcppArmadillo)


########################
### calculate kappa  ###
########################


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
Kappa.cal <- function(KC.layer, KC.X, KC.n, KC.p, KC.beta_1, KC.beta_2, KC.N, KC.seed){
  KC.n1 <- as.integer(KC.n/2)
  KC.n2 <- KC.n - KC.n1
  KC.done <- unlist(KC.layer) ## The set S
  set.seed(KC.seed) #### Add seed here
  KC.sample <- sample(1:KC.n, KC.n, replace = FALSE)
  KC.sinkNodes <- which(colMeans(KC.X) == KC.N)
  data.perm <- KC.X[KC.sample, ]
  
  KC.data1 <- data.perm[1:KC.n1, ]
  KC.data2 <- data.perm[-(1:KC.n1), ]
  KC.XBino1 <- sapply(1:KC.p, function(i) ConvBinomial(KC.data1[, i], KC.N))
  KC.XBino2 <- sapply(1:KC.p, function(i) ConvBinomial(KC.data2[, i], KC.N))
  
  if(length(KC.done) == 0){ ## A_0 for roots
    KC.set1 <- layer_A0(KC.data1, KC.beta_1, KC.beta_2)
    KC.set2 <- layer_A0(KC.data2, KC.beta_1, KC.beta_2)
  }else if( length(KC.done) >=1 & length(KC.done) < (KC.p-1)){  ## layer A_1 to A_{T-1}
    KC.set1 <- layer_At(KC.layer, KC.data1, KC.XBino1, KC.beta_1, KC.beta_2, KC.N)
    KC.set2 <- layer_At(KC.layer, KC.data2, KC.XBino2, KC.beta_1, KC.beta_2, KC.N)
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
Kappa <- function(Ka.layer, Ka.X, Ka.grid.t, Ka.beta_1, Ka.beta_2, Ka.N){
  Ka.n <- dim(Ka.X)[1]
  Ka.p <- dim(Ka.X)[2]
  Ka.done<-unlist(Ka.layer)
  #Ka.nodes <- seq(Ka.p)[-Ka.done]
  #Ka.sinkNodes <- which(colMeans(Ka.X) == Ka.N)
  # print(Ka.done)
  
  if(is.null(Ka.done) == T){
    Ka.left <- seq(Ka.p)
  }else{
    Ka.left <- seq(Ka.p)[-Ka.done]
  }
  
  Kappa.temp<-sapply(c(1,2,3,4,5), Kappa.cal, KC.layer=Ka.layer, KC.X=Ka.X, KC.n=Ka.n, KC.p=Ka.p, KC.beta_1=Ka.beta_1, KC.beta_2=Ka.beta_2, KC.N=Ka.N)
  kappa.temp1<-matrix(c(Kappa.temp[[1]],Kappa.temp[[3]],Kappa.temp[[5]],Kappa.temp[[7]],Kappa.temp[[9]]),5, length(Kappa.temp[[1]]),byrow = T)
  kappa.temp2<-matrix(c(Kappa.temp[[2]],Kappa.temp[[4]],Kappa.temp[[6]],Kappa.temp[[8]],Kappa.temp[[10]]), 5, length(Kappa.temp[[2]]),byrow = T)
  Ka.result<-sapply(Ka.grid.t, Kapp.compute, kap.re1=kappa.temp1, kap.re2=kappa.temp2, kap.left=Ka.left)
  return(Ka.result/5)
}




TuningPara <- function(Tune.layer, Tune.X, Tune.grid.t, Tune.alpha, Tune.beta_1, Tune.beta_2, Tune.N){
  kappa.vec <- Kappa(Tune.layer, Tune.X, Tune.grid.t, Tune.beta_1, Tune.beta_2, Tune.N)
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



##############################
## roots for the layer A_0  ##
##############################

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


########################################################
## calculate the ratio for the nodes in the layer A_t ##
########################################################

layer_At <-function(At.layer, At.X, At.XBino, At.beta_1, At.beta_2, At.N){
  At.n1 <- dim(At.XBino)[1]
  At.p <- dim(At.XBino)[2]
  At.Sdone <- unlist(At.layer)
  At.left <- seq(At.p)[-At.Sdone]
  At.n.left <- length(At.left)
  At.con.mean = At.con.var = c();
  for (k in 1:At.n.left) {
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
    w.Tt <- (At.beta_1 + At.beta_2 * con.mean.X)^(-1)
    At.Tt <- At.X[, At.left][, k] * w.Tt
    At.con.mean[k] <- con.prob * w.Tt * sum(At.Tt)
    
    sec.moment <- con.prob * w.Tt * sum(At.Tt^2)
    At.con.var[k] <- sec.moment - At.con.mean[k]^2
  }
  
  At.ratio <- At.con.var / At.con.mean
  At.final.result <- At.ratio / min(At.ratio)
  return(At.final.result) 
}



###################################################
######  the main function for the proposed method
###################################################

EDAG <- function(ED.X, ED.grid.t, ED.beta_1, ED.beta_2, ED.N){
  ptm <- proc.time() 
  ED.n <- dim(ED.X)[1]
  ED.p <- dim(ED.X)[2]
  ED.XBino <- sapply(1:ED.p, function(i) ConvBinomial(ED.X[, i], ED.N))
  ED.TF <- matrix(0, ED.n, ED.p)
  ED.done <- NULL
  ED.layer <- list()
  ED.X.step <- list()
  ED.ratio.t <- list()
  ED.thres <- list()
  ED.final.layer<-list()
  
  #######################################
  # Step1: layer and ordering
  # roots layer A_0
  ############## Node identification ####
  
  ED.A0.thresh <- TuningPara(ED.layer, ED.X, ED.grid.t, 0.1, ED.beta_1, ED.beta_2, ED.N)$thres
  ED.ratio.final<-layer_A0(ED.X, ED.beta_1, ED.beta_2)
  ED.layer[[1]] <- which(abs(ED.ratio.final - 1) <= ED.A0.thresh)   ## the index in the original variables
  #ED.layer[[1]] <- ED.nodes[ED.A0.index]     ## the nodes estimated in the layer A_0
  #ED.X.step[[1]] <- ED.nodes[-ED.A0.index]   ## the nodes remained to be estimated after the layer A_0
  ED.thres[[1]] <- ED.A0.thresh
  ED.ratio.t[[1]]<-ED.ratio.final
  
  #######################################
  ## identify nodes in layer A_t
  
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
          ED.AT.thresh <- TuningPara(ED.layer, ED.X, ED.grid.t, 0.1, ED.beta_1, ED.beta_2, ED.N)$thres
          ED.ATratio.final<-layer_At(ED.layer, ED.X, ED.XBino, ED.beta_1, ED.beta_2, ED.N)
          ED.layer[[lay]] <- ED.left[which(abs(ED.ATratio.final - 1) <= ED.AT.thresh)]
          ED.thres[[lay]] <- ED.AT.thresh
          ED.ratio.t[[lay]]<-ED.ATratio.final
          
          
          #########################################################################
          
        }else if(length(unlist(ED.layer)) == ED.p-1){
          ## sink node is in A_T 
          #ED.left <- ED.X.step[[lay - 1]][-unlist(ED.layer)]
          ED.left <- seq(ED.p)[-unlist(ED.layer)]
          ED.n.left <- length(ED.left)
          ED.layer[[length(ED.layer) + 1]] <- ED.left
          
        }else (length(unlist(ED.layer)) >= ED.p)
        cont <- FALSE
        #break
      }
    }
    
    ED.final.layer <- list(ED.layer, ED.ratio.t, ED.thres)  ##list 其中layer[[1]]表示topological layer
    
    ###############################################
    ##Step2:Sparse Estimation for Direction Edges
    ###############################################
    
    EDS.layernodes <- rev(ED.final.layer[[1]]) ### Selected layer structure with nodes in each layer
    EDS.layernum  <- length(ED.final.layer[[1]]) ### number of layer in the DAG
    #EDS.XBino <- ED.XBino
    EDS.result <- matrix(0, ED.p, ED.p)
    
    for (ii in 1:(EDS.layernum-1)) {   ###
      for (k in 1:length(EDS.layernodes[[ii]])) {
        EDS.now <- EDS.layernodes[[ii]][k]  ## each element in each layer ### the node as a child
        ## now turn to reconstruct the directed structure among SE.now and all the nodes in the upper layers
        
        #for (j in (ii+1):EDS.layernum) {
        
        #EDS.this <- EDS.layernodes[[j]]  ## each layer ### should be the nodes in the upper layer
        EDS.this <- potentialParents(EDS.layernodes, ii)
        
        if(length(EDS.this) > 1){
          EDS.nfold <- 5
          EDS.fold0 <- sample.int(sum(ED.XBino[, EDS.now]==0)) %% EDS.nfold
          EDS.fold1 <- sample.int(sum(ED.XBino[, EDS.now]==1)) %% EDS.nfold
          EDS.foldid <- numeric(length(ED.XBino[, EDS.now]))
          EDS.foldid[ED.XBino[, EDS.now]==0] <- EDS.fold0
          EDS.foldid[ED.XBino[, EDS.now]==1] <- EDS.fold1
          EDS.foldid <- EDS.foldid + 1
          EDS.fraction <- table(ED.XBino[, EDS.now])/length(ED.XBino[, EDS.now])
          EDS.weights <- 1 - EDS.fraction[as.character(ED.XBino[, EDS.now])]
          
          EDS.lassom <- glmnet::cv.glmnet(ED.XBino[, EDS.this], ED.XBino[, EDS.now], family = "binomial", foldid = EDS.foldid, weights = EDS.weights )
          EDS.bfit <- coefficients(EDS.lassom)[-1]
          for (mm in 1:length(EDS.this)) {
            if(EDS.bfit[mm] != 0)
              EDS.result[EDS.this[mm], EDS.now] <- 1
          }
        }else{
          EDS.lmod <- glm(ED.XBino[, EDS.now] ~ ED.XBino[, EDS.this], family = "binomial", control=list(maxit=100))
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






##############################################
### ODS(GLMLasso) in Park and Raskutti (2018)
##############################################

ODSGLMLasso <- function(X, beta_1, beta_2, N){
  ptm<-proc.time()
  n <- dim(X)[1]
  p <- dim(X)[2]
  XBino <- sapply(1:p, function(i) ConvBinomial(X[, i], N))
  
  ##########################
  #Step1: NormalizedGraph
  ##########################
  
  NormalizedGraph <- matrix(0, p, p)
  for (i in 1:p) {
    left <- seq(p)[-i]
    
    nfold <- 5
    fold0 <- sample.int(sum(XBino[, i]==0)) %% nfold
    fold1 <- sample.int(sum(XBino[, i]==1)) %% nfold
    foldid <- numeric(length(XBino[, i]))
    foldid[XBino[, i]==0] <- fold0
    foldid[XBino[, i]==1] <- fold1
    foldid <- foldid + 1
    fraction <- table(XBino[, i])/length(XBino[, i])
    weights <- 1 - fraction[as.character(XBino[, i])]
    
    lassom <- glmnet::cv.glmnet(XBino[, left], XBino[, i], family = "binomial", foldid = foldid, weights = weights)
    bfit <- coefficients(lassom)[-1]
    for (j in 1:length(bfit)) {
      if(bfit[j] != 0)
        NormalizedGraph[i, left[j]] <- 1
    }
  }
  ## return(NormalizedGraph)
  
  ####################
  #Step2: ordering
  ####################
  
  ODS.order <- c()
  #ODS.left <- list()
  
  ## first element of the ordering
  con.mean1 <- apply(X, 2, mean)
  con.var1  <- apply(X, 2, var)
  w1 <- (beta_1 + beta_2* con.mean1)^(-1)
  S1 <- w1^2 * con.var1 - w1 * con.mean1
  ODS.order[1] <- seq(p)[which.min(S1)]
  #ODS.left[[1]] <- ODS.nodes[-which.min(S1)]
  
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
      if(length(C[[k]]) == 0 ){  
        con.mean <- mean(X[, k])
        con.var  <- var(X[, k])
        
        w <- (beta_1 + beta_2 * con.mean)^(-1)
        kk <- which(left.order == k)
        S[kk] <- w^2 * con.var - w * con.mean
        
      }else{
        
        k_index <- length(C[[k]])
        
        glmfit <- glm(XBino[, k] ~ XBino[, C[[k]]], family = "binomial", control=list(maxit=100))
        glmfit.coef <- as.matrix(glmfit$coefficients)
        
        if(any(is.na(glmfit.coef))){
          na.index <- which(is.na(glmfit.coef) == TRUE)
          glmfit.coef[na.index] <- 0
        }
        
        glm_St <- cbind(XBino[, k], XBino[, C[[k]]])
        exp.term <- exp(glm_St %*% glmfit.coef)
        con.prob <- mean(exp.term^(XBino[, k1]) * (factorial(1)/(factorial(XBino[, k1])*factorial(1-XBino[, k1])))
                         / (1 + exp.term)^(XBino[, k1]))
        con.mean <-  con.prob * sum(X[, k1])
        sec.moment <- con.prob * sum(X[, k1]^2)
        con.var <- sec.moment - con.mean^2
        
        w <- (beta_1 + beta_2 * con.mean)^(-1)
        nx <- n*length(C[[k]])
        nC <- sum( nx * ifelse(nx >= c0*n, 1, 0) )   
        kk <- which(left.order == k)
        S[kk] <- sum( (w^2 * con.var - w * con.mean) * nx / nC ) 
      }   
    }
    
    ODS.order[m] <- left.order[which.min(S)]
    #ODS.left[[m]] <- ODS.left[[m-1]][-which.min(S)]
  }
  
  ODS.order[p] <- seq(p)[-ODS.order]
  ## return(ODS.order)
  
  ###########################
  ##Step3: directed edges
  ###########################
  
  ordering <- ODS.order
  rr <- rev(ordering)
  result <- matrix(0, p, p)
  
  for (ii in 1:(p-1)){
    now <- rr[ii]
    this <- sort(rr[(ii+1):p])
    if (length(this) > 1){
      nfold <- 5
      fold0 <- sample.int(sum(XBino[, now]==0)) %% nfold
      fold1 <- sample.int(sum(XBino[, now]==1)) %% nfold
      foldid <- numeric(length(XBino[, now]))
      foldid[XBino[, now]==0] <- fold0
      foldid[XBino[, now]==1] <- fold1
      foldid <- foldid + 1
      fraction <- table(XBino[, now])/length(XBino[, now])
      weights <- 1 - fraction[as.character(XBino[, now])]
      
      lassom <- glmnet::cv.glmnet(XBino[, this], XBino[, now], family = "binomial",foldid = foldid, weights = weights)
      bfit <- coefficients(lassom)[-1]
      
      for (jj in 1:length(this)){
        if( bfit[jj] != 0 )
          result[this[jj],now] <- 1
      }
    }else {
      # deal with the last two nodes
      lmod <- glm(XBino[, now] ~ XBino[, this], family = "binomial", control = list(maxit=100))
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

## N: Binomial parameter
do_one <- function(do.n, do.p, do.grid.t, do.T, do.graph=c('dense','sparse'), do.theta1, do.theta2, do.N, do.seed){
  do.result <- matrix(0, 6, 5) ##
  # Data generation
  #dat <- hubgraph(n, p)
  #do.dat <- get_DAGstru(do.n, do.p, do.at, do.T, do.d, do.fixdegree, do.theta1, do.theta2, do.N, do.seed) # N is the size for binomial
  do.dat <- get_DAGdata(do.n, do.p, do.T, do.graph, do.theta1, do.theta2, do.N, do.seed) # N is the size for binomial
  do.truth <- do.dat$DAG
  do.X <- do.dat$X
  
  
  # EDAG
  ptm <- proc.time()
  do.sresult <- EDAG(do.X, do.grid.t, 1, -1/do.N, do.N)
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
  sresult <- ODSGLMLasso(do.X, beta_1 = 1, beta_2 = -1/do.N, do.N)
  sx <- sresult$adj 
  do.result[1, 2] <- (proc.time() - ptm)[1]
  do.result[2, 2] <- hammingDistance(sx, do.truth)/(do.p*(do.p-1))     ##normalized hammdis
  do.result[3, 2] <- ifelse(sum(do.truth), sum(do.truth*sx)/sum(do.truth), 0) #recall
  do.result[4, 2] <- ifelse(sum(do.truth), sum(do.truth*sx)/(sum(sx)+1e-6), 0) # precision
  do.result[5, 2] <- 2*do.result[3, 2]*do.result[4,2]/(do.result[3, 2]+do.result[4,2]+1e-6)  # F1-score 
  do.result[6, 2] <- ifelse(sum(sx), 1 - sum(do.truth*sx)/(sum(sx)+1e-6), 0) #FDR
  
  
  #DirectLINGAM
  ptm <- proc.time()
  lx <- DirectLINGAM(do.X)
  do.result[1, 3] <- (proc.time() - ptm)[1]
  do.result[2, 3] <- hammingDistance(lx, do.truth)/(do.p*(do.p-1))     ##normalized hammdis
  do.result[3, 3] <- ifelse(sum(do.truth), sum(do.truth*lx)/sum(do.truth), 0) #recall
  do.result[4, 3] <- ifelse(sum(do.truth), sum(do.truth*lx)/(sum(lx)+1e-6), 0) # precision
  do.result[5, 3] <- 2*do.result[3, 3]*do.result[4,3]/(do.result[3, 3]+do.result[4,3]+1e-6)  # F1-score 
  do.result[6, 3] <- ifelse(sum(lx), 1 - sum(do.truth*lx)/(sum(lx)+1e-6), 0) #FDR
  
  
  # GES
  ptm <- proc.time()
  gx.log <- gesWrap(do.X)$Adj
  gx <- matrix(as.numeric(gx.log), do.p, do.p)
  do.result[1, 4] <- (proc.time() - ptm)[1]
  do.result[2, 4] <- hammingDistance(gx, do.truth)/(do.p*(do.p-1))     ##normalized hammdis
  do.result[3, 4] <- ifelse(sum(do.truth), sum(do.truth*gx)/sum(do.truth), 0) #recall
  do.result[4, 4] <- ifelse(sum(do.truth), sum(do.truth*gx)/(sum(gx)+1e-6), 0) # precision
  do.result[5, 4] <- 2*do.result[3, 4]*do.result[4,4]/(do.result[3, 4]+do.result[4,4]+1e-6)  # F1-score 
  do.result[6, 4] <- ifelse(sum(gx), 1 - sum(do.truth*gx)/(sum(gx)+1e-6), 0) #FDR
  
  
  #MMHC
  ptm <- proc.time()
  mmhc.cha <- mmhc(data.frame(do.X))$arcs
  mx <- charaToNume(mmhc.cha, do.p)
  do.result[1, 5] <- (proc.time() - ptm)[1]
  do.result[2, 5] <- hammingDistance(mx, do.truth)/(do.p*(do.p-1))     ##normalized hammdis
  do.result[3, 5] <- ifelse(sum(do.truth), sum(do.truth*mx)/sum(do.truth), 0) #recall
  do.result[4, 5] <- ifelse(sum(do.truth), sum(do.truth*mx)/(sum(mx)+1e-6), 0) # precision
  do.result[5, 5] <- 2*do.result[3, 5]*do.result[4, 5]/(do.result[3, 5]+do.result[4, 5]+1e-6)  # F1-score 
  do.result[6, 5] <- ifelse(sum(mx), 1 - sum(do.truth*mx)/(sum(mx)+1e-6), 0) #FDR
  
  #return(do.result)
  return(list(DAG=do.truth, EDAG=do.sresult, ODS=sresult, Perf=do.result))
}



####### repeat ################################


simu_replicate <- function(simr.n, simr.p, simr.times, simr.grid.t, simr.T, simr.graph=c('dense','sparse'), simr.theta1, simr.theta2, simr.N){
  simr.results <- list()
  simr.sd <- list()  ## standard deviation
  
  for (i in 1:simr.times) {
    print(i)
    #results[[i]] <- do_one(n = n, p = p, grid.t)
    simr.results[[i]] <- do_one(simr.n, simr.p, simr.grid.t, simr.T, simr.graph, simr.theta1, simr.theta2, simr.N, i)$Perf
  }   
  simr.result  <- Reduce('+', simr.results)
  Metric = simr.result/simr.times   ## mean of all the metrics 
  simr.sd <- apply(array(unlist(simr.results), c(6, 5, simr.times)), c(1,2), sd)  ## standard deviation using 3Darray
  
  return(list(Metrics = round(Metric, 2), SD = round(simr.sd, 2)))
  
}


###################################################
## generate data for random graph for binomial DAG  
###################################################



randomBinomialDag <- function(RBD.p, RBD.T, RBD.graph=c('dense','sparse'), RBD.seed){
  # T is the number of layer minus 1
  # d is the number of parents for each node except for the roots
  
  RBD.numberLayer <- RBD.T + 1
  
  ## the number of nodes in each layer is randomly generated
  RBD.layer <- list()
  set.seed(RBD.seed)
  RBD.layer[1] <- sample(RBD.p - RBD.T, 1, replace = FALSE)
  for (i in 2:RBD.T){
    RBD.done <- sum(unlist(RBD.layer))
    RBD.left <- RBD.numberLayer - length(RBD.layer) -1
    RBD.sample <- RBD.p - RBD.done - RBD.left
    set.seed(RBD.seed)
    RBD.layer[i] <- sample(RBD.sample, 1, replace = FALSE)
  }
  RBD.layer[RBD.numberLayer] <- RBD.p - sum(unlist(RBD.layer))
  
  
  ## sample does not work properly when choosing from sets with one element, we need to deal with it carefully
  
  ### generate node in A0 layer
  RBD.adjaMatrix <- matrix(0, RBD.p, RBD.p)
  RBD.layerNodes <- list()
  set.seed(RBD.seed)
  RBD.layerNodes[[1]] <- sample(RBD.p, RBD.layer[[1]], replace = FALSE)
  
  
  for (t in 2:RBD.numberLayer) {
    RBD.done <- unlist(RBD.layerNodes)
    RBD.left <- seq(RBD.p)[-RBD.done]
    if(t < RBD.numberLayer){
      set.seed(RBD.seed)
      RBD.layerNodes[[t]] <- sample(RBD.left, RBD.layer[[t]], replace = FALSE)
    }else{
      RBD.layerNodes[[t]] <- RBD.left
    }
    
    for (k in 1:RBD.layer[[t]]) {
      
      if(RBD.graph == 'sparse'){
        RBD.numberParents <- ifelse(length(RBD.done) == 1, 1, length(RBD.done)/2)
      }else{
        RBD.numberParents <- length(RBD.done)
      }
      
      RBD.Parents <- sample(x  = RBD.done, size = RBD.numberParents, replace = FALSE)
      RBD.child <- RBD.layerNodes[[t]][k]
      RBD.adjaMatrix[RBD.Parents, RBD.child] <- rep(1, RBD.numberParents)
    }
  }
  
  RBD.edges <- sum(RBD.adjaMatrix == 1) 
  return(list(DAG = RBD.adjaMatrix, layer = RBD.layer, Nodes = RBD.layerNodes, EdgesNumber = RBD.edges))
  
}



##############################
##### generate data ##########

get_DAGdata <- function(getD.n, getD.p, getD.T, getD.graph=c('dense',"sparse"), getD.theta1, getD.theta2, getD.N, getD.seed){
  # N is size for binomial parameter
  set.seed(getD.seed)
  getD.dat <- randomBinomialDag(getD.p, getD.T, getD.graph, getD.seed)
  getD.DAG <- getD.dat$DAG
  getD.layer <- getD.dat$layer
  getD.nodes <- getD.dat$Nodes
  getD.edges <- getD.dat$EdgesNumber
  set.seed(getD.seed)
  getD.theta <- runif(getD.p, getD.theta1, getD.theta2)
  getD.beta  <- matrix(runif(getD.p*getD.p, getD.theta1, getD.theta2), getD.p, getD.p)
  getD.X <- matrix(0, getD.n, getD.p)
  
  ## generate data for roots
  for (i in 1:getD.layer[[1]]) {
    getD.A0node <- getD.nodes[[1]]
    getD.root <- getD.A0node[i]
    #X[, root] <- rpois(n, theta[root])
    getD.g.root <- exp(getD.theta[i])/(1+exp(getD.theta[i])) ### probablity should be logit^{-1}
    set.seed(getD.seed)
    getD.X[, getD.root] <- rbinom(getD.n, size = getD.N, prob = getD.g.root) 
  }
  ## generate data for A1 layer
  for (k in 2:(getD.T+1)) {
    
    for (j in 1:getD.layer[[k]]) {
      getD.node <- getD.nodes[[k]]
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
  
  return(list(X = getD.X, DAG = getD.DAG, layer = getD.layer, nodes = getD.nodes, edges = getD.edges))
  
}











