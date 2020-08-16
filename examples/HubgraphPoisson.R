################################################
## Main Functions for Hub Graph for Poisson DAG
################################################


#################################################
### Tuning Parameter by Sun et al (2013 JMLR) ###
#################################################
## Helper Function
## Calculate the Kappa coefficient of Two Sets
agree.twosets <- function(set1, set2, left, p.tot){
  if(length(set1)+length(set2) ==0 || length(set1)+length(set2)==2*p.tot )
    kap <- -1 else{
      n11 <- length(intersect(set1, set2))
      n12 <- length(setdiff(set1,  intersect(set1, set2))) 
      n21 <- length(setdiff(set2,  intersect(set1, set2)))
      n22 <- length(intersect(setdiff(left, set1), setdiff(left, set2)))
      kap <- ( (n11+n22)/p.tot - ((n11+n12)*(n11+n21)+(n12+n22)*(n21+n22)) / (p.tot*p.tot) ) / 
        (1 -  ((n11+n12)*(n11+n21)+(n12+n22)*(n21+n22)) / (p.tot*p.tot) )
      
    }
  return(kap)
}



TuningPara <- function(layer, X,  num.split = 5, alpha = 0.1, beta_1 =1, beta_2 = 0){
  n <- dim(X)[1]
  p <- dim(X)[2]
  #num.grid <- 20
  #num.split <- 20
  n1 <- as.integer(n/2)
  n2 <- n - n1
  #alpha <- 0.1
  #grid.t <- 10^(seq(-3,2,by=0.05))
  s <- seq(0,60,by=1)
  grid.t <- 10^(-2+0.15*s)[-1]
  num.grid <- length(grid.t)
  
  done <- unlist(layer)
  kappa.mat <- matrix(0, num.split, num.grid)
  
  for (i.split in 1:num.split) {
    set.seed(i.split)
    for (i.grid in 1:num.grid) {
      
      data.perm <- X[sample(1:n, n, replace = FALSE), ]
      data1 <- data.perm[1:n1, ]
      data2 <- data.perm[-(1:n1), ]
      
      if(length(done) == 0){ ## A_0 for roots
        
        set1 <- layer_A0(layer, X = data1, thres = grid.t[i.grid], beta_1, beta_2 )
        set2 <- layer_A0(layer, X = data2, thres = grid.t[i.grid], beta_1, beta_2 )
        left <- seq(p)
        p.tot <- length(left)
        
      }else if( length(done) >=1 & length(done) < (p-1)){  ## layer A_1 to A_{T-1}
        
        set1 <- layer_At(layer, X = data1, thres = grid.t[i.grid], beta_1, beta_2 )
        set2 <- layer_At(layer, X = data2, thres = grid.t[i.grid], beta_1, beta_2 )
        ##出错原因在于set1和set2得到的是变量
        left <- seq(p)[-done]
        p.tot <- length(left)
        
      }
      
      kappa.mat[i.split, i.grid] <- agree.twosets(set1, set2, left, p.tot) 
      ###对于每一层而言计算kappa的p不同
    }
  }
  
  kappa.vec <- apply(kappa.mat, 2, mean)
  if(is.nan(kappa.vec) || is.infinite(kappa.vec)){
    kappa.vec <- kappa.vec[-which(is.nan(kappa.vec))]
    kappa.vec <- kappa.vec[-which(is.infinite(kappa.vec))]
  }
  
  i.opt.kappa <- grid.t[which( kappa.vec / max(kappa.vec) >= (1-alpha) )[1]]
  result <- list(alpha = alpha, kappa.values = kappa.vec, thres = i.opt.kappa)
  return(result)
}




## roots for the layer A_0
layer_A0 <- function(layer, X, thres, beta_1, beta_2){
  n <- dim(X)[1]
  p <- dim(X)[2]
  T <- matrix(0, n, p)
  #layer <- list()
  
  # roots layer A_0
  for (i in 1:p) {
    T[, i] <- (beta_1 + beta_2 * mean(X[,i]))^(-1) * X[,i]
  }
  ratio <- c()
  for (j in 1:p) {
    if(mean(T[,j]) == 0){ ## 避免分母为0造成ratio NaN
      ratio[j] <- (var(T[, j]) + 1e-6) / (mean(T[, j]) + 1e-6)
    }else{
      ratio[j] <- var(T[, j]) / mean(T[, j])
    }
  }
  
  ratio <- ratio / min(ratio)
  A0 <- which(abs(ratio - 1) <= thres)
  return(A0)
}



## nodes for other layer A_t
layer_At <-function(layer, X, thres, beta_1, beta_2){
  
  n1 <- dim(X)[1]
  p <- dim(X)[2]
  
  done <- unlist(layer)
  left <- seq(p)[-done]
  n.left <- length(left)
  para <- matrix(0, n.left, 1 + length(done)) # regression intercept
  #para_square <- matrix(0, n.left, 1 + length(done))
  for (k in 1:n.left) {
    para[k, ] <- glm(X[, left][, k] ~ X[, done], family = "poisson",control=list(maxit=100))$coefficients
    #para[k, ] <- zeroinfl(X[, left][, k] ~ X[, done], dist = "poisson")$coefficients$count
    #每行表示x_k对x_done的regress的参数值
    #para_square[k, ] <- glm(X[, left][, k]^2 ~ X[, done], family = "poisson",control=list(maxit=100))$coefficients
    
  }
  con.mean1 <- c()
  w  <- c()
  
  for (l in 1:n.left) {
    #mean.pre <- X[, done] * para[, 2:(1 + length(done))]
    con.mean1[l] <- mean( exp(para[l, 1] +  X[, done] * para[l, 2:(1 + length(done))] ))
    #con.var[l] <- mean( exp(para_square[l, 1] +  X[, done] * para_square[l, 2:(1 + length(done))])) - con.mean[l]^2
    w[l]  <- (beta_1 + beta_2 * con.mean1[l])^(-1)
  }
  Tt <- X[, left] * w
  
  con.mean <- c()
  con.var  <-  c()
  para.mean <- matrix(0, n.left, 1+length(done))
  para.var  <- matrix(0, n.left, 1+length(done))
  
  for (m in 1:n.left) {
    para.mean[m, ] <- glm(Tt[, m] ~ X[, done], family = "poisson",control=list(maxit=100))$coefficients 
    para.var[m, ]  <- glm(Tt[, m]^2 ~ X[, done], family = "poisson",control=list(maxit=100))$coefficients 
    #para.mean[m, ] <- zeroinfl(Tt[, m] ~ X[, done], dist = "poisson")$coefficients$count 
    #para.var[m, ]  <- zeroinfl(Tt[, m]^2 ~ X[, done], dist = "poisson")$coefficients$count 
    con.mean[m] <- mean( exp(para.mean[m, 1] +  X[, done] * para.mean[l, 2:(1 + length(done))] ))
    con.var[m]  <- mean( (exp(para.var[l, 1] +  X[, done] * para.var[l, 2:(1 + length(done))] ))^2 ) - con.mean[m]^2
  }
  
  
  
  ratio.t <- w * con.var / con.mean
  ratio.t <- ratio.t / min(ratio.t)
  At <- left[which(abs(ratio.t - 1) <= thres)]
  
  return(At) 
}




##########################
############  EDAG  ######
##########################

## Exponential DAG 
## First step: finding layer and topological ordering; Second Step: directed edges

EDAG <- function(X, beta_1, beta_2){
  ptm <- proc.time() 
  n <- dim(X)[1]
  p <- dim(X)[2]
  T <- matrix(0, n, p)
  done <- NULL
  layer <- list()
  ratio.t <- list()
  thres <- list()
  
  # Step1: layer and ordering
  # roots layer A_0
  for (i in 1:p) {
    T[, i] <- (beta_1 + beta_2 * mean(X[,i]))^(-1) * X[,i]
  }
  ratio <- c()
  for (j in 1:p) {
    if(mean(T[,j]) == 0){ ## 避免分母为0造成ratio NaN
      #ratio[j] <- (var(T[, j]) + 1e-6) / (mean(T[, j]) + 1e-6)
      ratio[j] <- var(T[, j]) / mean(T[, j])
    }else{
      ratio[j] <- var(T[, j]) / mean(T[, j])
    }
  }
  ratio.t[[1]] <- ratio
  ratio <- ratio / min(ratio)   ###scale
  thresh <- TuningPara(layer, X, num.split = 5, alpha = 0.1, beta_1 =1, beta_2 = 0)$thres
  layer[[1]] <- which(abs(ratio - 1) <= thresh)
  thres[[1]] <- thresh
  ##如果A_0 选不出 roots, 判断第二层就会报错,piosson regression fails,原因在于thres太小
  
  # layer A_t
  cont <- TRUE
  while( cont ){
    for (lay in 2:p) {
      if( length(unlist(layer)) < p-1){   
        
        print(lay)
        done <- unlist(layer)
        left <- seq(p)[-done]
        n.left <- length(left)
        para <- matrix(0, n.left, 1 + length(done)) # regression intercept
        # para_square <- matrix(0, n.left, 1 + length(done))
        for (k in 1:n.left) {
          para[k, ] <- glm(X[, left][, k] ~ X[, done], family = "poisson",control=list(maxit=100))$coefficients
          #para[k, ] <- zeroinfl(X[, left][, k] ~ X[, done], dist = "poisson")$coefficients$count
          #每行表示x_k对x_done的regress的参数值
          #para_square[k, ] <- glm(X[, left][, k]^2 ~ X[, done], family = "poisson",control=list(maxit=100))$coefficients
          
        }
        
        con.mean1 <- c()
        w  <- c()
        
        for (l in 1:n.left) {
          #mean.pre <- X[, done] * para[, 2:(1 + length(done))]
          con.mean1[l] <- mean( exp(para[l, 1] +  X[, done] * para[l, 2:(1 + length(done))] ))
          #con.var[l] <- mean( exp(para_square[l, 1] +  X[, done] * para_square[l, 2:(1 + length(done))])) - con.mean[l]^2
          w[l]  <- (beta_1 + beta_2 * con.mean1[l])^(-1)
        }
        Tt <- X[, left] * w
        #Tt.var <- Tt - apply(Tt, 2, mean)
        
        con.mean <- c()
        con.var  <-  c()
        para.mean <- matrix(0, n.left, 1+length(done))
        para.var  <- matrix(0, n.left, 1+length(done))
        
        for (m in 1:n.left) {
          para.mean[m, ] <- glm(Tt[, m] ~ X[, done], family = "poisson",control=list(maxit=100))$coefficients 
          para.var[m, ]  <- glm(Tt[, m]^2 ~ X[, done], family = "poisson",control=list(maxit=100))$coefficients 
          
          #para.mean[m, ] <- zeroinfl(Tt[, m] ~ X[, done], dist = "poisson")$coefficients$count 
          #para.var[m, ]  <- zeroinfl(Tt[, m]^2 ~ X[, done], dist = "poisson")$coefficients$count 
          con.mean[m] <- mean( exp(para.mean[m, 1] +  X[, done] * para.mean[l, 2:(1 + length(done))] ))
          con.var[m]  <- mean( (exp(para.var[l, 1] +  X[, done] * para.var[l, 2:(1 + length(done))] ))^2 ) - con.mean[m]^2
        }
        
        ratio.t[[lay]] <-  con.var / con.mean
        ratio <-  con.var / con.mean
        ratio.ts <- ratio / min(ratio)
        A <- TuningPara(layer, X, num.split = 5, alpha = 0.1, beta_1 =1, beta_2 = 0)
        thresh <- A$thres
        thres[[lay]] <- thresh
        layer[[lay]] <- left[which(abs(ratio.ts - 1) <= thresh)]
        
        
      }else if(length(unlist(layer)) == p-1){
        ## sink node is in A_T 
        
        done <- unlist(layer)
        left <- seq(p)[-done]
        n.left <- length(left)
        layer[[length(layer) + 1]] <- left
        
      }else (length(unlist(layer)) >= p)
      cont <- FALSE
    }
    
  }
  
  layer <- list(layer,ratio.t,thres)
  
  ##Step2:Sparse Estimation for Direction Edges
  rr <- rev(layer[[1]])
  T  <- length(layer[[1]])
  result <- matrix(0, p, p)
  for (ii in 1:(T-1)) {
    for (k in 1:length(rr[[ii]])) {
      now <- rr[[ii]][k]  ## each element in each layer
      
      for (j in (ii+1):T) {
        this <- rr[[j]]  ## each layer
        if(length(this) > 1){
          lassom <- glmnet::cv.glmnet(X[, this], X[, now], family = "poisson", nfolds = 5 )
          bfit <- coefficients(lassom)[-1]
          for (mm in 1:length(this)) {
            if(bfit[mm] != 0)
              result[this[mm], now] <- 1
          }
        }else{
          
          lmod <- glm(X[, now] ~ X[, this], family = "poisson", control=list(maxit=100))
          coeff <- as.vector(lmod$coefficients[2])
          if(coeff !=0 )
            result[this, now] <- 1
        }
      }
    }
  }
  #elapse <- proc.time() - ptm
  return(list(adj = result, OrderLayer = rev(rr)))
  
  
}


#EDAG(X, beta_1 = 1, beta_2 = 0)


###################################################################################### 
### L1 penalized likelihood regression for GLMs                         ##############
### for the generalized ODS algorithm, Park and Raskutti (2018 JMLR)    ##############
######################################################################################

### ODS(GLMLasso)

ODSGLMLasso <- function(X, beta_1, beta_2){
  ptm<-proc.time()
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  #Step1: NormalizedGraph
  NormalizedGraph <- matrix(0, p, p)
  for (i in 1:p) {
    left <- seq(p)[-i]
    lassom <- glmnet::cv.glmnet(X[, left], X[, i], family = "poisson", alpha = 1, nfolds = 5)
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
      neigh <- union(which( norgraph[k, ] != 0), which( norgraph[, k] != 0))     # k' neighborhood 是否要考虑第k列不为0，两者的union, YES!
      C[[k]] <- intersect(neigh, ODS.order)   # given parents
      
      # calculate conditinal mean/variance  
      if(length(C[[k]]) ==0 ){
        con.mean <- mean(X[, k])
        con.var  <- var(X[, k])
        
        w <- (beta_1 + beta_2 * con.mean)^(-1)
        kk <- which(left.order == k)
        S[kk] <- w^2 * con.var - w * con.mean
        
      }else{
        para.mean <- glm( X[, k] ~ X[, C[[k]]], family = "poisson", control = list(maxit=100))$coefficients 
        para.var  <- glm( X[, k]^2 ~ X[, C[[k]]], family = "poisson", control = list(maxit=100))$coefficients
        
        con.mean <- mean( exp(para.mean[1] + X[, C[[k]]] * para.mean[2:(1+length(C[[k]]))] )) 
        con.var  <- mean( exp(para.var[1] + X[, C[[k]]] * para.var[2:(1+length(C[[k]]))] )) 
        
        
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
      # variable selection
      lassom <- glmnet::cv.glmnet(X[, this], X[, now], family = "poisson", alpha = 1, nfolds = 5 )
      bfit <- coefficients(lassom)[-1]
      
      for (jj in 1:length(this)){
        if( bfit[jj] != 0 )
          result[this[jj],now] <- 1
      }
    } else {
      # deal with the last two nodes
      lmod <- glm(X[, now] ~ X[, this], family = "poisson", control = list(maxit=100))
      coeff <- as.vector(lmod$coefficients[2])
      if (coeff != 0) {
        result[this,now] <- 1
      }
    }
  }
  
  #result[lower.tri(result)] <- 0   # aviod cycles ??????如果同时出现i->j 和j->i应该怎么选择？
  #elapse <- proc.time()-ptm
  #return(list(adj=result,ordering=rev(rr), elapse))
  return(list(adj=result,ordering=rev(rr)))
}


#ODSGLMLasso(X, beta_1 = 1, beta_2 = 0)


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


do_one <- function(n, p){
  result <- matrix(0, 6, 2)
  # Data generation
  dat <- hubgraph(n, p)
  truth <- dat$truth
  X <- dat$X
 
  # EDAG
  ptm <- proc.time()
  sresult <- EDAG(X, beta_1 = 1, beta_2 = 0)
  sx <- sresult$adj
  result[1, 1] <- (proc.time() - ptm)[1]
  result[2, 1] <- hammingDistance(sx, truth)/(p*(p-1))  ##normalized hammdis
  #result[3, 1] <- Kendall(sapply(1:p, function(i){which(TO == i)}), sapply(1:p, function(i){which(unlist(sresult$OrderLayer) == i)}))$tau[1]
  result[3, 1] <- ifelse(sum(truth), sum(truth*sx)/sum(truth), 0) #recall
  result[4, 1] <- ifelse(sum(truth), sum(truth*sx)/sum(sx), 0) # precision
  result[5, 1] <- 2*result[3, 1]*result[4,1]/(result[3, 1]+result[4,1]+1e-6)  # F1 score
  result[6, 1] <- ifelse(sum(sx), 1 - sum(truth*sx)/sum(sx), 0) #FDR
  #sxtd <- sx
  
  
  # ODSGLMLasso
  ptm <- proc.time()
  sresult <- ODSGLMLasso(X, beta_1 = 1, beta_2 = 0)
  sx <- sresult$adj
  result[1, 2] <- (proc.time() - ptm)[1]
  result[2, 2] <- hammingDistance(sx, truth)/(p*(p-1))     ##normalized hammdis
  #result[3, 2] <- Kendall(sapply(1:p, function(i){which(TO == i)}), sapply(1:p, function(i){which(sresult$ordering == i)}))$tau[1]
  result[3, 2] <- ifelse(sum(truth), sum(truth*sx)/sum(truth), 0) #recall
  result[4, 2] <- ifelse(sum(truth), sum(truth*sx)/sum(sx), 0) # precision
  result[5, 2] <- 2*result[3, 2]*result[4,2]/(result[3, 2]+result[4,2]+1e-6)  # F1 score  # for NaN
  result[6, 2] <- ifelse(sum(sx), 1 - sum(truth*sx)/sum(sx), 0) #FDR
  
  return(result)
}




####### repeat

sim_run <- function(n, p, N, pa){
  # p: number of nodes
  # n: sample size
  # N: simulation times
  if(pa == T){
    cl <- makeCluster(detectCores()-1)
    registerDoSNOW(cl)
    print(detectCores())
    results <- foreach::foreach(i = 1:N, .export = ls(globalenv()), .packages = c("glmnet", "matrixcalc")) %dopar% {
      do_one(n, p)
    }
    stopCluster(cl)
  }else{
    pb <- txtProgressBar(min = 0, max = N, style = 3)
    result <- matrix(0, 6, 2)
    for (i in 1:N) {
      results <- foreach::foreach(i = 1:N, .export = ls(globalenv()),.packages = c("glmnet", "matrixcalc")) %do% {
        do_one(n, p)
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  result <- Reduce('+', results) #将多个数据框按照同一列merge
  #result <- rbind(result, apply(FUN = sd, MARGIN = 2, X = Reduce('rbind', lapply(results, function(x)x[1, ]))))
  #normalized hamming distance
  # 结果变为7*2矩阵
  
  return(result/N)
}




simulation_hubPoi <- function(n, p, N){
  #set.seed(1)
  #results <-
  # return hub graph
  cat('n=',n, "p=",p, "N=",N)
  result <- sim_run(n, p, N, T)
  sink(paste("p=",p,"n=",n,"_hubPoi.txt"))
  cat('Setting', 'p=',p, "n=",n,"N=",N, "\n")
  cat('\tEDAG\tODS\n')
  cat('Avg. time\t', round(result[1,],4),'\n')
  cat('Avg. Dist\t', round(result[2,],4),'\n')
  cat('Avg. Reca\t', round(result[3,],4),'\n')
  cat('Avg. Prec\t', round(result[4,],4),'\n')
  cat('Avg. fsco\t', round(result[5,],4),'\n')
  cat('Avg. FDR \t', round(result[6,],4),'\n')
  #cat('SD . Dist\t', round(result[7,],4),'\n')
  sink()
    
}

cat("done.\n")

#####################################
## generate data for hub graph (d=1)  
#####################################

hubgraph <- function(n, p){
  ordering <- seq(p)
  truth <- matrix(0,p,p)
  truth[1, 2:p] <- 1
  
  theta <- runif(p, 1, 3)
  beta  <- matrix(runif(p*p, -1, -0.5), p ,p)
  X <- matrix(0, n, p)
  X[, 1] <- rpois(n, theta[1])
  for (j in 2:p) {
    g <- exp(theta[j] + beta[j, 1] * X[, 1])
    X[, j] <- rpois(n, lambda = g)
  }
  
  cont <- TRUE
  while( cont ){
    if(sum(apply(X,2,function(x) all(x)==0)==TRUE)/p >= 0.5){ #如果产生数据X有一列为0, regenerate
      theta <- runif(p, 1, 3)
      beta  <- matrix(runif(p*p, -1, -0.5), p ,p)
      X <- matrix(0, n, p)
      X[, 1] <- rpois(n, theta[1])
      for (j in 2:p) {
        g <- exp(theta[j] + beta[j, 1] * X[, 1])
        X[, j] <- rpois(n, lambda = g)
      }
    }else{
      cont <- FALSE
    }    #(sum(apply(X,2,function(x) all(x)==0)==TRUE)/p < 0.5)
  }
  
  return(list(X = X, truth = truth, TO = ordering))
}




## performance for one experiment
do_one(n = 200, p = 100)
sim_run(n = 200, p = 5, N = 50, pa = T)

simulation_hubPoi(n = 200, p = 5, N = 50)
simulation_hubPoi(n = 300, p = 5, N = 50)
simulation_hubPoi(n = 500, p = 5, N = 50)

simulation_hubPoi(n = 200, p = 20, N = 50)
simulation_hubPoi(n = 300, p = 20, N = 50)
simulation_hubPoi(n = 500, p = 20, N = 50)

simulation_hubPoi(n = 200, p = 100, N = 50)
simulation_hubPoi(n = 300, p = 100, N = 50)
simulation_hubPoi(n = 500, p = 100, N = 50)
