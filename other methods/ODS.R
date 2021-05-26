####################################
### ODS(GLMLasso) for mixed graph ##
####################################

ODSGLMLasso <- function(X, beta1.pois, beta2.pois, beta1.bino, beta2.bino, bino.N, distri = c("pois","mix")){
  ptm<-proc.time()
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  if(distri == "mix"){
   
    IndexBino <- CheckIndex(X, bino.N)$IndexBino
    IndexPois <- CheckIndex(X, bino.N)$IndexPois
    XBino <- as.matrix(X[, IndexBino])
    XPois <- as.matrix(X[, IndexPois])
    
    if(length(IndexBino) != 0){
      X.all <- matrix(0, n*bino.N, p)
      X.all[, IndexBino] <- sapply(1:length(IndexBino), function(i) ConvBinomial(X[, IndexBino[i]], bino.N))
      
      if(length(IndexPois) != 0){
        
          X.all[, IndexPois] <- sapply(1:length(IndexPois), function(i) PoissonNew2(X[, IndexPois[i]], bino.N))
          
       }
    }else{
      X.all <- XPois
    } 
  }
  
  ##########################
  #Step1: NormalizedGraph
  ##########################
  
  NormalizedGraph <- matrix(0, p, p)
  for (i in 1:p) {
    left <- seq(p)[-i]
    
    
    if(distri == "mix"){
      
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
        
        
        lassomBino <- glmnet::cv.glmnet(X.all[, left], X.all[, i], family = "binomial", foldid = foldid, weights = weights, lambda=exp(seq(log(0.001), log(5), length.out=100)) )
        bfit <- coefficients(lassomBino)[-1]
        
      }else{

        lassomPois <- glmnet::cv.glmnet(X[, left], X[, i], family = "poisson", alpha = 1, nfolds = 5)
        bfit <- coefficients(lassomPois)[-1]
      }
      
    }else{

      set.seed(4)
      lassomPois <- glmnet::cv.glmnet(X[, left], X[, i], family = "poisson", alpha = 1, nfolds = 5)
      bfit <- coefficients(lassomPois)[-1]
      
    }
    

    for (j in 1:length(bfit)) {
      if(bfit[j] != 0)
        NormalizedGraph[i, left[j]] <- 1
    }
  }
  
  ###########################
  #####  Step2: ordering  ###
  ###########################
  
  ODS.order <- c()
  
  if(distri == "mix"){
    
    ratioBino <- RatioODS(XBino, beta1.bino, beta2.bino)
    ratioPois <- RatioODS(XPois, beta1.pois, beta2.pois)
    S1 <- rep(0,p)
    S1[IndexBino] <- ratioBino
    S1[IndexPois] <- ratioPois
    
  }else{
    
    S1 <- RatioODS(X, beta1.pois, beta2.pois)
  }
  ODS.order[1] <- seq(p)[which.min(S1)]
  

  
  # m-th element of the ordering m=2 to p-1
  # find condidate parents
  left.order <- seq(p)[-ODS.order]
  norgraph <- NormalizedGraph
  C <- list()
  c0 <- 0.005
  
  # calculate scores
  for (m in 2:(p-1)) {
    ODS.order <- ODS.order  #### renew the ordering
    left.order <- seq(p)[-ODS.order]
    S <- c()
    
    for (k1 in 1:length(left.order) ) { # k original variable
      k <- left.order[k1]
      neigh <- union(which( norgraph[k, ] != 0), which( norgraph[, k] != 0))   #k'neighborhood, union, YES!
      C[[k]] <- intersect(neigh, ODS.order)   # given parents
      kk <- which(left.order == k)
      
      # calculate conditinal mean/variance  
      if(length(C[[k]]) ==0 ){  
        
        if(distri == "mix"){
          
          if(k %in% IndexBino){
            S[kk] <- RatioODS(as.matrix(X[, k]), beta1.bino, beta2.bino)
          }else{
            S[kk] <- RatioODS(as.matrix(X[, k]), beta1.pois, beta2.pois)
          }
        }else{
          
          S[kk] <- RatioODS(as.matrix(X[, k]), beta1.pois, beta2.pois)
        }
        
        
        
      }else{
        
        if(distri == "mix"){
          
          if(k %in% IndexBino){
            
            glmfit <- glm(X.all[, k] ~ X.all[, C[[k]]], family = "binomial", control=list(maxit=100))
            glmfit.coef <- as.matrix(glmfit$coefficients)
            
            if(any(is.na(glmfit.coef))){
              na.index <- which(is.na(glmfit.coef) == TRUE)
              glmfit.coef[na.index] <- 0
            }
            
            ## add details for binomial coefficients
            if(any(glmfit.coef <= 0.5)){
              bi_index <- which(glmfit.coef <= 0.5)
              glmfit.coef[bi_index] <- 0
            }
            
            if(any(glmfit.coef >= 5)){
              bi_index <- which(glmfit.coef <= 5)
              glmfit.coef[bi_index] <- 0
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
          
        }else{
          
          glmPois <- glm( X[, k] ~ X[, C[[k]]], family = "poisson", control = list(maxit=100)) 
          glm2Pois <- glm( X[, k]^2 ~ X[, C[[k]]], family = "poisson", control = list(maxit=100))
          con.mean <- mean(exp(predict(glmPois)))
          sec.moment <- mean(exp(predict(glm2Pois)))
          con.var <- sec.moment - con.mean^2
          
          w <- (beta1.pois + beta2.pois * con.mean)^(-1)
          nx <- n*length(C[[k]])
          nC <- sum( nx * ifelse(nx >= c0*n, 1, 0) )   
          kk <- which(left.order == k)
          S[kk] <- sum( (w^2 * con.var - w * con.mean) * nx / nC ) 
        }

      } 
    }
    
    ODS.order[m] <- left.order[which.min(S)]
  }
  
  ODS.order[p] <- seq(p)[-ODS.order]
  
  
  #################################
  ####  Step3: directed edges  ####
  #################################
  
  ordering <- ODS.order
  rr <- rev(ordering)
  result <- matrix(0, p, p)
  
  for (ii in 1:(p-1)){
    now <- rr[ii]
    this <- sort(rr[(ii+1):p])
    if (length(this) > 1){
      
      if(distri == "mix"){
        
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
          
          lassomBino <- glmnet::cv.glmnet(X.all[, left], X.all[, now], family = "binomial", foldid = foldid, weights = weights, lambda=exp(seq(log(0.001), log(5), length.out=100)) )
          bfit <- coefficients(lassomBino)[-1]
          
        }else{
          lassomPois <- glmnet::cv.glmnet(X[, left], X[, now], family = "poisson", alpha = 1, nfolds = 5)
          bfit <- coefficients(lassomPois)[-1]
        }
        
      }else{
        
        lassomPois <- glmnet::cv.glmnet(X[, this], X[, now], family = "poisson", alpha = 1, nfolds = 5)
        bfit <- coefficients(lassomPois)[-1]
      }

      
      for (jj in 1:length(this)){
        if( bfit[jj] != 0 )
          result[this[jj],now] <- 1
      }
      
    }else {
      # deal with the last two nodes
      if(distri == "mix"){
        
        if(now %in% IndexBino){
          lmod <- glm(X.all[, now] ~ X.all[, this], family = "binomial", control = list(maxit=100))
          coeff <- as.vector(lmod$coefficients[2])
        }else{
          lmod <- glm(X[, now] ~ X[, this], family = "poisson", control = list(maxit=100))
          coeff <- as.vector(lmod$coefficients[2])
        }
        
      }else{
        
        lmod <- glm(X[, now] ~ X[, this], family = "poisson", control = list(maxit=100))
        coeff <- as.vector(lmod$coefficients[2])
      }
      
      if (coeff != 0) 
        result[this,now] <- 1
    }
  }
  
  return(list(DAG = result, ordering = rev(rr)))
}




RatioODS <- function(X, beta_1, beta_2){
  con.mean1 <- apply(X, 2, mean)
  con.var1  <- apply(X, 2, var)
  w1 <- (beta_1 + beta_2 * con.mean1)^(-1)
  S1 <- w1^2 * con.var1 - w1 * con.mean1
  
  return(S1)
}


