###########################################################
### To calculate the ratio for every node in each layer ###
###########################################################



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
  
  uR.ratio <- apply(uR.T, 2, var) / colMeans(uR.T)
  
  return(uR.ratio)
}



##################################
###  roots for the layer A_0   ###
##################################

layer_A0 <- function(A0.X, A0.beta1.pois, A0.beta2.pois, A0.beta1.bino, A0.beta2.bino, A0.bino.N, Distri = c("pois","mix")){
  
  A0.n <- dim(A0.X)[1]
  A0.p <- dim(A0.X)[2]
  
  if(Distri == "pois"){
    ## calcalating ratio for poission distribution
    
    A0.TF <- matrix(0, A0.n, A0.p)
    A0.TF <- t( (A0.beta1.pois + A0.beta2.pois *colMeans(A0.X))^(-1) ) %x% rep(1, A0.n) * A0.X
    A0.ratio <-  apply(A0.TF, 2, var) / colMeans(A0.TF)
    
  }else{
    ## calcalating ratio for mixed distributions
    
    A0.IndexPois <- CheckIndex(A0.X, A0.bino.N)$IndexPois
    A0.IndexBino <- CheckIndex(A0.X, A0.bino.N)$IndexBino
    
    A0.XPois <- as.matrix(A0.X[, A0.IndexPois])
    A0.XBino <- as.matrix(A0.X[, A0.IndexBino])
    
    A0.ratio <- rep(0, A0.p)
    A0.ratio[A0.IndexPois] <- unconRatio(A0.XPois, uR.beta_1=A0.beta1.pois, uR.beta_2=A0.beta2.pois)
    A0.ratio[A0.IndexBino] <- unconRatio(A0.XBino, uR.beta_1=A0.beta1.bino, uR.beta_2=A0.beta2.bino)
  }
  
  A0.final.ratio <- A0.ratio / min(A0.ratio)
  return(A0.final.ratio)
}



#####################################################
## calculate ratio for the nodes in the layer A_t ###
#####################################################



layer_At <-function(At.layer, At.X, At.XBino, At.beta1.pois, At.beta2.pois, At.beta1.bino, At.beta2.bino, At.bino.N, Distri = c("pois","mix"), At.IndexPois, At.IndexBino){
  
  At.n <- dim(At.X)[1]
  At.p <- dim(At.X)[2]
  At.Sdone <- unlist(At.layer)
  At.left <- seq(At.p)[-At.Sdone]
  At.n.left <- length(At.left)

  At.left.index <- sapply(1:At.n.left, function(i) To_dist(At.left[i], At.IndexBino, At.IndexPois))
  At.ratio <- sapply(1:At.n.left, function(i) At_single_ratio(At.X[, At.left][, i], At.X[, At.Sdone],
                    At.XBino[, At.left][, i], At.XBino[, At.Sdone],
                    At.beta1.pois, At.beta2.pois, At.beta1.bino, At.beta2.bino, At.bino.N,
                    At.IndexBino, At.IndexPois, Distri, At.left.index, i))
  
  At.final.ratio <- At.ratio / min(At.ratio)

  return(At.final.ratio) 
}





At_single_ratio <- function(AR.X, AR.St, AR.XBino, AR.XBino.St, AR.beta1.pois, AR.beta2.pois, AR.beta1.bino, AR.beta2.bino, AR.bino.N,  
                            AR.IndexBino, AR.IndexPois, AR.dist = c("pois", "mix"), AR.left.index, AR.X.index){
  
  AR.X.node = AR.left.index[AR.X.index]
  dim_x <- dim(as.matrix(AR.St))[2]
  
  if(AR.dist == "pois"){
    
      
      glmfit <- glm(AR.X ~ AR.St, family = "poisson", control=list(maxit=100))
      con.mean1Pois <- mean(predict(glmfit))
      Tt2Pois <- (AR.X)^2
      AR.con.mean <- mean(AR.X)    
      
      fitT2Pois <- glm(Tt2Pois ~ AR.St, family = "poisson", control=list(maxit=10000)) 
      sec.moment <- mean(exp(predict(fitT2Pois)))
      AR.con.var <- sec.moment - con.mean1Pois^2
    
  }else if(AR.dist == "mix"){
    
    if(AR.X.node == "bino"){

      glmfit <- glm(AR.XBino ~ AR.XBino.St, family = "binomial")
      glmfit.coef <- as.matrix(coef(glmfit)[-1])
        
      
      if(any(is.na(glmfit.coef))){
        na.index <- which(is.na(glmfit.coef) == TRUE)
        glmfit.coef[na.index] <- 0
      }
        
      ## Below is only for binomial distribution for computational simplicity and 
      # to avoid some extreme case leading to NA in computing ratio
      if(any(glmfit.coef <= 0.5)){
          bi_index <- which(glmfit.coef <= 0.5)
          glmfit.coef[bi_index] <- 0
      }
      
      if(any(glmfit.coef >= 5)){
        bi_index <- which(glmfit.coef >= 5)
        glmfit.coef[bi_index] <- 0
      }
   
      
      glm_St <- cbind(AR.XBino, AR.XBino.St)
      glmfit.coef <- rbind(1, glmfit.coef)
      exp.term <- exp(glm_St %*% glmfit.coef)
      con.prob <- mean(exp.term^(AR.XBino) * (factorial(1)/(factorial(AR.XBino)*factorial(1-AR.XBino))) / (1 + exp.term)^(AR.XBino))
      con.mean.X <-  con.prob * sum(AR.X)
      w.Tt <- (AR.beta1.bino + AR.beta2.bino * con.mean.X)^(-1)
      AR.Tt <- AR.X * w.Tt
      AR.con.mean <- con.prob * w.Tt * sum(AR.Tt)
      
      sec.moment <- con.prob * w.Tt * sum(AR.Tt^2)
      AR.con.var <- sec.moment - AR.con.mean^2
      
    }else if(AR.X.node == "pois"){
      

        glmfit <- glm(AR.X ~ AR.St, family = "poisson", control=list(maxit=100))
        con.mean1Pois <- mean(exp(predict(glmfit)))
        Tt2Pois <- (AR.X)^2
        AR.con.mean <- mean(AR.X)
        
        fitT2Pois <- glm(Tt2Pois ~ AR.St, family = "poisson", control=list(maxit=100))
        sec.moment <- mean(exp(predict(fitT2Pois)))
        AR.con.var <- sec.moment - con.mean1Pois^2
        
        
    }
  }
  
  AR.single.ratio <- AR.con.var / AR.con.mean
  return(AR.single.ratio) 
}





