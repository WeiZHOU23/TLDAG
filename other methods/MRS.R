########################################
### MRS algorithm in Park (2019 JMLR)
########################################



MRS <- function(X){
  ptm<-proc.time()
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  ####################################
  #Step1: ordering
  ####################################
  
  MRS.order <- c()
  coeff_list <- list()  ## record the non-zero coefficients in the ordering ooder for parents recovery
  
  ## first element of the ordering
  
  score1.nu <- apply(X^2, 2, mean)
  score1.de <- apply(X, 2, mean) + (apply(X, 2, mean))^2
  S1 <- score1.nu / score1.de
  MRS.order[1] <- seq(p)[which.min(S1)]
  
  # m-th element of the ordering m=2 to p-1
  # calculate scores
  for (m in 2:p) {
    MRS.order <- MRS.order  #### renew the ordering
    left.order <- seq(p)[-MRS.order]
    S <- c()
    coeff_list[[m]] <- list()   ### each ordering node corresponds to a list 
    
    for (k1 in 1:length(left.order) ) { 
      
      k <- left.order[k1]    # k  the original variable index
      scoreM.nu <- mean(X[, k]^2)
      
      if(length(MRS.order) == 1){
        
        glmPois_n1 <- glm( X[, k] ~ X[, MRS.order], family = "poisson", control = list(maxit=100))
        scoreM.de1 <- mean(exp(predict(glmPois_n1)))
        scoreM.de2 <- mean(exp(2 * predict(glmPois_n1)))
        scoreM.de <- scoreM.de1 + scoreM.de2
        
        S[k1] <- scoreM.nu / scoreM.de
        
        ## record the non-zero coefficients when applying the L1-GLM
        nonzero_coef_1 <- as.matrix(coefficients(glmPois_n1))[-1]
        nonzero_length_1 <- length(nonzero_coef_1)
        coeff_list[[m]][[k1]] <- c(nonzero_coef_1, rep(0, p - nonzero_length_1))
        
      }else{
        
        set.seed(k)
        glmPois <- glmnet::cv.glmnet(X[, MRS.order], X[, k], family = "poisson", alpha = 1, nfolds = 5)
        newX <- cbind(rep(1, n),  X[, MRS.order])
        scoreM.de1 <- mean(exp(newX %*% as.matrix(coefficients(glmPois)))) 
        scoreM.de2 <- mean(exp(2 * newX %*% as.matrix(coefficients(glmPois)))) 
        scoreM.de <- scoreM.de1 + scoreM.de2
        
        S[k1] <- scoreM.nu / scoreM.de
        
        
        ## record the non-zero coefficients when applying the L1-GLM
        nonzero_coef <- as.matrix(coefficients(glmPois))[-1]
        nonzero_length <- length(nonzero_coef)
        coeff_list[[m]][[k1]] <- c(nonzero_coef, rep(0, p - nonzero_length))
        
      }
      
      MRS.order[m] <- left.order[which.min(S)]
      
    }
    
  }
  
  
  ########################################
  ##Step2: directed edges recovery
  ########################################
  
  ordering <- MRS.order
  #rr <- rev(ordering)
  coeff.matrix <- matrix(0, p, p)
  result <- matrix(0, p, p)
  
  for (m in 2:p) {
    node <- ordering[m]
    order.ident <- ordering[1:(m-1)]
    left.oder <- seq(1:p)[-order.ident]
    k <- which(left.oder == node)
    coeff_value <- coeff_list[[m]][[k]]
    pa <- ordering[which(coeff_value > 1e-5)]
    
    result[pa, node] <- 1
  }
  
  
  return(list(DAG = result, ordering = ordering))
}


