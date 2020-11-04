###############################################
### estimate the order by DirectLINGAM
### Copyright in Nicola Gnecco (2019 arXiv)
###############################################

sourceCpp('direct_lingam_funcs.cpp')


direct_lingam_search <- function(dat){
  
  if (typeof(dat) != "double"){
    return(NA)
  }
  
  out <- tryCatch({
    #order <- .Call('direct_lingam_funcs_', dat)
    order <- direct_lingam_c(dat)         #######
    return(order)
  },
  error = function(e){
    order <- NA
    return(order)
  })
  
  return(out)
}




###########################################################################
### estimate the directed edges by least squares with the estimated order
###########################################################################

DirectLINGAM <- function(dat){
  
  p <- ncol(dat)
  K <- direct_lingam_search(dat)   ## estimated order
  DL.result <- matrix(0, p, p)
  
  for (i in 2:p) {
    for (j in 1:(i-1)) {
      
      y.index <- K[i]
      x.index <- K[c(1:(i-1))]
      x <- dat[, x.index]
      y <- dat[, y.index]
      mod <- lm(y~x)
      est.coef <- as.vector(mod$coefficients)[-1]  # remove the constant
      
      direct.index <- which(abs(est.coef) >= 0.05)
      child <- x.index[direct.index]
      
      DL.result[child, y.index] <- 1
      
    }
  }
  
  return(DL.result)  
}

