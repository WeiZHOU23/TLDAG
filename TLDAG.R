#########################################################
##### Main function for the proposed method: TLDAG  #####
#########################################################

TLDAG <- function(TD.X, TD.grid.t, TD.beta1.pois, TD.beta2.pois, TD.beta1.bino, TD.beta2.bino, TD.bino.N, TD.Distri = c("pois","mix")){
  
  ptm <- proc.time() 
  TD.n <- dim(TD.X)[1]
  TD.p <- dim(TD.X)[2]
  #TD.TF <- matrix(0, TD.n, TD.p)
  TD.done <- NULL
  TD.layer <- list()
  TD.ratio.t <- list()
  TD.thres <- list()
  TD.final.layer <- list()
  TD.time.list <- list()
  
  if(TD.Distri == "mix"){
    
    TD.IndexPois <- CheckIndex(TD.X, TD.bino.N)$IndexPois
    TD.IndexBino <- CheckIndex(TD.X, TD.bino.N)$IndexBino
    
    if(length(TD.IndexBino) != 0){
      TD.X.Binary <- matrix(0, TD.n*TD.bino.N, TD.p)
      TD.X.Binary[, TD.IndexBino] <- sapply(1:length(TD.IndexBino), function(i) ConvBinomial(TD.X[, TD.IndexBino[i]], TD.bino.N))
      
      
      if(length(TD.IndexPois) != 0){
    
          TD.X.Binary[, TD.IndexPois] <- sapply(1:length(TD.IndexPois), function(i) PoissonNew2(TD.X[, TD.IndexPois[i]], TD.bino.N))
      }
    }else{
      
      TD.X.Binary <- TD.X
    }
    
  }else{
    
    TD.IndexPois <- 0
    TD.IndexBino <- 0
    TD.X.Binary <- TD.X
  }


  ##################################################################################  
  ### Step1: layer and ordering
  ### roots layer A_0
  ############## Node identification ##################
  TD.A0.thresh <- TuningPara(TD.layer, TD.X, TD.grid.t, 0.1, TD.beta1.pois, TD.beta2.pois, TD.beta1.bino, TD.beta2.bino, TD.bino.N, TD.Distri, TD.IndexPois, TD.IndexBino)$thres
  TD.ratio.final<-layer_A0(TD.X, TD.beta1.pois, TD.beta2.pois, TD.beta1.bino, TD.beta2.bino, TD.bino.N, TD.Distri)
  TD.layer[[1]] <- which(abs(TD.ratio.final - 1) <= TD.A0.thresh)
  TD.thres[[1]] <- TD.A0.thresh
  TD.ratio.t[[1]]<-TD.ratio.final
  #######################################################################################
  # identify nodes in layer A_t
  if( length(unlist(TD.layer)) == TD.p){
    TDS.result <- matrix(0, TD.p, TD.p)
    TDS.layernodes <- unlist(TD.layer)
    
  }else{
    
    cont <- TRUE
    while( cont ){
      for (lay in 2:TD.p) {
        if( length(unlist(TD.layer)) < TD.p-1){   
          
          TD.left <- seq(TD.p)[-unlist(TD.layer)]
          #print(lay)
          ####################################################################
          TD.AT.thresh <- TuningPara(TD.layer, TD.X, TD.grid.t, 0.1, TD.beta1.pois, TD.beta2.pois, TD.beta1.bino, TD.beta2.bino, TD.bino.N, TD.Distri, TD.IndexPois, TD.IndexBino)$thres
          TD.thres[[lay]] <- TD.AT.thresh
          TD.AT.ratio.final<-layer_At(TD.layer, TD.X, TD.X.Binary, TD.beta1.pois, TD.beta2.pois, TD.beta1.bino, TD.beta2.bino, TD.bino.N, TD.Distri, TD.IndexPois, TD.IndexBino)
          TD.layer[[lay]] <- TD.left[which(abs(TD.AT.ratio.final - 1) <= TD.AT.thresh)]
          TD.ratio.t[[lay]]<-TD.AT.ratio.final
          #########################################################################
          
        }else if(length(unlist(TD.layer)) == TD.p-1){
          ## sink node is in A_T 
          TD.left <- seq(TD.p)[-unlist(TD.layer)]
          TD.n.left <- length(TD.left)
          TD.layer[[length(TD.layer) + 1]] <- TD.left
          
        }else (length(unlist(TD.layer)) >= TD.p)
        cont <- FALSE
        #break
      }
    }
    
    TD.final.layer <- list(TD.layer, TD.ratio.t, TD.thres)  
    
 
   
     #####################################################
     ### Step2:Sparse Estimation for Direction Edges   ###
     #####################################################

     TDS.layernodes <- rev(TD.final.layer[[1]]) ### Selected layer structure with nodes in each layer
     TDS.layernum  <- length(TD.final.layer[[1]]) ### number of layer in the DAG
     TDS.result <- matrix(0, TD.p, TD.p)


    node_parents <- list()

    for (ii in 1:(TDS.layernum-1)) {

      TDS.node.num <- length(TDS.layernodes[[ii]])

      node_parents[[ii]] <- sapply(1:length(TDS.layernodes[[ii]]), function(k) Parall_parents(TDS.layernodes[[ii]][k], potentialParents(TDS.layernodes, ii), pa.Distri = TD.Distri,
                                                                                              pa.IndexBino = TD.IndexBino, pa.X.Binary = TD.X.Binary, pa.X = TD.X, pa.P = TD.p) )

    }


    TDS.root <- TDS.layernodes[[TDS.layernum]]
    TDS.have.pa <- seq(1:TD.p)[-TDS.root]

    matrix_node_parents <- cbind_parents(node_parents, TD.p, TDS.have.pa)
    TDS.pa.index <- matrix_node_parents[1, ]
    TDS.result[, TDS.pa.index] <- sapply(1:length(TDS.have.pa), function(i) Trans_parents(matrix_node_parents[, i], TD.p))


   }
  # elapse2 <- proc.time() - ptm
  return(list(DAG = TDS.result, OrderLayer = rev(TDS.layernodes)))
}




## find the parents of one node
Parall_parents <- function(pa.now, pa.this, pa.Distri = c("mix","pois"), pa.IndexBino, pa.X.Binary, pa.X, pa.P){
  
  pa.result <- rep(0, pa.P+1)
  pa.result[1] <- pa.now
  # the first: now node, the last p variables: non-zero indicates parents nodes
  
  if(length(pa.this) > 1){
    
    if(pa.Distri == "mix"){
      
      if(pa.now %in% pa.IndexBino){
        pa.nfold <- 5
        pa.fold0 <- sample.int(sum(pa.X.Binary[, pa.now]==0)) %% pa.nfold
        pa.fold1 <- sample.int(sum(pa.X.Binary[, pa.now]==1)) %% pa.nfold
        pa.foldid <- numeric(length(pa.X.Binary[, pa.now]))
        pa.foldid[pa.X.Binary[, pa.now]==0] <- pa.fold0
        pa.foldid[pa.X.Binary[, pa.now]==1] <- pa.fold1
        pa.foldid <- pa.foldid + 1
        pa.fraction <- table(pa.X.Binary[, pa.now]) / length(pa.X.Binary[, pa.now])
        pa.weights <- 1 - pa.fraction[as.character(pa.X.Binary[, pa.now])]
        
        pa.lassom <- glmnet::cv.glmnet(pa.X.Binary[, pa.this], pa.X.Binary[, pa.now], family = "binomial", foldid = pa.foldid, weights = pa.weights, lambda=exp(seq(log(0.001), log(5), length.out=100)) )
        pa.bfit <- coefficients(pa.lassom)[-1]
      }else{
        pa.lassomPois <- glmnet::cv.glmnet(pa.X[, pa.this], pa.X[, pa.now], family = "poisson", nfolds = 5 )
        pa.bfit <- coefficients(pa.lassomPois)[-1]
      }
      
    }else{
      
      pa.lassom <- glmnet::cv.glmnet(pa.X[, pa.this], pa.X[, pa.now], family = "poisson", nfolds = 5)
      pa.bfit <- coefficients(pa.lassom)[-1]
    }
    
    
    for (mm in 1:length(pa.this)) {
      if(pa.bfit[mm] != 0){

        pa.result[pa.this[mm]+1] <- 1
      }
      
    }
    
  }else{
    
    if(pa.Distri == "mix"){
      
      if(pa.now %in% pa.IndexBino){
        pa.lmod <- glm(pa.X.Binary[, pa.now] ~ pa.X.Binary[, pa.this], family = "binomial", control=glm.control(maxit=50))
        pa.bfit <- as.vector(pa.lmod$coefficients[2])
      }else{
        lmodPois <- glm(pa.X[, pa.now] ~ pa.X[, pa.this], family = "poisson", control=list(maxit=100))
        pa.bfit <- as.vector(lmodPois$coefficients[2])
      }
      
    }else{
      
      pa.lmod <- glm(pa.X[, pa.now] ~ pa.X[, pa.this], family = "poisson", control=list(maxit=100))
      pa.bfit <- as.vector(pa.lmod$coefficients[2])
      
    }
    
    if(pa.bfit != 0){
      pa.result[pa.this+1] <- 1
    }
    
  }
  
  return(pa.result)
}




Trans_parents <- function(node_pa_col, TD.p){
  
  TDS.now <- node_pa_col[1]
  Trans_result <- node_pa_col[2:(TD.p+1)]
  
  return(Trans_result)
}



cbind_parents <- function(node_list, p, pa_num){
  
  if(length(node_list) == 1){
    cbind_re <- as.matrix(node_list[[1]])
    
  }else{
    
    cbind_re <- as.matrix(node_list[[1]])
    list_length <- length(node_list)
    
    for (i in 2:list_length) {
      cbind_re <- cbind(cbind_re, as.matrix(node_list[[i]]))
    }
  }
  
  
  
  return(cbind_re)
}






