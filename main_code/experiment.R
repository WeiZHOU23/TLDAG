##########################
######  Experiment  ######
##########################


GetLayer <- function(A){
  
  layer.result <- list()
  vector.result <- c()
  p <- ncol(A)
  root_TF <- sapply(1:p, function(i) {all(A[, i] == 0)})
  root_nodes <- which(root_TF)
  layer.result[[1]] <- root_nodes
  vector.result[root_nodes] <- 1
  
  cont <- TRUE
  while( cont ){
    for (lay in 2:p) {
      
      upper_nodes <- unlist(layer.result)
      A_update <- A[-upper_nodes, -upper_nodes]
      curr_p <- ncol(A_update)
      node_update <- seq(p)[-upper_nodes]
      
      if( length(upper_nodes) < p-1){   
        
        curr_TF_update <- sapply(1:curr_p, function(i) {all(A_update[, i] == 0)})
        curr_TF_nodes <- which(curr_TF_update)
        curr_nodes_index <- node_update[curr_TF_nodes]
        layer.result[[lay]] <- curr_nodes_index
        vector.result[curr_nodes_index] <- lay
        
      }else if(length(upper_nodes) == p-1){
        
        layer.result[[lay]] <- node_update
        vector.result[node_update] <- lay
        
      }else(length(unlist(layer.result)) >= p)
      cont <- FALSE
      #break
    }
  }
  
 return(list(Layer = layer.result, LayerVec = vector.result)) 
}




  
 

## performance for one experiment with one method
do.criteria <- function(do.truth, do.estimate, do.order, do.esti.order, do.time){
  
  do.p <- ncol(do.truth)
  
  do.result <- c()
  do.result[1] <- do.time
  do.result[2] <- hammingDistance(do.estimate, do.truth)/(do.p*(do.p-1))
  do.result[3] <- ifelse(sum(do.truth), sum(do.truth*do.estimate)/sum(do.truth), 0) #recall
  do.result[4] <- ifelse(sum(do.truth), sum(do.truth*do.estimate)/(sum(do.estimate)+1e-6), 0) # precision
  do.result[5] <- 2*do.result[3]*do.result[4] / (do.result[3]+do.result[4]+1e-6)  # F1 score
  do.result[6] <- ifelse(sum(do.estimate), 1 - sum(do.truth*do.estimate)/(sum(do.estimate)+1e-6), 0) #FDR
  
  ## to deal with the null graph case, make the standard deviation non-zero for computing the kendall tau 
  if(all(do.esti.order == 1) | all(do.order == 1)){
    do.result[7] <- 0
  }else{
    do.result[7] <- cor(do.order, do.esti.order, method = "kendall")  ## Kendall's Tau-b 
    ## kendall tau b
  }
  
  return(do.result)
}



  one_experiment <- function(one.n, one.p, one.grid.t, one.theta1.pois, one.theta2.pois, one.beta1.pois, one.beta2.pois, 
                            one.theta1.bino, one.theta2.bino, one.bino.N, one.distri = c("pois","bino", "mix"), one.graph = c("hub","ER", "BA"), one.BA.e, one.ER.pE, one.seed){


  ## data generation
  if(one.graph == "hub" && one.distri == "pois"){
    one.data <- hubPoisGraph(one.n, one.p, one.theta1.pois, one.theta2.pois, one.beta1.pois, one.beta2.pois, one.seed)
  }

  if(one.distri == "mix"){

    if(one.graph == "hub"){
      one.data <- hubMixedGraph(one.n, one.p, one.theta1.pois, one.theta2.pois, one.beta1.pois, one.beta2.pois, one.theta1.bino, one.theta2.bino, one.bino.N, one.seed)
    }
    if(one.graph == "ER"){
      one.data <- simul.data.random(one.n, one.p, one.BA.e, one.ER.pE, simul.model = one.graph,
                                    one.theta1.pois, one.theta2.pois, one.beta1.pois, one.beta2.pois, one.theta1.bino, one.theta2.bino, one.bino.N, one.seed)
    }
    if(one.graph == "BA"){
      one.data <- simul.data.random(one.n, one.p, one.BA.e, one.ER.pE, simul.model = one.graph,
                                    one.theta1.pois, one.theta2.pois, one.beta1.pois, one.beta2.pois, one.theta1.bino, one.theta2.bino, one.bino.N, one.seed)
    }

  }


  one.X <- one.data$X
  one.DAG <- one.data$DAG
  one.order <- GetLayer(one.DAG)$LayerVec
  one.index <- one.data$Distri


  ## the estimated DAG by the proposed methods
  ptm <- proc.time()
  one.TLDAG <- TLDAG(one.X, one.grid.t, 1, 0, 1, -1/one.bino.N, one.bino.N, one.distri)
  one.TL.order <- GetLayer(one.TLDAG$DAG)$LayerVec
  one.perf.TLDAG <- do.criteria(one.DAG, one.TLDAG$DAG, one.order, one.TL.order, (proc.time() - ptm)[1])
   
  one.result <- one.perf.TLDAG

  return(one.result)
}


  
  
  
## repeat the experiments 
## simu_num: the number of times to repeat the trails
simu_replicate <- function(simu_num, one.n, one.p, one.grid.t, one.theta1.pois, one.theta2.pois, one.beta1.pois, one.beta2.pois, 
                             one.theta1.bino, one.theta2.bino, one.bino.N, one.distri = c("pois", "mix"), one.graph = c("hub","ER", "BA"), 
                             one.BA.e, one.ER.pE){
  

    
    pb <- txtProgressBar(min = 0, max = simu_num, style = 3)
    result <- matrix(0, 1, 7)
    results <-
      lapply(1:simu_num, function(i){
        setTxtProgressBar(pb, i)
        one_experiment(one.n, one.p, one.grid.t, one.theta1.pois, one.theta2.pois, one.beta1.pois, one.beta2.pois, 
                       one.theta1.bino, one.theta2.bino, one.bino.N, one.distri, one.graph, one.BA.e, one.ER.pE, one.seed = i)
      })
    close(pb)
    
    
    result <- Reduce('+', results)
    result.mean <- result / simu_num
    result.se <- apply(array(unlist(results), c(1, 7, simu_num)), c(1,2), sd) / sqrt(simu_num)
    result.se <- as.vector(result.se)
    
    names(result.mean) <- c("Time", "HM", "Recall", "Precision", "F1-score", "FDR", "Kendall")
    names(result.se) <- c("Time", "HM", "Recall", "Precision", "F1-score", "FDR", "Kendall")
    
    return(list(TLDAG.MEAN = result.mean, TLDAG.SE = result.se))
  }
  






