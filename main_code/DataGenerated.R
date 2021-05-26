
####################################################
#### generate a random graph by ER and BA model  ###
####################################################
## p: the number of nodes
## n: the sample size for data generation
## BA.e: the number of edges to add in each time step
## ER.pE: connected probability 

simul.DAG.random <- function(n, p, BA.e, ER.pE, simul.model = c("ER","BA"), DAG.seed){
  
  if(simul.model == "ER"){
    
    stopifnot(ER.pE > 0, ER.pE < 1)
    A <- matrix(0, p, p)
    w <- which(lower.tri(A))
    set.seed(DAG.seed)
    A[w] <- rbinom(n = length(w), size=1, prob = ER.pE) * runif(n = length(w), min=0.1, max=1)
    
  }else{
    
    set.seed(DAG.seed)
    A <- as(igraph.to.graphNEL(barabasi.game(n = p, m = BA.e)), "matrix")
    w <- which(A != 0)
    set.seed(DAG.seed)
    A[w] <-  runif(n = length(w), min=0.1, max=1)
    
  }
  
  ## shuffling
  set.seed(DAG.seed)
  id.shuffle <- sample(1:p, p, replace = FALSE)
  A <- A[id.shuffle, id.shuffle]
  
  A.binary <- Value2Binary(A)
  A = A.binary
  
  #return(list(A = A.binary, G = as(t(A),"graphNEL")))
  return(A)
}



## transform the adjacency matrix with values into binary one in random graph
Value2Binary <- function(A){
  
  p <- ncol(A)
  TF_A <- sapply(1:p^2, function(i) ifelse(A[i] != 0, 1, 0))
  Binary_A <- matrix(TF_A, p, p, byrow = FALSE)
  
  return(Binary_A)
}


#######################################################
#### generate data in random graph for mixed graph ####
#######################################################

simul.data.random <- function(n, p, BA.e, ER.pE, simul.model = c("ER","BA"), 
                              theta1.pois, theta2.pois, beta1.pois, beta2.pois, 
                              theta1.bino, theta2.bino, bino.N, simul.seed){
  
  ### get the DAG structure 
  simul.dag <- simul.DAG.random(n, p, BA.e, ER.pE, simul.model, simul.seed)
  simul.A <- simul.dag
  #simul.A <- simul.dag$A
  #simul.graph <- simul.dag$G
  simul.layer <- GetLayer(simul.A)$Layer
  
  ### data generation
  simul.X <- matrix(0, n, p)
  simul.index <- list()
  simul.prob <- list()
  
  simul.X.init <- matrix(floor(bino.N*runif(n*p, min = 0, max = 1)), n, p)
  
  ## generate the parameters
  set.seed(simul.seed)
  theta.pois <- runif(p, theta1.pois, theta2.pois)
  set.seed(simul.seed)
  beta.pois <- matrix(runif(p*p, beta1.pois, beta2.pois), p ,p)
  
  set.seed(simul.seed)
  theta.bino <- runif(p, theta1.bino, theta2.bino)
  set.seed(simul.seed)
  beta.bino  <- matrix(runif(p*p, theta1.bino, theta2.bino), p ,p)
  
  
  ## generate data for root nodes
  simul.A1.nodes <- simul.layer[[1]]
  simul.A1.length <- length(simul.A1.nodes)
  set.seed(simul.seed)
  simul.A1.index <- sample(2, size = simul.A1.length, replace = TRUE)

  
  simul.index[[1]] <- simul.A1.index
  
  
  for (root in 1:simul.A1.length) {
    
    A1.root.index <- simul.A1.nodes[root]
    
    if(simul.A1.index[root] == 1){
      ## poission root
      A1.lambda <- exp(theta.pois[A1.root.index])
      set.seed(simul.seed)
      simul.X[, A1.root.index] <- rpois(n, lambda = A1.lambda)
      
    }else{
      ## binomial root
      A1.prob <- exp(theta.bino[A1.root.index]) / (1+exp(theta.bino[A1.root.index]))
      set.seed(simul.seed)
      simul.X[, A1.root.index] <- rbinom(n, size = bino.N, prob = A1.prob)
      
    }
  }
  
  
  ## generate data for A_t layer
  if(length(simul.layer) >= 2){
    
    for(i in 2:length(simul.layer)){
      
      simul.At.nodes <- simul.layer[[i]]
      simul.At.length <- length(simul.At.nodes)
      set.seed(simul.seed)
      simul.At.index <- sample(2, simul.At.length, replace = TRUE)
      simul.index[[i]] <- simul.At.index
      
      for (node in 1:simul.At.length){
        
        At.node.index <- simul.At.nodes[node]
        At.parent <- which(simul.A[, At.node.index] != 0)
        
        if(simul.At.index[node] == 1){
          ## poission nodes
          
          if(length(At.parent) == 1){
            
            At.lambda <- exp(theta.pois[At.node.index] + simul.X[, At.parent] * beta.pois[At.node.index, At.parent])
            set.seed(simul.seed)
            simul.X[, At.node.index] <- rpois(n, lambda = At.lambda)
            
          }else{
            
            At.lambda <- exp(theta.pois[At.node.index] + simul.X[, At.parent] %*% beta.pois[At.node.index, At.parent])
            set.seed(simul.seed)
            simul.X[, At.node.index] <- rpois(n, lambda = At.lambda)
            
          }
          
          
        }else{
          ## binomial nodes

          if(length(At.parent) == 1){
            At.deno <- exp(theta.bino[At.node.index] + simul.X[, At.parent] * beta.bino[At.node.index, At.parent])
            At.prob <- At.deno / (1+At.deno)
            set.seed(simul.seed)
            simul.X[, At.node.index] <- rbinom(n, size = bino.N, prob = At.prob)
            
          }else{
            At.deno <- exp(theta.bino[At.node.index] + simul.X[, At.parent] %*% beta.bino[At.node.index, At.parent])
            At.prob <- At.deno / (1+At.deno)
            set.seed(simul.seed)
            simul.X[, At.node.index] <- rbinom(n, size = bino.N, prob = At.prob)
            
          }
        }
      }
    }
  }
  
  
  
  Index.matrix <- matrix(0, p, 2)
  Index.matrix[, 1] <- unlist(simul.layer)
  Index.matrix[, 2] <- unlist(simul.index)
  
  all.IndexBino <- Index.matrix[which(Index.matrix[, 2] == 2), 1]
  
  #return(list(X = simul.X, DAG = simul.A, Layer = simul.layer, BProb = simul.prob, XBino = all.IndexBino, G = simul.graph, Distri = simul.index))
  return(list(X = simul.X, DAG = simul.A, Layer = simul.layer, BProb = simul.prob, XBino = all.IndexBino, Distri = simul.index))
  
}






#####################################
#### generate data for hub graph ####
#####################################


###  hub graph for poisson distribution 
hubPoisGraph <- function(n, p, theta1.pois, theta2.pois, beta1.pois, beta2.pois, h.seed){
  ordering <- seq(p)
  truth <- matrix(0, p, p)
  truth[1, 2:p] <- 1
  X <- matrix(0, n, p)
  
  set.seed(h.seed)
  theta <- runif(p, theta1.pois, theta2.pois)
  set.seed(h.seed)
  beta <- matrix(runif(p*p, beta1.pois, beta2.pois), p, p)
  
  set.seed(h.seed)
  X[, 1] <- rpois(n, lambda = exp(theta[1]))
  
  for (j in 2:p) {
    pois.lambda <- exp(theta[j] + X[, 1] * beta[j, 1])
    set.seed(h.seed)
    X[, j] <- rpois(n, lambda = pois.lambda)
  }
  
  return(list(X = X, DAG = truth, Distri = 1, Order = c(1, rep(2,p-1))))
}






### hub graph for mixed distributions
hubMixedGraph <- function(n, p, theta1.pois, theta2.pois, beta1.pois, beta2.pois, theta1.bino, theta2.bino, bino.N, hub.seed){
  
  ordering <- seq(p)
  hub.DAG <- matrix(0, p, p)
  hub.DAG[1, 2:p] <- 1
  hub.X <- matrix(0, n, p)
  
  ## generate the parameters
  set.seed(hub.seed)
  theta.pois <- runif(p, theta1.pois, theta2.pois)
  set.seed(hub.seed)
  beta.pois <- matrix(runif(p*p, beta1.pois, beta2.pois), p ,p)
  
  set.seed(hub.seed)
  theta.bino <- runif(p, theta1.bino, theta2.bino)
  set.seed(hub.seed)
  beta.bino  <- matrix(runif(p*p, theta1.bino, theta2.bino), p ,p)
  
  ## generate 
  set.seed(hub.seed)
  hubMix.index <- sample(2, p, replace = TRUE)
  
  
  ## generate date for the root node
  if(hubMix.index[1] == 1){
    ## poisson root node
    hub.lambda <- exp(theta.pois[1])
    set.seed(hub.seed)
    hub.X[, 1] <- rpois(n, lambda = hub.lambda)
    
  }else{
    ## binomial root node
    hub.prob <- exp(theta.bino[1]) / (1+exp(theta.bino[1]))
    set.seed(hub.seed)
    hub.X[, 1] <- rbinom(n, size = bino.N, prob = hub.prob)
    
  }
  
  
  for(i in 2:p){
    
    if(hubMix.index[i] == 1){
      ### poisson node
      pois.lambda <- exp(theta.pois[i] + hub.X[, 1] * beta.pois[i, 1])
      set.seed(hub.seed)
      hub.X[, i] <- rpois(n, lambda = pois.lambda)

    }else{
      ### binomial node
      prob.bino <- exp(theta.bino[i] + hub.X[, 1] * beta.bino[i, 1]) / (1+exp(theta.bino[i] + hub.X[, 1] * beta.bino[i, 1]))
      set.seed(hub.seed)
      hub.X[, i] <- rbinom(n, size = bino.N, prob = prob.bino)
    }
  }
  
  
  return(list(X = hub.X, DAG = hub.DAG, Distri = hubMix.index, Order = c(1, rep(2,p-1))))
}








  



  





