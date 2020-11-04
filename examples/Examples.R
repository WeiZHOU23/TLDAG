#########################################
####  All the experiments 
######################################


grid.t <- 10^(-2 + 0.15*seq(0, 60, by=1))

## the simulation result for n={200,500} and p={5,20,100} respectively


###########################################
### Example 1 for dense hub poisson graph 
###########################################
source("HubPoisson.R")

a25 <- simu_replicate(simr.n=200, simr.p=5, simu.times=50, simr.grid.t=grid.t)
a55 <- simu_replicate(simr.n=500, simr.p=5, simu.times=50, simr.grid.t=grid.t)

b25 <- simu_replicate(simr.n=200, simr.p=20, simu.times=50, simr.grid.t=grid.t)
b55 <- simu_replicate(simr.n=500, simr.p=20, simu.times=50, simr.grid.t=grid.t)

c25 <- simu_replicate(simr.n=200, simr.p=100, simu.times=50, simr.grid.t=grid.t)
c55 <- simu_replicate(simr.n=500, simr.p=100, simu.times=50, simr.grid.t=grid.t)




###########################################
### Example 2 for dense hub mixed graph 
###########################################
source("HubMixed.R")

a1<-simu_replicate(simr.n=200, simr.p=5, simr.N=4, simr.sparsity="dense", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                   simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)
a3<-simu_replicate(simr.n=500, simr.p=5, simr.N=4, simr.sparsity="dense", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                   simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.3)


b1 <- simu_replicate(simr.n=200, simr.p=20, simr.N=4, simr.sparsity="dense", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                     simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)
b3 <- simu_replicate(simr.n=500, simr.p=20, simr.N=4, simr.sparsity="dense", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                     simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)


c1 <- simu_replicate(simr.n=200, simr.p=100, simr.N=4, simr.sparsity="dense", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                     simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)
c3 <- simu_replicate(simr.n=500, simr.p=100, simr.N=4, simr.sparsity="dense", simr.times=50, simr.grid.t=grid.t, simr.theta1.pois=1, simr.theta2.pois=3, 
                     simr.theta1.pois.beta=0.1, simr.theta2.pois.beta=0.2, simr.theta1.bino=0.1, simr.theta2.bino=0.2)





###############################################
### Example 3 for dense random binomial graph 
###############################################
source("RanBino.R")

a1 <- simu_replicate(simr.n=200, simr.p=5, simr.times=50, simr.grid.t=grid.t, simr.T=2, simr.graph='dense', simr.theta1=0.1, simr.theta2=0.5, simr.N=4)
a3 <- simu_replicate(simr.n=500, simr.p=5, simr.times=50, simr.grid.t=grid.t, simr.T=2, simr.graph='dense', simr.theta1=0.1, simr.theta2=0.5, simr.N=4)

c1 <- simu_replicate(simr.n=200, simr.p=20, simr.times=50, simr.grid.t=grid.t, simr.T=2, simr.graph='dense', simr.theta1=0.1, simr.theta2=0.5, simr.N=4)
c3 <- simu_replicate(simr.n=500, simr.p=20, simr.times=50, simr.grid.t=grid.t, simr.T=2, simr.graph='dense', simr.theta1=0.1, simr.theta2=0.5, simr.N=4)


h1 <- simu_replicate(simr.n=200, simr.p=100, simr.times=50, simr.grid.t=grid.t, simr.T=3, simr.graph='dense', simr.theta1=0.01, simr.theta2=0.1, simr.N=4)
h3 <- simu_replicate(simr.n=500, simr.p=100, simr.times=50, simr.grid.t=grid.t, simr.T=3, simr.graph='dense', simr.theta1=0.01, simr.theta2=0.1, simr.N=4)



###############################################
### Example 4 for dense random mixed graph 
###############################################
source("RanMixed.R")


a1 <- simu_replicate(200, 5, 50, grid.t, 2, simr.graph = 'dense', 1, 3, 0.01, 0.02, 0.01, 0.02, 4)
a3 <- simu_replicate(500, 5, 50, grid.t, 2, simr.graph = 'dense', 1, 3, 0.01, 0.02, 0.01, 0.02, 4)

c1 <- simu_replicate(200, 20, 50, grid.t, 2, simr.graph = 'dense', 1, 3, 0.01, 0.02, 0.01, 0.02, 4)
c3 <- simu_replicate(500, 20, 50, grid.t, 2, simr.graph = 'dense', 1, 3, 0.01, 0.02, 0.01, 0.02, 4)

k1 <- simu_replicate(200, 100, 50, grid.t, 3, simr.graph = 'dense', 1, 3, 0.0001, 0.001, 0.0001, 0.001, 4)
k3 <- simu_replicate(500, 100, 50, grid.t, 3, simr.graph = 'dense', 1, 3, 0.0001, 0.001, 0.0001, 0.001, 4)




###############################################
### Example 5 for sparse random mixed graph 
###############################################
source("RanMixed.R")

a1s <- simu_replicate(200, 5, 50, grid.t, 2, simr.graph = 'sparse', 1, 3, 0.01, 0.02, 0.01, 0.02, 4)
a3s <- simu_replicate(500, 5, 50, grid.t, 2, simr.graph = 'sparse', 1, 3, 0.01, 0.02, 0.01, 0.02, 4)  

c1s <- simu_replicate(200, 20, 50, grid.t, 2, simr.graph = 'sparse', 1, 3, 0.01, 0.02, 0.01, 0.02, 4)
c3s <- simu_replicate(500, 20, 50, grid.t, 2, simr.graph = 'sparse', 1, 3, 0.01, 0.02, 0.01, 0.02, 4)

k1s <- simu_replicate(200, 100, 50, grid.t, 3, simr.graph = 'sparse', 1, 3, 0.0001, 0.005, 0.0001, 0.001, 4)
k3s <- simu_replicate(500, 100, 50, grid.t, 3, simr.graph = 'sparse', 1, 3, 0.0001, 0.005, 0.0001, 0.001, 4)
