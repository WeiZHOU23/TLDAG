getwd()
setwd()
# Load Packages and Programs
source("head.r")

## the grid set for stability selection for epsilon
grid.t <- 10^(-2 + 0.15*seq(0, 60, by=1))

###################################################################################################################
# Instruction for function"simu_replicate"
# simu_replicate(simu_num, one.n, one.p, one.grid.t, one.theta1.pois, one.theta2.pois, one.beta1.pois, one.beta2.pois, 
#    one.theta1.bino, one.theta2.bino, one.bino.N, one.distri, one.graph, one.BA.e, one.ER.pE)
###################################################################################################################
########## simu_num: number of repeat time;
########## one.n:    number of sample size;
########## one.p:    number of node size;
########## one.grid.t: The grid search vector for the tuning parameter epsilon_t;
########## one.theta1.pois, one.theta2.pois, one.beta1.pois, one.beta2.pois: the parameter for Poisson distribution;
########## one.theta1.bino, one.theta2.bino, one.bino.N:     the parameter  for Binomial distribution;
########## one.distri: the considering distibution, two choice:  pois and mix;
########## one.graph: the considering graph, three choice: hub, ER, BA;
########## one.BA.e, one.ER.pE: the number of edges to add for BA graph, and the probability of connecting an edge for ER graph
###################################################################################################################

## Example 1 for hub poisson graph when (p,n)=(5,200)
hubPoi_p5n200 <- simu_replicate(simu_num=50, one.n=200, one.p=5, one.grid.t=grid.t, one.theta1.pois=1, one.theta2.pois=3, one.beta1.pois=0.1, one.beta2.pois=0.5, 
                                 one.theta1.bino=0.5, one.theta2.bino=1, one.bino.N=4, one.distri = "pois", one.graph = "hub", one.BA.e=2, one.ER.pE=0.35)



## Example 2 for mixed hub graph when (p,n)=(5,200)
hubmix_p5n200 <- simu_replicate(simu_num=50, one.n=200, one.p=5, one.grid.t=grid.t, one.theta1.pois=1, one.theta2.pois=3, one.beta1.pois=0.1, one.beta2.pois=0.2, 
                                one.theta1.bino=0.1, one.theta2.bino=0.2, one.bino.N=4, one.distri = "mix", one.graph = "hub", one.BA.e=2, one.ER.pE=0.35)


## Example 3 for mixed ER graph when (p,n)=(5,200)
ER_p5n200 <- simu_replicate(simu_num=50, one.n=200, one.p=5, one.grid.t=grid.t, one.theta1.pois=1, one.theta2.pois=3, one.beta1.pois=0.01, one.beta2.pois=0.05, 
                            one.theta1.bino=0.01, one.theta2.bino=0.05, one.bino.N=4, one.distri = "mix", one.graph = "ER", one.BA.e=2, one.ER.pE=0.35)


## Example 4 for mixed BA graph when (p,n)=(5,200)
BA_p5n200 <- simu_replicate(simu_num=50, one.n=200, one.p=5, one.grid.t=grid.t, one.theta1.pois=1, one.theta2.pois=3, one.beta1.pois=0.01, one.beta2.pois=0.03, 
                            one.theta1.bino=0.01, one.theta2.bino=0.05, one.bino.N=4, one.distri = "mix", one.graph = "BA", one.BA.e=2, one.ER.pE=0.35)
