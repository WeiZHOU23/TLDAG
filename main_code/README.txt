The program files are for the implementation of our TLDAG method proposed in the paper 
"Efficient Learning of Quadratic Variance Function Directed Acyclic Graphs via Topological Layers".

A complete list of the programs:
DataGenerated.r

kappa_tune.r
ratio_cal.r
TLDAG.r
experiment.r
simulation.r


The main program is in simulation.r. It will load the functions and R packages, and perform a simulation 
trial of the paper. Moreover, rtools is also needed to download and the website is given in head.r. One can 
perform simulations with different settings by modifying the parameters given in the paper. Be careful for 
the parameter settings as they should satisfy the theoretical conditions such that the generated values do 
not blow up. Other files contain the support functions for method implementation and result evaluation.
