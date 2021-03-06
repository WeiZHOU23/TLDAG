# TLDAG
 Learning QVF-DAG via the topological layers
# Abstract


We study a special class of non-Gaussian DAG models, where the conditional variance of each node given its parents is a quadratic function of its conditional mean. Such a class of non-Gaussian DAG models are fairly flexible and admit many popular distributions as special cases, including Poisson, Binomial, Geometric, Exponential, and Gamma. To facilitate learning, we introduce a novel concept of topological layers, and develop an efficient DAG learning algorithm. It first identifies the topological layers in a hierarchical fashion and then reconstructs the directed edges between nodes in different layers, which requires much less computational cost than most existing algorithms in literature.

This repository maintains the code for this project

# The TLDAG package
The package contains the four cases (Poisson hub graph, mixed hub graph, ER random graph,
and BA random graph ). Among these methods, the ODS algorithm in Park and Raskutti (2018) and ODS algorithm in Park and Park (2019) is recoded here. For GES and MMHC, we refer to the package pcalg and bnlearn. As the Direct LiNGAM is conducted with 
https://github.com/cdt15/lingam with Matlab and coded in Gnecco et al (2019 arXiv) with Rcpp in https://arxiv.org/abs/1908.05097,
a slight modification version is maintained in this repository.



# Simulation Studies

The details of the simulation settings are specified in our paper, and an implementation for a toy example is given in the file "simulation.r".
Please be careful for installing the required packages before running the code.
