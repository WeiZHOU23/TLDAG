# TLDAG
 Learning DAG with EQ-DAG via the topological layers
# Abstract

We introduce a novel concept named topological layer of a DAG (TLDAG), which reformulates any DAG into
a unique topological structure. Then, we focus on a special class of directed acyclic graph
(DAG) models, where the conditional variance of each node given its parents is a quadratic
function of its conditional mean, and establish the identifiability via topological layers for
this class of DAG models. Note that the interested class of DAG models are fairly large that
many interesting classes of distributions satisfy this property, including Poisson, Binomial,
Geometric, Exponential, Gamma and many other distributions. Correspondingly, an efficient
learning method is proposed, which attempts to identify the topological layers in a hierarchical
fashion at first and then reconstruct the directed structures in a parallel fashion.

This repository maintains the code for this project

# The TLDAG package
The package contains the four cases (hub poisson graph, hub mixed graph for different sparsity, random binomial graph,
and random mixed graph for different sparsity). Among these methods, the ODS algorithm in Park and Raskutti (2018) is recoded
here. For GES and MMHC, we refer to the package pcalg and bnlearn. As the Direct LiNGAM is conducted with 
https://github.com/cdt15/lingam with Matlab and coded in Gnecco et al (2019 arXiv) with Rcpp in https://arxiv.org/abs/1908.05097,
a slightly modifief version is maintained in this repository.



# Simulation Studies

The details of the simulation settings are specified in our paper, also shown in the file "examples" 
for hub graphs and random graphs with n=\{200,500\} and p=\{5,20,100\} respectively.
