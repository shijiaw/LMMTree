# LMMTree
Linear mixed model with uncertainty in genetic similarity matrices

Summary
-------

Genome-wide association studies are often confounded by population stratification and structure. Linear mixed models (LMMs) are a powerful class of methods for uncovering genetic effects, while controlling for such confounding.  LMMs use the genetic similarity matrix as a fixed effect, and they assume that a true genetic similarity matrix is known. However, uncertainty about the phylogenetic structure of a study population may degrade the quality of LMM results. This may happen in bacterial studies in which the number of samples or loci are small, or in studies with low quality genotyping. In this work, we develop methods for linear mixed models in which the genetic similarity matrix is unknown and is derived from MCMC estimates of the phylogeny. We apply our model to a genome-wide association study of multidrug-resistance in tuberculosis, and illustrate our methods on simulated data.

- Simulation: a folder with code to reproduce results of our simulation study
- TB: code for MDR-GWAS. 
- Methods: `` LM'' denotes a linear regression model, `` LMM'' denotes a linear mixed effects model with empirical genetic similarity matrix (inner product of SNP sequences), `` LiMU'' denotes linear mixed model with uncertainty in genetic similarity matrices (using posterior samples to compute expected genetic similarity matrix). 

Simulation
------------


There are several files in folder simulation:

### For each Scenario (S)

- run`` treesim.R'' to simulate trees and sequences
- run ``phesim.R'' to simulate phenotype
- run ``mb.R'' to estimate posterior of phylogenetic trees
- use ``https://github.com/shijiaw/Expected-Genetic-Similarity-Matrices'' to get the estimated GSM with the thinnined posterior tree samples 
- run ``comparison.R'' for GWAS on simulated data set 
- ``pagg.R'' includes our code for plots

