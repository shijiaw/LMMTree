# LMMTree
Linear mixed model with uncertainty in genetic similarity matrices

Summary
-------

Genome-wide association studies are often confounded by population stratification and structure. Linear mixed models (LMMs) are a powerful class of methods for uncovering genetic effects, while controlling for such confounding.  LMMs use the genetic similarity matrix as a fixed effect, and they assume that a true genetic similarity matrix is known. However, uncertainty about the phylogenetic structure of a study population may degrade the quality of LMM results. This may happen in bacterial studies in which the number of samples or loci are small, or in studies with low quality genotyping. In this work, we develop methods for linear mixed models in which the genetic similarity matrix is unknown and is derived from MCMC estimates of the phylogeny. We apply our model to a genome-wide association study of multidrug-resistance in tuberculosis, and illustrate our methods on simulated data.

- Simulation: a folder with code to reproduce results of our simulation study
- TB: code for MDR-GWAS. 

Simulation
------------


There are several files in folder simulation:

### For each Scenario (S)

- run`` treesim.R'' to simulate trees and sequences
- run ``phesim.R'' to simulate phenotype
- run ``mb.R'' to estimate posterior of phylogenetic trees
- use ``https://github.com/shijiaw/Expected-Genetic-Similarity-Matrices'' to get the estimated GSM with the thinnined posterior tree samples 
- run ``comparison.R'' for GWAS on simulated data set. Method `LM' denotes linear regression model, `LMM' denotes linear mixed model with empirical genetic similarity matrix (inner product of SNP matrices), `LiMU' denotes linear mixed model with expected genetic similarity matrices (computed from posterior samples of phylogenies), `pyseer' denotes the ﬁxed eﬀects model with the genetic similarity matrix represented by MDS implemented in the pyseer software.
- ``pagg.R'' includes our code for plots


Software
------------

The package folder includes an R function that implement LiMU

- ``LiMU.R'' is the main function for LiMU method, which include Mr. Bayes run for posterior trees, computation of GSM and run LMMs via EMMA method.
- ``emma.R'' is the source code for EMMA method.
-``demo.R" shows an example for running LiMU with specific inputs.
- folder ``src'' includes java implementation for computation of expected GSM.
- folder ``jars'' includes jar files that is required for GSM computation.
- folder ``DataDemo'' includes data for ``demo.R", in which ``seq1.nex" is a .nex file includes DNA sequences for Mr. Bayes, ``y_1.csv" is phenotype data and ``SNP1.csv" is SNPs. 
- After running 'LiMU.R' you will get p-values of LiMU in folder 'Pvalue' (i.e. 'minuslog10P.csv' includes -log_10(p-value) for each site summarized by mean statistics, 'originalP.csv' includes all original p-values).


