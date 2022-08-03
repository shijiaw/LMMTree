# LMMTree
Linear mixed model with uncertainty in genetic similarity matrices

Summary
-------

Genome-wide association studies are often confounded by population stratification and structure. Linear mixed models (LMMs) are a powerful class of methods for uncovering genetic effects, while controlling for such confounding.  LMMs use the genetic similarity matrix as a fixed effect, and they assume that a true genetic similarity matrix is known. However, uncertainty about the phylogenetic structure of a study population may degrade the quality of LMM results. This may happen in bacterial studies in which the number of samples or loci are small, or in studies with low quality genotyping. In this work, we develop methods for linear mixed models in which the genetic similarity matrix is unknown and is derived from MCMC estimates of the phylogeny. We apply our model to a genome-wide association study of multidrug-resistance in tuberculosis, and illustrate our methods on simulated data.

This repository is organized as follows:

- `Simulation`: a folder with code to reproduce results of our simulation study (in R).
- `TB`: code for the experiment on multidrug-resistance in tuberculosis (in R).
- `src/GSM`: code for estimation of the kinship matrix and LMM effect sizes (in java).

Simulation
------------


There are several files in folder simulation:

### For each Scenario (S)

- Run `treesim.R` to simulate trees and sequences.
- Run `phesim.R` to simulate phenotype.
- Run `mb.R` to estimate posterior of phylogenetic trees.
- Use `https://github.com/shijiaw/Expected-Genetic-Similarity-Matrices` (Wang et al. 2021) to get estimate genetic similarity matrices with the thinnined posterior tree samples.
- Run `comparison.R` to conduct GWAS on the simulated data set. Method 'LM' denotes linear regression model, 'LMM' denotes linear mixed model with an empirical genetic similarity matrix (inner product of SNP matrices), 'LiMU' denotes linear mixed model with expected genetic similarity matrices (computed from posterior samples of phylogenies), 'pyseer' denotes the ﬁxed eﬀects model with the genetic similarity matrix represented by MDS implemented in the pyseer software.
- The script `pagg.R` includes code for plots.


Software
------------

The package folder includes an R function that implement LiMU, the code is organized as follows:

- `LiMU.R` is the main function for LiMU method, which includes Mr. Bayes run for posterior trees, computation of genetic similarity matrices, and LMMs through EMMA method.
- `emma.R` is the source code for EMMA method.
- `demo.R` provides an example for running LiMU with specific inputs.
- The folder `src` includes a java implementation for computation of the expected genetic similarity matrix.
- The folder `jars` includes required libraries (as jar files).
- The folder `DataDemo` includes data for `demo.R`, in which `seq1.nex` is a .nex file includes DNA sequences for Mr. Bayes, `y_1.csv` is example phenotype data, and `SNP1.csv` is example SNP data. 
- After running `LiMU.R`, we get p-values for LiMU output in the folder `Pvalue` (in particular, `minuslog10P.csv` includes -log_10(p-value) for each site summarized by mean statistics, and `originalP.csv` includes all original, non-transformed p-values).

Citations
-----------
S. Wang, S. Ge, C. Colijn, P. Biller, L. Wang, and L.T. Elliott. Estimating Genetic Similarity Matrices Using Phylogenies. 2021. Journal of Computational Biology. 28(6)

