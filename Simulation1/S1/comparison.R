rm(list=ls())
# simluate phenotype
library("phyclust")
library("MASS")
library(emma)
#library(nlme)
#library(lme4)


N = 50
M = 300
L = 500
P_LM <- matrix(NA, nr = N, nc = L)
P_LMM <- matrix(NA, nr = N, nc = L)
P_LMMT <- matrix(NA, nr = N, nc = L)

set.seed(321)

for(i in 1:50){
  readseq_dir <- paste('output_seq/seq',i ,'.txt', sep = '')
  seq_i <- do.call(rbind,lapply(strsplit(scan(readseq_dir, 
                                              what=""), ""), as.numeric))
  mean_mat <- apply(seq_i, 2, mean)
  sd_mat <- sqrt(mean_mat*(1-mean_mat))
  X <- apply(seq_i, 2, scale)
  phe_name <- paste('output_phenotype/y_', i, '.csv', sep = "")
  y <- unlist(read.csv(phe_name, sep = ",", header = FALSE))
  
  for(l1 in 1:L){
    out <- summary(lm(y~X[,l1]-1))
    # we store the pvalue in minus log form with 10 base
    P_LM[i, l1] <- -log10(out$coefficients[,4])
  }
  
  # LMM with GSM
  K <- X %*% t(X)/L
  rs <- emma.REML.t(matrix(y, nr = 1), t(seq_i), K)
  #rs2 <- emma.REML.t(matrix(y, nr = 1), t(X), K)
  P_LMM[i, ] <- -log10(rs$ps)
  
  # LMM with EGSM
  readK_dir <- paste('output_tree/K_',i ,'.csv', sep = '')
  #K_tree <- read.csv(readK_dir, sep = ",", header = FALSE)
  K_tree <- matrix(unlist(read.csv(readK_dir, sep = ",", header = FALSE)), nr = M)
  rsT <- emma.REML.t(matrix(y, nr = 1), t(seq_i), K_tree)
  P_LMMT[i, ] <- -log10(rsT$ps)
  print(i)
  
  #write.table(y, file = filename, row.names = FALSE, col.names = FALSE)
}
filename_LM = paste('Pvalue/LM.csv', sep = "")
write.table(P_LM, file = filename_LM, row.names = FALSE, col.names = FALSE)
filename_LMM = paste('Pvalue/LMM.csv', sep = "")
write.table(P_LMM, file = filename_LMM, row.names = FALSE, col.names = FALSE)
filename_LMMT = paste('Pvalue/LMMT.csv', sep = "")
write.table(P_LMMT, file = filename_LMMT, row.names = FALSE, col.names = FALSE)

