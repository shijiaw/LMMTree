rm(list=ls())
# simluate phenotype
library("phyclust")
library("MASS")
library(emma)
#library(nlme)
#library(lme4)


N = 1
M = 469
nunc <- 50
#L = 500
# P_LM <- matrix(NA, nr = N, nc = L)
# P_LMM <- matrix(NA, nr = N, nc = L)
# P_LMMT <- matrix(NA, nr = N, nc = L)
P_LM <- list()
P_LMM <- list()
P_LMMT <- list()
P_LMMTmat <- list()

set.seed(321)

#readseq_dir <- paste('output_SNP/SNP',i ,'.csv', sep = '')
y <- read.csv('MDR.csv', sep = ',', header = TRUE)[,-1]
seq <- read.csv('x.csv', sep = ',', header = TRUE)[,-1]

MAF <- rep(NA, ncol(seq))
for(i in 1:ncol(seq)){
  MAF[i] <- sum(seq[,i])/nrow(seq)
}

MAF002 <- which(MAF >= 0.02)
write.csv(MAF002, file = "MAF002.csv")

#mean_mat <- apply(seq, 2, mean)
#sd_mat <- sqrt(mean_mat*(1-mean_mat))
X <- apply(seq, 2, scale)
#phe_name <- paste('output_phenotype/y_', i, '.csv', sep = "")
#y <- unlist(read.csv(phe_name, sep = ",", header = FALSE))

L <- length(MAF002)
P_LMtemp <- rep(NA, L)
P_LMMtemp <- rep(NA, L)
P_LMMTtemp <- rep(NA, L)
P_LMMTmattemp <- matrix(NA, nunc, L)

for(l1 in 1:L){
  out <- summary(lm(y~X[,MAF002[l1]]-1))
  # we store the pvalue in minus log form with 10 base
  P_LMtemp[l1] <- -log10(out$coefficients[,4])
}

# LMM with GSM
K <- X %*% t(X)/L
rs <- emma.REML.t(matrix(y, nr = 1), t(seq[,MAF002]), K)
#rs2 <- emma.REML.t(matrix(y, nr = 1), t(X), K)
P_LMMtemp <- -log10(rs$ps)

# LMM with EGSM
for(index in 1:nunc){
  readK_dir <- paste('mb_tree/K_',index,'.csv', sep = '')
  #K_tree <- read.csv(readK_dir, sep = ",", header = FALSE)
  K_tree <- matrix(unlist(read.csv(readK_dir, sep = ",", header = FALSE)), nr = M)
  rsT <- emma.REML.t(matrix(y, nr = 1), t(seq[,MAF002]), K_tree)
  P_LMMTmattemp[index,] <- -log10(rsT$ps)
  print(index)
}





filename_LM = paste('Pvalue/LM.csv', sep = "")
write.table(P_LMtemp, file = filename_LM, row.names = FALSE, col.names = FALSE)
filename_LMM = paste('Pvalue/LMM.csv', sep = "")
write.table(as.vector(P_LMMtemp), file = filename_LMM, row.names = FALSE, col.names = FALSE)
filename_LMMT = paste('Pvalue/LMMT.csv', sep = "")
write.table(apply(P_LMMTmattemp, 2, median), file = filename_LMMT, row.names = FALSE, col.names = FALSE)
filename_LMMTun = paste('Pvalue/LMMTun.csv', sep = "")
write.table(matrix(P_LMMTmattemp, nr = 1), file = filename_LMMTun, row.names = FALSE, col.names = FALSE)

#for(i in 1:N){
#  filename_LMMT = paste('PvalueUncertainty/LMMTmat_', i,'.csv',sep = "")
#  write.table(matrix(P_LMMTmat[[i]], nr = 1), file = filename_LMMT, row.names = FALSE, col.names = FALSE)
#}
