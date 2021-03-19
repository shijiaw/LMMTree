#install.packages("seqinr")
rm(list = ls())
library("seqinr")
library(stringr)
ncrna <- read.fasta(file = "out475_snps.fa")
length(ncrna)
nc = 20976
rind <- read.csv("common_Index_fa.csv")[,-1]
#colind <- read.csv("MAF01.csv")[,-1]
singular <- read.csv("singular.csv")[,-1]

matrixDNA <- matrix(NA, nr = length(rind), nc = length(singular))
for(i in 1:length(rind)){
  temp <- ncrna[[i]][singular]
  matrixDNA[i, ] <- temp
  matrixDNA[i,which(matrixDNA[i,] == "a")] <- "A"
  matrixDNA[i,which(matrixDNA[i,] == "c")] <- "C"
  matrixDNA[i,which(matrixDNA[i,] == "g")] <- "G"
  matrixDNA[i,which(matrixDNA[i,] == "t")] <- "T"
  matrixDNA[i,which(matrixDNA[i,] == "n")] <- "-"
}


write.nexus.data(matrixDNA, 'file.nex')

