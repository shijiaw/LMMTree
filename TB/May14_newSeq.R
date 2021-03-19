library("seqinr")
#ncrna <- read.fasta(file = "1Nov.fa")

snps <- read.fasta(file = "out475_snps.fa")

length(snps)
nc = 20976
matrixDNA <- matrix(NA, nr = 471, nc = nc)
for(i in 1:471){
  matrixDNA[i, ] <- snps[[i]]
}

# Find the major allels
ma <- rep(NA, nc)
for(i in 1:nc){
  ma[i] <- names(which.max(table(matrixDNA[,i])))
}

SNPmatrix <- matrix(NA, nr = 471, nc = nc)
for(i in 1:nc){
  SNPmatrix[ ,i] <- 1-as.numeric(matrixDNA[,i] == ma[i])
}

stdSNPmatrix <- matrix(NA, nr = 471, nc = nc)
for(i in 1:nc){
  stdSNPmatrix[ ,i] <- (SNPmatrix[ ,i] - mean(SNPmatrix[ ,i]))/sd(SNPmatrix[ ,i])
}

commonSNP_Index <- read.csv("common_Index_fa.csv")[,-1]
write.csv(SNPmatrix[commonSNP_Index,], "SNP.csv")

