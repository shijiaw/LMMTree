##make data for pyseer
library(ggplot2)
library(ggpubr)
library("seqinr")
library(stringr)
library("phyclust")

snps <- read.fasta(file = "out475_snps.fa")

data <- read.csv("combined_metadata_470_strains_24032016.csv", header = TRUE)

commonnames <- intersect(snps$seqname,data$tip_labels2)

common_Index <- which(data$tip_labels2%in%commonnames)

common_Index_fa <- which(names(snps)%in%commonnames)

MDR <- 1-as.numeric(data$drug_susceptibility_result == 'MDR')
samples <- paste0("ID_", 1:length(MDR[common_Index]))
res <- MDR[common_Index]
pheno <- cbind(samples, res)
write.table(pheno, "resistance.pheno", row.names = FALSE, quote=FALSE, sep = '\t')

##genotype
#length(snps)
length(snps)
nc = 20976
matrixDNA_temp <- matrix(NA, nr = 471, nc = nc)
for(i in 1:471){
  matrixDNA_temp[i, ] <- snps$org.code[i,]
}
commonSNP_Index <- read.csv("common_Index_fa.csv")[,-1]
matrixDNA <- matrixDNA_temp[commonSNP_Index, ]

# Find the major allels
ma <- rep(NA, nc)
for(i in 1:nc){
  ma[i] <- names(which.max(table(matrixDNA[,i])))
}
NAindex <- which(ma == 'N')


SNPmatrix <- matrix(NA, nr = 467, nc = nc)
for(i in 1:nc){
  SNPmatrix[ ,i] <- 1-as.numeric(matrixDNA[,i] == ma[i])
}
SNPmatrix <- SNPmatrix[ ,-NAindex]

seq <- SNPmatrix
MAF <- rep(NA, ncol(seq))
for(i in 1:ncol(seq)){
  MAF[i] <- sum(seq[,i])/nrow(seq)
}

MAF005 <- which(MAF >= 0.005)

geno <- rbind(paste0("ID_", 1:length(MDR[samples])), t(SNPmatrix[,MAF005]))

write.table(cbind(c("gene",paste0("gene_", 1:ncol(SNPmatrix[,MAF005]))), geno), "SNP.Rtab", row.names = FALSE,col.names = FALSE, quote=FALSE, sep = '\t')

## VCF
VCFmat00 <- seq[,MAF005]
VCFmat <- VCFmat00 
for(i in 1:nrow(VCFmat00)){
  VCFmat[i,] <- paste0(VCFmat00[i,], "/", VCFmat00[i,])
}
VCFgeno <- t(VCFmat) 
VCFgenofmat <- matrix(NA, nr = nrow(VCFgeno)+1, nc = ncol(VCFgeno)+9)
nnr <- nrow(VCFgeno)
VCFgenofmat[,1] <- c('#CHROM', rep(11126, nnr))
VCFgenofmat[,2] <- c('POS', 1:nnr)
VCFgenofmat[,3] <- c('ID', 1:nnr)
VCFgenofmat[,4] <- c('REF', rep('G', nnr))
VCFgenofmat[,5] <- c('ALT', rep('A', nnr))
VCFgenofmat[,6] <- c('QUAL', rep('.', nnr))
VCFgenofmat[,7] <- c('FILTER', rep('.', nnr))
VCFgenofmat[,8] <- c('INFO', rep('PR', nnr))
VCFgenofmat[,9] <- c('FORMAT', rep('GT', nnr))
VCFgenofmat[,10:(ncol(VCFgeno)+9)] <- rbind(paste0("ID_", 1:length(MDR[samples])), t(VCFmat))

fileConn<-file("SNP.vcf")
writeLines(c("##fileformat=VCFv4.2","##fileDate=20170502","##source=PLINKv1.90","##contig=<ID=11126,length=5560>", 
             '##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'), fileConn)
write.table(VCFgenofmat, "SNP.vcf", row.names = FALSE,col.names = FALSE, quote=FALSE, sep = '\t', append = TRUE)
close(fileConn) 
#write.table(t(VCFmat[commonSNP_Index,MAF002]), "SNP.vcf", row.names = FALSE,col.names = FALSE, quote=FALSE, sep = '\t')


## compute GSM
mysub <- function(test){
  test.new <- test
  index1 <- unlist(gregexpr(pattern="\\[", test,perl = TRUE)) 
  index2 <- unlist(gregexpr(pattern="\\]", test,perl = TRUE)) 
  substr(test,index1[2], index2[2])
  M <- substr(test,1,(index2[1]+1))
  test.new <- str_replace(test.new,fixed(M),"")
  #test.new <- str_replace(test.new,";" ,"")
  for(i in 2:length(index1)){
    M <- substr(test, index1[i],index2[i])
    test.new <- str_replace(test.new,fixed(M),"")
  }
  return(list(tree=test.new))
}


M = 467
skipLines = (M + 5)
tree_name = paste("file.nex.t", sep = '')
tree.anc <- read.csv(tree_name, skip = skipLines, sep="\t", header = FALSE)
ntree <- nrow(tree.anc)
#thinned_index <- (ntree-nuncertanty):(ntree-1)
tree_name2 = paste("file.nex.con.tre", sep = '')
tips <- read.csv(tree_name2, skip = 5, sep="\t", header = FALSE)[1:M,3]
tree <- read.csv(tree_name2, skip = 2*M+10, sep="\t", header = FALSE)[1,]
test <- mysub(tree)
phylo_name = "phylo.newick"
mytree <- read.tree(text = test$tree)
mytree$tip.label <- paste("ID_", mytree$tip.label, sep = '') 
write.tree(mytree, file = phylo_name, digits = 12)



#sudo python scripts/phylogeny_distance.py --lmm phylo.newick > MDS_K1.tsv
#pyseer --lmm --phenotypes resistance.pheno --vcf SNP.vcf --similarity MDS_K1.tsv --min-af 0.00 --max-af 1.00 > res.txt

