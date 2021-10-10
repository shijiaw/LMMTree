##make data for pyseer
library(ggplot2)
library(ggpubr)
library("seqinr")
library(stringr)
library("phyclust")

for(i in 1:50){
  res <- unlist(read.csv(paste0("output_phenotype/y_", i, ".csv"), header = FALSE))
  samples <- paste0("ID_", 1:length(res))
  pheno <- cbind(samples, res)
  write.table(pheno, paste0("pyseerData/y_", i, ".pheno"), row.names = FALSE, quote=FALSE, sep = '\t')
  
  VCFmat00 <- read.csv(paste0("output_SNP/SNP", i, ".csv"), sep = ' ', header = FALSE)
  VCFmat <- VCFmat00 
  for(j in 1:nrow(VCFmat00)){
    VCFmat[j,] <- paste0(VCFmat00[j,], "/", VCFmat00[j,])
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
  VCFgenofmat[,10:(ncol(VCFgeno)+9)] <- rbind(paste0("ID_", 1:length(res)), t(VCFmat))
  
  fileConn<-file(paste0("pyseerData/SNP_", i,".vcf"))
  writeLines(c("##fileformat=VCFv4.2","##fileDate=20170502","##source=PLINKv1.90","##contig=<ID=11126,length=500>", 
               '##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">',
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'), fileConn)
  write.table(VCFgenofmat, paste0("pyseerData/SNP_", i,".vcf"), row.names = FALSE,col.names = FALSE, quote=FALSE, sep = '\t', append = TRUE)
  close(fileConn) 
}


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

for(i in 1:50){
  M = 30
  skipLines = (M + 5)
  tree_name2 = paste("output_seq/seq", i,".nex.con.tre", sep = '')
  tips <- read.csv(tree_name2, skip = 5, sep="\t", header = FALSE)[1:M,3]
  tree <- read.csv(tree_name2, skip = 2*M+10, sep="\t", header = FALSE)[1,]
  test <- mysub(tree)
  phylo_name = paste0("pyseerData/phylo", i,".newick")
  mytree <- read.tree(text = test$tree)
  mytree$tip.label <- paste("ID_", mytree$tip.label, sep = '') 
  write.tree(mytree, file = phylo_name, digits = 12)
  #command_tree = paste("sudo python scripts/phylogeny_distance.py --lmm pyseerData/phylo", i, ".newick > pyseerData/MDS_K", i, ".tsv", sep = '')
  #system(command_tree)
}

write("sudo python scripts/phylogeny_distance.py --lmm pyseerData/phylo1.newick > pyseerData/MDS_K1.tsv", file="k3.sh")
for(i in 2:50){
  command <- paste("sudo python scripts/phylogeny_distance.py --lmm pyseerData/phylo",i,".newick > pyseerData/MDS_K",i,".tsv", sep = "")
  write(command, file="k3.sh", append=T)
}

write("pyseer --lmm --phenotypes pyseerData/y_1.pheno --vcf pyseerData/SNP_1.vcf --similarity pyseerData/MDS_K1.tsv --min-af 0.00 --max-af 1.00 > pyseerData/result1.txt", sep = "", file="lmm.sh")
for(i in 2:50){
  command1 <- paste("pyseer --lmm --phenotypes pyseerData/y_", i, ".pheno --vcf pyseerData/SNP_", i, ".vcf --similarity pyseerData/MDS_K", i, ".tsv --min-af 0.00 --max-af 1.00 > pyseerData/result", i, ".txt", sep = "")
  write(command1, file="lmm.sh", append=T)
}


rm(list=ls())
# simluate phenotype
#library("phyclust")
#library("MASS")
#source("emma")

N = 50
M = 30
P_pyseer <- list()

set.seed(321)

for(i in 1:N){
  readseq_dir <- paste('output_SNP/SNP',i ,'.csv', sep = '')
  seq_i <- read.csv(readseq_dir, sep = ' ', header = FALSE)
  
  mean_mat <- apply(seq_i, 2, mean)
  sd_mat <- sqrt(mean_mat*(1-mean_mat))
  X <- apply(seq_i, 2, scale)
  phe_name <- paste('output_phenotype/y_', i, '.csv', sep = "")
  y <- unlist(read.csv(phe_name, sep = ",", header = FALSE))
  
  L <- ncol(X)
  P_pyseertemp <- rep(NA, L)

  for(l1 in 1:L){
    out <- summary(lm(y~X[,l1]-1))
    # we store the pvalue in minus log form with 10 base
    P_pyseertemp[l1] <- -log10(out$coefficients[,4])
  }
  
  results <- read.table(paste0("pyseerData/result", i, ".txt"), sep = "\t", header = TRUE)
  ppyseer <- -log10(results$filter.pvalue)
  final1 <- gsub("11126_", "", results$variant)
  final_index <- strtoi(gsub("_G_A", "", final1))
  P_pyseertemp[final_index] <- ppyseer
  
  P_pyseer[[i]] <- P_pyseertemp

}
filename_LM = paste('Pvalue/pyseer.csv', sep = "")
write.table(unlist(P_pyseer), file = filename_LM, row.names = FALSE, col.names = FALSE)

#pyseer --lmm --phenotypes resistance.pheno --pres SNP.Rtab --similarity MDS_K1.tsv --save-m mash_mds --max-dimensions 8 > res.txt
#sudo python scripts/phylogeny_distance.py --lmm phylo.newick > MDS_K1.tsv
#pyseer --phenotypes resistance.pheno --pres SNP.Rtab --distances MDS_K1.tsv --save-m mash_mds --max-dimensions 8 > resGWAS.txt
