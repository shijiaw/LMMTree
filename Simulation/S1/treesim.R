library("phyclust")
library(phangorn)
# simulate N datasets, each with one tree, M subjects, 
# simulate N datasets, each with sequences for M samples, length L (mutations),


N = 50
M = 30
L = 500


tree_Dir <- 'mkdir output_tree'
if (!file.exists('output_tree')){
  system('mkdir output_tree')
}

# store ACGT
seq_Dir <- 'mkdir output_seq'
if (!file.exists('output_seq')){
  system('mkdir output_seq')
}

# store 0, 1 SNPs
seq_Dir <- 'mkdir output_SNP'
if (!file.exists('output_SNP')){
  system('mkdir output_SNP')
}


DNA2SNP <- function(matrixDNA){
  # Find the major allels
  nc = ncol(matrixDNA)
  ma <- rep(NA, nc)
  for(i in 1:nc){
    ma[i] <- names(which.max(table(matrixDNA[,i])))
  }
  
  SNPmatrix <- matrix(NA, nr = nrow(matrixDNA), nc = nc)
  for(i in 1:nc){
    SNPmatrix[ ,i] <- 1-as.numeric(matrixDNA[,i] == ma[i])
  }
  
  index <- which(colMeans(SNPmatrix) > 0)
  
  output <- SNPmatrix[,index]
  return(output)
}


for(n in 1:N){
  command_tree = paste('./ms ',M,' 1 -T -seeds ',n,' 10 10 | tail +4 | grep -v //>output_tree/tree',n,'.newick', sep = '')
  system(command_tree)
  tree_name = paste("output_tree/tree", n, ".newick", sep = '')
  tree.anc <- read.tree(tree_name)
  leaf_order <- as.numeric(tree.anc$tip.label)
  tree.anc$tip.label <- paste("leaf_", tree.anc$tip.label, sep = "")
  write.tree(tree.anc, file = tree_name, digits = 12)
  seq <- simSeq(tree.anc, l = L, type="DNA", rate = 3)
  seq_name = paste("output_seq/seq", n, ".nex", sep = '')
  write.nexus.data(toupper(seq), file = seq_name)
  
  SNP <- DNA2SNP(as.character(seq[leaf_order,]))
  SNP_name = paste("output_SNP/SNP", n, ".csv", sep = '')
  write.table(SNP, file = SNP_name, row.names = FALSE, col.names = FALSE)
  
  #command_seq = paste('./ms ',M,' 1  -s ', L, ' -seeds ',n,' 10 10 | tail +7 | grep -v //>output_seq/seq',n,'.txt', sep = '')
  #system(command_seq)
}

