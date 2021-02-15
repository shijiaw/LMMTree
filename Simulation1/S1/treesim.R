library("phyclust")
# simulate N datasets, each with one tree, M subjects, 
# simulate N datasets, each with sequences for M samples, length L (mutations),


N = 50
M = 300
L = 500


tree_Dir <- 'mkdir output_tree'
if (!file.exists('output_tree')){
  system('mkdir output_tree')
}

seq_Dir <- 'mkdir output_seq'
if (!file.exists('output_seq')){
  system('mkdir output_seq')
}

for(n in 1:N){
  command_tree = paste('./ms ',M,' 1 -T -seeds ',n,' 10 10 | tail +4 | grep -v //>output_tree/tree',n,'.newick', sep = '')
  system(command_tree)
  command_seq = paste('./ms ',M,' 1  -s ', L, ' -seeds ',n,' 10 10 | tail +7 | grep -v //>output_seq/seq',n,'.txt', sep = '')
  system(command_seq)
}

#readChar(file = 'output_seq/seq1.txt', header = FALSE)
for(i in 1:N){
  readseq_dir <- paste('output_seq/seq',i ,'.txt', sep = '')
  seq_i <- do.call(rbind,lapply(strsplit(scan(readseq_dir, 
                                           what=""), ""), as.numeric))
  mean_mat <- apply(seq_i, 2, mean)
  sd_mat <- sqrt(mean_mat*(1-mean_mat))
  std_mat_i <- apply(seq_i, 2, scale)
  for(j in 1:L){
    std_mat_i[,j] <- (seq_i[,j]-mean_mat[j])/sd_mat[j]
  }
  K_i <- std_mat_i %*% t(std_mat_i)/L
  
  filename = paste('output_seq/K_', i, '.csv', sep = "")
  write.table(K_i, file = filename, row.names = FALSE, col.names = FALSE)
  
}  

# change tip name of trees
for(i in 1:N){
  tree_name = paste("output_tree/tree", i, ".newick", sep = '')
  tree.anc <- read.tree(tree_name)
  tree.anc$tip.label <- paste("leaf_", tree.anc$tip.label, sep = "")
  write.tree(tree.anc, file = tree_name, digits = 12)
}

reftree_name = "output_tree/reftree.newick"
tree.anc <- read.tree(tree_name)
tree.anc$tip.label <- paste("leaf_", 1:M, sep = "")
write.tree(tree.anc, file = reftree_name, digits = 12)
