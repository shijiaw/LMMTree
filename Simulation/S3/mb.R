library(stringr)
library("phyclust")
#library(phylotate)

mbtree_Dir <- 'mkdir mb_tree'
if (!file.exists('mb_tree')){
  system('mkdir mb_tree')
}

nuncertanty <- 100
M = 100
mysub <- function(test){
  test.new <- test
  index1 <- unlist(gregexpr(pattern="\\[", test,perl = TRUE)) 
  index2 <- unlist(gregexpr(pattern="\\]", test,perl = TRUE)) 
  substr(test,index1[1], index2[1])
  M <- substr(test,1,(index2[1]+1))
  test.new <- str_replace(test.new,fixed(M),"")
  #test.new <- str_replace(test.new,";" ,"")
  # for(i in 2:length(index1)){
  #   M <- substr(test, index1[i],index2[i])
  #   test.new <- str_replace(test.new,fixed(M),"")
  # }
  return(list(tree=test.new))
}


#generate mb files and run mb
for(n in 1:N){
  file <- paste("output_seq/seq_error", n, ".nex", sep = '')
  batchname <- paste("output_seq/batch", n, ".nex", sep = '')
  cat("begin mrbayes;",file=batchname,sep="\n")
  cat("set autoclose=yes nowarn=yes;",file=batchname, sep="\n", append=TRUE)
  cat(paste("execute", file,";"),file=batchname, sep="\n",append=TRUE)
  cat("prset brlenspr=clock:uniform;",file=batchname, sep="\n",append=TRUE)
  cat("mcmc nruns=1 nchains=1 ngen=500000 samplefreq=1000;",file=batchname, sep="\n",append=TRUE)
  cat("sumt burnin = 300000;",file=batchname, sep="\n",append=TRUE)
  cat("end;",file=batchname,append=TRUE)
  command_mb = paste("mb output_seq/batch", n, ".nex > log.txt", sep = '')
  system(command_mb)
  
  
  skipLines = (M + 5)
  tree_name = paste("output_seq/seq_error", n, ".nex.t", sep = '')
  tree.anc <- read.csv(tree_name, skip = skipLines, sep="\t", header = FALSE)
  ntree <- nrow(tree.anc)
  thinned_index <- (ntree-nuncertanty):(ntree-1)
  tree_name2 = paste("output_seq/seq_error", n, ".nex.con.tre", sep = '')
  tips <- read.csv(tree_name2, skip = 5, sep="\t", header = FALSE)[1:M,3]
  index <- 1
  for(i in thinned_index){
    test <- mysub(tree.anc[i,])
    mbtree_name = paste("mb_tree/tree", n, "_", index,".newick", sep = '')
    index <- index + 1
    mytree <- read.tree(text = test$tree)
    mytree$tip.label <- tips[as.numeric(mytree$tip.label)]
    write.tree(mytree, file = mbtree_name, digits = 12)
  }
  print(n)

}


