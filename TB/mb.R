library(stringr)
library("phyclust")
library("seqinr")
#library(phylotate)

mbtree_Dir <- 'mkdir mb_tree'
if (!file.exists('mb_tree')){
  system('mkdir mb_tree')
}

N=1
nuncertanty <- 50
M = 469
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
file <- paste("file.nex", sep = '')
batchname <- "batch.nex"
cat("begin mrbayes;",file=batchname,sep="\n")
cat("set autoclose=yes nowarn=yes;",file=batchname, sep="\n", append=TRUE)
cat(paste("execute", file,";"),file=batchname, sep="\n",append=TRUE)
cat("prset brlenspr=clock:uniform;",file=batchname, sep="\n",append=TRUE)
cat("mcmc nruns=1 nchains=1 ngen=1000000 samplefreq=5000;",file=batchname, sep="\n",append=TRUE)
cat("sumt burnin = 500000;",file=batchname, sep="\n",append=TRUE)
cat("end;",file="batch.nex",append=TRUE)
command_mb = paste("mb batch.nex > log.txt", sep = '')
system(command_mb)


skipLines = (M + 5)
tree_name = paste("file.nex.t", sep = '')
tree.anc <- read.csv(tree_name, skip = skipLines, sep="\t", header = FALSE)
ntree <- nrow(tree.anc)
thinned_index <- (ntree-nuncertanty):(ntree-1)
tree_name2 = paste("file.nex.con.tre", sep = '')
tips <- read.csv(tree_name2, skip = 5, sep="\t", header = FALSE)[1:M,3]
index <- 1
for(i in thinned_index){
  test <- mysub(tree.anc[i,])
  mbtree_name = paste("mb_tree/tree_", index,".newick", sep = '')
  index <- index + 1
  mytree <- read.tree(text = test$tree)
  #mytree$tip.label <- tips[as.numeric(mytree$tip.label)]
  mytree$tip.label <- paste("leaf_", mytree$tip.label, sep = '') 
  write.tree(mytree, file = mbtree_name, digits = 12)
}


