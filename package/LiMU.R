LiMU <- function(nunc, M, readseq_dir, file, phe_name){
  if (!file.exists('mb_tree')){
    system('mkdir mb_tree')
  }
  
  if (!file.exists('output_batch')){
    system('mkdir output_batch')
  }
  
  if (!file.exists('Pvalue')){
    system('mkdir Pvalue')
  }
  
  batchname <- "output_batch/batch.nex"
  command_mb = "mb output_batch/batch.nex > log.txt"
  tree_name = paste(file, ".t", sep = '')
  tree_name2 = paste(file, ".con.tre", sep = '')
  write.table(nunc, "M.txt", row.names = FALSE, col.names = FALSE)
  
  mysub <- function(test){
    test.new <- test
    index1 <- unlist(gregexpr(pattern="\\[", test,perl = TRUE)) 
    index2 <- unlist(gregexpr(pattern="\\]", test,perl = TRUE)) 
    substr(test,index1[1], index2[1])
    M <- substr(test,1,(index2[1]+1))
    test.new <- str_replace(test.new,fixed(M),"")
    return(list(tree=test.new))
  }
  
  #generate mb files and run mb
  cat("begin mrbayes;",file=batchname,sep="\n")
  cat("set autoclose=yes nowarn=yes;",file=batchname, sep="\n", append=TRUE)
  cat(paste("execute", file,";"),file=batchname, sep="\n",append=TRUE)
  cat("prset brlenspr=clock:uniform;",file=batchname, sep="\n",append=TRUE)
  cat("mcmc nruns=1 nchains=1 ngen=500000 samplefreq=1000;",file=batchname, sep="\n",append=TRUE)
  cat("sumt burnin = 300000;",file=batchname, sep="\n",append=TRUE)
  cat("end;",file=batchname,append=TRUE)
  system(command_mb)
  
  
  skipLines = (M + 5)
  
  tree.anc <- read.csv(tree_name, skip = skipLines, sep="\t", header = FALSE)
  ntree <- nrow(tree.anc)
  thinned_index <- (ntree-nunc):(ntree-1)
  
  tips <- read.csv(tree_name2, skip = 5, sep="\t", header = FALSE)[1:M,3]
  curwd <- getwd()
  setworkd <- paste(curwd, "/src", sep = '')
  write.table(getwd(), "dir.txt", row.names = FALSE, col.names = FALSE)
  index <- 1
  for(i in thinned_index){
    test <- mysub(tree.anc[i,])
    mbtree_name = paste("mb_tree/tree_", index,".newick", sep = '')
    index <- index + 1
    mytree <- read.tree(text = test$tree)
    mytree$tip.label <- tips[as.numeric(mytree$tip.label)]
    write.tree(mytree, file = mbtree_name, digits = 12)
    
  }
  
  setwd(setworkd)
  comp_javac <- "javac -cp ../jars/sfu-legacy-1.0.13.jar:../jars/apache-lucene.jar:../jars/colt-1.2.0.jar/:../jars/commons-math3-3.0.jar/:../jars/google-collections-1.0.jar:../jars/jama-1.0.3.jar:../jars/jblas-1.2.3.jar ./GSM/*.java"  
  javacomand <- "java -cp ../jars/*:. GSM/gsm"
  system(comp_javac)
  system(javacomand)
  
  setwd(curwd)
  
  seq_i <- read.csv(readseq_dir, sep = ' ', header = FALSE)
  mean_mat <- apply(seq_i, 2, mean)
  sd_mat <- sqrt(mean_mat*(1-mean_mat))
  X <- apply(seq_i, 2, scale)
  y <- unlist(read.csv(phe_name, sep = ",", header = FALSE))
  
  L <- ncol(X)
  P_LMMT <- rep(NA, L)
  P_LMMTmat <- matrix(NA, nunc, L)
  
  # LMM with EGSM
  for(index in 1:nunc){
    readK_dir <- paste('mb_tree/K_',index,'.csv', sep = '')
    K_tree <- matrix(unlist(read.csv(readK_dir, sep = ",", header = FALSE)), nr = M)
    rsT <- emma.REML.t(matrix(y, nr = 1), t(seq_i), K_tree)
    P_LMMTmat[index,] <- (rsT$ps)
  }
  
  P_LMMT <- -log10(apply(P_LMMTmat, 2, mean))
  
  write.csv(P_LMMT, "Pvalue/minuslog10P.csv")
  write.csv(P_LMMTmat, "Pvalue/originalP.csv")
}
