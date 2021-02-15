rm(list=ls())
# simluate phenotype
library("phyclust")
library("MASS")

N = 50
M = 300
L = 500
beta = 0.2
sigma_e = 0.5
sigma_b = 0.6

set.seed(321)

index <- sample.int(L, size = N, replace = TRUE)
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
  
  readK_dir <- paste('output_tree/K_',i ,'.csv', sep = '')
  K_tree <- read.csv(readK_dir, sep = ",", header = FALSE)
  
  b <- mvrnorm(n = 1, mu = rep(0, M), Sigma = K_tree)*sigma_b
  eps <- rnorm(1, mean = 0, sd = sigma_e)
  
  y = std_mat_i[,index[i]]*beta + b + eps
  
  filename = paste('output_phenotype/y_', i, '.csv', sep = "")
  write.table(y, file = filename, row.names = FALSE, col.names = FALSE)
}
write.table(index, file = 'index.csv', row.names = FALSE, col.names = FALSE)



