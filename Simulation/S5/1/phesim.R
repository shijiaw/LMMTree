rm(list=ls())
# simluate phenotype
library("phyclust")
library("MASS")

N = 50
M = 15
L = 100
beta = 0.1
sigma_e = 0.1
sigma_b = 0.1

set.seed(321)

phe_Dir <- 'mkdir output_phenotype'
if (!file.exists('output_phenotype')){
  system('mkdir output_phenotype')
}

index <- rep(NA, N)
for(i in 1:N){
  readseq_dir <- paste('output_SNP/SNP',i ,'.csv', sep = '')
  seq_i <- read.csv(readseq_dir, sep = ' ', header = FALSE)
  
  mean_mat <- apply(seq_i, 2, mean)
  sd_mat <- sqrt(mean_mat*(1-mean_mat))
  std_mat_i <- apply(seq_i, 2, scale)
  
  index[i] <- sample.int(ncol(seq_i), size = 1)
  
  readK_dir <- paste('output_tree/K_',i ,'.csv', sep = '')
  K_tree <- read.csv(readK_dir, sep = ",", header = FALSE)
  
  b <- mvrnorm(n = 1, mu = rep(0, M), Sigma = K_tree)*sigma_b
  eps <- rnorm(M, mean = 0, sd = sigma_e)
  
  y = std_mat_i[,index[i]]*beta + b + eps
  
  filename = paste('output_phenotype/y_', i, '.csv', sep = "")
  write.table(y, file = filename, row.names = FALSE, col.names = FALSE)
}
write.table(index, file = 'index.csv', row.names = FALSE, col.names = FALSE)



