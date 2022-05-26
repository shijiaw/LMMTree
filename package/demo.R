rm(list=ls())
library(stringr)
library(phyclust)
library(MASS)
source("emma.R")
source("LiMU.R")

nunc <- 100
M <- 15
readseq_dir <- 'DataDemo/SNP1.csv'
file <- "DataDemo/seq1.nex"
phe_name <- 'DataDemo/y_1.csv'

# run LiMU
run <- LiMU(nunc, M, readseq_dir, file, phe_name)


