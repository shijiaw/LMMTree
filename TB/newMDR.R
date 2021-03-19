# new MDR
#library(Biostrings)
library(ggplot2)
library(ggpubr)

#s = readDNAStringSet("1Nov.fa")
snps <- read.fasta(file = "out475_snps.fa")

data <- read.csv("combined_metadata_470_strains_24032016.csv", header = TRUE)

commonnames <- intersect(names(snps),data$tip_labels2)

common_Index <- which(data$tip_labels2%in%commonnames)

common_Index_fa <- which(names(snps)%in%commonnames)

write.csv(common_Index_fa, "common_Index_fa.csv")

N = length(names(snps))

Y <- 1-as.numeric(data$drug_susceptibility_result == 'MDR')
stdY <- (Y-mean(Y))/sd(Y)
finalY <- stdY[common_Index]
write.csv(finalY, "MDR.csv")
