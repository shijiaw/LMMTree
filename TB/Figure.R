library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(stringr)
library(seqinr)
library(dplyr)
#ncrna <- read.fasta(file = "1Nov.fa")

snps <- read.fasta(file = "out475_snps.fa")

nc = 20976
matrixDNA_temp <- matrix(NA, nr = 471, nc = nc)
for(i in 1:471){
  matrixDNA_temp[i, ] <- snps$org.code[i,]
}
commonSNP_Index <- read.csv("common_Index_fa.csv")[,-1]
matrixDNA <- matrixDNA_temp[commonSNP_Index, ]

# Find the major allels
ma <- rep(NA, nc)
for(i in 1:nc){
  ma[i] <- names(which.max(table(matrixDNA[,i])))
}
NAindex <- which(ma == 'N')



load("index_to_reference.Rdata")
alrefdf2 <- alrefdf2[-NAindex,]
singular <- read.csv("singular.csv")[,-1]
alrefdf2 <- alrefdf2[singular,]
MAF005 <- read.csv("MAF005.csv")[,-1]
BP <- alrefdf2[MAF005,2]


#linear regression with MAF = 0.02
p_LM <- read.csv("Pvalue/LM.csv", header = FALSE)[,1]
MAF005 <- read.csv("MAF005.csv")[,-1]
pvalue1_LM <- data.frame(p = p_LM, SNP_Index = BP, method = 'LM' )

#LMM with empirical GSM
#index_LMM <- MAF002
p_LMM <- read.csv("Pvalue/LMM.csv", header = FALSE)[,1]
data_LMM <- data.frame(p = p_LMM, SNP_Index = BP, method = "LMM")

#TreeLMM
p_tLMM <- read.csv("Pvalue/LMMT.csv", header = FALSE)[,1]
data_tLMM <- data.frame(p = p_tLMM, SNP_Index = BP, method = "LiMU")

#pyseer
results <- read.table("res.txt", sep = "\t", header = TRUE)
p_pyseer <- -log10(results$filter.pvalue)
data_pyseer <- data.frame(p = p_pyseer, SNP_Index = BP, method = "PYSEER")
#final1 <- gsub("11126_", "", results$variant[index])
#final <- strtoi(gsub("_G_A", "", final1))

data <- rbind(pvalue1_LM, data_pyseer, data_LMM, data_tLMM)

data$sizes = as.numeric(as.factor(data$method))
#data$sizes = rep(1:4, each = length(MAF002))
data$sizes[data$sizes == 1] = 1.1
data$sizes[data$sizes == 2] = 1.1
data$sizes[data$sizes == 3] = 1.1
data$sizes[data$sizes == 4] = 1.1
#data$shape = rep(c(3,1,2,4), each = length(MAF002))
g1 <- ggplot(data, aes( x=SNP_Index, y=p, color = method, type = method, shape = method)) +
  geom_point(size = data$sizes) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05/length(MAF005)), linetype="dashed", color = "#595959") +
  xlab("LM") +
  ylab("") +
  scale_color_manual(values = c('red', '#595959', 'blue', 'grey')) + scale_shape_manual(values = c(1,2,3,4))+
  xlab("Genome location") +
  ylab('-log10 p-value') +
  ylim(0, 20) + 
  xlim(0, max(data$SNP_Index)) +
  theme(legend.title = element_blank(), axis.title=element_text(size=12),axis.text=element_text(size=12), legend.position = c(0.917,0.87), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18)) +
   theme(legend.title = element_blank(),
        legend.spacing.y = unit(-1, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
   ggtitle('Association studies of multidrug-resistant tuberculosis') +
   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#print(g1)
gname = c("TBpvalue.pdf",sep="")  
pdf(gname,width=7,height=3, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(g1,
          ncol = 1, nrow = 1, common.legend = TRUE)

print(g1)
dev.off()

ref <- data_tLMM[which(data_tLMM[,1] > 6),2]

