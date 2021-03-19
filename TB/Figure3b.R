library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(stringr)

#linear regression with MAF = 0.02
p_LM <- read.csv("Pvalue/LM.csv", header = FALSE)[,1]
MAF002 <- read.csv("MAF002.csv")[,-1]
pvalue1_LM <- data.frame(p = p_LM, SNP_Index = MAF002, method = 'LM' )

#LMM with empirical GSM
#index_LMM <- MAF002
p_LMM <- read.csv("Pvalue/LMM.csv", header = FALSE)[,1]
data_LMM <- data.frame(p = p_LMM, SNP_Index = MAF002, method = "LMM")

#TreeLMM
p_tLMM <- read.csv("Pvalue/LMMT.csv", header = FALSE)[,1]
data_tLMM <- data.frame(p = p_tLMM, SNP_Index = MAF002, method = "TreeLMM")


data <- rbind(pvalue1_LM, data_LMM, data_tLMM)

data$sizes = as.numeric(as.factor(data$method))
data$sizes[data$sizes == 1] = 2.1
data$sizes[data$sizes == 2] = 1.5
data$sizes[data$sizes == 3] = 1.5
g1 <- ggplot(data, aes( x=SNP_Index, y=p, color = method, type = method, shape = method)) +
  geom_point(size = data$sizes) +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05/length(MAF002)), linetype="dashed", color = "#595959") +
  xlab("LM") +
  ylab("") +
  scale_color_manual(values = c('#595959', 'blue', 'red')) +
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

print(g1)
gname = c("TBpvalue.pdf",sep="")  
pdf(gname,width=7,height=3, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
ggarrange(g1,
          ncol = 1, nrow = 1, common.legend = TRUE)

print(g1)
dev.off()
