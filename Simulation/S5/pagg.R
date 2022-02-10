rm(list = ls())
library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(stringr)
library(xtable)

AUC0 <- read.csv("0/AUC0.csv")[,-1]
AUC1 <- read.csv("1/AUC0.csv")[,-1]
AUC2 <- read.csv("2/AUC0.csv")[,-1]
AUC3 <- read.csv("3/AUC0.csv")[,-1]
AUC4 <- read.csv("4/AUC0.csv")[,-1]
AUC5 <- read.csv("5/AUC0.csv")[,-1]
AUC6 <- read.csv("6/AUC0.csv")[,-1]
AUC7 <- read.csv("7/AUC0.csv")[,-1]
AUC8 <- read.csv("8/AUC0.csv")[,-1]
AUC9 <- read.csv("9/AUC0.csv")[,-1]
AUC10 <- read.csv("10/AUC0.csv")[,-1]
AUC <- c(AUC0, AUC1, AUC2, AUC3, AUC4, AUC5, AUC6, AUC7, AUC8, AUC9, AUC10)

df <- data.frame(method = rep(c("LM", "LMM", "LiMU", "pyseer"), 11),
                  sigma = rep(seq(0, 1, 0.1), each = 4),
                  AUC = AUC)
head(df)

p <- ggplot(data=df, aes(x=sigma, y=AUC, group=method)) +
  geom_line(linetype="dashed", aes(color=method), size=1.2)+
  geom_point(aes(color=method), size=3) + theme_classic() +
  labs(title="",x=expression(sigma[e]), y = "AUC")



gname = c("plotAUC.eps",sep="")  
postscript(gname,width=8,height=6,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
p
dev.off()



