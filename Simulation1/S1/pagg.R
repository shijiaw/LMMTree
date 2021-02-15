rm(list = ls())
library(ggplot2)
library(Rmisc) 
library(Sleuth2)
library(ggpubr)
library(stringr)
library(xtable)

N = 50
L = 500

pLMc = c()
pLMMc = c()
pLMMTc = c()
statusc = c()

#for(i in 1:rep){
  SNPdir <- paste("index.csv",sep = '')
  pvaluedir <- paste("Pvalue", "/LM.csv",sep = '')
  pLMMdir <- paste("Pvalue", "/LMM.csv",sep = '')
  pLMMTdir <- paste("Pvalue", "/LMMT.csv",sep = '')
  
  SNP <- as.vector(unlist(read.csv(SNPdir, sep = " ", header = FALSE)))
  pLM <- as.matrix(read.csv(pvaluedir, sep = " ", header = FALSE))
  pLMM <- as.matrix(read.csv(pLMMdir, sep = " ", header = FALSE))
  pLMMT <- as.matrix(read.csv(pLMMTdir, sep = " ", header = FALSE))
  
  status_mat <- matrix(0, nr = N, nc = L)
  for(i in 1:N){
    status_mat[i, SNP[i]] <- 1
  }
  
  # pLMM <- as.matrix(read.csv(pLMMdir, sep = ",", header = FALSE))
  # pLMMF <- as.matrix(read.csv(pLMMFdir, sep = ",", header = FALSE))
  # 
  # pLM <- matrix(NA, nr = nrow(pLMMF), nc = ncol(pLMMF))
  # pLMM <- pLMM[(1:nrow(pLMMF)),]
  # for(i in 1: nrow(pLMMF)){
  #   pLM[i,] <- unlist(pvalue[i,])[unlist(SNP[i,]+1)]
  # }
  
  status1 <- as.vector(status_mat)
  pLM1 <- as.vector(pLM)
  pLMM1 <- as.vector(pLMM)
  pLMMT1 <- as.vector(pLMMT)
  
  # pLMc = c(pLMc, pLM1)
  # pLMMc = c(pLMMc, pLMM1)
  # pLMMFc = c(pLMMFc, pLMMF1)
  #statusc = c(statusc, status1)
  
#}



pvData <- data.frame(pLM = pLM1, pLMMF =  pLMMT1 , pLMM = pLMM1, status = as.numeric(status1))
#write.csv(pvData, "p1.csv")
#colsum <- apply(pvData, 1, sum)
#index <- which(colsum > 0)
#pvData <- pvData[index, ]
#apply(pvData, 2, sum)

library(pROC)
f = pvData
model.lm = roc(response = f$status, predictor = 1 - exp(-log(10) * f$pLM))
model.tr = roc(response = f$status, predictor = 1 - exp(-log(10) * f$pLMMF))
model.lmm = roc(response = f$status, predictor = 1 - exp(-log(10) * f$pLMM))
#N = dim(f)[1]
df = data.frame(
  TPR = c(model.lm$sens, model.lmm$sens, model.tr$sens),
  FPR = 1 - c(model.lm$spec, model.lmm$spec, model.tr$spec),
  model = c(rep("LM", length(model.lm$sens)), rep("LMM", length(model.lmm$sens)), rep("TR", length(model.tr$sens)))
)

g = ggplot(df) + geom_line(aes(x = FPR, y = TPR, colour = model)) +
  labs(title = 'N = 300, M = 500') + theme_minimal() +
  xlim(c(0,1)) + ylim(c(0,1)) + geom_vline(xintercept=c(0.05), linetype="dotted")

gname = c("LMMROC1.eps",sep="")  
postscript(gname,width=4,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
g
dev.off()


lmind <- which.min(abs(1- model.lm$spec-0.05))
lmTPR <- model.lm$sens[lmind]

lmmind <- which.min(abs(1- model.lmm$spec-0.05))
lmmTPR <- model.lmm$sens[lmmind]

lmmTind <- which.min(abs(1- model.tr$spec-0.05))
lmmTTPR <- model.tr$sens[lmmTind]

print(c(lmTPR, lmmTPR, lmmTTPR))







