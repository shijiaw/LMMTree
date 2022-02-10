#plot
library("ggplot2")
index <- unlist(read.table("index.csv"))
setwd("./PvalueUncertainty")
list1 <- c()
list2 <- c()
for(j in 1:50){
  fn <- paste("LMMTmat_",j,".csv",sep = '')
  d1 <- as.matrix(read.csv(fn, sep = " ", header = FALSE))
  dd <- matrix(d1,nrow=100)
  list1 <- c(list1, dd[,index[j]])
  list2 <- c(list2, as.vector(dd[,-index[j]]))
  print(j)
}

pvalue1 <- data.frame(pv = list1, type = rep('effective', length(list1)))
pvalue2 <- data.frame(pv = list2, type = rep('non-effective', length(list2)))
pvalue = rbind(pvalue1, pvalue2)

p<-ggplot(pvalue, aes(x=type, y=pv, color=type)) +
  geom_boxplot()+ scale_fill_grey() + theme_classic()+theme(axis.title.x=element_blank(),
                                                            axis.text.x=element_blank(),
                                                            axis.ticks.x=element_blank()) +labs(y = expression(-log[10](P-value)))

gname = c("punc.eps",sep="")  
postscript(gname,width=5,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
p
dev.off()


#}