#plot
library("ggplot2")
setwd("./PvalueUncertainty")
#d1 <- read.table("LMMTmat_1csv")
j <- 1  # replicate 1
fn <- paste("LMMTmat_",j,".csv",sep = '')

d1 <- as.matrix(read.csv(fn, sep = " ", header = FALSE))

#d1 <- read.table(fn) #please note the data read-in for each file takes roughly 2 minutes.
dd <- matrix(d1,nrow=50) #2000 sites; 50 posterior samples. -log_10(p)
p.ci <- matrix(NA, nrow=3,ncol=2000)
rownames(p.ci) <- c("lower","CI","upper") # CI is used for the label; actually this column saves the mean of the posterior samples.
for(i in 1:2000){
  ci <- quantile(unlist(dd[,i]),probs=c(0.025,0.5,0.95))
  p.ci[,i] <- ci
}
plot(x=1:2000,y=p.ci[2,],type="l")

df <- as.data.frame(t(p.ci))
col.idx <- rep("black",2000)
col.idx[df$CI>(-log10(0.05/2000))] <- "red"
df2 <- cbind(df,Sites=c(1:2000),col.idx =col.idx )
head(df2)


my_ggplot <- ggplot(df2,aes(x = Sites, y = CI)) +geom_point(color=col.idx)
my_ggplot      


gname = c("punc2.eps",sep="")  
postscript(gname,width=5,height=3,horizontal = FALSE, onefile = FALSE, paper = "special")
par(mfrow=c(1,1),oma=c(0.2,1.5,0.2,1.5),mar=c(3,2,0.2,2),cex.axis=1,las=1,mgp=c(1,0.5,0),adj=0.5)
my_ggplot + geom_errorbar(aes(ymin = lower, ymax = upper),color=col.idx)
dev.off()


#}