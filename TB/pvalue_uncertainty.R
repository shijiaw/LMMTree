library("ggplot2")

d1 <- read.table("Pvalue/LMMTun.csv") #please note the data read-in for each file takes roughly 2 minutes.
dd <- matrix(d1,nrow=50) 
nn <- ncol(dd)
p.ci <- matrix(NA, nrow=3,ncol=nn)
rownames(p.ci) <- c("lower","CI","upper") # CI is used for the label; actually this column saves the mean of the posterior samples.
for(i in 1:nn){
  ci <- quantile(unlist(dd[,i]),probs=c(0.025,0.5,0.975))
  p.ci[,i] <- ci
}
#plot(x=1:500,y=p.ci[2,],type="l")

df <- as.data.frame(t(p.ci))
col.idx <- rep("black", nn)
col.idx[df$CI>(-log10(0.05/nn))] <- "red"
MAF002 <- read.csv("MAF002.csv")[,-1]
df2 <- cbind(df,Sites = MAF002,col.idx =col.idx )
head(df2)


my_ggplot <- ggplot(df2,aes(x = Sites, y = CI)) +geom_point(color=col.idx)
my_ggplot      
my_ggplot + geom_errorbar(aes(ymin = lower, ymax = upper),color=col.idx)
#}