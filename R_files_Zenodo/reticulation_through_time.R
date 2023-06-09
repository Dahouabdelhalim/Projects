#Written by: Dylan DeBaun

library(ape)
library(phytools)
library(ggplot2)
tree<-read.tree()
tree1<-force.ultrametric(tree)
tree1<-drop.tip(tree1,c(1:14,124:130))
x<-ltt(tree1,log.lineages = F)

bt<-branching.times(tree1)

ltt1 <- as.data.frame(cbind(x[["times"]]-22.17059,log(x[["ltt"]])))
ltt15 <- format(ltt1[-c(1:11),],scientific=F)
rtt <- read.csv("~/Desktop/rtt1.csv")
ltt15[1,] = as.numeric(c(-14.71835,2.639057))
ltt15[97,1] = as.numeric(0)
rtt[13,] = as.numeric(c(0,109,12,0,0,0))

#GEOM STEP Plot 
p <- ggplot()+geom_step(aes(x=as.numeric(ltt15[,1]), y = as.numeric(ltt15[,2])-as.numeric(min(ltt15[,2]))), size= 1.5,direction="vh")+  geom_step(aes(x=-rtt$year..my.,y=log(rtt$num_retic)),color="dark blue",size=1.5,direction="hv")+ scale_y_continuous(name="ln(Number of Lineages)",breaks = seq(0,2.5,0.5), labels = seq(2.56,5.06,0.5),sec.axis = sec_axis(trans=~.*1, name="ln(Number of Reticulations)"))
p
#p <- ggplot() + geom_point(aes(x=rtt$year.my.,y=rtt$log.num.),color="dark blue",size=2) +scale_x_continuous(breaks = seq(-15,0,5),name = "Time (mya)")              
p <- p + geom_vline(xintercept=-2.588) + geom_vline(xintercept=-5.332)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size=14),panel.border = element_rect(colour = "black", fill=NA, size=0.5))
p
p<- p +  geom_errorbarh(data = rtt, aes(xmin=-positive,xmax = -negative, y= log(num_retic)), alpha=0.4,linetype = "dashed",color="dark blue") +theme(legend.position = "none")  +xlab("Time (mya)") +scale_x_continuous(limits=c(-max(rtt$positive),0), expand = c(0.01, 0.02)) 
p
pdf("~/Desktop/Figure_LTTRTTLOGSCALE_VH.pdf",height = 10, width = 12)
plot(p)
dev.off()

#supplement figure
p <- ggplot()+  geom_step(aes(x=-rtt$year..my.,y=log(rtt$num_retic)),color="dark blue",size=1.5,direction="hv")+ scale_y_continuous(name="ln(Number of Lineages)",breaks = seq(0,2.5,0.5), labels = seq(2.56,5.06,0.5),sec.axis = sec_axis(trans=~.*1, name="ln(Number of Reticulations)"))
p <- p +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"),text = element_text(size=14),panel.border = element_rect(colour = "black", fill=NA, size=0.5))
p<- p +  geom_errorbarh(data = rtt, aes(xmin=-positive,xmax = -negative, y= log(num_retic)), alpha=0.4,linetype = "dashed",color="dark blue") +theme(legend.position = "none")  +xlab("Time (mya)") +scale_x_continuous(limits=c(-max(rtt$positive),0), expand = c(0.01, 0.02)) 
p
p + geom_line(aes(x=c(-rtt$year..my.[1],0),y=c(log(1),log(12))),col="red")
#calculation for the rate
log(12)/rtt$year..my.[1]
#0.1688305
