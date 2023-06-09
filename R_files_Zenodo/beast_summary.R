library(ape)
library(xlsx)
library(treeio)

# set working directory 
setwd("/path/to/BEAST/output/")


# load trees with phyloch package
phyClock<-phyloch::read.beast("clock_out_MCC.nex")
phyPIS<-phyloch::read.beast("phylo_out_MCC.nex")
phyR1<-phyloch::read.beast("r1_out_MCC.nex")
phyR2<-phyloch::read.beast("r2_out_MCC.nex")
phyR3<-phyloch::read.beast("r3_out_MCC.nex")
phyR4<-phyloch::read.beast("r4_out_MCC.nex")
phyR5<-phyloch::read.beast("r5_out_MCC.nex")
phyR6<-phyloch::read.beast("r6_out_MCC.nex")
phyR7<-phyloch::read.beast("r7_out_MCC.nex")
phyR8<-phyloch::read.beast("r8_out_MCC.nex")


# get node ages from each beast tree
Clock_Age<-phyClock$height
PIS_Age<-phyPIS$height
R1_Age<-phyR1$height
R2_Age<-phyR2$height
R3_Age<-phyR3$height
R4_Age<-phyR4$height
R5_Age<-phyR5$height
R6_Age<-phyR6$height
R7_Age<-phyR7$height
R8_Age<-phyR8$height


# put all the ages in a dataframe
t1<-data.frame(c(410:817), 
               Clock_Age, PIS_Age, R1_Age, R2_Age, R3_Age, R4_Age, R5_Age, R6_Age, R7_Age, R8_Age) 


names(t1)<-c("Node", "Clock_Age", "PIS_Age", "R1_Age", "R2_Age", "R3_Age", "R4_Age", "R5_Age", "R6_Age", "R7_Age", "R8_Age")


# calculate the mean age of each node from all the trees and add it to dataframe 
Mean_Age<-(Clock_Age + PIS_Age + R1_Age + R2_Age + R3_Age + R4_Age + R5_Age + R6_Age + R7_Age + R8_Age)/10

t1$Mean<-Mean_Age


# dataframe with only random trees
t2<-data.frame(c(410:817), R1_Age, R2_Age, R3_Age, R4_Age, R5_Age, R6_Age, R7_Age, R8_Age) 
names(t2)<-c("Node", "R1_Age", "R2_Age", "R3_Age", "R4_Age", "R5_Age", "R6_Age", "R7_Age", "R8_Age")
Mean_Age<-(R1_Age + R2_Age + R3_Age + R4_Age + R5_Age + R6_Age + R7_Age + R8_Age)/8
t2$Mean<-Mean_Age


# create consensus tree with mean node heights
phyMean<-phyClock
phyMean$height<-(phyClock$height + phyPIS$height + phyR1$height + phyR2$height + phyR3$height + phyR4$height + phyR5$height + phyR6$height + phyR7$height + phyR8$height)/10
phyMean$length<-(phyClock$length + phyPIS$length + phyR1$length + phyR2$length + phyR3$length + phyR4$length + phyR5$length + phyR6$length + phyR7$length + phyR8$length)/10
phyMean$edge.length<-(phyClock$edge.length + phyPIS$edge.length + phyR1$edge.length + phyR2$edge.length + phyR3$edge.length + phyR5$edge.length + phyR6$edge.length + phyR4$edge.length + phyR7$edge.length + phyR8$edge.length)/10


### option: plot all trees in one plot with consensus tree on top all others as bold & black 

# plot(phyR1, edge.color = "skyblue1", type="c", show.tip.label=FALSE)
# par(new=TRUE)
# plot(phyR2, edge.color = "skyblue2", type="c", show.tip.label=FALSE)
# par(new=TRUE)
# plot(phyR3, edge.color = "skyblue3", type="c", show.tip.label=FALSE)
# par(new=TRUE)
# plot(phyR4, edge.color = "skyblue4", type="c", show.tip.label=FALSE)
# par(new=TRUE)
# plot(phyR5, edge.color = "lightskyblue1", type="c", show.tip.label=FALSE)
# par(new=TRUE)
# plot(phyR6, edge.color = "lightskyblue2", type="c", show.tip.label=FALSE)
# par(new=TRUE)
# plot(phyR7, edge.color = "lightskyblue3", type="c", show.tip.label=FALSE)
# par(new=TRUE)
# plot(phyR8, edge.color = "lightskyblue4", type="c", show.tip.label=FALSE)
# par(new=TRUE)
# plot(phyClock, edge.color = "chocolate4", type="c", show.tip.label=FALSE)
# par(new=TRUE)
# plot(phyPIS, edge.color = "chartreuse4", type="c", show.tip.label=FALSE)
# par(new=TRUE)
# plot(phyMean, edge.color = "black", type="c", show.tip.label=FALSE)



# export files 

write.tree(phyClock, file = "out_clock_20p_MCC.tre")
write.tree(phyPIS, file = "out_phylo_20p_MCC.tre")
write.tree(phyR1, file = "out1_20p_MCC.tre")
write.tree(phyR2, file = "out2_20p_MCC.tre")
write.tree(phyR3, file = "out3_20p_MCC.tre")
write.tree(phyR4, file = "out4_20p_MCC.tre")
write.tree(phyR5, file = "out5_20p_MCC.tre")
write.tree(phyR6, file = "out6_20p_MCC.tre")
write.tree(phyR7, file = "out7_20p_MCC.tre")
write.tree(phyR8, file = "out8_20p_MCC.tre")
write.tree(phyMean, file = "out_mean_20p_MCC.tre")

write.xlsx(t1, file = "BEAST_summary.xlsx")
write.csv(t1, file = "BEAST_summary.csv")




