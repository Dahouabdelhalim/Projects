# Data Analysis and Plots for Rashidi and Ostrowski
#########

library(tidyverse)
library(cowplot)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

data = read.table('bacterial_migration.txt', header=TRUE)
d=as_tibble(data)
attach(d)

gram=data.frame(cbind(c('Pf','Kp','Bu','Ec','Ml','Bs','Ef','Sa'),c('N', 'N', 'N','N', 'P', 'P', 'P', 'P')))
names(gram)=c('species', 'gram_status')

# Take average across replicates
e=merge(d,gram, by.x="Against", by.y="species")
names(e)[7]="Gram_Ag"
means = e %>% group_by(Strain, Bacteria, Against, Gram, Gram_Ag) %>% summarise(mean=mean(NumCells), se=sd(NumCells)/sqrt(length(NumCells)))
means=mutate(means,bact_ag=paste(Bacteria, "(",Against, ")",sep=''))  #Create a label for Bacteria(Alternative Choice)

# Pair up the reciprocal bacterial pairs. Assign to an N(P) combo--e.g, Kp(Bs) and Bs(Kp) will both be listed as Kp(Bs)
# with this method, so that they get plotted together.
means2=mutate(means[means$Gram=='P',],bpairs=paste(Against, "  ",Bacteria, " ",sep=''))
means3=mutate(means[means$Gram=='N',],bpairs=paste(Bacteria, "  ",Against, " ",sep=''))
means4=rbind(means2,means3)

# Separate by Strains
qs32_paired=filter(means4,Strain=='QS32')
qs39_paired=filter(means4,Strain=='QS39')
qs40_paired=filter(means4,Strain=='QS40')

#########  MAKE PLOT SHOWING PAIRWISE OUTCOMES, BY STRAIN (=Fig. 1B)
my.cols=c('black', 'gray73')
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
limits <- aes(ymax = qs32_paired$mean + qs32_paired$se, ymin = qs32_paired$mean - qs32_paired$se)
p1 <- ggplot(data=qs32_paired,aes(x=reorder(bpairs,-qs32_paired$mean),mean,fill=Gram,colour=factor(Gram), alpha=0.7,width=0.5)) + 
  geom_bar(stat='identity',position=position_dodge(0.5),width=0.7) +
  geom_errorbar(limits, position=position_dodge(0.5),width=0.2,col='black') +
  xlab(NULL) + ylab("Number of \\nMigrating Amoebae") + ggtitle('QS32 (Texas)') +
  scale_colour_manual(qs32_paired$Gram, values=my.cols) + scale_fill_manual(values=cbPalette[c(4,1)])+ theme(plot.title=element_text(size=20)) +
  ylim(0,300)

limits <- aes(ymax = qs39_paired$mean + qs39_paired$se, ymin = qs39_paired$mean - qs39_paired$se)
p2 <- ggplot(data=qs39_paired,aes(x=reorder(bpairs,-qs32_paired$mean),mean,fill=Gram,colour=factor(Gram), alpha=0.7,width=0.5)) + 
  geom_bar(stat='identity',position=position_dodge(0.5),width=0.7) +
  geom_errorbar(limits, position=position_dodge(0.5),width=0.2,col='black') +
  xlab(NULL) + ylab("Number of \\nMigrating Amoebae") + ggtitle('QS39 (Tennessee)') +
  scale_colour_manual(qs39_paired$Gram, values=my.cols) + scale_fill_manual(values=cbPalette[c(2,1)])+ theme(plot.title=element_text(size=20)) +
  ylim(0,300)

limits <- aes(ymax = qs40_paired$mean + qs40_paired$se, ymin = qs40_paired$mean - qs40_paired$se)
p3 <- ggplot(data=qs40_paired,aes(x=reorder(bpairs,-qs32_paired$mean),mean,fill=Gram,colour=factor(Gram), alpha=0.7,width=0.5)) + geom_bar(stat='identity',position=position_dodge(0.5),width=0.7) +
  geom_errorbar(limits, position=position_dodge(0.5),width=0.2,col='black') +
  xlab('Gram-/Gram+ Choice') + ylab("Number of \\nMigrating Amoebae") + ggtitle('QS40 (Massachusetts)') +
  scale_colour_manual(qs40_paired$Gram, values=my.cols) + scale_fill_manual(values=cbPalette[c(5,1)])+ theme(plot.title=element_text(size=20)) +
  ylim(0,300) +

  guides(fill = FALSE) +
  guides(size = FALSE) +
  guides(color = FALSE)

outplot <- plot_grid(p1,p2,p3, nrow=3,ncol=1,labels=c("A", "B", "C"), label_size=24)
save_plot("Figure_1B.pdf", outplot, base_height=10, base_aspect_ratio=0.9)



##### MAKE PLOT WITH AVERAGE: Across all three strains (=Fig.S2) or across G+ vs G- (=Fig. 1C)
qs32=filter(means,Strain=='QS32')
qs39=filter(means,Strain=='QS39')
qs40=filter(means,Strain=='QS40')

# Plotting parameters
my.cols=c('black', 'gray90')
my.lwd=1

# Take means across bacteria[ag], then across bacteria, then across G+ vs G-
mpbag = means %>% group_by(bact_ag,Bacteria,Gram) %>% summarise(mean_per_bact_ag=mean(mean))# mbpag=means per bacteria(against) combo
mpb= mpbag %>% group_by(Bacteria,Gram) %>% summarise(mean_per_bact=mean(mean_per_bact_ag),se=sd(mean_per_bact_ag)/sqrt(length(mean_per_bact_ag))) 
mpg = mpb %>% group_by(Gram) %>% summarise(mean_per_gram=mean(mean_per_bact),se=sd(mean_per_bact)/sqrt(length(mean_per_bact)))

# Plot means per bacteria, sorted (=Fig. S2)
s<-ggplot(data=mpb) + geom_bar(mapping=aes(x=reorder(Bacteria,-mean_per_bact),y=mean_per_bact,fill=Bacteria,colour=Gram),stat='identity',width=0.7, size=my.lwd,alpha=0.8) +
  geom_errorbar(data=mpb, aes(x=reorder(Bacteria,-mean_per_bact),ymin=mean_per_bact-se, ymax=mean_per_bact+se),width=0.2, position=position_dodge(0.9)) +
  scale_colour_manual("Gram", qs32$Gram, values=my.cols) + scale_fill_manual(values=cbPalette) +
  xlab('Bacteria') + ylab("Mean Number of Migrating Amoebae") +
  labs(fill="Bacteria", colour="Gram")
  save_plot("Figure_S2.pdf",s,base_aspect_ratio=1.3)

# Plot means per Gram status (=Fig. 1C)
  cols=cbPalette[c(2,3)]
t<- ggplot(data=mpg) + 
  geom_bar(mapping=aes(x=reorder(Gram,-mean_per_gram),y=mean_per_gram,fill=Gram,colour='black'),stat='identity',width=0.7,size=my.lwd) +
  geom_errorbar(data=mpg, aes(x=reorder(Gram,-mean_per_gram),ymin=mean_per_gram-se, ymax=mean_per_gram+se),width=0.2, position=position_dodge(0.9)) +
  xlab('Gram Status') + ylab('Grand Mean \\n Number of Migrating Amoebae') +
  scale_fill_manual(mpg$Gram, values=cols) +
  scale_colour_manual(Gram, values=my.cols) +
  scale_x_discrete(labels=c("Negative", "Positive")) + ylim(c(0,200))+
guides(fill = FALSE) +
  guides(size = FALSE) +
  guides(color = FALSE)
save_plot("Figure_1C.pdf",t,base_aspect_ratio=1.1)


################################## CAMP MIGRATION ########################################

## cAMP migration plots, presenting PAIRED choice tests = Fig. 2B)  #####

library(tidyverse)
library(cowplot)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

raw=read.table('cAMP_migration.txt', header=TRUE)
raw2 = raw %>% group_by(strain,bacteria, against,block) %>% summarise(mean_cells=mean(num_cells), se=sd(num_cells)/sqrt(length(num_cells)))
raw3 = raw2 %>% group_by(strain,against, bacteria) %>% summarise(mean=mean(mean_cells), se1=sd(mean_cells)/sqrt(length(mean_cells))) 
camp = raw2 %>% group_by(strain,bacteria) %>% summarise(mean=mean(mean_cells), se1=sd(mean_cells)/sqrt(length(mean_cells))) 

qs32c = filter(camp,strain=='32')
qs39c = filter(camp,strain=='39')
qs40c = filter(camp,strain=='40')

# Ugly code, but accomplishes a re-order to reflect desired plot order
qs32c_pair=filter(raw3,strain=='32')
 row1=qs32c_pair %>% filter(bacteria=='nothing' & against=='N')
 row2=qs32c_pair %>% filter(bacteria=='N' & against=='nothing')
 row3=qs32c_pair %>% filter(bacteria=='N' & against=='WT')
 row4=qs32c_pair %>% filter(bacteria=='WT' & against=='N')
 row5=qs32c_pair %>% filter(bacteria=='WT' & against=='OE')
 row6=qs32c_pair %>% filter(bacteria=='OE' & against=='WT')
 row7=qs32c_pair %>% filter(bacteria=='N' & against=='OE')
 row8=qs32c_pair %>% filter(bacteria=='OE' & against=='N')
 qs32c_ord=rbind(row1,row2,row3,row4,row5,row6,row7,row8)
 
 qs32c_ord$treat=paste(qs32c_ord$bacteria, "(",qs32c_ord$against, ")",sep='')
 qs32c_ord$row=c(1,2,3,4,5,6,7,8)
 qs32c_ord$combo=c("None  Null", "None  Null","Null  WT","Null  WT","WT  OE", "WT  OE", "Null  OE", "Null  OE")
 
 qs39c_pair=filter(raw3,strain=='39')
 row1=qs39c_pair %>% filter(bacteria=='nothing' & against=='N')
 row2=qs39c_pair %>% filter(bacteria=='N' & against=='nothing')
 row3=qs39c_pair %>% filter(bacteria=='N' & against=='WT')
 row4=qs39c_pair %>% filter(bacteria=='WT' & against=='N')
 row5=qs39c_pair %>% filter(bacteria=='WT' & against=='OE')
 row6=qs39c_pair %>% filter(bacteria=='OE' & against=='WT')
 row7=qs39c_pair %>% filter(bacteria=='N' & against=='OE')
 row8=qs39c_pair %>% filter(bacteria=='OE' & against=='N')
 qs39c_ord=rbind(row1,row2,row3,row4,row5,row6,row7,row8)
 qs39c_ord$treat=paste(qs39c_ord$bacteria, "(",qs39c_ord$against, ")",sep='')
 qs39c_ord$row=c(1,2,3,4,5,6,7,8)
 qs39c_ord$combo=c("None  Null", "None  Null","Null  WT","Null  WT","WT  OE", "WT  OE", "Null  OE", "Null  OE")
 
 qs40c_pair=filter(raw3,strain=='40')
 row1=qs40c_pair %>% filter(bacteria=='nothing' & against=='N')
 row2=qs40c_pair %>% filter(bacteria=='N' & against=='nothing')
 row3=qs40c_pair %>% filter(bacteria=='N' & against=='WT')
 row4=qs40c_pair %>% filter(bacteria=='WT' & against=='N')
 row5=qs40c_pair %>% filter(bacteria=='WT' & against=='OE')
 row6=qs40c_pair %>% filter(bacteria=='OE' & against=='WT')
 row7=qs40c_pair %>% filter(bacteria=='N' & against=='OE')
 row8=qs40c_pair %>% filter(bacteria=='OE' & against=='N')
 qs40c_ord=rbind(row1,row2,row3,row4,row5,row6,row7,row8)
 qs40c_ord$treat=paste(qs40c_ord$bacteria, "(",qs40c_ord$against, ")",sep='')
 qs40c_ord$row=c(1,2,3,4,5,6,7,8)
 qs40c_ord$combo=c("None  Null", "None  Null","Null  WT","Null  WT","WT  OE", "WT  OE", "Null  OE", "Null  OE")
 
 my.labels=levels=c(expression(paste("Control    ", Delta, italic("cyaA"))), expression(paste(Delta, italic("cyaA"), "   WT")),expression(paste("  WT  ",italic(  cyaA^OE))),expression(paste(Delta, italic("cyaA  "), italic(cyaA^OE), sep='')))
 
 mig32=qs32c_ord
 mig39=qs39c_ord
 mig40=qs40c_ord
 
# Calculate errorbars 
limits.32 <- aes(ymax = mig32$mean + mig32$se1, ymin = mig32$mean - mig32$se1)
limits.39 <- aes(ymax = mig39$mean + mig39$se1, ymin = mig39$mean - mig39$se1)
limits.40 <- aes(ymax = mig40$mean + mig40$se1, ymin = mig40$mean - mig40$se1)

mig32$combo <- factor(mig32$combo, levels=c("None  Null","Null  WT","WT  OE","Null  OE"))
mig32$treat  <- factor(mig32$treat, levels=c("nothing(N)","N(nothing)","N(WT)","WT(N)","N(OE)", "OE(N)", "WT(OE)", "OE(WT)"))
mig39$combo <- factor(mig39$combo, levels=c("None  Null","Null  WT","WT  OE","Null  OE"))
mig39$treat  <- factor(mig39$treat, levels=c("nothing(N)","N(nothing)","N(WT)","WT(N)","N(OE)", "OE(N)", "WT(OE)", "OE(WT)"))

m1<-ggplot(mig32, aes(x=combo,y=mean,fill=treat)) + geom_bar(stat="identity",position="dodge")+
  geom_errorbar(limits.32, position=position_dodge(0.9),width=0.1,col='black') + 
  scale_x_discrete(labels=my.labels) + 
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold")) +
  labs(x='',y = "Number of \\nMigrating Amoebae") +  
  ggtitle('QS32 (Texas)') + guides(fill=FALSE) +theme(plot.title=element_text(size=20)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

m2<-ggplot(mig39, aes(x=combo,y=mean,fill=treat)) + geom_bar(stat="identity",position="dodge")+
  geom_errorbar(limits.39, position=position_dodge(0.9),width=0.1,col='black') + 
  scale_x_discrete(labels=my.labels) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold")) +
  labs(x='', y = "Number of \\nMigrating Amoebae") +  
  ggtitle('QS39 (Tennessee)') + guides(fill=FALSE)+theme(plot.title=element_text(size=20))

mig40$combo <- factor(mig40$combo, levels=c("None  Null","Null  WT","WT  OE","Null  OE"))
mig40$treat  <- factor(mig40$treat, levels=c("nothing(N)","N(nothing)","N(WT)","WT(N)","N(OE)", "OE(N)", "WT(OE)", "OE(WT)"))
m3<-ggplot(mig40, aes(x=combo,y=mean,fill=treat)) + geom_bar(stat="identity",position="dodge")+
  geom_errorbar(limits.40, position=position_dodge(0.9),width=0.1,col='black') + 
  scale_x_discrete(labels=my.labels) +
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold")) +
  labs(x = "Bacterial Pair", y = "Number of \\nMigrating Amoebae") +  
  ggtitle('QS40 (Massachusetts)') + guides(fill=FALSE)+theme(plot.title=element_text(size=20))

library(cowplot)
require(cowplot)
 outplot <- plot_grid(m1,m2,m3,nrow=3,ncol=1,labels=c("A", "B", "C"), label_size=24)
 save_plot("Figure_2B.pdf", outplot, base_height=11, base_aspect_ratio=0.8)
