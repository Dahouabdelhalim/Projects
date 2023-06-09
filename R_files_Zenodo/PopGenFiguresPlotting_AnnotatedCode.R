#### code for generating recombination rate plots

## Import PopgenStats
# this file contains loess values generated in r/MareyMaps, which is a GUI. In this case, i made a loess transform of my data with a span of 0.5 and a degree of 1.
# (outputs from r/MareyMap are presented in loess.txt)
# Population genetics statistics were generated using code from Simon Martin's Github https://github.com/simonhmartin/genomics_general
# Positions were transformed so that they could be plotted, using code by Steven Van Belleghem, found in Manhattan_SVB.r
stats <- read.table('PopgenStats.txt',header=TRUE)

# (Previously, in order to get the scaffolds into the correct orientation, I generated a list of scaffolds to flip
InvertList <- c("Sc0000009_pilon","Sc0000044_pilon","Sc0000098_pilon","Sc0000002_pilon","Sc0000094_pilon","Sc0000033_pilon","Sc0000073_pilon","Sc0000011_pilon","Sc0000029_pilon","Sc0000103_pilon","Sc0000086_pilon","Sc0000065_pilon","Sc0000030_pilon","Sc0000018_pilon","Sc0000102_pilon","Sc0000016_pilon","Sc0000025_pilon","Sc0000056_pilon","Sc0000048_pilon","Sc0000050_pilon","Sc0000051_pilon","Sc0000006_pilon","Sc0000027_pilon","Sc0000040_pilon","Sc0000062_pilon","Sc0000008_pilon","Sc0000064_pilon","Sc0000095_pilon","Sc0000015_pilon","Sc0000084_pilon","Sc0000021_pilon","Sc0000005_pilon","Sc0000072_pilon","Sc0000032_pilon","Sc0000058_pilon","Sc0000042_pilon","Sc0000039_pilon","Sc0000054_pilon","Sc0000113_pilon","Sc0000045_pilon","Sc0000026_pilon","Sc0000007_pilon","Sc0000037_pilon","Sc0000091_pilon","Sc0000077_pilon", "Sc0000046_pilon")

# this loop retrieves the end value for the scaffold, then flips every mid position to the correct orientation. Does not need to be run on contents of PopgenStats.txt)
for (i in InvertList){
  j <- max(stats$end[stats$scaffold==i])
 stats$mid[stats$scaffold==i] <- j - stats$mid[stats$scaffold==i]
}

## sort the object, allowing proper plotting of line graphs
stats <- stats[with(stats, order(chromosome, stats)),]

## make all negative values == 0

stats[stats < 0] <- 0
stats <- na.omit(stats)


###########################################################################
# plot of cM per Mb (loess) recombination values as a Manhattan line chart#
###########################################################################

# I put my files in a list. If I have multiple, it can loop and plot one underneath the other (adjust mfrow!)
comp=list(stats)
names=c("fileXX_name")
col=c("black")
par(mfrow=c(1,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)
top = 40
bot = 0.1

#check how many chromosomes you have
begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[31]/1000000

for (i in 1:length(comp)){
  # plot background
  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))

  # plot your stats
  par(new=TRUE)
  plot(comp[[i]]$chromPos/1000000,comp[[i]]$loess,type="l",pch=19, cex=.5,col=adjustcolor(col[i], alpha=1),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n",lwd=3)
  axis(2,cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("Pi")), cex=0.5,line = 0.8)
}

#label chromosomes
axis(1, at=chrom_coords[,5][1:31]/1000000, labels=(1:31),lwd=0, lwd.ticks=0)

# plot scale (you might have to adjust this to show up)
segments(c(390,390,400), c(250.07,250.0702,250.0702), c(400,390,400) ,c(250.07,250.0702,250.0702), lwd = 1)
text(395,300.086,labels = "10Mb", cex = 1)




##############################################
# plot of cM/Mb vs Fst graph, with trendline #
##############################################

# highlight values of loess greater than the 0.95 quantile 
stats$Percentile <- ifelse(stats$loess<quantile(stats$loess, .5), "F",
   ifelse(stats$Fst_pos_neg<quantile(stats$Fst_pos_neg, .95),"F", "T"))
highlight <- stats[stats$Percentile=='T',]

##### ggplot for plotting fst vs loess
cmMbfig <- ggplot(stats, aes(loess,Fst_pos_neg))+
  geom_point(alpha=0.05,col='black')+
  geom_smooth(,se=F,span=0.5)+
  geom_point(data=highlight,col='red',alpha=0.2)+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  xlab("Recombination Rate, cM/Mb")+
  ylab("Fst")+
  theme_classic()+
  theme(axis.title = element_text(size = 20))             

###################################################################
# Manhattan plot of Fst with SNPs highighted above coloured in red#
###################################################################
comp=list(stats)
names=c("Fst")
col=c("black")


par(mfrow=c(1,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

top = 1
bot = 0

#check how many chromosomes you have
begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[31]/1000000

for (i in 1:length(comp)){
  # plot background
  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))

  # plot your stats
  par(new=TRUE)
  plot(comp[[i]]$chromPos/1000000,comp[[i]]$Fst_pos_neg,type="p",pch=19, cex=.5,col=adjustcolor(col[i], alpha=1),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n",lwd=3)
  axis(2,cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("Pi")), cex=0.5,line = 0.8)
}

#put on the sig line
# plot chromosome labels
axis(1, at=chrom_coords[,5][1:32]/1000000, labels=(1:32),lwd=0, lwd.ticks=0)

# this adds the highlighted points to the Fst plot 
points(highlight$chromPos/1000000,highlight$Fst_pos_neg,type="p",pch=19, cex=.5,col=adjustcolor('red', alpha=1))

####################################
###### finalized Popgen Plots ######
####################################

##########
##pi_neg##
##########

comp=list(stats)
names=c("Fst")
col=c("black")


par(mfrow=c(1,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

top = 0.3
bot = 0

#check how many chromosomes you have
begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[31]/1000000

for (i in 1:length(comp)){
  # plot background
  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))

  # plot your stats
  par(new=TRUE)
  plot(comp[[i]]$chromPos/1000000,comp[[i]]$pi_neg,type="p",pch=19, cex=.5,col=adjustcolor(col[i], alpha=0.9),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n",lwd=3)
  axis(2,cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("Pi")), cex=0.5,line = 0.8)
}

# plot chromosome labels
axis(1, at=chrom_coords[,5][1:31]/1000000, labels=(1:31),lwd=0, lwd.ticks=0)

dev.copy(pdf,'pi_neg.pdf',width=20,height=3)
dev.off()

##########
##pi_pos##
##########

comp=list(stats)
names=c("Fst")
col=c("black")


par(mfrow=c(1,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

top = 0.3
bot = 0

#check how many chromosomes you have
begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[31]/1000000

for (i in 1:length(comp)){
  # plot background
  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))

  # plot your stats
  par(new=TRUE)
  plot(comp[[i]]$chromPos/1000000,comp[[i]]$pi_pos,type="p",pch=19, cex=.5,col=adjustcolor(col[i], alpha=0.9),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n",lwd=3)
  axis(2,cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("Pi")), cex=0.5,line = 0.8)
}

# plot chromosome labels
axis(1, at=chrom_coords[,5][1:31]/1000000, labels=(1:31),lwd=0, lwd.ticks=0)

dev.copy(pdf,'pi_pos.pdf',width=20,height=3)
dev.off()

###############
##dxy_pos_neg##
###############

comp=list(stats)
names=c("Fst")
col=c("black")


par(mfrow=c(1,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

top = 0.4
bot = 0

#check how many chromosomes you have
begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[31]/1000000

for (i in 1:length(comp)){
  # plot background
  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))

  # plot your stats
  par(new=TRUE)
  plot(comp[[i]]$chromPos/1000000,comp[[i]]$dxy_pos_neg,type="p",pch=19, cex=.5,col=adjustcolor(col[i], alpha=0.9),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n",lwd=3)
  axis(2,cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("Pi")), cex=0.5,line = 0.8)
}

# plot chromosome labels
axis(1, at=chrom_coords[,5][1:31]/1000000, labels=(1:31),lwd=0, lwd.ticks=0)

dev.copy(pdf,'dxy_pos_neg.pdf',width=20,height=3)
dev.off()



##finalised unhighlighted fst
comp=list(stats)
names=c("Fst")
col=c("black")


par(mfrow=c(1,1),mai=c(0.05,0.8,0.05,0.8), oma=c(2,0,1,0)+0)

top = 0.8
bot = 0

#check how many chromosomes you have
begin = chrom_coords$chromStarts[1]/1000000
end = chrom_coords$chromEnds[31]/1000000

for (i in 1:length(comp)){
  # plot background
  plot(0, pch = "",xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "n", main = "",axes = FALSE)
  rect(chrom_coords$chromStarts/1000000,rep(bot,length(chrom_coords$chromStarts)),chrom_coords$chromEnds/1000000,rep(top,length(chrom_coords$chromStarts)), col = c("gray95","gray90"), lwd = 0, border = c("gray95","gray90"))

  # plot your stats
  par(new=TRUE)
  plot(comp[[i]]$chromPos/1000000,comp[[i]]$Fst_pos_neg,type="p",pch=19, cex=.5,col=adjustcolor(col[i], alpha=0.9),xlim = c(begin,end), ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "",yaxt="n",xaxt="n",lwd=3)
  axis(2,cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("Pi")), cex=0.5,line = 0.8)
}

# plot chromosome labels
axis(1, at=chrom_coords[,5][1:31]/1000000, labels=(1:31),lwd=0, lwd.ticks=0)

dev.copy(pdf,'fst_noHighlight.pdf',width=20,height=3)
dev.off()


########################################
### Fst, autosomes vs sex chr plots ####
########################################

stats$SorA <- ifelse(stats$chromosome == 1, "Z Chromosome", "Autosomes")

# set negative values to zero
#stats$Fst_pos_neg[stats$Fst_pos_neg < 0] <- 0

ggplot(stats, aes(Fst_pos_neg, col=SorA,fill=SorA))+
  geom_density(alpha=0.2,n=64)+
  theme_classic()+
  theme(legend.text=element_text(size=11),legend.position=c(0.8,1), legend.justification=c('right','top'),legend.title=element_blank())+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
  annotate("text", x = 0.2, y=57.7, label = paste("mean Fst on Autosomes =",round(mean(stats$Fst_pos_neg[stats$SorA=="Autosomes"]),4)))+
  annotate("text", x = 0.2, y=54.5, label = paste("mean Fst on Z Chromosome =",round(mean(stats$Fst_pos_neg[stats$SorA=="Z Chromosome"]),4)))

########################################
###  Pi, autosomes vs sex chr plots ####
########################################

# (import jack function written by Steven Van Belleghem)
jack <- function(stat, w){

  start <- 0
  end <- w

  n <- 1

  D <- sum(stat)/length(stat)

  len <- as.integer(length(stat)/w)
  pseudo_stat <- vector(length = len)

  while(end <= length(stat)){

    jacked_stat <- sum(stat[-c(start:end)])/length(stat[-c(start:end)])

    pseudo_stat[n] <- D*len - jacked_stat*(len-1)

    start <- start + w
    end <- start + w
    n <- n + 1
  }
  return(pseudo_stat)
}

PIjacks <- list(pi_negA =jack(piA$pi_neg, 10), 
                pi_negZ = jack(piZ$pi_neg, 10), 
                pi_posA = jack(piA$pi_pos, 10),
                pi_posZ = jack(piZ$pi_pos, 10))

PIJpi_negA <- sum(PIjacks$pi_negA)/length(PIjacks$pi_negA)
PIsdpi_negA <- sd(PIjacks$pi_negA)
PIerrpi_negA <- sqrt(var(PIjacks$pi_negA)/length(PIjacks$pi_negA))

PIJpi_negZ <- sum(PIjacks$pi_negZ)/length(PIjacks$pi_negZ)
PIsdpi_negZ <- sd(PIjacks$pi_negZ)
PIerrpi_negZ <- sqrt(var(PIjacks$pi_negZ)/length(PIjacks$pi_negZ))

PIJpi_posA <- sum(PIjacks$pi_posA)/length(PIjacks$pi_posA)
PIsdpi_posA <- sd(PIjacks$pi_posA)
PIerrpi_posA <- sqrt(var(PIjacks$pi_posA)/length(PIjacks$pi_posA))

PIJpi_posZ <- sum(PIjacks$pi_posZ)/length(PIjacks$pi_posZ)
PIsdpi_posZ <- sd(PIjacks$pi_posZ)
PIerrpi_posZ <- sqrt(var(PIjacks$pi_posZ)/length(PIjacks$pi_posZ))

PInames <- c("negA","negZ","posA","posZ")
PIJ <- c(PIJpi_negA,PIJpi_negZ,PIJpi_posA,PIJpi_posZ)
PIsd <- c(PIsdpi_negA,PIsdpi_negZ,PIsdpi_posA,PIsdpi_posZ)
PIJT <- c(PIJpi_negA/PIJpi_negA, PIJpi_negZ/PIJpi_negA, PIJpi_posA/PIJpi_posA, PIJpi_posZ/PIJpi_posA)
PIsdT <- c(PIsdpi_negA/PIJpi_negA, PIsdpi_negZ/PIJpi_negA,PIsdpi_posA/PIJpi_negA, PIsdpi_posZ/PIJpi_negA)
PIStats <- data.frame(PInames,PIJ,PIsd,PIJT,PIsdT)

## plot in ggplot
PIStats$group <- c('Cph','Cph','Ceu','Ceu')
PIStats$chr <- c('A','Z','A','Z')

col <- c('#49AAE0','#a6d9f1', '#b10079','#e3a5c8')
ggplot(PIStats, aes(x=PInames, PIJT, fill=PInames))+
  geom_bar(stat='identity')+
  geom_errorbar(ymin=PIJT-PIsdT,ymax=PIJT+PIsdT)+
  ylim(0,1.2)+
  scale_fill_manual(values=col)+
  theme_classic()

#############################################################
### code for plotting PCAs, data generated wtih smartpca ####
#############################################################

## load packages
library(ggplot2)		# plotting
library(ggrepel)		# for making visually appealing labels
library(cowplot)		#for making faceted plots

## plot command
pcaDO <- read.table('HedgeapplePCA_DroppedOutliers.txt',h=T)
pcaDOZ <- read.table('HedgeapplePCA_DroppedOutliersZ.txt',h=T)

plot.pcaALL <- ggplot(pca,aes(PC1,PC2,col=group))+
  geom_point()+
  ggtitle("All individuals, autosomes")+
  scale_color_manual(labels=c('C. eurytheme','C. philodice'), values=c('orange','gold'))+
  geom_label_repel(aes(label = name,),colour='black',max.overlaps=200)+
  labs(x="PC1 2.631%", y="PC2 1.799%", col="Species")+
  theme_classic()+
  theme(legend.position='none')+
  theme(panel.background = element_rect(colour = "black", size=1, fill=NA),axis.text=element_text(size=20),axis.title=element_text(size=14,face="bold"))

plot.pcaDROP <- ggplot(pcaDO,aes(PC1,PC2,col=group))+
  geom_point()+
  ggtitle("Dropped, autosomes")+
  scale_color_manual(labels=c('C. eurytheme','C. philodice'), values=c('orange','gold'))+
  geom_label_repel(aes(label = name,),colour='black',max.overlaps=200)+
  labs(x="PC1 2.427%", y="PC2 1.032%", col="Species")+
  theme_classic()+
  theme(legend.position='none')+
  theme(panel.background = element_rect(colour = "black", size=1, fill=NA),axis.text=element_text(size=20),axis.title=element_text(size=14,face="bold"))

plot.pcaDROPZ <- ggplot(pcaDOZ,aes(PC1,PC2,col=group))+
  geom_point()+
  ggtitle("Dropped, Z chromosome")+
  scale_color_manual(labels=c('C. eurytheme','C. philodice'), values=c('orange','gold'))+
  geom_label_repel(aes(label = name,),colour='black',max.overlaps=200)+
  labs(x="PC1 9.614%", y="PC2 0.712%", col="Species")+
  theme_classic()+
  theme(legend.position='none')+
  theme(panel.background = element_rect(colour = "black", size=1, fill=NA),axis.text=element_text(size=20),axis.title=element_text(size=14,face="bold"))

plot_grid(plot.pcaALL, plot.pcaDROP,plot.pcaDROPZ)




