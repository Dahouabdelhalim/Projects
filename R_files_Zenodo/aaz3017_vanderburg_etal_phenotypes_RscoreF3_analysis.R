#Van der Burg et al.
#Genomic architecture of a genetically assimilated seasonal color pattern
#Supplemental script
#phenotype F3 cross analysis

library(ggplot2)
library(plyr)
library(lattice)


#plots frequency of mean wing color scores, separated for males and females.
#data.frames
scoresdf<-read.table("aaz3017_vanderburg_etal_Relative_freqscoreF3.txt",header=TRUE)
attach(scoresdf)
scoresdf


#backcross crxhz dataframe and plot

dfcrhzMF<-data.frame(
  group = rep(c("RED X HZ Males","RED X HZ Females"), each=13),
  x=rep(c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0),2),
  y=c(scoresdf$Rfreq[sex=="M"&mxf=="crXhz"],
      scoresdf$Rfreq[sex=="F"&mxf=="crXhz"])
)
dfcrhzMF

ggplot(dfcrhzMF,aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity",position = "identity" ) +
  scale_y_continuous(breaks=seq(-1,1,.1)) +
  scale_x_continuous(breaks=
        c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=c( "chocolate3","chocolate1"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom") +
  ggtitle("RED Male X HZ Female") +
  theme(plot.title=element_text(size=16, face="bold")) +
  xlab("Wing Score") + ylab("Frequency (%)") +
  theme(text = element_text(size=15))


#backcross hzxcr dataframe and plot
dfhzcrMF<-data.frame(
  group = rep(c("HZ X RED Males","HZ X RED Females"), each=13),
  x=rep(c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0),2),
  y=c(scoresdf$Rfreq[sex=="M"&mxf=="hzXcr"],
      scoresdf$Rfreq[sex=="F"&mxf=="hzXcr"])
)
dfhzcrMF

ggplot(dfhzcrMF,aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity",position = "identity" ) +
  scale_y_continuous(breaks=seq(-1,1,.1)) +
  scale_x_continuous(breaks=
                       c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=c( "chocolate3","chocolate1"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom") +
  ggtitle("HZ Male X RED Female") +
  theme(plot.title=element_text(size=16, face="bold")) +
  xlab("Wing Score") + ylab("Frequency (%)") +
  theme(text = element_text(size=15))


#heterozygote hzxhz dataframe and plot
dfhzhzMF<-data.frame(
  group = rep(c("HZ X HZ Males","HZ X HZ Females"), each=13),
  x=rep(c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0),2),
  y=c(scoresdf$Rfreq[sex=="M"&mxf=="hzXhz"],
      scoresdf$Rfreq[sex=="F"&mxf=="hzXhz"])
)
dfcrhzMF

ggplot(dfhzhzMF,aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity",position = "identity" ) +
  scale_y_continuous(breaks=seq(-1,1,.1)) +
  scale_x_continuous(breaks=
                       c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=c( "chocolate3","chocolate1"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom") +
  ggtitle("HZ Male X HZ Female") +
  theme(plot.title=element_text(size=16, face="bold")) +
  xlab("Wing Score") + ylab("Frequency (%)") +
  theme(text = element_text(size=15))


#F3 bcXbc dataframe and plot

dfbcbcMF<-data.frame(
  group = rep(c("F3 Males","F3 Females"), each=13),
  x=rep(c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0),2),
  y=c(scoresdf$Rfreq[sex=="M"&mxf=="bcXbc"],
      scoresdf$Rfreq[sex=="F"&mxf=="bcXbc"])
)
dfbcbcMF

ggplot(dfbcbcMF,aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity",position = "identity" ) +
  scale_y_continuous(breaks=seq(-1,1,.1)) +
  scale_x_continuous(breaks=
                       c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=c( "chocolate3","chocolate1"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom") +
  ggtitle("BackCross Male X BackCross Female") +
  theme(plot.title=element_text(size=16, face="bold")) +
  xlab("Wing Score") + ylab("Frequency (%)") +
  theme(text = element_text(size=15))

#F3 plXhz dataframe and plot

dfplhzMF<-data.frame(
  group = rep(c("PLAS X HZ Males","PLAS x HZ Females"), each=13),
  x=rep(c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0),2),
  y=c(scoresdf$Rfreq[sex=="M"&mxf=="plXhz"],
      scoresdf$Rfreq[sex=="F"&mxf=="plXhz"])
)
dfplhzMF

ggplot(dfplhzMF,aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity",position = "identity" ) +
  scale_y_continuous(breaks=seq(-1,1,.1)) +
  scale_x_continuous(breaks=
                       c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=c( "chocolate1","chocolate3"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom") +
  ggtitle("PLAS Male X HZ Female") +
  theme(plot.title=element_text(size=16, face="bold")) +
  xlab("Wing Score") + ylab("Frequency (%)") +
  theme(text = element_text(size=15))


#F1 cross plXcr 

dfplcrMF<-data.frame(
  group = rep(c("PLAS X RED Males","PLAS x RED Females"), each=13),
  x=rep(c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0),2),
  y=c(scoresdf$Rfreq[sex=="M"&mxf=="plXcr"],
      scoresdf$Rfreq[sex=="F"&mxf=="plXcr"])
)
dfplcrMF

ggplot(dfplcrMF,aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity",position = "identity" ) +
  scale_y_continuous(breaks=seq(-1,1,.1)) +
  scale_x_continuous(breaks=
                       c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=c( "chocolate1","chocolate3"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom") +
  ggtitle("PLAS Male X RED Female") +
  theme(plot.title=element_text(size=16, face="bold")) +
  xlab("Wing Score") + ylab("Frequency (%)") +
  theme(text = element_text(size=15))




#F1 cross crXpl 

dfcrplMF<-data.frame(
  group = rep(c("RED X PLAS Males","RED x PLAS Females"), each=13),
  x=rep(c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0),2),
  y=c(scoresdf$Rfreq[sex=="M"&mxf=="crXpl"],
      scoresdf$Rfreq[sex=="F"&mxf=="crXpl"])
)
dfcrplMF

ggplot(dfcrplMF,aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity",position = "identity" ) +
  scale_y_continuous(breaks=seq(-1,1,.1)) +
  scale_x_continuous(breaks=
                       c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=c( "chocolate1","chocolate3"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom") +
  ggtitle("RED Male X PLAS Female") +
  theme(plot.title=element_text(size=16, face="bold")) +
  xlab("Wing Score") + ylab("Frequency (%)") +
  theme(text = element_text(size=15))


#P cross plXpl 

dfplplMF<-data.frame(
  group = rep(c("PLAS X PLAS Males","PLAS x PLAS Females"), each=13),
  x=rep(c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0),2),
  y=c(scoresdf$Rfreq[sex=="M"&mxf=="plXpl"],
      scoresdf$Rfreq[sex=="F"&mxf=="plXpl"])
)
dfplplMF

ggplot(dfplplMF,aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity",position = "identity" ) +
  scale_y_continuous(breaks=seq(-1,1,.1)) +
  scale_x_continuous(breaks=
                       c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=c( "chocolate1","chocolate3"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom") +
  ggtitle("PLAS Male X PLAS Female") +
  theme(plot.title=element_text(size=16, face="bold")) +
  xlab("Wing Score") + ylab("Frequency (%)") +
  theme(text = element_text(size=15))

#P cross crXcr 

dfcrcrMF<-data.frame(
  group = rep(c("RED X RED Males","RED x RED Females"), each=13),
  x=rep(c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0),2),
  y=c(scoresdf$Rfreq[sex=="M"&mxf=="crXcr"],
      scoresdf$Rfreq[sex=="F"&mxf=="crXcr"])
)
dfcrcrMF

ggplot(dfcrcrMF,aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity",position = "identity" ) +
  scale_y_continuous(breaks=seq(-1,1,.1)) +
  scale_x_continuous(breaks=
                       c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=c( "chocolate1","chocolate3"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom") +
  ggtitle("RED Male X RED Female") +
  theme(plot.title=element_text(size=16, face="bold")) +
  xlab("Wing Score") + ylab("Frequency (%)") +
  theme(text = element_text(size=15))

replications(~.-ID,cr)

dfcrcrMF<-data.frame(
  group = rep(c("RED X RED Males","RED x RED Females"), each=13),
  x=rep(c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0),2),
  y=c(scoresdf$Rfreq[sex=="M"&mxf=="crXcr"],
      scoresdf$Rfreq[sex=="F"&mxf=="crXcr"])
)
dfcrcrMF

ggplot(dfcrcrMF,aes(x=x, y=y, fill=group)) +
  geom_bar(stat="identity",position = "identity" ) +
  scale_y_continuous(breaks=seq(-1,1,.1)) +
  scale_x_continuous(breaks=
                       c(1.0,1.3,1.7,2.0,2.3,2.7,3.0,3.3,3.7,4.0,4.3,4.7,5.0)) +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=c( "chocolate1","chocolate3"))+
  guides(fill=guide_legend(title=NULL)) +
  theme(legend.position="bottom") +
  ggtitle("RED Male X RED Female") +
  theme(plot.title=element_text(size=16, face="bold")) +
  xlab("Wing Score") + ylab("Frequency (%)") +
  theme(text = element_text(size=15))

count(cr$cross)
?table()
latmean<-data.frame(group=table(cr$meannum,cr$cross,cr$sex))
colnames(latmean)<-c('mean','cross','sex','Freq')
attach(latmean)
latmean

barchart(latmean$Freq[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]~
           latmean$mean[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]|
           latmean$sex[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]*
           latmean$cross[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"],
         horizontal = FALSE,xlab="score",ylab="Freq", main='MEAN')

barchart(latmean$Freq[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]~
           latmean$mean[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]|
           latmean$sex[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]*
           latmean$cross[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"],
         horizontal = FALSE,xlab="score",ylab="Freq", main='MEAN')

barchart(latmean$Freq[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]~
           latmean$mean[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]|
           latmean$sex[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]*
           latmean$cross[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"],
         horizontal = FALSE,xlab="score",ylab="Freq", main='MEAN')





lathu<-data.frame(group=table(cr$hwdist,cr$cross,cr$sex))
colnames(lathu)<-c('hu','cross','sex','Freq')
attach(lathu)
lathu
barchart(lathu$Freq[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]~
           lathu$hu[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]|
           lathu$sex[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]*
           lathu$cross[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"],
         horizontal = FALSE,xlab="score",ylab="Freq", main="HW UMBRAL")

barchart(lathu$Freq[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]~
           lathu$hu[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]|
           lathu$sex[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]*
           lathu$cross[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"],
         horizontal = FALSE,xlab="score",ylab="Freq", main="HW UMBRAL")

barchart(lathu$Freq[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]~
           lathu$hu[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]|
           lathu$sex[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]*
           lathu$cross[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"],
         horizontal = FALSE,xlab="score",ylab="Freq", main="HW UMBRAL")


lathm<-data.frame(group=table(cr$hwmed,cr$cross,cr$sex))
colnames(lathm)<-c('hm','cross','sex','Freq')
attach(lathm)
lathm
barchart(lathm$Freq[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]~
           lathm$hm[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]|
           lathm$sex[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]*
           lathm$cross[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"],
         horizontal = FALSE,xlab="score",ylab="Freq", main="HW MEDIAL")

barchart(lathm$Freq[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]~
           lathm$hm[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]|
           lathm$sex[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]*
           lathm$cross[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"],
         horizontal = FALSE,xlab="score",ylab="Freq", main="HW MEDIAL")

barchart(lathm$Freq[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]~
           lathm$hm[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]|
           lathm$sex[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]*
           lathm$cross[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"],
         horizontal = FALSE,xlab="score",ylab="Freq", main="HW MEDIAL")


latfw<-data.frame(group=table(cr$fw,cr$cross,cr$sex))
colnames(latfw)<-c('fw','cross','sex','Freq')
attach(latfw)
latfw
barchart(latfw$Freq[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]~
           latfw$fw[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]|
           latfw$sex[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"]*
           latfw$cross[cross=="f2-18"|cross=="f2-20"|cross=="f2-24"|cross=="f2-27"],
         horizontal = FALSE,xlab="score",ylab="Freq", main="FW")

barchart(latfw$Freq[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]~
           latfw$fw[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]|
           latfw$sex[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"]*
           latfw$cross[cross=="hz21"|cross=="hz23"|cross=="hz5"|cross=="hz6"|cross=="hz7"],
         horizontal = FALSE,xlab="score",ylab="Freq", main="FW")

barchart(latfw$Freq[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]~
           latfw$fw[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]|
           latfw$sex[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"]*
           latfw$cross[cross=="hz1"|cross=="hz13"|cross=="hz15"|cross=="hz16"|cross=="hz2"],
         horizontal = FALSE,xlab="score",ylab="Freq", main="FW")


