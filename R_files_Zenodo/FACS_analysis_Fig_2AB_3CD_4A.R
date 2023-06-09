## Flow cytometry data visualization from evolution experiments in Figure 2A-B
## additional analysis of the same evolution experiment shown in Fig. 3C-D and Figure 4A

## Flow cytometry plots of all 96-well plates containing populations shown in Figure 2A-B. 
## compare different 96-well plates (different days, same settings) medium B + C, (in different folders)
## frozen sample pinned 1:1000 into M9 buffer on ice --> FACS canto

#"IS" = delta IS; "30" = IT030 (IS wt)

# data: EE24.4/8/12 (evolution experiment 24, day 4/8/12) plates 1-3 (delta IS1C - medium C,B,A) + plates 7-10 (IT030 - medium C, B, A; controll plate medium E)
# evolution in Gal (A-1%, B-0.1%, C-0.01%) + 0.1% CAS; plate E = CAS only ctr.

#FlowJo: select all, export as scale values into new folder with naming scheme: "ScaleVal_EE24_12C30"
#use autogating in flow jo to gate the a single concise population of cells (by eye, same within plate, and similar (~60%) for all different plates)
#folder ...gated; files ...H3g

library("plyr")
library("ggplot2")
library("Rmisc")
#### data locations (see below): 
#for medium C data
#("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.12C30 etc. ")

#for medium B data
##("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190606_FACS/ScaleVal_EE2412B30")

## medium C - 12C30
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.12C30_gated")
### read in files in a loop#  # # # # # # # # # # # # # # # # # # # # # # # # # # # 
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_12C30=filelist

## medium C - 8C30 
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.8C30_gated")
### read in files in a loop
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_8C30=filelist

## medium C - 4C30 
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.4C30_gated")
### read in files in a loop
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_4C30=filelist

## medium C - 8CIS
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190902_EE24/ScaleVal_EE24.8CIS_gated")
### read in files in a loop
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_8CIS=filelist

## medium C - 12CIS
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190906_EE24/ScaleVal_EE24.12CIS_gated")
### read in files in a loop
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_12CIS=filelist

## medium B - 12B30
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190606_FACS/ScaleVal_EE24.12B30_gated")
### read in files in a loop
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_12B30=filelist

## medium B - 8B30
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190902_EE24/ScaleVal_EE24.8B30_gated")
### read in files in a loop
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_8B30=filelist

## medium B - 12BIS (delta ISC)
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190606_FACS/ScaleVal_EE24.12BIS_gated")
### read in files in a loop
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_12BIS=filelist



## medium A - 12A30 
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190923_EE24/ScaleVal_EE24.12A30_gated")
### read in files in a loop
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_12A30=filelist

## medium A - 12AIS
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190923_EE24/ScaleVal_EE24.12AIS_gated")
### read in files in a loop
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_12AIS=filelist

## medium E - 8E  #ancestral controls (only some wells, all look the same!)
setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.8E_gated")
### read in files in a loop
filelist = list.files(pattern = ".*.txt")
for (i in 1:length(filelist)){ 
  x <- read.csv(filelist[i], header=TRUE, sep="\\t") 
  assign(paste0(filelist[i]), x)
}#end loop
filelist_8E=filelist



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#C 30 typical plot (amplification gained in low galactose)
ggplot(EE24_12C30_F3g.txt,aes(log(FITC.H)))+
  geom_density(data=EE24_4C30_F3g.txt, aes(log(FITC.H), col="4"))+
  geom_density(data=EE24_8C30_F3g.txt, aes(log(FITC.H), col="8"))+
  geom_density(data=EE24_12C30_F3g.txt, aes(log(FITC.H), col="12"))

ggplot(EE24_12C30_F3g.txt,aes(log(Pacific.Blue.H)))+
  geom_density(data=EE24_4C30_F3g.txt, aes(log(Pacific.Blue.H), col="4"))+
  geom_density(data=EE24_8C30_F3g.txt, aes(log(Pacific.Blue.H), col="8"))+
  geom_density(data=EE24_12C30_F3g.txt, aes(log(Pacific.Blue.H), col="12"))#

#A 30 typical plot
ggplot(EE24_12A30_F2g.txt,aes(log(FITC.H)))+
  geom_density(data=EE24_8E_B3g.txt, aes(log(FITC.H), col="ancestral"))+
  geom_density(data=EE24_12A30_F2g.txt, aes(log(FITC.H), col="12"))

ggplot(EE24_12C30_F2g.txt,aes(log(Pacific.Blue.H)))+
  geom_density(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H), col="ancestral"))+
  geom_density(data=EE24_12A30_F2g.txt, aes(log(Pacific.Blue.H), col="12"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### LOOPS FOR EXPORTING VARIOUS PLOTS (ALL WELLS)   ### ### ###

## FOR C - * * 30 * * 
#log(YFP/CFP)- save ggplots for 96 well in a loop -filelist_8C30, 430, 12C30
for (i in 1:length(filelist_4C30)) {
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.4C30_gated")
  file4 <- read.csv(filelist_4C30[i], header=TRUE, sep="\\t") 
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.8C30_gated")
  file8 <- read.csv(filelist_8C30[i], header=TRUE, sep="\\t") 
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.12C30_gated")
  file12<- read.csv(filelist_12C30[i], header=TRUE, sep="\\t")
  p <- ggplot(file4,aes(log(FITC.H/Pacific.Blue.H)))+
    xlim(-3,4)+
    ylim(0,1.5)+
    geom_density(data=EE24_8E_B3g.txt,aes(log(FITC.H/Pacific.Blue.H),color="ancestral"))+
    geom_density(data=file4, aes(log(FITC.H/Pacific.Blue.H), col="4"))+
    geom_density(data=file8, aes(log(FITC.H/Pacific.Blue.H), col="8"))+
    geom_density(data=file12, aes(log(FITC.H/Pacific.Blue.H), col="12"))
  png(paste("plot_", filelist_8C30[i], "ratio.png", sep = ""), width=600, height=500, res=120) # start export
  print(p) 
  dev.off() # finish export
}




## FOR C - * * 30 * * 
#YFP - save ggplots for 96 well in a loop -filelist_8C30, 430, 12C30
for (i in 1:length(filelist_4C30)) {
  #setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.4C30_gated")
  #file4 <- read.csv(filelist_4C30[i], header=TRUE, sep="\\t") 
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.8C30_gated")
  file8 <- read.csv(filelist_8C30[i], header=TRUE, sep="\\t") 
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.12C30_gated")
  file12<- read.csv(filelist_12C30[i], header=TRUE, sep="\\t")
  p <- ggplot(file4,aes(log(FITC.H)))+
     xlim(2,7.5)+
    ylim(0,1.0)+
    geom_density(data=EE24_8E_B3g.txt,aes(log(FITC.H),color="ancestral"))+
  #  geom_density(data=file4, aes(log(FITC.H), col="4"))+
    geom_density(data=file8, aes(log(FITC.H), col="8"))+
    geom_density(data=file12, aes(log(FITC.H), col="12"))
  png(paste("plot_", filelist_8C30[i], "YFP.png", sep = ""), width=600, height=500, res=120) # start export
  print(p) 
  dev.off() # finish export
}


## FOR C - * * 30 * * 
#CFP - save ggplots for 96 well in a loop -filelist_8C30, 430, 12C30
for (i in 1:length(filelist_4C30)) {
  #setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.4C30_gated")
  file4 <- read.csv(filelist_4C30[i], header=TRUE, sep="\\t") 
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.8C30_gated")
  file8 <- read.csv(filelist_8C30[i], header=TRUE, sep="\\t") 
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.12C30_gated")
  file12<- read.csv(filelist_12C30[i], header=TRUE, sep="\\t")
  p <- ggplot(file4,aes(log(Pacific.Blue.H)))+
    xlim(2,7.5)+
    ylim(0,1.5)+
    geom_density(data=EE24_8E_B3g.txt,aes(log(Pacific.Blue.H),color="ancestral"))+
    geom_density(data=file4, aes(log(Pacific.Blue.H), col="4"))+
    geom_density(data=file8, aes(log(Pacific.Blue.H), col="8"))+
    geom_density(data=file12, aes(log(Pacific.Blue.H), col="12"))
  png(paste("plot_", filelist_8C30[i], "CFP.png", sep = ""), width=600, height=500, res=120) # start export
  print(p) 
  dev.off() # finish export
}
getwd()

## FOR C - * * 30 * * 
#YFP versus CFP - save ggplots for 96 well in a loop -filelist_8C30, 430, 12C30
for (i in 1:length(filelist_4C30)) {
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.4C30_gated")
  file4 <- read.csv(filelist_4C30[i], header=TRUE, sep="\\t") 
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.8C30_gated")
  file8 <- read.csv(filelist_8C30[i], header=TRUE, sep="\\t") 
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.12C30_gated")
  file12<- read.csv(filelist_12C30[i], header=TRUE, sep="\\t")  
  p <- ggplot(file4,aes(x=log(Pacific.Blue.H),y=log(FITC.H)))+
    xlim(2,8)+
    ylim(1,8)+
    geom_point(data=EE24_8E_B3g.txt,aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="ancestral"))+
   # geom_point(data=file4, aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="4"))+
    geom_point(data=file8, aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="8"),alpha=0.5)+
    geom_point(data=file12, aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="12"), alpha=0.2)+
  geom_point(data=EE24_8E_B3g.txt,aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="ancestral"), alpha=.05)
  png(paste("plot_", filelist_8C30[i], "both.png", sep = ""), width=600, height=500, res=120) # start export
  print(p) 
  dev.off() # finish export
}



## FOR A - * * 30 * * 
#YFP versus CFP - save ggplots for 96 well in a loop -filelist_8C30, 430, 12C30
for (i in 1:length(filelist_12A30)) {
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190923_EE24/ScaleVal_EE24.12A30_gated") 
  file12<- read.csv(filelist_12A30[i], header=TRUE, sep="\\t")  
  p <- ggplot(file12,aes(x=log(Pacific.Blue.H),y=log(FITC.H)),col="day12")+
    geom_point(data=file12)+
    xlim(2,8)+
    ylim(1,8)+
    geom_point(data=EE24_8E_B3g.txt,aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="ancestral"), alpha=0.2)
  png(paste("plot_", filelist_12A30[i], "both.png", sep = ""), width=600, height=500, res=120) # start export
  print(p) 
  dev.off() # finish export
}


## FOR A - * * IS * * 
#YFP versus CFP - save ggplots for 96 well in a loop -filelist_8C30, 430, 12C30
for (i in 1:length(filelist_12AIS)) {
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190923_EE24/ScaleVal_EE24.12AIS_gated") 
  file12<- read.csv(filelist_12AIS[i], header=TRUE, sep="\\t")  
  p <- ggplot(file12,aes(x=log(Pacific.Blue.H),y=log(FITC.H)),col="day12")+
    geom_point(data=file12)+
    xlim(2,8)+
    ylim(1,8)+
    geom_point(data=EE24_8E_B3g.txt,aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="ancestral"), alpha=0.2)
  png(paste("plot_", filelist_12A30[i], "both.png", sep = ""), width=600, height=500, res=120) # start export
  print(p) 
  dev.off() # finish export
}

## FOR B - * * 30 * *  - 8d
#YFP versus CFP - save ggplots for 96 well in a loop -filelist_8B30, well E1 (8B30) ancestral comparison (pin)

for (i in 1:length(filelist_8B30)) {
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190902_EE24/ScaleVal_EE24.8B30_gated")
  file8 <- read.csv(filelist_8B30[i], header=TRUE, sep="\\t") 
  p <- ggplot(file4,aes(x=log(Pacific.Blue.H),y=log(FITC.H)))+
    xlim(2,8)+
    ylim(1,8)+
    geom_point(data=file8, aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="8d"), alpha=0.2)+
    geom_point(data=EE24_8B30_D1g.txt,aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="ancestral"),alpha=0.05)+
    png(paste("plot_", filelist_8B30[i], "both.png", sep = ""), width=600, height=500, res=120) # start export
  print(p) 
  dev.off() # finish export
}


## FOR B - * * 30 * *   -12d
#YFP versus CFP - save ggplots for 96 well in a loop -filelist_8B30, well E1 (8B30) ancestral comparison (pin)
i=1
for (i in 1:length(filelist_12B30)) {
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190606_FACS/ScaleVal_EE24.12B30_gated")
  file12 <- read.csv(filelist_12B30[i], header=TRUE, sep="\\t") 
  p <- ggplot(file4,aes(x=log(Pacific.Blue.H),y=log(FITC.H)))+
    xlim(2,8)+
    ylim(1,8)+
    geom_point(data=file12, aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="12d"), alpha=0.2)+
    geom_point(data=EE24_12B30_D1g.txt,aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="ancestral"),alpha=0.05)+
    png(paste("plot_", filelist_8B30[i], "both.png", sep = ""), width=600, height=500, res=120) # start export
  print(p) 
  dev.off() # finish export
}


## FOR B - * * IS * *  - 8d
#YFP versus CFP - save ggplots for 96 well in a loop -filelist_8B30, well E1 (8B30) ancestral comparison (pin)

for (i in 1:length(filelist_12BIS)) {
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/20190606_FACS/ScaleVal_EE24.12BIS_gated")
  file12 <- read.csv(filelist_12BIS[i], header=TRUE, sep="\\t") 
  p <- ggplot(file4,aes(x=log(Pacific.Blue.H),y=log(FITC.H)))+
    xlim(2,8)+
    ylim(1,8)+
    geom_point(data=file12, aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="12d"), alpha=0.2)+
    geom_point(data=EE24_12BIS_H8g.txt,aes(x=log(Pacific.Blue.H),y=log(FITC.H), col="ancestral"),alpha=0.05)+
    png(paste("plot_", filelist_12BIS[i], "both.png", sep = ""), width=600, height=500, res=120) # start export
  print(p) 
  dev.off() # finish export
}


#multiplot
library("Rmisc")
#CFP - save ggplots for 96 well in a loop -filelist_8C30 
for (i in 1:length(filelist_4C30)) {
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.4C30_gated")
  file4 <- read.csv(filelist_4C30[i], header=TRUE, sep="\\t") 
  setwd("/Users/itomanek/Documents/promoter_evolution/experiments/FACS_data/EE24_20190805/ScaleVal_EE24.8C30_gated")
  file8 <- read.csv(filelist_8C30[i], header=TRUE, sep="\\t") 
  p1 <- ggplot(file4,aes(log(Pacific.Blue.H)))+
    xlim(1.5,6.5)+
    ylim(0,1.5)+
    geom_density(data=file4, aes(log(FITC.H), col="4"))+
    geom_density(data=file8, aes(log(FITC.H), col="8"))
  p2 <- ggplot(file4,aes(log(Pacific.Blue.H)))+
    xlim(2.5,6.5)+
    ylim(0,1.5)+
    geom_density(data=file4, aes(log(Pacific.Blue.H), col="4"))+
    geom_density(data=file8, aes(log(Pacific.Blue.H), col="8"))
  p<-multiplot(p1,p2,cols=1)
  png(paste("plot_", filelist_8C30[i], ".png", sep = ""), width=600, height=500, res=120) # start export
  print(p) 
  dev.off() # finish export
}



# YFP-only amplifications show tailed YFP distribution characteristic of amplification
## medium B - IS high YFP+ fraction for qPCR   #y...YFP+, c...CFP+
##wells A5y, B2y, G1y, B9c, F12a (looks ancestral in pin)
ggplot(EE24_12BIS_B12g.txt,aes(log(FITC.H)))+
  geom_density(data=EE24_12BIS_A5g.txt, aes(log(FITC.H), col="A5"))+
  geom_density(data=EE24_12BIS_B2g.txt, aes(log(FITC.H), col="B2"))+
  geom_density(data=EE24_12BIS_G1g.txt, aes(log(FITC.H), col="G1"))+
  geom_density(data=EE24_12BIS_B9g.txt, aes(log(FITC.H), col="B9cfp"))+
  geom_density(data=EE24_12BIS_F12g.txt, aes(log(FITC.H), fill="ancestral"),col=0, alpha=0.2)

# medium A
# medium A and medium E - Same YFP! 
plot(log(EE24_12A30_A2g.txt$Pacific.Blue.H),log(EE24_12A30_A2g.txt$FITC.H)) 
points(log(EE24_8E_B2g.txt$Pacific.Blue.H),log(EE24_8E_B2g.txt$FITC.H), col=" light blue") 







## @Figure 3 
#12B30 -A10,B2, B8, E5, H3 -- all mixed fraction pops (plot individually)
plot=ggplot(EE24_12B30_A10g.txt,aes(log(Pacific.Blue.H), log(FITC.H)))+
 xlim(1,8)+ylim(1,8)+
  scale_color_manual(values="black")+ theme(plot.title = element_text(size=10))+
  theme_bw()+ theme(legend.position = "none")+theme(axis.title.x = element_blank(),axis.title.y=element_blank())+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15)) 
  A10=plot+geom_point(data=EE24_12B30_A10g.txt, aes(log(Pacific.Blue.H),log(FITC.H), col="A10"), alpha=0.2)+ggtitle("A10")+
    geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)),col="grey", alpha=0.2)
  B2=plot+geom_point(data=EE24_12B30_B2g.txt, aes(log(Pacific.Blue.H),log(FITC.H), col="B2"), alpha=0.2)+ggtitle("B2")+
    geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)),col="grey", alpha=0.2)
  B8=plot+geom_point(data=EE24_12B30_B8g.txt, aes(log(Pacific.Blue.H),log(FITC.H), col="B8"), alpha=0.2)+ggtitle("B8")+
    geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)),col="grey", alpha=0.2)
  E5=plot+geom_point(data=EE24_12B30_E5g.txt, aes(log(Pacific.Blue.H),log(FITC.H), col="E5"), alpha=0.2)+ggtitle("E5")+
    geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)),col="grey", alpha=0.2)
  H3=plot+geom_point(data=EE24_12B30_H3g.txt, aes(log(Pacific.Blue.H),log(FITC.H), col="H3"), alpha=0.2)+ggtitle("H3")+
    geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)),col="grey", alpha=0.2)
multiplot(A10,B2,B8,E5,H3,cols=3)

#12A30 -D11,F3,G7,H1 -- all mixed fraction pops (plot individually)
plot=ggplot(EE24_12B30_A10g.txt,aes(log(Pacific.Blue.H), log(FITC.H)))+
  xlim(1,8)+ylim(1,8)+ theme(plot.title = element_text(size=10))
  scale_color_manual(values="black")+
  geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)),col="grey", alpha=0.2)+
  theme_bw()+ theme(legend.position = "none")+ theme(axis.title.x = element_blank(),axis.title.y=element_blank())+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15)) 
D11=plot+geom_point(data=EE24_12A30_D11g.txt, aes(log(Pacific.Blue.H),log(FITC.H), col="D11"), alpha=0.2)+ggtitle("D11")+
  geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)),col="grey", alpha=0.2)
F3=plot+geom_point(data=EE24_12A30_F3g.txt, aes(log(Pacific.Blue.H),log(FITC.H), col="F3"), alpha=0.2)+ggtitle("F3")+
  geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)),col="grey", alpha=0.2)
G7=plot+geom_point(data=EE24_12A30_G7g.txt, aes(log(Pacific.Blue.H),log(FITC.H), col="G7"), alpha=0.2)+ggtitle("G7")+
  geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)),col="grey", alpha=0.2)
H1=plot+geom_point(data=EE24_12A30_H1g.txt, aes(log(Pacific.Blue.H),log(FITC.H), col="H1"), alpha=0.2)+ggtitle("H1")+
  geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)),col="grey", alpha=0.2)
multiplot(D11,F3,G7,H1,cols=3)


## supplement 
#same with multiplot no legends and labels (add in inkscape)
p=ggplot(EE24_12B30_A10g.txt,aes(log(Pacific.Blue.H), log(FITC.H)))+
  xlim(1,8)+ylim(1,8)+  theme_bw()+ theme(legend.position = "none")+theme(axis.title.x=element_blank(),axis.title.y=element_blank())+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))+ theme(plot.title = element_text(size=10))
ctr=geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="grey", alpha=0.2)

a=p+ geom_point(data=EE24_12A30_D11g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr
b=p+geom_point(data=EE24_12A30_F3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr
c=p+geom_point(data=EE24_12A30_G7g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr
d=p+geom_point(data=EE24_12A30_H1g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr


#same with multiplot no legends and labels (add in inkscape)
ab=p+ geom_point(data=EE24_12B30_B2g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr
bb=p+geom_point(data=EE24_12B30_B8g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr
cb=p+geom_point(data=EE24_12B30_H3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr
db=p+geom_point(data=EE24_12B30_E5g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr
eb=p+geom_point(data=EE24_12B30_A10g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr
empty=p+geom_point(data=EE24_12B30_A10g.txt, aes(log(1),log(1), col="H1"), alpha=0.2)


#library("Rmisc")
multiplot(cols=3,a,b,c,d) #use 
multiplot(cols=3,ab,bb,cb,db,eb)#use 
multiplot(cols=5,ab,bb,cb,db,eb,a,b,c,d)
multiplot(cols=5,ab,bb,cb,db,eb,empty,a,b,c,d)


#@Figure 2C 
#plot with 3 FACS FRACTIONS
#same with multiplot no legends and labels (add in inkscape)
p=ggplot(EE24_12B30_A10g.txt,aes(log(Pacific.Blue.H), log(FITC.H)))+
  xlim(1,8)+ylim(1,8)+  theme_bw()+ theme(legend.position = "none")+theme(axis.title.x=element_blank(),axis.title.y=element_blank())+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))
ctr=geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), alpha=0.2,col="grey")

middle=p+ geom_point(data=EE24_12B30_A10g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr
right=p+geom_point(data=EE24_12A30_A5g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr
left=p+geom_point(data=EE24_12BIS_B2g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.2)+ctr

multiplot(cols=1,left,middle,right)

#plot medium C increased wells/mixed pops
p=ggplot(EE24_12B30_A10g.txt,aes(log(Pacific.Blue.H)))+
 theme_bw()+ theme(legend.position = "none")+theme(axis.title.x=element_blank(),axis.title.y=element_blank())+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))
  #xlim(1,8)+ylim(1,8)
ctr=geom_density(data=EE24_8E_B3g.txt, aes(log(FITC.H)), col="dark grey", alpha=0.2)
ctrCFP=geom_density(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H)), col="dark grey", alpha=0.2)
#mixed pop: B1
YFP_B1=p+ geom_density(data=EE24_4C30_B1g.txt,aes(log(FITC.H)),col="black")+ geom_density(data=EE24_12C30_B1g.txt,aes(log(FITC.H)),col="purple")+ctr+
  geom_density(data=EE24_8C30_B1g.txt,aes(log(FITC.H)),col="dark blue")
CFP_B1=p+geom_density(data=EE24_4C30_B1g.txt,aes(log(Pacific.Blue.H)),col="black")+geom_density(data=EE24_12C30_B1g.txt,aes(log(Pacific.Blue.H)),col="purple")+ctrCFP+
  geom_density(data=EE24_8C30_B1g.txt,aes(log(Pacific.Blue.H)),col="dark blue")
multiplot(YFP_B1,CFP_B1)  #270 250 size export
#nice amplification (as an illustration): B3
YFP_B3=p+ geom_density(data=EE24_4C30_B3g.txt,aes(log(FITC.H)),col="black")+ geom_density(data=EE24_12C30_B3g.txt,aes(log(FITC.H)),col="purple")+ctr+
  geom_density(data=EE24_8C30_B3g.txt,aes(log(FITC.H)),col="dark blue")
CFP_B3=p+geom_density(data=EE24_4C30_B3g.txt,aes(log(Pacific.Blue.H)),col="black")+geom_density(data=EE24_12C30_B3g.txt,aes(log(Pacific.Blue.H)),col="purple")+ctrCFP+
  geom_density(data=EE24_8C30_B3g.txt,aes(log(Pacific.Blue.H)),col="dark blue")
multiplot(YFP_B3,CFP_B3)  #270 250 size export
 
## scatter plot (this is a mixed pop)
p=ggplot(EE24_12C30_B1g.txt,aes(log(Pacific.Blue.H), log(FITC.H)))+
  xlim(1,8)+ylim(1,8)+  theme_bw()+ theme(legend.position = "none")+theme(axis.title.x=element_blank(),axis.title.y=element_blank())+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))
ctr=geom_point(data=EE24_8E_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="grey", alpha=0.3)
amp=geom_point(data=EE24_12C30_B3g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="black", alpha=0.7, shape=1)
B1=geom_point(data=EE24_12C30_B1g.txt, aes(log(Pacific.Blue.H),log(FITC.H)), col="magenta", alpha=0.7)
p+amp+B1+ctr #Figure 5a

multiplot((p+amp+ctr),(p+B1+ctr),col=2)  #270 250 size export


## @Figure 4D
### 12C30 check outliers C9&F6
# - plot medium C increased wells/mixed pops
p=ggplot(EE24_12C30_F6g.txt,aes(log(Pacific.Blue.H)))+
  theme_bw()+ theme(legend.position = "none")+theme(axis.title.x=element_blank(),axis.title.y=element_blank())+ #no grey background
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),text=element_text(size=15))
#xlim(1,8)+ylim(1,8)
ctr=geom_density(data=EE24_8E_B3g.txt, aes(log(FITC.H)), col="dark grey", alpha=0.2)
F6=geom_density(data=EE24_12C30_F6g.txt, aes(log(FITC.H)), col="blue", alpha=0.2)
C9=geom_density(data=EE24_12C30_C9g.txt, aes(log(FITC.H)), col="dark green", alpha=0.2)
F6_8=geom_density(data=EE24_8C30_F6g.txt, aes(log(FITC.H)), col="light blue", alpha=0.2)
C9_8=geom_density(data=EE24_8C30_C9g.txt, aes(log(FITC.H)), col=" green", alpha=0.2)
F6_4=geom_density(data=EE24_4C30_F6g.txt, aes(log(FITC.H)), col=" red", alpha=0.2)

#density plot YFP
p+ctr+F6+F6_8+F6_4 #increased yfp
p+ctr+C9+C9_8

#density plot YFP/CFP
ctr=geom_density(data=EE24_8E_B3g.txt, aes(log(FITC.H/Pacific.Blue.H)), col="dark grey", alpha=0.2)
F6_12=geom_density(data=EE24_12C30_F6g.txt, aes(log(FITC.H/Pacific.Blue.H)), col="blue", alpha=0.2)
F6_8=geom_density(data=EE24_8C30_F6g.txt, aes(log(FITC.H/Pacific.Blue.H)), col="light blue", alpha=0.2)
F6_4=geom_density(data=EE24_4C30_F6g.txt, aes(log(FITC.H/Pacific.Blue.H)), col=" red", alpha=0.2)

C9_12=geom_density(data=EE24_12C30_C9g.txt, aes(log(FITC.H/Pacific.Blue.H)), col="dark green", alpha=0.2)
C9_8=geom_density(data=EE24_8C30_C9g.txt, aes(log(FITC.H/Pacific.Blue.H)), col=" green", alpha=0.2)
p+ctr+F6_8+F6_12+F6_4
p+ctr+C9_8+C9_12  


