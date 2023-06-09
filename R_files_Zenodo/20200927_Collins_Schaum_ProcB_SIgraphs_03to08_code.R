##### Strategies paper ProcB submission  - SI Figures #####
# This file makes the Figures in the SI , apart from the experimental design overview figures and the 'regression to mean' Figure


library(ggplot2) 
library(plyr) 
library(nlme)
library(lme4)
library(reshape2)
library(MuMIn)
library(visreg)
library(corrplot)
library(latticeExtra)
#library(Hmisc) # should work for group.CI on old R versions
library(Rmisc) # group.CI on new R versions 
library(nlsMicrobio)

setwd("~/Dropbox/Collins Lab Shared Folder/'quorum' - transiently restored/Proc B/Data files as put onto dryad")

#### SI Figure 3 Biomass production in µgC per hour as a function of Net photosynthesis (NP) in µgC produced per µgC in the sample per hour for samples evolved at ambient CO2 levels (400ppm). #####

rm(list=ls()) #remove everything that may still be in memory 

#this is the combined growth rate and PS data, pooled for techreps (standard error is negligible, nearly doesn't even show on graph)

Si03<-read.csv("MonoCoStrategies_SI_Figure_03_data.csv")
str(Si03)
head(Si03)

# from raw data growth rate to biomass calculated assume spherical shape. C contents per cell obtained through spontanoues combustion and gas chromatography. Cell counts came from flow cytometry. 
#photosynthesis and resp data converted from O2 into C using * 32 * 0.75 * (12/44))*10^-4 (the 0.75 is a conversion factor arising from C:N ratio)
#as this is NP, this is taking into consideration that photosynthesis occurs 12hrs, but resp for 24 hours under a 12:12 light:dark cycle

plotSI03<-qplot(PS.avg,mue.avg, data=Si03, shape=lineage, colour=lineage,xlab="gC produced through photosynthesis per cell and hour",ylab="µ (biomass/hour)",xlim=c(0.01,0.072),ylim=c(0.01,0.072))+
  scale_shape_manual(values=c(1:16))+
  geom_errorbar(aes(ymin=mue.avg-mue.se, ymax=mue.avg+mue.se))+
  geom_errorbarh(aes(xmin=PS.avg-PS.se,xmax=PS.avg+PS.se))+
  theme_classic(base_size = 16, base_family = 'Helvetica')+
  scale_colour_ordinal()+
  facet_wrap(~lineage,ncol=2)+
  geom_smooth(aes(group=lineage), method="lm", se=F)+
  geom_point(position="jitter")+
  geom_abline(slope = 1, intercept = 0, lty = 2)
plotSI03 # this one for SI , might just want to clean up axes titles a bit (needs to state this is NP!)

#### SI Figure 4 Biomass production in µgC per hour as a function of Net photosynthesis (NP) in µgC produced per µgC in the sample per hour for samples evolved at elevated CO2 levels (1000ppm). #####
#data processing and conversion factors as above

Si04<-read.csv("MonoCoStrategies_SI_Figure_04_data.csv")
str(Si04)
head(Si04)

plotSI04<-qplot(PS.avg,mue.avg, data=Si04, shape=lineage, colour=lineage,xlab="gC produced through photosynthesis per cell and hour",ylab="µ (biomass/hour)",xlim=c(0.01,0.1),ylim=c(0.01,0.1))+
  scale_shape_manual(values=c(1:16))+
  geom_errorbar(aes(ymin=mue.avg-mue.se, ymax=mue.avg+mue.se))+
  geom_errorbarh(aes(xmin=PS.avg-PS.se,xmax=PS.avg+PS.se))+
  theme_classic(base_size = 16, base_family = 'Helvetica')+
  scale_colour_ordinal()+
  facet_wrap(~lineage,ncol=2)+
  geom_smooth(aes(group=lineage), method="lm", se=F)+
  geom_point(position="jitter")+
  geom_abline(slope = 1, intercept = 0, lty = 2)
plotSI04 # this one for SI , might just want to clean up axes titles a bit (needs to state this is NP!), also can tweak xlim and ylim a little to get a better view 

#####SI Figure 5| Selection pCO2 affects how much of carbon produced during net photosynthesis is not directly used for growth. #####
rm(list=ls()) #remove everything that may still be in memory 

#need carbon allocation data frame
alloc<-read.csv("carballoc400_1000ppm.csv")
alloc$ppm <- factor(alloc$ppm , levels = c('400ppm', '1000ppm'))

head(alloc)
alloc$react_to[alloc$react_to == "SPIKES"] <- "SPIKE"
alloc$back<-alloc$muecomp.avg*(alloc$mue.avg*24) #per day
alloc$react_to <- factor(alloc$react_to , levels = c('COCULT', 'SPIKE',"GFP"))

#not putting facets, this is all data pooled and just checking for effect of CO2 level
qplot(ppm,surplus_percent_PS, data=alloc,geom="boxplot") +        theme_classic(base_size = 16, base_family = 'Helvetica')+ geom_jitter(aes(shape=ecotype) )+ scale_shape_manual(values=c(1:20))

#####SI Figure 6| Relative importance of carbon allocation strategies in ambient and elevated pCO2 selected lineages in the indirect (ThinCert), perceived (Spike) or direct presence of lineages from the same species complex. #####
#uses same (carbon allocation data frame) as above, so don't need to clean work space
# fills the conceptual graph in Figure 1 of the main manuscript with data according to 
#A foldchange>1 and PS %>0 (extra PS put into react more)
#B foldchange>1 and PS %<0 (use C from storage for react more)
#C foldchange<1 and PS %>0 (extra PS NOT put into react more -put somewhere else?)
#C foldchange<1 and PS %<0 (c from storage needed for growth in mono, NOT put into react more -put somewhere else?)


react <- ddply(alloc, c("ppm","react_to" ), function(df) return(c(magnitude = mean(df$muecomp.avg), std=sd(df$muecomp.avg)/sqrt(6))))
#Now introduce ABCD as above in a new column 
alloc$quadrant<-ifelse(alloc$surplus_percent_PS>0&alloc$muecomp.avg>1,"A",ifelse(alloc$surplus_percent_PS<0&alloc$muecomp.avg>1,"B",ifelse(alloc$surplus_percent_PS>0&alloc$muecomp.avg<1,"C","D")))
sub400<-subset(alloc, ppm=="400ppm")
sub1000<-subset(alloc, ppm=="1000ppm")
mean(sub400$surplus_percent_PS)-mean(sub1000$surplus_percent_PS)
#now, do an ddply, where you count the #ecotypes out of #all ecotypes that are in one or the other category. Or do a treemap??
#this is basically a glorified pie chart.  - can chose not to display lineage/ecotype identity for clarity 
library(treemapify)

ggplot(alloc, aes(area =muecomp.avg , fill = quadrant, alpha=0.55,subgroup=ecotype)) +
  geom_treemap()+
  scale_fill_manual(values=c("lightgreen","red","purple","cornflowerblue"))+
  geom_treemap_subgroup_text(place = "centre", alpha = 0.5,  colour = "black", fontface = "italic", min.size = 0) +
  geom_treemap_subgroup_border() +
  theme_classic()+
  facet_grid(react_to~ppm)

#now make a means data frame 
alloc_av<-ddply(alloc,c("ecotype","ppm","quadrant","react_to"), function(df) return(c(PS.avg=mean(df$surplus_percent_PS),PS.se=sd(df$surplus_percent_PS)/sqrt(length(alloc)),mue.avg=mean(df$muecomp.avg),mue.se=sd(df$mmuecomp.avgue)/sqrt(length(alloc)))))
alloc_av <- ddply(alloc_av, .(ppm,react_to,quadrant ), mutate, numb.q = seq_along(quadrant))
ggplot(alloc_av, aes(area =numb.q , fill = quadrant, alpha=0.55, subgroup=ecotype)) +
  geom_treemap()+
  scale_fill_manual(values=c("lightgreen","red","purple","cornflowerblue"))+
  geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5,  colour = "black", fontface = "italic", min.size = 0) +
  #geom_treemap_subgroup_border() +
  theme_classic()+
  facet_grid(react_to~ppm)

#### SI Figure 7 - these data are stored with Schaum and Collins 2014 "Plasticity predicts evolution in a marine green alga" (also on data dryad) and are displayed graohically again in this SI so that readers don't lose it over having to flip between publications. #####

#### SI Figure 8| Growth rates (µ per day calculated from cells mL-1) for lineages selected at ambient and elevated pCO2 grown in mixed culture ####
#uses same data frame as above! 'carballoc'
#just uses slightly different facetting for visibility 
qplot(ppm,back, data=alloc,geom="boxplot", facets=react_to~., ylab="Growth rate µ day -1 in mixed culture", xlab="Selection environment", fill=ppm) +        
  theme_classic(base_size = 16, base_family = 'Helvetica') +
  scale_fill_manual("Selected",values=c("seagreen","coral"))+
  facet_wrap(~react_to,scales="free_y")+theme(legend.position="top")

