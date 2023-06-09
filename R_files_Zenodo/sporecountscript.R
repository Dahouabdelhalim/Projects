##Script for analysis of spore count data######################################
###########################################################################################################
###INITIALIZE PACKAGES AND SUCH############################################################################

rm(list=ls()) # clear memory and start fresh
#Load packages into memory (will need to be downloaded and installed first)
#This should only be necessary once per computer.
#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("dplyr")
#install.packages("readr")
#install.packages("lme4")
#install.packages("cowplot")
#install.packages("emmeans")
#install.packages("lmerTest")

library(ggplot2) #For plotting
library(reshape2)  #For converting data between wide and long forms
library(plyr)
library(dplyr)  #For grouping and manipulating data
library(readr)  #For importing CSV files
library("lme4")

#Set working directory
setwd("C:/Users/tlars/Google Drive/School/QSLab/ExperimentalEvolution/submission_121220(Evolution_Briefcommunications)/datafordryad/")

########################################################################################################
######################################Question 1########################################################
########################################################################################################


#Question 1:  does evolution in isolation of Burkholderia change its effects on Dictyostelium's spore production?
#For each strain, clade, and for all strains together, compare spore production when the ancestral Dicty is paired with Banc and when paired with any of the evolved lines
#To account for strains having different intrinsic spore production, must normalize each one by the spore production with no Burk present and express the others as changes in proportion to that
#Put another way, this is using spore counts to look at evolution in Burkholderia while leaving Dicty constant

#######################Import data#######################################################

#Import data
workingdata <- read_csv("sporecountdata.csv")

#Make sure everything that should be a factor is, and that factors are in the desired order that they will appear in the final graph
workingdata$Date<-factor(workingdata$Date)
workingdata$Dstrain<-factor(workingdata$Dstrain, levels=c("QS70","QS159","QS161","QS11","QS21","QS69"))
workingdata$Dline<-factor(workingdata$Dline, levels=c("None","Anc","E1","E2","E3"))
workingdata$Bstrain<-factor(workingdata$Bstrain, levels=c("bQS70","bQS159","bQS161","bQS11","bQS21","bQS69"))
workingdata$Bline<-factor(workingdata$Bline, levels=c("None","Anc","E1","E2","E3"))
workingdata$Dclade<-factor(workingdata$Dclade, levels=c("B1host","B2host"))
workingdata$Bclade<-factor(workingdata$Bclade, levels=c("B1","B2"))
workingdata$Dstatus<-factor(workingdata$Dstatus, levels=c("None","Anc","Evolved"))
workingdata$Bstatus<-factor(workingdata$Bstatus, levels=c("None","Anc","Evolved"))
workingdata$Note[which(is.na(workingdata$Note))]<-'none'


######################Calculations######################################################################################

workingdata$Average<-0 #initialize a new column for the average spore production
workingdata$Averageanc<-0
#For each strain, calculate the average of the rows with no Burk present and put it into the Average column for all rows matching that strain
for (strain in unique(workingdata$Dstrain)) {
  workingdata$Average[which(workingdata$Dstrain==strain)]<-mean(workingdata$Totalspores[which(workingdata$Dstrain==strain & workingdata$Bstatus=="None")])
}

#Calculate fold change by dividing Totalspores by Average
workingdata$Foldchange<-workingdata$Totalspores/workingdata$Average

#For each strain, calculate the average of the rows with the ancestor Burk present and put it into the Averageanc column for all rows matching that strain
for (strain in unique(workingdata$Dstrain)) {
  workingdata$Averageanc[which(workingdata$Dstrain==strain)]<-mean(workingdata$Totalspores[which(workingdata$Dstrain==strain & workingdata$Bstatus=="Anc")])
}

#Calculate fold change relative to Banc by dividing Totalspores by Averageanc
workingdata$Foldchangeanc<-workingdata$Totalspores/workingdata$Averageanc


#####################Graphing#########################################################
###################################################################

#remove the control lines
testdata<-workingdata[which(workingdata$Bstatus!="None"),]
#Calculate averages
testdata$Totalspores.repavg<-ave(testdata$Totalspores, testdata$Dstrain,testdata$Bstrain,testdata$Dline,testdata$Bline)
testdata$Foldchange.repavg<-ave(testdata$Foldchange, testdata$Dstrain,testdata$Bstrain,testdata$Dline,testdata$Bline)
testdata$Foldchangeanc.repavg<-ave(testdata$Foldchangeanc, testdata$Dstrain,testdata$Bstrain,testdata$Dline,testdata$Bline)

testdata2.anc<-testdata[which(testdata$Bstatus=="Anc"),]
testdata2.evo<-testdata[which(testdata$Bstatus!="Anc"),]
testdata2.evo$Totalspores.repavg<-ave(testdata2.evo$Totalspores, testdata2.evo$Dstrain,testdata2.evo$Bstrain,testdata2.evo$Dline,testdata2.evo$Bline)
testdata2.evo$Foldchange.repavg<-ave(testdata2.evo$Foldchange, testdata2.evo$Dstrain,testdata2.evo$Bstrain,testdata2.evo$Dline,testdata2.evo$Bline)
testdata2.evo$Foldchangeanc.repavg<-ave(testdata2.evo$Foldchangeanc, testdata2.evo$Dstrain,testdata2.evo$Bstrain,testdata2.evo$Dline,testdata2.evo$Bline)
testdata2.evo.collapsed<-testdata2.evo %>% distinct(Totalspores.repavg, .keep_all = TRUE)
testdata2.evo.collapsed$Totalspores<-testdata2.evo.collapsed$Totalspores.repavg
testdata2.evo.collapsed$Foldchange<-testdata2.evo.collapsed$Foldchange.repavg
testdata2.evo.collapsed$Foldchangeanc<-testdata2.evo.collapsed$Foldchangeanc.repavg

testdata3<-rbind(testdata2.anc,testdata2.evo.collapsed[,c(1:18)])

#Make graph
labels<-c(B1host="hosts of P. agricolaris",B2host="hosts of P.hayleyella")
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
xlabels<-c("+ancestral\\nParaburkholderia","+evolved\\nParaburkholderia")
ylabels<-c(0,expression("1x10"^8),expression("2x10"^8),expression("3x10"^8),expression("4x10"^8),expression("5x10"^8),expression("6x10"^8))
plot1e<-ggplot(data=testdata3, aes(x=Bstatus, y=Totalspores, fill=Bstrain, color=Bstrain))+
  geom_dotplot(binaxis='y',stackdir='center',position=position_dodge(width=.3), stackratio=1,dotsize=10000000,binwidth=1) +
  stat_summary(fun.y=mean,geom="line",aes(group=Bstrain, color=Bstrain),position=position_dodge(width=.3)) + 
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  scale_y_continuous(limits=c(0,600000000),breaks=seq(0,600000000, by=100000000), labels=ylabels) +
  scale_x_discrete(labels=xlabels) + 
  facet_grid(.~Dclade,switch='x',labeller=labeller(Dclade=labels)) +
  theme(plot.margin=margin(t=.1,l=.030,r=.030,b=.1, unit="npc"),
        axis.title.y=element_text(size=18,vjust=3),
        plot.title=element_text(size=16, hjust=5),
        axis.text.x=element_text(color="black",size=11),
        axis.text.y=element_text(color="black",size=12),
        strip.text.x=element_text(color="black",size=12),
        legend.title=element_blank(),
        legend.text=element_text(color="black",size=14),
        legend.key=element_rect(color="transparent", fill=alpha("red",0)),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_rect(color="black",fill="gray"),
        panel.background=(element_rect(color="black",fill="white"))) +
  labs(y=expression(paste(italic("D. discoideum")," spore production", sep="")),x="")
  #ggtitle("Evolution in isolation increases P. hayleyella's\\ntoxicity to D. discoideum hosts")
plot1e
ggsave("figure2.tiff", height=7, width=8)

##Supplementary figure 1 - to show absolute spore counts
absolutedata<-workingdata[which(workingdata$Bstatus!="Evolved"),]
xlabels<-c("P-","P+")
ylabels<-c(0,expression("1x10"^8),expression("2x10"^8),expression("3x10"^8),expression("4x10"^8),expression("5x10"^8),expression("6x10"^8),expression("7x10"^8),expression("8x10"^8))
plot1c<-ggplot(data=absolutedata, aes(x=Bstatus, y=Totalspores,fill=Dstrain,color=Dstrain))+
  geom_dotplot(binaxis='y',stackdir='center',stackratio=.5,dotsize=1500000000,binwidth=.01) +
  stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.25, color="black") +
  scale_y_continuous(limits=c(0,800000000),breaks=seq(0,800000000, by=100000000),labels=ylabels) +
  scale_x_discrete(labels=xlabels) + 
  facet_grid(.~Dstrain,switch='both',labeller=label_value) +
  theme(plot.margin=margin(t=.1,l=.030,r=.030,b=.1, unit="npc"),
        panel.spacing.x=unit(c(.01,.01,.05,.01,.01), "npc"), 
        axis.title.y=element_text(size=18),
        plot.title=element_text(size=14, hjust=.5),
        axis.text.x=element_text(color="black",size=12),
        axis.text.y=element_text(color="grey60",size=12),
        strip.text.x=element_text(color="black",size=12),
        legend.position='none',
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background=element_rect(color="black",fill="gray"),
        panel.background=(element_rect(color="black",fill="white"))) +
  labs(y=expression(paste(italic("D. discoideum")," spore production", sep="")),x="")
ggtitle("Absolute spore production in infected and uninfected Dictyostelium")
ggdraw(plot1c) + 
  draw_label("Hosts of", x=.25, y=.085, size=16) +
  draw_label("Hosts of", x=.67, y=.085, size=16) +
  draw_label("P. agricolaris", x=.41, y=.085, size=16, fontface="italic") +
  draw_label("P. hayleyella", x=.83, y=.085, size=16, fontface="italic") +
ggsave("suppfigure1_040520.png", height=7, width=7)


#####################Statistics#########################################################


newmodel<-lmer(Foldchange ~ Bstatus + (1|Bstrain),data=workingdata, REML=FALSE)
newmodel
coef(newmodel)
summary(newmodel)

newmodel2<-lmer(Foldchange ~ Bstatus*Bclade + (1|Bstrain),data=workingdata, REML=FALSE)
newmodel2
coef(newmodel2)
summary(newmodel2)

newmodel3<-lmer(Foldchange ~ Bstatus*Bclade + (1+Bstatus|Bstrain),data=workingdata, REML=FALSE)
summary(newmodel3)
coef(newmodel3)
summary(newmodel3)

newmodel4<-lmer(Foldchange ~ Bstatus, data=workingdata)


nullmodel<-lmer(Foldchange ~ 1 + (1|Bstrain),data=workingdata, REML=FALSE)
nullmodel2<-lmer(Foldchange ~ Bstatus+Bclade + (1|Bstrain), data=workingdata, REML=FALSE) #this is for checking model 2 against an additive version

#Test some models against one another
anova(nullmodel,newmodel) #including status vs excluding it
anova(nullmodel,newmodel2) #including status*clade vs excluding both
anova(nullmodel2,newmodel2) #including multiplicative interaction vs only additive
anova(newmodel,newmodel2) #Including clade vs excluding it
anova(newmodel2,newmodel3) #Random intercept or not
anova(newmodel3,newmodel4) #Skipping the random effects part entirely

#So the model you want is Bstatus*Bclade + (1|Bstrain) (newmodel2).  This suggests that clade and status are both important.