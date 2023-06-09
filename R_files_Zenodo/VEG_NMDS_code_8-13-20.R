library(extrafont)
library(vegan)
library(gridExtra)
library(ggplot2)
install.packages("pairwise.adonis") 
library(pairwise.adonis)

############################################ Early  2015 Veg NMDS  #######################################################

all15 = read.table("early15.csv",header=T,sep=",")
names(all15) ; dim(all15)

vegall15 = all15[,3:66] ; names(vegall15) 
names(vegall15) ; dim(vegall15)
colSums(vegall15) ; rowSums(vegall15)


#which distance matrix is best (largest value is better)
dummyVar = as.data.frame(all15$Treatment)
sort(rankindex(scale(dummyVar),vegall15), decreasing=T) # gow

# Default distance is bra
early15K2=metaMDS(vegall15, k=2, distance="gow")
early15K3=metaMDS(vegall15, k=3, distance="gow")

early15K2
early15K3
#high R2 is good here
stressplot(early15K2)
stressplot(early15K3)

#create dummy variables
treat = all15$Treatment ; treat
as.numeric(treat)

#NMDS customizable graph
MDS1 = early15K3$points[,1]
MDS2 = early15K3$points[,2]
NMDSplotE15 = data.frame(MDS1 = as.numeric(MDS1), MDS2 = as.numeric(MDS2), treat = as.factor(treat))

Early2015<-ggplot(NMDSplotE15, aes(x=MDS1, y=MDS2))+ theme_minimal() +
  scale_fill_identity()+stat_ellipse(aes(color = factor(treat)),lwd=2)+
  scale_color_manual(values=alpha(c("#000000", "#0072B2", "#009E73", 
                                    "#D55E00", "#56B4E9", "#CC79A7", "#F0E442"),0.80),labels=c("Vacant Lot", "Meadow", "Fescue Lawn","Flowering Lawn", "Grass Prairie", "Low Diversity Prairie", 
                                                                                               "High Diversity Prairie"))+
  theme(legend.position="right", legend.text = element_text(size=14,family="Calibri"),
        legend.title=element_blank(),legend.key.height=unit(1.3,"line"))+
  labs(title = "A. Early 2015 Vegetation\\n", x="\\nMDS1", y="MDS2\\n") +
  theme(plot.title=element_text(size=20,face="bold",family="Calibri"),
        axis.text=element_text(size=16,family="Calibri",lineheight=.5),
        axis.title=element_text(size=18,family="Calibri",lineheight=.5))
Early2015

#assumption tests

#there is NO difference among treatments AND season. 
bdvegall15 = betadisper(vegdist(vegall15,dist="gow"), treat)
#Visualize dispersion test (variance among groups)
boxplot(bdvegall15)

#Test variance among groups. Null Hypo: variance is equal 
anova(bdvegall15)
permutest(bdvegall15)
adonis2(vegdist(vegall15,"gow") ~ factor(treat), method= "gow")
pairwise.adonis(vegall15, factor(treat), sim.method='gower')


################################ LATE  2015 VEG NMDS  #######################################################

dataL15.0 = read.table("late15.csv",header=T,sep=",")
names(dataL15.0) ; dim(dataL15.0)
vegdataL15.0 = dataL15.0[,3:63] ; names(vegdataL15.0)
colSums(vegdataL15.0) ; rowSums(vegdataL15.0) # remove row 16

dataL15 = dataL15.0[-16,]
names(dataL15) ; dim(dataL15)
vegdataL15 = dataL15[,3:63] ; names(vegdataL15)
colSums(vegdataL15) ; rowSums(vegdataL15)

#figure out which distance matrix is best (largest value is better)
dummyVar = as.data.frame(dataL15$Treatment)
sort(rankindex(scale(dummyVar),vegdataL15), decreasing=T) # GOWER

# Default distance is bra
late15K2=metaMDS(vegdataL15, k=2, distance="gow")
late15K3=metaMDS(vegdataL15, k=3, distance="gow")

#your aim is to lessen the stress which you can check in the output of the following
late15K2
late15K3
#high R2 is good here
stressplot(late15K2)
stressplot(late15K3)

#create dummy variables
treat = dataL15$Treatment ; treat
as.numeric(treat)

#NMDS customizabe plot
MDS1 = late15K3$points[,1]
MDS2 = late15K3$points[,2]
NMDSplotL15 = data.frame(MDS1 = as.numeric(MDS1), MDS2 = as.numeric(MDS2), treat = as.factor(treat))

Late2015<-ggplot(NMDSplotL15, aes(x=MDS1, y=MDS2))+ theme_minimal() +
  scale_fill_identity()+stat_ellipse(aes(color = factor(treat)),lwd=2)+
  scale_color_manual(values=alpha(c("#000000", "#0072B2", "#009E73", 
                                    "#D55E00", "#56B4E9", "#CC79A7", "#F0E442"),0.80),labels=c("Vacant Lot", "Meadow", "Fescue Lawn","Flowering Lawn", "Grass Prairie", "Low Diversity Prairie", 
                                                                                               "High Diversity Prairie"))+
  theme(legend.position="right", legend.text = element_text(size=14,family="Calibri"),
        legend.title=element_blank(),legend.key.height=unit(1.3,"line"))+
  labs(title = "B. Late 2015 Vegetation\\n", x="\\nMDS1", y="MDS2\\n") +
  theme(plot.title=element_text(size=20,face="bold",family="Calibri"),
        axis.text=element_text(size=16,family="Calibri",lineheight=.5),
        axis.title=element_text(size=18,family="Calibri",lineheight=.5))
Late2015

#check assumptions

#there is NO difference among treatments AND season. 
bdvegdataL15 = betadisper(vegdist(vegdataL15,dist="gow"), treat)
#Visualize dispersion test (variance among groups)
boxplot(bdvegdataL15)

#Test variance among groups. Null Hypo: variance is equal 
anova(bdvegdataL15)
permutest(bdvegdataL15)
adonis2(vegdist(vegdataL15,"gow") ~ factor(treat), method= "gow")
pairwise.adonis(vegdataL15, factor(treat), sim.method="gow")

####################################   EARLY Veg 2016  #######################################################

all16.0 = read.table("early16.csv",header=T,sep=",")
names(all16.0) ; dim(all16.0)

vegall16.0 = all16.0[,3:52] ; names(vegall16.0) #
colSums(vegall16.0) ; rowSums(vegall16.0)
#remove empty rows 16 and 21 from the original data and start over

all16 = all16.0[-c(16,21),]
vegall16 = all16[,3:52] ; names(vegall16)
colSums(vegall16) ; rowSums(vegall16)

#figure out which distance matrix is best (larger value is better)
dummyVar = as.data.frame(all16$Treatment)
sort(rankindex(scale(dummyVar),vegall16), decreasing=T) # gow is 0.12, bra is 0.10
#chose default distance bra

early16K2=metaMDS(vegall16, k=2, distance="bra")
early16K3=metaMDS(vegall16, k=3, distance="bra")

#your aim is to lessen the stress which you can check in the output of the following
early16K2
early16K3
#high R2 is good here
stressplot(early16K2)
stressplot(early16K3)

#customizable NMDS plot
treat = all16$Treatment ; treat
as.numeric(treat)

MDS1 = early16K3$points[,1]
MDS2 = early16K3$points[,2]
NMDSplotE16 = data.frame(MDS1 = as.numeric(MDS1), MDS2 = as.numeric(MDS2), treat = as.factor(treat))

Early2016<-ggplot(NMDSplotE16, aes(x=MDS1, y=MDS2))+ theme_minimal() +
  scale_fill_identity()+stat_ellipse(aes(color = factor(treat)),lwd=2)+
  scale_color_manual(values=alpha(c("#000000", "#0072B2", "#009E73", 
                                    "#D55E00", "#56B4E9", "#CC79A7", "#F0E442"),0.80),labels=c("Vacant Lot", "Meadow", "Fescue Lawn","Flowering Lawn", "Grass Prairie", "Low Diversity Prairie", 
                                                                                               "High Diversity Prairie"))+
  theme(legend.position="right", legend.text = element_text(size=14,family="Calibri"),
        legend.title=element_blank(),legend.key.height=unit(1.3,"line"))+
  labs(title = "C. Early 2016 Vegetation\\n", x="\\nMDS1", y="MDS2\\n") +
  theme(plot.title=element_text(size=20,face="bold",family="Calibri"),
        axis.text=element_text(size=16,family="Calibri",lineheight=.5),
        axis.title=element_text(size=18,family="Calibri",lineheight=.5))

Early2016

#check assumptions

#there is NO difference among treatments AND season. 
bdvegall16 = betadisper(vegdist(vegall16,dist="bra"), treat)
#Visualize dispersion test (variance among groups)
boxplot(bdvegall16)

#Test variance among groups. Null Hypo: variance is equal
anova(bdvegall16)
permutest(bdvegall16)
adonis2(vegdist(vegall16,"bray") ~ factor(treat), method= "bray")
pairwise.adonis(vegall16, factor(treat), sim.method="bray")

######################################## LATE  2016 VEG NMDS #######################################################

dataL16.0 = read.table("late16.csv",header=T,sep=",")
names(dataL16.0) ; dim(dataL16.0)

vegdataL16.0 = dataL16.0[,3:60] ; names(vegdataL16.0) #
colSums(vegdataL16.0) ; rowSums(vegdataL16.0)
#remove empty rows 16 and 21 from the original data and start over

dataL16 = dataL16.0[-16,]
vegdataL16 = dataL16[,3:52] ; names(vegdataL16)
colSums(vegdataL16) ; rowSums(vegdataL16)

#figure out which distance matrix is best (largest value is better)
dummyVar = as.data.frame(dataL16$Treatment)
sort(rankindex(scale(dummyVar),vegdataL16), decreasing=T) # went with bra for comparison to early 2016

# Default distance is bra
late16K2=metaMDS(vegdataL16, k=2, distance="bra")
late16K3=metaMDS(vegdataL16, k=3, distance="bray")

#your aim is to lessen the stress which you can check in the output of the following
late16K2
late16K3
#high R2 is good here
stressplot(late16K2)
stressplot(late16K3)

treat = dataL16$Treatment ; treat
as.numeric(treat)
MDS1 = late16K3$points[,1]
MDS2 = late16K3$points[,2]
NMDSplotL16 = data.frame(MDS1 = as.numeric(MDS1), MDS2 = as.numeric(MDS2), treat = as.factor(treat))


Late2016<-ggplot(NMDSplotL16, aes(x=MDS1, y=MDS2))+ theme_minimal() +
  scale_fill_identity()+stat_ellipse(aes(color = factor(treat)),lwd=2)+
  scale_color_manual(values=alpha(c("#000000", "#0072B2", "#009E73", 
                                   "#D55E00", "#56B4E9", "#CC79A7", "#F0E442"),0.80),labels=c("Vacant Lot", "Meadow", "Fescue Lawn","Flowering Lawn", "Grass Prairie", "Low Diversity Prairie", 
                              "High Diversity Prairie"))+
  theme(legend.position="right", legend.text = element_text(size=14,family="Calibri"),
        legend.title=element_blank(),legend.key.height=unit(1.3,"line"))+
  labs(title = "D. Late 2016 Vegetation\\n", x="\\nMDS1", y="MDS2\\n") +
  theme(plot.title=element_text(size=20,face="bold",family="Calibri"),
        axis.text=element_text(size=16,family="Calibri",lineheight=.5),
        axis.title=element_text(size=18,family="Calibri",lineheight=.5))
                     
Late2016

#check assumptions late 2016     
#there is NO difference among treatments AND season. How about only about T1 vs T8?
bdvegdataL16 = betadisper(vegdist(vegdataL16,dist="bra"), treat)
#Visualize dispersion test (variance among groups)
boxplot(bdvegdataL16) # there is some variance among groups, therefore interpret with caution. significant differences may be in part due to this variance

#Test variance among groups. Null Hypo: variance is equal (proceed with adonis)
# p-val HIGHER than 0.05 then proceed. if p-val SMALLER than 0.05 then assumption is violated

anova(bdvegdataL16) #pvalue still higher than 0.05
permutest(bdvegdataL16)#pvalue still higher than 0.05
adonis2(vegdist(vegdataL16,"bra") ~ factor(treat), method= "bra")
pairwise.adonis(vegdataL16, factor(treat), sim.method="bray")

############################### compiled figure #################################

#all four graphs together
grid.arrange(Early2015,Late2015,Early2016,Late2016,ncol=2, nrow =2)
                    