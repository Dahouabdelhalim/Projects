#Cat brain size code 5.8.2021

#Libraries
library(NCmisc)
library(plyr)
library(readxl)
library(psych)
library(robust)
library(nlme)
library("pscl")
library(lme4) #Glmm
library(MASS)
library(car)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(multcomp)
library(MuMIn) #R2 calculation
library(piecewiseSEM) #BUT doesnt work on Gamma link log
library(r2glmm)
library(visreg)
library(dplyr)
library(doBy)
library(data.table)
library(GGally)
library(openxlsx)
library(ggplot2)
library(car)
library(readxl)
library(plyr)
library(stats)
source("/Users/raffaela/Documents/PhD_Classes/Stats_smallSamples20180827/diagnostic_fcns.r")
library(mclust)
library(doBy)
library(cowplot)
library(gridExtra)
library(gtable)
library(grid)
source("/Users/raffaela/Documents/PhD_Classes/Stats_smallSamples20180827/glmm_stability.r")
library(influence.ME)
library(factoextra)
library(olsrr)
library(parallel)
#install.packages("performance")
library(performance)
#install.packages("see")
library(see)
#install.packages("hbmem")
#library(hbmem)
detach("package:hbmem", unload = TRUE)
#install.packages("msm")
library(msm)
library(boot)
#install.packages("lmodel2")
library(lmodel2)
#install.packages("stargazer")
library(stargazer)
#install.packages("digitize")
library(digitize)
#install.packages("writexl")
library("writexl")
#install.packages("viridis") 
library("viridis")
#install.packages("ggalt")
library(ggalt)
#install.packages("kader")
library(kader)
#install.packages("geometry")
library(geometry)
###################################################################################################

#set working directory

setwd("/Users/raffaela/Documents/PhD_Projects/20170222_Cats/Paper/RSOS_submission")

#import data set
cat<- read_excel("RSOS_finalData_cat_data_whole_20200316.xlsx", 
                      col_types = c("text", "text", "text", 
                      "text", "text", "numeric", "numeric", 
                      "numeric", "numeric", "numeric", 
                      "numeric", "text", "text", "numeric", 
                      "numeric", "numeric", "text", "numeric", 
                      "numeric"))



#round palate length and skull length to first digit after comma
cat$ForamenMagnumMM<-round(x=cat$ForamenMagnumMM, digits=1)
cat$PalateLengthMM<-round(x=cat$PalateLengthMM, digits=1)

#convert brain measurements from g to cm3; 1ml=1.537g -> use in conversion; 1ml = 1cm3
cat$brainVolumeCM3 <- cat$BrainVolumeG/1.537

#exclude misidentified skull (GH 60.10) and skull with error in measurement (Z 2013.42.69)
#located in row number 1 and 17
cat2 <- cat[-c(1,17), ]

# save rounded file as CSV for further analysis
#CSV is more reliable to open correctly across platforms
write.csv(cat2, "RSOS_cat_brainsize_20210805.csv")

#import new csv file
Katze<- read.csv("RSOS_cat_brainsize_20210805.csv")


# C M 3  V O L U M E

fullcm3<- lm(brainVolumeCM3 ~ ForamenMagnumMM+ Hybrid_ID_LBQbased_simplified, 
              data=Katze)

nullcm3<- lm(brainVolumeCM3~ ForamenMagnumMM, 
              data=Katze)

residualPlot(fullcm3)
hist(resid(fullcm3))

#check residuals
diagnostics.plot(fullcm3)

#model summary
print(summary(fullcm3), corr=F)

#collinearity; model without random effects to check
detach("package:usdm", unload=TRUE)
library(car)

vif(fullcm3)

ID<- as.factor(Katze$Hybrid_ID_LBQbased_simplified)

plot(y=Katze$WeightG, x=ID)
#full null test
anova(nullcm3, fullcm3, test="Chisq")
overdisp.test(fullcm3)
#plot(fullcm3)
#plot model
plot(residuals(fullcm3))
summary(fullcm3)
summary(nullcm3)
hist(residuals(fullcm3), breaks=30)
summary(fullcm3)
Anova(fullcm3)

r2(fullcm3)
check_singularity(fullcm3)
check_model(fullcm3)
plot(compare_performance(fullcm3, nullcm3))
#check model stability; values under 2 are fine -> my model has 0.741
max(abs(dffits(fullcm3)))

#prep polygon clouds
cat_poly <-cbind(Katze$brainVolumeCM3,Katze$ForamenMagnumMM, Katze$Hybrid_ID_LBQbased_simplified, Katze$BrainVolumeG, Katze$PalateLengthMM)
cat_poly <- as.data.frame(cat_poly)

colnames(cat_poly)[1] <- "BrainVolumeCM3"
colnames(cat_poly)[2] <- "ForamenMagnumMMs"
colnames(cat_poly)[3] <- "HybridIDLBQbasedSimplified"
colnames(cat_poly)[4] <- "BrainVolumeGRAMS"
colnames(cat_poly)[5] <- "PalateLengthMMs"

cat_poly$BrainVolumeCM3<- as.numeric(as.character(cat_poly$BrainVolumeCM3))
cat_poly$ForamenMagnumMMs <- as.numeric(as.character(cat_poly$ForamenMagnumMMs))
cat_poly$HybridIDLBQbasedSimplified<- as.factor(cat_poly$HybridIDLBQbasedSimplified)
cat_poly$BrainVolumeGRAMS<- as.numeric(as.character(cat_poly$BrainVolumeGRAMS))
cat_poly$PalateLengthMMs<- as.numeric(as.character(cat_poly$PalateLengthMMs))

cat_poly2 <- na.omit(cat_poly)

#create polygon clouds
find_hull <- function(cat_poly2) cat_poly2[chull(cat_poly2$ForamenMagnumMMs, cat_poly2$BrainVolumeCM3), ]
hulls <- ddply(cat_poly2, "HybridIDLBQbasedSimplified", find_hull)

# convert ID to factor
cat_poly2$HybridIDLBQbasedSimplified<- as.factor(cat_poly2$HybridIDLBQbasedSimplified)

#graph
CM3_FM<- ggplot(data=cat_poly2)+
  geom_point(size=3,alpha=0.8, aes(y=BrainVolumeCM3, x=ForamenMagnumMMs, 
                                   shape=HybridIDLBQbasedSimplified, fill=HybridIDLBQbasedSimplified))+
  ylab(expression(paste(
    "cranial volume (",
    cm^3,
    ")", sep=""))) +
  xlab("basal skull length (mm)")+
  scale_shape_manual(name="cat groups", breaks= c("1", "2", "3", "4"), 
                     labels= c("domesticated", "hybrid", "lybica", "silvestris"),
                     values=c("1" = 21, "2"=22, "3"=23, "4"=24))+
  scale_fill_manual(name="cat groups", breaks= c("1", "2", "3", "4"), 
                    labels= c("domesticated", "hybrid", "lybica", "silvestris"), 
                    values=c("1"="#d28032", "2"="#445632", "3"="#95563f", "4"="#475f73"))+
  theme_bw()+
  theme(legend.position = "none")+
  geom_polygon(data = hulls, alpha = 0.2,aes(x=ForamenMagnumMMs, 
                                             y=BrainVolumeCM3, 
                                             fill=HybridIDLBQbasedSimplified))



# V T L  L E N G T H 

fullpalate2<- lm(PalateLengthMM~ ForamenMagnumMM+ Hybrid_ID_LBQbased_simplified, 
                 data=Katze)

nullpalate2<- lm(PalateLengthMM~ ForamenMagnumMM, 
                 data=Katze)

residualPlot(fullpalate2)
hist(resid(fullpalate2))

#check residuals
diagnostics.plot(fullpalate2)

#model summary
print(summary(fullpalate2), corr=F)

#collinearity; model without random effects to check
detach("package:usdm", unload=TRUE)
library(car)

vif(fullpalate2)

ID<- as.factor(Katze$Hybrid_ID_LBQbased_simplified)

plot(y=Katze$WeightG, x=ID)
#full null test
anova(nullpalate2, fullpalate2, test="Chisq")
overdisp.test(fullpalate2)
#plot(fullpalate)
#plot model
plot(residuals(fullpalate2))
summary(fullpalate2)
hist(residuals(fullpalate2), breaks=30)
summary(fullpalate2)
summary(nullpalate2)
Anova(fullpalate2)

r2(fullpalate2)
check_singularity(fullpalate2)
check_model(fullpalate2)
plot(compare_performance(fullpalate2, nullpalate2))
#model stability; results 0.64
max(abs(dffits(fullpalate2)))
    

find_hull <- function(cat_poly2) cat_poly2[chull(cat_poly2$ForamenMagnumMMs, cat_poly2$PalateLengthMMs), ]
hulls <- ddply(cat_poly2, "HybridIDLBQbasedSimplified", find_hull)

VTL_FM<- ggplot(data=cat_poly2)+
  geom_point(size=3,alpha=0.8, aes(y=PalateLengthMMs, x=ForamenMagnumMMs, 
                                   shape=HybridIDLBQbasedSimplified, fill=HybridIDLBQbasedSimplified))+
  ylab("palate length (mm)") +
  xlab("basal skull length (mm)")+
  scale_shape_manual(name="cat groups", breaks= c("1", "2", "3", "4"), 
                     labels= c("domesticated", "hybrid", "lybica", "silvestris"),
                     values=c("1" = 21, "2"=22, "3"=23, "4"=24))+
  scale_fill_manual(name="cat groups", breaks= c("1", "2", "3", "4"), 
                    labels= c("domesticated", "hybrid", "lybica", "silvestris"), 
                    values=c("1"="#d28032", "2"="#445632", "3"="#95563f", "4"="#475f73"))+
  theme_bw()+
  theme(legend.position = "none")+
  geom_polygon(data = hulls, alpha = 0.2,aes(x=ForamenMagnumMMs, 
                                             y=PalateLengthMMs, 
                                             fill=HybridIDLBQbasedSimplified))


# graph combo create legend dataset
catlegends <-read_excel("/Users/raffaela/Documents/PhD_Projects/20170222_Cats/Paper/Images_digitize/Legend_data_NOACTUALDATA_JUSTLEGEND.xlsx")

catlegends$species_simplified<-as.factor(catlegends$species_simplified)

legendplot<- ggplot(data=catlegends)+
  geom_point(size=3, alpha=0.8, aes(y='kubikWurzelKapazitaet(cm)', x='basallaenge(mm)',
                                    fill=species_simplified, shape=species_simplified))+
  ylab("palate length (mm)")+
  xlab("basal skull length (mm)")+
  scale_shape_manual(name="cat groups", breaks= c("domesticated", "hybrid", "lybica", "silvestris", "feral (hybrid?)"),
                     labels= c("domesticated", "hybrid", "lybica", "silvestris","feral (hybrid?)"),
                     values=c("domesticated" = 21, "hybrid"=22, "lybica"=23, "silvestris"=24, "feral (hybrid?)"=25))+
  scale_fill_manual(name="cat groups", breaks= c("domesticated", "hybrid", "lybica", "silvestris","feral (hybrid?)"),
                    labels= c("domesticated", "hybrid", "lybica", "silvestris","feral (hybrid?)"),
                    values=c("domesticated"="#d28032", "hybrid"="#445632", "lybica"="#95563f", "silvestris"="#475f73","feral (hybrid?)"="#445632"))+
  theme_bw()

legend_1 <- get_legend(legendplot)

title_1 <- ggdraw() + draw_label("cranial capacity and skull  measurements", fontface='bold')

p_2 <- plot_grid(VTL_FM, CM3_FM, 
                 labels = c("A", "B"),
                 nrow = 1, ncol= 2)


completeGraph <- plot_grid(title_1, p_2, ncol=1, rel_heights=c(0.1,1)) # rel_heights values control title margins

leg._plot<- plot_grid(legend_1,
                      nrow=1, ncol=1)

skullplot<- plot_grid(completeGraph, leg._plot,
                      ncol=2, rel_widths = c(1,0.2))

#ggsave(skullplot, filename = "/Users/raffaela/Documents/PhD_Projects/20170222_Cats/Paper/SkullGraph_20210331.tiff", dpi = 300, width = 10, height = 4)

############################################################################################
#                                                                                          #
#     D A T A  S E T  W I T H   E X C L U D E D  S K U L L S - NOT PRESENTED IN RESULTS    #
#                                                                                          #
############################################################################################

# C M 3  V O L U M E

fullcm3R<- lm(brainVolumeCM3 ~ ForamenMagnumMM+ Hybrid_ID_LBQbased_simplified, 
             data=cat)

nullcm3R<- lm(brainVolumeCM3~ ForamenMagnumMM, 
             data=cat)

residualPlot(fullcm3R)
hist(resid(fullcm3R))

#check residuals
diagnostics.plot(fullcm3R)

#model summary
print(summary(fullcm3R), corr=F)

#collinearity; model without random effects to check
detach("package:usdm", unload=TRUE)
library(car)

vif(fullcm3R)

ID<- as.factor(cat$Hybrid_ID_LBQbased_simplified)

plot(y=cat$WeightG, x=ID)
#full null test
anova(nullcm3R, fullcm3R, test="Chisq")
overdisp.test(fullcm3R)
#plot(fullcm3R)
#plot model
plot(residuals(fullcm3R))
summary(fullcm3R)
hist(residuals(fullcm3R), breaks=30)
summary(fullcm3R)
Anova(fullcm3R)

r2(fullcm3R)
check_singularity(fullcm3R)
check_model(fullcm3R)
plot(compare_performance(fullcm3R, nullcm3R))
#check model stability; values under 2 are fine -> my model has 1.47
max(abs(dffits(fullcm3R)))


# V T L  L E N G T H 

fullpalate2R<- lm(PalateLengthMM~ ForamenMagnumMM+ Hybrid_ID_LBQbased_simplified, 
                 data=cat)

nullpalate2R<- lm(PalateLengthMM~ ForamenMagnumMM, 
                 data=cat)

residualPlot(fullpalate2R)
hist(resid(fullpalate2R))

#check residuals
diagnostics.plot(fullpalate2R)

#model summary
print(summary(fullpalate2R), corr=F)

#collinearity; model without random effects to check
detach("package:usdm", unload=TRUE)
library(car)

vif(fullpalate2R)

ID<- as.factor(cat$Hybrid_ID_LBQbased_simplified)

plot(y=cat$WeightG, x=ID)
#full null test
anova(nullpalate2R, fullpalate2R, test="Chisq")
overdisp.test(fullpalate2R)
#plot(fullpalate)
#plot model
plot(residuals(fullpalate2R))
summary(fullpalate2R)
hist(residuals(fullpalate2R), breaks=30)
summary(fullpalate2R)
Anova(fullpalate2R)

r2(fullpalate2R)
check_singularity(fullpalate2R)
check_model(fullpalate2R)
plot(compare_performance(fullpalate2R, nullpalate2R))
#model stability; results 2.39
max(abs(dffits(fullpalate2R)))

###########################################
#                                         #
#     D I G I T I Z E  O L D  D A T A     #
#                                         #
###########################################

#####################################################################
#                                                                   #
#     Hirngrössenvariation im Felis silvestirs-Kreis -> Hemmer      #
#                                                                   #
#####################################################################



#set working directory
# setwd("/Users/raffaela/Documents/PhD_Projects/20170222_Cats/Paper/Images_digitize")

#read, calibrate and get image data
# cal = ReadAndCal("Hemmer_cm3_basal.png")

# This opens the jpeg in a plotting window and lets you 
# define points on the x and y axis. You must start by 
# clicking on the left-most x-axis point, then the right-most
# axis point, followed by the lower y-axis point and finally 
# the upper y-axis point. You don’t need to choose the end points
# of the axis, only two points on the axis that you know the 
# x or y value for. As you click on each of the 4 points, 
# the coordinates are saved in the object cal.

#next step is data point collection
# data.points = DigitData(col = 'red')

# Finally, you need to convert those raw x,y coordinates into the same 
# scale as the original graph. You do this by calling the Calibrate 
# function and feeding it your data.point list, the cal list that 
# contains your 4 control points from the first step, and then 4 numeric 
# values that represent the 4 original points you clicked on the x and 
# y axes. These values should be in the original scale of the figure 
# (i.e. read the values off the graph’s tick marks).

# #calibrate
# df_domesticated = Calibrate(data.points, cal, 70, 90, 2.8, 3.6)
# head(df_domesticated)
# 
# #repeat process for all other groups!
# cal = ReadAndCal("Hemmer_cm3_basal.png")
# data.points = DigitData(col = 'red')
# df_silvestris = Calibrate(data.points, cal, 70, 90, 2.8, 3.6)
# head(df_silvestris)
# 
# cal = ReadAndCal("Hemmer_cm3_basal.png")
# data.points = DigitData(col = 'red')
# df_hybrid = Calibrate(data.points, cal, 70, 90, 2.8, 3.6)
# head(df_hybrid)
# 
# cal = ReadAndCal("Hemmer_cm3_basal.png")
# data.points = DigitData(col = 'red')
# df_FelSilOrn = Calibrate(data.points, cal, 70, 90, 2.8, 3.6)
# head(df_FelSilOrn)
# 
# cal = ReadAndCal("Hemmer_cm3_basal.png")
# data.points = DigitData(col = 'red')
# df_FelSil = Calibrate(data.points, cal, 70, 90, 2.8, 3.6)
# head(df_FelSil)
# 
# cal = ReadAndCal("Hemmer_cm3_basal.png")
# data.points = DigitData(col = 'red')
# df_FelSilLyb = Calibrate(data.points, cal, 70, 90, 2.8, 3.6)
# head(df_FelSilLyb)
# 
# #add species name to extracted data frame and change columnn names
# df_domesticated["species"]="domesticated"
# names(df_domesticated)[names(df_domesticated) == "x"] <- "basallaenge(mm)"
# names(df_domesticated)[names(df_domesticated) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# df_silvestris["species"]="silvestris"
# names(df_silvestris)[names(df_silvestris) == "x"] <- "basallaenge(mm)"
# names(df_silvestris)[names(df_silvestris) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# df_hybrid["species"]="hybrid"
# names(df_hybrid)[names(df_hybrid) == "x"] <- "basallaenge(mm)"
# names(df_hybrid)[names(df_hybrid) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# df_FelSilOrn["species"]="FelSilOrn"
# names(df_FelSilOrn)[names(df_FelSilOrn) == "x"] <- "basallaenge(mm)"
# names(df_FelSilOrn)[names(df_FelSilOrn) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# df_FelSil["species"]="FelSilFalbkatze"
# names(df_FelSil)[names(df_FelSil) == "x"] <- "basallaenge(mm)"
# names(df_FelSil)[names(df_FelSil) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# df_FelSilLyb["species"]="FelSilLyb"
# names(df_FelSilLyb)[names(df_FelSilLyb) == "x"] <- "basallaenge(mm)"
# names(df_FelSilLyb)[names(df_FelSilLyb) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# #combine all data frames
# digitized_Hemmer_complete <- rbind(df_domesticated, df_silvestris,df_hybrid,df_FelSilOrn,df_FelSil,df_FelSilLyb)
# 
# write_xlsx(digitized_Hemmer_complete,"Hemmer_Basallaennge_CM3.xlsx")


############## second figure
# 
# cal = ReadAndCal("cm3_schaedellaenge_hemmer.png")
# data.points = DigitData(col = 'red')
# df_Hauskatze = Calibrate(data.points, cal, 80, 110, 2.8, 3.6)
# head(df_Hauskatze)
# 
# cal = ReadAndCal("cm3_schaedellaenge_hemmer.png")
# data.points = DigitData(col = 'red')
# df_EuroWild = Calibrate(data.points, cal, 80, 110, 2.8, 3.6)
# head(df_EuroWild)
# 
# cal = ReadAndCal("cm3_schaedellaenge_hemmer.png")
# data.points = DigitData(col = 'red')
# df_FELSilOrn = Calibrate(data.points, cal, 80, 110, 2.8, 3.6)
# head(df_FELSilOrn)
# 
# cal = ReadAndCal("cm3_schaedellaenge_hemmer.png")
# data.points = DigitData(col = 'red')
# df_FELSil = Calibrate(data.points, cal, 80, 110, 2.8, 3.6)
# head(df_FELSil)
# 
# cal = ReadAndCal("cm3_schaedellaenge_hemmer.png")
# data.points = DigitData(col = 'red')
# df_FELSilLyb = Calibrate(data.points, cal, 80, 110, 2.8, 3.6)
# head(df_FELSilLyb)
# 
# #add species name to extracted data frame and change columnn names
# df_Hauskatze["species"]="domesticated"
# names(df_Hauskatze)[names(df_Hauskatze) == "x"] <- "schaedellaenge(mm)"
# names(df_Hauskatze)[names(df_Hauskatze) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# df_EuroWild["species"]="silvestris"
# names(df_EuroWild)[names(df_EuroWild) == "x"] <- "schaedellaenge(mm)"
# names(df_EuroWild)[names(df_EuroWild) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# df_FELSilOrn["species"]="FelSilOrn"
# names(df_FELSilOrn)[names(df_FELSilOrn) == "x"] <- "schaedellaenge(mm)"
# names(df_FELSilOrn)[names(df_FELSilOrn) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# df_FELSil["species"]="FelSilFalbkatze"
# names(df_FELSil)[names(df_FELSil) == "x"] <- "schaedellaenge(mm)"
# names(df_FELSil)[names(df_FELSil) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# df_FELSilLyb["species"]="FelSilLyb"
# names(df_FELSilLyb)[names(df_FELSilLyb) == "x"] <- "schaedellaenge(mm)"
# names(df_FELSilLyb)[names(df_FELSilLyb) == "y"] <- "kubikWurzelKapazitaet(cm)"
# 
# #combine all data frames
# digitized_Hemmer2_complete <- rbind(df_Hauskatze, df_EuroWild,df_FELSilOrn,df_FELSil,df_FELSilLyb)
# 
# write_xlsx(digitized_Hemmer2_complete,"Hemmer_Schaedellaennge_CM3.xlsx")


###########################################################################################################################
#                                                                                                                         #
#     L'identification du Chat forestier d'Europe Felis s. silvestris par une méthode ostéomatrique -> Schauenberger      #
#                                                                                                                         #
###########################################################################################################################

# #read, calibrate and get image data
# cal = ReadAndCal("Schauenberger_Fig2_1969_rotated.png")
# data.points = DigitData(col = 'red')
# df_dom = Calibrate(data.points, cal, 75, 110, 20, 50)
# head(df_dom)
# 
# cal = ReadAndCal("Schauenberger_Fig2_1969_rotated.png")
# data.points = DigitData(col = 'red')
# df_sil = Calibrate(data.points, cal, 75, 110, 20, 50)
# head(df_sil)
# 
# 
# df_dom["species"]="domesticated"
# names(df_dom)[names(df_dom) == "x"] <- "Schaedellaenge(mm)"
# names(df_dom)[names(df_dom) == "y"] <- "Kranialvolumen(cm)"
# 
# df_sil["species"]="silvestris"
# names(df_sil)[names(df_sil) == "x"] <- "Schaedellaenge(mm)"
# names(df_sil)[names(df_sil) == "y"] <- "Kranialvolumen(cm)"
# 
# #combine all data frames
# digitized_Schauenberger_complete <- rbind(df_dom, df_sil)
# 
# write_xlsx(digitized_Schauenberger_complete,"Schauenberger_Schaedellaennge_CM3.xlsx")

#set working directory to digitized data
setwd("/Users/raffaela/Documents/PhD_Projects/20170222_Cats/Paper/Images_digitize")


#open all previously created and saved excel files
S_Schädel_CM3 <- read_excel("Schauenberger_Schaedellaennge_CM3.xlsx")
H_Schädel_CM3 <- read_excel("Hemmer_Schaedellaennge_CM3.xlsx")
H_Basal_CM3 <- read_excel("Hemmer_Basallaennge_CM3.xlsx")

#create polygon clouds
S_Schädel_CM3 <- na.omit(S_Schädel_CM3)
find_hull <- function(S_Schädel_CM3) S_Schädel_CM3[chull(S_Schädel_CM3$`Schaedellaenge(mm)`, S_Schädel_CM3$`Kranialvolumen(cm)`), ]
hulls <- ddply(S_Schädel_CM3, "species", find_hull)


S_Schädel_CM3_graph<- ggplot(data=S_Schädel_CM3)+
  geom_point(size=3, alpha=0.8, aes(x=`Schaedellaenge(mm)`, 
                 y=`Kranialvolumen(cm)`, 
                 shape=species,
                 fill=species))+
  ylab(expression(paste(
    "cranial volume (",
    cm^3,
    ")", sep=""))) +
  xlab("skull length (mm)")+
  scale_shape_manual(name="cat groups", breaks= c("domesticated",  "silvestris"), 
                     labels= c("domesticated", "silvestris"),
                     values=c("domesticated" = 21, "silvestris"=24))+
  scale_fill_manual(name="cat groups", breaks= c("domesticated", "silvestris"), 
                    labels= c("domesticated", "silvestris"), 
                    values=c("domesticated"="#d28032",  "silvestris"="#475f73"))+
 # geom_abline(size=1,colour="#d28032", lty="dashed",aes(intercept=5.72553, slope=0.38373))+ #domesticated
 # geom_abline(size=1,colour="#475f73", lty="dashed",aes(intercept=6.5176, slope=0.38373))+ #silvestris
  theme_bw()+
  theme(legend.position = "none")+
  annotate("text", x = 85, y = 48, label = "from Schauenberger 1969; Fig.2", fontface="bold", size=3)+
  geom_polygon(data = hulls, alpha = 0.2,aes(x=`Schaedellaenge(mm)`, 
                                             y=`Kranialvolumen(cm)`, 
                                             fill=species))


#create polygon clouds
H_Basal_CM3 <- na.omit(H_Basal_CM3)
find_hull <- function(H_Basal_CM3) H_Basal_CM3[chull(H_Basal_CM3$`basallaenge(mm)`, H_Basal_CM3$`kubikWurzelKapazitaet(cm)`), ]
hulls <- ddply(H_Basal_CM3, "species_simplified", find_hull)


H_Basal_CM3_graph<- ggplot(data=H_Basal_CM3)+
  geom_point(size=3, alpha=0.8,aes(x=H_Basal_CM3$`basallaenge(mm)`, 
                 y=H_Basal_CM3$`kubikWurzelKapazitaet(cm)`, 
                 shape=H_Basal_CM3$species_simplified,
                 fill=H_Basal_CM3$species_simplified))+
  ylab(expression(paste(
    "cranial volume (",
    sqrt(capacity,3),
    " (cm)", sep=""))) +
  xlab("basal skull length (mm)")+
  scale_shape_manual(name="cat groups", breaks= c("domesticated", "feral (hybrid?)", "lybica", "silvestris"), 
                     labels= c("domesticated", "feral (hybrid?)", "lybica", "silvestris"),
                     values=c("domesticated" = 21, "feral (hybrid?)"=25, "lybica"=23, "silvestris"=24))+
  scale_fill_manual(name="cat groups", breaks= c("domesticated", "feral (hybrid?)", "lybica", "silvestris"), 
                    labels= c("domesticated", "feral (hybrid?)", "lybica", "silvestris"), 
                    values=c("domesticated"="#d28032", "feral (hybrid?)"="#445632", "lybica"="#95563f", "silvestris"="#475f73"))+
#  geom_abline(size=1, colour="#d28032", lty="dashed",aes(intercept=-0.56000, slope=0.54334))+ #domesticated
#  geom_abline(size=1, colour="#445632", lty="dashed",aes(intercept=8.22899, slope=0.54334))+ #hybrid
#  geom_abline(size=1, colour="#95563f", lty="dashed",aes(intercept=3.45417, slope=0.54334))+ #ornata
#  geom_abline(size=1, colour="#475f73", lty="dashed",aes(intercept=18.00646, slope=0.54334))+ #silvestris
  theme_bw()+
  theme(legend.position = "none")+
  annotate("text", x = 71, y = 3.5, label = "from Hemmer 1972; Fig.2", fontface="bold", size=3)+
  geom_polygon(data = hulls, alpha = 0.2,aes(x=`basallaenge(mm)`, 
                                             y=`kubikWurzelKapazitaet(cm)`, 
                                             fill=species_simplified))


#create polygon clouds
H_Schädel_CM3 <- na.omit(H_Schädel_CM3)
find_hull <- function(H_Schädel_CM3) H_Schädel_CM3[chull(H_Schädel_CM3$`schaedellaenge(mm)`, H_Schädel_CM3$`kubikWurzelKapazitaet(cm)`), ]
hulls <- ddply(H_Schädel_CM3, "species_simplified", find_hull)


H_Schädel_CM3_graph<- ggplot(data=H_Schädel_CM3)+
  geom_point(size=3, alpha=0.8,aes(x=H_Schädel_CM3$`schaedellaenge(mm)`, 
                 y=H_Schädel_CM3$`kubikWurzelKapazitaet(cm)`, 
                 shape=H_Schädel_CM3$species_simplified,
                 fill=H_Schädel_CM3$species_simplified))+
  ylab(expression(paste(
    "cranial volume (",
    sqrt(capacity,3),
    " (cm)", sep=""))) +
  xlab("skull length (mm)")+
  scale_shape_manual(name="cat groups", breaks= c("domesticated", "hybrid", "lybica", "silvestris"), 
                     labels= c("domesticated", "hybrid", "lybica", "silvestris"),
                     values=c("domesticated" = 21, "hybrid"=22, "lybica"=23, "silvestris"=24))+
  scale_fill_manual(name="cat groups", breaks= c("domesticated", "hybrid", "lybica", "silvestris"), 
                    labels= c("omesticated", "hybrid", "lybica", "silvestris"), 
                    values=c("domesticated"="#d28032", "hybrid"="#445632", "lybica"="#95563f", "silvestris"="#475f73"))+
  #  geom_abline(size=1, colour="#d28032", lty="dashed",aes(intercept=-0.56000, slope=0.54334))+ #domesticated
  #  geom_abline(size=1, colour="#445632", lty="dashed",aes(intercept=8.22899, slope=0.54334))+ #hybrid
  #  geom_abline(size=1, colour="#95563f", lty="dashed",aes(intercept=3.45417, slope=0.54334))+ #ornata
  #  geom_abline(size=1, colour="#475f73", lty="dashed",aes(intercept=18.00646, slope=0.54334))+ #silvestris
  theme_bw()+
  theme(legend.position = "none")+
  annotate("text", x = 85, y = 3.5, label = "from Hemmer 1972; Fig.1", fontface="bold", size=3)+
  geom_polygon(data = hulls, alpha = 0.2,aes(x=`schaedellaenge(mm)`, 
                                             y=`kubikWurzelKapazitaet(cm)`, 
                                             fill=species_simplified))






# graph combo
legend_1 <- get_legend(legendplot)

title_1 <- ggdraw() + draw_label("cranial capacity and skull  measurements", fontface='bold')

p_2 <- plot_grid(VTL_FM, CM3_FM,
                 labels = c("A", "B"),
                 nrow = 1, ncol= 2)

p_3 <- plot_grid(S_Schädel_CM3_graph, H_Basal_CM3_graph,H_Schädel_CM3_graph,
                 labels = c("C", "D", "E"),
                 nrow = 1, ncol= 3)


completeGraph <- plot_grid( p_2, p_3, ncol=1, rel_heights=c(1,1)) # rel_heights values control title margins

leg._plot<- plot_grid(legend_1,
                      nrow=1, ncol=1)

skullplot<- plot_grid(completeGraph, leg._plot,
                      ncol=2, rel_widths = c(1,0.2))

#ggsave(skullplot, filename = "/Users/raffaela/Documents/PhD_Projects/20170222_Cats/Paper/SkullGraph_20210806_noRegressionlines.tiff", dpi = 300, width = 15, height = 8)


