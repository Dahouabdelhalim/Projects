
####################Setting directory###############################
getwd()
################### Data import ######################################
##### File name "analyses.csv" #######################################
bamboo = read.delim("analyses.csv",dec=".", sep=",",header=TRUE)
View(bamboo)
str(bamboo)

############### data tidying
bamboo = bamboo[, -25]  #### delete empty column 25,

############ Convert plot number and elevation into factor ###############
bamboo$ï..Plot_No = as.factor(bamboo$ï..Plot_No) ## as factor plot no
bamboo$Elevation=as.factor(bamboo$Elevation) ### as factor elevation
str(bamboo)

########## Convert the fresh weight integers to numeric #################
bamboo$S_0grW = as.numeric(bamboo$S_0grW)
bamboo$S_25grW = as.numeric(bamboo$S_25grW)
bamboo$S_50grW = as.numeric(bamboo$S_50grW)
bamboo$S_75grW = as.numeric(bamboo$S_75grW)
bamboo$Br_grW = as.numeric(bamboo$Br_grW)
bamboo$Fol_grW = as.numeric(bamboo$Fol_grW)
str(bamboo)


################ Computing  sub samples biomass from fresh to oven dry weight ratio #####
#1. Biomass of stem sub samples
bamboo$AGB_0 = (bamboo$Wt_stem*bamboo$S_0dryW)/(bamboo$S_0grW*4) ##Sub sample at 0 (basal) potion of the stem
bamboo$AGB_25 = (bamboo$Wt_stem*bamboo$S_25dryW)/(bamboo$S_25grW*4)# sub sample at 0.25 portion of the stem
bamboo$AGB_50 = (bamboo$Wt_stem*bamboo$S_50dryW)/(bamboo$S_50grW*4)# sub sample at 0.5 portion of the stem
bamboo$AGB_75 = (bamboo$Wt_stem*bamboo$S_75dryW)/(bamboo$S_75grW*4)# sub sample at 0.75 portion of the stem
bamboo$Stem_AGB = bamboo$AGB_0+bamboo$AGB_25+bamboo$AGB_50+bamboo$AGB_75 #TOTAL SUBSAMPLE STEM  biomass
View(bamboo)

#2. Biomass of bamboo BRANCHES sub sample #########
bamboo$branches_AGB = (bamboo$Wt_branch*bamboo$Br_dryW)/(bamboo$Br_grW)

#3. Biomass of bamboo FOLIAGE sub sample ##########
bamboo$Foliage_AGB = (bamboo$Wt_Foliage*bamboo$Fol_dryW)/(bamboo$Fol_grW)

#############  TOTAL SAMPLE BIOMASS ##############
bamboo$Total_AGB = bamboo$Stem_AGB+bamboo$branches_AGB+bamboo$Foliage_AGB


View(bamboo)


########################## COMPUTING ABOVEGROUND BIOMASS IN TONNES PER HACTARE ##############

1. ########### Calculating plot area in per hectare basis 
names(bamboo)[names(bamboo) == "Total_AGB"] <- "AGB" ###Renaming total AGB
names(bamboo)[names(bamboo) == "ï..Plot_No"] <- "Plot" ### Renaming plot


############ Creating a subset for the columns to be used in computation only #############
bamboo_select  =  subset(bamboo, select = c( "Plot", "Elevation", "Eco_type", "Mgt_option", "Species_name", "Dbh","Height","Age", "Stem_AGB","branches_AGB","Foliage_AGB","AGB"))
View(bamboo_select)

############## REGRESSION MODEL FITTING ######################################################
bamboo.del = na.omit(bamboo_select) # omitting all rows with missing variables (NA)
View(bamboo.del)
str(bamboo.del)

######################### Species_generalized models for both o.abyssinica and B.vulgaris  ####

Dbh  = bamboo.del$Dbh  #######Defining variables 
AGB = bamboo.del$AGB
Height = bamboo.del$Height
Plot = bamboo.del$Plot

n = bamboo.del$Dbh ####### Making a data frame for exponential model
m = bamboo.del$AGB
nm = data.frame(n,m)


windowsFonts(s = windowsFont("Times New Roman")) ##### font selection

library(basicTrendline)

########### Logarithmic model ###############################################################
trendline(Dbh, AGB, model="log2P", ePos.x= "bottomright", show.equation = TRUE,show.Rpvalue = TRUE, eSize = 1.2, eDigit =4, yhat = TRUE, CI.fill = FALSE, CI.color = "black", CI.lty = 2, linecolor = "blue", lwd = 3, pch = 19,col = "red", las = 1.3, cex=1.5,  cex.axis = 1.5, cex.lab =1.5, family = "s", bty = "l", main = "3", xlab = "Culm diameter (cm)",  ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))

############### Power equation ############################################################

trendline(Dbh, AGB, model = "power2P", ePos.x = "bottomright", show.equation = TRUE,show.Rpvalue = TRUE, eSize = 1.2, eDigit =4, yhat = TRUE, CI.fill = FALSE, CI.color = "black", CI.lty = 2, linecolor = "blue", lwd = 3, pch = 19,col = "red", las = 1.5, cex.axis = 1.5, cex.lab =1.5, family = "s", bty = "l", main = "3", xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))
par(mgp=c(2.3,1,0))

############## Exponential equation #########################################################
trendline(n, m, model="exp2P",ePos.x = "bottomright", show.equation = TRUE, show.Rpvalue = TRUE,  eSize = 1.2, eDigit =4, yhat = TRUE, CI.fill = FALSE, CI.color = "black", CI.lty = 2, linecolor = "blue", lwd = 3, pch = 19,col = "red", las = 1.3, cex.axis = 1.5, cex.lab =1.5, family = "s", bty = "l", main = "3", xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))

############## Quadratic equation ###########################################################
text(x=5,y=5, family="serif")
op = par(family = "serif")

trendline(Dbh, AGB,  model="line3P",ePos.x = "topleft", show.equation = TRUE, show.Rpvalue = TRUE, eSize = 1.5, eDigit =4, yhat = TRUE, CI.fill = FALSE,  CI.color = NA, linecolor = "black", lwd = 4, pch = 19, col = "black", tck = 0.03, las = 1.3, cex = 2, cex.axis = 1.8, cex.lab =2, family = "s", bty = "l", main = "3",plot(x, y), box(col = "black", lty = "solid", lwd = 4), xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))

minor.tick (nx = 4, ny = 2, tick.ratio =1)

par(mgp=c(2.4,1,0))


############################## Nonlinear regression model fitting for B.vulgaris############## ################## #########################################################################

bambu = bamboo.del[which(bamboo.del$Species_name == "B.vulgaris"), ] # selecting B.vulgaris data

Dbh  = bambu$Dbh #Defining varaibles
AGB = bambu$AGB
Height = bambu$Height
Plot = bamboo.del$Plot

x = bambu$Dbh ####### Making a data frame for exponential model
y = bambu$AGB
xy = data.frame(x,y)

windowsFonts(s = windowsFont("Times New Roman")) ##### font selection

########### Logarithmic model
text(x=5,y=5, family="serif")
op = par(family = "serif")

trendline(Dbh, AGB,  model="log2P",ePos.x = "topleft", show.equation = TRUE, show.Rpvalue = TRUE, eSize = 1.3, eDigit =4, yhat = TRUE, CI.fill = FALSE,  CI.color = NA, linecolor = "black", lwd = 2, pch = 19, col = "black", tck = 0.03, las = 1.3, cex = 1.5, cex.axis = 1.3, cex.lab =1.3, family = "s", bty = "l", main = "3",plot(x, y), box(col = "black", lty = "solid", lwd = 2.5), xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))

minor.tick (nx = 2, ny = 2, tick.ratio =-0.5)

par(mgp=c(2.4,1,0))

############### Power equation ###############################################################
trendline(Dbh, AGB, model = "power2P", ePos.x = "bottomright", show.equation = TRUE,show.Rpvalue = TRUE, eSize = 1.2, eDigit =4, yhat = TRUE, CI.fill = FALSE, CI.color = "black", CI.lty = 2, linecolor = "blue", lwd = 3, pch = 19,col = "red", las = 1.5, cex.axis = 1.5, cex.lab =1.5, family = "s", bty = "l", main = "3", xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))
par(mgp=c(2.3,1,0))


############## Exponential equation #########################################################
trendline(x, y, model="exp2P",ePos.x = "bottomright", show.equation = TRUE, show.Rpvalue = TRUE,  eSize = 1.2, eDigit =4, yhat = TRUE, CI.fill = FALSE, CI.color = "black", CI.lty = 2, linecolor = "blue", lwd = 3, pch = 19,col = "red", las = 1.3, cex.axis = 1.5, cex.lab =1.5, family = "s", bty = "l", main = "3", xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))


############################# Quadratic equation ###########################################

text(x=5,y=5, family="serif")
op = par(family = "serif")

trendline(Dbh, AGB,  model="line3P",ePos.x = "topleft", show.equation = FALSE, show.Rpvalue = TRUE, eSize = 1.3, eDigit =4, yhat = TRUE, CI.fill = FALSE,  CI.color = NA, linecolor = "black", lwd = 2, pch = 19, col = "black", tck = 0.03, las = 1.3, cex = 1.5, cex.axis = 1.3, cex.lab =1.3, family = "s", bty = "l", main = "3",plot(x, y), box(col = "black", lty = "solid", lwd = 2.5), xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))

minor.tick (nx = 2, ny = 2, tick.ratio =-0.5)

par(mgp=c(2.4,1,0))


##################################################################################################### Nonlinear Regression model fits for oxytenanthera abyssinica ##########################

oxyt = bamboo.del[which(bamboo.del$Species_name == "O.abyssinica"), ] #### select o.abyssinica dataset only
str(oxyt)

Dbh  = oxyt$Dbh ##############Defining variables
AGB = oxyt$AGB
Height = oxyt$Height

a = oxyt$Dbh ####### Making a data frame for exponential model
b = oxyt$AGB
ab = data.frame(x,y)


windowsFonts(s = windowsFont("Times New Roman")) ##### font selection

################ Logarithmic model ###################################################

text(x=5,y=5, family="serif")
op = par(family = "serif")

trendline(Dbh, AGB,  model="log2P",ePos.x = "topleft", show.equation = TRUE, show.Rpvalue = TRUE, eSize = 1.3, eDigit =4, yhat = TRUE, CI.fill = FALSE,  CI.color = NA, linecolor = "black", lwd = 2, pch = 19, col = "black", tck = 0.03, las = 1.3, cex = 1.5, cex.axis = 1.3, cex.lab =1.3, family = "s", bty = "l", main = "2.5",plot(x, y), box(col = "black", lty = "solid", lwd = 2), xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))

############### Power equation ###########################################################
text(x=5,y=5, family="serif")
op = par(family = "serif")

trendline(Dbh, AGB, model="log2P", ePos.x= "bottomright", show.equation = TRUE,show.Rpvalue = TRUE, eSize = 1.2, eDigit =4, yhat = TRUE, CI.fill = FALSE, CI.color = "black", CI.lty = 2, linecolor = "blue", lwd = 3, pch = 19,col = "red", las = 1.3, cex=1.5,  cex.axis = 1.5, cex.lab =1.5, family = "s", bty = "l", main = "3", xlab = "Culm diameter (cm)",  ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))

############### Power equation
trendline(Dbh, AGB, model = "power2P", ePos.x = "bottomright", show.equation = TRUE,show.Rpvalue = TRUE, eSize = 1.2, eDigit =4, yhat = TRUE, CI.fill = FALSE, CI.color = "black", CI.lty = 2, linecolor = "blue", lwd = 3, pch = 19,col = "red", las = 1.5, cex.axis = 1.5, cex.lab =1.5, family = "s", bty = "l", main = "3", xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))
par(mgp=c(2.3,1,0))

############## Exponential equation
trendline(a, b, model="exp2P",ePos.x = "bottomright", show.equation = TRUE, show.Rpvalue = TRUE,  eSize = 1.2, eDigit =4, yhat = TRUE, CI.fill = FALSE, CI.color = "black", CI.lty = 2, linecolor = "blue", lwd = 3, pch = 19,col = "red", las = 1.3, cex.axis = 1.5, cex.lab =1.5, family = "s", bty = "l", main = "3", xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))


############## Quadratic equation

text(x=5,y=5, family="serif")
op = par(family = "serif")

trendline(Dbh, AGB,  model="line3P",ePos.x = "topleft", show.equation = TRUE, show.Rpvalue = TRUE, eSize = 1.3, eDigit =4, yhat = TRUE, CI.fill = FALSE,  CI.color = NA, linecolor = "black", lwd = 2, pch = 19, col = "black", tck = 0.03, las = 1.3, cex = 1.5, cex.axis = 1.3, cex.lab =1.3, family = "s", bty = "l", main = "2.5",plot(x, y), box(col = "black", lty = "solid", lwd = 2), xlab = "Culm diameter (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))

minor.tick (nx = 2, ny = 2, tick.ratio =-0.5)

par(mgp=c(2.4,1,0))



############# Computing the plot area per hectare ###################################
Area = (3.14*5.65^2)/10000   #pie/radius(5.65^2m)/metres-ha conversion
bamboo$Area = rep(Area, length(bamboo$AGB))
bamboo$Area =  Area


############# Substituting Na's by the equations for calculating AGB tonnes per hectare (with each species selected model)

bamboo_select$AGB[is.na(bamboo_select$AGB)& bamboo_select$Species_name =="B.vulgaris"] = 
  0.2329*(bamboo_select$Dbh)^2 - 0.07558*bamboo_select$Dbh - 0.438

bamboo_select$AGB[is.na(bamboo_select$AGB)& bamboo_select$Species_name =="O.abyssinica"] = 
  -0.1091*(bamboo_select$Dbh)^2 + 4.481*bamboo_select$Dbh - 10.3

View(bamboo_select)

################### Calculating AGB in tonnes per hectare #######################################
bamboo_select$AGB_tha = (bamboo_select$AGB/bamboo$Area)/1000 ###BIOMASS/AREA (hectares)/KG-TONNES CONVERSION
AGB_tha = bamboo_select$AGB_tha 


###################### Convert elevation values into numeric for assigning gradients #############
bamboo_select$Elevation <-as.numeric(as.character(bamboo_select$Elevation, header = TRUE))
str(bamboo_select)

######################## Add a column "Gradient1" for allocating elevation gradients #############
#EACH TYPE OF SPECIES HAD  DIFFERENT MEDIUM AND HIGH ELEVATION LIMITS ##########################


Elevation = bamboo_select$Elevation
Gradient1 = bamboo_select$Elevation_gradient
Mgt_option = bamboo_select$Mgt_option


Gradient1 > 0
Gradient1[which(bamboo_select$Species_name == "B.vulgaris"& Elevation <= 900)] = "Low"
Gradient1[which(bamboo_select$Species_name == "O.abyssinica"& Elevation <= 900)] = "Low"
Gradient1[which(bamboo_select$Species_name == "B.vulgaris"& Elevation > 900 & Elevation <= 1400)] = "Medium"
Gradient1[which(bamboo_select$Species_name == "O.abyssinica"& Mgt_option == "Intensive" & Elevation > 900 & Elevation <= 1750)] = "Medium"
Gradient1[which(bamboo_select$Species_name == "O.abyssinica"& Mgt_option == "Extensive" & Elevation > 900 & Elevation <= 1900)] = "Medium"
Gradient1[which(bamboo_select$Species_name == "O.abyssinica"& Mgt_option == "Extensive" & Elevation > 1900)] = "High"
Gradient1[which(bamboo_select$Species_name == "O.abyssinica"& Mgt_option == "Intensive" & Elevation > 1750)] = "High"
Gradient1[which(bamboo_select$Species_name == "B.vulgaris" & Elevation > 1400)] = "High"
bamboo_select$Gradient1 = Gradient1


############################## CALCULATING TOTAL ABOVE GROUND CARBON PER HECTARE ################
###############################################################################################
bamboo_select
bamboo_select$AGC_tha = bamboo_select$AGB_tha * 0.5  ###AGB* an AGB-CARBON CONVERSION FACTOR
AGC_tha = bamboo_select$AGC_tha

############################## CALCULATING CULM DENSITY PER HECTARE ##############################
bamboo_select
bamboo_select$Area = (3.14*5.65*5.65)/10000
bamboo_select$Culm_density = 1/bamboo_select$Area
Culm_density = bamboo_select$Culm_density  

####################################################################################################################### COMPUTING CARBON OF 60 DESTRUCTIVELY SAMPLED CULMS ########################

bamboo.del$carbon = 0.5 *bamboo.del$AGB
Carbon  = bamboo.del$carbon

############################## ASSIGNING AGE CLASSES FOR THE 60 DESTRUCTIVELY SAMPLED CULMS ######
Age = bamboo.del$Age
bamboo.del$Age_class > 0
bamboo.del$Age_class[Age <= 2 ] = "I"
bamboo.del$Age_class[Age > 2 & Age <= 4] = "II"
bamboo.del$Age_class[Age > 4] = "III"
View(bamboo.del)

###############################################################################################
#########ASSIGNING ELEVATION GRADIENTS FOR THE 60 SAMPLED CULMS #############################

bamboo.del$Elevation <-as.numeric(as.character(bamboo.del$Elevation, header = TRUE))
str(bamboo.del)

Elevation = bamboo.del$Elevation
Gradient = bamboo.del$Elevation_gradient
Mgt_option = bamboo.del$Mgt_option


Gradient > 0
Gradient[which(bamboo.del$Species_name == "B.vulgaris"& Elevation <= 900)] = "Low"
Gradient[which(bamboo.del$Species_name == "O.abyssinica"& Elevation <= 900)] = "Low"
Gradient[which(bamboo.del$Species_name == "B.vulgaris"& Elevation > 900 & Elevation <= 1400)] = "Medium"
Gradient[which(bamboo.del$Species_name == "O.abyssinica"& Mgt_option == "Intensive" & Elevation > 900 & Elevation <= 1750)] = "Medium"
Gradient[which(bamboo.del$Species_name == "O.abyssinica"& Mgt_option == "Extensive" & Elevation > 900 & Elevation <= 1900)] = "Medium"
Gradient[which(bamboo.del$Species_name == "O.abyssinica"& Mgt_option== "Extensive" & Elevation > 1900)] = "High"
Gradient[which(bamboo.del$Species_name == "O.abyssinica"& Mgt_option== "Intensive" & Elevation > 1750)] = "High"
Gradient[which(bamboo.del$Species_name == "B.vulgaris" & Elevation > 1400)] = "High"
View(Gradient)
bamboo.del$Gradient = Gradient

View(bamboo.del)
##################### COMPUTING CARBON SEQUESTRATION RATE FOR THE SAMPLED CULMS ##################
bamboo.del$C_sequestration = bamboo.del$carbon/Age
C_sequestration = bamboo.del$C_sequestration

############################### Aggregate AGE CLASSES WITH THE SAMPLED COMPONENTS CARBON STOCK  ########################################

##### DEFINING COMPONENTS AGB
Stem_AGB = bamboo.del$Stem_AGB
branches_AGB = bamboo.del$branches_AGB
Foliage_AGB= bamboo.del$Foliage_AGB
aggregate(cbind(Stem_AGB, branches_AGB, Foliage_AGB) ~ Age_class, bamboo.del, sum)

#COMPUTING COMPONENTS CARBON SEQUESTRATION
bamboo.del$Stem.seq = (bamboo.del$Stem_AGB*0.5)/bamboo.del$Age  
bamboo.del$Branch.seq = (bamboo.del$branches_AGB*0.5)/bamboo.del$Age   
bamboo.del$Foliage.seq = (bamboo.del$Foliage_AGB*0.5)/bamboo.del$Age   
View(bamboo.del)

### Define components sequestration

Stem.seq = bamboo.del$Stem_AGB
Branch.seq = bamboo.del$branches_AGB
Foliage.seq = bamboo.del$Foliage_AGB
aggregate(cbind(Stem.seq, Branch.seq, Foliage.seq) ~ Age_class, bamboo.del, sum)
aggregate(cbind(C_sequestration) ~ Age_class, bamboo.del, sum)
mean(bamboo.del$C_sequestration)


###################### Percentage contribution of components AGB to the total AGB ############### #################################################################################################
bamboo.del$Perc_branch = (bamboo.del$branches_AGB / bamboo.del$AGB)*100
bamboo.del$perc_stem = (bamboo.del$Stem_AGB / bamboo.del$AGB)*100 
bamboo.del$perc_foliage = (bamboo.del$Foliage_AGB / bamboo.del$AGB)*100 
Perc_branch = bamboo.del$Perc_branch
perc_foliage = bamboo.del$perc_foliage

################# Aggregation whole data by plots #######################################
dat = aggregate(cbind(AGC_tha) ~ Plot + Mgt_option + Species_name + Gradient1, bamboo_select,sum)
dat = arrange(dat, Plot)

#### Displaying carbon stock or sequestration distribution in each the studied variables
qhpvt(dat, c("Gradient1","Mgt_option"), "Species_name", "mean(AGC_tha)", addTotal  = TRUE)
qhpvt(bamboo.del, c("Gradient","Mgt_option"), "Species_name", "mean(C_sequestration)", addTotal  = TRUE)
qhpvt(bamboo.del, c("Gradient","Mgt_option"), "Species_name", "mean(carbon)", addTotal  = TRUE)



###################################################################################################################################### ANOVA ################################################
####### Running factorial 3-way ANOVA
#data = "dat"
dat
anova2 <- lm(AGC_tha ~ Gradient1*Species_name*Mgt_option, data = dat) # TESTING VARIATION BETWEEN AND WITHI FACTORS AND HOW THEY INTERACT
anova(anova2)

ins.aov <- aov(AGC_tha ~ Gradient1 * Species_name * Mgt_option, dat = dat)
summary(ins.aov)

##### TUKEYS's HSD MEAN COMPARISON ##########################

TukeyHSD(ins.aov, "Gradient1") ###################dIFFERENCES ACROSS ELEVATION GRADIENT
HSD.test(anova2, "Gradient1", group = TRUE, console = TRUE)

TukeyHSD(ins.aov, "Species_name") #########VARIATION BETWEEN SPECIES
HSD.test(anova2, "Species_name", group = TRUE, console = TRUE)

TukeyHSD(ins.aov, "Mgt_option") ########### VARIATION BETWEEN MGT OPTIONS
HSD.test(anova2, "Mgt_opt", group = TRUE, console = TRUE)


######################### COMPARISON WITH OTHER BAMBOO BIOMASS MODELS #######################
bamboo.del ####data
View(bamboo.del)
###############################################################################################

bambu = bamboo.del[which(bamboo.del$Species_name == "B.vulgaris"), ] ###Subsetting for a B.vulgaris species-specific model

View(bambu)
bambu_select  =  subset(bambu, select = c( "Plot", "Species_name","Dbh", "AGB"))
View(bambu_select)
windowsFonts(s = windowsFont("Times New Roman"))

########### A scatter plot of Dbh against AGB
plot(bambu$Dbh,bambu$AGB, box(col = "black", lty = "solid", lwd = 2.5),tck = 0.03, minor.tick(nx = 2, ny = 2, tick.ratio = -0.6), pch =19, col ="black", cex= 1.5, cex.axis = 1.3, cex.lab = 1.3, family = "s", bty = "l", xlab = "Culm diameter-DBH (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))

text(x=5,y=5, family="serif")
op = par(family = "serif")

minor.tick (nx = 2, ny = 4, tick.ratio =-0.6)

par(mgp=c(2.5,1,0), las = 0)


dbh.seq <- seq(min(bambu$Dbh), max(bambu$Dbh), length.out = 100)


################## Adding lines of different species-specfic equations for B.vulgaris

agb1 = 0.2329*dbh.seq^2 - 0.07558*dbh.seq - 0.438
lines(dbh.seq, agb1, col = "black", lwd = 2.5, lty = "solid")
a =lm(dbh.seq ~ agb1)
summary(a)


agb =0.1776*dbh.seq^2 + 0.1881*dbh.seq + 1.837
lines(dbh.seq, agb, col = "black", lwd = 2.5, lty = "twodash")
g =lm(dbh.seq ~ agb)
summary(g)


agb2 = 0.09814*dbh.seq^2.36569 + 0.05216*dbh.seq^2.00483 + 0.03044*dbh.seq^1.74187
lines(dbh.seq, agb2, lty = "longdash", col = "black", lwd = 2.5)
f =lm(dbh.seq ~ agb2)
summary(f)

agb3 = (10^(2.281 + 2.149*log10(dbh.seq)))/1000
lines( dbh.seq, agb3, lty = "dotted", col = "black", lwd = 2.5)
e =lm(dbh.seq ~ agb3)
summary(e)


#####insert legend

legend("topleft",inset = c(0.015, -0.015), legend = c("B.vulgaris (Eqn. iii) "," B.vulgaris and O.abyssinica (Eqn. v)","B.procera (Huy et al., 2019b)", "B.vulgaris (Nath et al., 2009)"),
       col = c("black","black", "black", "black"), lty = c("solid","twodash","longdash","dotted"), lwd = c(2.5,2,2,2.5), bty = "n" , cex = 0.8)




###################### O. abbysinica #######################################################

oxyt = bamboo.del[which(bamboo.del$Species_name == "O.abyssinica"), ]
(oxyt)

oxyt_select  =  subset(oxyt, select = c( "Plot", "Species_name","Dbh", "AGB"))
View(oxyt_select)


windowsFonts(s = windowsFont("Times New Roman"))

plot(oxyt$Dbh,oxyt$AGB,box(col = "black", lty = "solid", lwd = 2.5), pch =19, col ="black", cex= 1.5, tck = 0.03, cex.axis = 1.3, cex.lab = 1.3, family = "s", bty = "l",  xlab = "Culm diameter-DBH (cm)", ylab = expression(paste("Biomass" ~ ("Kg"  ~ Culm^-1))))

text(x=5,y=5, family="serif")
op = par(family = "serif")

minor.tick (nx = 2, ny = 4, tick.ratio =-0.65)

par(mgp=c(2.5,1,0), las = 0)


dbh.seq <- seq(min(oxyt$Dbh), max(oxyt$Dbh), length.out = 100)


agb5 =-0.1091*dbh.seq^2 + 4.481*dbh.seq - 10.3
lines(dbh.seq,agb5, col = "black", lty = "solid", lwd = 2.5)
j =lm(dbh.seq ~ agb5)
summary(j)



agb7 =0.1776*dbh.seq^2 + 0.1881*dbh.seq + 1.837
lines(dbh.seq, agb7, col= "black", lty = "twodash", lwd = 2.5)
c =lm(dbh.seq ~ agb7)
summary(c)

agb6 = -0.2559*dbh.seq^2+2.8366*dbh.seq-3.9037
lines(dbh.seq, agb6, col = "black", lty = "longdash",lwd  = 2.5)
m =lm(dbh.seq ~ agb6)
summary(m)

agb8 = exp(-0.97+1.68*log(dbh.seq))
lines(dbh.seq,agb8, col = "black", lty = "dotted",lwd = 2.5)
k =lm(dbh.seq ~ agb8)
summary(k)


##gurmesa 2016

windowsFonts(s = windowsFont("Times New Roman"))
legend("topleft",inset = 0.03, legend = c("O.abyssinica (Eqn. iv)", " O.abyssinica and B. vulgaris (Eqn. v)","O.abyssinica (Darcha and Birhane, 2015)", "O.abyssinica (Gurmessa et al., 2016)"), col = c("black", "black", "black", "black"), lty = c("solid","twodash","longdash","dotted"), lwd = c(2.5,2.5,2.5,2.5), bty = "n", cex = 0.8)



