library(dplyr)
library(lme4)

frost_data_full <- read.csv("FrostDamage.csv", header = T)
frost_data<-frost_data_full  

# Summarize data by plot
plots<-group_by(frost_data,P)
damage <- summarise(plots, 
                    damage.mean=mean(LeafDamage, na.rm = TRUE),
                    BA = mean(PlotBasalArea, na.rm = TRUE),
                    T = mean(MinTemp, na.rm = TRUE)
)
damage<-damage[is.finite(damage$damage.mean),]

# least-squares fit of logistic function across plots
# First fit one-parameter model for comparison with logistic
a_start <- 0.5
FitFlat<-nls(damage.mean~1/(1+exp(-1*(a))),
             data=damage, 
             start=list(a=a_start))
# Fit two-parameter logistic model to the relationship between basal area and damage
a_start <- 1
b_start <-0.1
FitLogistic<-nls(damage.mean~1/(1+exp(-1*(a+b*BA))),
                                 data=damage, 
                                 start=list(a=a_start,b=b_start))
summary(FitLogistic)
anova1<-anova(FitLogistic,FitFlat)
#calculate r2 from ANOVA output
r2<-(anova1$'Res.Sum Sq'[2]-anova1$'Res.Sum Sq'[1])/(anova1$'Res.Sum Sq'[2])
print(r2)
plot(damage$BA,predict(FitLogistic))

# compile damage by functional type across plots
plots2<-group_by(frost_data,P,Type)
damage2 <- summarise(plots2, 
                    damage.mean=mean(LeafDamage, na.rm = TRUE),
                    BA = mean(PlotBasalArea, na.rm = TRUE),
                    T = mean(MinTemp, na.rm = TRUE)
)
damage2<-damage2[is.finite(damage2$damage.mean)&is.finite(damage2$Type),]

# fit logistic curve across plots. Parameter c provides test of a difference between savanna and forest spp. at the plot level
a_start <- 1
b_start <-0.1
c_start<-0
FitLogistic2<-nls(damage.mean~1/(1+exp(-1*(a+b*BA+c*(Type=="F")))),
                  data=damage2, 
                  start=list(a=a_start,b=b_start,c=c_start))

summary(FitLogistic2)
plot(damage2$BA,predict(FitLogistic2))


# eliminate species that have at least 10 individuals and some variation in damage
species_counts <- table(frost_data$Species.Code)
species_counts <- species_counts[species_counts>=10] # only use species with N>=10

#some statistics to help exclude species
species_acceptable <- data.frame(Species.Code = unique(frost_data$Species.Code),
                                 N_100 = numeric(length(unique(frost_data$Species.Code))),
                                 N_0 = numeric(length(unique(frost_data$Species.Code))),
                                 N_not = numeric(length(unique(frost_data$Species.Code))),
                                 prop_100 = numeric(length(unique(frost_data$Species.Code))),
                                 prop_0 = numeric(length(unique(frost_data$Species.Code))),
                                 prop_not = numeric(length(unique(frost_data$Species.Code))))

for(i in 1:nrow(species_acceptable)){
  code <- species_acceptable$Species.Code[i]
  geada_table <- table(frost_data[frost_data$Species.Code == code, "dieback.binary"])
  species_acceptable[i, "N_100"] <- ifelse("1" %in% names(geada_table), 
                                           geada_table[names(geada_table) == "1"],
                                           0)
  species_acceptable[i, "N_0"] <- ifelse("0" %in% names(geada_table), 
                                         geada_table[names(geada_table) == "0"],
                                         0)
#  species_acceptable[i, "N_not"] <- sum(geada_table) - species_acceptable[i, "N_100"] - 
#                                                      species_acceptable[i, "N_0"]
  species_acceptable[i, "prop_100"] <- species_acceptable[i, "N_100"] / sum(geada_table)
  species_acceptable[i, "prop_0"] <- species_acceptable[i, "N_0"] / sum(geada_table)
#  species_acceptable[i, "prop_not"] <- species_acceptable[i, "N_not"] / sum(geada_table)
  
}

#exclude species with N<10
frost_data <- frost_data[frost_data$Species.Code %in% names(species_counts)[-1], ]
#exclude species that have >98% with either 100% damage or 0 damage
frost_data <- frost_data[frost_data$Species.Code %in% 
                           species_acceptable[species_acceptable$N_100 >= 2 & species_acceptable$N_0 >= 2, "Species.Code"], ]
frost_data$Species.Code <- droplevels(frost_data$Species.Code)



# linear fit of leaf damage
#frost_temp.lm <- lm(frost_data$LeafDamage ~ frost_data$MinTemp)
#summary(frost_temp.lm)

# Logistic fit of binary damage using the full data (including species excluded above)
frost_temp.binom <- glm(dieback.binary ~ MinTemp, family = "binomial",
                       data = frost_data_full)
summary(frost_temp.binom)

# model to obtain estimate LT50
frost_temp.binom_mix <- glmer(dieback.binary ~ MinTemp + (1|Species.Code), family = "binomial",
                            data = frost_data)
  summary(frost_temp.binom_mix)
  coefs <- coef(frost_temp.binom_mix)$Species.Code
  LT50 <- coefs
  LT50$LT50 <- -(LT50$`(Intercept)`) / (LT50$MinTemp)
  View(LT50)
 
# test for difference in damage between savanna and forest species
frost_temp.binom_nest <- glmer(frost_data$dieback.binary ~ MinTemp + Type + (1|Species.Code), family = "binomial",
                               data = frost_data)
  summary(frost_temp.binom_nest)
  
# include height in model
  frost_temp.binom_nest_height <- glmer(frost_data$dieback.binary ~ MinTemp + Type + InitialHeight +(1|Species.Code), family = "binomial",
                                 data = frost_data)  
  summary(frost_temp.binom_nest_height)

  # Test of species type with plot and species as random factors
  frost_plot.binom_nest <- glmer(frost_data$dieback.binary ~ Type +(1|Species.Code)+(1|P), family = "binomial",
                                  data = frost_data)
  summary(frost_plot.binom_nest)
  

####################################
# Plot-by-plot analysis of frost damage ~ type and height
####################################
library(effects)

frost_data2 <- read.csv("FrostDamage.csv", header = T)
frost_data2$P <- as.factor(frost_data2$P)

#reduce data to extract plot basal area
plot_data <- frost_data2[!duplicated(frost_data2$P) & frost_data2$P %in% c(1:30), ]

models <- vector(mode="list", length=30) #list to store fitted models
predictions <- data.frame(P = c(1:30), #dataframe to store predictions from models for ref. diameter
                          BA = plot_data$PlotBasalArea,
                          predF = rep(0, length = 30),
                          predS = rep(0, length = 30))

for(i in 1:30){ #loop over all plots

  select_data <- frost_data2[frost_data2$P == i, ] #subset data for one plot
  
  
  flag <- FALSE
  if(length(table(select_data$LeafDamage)) == 1){
    flag <- TRUE #flag is used to fit models from plots with more than one value of geada;
                 #i.e. ignore the plots with no frost damage
  }
  
  
  if(flag == FALSE){ #only fit models for plots with variability in frost damage (not all == 0)
    mod <- lm(LeafDamage~InitialHeight*Type, select_data) #fit model
    models[[i]] <- mod #store model in list of models
    #plot(allEffects(mod)) #plot effects
    summary(mod)
    
    #get predictions for a 3-m tall tree of each type
    predicts <- predict(mod, newdata = data.frame(InitialHeight = rep(300,2),Type = c("F", "S")))
    #print(predicts)
    predictions[i, c(3, 4)] <- predicts #store predictions
  }

}

# fit logistic across plots
library(reshape2)
predictions2<- melt(predictions,id.vars = c("P","BA"),
                    variable.name = "Type",
                    value.name = "PredDamage")

predictions2$Type2 <- ifelse(predictions2$Type == "predF", 1, 0)
a_start <- 2
b_start <- -0.24
c_start<- -2.5
PlotLogistic<-nls(PredDamage~1/(1+exp(-1*(a+b*BA+c*(Type=="predF")))),
                  data=predictions2, 
                  start=list(a=a_start,b=b_start,c=c_start))
summary(PlotLogistic)
plot(damage2$BA,predict(PlotLogistic))

#############################################################
#             Phylogenetic analyses                         #
#############################################################
library('ape')
library(geiger)
library (phylolm)
tree1<-read.tree('FrostPhylogeny.txt')
plot(tree1)
Traits <- read.csv("FrostTraits.csv",header=T)

Traits<-Traits[is.finite(Traits$type),]
row.names(Traits)<-Traits[,1]
tree2<-multi2di(tree1, rantom =TRUE, tol = 0.0001)
is.ultrametric(tree2)
is.rooted(tree2)
is.binary.tree(tree2)
plot(tree2)
T50<-Traits$T50
Lat<-Traits$Lat
type<-Traits$type
Bio6<-Traits$Bio6
names(T50)<-Traits$species
names(type)<-Traits$species
names(Lat)<-Traits$species
names(Bio6)<-Traits$species
# non-phylogenetic ANOVA compare LT50 between savanna and forest species
a<-aov(T50~type)
summary(a)
# phylogenetic ANOVA to compare LT50 between savanna and forest species
aa<-aov.phylo(T50~type,tree2,nsim=10000)
summary(aa)
#non-phylogenetic ANOvA to compare southern latitude of savanna and forest species
b<-aov(Lat~type)
summary(b)
# phylogenetic ANOVA to compare southern latitude of savanna and forest species
bb<-aov.phylo(Lat~type,tree2,nsim=10000)
summary(bb)
#non-phylogenetic ANOvA to compare min temperature of range limit of savanna and forest species
c<-aov(Bio6~type)
summary(c)
# phylogenetic ANOVA to compare min temperature of range limit latitude of savanna and forest species
cc<-aov.phylo(Bio6~type,tree2,nsim=100000)
summary(cc)

damage<-damage[is.finite(damage$damage.mean),]

# non phylogenetic regression to test for correlation between LT50 and Latitude
d<-lm(Traits$Lat~Traits$T50)
summary(d)
# phylogenetic regression to test for correlation between LT50 and Latitude
dd<-phylolm(T50~Lat, data = Traits, tree2, model = c("BM"), lower.bound = NULL, upper.bound = NULL, starting.value = NULL, measurement_error = FALSE, boot=0,full.matrix = TRUE)
summary(dd)

# non phylogenetic regression to test for correlation between LT50 and Latitude with type included
e<-lm(Lat~T50+type, data = Traits)
summary(e)
# phylogenetic regression to test for correlation between LT50 and Latitude
ee<-phylolm(Lat~T50+type, data = Traits, tree2, model = c("BM"), lower.bound = NULL, upper.bound = NULL, starting.value = NULL, measurement_error = FALSE, boot=0,full.matrix = TRUE)
summary(ee)

# non phylogenetic regression to test for correlation between LT50 and Bio6
f<-lm(Traits$Bio6~Traits$T50)
summary(f)
# phylogenetic regression to test for correlation between LT50 and Bio6
ff<-phylolm(Bio6~Traits$T50, data = Traits, tree2, model = c("BM"), lower.bound = NULL, upper.bound = NULL, starting.value = NULL, measurement_error = FALSE, boot=0,full.matrix = TRUE)
summary(ff)

# non phylogenetic regression to test for correlation between LT50 and Latitude with type included
e<-lm(Bio6~T50+type, data = Traits)
summary(e)
# phylogenetic regression to test for correlation between LT50 and Latitude
ee<-phylolm(Bio6~T50+type, data = Traits, tree2, model = c("BM"), lower.bound = NULL, upper.bound = NULL, starting.value = NULL, measurement_error = FALSE, boot=0,full.matrix = TRUE)
summary(ee)



e<-phylolm(T50~Bio6, data = Traits, tree2, model = c("BM"), lower.bound = NULL, upper.bound = NULL, starting.value = NULL, measurement_error = FALSE, boot=0,full.matrix = TRUE)
summary(e)

