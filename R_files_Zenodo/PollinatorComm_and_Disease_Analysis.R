# Full Code for "Pollinator community species richness dilutes prevalence of 
# multiple viruses within multiple host species"

# code written by: Michelle L Fearon

# This code includes all analyses, figures, and many of the tables included in the main text and Appendix S1 
# (each is labeled with it's corresponding figure or table label in the text)
# Additionally this code produces a matrix of the relative abundance of all species collected at each site (line 970)


# set the working directory
setwd("C:/Users/mlfea/OneDrive/Documents/Projects/Analyses/PollinatorComm")


#### Abbreviations used in the code below:

## Bee host species
# APME = APIS = Apis = Apis mellifera (honey bees)
# BOIM = BOMB = Bombus = Bombus impatiens (Easter bumblebee)
# LASI = Lasi = Lasioglossum spp. (Halictid sweat bee)
# PEPR = PEPO = Pepo = Eucera pruinosa (previously named Peponapis pruinosa) (squash bee)

## Viruses
# DWV = deformed wing virus
# BQCV = black queen cell virus
# SBV = sacbrood virus


# For additional abbreviations for each species found in the 'comm" matrix of each species relative abundance per site
# please refer to SpeciesCodes.xlsx in the associated DRYAD data.


###########################################
###### GLMM global models for virus prevalence of BQCV, DWV, and SBV in four host species
###########################################


# load binary disease data with additional covariates
disease <- read.csv("Disease_Pos_global_publish.csv", stringsAsFactors = F)

# samples with low RNA concentration and failure to amplify 18S gene were pulled out in excel prior to analsis (explained in MS)
# code for caluclated pollinator community metrics (e.g. richness, abundance, etc.) used in this data set can be found below starting on line 970

head(disease)
dim(disease)



# Function to test for the Variance Inflation Factor, want it to be <6 for all main factors
# vif.mer function developed by Austin F. Frank https://raw.github.com/aufrank/R-hacks/master/mer-utils.R
vif.mer<-function(fit){   ## adapted from rms:vif
  v<-vcov(fit)
  nam<-names(fixef(fit)) ## exclude intercepts
  ns<-sum(1*(nam=="Intercept"|nam=="Intercept"))
  if(ns>0){
    v<-v[-(1:ns),-(1:ns),drop=FALSE]
    nam<-nam[-(1:ns)]
  }
  d<-diag(v)^(0.5)
  v<-diag(solve(v/(d%o%d)))
  names(v)<-nam
  v
}

# Transformation functions to use within models developed by Austin F. Frank https://github.com/aufrank/R-hacks/blob/master/regression-utils.R
# Z transform all main factors (scale and centers variables)
z. <- function (x) scale(x)

# log transform all abundance factors
l. <- function (x) log(x)


# overdispersion test for glmer models
# if overdispersed in a glm --> use quasipoisson and an F test to compare models
# if overdispersed in a glmer --> add ID as a random effect

# Overdispersion parameter estimation function from Ben Bolker: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor
# Comment additions from http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model)) - model.df
  # extracts the Pearson residuals
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}



# load libraries
library(MuMIn) 
library(car)
library(lme4)
library(ggeffects)
library(DHARMa)
library(ggplot2)



# modify factors to correct class
disease$Year <- as.factor(disease$Year)
disease$Site <- as.factor(disease$Site)
disease$VisitNumber <- as.factor(disease$VisitNumber)
disease$Transect.ID <- as.factor(disease$Transect.ID)
disease$Type <- as.factor(disease$Type)
disease$Virus_Type <- as.factor(disease$Virus_Type)
disease$Species <- as.factor(disease$Species)
disease$Genus <- as.factor(disease$Genus)
disease$Date_Collected <- as.factor(disease$Date.Collected)
disease$NetorPan <- as.factor(disease$NetorPan)



# Analysis in Appendix S1 Table S2
# Test of RNA concentration or collection method are correlated with virus presence/absence
###############################

# samples with "too high" RNA were beyond the range of detection and above 600 ng/ul.
# set "too high" to 600 as a conservative estimate
disease$RNA.conc.[disease$RNA.conc. == "too high"] <- 600  

# samples that were listed as "too low" were below the range of detection, and about 1 ng/ul
disease$RNA.conc.[disease$RNA.conc. == "Sample too low"] <- 1
disease$RNA.conc. <- as.numeric(disease$RNA.conc.)

# mean and range of RNA concentrations
mean(disease$RNA.conc., na.rm = T)
range(disease$RNA.conc., na.rm = T)


# check of virus differences based on RNA conc or collection method
mod <- glmer(Virus ~ NetorPan + z.(RNA.conc.) + (1|Virus_Type) + (1|Genus) + (1|VisitNumber:Site), family = binomial, data = disease)
summary(mod)  # neither RNA concentration or collection method are significant
Anova(mod)




# Scale variables
disease$Richness_z <- as.numeric(z.(disease$Richness))
disease$EstRich_Asyp_z <- as.numeric(z.(disease$EstRich_Asyp))
disease$EstRich_46_z <- as.numeric(z.(disease$EstRich_46))
disease$EstRich_183_z <- as.numeric(z.(disease$EstRich_183))
disease$EstRich_338_z <- as.numeric(z.(disease$EstRich_338))
disease$Abundance_z <- as.numeric(z.(l.(disease$Abundance)))  # all abundance variables are natural log and z transformed
disease$APME_ABUND_z <- as.numeric(z.(l.(disease$APME_ABUND)))
disease$BOIM_ABUND_z <- as.numeric(z.(l.(disease$BOIM_ABUND)))
disease$LASI_ABUND_z <- as.numeric(z.(l.(disease$LASI_ABUND)))
disease$PEPR_ABUND_z <- as.numeric(z.(l.(disease$PEPR_ABUND+1))) # one site had zero PEPR, so needed to add +1 to prevent inf in the log transformation
summary(disease)

#Pearson correlation
cor(z.(disease$Richness), z.(disease$Abundance))
cor(z.(disease$Richness), z.(disease$APME_ABUND))
cor(z.(disease$Richness), z.(disease$BOIM_ABUND))




### Global model of virus prevalence
####################################

# Ran models with no interactions, and all combinations of two-way, three-way, and four-way interactions
# Then used a model selection table to find the top models
# Residuals of all models were checked for spatial autocorrelation


#### model with no interactions
mod1 <- glmer(Virus ~ Abundance_z + Richness_z + Virus_Type + Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
# increase the number of iterations to avoid convergence issue
ss_1 <- getME(mod1,c("theta","fixef"))
mod1_2 <- update(mod1,start=ss_1,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod1_2)
overdisp_fun(mod1_2)


##### two-way interaction models for all pairs of main effects
mod2a <- glmer(Virus ~ Abundance_z + Richness_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
# increase the number of iterations to avoid convergence issue
ss_2a <- getME(mod2a,c("theta","fixef"))
mod2a_2 <- update(mod2a,start=ss_2a,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2a_2)

mod2b <- glmer(Virus ~ Abundance_z + Richness_z * Virus_Type + Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_2b <- getME(mod2b,c("theta","fixef"))
mod2b_2 <- update(mod2b,start=ss_2b,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2b_2)

mod2c <- glmer(Virus ~ Abundance_z * Richness_z + Virus_Type + Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_2c <- getME(mod2c,c("theta","fixef"))
mod2c_2 <- update(mod2c,start=ss_2c,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2c_2)

mod2d <- glmer(Virus ~ Richness_z + Abundance_z * Virus_Type + Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_2d <- getME(mod2d,c("theta","fixef"))
mod2d_2 <- update(mod2d,start=ss_2d,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2d_2)

mod2e <- glmer(Virus ~ Abundance_z + Richness_z * Genus + Virus_Type + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_2e <- getME(mod2e,c("theta","fixef"))
mod2e_2 <- update(mod2e,start=ss_2e,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2e_2)

mod2f <- glmer(Virus ~ Richness_z + Abundance_z * Genus + Virus_Type + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_2f <- getME(mod2f,c("theta","fixef"))
mod2f_2 <- update(mod2f,start=ss_2f,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2f_2)




#### three-way interactions for all sets of 3 variables
mod3a <- glmer(Virus ~ Abundance_z + Richness_z * Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_3a <- getME(mod3a,c("theta","fixef"))
mod3a_2 <- update(mod3a,start=ss_3a,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod3a_2)

mod3b <- glmer(Virus ~ Richness_z + Abundance_z * Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_3b <- getME(mod3b,c("theta","fixef"))
mod3b_2 <- update(mod3b,start=ss_3b,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod3b)

mod3c <- glmer(Virus ~ Richness_z * Abundance_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_3c <- getME(mod3c,c("theta","fixef"))
mod3c_2 <- update(mod3c,start=ss_3c,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod3c_2)

mod3d <- glmer(Virus ~ Richness_z * Abundance_z * Virus_Type + Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_3d <- getME(mod3d,c("theta","fixef"))
mod3d_2 <- update(mod3d,start=ss_3d,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod3d_2)

mod3e <- glmer(Virus ~ Richness_z * Abundance_z * Genus + Virus_Type + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_3e <- getME(mod3e,c("theta","fixef"))
mod3e_2 <- update(mod3e,start=ss_3e,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod3e_2)





##### four way interaction model
mod4 <- glmer(Virus ~ Abundance_z * Richness_z * Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease) # convergence issue
ss_4 <- getME(mod4,c("theta","fixef"))
mod4_2 <- update(mod4,start=ss_4,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod4_2)



### MODEL SELECTION TABLE (Appendix S1 Table S8)
out.put_all <- model.sel(mod1_2, mod2a_2, mod2b_2, mod2c_2, mod2d_2, mod2e_2, mod2f_2, mod3a_2, mod3b_2, mod3c_2, mod3d_2, mod3e_2, mod4_2)
write.csv(out.put_all, "Model Selection Table_all.csv")
# Mod2a_2 and Mod3a_2 are the best (described as Mod 2a and Mod 3a in the text)

### Subset of model selection Table 1 in Manuscript
out.put_subset <- model.sel(mod2a_2, mod3a_2, mod3b_2, mod3c_2)
write.csv(out.put_subset, "Model Selection Table_subset.csv")



# Tests of assumptions for Mod 2a_2
plot(mod2a_2)
qqnorm(resid(mod2a_2))
qqline(resid(mod2a_2))

#Test for overdispersion
overdisp_fun(mod2a_2)
print(mod2a_2, correlation = T)

# variance-covariance matrix for the model
vcov(mod2a_2)

# check if singular
isSingular(mod2a_2)





# TABLE 3 in Manuscript: Model output for Mod 2a
Anova(mod2a_2)


# Appendix S1 Table S9: Model output for Mod 3a
Anova(mod3a_2)


# Appendix S1 Table S10: VIF for main factors in Mod 2a
vif.mer(mod2a_2)



# Appendix S1 Table S11: 
# model with richness, species specific abundances and virus * genus interaction
mod5 <- glmer(Virus ~ Richness_z + APME_ABUND_z + BOIM_ABUND_z + LASI_ABUND_z + PEPR_ABUND_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_5 <- getME(mod5,c("theta","fixef"))
mod5_2 <- update(mod5,start=ss_5,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

summary(mod5_2) # none of the species specifc abundances are significant
vif.mer(mod5_2) 
Anova(mod5_2) # used for Table S11



# check models with each host species' abundance separately - SAME RESULTS AS ABOVE (none of the abundance factors are significant)
# Apis abundance
mod5a <- glmer(Virus ~ Richness_z + APME_ABUND_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_5a <- getME(mod5a,c("theta","fixef"))
mod5a_2 <- update(mod5a,start=ss_5a,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

summary(mod5a_2)
vif.mer(mod5a_2) 

# bombus abundance
mod5b <- glmer(Virus ~ Richness_z + BOIM_ABUND_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_5b <- getME(mod5b,c("theta","fixef"))
mod5b_2 <- update(mod5b,start=ss_5b,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

summary(mod5b_2)
vif.mer(mod5b_2) 


# LASI abundance
mod5c <- glmer(Virus ~ Richness_z + LASI_ABUND_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_5c <- getME(mod5c,c("theta","fixef"))
mod5c_2 <- update(mod5c,start=ss_5c,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

summary(mod5c_2)
vif.mer(mod5c_2) 


# PEPO abundance
mod5d <- glmer(Virus ~ Richness_z + PEPR_ABUND_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_5d <- getME(mod5d,c("theta","fixef"))
mod5d_2 <- update(mod5d,start=ss_5d,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

summary(mod5d_2)
vif.mer(mod5d_2) 





# Appendix S1 Table S14:
# model with asymptote estimated species richness and virus * genus interaction
mod6 <- glmer(Virus ~ EstRich_Asyp_z + Abundance_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_6 <- getME(mod6,c("theta","fixef"))
mod6_2 <- update(mod6,start=ss_6,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

summary(mod6_2)
vif.mer(mod6_2) 
Anova(mod6_2)  # used for Table S14


# Appendix S1 Table S15:
# model with estimated species richness at 46 individuals and virus * genus interaction
# 46 is the lowest number of individuals collected at a site
mod7 <- glmer(Virus ~ EstRich_46_z + Abundance_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_7 <- getME(mod7,c("theta","fixef"))
mod7_2 <- update(mod7,start=ss_7,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

summary(mod7_2)
vif.mer(mod7_2) 
Anova(mod7_2) # used for Table S15


# [not shown in text]

# model with estimated species richness at 183 individuals and virus * genus interaction
# 183 is the second lowest number of individuals collected at a site
mod8 <- glmer(Virus ~ EstRich_183_z + Abundance_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_8 <- getME(mod8,c("theta","fixef"))
mod8_2 <- update(mod8,start=ss_8,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

summary(mod8_2)
vif.mer(mod8_2) 


# model with estimated species richness at 338 individuals and virus * genus interaction
# 338 is the mean number of individuals collected at each site
mod9 <- glmer(Virus ~ EstRich_338_z + Abundance_z + Virus_Type * Genus + (1|VisitNumber:Site) + (1|Individual.ID), family = binomial, data = disease)
ss_9 <- getME(mod9,c("theta","fixef"))
mod9_2 <- update(mod9,start=ss_9,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

summary(mod9_2)
vif.mer(mod9_2) 





###### Moran's I spatial autocorrelation tests, stats for Table S12 in Appendix S1
##################################
## Appears in Appendix S1 Table S12

# none of the tests of spatial autocorrelation were significant

#Test for spatial autocorrelation
mod1_resid <- simulateResiduals(mod1)
mod1_resid2 <- recalculateResiduals(mod1_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod1_resid2, unique(disease$Long), unique(disease$Lat))

mod2a_resid <- simulateResiduals(mod2a_2, n = 5000, integerResponse = TRUE) 
mod2a_resid2 <- recalculateResiduals(mod2a_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod2a_resid2, unique(disease$Long), unique(disease$Lat))

mod2b_resid <- simulateResiduals(mod2b_2, n = 5000, integerResponse = TRUE) 
mod2b_resid2 <- recalculateResiduals(mod2b_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod2b_resid2, unique(disease$Long), unique(disease$Lat))

mod2c_resid <- simulateResiduals(mod2c_2, n = 5000, integerResponse = TRUE) 
mod2c_resid2 <- recalculateResiduals(mod2c_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod2c_resid2, unique(disease$Long), unique(disease$Lat))

mod2d_resid <- simulateResiduals(mod2d_2, n = 5000, integerResponse = TRUE) 
mod2d_resid2 <- recalculateResiduals(mod2d_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod2d_resid2, unique(disease$Long), unique(disease$Lat))

mod2e_resid <- simulateResiduals(mod2e_2, n = 5000, integerResponse = TRUE) 
mod2e_resid2 <- recalculateResiduals(mod2e_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod2e_resid2, unique(disease$Long), unique(disease$Lat))

mod2f_resid <- simulateResiduals(mod2f_2, n = 5000, integerResponse = TRUE) 
mod2f_resid2 <- recalculateResiduals(mod2f_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod2f_resid2, unique(disease$Long), unique(disease$Lat))


mod3a_resid <- simulateResiduals(mod3a_2, n = 5000, integerResponse = TRUE)
mod3a_resid2 <- recalculateResiduals(mod3a_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod3a_resid2, unique(disease$Long), unique(disease$Lat))

mod3b_resid <- simulateResiduals(mod3b_2, n = 5000, integerResponse = TRUE)
mod3b_resid2 <- recalculateResiduals(mod3b_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod3b_resid2, unique(disease$Long), unique(disease$Lat))

mod3c_resid <- simulateResiduals(mod3c_2, n = 5000, integerResponse = TRUE)
mod3c_resid2 <- recalculateResiduals(mod3c_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod3c_resid2, unique(disease$Long), unique(disease$Lat))

mod3d_resid <- simulateResiduals(mod3d_2, n = 5000, integerResponse = TRUE)
mod3d_resid2 <- recalculateResiduals(mod3d_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod3d_resid2, unique(disease$Long), unique(disease$Lat))

mod3e_resid <- simulateResiduals(mod3e_2, n = 5000, integerResponse = TRUE)
mod3e_resid2 <- recalculateResiduals(mod3e_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod3e_resid2, unique(disease$Long), unique(disease$Lat))



mod4_resid <- simulateResiduals(mod4_2, n = 5000, integerResponse = TRUE)
mod4_resid2 <- recalculateResiduals(mod4_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod4_resid2, unique(disease$Long), unique(disease$Lat))


mod5_resid <- simulateResiduals(mod5_2, n = 1000)
mod5_resid2 <- recalculateResiduals(mod5_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod5_resid2, unique(disease$Long), unique(disease$Lat))


mod6_resid <- simulateResiduals(mod6_2, n = 1000)
mod6_resid2 <- recalculateResiduals(mod6_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod6_resid2, unique(disease$Long), unique(disease$Lat))


mod7_resid <- simulateResiduals(mod7_2, n = 1000)
mod7_resid2 <- recalculateResiduals(mod7_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod7_resid2, unique(disease$Long), unique(disease$Lat))








# Load libraries
library(grid)
library(gridExtra)
library(reshape)
library(ggiraph)
library(ggiraphExtra)
library(dplyr)
library(moonBook)
library(ggeffects)
library(colorBlindness)
library(emmeans)


# function to make a grid plot with a shared legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}



# FIGURE 1 in Manuscript, Appendix S1 Table S17 and S18
#################################
# Plot of virus prevalence by host species and virus


# calculates the pairwise tests for each genus within each virus type combination
a <- emmeans(mod2a_2, specs = pairwise ~ Genus|Virus_Type, type = "response")
a



# data frame of pairwise tests and predicted prevalence
PrevByHost <- a$emmeans %>%
  as.data.frame()
PrevByHost

# Appendix S1 Table S17 
write.csv(PrevByHost, file = "Predicted_Prevalence_By_Host_and_Virus.csv", quote = F)


pairwise <- a$contrasts %>%
  summary(infer = TRUE) %>%
  as.data.frame()
pairwise

# Appendix S1 Table S18 
write.csv(pairwise, file = "Pairwise_Tests_Prevalence_By_Host_and_Virus.csv", quote = F)





# color blind palette (from https://zenodo.org/record/3381072#.X8VKS2hKjZt)
Tol_muted <- c('#882255','#44AA99', '#117733', '#332288', '#DDCC77','#CC6677', '#88CCEE', '#AA4499', '#DDDDDD', '#999933')
        # Goedhart, Joachim. (2019, August 29). Material related to the blog "Dataviz with Flying Colors". Zenodo. http://doi.org/10.5281/zenodo.3381072



# FIGURE 1 in Manuscript
me_inter <- ggpredict(mod2a_2, c("Virus_Type", "Genus"))
prev_host <- plot(me_inter, dodge = 0.4) +
  labs(x = "Virus", y = "Virus prevalence", title = NULL) +
  ylim(0,1) +
  scale_color_manual(values = Tol_muted, name = "Species", labels = c(bquote(italic("Apis mellifera")), bquote(italic("Bombus impatiens")), 
                                                                    bquote(italic("Lasioglossum")~ "spp."), bquote(italic("Eucera pruinosa")))) +
  scale_shape_discrete(name = "Species", labels = c(bquote(italic("Apis mellifera")), bquote(italic("Bombus impatiens")), 
                                                  bquote(italic("Lasioglossum")~ "spp."), bquote(italic("Eucera pruinosa")))) +
  geom_text(aes(y = predicted + 0.13), label = c("a","b","c","c","a", "b", "c", "c", "a", "a", "b", "b"), position = position_dodge(width = 0.4), show.legend = F, size = 9/.pt) +
  theme_classic() +
  theme(legend.position = c(0.75,0.8), legend.box.background = element_rect(), axis.text = element_text(size = 9, color = "black"), axis.title = element_text(size = 11, color = "black"), 
        legend.text = element_text(size = 7), legend.title = element_text(size = 10))
print(prev_host)
ggsave("Fig1.tiff", plot = prev_host, dpi = 600, width = 3, height = 3.5, units = "in", compression="lzw")

# color blind check
cvdPlot(prev_host)



# Note: FIGURE 2 code is below with the compilation of the pollinator community data (line 1237)



# FIGURE 3A in Manuscript

# Marginal effects from best model with Richness on the x axis
me_rich <- ggpredict(mod2a_2, c("Richness_z", "Genus", "Virus_Type"))
plot(me_rich, add.data = F, jitter = 0.05) # initial plot (aim to replicate this)

# recalculate original richness 
me_rich$x # scaled richness
rich_mean <- mean(disease$Richness) # mean of original richness from disease data set
rich_sd <- sd(disease$Richness) # sd of original richness from disease data set

me_rich$Rich <- t((t(me_rich$x) * rich_sd) + rich_mean)

# plot of Richness vs virus prevalence by virus and by species (with original richness scale)
virus_rich <- ggplot(me_rich, aes(x = Rich, y = predicted)) +
  scale_fill_manual(values = Tol_muted, name = "Species", labels = c(bquote(italic("Apis mellifera")), bquote(italic("Bombus impatiens")), 
                                                                   bquote(italic("Lasioglossum")~ "spp."), bquote(italic("Eucera pruinosa")))) +
  scale_color_manual(values = Tol_muted, name = "Species", labels = c(bquote(italic("Apis mellifera")), bquote(italic("Bombus impatiens")), 
                                                                    bquote(italic("Lasioglossum")~ "spp."), bquote(italic("Eucera pruinosa")))) +
  labs(color = "Species", tag = "a") +
  geom_line(aes(color = group), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, outline.type = NULL) +
  labs(x = "Species Richness", y = "Virus Prevalence") +
  facet_wrap(~facet) +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=10, color = "black"))+
  theme(axis.title = element_text(size=12, color="black"))
print(virus_rich)
#ggsave("VirusPrev_Richness_predict.png", plot = virus_rich, dpi = 300, width = 18, height = 8, units = "cm")




# FIGURE 3B in Manuscript

# Marginal effects from best model with Abundance on the x axis
me_abund <- ggpredict(mod2a_2, c("Abundance_z", "Genus", "Virus_Type"))
plot(me_abund, add.data = F, jitter = 0.05) # initial plot (aim to replicate this)

# recalculate original Abundance 
me_abund$x # scaled richness
abund_mean <- mean(l.(disease$Abundance)) # mean of original richness from disease data set
abund_sd <- sd(l.(disease$Abundance)) # sd of original richness from disease data set
me_abund$Abund_log <- t((t(me_abund$x) * abund_sd) + abund_mean)
me_abund$Abund <- exp(me_abund$Abund_log)

range(me_abund$Abund)
cbind(disease$Site,disease$Abundance_z, disease$Abundance)

# plot of Richness vs virus prevalence by virus and by species (with original richness scale)
virus_abund <- ggplot(me_abund, aes(x = Abund, y = predicted)) +
  scale_fill_manual(values = Tol_muted, name = "Species", labels = c(bquote(italic("Apis mellifera")), bquote(italic("Bombus impatiens")), 
                                                                   bquote(italic("Lasioglossum")~ "spp."), bquote(italic("Eucera pruinosa")))) +
  scale_color_manual(values = Tol_muted, name = "Species", labels = c(bquote(italic("Apis mellifera")), bquote(italic("Bombus impatiens")), 
                                                                    bquote(italic("Lasioglossum")~ "spp."), bquote(italic("Eucera pruinosa")))) +
  labs(color = "Species", tag = "b") +
  geom_line(aes(color = group), size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, outline.type = NULL) +
  labs(x = "Total Pollinator Abundance", y = "Virus Prevalence") +
  facet_wrap(~facet) +
  ylim(0,1) +
  xlim(46,800) +
  scale_x_log10(breaks = c(50,100,300,600)) +
  theme_classic() +
  theme(axis.text = element_text(size=10, color = "black"))+
  theme(axis.title = element_text(size=12, color="black"))
print(virus_abund)
#ggsave("VirusPrev_Abundance_predict.png", plot = virus_abund, dpi = 300, width = 18, height = 8, units = "cm")

# color blind check
cvdPlot(virus_rich)
cvdPlot(virus_abund)



### FULL FIGURE 3 in Manuscript
# Grid of both richness and abundance plots
virus_rich_abund <- grid_arrange_shared_legend(virus_rich, virus_abund, nrow=2, ncol = 1, position = "right")
ggsave("Fig3.tiff", plot = virus_rich_abund, dpi = 600, width = 6, height = 6, units = "in", compression="lzw")











################################
# INFECTION PREVALENCE CALCULATIONS, ANALYESES, AND TABLES
################################


# load epiR library
library(epiR)



######## VIRUS INFECTION PREVALENCE CALCULATIONS 
## (based on the presence of the negative strand, indicating active replication in host)



# NEGATIVE STRAND DWV
###############################
# load binary disease data with additional covariates
negDWV <- read.csv("NegStdDWV_publish.csv", stringsAsFactors = F)
head(negDWV)
dim(negDWV)


### total number of negative strand positive samples for each virus
sum(negDWV$NegDWV)  #27 samples were positive for the DWV negative strand


# DWV Negative Strand prevalence calculations
GENUS_DWV <- sort(unique(negDWV$Genus))

# make a matrix of apparent prevalences and lower and upper 95% confidence intervals (in percents)
DWVnegprev <- matrix(NA, nrow = 4, ncol = 5)
rownames(DWVnegprev) <- GENUS_DWV
colnames(DWVnegprev) <- c("DWV", "DWV_L", "DWV_U", "N_positive", "N_tested")


# loop through each genus, and calculate the apparent prevalence for all three viruses
for (i in 1:4){
  # subset data for each genus
  temp <- negDWV[ negDWV$Genus == GENUS_DWV[i], ]
  # calculate the number of infected ind for each virus (store in N_infected matrix)
  infectedDWV <- sum(temp$NegDWV)
  DWVnegprev[ i , "N_positive"] <- infectedDWV
  # calculate the total number tested
  tested <- sum(temp$S18Retest)
  DWVnegprev[ i , "N_tested"] <- tested
  # unlist prevalence estimates to make easier to work with
  DWV <- unlist(epi.prev(pos = infectedDWV, tested = tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))
  
  # store virus prevalence and lower and upper confidence interval in prev matrix
  DWVnegprev[ i, "DWV"] <- DWV[1]
  DWVnegprev[ i, "DWV_L"] <- DWV[2]
  DWVnegprev[ i, "DWV_U"] <- DWV[3]
}

#check matrix and export as CSV
DWVnegprev
write.csv(DWVnegprev, file = "DWVNegStrdPrevBySpecies.csv", quote = F)



# NEGATIVE STRAND BQCV
###############################
# load binary disease data with additional covariates
negBQCV <- read.csv("NegStdBQCV_publish.csv", stringsAsFactors = F)
head(negBQCV)
dim(negBQCV)


### total number of negative strand positive samples for each virus
sum(negBQCV$NegBQCV)  #41 samples were positive for the BQCV negative strand


# BQCV Negative Strand prevalence calculations
GENUS_BQCV <- sort(unique(negBQCV$Genus))

# make a matrix of apparent prevalences and lower and upper 95% confidence intervals (in percents)
BQCVnegprev <- matrix(NA, nrow = 4, ncol = 5)
rownames(BQCVnegprev) <- GENUS_BQCV
colnames(BQCVnegprev) <- c("BQCV", "BQCV_L", "BQCV_U", "N_positive", "N_tested")



# loop through each genus, and calculate the apparent prevalence for all three viruses
for (i in 1:4){
  # subset data for each genus
  temp <- negBQCV[ negBQCV$Genus == GENUS_BQCV[i], ]
  # calculate the number of infected ind for each virus (store in N_infected matrix)
  infectedBQCV <- sum(temp$NegBQCV)
  BQCVnegprev[ i , "N_positive"] <- infectedBQCV
  # calculate the total number tested
  tested <- sum(temp$S18Retest)
  BQCVnegprev[ i , "N_tested"] <- tested
  # unlist prevalence estimates to make easier to work with
  BQCV <- unlist(epi.prev(pos = infectedBQCV, tested = tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))
  
  # store virus prevalence and lower and upper confidence interval in prev matrix
  BQCVnegprev[ i, "BQCV"] <- BQCV[1]
  BQCVnegprev[ i, "BQCV_L"] <- BQCV[2]
  BQCVnegprev[ i, "BQCV_U"] <- BQCV[3]
}

#check matrix and export as CSV
BQCVnegprev
write.csv(BQCVnegprev, file = "BQCVNegStrdPrevBySpecies.csv", quote = F)



# NEGATIVE STRAND SBV
###############################
# load binary disease data with additional covariates
negSBV <- read.csv("NegStdSBV_publish.csv", stringsAsFactors = F)
head(negSBV)
dim(negSBV)

negSBV$S18Retest

### total number of negative strand positive samples for each virus
sum(negSBV$NegSBV)  #48 samples were positive for the SBV negative strand


# SBV Negative Strand prevalence calculations
GENUS_SBV <- sort(unique(negSBV$Genus))

# make a matrix of apparent prevalences and lower and upper 95% confidence intervals (in percents)
SBVnegprev <- matrix(NA, nrow = 4, ncol = 5)
rownames(SBVnegprev) <- GENUS_SBV
colnames(SBVnegprev) <- c("SBV", "SBV_L", "SBV_U", "N_positive", "N_tested")


# loop through each genus, and calculate the apparent prevalence for all three viruses
for (i in 1:4){
  # subset data for each genus
  temp <- negSBV[ negSBV$Genus == GENUS_SBV[i], ]
  # calculate the number of infected ind for each virus (store in N_infected matrix)
  infectedSBV <- sum(temp$NegSBV)
  SBVnegprev[ i , "N_positive"] <- infectedSBV
  # calculate the total number tested
  tested <- sum(temp$S18Retest)
  SBVnegprev[ i , "N_tested"] <- tested
  # unlist prevalence estimates to make easier to work with
  SBV <- unlist(epi.prev(pos = infectedSBV, tested = tested, method = "blaker", conf.level = 0.95, sp = 1, se = 0.95))
  
  # store virus prevalence and lower and upper confidence interval in prev matrix
  SBVnegprev[ i, "SBV"] <- SBV[1]
  SBVnegprev[ i, "SBV_L"] <- SBV[2]
  SBVnegprev[ i, "SBV_U"] <- SBV[3]
}

#check matrix and export as CSV
SBVnegprev
write.csv(SBVnegprev, file = "SBVNegStrdPrevBySpecies.csv", quote = F)




### TABLE 2 Infection prevalence and 95% confidence interval data in manuscript, and 
# Appendix S1: TABLE S7 with numbers of each host species tested for each virus
###########################
DWVnegprev
BQCVnegprev
SBVnegprev





###### Pairwise Chi-squared tests of INFECTION prevalence (2 proportions tests)
#######################################
# Based on the presence of the negative strand in virus-positive samples


##### Appendix S1, Table S19 shows the p-values from these tests



####### Bonferoni Correction for Multiple comparisons (6 comparisons among 4 hosts per virus)
alpha <- 0.05/6
# threshold for significance is alpha = 0.008333

# P-values recorded in TABLE S19 in Appendix S1

##### Tests of two proportions for DWV infection prevalence

#Apis vs Bombus DWV
infected <- c(DWVnegprev[ "APIS", "N_positive"], DWVnegprev[ "BOMB", "N_positive"])
total <- c(DWVnegprev[ "APIS", "N_tested"], DWVnegprev[ "BOMB", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Apis vs Lasioglossum DWV
infected <- c(DWVnegprev[ "APIS", "N_positive"], DWVnegprev[ "LASI", "N_positive"])
total <- c(DWVnegprev[ "APIS", "N_tested"], DWVnegprev[ "LASI", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Apis vs Peponapis DWV
infected <- c(DWVnegprev[ "APIS", "N_positive"], DWVnegprev[ "PEPO", "N_positive"])
total <- c(DWVnegprev[ "APIS", "N_tested"], DWVnegprev[ "PEPO", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Bombus vs Lasioglossum DWV
infected <- c(DWVnegprev[ "BOMB", "N_positive"], DWVnegprev[ "LASI", "N_positive"])
total <- c(DWVnegprev[ "BOMB", "N_tested"], DWVnegprev[ "LASI", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Bombus vs Peponapis DWV
infected <- c(DWVnegprev[ "BOMB", "N_positive"], DWVnegprev[ "PEPO", "N_positive"])
total <- c(DWVnegprev[ "BOMB", "N_tested"], DWVnegprev[ "PEPO", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Lasioglossum vs Peponapis DWV
infected <- c(DWVnegprev[ "LASI", "N_positive"], DWVnegprev[ "PEPO", "N_positive"])
total <- c(DWVnegprev[ "LASI", "N_tested"], DWVnegprev[ "PEPO", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)



##### Tests of two proportions for BQCV infection prevalence

#Apis vs Bombus BQCV
infected <- c(BQCVnegprev[ "APIS", "N_positive"], BQCVnegprev[ "BOMB", "N_positive"])
total <- c(BQCVnegprev[ "APIS", "N_tested"], BQCVnegprev[ "BOMB", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Apis vs Lasioglossum BQCV
infected <- c(BQCVnegprev[ "APIS", "N_positive"], BQCVnegprev[ "LASI", "N_positive"])
total <- c(BQCVnegprev[ "APIS", "N_tested"], BQCVnegprev[ "LASI", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Apis vs Peponapis BQCV
infected <- c(BQCVnegprev[ "APIS", "N_positive"], BQCVnegprev[ "PEPO", "N_positive"])
total <- c(BQCVnegprev[ "APIS", "N_tested"], BQCVnegprev[ "PEPO", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Bombus vs Lasioglossum BQCV
infected <- c(BQCVnegprev[ "BOMB", "N_positive"], BQCVnegprev[ "LASI", "N_positive"])
total <- c(BQCVnegprev[ "BOMB", "N_tested"], BQCVnegprev[ "LASI", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Bombus vs Peponapis BQCV
infected <- c(BQCVnegprev[ "BOMB", "N_positive"], BQCVnegprev[ "PEPO", "N_positive"])
total <- c(BQCVnegprev[ "BOMB", "N_tested"], BQCVnegprev[ "PEPO", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Lasioglossum vs Peponapis BQCV
infected <- c(BQCVnegprev[ "LASI", "N_positive"], BQCVnegprev[ "PEPO", "N_positive"])
total <- c(BQCVnegprev[ "LASI", "N_tested"], BQCVnegprev[ "PEPO", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)


##### Tests of two proportions for SBV infection prevalence

#Apis vs Bombus SBV
infected <- c(SBVnegprev[ "APIS", "N_positive"], SBVnegprev[ "BOMB", "N_positive"])
total <- c(SBVnegprev[ "APIS", "N_tested"], SBVnegprev[ "BOMB", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Apis vs Lasioglossum SBV
infected <- c(SBVnegprev[ "APIS", "N_positive"], SBVnegprev[ "LASI", "N_positive"])
total <- c(SBVnegprev[ "APIS", "N_tested"], SBVnegprev[ "LASI", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Apis vs Peponapis SBV
infected <- c(SBVnegprev[ "APIS", "N_positive"], SBVnegprev[ "PEPO", "N_positive"])
total <- c(SBVnegprev[ "APIS", "N_tested"], SBVnegprev[ "PEPO", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Bombus vs Lasioglossum SBV
infected <- c(SBVnegprev[ "BOMB", "N_positive"], SBVnegprev[ "LASI", "N_positive"])
total <- c(SBVnegprev[ "BOMB", "N_tested"], SBVnegprev[ "LASI", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Bombus vs Peponapis SBV
infected <- c(SBVnegprev[ "BOMB", "N_positive"], SBVnegprev[ "PEPO", "N_positive"])
total <- c(SBVnegprev[ "BOMB", "N_tested"], SBVnegprev[ "PEPO", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)

#Lasioglossum vs Peponapis SBV
infected <- c(SBVnegprev[ "LASI", "N_positive"], SBVnegprev[ "PEPO", "N_positive"])
total <- c(SBVnegprev[ "LASI", "N_tested"], SBVnegprev[ "PEPO", "N_tested"])
prop.test(infected,total, alternative=c("two.sided"), conf.level = 0.95)






################################
# POLLINATOR COMMUNITY CALCULATIONS, ANALYSES, FIGURES, AND TABLES
################################

# Read in Pollinator Community data
data <- read.csv("PollinatorComm2015_2016_publish.csv", sep=",", header=T, stringsAsFactors=F)

# check data set
summary(data)
head(data)
dim(data)


# load vegan package
library(vegan)

# 2015 total = 1490
sum(data$Year == "2015")

# 2016 total = 3247
sum(data$Year == "2016")

# number of different species = 127
length(unique(data$Code))
spp_codes <- sort(unique(data$Code))
table(data$Code)

# number of different genera = 80 (78 when "unknown" and "unknown wasp" excluded)
length(unique(data$Genus))
ugenera <- sort(unique(data$Genus))
table(data$Genus)

# number of different families = 18
length(unique(data$Family))

table(data$Family)


# vector of unique field sites = 14 sites
sites <- unique(data$Site)
length(sites)

sitevisit <- unique(data$SiteVisit)



##### Create a matrix of site x species
# where each element is number of catpures
# of species at a given site

# make an empty matrix (14 rows x 127 columns)
comm <- matrix(0, nrow = length(sites), ncol = length(spp_codes))
dim(comm)

# add a rownames to comat matrix

rownames(comm) <- sites
colnames(comm) <- spp_codes
comm

# loop over sites
# use table to get occurence counts
# use names to assign to community matrix


for (zz in sites){
  
  counts <- table(data$Code[ data$Site == zz ])
  comm[ as.character(zz), names(counts) ] <- counts
  
}

comm
sum(comm)  # should equal 4737


# export data to Community Matrix file
write.csv(comm, file = "CommunityMatrix.csv", quote = F)








### TABLE S13 in Supplementary Material (first 4 columns: sites, richness, abundance, simpson index)
#####################################
# remainder of Table 13 is generated below at lines and 1285 and 1301


## Create a summary data set with calcualtions for species richness, total abundance, and simpson diversity index

# for each site calculate:
# species richness
# total abundance
# simpson diversity index (1-D)


# create empty matrix to store data
summary <- matrix(NA, nrow = length(sites), ncol = 5)
rownames(summary) <- sites
colnames(summary) <- c("richness", "abundance", "simpson", "invsimpson", "shannon")
summary

# mean species richness, abundance and simpson diversity across all sites
mean(summary[ , "richness"])
mean(summary[ , "abundance"])
mean(summary[ , "simpson"])

# practice calculating for one site
richness <- sum(comm[ "BB", ] > 0)
abundance <- sum(comm[ "BB", ])
simpson <- diversity(comm[ "BB", ], index = "simpson")


# loop through each site

for (x in sites){
  summary[x, "richness"] <- sum(comm[ x, ] > 0)
  summary[x, "abundance"] <- sum(comm[ x, ])
  summary[x, "simpson"] <- diversity(comm[ x, ], index = "simpson")
  summary[x, "invsimpson"] <- diversity(comm[x, ], index = "invsimpson")
  summary[x, "shannon"] <- diversity(comm[ x, ], index = "shannon")
  
}

summary

summary <- data.frame(summary)
summary$sites <- sites
summary <- summary %>% arrange(sites) %>% select(sites, richness:shannon)

# export data to to CommSummary file
write.csv(summary, file = "CommSummary15_16.csv", quote = F)







### FIGURE 2 in Mansucript
####################################
###  repeat creation of community matrix, but with genera rather than species

# make an empty matrix (14 rows x 80 columns)
genera_comm <- matrix(0, nrow = length(sites), ncol = length(ugenera))
dim(genera_comm)


# add a rownames to comat matrix

rownames(genera_comm) <- sites
colnames(genera_comm) <- ugenera
genera_comm

# loop through site and count # in each genera

for (yy in sites){
  
  gencounts <- table(data$Genus[ data$Site == yy ])
  genera_comm[ as.character(yy), names(gencounts) ] <- gencounts
  
}

genera_comm
sum(genera_comm) # should equal 4737

# export data to Community Matrix file
write.csv(genera_comm, file = "Genera_CommunityMatrix.csv", quote = F)


##### Make a stacked bar plot with relative abundance of each genera
# group genera with less than 15 individuals into "other bee" and "other wasp" to reduce clutter


# find the maximum number of occurences for each genus 
max_num <- rep(NA, length(ugenera))
names(max_num) <- ugenera

for (i in ugenera){
  temp <- genera_comm[ ,i]
  max_num[i] <- max(as.numeric(temp))
}
max_num

# sort so that highest num ind is first
max_num <- sort(max_num, decreasing = T)

# find genera that have >20 ind for at least one site
x <- max_num[max_num > 20]
main <- names(x)

# find all other genera with < 20 individuals at any site (rarer genera)
y <- max_num[max_num < 20]
others <- names(y)

# make a data set with the genera with >20 
updated_genera <- genera_comm[ , main]


# sum all other genera and create a column in updated_genera matrix with sum named "Other"
Other_Genera <- matrix(NA, nrow = length(sites), ncol = 1)
rownames(Other_Genera) <- sites
colnames(Other_Genera) <- "Other Genera"

for (ii in sites){
  Other_Genera[ii, ] <- sum(genera_comm[ ii, others])
}

# put the two data frames together to get the main genera with > 20 and a group of all other genera combined
updated_genera <- as.data.frame(cbind(updated_genera, Other_Genera))

total_abundance <- rowSums(updated_genera)
totals <- as.data.frame(cbind(sites, total_abundance))

######## Stacked bar plot with relative abundance of each genera at each site
library(tidyr)


# colors to be used in the stacked bar plot
#Update order of colors to match colors with appropriate host species
Tol_muted2 <- c('#882255','#DDCC77', '#44AA99', '#117733','#88CCEE', '#332288', '#DDDDDD', '#AA4499', '#CC6677', '#999933')


# add sites to data frame
updated_genera$sites <- sites


# pivot data longer to get a column for genus and column with relative abundance
updated_genera <- updated_genera %>%
  pivot_longer(cols = -sites, names_to = "genus", values_to = "abundance")


# add a coulmn with the total abundance and a height that will be used for species richness labels
sum_genera <- left_join(updated_genera, totals, by= "sites")
sum_genera$total_abundance <- as.numeric(sum_genera$total_abundance)
sum_genera <- sum_genera %>% 
  arrange(total_abundance) # order by lowest to highest total abundance

# convert genus to a factor
sum_genera$genus <- as.factor(sum_genera$genus)
str(sum_genera)

# species richness for each site (to be used as labels in the plot)
richness <- summary %>% 
  arrange(abundance) %>%
  select(sites, richness)

sum_genera <- left_join(sum_genera, richness, by = "sites")
View(sum_genera)



# ggplot stacked bar plot (Figure 2)
p <- ggplot(sum_genera, aes(x = reorder(sites, abundance), y = abundance, fill = genus)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(name = "Genus", values = Tol_muted2,
                    breaks = c("Apis", "Bombus", "Lasioglossum", "Peponapis", "Augochlora", "Vespula", "Other Genera"),
                    labels = c(bquote(italic("Apis")), bquote(italic("Bombus")), bquote(italic("Lasioglossum")), bquote(italic("Eucera")),
                               bquote(italic("Augochlora")), bquote(italic("Vespula")), bquote("Other Genera"))) +
  labs(x = "Sites", y = "Pollinator abundance") +
  geom_text(aes(y = total_abundance + 20, label = richness), size = 10/.pt, color = "black") +
  ylim(0,800) +
  theme_classic() +
  theme(axis.text = element_text(size=8, color = "black"), axis.title = element_text(size=10, color="black"),
        legend.title =  element_text(size=10, color="black"), legend.text = element_text(size=8, color = "black"))
print(p)
ggsave("Fig2.tiff", plot = p, dpi = 600, width = 6, height = 6, units = "in", compression="lzw")


# color blind check
cvdPlot(p)






#### Estimate the asymptotes of the species rarefaction curves
####################################
### Appendix S1 Figure S2 and remainder of Table S13

# load libraries
library(iNEXT)


# Making a list of species abundances found at each site (must be a list to work)
sites
BB <- comm["BB",(comm["BB",] > 0)]
BC <- comm["BC",(comm["BC",] > 0)]
BJ <- comm["BJ",(comm["BJ",] > 0)]
BP <- comm["BP",(comm["BP",] > 0)]
E <- comm["E",(comm["E",] > 0)]
G <- comm["G",(comm["G",] > 0)]
GT <- comm["GT",(comm["GT",] > 0)]
K <- comm["K",(comm["K",] > 0)]
PL <- comm["PL",(comm["PL",] > 0)]
PR <- comm["PR",(comm["PR",] > 0)]
S <- comm["S",(comm["S",] > 0)]
SP <- comm["SP",(comm["SP",] > 0)]
TA <- comm["T",(comm["T",] > 0)]
W <- comm["W",(comm["W",] > 0)]


# List with species from each site divided up
bees <- list("BB" = BB, "BC" = BC,"BJ" = BJ, "BP" = BP, "E" = E, "G" = G, "GT" = GT, "K" = K,
             "PL" = PL, "PR" = PR, "S" = S, "SP" = SP, "T" = TA, "W" = W)

mean(summary[ , "abundance"])

#rarefaction at asymptote
xx <- iNEXT(bees, q = 0, datatype = "abundance")
xx$AsyEst

#rarefaction at minimum abundance collected (46 individuals)
xx2 <- iNEXT(bees, q = 0, datatype = "abundance", endpoint = 46)
xx2$AsyEst  # Species richness data used in Table S13 (middle 3 columns)

#rarefaction at second minimum abundance collected
xx3 <- iNEXT(bees, q = 0, datatype = "abundance", endpoint = 183)
xx3$AsyEst

#rarefaction at average abundance collected
xx4 <- iNEXT(bees, q = 0, datatype = "abundance", endpoint = 338)
xx4$AsyEst

#estimator of total bee richness based on asymptote
BeeRichness <- ChaoRichness(bees, datatype = "abundance", conf = 0.95)


# Appendix S1 TABLE S13 in Supplementary Material (last 4 columns)
####################################
# export data to CSV
write.csv(BeeRichness, file = "Estimated Bee Richness.csv", quote = F)





# Appendix S1 FIGURE S2 in Supplementary Material
###############

# added colors to the color blind palette to get to 14
Tol_muted3 <- c('#882255','#3399FF', '#117733', '#332288', '#DDCC77','#CC6677', '#88CCEE', '#AA4499', '#999999', "#000000", '#EE8866', "#E69F00", "#009E73", '#999933')

# shapes to use in the figure
shapes <- c(1:14)


rarefaction <- ggiNEXT(xx, type = 1, se = FALSE, grey = FALSE) +
  scale_color_manual(values = Tol_muted3, name = "Sites") +
  scale_shape_manual(values = shapes, name = "Sites") +
  theme(axis.text = element_text(size=10, color = "black"), axis.title = element_text(size=12, color="black"), 
        legend.text = element_text(size = 10, color = "black"), legend.position = "bottom")
print(rarefaction)
ggsave("FigS2.tiff", plot = rarefaction, dpi = 600, width = 8, height = 5, units = "in", compression="lzw")

# color blind check
cvdPlot(rarefaction)





###########################################
# TEST OF COMMUNITY NESTEDNESS
###########################################

#convert community matrix to a data frame
comm <- as.data.frame(comm)


# test of nestedness temperature against a null community model
oecosimu(comm, nestedtemp, method = "r00")



x <- nestedtemp(comm)




###### Appendix S1 Figure S1
tiff("FigS1.tiff", units = "in", res = 600, height = 5, width = 9, compression = "lzw")
#windows(height = 10, width = 20)
par(oma = c(1,1,1,1))
par(mar = c(0.5,2.5,2.5,0.5))
#plot.new()
plot(x, kind = "incidence", col= c("white", "black"),  names = F)
mtext(text = "Species", side = 3, line = 2.2, cex = 0.85)
mtext(text = "Sites", side = 2, line = 2, cex = 0.85)
axis(side = 3, at = seq(0, 1, len = ncol(x$u)), labels = colnames(x$u), cex.axis = 0.4, las = 2)
axis(side = 2, at = seq(1, 0, len = nrow(x$u)), labels = rownames(x$u), cex.axis = 0.75)
dev.off()



# additional nestedness plot with heat map colors
plot(x, kind = "temperature", col=rev(heat.colors(100)),  names = T)






#################################
### NMDS FIGURES AND ANALYSES FOR COMMUNITY COMPOSITION VS VIRUS PREVALENCE
#################################

# add sites as a column to use for reordering dataframe
comm2 <- as.data.frame(comm)
comm2$sites <- rownames(comm2)



#re-order sites in alphabetical order
comm2 <- arrange(comm2, sites)
rownames(comm2) <- comm2$sites # add back the row names
comm2 <- select(comm2, -sites) #remove sites column for NMDS

# read in file with virus prevalence for each species and virus
viruses <- read.csv("VirusPrevbySpp+Site.csv", sep = ",", stringsAsFactors = F)

# fix NA's for site K Pepo prevalence (no Pepo recorded at K, therefore no prevalence) Changed to zeros to allow models to work
viruses[ 8, "PepoDWV"] <- 0
viruses[ 8, "PepoBQCV"] <- 0
viruses[ 8, "PepoSBV"] <- 0

# NMDS plot with only sites 
res <- metaMDS(comm2, distance = "bray", k=2)
stressplot(res)
res$stress


## Appendix S1 FIGURE S3 in Supplementary Material
# NMDS of SITES only 
tiff("FigS3.tiff", units = "in", res = 600, height = 5, width = 5, compression = "lzw")
#windows(height = 7, width = 7)
#plot.new()
par(oma = c(1,1,1,1))
par(mar = c(4,4,0.5,0.5))
ordiplot(res, type = "none", display = "sites")
orditorp(res, display = "sites", cex = 0.95)
text(x = 0.6, y = -0.9, labels = "Stress = 0.13", font = 2, col = "black", cex = 1.2)
dev.off()

# extract x and y coordinates from MDS plot in a new data frame
MDS_xy <- data.frame(res$points)

# figure label function for base R figures
    # from https://www.r-bloggers.com/2017/03/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/
fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}




### Appendix S1 FIGURE S4 in Manuscript
### NMDS ANALYSES FOR COMMUNITY COMPOSITION VS VIRUS PREVALENCE
# Plot of APIS viruses and NMDS
tiff("FigS4_AB.tiff", units = "in", res = 600, height = 6.5, width = 9, compression = "lzw")
#windows(height = 8, width = 12)
xx <- c(1:6)
lmat <- matrix( xx, nrow = 2, ncol = 3, byrow = T)
layout(lmat)   # function that tells R how to divide up the plot window
layout.show(6)
par(oma = c(1,2,1,1))
par(mar = c(3,3,3,1))
##### Appendix S1 FIGURE s4 A
### NMDS plots of Apis with DWV, BQCV, and SBV

# Plot with sites and overlayed with contour lines for apis DWV prevalence
#plot.new()
#ordiplot(res, type = "none", display = "sites")

ordisurf(res, viruses$ApisDWV, col= "gray34", plot = T, display = "sites", labcex = 0.8, cex = 0.2, nlevels = 6, main = "Apis DWV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9)
A_Dwv <- envfit(res ~ ApisDWV, data=viruses)
fig_label("a", region = "device", pos = "topleft", cex = 2)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)


# comm stats for Apis DWV prevalence
ordi_ApisDWV <- ordisurf(res,viruses$ApisDWV,plot = FALSE)
summary(ordi_ApisDWV)
text(x = 0.6, y = 0.9, labels = "Stress = 0.13", font = 2, col = "black", cex = 1.1)
text(x = 0.6, y = -0.65, labels = "F=0.55, p=0.12", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "= 0.28")), col = "black", cex = 1.1)


# Plot with sites and overlayed with contour lines for apis BQCV prevalence
#windows(height = 7, width = 7)
#plot.new()
#ordiplot(res, type = "n", display = "sites")
ordisurf(res, viruses$ApisBQCV, col= "red", display = "sites", labcex = 0.8, cex = 0.2, nlevels = 6, main = "Apis BQCV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9, air = 1)
A_Bqcv <- envfit(res ~ ApisBQCV, data=viruses)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# comm stats for Apis BQCV prevalence
ordi_ApisBQCV <- ordisurf(res,viruses$ApisBQCV,plot = FALSE)
summary(ordi_ApisBQCV)
text(x = 0.6, y = -0.65, labels = "F=0.07, p=0.3", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "= 0.04")), col = "black", cex = 1.1)

# Plot with sites and overlayed with contour lines for apis SBV prevalence
#windows(height = 7, width = 7)
#plot.new()
#ordiplot(res, type = "n", display = "sites")
ordisurf(res, viruses$ApisSBV, col= "blue", display = "sites", labcex = 0.8, cex = 0.2, nlevels = 6, main = "Apis SBV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9)
A_Sbv <- envfit(res ~ ApisSBV, data=viruses)
plot(A_Sbv, col = "black")
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# comm stats for Apis SBV prevalence
ordi_ApisSBV <- ordisurf(res, viruses$ApisSBV,plot = FALSE)
summary(ordi_ApisSBV)
text(x = 0.6, y = -0.65, labels = "F=2.98, p=0.03", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "= 0.68")), col = "black", cex = 1.1)


##### Appendix S1 FIGURE S4 B
### NMDS plots of Bombus with DWV, BQCV, and SBV

# Plot with sites and overlayed with contour lines for bombus DWV prevalence
ordisurf(res, viruses$BombusDWV, col= "gray34", plot = T, display = "sites", labcex = 0.8, cex = 0.2, nlevels = 6, main = "Bombus DWV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9)
B_Dwv <- envfit(res ~ BombusDWV, data=viruses)
fig_label("b", region = "device", pos = "left", cex = 2)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# modified comm stats for Bombus DWV prevalence
ordi_BombusDWV <- ordisurf(res,viruses$BombusDWV,plot = FALSE)
summary(ordi_BombusDWV)
text(x = 0.6, y = -0.65, labels = "F<0.01, p=0.8", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "<0.01")), col = "black", cex = 1.1)



# Plot with sites and overlayed with contour lines for bombus BQCV prevalence
ordisurf(res, viruses$BombusBQCV, col= "red", display = "sites", labcex = 0.8, cex = 0.2, nlevels = 5, main = "Bombus BQCV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9, air = 1)
B_Bqcv <- envfit(res ~ BombusBQCV, data=viruses)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# modified comm stats for Bombus BQCV prevalence
ordi_BombusBQCV <- ordisurf(res,viruses$BombusBQCV,plot = FALSE)
summary(ordi_BombusBQCV)
text(x = 0.6, y = -0.65, labels = "F=0.09, p=0.33", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "= 0.06")), col = "black", cex = 1.1)


# Plot with sites and overlayed with contour lines for Bombus SBV prevalence
ordisurf(res, viruses$BombusSBV, col= "blue", display = "sites", labcex = 0.8, cex = 0.2, main = "Bombus SBV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9)
B_Sbv <- envfit(res ~ BombusSBV, data=viruses)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# modified comm stats for Bombus SBV prevalence
ordi_BombusSBV <- ordisurf(res,viruses$BombusSBV,plot = FALSE)
summary(ordi_BombusSBV)
text(x = 0.6, y = -0.65, labels = "F=0.04, p=0.34", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "= 0.03")), col = "black", cex = 1.1)
dev.off()




##### Appendix S1 FIGURE s4 C
### NMDS plots of Lasi with DWV, BQCV, and SBV
tiff("FigS4_CD.tiff", units = "in", res = 600, height = 6.5, width = 9, compression = "lzw")
#windows(height = 8, width = 12)
xx <- c(1:6)
lmat <- matrix( xx, nrow = 2, ncol = 3, byrow = T)
layout(lmat)   # function that tells R how to divide up the plot window
layout.show(6)
par(oma = c(1,2,1,1))
par(mar = c(3,3,3,1))

# Plot with sites and overlayed with contour lines for Lasi DWV prevalence
ordisurf(res, viruses$LasiDWV, col= "gray34", plot = T, display = "sites", labcex = 0.8, cex = 0.2, nlevels = 5, main = "Lasioglossum DWV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9)
L_Dwv <- envfit(res ~ LasiDWV, data=viruses)
fig_label("c", region = "device", pos = "topleft", cex = 2)
plot(L_Dwv, col = "black", cex = 0.9)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# modified comm stats for Lasi DWV prevalence
ordi_LasiDWV <- ordisurf(res,viruses$LasiDWV,plot = FALSE)
summary(ordi_LasiDWV)
text(x = 0.6, y = -0.65, labels = "F=2.17, p=0.02", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "= 0.60")), col = "black", cex = 1.1)


# Plot with sites and overlayed with contour lines for Lasi BQCV prevalence
ordisurf(res, viruses$LasiBQCV, col= "red", display = "sites", labcex = 0.8, cex = 0.2, nlevels = 6, main = "Lasioglossum BQCV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9, air = 1)
L_Bqcv <- envfit(res ~ LasiBQCV, data=viruses)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# modified comm stats for Lasi BQCV prevalence
ordi_LasiBQCV <- ordisurf(res,viruses$LasiBQCV,plot = FALSE)
summary(ordi_LasiBQCV)
text(x = 0.6, y = -0.65, labels = "F=0.38, p=0.11", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "= 0.21")), col = "black", cex = 1.1)


# Plot with sites and overlayed with contour lines for Lasi SBV prevalence
ordisurf(res, viruses$LasiSBV, col= "blue", display = "sites", labcex = 0.8, cex = 0.2, main = "Lasioglossum SBV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9)
L_Sbv <- envfit(res ~ LasiSBV, data=viruses)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# modified comm stats for Lasi SBV prevalence
ordi_LasiSBV <- ordisurf(res,viruses$LasiSBV,plot = FALSE)
summary(ordi_LasiSBV)
text(x = 0.65, y = -0.65, labels = "F<0.01, p=0.45", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "< 0.01")), col = "black", cex = 1.1)


##### Appendix S1 FIGURE s4 C
### NMDS plots of Pepo with DWV, BQCV, and SBV


# Plot with sites and overlayed with contour lines for Pepo DWV prevalence
ordisurf(res, viruses$PepoDWV, col= "gray34", plot = T, display = "sites", labcex = 0.8, cex = 0.2, nlevels = 5, main = "Eucera DWV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9)
P_Dwv <- envfit(res ~ PepoDWV, data=viruses)
fig_label("d", region = "device", pos = "left", cex = 2)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# modified comm stats for Pepo DWV prevalence
ordi_PepoDWV <- ordisurf(res,viruses$PepoDWV,plot = FALSE)
summary(ordi_PepoDWV)
text(x = 0.6, y = -0.65, labels = "F=0.28, p=0.15", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "= 0.16")), col = "black", cex = 1.1)


# Plot with sites and overlayed with contour lines for Pepo BQCV prevalence
ordisurf(res, viruses$PepoBQCV, col= "red", display = "sites", labcex = 0.8, cex = 0.2, nlevels = 5, main = "Eucera BQCV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9, air = 1)
P_Bqcv <- envfit(res ~ PepoBQCV, data=viruses)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# modified comm stats for Pepo BQCV prevalence
ordi_PepoBQCV <- ordisurf(res,viruses$PepoBQCV,plot = FALSE)
summary(ordi_PepoBQCV)
text(x = 0.6, y = -0.65, labels = "F<0.01, p=0.79", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "< 0.01")), col = "black", cex = 1.1)


# Plot with sites and overlayed with contour lines for Pepo SBV prevalence
ordisurf(res, viruses$PepoSBV, col= "blue", display = "sites", labcex = 0.8, cex = 0.2, nlevels = 6, main = "Eucera SBV Prevalence", cex.axis = 1.3, ann = F)
orditorp(res, display = "sites", cex = 0.9)
P_Sbv <- envfit(res ~ PepoSBV, data=viruses)
mtext("NMDS1", side = 1, cex = 0.9, line = 2.5)
mtext("NMDS2", side = 2, cex = 0.9, line = 2.5)

# modified comm stats for Pepo SBV prevalence
ordi_PepoSBV <- ordisurf(res,viruses$PepoSBV,plot = FALSE)
summary(ordi_PepoSBV)
text(x = 0.6, y = -0.65, labels = "F=0.006, p=0.41", col = "black", cex = 1.1)
text(x = 0.6, y = -0.78, labels = expression(paste("Adj.", R^2, "= 0.004")), col = "black", cex = 1.1)
dev.off()




