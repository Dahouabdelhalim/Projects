####################################################################################################################################
####################################################################################################################################
# Charles D. Waters
# NOAA Alaska Fisheries Science Center
# Auke Bay Laboratories
# 17109 Point Lena Loop Road
# Juneau, AK 99801
# Email: Charlie.Waters@noaa.gov
####################################################################################################################################
####################################################################################################################################

# Code to quantify the effect of inbreeding coefficient, F, from PLINK 
# on each of eight fitness-related traits using linear mixed models 

####################################################################################################################################
####################################################################################################################################

# Load packages that will be needed
install.packages("ggplot2", dependencies = T)
install.packages("plyr", dependencies = T)
install.packages("GGally", dependencies = T)
install.packages("car", dependencies = T)
install.packages("MASS", dependencies = T)
install.packages("stargazer", dependencies = T)
install.packages("pwr", dependencies = T)
install.packages("gridExtra", dependencies = T)
install.packages("RcmdrMisc", dependencies = T)
install.packages("openxlsx", dependencies = T)
library(ggplot2)
library(plyr)
library(GGally)
library(car)
library(MASS)
library(stargazer)
library(pwr)
library(gridExtra)
library(RcmdrMisc)
library(openxlsx)

# Read in empirical estimates of F obtained from PLINK. The estimates were obtained using allele freqs
# from the P1 Founders, 1100 unlinked loci, and all 452 study individuals. 

F_Plink_1100 <- read.csv("Data file for effects of F on phenotypes.csv", header=T)

cor(F_Plink_1100$Plink_het_1100, F_Plink_1100$Plink_het_6805)
plot(F_Plink_1100$Plink_het_1100, F_Plink_1100$Plink_het_6805)

# Make sure that some covariates are factors
F_Plink_1100$Line <- factor(F_Plink_1100$Line)
F_Plink_1100$Sex <- factor(F_Plink_1100$Sex) 
F_Plink_1100$Gen <- factor(F_Plink_1100$Gen) 
F_Plink_1100$Age <- factor(F_Plink_1100$Age)

# Add a column to the data containing F^2, to test for non-linear (i.e. quadratic) effect of F on fitness
F_Plink_1100$Fsquared <- F_Plink_1100$Plink_het_1100^2

# Look at distributions of samples within levels of each covariate
table(F_Plink_1100$Age) # Age 4 fish comprise 90% of all of the samples in this data set (for which age is known)
table(F_Plink_1100$Sex)
table(F_Plink_1100$Gen)
table(F_Plink_1100$Line) 

table(F_Plink_1100$Line, F_Plink_1100$Age)
table(F_Plink_1100$Line, F_Plink_1100$Gen) # P1 Founders are only present in 1998, not the subsequent generations 

# Given that there is an unequal distribution of ages by sex and year, and because Age 4 fish comprise 90% of the samples,
# we will exclude Age 3 and Age 5 fish from these regression analyses. We will also exclude the P1 Founders to achieve a balanced design,
# since the Founders are only present in 1998 while 1998 also lacks INT and SEG fish (we therefore cannot separate 
# the effect of the P1 founders from the year/generation effect of 1998).

F_Plink_1100 <- F_Plink_1100[which(F_Plink_1100$Age==4 & F_Plink_1100$Line!="P1"),]
F_Plink_1100$Gen <- factor(F_Plink_1100$Gen, levels=c("F1", "F2", "F3", "F4"))


####################################################################################################################################
####################################################################################################################################

# Then visually check for correlations among covariates as recommended by Zuur et al. 2010

# Select covariates to explore
covariates <- F_Plink_1100[,c("Line","Gen","Sex", "Plink_het_1100")]

# Now visually examine for collinearity amongst our covariates, as recommended by Zuur et al. 2010
ggpairs(covariates)

# The resulting plot shows some imbalances within predictors
# but there does not appear to be any large correlations between predictors, which is good. 

####################################################################################################################################
####################################################################################################################################

# Now standardize (mean center and scale by standard deviation) the numeric variables (both response and explanatory) 
# so that direct comparisons can be made between them. I'm going to standardize each one individually so that I don't get confused.

F_Plink_1100$scaled_fecundity <- scale(F_Plink_1100$bias_corrected_fecundity, center=T, scale=T) #I confirmed manually that this function correctly standardizes the values
F_Plink_1100$scaled_repro_effort <- scale(F_Plink_1100$repro_effort, center=T, scale=T) 
F_Plink_1100$scaled_roza_day <- scale(F_Plink_1100$roza_day, center=T, scale=T) 
F_Plink_1100$scaled_spawn_day <- scale(F_Plink_1100$spawn_day, center=T, scale=T) 
F_Plink_1100$scaled_length <- scale(F_Plink_1100$forklength, center=T, scale=T) 
F_Plink_1100$scaled_weight <- scale(F_Plink_1100$weight, center=T, scale=T) 
F_Plink_1100$scaled_condition <- scale(F_Plink_1100$Condition_factor, center=T, scale=T) 
F_Plink_1100$scaled_DGC <- scale(F_Plink_1100$DGC, center=T, scale=T) 
F_Plink_1100$scaled_F <- scale(F_Plink_1100$Plink_het_1100, center=T, scale=T) 
F_Plink_1100$scaled_F_squared <- scale(F_Plink_1100$Fsquared)

####################################################################################################################################
####################################################################################################################################

# We are now ready to plot level of inbreeding vs the eight fitness measures, and quantify the effect of F on the traits using linear models. 
par(mar=c(6,6,4,4))

####################################################################################################################################
####################################################################################################################################

# Fecundity, after correcting for residual ovarian fluid that remained in the egg mass (Knudsen et al. 2008)

# First examine the phenotypic distribution
hist(F_Plink_1100$scaled_fecundity) # normally distributed, for the most part

# Now plot the trait vs inbreeding coefficient to visually assess effects
plot(F_Plink_1100$Plink_het_1100, F_Plink_1100$scaled_fecundity,xlab="Inbreeding Coefficient",ylab="Fecundity", cex.lab=2, cex.axis=2,pch=19,cex=1.5)

# Construct linear model containing all possible covariates (i.e. full model)
IBD_fecundity.full_scaled <- lm(scaled_fecundity~scaled_F + scaled_F_squared + Gen + Line + scaled_F:Gen + scaled_F:Line, data=F_Plink_1100) #Omit sex as a covariate since this trait is for females only
summary(IBD_fecundity.full_scaled)
vif(IBD_fecundity.full_scaled) #VIF less than 3 for all terms (look at GVIF^(1/2*df))

# Now perform backward model selection to identify an optimal model (based on AIC) using the stepAIC function from the MASS package
# However, we want to retain the main effect of inbreeding coefficient regardless of its significance
# because we are interested in its specific effect.

IBD_fecundity_step_scaled <- stepAIC(IBD_fecundity.full_scaled, scope=list(upper=IBD_fecundity.full_scaled, lower=~scaled_F), direction="backward")
IBD_fecundity_step_scaled$anova  # The optimal model includes F, Gen, and Line

# Now run the optimal model
IBD_fecundity.final_scaled <- lm(scaled_fecundity~scaled_F + Gen + Line, data=F_Plink_1100)
summary(IBD_fecundity.final_scaled) 
vif(IBD_fecundity.final_scaled) # VIF close to 1 for all variables, which is good

####################################################################################################################################
####################################################################################################################################

# Reproductive effort

# First examine the phenotypic distribution
hist(F_Plink_1100$scaled_repro_effort) # somewhat skewed normal distribution

# Now plot the trait vs inbreeding coefficient
plot(F_Plink_1100$Plink_het_1100, F_Plink_1100$scaled_repro_effort,xlab="Inbreeding Coefficient",ylab="Reproductive Effort", cex.lab=2, cex.axis=2,pch=19,cex=1.5)

# Construct linear model containing all possible covariates (i.e. full model)
IBD_repro.full_scaled <- lm(scaled_repro_effort~scaled_F + scaled_F_squared + Gen + Line + scaled_F:Gen + scaled_F:Line, data=F_Plink_1100) #Omit sex as a covariate since this trait is for females only
summary(IBD_repro.full_scaled) 
vif(IBD_repro.full_scaled) #VIF less than 3 for all terms, which is ok (look at GVIF^(1/2*df))

# Now perform backward model selection to identify an optimal model (based on AIC) using the stepAIC function from the MASS package
# However, we want to retain the main effect of inbreeding coefficient regardless of its significance
# because we are interested in its specific effect.

IBD_repro_step_scaled <- stepAIC(IBD_repro.full_scaled, scope=list(upper=IBD_repro.full_scaled, lower=~scaled_F), direction="backward")
IBD_repro_step_scaled$anova  # The optimal model includes F and Gen

# Now run the optimal model
IBD_repro.final_scaled <- lm(scaled_repro_effort~scaled_F + Gen, data=F_Plink_1100)
summary(IBD_repro.final_scaled) 
vif(IBD_repro.final_scaled) # VIF close to 1 for all variables, which is good

####################################################################################################################################
####################################################################################################################################

# Return timing

# First examine the phenotypic distribution
hist(F_Plink_1100$scaled_roza_day) # somewhat skewed normal distribution

# Now plot the trait vs inbreeding coefficient
plot(F_Plink_1100$Plink_het_1100, F_Plink_1100$scaled_roza_day,xlab="Inbreeding Coefficient",ylab="Return Timing", cex.lab=2, cex.axis=2,pch=19,cex=1.5)

# scaled variables
IBD_return.full_scaled <- lm(scaled_roza_day~scaled_F + scaled_F_squared + Gen + Line + Sex + scaled_F:Gen + scaled_F:Line + scaled_F:Sex, data=F_Plink_1100) 
summary(IBD_return.full_scaled) 
vif(IBD_return.full_scaled) #VIF equal to or less than 3 for all terms (look at GVIF^(1/2*df))

# Now perform backward model selection to identify an optimal model (based on AIC) using the stepAIC function from the MASS package
# However, we want to retain the main effect of inbreeding coefficient regardless of its significance
# because we are interested in its specific effect.

IBD_return_step_scaled <- stepAIC(IBD_return.full_scaled, scope=list(upper=IBD_return.full_scaled, lower=~scaled_F), direction="backward")
IBD_return_step_scaled$anova  # The optimal model includes F, Gen, and Sex

# Now run the optimal model
IBD_return.final_scaled <- lm(scaled_roza_day~scaled_F + Gen + Sex, data=F_Plink_1100)
summary(IBD_return.final_scaled) 
vif(IBD_return.final_scaled) #VIF less than 2, which is still good.

####################################################################################################################################
####################################################################################################################################

# Spawn timing

# First examine the phenotypic distribution
hist(F_Plink_1100$scaled_spawn_day) # somewhat normally distribution

# Now plot the trait vs inbreeding coefficient
plot(F_Plink_1100$Plink_het_1100, F_Plink_1100$scaled_spawn_day,xlab="Inbreeding Coefficient",ylab="Spawn Timing", cex.lab=2, cex.axis=2,pch=19,cex=1.5)

# Construct linear model containing all possible covariates (i.e. full model)
IBD_spawn.full_scaled <- lm(scaled_spawn_day~scaled_F + scaled_F_squared + Gen + Line + Sex + scaled_F:Gen + scaled_F:Line + scaled_F:Sex, data=F_Plink_1100) 
summary(IBD_spawn.full_scaled) 
vif(IBD_spawn.full_scaled) #VIF less than 3.1 for all terms (look at GVIF^(1/2*df))

# Now perform backward model selection to identify an optimal model (based on AIC) using the stepAIC function from the MASS package
# However, we want to retain the main effect of inbreeding coefficient regardless of its significance
# because we are interested in its specific effect.

IBD_spawn_step_scaled <- stepAIC(IBD_spawn.full_scaled, scope=list(upper=IBD_spawn.full_scaled, lower=~scaled_F), direction="backward")
IBD_spawn_step_scaled$anova  # The optimal model includes F, Gen, Line, and F:Line

# Now run the optimal model
IBD_spawn.final_scaled <- lm(scaled_spawn_day~scaled_F + Gen + Line + scaled_F:Line, data=F_Plink_1100)
summary(IBD_spawn.final_scaled) 
vif(IBD_spawn.final_scaled) #VIF less than 2, which is still good.

####################################################################################################################################
####################################################################################################################################

# Fork length at Roza

# First examine the phenotypic distribution
hist(F_Plink_1100$scaled_length) # somewhat normally distribution

# Now plot the trait vs inbreeding coefficient
plot(F_Plink_1100$Plink_het_1100, F_Plink_1100$scaled_length,xlab="Inbreeding Coefficient",ylab="Forklength", cex.lab=2, cex.axis=2,pch=19,cex=1.5)

# Construct linear model containing all possible covariates (i.e. full model)
IBD_length.full_scaled <- lm(scaled_length~scaled_F + scaled_F_squared + Gen + Line + Sex + scaled_F:Gen + scaled_F:Line + scaled_F:Sex, data=F_Plink_1100) 
summary(IBD_length.full_scaled) 
vif(IBD_length.full_scaled) #VIF less than 3.1 for all terms (look at GVIF^(1/2*df))

# Now perform backward model selection to identify an optimal model (based on AIC) using the stepAIC function from the MASS package
# However, we want to retain the main effect of inbreeding coefficient regardless of its significance
# because we are interested in its specific effect.

IBD_length_step_scaled <- stepAIC(IBD_length.full_scaled, scope=list(upper=IBD_length.full_scaled, lower=~scaled_F), direction="backward")
IBD_length_step_scaled$anova  # The optimal model includes F, Gen, Line, and Sex

# Now run the optimal model
IBD_length.final_scaled <- lm(scaled_length~scaled_F + Gen + Line + Sex, data=F_Plink_1100)
summary(IBD_length.final_scaled) 
vif(IBD_length.final_scaled) # VIF close to 1 for all variables, which is good

####################################################################################################################################
####################################################################################################################################

# Weight at Roza

# First examine the phenotypic distribution
hist(F_Plink_1100$scaled_weight) # largely normal distribution

# Now plot the trait vs inbreeding coefficient
plot(F_Plink_1100$Plink_het_1100, F_Plink_1100$scaled_weight,xlab="Inbreeding Coefficient",ylab="Weight", cex.lab=2, cex.axis=2,pch=19,cex=1.5)

# Construct linear model containing all possible covariates (i.e. full model)
IBD_weight.full_scaled <- lm(scaled_weight~scaled_F + scaled_F_squared + Gen + Line + Sex + scaled_F:Gen + scaled_F:Line + scaled_F:Sex, data=F_Plink_1100) 
summary(IBD_weight.full_scaled) 
vif(IBD_weight.full_scaled) #VIF less than 3.1 for all terms (look at GVIF^(1/2*df))

# Now perform backward model selection to identify an optimal model (based on AIC) using the stepAIC function from the MASS package
# However, we want to retain the main effect of inbreeding coefficient regardless of its significance
# because we are interested in its specific effect.

IBD_weight_step_scaled <- stepAIC(IBD_weight.full_scaled, scope=list(upper=IBD_weight.full_scaled, lower=~scaled_F), direction="backward")
IBD_weight_step_scaled$anova  # The optimal model includes F, Gen, Line, and Sex

# Now run the optimal model
IBD_weight.final_scaled <- lm(scaled_weight~scaled_F + Gen + Line + Sex, data=F_Plink_1100)
summary(IBD_weight.final_scaled) 
vif(IBD_weight.final_scaled) # VIF is 2 or lower, which is good.

####################################################################################################################################
####################################################################################################################################

# Condition factor

# First examine the phenotypic distribution
hist(F_Plink_1100$scaled_condition) # mostly normally distributed

# Now plot the trait vs inbreeding coefficient
plot(F_Plink_1100$Plink_het_1100, F_Plink_1100$scaled_condition,xlab="Inbreeding Coefficient",ylab="Condition Factor", cex.lab=2, cex.axis=2,pch=19,cex=1.5)

# Construct linear model containing all possible covariates (i.e. full model)

IBD_condition.full_scaled <- lm(scaled_condition~scaled_F + scaled_F_squared + Gen + Line + Sex + scaled_F:Gen + scaled_F:Line + scaled_F:Sex, data=F_Plink_1100) 
summary(IBD_condition.full_scaled) 
vif(IBD_condition.full_scaled) #VIF less than 3.1 for all terms (look at GVIF^(1/2*df))

# Now perform backward model selection to identify an optimal model (based on AIC) using the stepAIC function from the MASS package
# However, we want to retain the main effect of inbreeding coefficient regardless of its significance
# because we are interested in its specific effect.

IBD_condition_step_scaled <- stepAIC(IBD_condition.full_scaled, scope=list(upper=IBD_condition.full_scaled, lower=~scaled_F), direction="backward")
IBD_condition_step_scaled$anova  # The optimal model includes F, F squared, Gen, and Line

# Now run the optimal model
IBD_condition.final_scaled <- lm(scaled_condition~scaled_F + scaled_F_squared + Gen + Line, data=F_Plink_1100)
summary(IBD_condition.final_scaled) 
vif(IBD_condition.final_scaled) # VIF is 2 or lower, which is good.

####################################################################################################################################
####################################################################################################################################

# Daily growth coefficient

# First examine the phenotypic distribution
hist(F_Plink_1100$scaled_DGC) # somewhat skewed distribution

# Now plot the trait vs inbreeding coefficient
plot(F_Plink_1100$Plink_het_1100, F_Plink_1100$scaled_DGC,xlab="Inbreeding Coefficient",ylab="DGC", cex.lab=2, cex.axis=2,pch=19,cex=1.5)

# Construct linear model containing all possible covariates (i.e. full model)

IBD_DGC.full_scaled <- lm(scaled_DGC~scaled_F + scaled_F_squared + Gen + Line + Sex + scaled_F:Gen + scaled_F:Line + scaled_F:Sex, data=F_Plink_1100) 
summary(IBD_DGC.full_scaled) 
vif(IBD_DGC.full_scaled) #VIF less than 3 for all terms (look at GVIF^(1/2*df))

# Now perform backward model selection to identify an optimal model (based on AIC) using the stepAIC function from the MASS package
# However, we want to retain the main effect of inbreeding coefficient regardless of its significance
# because we are interested in its specific effect.

IBD_DGC_step_scaled <- stepAIC(IBD_DGC.full_scaled, scope=list(upper=IBD_DGC.full_scaled, lower=~scaled_F), direction="backward")
IBD_DGC_step_scaled$anova  # The optimal model F, Gen, Line, Sex, and F:Sex

# Now run the optimal model
IBD_DGC.final_scaled <- lm(scaled_DGC~scaled_F + Gen + Line + Sex + scaled_F:Sex, data=F_Plink_1100)
summary(IBD_DGC.final_scaled) 
vif(IBD_DGC.final_scaled) #VIF's are close to 1, which is good.

####################################################################################################################################
####################################################################################################################################

# Now plot results for the final models, broken into two plots to make it visually easier
install.packages('dotwhisker')
library(dotwhisker)
library(broom)
library(dplyr)


IBD_traits1 <- dwplot(list(IBD_fecundity.final_scaled, IBD_repro.final_scaled, IBD_return.final_scaled, IBD_spawn.final_scaled),
                     show_intercept=FALSE, vline=geom_vline(xintercept=0, colour="grey60", linetype=2, size=1.5),
                     line_args=list(size=1.75), whisker_args = list(size=1.75), dot_args = list(size=3, shape=19), dodge_size = 0.3) +
  theme_classic() + xlab(expression(paste("Coefficient"))) +
  theme(axis.text=element_text(size=16, color="black"),axis.title=element_text(size=16),
        legend.text=element_text(size=12), legend.title=element_text(size=14),
        legend.position = "bottom",plot.margin = unit(c(0.55,0.75,0.25,0.1), "cm")) +
  guides(colour=guide_legend(title="")) +
  scale_colour_manual(values=c("black", "blue", "darkgreen", "darkcyan"), breaks = c("Model 1", "Model 2", "Model 3", "Model 4"), 
                      labels = c("Fecundity", "Repro. Effort", "Return Timing","Spawn Timing")) + 
  scale_y_discrete(breaks=c("scaled_F", "GenF2", "GenF3", "GenF4", "LineSEG", "Sexmale","scaled_F:LineSEG"),
                   labels=c("F", bquote('F'[2]~'Gen'), bquote('F'[3]~'Gen'), bquote('F'[4]~'Gen'), "SEG Line", "Sex", "F:SEG Line"))

IBD_traits2 <- dwplot(list(IBD_condition.final_scaled, IBD_length.final_scaled, IBD_weight.final_scaled, IBD_DGC.final_scaled),
                      show_intercept=FALSE, vline=geom_vline(xintercept=0, colour="grey60", linetype=2, size=1.5),
                      line_args=list(size=1.75), whisker_args = list(size=1.75), dot_args = list(size=3, shape=19), dodge_size = 0.4) +
  theme_classic() + xlab(expression(paste("Coefficient"))) +
  theme(axis.text=element_text(size=16, color="black"),axis.title=element_text(size=16),
        legend.text=element_text(size=12), legend.title=element_text(size=14),
        legend.position = "bottom",plot.margin = unit(c(0.55,0.75,0.25,0.1), "cm")) +
  guides(colour=guide_legend(title="")) +
  scale_x_continuous(name="Coefficient",limits=c(-2, 2), breaks=c(-1.5,-1,-0.5, 0, 0.5, 1, 1.5), labels=c('-1.5', '-1.0', '-0.5','0.0','0.5','1.0','1.5')) +
  scale_colour_manual(values=c("gray14", "coral4", "dodgerblue3", "purple4"), breaks = c("Model 1", "Model 2", "Model 3", "Model 4"), 
                      labels = c("Condition Factor","Fork Length", "Weight","DGC")) + 
  scale_y_discrete(breaks=c("scaled_F", "scaled_F_squared", "GenF2", "GenF3", "GenF4", "LineSEG", "Sexmale","scaled_F:Sexmale"),
                   labels=c("F", bquote('F'^2), bquote('F'[2]~'Gen'), bquote('F'[3]~'Gen'), bquote('F'[4]~'Gen'), "SEG Line", "Sex", "F:Sex"))

## Combine both plots onto a single plot
pdf(file = "IBD_model_coefficients.pdf", width=18, height=11) 
grid.arrange(IBD_traits1, IBD_traits2, nrow=1, ncol=2)
dev.off()

####################################################################################################################################
####################################################################################################################################

# Now calculate post-hoc power for each explanatory variable in the FINAL MODELS (after model selection) using the pwr.f2.test
# function from the 'pwr' package. We set the numerator degrees of freedom (u) = 1 for each coefficient. The denominator degrees
# of freedom (v) is calculated from the sample sizes for the group of fish from which the coefficients were estimated. 
# The effect size, f2, must be positive for the function; these values are the regression coefficients themselves.
# We then use an alpha=0.05 for all calculations and determine the power for each coefficient.

####################################################################################################################################
# Fecundity (females only)

# Inbreeding coefficient, F. The coefficient was estimated from all fish analyzed for this trait (the regression coefficient 
# doesn't change regardless of the reference level used for the intercept in the regression). Thus, the 
# degrees of freedom for the denominator (v) were determined from this sample size (n=204). 

pwr.f2.test(u=1, v=202, f2=0.004, sig.level = 0.05, power = NULL) # pwr=0.15, not sensitive to +/- 5 for the v parameter

# Coefficients for the F2, F3, and F4 generations. These coefficients reflect the difference between each generation and the 
# F1 generation, which is included in the intercept (i.e. reference level). The coefficients were estimated using all fish from
# each generation (females only for this trait but both INT and sEG lines): n=52, n=59, n=70, and n=23 for the F1, F2, F3, and F4 generations,
# respectively, so the degrees of freedom for the denominator (v) were determined from these sample sizes. 

pwr.f2.test(u=1, v=109, f2=0.775, sig.level = 0.05, power = NULL) # F2 gen pwr=1.00

pwr.f2.test(u=1, v=120, f2=0.160, sig.level = 0.05, power = NULL) # F3 gen pwr=0.99

pwr.f2.test(u=1, v=73, f2=0.281, sig.level = 0.05, power = NULL) # F4 gen pwr=0.99

# Hatchery line. This coefficient was estimated by comparing all INT fish analzed for this trait (n=95) to 
# all SEG fish analyzed for this trait (n=109), so the degrees of freedom for the denominator (v) were determined
# from these sample sizes.

pwr.f2.test(u=1, v=202, f2=0.212, sig.level = 0.05, power = NULL) # pwr=1.00

####################################################################################################################################
# Reproductive Effort (females only)

# Inbreeding coefficient, F. The coefficient was estimated from all fish analyzed for this trait (the regression coefficient 
# doesn't change regardless of the reference level used for the intercept in the regression). Thus, the 
# degrees of freedom for the denominator (v) were determined from this sample size (n=201). 

pwr.f2.test(u=1, v=199, f2=0.057, sig.level = 0.05, power = NULL) # pwr=0.92, not sensitive to +/- 5 for the v parameter

# Coefficients for the F2, F3, and F4 generations. These coefficients reflect the difference between each generation and the 
# F1 generation, which is included in the intercept (i.e. reference level). The coefficients were estimated using all fish from
# each generation (females only for this trait but both INT and sEG lines): n=50, n=58, n=70, and n=23 for the F1, F2, F3, and F4 generations,
# respectively, so the degrees of freedom for the denominator (v) were determined from these sample sizes. 

pwr.f2.test(u=1, v=106, f2=0.151, sig.level = 0.05, power = NULL) # F2 gen pwr=0.98

pwr.f2.test(u=1, v=118, f2=0.428, sig.level = 0.05, power = NULL) # F3 gen pwr=1.00

pwr.f2.test(u=1, v=71, f2=0.512, sig.level = 0.05, power = NULL) # F4 gen pwr=1.00

####################################################################################################################################
# Return timing

# Inbreeding coefficient, F. The coefficient was estimated from all fish analyzed for this trait (the regression coefficient 
# doesn't change regardless of the reference level used for the intercept in the regression). Thus, the 
# degrees of freedom for the denominator (v) were determined from this sample size (n=336). 

pwr.f2.test(u=1, v=334, f2=0.006, sig.level = 0.05, power = NULL) # pwr=0.29, not sensitive to +/- 5 for the v parameter

# Coefficients for the F2, F3, and F4 generations. These coefficients reflect the difference between each generation and the 
# F1 generation, which is included in the intercept (i.e. reference level). The coefficients were estimated using all fish from
# each generation (both sexes and both lines): n=86, n=96, n=112, and n=42 for the F1, F2, F3, and F4 generations,
# respectively, so the degrees of freedom for the denominator (v) were determined from these sample sizes. 

pwr.f2.test(u=1, v=180, f2=0.299, sig.level = 0.05, power = NULL) # F2 gen pwr=1.00

pwr.f2.test(u=1, v=196, f2=0.181, sig.level = 0.05, power = NULL) # F3 gen pwr=1.00

pwr.f2.test(u=1, v=126, f2=0.190, sig.level = 0.05, power = NULL) # F4 gen pwr=1.00

# Sex. This coefficient was estimated by comparing all females to all males (across both lines and all generations, n= 
# 207 and 129 respectively), so the degrees of freedom for the denominator (v) were determined from these sample sizes.

pwr.f2.test(u=1, v=334, f2=0.217, sig.level = 0.05, power = NULL) # pwr=1.00

####################################################################################################################################
# Spawn timing

# Inbreeding coefficient, F. The inclusion of the F x Line interaction means that 
# the regression coefficient for F was estimated from all fish in the integrated line ONLY. This is evident because
# the regression coefficient for F changes when you change the reference level from INT to SEG. Thus, the 
# degrees of freedom for the denominator (v) were determined from the number of INT fish in the analysis (n=176). 

pwr.f2.test(u=1, v=174, f2=0.250, sig.level = 0.05, power = NULL) # pwr=1.00, not sensitive to +/- 5 for the v parameter

# Coefficients for the F2, F3, and F4 generations. These coefficients reflect the difference between each generation and the 
# F1 generation, which is included in the intercept (i.e. reference level). The coefficients were estimated using all fish from
# each generation (both sexes and both lines): n=84, n=96, n=112, and n=42 for the F1, F2, F3, and F4 generations,
# respectively, so the degrees of freedom for the denominator (v) were determined from these sample sizes. 

pwr.f2.test(u=1, v=178, f2=0.278, sig.level = 0.05, power = NULL) # F2 gen pwr=1.00

pwr.f2.test(u=1, v=194, f2=0.542, sig.level = 0.05, power = NULL) # F3 gen pwr=1.00

pwr.f2.test(u=1, v=124, f2=0.743, sig.level = 0.05, power = NULL) # F4 gen pwr=1.00

# Hatchery line. This coefficient was estimated by comparing all INT fish analzed for this trait (n=176) to 
# all SEG fish analyzed for this trait (n=158), so the degrees of freedom for the denominator (v) were determined
# from these sample sizes.

pwr.f2.test(u=1, v=332, f2=0.929, sig.level = 0.05, power = NULL) # pwr=1.00

# Interaction of inbreeding coefficient, F, and hatchery line. The inclusion of this term indicates that the effect of F differs
# between the hatchery lines. So, this regression coefficient simply refers to the effect of F in the SEG line, while the
# regression coefficient for F (above) refers to the effect of F in the INT line. This interpretation is evident 
# when you change the reference level of the LINE variable from INT to SEG. Thus, the degrees of freedom for the 
# denominator (v) were determined from the number of SEG fish in the analysis (n=158). 

# NOTE that, to calculate actual effect size (f2), I add the regression coefficient to that of the INT line in order
# to calculate the overall effect of F in the SEG line.

# Effect size for SEG line: 0.250 - 0.297 = -0.047
pwr.f2.test(u=1, v=156, f2=0.047, sig.level = 0.05, power = NULL) # SEG pwr=0.77

h####################################################################################################################################
# Fork length

# Inbreeding coefficient, F. The coefficient was estimated from all fish analyzed for this trait (the regression coefficient 
# doesn't change regardless of the reference level used for the intercept in the regression). Thus, the 
# degrees of freedom for the denominator (v) were determined from this sample size (n=335). 

pwr.f2.test(u=1, v=333, f2=0.064, sig.level = 0.05, power = NULL) # pwr=1, not sensitive to +/- 5 for the v parameter

# Coefficients for the F2, F3, and F4 generations. These coefficients reflect the difference between each generation and the 
# F1 generation, which is included in the intercept (i.e. reference level). The coefficients were estimated using all fish from
# each generation (both sexes and both lines): n=86, n=96, n=111, and n=42 for the F1, F2, F3, and F4 generations,
# respectively, so the degrees of freedom for the denominator (v) were determined from these sample sizes. 

pwr.f2.test(u=1, v=180, f2=0.618, sig.level = 0.05, power = NULL) # F2 gen pwr=1.00

pwr.f2.test(u=1, v=195, f2=0.234, sig.level = 0.05, power = NULL) # F3 gen pwr=1.00

pwr.f2.test(u=1, v=126, f2=0.157, sig.level = 0.05, power = NULL) # F4 gen pwr=0.99

# Hatchery line. This coefficient was estimated by comparing all INT fish analzed for this trait (n=177) to 
# all SEG fish analyzed for this trait (n=158), so the degrees of freedom for the denominator (v) were determined
# from these sample sizes.

pwr.f2.test(u=1, v=333, f2=0.392, sig.level = 0.05, power = NULL) # pwr=1.00

# Sex. This coefficient was estimated by comparing all females to all males (across both lines and all generations, n= 
# 206 and 129 respectively), so the degrees of freedom for the denominator (v) were determined from these sample sizes.

pwr.f2.test(u=1, v=333, f2=0.205, sig.level = 0.05, power = NULL) # pwr=1.00

####################################################################################################################################
# Weight

# Inbreeding coefficient, F. The inclusion of the F x Gen interaction means that 
# the regression coefficient for F was estimated from all fish in the F1 generation line ONLY. This is evident because
# the regression coefficient for F changes when you change the reference level from the F1 to F2, F3, or F4 generations. Thus, the 
# degrees of freedom for the denominator (v) were determined from the number of F1 fish in the analysis (n=86). 

pwr.f2.test(u=1, v=84, f2=0.077, sig.level = 0.05, power = NULL) # pwr=0.72, not sensitive to +/- 5 for the v parameter

# Coefficients for the F2, F3, and F4 generations. These coefficients reflect the difference between each generation and the 
# F1 generation, which is included in the intercept (i.e. reference level). The coefficients were estimated using all fish from
# each generation (both sexes and both lines): n=86, n=96, n=112, and n=42 for the F1, F2, F3, and F4 generations,
# respectively, so the degrees of freedom for the denominator (v) were determined from these sample sizes. 

pwr.f2.test(u=1, v=180, f2=0.743, sig.level = 0.05, power = NULL) # F2 gen pwr=1.00

pwr.f2.test(u=1, v=196, f2=0.164, sig.level = 0.05, power = NULL) # F3 gen pwr=1.00

pwr.f2.test(u=1, v=126, f2=0.855, sig.level = 0.05, power = NULL) # F4 gen pwr=1.00

# Hatchery line. This coefficient was estimated by comparing all INT fish analzed for this trait (n=177) to 
# all SEG fish analyzed for this trait (n=159), so the degrees of freedom for the denominator (v) were determined
# from these sample sizes.

pwr.f2.test(u=1, v=334, f2=0.437, sig.level = 0.05, power = NULL) # pwr=1.00

# Sex. This coefficient was estimated by comparing all females to all males (across both lines and all generations, n= 
# 207 and 129 respectively), so the degrees of freedom for the denominator (v) were determined from these sample sizes.

pwr.f2.test(u=1, v=334, f2=0.207, sig.level = 0.05, power = NULL) # pwr=1.00

####################################################################################################################################
# Condition Factor

# Inbreeding coefficient, F. The inclusion of the F x Gen interaction means that 
# the regression coefficient for F was estimated from all fish in the F1 generation line ONLY. This is evident because
# the regression coefficient for F changes when you change the reference level from the F1 to F2, F3, or F4 generations. Thus, the 
# degrees of freedom for the denominator (v) were determined from the number of F1 fish in the analysis (n=86). 

pwr.f2.test(u=1, v=84, f2=0.053, sig.level = 0.05, power = NULL) # pwr=0.98, not sensitive to +/- 5 for the v parameter

# Inbreeding coefficient squared (F^2). The coefficient was also estimated from all fish analyzed for this trait (the regression coefficient 
# doesn't change regardless of the reference level used for the intercept in the regression). Thus, the 
# degrees of freedom for the denominator (v) were determined from this sample size (n=335). 

pwr.f2.test(u=1, v=333, f2=0.094, sig.level = 0.05, power = NULL) # pwr=1.00

# Coefficients for the F2, F3, and F4 generations. These coefficients reflect the difference between each generation and the 
# F1 generation, which is included in the intercept (i.e. reference level). The coefficients were estimated using all fish from
# each generation (both sexes and both lines): n=86, n=96, n=111, and n=42 for the F1, F2, F3, and F4 generations,
# respectively, so the degrees of freedom for the denominator (v) were determined from these sample sizes. 

pwr.f2.test(u=1, v=180, f2=0.626, sig.level = 0.05, power = NULL) # F2 gen pwr=1.00

pwr.f2.test(u=1, v=195, f2=0.061, sig.level = 0.05, power = NULL) # F3 gen pwr=0.93

pwr.f2.test(u=1, v=126, f2=1.672, sig.level = 0.05, power = NULL) # F4 gen pwr=1.00

# Hatchery line. This coefficient was estimated by comparing all INT fish analzed for this trait (n=177) to 
# all SEG fish analyzed for this trait (n=158), so the degrees of freedom for the denominator (v) were determined
# from these sample sizes.

pwr.f2.test(u=1, v=333, f2=0.185, sig.level = 0.05, power = NULL) # pwr=1.00

####################################################################################################################################
# Daily growth coefficient

# Inbreeding coefficient, F. The inclusion of the F x Sex interaction means that 
# the regression coefficient for F was estimated from all female fish ONLY (reference level). This is evident because
# the regression coefficient for F changes when you change the reference level from female to male. Thus, the 
# degrees of freedom for the denominator (v) were determined from the number of female fish in the analysis (n=201). 

pwr.f2.test(u=1, v=199, f2=0.013, sig.level = 0.05, power = NULL) # pwr=0.36, not sensitive to +/- 5 for the v parameter

# Coefficients for the F2, F3, and F4 generations. These coefficients reflect the difference between each generation and the 
# F1 generation, which is included in the intercept (i.e. reference level). The coefficients were estimated using all fish from
# each generation (both sexes and both lines): n=84, n=94, n=108, and n=42 for the F1, F2, F3, and F4 generations,
# respectively, so the degrees of freedom for the denominator (v) were determined from these sample sizes. 

pwr.f2.test(u=1, v=176, f2=0.263, sig.level = 0.05, power = NULL) # F2 gen pwr=1.00

pwr.f2.test(u=1, v=190, f2=0.690, sig.level = 0.05, power = NULL) # F3 gen pwr=1.00

pwr.f2.test(u=1, v=124, f2=1.554, sig.level = 0.05, power = NULL) # F4 gen pwr=1.00

# Hatchery line. This coefficient was estimated by comparing all INT fish analzed for this trait (n=175) to 
# all SEG fish analyzed for this trait (n=153), so the degrees of freedom for the denominator (v) were determined
# from these sample sizes.

pwr.f2.test(u=1, v=326, f2=0.212, sig.level = 0.05, power = NULL) # pwr=1.00

# Sex. This coefficient was estimated by comparing all females to all males (across both lines and all generations, n= 
# 201 and 127 respectively), so the degrees of freedom for the denominator (v) were determined from these sample sizes.

pwr.f2.test(u=1, v=326, f2=0.625, sig.level = 0.05, power = NULL) # pwr=1.00

# Interaction of inbreeding coefficient, F, and sex. The inclusion of this term indicates that the effect of F differs
# between the sexes. So, this regression coefficient simply refers to the effect of F for males, while the
# regression coefficient for F (above) refers to the effect of F for females. This interpretation is evident 
# when you change the reference level of the SEX variable from female to male. Thus, the degrees of freedom for the 
# denominator (v) were determined from the number of male fish in the analysis (n=127). 

# NOTE that, to calculate actual effect size (f2), I add the regression coefficient to that of the females in order to calculate the 
# overall effect of F in the males.

# Effect size for males: 0.013 - 0.135 = -0.122
pwr.f2.test(u=1, v=125, f2=0.122, sig.level = 0.05, power = NULL) # male pwr=0.97
