rm(list=ls())

library(vegan)
library(lme4)  

data <- read.csv("Negev_bat_data.csv", header=TRUE)

attach(data)

# specifying group level variables (group effect term)
subject<-factor(site)
season <-factor(season)
water<-factor(water)
type <-factor(type)
lights <-factor (lights)
desert <-factor(desert)


## test for correlations between all continuous variables
# significant: Tmax~Tmin; Pond length~width; 

autocor1 <- lm(log_length~log_area)
summary(autocor1)


par(mfrow=c(2,3))

## look for associations between variables

plot(log_activity~log_depth, pch=16, col= "red ",cex=2, cex.main=2,
     cex.lab=1.5, ylab= "Log bat activity", xlab= "Log pond depth")

plot(species~altitude, pch=16, col= "black", cex=1.5, ylab= " ", xlab= "Altitude (m)")

plot(log_activity~log_length, pch=16, col= "black", cex=1.5, ylab= " ", xlab= "Log pond length")

plot(species~Tmin, pch=16, col= "black", cex=1.5, ylab= "Summer species richness", xlab= "Minimum temperature")

boxplot(log_pk~water, ylab= "Activity p. kuhlii", xlab= "Water presence")

boxplot(species~water, ylab= " ", xlab= "Water presence")

boxplot(species~desert, ylab= "Species richness", xlab= "Desert")

## competitors figures:
par(mfrow=c(3,2))
boxplot(PK_log ~water, ylab= "Activity P. kuhlii", xlab= "Water presence")
plot(PK_log~log_depth, pch=16, col= "black", cex=1.5, ylab= " ", xlab= "Log pond depth")
boxplot(PB_log~water, ylab= "Activity H. ariel", xlab= "Water presence")
plot(PB_log~Tmin, pch=16, col= "black", cex=1.5, ylab= " ", xlab= "Minimum temperatures")
boxplot(PR_log~desert, ylab= "Activity P. rueppellii", xlab= "Desert")
plot(PR_log~log_length, pch=16, col= "black", cex=1.5, ylab= " ", xlab= "Log pond length")


## test for normality

shapiro.test(species)
hist(species,freq=F,col="grey")


## Generalised mixed effect models for species richness:

glm1 <- glmer(species ~ log_length +(1 | subject), family=poisson)
summary(glm1)
r.squaredGLMM(glm1)

glm2 <- glmer(species ~ water + lights +(1 | subject), family=poisson)
summary(glm2)

glm3 <- glmer(species ~ veg_index +(1 | subject), family=poisson)
summary(glm3)

glm4 <- glmer(species ~ scale(altitude) +(1 | subject), family=poisson)
summary(glm4)

glm5 <- glmer(species ~ Tmin +(1 | subject), family=poisson)
summary(glm5)

glm6 <- glmer(species ~ desert + water +(1 | subject), family=poisson)
summary(glm6)

glm7 <- glmer(species ~ water + scale(altitude) +(1 | subject), family=poisson)
summary(glm7)

glm8 <- glmer(species ~ water + Tmin +(1 | subject), family=poisson)
summary(glm8)

glm9 <- glmer(species ~ desert + water + log_length + scale(altitude) +(1 | subject), family=poisson)
summary(glm9)

anova (glmm2, glmm10)

## Generalised mixed effect models with negative binomial distribution for activity (glmer.nb):

glmm1 <- glmer.nb(PB_log ~ desert +(1 | subject))
summary(glmm1)
r.squaredLR(glmm1)
r.squaredGLMM(glmm1)

glmm10 <- glmer.nb(Pk_log ~ desert + Tmin +(1 | subject))
summary(glmm10)

glmm100 <- glmer.nb(Pk_log ~ log_length +(1 | subject))
summary(glmm100)

glmm2 <- glmer.nb(Pk_log ~ Tmin +(1 | subject))
summary(glmm2)

glmm3 <- glmer.nb(Pk_log ~ scale(altitude) +(1 | subject))
summary(glmm3)

glmm30 <- glmer.nb(Pk_log ~ water +(1 | subject))
summary(glmm30)

glmm4 <- glmer.nb(Pk_log ~ scale(altitude) + Tmin +(1 | subject))
summary(glmm4)

glmm200 <- glmer.nb(Pk_log ~ water + Tmin +(1 | subject))
summary(glmm200)

glmm201 <- glmer.nb(Pk_log ~ water + scale(altitude) +(1 | subject))
summary(glmm201)

glmm21 <- glmer.nb(Pk_log ~ water + scale(altitude) + Tmin + (1 | subject))
summary(glmm21)

glmm22 <- glmer.nb(Pk_log ~ desert + water +(1 | subject))
summary(glmm22)

glmm23 <- glmer.nb(Pk_log ~ desert + water + Tmin +(1 | subject))
summary(glmm23)

glmm24 <- glmer.nb(Pk_log ~ water + scale(altitude) + Tmin + log_length +(1 | subject))
summary(glmm24)

glmm_pr <- glmer.nb(Pk_log ~ water + scale(altitude) + Tmin + PR_log + (1 | subject))
summary(glmm_pr)

glmm_pr2 <- glmer.nb(Pk_log ~ PR_log + (1 | subject))
summary(glmm_pr2)

glmm_pk <- glmer.nb(Pk_log ~ water + scale(altitude) + Tmin + PK_log + (1 | subject))
summary(glmm_pk)

glmm_pk2 <- glmer.nb(Pk_log ~ PK_log + (1 | subject))
summary(glmm_pk2)


glmm5 <- glmer.nb(Pk_log ~ desert + scale(altitude) + Tmin +(1 | subject))
summary(glmm5)


glmm6 <- glmer.nb(activity ~ season + water + log_depth +(1 | subject))
summary(glmm6)

glmm7 <- glmer.nb(activity ~ season + water + log_alt +(1 | subject))
summary(glmm7)

glmm8 <- glmer.nb(activity ~ season + water + log_alt + Tmin +(1 | subject))
summary(glmm8)

glmm9 <- glmer.nb(activity ~ season + water + log_depth + log_alt + Tmax +(1 | subject))
summary(glmm9)

glmm10 <- glmer.nb(activity ~ season + log_depth + log_alt +(1 | subject))
summary(glmm10)

glmm11 <- glmer.nb(activity ~ Tmin + (1 | subject))
summary(glmm11)

glmm12 <- glmer.nb(activity ~ water + log_depth + Tmax + (1 | subject))
summary(glmm12)

glmm13 <- glmer.nb(activity ~ water + Tmax + (1 | subject))
summary(glmm13)

glmm14 <- glmer.nb(activity ~ water  + log_alt + (1 | subject))
summary(glmm14)

glmm15 <- glmer.nb(activity ~ water + log_depth + Tmax + (1 | season:subject))
summary(glmm15)

anova (glmm4, glmm2)

