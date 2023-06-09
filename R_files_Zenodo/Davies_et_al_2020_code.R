###-----------------------------------------------------------------------------------------------###
### Davies et al., 2020 Full code                                                                 ###
### Reference: Davies, Lewis & Dougherty (2020) A meta-analysis of factors influencing the        ###
### strength of mate choice copying in animals. Behavioural Ecology.                              ###
### Author: Liam Dougherty, University of Liverpool [liam.dougherty@liv.ac.uk]                    ###
### Date: 05/05/2020                                                                              ###
###-----------------------------------------------------------------------------------------------###


#-------------------------#
# 1. Setup                #
#-------------------------#

library(ape)
library(rotl)
library(metafor)
library(MCMCglmm)

data <- read.delim("MCC_meta_analysis_final.txt", h=T) 
data$Study_no <- as.factor(data$Study_no) 
data$species2 <- as.factor(data$Species_latin)
data$Species_latin <- gsub(" ", "_", data$Species_latin) # Remove spaces between binomial names
data$Species_latin <- as.factor(data$Species_latin) 
data$animal <- as.factor(data$Species_latin) # MCMCglmm model needs column called 'animal' for phylogeny
data$Effect_size_no <- as.factor(data$Effect_size_no)
data$obs <- as.factor(data$Effect_size_no)
precision <- sqrt(1/data$Variance) # Precision= inverse standard error
data[,"precision"] <- precision 

tree2 <- read.tree("tree1_edited.tre") # Import tree
tree2_grafen <- compute.brlen(tree2, method="Grafen", power=1)
matrix <- vcv(tree2_grafen, cor=TRUE, model="Brownian")

load("meta2b.rda") # Load MCMCglmm model to save time


#----------------------------------#
# 2. Creating phylogenetic tree    #
#----------------------------------#

# [Only need to do this once]

taxa1 <- levels(data$Species_latin) 
resolved_names1 <- tnrs_match_names(taxa1,context_name = "Animals") # Check OTL for species names

# Plot tree using TOL data
tree1 <- tol_induced_subtree(ott_ids = resolved_names1$ott_id) 
tree1$tip.label <- strip_ott_ids(tree1$tip.label) # Remove ott IDs for presentation
tree1$node.label <- NULL # Remove node labels (these might be a problem later)

write.nexus(tree1, file="tree1.nex")


#-------------------------------------#
# 3. Meta-analysis- overall models    #
#-------------------------------------#

# Simplest model- no random effects
meta1 <- rma.uni(Hedges_d_directional, Variance, data= data, method= "REML") 

# Funnel plot
funnel(meta1, yaxis="seinv", xlab="Effect size (Hedges' d)", steps=11, digits=1, ylim=c(1, 11), xlim=c(-2.2, 5), 
       back="white", shade="white", hlines="white", pch=21, col=rgb(0,0,0, max=255), bg=rgb(255,255,255, max=255, alpha=150))
abline(v=0, lty=2) 

### Overall model- including random effects
meta2 <- rma.mv(Hedges_d_directional, Variance, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                R= list(Species_latin = matrix), data= data, method= "REML")
summary(meta2)


### Heteogeneity (I^2)
W <- diag(1/data$Variance)
X <- model.matrix(meta2)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W

100 * sum(meta2$sigma2) / (sum(meta2$sigma2) + 
                             (meta2$k-meta2$p)/sum(diag(P))) # Total I^2
100 * meta2$sigma2 / (sum(meta2$sigma2) + (meta2$k-meta2$p)/sum(diag(P))) # I^2 for each random factor


### MCMCglmm model
prior3 <- list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)))
meta2b <- MCMCglmm(Hedges_d_directional~1, random= ~Study_no + Species_latin + animal, mev=data$Variance, data=data,
                   nitt=300000, thin=50, burnin=200000, pr=TRUE, prior=prior3, pedigree=tree2_grafen)


### Meta-analytic residuals [already present in datafile]
meta2b$Random$formula<-update(meta2b$Random$formula, ~.+leg(mev, -1, FALSE):units)
fitted <- predict(meta2b, marginal=~leg(mev, -1, FALSE):units) # Only works using R v3.2 or earlier
residuals <- data$Hedges_d_directional - fitted # raw data - predictions = meta-analytic residuals
zresiduals <- residuals*precision
data$zresiduals <- zresiduals
write.table(data, file= "newdata.txt", sep="\\t") # Write new data file

# Funnel plot
meta_resid <- rma.uni(zresiduals, Variance, data= data, method= "REML") 
funnel(meta_resid, yaxis="seinv", xlab="Effect size (Hedges' d)", steps=11, digits=1, ylim=c(1, 11),
       back="white", shade="white", hlines="white", pch=21, col=rgb(0,0,0, max=255), bg=rgb(255,255,255, max=255, alpha=150))


#-----------------------------#
# 4. Effect of moderators     #
#-----------------------------#

### Single-factor meta-regression models
meta_group <- rma.mv(Hedges_d_directional, Variance, mod= ~ Group, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")
meta_sex <- rma.mv(Hedges_d_directional, Variance, mod= ~ Choosing_sex, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")
meta_design <- rma.mv(Hedges_d_directional, Variance, mod= ~ Design_type, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")
meta_pairing <- rma.mv(Hedges_d_directional, Variance, mod= ~ Model_pairing, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")
meta_demo <- rma.mv(Hedges_d_directional, Variance, mod= ~ Demonstration_type, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")
meta_behaviour <- rma.mv(Hedges_d_directional, Variance, mod= ~ Behaviour_measured, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, 
                     ~ 1|obs), R= list(Species_latin = matrix), data= data, method= "REML")
meta_copy_type <- rma.mv(Hedges_d_directional, Variance, mod= ~ Generalised_or_individual, random= list(~ 1|Species_latin, ~ 1|species2, 
                     ~ 1|Study_no, ~ 1|obs), R= list(Species_latin = matrix), data= data, method= "REML")
meta_birth <- rma.mv(Hedges_d_directional, Variance, mod= ~ Animal_born, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                         R= list(Species_latin = matrix), data= data, method= "REML")

summary(meta_group)
summary(meta_sex)
summary(meta_design)
summary(meta_pairing)
summary(meta_demo)
summary(meta_behaviour)
summary(meta_copy_type)
summary(meta_birth)


### To determine means and 95% CIs for each factor level- run same model as before but remove intercept
meta_group2 <- rma.mv(Hedges_d_directional, Variance, mod= ~ Group -1, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")
meta_sex2 <- rma.mv(Hedges_d_directional, Variance, mod= ~ Choosing_sex -1, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")
meta_design2 <- rma.mv(Hedges_d_directional, Variance, mod= ~ Design_type -1, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")
meta_pairing2 <- rma.mv(Hedges_d_directional, Variance, mod= ~ Model_pairing -1, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")
meta_demo2 <- rma.mv(Hedges_d_directional, Variance, mod= ~ Demonstration_type -1, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")
meta_behaviour2 <- rma.mv(Hedges_d_directional, Variance, mod= ~ Behaviour_measured -1, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, 
                     ~ 1|obs), R= list(Species_latin = matrix), data= data, method= "REML")
meta_copy_type2 <- rma.mv(Hedges_d_directional, Variance, mod= ~ Generalised_or_individual -1, random= list(~ 1|Species_latin, ~ 1|species2, 
                   ~ 1|Study_no, ~ 1|obs), R= list(Species_latin = matrix), data= data, method= "REML")
meta_birth2 <- rma.mv(Hedges_d_directional, Variance, mod= ~ Animal_born -1, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                     R= list(Species_latin = matrix), data= data, method= "REML")

summary(meta_group2)
summary(meta_sex2)
summary(meta_design2)
summary(meta_pairing2)
summary(meta_demo2)
summary(meta_behaviour2)
summary(meta_copy_type2)
summary(meta_birth2) 


### Mating rate- reduced dataset
mating_data <- subset(data, Rate_of_multiple_mating=="High" | Rate_of_multiple_mating=="Low")

# Getting pruned tree
species_mating <- unique(mating_data$Species_latin) # List of species names
tree_mating <- drop.tip(tree2_grafen, tree2_grafen$tip.label[-match(species_mating, tree2_grafen$tip.label)])
matrix2 <- vcv(tree_mating, cor=TRUE, model="Brownian")

# Meta-regressions
meta_mating_rate <- rma.mv(Hedges_d_directional, Variance, mod= ~ Rate_of_multiple_mating, random= list(~ 1|Species_latin, ~ 1|species2, 
                    ~ 1|Study_no, ~ 1|obs), R= list(Species_latin = matrix), data= mating_data, method= "REML")
summary(meta_mating_rate)

meta_mating_rate2 <- rma.mv(Hedges_d_directional, Variance, mod= ~ Rate_of_multiple_mating -1, random= list(~ 1|Species_latin, ~ 1|species2, 
                     ~ 1|Study_no, ~ 1|obs), R= list(Species_latin = matrix), data= mating_data, method= "REML")
summary(meta_mating_rate2)


#-------------------------#
# 5. Publication bias     #
#-------------------------#

### Effect size over time
meta_year <- rma.mv(Hedges_d_directional, Variance, mod= ~ Year, random= list(~ 1|Species_latin, ~ 1|species2, ~ 1|Study_no, ~ 1|obs), 
                    R= list(Species_latin = matrix), data= data, method= "REML")
summary(meta_year)

# Bubble plot
preds <- predict(meta_year, newmods=c(1992:2018))
data$size <- (1 / sqrt(data$Variance))/10 

plot(NA, NA, xlim=c(1992,2018), ylim=c(-2.6,4.5),
     xlab="Publication year", ylab="Effect size (Hedges' d)",
     las=1, bty="l")
symbols(data$Year, data$Hedges_d_directional, circles=data$size, inches=FALSE, add=TRUE, 
        fg= rgb(0,0,0, max=255, alpha=150), pch=21, col=rgb(0,0,0, max=255), bg=rgb(255,255,255, max=255, alpha=150))
lines(1992:2018, preds$pred)
lines(1992:2018, preds$ci.lb, lty="dashed")
lines(1992:2018, preds$ci.ub, lty="dashed")
abline(h=0, lty="dotted")


### Trim and fill

# Raw data
trimfill_left <- trimfill(meta1, side="left") 
summary(trimfill_left) 

trimfill_right <- trimfill(meta1, side="right")
summary(trimfill_right)
funnel(trimfill_right, yaxis="seinv", xlab="Effect size (Hedges' d)", steps=11, digits=1, ylim=c(1, 11), xlim=c(-2.2, 5), 
       back="white", shade="white", hlines="white", pch=21, col=rgb(0,0,0, max=255), bg=rgb(255,255,255, max=255, alpha=150), 
       pch.fill=16)
abline(v=0, lty=2)

# Residuals
trimfill_left2 <- trimfill(meta_resid, side="left") 
summary(trimfill_left2) 

trimfill_right2 <- trimfill(meta_resid, side="right")
summary(trimfill_right2)
funnel(trimfill_right2, yaxis="seinv", xlab="Residual effect size", steps=11, digits=1, ylim=c(1, 11),back="white", shade="white", 
       hlines="white", pch=21, col=rgb(0,0,0, max=255), bg=rgb(255,255,255, max=255, alpha=150), pch.fill=16)
abline(v=0, lty=2)


### Egger's regression

# Raw data
egger <- lm(precision ~ Hedges_d_directional, data= data)
summary(egger)
plot(data$Hedges_d_directional, precision, xlab="Effect size (Hedges' d)", ylab="Inverse standard error", ylim=c(0,10),
      pch=21, col=rgb(0,0,0, max=255), bg=rgb(255,255,255, max=255, alpha=150)) 
abline(egger)

# Residuals
egger2 <- lm(precision ~ zresiduals, data= data)
summary(egger2)
plot(data$zresiduals, data$precision, xlab="Residual effect size", ylab="Inverse standard error", ylim=c(0,10), 
     pch=21, col=rgb(0,0,0, max=255), bg=rgb(255,255,255, max=255, alpha=150)) 
abline(egger2)



##########################################################################################################
##########################################################################################################