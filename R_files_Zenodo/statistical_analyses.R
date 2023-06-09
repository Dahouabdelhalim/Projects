library(phytools)
library(phangorn)
library(ks)
library(ape)
library(geiger)
library(nlme)
library(effects)
library(car)
library(caper)
library(MuMIn)
library(corrplot)
library(geomorph)

#-----------#
#import tree#
#-----------#
small_tree<-read.nexus(file = "path/small_phylogeny.nex")

#-------------------------#
# reading in trait data  #
#-------------------------#
mam <- read.csv("path/trait_data.csv", fileEncoding = "UTF-8", stringsAsFactors = F)

#-----------------------------#
# trim by species of interest #
#-----------------------------#
grinnell_mammals <- c("Ammospermophilus_leucurus","Chaetodipus_californicus","Chaetodipus_fallax","Chaetodipus_formosus","Chaetodipus_penicillatus","Dipodomys_deserti","Dipodomys_merriami","Dipodomys_microps","Dipodomys_panamintinus","Dipodomys_ordii","Microtus_californicus","Microtus_longicaudus","Neotoma_lepida","Neotoma_macrotis","Onychomys_torridus","Otospermophilus_beecheyi","Perognathus_longimembris","Peromyscus_boylii","Peromyscus_crinitus","Peromyscus_eremicus","Peromyscus_maniculatus","Peromyscus_truei","Reithrodontomys_megalotis","Neotamias_merriami","Neotamias_panamintinus")
grinnell_mammals_traits <- mam[mam$Binomial.1.2 %in% grinnell_mammals, ]

#---------#
# pruning #
#---------#
tips <- grinnell_mammals_traits$Binomial.1.2
pruned_small<-lapply(small_tree,function(t,tips)drop.tip(t,setdiff(t$tip.label,tips)),tips=tips)
class(pruned_small)<-"multiPhylo"

#----------#
# avg tree #
#----------#
ave.tree <- midpoint(ls.consensus(pruned_small))

#read df
df <- read.csv("path/pelage_properties_spp.csv")

#check names
row.names(df) <- c("Ammospermophilus_leucurus","Chaetodipus_californicus","Chaetodipus_fallax","Chaetodipus_formosus","Chaetodipus_penicillatus","Dipodomys_deserti","Dipodomys_merriami","Dipodomys_microps","Dipodomys_panamintinus","Dipodomys_ordii","Microtus_californicus","Microtus_longicaudus","Neotoma_lepida","Neotoma_macrotis","Onychomys_torridus","Otospermophilus_beecheyi","Perognathus_longimembris","Peromyscus_boylii","Peromyscus_crinitus","Peromyscus_eremicus","Peromyscus_maniculatus","Peromyscus_truei","Reithrodontomys_megalotis","Neotamias_merriami","Neotamias_panamintinus")
name.check(ave.tree, df) #checks to see that df matches the tree

#################################
#    phylogenetic analyses      #
#################################

df$log_mass_scaled <- scale(df$log_mass, center=TRUE,scale=TRUE)
df$insulation_scaled <- scale(df$insulation, center=TRUE,scale=TRUE)
df$depth_dorsal_scaled <- scale(df$depth_dorsal, center=TRUE,scale=TRUE)
df$length_dorsal_scaled <- scale(df$length_dorsal, center=TRUE,scale=TRUE)

#final model for conductivity
model4<-pgls(conductivity ~ log_mass_scaled + length_dorsal_scaled + activity + habitat_peterson2, data = comparative_data_object, lambda ='ML')
summary(model4)
modaov <- procD.pgls(conductivity ~ log_mass_scaled + length_dorsal_scaled + activity + habitat_peterson2, phy = ave.tree, data = df,iter = 999, SS.type = c("II"))
anova(modaov)

#now insulation
model1<-pgls(insulation ~ habitat_peterson2 + activity + depth_dorsal_scaled + log_mass_scaled + length_dorsal_scaled, data = comparative_data_object, lambda ='ML')
summary(model1)
modaov <- procD.pgls(insulation ~ habitat_peterson2 + activity + depth_dorsal_scaled + log_mass_scaled + length_dorsal_scaled, phy = ave.tree, data = df,iter = 999, SS.type = c("II"))
anova(modaov)

#now let's calculate the residuals for insulation
mod1<-pgls(insulation ~ habitat_peterson2 + activity + log_mass_scaled + length_dorsal_scaled, data = comparative_data_object, lambda ='ML')
mod2<-pgls(depth_dorsal_scaled ~ habitat_peterson2 + activity + log_mass_scaled + length_dorsal_scaled, data = comparative_data_object, lambda ='ML')
resid(mod1)
resid(mod2)

df <- cbind(resid(mod1),resid(mod2))
plot(df)
write.csv(df,file = "Desktop/insulation_depth.csv")

#now let's calculate the residuals for conductivity
mod1<-pgls(conductivity ~ log_mass_scaled + activity + habitat_peterson2, data = comparative_data_object, lambda ='ML')
mod2<-pgls(length_dorsal_scaled ~ log_mass_scaled + activity + habitat_peterson2, data = comparative_data_object, lambda ='ML')
resid(mod1)
resid(mod2)

df <- cbind(resid(mod1),resid(mod2))
plot(df)
write.csv(df,file = "Desktop/conductivity_length.csv")

mod1<-pgls(conductivity ~ length_dorsal_scaled + habitat_peterson2, data = comparative_data_object, lambda ='ML')

##############################
#run PGLS for all 1000 trees #
##############################

#read df
df <- read.csv("path/pelage_properties_spp.csv")

#check names
row.names(df) <- c("Ammospermophilus_leucurus","Chaetodipus_californicus","Chaetodipus_fallax","Chaetodipus_formosus","Chaetodipus_penicillatus","Dipodomys_deserti","Dipodomys_merriami","Dipodomys_microps","Dipodomys_panamintinus","Dipodomys_ordii","Microtus_californicus","Microtus_longicaudus","Neotoma_lepida","Neotoma_macrotis","Onychomys_torridus","Otospermophilus_beecheyi","Perognathus_longimembris","Peromyscus_boylii","Peromyscus_crinitus","Peromyscus_eremicus","Peromyscus_maniculatus","Peromyscus_truei","Reithrodontomys_megalotis","Neotamias_merriami","Neotamias_panamintinus")
name.check(ave.tree, df) #checks to see that df matches the tree

#comparative data object
comparative_data_object <- comparative.data(phy = ave.tree, data = df, names.col = spp, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

uncertainty <- function(df,trees){
  params = c()
  for (i in seq(1000))
  {
    comparative_data_object <- comparative.data(phy = trees[[i]], data = df, names.col = spp, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
    mod<-try(pgls(conductivity ~ log_mass_scaled + length_dorsal_scaled + activity + habitat_peterson2, data = comparative_data_object, lambda ='ML'),silent = TRUE)
    if (mod == "Error in pgls(conductivity ~ depth_dorsal_scaled + habitat: \\n  Problem with optim:52ERROR: ABNORMAL_TERMINATION_IN_LNSRCH\\n"){
    }
    else {params = rbind(params, data.frame(row.names = NULL, tree = i, logLik = mod$model$log.lik,
                                            AIC = mod$aicc, lambda = mod$param[2], lambda_bounds_low = mod$param.CI$lambda$bounds.val[1], 
                                            lambda_bounds_high = mod$param.CI$lambda$bounds.val[2], lambda_p_low = mod$param.CI$lambda$bounds.p[1],
                                            lambda_p_high = mod$param.CI$lambda$bounds.p[2], Int = mod$model$coef[1], Int_err = mod$sterr[1],
                                            beta = mod$model$coef[2], beta_err = mod$sterr[2]))
    }
  }
  return(params)
}
results <- uncertainty(df,pruned_small)
median(results$lambda)
median(results$lambda_p_low)
median(results$lambda_p_high)

"> median(results$lambda)
[1] 1e-06
> median(results$lambda_p_low)
[1] 1
> median(results$lambda_p_high)
[1] 0.01205838"

#now insulation
uncertainty <- function(df,trees){
  params = c()
  for (i in seq(1000))
  {
    comparative_data_object <- comparative.data(phy = trees[[i]], data = df, names.col = spp, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
    mod<-try(pgls(insulation ~ log_mass_scaled + length_dorsal_scaled + depth_dorsal_scaled + activity + habitat_peterson2, data = comparative_data_object, lambda ='ML'),silent = TRUE)
    if (mod == "Error in pgls(conductivity ~ depth_dorsal_scaled + habitat: \\n  Problem with optim:52ERROR: ABNORMAL_TERMINATION_IN_LNSRCH\\n"){
    }
    else {params = rbind(params, data.frame(row.names = NULL, tree = i, logLik = mod$model$log.lik,
                                            AIC = mod$aicc, lambda = mod$param[2], lambda_bounds_low = mod$param.CI$lambda$bounds.val[1], 
                                            lambda_bounds_high = mod$param.CI$lambda$bounds.val[2], lambda_p_low = mod$param.CI$lambda$bounds.p[1],
                                            lambda_p_high = mod$param.CI$lambda$bounds.p[2], Int = mod$model$coef[1], Int_err = mod$sterr[1],
                                            beta = mod$model$coef[2], beta_err = mod$sterr[2]))
    }
  }
  return(params)
}
results <- uncertainty(df,pruned_small)
median(results$lambda)
median(results$lambda_p_low)
median(results$lambda_p_high)
"> median(results$lambda)
[1] 1e-06
> median(results$lambda_p_low)
[1] 1
> median(results$lambda_p_high)
[1] 0.0004420929"

#############################
# non-phylogenetic approach #
#############################

#read data
data <- read.csv("path/conductance_nonphylo.csv")

#scale and center
data$log_mass_scaled <- scale(data$log_mass, center=TRUE,scale=TRUE)
data$insulation_scaled <- scale(data$insulation, center=TRUE,scale=TRUE)
data$depth_dorsal_scaled <- scale(data$dorsal_fur_depth, center=TRUE,scale=TRUE)
data$length_dorsal_scaled <- scale(data$dorsal_fur_length, center=TRUE,scale=TRUE)

#scaled analyses
mod1 <- lm(conductivity ~ log_mass_scaled + length_dorsal_scaled + activity + habitat, data=data)
Anova(mod1,type="2")

mod2 <- lm(insulation ~ log_mass_scaled + length_dorsal_scaled + depth_dorsal_scaled + activity + habitat, data=data)
Anova(mod2,type="2")

#non-scaled analyses
mod1 <- lm(conductivity ~ log_mass_scaled + dorsal_fur_length + activity + habitat, data=data)
Anova(mod1,type="2")

mod2 <- lm(insulation ~ log_mass_scaled + length_dorsal_scaled + dorsal_fur_depth + activity + habitat, data=data)
Anova(mod2,type="2")

#test for an association between pelage thickness and age of specimen
data$years_scaled <- scale(data$years_in_museum, center=TRUE,scale=TRUE)
mod_depth <- lmer(dorsal_fur_depth ~ years_scaled + (1|spp),data=data)
summary(mod_depth)
anova(mod_depth,type="2")
sjstats::omega_sq(mod_depth)

#test for an association between conductivity and age of specimen, accounting for significant effects
mod1 <- lm(conductivity ~ dorsal_fur_length + habitat +years_scaled, data=data)
Anova(mod1,type="2")
sjstats::omega_sq(mod1)

#test for an association between insulation and age of specimen, accounting for significant effects
mod1 <- lm(insulation ~ dorsal_fur_depth + habitat + years_scaled, data=data)
Anova(mod1,type="2")
sjstats::omega_sq(mod1)

#test for an association between insulation and seasonality, accounting for significant effects
mod1 <- lm(insulation ~ dorsal_fur_depth + habitat + seasonality, data=data)
Anova(mod1,type="2")
sjstats::omega_sq(mod1)

#test for an association between insulation and seasonality, accounting for significant effects
mod1 <- lm(conductivity ~ dorsal_fur_length + habitat + seasonality, data=data)
Anova(mod1,type="2")
sjstats::omega_sq(mod1)

#Mann Whitney for fresh vs dried
data <- read.csv("path/fresh_vs_dried.csv")
shapiro.test(data$conductance_wet)
t.test(conductance_wet ~ trt, paired = TRUE, alternative = "less",data=data)

#t.test for heat savings, nocturnal
data <- read.csv("path/heating_costs_diurnal.csv")
shapiro.test(data$Qgen_heat_sum)
t.test(Qgen_heat_sum ~ species, paired = TRUE,data=data)

#t.test for heat savings, diurnal
data <- read.csv("path/heating_costs_nocturnal.csv")
shapiro.test(data$Qgen_heat_sum)
t.test(Qgen_heat_sum ~ species, paired = TRUE,data=data)


