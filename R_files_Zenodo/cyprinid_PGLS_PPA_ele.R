setwd()#set working directory to file source

#load packages
library(car)
library(MASS)
library(ggplot2)
library(phytools)
library(geiger)
library(phylolm)
library(phylobase)
library(nlme)
library(MuMIn)
library(phylopath)
library(phytools)
library(lme4)
library(lmerTest)
library(lsmeans)
library(plyr)
library(letsR)
library(phylosignal)

#upload dataset
cyp <- read.csv('cyp.SSD.CTmax.csv')
rownames(cyp) <- cyp$binom

#upload pruned phylogeny
cyp.phylo <- read.newick('phylo.for.cyp.SSD.CTmax.tre')
summary(cyp.phylo)


##PGLS to test for correlations between CTmax, SSD, and mating system

## compare phylogenetic transformations. 
##Model with pagel's lambda starting at 1.0 does not converge, produces an error.
ct.ssd.m.pagel1 <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(1, phy = cyp.phylo), data = cyp, method = 'ML', na.action = na.fail) #pagel's lambda starting value set to 1

#Here, one pagel's lambda model has a starting value of 0.7 and another pagel's lambda model a starting value of 0.4, but they produce different results (look at AIC scores)
#This suggests that we should probably avoid using a model that separately estimates pagel's lambda for this analysis because the model does not converge.
#Build models with Brownian motion, pagel's l starting at .7 and .4, and no phylogeny.
ct.ssd.m.mod00a <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(1, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #brownian motion
ct.ssd.m.mod00b <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(0.7, phy = cyp.phylo), data = cyp, method = 'ML', na.action = na.fail) #pagel's lambda starting value set to .7
ct.ssd.m.mod00c <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(0.4, phy = cyp.phylo), data = cyp, method = 'ML', na.action = na.fail) #pagel's lambda starting value set to .4
ct.ssd.m.mod00d <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(0, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #star phylogeny

#compare using AIC.
#no phylogeny and brownian motion are pretty similar, but no phylogeny has stronger support. Pagel's lambda is unstable, different starting lamdba values give different results. Just stick with brownian motion or no phylogeny.
AIC(ct.ssd.m.mod00a, ct.ssd.m.mod00b, ct.ssd.m.mod00c, ct.ssd.m.mod00d)
#                df      AIC
#ct.ssd.m.mod00a  5 132.3289
#ct.ssd.m.mod00b  6 111.0244
#ct.ssd.m.mod00c  6 115.5364
#ct.ssd.m.mod00d  5 129.5556

#Use likelihood ratio tests to determine significance of each effect in the model.
#First, let's do a model with no phylogeny (i.e. the best-supported model)
#look at effect of mean size 
ct.ssd.m.mod01a <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(0, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #no phylogeny
ct.ssd.m.mod01b <- gls(CTmax ~ log.msl.fsl + spawing.mode, correlation = corPagel(0, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #no phylogeny
anova(ct.ssd.m.mod01a, ct.ssd.m.mod01b) 

#significant effect of mean size on ctmax
#                Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#ct.ssd.m.mod01a     1  5 129.5556 135.6500 -59.77781                        
#ct.ssd.m.mod01b     2  4 131.4577 136.3332 -61.72886 1 vs 2 3.902088  0.0482

#look at effect of mating system
ct.ssd.m.mod01a <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(0, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #no phylogeny
ct.ssd.m.mod01c <- gls(CTmax ~ log.msl.fsl + mean.size, correlation = corPagel(0, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #no phylogeny
anova(ct.ssd.m.mod01a, ct.ssd.m.mod01c) 

#significant correlation of mating system and ctmax
#                Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#ct.ssd.m.mod01a     1  5 129.5556 135.6500 -59.77781                        
#ct.ssd.m.mod01c     2  4 131.5749 136.4504 -61.78747 1 vs 2 4.019309   0.0450

# look at effect of ssd
ct.ssd.m.mod01a <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(0, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #no phylogeny
ct.ssd.m.mod01d <- gls(CTmax ~ spawing.mode + mean.size, correlation = corPagel(0, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #no phylogeny
anova(ct.ssd.m.mod01a, ct.ssd.m.mod01d) 

#significant correlation of ssd and ctmax
#                Model df      AIC      BIC    logLik   Test   L.Ratio p-value
#ct.ssd.m.mod01a     1  5 129.5556 135.6500 -59.77781                         
#ct.ssd.m.mod01d     2  4 128.3780 133.2535 -60.18900 1 vs 2 0.8223758  0.3645


## now let's do a model with Brownian motion

ct.ssd.m.mod02a <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(1, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #brownian motion
ct.ssd.m.mod02b <- gls(CTmax ~ log.msl.fsl + spawing.mode, correlation = corPagel(1, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #brownian motion
anova(ct.ssd.m.mod02a, ct.ssd.m.mod02b) 

#No effect of mean size on ctmax
#                Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#ct.ssd.m.mod02a     1  5 132.3289 138.4233 -61.16445                        
#ct.ssd.m.mod02b     2  4 133.3123 138.1878 -62.65614 1 vs 2 2.983375  0.0841

ct.ssd.m.mod02a <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(1, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #brownian motion
ct.ssd.m.mod02c <- gls(CTmax ~ log.msl.fsl + mean.size, correlation = corPagel(1, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #brownian motion
anova(ct.ssd.m.mod02a, ct.ssd.m.mod02c) 

#no effect of mating system and ctmax
#                Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#ct.ssd.m.mod02a     1  5 132.3289 138.4233 -61.16445                        
#ct.ssd.m.mod02c     2  4 133.8338 138.7093 -62.91691 1 vs 2 3.504928  0.0612

ct.ssd.m.mod02a <- gls(CTmax ~ log.msl.fsl + spawing.mode + mean.size, correlation = corPagel(1, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #brownian motion
ct.ssd.m.mod02d <- gls(CTmax ~ spawing.mode + mean.size, correlation = corPagel(1, phy = cyp.phylo, fixed = TRUE), data = cyp, method = 'ML', na.action = na.fail) #brownian motion
anova(ct.ssd.m.mod02a, ct.ssd.m.mod02d) 

#no correlation of ssd and ctmax
#                Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#ct.ssd.m.mod02a     1  5 132.3289 138.4233 -61.16445                        
#ct.ssd.m.mod02d     2  4 131.6749 136.5504 -61.83743 1 vs 2 1.345966   0.246


###So, higher ctmax is correlated with pair spawning mating systems and smaller body size. Results are qualitatively the same between no phylogeny and Brownian motion, but effects are not significant under Brownian motion.


###figure to show the correlation between ctmax and mating system evolution
## need to make a new datafrom for the graph
cyp.traits <- data.frame(cyp[,c(4,2)])
rownames(cyp.traits) <- cyp$binom
colnames(cyp.traits) <- c('spawing.mode', 'CTmax')
cyp.traits

# determine how wide to make the limits on the axis
range(cyp$CTmax) # 28.37 39.30

# make a special data object with both the phylogeny and the data in it. Necessary for making the joint graph 
gphy4 <-phylo4d(cyp.phylo, cyp.traits)
gphy4

# make separate vectors for colors and shapes related to spawning mode
spawn <- ifelse(gphy4@data$spawing.mode == 'group', 'black', 'purple')
spawn2 <- ifelse(gphy4@data$spawing.mode == 'group', 19, 15)

# plot 
dotplot(gphy4, trait = 'CTmax', dot.col = spawn, tip.col = spawn, center = FALSE, scale = FALSE, data.xlim = c(28.0, 40.0), tree.ratio = 0.3, grid.horizontal = FALSE, grid.vertical = FALSE, show.box = TRUE, trait.bg.col = 'white', dot.pch = spawn2, dot.cex = 2)

#save as pdf
pdf('cyprinids.ctmax.mating.pdf', width = 5, height = 5)
dotplot(gphy4, trait = 'CTmax', dot.col = spawn, tip.col = spawn, center = FALSE, scale = FALSE, data.xlim = c(28.0, 40.0), tree.ratio = 0.3, grid.horizontal = FALSE, grid.vertical = FALSE, show.box = TRUE, trait.bg.col = 'white', dot.pch = spawn2, dot.cex = 2)
dev.off()






##### To determine most plausible causal pathways, do PPA

### Phylogenetic path analysis
###
# Define the model set with hypotheses. 
# Well established that sexual selection drives SSD in this system, so no reason to reverse the order on that. 
# Also makes sense that CTmax would only affect SSD indirectly by limiting species mean body size or indirectly by mediating the opportunity for sexual selection, or that SSD affects CTmax indirectly through mean body size, so no direct path between SSD and CTmax.
# For hypotheses showing that CTmax either constraints or compensates for trait evolution, and that SSD either drives or is constrained by the evolution of mean body size.
cyp.mods <- define_model_set(
  #ctmax compensates for size evolution, ssd drives evolution of mean size
  therm.comp.ssd.drive = c(CTmax ~ spawing.mode + mean.size, mean.size ~ log.msl.fsl),
  #ctmax constrains size evolution, mean size constrains ssd evolution
  therm.const.ssd.drive = c(spawing.mode ~ CTmax, mean.size ~ log.msl.fsl + CTmax),
  #ctmax compensates for size evolution, ssd drives evolution of mean size
  therm.comp.size.const = c(CTmax ~ spawing.mode + mean.size, log.msl.fsl ~ mean.size),
  #ctmax constrains size evolution, mean size constrains ssd evolution
  therm.const.size.const = c(spawing.mode ~ CTmax, mean.size ~ CTmax, log.msl.fsl ~ mean.size),
  .common = c(log.msl.fsl ~ spawing.mode)
)

# look at the models to make sure we have correctly specified the right pathways. 
plot_model_set(cyp.mods)

# now we run the path analysis. Exploratory analyses suggested Pagel's lambda was unstable, let's just stick with Brownian motion because you have to specify a model of evolution and Brownian motion had relatively high support in the PGLS. 
cyp.result <- phylo_path(cyp.mods, data = cyp, tree = cyp.phylo, model = 'BM')
cyp.sum <- summary(cyp.result)
cyp.sum
#                                        model k q    C      p CICc delta_CICc      l      w
#therm.const.ssd.drive   therm.const.ssd.drive 2 8 2.26 0.6889 27.3       0.00 1.0000 0.7710
#therm.comp.ssd.drive     therm.comp.ssd.drive 2 8 5.42 0.2469 30.4       3.16 0.2055 0.1585
#therm.comp.size.const   therm.comp.size.const 2 8 8.11 0.0876 33.1       5.86 0.0535 0.0412
#therm.const.size.const therm.const.size.const 2 8 8.79 0.0665 33.8       6.54 0.0381 0.0293


#one model has high akaike weight (0.77) and no other models have delta cCIC < 2. The best model is the hypothesis that ctmax constrains mating system and mean size evolution, and SSD drives the evolution of mean size. 
#bootstrap the best model with 500 iterations to get confidence intervals
cyp.best <- best(cyp.result, boot = 500)
#Take a look at the best model with path coefficients
plot(cyp.best)
#make it a little prettier
plot(cyp.best, text_size = 4, type = 'color', edge_width = 0.75, box_x = 20, box_y = 15, flip_x = TRUE, show.legend = FALSE)

#See which pathways show significant effects (confidence intervals don't overlap zero)
cyp.best.coef <- coef_plot(cyp.best, error_bar = "ci", order_by = "strength") + ggplot2::theme_classic()
pdf('cyp.best.coef.pdf', width = 6.5, height = 3)
cyp.best.coef
dev.off()

#values
cyp.best
#$coef
#             CTmax spawing.mode log.msl.fsl  mean.size
#CTmax            0      1.07047   0.0000000 -0.3829876
#spawing.mode     0      0.00000   0.5408466  0.0000000
#log.msl.fsl      0      0.00000   0.0000000  0.3767828
#mean.size        0      0.00000   0.0000000  0.0000000

#$se
#             CTmax spawing.mode log.msl.fsl mean.size
#CTmax            0     0.579264   0.0000000 0.1411369
#spawing.mode     0     0.000000   0.4128482 0.0000000
#log.msl.fsl      0     0.000000   0.0000000 0.1663071
#mean.size        0     0.000000   0.0000000 0.0000000

#$lower
#             CTmax spawing.mode log.msl.fsl   mean.size
#CTmax            0    0.0663799   0.0000000 -0.65811146
#spawing.mode     0    0.0000000  -0.2476349  0.00000000
#log.msl.fsl      0    0.0000000   0.0000000  0.09434247
#mean.size        0    0.0000000   0.0000000  0.00000000

#$upper
#             CTmax spawing.mode log.msl.fsl  mean.size
#CTmax            0     2.555531    0.000000 -0.1277796
#spawing.mode     0     0.000000    1.323128  0.0000000
#log.msl.fsl      0     0.000000    0.000000  0.6820516
#mean.size        0     0.000000    0.000000  0.0000000


