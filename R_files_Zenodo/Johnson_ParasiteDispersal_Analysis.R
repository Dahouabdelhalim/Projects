### Analysis for Johnson et al. "The cost of travel: How dispersal ability limits local adaptation in host-parasite interactions
### 7/7/2020

####### LIBRARIES ######
library(tidyverse)
library(lme4)
library(MASS)
library(multcomp)

################################################
####### Pseudacris regilla (host) analysis  ####
################################################
exp.dat.psre <- read.csv("data/Experimental_Crosses_PSRE.csv")
exp.dat.psre %>% mutate(
  fhost = as.factor(Code),
  fobs = as.factor(Observation),
  fhost_source = as.factor(Host_Pond),
  fparas_source = as.factor(Parasite_Source)
) -> exp.dat.psre

exp.dat.psre %>% filter(!is.na(SVL)) -> exp.dat.psre


##### Run GLMM ################################
#### Define best random effect structure
# host as RE
home1 <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ scale(SVL) + 
                 Parasite*scale(dist_host_paras) + 
                 (1 | fhost_source) + (1 | fparas_source) + (1 | fhost),
               family = binomial, data = exp.dat.psre,
               control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# obs as RE
home2 <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ scale(SVL) + 
                 Parasite*scale(dist_host_paras) + 
                 (1 | fhost_source) + (1 | fparas_source) + (1 | fobs),
               family = binomial, data = exp.dat.psre,
               control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# both as RE
home3 <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ scale(SVL) + 
                 Parasite*scale(dist_host_paras) + 
                 (1 | fhost_source) + (1 | fparas_source) + (1|fobs) + (1|fhost),
               family = binomial, data = exp.dat.psre,
               control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
AIC(home1, home2, home3) # using obs as RE is better

#### explore model output using likelihood ratio tests
summary(home2)

# test interaction significance using LRT
home2_red <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ scale(SVL) + 
                     Parasite + scale(dist_host_paras) + 
                     (1 | fhost_source) + (1 | fparas_source) + (1 | fobs),
                   family = binomial, data = exp.dat.psre,
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(home2, home2_red)

# test significance of distance for mano only using reduced model
# need to manually remove distance for just mano
exp.dat.psre$ScaledDist <- scale(exp.dat.psre$dist_host_paras)
exp.dat.psre %>% mutate(ribDist = case_when(
  Parasite == "Mano"~0,
  Parasite=="Rib"~ScaledDist
)) -> exp.dat.psre

home2_red <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ scale(SVL) + 
                     Parasite + ribDist+
                     (1 | fhost_source) + (1 | fparas_source) + (1 | fobs),
                   family = binomial, data = exp.dat.psre,
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(home2_red)
anova(home2, home2_red)

# re-level parasites to get comparisons
exp.dat.psre$rParasite <- exp.dat.psre$Parasite
exp.dat.psre$rParasite <- factor(exp.dat.psre$rParasite, levels = c("Rib", "Mano"))
home2R <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ scale(SVL) + 
                  rParasite*scale(dist_host_paras) + 
                  (1 | fhost_source) + (1 | fparas_source) + (1 | fobs),
                family = binomial, data = exp.dat.psre,
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(home2R)

# test significance of SVL using LRT
home2_nosvl <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ 
                       Parasite*scale(dist_host_paras) + 
                       (1 | fhost_source) + (1 | fparas_source) + (1 | fobs),
                     family = binomial, data = exp.dat.psre,
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(home2, home2_nosvl )

# test whether developmental stage is significant
exp.dat.psre$Dstage <- as.numeric(exp.dat.psre$Stage)
home2_dev <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ scale(SVL)+
                     Parasite*scale(dist_host_paras) + scale(Dstage)+
                     (1 | fhost_source) + (1 | fparas_source) + (1 | fobs),
                   family = binomial, data = exp.dat.psre %>% filter(!is.na(Dstage)),
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(home2_dev)
home_noDev <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ scale(SVL) + 
                    Parasite*scale(dist_host_paras) + 
                    (1 | fhost_source) + (1 | fparas_source) + (1 | fobs),
                  family = binomial, data = exp.dat.psre %>% filter(!is.na(Dstage)),
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(home_noDev, home2_dev)

# test whether parasite species is significant
home_nosp <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ scale(SVL) + 
                     scale(dist_host_paras) + 
                     (1 | fhost_source) + (1 | fparas_source) + (1 | fobs),
                   family = binomial, data = exp.dat.psre ,
                   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
anova(home2, home_nosp)

##### Variance partitioning (RE only) ##########
exp.dat.psre.rib <- subset(exp.dat.psre, Parasite == "Rib")
homeRE.rib <- glmer(cbind(Paras_Inf, Paras_Uninf) ~ 
                      (1|fhost_source) + (1|fparas_source) + (1|fobs),
                    family = binomial,
                    data = exp.dat.psre.rib,
                    control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

randoms <- lme4::ranef(homeRE.rib, condVar = TRUE)
qq <- attr(ranef(homeRE.rib, condVar=TRUE)[[1]], "postVar")
qq.host <- attr(ranef(homeRE.rib, condVar=TRUE)$fhost_source, "postVar")
qq.para <- attr(ranef(homeRE.rib, condVar=TRUE)$fparas_source, "postVar")
rand.interc.host <- randoms$fhost_source
rand.interc.para <- randoms$fparas_source
df.host <- data.frame(Intercepts = randoms$fhost_source[,1],
                      sd.interc = 2*sqrt(qq.host[,,1:length(qq.host)]),
                      lev.names = rownames(rand.interc.host),
                      who = "Host Population")
df.para <- data.frame(Intercepts = randoms$fparas_source[,1],
                      sd.interc = 2*sqrt(qq.para[,,1:length(qq.para)]),
                      lev.names = rownames(rand.interc.para),
                      who = "Parasite Population")
df.ranef <- rbind(df.host, df.para)
df.ranef$raneflab <- paste(df.ranef$lev.names, df.ranef$who, sep = "_")
df.ranef$raneflab <- factor(df.ranef$raneflab, levels = df.ranef$raneflab[order(df.ranef$who, df.ranef$Intercepts)])
df.ranef$lev.names <- as.factor(df.ranef$lev.names)
levels(df.ranef$lev.names) <- c("BLA", "SLA", "OHL", "FLOR", "KAM", "MAL", "MQL", "FLOR", "ORP", "PEN", "PRP", "SWMP")
# variance partitioning 
summary(glmer(cbind(Paras_Inf, Paras_Uninf)~(1|fparas_source)+(1|fhost_source), data = exp.dat.psre.rib, family = "binomial"))
# prop var host source: 
(0.27)/(0.27+0.24)
(0.24)/(0.27+0.24)
##### Fig 5 (Random effects) ###################
ggplot(df.ranef, aes(raneflab, Intercepts, shape = who)) +
  geom_hline(yintercept = 0, linetype = 2, color = "darkgrey")+
  geom_errorbar(aes(ymin = Intercepts - sd.interc, ymax = Intercepts + sd.interc, color = who), width = 0)+
  geom_point(aes(size = 1, color = who))+
  guides(size = FALSE, shape = FALSE)+
  scale_x_discrete(name = "Population", labels = df.ranef$lev.names, breaks = df.ranef$raneflab)+
  scale_y_continuous(name = "Random effect", limits = c(-1.5, 1.5))+
  coord_flip()+
  scale_color_manual(name = "", values = c("forestgreen", "lightblue3"))+
  theme_classic()+
  theme(legend.position = "top")
##### Fig 3 (GLMM) #######
# plot final model 
predict.data <- data.frame(expand.grid(SVL = mean(exp.dat.psre$SVL), 
                                       dist_host_paras = seq(min(exp.dat.psre$dist_host_paras), max(exp.dat.psre$dist_host_paras), by = 10), 
                                       fparas_source = "fake",
                                       Parasite = c("Rib", "Mano"), 
                                       fhost_source = "fake", fobs = "fake"))
# predict over that fake data to get best fit lines
predict.data$predpar <- predict(home2, newdata = predict.data, allow.new.levels = TRUE, type = "response")

# use bootstrapping to get confidence interval (takes a long time -- used nsim = 999 for pub figure)
bs <- bootMer(home2, nsim = 99, function(x) predict(x, newdata = predict.data, re.form=NA, allow.new.levels = TRUE, type = "response"))
predict.data$CI.lower = apply(bs$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))
predict.data$CI.upper = apply(bs$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))

# extract random effects from model so we can plot partial residuals instead of raw data
r.e.parasite <- data.frame("ParasitePopulation" = rownames(ranef(home2)$fparas_source), "RandomEffect" = ranef(home2)$fparas_source[,1])
r.e.host <- data.frame("HostPopulation" = rownames(ranef(home2)$fhost_source), "RandomEffect" = ranef(home2)$fhost_source[,1])
# these are the predicted log-odds for each population
# correct each observation for the host source and parasite source so we are plotting residuals
exp.dat.psre$residProp <- NA

for(i in 1:nrow(exp.dat.psre)){ # for each individual observation
  host.re <- r.e.host$RandomEffect[which(r.e.host$HostPopulation == as.character(exp.dat.psre$fhost_source[i]))] # extract the random effect for the host pop
  parasite.re <- r.e.parasite$RandomEffect[which(as.character(r.e.parasite$ParasitePopulation) == as.character(exp.dat.psre$fparas_source[i]))] # extract the random effect for the parasite pop
  adj <- host.re + parasite.re # add together the two re and exponentiate
  exp.dat.psre$residProp[i] <- plogis(qlogis(exp.dat.psre$Paras_Inf[i]/exp.dat.psre$Dosage[i]) - adj)
}
# RIB
polygon.x.r <- predict.data %>% filter(Parasite == "Rib") %>% pull(dist_host_paras) 
ci.lower.plot.r <- predict.data %>% filter(Parasite == "Rib") %>% pull(CI.lower)
ci.upper.plot.r <- predict.data %>% filter(Parasite == "Rib") %>% pull(CI.upper) 
polygon.y.r <- c(ci.lower.plot.r, rev(ci.upper.plot.r))

plot(residProp ~ dist_host_paras, data =  subset(exp.dat.psre, Parasite == "Rib"), pch = 16, col = "lightblue", xlab = "Geographic Distance (km)", ylab = "Parasite Infection Success", ylim = c(0,1), type = "n")
polygon(x=c(polygon.x.r, rev(polygon.x.r)), y=polygon.y.r, col='lightblue', border=NA)
# add the best fit line
points(predpar~dist_host_paras, subset(predict.data, Parasite == "Rib"), type = "l", lty = 1, col = "lightblue4", lwd = 2)
# add the lower CI
points(CI.lower~dist_host_paras, subset(predict.data, Parasite == "Rib"), type = "l", lty = 2, col = "lightblue4")
# add the upper CI
points(CI.upper~dist_host_paras, subset(predict.data, Parasite == "Rib"), type = "l", lty = 2, col = "lightblue4")
# add raw data  
points(residProp ~ dist_host_paras, data =  subset(exp.dat.psre, Parasite == "Rib"), pch = 16, col = "lightblue4")

# MANO
polygon.x.m <- predict.data %>% filter(Parasite == "Mano") %>% pull(dist_host_paras) 
ci.lower.plot.m <- predict.data %>% filter(Parasite == "Mano") %>% pull(CI.lower)
ci.upper.plot.m <- predict.data %>% filter(Parasite == "Mano") %>% pull(CI.upper) 
polygon.y.m <- c(ci.lower.plot.m, rev(ci.upper.plot.m))
plot(residProp ~ dist_host_paras, data =  subset(exp.dat.psre, Parasite == "Mano"), pch = 16, ylim = c(0,1), type = "n", col = "indianred", xlab = "Geographic Distance (km)", ylab = "Parasite Infection Success")
polygon(x=c(polygon.x.m, rev(polygon.x.m)), y=polygon.y.m, col='indianred1', border=NA)
# add the best fit line
points(predpar~dist_host_paras, subset(predict.data, Parasite == "Mano"), type = "l", lty = 1, col = "indianred4", lwd = 2)
# add the lower CI
points(CI.lower~dist_host_paras, subset(predict.data, Parasite == "Mano"), type = "l", lty = 2, col = "indianred4")
# add the upper CI
points(CI.upper~dist_host_paras, subset(predict.data, Parasite == "Mano"), type = "l", lty = 2, col = "indianred4")
points(residProp ~ dist_host_paras, data =  subset(exp.dat.psre, Parasite == "Mano"), pch = 16, col = "indianred4")


################################################
####### Genetic data ###########################
################################################
pg <- read.csv("data/Parasite_Genetic_Dist.csv") # Parasite genetic distance
NeiD <- read.csv("data/neisD_frog.csv") # PSRE genetic distance (Nei's D)
FST <- read.csv("data/MeanFST_frog.csv") # PSRE genetic distance (FST)
GeogDist <- read.csv("data/Geog_Dist.csv") # Geographic Distances
pcadat <- read.csv("data/PSRE_PCA_scores.csv") # PSRE PCA scores from genetic clustering

# format as matrix
GeogDist %>% column_to_rownames("X") %>% as.matrix() -> GeogDist
NeiD %>% column_to_rownames("X") %>% as.matrix() -> NeiD
FST %>% column_to_rownames("X") %>% as.matrix() -> FST

# make sure dimensions match and they're sorted same way 
GeogDist <- subset(GeogDist, rownames(GeogDist) %in% rownames(NeiD))
GeogDist <- GeogDist[, colnames(GeogDist) %in% colnames(NeiD)]
GeogDist <- GeogDist[sort(rownames(GeogDist)),sort(colnames(GeogDist))]
NeiD <- NeiD[sort(rownames(NeiD)), sort(colnames(NeiD))]
FST <- FST[sort(rownames(FST)), sort(colnames(FST))]

# mantel test: Host Dissimilarity vs. Geographic Distance
mantel.test(m1 = GeogDist, m2 =  NeiD) 
mantel.test(m1 = GeogDist, m2 =  FST) 

# create df for plotting
# remove redundancies in matrix by getting rid of lower triangle
GeogDist[lower.tri(GeogDist)] <- NA
GeoNei <- melt(GeogDist)
NeiD[lower.tri(NeiD)] <- NA
GeoNei <- cbind(GeoNei, melt(NeiD))
GeoNei <- GeoNei[, c(1,2,3,6)]
colnames(GeoNei)  <- c("Pop1", "Pop2", "GeoD", "NeiD")

FST[lower.tri(FST)] <- NA
GeoFST <- cbind(melt(GeogDist), melt(FST))
GeoFST <- GeoFST[, c(1,2,3,6)]
colnames(GeoFST)  <- c("Pop1", "Pop2", "GeoD", "MeanFST")

##### Fig 4 (Genetic dissimilarity) ###########
# hosts isolation by distance
plot(NeiD~GeoD, GeoNei, pch = 16, xlab = "Geographic distance (km)", ylab = "Genetic dissimilarity (Nei's D)", col = "#568D3F", cex.axis = .8, cex = 1.5)
plot(MeanFST~GeoD, GeoFST, pch = 16, xlab = "Geographic distance (km)", ylab = expression('Genetic dissimilarity (F'[ST]*')'), col = "#568D3F", cex.axis = .8, cex = 1.5)
# host genetic clustering
cbbPalette <- c("darkgoldenrod1", "firebrick1", "darkmagenta", "dodgerblue4", "darkorange", "lightslateblue", "mediumturquoise")
pcameans <- pcadat %>% group_by(Label) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2)) %>% data.frame()
palette(cbbPalette)

plot(PC2~PC1, pcadat, col = as.factor(Label), pch = 1, cex = 1, xlab = "PC1 (25% variance explained)",
     ylab = "PC2 (13% variance explained)" )
points(PC2~PC1, pcameans, bg = as.factor(Label), pch = 21, cex = 1.5, ylim = c(-8,13), xlim = c(-7,16))
legend(legend = unique(pcadat$Label), x = 0, y = 12, pt.bg = cbbPalette, pch = 21,  ncol = 2, cex = .8, bty = "n")

# parasite isolation by distance
plot(meanD~jitter(GeogDist, .1), col = "lightblue4", subset(pg, Parasite=="Rib"), ylab = "Sequence dissimilarity (% bp)", xlab = "Geographic Distance (km)", pch = 16, xlim = c(-10, 1000), ylim = c(-.01, .15), cex.axis = 0.8, cex = 1.5)
points(meanD~jitter(GeogDist, .1), col = "indianred", subset(pg, Parasite=="Mano"), pch = 16, cex = 1.5)



################################################
####### Other host species analysis ###########
################################################
##### Run GLMM ######
allspprib <- read.csv("data/Experimental_Crosses_AllSpp.csv")
allspprib$Species <- as.factor(allspprib$Species)
allspprib$Species <- relevel(allspprib$Species, ref = "PSRE")

allspglm <- glmer(cbind(Paras_Inf, Paras_Uninf)~(1|fhost_source)+(1|fparas_source)+scale(SVL)+
                    scale(dist_host_paras)*Species, data = allspprib, family = binomial)
summary(allspglm)

# no distance
allspglm2 <- glmer(cbind(Paras_Inf, Paras_Uninf)~(1|fhost_source)+(1|fparas_source)+scale(SVL)+
                    Species, data = allspprib, family = binomial)

summary(allspglm2)

# effect of distance (not significant)
anova(allspglm, allspglm2)

summary(glht(allspglm, mcp(Species="Tukey")))
summary(glht(allspglm2, mcp(Species="Tukey")))

predict.species <- expand.grid(Species = c("BUBO", "PSRE", "TATO"), fhost_source = "fake", fparas_source = "fake",
                        SVL = mean(allspprib$SVL, na.rm=TRUE), dist_host_paras = mean(allspprib$dist_host_paras, na.rm=TRUE))

predict(allspglm2, newdata = predict.species, type = "response", allow.new.levels=TRUE)
predict.species$Species
(0.67-0.56)/0.56
(0.67-0.51)/0.51

