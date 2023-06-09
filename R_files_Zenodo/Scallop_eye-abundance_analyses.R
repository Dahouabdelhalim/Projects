### Variation in eye abundance among scallops reveals ontogenetic
### and evolutionary convergence associated with life habits

# Jorge A Audino, Dean C Adams, and Jeanne M Serb

library(ape)
library(phytools)
library(geiger)
library(ggplot2)
library(RRPP)
library(geomorph)
library(OUwie)
library(dplyr)
library(hrbrthemes)
library(viridis)

### 1 ###
### DOES EYE COUNT VARY AMONG SPECIES? ###

# Read table for 37 species (total of 324 observations)
all.data <- read.table("Scallop_eyecount_shell_height.txt", sep = "\\t", header = T)
eye.count <- log(all.data[,2]) # 
Species <- as.factor(all.data[,1])
### ANOVA via RRPP
LM <- lm.rrpp(eye.count~Species, iter = 9999, print.progress = FALSE)
anova(LM)$table
### Supplementary material ###
# Visual inspection of the linear model
pdf("Supplementary_FigS1a.pdf", width = 5, height = 4)
par(mfrow=c(1,1))
hist(LM$LM$residuals, xlab = "Linear model residuals", main = "ANOVA eye count ~ species",
     cex.lab = 1.2, cex.main = 1.4)
dev.off()
pdf("Supplementary_FigS1b.pdf", width = 7, height = 7)
par(mfrow=c(2,2))
plot(LM, cex = 1.5)
dev.off()

# Boxplot of eye abundance per species
eyes.species <- data.frame(all.data[,2], Species)
names(eyes.species) <- c("Eyes","Species")
eyes.species <- eyes.species[c("Species", "Eyes")]
eyes.species.plot <- ggplot(eyes.species, aes(x=Species, y=Eyes)) +
  geom_boxplot(outlier.shape=21)
eyes.species.plot.full <- eyes.species.plot +
  theme_gray(base_size = 11) +
  xlab("Species of scallops") + 
  theme(axis.title.x=element_text(angle=0, size=16)) +
  theme(axis.title.x=element_text(vjust=-0.50)) +
  theme(axis.text.x=element_text(angle=90, size=7, face = "italic")) +
  ylab("Eye abundance") +
  theme(axis.title.y=element_text(angle=90, size=16)) +
  theme(axis.title.y=element_text(vjust=1)) +
  theme(axis.text.y=element_text(angle=0, size=10)) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  guides(fill=FALSE)
eyes.species.plot.full
pdf("Supplementary_Fig4.pdf", width = 7, height = 6)
eyes.species.plot.full
dev.off()

### 2 ###
### ALLOMETRY ###
## Testing whether eye count varies with shell height (~ body size)
# Read data for 5 species
allometry.data <- read.table("Scallop_allometry.txt", sep = "\\t", header = T)
dat <- log(allometry.data[,1:2])
sp <- as.factor(allometry.data[,3])
## Fit Linear Models (LM) for regression analysis
# Aequipecten glyptus
Aglyptus <- log(allometry.data[(1:18),(1:2)])
Aglyptus.model <- lm.rrpp(Aglyptus$eyecount~Aglyptus$height, iter = 9999)
Aglyptus.coef <- coef(Aglyptus.model)
anova(Aglyptus.model)$table
# Amusium pleuronectes
Apleuro <- log(allometry.data[(19:54),(1:2)])
Apleuro.model <- lm.rrpp(Apleuro$eyecount ~Apleuro$height, iter = 9999)
Apleuro.coef <- coef(Apleuro.model)
anova(Apleuro.model)$table
# Chlamys behringiana
Cbehr <- log(allometry.data[(55:74),(1:2)])
Cbehr.model <- lm.rrpp(Cbehr$eyecount~Cbehr$height, iter = 9999)
Cbehr.coef <- coef(Cbehr.model)
anova(Cbehr.model)$table
# Ylistrum balloti
Ybal <- log(allometry.data[(75:119),(1:2)])
Ybal.model <- lm.rrpp(Ybal$eyecount~Ybal$height, iter = 9999)
Ybal.coef <- coef(Ybal.model)
anova(Ybal.model)$table
# Placopecten magellanicus
Pmag <- log(allometry.data[(120:167),(1:2)])
Pmag.model <- lm.rrpp(Pmag$eyecount~Pmag$height, iter = 9999)
Pmag.coef <- coef(Pmag.model)
anova(Pmag.model)$table

# Plot regression lines for each species
pdf("Fig2a.pdf", width = 7, height = 6)
par(mar = c(5, 5, 4, 3) + 0.1)
plot(dat$height, dat$eyecount, pch = 21, bg = as.numeric(sp),
     xlab = "log shell height", ylab = "log eye abundance", cex.lab = 1.6,
     main = "Static allometry", cex.main = 1.5, cex = 1.5)
abline(Aglyptus.coef, col = "black", lwd = 2)
abline(Apleuro.coef, col = "#DF536B", lwd = 2)
abline(Cbehr.coef, col = "#61D04F", lwd = 2)
abline(Ybal.coef, col = "#28E2E5", lwd = 2)
abline(Pmag.coef, col = "#2297E6", lwd = 2)
legend("bottomright", c("Aequipecten glyptus","Amusium pleuronectes",
                        "Chlamys behringiana","Ylistrum balloti",
                        "Placopecten magellanicus"),
       pch = 21, pt.bg = c("black", "#DF536B", "#61D04F",
                           "#28E2E5", "#2297E6"), cex = 1)
dev.off()

### ANCOVA
# Testing differences in allometry among species
ancova.model <- lm.rrpp(dat$eyecount~dat$height*sp, iter = 9999)
common.slope.model <- lm.rrpp(dat$eyecount~dat$height+sp, iter = 9999)
anova(ancova.model)$table # interaction term significant
#
PW <- pairwise(ancova.model, fit.null = common.slope.model,
               groups = sp, covariate = dat$height)
summary(PW, test.type = "DL", show.vectors = T)
# differences between vector lengths are significant

### Supplementary material ###
# Visual inspection of the linear model
pdf("Supplementary_FigS2a.pdf", width = 5, height = 4)
par(mfrow=c(1,1))
hist(ancova.model$LM$residuals, xlab = "Linear model residuals", main = "ANCOVA eye count~shell height*species",
     cex.lab = 1.2, cex.main = 1.2)
dev.off()
pdf("Supplementary_FigS2b.pdf", width = 7, height = 7)
par(mfrow=c(2,2))
plot(ancova.model, cex = 1.5)
dev.off()


### 3 ###
### ASSOCIATIONS WHILE ACCOUNTING FOR PHYLOGENETIC RELATIONSHIPS ###
# Read data for 33 species species of Pectinidae
mydata <- read.table("eyecount_means.txt", sep = "\\t", header = T)
species <- mydata[,1]
log.eyes <- log(mydata[,2]) # means of eye abundance
log.height <- log(mydata[,3]) # means of shell heights
lifestyle <- as.factor(mydata[,4]) # habits of life
# Organize data and tree for 31 species (no shell heights for A. purpuratus and B. vexillum)
keep <- complete.cases(log.eyes, log.height, lifestyle)
scallops <- data.frame(log.eyes = log.eyes[keep], log.height = log.height[keep],
                       lifestyle = lifestyle[keep])
rownames(scallops) <- mydata$species[keep]
# Read phylogenetic MCC tree
scallop.tree <- read.nexus(file = "Scallop_MCC.tre")
taxa <- name.check(scallop.tree, scallops)
pruned.tree <- drop.tip(scallop.tree, tip=taxa$tree_not_data)
name.check(pruned.tree, scallops) # matching tips labels and dataset
# Phylogenetic covariance Matrix
phycov <- vcv(pruned.tree)

## Phylogenetic ANCOVA via Generalized Least Squares (GLS)
# Eye abundance ~ shell height * habits of life
GLS1 <- lm.rrpp(log.eyes~log.height*lifestyle, data = scallops,
                Cov = phycov, iter = 9999)
anova(GLS1)$table
# Interaction term not significant (= similar slopes among groups)
## Eye abundance ~ shell height + habits of life (= common slopes model)
GLS2 <- lm.rrpp(log.eyes~log.height+lifestyle,data = scallops,
                Cov = phycov, iter = 9999)
anova(GLS2)$table
# Both variables significantly explain variation in eye abundance
# Body size contributes to explain 21% of variation in eye abundance
# Lifestyle contributes to explain 52% of variation in eye abundance

### Supplementary material ###
# Visual inspection of the linear model
pdf("Supplementary_FigS3a.pdf", width = 5, height = 4)
par(mfrow=c(1,1))
hist(GLS2$LM$gls.residuals, xlab = "Linear model residuals", main = "Phylogenetic ANCOVA",
     cex.lab = 1.2, cex.main = 1.2)
dev.off()
pdf("Supplementary_FigS3b.pdf", width = 7, height = 7)
par(mfrow=c(2,2))
plot(GLS2, cex = 1.5)
dev.off()

# Plot regression for evolutionary allometry
pdf("Fig2b.pdf", width = 7, height = 6)
par(mar = c(5, 5, 4, 3) + 0.1)
plot(log.eyes~log.height, pch = 19, xlab = "log shell height",
     ylab = "log eye abundance", main = "Evolutionary allometry",
     cex.main = 1.5, cex.lab = 1.6, cex = 1.5)
abline(coef(lm(log.eyes~log.height)), col = "black", lwd = 2)
dev.off()

# Pairwise comparisons between habits of life
pw <- pairwise(GLS2, groups = scallops$lifestyle)
test.pw <- summary(pw, test.type = "dist"); test.pw
# significant comparisons: byssal-attach:cement, byssal-attach:recess,
# cement:free-living, cement:glide, glide:recess

# Plot eye abundance per habit of life
eyes.lifestyle <- data.frame(mydata[,2], lifestyle)
eyes.lifestyle <- eyes.lifestyle[c("lifestyle", "mydata...2.")]
eyes.lifestyle.plot <- ggplot(eyes.lifestyle, aes(x=lifestyle, y=mydata...2., fill = lifestyle)) +
  geom_boxplot()
eyes.lifestyle.plot.full <- eyes.lifestyle.plot +
  scale_fill_manual(values=c("#99FF99", "#FF6633", "#FFFF66", "#666FFF", "#CC99FF")) +
  theme_gray(base_size = 11) +
  xlab("Life habits") + 
  theme(axis.title.x=element_text(angle=0, size=18)) +
  theme(axis.title.x=element_text(vjust=-1)) +
  theme(axis.text.x=element_text(angle=0, size=14)) +
  ylab("Eye abundance per species") +
  theme(axis.title.y=element_text(angle=90, size=18)) +
  theme(axis.title.y=element_text(vjust=1)) +
  theme(axis.text.y=element_text(angle=0, size=12)) +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  guides(fill=FALSE)
eyes.lifestyle.plot.full
pdf("Fig1b.pdf", width = 7, height = 5)
par(mar = c(6, 5, 2, 2) + 0.1)
eyes.lifestyle.plot.full
dev.off()

### 4 ###
### TRAIT EVOLUTION ###
eyes <- mydata[,2] # eye count (means per species)
eyes <- setNames(eyes, mydata[,1]) # set species names to means
lifestyles <- mydata[,4] # habits of life
lifestyles <-setNames(lifestyles, mydata[,1])
data.taxa <- as.matrix(mydata[,2])
row.names(data.taxa) <- mydata[,1]
taxa2 <- name.check(scallop.tree, data.taxa)
tree <- drop.tip(scallop.tree, tip = taxa2$tree_not_data)
name.check(tree, data.taxa)

### 4.1. CONTINUOUS TRAIT EVOLUTION: Eye abundance
# Estimate ancestral values and plot mapped tree
eye.contMap <- contMap(tree, eyes, plot = FALSE)
colfunc <- colorRampPalette(c("black", "white", "#666FFF"))
eye.contMap <- setMap(eye.contMap, colors = c(colfunc(5)))
pdf("Fig3a.pdf", width = 3, height = 4)
plot(eye.contMap, fsize = c(0.6,0.6),
     leg.txt = "Eye abundance", ftype = "off")
dev.off()

### 4.3. DISCRETE TRAIT EVOLUTION: scallop habits of life
## Fit data to the phylogeny under three evolutionary models:
# equal rates (ER), symmetric (SYM), all rates different (ARD)
discrete.models <- list("ER", "SYM", "ARD")
models.output <- lapply(discrete.models, function(x) {
  fitDiscrete(tree, lifestyles, model = x)
}) 
# Identify the optimal model of discrete trait evolution
AICs <- lapply(1:length(models.output), function(x) {models.output[[x]]$opt$aic}) # AIC values
AICs # Best fit with "ER" model

## Estimate ancestral states under the optimal model "ER"
# (1) Maximum likelihood framework
fitER <- ace(lifestyles, tree, model = "ER", type = "discrete")
pdf("Supplementary_Fig8.pdf", width = 4, height = 4)
plotTree(tree, fsize = 0.7, ftype = "i")
cols <- setNames(c("#99FF99", "#FF6633", "#FFFF66", "#666FFF", "#CC99FF"), sort(unique(lifestyles)))
habit.mode<-as.factor(setNames(mydata[,4],mydata[,1]))
tiplabels(pie=to.matrix(habit.mode[tree$tip.label], sort(unique(lifestyles))),
          piecol = cols,cex=0.45)
nodelabels(node = 1:tree$Nnode+Ntip(tree),
           pie = fitER$lik.anc, piecol = cols, cex = 0.7)
add.simmap.legend(colors = cols, fsize = 0.8,
                  prompt = FALSE, x = 1, y = 9)
dev.off()
# (2) Bayesian stochastic character mapping (SIMMAP)
# accounting for ancestral state uncertainty
# accounting for phylogenetic uncertainty
posterior.trees <- read.nexus(file = "Scallop_posterior_trees.trees") #10,001 trees
sampled.trees <- sample(posterior.trees, size = 100) # sampling 100 trees
trees <- drop.tip.multiPhylo(sampled.trees, tip = "Spondylus_nicobaricus")
# Perform 10 simulations for ancestral states in each of 100 sampled trees
simmap.habits <- make.simmap(trees, lifestyles, model = "ER", nsim = 10, message = F)
sum.simmap.habits <- summary(simmap.habits, check.equal = T, message = F)
# plot(sum.simmap.habits, colors = cols, fsize = 0.6, ftype = "i", direction="leftwards")
# rotate tree to match contMap tree
pdf("Fig3b.pdf", width = 3, height = 4)
sum.simmap.habits$ref.tree <- rotateConstr(sum.simmap.habits$ref.tree, tree$tip.label)
plot(sum.simmap.habits, fsize = 0.6, ftype = "i", direction="leftwards",
     colors = cols, cex = 0.65)
add.simmap.legend(colors = cols, fsize = 0.7,
                  prompt = FALSE, x = 0.45, y = 6)
dev.off()

### 4.4. COMPARE EVOLUTIONARY MODELS OF EVOLUTION USING OUwie
# Trait: eye abundance
# Selective regimes: lifestyles
OUwie.data <- data.frame(mydata[,1], mydata[,4], mydata[,2])
trees.OUwie <- sample(posterior.trees, size = 100) # sample 100 trees
trees.OUwie <- drop.tip.multiPhylo(trees.OUwie, tip = "Spondylus_nicobaricus")
# Use 1000 simmap trees for analysis (10 simulations on each sampled tree)
simmap.OUwie <- make.simmap(trees.OUwie, lifestyle, model = "ER", nsim = 10, message = F)
# Single-rate Brownian motion (BM1)
# Brownian motion with different rate parameters for each state (BMS)
# Ornstein-Uhlenbeck model with a single optimum for all species (OU1)
# Ornstein-Uhlenbeck model with different state means
# and single alpha and sigma^2 acting all selective regimes (OUM)

# Fitting BM, BMS, OU1 & OUM for each of 1,000 SIMMAP trees
OUwie.run <- vector("list", 1000)
for (i in 1:1000){OUwie.run[[i]] <- vector("list", 4)}
OUwie.models <- list("BM1", "BMS", "OU1", "OUM")
for (i in 1:1000){
  OUwie.run[[i]] <- lapply(1:length(OUwie.models), function(x) {
    OUwie(simmap.OUwie[[i]], OUwie.data, model = OUwie.models[[x]], ub = 500,
          algorithm = "invert", simmap.tree = TRUE, diagn = T)})
}
#Extract key information from model results
BM1 <- matrix(NA, nrow = 1000, ncol = 4)
BMS <- matrix(NA, nrow = 1000, ncol = 12)
OU1 <- matrix(NA, nrow = 1000, ncol = 6)
OUM <- matrix(NA, nrow = 1000, ncol = 10)
colnames(BM1) <- c("AICc", "eigval", "sigmasq", "theta")
colnames(BMS) <- c("AICc", "eigval1", "eigval2","eigval3", "eigval4", 
                   "eigval5", "sigmasq1", "sigmasq2", "sigmasq3",
                   "sigmasq4", "sigmasq5", "theta")
colnames(OU1) <- c("AICc","eigval1","eigval2","alpha","sigmasq","theta")
colnames(OUM) <- c("AICc", "eigval1", "eigval2","alpha","sigmasq",
                   "theta1","theta2","theta3","theta4","theta5" )
# Organize that information in matrices for easier comparisons
for (i in 1:1000) {
  BM1[,1][[i]] <- OUwie.run[[i]][[1]]$AICc #AICc BM1
  BM1[,2][[i]] <- OUwie.run[[i]][[1]]$eigval # eigenvalues
  BM1[,3][[i]] <- OUwie.run[[i]][[1]]$solution[2]# sigma.sq
  BM1[,4][[i]] <- OUwie.run[[i]][[1]]$theta[1,1]# optimum (theta)
  BMS[,1][[i]] <- OUwie.run[[i]][[2]]$AICc #AICc BMS
  BMS[,2][[i]] <- OUwie.run[[i]][[2]]$eigval[1] # eigenvalues
  BMS[,3][[i]] <- OUwie.run[[i]][[2]]$eigval[2]
  BMS[,4][[i]] <- OUwie.run[[i]][[2]]$eigval[3]
  BMS[,5][[i]] <- OUwie.run[[i]][[2]]$eigval[4]
  BMS[,6][[i]] <- OUwie.run[[i]][[2]]$eigval[5]
  BMS[,7][[i]] <- OUwie.run[[i]][[2]]$solution[2]# sigma.sq byssal-attach
  BMS[,8][[i]] <- OUwie.run[[i]][[2]]$solution[4]# sigma.sq cement
  BMS[,9][[i]] <- OUwie.run[[i]][[2]]$solution[6]# sigma.sq free-living
  BMS[,10][[i]] <- OUwie.run[[i]][[2]]$solution[8]# sigma.sq glide
  BMS[,11][[i]] <- OUwie.run[[i]][[2]]$solution[10]# sigma.sq recess
  BMS[,12][[i]] <- OUwie.run[[i]][[2]]$theta[1,1]# optimum (theta)
  OU1[,1][[i]] <- OUwie.run[[i]][[3]]$AICc # AICc OU1
  OU1[,2][[i]] <- OUwie.run[[i]][[3]]$eigval[1] # eigenvalues
  OU1[,3][[i]] <- OUwie.run[[i]][[3]]$eigval[2]
  OU1[,4][[i]] <- OUwie.run[[i]][[3]]$solution[1] # alpha
  OU1[,5][[i]] <- OUwie.run[[i]][[3]]$solution[2] # sigma.sq
  OU1[,6][[i]] <- OUwie.run[[i]][[3]]$theta[1] # optimum (theta)
  OUM[,1][[i]] <- OUwie.run[[i]][[4]]$AICc # AICc OUM
  OUM[,2][[i]] <- OUwie.run[[i]][[4]]$eigval[1] # eigenvalues
  OUM[,3][[i]] <- OUwie.run[[i]][[4]]$eigval[2]
  OUM[,4][[i]] <- OUwie.run[[i]][[4]]$solution[1] # alpha
  OUM[,5][[i]] <- OUwie.run[[i]][[4]]$solution[2] # sigma.sq
  OUM[,6][[i]] <- OUwie.run[[i]][[4]]$theta[1,1]# theta byssal-attach
  OUM[,7][[i]] <- OUwie.run[[i]][[4]]$theta[2,1]# theta cement
  OUM[,8][[i]] <- OUwie.run[[i]][[4]]$theta[3,1]# theta gree-living
  OUM[,9][[i]] <- OUwie.run[[i]][[4]]$theta[4,1]# theta glide
  OUM[,10][[i]] <- OUwie.run[[i]][[4]]$theta[5,1]# theta recess
}

### Inspection of model parameters and estimates
# Eigenvalues produced by each model
eigvals <- cbind(BM1[,2],BMS[,2],BMS[,3],BMS[,4], BMS[,5], BMS[,6],
                 OU1[,2],OU1[,3],OUM[,2],OUM[,3])
colnames(eigvals) <- c("BM1","BMS","BMS","BMS","BMS","BMS","OU1","OU1","OUM","OUM")
# Negative eigenvalues: model failed in obtaining reliable parameter estimates
eigvals.failed <- which(eigvals < 0, arr.ind = TRUE)
eigvals.clean <- eigvals[-c(eigvals.failed[,1]),] # discarding models

# Theta: optima (eye abundance)
theta <- cbind(BM1[,4],BMS[,12],OU1[,6],OUM[,6],OUM[,7],OUM[,8],OUM[,9],OUM[,10])
colnames(theta) <- c("BM1","BMS","OU1","OUM","OUM","OUM","OUM","OUM")
# Check for unrealistic optima (theta < 0)
theta.failed <- which(theta < 0, arr.ind = TRUE)
theta.clean1 <- theta[-c(theta.failed[,1]),] # discarding models
theta.failed2 <- which(theta.clean1 > 200, arr.ind = TRUE)
theta.clean <- theta.clean1[-c(theta.failed2[,1]),]
theta.means <- apply(theta.clean, 2, mean)
# Data frame for Violin plot
theta.df <- data.frame(
  model = c( rep("BM1",907), rep("BMS",907),rep("OU1",907),
             rep("OUM1",907),rep("OUM2",907), rep("OUM3",907),
             rep("OUM4",907), rep("OUM5",907)),
  theta = c( theta.clean[,1], theta.clean[,2], theta.clean[,3], theta.clean[,4],
                theta.clean[,5], theta.clean[,6], theta.clean[,7], theta.clean[,8]))
# Sample size and summary statistics for boxplot
sample_theta = theta.df %>% group_by(model) %>% summarize(num=n())
# Plot Violin with Boxplot
pdf("Fig4d.pdf", width = 5, height = 5)
par(mar = c(5,5,3,3))
theta.df %>%
  left_join(sample_theta) %>%
  mutate(myaxis = paste0(model)) %>%
  ggplot( aes(x=model, y=theta, fill=model)) +
  scale_fill_manual(values=c(c("#666666", "#999999", "#CCCCCC"), c(cols))) +
  geom_violin(bw=2) +
  geom_boxplot(width=0.5, color="grey", alpha=0.1) +
  theme_gray(base_size = 11) +
  guides(fill=FALSE)
dev.off()

# Sigma.sq: rate of evolution under different models
sigma.sq <- cbind(BM1[,3],BMS[,7],BMS[,8],BMS[,9], BMS[,10], BMS[,11],
                  OU1[,5],OUM[,5])
colnames(sigma.sq) <- c("BM1","BMS","BMS","BMS","BMS","BMS","OU1","OUM")
sigma.sq.means <- apply(sigma.sq, 2, mean); sigma.sq.means
# Data frame for Violin plot
sigma.df <- data.frame(
  model = c( rep("BM1",1000), rep("BMS1",1000),rep("BMS2",1000),
             rep("BMS3",1000),rep("BMS4",1000), rep("BMS5",1000),
             rep("OU1",1000), rep("OUM",1000)),
  sigma.sq = c( sigma.sq[,1], sigma.sq[,2], sigma.sq[,3], sigma.sq[,4],
                sigma.sq[,5], sigma.sq[,6], sigma.sq[,7], sigma.sq[,8])
)
# Sample size and summary statistics for boxplot
sample_sigma = sigma.df %>% group_by(model) %>% summarize(num=n())
# Plot Violin with Boxplot
pdf("Fig4c.pdf", width = 5, height = 5)
par(mar = c(5,5,3,3))
sigma.df %>%
  left_join(sample_sigma) %>%
  mutate(myaxis = paste0(model)) %>%
  ggplot( aes(x=model, y=sigma.sq, fill=model)) +
  scale_fill_manual(values=c("#666666", cols, "#CCCCCC", "#FFFFFF")) +
  geom_violin() +
  geom_boxplot(width=0.5, color="grey", alpha=0.1) +
  theme_gray(base_size = 11) +
  guides(fill=FALSE)
dev.off()

# alpha: selective force (estimated from dataset)
alpha <- cbind(OU1[,4],OUM[,4])
colnames(alpha) <- c("OU1","OUM")
alpha.means <- apply(alpha, 2, mean); alpha.means
pdf("Supplementary_FigS5.pdf", width = 5, height = 5)
par(mar = c(5,5,3,3))
boxplot(alpha, ylab = expression(paste("Selective force (",alpha,")")),
        outline = FALSE, xlab = "Ornstein-Uhlenbeck models", cex.lab = 1.3,
        col = c("#CCCCCC", "#FFFFFF"))
dev.off()

### Identify the optimal model
# Organize AICc scores for each reconstruction
AICc <- cbind(BM1[,1],BMS[,1],OU1[,1],OUM[,1])
colnames(AICc) <- c("BM1", "BMS", "OU1", "OUM")
# discard reconstructions with failed parameters (eigval & theta)
to.remove <- c(theta.failed[,1], eigvals.failed[,1])
AICc.clean <- AICc[-c(to.remove),] # 874 reconstructions
AICc.means <- apply(AICc.clean, 2, mean); AICc.means
# lower scores, better fit
# Plot performance of each evolutionary model based on AICc scores
# Data frame for Violin plot
AIC.df <- data.frame(
  model = c( rep("BM1",874), rep("BMS",874), rep("OU1",874), rep("OUM",874)),
  AICc.score = c( AICc.clean[,1], AICc.clean[,2], AICc.clean[,3], AICc.clean[,4])
)
# Sample size and summary statistics for boxplot
sample_AIC = AIC.df %>% group_by(model) %>% summarize(num=n())
# Plot Violin with Boxplot
pdf("Fig4a.pdf", width = 5, height = 5)
par(mar = c(5,5,3,3))
AIC.df %>%
  left_join(sample_AIC) %>%
  mutate(myaxis = paste0(model)) %>%
  ggplot( aes(x=model, y=AICc.score, fill=model)) +
  scale_fill_manual(values=c("#666666", "#999999", "#CCCCCC", "#FFFFFF")) +
  geom_violin() +
  geom_boxplot(width=0.5, color="grey", alpha=0.1) +
  theme_gray(base_size = 11) +
  guides(fill=FALSE)
dev.off()

# Select the optimal model
optimal.models <- lapply(1:length(AICc.clean[,1]), function(x) {sort(AICc.clean[x,])})
delta <- lapply(1:length(optimal.models), function(x) {optimal.models[[x]][2]-optimal.models[[x]][1]})
keep <- which(delta > 2, arr.ind = TRUE) # keeping delta AIcc > 2
significant.dif <- optimal.models[keep] 
best.fit <- lapply(1:length(significant.dif), function(x) {significant.dif[[x]][1]})
best.fit
# 556 out of 874 reconstructions have delta AICc > 2
length(grep("OUM", best.fit)) # optimal performance of OUM 426 times
length(grep("OU1", best.fit)) # optimal performance of OU1 130 times

# AICc weights
AIC.weights <- lapply(1:dim(AICc.clean)[1], function(x) {aicw(t(AICc.clean)[,x])})
AICw <- matrix(NA, nrow = length(AIC.weights), ncol = 4)
colnames(AICw) <- c("BM1", "BMS", "OU1", "OUM")
for (i in 1:length(AIC.weights)){
  AICw[i,] <- AIC.weights[[i]][,3]
}
AICw.means <- apply(AICw, 2, mean); AICw.means
#plot AICc weights
#Data frame for Violin plot
AICw.df <- data.frame(
  model = c( rep("BM1",874), rep("BMS",874), rep("OU1",874), rep("OUM",874)),
  AIC.weights = c( AICw[,1], AICw[,2], AICw[,3], AICw[,4])
)
# Sample size and summary statistics for boxplot
sample_AICw = AICw.df %>% group_by(model) %>% summarize(num=n())
# Plot Violin with Boxplot
pdf("Fig4b.pdf", width = 5, height = 5)
par(mar = c(5,5,3,3))
AICw.df %>%
  left_join(sample_AICw) %>%
  mutate(myaxis = paste0(model)) %>%
  ggplot( aes(x=model, y=AIC.weights, fill=model)) +
  scale_fill_manual(values=c("#666666", "#999999", "#CCCCCC", "#FFFFFF")) +
  geom_violin(bw=0.04) +
  geom_boxplot(width=0.5, color="grey", alpha=0.1) +
  theme_gray(base_size = 11) +
  guides(fill=FALSE)
dev.off()

sorted.weights <- lapply(1:length(AIC.weights), function(x) {sort.list(AIC.weights[[x]][,3], decreasing = T)})
largest.weights <- lapply(1:length(sorted.weights), function(x) {sorted.weights[[x]][1]})
length(grep("3", largest.weights)) # OU1 had greater support 265 times
length(grep("4", largest.weights)) # OUM had greater support 609 times
# Overall, the optimal model is OUM based on AICc and Akaike weights