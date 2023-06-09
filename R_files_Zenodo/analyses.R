#Final (hopefully!) analyses for manuscript with full annotation
# =============================================================================
# 0) Import data.
# 1) Methodological exploration
  # 1a) How well do the traits (normalized as ratios) predict function?
  # 1b) How well does PC1 and PC2 (and rest of PCs as an aside) on raw data predict function?
  # 1c) How well does PC1 of size corrected data predict function?
# 2) Parallelism exploration
  # 2a) Plots of parallelism in 3 functional traits
  # 2b) Effect size of habitat term: parallelism in 3 functional measures
  # 2c) Effect size of habitat term: parallelism in PC2 of the component traits for 3 functional measures
# 3) Background stats to round out paper
  # 3a) Show that component morphology and function varies significantly among populations using MANOVA
  # 3b) Show that there are within pair habitat effects.
  # 3c) Revision1: histogram of eta-sq values for component traits subsumed in histogram of all traits, to show that divergence likely under selection
  # 3d) Revision1: P-matrices of LR, SI, and KT component traits to estimate dimensionality
#=============================================================================

# 0) Import data.
# =============================================================================
#RUN THIS AFTER RUNNING main.cole.R
cleandatafinal <- cleandata.sizecorr.PCs
names(cleandatafinal)
# =============================================================================

# 1a) How well do the traits (normalized as ratios) predict function?
# =============================================================================
# PCA on non-size corrected, log-transformed data
data.1a <- cleandatafinal[ , c(1:4, 7:20, 48, 49, 52)]
#Replace NAs with column means.
for (col.1a in 5:ncol(data.1a)){
  m.t <- mean(data.1a[ , col.1a], na.rm = TRUE)
  rows.to.replace <- which(is.na(data.1a[ , col.1a]))
  data.1a[rows.to.replace , col.1a ] <- m.t
}

#**ratio and model on opening lever ratio**
data.1a.LR <- (data.1a[ , c(6, 8)])
  names = c("In.div.In", "Out.div.In") # Put your names for new columns in here
data.1a.LR.ratios.t = f.divide(dataset = data.1a.LR, range = c(1:2), index= 1, names = names)
  data.1a.LR.ratios <- cbind(data.1a[ , c(2:4, 21)], data.1a.LR.ratios.t)
model.LR.1a <- lm(data.1a.LR.ratios$OutOpeningInLeverRatio ~ data.1a.LR.ratios$Out.div.In) # + data.1a.LR.ratios$habitat + data.1a.LR.ratios$watershed)
summary(model.LR.1a)

#**ratio and model on SI**
data.1a.SI <- data.1a[ , c(9:13)]
names = c("Buccal.div.EpaxHt", "EpaxAr.div.EpaxHt", "Gape.div.EpaxHt", "OutL.div.EpaxHt", "EpaxHt.div.EpaxHt")
data.1a.SI.ratios.t = f.divide(dataset = data.1a.SI, range = c(1:5), index= 5, names = names)
data.1a.SI.ratios <- cbind(data.1a[ , c(2:4,19)], data.1a.SI.ratios.t)
model.SI.1a <- lm(data.1a.SI.ratios$SI ~ data.1a.SI.ratios$Buccal.div.EpaxHt + data.1a.SI.ratios$EpaxAr.div.EpaxHt + data.1a.SI.ratios$Gape.div.EpaxHt + data.1a.SI.ratios$OutL.div.EpaxHt) # + data.1a.SI.ratios$habitat + data.1a.SI.ratios$watershed)
summary(model.SI.1a)

#**ratio and model on 4bar**
data.1a.4bar <- (data.1a[ , c(14:18)])
  names = c("Fixed.div.Fixed", "Coupler.div.Fixed", "Input.div.Fixed", "Output.div.Fixed", "Diagonal.div.Fixed") # Put your names for new columns in here #LAC helped me with this
data.1a.4bar.ratios.t = f.divide(dataset = data.1a.4bar, range = c(1:5), index= 1, names = names)
data.1a.4bar.ratios <- cbind(data.1a[ , c(2:4, 20)], data.1a.4bar.ratios.t)
model.4bar.1a <- lm(data.1a.4bar.ratios$KTcoeff ~ + data.1a.4bar.ratios$Coupler.div.Fixed + data.1a.4bar.ratios$Input.div.Fixed + data.1a.4bar.ratios$Output.div.Fixed + data.1a.4bar.ratios$Diagonal.div.Fixed) # + data.1a.4bar.ratios$habitat + data.1a.4bar.ratios$watershed)
summary(model.4bar.1a)
# =============================================================================


# 1b) How well does PC1 and PC2 (and rest of PCs as an aside) on raw data predict function?
# =============================================================================
# PCA on non-size corrected, log-transformed data
data.1b <- cleandatafinal[ , c(1:4, 7:20, 48, 49, 52)]
#Replace NAs with column means.
for (col.1b in 5:ncol(data.1b)){
  m.t <- mean(data.1b[ , col.1b], na.rm = TRUE)
  rows.to.replace <- which(is.na(data.1b[ , col.1b]))
  data.1b[rows.to.replace , col.1b ] <- m.t
}

#**PCA and model on opening lever ratio**
pc.LR.1b <- prcomp(data.1b[ , c( 6, 8)] , scale = TRUE) #size, out lever, opening inlever
pc.LR.1b.data <- cbind(data.1b[ , c(2:4, 5, 6, 8, 21)], pc.LR.1b$x)
model.LR.1b.full <- lm(pc.LR.1b.data$OutOpeningInLeverRatio ~ pc.LR.1b.data$PC1 + pc.LR.1b.data$PC2 ) # + pc.LR.1b.data$habitat + pc.LR.1b.data$watershed)
summary(model.LR.1b.full)

#**PCA and model on SI**
pc.SI.1b <- prcomp(data.1b[ , c(9, 10, 11, 12, 13)] , scale = TRUE) #Should 14 be included?? Looks to be a KT component #
#NB: consider sqrt(epaxial area)
pc.SI.1b.data <- cbind(data.1b[ , c(2:4, 5, 9:14, 19)], pc.SI.1b$x)
model.SI.1b.full <- lm(pc.SI.1b.data$SI ~ pc.SI.1b.data$PC1 + pc.SI.1b.data$PC2)# + pc.SI.1b.data$PC3 + pc.SI.1b.data$PC4 + pc.SI.1b.data$PC5) # pc.SI.1b.data$habitat + pc.SI.1b.data$watershed)
summary(model.SI.1b.full)

#**PCA and model on 4bar**
pc.4bar.1b <- prcomp(data.1b[ , c(14:18)] , scale = TRUE)
pc.4bar.1b.data <- cbind(data.1b[ , c(2:4, 5, 14:18, 20)], pc.4bar.1b$x)
model.4bar.1b.full <- lm(pc.4bar.1b.data$KTcoeff ~ pc.4bar.1b.data$PC1 + pc.4bar.1b.data$PC2) # + pc.4bar.1b.data$PC3 + pc.4bar.1b.data$PC4 + pc.4bar.1b.data$PC5) # + pc.4bar.1b.data$habitat + pc.4bar.1b.data$watershed)
summary(model.4bar.1b.full)
# =============================================================================


# 1c) How well does PC1 of size corrected data predict function?
# =============================================================================
# PCA on size corrected, log-transformed data
data.1c <- cleandatafinal[ , c(1:4, 21, 35:49, 52)]
#Replace NAs with column means.
for (col.1c in 5:ncol(data.1c)){
  m.t <- mean(data.1c[ , col.1c], na.rm = TRUE)
  rows.to.replace <- which(is.na(data.1c[ , col.1c]))
  data.1c[rows.to.replace , col.1c ] <- m.t 
}
#**PCA and model on opening lever ratio**
pc.LR.1c <- prcomp(data.1c[ , c(6, 8)] , scale = TRUE) #opening in lever.sc, outlever.sc
pc.LR.1c.data <- cbind(data.1c[ , c(2:4, 5, 6, 8, 21)], pc.LR.1c$x)
model.LR.1c.full <- lm(pc.LR.1c.data$OutOpeningInLeverRatio ~ pc.LR.1c.data$PC1 + pc.LR.1c.data$PC2) # + pc.LR.1c.data$habitat + pc.LR.1c.data$watershed)
summary(model.LR.1c.full)

#**PCA and model on SI**
pc.SI.1c <- prcomp(data.1c[ , c(9:13)] , scale = TRUE) #NB: consider sqrt(epaxial area)
pc.SI.1c.data <- cbind(data.1c[ , c(2:4, 5, 9:13, 19)], pc.SI.1c$x)
model.SI.1c.full <- lm(pc.SI.1c.data$SI ~ pc.SI.1c.data$PC1 + pc.SI.1c.data$PC2) # + pc.SI.1c.data$PC3 + pc.SI.1c.data$PC4 + pc.SI.1c.data$PC5) #+ pc.SI.1c.data$habitat + pc.SI.1c.data$watershed) #use PC axes that explain more than 1% of the variance
summary(model.SI.1c.full)

#**PCA and model on 4bar**
pc.4bar.1c <- prcomp(data.1c[ , c(14:18)] , scale = TRUE)
pc.4bar.1c.data <- cbind(data.1c[ , c(2:4, 5, 14:18, 20)], pc.4bar.1c$x)
model.4bar.1c.full <- lm(pc.4bar.1c.data$KTcoeff ~ pc.4bar.1c.data$PC1 + pc.4bar.1c.data$PC2) # + pc.4bar.1c.data$PC3 + pc.4bar.1c.data$PC4 + pc.4bar.1c.data$PC5) # + pc.4bar.1c.data$habitat + pc.4bar.1c.data$watershed) #use PC axes that explain more than 1% of the variance
summary(model.4bar.1c.full)
# =============================================================================


# 2a) Plots of parallelism in 3 functional traits
# =============================================================================
#plotting function
p.plot.parallelism <- function(data.to.plot, sites.to.plot, y.lab) {
means <- tapply(data.to.plot, sites.to.plot, mean, na.rm = TRUE)
sds <- tapply(data.to.plot, sites.to.plot, sd, na.rm = TRUE)
Ns <- summary(factor(sites.to.plot))
lowerCI <- means - sds/sqrt(Ns)
upperCI <- means + sds/sqrt(Ns)
xvals <- seq(1:32)
colors <- rep(c("blue", "dark green"), 16, each = TRUE)
datatoplot <- data.frame(means, sds, Ns, lowerCI, upperCI, xvals, colors)
difference <- datatoplot$means[seq(1,31, by = 2)] - datatoplot$means[seq(2,32, by = 2)]
datatoplot$difference <- rep(difference, times = 1, each = 2)
datatoplot <- datatoplot[order(datatoplot$difference),]
datatoplot$neworder <- seq(1:32)
plot(means ~ neworder, datatoplot, col = as.character(datatoplot$colors), pch = 16, cex = 1.5, axes = F, xlab = "", ylab = y.lab, ylim = c(min(datatoplot$lowerCI), max(datatoplot$upperCI)))  
box()
axis(2)
axis1names.abbr <- unique(substr(rownames(datatoplot), 1, 3))
#axis1names.final <- c("Moore" ,"Kennedy" ,"Pachena" ,"Frederick" ,"Misty" ,"Thiemer" ,"Swan", "Northy" ,"Comida" ,"Roberts", "Joe", "Beaver" ,"VillageBay", "Pye", "Muchalat" ,"Boot")
axis(1, at = seq(1.5, 31.5, by = 2), labels = axis1names.abbr, las = 2)
for(i in seq(1, 31, by = 2)){lines(c(i, i+1), c(datatoplot$means[i], datatoplot$means[i+1]))}
arrows(datatoplot$neworder, datatoplot$lowerCI, datatoplot$neworder, datatoplot$upperCI, angle = 90, length = 0.04, code = 3, col = as.character(datatoplot$colors))
}

#Plots
data.forplotting <- cleandatafinal[ , c(4, 48, 49, 52)]
#Replaced NAs with column means.
for (col.plot in 2:ncol(data.forplotting)){
  m.t <- mean(data.forplotting[ , col.plot], na.rm = TRUE)
  rows.to.replace <- which(is.na(data.forplotting[ , col.plot]))
  data.forplotting[rows.to.replace , col.plot ] <- m.t
}
par(mfrow = c(3, 1), mar = c(3, 4, 0, 1))
#LR
data.to.plot.LR <- data.forplotting$OutOpeningInLeverRatio ; sites.to.plot.LR <- substr(data.forplotting$fishID.univ, 1, 4) ; y.lab.LR <- "Lever Ratio"
p.plot.parallelism(data.to.plot.LR, sites.to.plot.LR, y.lab.LR)
#SI
data.to.plot.SI <- data.forplotting$SI ; sites.to.plot.SI <- substr(data.forplotting$fishID.univ, 1, 4) ; y.lab.SI <- "Suction Index"
p.plot.parallelism(data.to.plot.SI, sites.to.plot.SI, y.lab.SI)
#KT
data.to.plot.KT <- data.forplotting$KTcoeff ; sites.to.plot.KT <- substr(data.forplotting$fishID.univ, 1, 4) ; y.lab.KT <- "Opercular 4-bar KT"
p.plot.parallelism(data.to.plot.KT, sites.to.plot.KT, y.lab.KT)

# =============================================================================

# 2b) Effect size of habitat term: parallelism in 3 functional measures
# =============================================================================
data.2b <- cleandatafinal[ , c(4, 48, 49, 52)]
#Replace NAs with column means.
for (col.2b in 2:ncol(data.2b)){
  m.t <- mean(data.2b[ , col.2b], na.rm = TRUE)
  rows.to.replace <- which(is.na(data.2b[ , col.2b]))
  data.2b[rows.to.replace , col.2b ] <- m.t
}
data.2b$habitat <- substr(data.2b$fishID.univ, 4, 4) ; data.2b$watershed <- substr(data.2b$fishID.univ, 1, 3)
#LR
anova(model.2b.LR <- lm(data.2b$OutOpeningInLeverRatio ~ data.2b$habitat * data.2b$watershed))
  EtaSq(model.2b.LR)  
#SI
anova(model.2b.SI <- lm(data.2b$SI ~ data.2b$habitat * data.2b$watershed))
  EtaSq(model.2b.SI)    
#KT
anova(model.2b.KT <- lm(data.2b$KTcoeff ~ data.2b$habitat * data.2b$watershed))
  EtaSq(model.2b.KT)    
# =============================================================================
  
# 2c) Effect size of habitat term: parallelism in PC2 of the component traits for 3 functional measures
# =============================================================================
#LR
anova(model.2c.LR <- lm(pc.LR.1b.data$PC2 ~ pc.LR.1b.data$habitat * pc.LR.1b.data$watershed))
  EtaSq(model.2c.LR) 
  #for revisions Table 1b
  anova(model.2cc.LR <- lm(data.1a.LR.ratios$Out.div.In ~ data.1a.LR.ratios$habitat * data.1a.LR.ratios$watershed))
  EtaSq(model.2cc.LR)
#SI
anova(model.2c.SI <- lm(pc.SI.1b.data$PC2 ~ pc.SI.1b.data$habitat * pc.SI.1b.data$watershed))
  EtaSq(model.2c.SI)
  #for revisions Table 1b
  model.2cc.SI <- manova(as.matrix(data.1a.SI.ratios[, c(10:13)]) ~ data.1a.SI.ratios$habitat * data.1a.SI.ratios$watershed)
  EtaSq(model.2cc.SI)
#KT
anova(model.2c.KT <- lm(pc.4bar.1b.data$PC2 ~ pc.4bar.1b.data$habitat * pc.4bar.1b.data$watershed))
  EtaSq(model.2c.KT)
  #for revisions Table 1b
  anova(model.2cc.KT <- manova(as.matrix(data.1a.4bar.ratios[ , c(11:14)]) ~ data.1a.4bar.ratios$habitat * data.1a.4bar.ratios$watershed))
  EtaSq(model.2cc.KT)
# =============================================================================

# 2d) Plot of effect size for trait against effect size for function
# =============================================================================
table.to.plot <- as.data.frame(matrix(NA, nrow = 3, ncol = 2)); colnames(table.to.plot) <- c("HabEtaFunction", "HabEtaTrait") ; rownames(table.to.plot) <- c("LR", "SI", "KT")
  table.to.plot[1,1] <- round(EtaSq(model.2b.LR)[1,1], 3)
  table.to.plot[1,2] <- round(EtaSq(model.2c.LR)[1,1], 3)
  table.to.plot[2,1] <- round(EtaSq(model.2b.SI)[1,1], 3)
  table.to.plot[2,2] <- round(EtaSq(model.2c.SI)[1,1], 3)
  table.to.plot[3,1] <- round(EtaSq(model.2b.KT)[1,1], 3)
  table.to.plot[3,2] <- round(EtaSq(model.2b.KT)[1,1], 3)

par(mfrow = c(1, 1), mar = c(4, 4, 1, 1))
plot(table.to.plot$HabEtaTrait ~ table.to.plot$HabEtaFunction, pch = 16, cex = 3, xlab = "Effect Size (Eta2) of Habitat on Component Traits" , ylab = "Effect Size (Eta2) of Habitat on Functional Score", xlim = c(0, 0.05), ylim = c(0, 0.05))
text(0.006, 0.00, "KT: Most Complex") ; text(0.046, 0.0395, "LR: Least Complex") ; text(0.0185, 0.01, "SI")
abline(a = 0, b = 1, lty = 2)
#segments(x0 = 0.000, y0 = 0.000, x1 = 0.017, y1= 0.01, lty = 2) ; segments(x0 = 0.000, y0 = 0.000, x1 = 0.017, y1= 0.01, lty = 2)
# =============================================================================  
  
# 3) Background stats to round out paper
#=============================================================================
# 3a) Show that component morphology and function varies significantly among populations using MANOVA
# =============================================================================
#LR
data.LR <- data.1a.LR.ratios
comp.data.cols <- c(5,6) ; func.data.cols <- 4 ; population <- substr(data.LR$fishID.univ, 1, 4)
manova.LR.comp <- manova(as.matrix(data.LR[, comp.data.cols]) ~ population) ; summary(manova.LR.comp, test = "Pillai")
anova.LR.func <- aov(as.matrix(data.LR[, func.data.cols]) ~ population) ; summary(anova.LR.func)
EtaSq(manova.LR.comp); EtaSq(anova.LR.func)
#SI
data.SI <- data.1a.SI.ratios
comp.data.cols <- c(5:9) ; func.data.cols <- 4 ; population <- substr(data.SI$fishID.univ, 1, 4)
manova.SI.comp <- manova(as.matrix(data.SI[, comp.data.cols]) ~ population) ; summary(manova.SI.comp, test = "Pillai")
anova.SI.func <- aov(as.matrix(data.SI[, func.data.cols]) ~ population) ; summary(anova.SI.func)
EtaSq(manova.SI.comp); EtaSq(anova.SI.func)
#KT
data.KT <- data.1a.4bar.ratios
comp.data.cols <- c(5:9) ; func.data.cols <- 4 ; population <- substr(data.KT$fishID.univ, 1, 4)
manova.KT.comp <- manova(as.matrix(data.KT[, comp.data.cols]) ~ population) ; summary(manova.KT.comp, test = "Pillai")
anova.KT.func <- aov(as.matrix(data.KT[, func.data.cols]) ~ population) ; summary(anova.KT.func)
EtaSq(manova.KT.comp); EtaSq(anova.KT.func)
# =============================================================================

# 3b) Show that there are within pair habitat effects.
# =============================================================================
data.LR <- data.1a.LR.ratios; comp.data.cols.LR <- c(5,6)
data.SI <- data.1a.SI.ratios; comp.data.cols.SI <- c(5:9)
data.KT <- data.1a.4bar.ratios; comp.data.cols.KT <- c(5:9) 
func.data.cols <- 4 
output.comp <- data.frame(matrix(NA, nrow = 16, ncol = 3)) ; 
rownames(output.comp) <- unique(substr(data.LR$fishID.univ, 1, 3)); colnames(output.comp) <- c("LR", "SI", "KT")
output.func <- data.frame(matrix(NA, nrow = 16, ncol = 3)) ; 
rownames(output.func) <- unique(substr(data.LR$fishID.univ, 1, 3)); colnames(output.func) <- c("LR", "SI", "KT")
wshds <- unique(substr(data.LR$fishID.univ, 1, 3))

#Components
#Pillai and P
for (i in 1:16){
  #LR
  model1.LR.t <- manova(as.matrix(data.LR[substr(data.LR$fishID.univ, 1, 3) == wshds[i] , comp.data.cols.LR]) ~ data.LR$habitat[substr(data.LR$fishID.univ, 1, 3) == wshds[i]])
  output.comp[i, 1] <- paste("/", round(summary(model1.LR.t, test = "Pillai")$stats[1,2],2) , " (", round(summary(model1.LR.t, test = "Pillai")$stats[1,6],2) , ")" , sep = "")
  #SI
  model1.SI.t <- manova(as.matrix(data.SI[substr(data.SI$fishID.univ, 1, 3) == wshds[i] , comp.data.cols.SI]) ~ data.SI$habitat[substr(data.SI$fishID.univ, 1, 3) == wshds[i]])
  output.comp[i, 2] <- paste("/", round(summary(model1.SI.t, test = "Pillai")$stats[1,2],2), " (", round(summary(model1.SI.t, test = "Pillai")$stats[1,6],2) , ")" , sep = "")
  #KT
  model1.KT.t <- manova(as.matrix(data.KT[substr(data.KT$fishID.univ, 1, 3) == wshds[i] , comp.data.cols.KT]) ~ data.KT$habitat[substr(data.KT$fishID.univ, 1, 3) == wshds[i]])
  output.comp[i, 3] <- paste("/", round(summary(model1.KT.t, test = "Pillai")$stats[1,2],2), " (", round(summary(model1.KT.t, test = "Pillai")$stats[1,6],2) , ")" , sep = "")
} 

#Function
#F and P
for (i in 1:16){
  #LR
  model2.LR.t <- aov(as.matrix(data.LR[substr(data.LR$fishID.univ, 1, 3) == wshds[i] , func.data.cols]) ~ data.LR$habitat[substr(data.LR$fishID.univ, 1, 3) == wshds[i]])
  output.func[i, 1] <- paste("/", round(unlist(summary(model2.LR.t))[7], 2) , " (", round(unlist(summary(model2.LR.t))[9], 2) , ")" , sep = "")
  #SI
  model2.SI.t <- aov(as.matrix(data.SI[substr(data.SI$fishID.univ, 1, 3) == wshds[i] , func.data.cols]) ~ data.SI$habitat[substr(data.SI$fishID.univ, 1, 3) == wshds[i]])
  output.func[i, 2] <- paste("/", round(unlist(summary(model2.SI.t))[7], 2) , " (", round(unlist(summary(model2.SI.t))[9], 2) , ")" , sep = "")  
  #KT
  model2.KT.t <- aov(as.matrix(data.KT[substr(data.KT$fishID.univ, 1, 3) == wshds[i] , func.data.cols]) ~ data.KT$habitat[substr(data.KT$fishID.univ, 1, 3) == wshds[i]])
  output.func[i, 3] <- paste("/", round(unlist(summary(model2.KT.t))[7], 2) , " (", round(unlist(summary(model2.KT.t))[9], 2) , ")" , sep = "")
}
#=============================================================================


# 3c) histogram of eta-sq values for all traits from Nature EE with arrows added for our observed EtaSq values, to show that divergence likely under selection
#=============================================================================
eta.from.natureEE <- read.csv("/Users/yestuart/Desktop/lake.stream.github.analysis/project.1.parallelism/proj.1.output.data/05a.etas.morpho.linear.model.csv", header = TRUE, stringsAsFactors = FALSE)
data.eta <- eta.from.natureEE[c(75, 97:140, 141:181), c(1:4)]
data.eta$trait.type <- c(rep("non-trophic", 16), rep("SI", 2), rep("non-trophic", 19), rep("LR", 2), rep("SI", 3), rep("KT", 5), rep("non-trophic", 39))

par(mar = c(4, 4, 1, 1), mfrow = c(2, 2))
#adding habitat eta squared arrows for function
hist(data.eta$etaSq.habitat, xlim = c(0, 0.25), col = "gray90", ylim = c(0,65), main = "", xlab = "Function: Habitat Effect Size")
arrows(0.044, 55, 0.044, 30 , length = 0.15, angle = 20)
text(0.045, 60, "LR")
arrows(0.017, 60, 0.017, 30 , length = 0.15, angle = 20)
text(0.017, 64, "SI")
arrows(0.00, 60, 0.00, 30, length = 0.15, angle = 20)
text(0.000, 64, "KT")
text(.24, 60, "(A)")
  sum(data.eta$etaSq.habitat >  0.044)/length(data.eta$etaSq.habitat) # LR: 0.267
  sum(data.eta$etaSq.habitat >  0.017)/length(data.eta$etaSq.habitat) # SI: 0.407
  sum(data.eta$etaSq.habitat >  0.000)/length(data.eta$etaSq.habitat) # KT: 1

#adding interaction eta squared arrows for function
hist(data.eta$etaSq.habitat.x.watershed, col = "gray90", ylim = c(0,35), xlim = c(0, 0.35), main = "", xlab = "Function: Interaction Effect Size")
arrows(0.069, 33, 0.069, 22 , length = 0.15, angle = 20)
text(0.06, 35, "LR")
arrows(0.285, 32, 0.285, 20 , length = 0.15, angle = 20)
text(0.285, 35, "SI")
arrows(0.073, 30, 0.073, 18, length = 0.15, angle = 20)
text(0.09, 31, "KT")
text(.34, 35, "(B)")
  sum(data.eta$etaSq.habitat.x.watershed>0.069)/length(data.eta$etaSq.habitat.x.watershed) # LR: 0.860
  sum(data.eta$etaSq.habitat.x.watershed >  0.285)/length(data.eta$etaSq.habitat.x.watershed) # SI: 0.012
  sum(data.eta$etaSq.habitat.x.watershed >  0.073)/length(data.eta$etaSq.habitat.x.watershed) # KT: .826

  #adding habitat eta squared arrows for component traits
hist(data.eta$etaSq.habitat, xlim = c(0, 0.25), col = "gray90", ylim = c(0,65), main = "", xlab = "Traits: Habitat Effect Size")
arrows(0.044, 55, 0.044, 30 , length = 0.15, angle = 20)
text(0.045, 60, "LR")
arrows(0.035, 60, 0.035, 30 , length = 0.15, angle = 20)
text(0.034, 64, "SI")
arrows(0.00, 60, 0.00, 30, length = 0.15, angle = 20)
text(0.000, 64, "KT")
text(0.24, 60, "(C)")
  sum(data.eta$etaSq.habitat >  0.044)/length(data.eta$etaSq.habitat) # LR: 0.267
  sum(data.eta$etaSq.habitat >  0.035)/length(data.eta$etaSq.habitat) # SI: 0.291
  sum(data.eta$etaSq.habitat >  0.000)/length(data.eta$etaSq.habitat) # KT: 1

#adding interaction eta squared arrows for component traits
hist(data.eta$etaSq.habitat.x.watershed, col = "gray90", ylim = c(0,35), xlim = c(0, 0.35), main = "", xlab = "Traits: Interaction Effect Size")
arrows(0.069, 32, 0.069, 22 , length = 0.15, angle = 20)
text(0.09, 30, "LR")
arrows(0.066, 34, 0.066, 25 , length = 0.15, angle = 20)
text(0.066, 35, "SI")
arrows(0.061, 30, 0.061, 20, length = 0.15, angle = 20)
text(0.05, 30, "KT")
text(.34, 35, "(D)")
  sum(data.eta$etaSq.habitat.x.watershed>0.069)/length(data.eta$etaSq.habitat.x.watershed) # LR: 0.860
  sum(data.eta$etaSq.habitat.x.watershed >  0.066)/length(data.eta$etaSq.habitat.x.watershed) # SI: 0.872
  sum(data.eta$etaSq.habitat.x.watershed >  0.061)/length(data.eta$etaSq.habitat.x.watershed) # KT: 0.884
#=============================================================================

# 3d) P-matrices of LR, SI, and KT component traits to estimate dimensionality
#=============================================================================
p.data <- data.1a[ , c(2:5, 7:18)]
#Replace NAs with column means.
for (col.1a in 4:ncol(p.data)){
  m.t <- mean(p.data[ , col.1a], na.rm = TRUE)
  rows.to.replace <- which(is.na(p.data[ , col.1a]))
  p.data[rows.to.replace , col.1a ] <- m.t
}

#Mean-scale
cols.to.transform <- c(4:16)
trait.means <- as.data.frame(apply(p.data[,cols.to.transform], 2, mean, na.rm = TRUE))
trait <- rownames(trait.means) ; trait.means <- data.frame(trait, trait.means) ; colnames(trait.means)[2] <- "trait.mean"
mean.scaled.data <- as.data.frame(matrix(NA, nrow = length(p.data$fishID.univ), ncol = length(p.data))) ; colnames(mean.scaled.data) <- colnames(p.data) 
#trait.means$trait[1:33] == colnames(mean.scaled.data)[2:34]
for (i in 2:length(trait.means$trait.mean)) mean.scaled.data[,i+3] <- p.data[,(i+3)] / trait.means$trait.mean[i]
mean.scaled.data$fishID.univ <- p.data$fishID.univ
mean.scaled.data$habitat <- substr(mean.scaled.data$fishID.univ, 4, 4)
mean.scaled.data$watershed <- substr(mean.scaled.data$fishID.univ, 1,3)
mean.scaled.data$standard.length.from.calipers <- p.data$standard.length.from.calipers

#Save residuals against SL
resid.for.Pmatrix3 <- as.data.frame(matrix(NA, nrow = length(mean.scaled.data$fishID.univ), ncol = length(mean.scaled.data)))
colnames(resid.for.Pmatrix3) <- paste(colnames(mean.scaled.data), ".resid", sep = "") 
colnames(resid.for.Pmatrix3)[c(1, 2, 3, 4)] <- c("habitat", "watershed", "fishID.univ", "Standard.Length.mm")
#linear model to size correct
for (i in 5:length(mean.scaled.data)){
  outcome <- lm(mean.scaled.data[ , i] ~ mean.scaled.data$standard.length.from.calipers, na.action = na.exclude)
  resid.for.Pmatrix3[ , i] <- residuals(outcome)
}
resid.for.Pmatrix3$fishID.univ <- mean.scaled.data$fishID.univ 
resid.for.Pmatrix3$Standard.Length.mm <- mean.scaled.data$standard.length.from.calipers
resid.for.Pmatrix3$habitat <- mean.scaled.data$habitat
resid.for.Pmatrix3$watershed <- mean.scaled.data$watershed
allfishdat <- resid.for.Pmatrix3

##Calculate PC scores and P matrix for LR##
p.LR <- allfishdat[ ,c(1:4, 5:6)]
PCA.LR <- prcomp(p.LR[, 5:6], scale = TRUE)
  PC.scores.LR <- data.frame(p.LR$fishID.univ, PCA.LR$x); names(PC.scores.LR)[1] <- "fishID.univ"
  screeplot(PCA.LR)
  PCA.LR$rotation
  summary(PCA.LR)
vcm.LR <- round(cov(PC.scores.LR[, c(2:3)]), 4)
  eigenvalue.LR <- unlist(eigen(vcm.LR)[1])
  matrix.size.LR <- sum(eigenvalue.LR)
  matrix.dimensions.LR <- sum(eigenvalue.LR/eigenvalue.LR[1])

##Calculate PC scores and P matrix for SI##
p.SI <- allfishdat[ ,c(1:4, 7:11)]
PCA.SI <- prcomp(p.SI[, c(5:9)], scale = TRUE)
  PC.scores.SI <- data.frame(p.SI$fishID.univ, PCA.SI$x); names(PC.scores.SI)[1] <- "fishID.univ"
  screeplot(PCA.SI)
  PCA.SI$rotation
  summary(PCA.SI)
vcm.SI <- round(cov(PC.scores.SI[, c(2:6)]), 4)
  eigenvalue.SI <- unlist(eigen(vcm.SI)[1])
  matrix.size.SI <- sum(eigenvalue.SI)
  matrix.dimensions.SI <- sum(eigenvalue.SI/eigenvalue.SI[1])

##Calculate PC scores and P matrix for KT##
p.KT <- allfishdat[ ,c(1:4, 12:16)]
PCA.KT <- prcomp(p.KT[, c(5:9)], scale = TRUE)
  PC.scores.KT <- data.frame(p.KT$fishID.univ, PCA.KT$x); names(PC.scores.KT)[1] <- "fishID.univ"
  screeplot(PCA.KT)
  PCA.KT$rotation
  summary(PCA.KT)
vcm.KT <- round(cov(PC.scores.KT[, c(2:6)]), 4)
  eigenvalue.KT <- unlist(eigen(vcm.KT)[1])
  matrix.size.KT <- sum(eigenvalue.KT)
  matrix.dimensions.KT <- sum(eigenvalue.KT/eigenvalue.KT[1])
#=============================================================================

# 3d) Incorporating environmental data
#=============================================================================
#traits or function ~ PC1 of diet data + habitat * watershed
stomach <- read.csv(paste("/Users/yestuart/Desktop/lake.stream.github.analysis/cleaned.data/03.dietCntnt.csv", sep = ""), stringsAsFactors = FALSE)
  stomach.t <- stomach[-(which(substr(stomach$fishID.univ, 4, 4) == "M")), ]
  names(stomach.t) 
  stomach.tt <- stomach.t[, c(1, 2:31, 33:45, 47:58)]
  str(stomach.tt)
  
  #convert percentages to present/absent
  f.colname.number(stomach.tt)
  stomach.PA <- stomach.tt[,c(1, 2:56)]
  for (stomPAcols in 2:length(colnames(stomach.PA))) {
    x.t <- which(stomach.PA[,stomPAcols] > 0)
    stomach.PA[x.t, stomPAcols] <- 1
  }
  colnames(stomach.PA) <- paste(colnames(stomach.PA), "PA", sep = ".")
  
  #z-transform
  stomach.PA.z <- as.data.frame(apply(stomach.PA[2:length(stomach.PA)], 2, f.zscore))
  f.colname.number(stomach.PA.z)
  levels(as.factor(stomach.tt$Psychodidae)) #0; only present in marine fish [22]
  levels(as.factor(stomach.tt$Nymphomyiidae)) #0 ; only present in marine fish [44]
  levels(as.factor(stomach.tt$Isopods)) #0 ; only present in marine fish [50]
  levels(as.factor(stomach.tt$crustacean.part)) #0 ; only present in marine fish. [53]
  stomach.PA.z <- stomach.PA.z[ , -c(22, 44, 50, 53)]
  colnames(stomach.PA.z) <- paste(colnames(stomach.PA.z), "zscore", sep = ".")
  stomach.PA.z <- cbind(stomach.PA$fishID.univ, stomach.PA.z)
  colnames(stomach.PA.z)[1] <- "fishID.univ"
  
  #principal component analysis
    #remove NAs
  stomach.PA.z.noNA <- na.omit(stomach.PA.z)
  stomach.pca <- prcomp(stomach.PA.z.noNA[ , 2:length(stomach.PA.z.noNA)])
  summary(stomach.pca) ; plot(stomach.pca)
  stomach.pca.scores <- as.data.frame(stomach.pca$x)
  plot(stomach.pca.scores$PC1 ~ stomach.pca.scores$PC2, xlim = c(-5,2)) #some outliers in PC2
  stomach.pca.scores.f <- data.frame(stomach.PA.z.noNA$fishID.univ, stomach.pca.scores) ; colnames(stomach.pca.scores.f)[1] <- "fishID.univ"
  
  #Calculate means for PCs
  stomach.pca.pop.means <- as.data.frame(round(apply(stomach.pca.scores.f[ , 2:length(stomach.pca.scores.f)] , 2, tapply, substr(stomach.pca.scores.f$fishID.univ, 1, 4), mean, na.rm = TRUE), 4))
  stomach.pca.pop.means$site <- rownames(stomach.pca.pop.means)
  
  #Calculate means for component traits, function. Component trait values are size corrected, and taken from allfishdat in the previous section.
  component.pop.means <- as.data.frame(round(apply(allfishdat[ , 4:length(allfishdat)] , 2, tapply, substr(allfishdat$fishID.univ, 1, 4), mean, na.rm = TRUE), 4))
    component.pop.means$site <- rownames(component.pop.means)
  function.pop.means <- as.data.frame(round(apply(cleandatafinal[ , c(52, 48, 49)] , 2, tapply, substr(cleandatafinal$fishID.univ, 1, 4), mean, na.rm = TRUE), 4))
    function.pop.means$site <- rownames(function.pop.means)
  
  #function models
    #LR
    func.x.env.model.full.LR <- lm(function.pop.means$OutOpeningInLeverRatio ~ stomach.pca.pop.means$PC1 + substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(func.x.env.model.full.LR) ; anova(func.x.env.model.full.LR)
    func.x.env.model.reduced.LR <- lm(function.pop.means$OutOpeningInLeverRatio ~ substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(func.x.env.model.reduced.LR) ; anova(func.x.env.model.reduced.LR)
      anova(func.x.env.model.reduced.LR, func.x.env.model.full.LR) # model with PC1 for gut contents not significantly better than reduced model
    #SI
    func.x.env.model.full.SI <- lm(function.pop.means$SI ~ stomach.pca.pop.means$PC1 + substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(func.x.env.model.full.SI) ; anova(func.x.env.model.full.SI)
    func.x.env.model.reduced.SI <- lm(function.pop.means$SI ~ substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(func.x.env.model.reduced.SI) ; anova(func.x.env.model.reduced.SI)
      anova(func.x.env.model.full.SI, func.x.env.model.reduced.SI) # model with PC1 for gut contents not significantly better than reduced model
    #KT
    func.x.env.model.full.KT <- lm(function.pop.means$KT ~ stomach.pca.pop.means$PC1 + substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(func.x.env.model.full.KT) ; anova(func.x.env.model.full.KT)
    func.x.env.model.reduced.KT <- lm(function.pop.means$KT ~ substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(func.x.env.model.reduced.KT) ; anova(func.x.env.model.reduced.KT)
    anova(func.x.env.model.full.KT, func.x.env.model.reduced.KT) # marginally signficant

  #component trait models
    #LR
      comp.x.env.model.full.LR <- manova(as.matrix(component.pop.means[ , c(2:3)]) ~ stomach.pca.pop.means$PC1 + substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(comp.x.env.model.full.LR, test = "Pillai")
      comp.x.env.model.reduced.LR <- manova(as.matrix(component.pop.means[ , c(2:3)]) ~ substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(comp.x.env.model.reduced.LR, test = "Pillai")
      anova(comp.x.env.model.full.LR, comp.x.env.model.reduced.LR) # full model not significantly better
    #SI
    comp.x.env.model.full.SI <- manova(as.matrix(component.pop.means[ , c(4:8)]) ~ stomach.pca.pop.means$PC1 + substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(comp.x.env.model.full.SI, test = "Pillai")
      comp.x.env.model.reduced.SI <- manova(as.matrix(component.pop.means[ , c(4:8)]) ~ substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(comp.x.env.model.reduced.SI, test = "Pillai")
      anova(comp.x.env.model.full.SI, comp.x.env.model.reduced.SI) # full model not significantly better
    #KT
    comp.x.env.model.full.KT <- manova(as.matrix(component.pop.means[ , c(9:13)]) ~ stomach.pca.pop.means$PC1 + substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(comp.x.env.model.full.KT, test = "Pillai")
    comp.x.env.model.reduced.KT <- manova(as.matrix(component.pop.means[ , c(9:13)]) ~ substr(stomach.pca.pop.means$site, 4,4) + substr(stomach.pca.pop.means$site, 1, 3)) ; summary(comp.x.env.model.reduced.KT, test = "Pillai")
      anova(comp.x.env.model.full.KT, comp.x.env.model.reduced.KT) # marginally significantly better
#=============================================================================
  