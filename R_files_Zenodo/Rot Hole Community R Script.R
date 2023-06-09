## SUPPLEMENTARY MATERIAL WINDSOR ET AL. COMMUNITY STRUCTURE OF ROT HOLE INVERTEBRATE COMMUNITIES...## 

# DEVELOPED IN R BY F. M. WINDSOR (17/01/2018)
# DATA EXPLORATION REMOVED FOR CONCISENESS, ONLY FINAL MODELS PROVIDED
# FOR ENQUIRIES PLEASE CONTACT windsorfm@cardiff.ac.uk 

######################################################################################################
############################################### SET-UP ###############################################
######################################################################################################

rm(list=ls())
setwd("D:/")
library(cooccur); library(vegan); library(betapart); library(vegetarian); library(mvabund)
library(randomForest); library(ggplot2); library(dplyr); library(bipartite)

######################################################################################################
############################################# DATA INPUT #############################################
######################################################################################################

dframe1 <- read.csv("community_family.csv", header = T, row.names = 1)
dframe2 <- read.csv("environment.csv")

######################################################################################################
########################################### DATA ANALYSIS ############################################
######################################################################################################

######################################## COMMUNITY STRUCTURE #########################################

#### COMMUNITY DESCRIPTIONS ####

richness <- specnumber(dframe1)
abundance <- rowSums(dframe1)
sum(abundance)

dframe1a <- cbind(dframe1, richness)
dframe1a <- cbind(dframe1, abundance)

## Alpha diversity ## 
alpha <- diversity(dframe1, index = 'simpson') #Calculate the alpha diversity for each site
mean(alpha)
sd(alpha)
sd(alpha)/sqrt(length(alpha))
dframe1a <- cbind(dframe1a, alpha) #Bind the alpha diversity to the existing data
dframe2$diversity <- alpha #Bind diversity to environmental data

kruskal.test(diversity ~ Site, data = dframe2)

## Beta diversity ##
beta.multi.abund(dframe1, index.family = "bray") #1-Turnover, 2-nestedness, 3-beta-diversity
beta <- beta.pair.abund(dframe1)


## Gamma diversity ## 
d(dframe1, lev = 'gamma', q = 2, boot = TRUE) #q = 2 means inverse simpson's
d(dframe1, lev = 'beta', q = 2, boot = TRUE)
d(dframe1, lev = 'alpha', q = 2, boot = TRUE) # Matches relatively well to invsimpson above


#### SPECIES COOCCURRENCE ####

dframe1_binomial <- data.frame(t((dframe1>0)*1L)) #Abundance to binomial and transpose the data t()

cooccur <- cooccur(dframe1_binomial, type = "spp_site", spp_names = TRUE, thresh = TRUE, eff_matrix = T)
summary(cooccur)
plot(cooccur)
obs.v.exp(cooccur)
pair.profile(cooccur)
pair.attributes(cooccur)

ptab <- cooccur$results
test <- cooccur$spp_key
ptab$signs <- ifelse(ptab$p_gt>=0.05,0,1) + ifelse(ptab$p_lt>=0.05,0,-1)
exp_cooccur <- ptab$exp_cooccur
obs_cooccur <- ptab$obs_cooccur
signs <- ptab$signs

plot <- ggplot(aes(x=exp_cooccur, y=obs_cooccur), data = ptab) + theme_bw() + 
  geom_point(aes(fill = factor(signs, levels = c(-1, 0, 1))), colour = "black", pch = 21, size = 4) +
  scale_fill_manual(values = c("#FFCC66","dark gray","light blue"), name = "Probability of \\nco-occurrence", 
                    labels = c("Negative", "Random", "Positive"), drop = FALSE) + 
  theme(plot.title = element_text(vjust = 2,size = 20, face = "bold"), legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12), axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "black"),axis.text.x=element_text(hjust = 0, vjust = 1), 
        legend.position = c(0.16,0.82), legend.background = element_rect(colour = "black", size =0.5),
        panel.border = element_rect(colour = "black", fill = NA, size =0.5)) +
  xlab("Expected co-occurrences") + 
  ylab("Observed co-occurrences") + 
  geom_abline(color="black", size = 0.8, linetype = "dashed") + 
  xlim(0,22.5) + 
  ylim(0,22.5)
plot


#### EcoSIM - C-scores #### 

cooccur_model <- EcoSimR::cooc_null_model(dframe1_binomial, algo = "sim9", nReps = 1000, burn_in = 500)
plot(cooccur_model, type = "cooc")


#### COMMUNITY STRUCTURE ####

NMDS <- metaMDS(dframe1, trymax = 100, distance = "jaccard", k = 2)

data.scores <- as.data.frame(scores(NMDS))  
data.scores$Sample <- rownames(data.scores)

merge <- merge(data.scores, dframe2, by = "Sample", all.x = T, all.y = T)

plot <- merge(merge, aggregate(cbind(mean.x=NMDS1, mean.y=NMDS2)~Site, merge, mean), by='Site')

plot1 <- ggplot(plot) +  
  geom_segment(aes(x=mean.x, y=mean.y, xend=NMDS1, yend=NMDS2, colour = Site)) + 
  geom_point(aes(x=mean.x, y=mean.y, fill = Site), pch = 21, cex = 5) +
  geom_point(aes(NMDS1, NMDS2, fill = Site), pch = 21, cex = 4) + 
  theme_bw() + 
  xlim(-0.75,0.75) +
  ylim(-0.75,0.75) +
  xlab("NMDS1") + 
  ylab("NMDS2") +
  theme(plot.title = element_text(vjust = 2,size = 20, face = "bold"), legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12), axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "black"),axis.text.x=element_text(hjust = 0, vjust = 1), 
        legend.position = c(0.16,0.82), legend.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), plot.margin = unit(c(1, 1, 2, 1), "lines")) 
plot1


#### SUMMARISE AXES #### 

envdata <- select(dframe2, Water.Content, Tree.Diameter, Rot.Hole.Size, Density, 
                  Water.Potential)

envdata <- scale(envdata)

pca <- princomp(envdata)
summary(pca) #Axis 1 = Soft and wet to dry and hard; Axis 2 = Large to small features
plot(pca) 
biplot(pca)
print(pca$loadings)
scores(pca, display = "species")

dframe2$C1 <- pca$scores[,1] #33% of variation
dframe2$C2 <- pca$scores[,2] #25% of variation

pcs <- select(dframe2, C1, C2, Sample)
pcs$SH <- dframe2$diversity

pca_data <- data.frame(pca$scores[,1:2])
pca_arrows <- as.data.frame(pca$loadings[,1:2])
pca_arrows$labels <- c("Water content", "Tree diameter",
                       "Rot hole size", "Density", 
                       "Water Potential")

n <- nrow(pca_arrows)
new_data <- data.frame(X = c(rep(0,n), rep(pca_arrows$Comp.1)),
                       Y = c(rep(0,n), rep(pca_arrows$Comp.2)))
new_data$group <- as.factor(rep(1:n, times = 2))

new_data1 <- new_data[c(1,5,6,10),]
new_data2 <- new_data[-c(1,5,6,10),]

pca_arrows$Comp.1 <- pca_arrows$Comp.1 + c(-0.15,0.23,-0.1,0,0)
pca_arrows$Comp.2 <- pca_arrows$Comp.2 + c(0.05,0,0.01,0.05,0.025)


plot <- ggplot(aes(x = Comp.1, y = Comp.2), data = pca_data) + 
  geom_point(size = 4, pch = 21, fill = "grey50") +
  geom_line(aes(x = X*4, y = Y*4, group = group), data = new_data1, arrow = arrow(type = "closed", ends = "last", length = unit(3, "mm")), size = 0.8) +
  geom_line(aes(x = X*4, y = Y*4, group = group), data = new_data2, arrow = arrow(type = "closed", ends = "first", length = unit(3, "mm")), size = 0.8) +
  geom_text(aes(x = Comp.1*4.1, y = Comp.2*4.1, label = labels), data = pca_arrows) +
  theme_bw() + 
  theme(plot.title = element_text(vjust = 2,size = 20, face = "bold"), legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12), axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "black"),axis.text.x=element_text(hjust = 0, vjust = 1), 
        legend.position = c(0.16,0.82), legend.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), plot.margin = unit(c(1, 1, 2, 1), "lines")) + 
  xlab("PC1 (33.6%)") + 
  ylab("PC2 (24.8%)")
plot


#### COMMUNITY NESTEDNESS ####
dframe1_nest <- data.frame((dframe1>0)*1L)

nest <- nestedtemp(dframe1_nest)

binmatnest <- nestedness(dframe1_nest, null.models = T, n.nulls = 100, popsize = 30, n.ind = 7, n.gen = 2000, binmatnestout=FALSE)
summary(binmatnest)

data <- nest$comm
data <- t(data[nrow(data):1,])

dframe <- data.frame(x = seq(0,50), y = seq(0,50))
dframe$y <- (1-nest$smooth$y * 33)
dframe$x <- (nest$smooth$x * 66)

melt_data <- melt(data)
melt_data$value <- as.factor(melt_data$value)
melt_data$Var2 <- as.numeric(melt_data$Var2)

plot <- ggplot(aes(x = Var1, y = Var2-33, fill = value), data = melt_data) + 
  geom_tile() + 
  scale_fill_manual(values = c("white","grey30"), labels = c("Absent", "Present"), name = c("Species Occurrence")) +
  geom_line(aes(x = x+1, y = y), data = dframe, inherit.aes = F, size = 0.8, linetype = "dashed") + 
  coord_cartesian(expand = F) + 
  theme_bw() + 
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12), axis.title = element_text(size = 14), 
        axis.text = element_blank(), legend.position = c(0.8, 0.15), legend.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), plot.margin = unit(c(1, 1, 2, 1), "lines"), axis.ticks = element_blank()) +
  ylab("Sample sites") + 
  xlab("Species")
plot

nested <- data.frame(nest$r)
nested$Sample <- rownames(nested)
nested <- merge(nested,dframe2, by = "Sample")

model1 <- glm(nest.r ~ C1 + C2, family = gaussian (link = "identity"),
              data = nested)
model0 <- glm(nest.r ~ 1, family = gaussian (link = "identity"),
              data = nested)
summary.lm(model1)
anova(model0, model1, test = "Chi")
plot(model1)

preds <- predict(model1, type = "response", se.fit = T)
nested$fit <- preds$fit
nested$se.fit <- preds$se.fit

plot <- ggplot(aes(x = C1, y = nest.r, fill = Site), data = nested) +
  geom_point(pch = 21, size = 4) +
  geom_line(aes(y = fit, x = C1, colour = Site), size = 0.8, linetype = "dashed") +
  #geom_line(aes(y = fit-se.fit, x = C1, colour = Site), size = 1, lty = "dashed") +
  #geom_line(aes(y = fit+se.fit, x = C1, colour = Site), size = 1, lty = "dashed") +
  theme_bw() + 
  theme(plot.title = element_text(vjust = 2,size = 20, face = "bold"), legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12), axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "black"),axis.text.x=element_text(hjust = 0, vjust = 1), 
        legend.position = c(0.16,0.82), legend.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size =0.5), plot.margin = unit(c(1, 1, 2, 1), "lines")) +
  ylab("Nestedness temperature (0-1)") + 
  xlab("PC1\\nWet and low density                        Dry and high density") + 
  scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2), limits = c(-4.05,2))
plot



####################################### ENVIRONMENTAL FACTORS ########################################

#### ORDINATION TECHNIQUES ####

dca <- decorana(dframe1, iweigh = 1, ira = 0) #iweigh = 1/5th downweighting for rare taxa
summary(dca)                                  #ira = use detrended 
plot(dca) #As axes are generally below 2.5 use RDA rather than CCA (not a good enough gradient)

# Select variables of interest for further analysis
dframe2a <- select(dframe2, Tree.Diameter, Rot.Hole.Size, Density, Water.Content)
row.names(dframe2a) <- dframe2$Sample

#Need to rescale some of the variables (Tree.Diameter and Rot.Hole.Size)
dframe2a$Tree.Diameter <- rescale(dframe2a$Tree.Diameter)
dframe2a$Rot.Hole.Size <- rescale(dframe2a$Rot.Hole.Size)

rda <- rda(X = dframe1, Y = dframe2a)
summary(rda)
plot(rda)

test <- dca$rproj
dframe2$dca2 <- test[,2]

plot(dframe2$Water.Content,dframe2$dca2)


#### MULTIVARIATE GLMS ####
spp <-mvabund(as.matrix(dframe1))
c1<- dframe2$C1
c2<- dframe2$C2
Site <- dframe2$Site

model1 <- manyglm(as.matrix(dframe1) ~ Site + C1,
                  family = "poisson")


aov.many <- anova.manyglm(model1, p.uni = "adjusted") # 1 min 47 secs
print(aov.many)
print.manyglm(model1)
summary.manyglm(model1, symbolic.cor = TRUE, show.est = TRUE)
residuals.manyglm(model1)
plot(model1)
predict(model1, p.uni = "adjusted")

best.r.sq(spp ~ Site + c1)



