###################
### Preparation ###
###################

# libraries
#####
library(psych) # for PCA
library(lme4) # for (G)LMM
library(lmerTest) # for (G)LMM
library(ggplot2) # for plots
library(ggbeeswarm) # for plots
library(ggpubr) # for plots
library(FactoMineR) # for PCA
library(factoextra) # for PCA
library(MASS) # for BoxCox transformation
#####

# data: load, check, prepare, relevel, reorder
#####
### load data
data <- read.table("behavioural_data.txt", header=T,stringsAsFactors=TRUE)
View(data)

### create data needed
# EMP
data$EPM.rel.OA.entries = data$EPM.entries.OA/(data$EPM.entries.CA+data$EPM.entries.OA) # calculate relative number of open arm entries for EPM
data$EPM.rel.OA.duration = data$EPM.duration.OA/(data$EPM.duration.CA+data$EPM.duration.OA) # calculate relative time spent on open arms on EPM
# FET
data$FET.duration <- 15*60-data$FET.outside.time # calculate time spent in the FET by substracting time spent outside from max. observation time (900s)
# LM
data$LM.rel.diff.dur <- (data$LM.duration2/data$LM.duration1) # calculate relative difference between rounds
data$LM.rel.diff.dist <- (data$LM.distance2/data$LM.distance1) # calculate relative difference between rounds
data$LM.rel.diff.mist <- (data$LM.mistakes2/data$LM.mistakes1) # calculate relative difference between rounds
# centering
data$genotypeC = as.numeric(data$genotype=="B6D2F1N") - 0.5
data$environmentC = as.numeric(data$environment=="complex") - 0.5
data$TSC = as.numeric(data$TS=="trained") - 0.5

### relevel
data <- within(data, genotype <- relevel(genotype, ref = "C57BL/6J"))
data <- within(data, environment <- relevel(environment, ref = "scarce"))
data <- within(data, TS <- relevel(TS, ref = "non-trained"))

### re-order columns in dataset (for PCA)
col_order <- c("EPM.distance","EPM.rel.OA.entries","EPM.rel.OA.duration",
               "OFT.distance","OFT.entries","OFT.duration.center", 
               "FET.duration","FET.entries","FET.distance","FET.latency",        
               "LM.rel.diff.dur","LM.rel.diff.dist","LM.rel.diff.mist",
               "ID","cage","genotypeC","environmentC","TSC",
               "EPM.entries.OA","EPM.duration.OA","EPM.entries.CA","EPM.duration.CA",
               "FET.outside.time",
               "LM.duration1","LM.distance1","LM.mistakes1","LM.duration2","LM.distance2","LM.mistakes2",
               "genotype","environment","TS")
data <- data[, col_order]
#####



############################
### Statistical analysis ###
############################

# Supplementary Table 2: Spatial learning in the labyrinth maze (LM). 
#####
wilcox.test(data$LM.duration1,data$LM.duration2,paired=T) 
wilcox.test(data$LM.distance1,data$LM.distance2,paired=T)
wilcox.test(data$LM.mistakes1,data$LM.mistakes2,paired=T)
#####

# Supplementary Table 4: Statistical analysis of CJB test and behavioral test battery
# behavioral part only
#####
### check suitability for PCA
# Kaiser-Meyer-Olkin Kriterium (KMO)
KMO(na.omit(data[,c(1:13)])) # library(psych) 
# Bartlett test of sphericity (BTS)
Cor <- cor(na.omit(data[,c(1:13)]))
cortest.bartlett(Cor, n = 36) # library(psych)

### PCA
PCA1 <- prcomp(na.omit(data[,c(1:13)]), center = TRUE,scale. = TRUE)
summary(PCA1) # importance of components
PCA1$rotation # extract loadings

### GLMM
# prepare dataset
data_reduced <- na.omit(data)
data_reduced_withPC <- cbind(data_reduced, PCA1$x)

# model function
model <- function(x){
  lmer(x~G*E+TS+(1|cage))
}

G <- data_reduced_withPC$genotypeC
E <- data_reduced_withPC$environmentC
TS <- data_reduced_withPC$TSC
cage <- data_reduced_withPC$cage

# model PC1
lmer_PC1 <- model(data_reduced_withPC$PC1) # library(lme4)
anova(lmer_PC1) # library(lmerTest)
#summary(lmer_PC1)

# model PC2, with square root transformation of PC2
b <- sqrt(data_reduced_withPC$PC2+abs(min(data_reduced_withPC$PC2)))
lmer_PC2transf <-model(b) # library(lme4)
anova(lmer_PC2transf) # library(lmerTest)
#summary(lmer_PC2transf)

### Full model test (comparing final model and null model)
# PC1
PC1_mod_final <- lmer(PC1 ~ genotypeC * environmentC + TSC + (1|cage), 
                     data=data_reduced_withPC, REML = FALSE)
PC1_mod_null <- lmer(PC1 ~  1 + (1|cage), 
                    data=data_reduced_withPC, REML = FALSE)
anova(PC1_mod_final, PC1_mod_null)

# PC2
PC2_mod_final <- lmer(b ~ genotypeC * environmentC + TSC + (1|cage), 
                     data=data_reduced_withPC, REML = FALSE)
PC2_mod_null <- lmer(b ~  1 + (1|cage), 
                    data=data_reduced_withPC, REML = FALSE)
anova(PC2_mod_final, PC2_mod_null)
#####

# Supplementary Table 7: Statistical analysis of behavioral test battery without PCA
#####
# model function
model <- function(x){
  lmer(x~G*E+TS+(1|cage))
}
G <- data$genotypeC
E <- data$environmentC
TS <- data$TSC
cage <- data$cage

### EPM
# entries
x <- data$EPM.rel.OA.entries
EPM.entr <- model(x)
table1=as.data.frame(round(coef(summary(EPM.entr)),3))
table2=round(as.data.frame(anova(EPM.entr)),3)
# duration
x <- data$EPM.rel.OA.duration
EPM.dur <- model(x)
table1=as.data.frame(round(coef(summary(EPM.dur)),3))
table2=round(as.data.frame(anova(EPM.dur)),3)
# distance with BoxCox transformation
x <- data$EPM.distance                     # library(MASS)
Box = boxcox(x ~ 1,lambda = seq(-6,6,0.1)) # BoxCox transformation
Cox = data.frame(Box$x, Box$y)             # Create a data frame with the results
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]  # Order the new data frame by decreasing y
Cox2[1,]                                   # Display the lambda with the greatest log likelihood
lambda = Cox2[1, "Box.x"]                  # Extract that lambda
x.BC = (x ^ lambda - 1)/lambda             # Transform the original data
EPM.dist.BC <- model(x.BC)
table1=as.data.frame(round(coef(summary(EPM.dist.BC)),3))
table2=round(as.data.frame(anova(EPM.dist.BC)),3)

### OFT
# entries
x <- data$OFT.entries
OFT.entry <- model(x)
table1=as.data.frame(round(coef(summary(OFT.entry)),3))
table2=round(as.data.frame(anova(OFT.entry)),3)
# duration
x <- data$OFT.duration.center
OFT.dur <- model(x)
table1=as.data.frame(round(coef(summary(OFT.dur)),3))
table2=round(as.data.frame(anova(OFT.dur)),3)
# distance
x <- data$OFT.distance
OFT.dist <- model(x)
table1=as.data.frame(round(coef(summary(OFT.dist)),3))
table2=round(as.data.frame(anova(OFT.dist)),3)

### FET
# entries
x <- data$FET.entries
FET.entry <- model(x)
table1=as.data.frame(round(coef(summary(FET.entry)),3))
table2=round(as.data.frame(anova(FET.entry)),3)
# latency with log transformation
x <- data$FET.latency
FET.late.log <- model(log(x))
table1=as.data.frame(round(coef(summary(FET.late.log)),3))
table2=round(as.data.frame(anova(FET.late.log)),3)
# duration
x <- data$FET.duration
FET.dur <- model(x)
table1=as.data.frame(round(coef(summary(FET.dur)),3))
table2=round(as.data.frame(anova(FET.dur)),3)
# distance
x <- data$FET.distance
FET.dist <- model(x)
table1=as.data.frame(round(coef(summary(FET.dist)),3))
table2=round(as.data.frame(anova(FET.dist)),3)

### LM
# mistakes
x <- data$LM.rel.diff.mist
LM.mist <- model(x)
table1=as.data.frame(round(coef(summary(LM.mist)),3))
table2=round(as.data.frame(anova(LM.mist)),3)
# distance with log transformation
x <- data$LM.rel.diff.dist
LM.dist.log <- model(log(x))
table1=as.data.frame(round(coef(summary(LM.dist.log)),3))
table2=round(as.data.frame(anova(LM.dist.log)),3)
# duration with sqrt transformation
x <- data$LM.rel.diff.dur
LM.dur.sqrt <- model(sqrt(x))
table1=as.data.frame(round(coef(summary(LM.dur.sqrt)),3))
table2=round(as.data.frame(anova(LM.dur.sqrt)),3)
#####




###############
### Figures ###
###############

# Figure 5: Influence of genotype and environment on anxiety-like behavior (PC1) and spatial learning (PC2)
#####
Figure5a <- ggplot(data= data_reduced_withPC,aes(x=genotype, y=PC1)) + # library(ggplot2)
  geom_beeswarm(aes(colour = environment, shape=environment), # library(ggbeeswarm)
                dodge.width = 0.5,size=2,alpha=0.5,cex = 2.2) +
  stat_summary(aes(colour = environment),
               fun=median,geom="crossbar", 
               position=position_dodge(width=0.5), 
               width=0.3,size=0.3,fatten=3)+
  stat_summary(aes(colour = environment),
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  annotate("text", x = 1.5, y = 5, label = "Genotype: p = 0.001\\nEnvironment: p = 0.02\\nGxE: p = 0.72",
           hjust = 0.5, size = 3) +
  ylim(-3.5,6.5) +
  ylab("Anxiety (PC1)") +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),values = c("#fc7e4c", "#799df2")) + 
  scale_shape_manual(name="Environment",breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6J" = "C57", "B6D2F1N" = "F1"),expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x = element_blank(),axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = .5,vjust = .0,face = "plain",size = 12),
        strip.text = element_text(size = 12)) +
  labs(tag="a)") 

Figure5b <- ggplot(data= data_reduced_withPC,aes(x=genotype, y=PC2)) + # library(ggplot2)
  geom_beeswarm(aes(colour = environment, shape=environment), # library(ggbeeswarm)
                dodge.width = 0.5,size=2,alpha=0.5,cex = 2.2) +
  stat_summary(aes(colour = environment),fun=median,geom="crossbar", 
               position=position_dodge(width=0.5),width=0.3,size=0.3,fatten=3)+
  stat_summary(aes(colour = environment),
               fun.min = function(z) {quantile(z,0.25)},
               fun.max = function(z) {quantile(z,0.75)},
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  annotate("text", x = 1.5, y = 5, label = "Genotype: p = 0.25\\nEnvironment: p = 0.81\\nGxE: p = 0.55",
           hjust = 0.5, size = 3) +
  ylim(-3.5,6.5) +
  ylab("Spatial learning (PC2)") +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),values = c("#fc7e4c", "#799df2")) + 
  scale_shape_manual(name="Environment",breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6J" = "C57", "B6D2F1N" = "F1"),expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x = element_blank(),axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = .5,vjust = .0,face = "plain",size = 12),
        strip.text = element_text(size = 12)) +
  labs(tag="b)") 
#####

# Supplementary Figure 1: Scree plot from PCA on 13 behavioral parameters from the battery of behavioral tests 
#####
PCA2<-PCA(na.omit(data[,c(1:13)]),ncp=15, scale.unit=TRUE,graph=F) # library(FactoMineR)

screeplot_PCA2 <- fviz_eig(PCA2,addlabels = TRUE,ncp=10,linecolor ="red",
                           ggtheme=theme_classic(base_size=18)) # library(factoextra)
screeplot_PCA2+labs(title = "",x = "Principal Components",y = "% of variances explained")
#####

# Supplementary Figure 4: State anxiety.
#####
Suppl.Fig.4a <- ggplot(data=data,aes(x=genotype, y=EPM.rel.OA.entries)) +
  geom_beeswarm(aes(colour = environment, shape=environment),
                dodge.width = 0.5,size=2,alpha=0.6,cex = 2.2) +
  stat_summary(aes(colour = environment),fun=mean,geom="crossbar", 
               position=position_dodge(width=0.5), width=0.3, size=0.3)+
  stat_summary(aes(colour = environment),fun.data = mean_sd, 
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  scale_y_continuous(name="EPM: relative open arm entries",breaks=seq(0,1,0.2),limits=c(0,1)) +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),values = c("#fc7e4c", "#799df2")) +
  scale_shape_manual(name="Environment",breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6J" = "C57", "B6D2F1N" = "F1"),expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x = element_blank(),axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = .5,vjust = .0,face = "plain",size = 12),
        strip.text = element_text(size = 12))+
  annotate("text", x = 1.5, y = 0.9,
           label = "Genotype: p < 0.0001\\nEnvironment: p = 0.031\\nGxE: p = 0.016",hjust = 0.5, size = 3) +
  labs(tag="a)")

Suppl.Fig.4b <- ggplot(data=data,aes(x=genotype, y=EPM.rel.OA.duration)) +
  geom_beeswarm(aes(colour = environment, shape=environment),
                dodge.width = 0.5,size=2,alpha=0.6,cex = 2.2) +
  stat_summary(aes(colour = environment),fun=mean,geom="crossbar", 
               position=position_dodge(width=0.5), width=0.3, size=0.3)+
  stat_summary(aes(colour = environment),fun.data = mean_sd, 
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  scale_y_continuous(name="EPM: relative time on open arms",breaks=seq(0,1,0.2),limits=c(0,1)) +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),values = c("#fc7e4c", "#799df2")) +
  scale_shape_manual(name="Environment",breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6J" = "C57", "B6D2F1N" = "F1"),expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x = element_blank(),axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = .5,vjust = .0,face = "plain",size = 12),
        strip.text = element_text(size = 12))+
  annotate("text", x = 1.5, y = 0.9,
           label = "Genotype: p < 0.0001\\nEnvironment: p = 0.006\\nGxE: p = 0.016",hjust = 0.5, size = 3) +
  labs(tag="b)")

Suppl.Fig.4c <- ggplot(data=data,aes(x=genotype, y=OFT.entries)) +
  geom_beeswarm(aes(colour = environment, shape=environment),
                dodge.width = 0.5,size=2,alpha=0.6,cex = 2.2) +
  stat_summary(aes(colour = environment),fun=mean,geom="crossbar", 
               position=position_dodge(width=0.5), width=0.3, size=0.3)+
  stat_summary(aes(colour = environment),fun.data = mean_sd, 
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  scale_y_continuous(name="OFT: entries into the center zone",breaks=seq(0,25,5),limits=c(0,25)) +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),values = c("#fc7e4c", "#799df2")) +
  scale_shape_manual(name="Environment",breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6J" = "C57", "B6D2F1N" = "F1"),expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x = element_blank(),axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = .5,vjust = .0,face = "plain",size = 12),
        strip.text = element_text(size = 12))+
  annotate("text", x = 1.5, y = 22, label = "Genotype: p = 0.093\\nEnvironment: p = 0.523\\nGxE: p = 0.571",
           hjust = 0.5, size = 3) +
  labs(tag="c)")

Suppl.Fig.4d <- ggplot(data=data,aes(x=genotype, y=OFT.duration.center)) +
  geom_beeswarm(aes(colour = environment, shape=environment),
                dodge.width = 0.5,size=2,alpha=0.6,cex = 2.2) +
  stat_summary(aes(colour = environment),fun=mean,geom="crossbar", 
               position=position_dodge(width=0.5), width=0.3, size=0.3)+
  stat_summary(aes(colour = environment),fun.data = mean_sd, 
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  scale_y_continuous(name="OFT: time spent in center [s]",breaks=seq(0,40,5),limits=c(0,40)) +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),values = c("#fc7e4c", "#799df2")) +
  scale_shape_manual(name="Environment",breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6J" = "C57", "B6D2F1N" = "F1"),expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x = element_blank(),axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = .5,vjust = .0,face = "plain",size = 12),
        strip.text = element_text(size = 12))+
  annotate("text", x = 1.5, y = 36, label = "Genotype: p = 0.833\\nEnvironment: p = 0.348\\nGxE: p = 0.796",
           hjust = 0.5, size = 3) +
  labs(tag="d)")
#####

# Supplementary Figure 5: Trait anxiety
#####
Suppl.Fig.5a <- ggplot(data=data,aes(x=genotype, y=FET.latency)) +
  geom_beeswarm(aes(colour = environment, shape=environment),
                dodge.width = 0.5,size=2,alpha=0.6,cex = 2.2) +
  stat_summary(aes(colour = environment),fun=mean,geom="crossbar", 
               position=position_dodge(width=0.5), width=0.3, size=0.3)+
  stat_summary(aes(colour = environment),fun.data = mean_sd, 
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  scale_y_continuous(name="FET: latency to enter the arena [s]",breaks=seq(0,1000,150),limits=c(-100,1000)) +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),values = c("#fc7e4c", "#799df2")) +
  scale_shape_manual(name="Environment",breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6J" = "C57", "B6D2F1N" = "F1"),expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x = element_blank(),axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = .5,vjust = .0,face = "plain",size = 12),
        strip.text = element_text(size = 12))+
  annotate("text", x = 1.5, y = 850, label = "Genotype: p = 0.04\\nEnvironment: p = 0.004\\nGxE: p = 0.793",
           hjust = 0.5, size = 3) +
  labs(tag="a)")

Suppl.Fig.5b <- ggplot(data=data,aes(x=genotype, y=FET.entries)) +
  geom_beeswarm(aes(colour = environment, shape=environment),
                dodge.width = 0.5,size=2,alpha=0.6,cex = 2.2) +
  stat_summary(aes(colour = environment),fun=mean,geom="crossbar", 
               position=position_dodge(width=0.5), width=0.3, size=0.3)+
  stat_summary(aes(colour = environment),fun.data = mean_sd, 
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  scale_y_continuous(name="FET: entries into the arena",breaks=seq(0,35,5),limits=c(-5,35)) +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),values = c("#fc7e4c", "#799df2")) +
  scale_shape_manual(name="Environment",breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6J" = "C57", "B6D2F1N" = "F1"),expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x = element_blank(),axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = .5,vjust = .0,face = "plain",size = 12),
        strip.text = element_text(size = 12))+
  annotate("text", x = 1.5, y = 30, label = "Genotype: p = 0.369\\nEnvironment: p < 0.0001\\nGxE: p = 0.875",
           hjust = 0.5, size = 3) +
  labs(tag="b)")
#####

# Supplementary Figure 6: Spatial learning
#####
Suppl.Fig.6a <- ggplot(data=data,aes(x=genotype, y=LM.rel.diff.mist)) +
  geom_beeswarm(aes(colour = environment, shape=environment),
                dodge.width = 0.5,size=2,alpha=0.6,cex = 2.2) +
  stat_summary(aes(colour = environment),fun=mean,geom="crossbar", 
               position=position_dodge(width=0.5), width=0.3, size=0.3)+
  stat_summary(aes(colour = environment),fun.data = mean_sd, 
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  scale_y_continuous(name="LM: relative difference in mistakes",breaks=seq(-0.2,2.6,0.2),limits=c(-0.2,2.6)) +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),values = c("#fc7e4c", "#799df2")) +
  scale_shape_manual(name="Environment",breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6J" = "C57", "B6D2F1N" = "F1"),expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x = element_blank(),axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = .5,vjust = .0,face = "plain",size = 12),
        strip.text = element_text(size = 12))+
  annotate("text", x = 1.5, y = 2.3, label = "Genotype: p = 0.065\\nEnvironment: p = 0.910\\nGxE: p = 0.571",
           hjust = 0.5, size = 3) +
  labs(tag="a)")

Suppl.Fig.6b <- ggplot(data=data,aes(x=genotype, y=LM.rel.diff.dur)) +
  geom_beeswarm(aes(colour = environment, shape=environment),
                dodge.width = 0.5,size=2,alpha=0.6,cex = 2.2) +
  stat_summary(aes(colour = environment),fun=mean,geom="crossbar", 
               position=position_dodge(width=0.5), width=0.3, size=0.3)+
  stat_summary(aes(colour = environment),fun.data = mean_sd, 
               geom = "errorbar", width=0.2, position=position_dodge(width=0.5)) +
  scale_y_continuous(name="LM: relative difference in time to reach exit",breaks=seq(-0.2,2.6,0.2),limits=c(-0.2,2.6)) +
  xlab("Genotype") +
  scale_colour_manual(name="Environment",breaks=c("scarce", "complex"),
                      labels=c("Scarce", "Complex"),values = c("#fc7e4c", "#799df2")) +
  scale_shape_manual(name="Environment",breaks=c("scarce", "complex"),
                     labels=c("Scarce", "Complex"),values= c("scarce"=18, "complex"=20)) +
  scale_x_discrete(labels=c("C57BL/6J" = "C57", "B6D2F1N" = "F1"),expand = c(0,0.3)) +
  theme_classic(base_size=12) +
  theme(legend.position = "bottom",legend.title = element_blank(),
        axis.ticks.x = element_blank(),axis.line.x=element_blank(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(hjust = .5,vjust = .0,face = "plain",size = 12),
        strip.text = element_text(size = 12))+
  annotate("text", x = 1.5, y = 2.3, label = "Genotype: p = 0.05\\nEnvironment: p = 0.719\\nGxE: p = 0.621",
           hjust = 0.5, size = 3) +
  labs(tag="b)")
#####