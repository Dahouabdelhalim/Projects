#!/usr/bin/env R
# Authors: Barbara Diez Rodriguez, Karen Kloth, Benedicte R. Albrectsen (2021/07/01)
# plot figures



#libraries ####

library(ggplot2)
library(viridisLite)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(gridGraphics)

#Datasets required ####

  #Aphid_fitness
  #Aphid_development
  #Tannin_concentrations
  #Tannin_induction
  #Aspen_lines_timebins_0-8h.csv
  #Aspen_lines_variables_0-8h_plsda.csv
  



#Figure 1 ####

#Timebins
# open timebin data
y <- read.csv("path/to/Aspen_lines_timebins_0-8h.csv")
y$wave <- factor(y$wave)

# calculate mean percentage per group per hour per waveform
outp <- numeric(0)
for (i in as.numeric(levels(y$wave))) { # for each wave
  for (j in 1:8) { #for each hour
  df1 <- numeric(0)
  df1 <- subset(y, y$wave == i)
  df2 <- numeric(0)
  df2 <- subset(df1, df1$hour == j)
  df3 <- numeric(0)
  #df3 <- round(tapply(df2$perc, df2$group, mean), 1)
  df3 <- c(round(tapply(df2$perc, df2$group, mean), 1), as.character(i), as.character(j))
  outp <- rbind(outp, df3)
  }
}
outp <- suppressWarnings(data.frame(outp))
View(outp)
colnames(outp)[6:7] <- c("wave", "hour")

####
library(reshape)
z <-melt(outp,id.vars=c("wave", "hour"),
          variable_name = "GT",
          value.name = "perc")
colnames(z)[4] <- "perc"
colnames(z)[3] <- "group"
# Plot --------------------------------
library(ggplot2)
library(RColorBrewer)

# correct data class
z$perc <- as.numeric(z$perc)
z$hour <- as.numeric(z$hour)

# remove wave 8 (pd) and wave 3 (not existing)
z <- subset(z, z$wave != 3)
z <- subset(z, z$wave != 8)
z$wave <- factor(z$wave)

######################
# subset per group
z50 <- subset(z, z$group == "GT50")
z60 <- subset(z, z$group == "GT60")
z69 <- subset(z, z$group == "GT69")
z72 <- subset(z, z$group == "GT72")
z79 <- subset(z, z$group == "GT79")

# GT50----------------------------------------
# calculate proportion (due to mean values, total is not exactly 100%)
outp <- numeric(0)
for (i in 1:8) { #for each hour
  df1 <- numeric(0)
  df1 <- subset(z50, z50$hour == i)
  df2 <- numeric(0)
  df2 <- 100*df1$perc / sum(df1$perc)
  df3 <- numeric(0)
  df3 <- cbind(df3, df2, as.character(i))
  outp <- rbind(outp, df3)
}
outp <- suppressWarnings(data.frame(outp))
df1$wave <- factor(df1$wave, levels= c("1", "2", "4", "5", "6", "7", "11"), ordered= T)
outp$wave <- rep(levels(df1$wave), 8)
zz <- numeric(0)
zz <- outp[,c(2,3,1)]
colnames(zz) <- c("Hour", "Wave", "Percentage")
zz$Activity <- ifelse(zz$Wave == 1, "non probing",
                        ifelse(zz$Wave == 2, "pathway",
                               ifelse(zz$Wave == 4, "salivation",
                                      ifelse(zz$Wave == 5, "phloem feeding",
                                             ifelse(zz$Wave == 6, "penetration difficulties",
                                                    ifelse(zz$Wave == 7, "xylem ingestion",
                                                           "Rpd"))))))
# correct data class
zz$Percentage <- as.numeric(as.character(zz$Percentage))
zz$Hour <- as.numeric(as.character(zz$Hour))
zz$Activity <- factor(zz$Activity)
zz$Activity <- factor(zz$Activity, levels = levels(zz$Activity)[c(1,3,2,5,6,4,7)])

ggplot(zz, aes(x=Hour, y=Percentage, fill=Activity)) +
  geom_area(colour="black", size=.2, alpha=.4) +
  #scale_fill_manual(values= cols) +
  scale_fill_brewer(palette="Paired", breaks=levels(zz$Activity)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank())
ggsave("GT50_timebudgetplot.pdf", width = 4.5, height = 2.5, units = "in")


# GT60----------------------------------------
# calculate proportion (due to mean values, total is not exactly 100%)
outp <- numeric(0)
for (i in 1:8) { #for each hour
  df1 <- numeric(0)
  df1 <- subset(z60, z60$hour == i)
  df2 <- numeric(0)
  df2 <- 100*df1$perc / sum(df1$perc)
  df3 <- numeric(0)
  df3 <- cbind(df3, df2, as.character(i))
  outp <- rbind(outp, df3)
}
outp <- suppressWarnings(data.frame(outp))
df1$wave <- factor(df1$wave, levels= c("1", "2", "4", "5", "6", "7", "11"), ordered= T)
outp$wave <- rep(levels(df1$wave), 8)
zz <- numeric(0)
zz <- outp[,c(2,3,1)]
colnames(zz) <- c("Hour", "Wave", "Percentage")
zz$Activity <- ifelse(zz$Wave == 1, "non probing",
                      ifelse(zz$Wave == 2, "pathway",
                             ifelse(zz$Wave == 4, "salivation",
                                    ifelse(zz$Wave == 5, "phloem feeding",
                                           ifelse(zz$Wave == 6, "penetration difficulties",
                                                  ifelse(zz$Wave == 7, "xylem ingestion",
                                                         "Rpd"))))))
# correct data class
zz$Percentage <- as.numeric(as.character(zz$Percentage))
zz$Hour <- as.numeric(as.character(zz$Hour))
zz$Activity <- factor(zz$Activity)
zz$Activity <- factor(zz$Activity, levels = levels(zz$Activity)[c(1,3,2,5,6,4,7)])

ggplot(zz, aes(x=Hour, y=Percentage, fill=Activity)) +
  geom_area(colour="black", size=.2, alpha=.4) +
  #scale_fill_manual(values= cols) +
  scale_fill_brewer(palette="Paired", breaks=levels(zz$Activity)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank())
ggsave("GT60_timebudgetplot.pdf", width = 4.5, height = 2.5, units = "in")

# GT69----------------------------------------
# calculate proportion (due to mean values, total is not exactly 100%)
outp <- numeric(0)
for (i in 1:8) { #for each hour
  df1 <- numeric(0)
  df1 <- subset(z69, z69$hour == i)
  df2 <- numeric(0)
  df2 <- 100*df1$perc / sum(df1$perc)
  df3 <- numeric(0)
  df3 <- cbind(df3, df2, as.character(i))
  outp <- rbind(outp, df3)
}
outp <- suppressWarnings(data.frame(outp))
df1$wave <- factor(df1$wave, levels= c("1", "2", "4", "5", "6", "7", "11"), ordered= T)
outp$wave <- rep(levels(df1$wave), 8)
zz <- numeric(0)
zz <- outp[,c(2,3,1)]
colnames(zz) <- c("Hour", "Wave", "Percentage")
zz$Activity <- ifelse(zz$Wave == 1, "non probing",
                      ifelse(zz$Wave == 2, "pathway",
                             ifelse(zz$Wave == 4, "salivation",
                                    ifelse(zz$Wave == 5, "phloem feeding",
                                           ifelse(zz$Wave == 6, "penetration difficulties",
                                                  ifelse(zz$Wave == 7, "xylem ingestion",
                                                         "Rpd"))))))
# correct data class
zz$Percentage <- as.numeric(as.character(zz$Percentage))
zz$Hour <- as.numeric(as.character(zz$Hour))
zz$Activity <- factor(zz$Activity)
zz$Activity <- factor(zz$Activity, levels = levels(zz$Activity)[c(1,3,2,5,6,4,7)])

ggplot(zz, aes(x=Hour, y=Percentage, fill=Activity)) +
  geom_area(colour="black", size=.2, alpha=.4) +
  #scale_fill_manual(values= cols) +
  scale_fill_brewer(palette="Paired", breaks=levels(zz$Activity)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank())
ggsave("GT69_timebudgetplot.pdf", width = 4.5, height = 2.5, units = "in")

# GT72----------------------------------------
# calculate proportion (due to mean values, total is not exactly 100%)
outp <- numeric(0)
for (i in 1:8) { #for each hour
  df1 <- numeric(0)
  df1 <- subset(z72, z72$hour == i)
  df2 <- numeric(0)
  df2 <- 100*df1$perc / sum(df1$perc)
  df3 <- numeric(0)
  df3 <- cbind(df3, df2, as.character(i))
  outp <- rbind(outp, df3)
}
outp <- suppressWarnings(data.frame(outp))
df1$wave <- factor(df1$wave, levels= c("1", "2", "4", "5", "6", "7", "11"), ordered= T)
outp$wave <- rep(levels(df1$wave), 8)
zz <- numeric(0)
zz <- outp[,c(2,3,1)]
colnames(zz) <- c("Hour", "Wave", "Percentage")
zz$Activity <- ifelse(zz$Wave == 1, "non probing",
                      ifelse(zz$Wave == 2, "pathway",
                             ifelse(zz$Wave == 4, "salivation",
                                    ifelse(zz$Wave == 5, "phloem feeding",
                                           ifelse(zz$Wave == 6, "penetration difficulties",
                                                  ifelse(zz$Wave == 7, "xylem ingestion",
                                                         "Rpd"))))))
# correct data class
zz$Percentage <- as.numeric(as.character(zz$Percentage))
zz$Hour <- as.numeric(as.character(zz$Hour))
zz$Activity <- factor(zz$Activity)
zz$Activity <- factor(zz$Activity, levels = levels(zz$Activity)[c(1,3,2,5,6,4,7)])

ggplot(zz, aes(x=Hour, y=Percentage, fill=Activity)) +
  geom_area(colour="black", size=.2, alpha=.4) +
  #scale_fill_manual(values= cols) +
  scale_fill_brewer(palette="Paired", breaks=levels(zz$Activity)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank())
ggsave("GT72_timebudgetplot.pdf", width = 4.5, height = 2.5, units = "in")

# GT79----------------------------------------
# calculate proportion (due to mean values, total is not exactly 100%)
outp <- numeric(0)
for (i in 1:8) { #for each hour
  df1 <- numeric(0)
  df1 <- subset(z79, z79$hour == i)
  df2 <- numeric(0)
  df2 <- 100*df1$perc / sum(df1$perc)
  df3 <- numeric(0)
  df3 <- cbind(df3, df2, as.character(i))
  outp <- rbind(outp, df3)
}
outp <- suppressWarnings(data.frame(outp))
df1$wave <- factor(df1$wave, levels= c("1", "2", "4", "5", "6", "7", "11"), ordered= T)
outp$wave <- rep(levels(df1$wave), 8)
zz <- numeric(0)
zz <- outp[,c(2,3,1)]
colnames(zz) <- c("Hour", "Wave", "Percentage")
zz$Activity <- ifelse(zz$Wave == 1, "non probing",
                      ifelse(zz$Wave == 2, "pathway",
                             ifelse(zz$Wave == 4, "salivation",
                                    ifelse(zz$Wave == 5, "phloem feeding",
                                           ifelse(zz$Wave == 6, "penetration difficulties",
                                                  ifelse(zz$Wave == 7, "xylem ingestion",
                                                         "Rpd"))))))
# correct data class
zz$Percentage <- as.numeric(as.character(zz$Percentage))
zz$Hour <- as.numeric(as.character(zz$Hour))
zz$Activity <- factor(zz$Activity)
zz$Activity <- factor(zz$Activity, levels = levels(zz$Activity)[c(1,3,2,5,6,4,7)])

ggplot(zz, aes(x=Hour, y=Percentage, fill=Activity)) +
  geom_area(colour="black", size=.2, alpha=.4) +
  #scale_fill_manual(values= cols) +
  scale_fill_brewer(palette="Paired", breaks=levels(zz$Activity)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(panel.border = element_blank()) +
  theme(axis.ticks.y = element_blank(), axis.ticks.x = element_blank())
ggsave("GT79_timebudgetplot.pdf", width = 4.5, height = 2.5, units = "in")

#PLSDA
library(mixOmics)
# Reference: http://mixomics.org/methods/pls-da/

###########################

x <- read.csv("path/to/Aspen_lines_variables_0-8h_plsda.csv")
x[is.na(x)] <- 0
ncol(x)

#Names
x$tannin <- factor(x$tannin)
levels(x$tannin) <- c("High tannin", "Low tannin")
colnames(x) <- c("tannin",   "id_rec",   "file",     "group",    
                 "non probing (total)",   "pathway (total)",
                 "pathway (mean)",   "rate.pd",  "salivation (latency)", "salivation (total)",   
                 "salivation (mean)", 
                 "salivation (max)", "sum.e2",   "sum.e2s",  "mean.e2",  "max.e2",   
                 "penetration difficulties (total)",    
                 "xylem ingestion (total)", "Rpd (total)",  "mean.w11")

# Define X (EPG variables (standard +1 to avoid zero values)) and Y (tannin class)
X <- x[,5:ncol(x)]+1
Y <- x[,1]

# Preliminary PLSDA model
splsda.x <- splsda(X, Y, ncomp = 5, logratio = 'CLR')

set.seed(2543) # for reproducibility here, only when the `cpus' argument is not used
perf.plsda <- perf(splsda.x, validation = "Mfold", folds = 5, 
                   progressBar = FALSE, auc = TRUE, nrepeat = 10) 
plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

# Chose smallest number of components with lowest error rate (ncomp = 2)

list.keepX <- c(10, 10) #= minimum number of variables to use for each component
set.seed(2543) # random seed
tune.splsda <- tune.splsda(X, Y, ncomp = 2, validation = 'Mfold', folds = 5, logratio = 'CLR',
                           progressBar = FALSE, dist = 'max.dist',
                           test.keepX = list.keepX, nrepeat = 50) # tuning parameters: number of components, number of variables from each components
choice.ncomp <- tune.splsda$choice.ncomp$ncomp
choice.keepX <- tune.splsda$choice.keepX[1:choice.ncomp]

# Build PLSDA model again using optimized parameters
splsda.RT.res <- splsda(X, Y, ncomp = choice.ncomp, keepX = choice.keepX, logratio = 'CLR') 

# plot
plotIndiv(splsda.RT.res, comp=c(1,2),ind.names = FALSE, legend = TRUE, ellipse = TRUE, title = 'sPLS-DA')
plotLoadings(splsda.RT.res, comp = 1, title = 'Loadings on comp 1', contrib = 'max', method = 'mean')
plotLoadings(splsda.RT.res, comp = 2, title = 'Loadings on comp 2', contrib = 'max', method = 'mean')

## Change colours
# Select colours:
library(RColorBrewer)
display.brewer.all()
display.brewer.pal(10,"Paired")
brewer.pal(10,"Paired")

# Plot
x <- plotLoadings(splsda.RT.res, comp = 1, title = 'Loadings on comp 1', 
             contrib = 'max', method = 'mean')
colnames(x$df)
x$df[,8] <- c(rep("#A6CEE3", 36), rep("#FDBF6F", 21))
plot(x)

#Figure 2 ####

#create boxplot with number of aphids in field and greenhouse conditions

aphids <- read.csv("path/to/Aphid_fitnes.csv", header = TRUE) #aphid fitness data
aphids$GT <- as.factor(aphids$GT)
aphids <- subset(aphids, aphids$aphid_treatment != "control")
aphids$env <- factor(aphids$env,levels= c("lab", "field"))

box1 <- ggplot(aphids, aes(x=GT, y=total_aphids, fill = env)) + geom_boxplot(position=position_dodge(0.7), width=0.5) +
  labs(title="", 
       x="SwAsp genotype",
       y="Total number aphids",
       tag = "B") + 
  theme_classic(base_size = 15) + 
  theme(axis.title.x = element_text(margin = margin(t=20), size = 15)) +
  theme(axis.title.y = element_text(margin = margin(r=20), size = 15)) +
  theme(plot.title = element_text(margin = margin(b=25))) +
  theme(aspect.ratio = 1) +
  scale_fill_manual(name   = 'Environment',
                    breaks = c('lab', 'field'), 
                    values = c("#E68226", "Darkcyan"),
                    labels = c('Greenhouse', 'Field')) +
  coord_cartesian(ylim = c(-3, max(aphids$total_aphids) + 3)) +
  theme(legend.position = "top") 
box1





#Create dataset for growth curve plot

x <- read.csv("path/to/Aphid_development.csv", header = TRUE) #growth curves dataset

x$GT <- factor(x$GT)
x$day <- factor(x$day)
x$age <- factor(x$age)

## calculate mean and standard error per genotype
gt=levels(x$GT) 
dy=levels(x$days_offspring)
ag=levels(x$age)

#mean
df=numeric(0)
df2=numeric(0)
df3=numeric(0)
mn=matrix(nr=0,nc=1)

for (i in gt) {		
  df <- subset(x,x$GT==i)
  for (j in ag) {
    df2 <- subset(df,df$age==j)
    df3 <- mean(df2[,5]) #col5=nymphs 
    mn <- c(mn,df3[1])
  }
}

# standard error
df=numeric(0)
df2=numeric(0)
df3=numeric(0)
df4=numeric(0)
se=matrix(nr=0,nc=1)
rep=matrix(nr=0,nc=1)

for (i in gt) {		
  df <- subset(x,x$GT==i)
  for (j in ag) {
    df2 <- subset(df,df$age==j)
    df3 <- sd(df2[,5]) #col5=nymphs 
    df4 <- df3/sqrt(nrow(df2))
    se <- c(se,df4)
    rep <- c(rep,nrow(df2))
  }
}

# merge in dataframe
z <- data.frame(mn)
z$se <- se
z$rep <- rep
maxage <- max(as.numeric(levels(x$age)[x$age]))
z$GT <- c(rep(gt[1],maxage),rep(gt[2],maxage),rep(gt[3],maxage),rep(gt[4],maxage),rep(gt[5],maxage))
z$age <- rep(1:maxage,5)
z
z <- z[,c(4,5,3,1,2)]

z2 <- subset(z,z$rep>=3)
z2 <- subset(z2,z2$age>=10)
z2 <- subset(z2,z2$age<=22)
z2$GT <- as.factor(z2$GT)
z2$age <- as.factor(z2$age)

#create plot 

curve <- ggplot(z2, aes(x=age, y=mn, group=GT)) +
  geom_line(aes(linetype=GT), size = 1, color = "black") +
  geom_point(aes(shape = GT), size = 3, color = "black") +
  labs(title="", 
       x="Days",
       y="Number of nymphs",
       tag = "A") + 
  theme_classic(base_size = 15) + 
  theme(axis.title.x = element_text(margin = margin(t=20))) +
  theme(axis.title.y = element_text(margin = margin(r=20))) +
  theme(plot.title = element_text(margin = margin(b=25))) +
  scale_x_discrete(breaks=c("10","15","20")) +
  theme(axis.text.x = element_text(size= 15)) +
  theme(axis.text.y = element_text(size= 15))  +
  theme(aspect.ratio = 1) +
  theme(legend.position = "top") +
  scale_shape_manual(name = "SwAsp genotype",
                     values = c(16,17,15,3,7)) +
  guides(linetype = "none")
curve

fig2 <- curve | box1

ggsave("fig2.pdf", width = 30, height = 20, units = "cm", device='pdf', dpi=600)
ggsave("fig2.tiff", width = 30, height = 20, units = "cm", device='tiff', dpi=600)

#Figure 3 ####

#tannins <- read.csv2("path/to/Tannin_concentrations.csv, header = TRUE, dec = ".") #tannin induction data with aphid numbers
tannins$GT <- as.factor(tannins$GT)
g <- subset(tannins, tannins$env == "lab" & tannins$induction == "local") 
f <- subset(tannins, tannins$env == "field" & tannins$induction == "local")



#Tannins and aphids in field conditions

t_field <- ggplot(f, aes(x=induced_tannins, y=total_aphids, group=GT, color = GT)) +
  geom_point(aes(shape = GT), size = 3)  +
  labs(title="", 
       x="\\u0394 Condensed tannins mg/g DW",
       y="Number of aphids",
       tag = "C",
       color = "GT") + 
  theme_classic(base_size = 15) + 
  theme(axis.title.x = element_text(margin = margin(t=20))) +
  theme(axis.title.y = element_text(margin = margin(r=20))) +
  theme(plot.title = element_text(margin = margin(b=25))) +
  theme(aspect.ratio = 1) +
  geom_smooth(method='lm', se = FALSE) +
  scale_color_manual(name = "SwAsp \\ngenotype",
                     values=c("#B35806", "#ED7307","#F4B967", "#58298E")) +
  scale_shape_manual(name ="SwAsp \\ngenotype",
                     values = c(16,17,18,3)) +
  theme(axis.text.x = element_text(size= 15)) +
  theme(axis.text.y = element_text(size= 15)) + ylim(0,70) + theme(legend.position = "none") 
t_field

#Tannins and aphids in greenhouse conditions

t_lab <- ggplot(g, aes(x=induced_tannins, y=total_aphids, group=GT, color = GT)) +
  geom_point(aes(shape = GT), size = 3)  +
  labs(title="", 
       x="\\u0394 Condensed tannins mg/g DW",
       y="Number of aphids",
       tag ="B",
       color = "GT") + 
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_text(margin = margin(t=20))) +
  theme(axis.title.y = element_text(margin = margin(r=20))) +
  theme(plot.title = element_text(margin = margin(b=25))) +
  geom_smooth(method='lm', se = FALSE) +
  theme(aspect.ratio = 1) +
  scale_color_manual(name= "SwAsp \\ngenotype",
                     values = c("#B35806", "#ED7307", "#9287C0", "#58298E", "#31174F")) +
  scale_shape_manual(name ="SwAsp \\ngenotype",
                     values = c(16,17,15,3,7)) +
  theme(axis.text.x = element_text(size= 15)) +
  theme(axis.text.y = element_text(size= 15)) +
  ylim(0,170) + theme(legend.position = "none")
t_lab


#seasonal tannins ###


season <- read.csv("path/to/Tannin_seasonal_changes.csv", header = T) #seasonal tannins
season$GT <- as.factor(season$GT)
season$N_treatment <- as.factor(season$N_treatment)
season$month <- factor(season$month, levels = c("july", "aug", "sep"))

p2 <- ggplot(season, aes(x=month, y=tannin_concentration, color = GT, 
                         group = interaction(N_treatment, GT))) +
  geom_line(aes(linetype = N_treatment), size = 1) +
  geom_point(aes(shape = GT), size = 3) +
  labs(title="", 
       x="",
       y="Condensed tannins mg/g DW",
       tag = "A",
       color = "SwAsp genotype",
       linetype = "Treatment") + 
  theme_classic(base_size = 15) +
  theme(aspect.ratio = 2/3) + 
  theme(axis.title.x = element_text(margin = margin(t=20))) +
  theme(axis.title.y = element_text(margin = margin(r=20))) +
  theme(plot.title = element_text(margin = margin(b=25))) +
  scale_x_discrete(labels= c("July", "August", "September")) +
  theme(axis.text.x = element_text(size= 15)) +
  theme(axis.text.y = element_text(size= 15))  +
  scale_color_manual(values = c("#B35806", "#ED7307","#F4B967", "#58298E")) +
  scale_shape_manual(name ="SwAsp genotype",
                     values = c(16,17,18,3)) +
  scale_linetype_manual(values = c("solid", "dotted"),
                        labels= c("Not fertilized", "Fertilized")) +
  theme(legend.position = "none", legend.direction="vertical") + guides(shape = "none", color = "none")
p2

#Create mock plots to extract full common legends
tannins$GT <-as.factor(tannins$GT)
p1 <- ggplot(tannins, aes(x=induced_tannins, y=total_aphids, group=GT, color = GT)) +
  geom_point(aes(shape = GT), size = 3)  +
  theme_classic(base_size = 25) +
  theme(aspect.ratio = 1) +
  scale_color_manual(name= "SwAsp \\ngenotype",
                     values = c("#B35806", "#ED7307","#F4B967", "#9287C0", "#58298E", "#31174F")) +
  scale_shape_manual(name ="SwAsp \\ngenotype",
                     values = c(16,17,18,15,3,7)) +
  theme(legend.position = "top", legend.direction="vertical") +
  guides(color=guide_legend(ncol=2))
p1

p5 <- ggplot(season, aes(x=month, y=tannin_concentration, color = GT, 
                         group = interaction(N_treatment, GT))) +
  geom_line(aes(linetype = N_treatment), size = 1) +
  geom_point(aes(shape = GT), size = 3) +
  labs(title="", 
       x="",
       y="Condensed tannins mg/g DW",
       tag = "A",
       color = "SwAsp genotype",
       linetype = "Treatment") + 
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 2/3) + 
  theme(axis.title.x = element_text(margin = margin(t=20))) +
  theme(axis.title.y = element_text(margin = margin(r=20))) +
  theme(plot.title = element_text(margin = margin(b=25))) +
  scale_x_discrete(labels= c("July", "August", "September")) +
  theme(axis.text.x = element_text(size= 15)) +
  theme(axis.text.y = element_text(size= 15))  +
  scale_color_manual(values = c("#B35806", "#ED7307","#F4B967", "#58298E")) +
  scale_shape_manual(name ="SwAsp genotype",
                     values = c(16,17,18,3)) +
  scale_linetype_manual(values = c("solid", "dotted"),
                        labels= c("Not fertilized", "Fertilized")) +
  theme(legend.position = "right", legend.direction="vertical") + guides(shape = "none", color = "none")
p5


leg1 <- get_legend(p1) #extract legend from mock plot p1
as_ggplot(leg1) 
leg2 <- get_legend(p5) #extract legend from mock plot p5
as_ggplot(leg2)

# create common legend
legend <- grid.arrange(leg1, ncol = 1)

#create final plot

fig3 <- p2 + wrap_elements(panel = legend, clip = T) + t_lab + t_field + plot_layout(ncol = 2, nrow = 2)

tiff("fig3_tannins_aphids_last_difflimit.tiff", width = 7000, height = 7000, res = 600) #save as tiff
fig3
dev.off()

ggsave("fig3.pdf", width = 25, height = 25, units = "cm", device=cairo_pdf, dpi=600) #save as pdf


#Figure S5

tannins$GT <- as.factor(tannins$GT)
tannins$induction_type <- as.factor(tannins$induction_type)
ind <- subset(tannins, tannins$env == "lab")

box2 <- ggplot(ind, aes(x=GT, y=tannin_concentration, fill = induction_type, alpha = 0.5)) + 
  stat_summary(
    fun = median,
    geom = 'line',
    size = 1,
    aes(group = induction_type, colour = induction_type),
    position = position_dodge(width = 0.50) #this has to be added
  ) + guides(color = "none", alpha = "none") +
  scale_color_manual(name   = 'Type of induction',
                     breaks = c('control', 'local', 'systemic'), 
                     values = c("#b2e880","#00a299", "#006c92"),
                     labels = c('Control', 'Local', 'Systemic')) +
  geom_boxplot(position=position_dodge(0.7),width=0.5) +
  labs(title="", 
       x="SwAsp genotype",
       y="Condensed tannins mg/g DW") + 
  theme_classic(base_size = 15) + 
  theme(axis.title.x = element_text(margin = margin(t=20), size = 15)) +
  theme(axis.title.y = element_text(margin = margin(r=20), size = 15)) +
  theme(plot.title = element_text(margin = margin(b=25))) +
  #theme(aspect.ratio = 1) +
  scale_fill_manual(name   = 'Type of induction',
                    breaks = c('control', 'local', 'systemic'), 
                    values = c("#b2e880","#00a299", "#006c92"),
                    labels = c('Control', 'Local', 'Systemic')) +
  theme(legend.position = "top") +
  theme(legend.key.size = unit(1.5, 'cm')) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.5)))
box2
ggsave("figS5.pdf", width = 25, height = 25, units = "cm", device='pdf', dpi=600)
