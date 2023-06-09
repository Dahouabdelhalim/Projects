########################################################################
###################  Phenomics and genomics for BYD ####################
###################             Paula Silva         ####################
########################################################################

###################
###   Figure 1  ###
###################

library(ggplot2)
library(dplyr)
library(viridis)
library(scales)
library(patchwork)

###### Load BLUEs
blue16 <- read.csv('output/blues/blues_ayn16.csv', header = T, check.names = F)
blue16 = blue16 %>% mutate(yield_red = yield/1000)
blue17 <- read.csv('output/blues/blues_ayn17.csv', header = T, check.names = F)
blue17 = blue17 %>% mutate(yield_red = yield/1000)
blue18 <- read.csv('output/blues/blues_ayn18.csv', header = T, check.names = F)
blue18 = blue18 %>% mutate(yield_red = yield/1000)
blue19 <- read.csv('output/blues/blues_ayn19.csv', header = T, check.names = F)
blue19 = blue19 %>% mutate(yield_red = yield/1000)
blue20 <- read.csv('output/blues/blues_ayn20.csv', header = T, check.names = F)
blue20 = blue20 %>% mutate(yield_red = yield/1000)

blue16 %>% group_by(trt) %>% summarise_if(is.numeric, mean)
data16 = blue16 %>% group_by(trt, entry) %>% summarise_if(is.numeric, mean)

p1=blue16 %>%
  ggplot(aes(fill=trt, x=byd)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = "2015-2016", x = "BYD severity (%)",
       y = "Density") + xlim (0,70) + ylim (0,1) +
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = c(0.7,0.8), legend.title = element_blank()) +
  geom_vline(xintercept = 5.42, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 20.8, colour = "#21908CFF", linetype = "dashed") 

blue17 %>% group_by(trt) %>% summarise_if(is.numeric, mean)
data17 = blue17 %>% group_by(trt, entry) %>% summarise_if(is.numeric, mean)

p2=blue17 %>%
  ggplot(aes(fill=trt, x=byd)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = "2016-2017", x = "BYD severity (%)",
       y = "") + xlim (0,70) + ylim (0,1) + 
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 6.06, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 27.4, colour = "#21908CFF", linetype = "dashed") 

blue18 %>% group_by(trt) %>% summarise_if(is.numeric, mean)
data18 = blue18 %>% group_by(trt, entry) %>% summarise_if(is.numeric, mean)

p3=blue18 %>%
  ggplot(aes(fill=trt, x=byd)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = "2017-2018", x = "BYD severity (%)",
       y = "") + xlim (0,70) + ylim (0,1) + 
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 6.06, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 15.4, colour = "#21908CFF", linetype = "dashed") 

blue19 %>% group_by(trt) %>% summarise_if(is.numeric, mean)
data19 = blue19 %>% group_by(trt, entry) %>% summarise_if(is.numeric, mean)

p4=blue19 %>%
  ggplot(aes(fill=trt, x=byd)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = "2018-2019", x = "BYD severity (%)",
       y = "") + xlim (0,70) + ylim (0,1) + 
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 3.2, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 6.50, colour = "#21908CFF", linetype = "dashed") 

blue20 %>% group_by(trt) %>% summarise_if(is.numeric, mean)
data20 = blue20 %>% group_by(trt, entry) %>% summarise_if(is.numeric, mean)

p5=blue20 %>%
  ggplot(aes(fill=trt, x=byd)) + 
  geom_density(position="dodge", alpha=0.5,) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = "2019-2020", x = "BYD severity (%)",
       y = "") + xlim (0,70) + ylim (0,1) + 
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 0.42, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 7.02, colour = "#21908CFF", linetype = "dashed")

blue16 %>% group_by(trt) %>% summarise_if(is.numeric, mean)

q1=blue16 %>%
  ggplot(aes(fill=trt, x=yield_red)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF"))+
  labs(title = "", x = "GY (tons/ha)",
       y = "Density") +  ylim(0,1) + xlim(1,8) +
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 3.4, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 3, colour = "#21908CFF", linetype = "dashed")

blue17 %>% group_by(trt) %>% summarise_if(is.numeric, mean)

q2=blue17 %>%
  ggplot(aes(fill=trt, x=yield_red)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF"))+
  labs(title = "", x = "GY (tons/ha)",
       y = "") +  ylim(0,1) + xlim(1,8) + 
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 6.3, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 4.4, colour = "#21908CFF", linetype = "dashed")

blue18 %>% group_by(trt) %>% summarise_if(is.numeric, mean)

q3=blue18 %>%
  ggplot(aes(fill=trt, x=yield_red)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = "", x = "GY (tons/ha)",
       y = "") +  ylim(0,1) + xlim(1,8) +
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 3.9, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 3.7, colour = "#21908CFF", linetype = "dashed")

blue19 %>% group_by(trt) %>% summarise_if(is.numeric, mean)

q4=blue19 %>%
  ggplot(aes(fill=trt, x=yield_red)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = "", x = "GY (tons/ha)",
       y = "") +  ylim(0,1) + xlim(1,8) +
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 5.6, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 5, colour = "#21908CFF", linetype = "dashed")

blue20 %>% group_by(trt) %>% summarise_if(is.numeric, mean)

q5=blue20 %>%
  ggplot(aes(fill=trt, x=yield_red)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = "", x = "GY (tons/ha)",
       y = "") + ylim(0,1) + xlim(1,8) +
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 5.9, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 4.4, colour = "#21908CFF", linetype = "dashed")

blue17 %>% group_by(trt) %>% summarise_if(is.numeric, mean)

m2=blue17 %>%
  ggplot(aes(fill=trt, x=ptht)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF"))+
  labs(title = "", x = expression('PTHT'[M]~(m)),
       y = "Density")  + ylim(0,12.5) + xlim(0.5,1.2) + 
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 0.99, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 0.94, colour = "#21908CFF", linetype = "dashed")

blue18 %>% group_by(trt) %>% summarise_if(is.numeric, mean)

m3=blue18 %>%
  ggplot(aes(fill=trt, x=ptht)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = "", x = expression('PTHT'[M]~(m)),
       y = "") + ylim(0,12.5) + xlim(0.5,1.2) +
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 0.67, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 0.65, colour = "#21908CFF", linetype = "dashed")

blue19 %>% group_by(trt) %>% summarise_if(is.numeric, mean)

m4=blue19 %>%
  ggplot(aes(fill=trt, x=ptht)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = "", x = expression('PTHT'[M]~(m)),
       y = "") + ylim(0,12.5) + xlim(0.5,1.2) +
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 0.982, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 0.97, colour = "#21908CFF", linetype = "dashed")

blue20 %>% group_by(trt) %>% summarise_if(is.numeric, mean)

m5=blue20 %>%
  ggplot(aes(fill=trt, x=ptht)) + 
  geom_density(position="dodge", alpha=0.5) +
  scale_fill_manual(values = c("purple","#21908CFF")) +
  labs(title = , x = expression('PTHT'[M]~(m)),
       y = "") + ylim(0,12.5) + xlim(0.5,1.2) +
  theme_set(theme_bw()) +
  theme(plot.title = element_text(size = 14),
        text = element_text(size = 12),
        legend.position = 'none') +
  geom_vline(xintercept = 0.83, colour = "purple", linetype = "dashed") +
  geom_vline(xintercept = 0.805, colour = "#21908CFF", linetype = "dashed")

pdf('Figure_1.pdf', width = 15, height = 12)
patchwork=p1+p2+p3+p4+p5+plot_spacer()+m2+m3+m4+m5+q1+q2+q3+q4+q5
patchwork + plot_annotation(tag_levels = 'A') +
  plot_layout(ncol = 5, nrow = 4, heights = c(1,1,1,1))
dev.off()

###################
###   Figure 2  ###
###################

h2 <- read.csv('heritability.csv', header = T, check.names = F)

x14.text <- element_text(size = 14)
x12.text <- element_text(size = 12)

pdf('Figure_2.pdf', width = 7, height = 6)
h2 %>% mutate(Trait = factor(Trait, levels=c("BYD","PTHTM","GY"))) %>%
  ggplot(aes(Season, H2, fill = factor(Trait))) +
  theme_bw() + facet_grid(Treatment ~ ., scales = "free_y") + 
  geom_bar(stat = "identity", position = position_dodge(width = 1), alpha=0.7) +
  scale_fill_manual(values = c("#FDE725FF","#21908CFF","purple"), name = "Trait", 
                    breaks = c("BYD", "GY","PTHTM"), labels = c("BYD", "GY", expression(PTHT[M]))) + 
  labs(title="Broad-sense Heritability",x = "Field Season", y = "Heritability") + 
  theme(axis.title = x14.text, axis.text =x12.text, strip.text = element_text(size = 14)) + ylim(0,1) + 
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype="dashed", 
             color = "grey50")
dev.off()

###################
###   Figure 3  ###
###################

library(rrBLUP)
library(ggplot2)
library(patchwork)

### Load numeric geno file
hmp01 <- read.delim('F_MAF0.01_Miss85_Het15-ayn.byd0.hmp.txt', header = T, check.names = F)
colnames(hmp01)
geno = t(as.matrix(hmp01[, 16:ncol(hmp01)])) # only genotypes
geno[1:15,1:5]
A = A.mat(geno, impute.method="mean", return.imputed = T)
e = eigen(A$A)

x=e$vectors[,1]
y=e$vectors[,2]
z=e$vectors[,3]

passport <- read.csv('passport.data.csv', header = T, check.names = F)
df=as.data.frame(cbind(x,y,z,passport))
df$year=as.factor(df$year)
df$group=as.factor(df$group)
df$predict=as.factor(df$predict)

cbp2 <- c('gray70', 'purple')

df <- df[order(df$group),]
p1=ggplot(df,aes(x=x, y=y)) + 
  geom_point(aes(shape=group,color=group),cex=1.8,stroke=1.3, alpha=0.5) +
  scale_shape_manual(values=c(19,17), labels=c("Breeding Line", "Cultivar")) +
  scale_color_manual(values=cbp2, labels=c("Breeding Line", "Cultivar")) +
  ggrepel::geom_text_repel(data = df, point.padding = 0.5, aes(label = ifelse(group=='cultivar',as.character(entry),''))) +
  scale_x_continuous(name = paste("PC 1 (", round(e$values[1]/sum(e$values), digits=2)*100, "%)", sep = ''), limits = c(-0.25,0.25)) +
  scale_y_continuous(name = paste("PC 2 (", round(e$values[2]/sum(e$values), digits=2)*100, "%)", sep = ''), limits = c(-0.25,0.25)) +
  theme_bw() + theme(legend.position=c(0.73,0.9), legend.title=element_blank(), legend.text = element_text(size = 14))

df <- df[order(df$predict),]
p2=ggplot(df,aes(x=x, y=y)) + 
  geom_point(aes(shape=predict,color=predict),cex=1.8,stroke=1.3, alpha=0.5) +
  scale_shape_manual(values=c(19,17), labels=c("Bdv2 -", "Bdv2 +")) +
  scale_color_manual(values=cbp3, labels=c("Bdv2 -", "Bdv2 +")) +
  ggrepel::geom_text_repel(data = df, point.padding = 0.5, aes(label = ifelse(group=='cultivar',as.character(entry),''))) +
  scale_x_continuous(name = paste("PC 1 (", round(e$values[1]/sum(e$values), digits=2)*100, "%)", sep = ''), limits = c(-0.25,0.25)) +
  scale_y_continuous(name = paste("PC 2 (", round(e$values[2]/sum(e$values), digits=2)*100, "%)", sep = ''), limits = c(-0.25,0.25)) +
  theme_bw() + theme(legend.position=c(0.81,0.9), legend.title=element_blank(), legend.text = element_text(size = 14))

df <- df[order(df$BYD_BLUE),]
p3=ggplot(df,aes(x=x, y=y)) + 
  geom_point(aes(color=BYD_BLUE),cex=1.8,stroke=1.3, alpha=0.5) +
  scale_color_gradient2(midpoint = mean(df$BYD_BLUE)-0.6, low = "#FDE725FF", mid = "#21908CFF" , high = "purple") +
  ggrepel::geom_text_repel(data = df, point.padding = 0.5, aes(label = ifelse(group=='cultivar',as.character(entry),''))) +
  scale_x_continuous(name = paste("PC 1 (", round(e$values[1]/sum(e$values), digits=2)*100, "%)", sep = ''), limits = c(-0.25,0.25)) +
  scale_y_continuous(name = paste("PC 2 (", round(e$values[2]/sum(e$values), digits=2)*100, "%)", sep = ''), limits = c(-0.25,0.25)) +
  theme_bw() + theme(legend.position=c(0.61,0.9), legend.text = element_text(size = 14), legend.direction="horizontal")

pdf('Figure_3.pdf', width = 14, height = 12)
patchwork=p1+p2+p3
patchwork + plot_annotation(tag_levels = 'A') +
  plot_layout(ncol = 3, nrow = 3, heights = c(1,1,1))
dev.off()

###################
###   Figure 4  ###
###################

#########
## GWAS
#########
if (!requireNamespace("BiocManager", quietly = TRUE))
  +     install.packages("BiocManager")
BiocManager::install(c("multtest","gplots", "LDheatmap", "genetics", "ape", "EMMREML", "scatterplot3d"))
BiocManager::install("Matrix")
library(multtest) 
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(scatterplot3d)
library(compiler)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

# read pheno file
myY = read.csv('blups.csv', header = T)
myY = myY[,-c(2)]
myY$bdv2 = as.numeric(myY$bdv2)

# read hap file
myG = hmp01
myKI = read.csv('GAPIT.Kin.VanRaden.csv', head = F) # kinship matrix
myCV = read.csv('GAPIT.PCA.csv',head = T) # principal components

# run gwas with two methods
myGAPIT = GAPIT(Y = myY, G = myG, KI = myKI, CV = myCV,
                model=c("GLM", "MLM"), Multiple_analysis=TRUE)

###### Re-doing MHT plots with GLM results
library(CMplot)
CMplot(Pmap,plot.type='m', col=c('gray20','gray50','gray80'), 
       cex.axis=1.4, cex.lab=1.8, pch=21, LOG10=T, threshold=0.000024041987884926, threshold.col='gray20', 
       signal.col=c('purple'), signal.pch=1, signal.cex = 1.2,width=25,ylim=c(0,25),
       chr.den.col=NULL, file="pdf",memo="",dpi=50,file.output=TRUE,verbose=TRUE)

###################
###   Figure 5  ###
###################

library(ggpubr)
library(tidyverse)
library(patchwork)

pal = c("#21908CFF","purple","#21908CFF","purple")
x12.text <- element_text(size = 12)
x14.text <- element_text(size = 14)
x16.text <- element_text(size = 16)

p1=ggplot(data, aes(x = Bdv2, y = pheno_value, fill = Bdv2, alpha = 0.4)) + 
  geom_boxplot() + scale_fill_manual(values = c("#440154FF","#35B779FF")) +
  facet_grid(trait ~ ., scales="free") +
  geom_point(pch = 21, position = position_jitterdodge(), color='gray30',alpha=0.4) +
  theme_bw() + theme(axis.title = x16.text, axis.text =x14.text, strip.text = element_text(size = 16)) +
  theme(legend.position="none",
        plot.title = element_text(size=12), axis.title = element_text(size=12),
        axis.text =  element_text(size=12)) +
  ylab("") + xlab('') +
  stat_compare_means(method = "anova", label.y.npc = 0.95,label.x = 1.8) +
  scale_x_discrete(breaks=c("Bdv2 -","Bdv2 +"),labels=c("Bdv2 -","Bdv2 +"))

dat=data[complete.cases(data[ ,3]),]
dat$QTL_5A <- factor(dat$QTL_5A, levels = c("ATCTTATAGG","GCACCGAGCA"),
                     labels = c("HAP 1","HAP 2"))

p2=ggplot(dat, aes(x = QTL_5A, y = pheno_value, fill = QTL_5A, alpha = 0.4)) + 
  geom_boxplot() + scale_fill_manual(values = c("#3B528BFF","#21908CFF")) +
  facet_grid(trait ~ ., scales="free") + 
  geom_point(pch = 21, position = position_jitterdodge(), color='gray30',alpha=0.4) +
  theme_bw() + theme(axis.title = x16.text, axis.text =x14.text, strip.text = element_text(size = 16)) +
  theme(legend.position="none",
        plot.title = element_text(size=12), axis.title = element_text(size=12),
        axis.text =  element_text(size=12)) +
  ylab("") + xlab('') + 
  stat_compare_means(method = "anova", label.y=55,label.x = 1.9)

dat$Bdv2 <- factor(dat$Bdv2, levels = c("Bdv2 FALSE","Bdv2 TRUE"),
                   labels = c("Bdv2 -","Bdv2 +"))

p3=ggplot(dat, aes(x = QTL_5A, y = pheno_value, fill = QTL_5A, alpha = 0.4)) + 
  geom_boxplot() + scale_fill_manual(values = c("#7AD151FF","#FDE725FF")) +
  facet_grid(trait ~ Bdv2, scales="free") + 
  geom_point(pch = 21, position = position_jitterdodge(), color='gray30',alpha=0.4) +
  theme_bw() + theme(axis.title = x14.text, axis.text =x12.text, 
                     strip.text = element_text(size = 14)) +
  theme(legend.position="none",
        plot.title = element_text(size=14), axis.title = element_text(size=14),
        axis.text =  element_text(size=12)) +
  ylab("") + xlab('')

pdf('Figure_5.pdf', width = 14, height = 7)
patchwork=p1+p2+p3
patchwork + plot_annotation(tag_levels = 'A') +
  plot_layout(ncol = 3, nrow = 1)
dev.off()

###################
###   Figure 6  ###
###################

library(ggplot2)
library(viridis)

gs <- read.csv('predictions.csv', header = T, check.names = F)

pdf('Figure_6.pdf', width = 7, height = 6)
gs %>%   
  mutate(Populations = factor(Populations, levels=c("Exclude 2016 (294, 63, 3)",
                                                    "Exclude 2017 (308, 49, 0)",
                                                    "Exclude 2018 (283, 74, 28)",
                                                    "Exclude 2019 (279, 78, 1)",
                                                    "Exclude 2020 (264, 93, 1)",
                                                    "Exclude 2016+2017 (245, 112, 3)",
                                                    "Exclude 2016+2018 (220, 137, 31)",
                                                    "Exclude 2016+2019 (216, 141, 4)",
                                                    "Exclude 2016+2020 (201, 156, 4)",
                                                    "Exclude 2017+2018 (234, 123, 28)",
                                                    "Exclude 2017+2019 (230, 127, 1)",
                                                    "Exclude 2017+2020 (215, 142, 1)",
                                                    "Exclude 2018+2019 (205, 152, 29)",
                                                    "Exclude 2018+2020 (190, 167, 29)",
                                                    "Exclude 2019+2020 (186, 171, 2)"))) %>% 
  mutate(trait = factor(trait, levels=c("BYD","BYD + Bdv2","PTHTM","PTHTM + Bdv2","GYRD","GYRD + Bdv2"))) %>%
  ggplot(aes(x=trait, y=Populations, fill=Prediction_Ability, label=Prediction_Ability)) +
  geom_tile() + xlab('Trait') + ylab('Training Population') + theme_bw() +
  scale_fill_viridis(discrete = F, alpha=0.6, begin = '0', end = '1', option = "D") +
  geom_text() +
  theme(legend.position="bottom",
        plot.title = element_text(size=12), axis.title = element_text(size=12),
        axis.text =  element_text(size=9)) + labs(fill = "Predictive Ability") +
  scale_x_discrete('Trait', breaks=c("BYD","BYD + Bdv2","PTHTM","PTHTM + Bdv2","GYRD","GYRD + Bdv2    "),
                   labels=c("BYD", "BYD + Bdv2", expression(PTHT[M]), expression(PTHT[M]~+~Bdv2), 
                            expression(GY[RD]),expression(GY[RD]~+~Bdv2))) +
  geom_vline(xintercept = c(2.5,4.5), linetype="dashed", 
             color = "grey40")
dev.off()



