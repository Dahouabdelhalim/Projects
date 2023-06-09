#========read me=======
#This file is the final code for generating the analyses and figures
#for the paper: Tang, T., Zhang, N., Bongers, F.J., Staab, M., Schuldt, A., Fornoff, F., et al. (2022). 
#Tree species and genetic diversity increase productivity via functional diversity and trophic feedbacks. 
#eLife, 11, e78703.

#Author:TingTang, the data and metadata can be find in Dryad:
#http://dx.doi.org/10.5061/dryad.gf1vhhmqx
#All the figures may adjusted by Illustrator applications, but only for format transformation or layout


#=======Environment setup=======
#the analyses were done in R4.0.5
setwd("F:/IBCAS/MyProjects/GD_SD_data/Final data and code")
options(digits=8)                     # set number of post-decimal digits in output to 8
# increased to show more digits if required
rm(list=ls())                         # clear workspace
library(readxl)

#packages for part I analyses
library(lmerTest)
library(lme4)                         # not really needed since loaded by lmerTest
library(nlme)

##packages for plotting figures
library(ggplot2)
library(ggsignif)
library(dplyr)
library(grid)
d<- read.csv("Tang_et_al_eLife_Data.csv")
## Structure of data set:
str(d)
names(d)
d$SPDIV <- factor(d$SP.div);nlevels(d$SPDIV)
d$SFDIV <- factor(d$SF.div);nlevels(d$SFDIV)
d$SPSF<-factor(d$SP.SF);nlevels(d$SPSF)
d$PLOT <- factor(d$PLOT);nlevels(d$PLOT)
d$SUBPL <- factor(d$SUBPL);nlevels(d$SUBPL)
d$SFMON <- factor(as.numeric(d$SPSF=="1.4"));nlevels(d$SFMON)
d$SFMIX <- factor(as.numeric(d$SPSF=="4.4"));nlevels(d$SFMIX)

d$div <- as.numeric(d$SPSF)
tapply(d$div,list(d$SPSF),length)
d$bm <- d$Tree_productivity
d$hbv <- d$Herbivory
d$fung <- d$Fungal_div
#transformations (after checking the distribution of residual values from linear models)
d$anghbv <- asin(sqrt(d$hbv/100))

#note: during revision, we add FDis_indiv which is functional diversity calculated 
#by the individual trait value, the calculation detail see the script:
#Tang_et_al_eLife_Rcode_FD_CWM_Calculation.R
#the descriptionof methods see the "Method" of the paper
#"FDis_mean" indicate the functional diversity calculated by seed family means

#========Fig 1=======
#concept figure, plotted by microsoft powerpoint

#=========Fig 2 and fig 2 supplement=========
#========Fig 2a=======
d.all <- d
d.all$SPSF <- as.factor(d.all$SPSF )
(p1.3<- ggplot(d.all, aes(x = SPDIV, y = bm)) +
    geom_boxplot(outlier.shape = NA,size = 1, aes(color = SFDIV))+
    geom_jitter(aes(color = SFDIV),
                size = 2, position = position_jitterdodge(0.4), alpha = .7)+
    # stat_summary(geom ="errorbar",
    #              fun.min = function(x) mean(x) - sd(x),
    #              fun.max = function(x) mean(x) + sd(x),
    #              width=0.2,size=1.5)+
    # stat_summary(geom = "point",
    #              fun = "mean",
    #              size = 5)+
    # geom_signif(y_position=c(50, 60), xmin=c(0.81, 1.81), xmax=c(1.19, 2.19),
    #             annotation=c("n.s.", "n.s."), tip_length=0.04) +
    geom_signif(y_position= 60, xmin=1, xmax=2,
                annotation=c("*"), tip_length=0.04)+
    scale_color_manual(values =c("#698D3F",
                                 #90BE6D",
                                 #"#5B8E7D"
                                 #, "#EAC81F",
                                 "#EA9010"
    ))+
    ylim(0,105)+
    #scale_shape_manual(values =c(21,16,24,17))+
    labs(title = "",
         x = "Species richness",
         y = expression('Community productivity ('~'Mg'~""~ ha^-1 ~')'))+
    guides(color = "none")+
    theme_bw()+
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      panel.grid = element_blank())
)

(p1.2<- ggplot(d.all, aes(x = SFDIV, y = bm)) +
    geom_boxplot(outlier.shape = NA)+
    # geom_signif(y_position=60, xmin=1, xmax=2,vjust=-0.2,
    #             annotation=c("n.s."))+
    #scale_shape_manual(values =c(0, 1))+
    labs(title = "",
         x = "Genetic richness",
         y = "")+
    ylim(0,60)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.margin = unit(rep(-1,4),"lines")
    )
  #+guides(color=F, shape = F)
)

library(grid)
pdf("Figure output/fig 2a.pdf", width = 6, height = 6)
subvp <- viewport(width = 0.16, height = 0.24, x=0.24, y=0.82)
p1.3
print(p1.2, vp = subvp)
dev.off()

#======figure 2b ~ fdis_mean===========
(p2.5.1<- ggplot(d.all, aes(x = SPDIV, y = FDis_mean)) +
   geom_boxplot(outlier.shape = NA,size =1,aes(color = SFDIV))+
   geom_jitter(aes(color = SFDIV),
               size = 2, position = position_jitterdodge(0.4), alpha = .7)+
   geom_signif(y_position= 3.3, xmin=1, xmax=2,
               annotation=c("***"), tip_length=0.04)+
   scale_color_manual(values =c("#698D3F",
                                "#EA9010"
   ))+
   ylim(0,4)+
   #scale_shape_manual(values =c(21,16,24,17))+
   labs(title = "",
        x = "",
        y = "Tree functional diversity")+
   
   theme_bw()+
   guides(color= "none")+
   theme(
     axis.title = element_text(size = 18),
     axis.text = element_text(size = 14),
     panel.grid = element_blank())
)

(p2.5.2<- ggplot(d.all, aes(x = SFDIV, y = FDis_mean)) +
    geom_boxplot(outlier.shape = NA)+
    #geom_signif(y_position=2.9, xmin=1, xmax=2,
    #annotation=c("+"))+
    #scale_shape_manual(values =c(0, 1))+
    labs(title = "",
         x = "Genetic richness",
         y = "")+
    ylim(0.5,3.2)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.margin = unit(rep(-1,4),"lines")
    )
  #+guides(color=F, shape = F)
)

library(grid)
pdf("Figure output/fig 2b.pdf.pdf", width = 6, height = 6)
subvp <- viewport(width = 0.16, height = 0.24, x=0.20, y=0.82)
p2.5.1
print(p2.5.2, vp = subvp)
dev.off()

#======figure 2 supplement ~ fdis_indiv===========
(p2.5.1<- ggplot(d.all, aes(x = SPDIV, y = FDis_indiv)) +
   geom_boxplot(outlier.shape = NA,size =1,aes(color = SFDIV))+
   geom_jitter(aes(color = SFDIV),
               size = 2, position = position_jitterdodge(0.4), alpha = .7)+
   geom_signif(y_position= 3.3, xmin=1, xmax=2,
               annotation=c("***"), tip_length=0.04)+
   scale_color_manual(values =c("#698D3F",
                                "#EA9010"
   ))+
   ylim(0,4)+
   #scale_shape_manual(values =c(21,16,24,17))+
   labs(title = "",
        x = "",
        y = "Tree functional diversity")+
   
   theme_bw()+
   guides(color= "none")+
   theme(
     axis.title = element_text(size = 18),
     axis.text = element_text(size = 14),
     panel.grid = element_blank())
)

(p2.5.2<- ggplot(d.all, aes(x = SFDIV, y = FDis_indiv)) +
    geom_boxplot(outlier.shape = NA)+
    #geom_signif(y_position=2.9, xmin=1, xmax=2,
    #annotation=c("+"))+
    #scale_shape_manual(values =c(0, 1))+
    labs(title = "",
         x = "Genetic richness",
         y = "")+
    ylim(0.5,3.2)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.margin = unit(rep(-1,4),"lines")
    )
  #+guides(color=F, shape = F)
)

library(grid)
pdf("Figure output/fig 2b_supplement.pdf.pdf", width = 6, height = 6)
subvp <- viewport(width = 0.16, height = 0.24, x=0.20, y=0.82)
p2.5.1
print(p2.5.2, vp = subvp)
dev.off()
#=======figure 2c ~hebivory=========

d$anghbv <- asin(sqrt(d$hbv/100))
d.all <- d

breaks.value <- c(0.1, 0.2, 0.3, 0.4, 0.5)
label.value <- (sin(breaks.value))^2
label.value <- round(label.value, digits = 2)

(p2.2.1<- ggplot(d.all, aes(x = SPDIV, y = anghbv)) +
    geom_boxplot(outlier.shape = NA,size =1, aes(color = SFDIV))+
    geom_jitter(aes(color = SFDIV),
                size = 2, position = position_jitterdodge(0.4), alpha = .7)+
    # geom_signif(y_position= 0.75, xmin=1, xmax=2,
    #             annotation=c("***"), tip_length=0.04)+
    scale_color_manual(values =c("#698D3F",
                                 "#EA9010"))+
    labs(title = "",
         x = "Species richness",
         y = "Herbivory (angular scale)")+
    #ylim(0,27)+
    scale_y_continuous(breaks = breaks.value, labels = label.value, limits = c(0.08,0.62))+
    theme_bw()+
    guides(color="none")+
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      panel.grid = element_blank())
)

(p2.2.2<- ggplot(d.all, aes(x = SFDIV, y = anghbv)) +
    geom_boxplot(outlier.shape = NA)+
    geom_signif(y_position=0.45, xmin=1, xmax=2,annotation=c("*"))+
    #scale_shape_manual(values =c(0, 1))+
    labs(title = "",
         x = "Genetic richness",
         y = "")+
    scale_y_continuous(breaks = breaks.value, labels = label.value, limits = c(0.15,0.5))+
    #ylim(0,18)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4),"lines"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
    )
  #+guides(color=F, shape = F)
)

library(grid)
pdf("Figure output/fig 2c.pdf", width = 6, height = 6)
subvp <- viewport(width = 0.16, height = 0.24, x=0.22, y=0.82)
p2.2.1
print(p2.2.2, vp = subvp)
dev.off()


#=========figure 2d ~ fungal=======
names(d)
(p2.3.1<- ggplot(d.all, aes(x = SPDIV, y =fung)) +
    geom_boxplot(outlier.shape = NA,size =1, aes(color = SFDIV))+
    geom_jitter(aes(color = SFDIV),
                size = 2, position = position_jitterdodge(0.4), alpha = .7)+
    # geom_signif(y_position= 0.75, xmin=1, xmax=2,
    #             annotation=c("***"), tip_length=0.04)+
    scale_color_manual(values =c("#698D3F",
                                 "#EA9010"))+
    labs(title = "",
         x = "Species richness",
         y = "Soil fungal diversity (Chao1 index)")+
    ylim(450,1350)+
    theme_bw()+
    guides(color="none")+
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      panel.grid = element_blank())
)

(p2.3.2<- ggplot(d.all, aes(x = SFDIV, y = fung)) +
    geom_boxplot(outlier.shape = NA)+
    #geom_signif(y_position=1100, xmin=1, xmax=2,annotation=c("**"))+
    #scale_shape_manual(values =c(0, 1))+
    labs(title = "",
         x = "Genetic richness",
         y = "")+
    ylim(450,1100)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.margin = unit(rep(-1,4),"lines")
    )
  #+guides(color=F, shape = F)
)
pdf("Figure output/fig 2d.pdf", width = 6, height = 6)
subvp <- viewport(width = 0.16, height = 0.24, x=0.225, y=0.82)
p2.3.1
print(p2.3.2, vp = subvp)
dev.off()


#=========figure 2E ~ CWM R1=======
names(d)
(p2.4.1<- ggplot(d.all, aes(x = SPDIV, y =cwm_rc1)) +
    geom_boxplot(outlier.shape = NA,size =1, aes(color = SFDIV))+
    geom_jitter(aes(color = SFDIV),
                size = 2, position = position_jitterdodge(0.4), alpha = .7)+
    # geom_signif(y_position= 0.75, xmin=1, xmax=2,
    #             annotation=c("***"), tip_length=0.04)+
    scale_color_manual(values =c("#698D3F",
                                 "#EA9010"))+
    labs(title = "",
         x = "Species richness",
         y = "Community weighted mean (RC1)")+
    ylim(-1,3.7)+
    theme_bw()+
    guides(color="none")+
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      panel.grid = element_blank())
)

(p2.4.2<- ggplot(d.all, aes(x = SFDIV, y = cwm_rc1)) +
    geom_boxplot(outlier.shape = NA)+
    #geom_signif(y_position=1100, xmin=1, xmax=2,annotation=c("**"))+
    #scale_shape_manual(values =c(0, 1))+
    labs(title = "",
         x = "Genetic richness",
         y = "")+
    ylim(-0.5,0.5)+
    #scale_y_continuous(breaks = c(-0.4,0,0.4))+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.margin = unit(rep(-1,4),"lines")
    )
  #+guides(color=F, shape = F)
)
pdf("Figure output/fig 2e.pdf", width = 6, height = 6)
subvp <- viewport(width = 0.16, height = 0.24, x=0.225, y=0.82)
p2.4.1
print(p2.4.2, vp = subvp)
dev.off()

#=========figure 2f ~ CWM R2=======
names(d)
(p2.5.1<- ggplot(d.all, aes(x = SPDIV, y =cwm_rc2)) +
    geom_boxplot(outlier.shape = NA,size =1, aes(color = SFDIV))+
    geom_jitter(aes(color = SFDIV),
                size = 2, position = position_jitterdodge(0.4), alpha = .7)+
    # geom_signif(y_position= 0.75, xmin=1, xmax=2,
    #             annotation=c("***"), tip_length=0.04)+
    scale_color_manual(values =c("#698D3F",
                                 "#EA9010"))+
    labs(title = "",
         x = "Species richness",
         y = "Community weighted mean (RC1)")+
    ylim(-1,2.8)+
    theme_bw()+
    guides(color="none")+
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      panel.grid = element_blank())
)

(p2.5.2<- ggplot(d.all, aes(x = SFDIV, y = cwm_rc2)) +
    geom_boxplot(outlier.shape = NA)+
    #geom_signif(y_position=1100, xmin=1, xmax=2,annotation=c("**"))+
    #scale_shape_manual(values =c(0, 1))+
    labs(title = "",
         x = "Genetic richness",
         y = "")+
    ylim(-0.5,0.5)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      plot.margin = unit(rep(-1,4),"lines")
    )
  #+guides(color=F, shape = F)
)
pdf("Figure output/fig 2f.pdf", width = 6, height = 6)
subvp <- viewport(width = 0.16, height = 0.24, x=0.225, y=0.82)
p2.5.1
print(p2.5.2, vp = subvp)
dev.off()


#===========figure 3==========
string(d)
#==========fig 3a===========
aov5 <- lm(terms(bm~FDis_mean+SPDIV+SFMON+SFMIX+
                   FDis_mean:(SPDIV+SFMON+SFMIX),keep.order=T),data=d)
anova(aov5)
d$pre.a <- predict(aov5, newdata = d)

( p3.1<- ggplot(data=d,aes(x = FDis_mean, y= pre.a)) +
    scale_color_manual(values =c(
      "#698D3F",
      "#EA9010"
    ))+
    geom_point(data = d,
               aes(x = FDis_mean, y= bm,color = SFDIV,shape = SPDIV )
    )+
    geom_line(aes(color = SFDIV, linetype = SPDIV ), size = 1 )+
    
    scale_shape_manual(values = c(1,19))+
    scale_linetype_manual(values = c("dashed","solid"))+
    # geom_smooth(method = "lm", se = F, aes(color =SPDIV))+
    #stat_summary(geom = "line", fun = "mean",aes(color = SPDIV ))+
    
    labs(title = "",
         x = "Tree functional diversity",
         y =  'Community productivity ('~'Mg'~""~ ha^-1 ~')')+
    # scale_x_continuous(breaks = breaksx.value, labels = labelx.value)+
    # scale_y_continuous(breaks = breaksy.value, labels = labely.value)+
    theme_bw()+
    guides(color = "none", shape = "none", linetype ="none")+
    theme(
      panel.grid = element_blank())
) 

#==========fig 3a supplement===========
aov5 <- lm(terms(bm~FDis_indiv+SPDIV+SFMON+SFMIX+
                   FDis_indiv:(SPDIV+SFMON+SFMIX),keep.order=T),data=d)
anova(aov5)
d$pre.a <- predict(aov5, newdata = d)
# Response: bm
# Df   Sum Sq  Mean Sq  F value     Pr(>F)    
# FDis        1  2746.29 2746.286 15.14400 0.00019739 ***
#   SPDIV     1   507.15  507.152  2.79662 0.09814075 .  
# SFMON       1   199.43  199.430  1.09973 0.29729851    
# SFMIX       1   540.23  540.226  2.97900 0.08798495 .  
# FDis:SPDIV  1   501.46  501.462  2.76524 0.10001565    
# FDis:SFMIX  1   181.53  181.535  1.00105 0.31989689    
# Residuals  85 15414.31  181.345                        
# ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#revision 1:
# Analysis of Variance Table
# 
# Response: bm
# Df   Sum Sq  Mean Sq  F value     Pr(>F)    
# FDis        1   977.00  977.002  5.36199 0.02301833 *  
# SPDIV       1  2248.13 2248.135 12.33825 0.00071704 ***
# SFMON       1     0.22    0.217  0.00119 0.97252433    
# SFMIX       1   364.84  364.841  2.00232 0.16075474    
# FDis:SPDIV  1   572.05  572.052  3.13954 0.08004149 .  
# FDis:SFMON  1   621.04  621.040  3.40840 0.06838755 .  
# FDis:SFMIX  1     1.59    1.587  0.00871 0.92587428    
# Residuals  84 15305.53  182.209                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

( p3.1.1<- ggplot(data=d,aes(x = FDis, y= pre.a)) +
    scale_color_manual(values =c(
      "#698D3F",
      "#EA9010"
    ))+
    geom_point(data = d,
               aes(x = FDis_indiv, y= bm,color = SFDIV,shape = SPDIV )
    )+
    geom_line(aes(color = SFDIV, linetype = SPDIV ), size = 1 )+
    
    scale_shape_manual(values = c(1,19))+
    scale_linetype_manual(values = c("dashed","solid"))+
    # geom_smooth(method = "lm", se = F, aes(color =SPDIV))+
    #stat_summary(geom = "line", fun = "mean",aes(color = SPDIV ))+
    
    labs(title = "",
         x = "Tree functional diversity",
         y =  'Community productivity ('~'Mg'~""~ ha^-1 ~')')+
    # scale_x_continuous(breaks = breaksx.value, labels = labelx.value)+
    # scale_y_continuous(breaks = breaksy.value, labels = labely.value)+
    theme_bw()+
    guides(color = "none", shape = "none", linetype ="none")+
    theme(
      panel.grid = element_blank())
) 
pdf("Figure output/fig 3a_supplement.pdf", width = 6, height = 8)
p3.1.1
dev.off()

#==========fig 3b=========
aov5 <- lm(terms(bm~anghbv+SPDIV+SFMON+SFMIX+
                   anghbv:(SPDIV+SFMON+SFMIX),keep.order=T),data=d)

anova(aov5)

d$pre.b <- predict(aov5, newdata = d)
# Response: bm
# Df   Sum Sq  Mean Sq  F value     Pr(>F)    
# anghbv        1   494.98  494.982  2.80047   0.097957 .  
# SPDIV         1  3030.70 3030.699 17.14683 8.1967e-05 ***
#   SFMON         1     2.15    2.150  0.01217   0.912433    
# SFMIX         1   797.92  797.920  4.51440   0.036549 *  
#   anghbv:SPDIV  1    10.00    9.998  0.05657   0.812588    
# anghbv:SFMON  1   907.61  907.613  5.13501   0.026019 *  
#   anghbv:SFMIX  1     0.05    0.052  0.00029   0.986354    
# Residuals    84 14846.98  176.750                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

breaks.value <- c(0.1, 0.2, 0.3, 0.4, 0.5)
label.value <- (sin(breaks.value))^2
label.value <- round(label.value, digits = 2)

( p3.2<- ggplot(data=d,aes(x = anghbv, y= pre.b)) +
    scale_color_manual(values =c(
      "#698D3F",
      "#EA9010"
    ))+
    scale_shape_manual(values = c(1,19))+
    scale_linetype_manual(values = c("dashed","solid"))+
    geom_point(data = d,
               aes(x = anghbv, y= bm,color = SFDIV,shape = SPDIV )
    )+
    geom_line(aes(color = SFDIV, linetype = SPDIV ), size = 1 )+
    
    labs(title = "",
         x = "Herbivory leaf damage (anguar scale)",
         y =  'Community productivity ('~'Mg'~""~ ha^-1 ~')')+
    scale_x_continuous(breaks = breaks.value, labels = label.value)+
    theme_bw()+
    guides(color = "none", shape = "none", linetype ="none")+
    theme(
      panel.grid = element_blank())
) 


#==========figure 3c=============
aov5 <- lm(terms(bm~fung+SPDIV+SFMON+SFMIX+
                   fung:(SPDIV+SFMON+SFMIX),keep.order=T),data=d)
anova(aov5)

# Analysis of Variance Table
# 
# Response: bm
# Df   Sum Sq  Mean Sq  F value     Pr(>F)    
# fung        1  1097.60 1097.605  6.33698 0.01372408 *  
# SPDIV       1  2592.93 2592.925 14.97016 0.00021497 ***
# SFMON       1    19.18   19.182  0.11074 0.74012682    
# SFMIX       1  1017.44 1017.438  5.87414 0.01751172 *  
# fung:SPDIV  1     6.27    6.273  0.03622 0.84952925    
# fung:SFMON  1    55.64   55.642  0.32125 0.57237123    
# fung:SFMIX  1   752.01  752.009  4.34170 0.04022899 *  
# Residuals  84 14549.32  173.206                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
d$pre.c <- predict(aov5, newdata = d)

( p3.3<- ggplot(data=d,aes(x = fung, y= pre.c)) +
    scale_color_manual(values =c(
      "#698D3F",
      "#EA9010"
    ))+
    scale_shape_manual(values = c(1,19))+
    geom_point(data = d,
               aes(x = fung, y= bm,color = SFDIV,shape = SPDIV )
    )+
    geom_line(aes(color = SFDIV, linetype = SPDIV ), size = 1 )+
    scale_linetype_manual(values = c("dashed","solid"))+
    labs(title = "",
         x = "Soil fungal diversity (Chao1 index)",
         y =  'Community productivity ('~'Mg'~""~ ha^-1 ~')')+
    scale_x_continuous(breaks = c(400, 600, 800, 1000), limits = c(390,1050))+
    theme_bw()+
    guides(color = "none", shape = "none", linetype ="none")+
    theme(
      panel.grid = element_blank())
) 

#fig 3 with 3 panels output
library(grid)
pdf("Figure output/figure 3_panel.pdf", width =6, height = 8)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
vplayout <- function(x, y)
  viewport(layout.pos.row = x, layout.pos.col = y)
print(p3.1, vp = vplayout(1, 1))
print(p3.2, vp = vplayout(1, 2))
print(p3.3, vp = vplayout(1, 3))
dev.off()

#=========SEM analyses=========
library(piecewiseSEM)
#======fig 4 and fig 5: for fdis_mean==========
#1) fig 4: global model
piece.fit.1 <- list(lm(FDis_mean ~ SP.div+SF.div,
                       data = d),
                    lm(anghbv ~ 
                         SP.div+
                         SF.div+
                         FDis_mean
                       ,
                       data = d),
                    lm(fung ~ 
                         SP.div+
                         SF.div
                       +FDis_mean,
                       data = d),
                    lm(bm ~ 
                         #SP.div+
                         SF.div+
                         anghbv+
                         fung
                       +FDis_mean
                       ,
                       data = d)
                    
)

fit.1 <- as.psem(piece.fit.1)
summary(fit.1)

#2)fig 5: multi-groups SEM

d$f.sp.div <- as.factor(d$SP.div)
multi.fit<- list(lm(FDis_mean ~ SF.div,
                    data = d),
                 lm(anghbv ~ 
                      SF.div
                    +FDis_mean
                    ,
                    data = d),
                 lm(fung ~ 
                      SF.div
                    #+FDis_mean,
                    ,data = d),
                 lm(bm ~ 
                      SF.div+
                      anghbv+
                      fung
                    +FDis_mean
                    ,
                    data = d)
                 
)

multi.1 <- as.psem(multi.fit)
summary(multi.1)
piece.fit.2 <- multigroup(multi.1, group = "SPDIV")
piece.fit.2


#==========fig 4 supplement and supplement fig 5: for fdis_indiv=======
#just replace fdis as individual level
#1) global model: fig 4 supplement
piece.fit.1 <- list(lm(FDis_indiv ~ SP.div+SF.div,
                       data = d),
                    lm(anghbv ~ 
                         SP.div+
                         SF.div+
                         FDis_indiv
                       ,
                       data = d),
                    lm(fung ~ 
                         SP.div+
                         SF.div
                       +FDis_indiv,
                       data = d),
                    lm(bm ~ 
                         #SP.div+
                         SF.div+
                         anghbv+
                         fung
                       +FDis_indiv
                       ,
                       data = d)
                    
)


fit.1 <- as.psem(piece.fit.1)
summary(fit.1)

#2) multi-groups model: fig 5 supplement
d$f.sp.div <- as.factor(d$SP.div)
multi.fit<- list(lm(FDis_indiv ~ SF.div,
                    data = d),
                 lm(anghbv ~ 
                      SF.div
                    +FDis_indiv
                    ,
                    data = d),
                 lm(fung ~ 
                      SF.div
                    #+FDis_indiv,
                    ,data = d),
                 lm(bm ~ 
                      SF.div+
                      anghbv+
                      fung
                    +FDis_indiv
                    ,
                    data = d)
                 
)



multi.1 <- as.psem(multi.fit)
summary(multi.1)
piece.fit.2 <- multigroup(multi.1, group = "SPDIV")
piece.fit.2

#note: the sem model plotted by Mircosoft powerpoint


