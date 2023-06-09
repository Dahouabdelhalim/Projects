######################################################################################################
#                                                                                                    #
# Script for 'Fear of sex: Sexual conflict exposed as avoidance in a parthenogenetic invertebrate'   #
# Marcus Lee, Carlota Solano-Udina & Lars-Anders Hansson 2021                                        #
#                                                                                                    #
###################################################################################################### 

dd <- read.csv('FoS_dataset_2021.csv')
dd.sum <- read.csv('FoS_dataset_2021_sum.csv')

suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(car))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lsmeans))
suppressPackageStartupMessages(library(nlme))
suppressPackageStartupMessages(library(arm))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(effects))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Rmisc))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(glmmTMB))
suppressPackageStartupMessages(library(bbmle))

dd.sum$sex_comb[dd.sum$Treatment=='MM'] <- 'same'
dd.sum$sex_comb[dd.sum$Treatment=='FM'] <- 'opposite'
dd.sum$sex_comb[dd.sum$Treatment=='FF'] <- 'same'

#### Q1) Does speed differ between sexes interacting with the treatment? ####

speedmod <- lmer(dd.sum$Median_speed ~ dd.sum$Size + dd.sum$sex_comb * dd.sum$Sex +  (1 | dd.sum$Recording), REML=F)
summary(speedmod)
dd.sum$speedres <- residuals(speedmod)
dd.sum$sex_comb <- as.factor(dd.sum$sex_comb)
leveneTest(dd.sum$speedres~dd.sum$Sex*dd.sum$sex_comb)

# heteroscedastic data therefore use lme
dd.sum$com <- with(dd.sum, interaction(Sex,  sex_comb))
speedmod2 <- lme(Median_speed ~ Size + sex_comb * Sex, data = dd.sum, random = ~1 | Recording, weights = varIdent(com), method='ML')
summary(speedmod2)
anova.lme(speedmod2, type='marginal')
emmeans(speedmod2, list(pairwise~Sex:sex_comb), type='response')

smeans <- lsmeans(speedmod2, list(pairwise~Sex:sex_comb), adjust = "tukey")
smeans
smeans <- as.data.frame(smeans$`lsmeans of Sex, sex_comb`)
smeans

plot(resid(speedmod2),dd.sum$Median_speed)
qqmath(speedmod, id=0.05)
plot(speedmod)
shapiro.test(residuals(speedmod2))

#### Q2) Does depth differ between sexes interacting with the treatment? ####

# As the water surface and bottom of the aquarium create a bounded conitinous variable convert to a ratio
dd.sum$dperc <- dd.sum$Median_Z/750
dd.sum$depercden <- 750
head(dd.sum$dperc)

# fit a beta regression model with random effects
glm.1 <- glmmTMB(dperc ~ Sex*sex_comb + (1 | Recording), data = dd.sum, family = list(family = "beta", link = "logit"))
summary(glm.1)

# allow the phi to vary among groups
glm.2 <- update(glm.1, dispformula = ~ Sex * sex_comb)
summary(glm.2)

glm.3 <- update(glm.1, dispformula = ~ Sex )
summary(glm.3)

glm.4 <- update(glm.1, dispformula = ~ sex_comb)
summary(glm.4)

# determine which structure to use with with relation to phi 
AICtab(glm.1, glm.2, glm.3, glm.4) # Sex * sex_comb
summary(glm.2)

# check pairwise comparisons
emmeans(glm.2, list(pairwise~Sex:sex_comb), type='response')

# plot fitted vs residuals
plot(resid(glm.2) ~ fitted(glm.2))

# save that data and back transform the ratio results
dmeans2 <- emmeans(glm.2, list(pairwise~Sex:sex_comb), type='response')
dmeans2 <- as.data.frame(dmeans2$`emmeans of Sex, sex_comb`)
dmeans2$back <- c(dmeans2$response*750)
dmeans2$backres <- c(dmeans2$SE*750)
dmeans2

#### Q3) Do horizontal movements differ between sexes or treatments? ####

#### create a dataframe including the horzontal movements and the turing ratios
needed.function <- function(x, y, z, lag){
  # It takes the coordinates a single track as vectors x, y and z and a frame
  # interval (lag) to use to calculate distances. It returns the average values
  # for the given track. Normally used by external functions.
  
  ngdr <- numeric()
  hndr <- numeric()
  gd <- numeric()
  
  for (i in 1:(length(x) - lag))
  {
    x.dis <- diff(x[i:(i+lag)])
    y.dis <- diff(y[i:(i+lag)])
    z.dis <- diff(z[i:(i+lag)])
    gross.dis <- sum(sqrt(x.dis^2 + y.dis^2 + z.dis^2))
    xy.dis <- sqrt((x[i] - x[i+lag])^2 + (y[i] - y[i+lag])^2)
    net.dis <- sqrt((x[i] - x[i+lag])^2 + (y[i] - y[i+lag])^2 + (z[i] - z[i+lag])^2)
    ngdr <- c(ngdr, net.dis / gross.dis)
    hndr <- c(hndr, xy.dis / gross.dis)
    gd <- c(gd, gross.dis)
  }
  
  # Please NOTE the na.rm=T below!
  return(data.frame(NGDR=mean(ngdr, na.rm=T), HNDR=mean(hndr, na.rm=T), GD=mean(gd, na.rm=T)))
}
get.average.ngdr.hndr <- function(df){
  # It takes the original data frame as input and returns a new data frame with 
  # calculated values for each track. NOTE: Sensitive to NA values in XYZ!!!
  
  # NOTE: Uncomment lines below to interpolate missing positions beforeheand
  # df <- ddply(df, .(Recording, Track), transform,
  #             X = stats::approx(x=Frame, y=X, xout=Frame, rule=2)$y,
  #             Y = stats::approx(x=Frame, y=Y, xout=Frame, rule=2)$y,
  #             Z = stats::approx(x=Frame, y=Z, xout=Frame, rule=2)$y)
  
  new.df <- ddply(df, .(Recording, Treatment, Track, Individual, Sex, Size), summarise,
                  NGDR = needed.function(X, Y, Z, lag=6)$NGDR, # NOTE lag!!!
                  HNDR = needed.function(X, Y, Z, lag=6)$HNDR,
                  GD = needed.function(X, Y, Z, lag=6)$GD) # NOTE lag!
  
  return(new.df)
}

# create the dataset
dd2 <- get.average.ngdr.hndr(dd)

dd2$sex_comb[dd2$Treatment=='MM'] <- 'same'
dd2$sex_comb[dd2$Treatment=='FM'] <- 'opposite'
dd2$sex_comb[dd2$Treatment=='FF'] <- 'same'
dd2$sex_comb <- as.factor(dd2$sex_comb)
str(dd2)

head(dd2)
Hmod1 <- glmmTMB(HNDR ~ Sex*sex_comb +  (1 | Recording), family = list(family = "beta", link = "logit"),  data = dd2)
summary(Hmod1)

Hmod2 <- update(Hmod1, dispformula = ~ Sex * sex_comb)
summary(Hmod2)

Hmod3 <- update(Hmod1, dispformula = ~ Sex )
summary(Hmod3)

Hmod4 <- update(Hmod1, dispformula = ~ sex_comb)
summary(Hmod4)

# determine which structure to use with with relation to phi 
AICtab(Hmod1, Hmod2, Hmod3, Hmod4) # sex_comb

emmeans(Hmod4, list(pairwise~Sex:sex_comb), type='response')
summary(Hmod4)

hmeans2 <- emmeans(Hmod4, list(pairwise~Sex:sex_comb), type='response')
hmeans2 <- as.data.frame(hmeans2$`emmeans of Sex, sex_comb`)

#### Q4) Does tortuosity of the swimming path differ between sexes or treatments? ####

Nmod1 <- glmmTMB(NGDR ~ Sex*sex_comb +  (1 | Recording), family = list(family = "beta", link = "logit"),  data = dd2)
summary(Nmod1)

Nmod2 <- update(Nmod1, dispformula = ~Sex*sex_comb)
summary(Nmod2)

Nmod3 <- update(Nmod1, dispformula = ~Sex)
summary(Nmod3)

Nmod4 <- update(Nmod1, dispformula = ~sex_comb)
summary(Nmod4)

# determine which structure to use with with relation to phi
AICtab(Nmod1,Nmod2,Nmod3,Nmod4) # S
summary(Nmod4)

nmeans2 <- emmeans(Nmod4, list(pairwise~Sex:sex_comb), type='response')
nmeans2 <- as.data.frame(nmeans2$`emmeans of Sex, sex_comb`)

emmeans(Nmod4, list(pairwise~Sex:sex_comb), type='response', adjust = "tukey")


#######################################################
#                                                     #
# Figure for the manuscript                           #
#                                                     #
####################################################### 

pd <- position_dodge(0.4)


# figure for speed
sp <- ggplot(data = smeans) +
  aes(x=sex_comb, y=lsmean, fill=Sex)+
  geom_jitter(mapping = aes(x=sex_comb, y=Median_speed, col=Sex, shape = Sex),
              data = dd.sum,
              alpha=0.4,
              stat = "identity",
              position = position_jitterdodge(jitter.width= 0.1, jitter.height=0.1, dodge.width = 0.4),
              size=2,
              na.rm = FALSE,
              show.legend = NA,
              inherit.aes = TRUE)+
  geom_errorbar(aes(ymin=lsmean-2*SE, ymax=lsmean+2*SE), width=.05, position = pd) +
  theme_classic() +
  xlab('Conspecific Sex') + ylab('Average speed (mm'~s^-1*')') +
  theme(legend.position="none") +
  geom_point(position = pd, size = 4, shape = c(21,24,21,24))



sp1 <- sp + scale_fill_manual(values=c('#F8766D','#00BFC4','#F8766D', '#00BFC4')) +
  scale_colour_manual(values=c('#F8766D','#00BFC4','#F8766D', '#00BFC4')) +
  ylim(0,25) 



# figure for depth

vp <- ggplot(data = dmeans2) +
  aes(x=sex_comb, y=back, fill=Sex)+
  geom_hline(yintercept=0, linetype="dashed", color = "black", alpha=0.5)+
  geom_hline(yintercept=750, linetype="dashed", color = "black", alpha=0.5)+
  geom_jitter(mapping = aes(x=sex_comb, y=Median_Z, col=Sex, shape = Sex),
              data = dd.sum,
              alpha=0.4,
              stat = "identity",
              position = position_jitterdodge(jitter.width= 0.1, jitter.height=0.1, dodge.width = 0.4),
              size=2,
              na.rm = FALSE,
              show.legend = NA,
              inherit.aes = TRUE)+
  geom_errorbar(aes(ymin=back-2*backres, ymax=back+2*backres), width=.05, position = pd) +
  theme_classic() +
  scale_y_reverse(name = "Water Depth (mm)", breaks=c(0,150,300,450,600,750)) +
  xlab('Conspecific Sex') + ylab('Average depth (mm)') +
  theme(legend.position="none") +
  geom_point(position = pd, size = 4, shape = c(21,24,21,24))



vp1 <- vp + scale_fill_manual(values=c('#F8766D','#00BFC4','#F8766D', '#00BFC4')) +
  scale_colour_manual(values=c('#F8766D','#00BFC4','#F8766D', '#00BFC4')) 
   


# Figure for HDNR

hndrdd <- summarySE(data = dd2, measurevar = 'HNDR', groupvars = c('sex_comb','Sex'))

hp <- ggplot(data = hndrdd) +
  aes(x=sex_comb, y=HNDR, fill=Sex)+
  geom_jitter(mapping = aes(x=sex_comb, y=HNDR, col=Sex, shape = Sex),
              data = dd2,
              alpha=0.4,
              stat = "identity",
              position = position_jitterdodge(jitter.width= 0.1, jitter.height=0.1, dodge.width = 0.4),
              size=2,
              na.rm = FALSE,
              show.legend = NA,
              inherit.aes = TRUE)+
  geom_errorbar(aes(ymin=HNDR-2*se, ymax=HNDR+2*se), width=.05, position = pd) +
  theme_classic() +
  xlab('Conspecific Sex') + ylab('HNDR') +
  theme(legend.position="none") +
  geom_point(position = pd, size = 4, shape = c(21,24,21,24))
  



hp1 <- hp + scale_fill_manual(values=c('#F8766D','#00BFC4','#F8766D', '#00BFC4')) +
  scale_colour_manual(values=c('#F8766D','#00BFC4','#F8766D', '#00BFC4')) +
  ylim(0,1) 



# Figure for NGDR

ngdrdd <- summarySE(data = dd2, measurevar = 'NGDR', groupvars = c('sex_comb','Sex'))

np <- ggplot(data = ngdrdd) +
  aes(x=sex_comb, y=NGDR, fill=Sex, shape = Sex)+
  geom_jitter(mapping = aes(x=sex_comb, y=NGDR, col=Sex, shape = Sex),
              data = dd2,
              alpha=0.4,
              stat = "identity",
              position = position_jitterdodge(jitter.width= 0.1, jitter.height=0.1, dodge.width = 0.4),
              size=2,
              na.rm = FALSE,
              show.legend = NA,
              inherit.aes = TRUE) +
  geom_errorbar(aes(ymin=NGDR-2*se, ymax=NGDR+2*se), width=.05, position = pd) +
  theme_classic() +
  xlab('Conspecific Sex') + ylab('NGDR') +
  theme(legend.position = c(0.8, 0.15)) + 
  geom_point(position = pd, size = 4, shape = c(21,24,21,24)) 

np1 <- np + scale_fill_manual(values=c('#F8766D','#00BFC4','#F8766D', '#00BFC4')) +
  scale_colour_manual(values=c('#F8766D','#00BFC4','#F8766D', '#00BFC4'), guide=F) + 
  scale_shape_manual(values=c(Female=21,Male=24), labels=c(Female= 'Female',Male='Male'), guide = F)+
  ylim(0,1) +
  guides(fill = guide_legend('Sex', override.aes = list(shape = c(21,24), colour = c('#F8766D','#00BFC4'), fill = c('#F8766D','#00BFC4'))))


# plot all together

y <- c('a) ', 'b) ', 'c)  ', 'd)  ')
fig1.eps <- plot_grid(sp1, vp1, hp1, np1, nrow=2, ncol=2, labels = y, label_x = 0.15)
fig1.eps

ggsave2("Fig1.eps", plot = fig1.eps, width = 174, height = 200, units = "mm", device=cairo_ps)


#######################################################
#                                                     #
# Analysis and figures for supplementary material     #
#                                                     #
####################################################### 

# Predation analysis 

predd <- read.csv('FoS_Predation_Data_2021.csv')

median_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),sd=sd(x[[col]], na.rm=TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func,varname)
  return(data_sum)
}

predd.sum <- median_summary(predd,varname="Z",groupnames=c("Treatment", 'Individual'))
colnames(predd.sum) <- c( 'Treatment', 'Individual', 'Depth', 'sd_Depth')

t.test(predd.sum$Depth~predd.sum$Treatment)

predd.sum$dperc <- predd.sum$Depth/750
predd.sum$depercden <- 750


predmod1 <- glmmTMB(dperc ~ Treatment, family = list(family = "beta", link = "logit"),  data = predd.sum)
summary(predmod1)
pmeans1 <- as.data.frame(emmeans(predmod1, 'Treatment' , type='response', adjust = "sidak"))
pmeans1$back <- pmeans1$response*750
pmeans1$backSE <- pmeans1$SE*750
head(pmeans1)

# Predation figures 

pd <- position_dodge(0.4)
lab <- c('Control /  \\u2640\\u2640', 'Predation /  \\u2640\\u2642')

pe <- position_dodge2(0.7)

p2 <- ggplot(data=pmeans1) +
  aes(x=Treatment, y=back, fill=Treatment)+
  geom_errorbar(aes(ymin=back-backSE*2, ymax=back+backSE*2), width=.05, position = pd) +
  theme_classic() +
  xlab('Treatment') + ylab('Average Depth (mm)') +
  theme(legend.position="none") +
  geom_hline(yintercept=0, linetype="dashed", color = "black", alpha=0.5)+
  geom_hline(yintercept=750, linetype="dashed", color = "black", alpha=0.5)+
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_text(size=14))+
  geom_point( size = 6, shape=c(15,18), alpha=0.65)+
  scale_y_reverse(name = "Water Depth (mm)", breaks=c(0,150,300,450,600,750)) +
  geom_jitter(mapping = aes(x=Treatment, y=Depth),
              shape = c(15,15,15,15,15,15,18,18,18,18,18,18),
              data = predd.sum,
              colour = 'dark grey',
              stat = "identity",
              alpha = 0.5,
              position = position_jitterdodge(jitter.width= 0.3, jitter.height=0, dodge.width = 0.4),
              size=3,
              na.rm = FALSE,
              show.legend = NA,
              inherit.aes = TRUE)
p2

dmeans3 <- dmeans2[1,]
dmeans3 <- rbind(dmeans3, dmeans2[3,])
dmeans3

dmeans3$Treatment[dmeans3$sex_comb=='same'] <- 'Control'
dmeans3$Treatment[dmeans3$sex_comb=='opposite'] <- 'Predation'


p2 + ylim(750,0) +
  geom_errorbar(aes(ymin=back-backres, ymax=back+backres, col=Treatment), width=.05, alpha=0.5, position = position_nudge(x = 0.2, y = 0), data = dmeans3) +
  geom_point(mapping = aes(x=Treatment, y=back, shape=Treatment, col=Treatment), size = 6, alpha=0.8, data=dmeans3, position = position_nudge(x = 0.2, y = 0)) +
  scale_x_discrete(labels = lab)


p2 <- p2 + ylim(750,0)

dmeans3$Treatment <- factor(dmeans3$Treatment)

p3 <- ggplot()+
  geom_point(mapping = aes(x=Treatment, y=back, shape=Treatment, col=c('red', 'red')), size = 6, alpha=0.8, data=dmeans3) +
  geom_errorbar(aes(ymin=back-backres, ymax=back+backres, x=Treatment, col=c('red', 'red')), width=.05,data = dmeans3)+
  geom_hline(yintercept=0, linetype="dashed", color = "black", alpha=0.5)+
  geom_hline(yintercept=750, linetype="dashed", color = "black", alpha=0.5)+
  scale_y_reverse(name = "Water Depth (mm)", breaks=c(0,150,300,450,600,750)) +
  theme_classic() +
  theme(legend.position="none") +
  theme(axis.line.y = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(axis.text.x=element_text(face="bold", size=12))+
  theme(axis.text.y=element_blank())+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())
 

lab1 <- c('\\u2640\\u2640', '\\u2640\\u2642')

p3 <- p3 + ylim(750,0) + scale_x_discrete(labels = lab1)
 
FigS2.eps <- plot_grid(p2,p3, nrow=1, ncol=2)
FigS2.eps

#ggsave2("FigS2.eps", plot = FigS2.eps, width = 174, height = 200, units = "mm", device=cairo_ps)

# Figure to display the size and speed relationship

sp2 <- ggplot(data = dd.sum, 
              aes(x = Size,
                  y = Median_speed,
                  colour=Sex,
                  linetype=sex_comb, shape=sex_comb)) + 
  geom_smooth(method=lm) + 
  geom_point(size = 2) +
  xlab("Size (mm)") + ylab('Average Speed mm'~s^-1*'')+
  theme_classic()

FigS1.eps <- sp2 + scale_linetype_manual(values=c('solid','dotted')) + scale_shape_manual(values = c(17, 19)) + labs(linetype = "Conspecific sex") + labs(shape = "Conspecific sex")
FigS1.eps

ggsave2("FigS1.eps", plot = FigS1.eps, width = 174, height = 100, units = "mm", device=cairo_ps)

