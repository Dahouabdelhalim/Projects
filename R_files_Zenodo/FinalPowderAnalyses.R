##powder spread, eye score, and peck rate analyses for Hawley et al. High virulence is associated with pathogen spreadability in a songbird-bacterial system

rm(list=ls())

#load various libraries- may need to install packages
library(effects)
library(emmeans)
library(multcompView)
library(multcomp)
library(ggplot2)
library(car)
library(dplyr)
library(data.table)

#load data file for powder scores and eye scores
library(readr)
Aviary_Contact_Rate_F17_Rformat <- read_csv("DHawley/Aviary_Contact_Rate_F17_Rformat.csv")
View(Aviary_Contact_Rate_F17_Rformat)
aviary<-Aviary_Contact_Rate_F17_Rformat

##first make sure R is treating status (index bird versus flockmate) as a factor
aviary$status = as.factor(aviary$status)
unique(aviary$status)

unique (aviary$treatment)

#relevel treatment to switch low and high virulence treatment to non-alpha order, so control is baseline intercept
aviary <- mutate(aviary, 
                 treatment=factor (treatment, 
                                   levels=c("CON", "LOW", "HIGH")))
levels(aviary$treatment)

#subset to flockmates only
aviary.fm=subset(aviary,status=="FLOCKMATE")

#subset to index birds only
aviary.in=subset(aviary,status=="INDEX")

#make sure DPI (days post-inoculation) is reading correctly
unique(aviary.in$pid)

##subset to days post-inoculation DPI < 17 before the single index bird died (1170) so that eye score max is over an equivalent number of days for all index birds (first 14 days)
Aviary.InDPI=subset(aviary.in, pid < 17)

##calculate max eye score for each index bird in the first two weeks post-inoculation
Index_max = Aviary.InDPI %>% 
  group_by(treatment, band_ID) %>% 
  summarize(max_eyescore=max(eyescore_tot,na.rm=TRUE), sum_eyescore=sum(eyescore_tot,na.rm=TRUE), max_MGload=max(log_MGload,na.rm=TRUE), sum_MGload=sum(log_MGload,na.rm=TRUE))
Index_max

##summarize max and sum of eye scores and pathogen loads across treatment for index birds
library(dplyr)
group_by(Index_max, treatment) %>%
  summarise(
    mean = mean(max_eyescore, na.rm = TRUE),
    sd = sd(max_eyescore, na.rm = TRUE),
)

#now we begin data analysis

##first let's look at MAX eye score in index birds based on treatment... use Kruskal-Wallis since small sample size and data do not meet assumptions of normality, etc.
kruskal.test(max_eyescore ~ treatment, data = Index_max)

##results are similar if I use sum eye score rather than max
kruskal.test(sum_eyescore ~ treatment, data = Index_max)

##now post-hoc tests; note that I had to load FSA package manually for post-hoc tests
library(dunn.test)
dunn.test(Index_max$max_eyescore, g=Index_max$treatment, method="bonferroni", kw=TRUE, label=TRUE, wrap=FALSE, table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05, altp=FALSE)

##now look at pathogen load (max and sum)- max is not quite statistically significant at p=0.066- reported in supplement as not directly relevant to this study
kruskal.test(max_MGload ~ treatment, data = Index_max)
kruskal.test(sum_MGload ~ treatment, data = Index_max)

##change level names of Virulence Treatments for graphing purposes
levels(Index_max$treatment)<-c("Control", "Low", "High")
levels(Index_max$treatment)

##graphing raw data for eye score for manuscript- this is Figure 3
Fig1=ggplot(data=Index_max,aes(x=treatment,y=max_eyescore))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(), alpha = 0.3)+ 
  ylab("Maximum Eye Score")+ #change y axes
  xlab("Virulence Treatment")+
  theme_bw()+
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=18),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(face="italic"))
Fig1

##now the main analysis of whether powder score, the proxy for spreadability, differs by treatment; we collected powder scores twice but here we are using round 1 when all index birds were still alive
##use ordinal regression given nature of powder score data (discrete, scores ranged from 0 to 2 for flockmates)
  
  #first subset to powder round 1 
  aviary.fm1=subset(aviary.fm, powder_round==1)
  
  ##given that the powder scores only span 0,1,2 and will be treated as ordinal, need to make Total_eye a factor
  aviary.fm1$Total_eye_powder = as.factor(aviary.fm1$Total_eye_powder)
  unique(aviary.fm1$Total_eye_powder)
  
  library(ordinal)
  mod2 = clmm(Total_eye_powder~treatment + (1|group), data =aviary.fm1)
  mod3 = clmm(Total_eye_powder~treatment + sex + (1|group), data =aviary.fm1)
  summary(mod2)
  summary(mod3)
  anova(mod2,mod3) ##test whether sex (mod3) significantly improves model in LR test - it does not
  emmeans(mod2, pairwise~treatment, type="response")
  ##mod2 results are what was reported in manuscript
  

  ##change level names of Virulence Treatments for graphing purposes
  levels(aviary.fm1$treatment)<-c("Control", "Low", "High")
  levels(aviary.fm1$treatment)

  library(plyr)
  library(dplyr)

##summarize data for stacked bar chart - calculate proportion of each flockmate powder score (0, 1, or 2) for each treatment
  ##first make sure R treats Powder score as factor for graph purposes
  aviary.fm1$Total_eye_powder = as.factor(aviary.fm1$Total_eye_powder)
  unique(aviary.fm1$Total_eye_powder)
  aviary.fm1.low=subset(aviary.fm1, treatment=="LOW")
  aviary.fm1.con=subset(aviary.fm1, treatment=="CON")
  aviary.fm1.high=subset(aviary.fm1, treatment=="HIGH")

  library(dplyr)
  
  ##generate proportions for each flock for each bird
  table(aviary.fm1.low$Total_eye_powder)
  table(aviary.fm1.con$Total_eye_powder)
  table(aviary.fm1.high$Total_eye_powder)
  
  ##input data file with proportions of each score per flock
  library(readxl)
  PropBirdsPowder <- read_excel("Desktop/PropBirdsPowder.xlsx")
  View(PropBirdsPowder)
  
  #relevel treatment to switch low and high virulence treatment to non-alpha order
 PropBirdsPowder <- mutate(PropBirdsPowder, 
                   Treatment=factor (Treatment, 
                                     levels=c("Control", "Low", "High")))
  levels(PropBirdsPowder$Treatment)
  
  #make stackedbargraph- THIS WAS THE FINAL FIGURE 2 FOR POWDER SCORE
  ggplot(data=PropBirdsPowder, aes(x = Treatment, y = Proportion, fill = Score)) + 
    geom_bar(stat = "identity", color = "black") +
    scale_fill_brewer(labels = c("0", "1", "2"))+
    guides(fill = guide_legend(title = "Powder Score"))+
    ylab("Proportion Flockmates")+ #change y axes
    xlab("Virulence Treatment")+
    theme(axis.title=element_text(size=22),
          axis.text=element_text(size=16),
          panel.grid = element_blank(), 
          legend.text = element_text(size=17),
          legend.title = element_text (size=17),
          axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
          axis.line=element_line(),
          legend.box.background = element_rect(colour = "black"),
          axis.text.x = element_text(face="italic"))
 
   
#now I'm going to load the video peck data from a separate file - these video data are for examining peck rate for index birds only while at feeder
  library(readr)
  Aviary_Peck_Data_R <- read_csv("Desktop/Aviary_Peck_Data_R.csv")
  View(Aviary_Peck_Data_R)

pecks<-Aviary_Peck_Data_R
  
unique(pecks$Status)


#subset to index birds only
pecks.in=subset(pecks,Status=="INDEX")
unique(pecks.in$Treatment)

#relevel to switch low and high to non-alpha order
pecks.in <- mutate(pecks.in, 
                 Treatment=factor (Treatment, 
                                   levels=c("CONTROL", "LOW", "HIGH")))
levels(pecks.in$Treatment)

##make sure it is treating Bird as factor
pecks.in$Bird = as.factor(pecks.in$Bird)
unique(pecks.in$Bird)

##Look at data with simple histogram
PeckData <- pecks.in$Pecks_per_sec
hist(PeckData)


##first will try a linear mixed model on "Pecks per sec" even though the data prob aren't good fit for linear model
lm_pecks<-lmer(Pecks_per_sec~Treatment + (1|Bird), data = pecks.in)
summary(lm_pecks)
plot(lm_pecks)
resid(lm_pecks)
hist(resid(lm_pecks))
shapiro.test(resid(lm_pecks))
##shapiro test says that the data do not meet assumptions of normality


##Because data do not quite meet assumptions of linear models, try transforming peck data
pecks.in = pecks.in %>% 
  mutate(pecks_asin = asin(sqrt(Pecks_per_sec)))->pecks.in
pecks.in

##now redo lm with transformed data 
lm2_pecks<-lmer(pecks_asin~Treatment + (1|Bird), data = pecks.in)
summary(lm2_pecks)
Anova(lm2_pecks)
plot(lm2_pecks)
resid(lm2_pecks)
hist(resid(lm2_pecks))
shapiro.test(resid(lm2_pecks))
emmeans(lm2_pecks, pairwise~Treatment, type="response")
##S-W normality test got even more significant so transformation seemed to make things worse


###make a new variable of pecks_nonzero so that I can try a GLMM with a gamma distribution
pecks.in %>%
  mutate(pecks_nozero = Pecks_per_sec + 0.0000001)->pecks.in
hist(pecks.in$pecks_nozero)

##trying gamma distribution in a GLM - this seems like the best model to use
glm_pecks=glmer(pecks_nozero~Treatment + (1|Bird), data = pecks.in, family = Gamma)
summary(glm_pecks)
Anova(glm_pecks)
emmeans(glm_pecks, pairwise~Treatment, type="response")
resid(glm_pecks)
hist(resid(glm_pecks))
plot(glm_pecks)
##results are similar to LMM results but not breaking assumptions 

##change level names of Virulence Treatments for graphing purposes
levels(pecks.in$Treatment)<-c("Control", "Low", "High")
levels(pecks.in$Treatment)

install.packages("viridis")

##just plotting the raw data for pecks per second - this is the final figure in the Manuscript
Fig4=ggplot(data=pecks.in,aes(x=Treatment,y=Pecks_per_sec))+
  geom_boxplot()+
  geom_jitter(aes(color=Bird), width = 0.1, height = 0.1)+
  ##scale_colour_viridis(discrete = T)+
  ylab("Pecks at Food per Second")+ #change y axes
  xlab("Virulence Treatment")+
  theme_bw()+
  theme(axis.title=element_text(size=24),
        axis.text=element_text(size=16),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        axis.text.x = element_text(face="italic"),
        axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 0, l = 0)),
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.background = element_blank(),
        legend.key=element_rect(fill="white",color="white"))
Fig4


