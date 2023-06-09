## plot/analysis for digging termites (2D experiments)
## The analysis is a post-hoc test arising from an observation of the data

## packages
library(ggplot2)
library(car)
library(lme4)
library(multcomp)

## data
d <- "data of 2D tunneling experiments"

## num of NumOfTunnelFaces
ggplot(d, aes(x=Species, y=NumOfTunnelFaces))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7)+
  theme_bw() + theme(aspect.ratio = 1) + ylim(1,8) +
  stat_summary(fun.y=mean, geom="point", color="red")
tapply(d$NumOfTunnelFaces, d$species, mean)
tapply(d$NumOfTunnelFaces, d$species, median)
ggsave("2d.pdf", width=3, height = 3)


r <- glm(NumOfTunnelFaces~species, family=poisson(link = "log"), data=d)
Anova(r)
#Response: NumOfTunnelFaces
#LR Chisq Df Pr(>Chisq)   
#species   11.568  2   0.003077 **
multicomparison<-glht(r,linfct=mcp(species="Tukey"))
summary(multicomparison)
#Simultaneous Tests for General Linear Hypotheses
#Multiple Comparisons of Means: Tukey Contrasts
#Fit: glm(formula = NumOfTunnelFaces ~ species, family = poisson(link = "log"), data = d)
#Linear Hypotheses:
#  Estimate Std. Error z value Pr(>|z|)   
#Reticuli - Paraneo == 0  0.09963    0.26494   0.376  0.92468   
#Hetero - Paraneo == 0    0.70263    0.23363   3.007  0.00731 **
#  Hetero - Reticuli == 0   0.60300    0.23097   2.611  0.02447 * 

r <- glm(NumOfTunnelFaces~colony, family=poisson(link = "log"), data=d[d$species=="Paraneo",])
Anova(r)
#LR Chisq Df Pr(>Chisq)
#colony   1.3478  2     0.5097
r <- glm(NumOfTunnelFaces~colony, family=poisson(link = "log"), data=d[d$species=="Hetero",])
Anova(r)
#LR Chisq Df Pr(>Chisq)
#colony   2.0493  2     0.3589
r <- glm(NumOfTunnelFaces~colony, family=poisson(link = "log"), data=d[d$species=="Reticuli",])
Anova(r)
#LR Chisq Df Pr(>Chisq)
#colony   4.1433  2      0.126
#

## Time for excavation
library("survminer")
library(survival)
df<-survfit(Surv(hour,cens)~species, type = "kaplan-meier", data=d)
df <- transform(df, species= factor(species, levels = c("Paraneo", "Hetero", "Reticuli")))
ggsurvplot(fit = df,
           pval = F, pval.method = TRUE,
           risk.table = F, conf.int = FALSE,
           ncensor.plot = FALSE, size = 1, linetype = 1:3,
           legend.title = "Species",
           title="Spending time at the front",  xlim = c(0,24),
           xlab="Time (sec)", ggtheme = theme_bw()  + theme(aspect.ratio = 0.6))

ggsave("2Dgrowth.pdf", width=6, height=3)

