#Code for Samojedny et al. "Specific leaf area is lower on 
#ultramafic than on neighbouring non-ultramafic soils"

library(car) #for Anova function
library(multcomp) #preplanned contrasts
library(dplyr) #data restructure
library(relaimpo) # % variance explained by model predictors
library(ggeffects) #for predicting model effects for plotting

setwd("~/Dropbox/serpentine_sla")
df <- read.csv("data_sla_means.csv",header=T)

#####################################
# Q1: Is there an effect of soil on SLA, 
# and does this differ across regions?
#####################################

# ANOVA procedure on original data
aov.org <- aov(SLA ~ Soil * Region, data = df,
  contrasts = list(Soil = 'contr.sum',Region = 'contr.sum'))
Anova(aov.org, type = 'III')
shapiro.test(residuals(aov.org)) #normality assumption broken
leveneTest(SLA ~ Soil * Region, df) #homogeneity assumption broken

# ANOVA procedure on rank-transformed data
aov.rnk <- aov(rank(SLA) ~ Soil * Region, data = df,
  contrasts = list(Soil = 'contr.sum',Region = 'contr.sum'))
Anova(aov.rnk, type = 'III')

# examine diagnostic plots on rank-transformed data
par(mfrow=c(2,2))
plot(aov.rnk, main = "Rank-Transformed") #looks okay

#pre-planned contrasts 
df$soil.region <- interaction(df$Soil, df$Region) 
mod <- aov(rank(SLA) ~ soil.region, data = df, contrasts = list(soil.region = 'contr.sum'))
Ca_U_NU <- c(-1,1,0,0,0,0,0,0,0,0) 
Co_U_NU <- c(0,0,-1,1,0,0,0,0,0,0)
Le_U_NU <- c(0,0,0,0,-1,1,0,0,0,0)
Pr_U_NU <- c(0,0,0,0,0,0,-1,1,0,0)
So_U_NU <- c(0,0,0,0,0,0,0,0,-1,1)
cntrMat <- rbind(Ca_U_NU,Co_U_NU,Le_U_NU,Pr_U_NU,So_U_NU)
summary(glht(mod, linfct=mcp(soil.region=cntrMat), alternative="two.sided"), test=adjusted("BH"))

#####################################
# Q2: Is there an effect of soil, clim,
# and their interaction on SLA?
#####################################

df.sum <- df %>% group_by(Region, Soil) %>%
  summarise(
    n = n(),
    mean = mean(SLA),
    sd = sd(SLA),
    Q1 = quantile(SLA,probs=.25),
    median = median(SLA),
    Q3 = quantile(SLA,.75),
    min = min(SLA),
    max= max(SLA),
    se = sd / sqrt(n)
  )
df.sum<- as.data.frame(df.sum)

#append climate data
df.clim <- read.csv("data_clim_region_means.csv")
df.sum <- merge(df.sum, df.clim, by=c("Region","Soil"))

mod<-lm(mean ~ AnnPrecip*Soil, dat=df.sum)
Anova(mod)
par(mfrow = c(2, 2))
plot(mod)
shapiro.test(mod$residuals) #normality met
calc.relimp(mod) # % variance explained by model predictors

mod<-lm(mean ~ AnnTemp*Soil, dat=df.sum)
Anova(mod)
plot(mod)
shapiro.test(mod$residuals) #normality met
calc.relimp(mod) # % variance explained by model predictors

summary(lm(AnnTemp ~ AnnPrecip, dat=df.sum))
plot(AnnTemp ~ AnnPrecip, dat=df.sum)

#####################################
# Figure 2
#####################################
df.sum$symbol <- c(21,21,22,22,23,23,24,24,25,25) #for making figure2
df.sum$col <- c('#737373','#C9B826','#737373','#C9B826','#737373','#C9B826','#737373','#C9B826','#737373','#C9B826')

tiff(file="Fig2.tiff", width=9,height=4, units = 'in', res = 600)
par(mfrow=c(1,2), mar=c(5,4,2,4))

# x axis: Annual Temp
mod<-lm(mean ~ AnnTemp*Soil, dat=df.sum)
par(mar=c(5,10,2,0))
plot(mean ~ AnnTemp, dat=df.sum,  col=ifelse(df.sum$Soil == "NU",'#737373','#C9B826'),
     pch=df.sum$symbol, bg=df.sum$col, ylim=c(0,40),ylab="", xlab="")
mtext("SLA",side=2,line=2.2)
mtext("Annual Temperature (C)",side=1,line=2.2)
tmp<-ggpredict(mod, c("AnnTemp","Soil"))
tmp <- as.data.frame(tmp)
tmp.NU <- subset(tmp, group=="NU")
tmp.U <- subset(tmp, group=="U")
lines(tmp.NU$predicted~tmp.NU$x, col='#737373', lwd=2, lty=1)
polygon(c(tmp.NU$x, rev(tmp.NU$x)),c(tmp.NU$conf.high,rev(tmp.NU$conf.low)),col='#73737326',border=NA)
lines(tmp.U$predicted~tmp.U$x, col='#C9B826', lwd=2, lty=1)
polygon(c(tmp.U$x, rev(tmp.U$x)),c(tmp.U$conf.high,rev(tmp.U$conf.low)),col='#C9B82626',border=NA)
arrows(df.sum$AnnTemp,df.sum$mean-df.sum$se,df.sum$AnnTemp,df.sum$mean+df.sum$se, col="black",code=3, length=0.0, angle = 90, lwd=0.3,cex=0.5)

# x axis: Annual Precip
mod<-lm(mean ~ AnnPrecip*Soil, dat=df.sum)
par(mar=c(5,4,2,6))
plot(mean ~ AnnPrecip, dat=df.sum,  col=ifelse(df.sum$Soil == "NU",'#737373','#C9B826'),
     pch=df.sum$symbol, bg=df.sum$col,ylim=c(0,40),ylab="", xlab="")
mtext("SLA",side=2,line=2.2)
mtext("Annual Precipitation (mm)",side=1,line=2.2)
tmp<-ggpredict(mod, c("AnnPrecip","Soil"))
tmp <- as.data.frame(tmp)
tmp.NU <- subset(tmp, group=="NU")
tmp.U <- subset(tmp, group=="U")
lines(tmp.NU$predicted~tmp.NU$x, col='#737373', lwd=2, lty=1)
polygon(c(tmp.NU$x, rev(tmp.NU$x)),c(tmp.NU$conf.high,rev(tmp.NU$conf.low)),col='#73737324',border=NA)
lines(tmp.U$predicted~tmp.U$x, col='#C9B826', lwd=2, lty=1)
polygon(c(tmp.U$x, rev(tmp.U$x)),c(tmp.U$conf.high,rev(tmp.U$conf.low)),col='#C9B82624',border=NA)
arrows(df.sum$AnnPrecip,df.sum$mean-df.sum$se,df.sum$AnnPrecip,df.sum$mean+df.sum$se, col="black",code=3, length=0.0, angle = 90, lwd=0.3,cex=0.5)
legend("topleft", c("NU", "U","","Lesbos", "South Africa", "California", "Costa Rica","Puerto Rico"), pch=c(19,19,19,23,25,21, 22, 24), col=c('#737373','#C9B826',"white","darkgrey","darkgrey","darkgrey","darkgrey","darkgrey"),
       inset=c(1,0), xpd=TRUE, bty="n", cex=0.9
)
dev.off()

