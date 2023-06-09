### Load data file created by NYData.R
data <- read.csv("NYdata.csv")

### Convert NaN to NA
data$TotAnts <- ifelse(is.nan(data$TotAnts), NA, data$TotAnts)
data$TotMass <- ifelse(is.nan(data$TotMass), NA, data$TotMass)
data$site.plant <- paste(data$site, data$plant)

### Load required library
library(lme4)

### Analysis corresponding with Figure 1 and Table 1
m0 <- glmer(oviposit~site+minday+TotFem+TotAnts+(1|site.plant)+(1|ID),
            data=data, family="binomial",
            control=glmerControl(optimizer="bobyqa"))
m1 <- glmer(oviposit~site+minday+TotFem+I(TotAnts/TotFem)+(1|site.plant)+(1|ID),
            data=data, family="binomial",
            control=glmerControl(optimizer="bobyqa"))

aics <- AIC(m0,m1)
aics
aics[2,2]-aics[1,2]

### Results presented in Table 1
summary(m1)

### Analysis corresponding to Table 3
m2 <- glmer(oviposit~minday+TotFem+I(TotAnts/TotFem)*spp+
                (1|site.plant)+(1|ID),
            data=data, family="binomial", subset=site=="A",
            control=glmerControl(optimizer="bobyqa"))
m3 <- glmer(oviposit~minday+TotFem+I(TotAnts/TotFem)+
                (1|site.plant)+(1|ID),
            data=data, family="binomial", subset=site=="A",
            control=glmerControl(optimizer="bobyqa"))

anova(m2,m3)

### Results presented in Table 3
summary(m2)

### Figure 1
## Load required library
library(sciplot)

mn.day <- mean(data$minday, na.rm=TRUE)
mn.ants <- mean(data$TotAnts, na.rm=TRUE)
mn.fem <- mean(data$TotFem, na.rm=TRUE)
parms <- fixef(m1)
plot.fem <- function(x) {
    logodds <- parms["(Intercept)"]+parms["minday"]*mn.day+
        parms["I(TotAnts/TotFem)"]*mn.ants/mn.fem+parms["TotFem"]*x
  odds <- exp(logodds)
  return(odds/(1+odds))
}


X11(width=6, height=5)
par(mai=c(1,1,0.25,0.25))
lineplot.CI(x.factor=floor(TotFem), response=as.numeric(oviposit),
            fun=mean, type="p", x.cont=TRUE, data=data,
            ylab="P(Oviposition)",
            xlab="# Females", cex.lab=1.5)
curve(plot.fem(x), from=1, to=13, add=TRUE)
