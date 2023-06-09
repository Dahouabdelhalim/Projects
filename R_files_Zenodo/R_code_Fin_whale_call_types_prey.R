setwd("E:/COPY PhD AZORES/2_TEMPORAL PRESENCE & PREY/2_ACOUSTICS AND PREY/FIN WHALE ANALYSIS/R")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(MASS)
library(datasets)
library(RColorBrewer)
library(tibble)
library(ggpubr)
library(glmm)
library(dsmodels)
library(interactions)
library(ggthemes)
library(modEvA)
library(jtools)
##FIGURES ========================================================================================
prey<-read.csv(file="all_prey.csv", header=TRUE, dec= '.')
prey$month <- factor(prey$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", 
                                            "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
prey$Year<-as.factor(prey$Year)

colors<-c("green", "black", "red","yellow", "grey")
df <-data.frame(x = factor(c("Jan", "Feb", "Mar", "Apr", "May", 
                             "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
                ,colour = factor(c("1", "1", "2", "2", "2", 
                                   "3", "3", "3", "4", "4", "4", "1")))

Fig2A<-ggplot(prey, aes(x=month,  y=zoo, fill=Year)) + 
  geom_boxplot(outlier.size=0.5) + 
  theme_bw() + 
  theme(plot.margin = unit(c(1,2.1,0,2), "cm")) + 
  theme(panel.grid.minor = element_blank(), axis.title = element_text(face="bold", size = 16),
        axis.text.x =element_text(colour="black", size = 22),
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour="black", size = 22),
        axis.text.y = element_text(colour="black", size = 22),
        strip.text.y = element_text(size = 16), legend.text=element_text(size=22), 
        legend.title=element_text(size = 22), legend.position="top",
        strip.background = element_rect(colour="black", fill="grey96")) +
  scale_x_discrete() +scale_fill_manual(values = c("brown4", "brown2", "salmon2","tan1", "wheat2"))+
  scale_y_continuous(name = expression("Zooplankton biomass(gWW/m " ^ " 2" )) 
Fig2A

fin_my<-read.csv(file="fin_month_year.csv", header=TRUE, dec= '.')
fin_my$month <- factor(fin_my$month, levels = c("Jan", "Feb", "Mar", "Apr", "May", 
                                                "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
fin_my$Year<-as.factor(fin_my$Year)
call_names <- c(X20="20-Hz note", X40="40-Hz call")

Fig2BC<-ggplot(fin_my, aes(x=month, y=rate, fill=Year)) + 
  geom_bar(width=0.5, stat="identity", position=position_dodge(), color="black") +
  theme_few() + 
  theme(plot.margin = unit(c(1,2.1,0,2), "cm")) + 
  theme(panel.grid.minor = element_blank(), axis.title = element_text(face="bold", size = 16),
        axis.text.x =element_text(colour="black", size = 22),
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour="black", size = 22),
        axis.text.y = element_text(colour="black", size = 22),
        legend.text=element_text(size=22), strip.text.y = element_text(size = 22),
        legend.title=element_text(size = 22), legend.position="top",
        strip.background = element_rect(colour="black", fill="white")) + 
  scale_x_discrete()+ facet_grid(Type~., scales="free", 
                                 labeller = as_labeller(call_names))  +
  scale_y_continuous(name = "Call rate (number of calls/hour)") + 
  scale_fill_manual(values = c("blue4", "blue", "dodgerblue1","lightskyblue", "azure1"))
Fig2BC

#STATISTICAL MODELS AND FIGURES===================================================================
#Database
fin_w_ns<-read.csv(file="weekly_fin_nosum.csv", header=TRUE, dec= '.')
summary(fin_w_ns)
fin_w_ns$year1<-as.factor(fin_w_ns$year1)
fin_w_ns$year<-as.factor(fin_w_ns$year)
fin_w_ns$X20_rate<-round(fin_w_ns$X20_rate)
fin_w_ns$X40_rate<-round(fin_w_ns$X40_rate)

#Model selection for the 20-Hz call-----------------------------------------------------
F20_1p<-glm(X20_rate~ zoo*season + year, family=poisson, data=fin_w_ns)
F20_2p<-glm(X20_rate~zoo + season + year, family=poisson, data=fin_w_ns)
F20_3p<-glm(X20_rate~zoo*season, family=poisson, data=fin_w_ns)
F20_4p<-glm(X20_rate~zoo + year, family=poisson, data=fin_w_ns)
F20_5p<-glm(X20_rate~zoo + season, family=poisson, data=fin_w_ns)
F20_6p<-glm(X20_rate~zoo, family=poisson, data=fin_w_ns)
F20_7p<-glm(X20_rate~season + year, family=poisson, data=fin_w_ns)
F20_8p<-glm(X20_rate~season, family=poisson, data=fin_w_ns)
F20_9p<-glm(X20_rate~year, family=poisson, data=fin_w_ns)

#QAIC table (Table S1)
dfun <- function(object) {
  with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
}

library(MuMIn); packageVersion("MuMIn")
x.quasipoisson <- function(...) {
  res <- quasipoisson(...)
  res$aic <- poisson(...)$aic
  res
}
F20_2p <- update(F20_1p,family="x.quasipoisson",
                 na.action=na.fail)
(gg <- dredge(F20_2p,rank="QAIC", chat=dfun(F20_1p)))


#Best selected model outputs
F20_7P<-glm(X20_rate~season_spr + year + X20_rate_lag,family=quasipoisson, data=fin_w_ns)
anova(F20_7P, test="F") #(Table 1)
summary(F20_7P) #(Table S2)
hist(resid(F20_7P))
par(mfrow=c(3,2))
plot(F20_7P)
acf(resid(F20_7P))
pacf(resid(F20_7P), ylim=c(-1,1))
Dsquared(F20_7P)

#CI and dispersion parameter (Table S2)
equation<-X20_rate~season_win + year + X20_rate_lag
F20_7P_1<-glm(equation,family=quasipoisson, data=fin_w_ns)
summary(F20_7P_1)
confint(F20_7P_1,level=0.95)
B<-1000
n<-nrow(fin_w_ns)
dispersao<-summary(F20_7P_1)$dispersion
coefs.dispersao<-rep(0,B)

muhat <- predict(F20_7P_1, type = "response")
dados<-fin_w_ns

# Parametric bootstrapping of the dispersion parameter
rqpois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}
for (b in 1:B) {
  dados$X20_rate <- rqpois(n, muhat,dispersao)
  bootGlm <- glm(equation, family = quasipoisson, data = dados)
  coefs.dispersao[b] <- summary(bootGlm)$dispersion
  print(summary(bootGlm)$dispersion)
}

# Mean of the estimates of the dispersion parameter based on the parametric bootstrap
mean.boot.dispersion<-print(mean(coefs.dispersao))

# Standard error of the estimates of the dispersion parameter based on the parametric bootstrap
se.dispersion<-print(sd(coefs.dispersao)/sqrt(n))

# Confidence interval (95%) for the dispersion parameter, based on the bootstrapped stand. error
CI.dispersion.lower<-print(dispersao-1.96*se.dispersion)
CI.dispersion.upper<-print(dispersao+1.96*se.dispersion)
Range<-CI.dispersion.upper-CI.dispersion.lower

# Confidence interval (95%) for the dispersion parameter, based on the empirical distribution of the bootstrapped dispersion parameter
quantile(coefs.dispersao, prob=c(0.025, 0.975))
Lower.limit<-unname(quantile(coefs.dispersao, prob=0.025))
Upper.limit<-unname(quantile(coefs.dispersao, prob=0.975))
Range_1<-Upper.limit - Lower.limit

#Figure 3A
s<-effect_plot(F20_7P, pred = season_spr, interval = TRUE, plot.points = TRUE,  point.color = "blue") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.title = element_text(face="bold", size = 22),
        strip.text.y = element_text(size = 22), legend.text=element_text(size=22), 
        axis.text.x =element_text(colour="black", size = 22),
        axis.text.y =element_text(colour="black", size = 22))+
  scale_y_continuous(name = "20-Hz call rate") 
s
#Figure 3B
y<-effect_plot(F20_7P, pred = year, interval = TRUE, plot.points = TRUE,  point.color = "blue") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.title = element_text(face="bold", size = 22),
        strip.text.y = element_text(size = 22), legend.text=element_text(size=22), 
        axis.text.x =element_text(colour="black", size = 22),
        axis.text.y =element_text(colour="black", size = 22), 
        axis.title.y = element_blank())+
  scale_y_continuous(name = "20-Hz call rate") 
y
grid.arrange(s,y,  ncol=2)

#Model selection for the 40-Hz call rate ----------------------------------------------------------------------
F40_1p<-glm(X40_rate~ zoo*season + year + season, family=poisson, data=fin_w_ns)
F40_2p<-glm(X40_rate~zoo + season + year, family=poisson, data=fin_w_ns)
F40_3p<-glm(X40_rate~zoo*season, family=poisson, data=fin_w_ns)
F40_4p<-glm(X40_rate~zoo + year, family=poisson, data=fin_w_ns)
F40_5p<-glm(X40_rate~zoo + season, family=poisson, data=fin_w_ns)
F40_6p<-glm(X40_rate~zoo, family=poisson, data=fin_w_ns)
F40_7p<-glm(X40_rate~season + year, family=poisson, data=fin_w_ns)
F40_8p<-glm(X40_rate~season, family=poisson, data=fin_w_ns)
F40_9p<-glm(X40_rate~year, family=poisson, data=fin_w_ns)

#Table S1
library(MuMIn); packageVersion("MuMIn")
x.quasipoisson <- function(...) {
  res <- quasipoisson(...)
  res$aic <- poisson(...)$aic
  res
}
F40_2p <- update(F40_1p,family="x.quasipoisson",
                 na.action=na.fail)
(gg <- dredge(F40_2p,rank="QAIC", chat=dfun(F40_1p)))

#Best selected model
F40_6P<-glm(X40_rate~zoo, family=quasipoisson,data=fin_w_ns)
anova(F40_6P, test="F") #Table 2
summary(F40_6P) #Table S2
hist(resid(F40_6P))
par(mfrow=c(3,2))
plot(F40_6P)
acf(resid(F40_6P))
pacf(resid(F40_6P))
Dsquared(F40_6P)

#CI and dispersion parameter (Table S2)
equation1<-X40_rate~zoo
F40_6P_1<-glm(equation1,family=quasipoisson, data=fin_w_ns)
summary(F40_6P_1)
confint(F40_6P_1,level=0.95)
B<-1000
n<-nrow(fin_w_ns)
dispersao<-summary(F40_6P_1)$dispersion
coefs.dispersao<-rep(0,B)

muhat1 <- predict(F40_6P_1, type = "response")
dados<-fin_w_ns

# Parametric bootstrapping of the dispersion parameter
# Generates a QuasiPoisson distribution based on the estimated parameters obtained from the original datset
rqpois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta-1))
}
for (b in 1:B) {
  dados$X40_rate <- rqpois(n, muhat1,dispersao)
  bootGlm <- glm(equation1, family = quasipoisson, data = dados)
  coefs.dispersao[b] <- summary(bootGlm)$dispersion
  print(summary(bootGlm)$dispersion)
}

# Mean of the estimates of the dispersion parameter based on the parametric bootstrap
mean.boot.dispersion<-print(mean(coefs.dispersao))

# Standard error of the estimates of the dispersion parameter based on the parametric bootstrap
se.dispersion<-print(sd(coefs.dispersao)/sqrt(n))

# Confidence interval (95%) for the dispersion parameter, based on the bootstrapped stand. error
CI.dispersion.lower<-print(dispersao-1.96*se.dispersion)
CI.dispersion.upper<-print(dispersao+1.96*se.dispersion)
Range<-CI.dispersion.upper-CI.dispersion.lower

# Confidence interval (95%) for the dispersion parameter, based on the empirical distribution of the bootstrapped dispersion parameter
quantile(coefs.dispersao, prob=c(0.025, 0.975))
Lower.limit<-unname(quantile(coefs.dispersao, prob=0.025))
Upper.limit<-unname(quantile(coefs.dispersao, prob=0.975))
Range_1<-Upper.limit - Lower.limit

#Figure 4
fin_w_ns<-read.csv(file="weekly_fin_nosum.csv", header=TRUE, dec= '.')
F40_6P<-glm(X40_rate~zoo, family=quasipoisson,data=fin_w_ns)
z<-effect_plot(F40_6P, pred = zoo, interval = TRUE, plot.points = TRUE, point.color = "blue") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.title = element_text(face="bold", size = 22),
        strip.text.y = element_text(size = 22), legend.text=element_text(size=22), 
        axis.text.x =element_text(colour="black", size = 22),
        axis.text.y =element_text(colour="black", size = 22))+
  scale_y_continuous(name = "40-Hz call rate index") +
  scale_x_continuous(name="Zooplankton biomass (gWW/m^2)")
z


