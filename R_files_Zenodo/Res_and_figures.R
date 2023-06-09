library(tidyverse)
library(viridis)
library(grid)
library(rje)
library(colorRamps)
library(mgcv)
library(visreg)
library(cowplot)  
library(GGally)

load("ssdat_28062020.Rdata")

##### compute deviations #####
# population abundance
ssdat$deviationAb <- as.vector(scale(ssdat$TREND, center=F) - scale(ssdat$meanAB)) 
# species richness
ssdat$deviation <- as.vector(scale(ssdat$slp_ric, center=F) - scale(ssdat$meanRIC)) 

##### descriptive deviations ######
ssdat$groups <- "T+A+"
ssdat$groups[scale(ssdat$TREND, center=F) > 0 & scale(ssdat$meanAB) <0 ] <- "T+A-"
ssdat$groups[scale(ssdat$TREND, center=F) < 0 & scale(ssdat$meanAB) <0 ] <- "T-A-"
ssdat$groups[scale(ssdat$TREND, center=F) > 0 & scale(ssdat$meanAB) >0 ] <- "T-A+"

ggplot(ssdat, aes(x=scale(meanAB), y=scale(TREND, center=F), fill=groups))+
  geom_point(shape=21,aes(size=1/SD), alpha=0.75)+
  # geom_abline(slope=1, intercept=0, size = 1, col= "red")+
  scale_fill_viridis_d("space-time \\n deviation", guide=F)+
  theme_classic()+
  guides(size="none")+
  stat_smooth(method='lm', se=T, linetype=1, col="black", fill = "white", alpha = 0.75)+
  xlab("relative population size")+
  ylab("relative population trend")+
  coord_fixed(ratio=1)+
  theme(legend.position = "bottom")

scatter<-
  ggplot(ssdat, aes(x=scale(meanAB), y=scale(TREND, center=F), fill=deviationAb))+
  geom_point(shape=21,aes(size=1/(SD*2)), alpha=0.75)+
  geom_abline(slope=1, intercept=0, size = 1, col= "red", linetype = 2)+
  scale_fill_viridis("space-time \\n deviation", guide=F)+
  theme_classic()+
  scale_size_continuous(range = c(0.1,3))+
  guides(size="none")+
  stat_smooth(method='lm', se=T, linetype=1, col="black", fill = "white", alpha = 0.75)+
  xlab("relative population size")+
  ylab("relative population trend")+
  coord_fixed(ratio=1)+
  theme(legend.position = "bottom")

distribution<-
  ggplot(ssdat, aes(x=deviationAb, fill=..x..))+
  geom_histogram()+
  geom_vline(xintercept = 0, col = "red", linetype = 2, size = 1)+
  theme_classic()+
  scale_fill_viridis(guide = F)+
  xlab("space-time difference")



png(file = 'fig2.png',width = 2000, height=1000, res =250)
plot_grid(scatter, distribution, axis = "b", labels = "auto")
dev.off()  
pdf(file = 'fig2.pdf',width = 16, height=8)
plot_grid(scatter, distribution, axis = "b", labels = "auto")
dev.off() 


map<-
  ggplot(ssdat, aes(x=X, y=Y, fill=deviationAb))+
  geom_point(shape=21, aes(size=1/SD), alpha=0.8)+
  scale_fill_viridis(guide = F)+
  theme_classic()+
  scale_size(range = c(0.1, 3))+
  xlab("Longitude")+
  ylab("Latitude")+
  coord_fixed(ratio=1)+
  guides(size="none", fill = "none")
  
mod_map <- gamm(deviation ~ s(X, Y), random = list(BIOGEO = ~1), data = na.omit(ssdat))
plot(mod_map$gam)
summary(mod_map$gam)

Longitude  <-
  ggplot(data = ssdat, aes(x = X, y = deviationAb, fill = deviation))+
  geom_point(shape=21,aes(size=1/SD), alpha=0.8)+
  geom_smooth(method = "lm", col="black", fill = "white", alpha = 0.8)+
  # scale_fill_cmocean(name = "tarn", limits = c(-4.1, 4.1))+
  scale_fill_viridis(guide = F)+
  theme_classic()+
  xlab("Longitude")+
  ylab("space-time deviation")+
  scale_size(guide = 'none', range =c(0.1,3))+
  guides(size="none", fill = "none")


mod_longitude <- lme(deviation ~ scale(X), random = ~1|BIOGEO, data = na.omit(ssdat))
summary(mod_longitude)

Latitude <-
  ggplot(data = ssdat, aes(x = Y, y = deviationAb, fill = deviation))+
  geom_point(shape=21,aes(size=1/SD), alpha=0.8)+
  geom_smooth(method = "lm", col="black", fill = "white", alpha = 0.8)+
  # scale_fill_cmocean(name = "tarn", limits = c(-4.1, 4.1))+
  scale_fill_viridis(guide = F)+
  theme_classic()+
  xlab("Latitude")+
    ylab("")+
    guides(size="none", fill = "none")+
  scale_size(guide = 'none', range =c(0.1,3))

mod_latitude <- lme(deviation ~ scale(Y), random = ~1|BIOGEO, data = na.omit(ssdat))
visreg(mod_latitude)
summary(mod_latitude)

png(file = 'fig3.png',width = 2400, height=800, res = 200)
plot_grid(map, Latitude, Longitude, nrow = 1, axis = "b", labels = "auto")
dev.off()  
pdf(file = 'fig3.pdf',width = 24, height=8)
plot_grid(map, Latitude, Longitude, nrow = 1, axis = "b", labels = "auto")
dev.off()  


## check the sites with significant temporal trends
ssdat$CHECK<-"YES" 
ssdat$CHECK[abs(ssdat$TREND)- (1.96*ssdat$SD)<0]<-"NO"


(nrow(ssdat[ssdat$CHECK=="YES",]) / nrow(ssdat))*100
ssdat$CHECK2<-"YES" 
ssdat$CHECK2[abs(ssdat$TREND) < 0.001]<-"NO"


ggplot(ssdat, aes(x=scale(meanAB), y=scale(TREND, center=F), fill=CHECK))+
  geom_point(shape=21,aes(size=1/SD), alpha=0.4)+
  geom_abline(slope=1, intercept=0)+
  # scale_fill_gradientn(colours=cubeHelix(30))+
  theme_classic()+
  guides(size="none")+
  xlab("relative abundance")+
  ylab("relative temporal trend")


#### correlation between predictors####
# Correlation plot
ggcorr(ssdat[, c("HII", "Natural",
          "colExtRatio", "completeness"
)], palette = "RdBu", label = TRUE, method = c("pairwise", "kendall"), digits =3)

png("SOURCE SINK/figures/final figures/Fig.S2_predictors_correlations.png", width = 1600, height= 1600, res = 150)
ggpairs(ssdat[, c("HII", "Natural",
                  "colExtRatio", "completeness")], title = "", axisLabels = "show")
dev.off()

#### correlation between deviations and predictors####

### colonization extinction ####
colExt<-
  ggplot(ssdat[ssdat$colExtRatio >0.9,], aes(y=deviationAb, x=colExtRatio, weight=1/SD))+
  geom_point(shape=21,aes(size=1/SD, fill = deviation), alpha=0.8)+
  geom_smooth(method = "lm", col="black", fill = "lightgrey", alpha = 0.8)+
  theme(legend.position="none")+
  labs(x="colonisation extinction ratio", y="space-time deviation")+
  theme_classic()+
  scale_fill_viridis(guide = F)+
  # facet_grid(~groups)+
  scale_size(guide = 'none', range =c(0.1,3))
  # annotate("text", x=0.9, y=4, label="(b)", size=6)

  
mod_colExt <- gam(DebtCredit ~ colExtRatio +  s(X,Y), data=na.omit(ssdat), na.action="na.fail")
summary(mod_colExt)
# colExt<-
visreg(mod_colExt, 'colExtRatio', gg=T)+
  geom_point(  shape=21, alpha=0.8)+
  stat_smooth(method='lm', se=F)+
  theme(legend.position="none")+
  labs(x="colonisation extinction ratio", y="Disequilibrium Value")+
  theme_classic()+
  scale_size(guide = 'none')
  # annotate("text", x=0.75, y=4, label="(b)", size=6)
dr <- dredge(mod)

#### completeness of communities #####
completeness <- 
ggplot(ssdat, aes(y=deviationAb, x=completeness, weight=1/SD))+
  geom_point(shape=21,aes(size=1/SD, fill = deviation), alpha=0.8)+
  geom_smooth(method = "lm", col="black", fill = "lightgrey", alpha = 0.8)+
  theme(legend.position="none")+
  labs(x="community completeness", y="")+
  theme_classic()+
  scale_fill_viridis(guide = F)+
  # facet_grid(~groups)+
  scale_size(guide = 'none', range =c(0.1,3))


# mod<-lm(DebtCredit ~ 1 + completeness, weight=1/SD, data=ssdat, na.action=na.fail)
mod_dark<-gam(deviationAb ~ completeness + s(X,Y), data=na.omit(ssdat), na.action='na.fail')
summary(mod_dark)

dark<-
visreg(mod_dark, "completeness", gg=T)+
geom_point()+
stat_smooth(method='lm', se=F)+
#   scale_y_continuous(trans="log1p", breaks=c(0,10,10,50,100, 200))+
#   scale_x_continuous(trans="log1p")+
theme(legend.position="none")+
labs(x="Completeness Index", y="Disequilibrium Value")+
theme_classic()+
scale_size(guide = 'none')+
annotate("text", x=1.3, y=4, label="(a)", size=6)


### human influence ####
hii<-
  ggplot(ssdat, aes(y=deviationAb, x=HII, weight=1/SD))+
  geom_point(shape=21,aes(size=1/SD, fill = deviation), alpha=0.8)+
  geom_smooth(method = "lm", col="black", fill = "lightgrey", alpha = 0.8)+
    theme(legend.position="none")+
    labs(x="human influence index", y="space-time deviation")+
    theme_classic()+
    scale_fill_viridis(guide = F)+
    # facet_grid(~groups)+
  scale_size(guide = 'none', range =c(0.1,3))
  # annotate("text", x=0.9, y=4, label="(b)", size=6)

mod<-gam(deviationAb ~ scale(HII) +  s(X,Y), weight=1/SD, data=ssdat, na.action=na.fail)
summary(mod)
visreg(mod)


natural <- 
  ggplot(ssdat, aes(y=deviationAb, x=Natural, weight=1/SD))+
  geom_point(shape=21,aes(size=1/SD, fill = deviation), alpha=0.8)+
  geom_smooth(method = "lm", col="black", fill ="lightgrey", alpha = 0.8)+
  theme(legend.position="none")+
  labs(x="wilderness", y="")+
  theme_classic()+
  scale_fill_viridis(guide = F)+
  # facet_grid(~groups)+
  scale_size(guide = 'none', range =c(0.1,3))
# annotate("text", x=0.9, y=

mod<-gam(deviationAb ~ 1 + scale(Natural) +  s(X,Y),weight=1/SD,data=ssdat, na.action=na.fail)
summary(mod)

ggplot(ssdat, aes(y=deviationAb, x=Natural, weight=1/SD))+
  geom_point(shape=21,aes(size=1/SD, fill = TREND), alpha=0.8)+
  geom_smooth(method = "lm", col="black", fill ="lightgrey", alpha = 0.8)+
  theme(legend.position="none")+
  labs(x="wilderness", y="")+
  theme_classic()+
  scale_fill_viridis(guide = F)+
  # facet_grid(~groups)+
  scale_size(guide = 'none')
# annotate("text", x=0.9, y=


ggplot(ssdat, aes(y=TREND, x=Agricultural, weight=1/SD))+
  geom_point(shape=21,aes(size=1/SD, fill = deviation), alpha=0.8)+
  geom_smooth(method = "lm", col="black", fill ="lightgrey", alpha = 0.8)+
  theme(legend.position="none")+
  labs(x="wilderness", y="")+
  theme_classic()+
  scale_fill_viridis(guide = F)+
  # facet_grid(~groups)+
  scale_size(guide = 'none')
# annotate("text", x=0.9, y=

mod<-gamm(deviationAb ~ 1 + scale(Natural) +  s(X,Y),weight=1/SD, random=list(BIOGEO=~1), data=na.omit(ssdat), na.action=na.pass)
summary(mod$gam)

# png(file = 'fig4.png',width = 1600, height= 1600, res =200)
plot_grid(hii, natural, colExt, completeness,  nrow  = 2, axis = "b, r", labels = "auto")
# dev.off()  
# pdf(file = 'fig4.pdf',width = 16, height=16)
plot_grid(hii, natural, colExt, completeness, nrow  = 2, axis = "b, r", labels = "auto")
# dev.off() 


####### figure 4 ########
load("SOURCE SINK/ssdat_10062016.Rdata")
mod_ssdat <- ssdat[-c(919, 1139, 1775), c("DebtCredit", "SD","completeness", "turnOver", "colExtRatio", "Natural", "HII", "BIOGEO", "X", "Y", "alpha")]
mod_ssdat$DebtCredit <- as.vector(mod_ssdat$DebtCredit)

## completeness ------------------------------------------------------------------------------------
mod_dark<-gamm(deviationAb ~ 1 + completeness +alpha + s(X,Y), random=list(BIOGEO=~1), data=mod_ssdat, na.action='na.fail')
mod_dark$gam$data <- mod_ssdat
summary(mod_dark$gam)

dark<-
  visreg(mod_dark$gam,"completeness",points=list(size=(1/mod_ssdat$SD)/500, pch=21, col="black", alpha=0.2), gg=T)+
  theme(legend.position="none")+
  labs(x="Completeness Index", y="Disequilibrium Value")+
  theme_classic()+
  scale_size(guide = 'none')+
  annotate("text", x=1.3, y=4, label="(a)", size=6)

## colonization extinction ------------------------------------------------------------------------------------
mod_colExt <- gamm(deviationAb ~ 1   + s(X,Y) +  colExtRatio, data=mod_ssdat,random=list(BIOGEO=~1), na.action="na.fail")
summary(mod_colExt$gam)
mod_colExt$gam$data <- mod_ssdat

colExt<-
  visreg(mod_colExt$gam, 'colExtRatio', points=list(size=(1/mod_ssdat$SD)/500, pch=21, col="black", alpha=0.2), gg=T)+ 
  theme(legend.position="none")+
  labs(x="Colonisation-Extinction ratio", y="Disequilibrium Value")+
  theme_classic()+
  scale_size(guide = 'none')+
  annotate("text", x=0.75, y=4, label="(b)", size=6)



## turnOver ----------------------------------------------------------------------------------------------
mod_turnOver<-gamm(deviationAb ~ 1 +  scale(turnOver) + s(X,Y) , data=mod_ssdat, random=list(BIOGEO=~1), na.action=na.fail)
summary(mod_turnOver$gam)
mod_turnOver$gam$data <- mod_ssdat

turnOver<-
  visreg(mod_turnOver$gam, "turnOver", points=list(size=(1/mod_ssdat$SD)/500, pch=21, col="black", alpha=0.2),gg=T)+
  theme(legend.position="none")+
  labs(x="Turn-over", y="Disequilibrium Value")+
  theme_classic()+
  scale_size(guide = 'none')+
  annotate("text", x=1.35, y=4, label="(c)", size=6)


## Natural areas ------------------------------------------------------------------------------------
mod_nat<-gamm(deviationAb~ 1 + Natural + s(X,Y), data=mod_ssdat, random=list(BIOGEO=~1),na.action=na.fail)
summary(mod_nat$gam)
mod_nat$gam$data = mod_ssdat

Natural <- 
  visreg(mod_nat$gam, "Natural", points=list(size=(1/mod_ssdat$SD)/500, pch=21, col="black", alpha=0.2), gg=T)+ 
  theme(legend.position="none")+
  labs(x="% Natural areas", y="Debt-Credit value")+
  theme_classic()+
  scale_size(guide = 'none')+
  annotate("text", x=0, y=4, label="(d)", size=6)


## Human Influence Index ------------------------------------------------------------------------------------
mod_hii<-gamm(deviationAb ~ 1+ HII +  s(X,Y), random=list(BIOGEO=~1), data=mod_ssdat, na.action=na.fail)
summary(mod_hii$gam)
mod_hii$gam$data <- mod_ssdat

hii<-
  visreg(mod_hii$gam, "HII", points=list(size=(1/mod_ssdat$SD)/500, pch=21, col="black", alpha=0.2), gg=T)+
  # geom_point()+
  # stat_smooth(method='lm', se=F)+
  #   scale_y_continuous(trans="log1p", breaks=c(0,10,10,50,100, 200))+
  #   scale_x_continuous(trans="log1p")+
  theme(legend.position="none")+
  labs(x="Human Influence Index", y="Disequilibrium Value")+
  theme_classic()+
  scale_size(guide = 'none')+
  annotate("text", x=15, y=4, label="(e)", size=6)+
  
  
  # pdf("fig4.pdf",width=4, height=6)
colExt
turnOver
dark 
Natural
hii
# dev.off()


library(MuMIn)
AIC_selection <- dredge(mod_all)
AIC_selection_rounded <- round(AIC_selection[1:5,-7],3)
write.table(AIC_selection_rounded, file="SOURCE SINK/figures/table_multivariate_AIC_review.csv", row.names = F, sep=";")
average_model <- model.avg(AIC_selection, subset = delta < 4, rank=AIC)





