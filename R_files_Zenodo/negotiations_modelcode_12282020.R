library(lme4)
library(lmerTest)
library(ggplot2)
library(gridExtra)
library(blmeco) 

# Import spreadsheet for egg counts per clutch
clutch_laid<-read.csv('clutch_laid.csv')


# Import spreadsheet for egg counts per month
month_laid<-read.csv('month_laid.csv')


# Import spreadsheet for embryo survival and embryos maturing per clutch
embryo_surv<-read.csv('embryo_surv.csv')


# Import spreadsheet for embryos maturing per month
embryo_mat<-read.csv('embryo_mat.csv')


# Import spreadsheet for parental care
Pcare<-read.csv('Pcare.csv')



## model for number of eggs laid per clutch with scaled FSL as covariate
mod.eggs1<-glmer(Laid~Treatment*as.factor(Month)+scale(FSL)+(1|ID),data=na.omit(clutch_laid),family="poisson")
dispersion_glmer(mod.eggs1) # model is overdispersed
# correct for overdispersion using an observation level random effect (OLRE)
mod.eggs2<-glmer(Laid~Treatment*as.factor(Month)+scale(FSL)+(1|ID)+(1|OLRE),data=na.omit(clutch_laid),family="poisson")
dispersion_glmer(mod.eggs2) # model is not overdispersed
summary(mod.eggs2)

# generate predicted values from the model for plots
clutch_laid<-na.omit(clutch_laid)
clutch_laid$predictlaid<-predict(mod.eggs2,type='response')



# model for number of eggs laid per month with scaled FSL as a covariate
mod.eggs3<-glmer(Laid~Treatment*as.factor(Month)+scale(FSL)+(1|ID),data=na.omit(month_laid),family=poisson)
dispersion_glmer(mod.eggs3) # model is overdispersed
# correct for overdispersion using an observation level random effect (OLRE)
mod.eggs4<-glmer(Laid~Treatment*as.factor(Month)+scale(FSL)+(1|OLRE)+(1|ID),data=na.omit(month_laid),family=poisson)
dispersion_glmer(mod.eggs4)
summary(mod.eggs4)
month_laid<-na.omit(month_laid)
month_laid$predictlaid<-predict(mod.eggs4,type='response')
plot(mod.eggs4)


# Models for female and male time tending with and without eggs laid as covariate
# fit model with eggs laid
mod.Mtend<-lmer(M.tend~Month*Treatment+Age+Laid+(1|ID),data=na.omit(Pcare))
mod.Ftend<-lmer(F.tend~Month*Treatment+Age+Laid+(1|ID),data=na.omit(Pcare))

summary(mod.Mtend)

summary(mod.Ftend)

Pcare<-na.omit(Pcare)


# predicted values from the models
Pcare$predictMtend<-predict(mod.Mtend)
Pcare$predictFtend<-predict(mod.Ftend)

## diagnostic plots male time tending
plot(mod.Mtend)
# examine normality of residuals
qqnorm(resid(mod.Mtend))
qqline(resid(mod.Mtend))
# examine homogeneity of variances
plot(mod.Mtend, resid(.,type="pearson")~fitted(.)|Treatment, abline=0,pch=20)
plot(mod.Mtend, resid(.,type="pearson")~fitted(.)|ID, abline=0,pch=20)

## diagnostic plots female time tending
plot(mod.Ftend)
# examine normality of residuals
qqnorm(resid(mod.Ftend))
qqline(resid(mod.Ftend))
# examine homogeneity of variances
plot(mod.Ftend, resid(.,type="pearson")~fitted(.)|Treatment, abline=0,pch=20)
plot(mod.Ftend, resid(.,type="pearson")~fitted(.)|ID, abline=0,pch=20)



## proportion of eggs laid that survived to hatching (embryo survival) with male time tending as covariate
mod_surv<-glmer(Prophatch~Treatment*as.factor(Month)+M.tend+(1|ID),data=na.omit(embryo_surv),family="binomial")
summary(mod_surv)
qqnorm(resid(mod_surv))
qqline(resid(mod_surv))
plot(mod_surv, resid(.,type="pearson")~fitted(.)|Treatment, abline=0,pch=20)
plot(mod_surv, resid(.,type="pearson")~fitted(.)|ID, abline=0,pch=20)


# predicted values from the model
embryo_surv<-na.omit(embryo_surv)
embryo_surv$predictprophatch<-predict(mod_surv,type='response')


## Embryos maturing per clutch with and without male tending and eggs laid
mod_mat1<-glmer(Hatched~Treatment*as.factor(Month)+scale(Laid)+scale(M.tend)+(1|ID)+(1|OLRE),data=na.omit(embryo_surv),family=poisson)
summary(mod_mat1)
dispersion_glmer(mod_mat1) # model is overdispersed, correct with observation level random effect

mod_mat2<-glmer(Hatched~Treatment*as.factor(Month)+(1|ID)+(1|OLRE),data=na.omit(embryo_surv),family=poisson)
dispersion_glmer(mod_mat2)
summary(mod_mat2)
# predicted values from the model
embryo_surv$predicthatch<-predict(mod_mat2,type='response')

## Number of embryos hatched per month with and without male tending and eggs laid
mod_mat3<-glmer(Hatched~Treatment*as.factor(Month)+scale(Laid)+scale(M.tend)+(1|ID),data=na.omit(embryo_mat),family=poisson)
dispersion_glmer(mod_mat3) # model is overdispersed, correct with observation level random effect

mod_mat4<-glmer(Hatched~Treatment*as.factor(Month)+scale(Laid)+scale(M.tend)+(1|ID)+(1|OLRE),data=na.omit(embryo_mat),family=poisson)
dispersion_glmer(mod_mat4)
summary(mod_mat4)

# predicted values from the model
embryo_mat<-na.omit(embryo_mat)
embryo_mat$predicthatch<-predict(mod_mat4,type='response')

## Figure 1 - Female Effort
laid2<-clutch_laid
laid2$Treatment <- factor(laid2$Treatment, levels = c("C", "F","M"),
                          labels = c("Control", "FH","MH"))
laid3<-month_laid
laid3$Treatment <- factor(laid3$Treatment, levels = c("C", "F","M"),
                          labels = c("Control", "FH","MH"))
Pcare$Treatment <- factor(Pcare$Treatment, levels = c("C", "F","M"),
                          labels = c("Control", "FH","MH"))
embryo_surv$Treatment <- factor(embryo_surv$Treatment, levels = c("C", "F","M"),
                                    labels = c("Control", "FH","MH"))

hatch2<-embryo_mat
hatch2$Treatment <- factor(hatch2$Treatment, levels = c("C", "F","M"),
                           labels = c("Control", "FH","MH"))


(f1<-ggplot(laid2, aes(x = as.factor(Month), y = Laid, lty=Treatment)) +
    #geom_point(aes(y = Laid,shape= Treatment),size=2,alpha=0.2,position=position_dodge(width=0.5))+
    #geom_point(aes(y = Laid,shape= Treatment),size=2,position=position_jitter(h=0,w=0.2),alpha=0.2) +
    geom_boxplot(fill=NA,colour=alpha("black",0.4),width=0.8,size=0.2)+
    #geom_pointrange(aes(ymin=predictEggslower,ymax=predictEggsupper,shape=Treatment,position='dodge'),alpha=0.2)+
    geom_smooth(size = 1,method="lm",se=F, col="black",aes(group=Treatment,Month,predictlaid)) +
    #geom_line(aes(y=meaneggs,group=Treatment),size=.75,position=position_dodge(width=0.2))+
    facet_grid(.~Treatment)+
    labs(x = "Manipulation",
         y = "Number of Eggs Laid per Clutch",tag = "(a)")+ scale_color_discrete(name = "Manipulation")+ 
    scale_x_discrete(labels = c('pre','post'))+theme_classic())


(f2<-ggplot(laid3, aes(x = as.factor(Month), y = Laid, lty=Treatment)) +
    #geom_point(aes(y = Laid,shape= Treatment),size=2,alpha=0.2,position=position_dodge(width=0.5))+
    #geom_point(aes(y = Laid,shape= Treatment),size=2,position=position_jitter(h=0,w=0.2),alpha=0.2) +
    geom_boxplot(fill=NA,colour=alpha("black",0.4),width=0.8,size=0.2)+
    #geom_pointrange(aes(ymin=predictEggslower,ymax=predictEggsupper,shape=Treatment,position='dodge'),alpha=0.2)+
    geom_smooth(size = 1,method="lm",se=F, col="black",aes(group=Treatment,Month,predictlaid)) +
    #geom_line(aes(y=meaneggs,group=Treatment),size=.75,position=position_dodge(width=0.2))+
    facet_grid(.~Treatment)+
    labs(x = "Manipulation",
         y = "Number of Eggs Laid per Month",tag = "(b)")+ scale_color_discrete(name = "Manipulation")+ 
    scale_x_discrete(labels = c('pre','post'))+theme_classic())

(f3<-ggplot(Pcare, aes(x = as.factor(Month), y = F.tend, lty=Treatment)) +
    #geom_point(aes(y = Laid,shape= Treatment),size=2,alpha=0.2,position=position_dodge(width=0.5))+
    #geom_point(aes(y = Laid,shape= Treatment),size=2,position=position_jitter(h=0,w=0.2),alpha=0.2) +
    geom_boxplot(fill=NA,colour=alpha("black",0.4),width=0.8,size=0.2)+
    #geom_pointrange(aes(ymin=predictEggslower,ymax=predictEggsupper,shape=Treatment,position='dodge'),alpha=0.2)+
    geom_smooth(size = 1,method="lm",se=F, col="black",aes(group=Treatment,Month,predictFtend)) +
    #geom_line(aes(y=meaneggs,group=Treatment),size=.75,position=position_dodge(width=0.2))+
    facet_grid(.~Treatment)+
    labs(x = "Manipulation",
         y = "Female Time Spent Tending (mins)",tag = "(c)")+ scale_color_discrete(name = "Manipulation")+ 
    scale_x_discrete(labels = c('pre','post'))+theme_classic())

grid.arrange(f1,f2,f3,ncol=1)



(p1<-ggplot(Pcare, aes(x = as.factor(Month), y = M.tend, lty=Treatment)) +
    #geom_point(aes(y = Laid,shape= Treatment),size=2,alpha=0.2,position=position_dodge(width=0.5))+
    #geom_point(aes(y = Laid,shape= Treatment),size=2,position=position_jitter(h=0,w=0.2),alpha=0.2) +
    geom_boxplot(fill=NA,colour=alpha("black",0.4),width=0.8,size=0.2)+
    #geom_pointrange(aes(ymin=predictEggslower,ymax=predictEggsupper,shape=Treatment,position='dodge'),alpha=0.2)+
    geom_smooth(size = 1,method="lm",se=F, col="black",aes(group=Treatment,Month,predictMtend)) +
    #geom_line(aes(y=meaneggs,group=Treatment),size=.75,position=position_dodge(width=0.2))+
    facet_grid(.~Treatment)+
    labs(x = "Manipulation",
         y = "Male Time Spent Tending (mins)",tag = "(a)")+ scale_color_discrete(name = "Manipulation")+ 
    scale_x_discrete(labels = c('pre','post'))+theme_classic())

(p2<-ggplot(embryo_surv, aes(x = as.factor(Month), y = Prophatch, lty=Treatment)) +
    #geom_point(aes(y = Laid,shape= Treatment),size=2,alpha=0.2,position=position_dodge(width=0.5))+
    #geom_point(aes(y = Laid,shape= Treatment),size=2,position=position_jitter(h=0,w=0.2),alpha=0.2) +
    geom_boxplot(fill=NA,colour=alpha("black",0.4),width=0.8,size=0.2)+
    #geom_pointrange(aes(ymin=predictEggslower,ymax=predictEggsupper,shape=Treatment,position='dodge'),alpha=0.2)+
    geom_smooth(size = 1,method="lm",se=F, col="black",aes(group=Treatment,Month,predictprophatch)) +
    #geom_line(aes(y=meaneggs,group=Treatment),size=.75,position=position_dodge(width=0.2))+
    facet_grid(.~Treatment)+
    labs(x = "Manipulation",
         y = "Proportion of Embryos Surviving per Clutch",tag = "(b)")+ scale_color_discrete(name = "Manipulation")+ 
    scale_x_discrete(labels = c('pre','post'))+theme_classic())

grid.arrange(p1,p2,nrow=2)

embryo_mat$Treatment <- factor(embryo_mat$Treatment, levels = c("C", "F","M"),
                                      labels = c("Control", "FH","MH"))
(g1<-ggplot(embryo_surv, aes(x = as.factor(Month), y = Hatched, lty=Treatment)) +
    #geom_point(aes(y = Laid,shape= Treatment),size=2,alpha=0.2,position=position_dodge(width=0.5))+
    #geom_point(aes(y = Laid,shape= Treatment),size=2,position=position_jitter(h=0,w=0.2),alpha=0.2) +
    geom_boxplot(fill=NA,colour=alpha("black",0.4),width=0.8,size=0.2)+
    #geom_pointrange(aes(ymin=predictEggslower,ymax=predictEggsupper,shape=Treatment,position='dodge'),alpha=0.2)+
    geom_smooth(size = 1,method="lm",se=F, col="black",aes(group=Treatment,Month,predicthatch)) +
    #geom_line(aes(y=meaneggs,group=Treatment),size=.75,position=position_dodge(width=0.2))+
    facet_grid(.~Treatment)+
    labs(x = "Manipulation",
         y = "Number of Embryos Maturing per Clutch",tag = "(a)")+ scale_color_discrete(name = "Manipulation")+ 
    scale_x_discrete(labels = c('pre','post'))+theme_classic())

(g2<-ggplot(embryo_mat, aes(x = as.factor(Month), y = Hatched, lty=Treatment)) +
    #geom_point(aes(y = Laid,shape= Treatment),size=2,alpha=0.2,position=position_dodge(width=0.5))+
    #geom_point(aes(y = Laid,shape= Treatment),size=2,position=position_jitter(h=0,w=0.2),alpha=0.2) +
    geom_boxplot(fill=NA,colour=alpha("black",0.4),width=0.8,size=0.2)+
    #geom_pointrange(aes(ymin=predictEggslower,ymax=predictEggsupper,shape=Treatment,position='dodge'),alpha=0.2)+
    geom_smooth(size = 1,method="lm",se=F, col="black",aes(group=Treatment,Month,predicthatch)) +
    #geom_line(aes(y=meaneggs,group=Treatment),size=.75,position=position_dodge(width=0.2))+
    facet_grid(.~Treatment)+
    labs(x = "Manipulation",
         y = "Number of Embryos Maturing per Month",tag = "(b)")+ scale_color_discrete(name = "Manipulation")+ 
    scale_x_discrete(labels = c('pre','post'))+theme_classic())
grid.arrange(g1,g2,nrow=2)

