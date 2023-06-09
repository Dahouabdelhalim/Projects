library(glmmTMB)
library(ggplot2)
library(sjstats)
library(MuMIn)
library(dplyr)
library(ggbeeswarm)
library(rlist)

###

#Analysis of the molerat activity
#Compute the full GLMM (table 1), the minimal GLMM (supplementary table 2), figure 1, figure 2 and supplementary figure 5

###

###set working directory
setwd('C:/Users/...')
###load data
table=read.csv('DMR_activity_table.csv')
###function to cmpute dispertion of a poisson model
dispfun <- function(m) {
  r <- residuals(m,type="pearson")
  n <- df.residual(m)
  dsq <- sum(r^2)
  c(dsq=dsq,n=n,disp=dsq/n)
}

###scaling continuous variables
table$Weight.s = scale(table$Weight)
table$Weight.s[table$Sex == 'Female'] = scale(table$Weight[table$Sex == 'Female'])
table$Weight.s[table$Sex == 'Male'] = scale(table$Weight[table$Sex == 'Male'])
table$activity.log = log(table$activity+1)
table$rainfall.s = scale(table$rainfall)
table$colonySize.s = scale(table$colonySize)
table$rainfall.log=log(table$rainfall+1)
table$rainfall.log = scale(table$rainfall.log)
table$activity.jitt = abs(jitter(table$activity))

#### Compute the full GLMM ####
mod0 <-glmmTMB(activity~(sin(year1period)+cos(year1period))+breeder+Sex*breeder+day+Weight.s*Sex+rainfall.log+breeder*colonySize.s
               +(1|Animal)
               +(1|ReadingID),
               family = 'poisson',           
               data=table)              
summary(mod0)

###r squared
MuMIn::r.squaredGLMM(mod0)

###dispersion factor of the poisson model
disp = dispfun(mod0)
print(paste0('dispersion: ',disp[3]))




#### Seasonality graph (figure 2) ####
theme2 = theme(legend.text = element_blank(),panel.grid.major = element_blank(),legend.position = "none", panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white", colour = "black",size = 2, linetype = "solid"),
        axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20,vjust=-0.2),plot.title = element_text(size=5,face = "bold"), axis.text.x = element_text(size = 18,angle=90, hjust=0.5),axis.text.y = element_text(size = 18))
###Create a table with the predicted results by the full GLMM
dat = table
minrain = min(table$rainfall.log)
preddat = expand.grid(year1period=seq(0,2*pi,0.1),rainfall.log = minrain,Sex=dat$Sex[1],nbReaders=1,breeder=0,day=dat$day[1],Weight.s = 0,colonySize.s = 0,Animal = dat$Animal[1],ReadingID=dat$ReadingID[1])
preddat$activity <- predict(mod0, newdata = preddat, level = 0, reform = NA)

###Computing 95CI
mm <- model.matrix(terms(mod0),preddat)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod0)[[1]],mm))
cmult=1.96
tvar1 <- pvar1+VarCorr(mod0)[[1]]$Animal[1]  
preddat <- data.frame(
  preddat
  , plo = preddat$activity-cmult*sqrt(pvar1)
  , phi = preddat$activity+cmult*sqrt(pvar1)
  , tlo = preddat$activity-cmult*sqrt(tvar1)
  , thi = preddat$activity+cmult*sqrt(tvar1)
)

preddat = preddat %>% mutate(plo=exp(plo))
preddat = preddat %>% mutate(phi=exp(phi))
preddat = preddat %>% mutate(tlo=exp(tlo))
preddat = preddat %>% mutate(thi=exp(thi))
preddat = preddat %>% mutate(activity=exp(activity))
###plotting
figure2=ggplot(preddat, aes(x=year1period,y=activity))+  geom_point(data=table,aes(x=year1period,y=activity.jitt),shape=4,size=2,col='grey')+geom_ribbon(preddat,mapping = aes(x=year1period,ymin = plo, ymax = phi,alpha=0.5),col='black',fill='grey',alpha=0.5)+ylim(c(0,15))+
  geom_line(size=2,col='black')+xlab('Time of the year')+ylab('Activity score')+theme2+scale_x_continuous(breaks = seq(0,pi*1.9,pi/3), labels=c('Jan','Mar','May','Jul','Sep','Nov'))
print(figure2)


#### Breeder vs non-breeder model prediction (supplementary table 5) ####

theme2 = theme(panel.grid.major = element_blank(),legend.position = "none", panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white", colour = "black",size = 2, linetype = "solid"),
               axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20),plot.title = element_text(size=5,face = "bold"), axis.text.x = element_text(size = 18,angle=0, hjust=0.5,vjust=0.5),axis.text.y = element_text(size = 18))
###Create a table with the predicted results by the full GLMM
preddat = expand.grid(year1period=0.5*pi,rainfall.log = minrain,Sex=c('Female','Male'),nbReaders=1,breeder=c(0,1),day=dat$day[1],colonySize.s = 0,Weight.s = dat$Weight.s[1],Animal = 'BO1F003',ReadingID=dat$ReadingID[1])
preddat$activity <- predict(mod0, newdata = preddat, reform = NA)

###Computing 95CI
mm <- model.matrix(terms(mod0),preddat)
pvar1 <- diag(mm %*% tcrossprod(vcov(mod0)[[1]],mm))
cmult=1.96
tvar1 <- pvar1+VarCorr(mod0)[[1]]$Animal[1]  
preddat <- data.frame(
  preddat
  , plo = preddat$activity-cmult*sqrt(pvar1)
  , phi = preddat$activity+cmult*sqrt(pvar1)
  , tlo = preddat$activity-cmult*sqrt(tvar1)
  , thi = preddat$activity+cmult*sqrt(tvar1)
)

preddat = preddat %>% mutate(plo=exp(plo))
preddat = preddat %>% mutate(phi=exp(phi))
preddat = preddat %>% mutate(tlo=exp(tlo))
preddat = preddat %>% mutate(thi=exp(thi))
preddat = preddat %>% mutate(activity=exp(activity))
###plotting
suppfigure_5 = ggplot(preddat,aes(x=as.factor(breeder),y=activity,col = Sex)) +scale_x_discrete(labels=c('non-breeder','breeder'))+ geom_point(size=3,position=position_dodge(.2)) +
  ylim(c(0,6))+geom_errorbar(ymin = preddat$plo,ymax =preddat$phi,width=0.1,position=position_dodge(.2))+theme2+ylab('Activity score')+xlab('')+scale_color_manual(values=c('black','darkgrey'))
print(suppfigure_5)

#### Breeder vs non-breeder raw data + boxplot (figure 1) ####
theme2 = theme(legend.text = element_text(size=20),legend.title = element_text(size=20),legend.position = c(0.85,0.9),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white", colour = "black",size = 2, linetype = "solid"),
               axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20,vjust=-0.2),plot.title = element_text(size=5,face = "bold"), axis.text.x = element_text(size = 18,angle=0, hjust=0.5),axis.text.y = element_text(size = 18))

figure1=ggplot(table,mapping=aes(x=as.factor(breeder),y=activity.jitt,col = Sex)) +scale_x_discrete(labels=c('non-breeder','breeder'))+
    geom_boxplot(position=position_dodge(.9),outlier.shape = NA)+geom_quasirandom(size=2,shape=4,dodge.width=.9)+
    ylim(c(0,15))+theme2+
    ylab('Activity score')+xlab('')+scale_color_manual(values=c('black','darkgrey'))
print(figure1)

#### minimal GLMM model (supplementary table 2) ####
minimal_model <-glmmTMB(activity~Sex+breeder+colonySize.s
                        +(1|Animal)
                        +(1|ReadingID),
                        family = 'poisson',           
                        data=table)              
summary(minimal_model)

###r squared
MuMIn::r.squaredGLMM(minimal_model)

###dispersion factor of the poisson model
disp = dispfun(minimal_model)
print(paste0('dispersion: ',disp[3]))

#### printing figures to pdf ####
pdf('figure1.pdf')
figure1
dev.off()
pdf('figure2.pdf')
figure2
dev.off()
pdf('supplementaryFigure7.pdf')
suppfigure_5
dev.off()