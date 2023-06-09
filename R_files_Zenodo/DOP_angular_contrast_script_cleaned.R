library(nlme)
library(lme4)
library(car)
library(effects)
library(lsmeans)
library(pbkrtest)
library(effects)
#install.packages("ggplot2")
library("r2glmm")
library(ggplot2)
library(MuMIn)
library(semEff)
library(RRPP)
#source("PRESS-stats/model_fit_stats.R")
#source("PRESS-stats/pred_r_squared.R")
#source("PRESS-stats/PRESS.R")
citation("RRPP")

###A note on R2 values in mixed effects models, synthesizing what I have learned. There are various R2 values
#for linear models, which address various deficiencies in the R2. As I understand it, the regular R2 is 
#insufficient for complex models, as it increases with increasing model fit without penalizing for overfitting,
#also for mixed effect models, the introduction of the Random effects introduces new sources of variance which 
#are not accounted for by regular or pseudo R2. the marginal R2 addresses this by estimating the proportion of 
#variance explained by the fixed effects alone, and the conditional R2 includes both the variance explained by 
#fixed and random effects.  However, none of these address the overfitting issue yet. This is where the predicted
#R2 is handy: as it uses cross-validation (subsetting the data and re-fitting the model) to determine how well 
#the model fits given different data (in this case a subset of your data), which should drastically penalize 
#models that are fit to the shape of the data rather than shape of the actual effect you are modeling.

#in this script I use two packages to calculate and report R2 values, the r.squaredGLMM function of 'MuMIn' for
#conditional and marginal R2, and the R2 function of the semEFF package for the predicted and adjusted R2

# setwd("C:/Users/richa/Dropbox/Harvard/Research/Gary polarization manuscript work/NY_CCT_specimen_model_polDATA")
# setwd("C:/Users/Astaroph/Dropbox/Harvard/Research/Gary polarization manuscript work/NY_CCT_specimen_model_polDATA")
dop=read.csv('NY_CCT_DOP_data_R_ready.csv')

##This dataset contains 2 new females that Cheng-Chia provided, but he notes that
##they are much less bright compared  to the other two, and indeed, they are outliers
##so I believe I need to exclude them...
#dop=read.csv('NY_CCT_DOP_data_R_ready_9_10_20_updated.csv')

#Ind=subset(dop,Individual=="F5")

I20=subset(dop,I_angle==20)
I30=subset(dop,I_angle==30)
I40=subset(dop,I_angle==40)
I50=subset(dop,I_angle==50)
I60=subset(dop,I_angle==60)
I70=subset(dop,I_angle==70)
mean_se <- function(x, mult = 1) {  
  x <- na.omit(x)
  se <- mult * sqrt(var(x) / length(x))
  mean <- mean(x)
  data.frame(y = mean, ymin = mean - se, ymax = mean + se)
}

library(ggplot2)
#Figures 4E (UV) and 4F (Blue)
violin_1=ggplot(
  data = dop
  , mapping = aes(
    y =DOP_UV
    , x=(Angular_contrast)
    , colour=Sex
  )
)+
  labs(
    x = 'Angular contrast'
    ,title= 'DOP-Blue by overall angular contrast in viewing and incident angle'
    , y = 'DOP-Blue'
  )+
  # geom_point(alpha=0.3,
  #             mapping = aes(x=(Angular_contrast),
  #                           colour =Sex
  #             )
  #             , size=1.8
  # )+
#  geom_violin(alpha=.4, draw_quantiles = c(0.25, 0.5, 0.75),width=1,
#              mapping=aes(fill=factor(Sex),colour=factor(Sex)))+
  stat_summary(fun.data=mean_se, geom="ribbon", alpha=0.15,mapping=aes(fill=Sex))+
  stat_summary(fun.data=mean_se, geom="line", alpha=0.95,size=1,mapping=aes(color=Sex))+
  # stat_summary(fun.data=mean_se, geom="line", alpha=0.35,size=1,linetype='dotted',mapping=aes(group=Individual))+
    # geom_smooth(method="loess",alpha=.4, se =TRUE)+
  facet_grid(.~.)
violin_1
violin_1_2=violin_1+scale_color_manual(values=c("#fdae61","#2c7bb6"))+
  scale_fill_manual(values=c("#fdae61","#2c7bb6"))+theme_bw()
violin_1_2

violin_1_3=violin_1_2

#violin_1_3=violin_1_2+scale_y_log10()
violin_1_3
violin_1_4=violin_1_3
violin_1_4
#Resizing and moving the text in the graph
violin_1_5=violin_1_4+theme(legend.title = element_text(colour="Black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 12, face = "bold"))+
  theme(axis.text=element_text(colour="Black", size=10, face="plain"))+
  theme(axis.title=element_text(colour="Black", size=12, face="bold"))+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 10, face = "plain"))
violin_1_5
library(grid)
violin_1_6=violin_1_5+theme(panel.spacing = unit(0.5, "lines"),
                            legend.position="right")
violin_1_6

#remove gridlines and top/right borders and facets
violin_1_7=violin_1_6+ theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))
#theme(strip.text=element_blank())+
#theme(strip.background=element_blank())
violin_1_7

str(dop)

###Trying various transformations and polynomial models, all with random intercepts only due to
##replication issues
model1<-lmer(DOP_Blue~Sex*Angular_contrast+
               (1|Individual),data=dop,REML=FALSE)
summary(model1)
Anova(model1)
plot(model1)
##fitted vs residuals shows clear heteroskedascity and deviance at the extremes
r2beta(model1)
r.squaredGLMM(model1)

model2<-lmer(sqrt(DOP_Blue)~Sex*Angular_contrast+
               (1|Individual),data=dop,REML=FALSE)
summary(model2)
Anova(model2)
plot(model2)
## heteroskedascity at the lower fitted values is improved but strong deviance in higher fitted values

r2beta(model2)
r.squaredGLMM(model2)

model3<-lmer((DOP_Blue)^2~Sex*Angular_contrast+
               (1|Individual),data=dop,REML=FALSE)
summary(model3)
Anova(model3)
plot(model3)
## extreme heteroskedascity and deviance at the extremes, very bad fit

r2beta(model3)
r.squaredGLMM(model3)

model4<-lmer(log2(DOP_Blue)~Sex*Angular_contrast+
               (1|Individual),data=dop,REML=FALSE)
summary(model4)
Anova(model4)
plot(model4)
## heteroskedascity is fine but strong deviance at the extremes

r2beta(model4)
r.squaredGLMM(model4)

anova(model1,model2,model3,model4)
####sqrt Y transformation clearly better via AIC, use sqrt

##Now trying to fit polynomial models using the sqrt y transform
##Apparently when fitting polynomials, using the poly() function fits 'orthogonal polynomials'
##which accounts for the correlation better than using var+var^2+var^3
model5<-lmer(sqrt(DOP_Blue)~Sex*poly(Angular_contrast,2)+
               (1|Individual),data=dop,REML=FALSE)
summary(model5)
Anova(model5)
plot(model5)
##a little  heteroskedascity, deviance is OK but not amazing
r2beta(model5)
r.squaredGLMM(model5)


model6<-lmer(sqrt(DOP_Blue)~Sex*poly(Angular_contrast,3)+
               (1|Individual),data=dop,REML=FALSE)
summary(model6)
Anova(model6)
plot(model6)
##fitted vs residuals looks ok, traces of heteroskedascitiy and deviance, but acceptable
r2beta(model6)
r.squaredGLMM(model6)


anova(model5,model6)

##3rd degree polynomial a much better fit compared to 2nd degree. Use 3rd degree
##However, fitting 7 fixed effects and a random effect on 11 individuals is dicey.
##Will plan to report a full conditional and marginal R2 table in the supplemental to detail the
##model selection process. 
############################
##Making a saveable model details dataframe
#using the r.squaredGLMM function of the MuMIn package to get marginal and conditional R2
r2_1=r.squaredGLMM(model1)
r2_2=r.squaredGLMM(model2)
r2_3=r.squaredGLMM(model3)
r2_4=r.squaredGLMM(model4)
r2_5=r.squaredGLMM(model5)
r2_6=r.squaredGLMM(model6)

#Using the R2 function of the semEFF package to get the R2, adjusted R2 and predicted R2
R2_1=R2(model1)
R2_2=R2(model2)
R2_3=R2(model3)
R2_4=R2(model4)
R2_5=R2(model5)
R2_6=R2(model6)

R2a=data.frame(rbind(r2_1,r2_2,r2_3,r2_4,r2_5,r2_6))
R2a

R2b=data.frame(rbind(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6))
R2b

table=anova(model1,model2,model3,model4,model5,model6)

##Incredibly annoying workaround I had to hack out to call the model formula in a way that could be printed
form1=paste(toString(formula(model1)[2]),toString(formula(model1)[1]),toString(formula(model1)[3]))
form2=paste(toString(formula(model2)[2]),toString(formula(model2)[1]),toString(formula(model2)[3]))
form3=paste(toString(formula(model3)[2]),toString(formula(model3)[1]),toString(formula(model3)[3]))
form4=paste(toString(formula(model4)[2]),toString(formula(model4)[1]),toString(formula(model4)[3]))
form5=paste(toString(formula(model5)[2]),toString(formula(model5)[1]),toString(formula(model5)[3]))
form6=paste(toString(formula(model6)[2]),toString(formula(model6)[1]),toString(formula(model6)[3]))
formtable=data.frame(Formula=rbind(form1,form2,form3,form4,form5,form6))

formtable
modeltable=data.frame(Model=c("model1","model2","model3","model4","model5","model6"),
                        Func=c("y=x","sqrt(y)=x","(y)^2=x","log2(y)=x","sqrt(y)=x+X^2","sqrt(y)=x+x^2+x^3"),
                        Formula=formtable$Formula,
                        R2=R2b$X.r.squared.,Marg_R2=R2a$R2m,Cond_R2=R2a$R2c,Adj_R2=R2b$X.adj.r.squared.,
                        Pred_R2=R2b$X.pred.r.squared,
                        Df=table[1],AIC=table[2],BIC=table[3],loglik=table[4])
modeltable
# write.csv(modeltable,"model_testing_details_Blue_DOP_vs_Angular_contrast.csv")


##Now proceed with significance testing via pbkrtest


#dropping the 1st order interaction

####First, I am adding 3 columns to the dop dataframe that are the product of the poly() function
##so that I can more easily drop factors and their interactions for pbkrtest
##Strictly speaking, columns 1, 2, and 3 are the orthogonal additive contributions of linear(1)
##quadratic(2) and 3rd degree components (3) of a 3rd degree polynomial function of angular_contrast
dop$AngCon1=poly(dop$Angular_contrast,3)[,1]
dop$AngCon2=poly(dop$Angular_contrast,3)[,2]
dop$AngCon3=poly(dop$Angular_contrast,3)[,3]


##Now I am making a new model that is the same as model6 above, but with the new factors
model7=lmer(sqrt(DOP_Blue)~AngCon1*Sex+AngCon2*Sex+AngCon3*Sex+
                (1|Individual),data=dop,REML=FALSE)
summary(model6)
summary(model7)
Anova(model7)

#now dropping each main effect
model7.1=update(model7,~.-Sex:AngCon3)
model7.2=update(model7,~.-Sex:AngCon2)
model7.3=update(model7,~.-Sex:AngCon1)

model7.4=lmer(sqrt(DOP_Blue)~AngCon1+AngCon2+AngCon3+Sex+
              (1|Individual),data=dop,REML=FALSE)

model7.5=update(model7.4,~.-AngCon3)
model7.6=update(model7.4,~.-AngCon2)
model7.7=update(model7.4,~.-AngCon1)
model7.8=update(model7.4,~.-Sex)

#enabling parallelization
library(parallel)
(nc <- detectCores())
## Create clusters
cl <- makeCluster(rep("localhost", 10))

#testing each model against the next largest model
A_1=PBmodcomp(model7,model7.1,nsim=10000,cl=cl)
A_2=PBmodcomp(model7,model7.2,nsim=10000,cl=cl)
A_3=PBmodcomp(model7,model7.3,nsim=10000,cl=cl)


B_1=PBmodcomp(model7.4,model7.5,nsim=10000,cl=cl)
B_2=PBmodcomp(model7.4,model7.6,nsim=10000,cl=cl)
B_3=PBmodcomp(model7.4,model7.7,nsim=10000,cl=cl)
B_4=PBmodcomp(model7.4,model7.8,nsim=10000,cl=cl)

stopCluster(cl)
options(scipen=999999999)
mctab=data.frame(rbind(A_1$test,A_2$test,A_3$test,B_1$test,B_2$test,B_3$test,B_4$test))
mctab
warnings()
modcomptable=data.frame(Model=c("model7.1","model7.1","model7.2","model7.2","model7.3",
                                "model7.3","model7.5","model7.5","model7.6","model7.6","model7.7","model7.7","model7.8","model7.8"),
                        Test=c("LRT","PBtest","LRT","PBtest","LRT","PBtest","LRT","PBtest",
                               "LRT","PBtest","LRT","PBtest","LRT","PBtest"),
                        Statistic=mctab$stat,P_value=mctab$p.value,
                        Factor=c("Sex:AngCon3","Sex:AngCon3",
                                 "Sex:AngCon2","Sex:AngCon2","Sex:AngCon1","Sex:AngCon1","AngCon3","AngCon3","AngCon2","AngCon2","AngCon1","AngCon1","Sex","Sex"),
                        Significance=c("ns","ns","*","*","***","***","***","***","***","***","***","***","ns","ns"))
modcomptable
write.csv(modcomptable,"PBmodcomp_factor_significance_Blue_DOP_vs_Angular_contrast_POLY3.csv")



################################################################################################
##############################################################################
#NOW DOING THE UV ANALYSIS

model1<-lmer(DOP_UV~Sex*Angular_contrast+
               (1|Individual),data=dop,REML=FALSE)
summary(model1)
Anova(model1)
plot(model1)
##fitted vs residuals shows clear heteroskedascity and deviance at the extremes
r2beta(model1)
r.squaredGLMM(model1)

model2<-lmer(sqrt(DOP_UV)~Sex*Angular_contrast+
               (1|Individual),data=dop,REML=FALSE)
summary(model2)
Anova(model2)
plot(model2)
## heteroskedascity at the lower fitted values is improved but strong deviance in higher fitted values

r2beta(model2)
r.squaredGLMM(model2)

model3<-lmer((DOP_UV)^2~Sex*Angular_contrast+
               (1|Individual),data=dop,REML=FALSE)
summary(model3)
Anova(model3)
plot(model3)
## extreme heteroskedascity and deviance at the extremes, very bad fit

r2beta(model3)
r.squaredGLMM(model3)

model4<-lmer(log2(DOP_UV)~Sex*Angular_contrast+
               (1|Individual),data=dop,REML=FALSE)
summary(model4)
Anova(model4)
plot(model4)
## heteroskedascity is fine but strong deviance at the extremes

r2beta(model4)
r.squaredGLMM(model4)

anova(model1,model2,model3,model4)
####sqrt Y transformation clearly better via AIC, use sqrt

##Now trying to fit polynomial models using the sqrt y transform
##Apparently when fitting polynomials, using the poly() function fits 'orthogonal polynomials'
##which accounts for the correlation better than using var+var^2+var^3
model5<-lmer(sqrt(DOP_UV)~Sex*poly(Angular_contrast,2)+
               (1|Individual),data=dop,REML=FALSE)
summary(model5)
Anova(model5)
plot(model5)
##a little  heteroskedascity, deviance is OK but not amazing
r2beta(model5)
r.squaredGLMM(model5)


model6<-lmer(sqrt(DOP_UV)~Sex*poly(Angular_contrast,3)+
               (1|Individual),data=dop,REML=FALSE)
summary(model6)
Anova(model6)
plot(model6)
##fitted vs residuals looks barely acceptable, some heteroskedascitiy and deviance, but acceptable
r2beta(model6)
r.squaredGLMM(model6)


anova(model5,model6)
##3rd degree polynomial NOT a much better fit compared to 2nd degree. Use 2nd degree
##However, fitting 5 fixed effects and a random effect on 11 individuals is dicey.
##Will plan to report a full conditional and marginal R2 table in the supplemental to detail the
##model selection process. 

##Making a saveable model details dataframe 
#using the r.squaredGLMM function of the MuMIn package to get marginal and conditional R2
r2_1=r.squaredGLMM(model1)
r2_2=r.squaredGLMM(model2)
r2_3=r.squaredGLMM(model3)
r2_4=r.squaredGLMM(model4)
r2_5=r.squaredGLMM(model5)
r2_6=r.squaredGLMM(model6)

#Using the R2 function of the semEFF package to get the R2, adjusted R2 and predicted R2
R2_1=R2(model1)
R2_2=R2(model2)
R2_3=R2(model3)
R2_4=R2(model4)
R2_5=R2(model5)
R2_6=R2(model6)

R2a=data.frame(rbind(r2_1,r2_2,r2_3,r2_4,r2_5,r2_6))
R2a

R2b=data.frame(rbind(R2_1,R2_2,R2_3,R2_4,R2_5,R2_6))
R2b

table=anova(model1,model2,model3,model4,model5,model6)

##Incredibly annoying workaround I had to hack out to call the model formula in a way that could be printed
form1=paste(toString(formula(model1)[2]),toString(formula(model1)[1]),toString(formula(model1)[3]))
form2=paste(toString(formula(model2)[2]),toString(formula(model2)[1]),toString(formula(model2)[3]))
form3=paste(toString(formula(model3)[2]),toString(formula(model3)[1]),toString(formula(model3)[3]))
form4=paste(toString(formula(model4)[2]),toString(formula(model4)[1]),toString(formula(model4)[3]))
form5=paste(toString(formula(model5)[2]),toString(formula(model5)[1]),toString(formula(model5)[3]))
form6=paste(toString(formula(model6)[2]),toString(formula(model6)[1]),toString(formula(model6)[3]))
formtable=data.frame(Formula=rbind(form1,form2,form3,form4,form5,form6))

formtable
modeltable=data.frame(Model=c("model1","model2","model3","model4","model5","model6"),
                      Func=c("y=x","sqrt(y)=x","(y)^2=x","log2(y)=x","sqrt(y)=x+X^2","sqrt(y)=x+x^2+x^3"),
                      Formula=formtable$Formula,
                      R2=R2b$X.r.squared.,Marg_R2=R2a$R2m,Cond_R2=R2a$R2c,Adj_R2=R2b$X.adj.r.squared.,
                      Pred_R2=R2b$X.pred.r.squared,
                      Df=table[1],AIC=table[2],BIC=table[3],loglik=table[4])
modeltable
write.csv(modeltable,"model_testing_details_UV_DOP_vs_Angular_contrast.csv")




##
##Now proceed with significance testing via pbkrtest

poly(dop$Angular_contrast,2)[,2]

####I will only use the first two Angular contrast columns (linear and quadratic) that 
##I added using the poly() function so that I can more easily drop factors and their interactions
### for pbkrtest
##Strictly speaking, columns 1 and 2, are the orthogonal additive contributions of linear(1)
## and quadratic(2) contributions of a 2nd  degree polynomial function of angular_contrast

##Now I am making a new model that is the same as model5 above, but with the new factors
model7=lmer(sqrt(DOP_UV)~AngCon1*Sex+AngCon2*Sex+
              (1|Individual),data=dop,REML=FALSE)
summary(model5)
summary(model7)
Anova(model7)

#now dropping each main effect
model7.1=update(model7,~.-Sex:AngCon2)
model7.2=update(model7,~.-Sex:AngCon1)

model7.3=lmer(sqrt(DOP_UV)~AngCon1+AngCon2+Sex+
                (1|Individual),data=dop,REML=FALSE)

model7.5=update(model7.3,~.-AngCon2)
model7.6=update(model7.3,~.-AngCon1)
model7.7=update(model7.3,~.-Sex)

Anova(model7.3)
#enabling parallelization
library(parallel)
(nc <- detectCores())
## Create clusters
cl <- makeCluster(rep("localhost", 10))

#testing each model against the next largest model
A_1=PBmodcomp(model7,model7.1,nsim=10000,cl=cl)
A_2=PBmodcomp(model7,model7.2,nsim=10000,cl=cl)


B_1=PBmodcomp(model7.3,model7.5,nsim=10000,cl=cl)
B_2=PBmodcomp(model7.3,model7.6,nsim=10000,cl=cl)
B_3=PBmodcomp(model7.3,model7.7,nsim=10000,cl=cl)

stopCluster(cl)
options(scipen=999999999)
mctab=data.frame(rbind(A_1$test,A_2$test,B_1$test,B_2$test,B_3$test))
mctab
warnings()
modcomptable=data.frame(Model=c("model7.1","model7.1","model7.2","model7.2","model7.5","model7.5","model7.6","model7.6","model7.7","model7.7"),
                        Test=c("LRT","PBtest","LRT","PBtest","LRT","PBtest","LRT","PBtest",
                               "LRT","PBtest"),
                        Statistic=mctab$stat,P_value=mctab$p.value,
                        Factor=c("Sex:AngCon2","Sex:AngCon2","Sex:AngCon1","Sex:AngCon1","AngCon2","AngCon2","AngCon1","AngCon1","Sex","Sex"),
                        Significance=c("ns","ns","ns","ns","***","***","***","***","ns","ns"))
modcomptable
write.csv(modcomptable,"PBmodcomp_factor_significance_UV_DOP_vs_Angular_contrast_POLY2.csv")


#######################################
#code for dataset with robots and specimens
dop=read.csv('NY_CCT_DOP_data_robots_and_specimens_R_ready.csv')

dop$Sex <- factor(dop$Sex,levels = c("F", "M", "HNPB_depol", "HNPB_pol","3M_Blue_R374_depol","3M_Blue_R374_pol"))

dop2=subset(dop,Comment=="Fine")
dop3=subset(dop, Sex!='3M_Blue_R374_depol')
dop3=subset(dop3, Sex!='HNPB_depol')

library(ggplot2)

##Figure 3D (UV) and E (Blue)
violin_1=ggplot(
  data = dop
  , mapping = aes(
    y =(DOP_Blue)
    , x=(Angular_contrast)
    , colour=Sex
  )
)+
  labs(
    x = 'Angular_contrast'
    ,title= 'DOP-Blue by Angular_contrast in viewing and incident angle'
    , y = 'DOP-Blue'
  )+
  geom_point(alpha=0.5,
             mapping = aes(x=(Angular_contrast),
                           colour =Sex
             ), size=1.7
             
  )+
  # geom_line(alpha=0.3,
  #           mapping = aes(x=(Angular_contrast),
  #                         colour =Sex
  #           )
  #           , size=1.8
  # )+
  #  geom_violin(alpha=.4, draw_quantiles = c(0.25, 0.5, 0.75),width=1,
  #              mapping=aes(fill=factor(Sex),colour=factor(Sex)))+
  # stat_summary(fun.data=mean_se, geom="ribbon", alpha=0.15,mapping=aes(fill=Sex))+
  # stat_summary(fun.data=mean_se, geom="line", alpha=0.95,size=1,mapping=aes(color=Sex))+
  # stat_summary(fun.data=mean_se, geom="line", alpha=0.35,size=1,linetype='dotted',mapping=aes(group=Individual))+
  geom_smooth(method="loess",alpha=.1, se =TRUE,aes(fill=Sex))+
  geom_hline(yintercept=0,linetype='dashed')+
  facet_grid(.~.)
violin_1
violin_1_2=violin_1+scale_color_manual(values=c("#fdae61","#2c7bb6","#9E71A8","#632770","#669999","#226666"))+
  scale_fill_manual(values=c("#fdae61","#2c7bb6","#9E71A8","#632770","#669999","#226666"))+theme_bw()
violin_1_2
violin_1_3=violin_1_2

# violin_1_3=violin_1_2+scale_y_continuous(limits=c(-1.25,1.25))
violin_1_3
violin_1_4=violin_1_3
violin_1_4
#Resizing and moving the text in the graph
violin_1_5=violin_1_4+theme(legend.title = element_text(colour="Black", size=12, face="bold"))+
  theme(legend.text = element_text(colour="Black", size = 12, face = "bold"))+
  theme(axis.text=element_text(colour="Black", size=10, face="plain"))+
  theme(axis.title=element_text(colour="Black", size=12, face="bold"))+
  theme(axis.title.y=element_text(vjust=1))+
  theme(strip.text=element_text(colour="Black", size = 10, face = "plain"))
violin_1_5
library(grid)
violin_1_6=violin_1_5+theme(panel.spacing = unit(0.5, "lines"),
                            legend.position="right")
violin_1_6

#remove gridlines and top/right borders and facets
violin_1_7=violin_1_6+ theme(axis.line = element_line(colour = "black"),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.border=element_blank())+
  theme(axis.line = element_line(color = 'black'))
#theme(strip.text=element_blank())+
#theme(strip.background=element_blank())
violin_1_7





