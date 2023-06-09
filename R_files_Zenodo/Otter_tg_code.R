setwd("~/Dropbox/Chapter 1/Writeup/Submission Version/RSOS/Dryad")
source("Data Prep.R")

# UNIVARIATE LOGISTIC MODELS ----------------------------------------------

varnames <- c("AC4","AC5","lengthper10cm","Dev10","Frm10","Grz1","For10","Wet1","Scrub10","Grass10","Othnat10","PD10","log2(PD)","PD5cat","PD4cat","HD10","log2(HD)","HD5cat","HD4cat","RD10","log2(zero.to.min(RD))","RD5cat","RD4cat","AREA2","imp10", "imp","Dev10")

data <- subset(datastore,SEX=="M")
results.univ.M <- univglms(y="toxo",varnames=varnames,data=datastore)
write.csv(results.univ.M,"results_univ_M.csv")
data <- subset(datastore,SEX=="F")
results.univ.F <- univglms(y="toxo",varnames=varnames,data=data)
write.csv(results.univ.F,"results_univ_F.csv")

# Combined Univariate table

varnames <- c("AC4","AC5","lengthper10cm","SEX","Dev10","Frm10","Grz1","For10","Wet1","Scrub10","Grass10","Othnat10","imp","log2(PD)","log2(HD)","log2(zero.to.min(RD))","AREA2")

data <- datastore
results.univ.both <- univglms("toxo",varnames=varnames,data=data)
write.csv(results.univ.both,"results_univ_both.csv")

# MULTIVARIATE MIXED EFFECTS MODELS ---------------------------------------
require(lme4)
glmer_landuse1 <- glmer(toxo~Dev10+Frm10+Grz1+(1|AREA2),binomial,data=data)
summary(glmer_landuse1)

glmer_landuse2 <- glmer(toxo~Dev10+Frm10+(1|AREA2),binomial,data=datastore)
summary(glmer_landuse2)

glmer_landuse3 <- glmer(toxo~Dev10*Frm10+Grz1+(1|AREA2),binomial,data=datastore)
summary(glmer_landuse3)

glmer_landuse4 <- glmer(toxo~Dev10+Frm10*Grz1+(1|AREA2),binomial,data=datastore)
summary(glmer_landuse4)

glmer_landuse5 <- glmer(toxo~SEX+AC4+Dev10+(1|AREA2),binomial,data=datastore)
summary(glmer_landuse5)

glmer1 <- glmer(toxo~SEX+AC4+(1|AREA2),binomial,data=datastore)
summary(glmer1)

glmer2 <- glmer(toxo~SEX+AC4+Frm10+(1|AREA2),binomial,data=datastore)
summary(glmer2)

glmer3 <- glmer(toxo~SEX+AC4+Dev10+Frm10+(1|AREA2),binomial,data=datastore)
summary(glmer3)

glmer4 <- glmer(toxo~SEX+AC4+Dev10*Frm10+(1|AREA2),binomial,data=datastore,nAGQ=4)
summary(glmer4)

glmer5 <- glmer(toxo~SEX+AC4+PD4cat+(1|AREA2),binomial,data=datastore)
summary(glmer5)

glmer6 <- glmer(toxo~SEX+AC4+log2(PD)+(1|AREA2),binomial,data=datastore)
summary(glmer6)

glmer7 <- glmer(toxo~SEX+AC4+HD4cat+(1|AREA2),binomial,data=datastore,nAGQ=4)
summary(glmer7)

glmer8 <- glmer(toxo~SEX+AC4+log2(HD)+(1|AREA2),binomial,data=datastore)
summary(glmer8)

glmer9 <- glmer(toxo~SEX*log2(HD)+AC4+(1|AREA2),binomial,data=datastore)
summary(glmer9)

# glmer11 <- glmer(toxo~SEX+AC4+Dev10+log2(HD)+(1|AREA2),binomial,data=datastore,nAGQ=25)
# summary(glmer11) Dev10 and HD are too correlated

glmer12 <- glmer(toxo~SEX+AC4+Frm10+log2(HD)+(1|AREA2),binomial,data=datastore)
summary(glmer12)

AIC(glmer6);AIC(glmer8);AIC(glmer12);

#glmer12 <- glmer(toxo~SEX+Dev10+log2(HD)+AC4+(1|AREA2),binomial,data=datastore
#summary(glmer12) # Don't use this model, as it conditions on an intermediate

glmer13 <- glmer(toxo~SEX+AC4+Devall+(1|AREA2),binomial,data=datastore)
summary(glmer13)

glmer14 <- glmer(toxo~SEX+AC4+Devlowplus+(1|AREA2),binomial,data=datastore)
summary(glmer14)

glmer15 <- glmer(toxo~SEX+AC4+Devmedplus+(1|AREA2),binomial,data=datastore)
summary(glmer15)

glmer16 <- glmer(toxo~SEX+AC4+Devhigh+(1|AREA2),binomial,data=datastore,nAGQ=4)
summary(glmer16)

# glmer17 <- glmer(toxo~SEX+AC4+Dev10*Frm10+Wet1+For10+(1|AREA2),binomial,data=datastore,nAGQ=25)
# summary(glmer17)

glmer18 <- glmer(toxo~SEX+AC4+Wet1+(1|AREA2),binomial,data=datastore,nAGQ=4)
summary(glmer18)

# glmer19 <- glmer(toxo~SEX+AC4+Scrub10+(1|AREA2),binomial,data=datastore,nAGQ=8)
# summary(glmer19) # NB: does not converge

glmer20 <- glmer(toxo~SEX+AC4+Grass10+(1|AREA2),binomial,data=datastore,nAGQ=4)
summary(glmer20)

glmer21 <- glmer(toxo~SEX+AC4+Dev10+(1|AREA2),binomial,data=datastore)
summary(glmer21)

glmer22 <- glmer(toxo~SEX+AC4+imp10+(1|AREA2),binomial,data=datastore)
summary(glmer22)

glmer23 <- glmer(toxo~SEX+AC4+Grz1+(1|AREA2),binomial,data=datastore)
summary(glmer23)

glmer24 <- glmer(toxo~SEX+AC4+For10+(1|AREA2),binomial,data=datastore)
summary(glmer24)

# glmer21 <- glmer(toxo~SEX+AC4+Othnat10+(1|AREA2),binomial,data=datastore,nAGQ=8)
# summary(glmer21) # NB: does not converge

# Model selection table ---------------------------------------------------

modlist <- list(glmer1,glmer2,glmer3,glmer4,glmer8,glmer12,glmer18,glmer_landuse5)
table1 <- model.sel.tab(modlist)
library(AICcmodavg)
names <- rep("NULL",length(modlist))
for(i in 1:length(modlist)){
  names[i] <- as.character(formula(modlist[[i]]))[3]
}
table2 <- aictab(modlist,modnames=names)

write.csv(table1,"Model Selection Table1.csv",row.names=F)
write.csv(table2,"Model Selection Table2.csv",row.names=F)

# Multivariate Tables ------------------------------------------------------

table <- multiv.tab.glmer(glmer6)
write.csv(table,"multiv_glmer6.csv")
table <- multiv.tab.glmer(glmer8)
write.csv(table,"multiv_glmer8.csv")

# Diet models -------------------------------------------------------------

# Minimum number of foraging bouts
minbouts=1
# Minimum number of dives
mindives=29

nott=length(datastore$OTTERNO)
data$N_bouts[which(is.na(data$N_bouts))]<-0
data$N_dives[which(is.na(data$N_dives))]<-0
data$hasdiet <- rep(T,length(datastore$OTTERNO))
for(i in 1:nott){
  if(data$N_bouts[i]<minbouts){data$hasdiet[i]<-F}
  if(data$N_dives[i]<mindives){data$hasdiet[i]<-F}
}

diet <- subset(data,data$hasdiet)
diet$BSMB <- 0
diet$BSMB[which(diet$AREA2=="BIGS")] <- 1
diet$BSMB[which(diet$AREA2=="MPEN")] <- 1
diet$BSMB[which(diet$AREA2=="MBAY")] <- 1
diet <- subset(diet,BSMB==1)

diet$AREA2 <- factor(diet$AREA2)

diet$snail10plus <- 0
diet$snail10plus[diet$snail>.1] <- 1

diet$snail5plus <- 0
diet$snail5plus[diet$snail>.05] <- 1
diet$AGE2 <- (diet$AGE)^2

glm_diet1 <- glm(toxo~SEX+AGE+snail,binomial,diet)
summary(glm_diet1)

glm_diet2 <- glm(toxo~SEX+snail,binomial,diet)
summary(glm_diet2)

glm_diet3 <- glm(toxo~SEX+snail5plus,binomial,diet)
summary(glm_diet3)

glm_diet4 <- glm(toxo~SEX+snail,binomial,diet)
summary(glm_diet4)

glm_diet5 <- glm(toxo~SEX+AGE+snail10plus,binomial,diet)
summary(glm_diet5)

glm_diet6 <- glm(toxo~AC5+snail10plus,binomial,diet)
summary(glm_diet6)

glm_diet7 <- glm(toxo~SEX+snail10plus,family=binomial,data=diet)
summary(glm_diet7)

glm_diet8 <- glm(toxo~SEX+log2(HD)+snail10plus,binomial,diet)
summary(glm_diet8)

glm_diet9 <- glm(toxo~SEX+imp+snail,binomial,diet)
summary(glm_diet9)

glm_diet10 <- glm(toxo~SEX+Dev10+snail,binomial,diet)
summary(glm_diet10)

glm_diet11 <- glm(toxo~Frm10,binomial,diet)
summary(glm_diet11)

modlist_diet <- list(glm_diet1,glm_diet2,glm_diet3,glm_diet4,glm_diet5,glm_diet6,glm_diet7,glm_diet8,glm_diet9,glm_diet10,glm_diet11)
model.sel.tab(modlist_diet)

write.csv(multiv.tab.glm(glm_diet7),"results_glm_diet.csv")

#### Model-based predictions ####

# Model-averaging prediction

predict.data <- read.csv("Predict_data.csv",header=T)
predict.data$AREA2[predict.data$AREA2=="PBSS"]<-"SNLO"

predict.data$Devall <- with(predict.data,rowSums(cbind(NLCDV21,NLCDV22,NLCDV23,NLCDV24)))
predict.data$Dev10 <- 10*(predict.data$Devall)
predict.data$Frm10 <- 10*(predict.data$NLCDV82)
predict.data$Grz1 <- 100*predict.data$NLCDV81
predict.data$For10 <- 10*(predict.data$NLCDV41+predict.data$NLCDV42)
predict.data$Wet1 <- 100*(predict.data$NLCDV90+predict.data$NLCDV95)
predict.data$Scrub10 <- 10*(predict.data$NLCDV52)
predict.data$Grass10 <- 10*(predict.data$NLCDV71)
predict.data$Othnat10 <- predict.data$Grass10+predict.data$Scrub10

newdata <- subset(predict.data,select=c("ATOS","AREA2","HD","Dev10","Frm10","Wet1","Devall"))
newdata$AC4 <- factor("A",levels=levels(data$AC4))
newdata$SEX <- factor("F",levels=levels(data$SEX))
output2f <- newdata
require(MuMIn)
avg <- model.avg(modlist)
output2f$P <- predict(avg,newdata,backtransform=T)
write.csv(output2f,"Predicted_Favg.csv")

newdata$SEX <- factor("M",levels=levels(data$SEX))
output2m <- newdata
output2m$P <- predict(avg,newdata,backtransform=T)
write.csv(output2m,"Predicted_Mavg.csv")


# GAM with non-parametric smoother - can give confidence bounds
# Spatially smoothed estimate of "observed prevalence"
require(mgcv)
calif_main <- subset(datastore,ATOS<1500)
atos.gam <- gam(toxo~s(ATOS),data=calif_main,binomial)

atos <- as.numeric(calif_main$ATOS)
atos_all <- seq(min(atos),max(atos))
newdat<-atos_all

pred <- predict(atos.gam,data.frame(ATOS=newdat),type="link",se=T)
gampred <- exp(pred$fit) / (1+exp(pred$fit))

plot(atos_all, gampred,type = "l", ylim = c(0, 1),xlab = "ATOS (km)", ylab = "Probability TG antibody positive")
lines(atos_all, exp(pred$fit+1.96*pred$se.fit) / (1 + exp(pred$fit + 1.96 * pred$se.fit)), lty = 2)
lines(atos_all, exp(pred$fit-1.96*pred$se.fit) / (1 + exp(pred$fit - 1.96 * pred$se.fit)), lty = 2)

# 1D kernel density estimate of prevalence

atos_samples <- calif_main$ATOS
nsamples <- length(atos_samples)
atos_positives <- calif_main$ATOS[calif_main$toxo==1]
npos <- length(atos_positives)

kd_samples <- density(atos_samples,adjust=1,from=min(atos_all),to=max(atos_all))
kd_positives <- density(atos_positives,adjust=1,from=min(atos_all),to=max(atos_all))

smooth_samples <- approx(x=kd_samples$x,y=kd_samples$y,xout=atos_all)
smooth_samples$y <- smooth_samples$y*length(atos_samples)

smooth_positives <- approx(x=kd_positives$x,y=kd_positives$y,xout=atos_all)
smooth_positives$y <- smooth_positives$y*length(atos_positives)

smooth_prevalence <- smooth_samples
smooth_prevalence$y <- smooth_positives$y/smooth_samples$y

plot(smooth_samples,type="l")
lines(smooth_positives,col="red")
plot(smooth_prevalence,type="l")

# Output results
require(xlsx)
est.prev <- data.frame(ATOS=atos_all,GAM=gampred,N=smooth_samples$y)
write.xlsx(est.prev,"Estimated_prev.xlsx",row.names=F)

# Estimate sampling density to gret out areas with inadequate sampling

kd_samplesize <- density(calif_main$ATOS,bw=10,from=min(atos_all),to=max(atos_all))
smooth_samplesize <- approx(x=kd_samplesize$x,y=kd_samplesize$y,xout=atos_all)
smooth_samplesize$y <- smooth_samplesize$y*length(calif_main$ATOS)

# Sex-specific GAM

atos.gam2 <- gam(toxo~SEX+s(ATOS),data=calif_main,binomial)
newdat<-data.frame(ATOS=atos_all,SEX="F")
pred.sex <- predict(atos.gam2,newdat,type="link",se=T)
gampred.sex <- exp(pred.sex$fit) / (1+exp(pred.sex$fit))

plot(atos_all, gampred.sex,type = "l", ylim = c(0, 1),xlab = "ATOS (km)", ylab = "Probability TG antibody positive")
lines(atos_all, exp(pred.sex$fit+1.96*pred.sex$se.fit) / (1 + exp(pred.sex$fit + 1.96 * pred.sex$se.fit)), lty = 2)
lines(atos_all, exp(pred.sex$fit-1.96*pred.sex$se.fit) / (1 + exp(pred.sex$fit - 1.96 * pred.sex$se.fit)), lty = 2)

# Output results
require(xlsx)
est.prev.F <- data.frame(ATOS=atos_all,GAM=gampred.sex,N=smooth_samplesize$y)
est.prev.F$N05 <- 1;est.prev.F$N05[which(est.prev.F$N<0.05)]<-0
est.prev.F$N10 <- 1;est.prev.F$N10[which(est.prev.F$N<0.10)]<-0
est.prev.F$GAM05 <- est.prev.F$GAM*est.prev.F$N05
est.prev.F$GAM10 <- est.prev.F$GAM*est.prev.F$N10

write.xlsx(est.prev.F,"Estimated_prev_F.xlsx",row.names=F)

# Table 2

sumtab <- read.csv("WS_data_summary.csv",header=T)
output <- sumtab
otdat <- datastore
X <- tapply(otdat$TGPOS,otdat$AREA2,FUN="sum",na.rm=T)
N <- tapply(otdat$TGPOS,otdat$AREA2,FUN="length")
sortorder <- order(names(X))
nsites <- length(output[,1])
output$X <- X[sortorder]
output$N <- N[sortorder]
output$P <- output$X/output$N
ci <- prop.ci(output$X,output$N)
output$cilo <- ci[,1];output$cihi <- ci[,2]
output$state <- c("CA","BC","AK","CA","CA","CA","BC","CA","CA","CA","WA","AK","AK")
list <- list(RD=output$RD,PD=output$PD,HD=output$HD,Forest=output$Forest,Wetland=output$Wetland,Developed=output$Developed,DevLowPlus=output$DevLowPlus,Impervious=output$Impervious,Cropping=output$Cropping,Pasture=output$Pasture,Anthropic=output$Anthropic)
res <- spearman.tab(list,output$P)
write.csv(res,"Table2.csv",row.names=F)
