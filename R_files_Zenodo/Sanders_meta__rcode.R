### Load packages
library("metafor") # for effect size estimation 
library("MCMCglmm") # for meta-analysis
library("dplyr") # for data manipulation
library ("dmetar") # for p-curve
library ("ggplot2") # for the figure

### Read in data
ALANdat1<-Sanders_et_al_datafile #1304 observations

### subset data (excluding not used measures, see methods)
ALANdat=subset(ALANdat1, Delete %in% c("Keep")) %>% droplevels() #1109 observations (198 removed)
###Check for subcategories that have been included
levels(ALANdat$Sub.Category)

###estimate effect sizes 
dat.n<-escalc(measure= "SMDH", n1i = Experimental_N, n2i = Control_N, m1i = Experimental_M, m2i = Control_M, sd1i = Experimental_SD, sd2i = Control_SD, data = ALANdat, append = TRUE)
### adjust for sign of effect size for fitness and physiology (see methods)
dat.n$yi=dat.n$yi*dat.n$Sign
ALAN<- dat.n

###set parameters for MCMCglmm 
nitt=15000;burnin=5000;thin=5

#We use flat priors for the residuals and
#for study as random factor. The random term idh(SE):units was fixed to
#one so that all measurement errors could be considered as independent of
#each other
#this prior is the inverse-Gamma priors (V = 1, nu = 0.002)
prior=list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, fix=1), G2=list(V=1, nu=0.002)))

# Meta- analysis by subcategories
#Physiology
p1<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Gene.expression"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
p2<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Hormone"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
p3<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Immune.response"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
p4<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Stress.response"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
p5<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Gland.structure"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
#Phenology     
m2<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Category=="Phenology"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
#Activity      
n1<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Activity.offset"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
n2<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Activity.onset"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
n3<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="diurnAct"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
n4<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="noctAct"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
# LHT
l1<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Seafinding.turtles"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
l2<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Predation.risk"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
l3<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Size"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
l4<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Cognition"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
l5<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Feeding"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
l6<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Predation"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
l7<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Reproductive.output"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
#Community     
c1<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Abundance"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
c2<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Bat.activity"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions
c3<-summary(MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=subset(ALAN, Sub.Category=="Diversity"),verbose = F, prior=prior, nitt=nitt,burnin=burnin,thin=thin))$solutions


#check for auto-correlation and the plots
#for variance and random components. This is an example and has been done for all models including the ones checking for publication bias

plot(m1$Sol)
plot(m1$VCV)

autocorr(m1$Sol[,1])
autocorr(m1$VCV[,1])

######## organise data for forest plot 
Forestpl.subcat=as.data.frame(rbind(n1,n2,n3,n4,c1,c2,c3,m2,p1,p2,p3,p4,p5,l1,l2,l3,l4,l5,l6,l7))
colnames(Forestpl.subcat)[c(2,3)]<-c("lCI","uCI")
Forestpl.subcat
Forestpl.subcat$category=c("Activity.offset","Activity.onset","diurnAct","noctAct","Phenology","Abundance" ,"Bat.activity","Diversity","Gene.expression","Hormone",
                           "Immune.response","Stress.response","Gland.structure","Seafinding.turtles" ,"Predation.risk","Size","Cognition","Feeding","Predation","Reproductive.output")
Forestpl.subcat$group.colour=c(rep("#FDAE61",4),"gray64",rep("#D7191C",3), rep("#2C7BB6",5),rep("#ABD9E9",7))
Forestpl.subcat$category<-factor(Forestpl.subcat$category, levels=c("Gene.expression","Hormone","Immune.response",
                                                                    "Stress.response","Gland.structure","Phenology",
                                                                    "Seafinding.turtles" ,"Predation.risk","Size","Cognition",
                                                                    "Feeding","Predation","Reproductive.output","Activity.offset",
                                                                    "Activity.onset","diurnAct","noctAct","Abundance" ,"Bat.activity","Diversity"))

#plot Figure 2, subcategory forest plot
str(Forestpl.subcat)
label <- Forestpl.subcat$category
mean  <- Forestpl.subcat$post.mean
lower <- Forestpl.subcat$lCI
upper <- Forestpl.subcat$uCI
group <- Forestpl.subcat$category
group.col <- Forestpl.subcat$group.colour
df <- data.frame(label,group,group.col, mean, lower, upper)
df$label<-ordered(df$label, levels=c("Diversity","Bat.activity","Abundance","noctAct","diurnAct",
                                     "Activity.onset","Activity.offset","Reproductive.output","Predation",
                                     "Feeding","Cognition","Size","Predation.risk","Seafinding.turtles",
                                     "Phenology","Gland.structure","Stress.response","Immune.response",
                                     "Hormone","Gene.expression"))

fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(col=group.col,size=1) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  theme(axis.title.x = element_text(size=30),axis.text.x  = element_text(size=20))+
  theme(axis.title.y = element_text(size=30),axis.text.y = element_text(size=20))+
  xlab("")+ ylab("Hedges'd (95% Credible Interval)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.background = element_rect(fill = "white", colour = "#6D9EC1",
                                        size = 2, linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey90")) 

print(fp)

##Check for publication bias 
#Gene expression
ALAN_gen=filter (ALAN, Sub.Category == "Gene.expression")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_gen,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_gen$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_gen$Precision<-sqrt(1/ALAN_gen$vi);ALAN_gen$MAR<-ALAN_gen$yi-ALAN_gen$Prediction ;ALAN_gen$zMAR<-ALAN_gen$MAR*ALAN_gen$Precision
pcurve(metagen(yi,vi, data = ALAN_gen))
plot(1/ALAN_gen$vi,ALAN_gen$zMAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_gen,nitt=nitt, thin=thin, burnin=burnin))

#Hormone
ALAN_ho=filter (ALAN, Sub.Category == "Hormone")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_ho,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_ho$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_ho$Precision<-sqrt(1/ALAN_ho$vi);ALAN_ho$MAR<-ALAN_ho$yi-ALAN_ho$Prediction ;ALAN_ho$zMAR<-ALAN_ho$MAR*ALAN_ho$Precision
plot(1/ALAN_ho$vi,ALAN_ho$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
pcurve(metagen(yi,vi, data = ALAN_ho))
plot(sqrt(1/ALAN_ho$vi),ALAN_ho$zMAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_ho,nitt=nitt, thin=thin, burnin=burnin))

#Immune response 
ALAN_imm=filter (ALAN, Sub.Category == "Immune.response")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_imm,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_imm$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_imm$Precision<-sqrt(1/ALAN_imm$vi);ALAN_imm$MAR<-ALAN_imm$yi-ALAN_imm$Prediction ;ALAN_imm$zMAR<-ALAN_imm$MAR*ALAN_imm$Precision
pcurve(metagen(yi,vi, data = ALAN_imm))
plot(1/ALAN_imm$vi,ALAN_imm$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_imm,nitt=nitt, thin=thin, burnin=burnin))

#Stress response
ALAN_stress=filter (ALAN, Sub.Category == "Stress.response")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_stress,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_stress$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_stress$Precision<-sqrt(1/ALAN_stress$vi);ALAN_stress$MAR<-ALAN_stress$yi-ALAN_stress$Prediction ;ALAN_stress$zMAR<-ALAN_stress$MAR*ALAN_stress$Precision
pcurve(metagen(yi,vi, data = ALAN_stress))
plot(1/ALAN_stress$vi,ALAN_stress$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_stress,nitt=nitt, thin=thin, burnin=burnin))

# Gland structure
ALAN_gs=filter (ALAN, Sub.Category == "Gland.structure")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_gs,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_gs$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_gs$Precision<-sqrt(1/ALAN_gs$vi);ALAN_gs$MAR<-ALAN_gs$yi-ALAN_gs$Prediction ;ALAN_gs$zMAR<-ALAN_gs$MAR*ALAN_gs$Precision
pcurve(metagen(yi,vi, data = ALAN_gs))
plot(1/ALAN_gs$vi,ALAN_gs$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_gs,nitt=nitt, thin=thin, burnin=burnin))

###Phenology
ALAN_phe=filter (ALAN, Sub.Category == "Phenology")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_phe,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_phe$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_phe$Precision<-sqrt(1/ALAN_phe$vi);ALAN_phe$MAR<-ALAN_phe$yi-ALAN_phe$Prediction ;ALAN_phe$zMAR<-ALAN_phe$MAR*ALAN_phe$Precision
pcurve(metagen(yi,vi, data = ALAN_phe))
plot(1/ALAN_phe$vi,ALAN_phe$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_phe,nitt=nitt, thin=thin, burnin=burnin))

#Seafinding turtles
ALAN_st=filter (ALAN, Sub.Category == "Seafinding.turtles")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_st,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_st$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_st$Precision<-sqrt(1/ALAN_st$vi);ALAN_st$MAR<-ALAN_st$yi-ALAN_st$Prediction ;ALAN_st$zMAR<-ALAN_st$MAR*ALAN_st$Precision
pcurve(metagen(yi,vi, data = ALAN_st))
plot(1/ALAN_st$vi,ALAN_st$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_st,nitt=nitt, thin=thin, burnin=burnin))

#Predation risk
ALAN_prr=filter (ALAN, Sub.Category == "Predation.risk")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_prr,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_prr$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_prr$Precision<-sqrt(1/ALAN_prr$vi);ALAN_prr$MAR<-ALAN_prr$yi-ALAN_prr$Prediction ;ALAN_prr$zMAR<-ALAN_prr$MAR*ALAN_prr$Precision
pcurve(metagen(yi,vi, data = ALAN_prr))
plot(1/ALAN_prr$vi,ALAN_prr$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_prr,nitt=nitt, thin=thin, burnin=burnin))

#Size
ALAN_siz=filter (ALAN, Sub.Category == "Size")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_siz,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_siz$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_siz$Precision<-sqrt(1/ALAN_siz$vi);ALAN_siz$MAR<-ALAN_siz$yi-ALAN_siz$Prediction ;ALAN_siz$zMAR<-ALAN_siz$MAR*ALAN_siz$Precision
pcurve(metagen(yi,vi, data = ALAN_siz))
plot(1/ALAN_siz$vi,ALAN_siz$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_siz,nitt=nitt, thin=thin, burnin=burnin))

#Cognition 
ALAN_cog=filter (ALAN, Sub.Category == "Cognition")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_cog,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_cog$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_cog$Precision<-sqrt(1/ALAN_cog$vi);ALAN_cog$MAR<-ALAN_cog$yi-ALAN_cog$Prediction ;ALAN_cog$zMAR<-ALAN_cog$MAR*ALAN_cog$Precision
pcurve(metagen(yi,vi, data = ALAN_cog))
plot(1/ALAN_cog$vi,ALAN_cog$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_cog,nitt=nitt, thin=thin, burnin=burnin))

#Feeding
ALAN_fee=filter (ALAN, Sub.Category == "Feeding")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_fee,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_fee$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_fee$Precision<-sqrt(1/ALAN_fee$vi);ALAN_fee$MAR<-ALAN_fee$yi-ALAN_fee$Prediction ;ALAN_fee$zMAR<-ALAN_fee$MAR*ALAN_fee$Precision
pcurve(metagen(yi,vi, data = ALAN_fee))
plot(1/ALAN_fee$vi,ALAN_fee$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_fee,nitt=nitt, thin=thin, burnin=burnin))

#Predation
ALAN_pred=filter (ALAN, Sub.Category == "Predation")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_pred,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_pred$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_pred$Precision<-sqrt(1/ALAN_pred$vi);ALAN_pred$MAR<-ALAN_pred$yi-ALAN_pred$Prediction ;ALAN_pred$zMAR<-ALAN_pred$MAR*ALAN_pred$Precision
pcurve(metagen(yi,vi, data = ALAN_pred))
plot(1/ALAN_pred$vi,ALAN_pred$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_pred,nitt=nitt, thin=thin, burnin=burnin))

#Offspring
ALAN_Off=filter (ALAN, Sub.Category == "Reproductive.output")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_Off,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_Off$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_Off$Precision<-sqrt(1/ALAN_Off$vi);ALAN_Off$MAR<-ALAN_Off$yi-ALAN_Off$Prediction ;ALAN_Off$zMAR<-ALAN_Off$MAR*ALAN_Off$Precision
pcurve(metagen(yi,vi, data = ALAN_Off))
plot(1/ALAN_Off$vi,ALAN_Off$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_Off,nitt=nitt, thin=thin, burnin=burnin))

#Activity offset 
ALAN_Aof=filter (ALAN, Sub.Category == "Activity.offset")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_Aof,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_Aof$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_Aof$Precision<-sqrt(1/ALAN_Aof$vi);ALAN_Aof$MAR<-ALAN_Aof$yi-ALAN_Aof$Prediction ;ALAN_Aof$zMAR<-ALAN_Aof$MAR*ALAN_Aof$Precision
pcurve(metagen(yi,vi, data = ALAN_Aof))
plot(1/ALAN_Aof$vi,ALAN_Aof$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_Aof,nitt=nitt, thin=thin, burnin=burnin))

#Activity onset 
ALAN_Aon=filter (ALAN, Sub.Category == "Activity.onset")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_Aon,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_Aon$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_Aon$Precision<-sqrt(1/ALAN_Aon$vi);ALAN_Aon$MAR<-ALAN_Aon$yi-ALAN_Aon$Prediction ;ALAN_Aon$zMAR<-ALAN_Aon$MAR*ALAN_Aon$Precision
pcurve(metagen(yi,vi, data = ALAN_Aon))
plot(1/ALAN_Aon$vi,ALAN_Aon$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_Aon,nitt=nitt, thin=thin, burnin=burnin))

#dirunal Activity 
ALAN_dA=filter (ALAN, Sub.Category == "diurnAct")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_dA,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_dA$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_dA$Precision<-sqrt(1/ALAN_dA$vi);ALAN_dA$MAR<-ALAN_dA$yi-ALAN_dA$Prediction ;ALAN_dA$zMAR<-ALAN_dA$MAR*ALAN_dA$Precision
pcurve(metagen(yi,vi, data = ALAN_dA))
plot(1/ALAN_dA$vi,ALAN_dA$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_dA,nitt=nitt, thin=thin, burnin=burnin))

# nocturnal Activity 
ALAN_nA=filter (ALAN, Sub.Category == "noctAct")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_nA,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_nA$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_nA$Precision<-sqrt(1/ALAN_nA$vi);ALAN_nA$MAR<-ALAN_nA$yi-ALAN_nA$Prediction ;ALAN_nA$zMAR<-ALAN_nA$MAR*ALAN_nA$Precision
pcurve(metagen(yi,vi, data = ALAN_nA))
plot(1/ALAN_nA$vi,ALAN_nA$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_nA,nitt=nitt, thin=thin, burnin=burnin))

#Abundance
ALAN_Ab=filter (ALAN, Sub.Category == "Abundance")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_Ab,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_Ab$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_Ab$Precision<-sqrt(1/ALAN_Ab$vi);ALAN_Ab$MAR<-ALAN_Ab$yi-ALAN_Ab$Prediction ;ALAN_Ab$zMAR<-ALAN_Ab$MAR*ALAN_Ab$Precision
pcurve(metagen(yi,vi, data = ALAN_Ab))
plot(1/ALAN_Ab$vi,ALAN_Ab$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_Ab,nitt=nitt, thin=thin, burnin=burnin))

#Bat Activity 
ALAN_Ba=filter (ALAN, Sub.Category == "Bat.activity")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_Ba,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_Ba$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_Ba$Precision<-sqrt(1/ALAN_Ba$vi);ALAN_Ba$MAR<-ALAN_Ba$yi-ALAN_Ba$Prediction ;ALAN_Ba$zMAR<-ALAN_Ba$MAR*ALAN_Ba$Precision
pcurve(metagen(yi,vi, data = ALAN_Ba))
Bat<-plot(1/ALAN_Ba$vi,ALAN_Ba$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_Ba,nitt=nitt, thin=thin, burnin=burnin))

#Diversity 
ALAN_Dv=filter (ALAN, Sub.Category == "Diversity")
m1<-MCMCglmm(yi~1, random=~us(sqrt(vi)):units+Study, data=ALAN_Dv,verbose = F, pr=T, prior=prior, nitt=nitt,burnin=burnin,thin=thin)
ALAN_Dv$Prediction<-predict(m1, marginal=~us(sqrt(vi)):units, type = "terms");ALAN_Dv$Precision<-sqrt(1/ALAN_Dv$vi);ALAN_Dv$MAR<-ALAN_Dv$yi-ALAN_Dv$Prediction ;ALAN_Dv$zMAR<-ALAN_Dv$MAR*ALAN_Dv$Precision
pcurve(metagen(yi,vi, data = ALAN_Dv))
plot(1/ALAN_Dv$vi,ALAN_Dv$MAR,pch=19,cex=1.25,cex.axis=1.25,cex.lab=1.5, ylab="",xlab="",col=rgb(58/255,95/255,205/255, 0.5));abline(a=0,b=0, lwd=1, lty=2)
summary(MCMCglmm(zMAR~Precision,family="gaussian",verbose=FALSE,data=ALAN_Dv,nitt=nitt, thin=thin, burnin=burnin))


