library(ape); library(phytools); library(MCMCglmm)
Mult=1;NITT=260000;THIN=200;BURN=60000

Tree<-read.newick("full.Tree.txt")
Tree$node.label<-NULL; Tree<-collapse.singles(Tree); inv.phylo<-inverseA(Tree,nodes="TIPS",scale=TRUE)

#### Prior Set ####
#Prior 2 fixes (effectively) the estimates for the coefficients to the a priori value
prior.default<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

prior.offset2<-list(B=list(mu=matrix(c(0,-.25),2),V=diag(2)*1e7),G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
prior.offset2$B$V[2,2]<-1e-7

#### Model Comparison Runs ####
m1.lym<-MCMCglmm(ln10lympho~1,random=~phylo,
                 ginverse=list(phylo=inv.phylo$Ainv),
                 nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult,
                 family="gaussian", prior=prior.default, data=CompositeData3,  verbose=F)
save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/DownsScaling/offset1.RData")

m2.lym<-MCMCglmm(ln10lympho~ln10Mass,random=~phylo,
                 ginverse=list(phylo=inv.phylo$Ainv),
                 nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult,
                 family="gaussian", prior=prior.default, data=CompositeData3,  verbose=F)
save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/DownsScaling/offset2.RData")

m3off.lym<-MCMCglmm(ln10lympho~ln10Mass,random=~phylo,
                    ginverse=list(phylo=inv.phylo$Ainv),
                    nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult,
                    family="gaussian",prior=prior.offset2, data=CompositeData3, verbose=F)
save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/DownsScaling/offset3.RData")

m1.neu<-MCMCglmm(ln10neutro~1,random=~phylo,
                 ginverse=list(phylo=inv.phylo$Ainv),
                 nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult,
                 family="gaussian", prior=prior.default, data=CompositeData3,  verbose=F)
save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/DownsScaling/offset5.RData")

m2.neu<-MCMCglmm(ln10neutro~ln10Mass,random=~phylo,
                 ginverse=list(phylo=inv.phylo$Ainv),
                 nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult,
                 family="gaussian", prior=prior.default, data=CompositeData3,  verbose=F)
save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/DownsScaling/offset6.RData")

m3off.neu<-MCMCglmm(ln10neutro~ln10Mass,random=~phylo,
                    ginverse=list(phylo=inv.phylo$Ainv),
                    nitt=NITT*Mult,thin=THIN*Mult,burnin=BURN*Mult,
                    family="gaussian",prior=prior.offset2, data=CompositeData3, verbose=F)
save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/DownsScaling/offset7.RData")

#### Full Models ####

m5.lym<-MCMCglmm(ln10lympho~ln10Mass+as.numeric(MaxReproEffort)+as.numeric(MaxLong)+
                   as.numeric(Diet.Type.Trophic.Level)+as.numeric(Sociality)+ln10Mass:as.numeric(MaxReproEffort)+
                   ln10Mass:as.numeric(MaxLong)+ln10Mass:as.numeric(Diet.Type.Trophic.Level)+ln10Mass:as.numeric(Sociality),random=~phylo,
                 family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),prior=prior.default,
                 data=CompositeData3, nitt=260000, thin=200, burnin=60000, verbose=F)
save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/DownsScaling/phyloNEW.RData")

m5.neu<-MCMCglmm(ln10neutro~ln10Mass+as.numeric(MaxReproEffort)+as.numeric(MaxLong)+
                   as.numeric(Diet.Type.Trophic.Level)+as.numeric(Sociality)+ln10Mass:as.numeric(MaxReproEffort)+
                   ln10Mass:as.numeric(MaxLong)+ln10Mass:as.numeric(Diet.Type.Trophic.Level)+ln10Mass:as.numeric(Sociality), random=~phylo,
                 family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),prior=prior.default,
                 data=CompositeData3, nitt=260000, thin=200, burnin=60000, verbose=F)
save.image("C:/Users/Dochtermann/Dropbox/Working/Projects/DownsScaling/phyloNEW.RData")

#### Model Summaries ####
summary(m1.lym)
summary(m2.lym)
summary(m3off.lym)
summary(m5.lym)

summary(m1.neu)
summary(m2.neu)
summary(m3off.neu)
summary(m5.neu)

#### DIC Values ####
m1.lym$DIC
m2.lym$DIC
m3off.lym$DIC
m5.lym$DIC

m1.neu$DIC
m2.neu$DIC
m3off.neu$DIC
m5.neu$DIC

-281.076-m1.lym$DIC
-281.076-m2.lym$DIC
-281.076-m3off.lym$DIC

-339.1428-m1.neu$DIC
-339.1428-m2.neu$DIC
-339.1428-m3off.neu$DIC

#### Pagel and R2 calculations ####
#Lym
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(m2.lym$Sol[i,] %*% t(m2.lym$X)))
  vmVarF[i]<-Var}
v.lym.F<-vmVarF

R2m<-vmVarF/(vmVarF+m2.lym$VCV[,1]+m2.lym$VCV[,2])
Pagel<-m2.lym$VCV[,1]/(m2.lym$VCV[,1]+m2.lym$VCV[,2])
unadjPagel<-m2.lym$VCV[,1]/(vmVarF+m2.lym$VCV[,1]+m2.lym$VCV[,2])

posterior.mode(R2m)
HPDinterval(R2m) 
posterior.mode(unadjPagel)
HPDinterval(unadjPagel)  
posterior.mode(Pagel)
HPDinterval(Pagel) 

vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(m5.lym$Sol[i,] %*% t(m5.lym$X)))
  vmVarF[i]<-Var}
v.lym.F<-vmVarF

R2m<-vmVarF/(vmVarF+m5.lym$VCV[,1]+m5.lym$VCV[,2])
Pagel<-m5.lym$VCV[,1]/(m5.lym$VCV[,1]+m5.lym$VCV[,2])
unadjPagel<-m5.lym$VCV[,1]/(vmVarF+m5.lym$VCV[,1]+m5.lym$VCV[,2])

posterior.mode(R2m)
HPDinterval(R2m) 
posterior.mode(unadjPagel)
HPDinterval(unadjPagel)  
posterior.mode(Pagel)
HPDinterval(Pagel) 

#Neu
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(m2.neu$Sol[i,] %*% t(m2.neu$X)))
  vmVarF[i]<-Var}
v.neu.F<-vmVarF

R2m<-vmVarF/(vmVarF+m2.neu$VCV[,1]+m2.neu$VCV[,2])
Pagel<-m2.neu$VCV[,1]/(m2.neu$VCV[,1]+m2.neu$VCV[,2])
unadjPagel<-m2.neu$VCV[,1]/(vmVarF+m2.neu$VCV[,1]+m2.neu$VCV[,2])

posterior.mode(R2m)
HPDinterval(R2m) 
posterior.mode(unadjPagel)
HPDinterval(unadjPagel)  
posterior.mode(Pagel)
HPDinterval(Pagel) 

vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(m5.neu$Sol[i,] %*% t(m5.neu$X)))
  vmVarF[i]<-Var}
v.neu.F<-vmVarF

R2m<-vmVarF/(vmVarF+m5.neu$VCV[,1]+m5.neu$VCV[,2])
Pagel<-m5.neu$VCV[,1]/(m5.neu$VCV[,1]+m5.neu$VCV[,2])
unadjPagel<-m5.neu$VCV[,1]/(vmVarF+m5.neu$VCV[,1]+m5.neu$VCV[,2])

posterior.mode(R2m)
HPDinterval(R2m) 
posterior.mode(unadjPagel)
HPDinterval(unadjPagel)  
posterior.mode(Pagel)
HPDinterval(Pagel) 

#### Stacked Barplot ####
lym.res<-posterior.mode(m2.lym$VCV[,2])
neu.res<-posterior.mode(m2.neu$VCV[,2])
lym.F<-posterior.mode(v.lym.F)
neu.F<-posterior.mode(v.neu.F)
lym.lambda<-posterior.mode(m2.lym$VCV[,1])
neu.lambda<-posterior.mode(m2.neu$VCV[,1])
lym.T<-lym.res+lym.F+lym.lambda
neu.T<-neu.res+neu.F+neu.lambda

prop.var<-matrix(c(lym.res/lym.T,neu.res/neu.T,
                   lym.F/lym.T,neu.F/neu.T,
                   lym.lambda/lym.T,neu.lambda/neu.T),
                 3,2,byrow=T)

colnames(prop.var)<-c("Lymphocytes","Neutrophils")

labels3<-c(expression("V"["R"]),expression("V"["F"]),expression("V"[lambda]))
labels4<-c(expression(atop('Lymphocytes','log'[10]*'(10'^9*'cells/L)')),           expression(atop('Neutrophils','log'[10]*'(10'^9*'cells/L)')))
#tiff('phylo.controlled.barNEW.tiff',res=600,width=7.5,height=7.5,units='in')
par(mfrow=c(1,1),mar=c(7.1,7.1,1,1))
barplot(prop.var, ylab="Proportion of Variation Explained",
        xaxt="n", col=c("black","gray","white"),
        ylim=c(0,1),legend = labels3,cex.lab=1.5,cex.axis=1.25,cex.names=1.5) 
axis(1,at=c(0.75,2),labels=labels4, line=2,cex.lab=2,lty=0,cex.axis=1.5)
#dev.off()


#Estimate fixed effects seperately for Mass + rest of fixed effects
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(m5.lym$Sol[i,c(1,3:10)] %*% t(m5.lym$X[,c(1,3:10)])))
  vmVarF[i]<-Var}
v.lym.NoMass<-as.mcmc(vmVarF)

vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(m5.lym$Sol[i,2] %*% t(m5.lym$X[,2])))
  vmVarF[i]<-Var}
v.lym.Mass<-as.mcmc(vmVarF)

vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(m5.neu$Sol[i,c(1,3:10)] %*% t(m5.neu$X[,c(1,3:10)])))
  vmVarF[i]<-Var}
v.neu.NoMass<-as.mcmc(vmVarF)

vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(m5.neu$Sol[i,2] %*% t(m5.neu$X[,2])))
  vmVarF[i]<-Var}
v.neu.Mass<-as.mcmc(vmVarF)

lym.res<-posterior.mode(m5.lym$VCV[,2])
lym.lambda<-posterior.mode(m5.lym$VCV[,1])
lym.FMass<-posterior.mode(v.lym.Mass)
lym.FNoMass<-posterior.mode(v.lym.NoMass)
lym.T2<-lym.res+lym.FMass+lym.FNoMass+lym.lambda

neu.res<-posterior.mode(m5.neu$VCV[,2])
neu.lambda<-posterior.mode(m5.neu$VCV[,1])
neu.FMass<-posterior.mode(v.neu.Mass)
neu.FNoMass<-posterior.mode(v.neu.NoMass)
neu.T2<-neu.res+neu.FMass+neu.FNoMass+neu.lambda

prop.var2<-matrix(c(lym.res/lym.T2,neu.res/neu.T2,
                   lym.FNoMass/lym.T2,neu.FNoMass/neu.T2,
                   lym.FMass/lym.T2,neu.FMass/neu.T2,
                   lym.lambda/lym.T2,neu.lambda/neu.T2),
                 4,2,byrow=T)

labels3a<-c(expression("V"["R"]),expression("V"["F, other"]),expression("V"["mass"]),expression("V"[lambda]))

#tiff('phylo.controlled.barFULL.tiff',res=600,width=7.5,height=7.5,units='in')
pdf('phylo.controlled.bar.pdf',useDingbats=F)
par(mfrow=c(1,1),mar=c(7.1,7.1,1,1))
barplot(prop.var2, ylab="Proportion of Variation Explained",
        xaxt="n", col=c("black","gray50","gray80","white"),args.legend =list(bg="white"),
        ylim=c(0,1),legend = labels3a,cex.lab=1.5,cex.axis=1.25,cex.names=1.5) 
axis(1,at=c(0.75,2),labels=labels4, line=2,cex.lab=2,lty=0,cex.axis=1.5)
#dev.off()

#### Line plots ####

newdat<-data.frame(x=seq(0,10,length=20))
mm<-model.matrix(~x,newdat)

v.lym<- var(m2.lym$Sol)
v.neu<- var(m2.neu$Sol)
fixef.lym<-posterior.mode(m2.lym$Sol)
fixef.neu<-posterior.mode(m2.neu$Sol)

newdat$lym<-mm%*%fixef.lym
newdat$neu<-mm%*%fixef.neu

pvar.lym<-diag(mm %*% tcrossprod(cov(m2.lym$Sol),mm))
pvar.neu<-diag(mm %*% tcrossprod(cov(m2.neu$Sol),mm))

newdat <- data.frame(
  newdat
  , plo.lym = newdat$lym-1.96*sqrt(pvar.lym)
  , phi.lym = newdat$lym+1.96*sqrt(pvar.lym)
  , plo.neu = newdat$neu-1.96*sqrt(pvar.neu)
  , phi.neu = newdat$neu+1.96*sqrt(pvar.neu)
)

polygon.x <- c(newdat$x, rev(newdat$x))
polygon.lym.y <- c(newdat$plo.lym, rev(newdat$phi.lym))
polygon.neu.y <- c(newdat$plo.neu, rev(newdat$phi.neu))

#tiff('phylo.controlled.sNEWW.tiff',res=600,width=7.5,height=15,units='in')
#something weird is going on with the sizing for mtext and the plots themselves, had to eyeball it
pdf('phylo.controlled.line1.pdf',height=21,width=7,useDingbats=F)
par(mfrow=c(3,1),mar=c(8,8,2,2))
par(pty="s") 
plot(ln10lympho~ln10neutro,data=CompositeData3,
     cex=2,cex.lab=2.25,cex.axis=1.75,
     xlab=" ",
     ylab=expression(atop('Lymphocytes','log'[10]*'(10'^9*'cells/L)')))
mtext(expression(atop('Neutrophils','log'[10]*'(10'^9*'cells/L)')), side=1, line=7,cex=1.5)

plot(ln10lympho~ln10Mass,data=CompositeData3, xlab=" ",
     ylab=expression(atop('Lymphocytes','log'[10]*'(10'^9*'cells/L)')),cex=2.25,cex.lab=2,cex.axis=1.75)
polygon(x=polygon.x, y=polygon.lym.y, col=adjustcolor("darkgray", alpha.f=0.4), border=NA)
lines(newdat$x,newdat$lym,col="black",lty=2,lwd=2)
mtext(expression(atop('Mass','log'[10]*'(g)')), side=1, line=7,cex=1.5)

plot(ln10neutro~ln10Mass,data=CompositeData3, xlab=" ",
     ylab=expression(atop('Neutrophils','log'[10]*'(10'^9*'cells/L)')),cex=2.25,cex.lab=2,cex.axis=1.75)
polygon(x=polygon.x, y=polygon.neu.y, col=adjustcolor("darkgray", alpha.f=0.4), border=NA)
lines(newdat$x,newdat$neu,col="black",lty=2,lwd=2)
mtext(expression(atop('Mass','log'[10]*'(g)')), side=1, line=7,cex=1.5)
#dev.off()