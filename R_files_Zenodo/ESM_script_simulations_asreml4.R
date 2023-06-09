#Open pedigree and phenotypic data
Young<-read.csv("Young_pheno.csv")
Pedigree<-read.csv("Pedigree.csv")

library(asreml)
ainv<-ainverse(Pedigree[,c(1,3,2)])

library(nadiv)
ginvD <- makeD(Pedigree)$Dinv
attr(ginvD,'INVERSE')=TRUE

##Prepare data for the models
Young$RearID<-as.factor(Young$RearID)
Young$Year<-as.factor(Young$Year)
Young$Mother<-as.factor(Young$Mother)
Young$Ring<-as.factor(Young$Ring)
Young$Ring2<-as.factor(Young$Ring2)

#Prepare output tables
PvalsDom<-vector()#stores p values of the test: is dominance statistically different from simulated values
PvalsDom2<-vector()#stores p values of the test:is dominance statistically different from simulated values
Output<-data.frame( "VM"=NA, "VE"=NA, "VA"=NA, "VD"=NA, "VR"=NA)#Store model1 output
Output3<-data.frame("VM"=NA, "VE"=NA, "VA"=NA, "VR"=NA)#Store model3 output
Output4<-data.frame( "VE"=NA, "VA"=NA, "VD"=NA, "VR"=NA)#Store model4 output
Output5<-data.frame( "VE"=NA, "VA"=NA, "VR"=NA)#Store model5 output

RealVals<-data.frame( "VM"=NA, "VE"=NA, "VA"=NA, "VD"=NA, "VR"=NA)#Store real values

for (i in 1:1){
  
  #Generate breeding values
  A <- makeA(Pedigree)
  Gmat <- matrix(0.20, 1, 1)
  breedingValues <- grfx(n = 10946, G = Gmat, incidence = A, saveIncidence = FALSE)
  RealVals[i,"VA"]<-var(breedingValues)
  
  #Generate dominance values
  D <- makeD(Pedigree)$D
  Dmat <- matrix(0.03, 1, 1)
  dominanceValues <- grfx(n = 10946, G = Dmat, incidence = D, saveIncidence = FALSE)
  RealVals[i,"VD"]<-var(dominanceValues)
  
  #Generate random effects values
  Emat <- matrix(0.30, 1, 1)
  cfxE <- drfx(G = Emat, fac = "RearID", dataf = Young)
  
  Rmat <- matrix(0.44, 1, 1)
  cfxR <- drfx(G = Rmat, fac = "Ring", dataf = Young)
  
  Mmat <- matrix(0.03, 1, 1)
  cfxM <- drfx(G = Mmat, fac = "Mother", dataf = Young)
  
  RealVals[i,"VE"]<-var(cfxE$fx)
  RealVals[i,"VR"]<-var(cfxR$fx)
  RealVals[i,"VM"]<-var(cfxM$fx)
  
  #Calculate phenotypic data
  Genotype<-as.data.frame(dominanceValues+breedingValues)
  Genotype<-data.frame("Gval"=Genotype$V1,"Ring"=as.factor(as.character(Pedigree$id)))
  Genotype<-droplevels(Genotype[Genotype$Ring%in% Young$Ring,])
  
  RandEffects<-as.vector(cfxE$fx)+ as.vector(cfxR$fx)+ as.vector(cfxM$fx)
  Young$RandEffects<-RandEffects
  Youngb<-merge(Young,Genotype)
  
  Youngb$Phenotype<-Youngb$Gval+ Youngb$RandEffects
  Youngb$Phenotype<-ifelse(is.na (Youngb$BR)==TRUE,NA,Youngb$Phenotype)
  
  #Run the "Full model"
  modeltest<-asreml(fixed=Phenotype~ 1 
                    , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD) + Mother
                    , data=Youngb
                    , na.action = na.method(x="include",y="omit" ))
  
  #Constrain dominance to be equal to the simulated value
  sv<-asreml(fixed=Phenotype~ 1 
             , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD) + Mother
             , data=Youngb
             , na.action = na.method(x="include",y="omit" ),start.values=T)
  gam<-sv$vparameters.table
  gam$Value[3]<-var(dominanceValues)
  gam$Constraint[3]<-"F"
  
  #Run the constrained full model
  modeltest2<-asreml(fixed=Phenotype~ 1 
                     , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD) + Mother
                     , data=Youngb
                     , na.action = na.method(x="include",y="omit"),G.param=gam)
  
  PvalsDom[i]<-1-pchisq(2*abs(modeltest$loglik-modeltest2$loglik),1)#Test if dominance differs from the simulated value
  
  #Remove dominance
  modeltest3<-asreml(fixed=Phenotype~ 1 
                     , random= ~vm(Ring,ainv)+ RearID  + Mother
                     , data=Youngb
                     , na.action = na.method(x="include",y="omit" ))
  
  PvalsDom2[i]<-1-pchisq(2*abs(modeltest$loglik-modeltest3$loglik),1)#Test if dominance differs from zero (1df)
  
  #Remove maternal effects
  modeltest4<-asreml(fixed=Phenotype~ 1 
                     , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD)
                     , data=Youngb
                     , na.action = na.method(x="include",y="omit" ))
  
  #Remove both maternal effects and dominance
  modeltest5<-asreml(fixed=Phenotype~ 1 
                     , random= ~vm(Ring,ainv)+ RearID 
                     , data=Youngb
                     , na.action = na.method(x="include",y="omit" ))
  
  Output[i,1:5]<-summary(modeltest)$varcomp[,1]
  Output3[i,1:4]<-summary(modeltest3)$varcomp[,1]
  Output4[i,1:4]<-summary(modeltest4)$varcomp[,1]
  Output5[i,1:3]<-summary(modeltest5)$varcomp[,1]

}


length(PvalsDom[PvalsDom<0.05])/length(PvalsDom)#% divergence between dominance variance and simulated value
length(PvalsDom2[(PvalsDom2/2)<0.05])/length(PvalsDom2)#% probability to detect dominance variance 


#Plot the results
Output[,6:10]<-NA
names(Output)[6:10]<-paste(names(Output)[1:5],"b",sep="")
Output[,6:10]<-Output[,1:5]-RealVals[1:5]

Output3[,5:8]<-Output3[,1:4]-RealVals[c(1:3,5)]
names(Output3)[5:8]<-paste(names(Output3)[1:4],"b",sep="")

Output4[,5:8]<-Output4[,1:4]-RealVals[c(2:5)]
names(Output4)[5:8]<-paste(names(Output4)[1:4],"b",sep="")

Output5[,4:6]<-Output5[,1:3]-RealVals[c(2,3,5)]
names(Output5)[4:6]<-paste(names(Output5)[1:3],"b",sep="")


#Plot mean estimates full model
estimates<-apply(Output[,c("VAb","VDb","VMb","VEb","VRb")],2,median)
CIinf<-as.numeric(c(quantile(Output$VAb,probs=c(0.025)),quantile(Output$VDb,probs=c(0.025)),
                    quantile(Output$VMb,probs=c(0.025)),quantile(Output$VEb,probs=c(0.025)),quantile(Output$VRb,probs=c(0.025))))

CIsup<-as.numeric(c(quantile(Output$VAb,probs=c(0.975)),quantile(Output$VDb,probs=c(0.975)),
                    quantile(Output$VMb,probs=c(0.975)),quantile(Output$VEb,probs=c(0.975)),quantile(Output$VRb,probs=c(0.975))))

#Plot mean estimates model w/o dominance
estimates3<-apply(Output3[,c("VAb","VMb","VEb","VRb")],2,median)
CIinf3<-as.numeric(c(quantile(Output3$VAb,probs=c(0.025)),quantile(Output3$VMb,probs=c(0.025)),
                     quantile(Output3$VEb,probs=c(0.025)),quantile(Output3$VRb,probs=c(0.025))))
CIsup3<-as.numeric(c(quantile(Output3$VAb,probs=c(0.975)),quantile(Output3$VMb,probs=c(0.975)),
                     quantile(Output3$VEb,probs=c(0.975)),quantile(Output3$VRb,probs=c(0.975))))

#Plot mean estimates model w/o mother
estimates4<-apply(Output4[,c("VAb","VDb","VEb","VRb")],2,median)
CIinf4<-as.numeric(c(quantile(Output4$VAb,probs=c(0.025)),quantile(Output4$VDb,probs=c(0.025)),
                     quantile(Output4$VEb,probs=c(0.025)),quantile(Output4$VRb,probs=c(0.025))))
CIsup4<-as.numeric(c(quantile(Output4$VAb,probs=c(0.975)),quantile(Output4$VDb,probs=c(0.975)),
                     quantile(Output4$VEb,probs=c(0.975)),quantile(Output4$VRb,probs=c(0.975))))

#Plot mean estimates model w/o mother and dominance
estimates5<-apply(Output5[,c("VAb","VEb","VRb")],2,median)
CIinf5<-as.numeric(c(quantile(Output5$VAb,probs=c(0.025)),quantile(Output5$VEb,probs=c(0.025)),
                     quantile(Output5$VRb,probs=c(0.025))))
CIsup5<-as.numeric(c(quantile(Output5$VAb,probs=c(0.975)),quantile(Output5$VEb,probs=c(0.975)),
                     quantile(Output5$VRb,probs=c(0.975))))


####Calculate heritability

#Full model
Output$h2<-Output$VA/rowSums(Output[,c("VA","VD","VM","VE","VR")])
h2_full<-median(Output$h2)
h2_full_CIinf<-as.numeric(quantile(Output$h2,probs=c(0.025)))
h2_full_CIsup<-as.numeric(quantile(Output$h2,probs=c(0.975)))


# Model _D
Output3$h2<-Output3$VA/rowSums(Output3[,c("VA","VM","VE","VR")])
h2_D<-median(Output3$h2)
h2_D_CIinf<-as.numeric(quantile(Output3$h2,probs=c(0.025)))
h2_D_CIsup<-as.numeric(quantile(Output3$h2,probs=c(0.975)))

# Model _M
Output4$h2<-Output4$VA/rowSums(Output4[,c("VA","VD","VE","VR")])
h2_M<-median(Output4$h2)
h2_M_CIinf<-as.numeric(quantile(Output4$h2,probs=c(0.025)))
h2_M_CIsup<-as.numeric(quantile(Output4$h2,probs=c(0.975)))

# Model _D_M
Output5$h2<-Output5$VA/rowSums(Output5[,c("VA","VE","VR")])
h2_DM<-median(Output5$h2)
h2_DM_CIinf<-as.numeric(quantile(Output5$h2,probs=c(0.025)))
h2_DM_CIsup<-as.numeric(quantile(Output5$h2,probs=c(0.975)))

#Combine all estimates
estimates.tot<-c(estimates,estimates3,estimates4,estimates5)
CIinf.tot<-c(CIinf,CIinf3,CIinf4,CIinf5)
CIsup.tot<-c(CIsup,CIsup3,CIsup4,CIsup5)
culr<-c(c("#85D4E3","#F4B5BD", "#9C964A", "#CDC08C",  "#FAD77B" ),
        c("#85D4E3", "#9C964A", "#CDC08C",  "#FAD77B" ),
        c("#85D4E3","#F4B5BD",  "#CDC08C",  "#FAD77B" ),
        c("#85D4E3", "#CDC08C",  "#FAD77B"))
x<- c(seq(0.5,2.5,0.5),seq(4.5,6,0.5),seq(8,9.5,0.5), seq(11.5,12.5,0.5))

#Combine all estimates of heritability
estimates.h2.tot<-c(h2_full,h2_D,h2_M,h2_DM)
CIinf.h2.tot<-c(h2_full_CIinf,h2_D_CIinf,h2_M_CIinf,h2_DM_CIinf)
CIsup.h2.tot<-c(h2_full_CIsup,h2_D_CIsup,h2_M_CIsup,h2_DM_CIsup)
x2<- c(1.5,5.25,8.75,12)

RealVals$h2<-RealVals$VA/rowSums(RealVals[,c("VA","VD","VM","VE","VR")])
mean(RealVals$h2)

#Plot estimates
par(mfrow=c(3,1),oma=c(0,0,0,0),mar=c(0,6,3,6))
layout(matrix(c(1,1,2), nrow = 3, ncol = 1, byrow = TRUE))
plot(x,estimates.tot,pch=21,xaxt="n",xlim=c(0,13),ylim=c(min(CIinf.tot),max(CIsup.tot)),xlab=" ",ylab="Error in variance",main="",cex=2,bg=culr,col=NULL,cex.axis=1.5,cex.lab=1.5)
arrows(x,CIinf.tot,x,CIsup.tot,angle=90,code=3,length=0.1,lwd=2,col=culr)
abline(0,0,col="black",lty=2,lwd=2)

#add legends
par(xpd=TRUE)
legend(x="top",c("VA","VD","VM","VCE","VR"),bty="n",
       fill=c("#85D4E3","#F4B5BD", "#9C964A", "#CDC08C",  "#FAD77B"),
       cex=1.5,pt.cex=1.5,ncol = 5,x.intersp = 0.2,inset = -0.12)

#Plot heritabilities
par(new=T)
par(mfrow=c(3,1),oma=c(2,0,1,0),mgp=c(3,1,0),mar=c(2,6,1,6))

plot(x2,estimates.h2.tot,pch=21,xaxt="n",xlim=c(0,13),ylim=c(min(CIinf.h2.tot),max(CIsup.h2.tot)),xlab=" ",ylab="Heritability",main="",cex=2,cex.axis=1.5,cex.lab=1.5,bg="black",col=NULL)
arrows(x2,CIinf.h2.tot,x2,CIsup.h2.tot,angle=90,code=3,length=0.1,lwd=2,col="black")
axis(side=1,at=x2,labels=c("Full model","Model-D","Model-M", "Model-D-M"),cex.axis=1.5,tick=FALSE)
par(xpd=FALSE)
abline(mean(RealVals$h2),0,col="black",lty=2,lwd=2)

#Look at correlations between variance components
round(cor(Output),2)

