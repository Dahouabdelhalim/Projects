#Open pedigree and phenotypic data
Young<-read.csv("Young_pheno.csv")
Pedigree<-read.csv("Pedigree.csv")

library(asreml)
ainv<-ainverse(Pedigree[,c(1,3,2)])

library(nadiv)
ginvD <- makeD(Pedigree)$Dinv
attr(ginvD,'INVERSE')=TRUE

#Prepare data for the models
Young$RearID<-as.factor(Young$RearID)
Young$Year<-as.factor(Young$Year)
Young$Mother<-as.factor(Young$Mother)
Young$Ring<-as.factor(Young$Ring)
Young$Ring2<-as.factor(Young$Ring2)

#Model for HA 
modelHA<-asreml(fixed=Agg_d16~ 1 + Sex  + handlingorder_d16
                , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD)+ Year + obs_d16+ Mother
                , data=Young 
                , na.action = na.method(x="include",y="omit" ))

summary(modelHA)$varcomp
wald.asreml(modelHA,ssType = "conditional", denDF="numeric")#Test significance of fixed effects
summary(modelHA,all=TRUE)$coef.fixed#Coefficients of fixed effects and their SE 

#Calculate variance proportions and their SEs using vpredict
vpredict(modelHA,Mother~V3/(V3+V4+V5+V6+V7))
vpredict(modelHA,RearID~V4/(V3+V4+V5+V6+V7))
vpredict(modelHA,Dom~V6/(V3+V4+V5+V6+V7))
vpredict(modelHA,Ped~V5/(V3+V4+V5+V6+V7))
vpredict(modelHA,Res~V7/(V3+V4+V5+V6+V7))
vpredict(modelHA,DomG~V6/(V5+V6))


#Model for BR 
modelBR<-asreml(fixed=BR~ 1 + Sex  + handlingorder_d16
                , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD)+ Year + obs_d16+ Mother
                , data=Young 
                , na.action = na.method(x="include",y="omit" ))

summary(modelBR)$varcomp

wald.asreml(modelBR,ssType = "conditional", denDF="numeric")$wald#Test significance of fixed effects
summary(modelBR,all=TRUE)$coef.fixed#Coefficients of fixed effects and their SE 

#Calculate variance proportions and their SEs using vpredict 
vpredict(modelBR,Mother~V3/(V3+V4+V5+V6+V7))
vpredict(modelBR,RearID~V4/(V3+V4+V5+V6+V7))
vpredict(modelBR,Dom~V6/(V3+V4+V5+V6+V7))
vpredict(modelBR,Ped~V5/(V3+V4+V5+V6+V7))
vpredict(modelBR,Res~V7/(V3+V4+V5+V6+V7))
vpredict(modelBR,DomG~V6/(V5+V6))

#Run model without dominance
modelBR.nodo<-asreml(fixed=BR~ 1 + Sex  + handlingorder_d16
                     , random= ~vm(Ring,ainv)+ RearID + Year + obs_d16+ Mother
                     , data=Young 
                     , na.action = na.method(x="include",y="omit" ))

1-pchisq(2*(modelBR$loglik-modelBR.nodo$loglik),1)#LRT test of dominance (1df)

#Model for Docility 
modelDoc<-asreml(fixed=docility~ 1 + Sex  + handlingorder_d16
                , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD)+ Year + obs_d16+ Mother
                , data=Young 
                , na.action = na.method(x="include",y="omit" ))


summary(modelDoc)$varcomp

wald.asreml(modelDoc,ssType = "conditional", denDF="numeric")#Test significance of fixed effects
summary(modelDoc,all=TRUE)$coef.fixed#Coefficients of fixed effects and their SE 

#Calculate variance proportions and their SEs using vpredict
vpredict(modelDoc,Mother~V3/(V3+V4+V5+V6+V7))
vpredict(modelDoc,RearID~V4/(V3+V4+V5+V6+V7))
vpredict(modelDoc,Dom~V6/(V3+V4+V5+V6+V7))
vpredict(modelDoc,Ped~V5/(V3+V4+V5+V6+V7))
vpredict(modelDoc,Res~V7/(V3+V4+V5+V6+V7))
vpredict(modelDoc,DomG~V6/(V5+V6))

#Run model without dominance
modelDoc.nodo<-asreml(fixed=docility~ 1 + Sex  + handlingorder_d16
                     , random= ~vm(Ring,ainv)+ RearID + Year + obs_d16+ Mother
                     , data=Young 
                     , na.action = na.method(x="include",y="omit" ))

1-pchisq(2*(modelDoc$loglik-modelDoc.nodo$loglik),1)#LRT test of dominance (1df)


#Model for tarsus length

modelTarsus<-asreml(fixed=Tarsus~ 1 + Sex  
                , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD)+ Year + obs_d16+ Mother
                , data=Young 
                , na.action = na.method(x="include",y="omit" ))

summary(modelTarsus)$varcomp
wald.asreml(modelTarsus,ssType = "conditional", denDF="numeric")#Test significance of fixed effects
summary(modelTarsus,all=TRUE)$coef.fixed#Coefficients of fixed effects and their SE 

#Calculate variance proportions and their SEs using vpredict
vpredict(modelTarsus,Mother~V3/(V3+V4+V5+V6+V7))
vpredict(modelTarsus,RearID~V4/(V3+V4+V5+V6+V7))
vpredict(modelTarsus,Dom~V6/(V3+V4+V5+V6+V7))
vpredict(modelTarsus,Ped~V5/(V3+V4+V5+V6+V7))
vpredict(modelTarsus,Res~V7/(V3+V4+V5+V6+V7))
vpredict(modelTarsus,DomG~V6/(V5+V6))#

#Run model without dominance
modelTarsus.nodo<-asreml(fixed=Tarsus~ 1 + Sex  
                         , random= ~vm(Ring,ainv)+ RearID + Year + obs_d16+ Mother
                         , data=Young 
                         , na.action = na.method(x="include",y="omit" ))

1-pchisq(2*(modelTarsus$loglik-modelTarsus.nodo$loglik),1)#LRT test of dominance (1df)

#Model for body mass
modelW<-asreml(fixed=Weight_d16~ 1 + Sex   + Tarsus
               , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD)+ Year + obs_d16+ Mother
               , data=Young 
               , na.action = na.method(x="include",y="omit" ))

summary(modelW)$varcomp
wald.asreml(modelW,ssType = "conditional", denDF="numeric")#Test significance of fixed effects
summary(modelW,all=TRUE)$coef.fixed#Coefficients of fixed effects and their SE 

#Calculate variance proportions and their SEs using vpredict
vpredict(modelW,Mother~V3/(V3+V4+V5+V6+V7))
vpredict(modelW,RearID~V4/(V3+V4+V5+V6+V7))
vpredict(modelW,Dom~V6/(V3+V4+V5+V6+V7))
vpredict(modelW,Ped~V5/(V3+V4+V5+V6+V7))
vpredict(modelW,Res~V7/(V3+V4+V5+V6+V7))
vpredict(modelW,DomG~V6/(V5+V6))#8.4%

#Run model without dominance
modelW.nodo<-asreml(fixed=Weight_d16~ 1 + Sex  +Tarsus
                    , random= ~vm(Ring,ainv)+ RearID +  Year + obs_d16+ Mother
                    , data=Young 
                    , na.action = na.method(x="include",y="omit" ))

1-pchisq(2*(modelW$loglik-modelW.nodo$loglik),1)#LRT test of dominance (1df)


#Model for wing length
modelWing<-asreml(fixed=Wing~ 1 + Sex  
                  , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD)+ Year + obs_d16+ Mother
                  , data=Young 
                  , na.action = na.method(x="include",y="omit" ))

summary(modelWing)$varcomp
wald.asreml(modelWing,ssType = "conditional", denDF="numeric")#Test significance of fixed effects
summary(modelWing,all=TRUE)$coef.fixed#Coefficients of fixed effects and their SE 

#Calculate variance proportions and their SEs using vpredict
vpredict(modelWing,Mother~V3/(V3+V4+V5+V6+V7))
vpredict(modelWing,RearID~V4/(V3+V4+V5+V6+V7))
vpredict(modelWing,Dom~V6/(V3+V4+V5+V6+V7))
vpredict(modelWing,Ped~V5/(V3+V4+V5+V6+V7))
vpredict(modelWing,Res~V7/(V3+V4+V5+V6+V7))
vpredict(modelWing,DomG~V6/(V5+V6))

#Run model without dominance
modelWing.nodo<-asreml(fixed=Wing~ 1 + Sex
                       , random= ~vm(Ring,ainv)+ RearID + Year + obs_d16+ Mother
                       , data=Young 
                       , na.action = na.method(x="include",y="omit" ))

1-pchisq(2*(modelWing$loglik-modelWing.nodo$loglik),1)#LRT test of dominance (1df)


#Model for head-bill length
modelHead<-asreml(fixed=Head~ 1 + Sex  
                  , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD)+ Year + obs_d16+ Mother
                  , data=Young 
                  , na.action = na.method(x="include",y="omit" ))

summary(modelHead)$varcomp
wald.asreml(modelHead,ssType = "conditional", denDF="numeric")#Test significance of fixed effects
summary(modelHead,all=TRUE)$coef.fixed#Coefficients of fixed effects and their SE 

#Calculate variance proportions and their SEs using vpredict
vpredict(modelHead,Mother~V3/(V3+V4+V5+V6+V7))
vpredict(modelHead,RearID~V4/(V3+V4+V5+V6+V7))
vpredict(modelHead,Dom~V6/(V3+V4+V5+V6+V7))
vpredict(modelHead,Ped~V5/(V3+V4+V5+V6+V7))
vpredict(modelHead,Res~V7/(V3+V4+V5+V6+V7))
vpredict(modelHead,DomG~V6/(V5+V6))

#Run model without dominance
modelHead.nodo<-asreml(fixed=Head~ 1 + Sex  
                       , random= ~vm(Ring,ainv)+ RearID + Year + obs_d16+ Mother
                       , data=Young 
                       , na.action = na.method(x="include",y="omit" ))

1-pchisq(2*(modelHead$loglik-modelHead.nodo$loglik),1)#LRT test of dominance (1df)

#Model for tail length
modelTail<-asreml(fixed=Tail~ 1 + Sex  
                  , random= ~vm(Ring,ainv)+ RearID + vm(Ring2,ginvD)+ Year + obs_d16+ Mother
                  , data=Young 
                  , na.action = na.method(x="include",y="omit" ))

summary(modelTail)$varcomp
wald.asreml(modelTail,ssType = "conditional", denDF="numeric")#Test significance of fixed effects
summary(modelTail,all=TRUE)$coef.fixed#Coefficients of fixed effects and their SE 

#Calculate variance proportions and their SEs using vpredict
vpredict(modelTail,Mother~V3/(V3+V4+V5+V6+V7))
vpredict(modelTail,RearID~V4/(V3+V4+V5+V6+V7))
vpredict(modelTail,Dom~V6/(V3+V4+V5+V6+V7))
vpredict(modelTail,Ped~V5/(V3+V4+V5+V6+V7))
vpredict(modelTail,Res~V7/(V3+V4+V5+V6+V7))
vpredict(modelTail,DomG~V6/(V5+V6))

#Plot variance estimates and their 95%CI for each trait

Variances<-matrix(NA,56,4)
Variances[,1]<-rep(c("HA","BR","Docility","Tarsus","Mass","Wing","Head","Tail"),each=7)
Variances[,2]<-rep(c("VO","VY","VM","VCE","VA","VD","VR"),8)
colnames(Variances)<-c("Trait","Level","Estimate","SE")
Variances<-as.data.frame(Variances)
Variances$Estimate<-c(summary(modelHA)$varcomp[,1],summary(modelBR)$varcomp[,1],summary(modelDoc)$varcomp[,1],
                      summary(modelTarsus)$varcomp[,1],summary(modelW)$varcomp[,1],summary(modelWing)$varcomp[,1],
                      summary(modelHead)$varcomp[,1],summary(modelTail)$varcomp[,1])
Variances$SE<-c(summary(modelHA)$varcomp[,2],summary(modelBR)$varcomp[,2],summary(modelDoc)$varcomp[,2],
                summary(modelTarsus)$varcomp[,2],summary(modelW)$varcomp[,2],summary(modelWing)$varcomp[,2],
                summary(modelHead)$varcomp[,2],summary(modelTail)$varcomp[,2])

Variances$CIsup<-Variances$Estimate+1.96*Variances$SE
Variances$CIinf<-Variances$Estimate-1.96*Variances$SE

library(wesanderson)

par(mfrow=c(2,4))
for (i in as.character(unique(Variances$Trait))){
  
  Vartrait<-Variances[Variances$Trait==i,]
  Vartrait<-Vartrait[c(5,6,3,4,7,2,1),]
  
  plot(Vartrait$Estimate, xlab=" ", ylim=c(min(Vartrait$CIinf,na.rm=T),max(Vartrait$CIsup,na.rm=T)),main=i,cex.main=2,
       ylab=" ",pch=21,cex=2,cex.axis=1.2,cex.lab=1.5, xaxt="n",
       bg=c(wes_palette("Moonrise3", 5, type = "discrete"),wes_palette("Moonrise2", 2, type = "discrete")),col=NULL)
  abline(0,0,lty=2)
  arrows(1:7,Vartrait$CIinf,1:7,Vartrait$CIsup,
         code=3,angle=90,length=0.1,col= c(wes_palette("Moonrise3", 5, type = "discrete"),wes_palette("Moonrise2", 2, type = "discrete")),lwd=2)
  axis(side=1,at=c(1:7),labels=c("VA","VD","VM","VCE","VR","VY","VO"),tick=TRUE)
  
}


#Compile estimates of variance ratios (V/VP) for each trait
estimatesHA<-c(as.numeric(vpredict(modelHA,Mother~V3/(V3+V4+V5+V6+V7))[1]),
               as.numeric(vpredict(modelHA,RearID~V4/(V3+V4+V5+V6+V7))[1]),
               as.numeric(vpredict(modelHA,Dom~V6/(V3+V4+V5+V6+V7))[1]),
               as.numeric(vpredict(modelHA,Ped~V5/(V3+V4+V5+V6+V7))[1]),
               as.numeric(vpredict(modelHA,Res~V7/(V3+V4+V5+V6+V7))[1]))

estimatesBR<-c(as.numeric(vpredict(modelBR,Mother~V3/(V3+V4+V5+V6+V7))[1]),
               as.numeric(vpredict(modelBR,RearID~V4/(V3+V4+V5+V6+V7))[1]),
               as.numeric(vpredict(modelBR,Dom~V6/(V3+V4+V5+V6+V7))[1]),
               as.numeric(vpredict(modelBR,Ped~V5/(V3+V4+V5+V6+V7))[1]),
               as.numeric(vpredict(modelBR,Res~V7/(V3+V4+V5+V6+V7))[1]))

estimatesDoc<-c(as.numeric(vpredict(modelDoc,Mother~V3/(V3+V4+V5+V6+V7))[1]),
                as.numeric(vpredict(modelDoc,RearID~V4/(V3+V4+V5+V6+V7))[1]),
                as.numeric(vpredict(modelDoc,Dom~V6/(V3+V4+V5+V6+V7))[1]),
                as.numeric(vpredict(modelDoc,Ped~V5/(V3+V4+V5+V6+V7))[1]),
                as.numeric(vpredict(modelDoc,Res~V7/(V3+V4+V5+V6+V7))[1]))

estimatesTars<-c(as.numeric(vpredict(modelTarsus,Mother~V3/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelTarsus,RearID~V4/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelTarsus,Dom~V6/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelTarsus,Ped~V5/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelTarsus,Res~V7/(V3+V4+V5+V6+V7))[1]))

estimatesWeight<-c(as.numeric(vpredict(modelW,Mother~V3/(V3+V4+V5+V6+V7))[1]),
                   as.numeric(vpredict(modelW,RearID~V4/(V3+V4+V5+V6+V7))[1]),
                   as.numeric(vpredict(modelW,Dom~V6/(V3+V4+V5+V6+V7))[1]),
                   as.numeric(vpredict(modelW,Ped~V5/(V3+V4+V5+V6+V7))[1]),
                   as.numeric(vpredict(modelW,Res~V7/(V3+V4+V5+V6+V7))[1]))

estimatesWing<-c(as.numeric(vpredict(modelWing,Mother~V3/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelWing,RearID~V4/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelWing,Dom~V6/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelWing,Ped~V5/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelWing,Res~V7/(V3+V4+V5+V6+V7))[1]))

estimatesHead<-c(as.numeric(vpredict(modelHead,Mother~V3/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelHead,RearID~V4/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelHead,Dom~V6/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelHead,Ped~V5/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelHead,Res~V7/(V3+V4+V5+V6+V7))[1]))

estimatesTail<-c(as.numeric(vpredict(modelTail,Mother~V3/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelTail,RearID~V4/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelTail,Dom~V6/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelTail,Ped~V5/(V3+V4+V5+V6+V7))[1]),
                 as.numeric(vpredict(modelTail,Res~V7/(V3+V4+V5+V6+V7))[1]))

estimates<-c(estimatesHA,estimatesBR,estimatesDoc,estimatesTars,estimatesWeight,estimatesWing,estimatesHead,estimatesTail)
estimates<-matrix(estimates,5,8,byrow = FALSE)

colnames(estimates)<-c("HA","BR","Docility","Tarsus","Mass","Wing","Head","Tail")
rownames(estimates)<-c("VM","VCE","VD","VA","VR")
estimates<-estimates[c("VA","VD","VM","VCE","VR"),]

#Plot variance proportions: stacked barplot
opar = par(oma = c(0,0,2,0))
barplot(estimates,beside=F,col=wes_palette("Moonrise3", 5, type = "discrete"))
opar =par(oma = c(0,0,2,0), mar = c(0,0,0,0), new = TRUE)
legend(x = "top", legend = rownames(estimates), fill = wes_palette("Moonrise3", 5, type = "discrete"), bty = "n", ncol = 5, inset = -0.12) 
par(opar) # reset par 


