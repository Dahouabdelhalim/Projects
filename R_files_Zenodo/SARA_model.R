
## =============================================================================
# Taghipoor Masoomeh1, Giger-Reverdin Sylvie1 (2020)
# UMR Modélisation Systémique Appliquée aux Ruminants (MoSAR), INRAE,
# AgroParisTech, Université Paris-Saclay, 75005 Paris, France

# Characterization of goats' response to a high concentrate diet: 
#   several time scales, one global index
# 
## =============================================================================          

rm(list = ls())
dev.off()

library(lattice)
library(plotrix)
library(ggplot2)
library(corrplot)#graphique de la correlation
library(Hmisc)     
library(lme4)
library(car)
library(emmeans)

# setwd("")

load("DataSG.Rdata")
load("SyntheticVar.Rdata")


#============ Graphs  ===============
#Two scales 

# # 
xyplot(Hplus*10^7 ~ Day2|GoatNumber,DataSG,
       scales=list(x=list(rot=45)),
       index.cond=list(c(6,7,8,3,4,5,1,2)),
       type=c("p","l"),pch=16,
       grid = TRUE,
       group = Time,auto.key = list(
         title="Post prandial time"
         , columns = 5
         , corner = c(1, 0)
         , space="top"))


   #============  Corr analysis  ===============

  mydata2 = SyntheticVar
  param = mydata2[,c(4:5,7:13)]
  colnames(param)= c("v0", "A","R","dur","AmpAC","vlast","a","b","c")
      
# Nul hyp: Gaussian distribution
  for(i in 1:length(param)){
    ist = shapiro.test(param[,i])
    print(ist$p.value)
  }
      
#normality is not respected , Spearman correaltion test
  cor.param = cor(param, method="spearman")
  sign.cor =  rcorr(as.matrix(param), type=c("spearman"))
  cor.pvalue = sign.cor$P
      
  corrplot(cor.param, type="upper", order="original", tl.col="black", tl.srt=45)
 
  #============== ================ ================ ================
  #======================  Mixed Model  
  #============== ================ ================ ================

# Definition of model parameters
      
  # str(mydata2)  
  ID = unique(mydata2$GoatNumber)
  mydata2$Day=as.numeric(mydata2$Day)
  
#long-term effect: week
  week = c("s2","s2","s3","s3","s3","s3","s4","s5","s6")
  mydata2$week = as.factor(rep(week,length(ID)))
      
#hort term effect: dayS3
  mydata2$Day.short = mydata2$Day - 2 #day before challenge <=0
  mydata2$dayS3 = mydata2$Day.short
  mydata2$dayS3[mydata2$week %in% c("s2","s4","s5","s6") ] = 0 

#Shapiro test for normality
# replace var by precovery, var0 and Amplitude
      
  var=mydata2$Amplitude     
  shapiro.test(var)
  power=powerTransform(var, family="bcPower") # pos value
  hist(var^(power$lambda))
  param = var^(power$lambda)

        mod0 = lmer(param~mydata2$dayS3 + (1|mydata2$GoatNumber),
                    REML = FALSE, na.action = na.omit) #Short term
        mod1 = lmer(param~mydata2$dayS3+mydata2$week+(1|mydata2$GoatNumber),
                    REML = FALSE, na.action = na.omit) #both
        mod2 = lmer( param ~  mydata2$week+(1|mydata2$GoatNumber),
                     REML = FALSE, na.action = na.omit) #long term
        
  anova(mod0,mod1)
  anova(mod2,mod1)
  # summary(mod2)
  emmeans(mod1, list(pairwise ~ week), adjust = "tukey") #weeks effect


# Normality and heterosedasticity of results
   plot(mod2,col=mydata2$GoatNumber,pch=param,id=0.05)
   qqmath(mod2,col=mydata2$GoatNumber,pch=param)   
   hist(residuals(mod2))
   shapiro.test(residuals(mod2))
# 

   
   #============== ================ ================ ================
    #==============  Quantification of goats response  ================
   #============== ================ ================ ================ 
 # Euclidien distance
     eucdis=function(v0,A){
       d=sqrt((v0[2]-v0[1])^2+(A[2]-A[1])^2)
       return(d) }

 
 d=c()
 dist = data.frame(Goat = double(),
                   d12 = double(),
                   d23 = double(),
                   d34 = double(),
                   d45 = double(),
                   d56 = double(),
                   d67 = double(),
                   d78 = double(),
                   d89 = double(),
                   stringsAsFactors=FALSE) 
 for(j in c(1:length(ID))){
   v=1e7*mydata2$var0[mydata2$GoatNumber==ID[j]]
   a=1e7*mydata2$Amplitude[mydata2$GoatNumber==ID[j]]
   for(i in c(1:8)){
     v1=v[i]
     v2=v[i+1]
     a1=a[i]
     a2=a[i+1]
     if(v1<v2 & a1<a2){weight=-2
     }else if(v1<v2 &a1>a2){weight=1
     }else if(v1>v2 &a1<a2){weight=-1
     }else if(v1>v2 &a1>a2){weight=2}
     d[i]=weight*eucdis(c(v1,v2),c(a1,a2))
   }
   dist[j,]=c(as.character(ID[j]),d)
 }
 
 #Non-weighted metric
 dist[, c(2:9)]=sapply(dist[, c(2:9)], as.numeric)
 dist$Index=rowSums(dist[,c(3:9)])
 
 #Weighted Metric
 dist1=dist
 dist1[, c(2)]=dist[,c(2)]*1/3
 dist1[, c(7)]=dist[,c(7)]*1/4
 dist1[, c(8)]=dist[,c(8)]*1/6
 dist1[, c(9)]=dist[,c(9)]*1/8
 
 dist1$Index=rowSums(dist1[,c(3:9)])
 dist1[,c(2:10)]=round(dist1[,c(2:10)],2)
 
 
 #==============  Plot of response trajectory  ================
 res=mydata2
 res$dayfac=as.factor(res$Day)
 
 ID = unique(as.factor(res$GoatNumber))
 
 
 plot( res$var0[res$GoatNumber==ID[1]]*10^7, res$Amplitude[res$GoatNumber==ID[1]]*10^7,type="l", col= 1,
       ylab=("A"), xlab=("v0"))
 text( res$var0[res$GoatNumber==ID[1]]*10^7, res$Amplitude[res$GoatNumber==ID[1]]*10^7,  res$dayfac,
       cex=2, pos=3,col="blue") 
 
 par(mfrow=c(2,4))
 for(i in 1:8){
   
   plot( res$var0[res$GoatNumber==ID[i]]*10^7, res$Amplitude[res$GoatNumber==ID[i]]*10^7, 
         col=1,type="l",ylab=("A"), xlab=("v0"),xlim=c(0,15),ylim=c(0,27),
         main=c(as.character(ID[i])),cex.main=2,cex.lab=1.5, lty=2)
   text(res$var0[res$GoatNumber==ID[i]]*10^7, res$Amplitude[res$GoatNumber==ID[i]]*10^7,  res$dayfac,
        cex=1.2, pos=3,col="blue" ) 
 }
 
 #=====
 
 