#### 2015 copepod fungus co-occurrence ###

## 4-1-20

library(MASS)     

### COPEPOD and FUNGUS ###
# N = 15 sites * 6 colonies = 90
15*6 # 90
8+34+4+44 # 90 
# 8 with both; 34 with copepod only; 4 with fungus only; 44 neither

infect <- as.table(rbind(c(8,34), c(4,44))) 
dimnames(infect) <- list(Cop = c("Y","N"),
                         Fung = c("Y","N"))
infect
(Xsq <- chisq.test(infect))  # p= 0.2376; not significantly different from random co-occurrence; Chisq = 1.39
Xsq$observed   
Xsq$expected   

### PLOT ###

par(mfrow=c(1,1))
# http://damirah.cz/r-barplot-with-error-bars/ 
means_chi<-c(4,6.4,34,36.4,8,5.6,44,41.6) # F, Fexp, C, Cexp, Both, Bothexp,N, N exp, 
mp_chi<-barplot(means_chi,axes=FALSE,axisnames=FALSE,ylim=c(0,50),col=c("orchid","black","pink","black","purple","black","forestgreen","black"),density = c(100),cex.lab=1, xlab="Infection status",ylab="Frequency (# of coloneis)")
# add axes
axis(1,labels=c("F_obs", "F_exp","C_obs", "C_exp",  "FC_obs", "FC_exp", "N_obs", "N_exp"), at=mp_chi,cex.lab=.5)
axis(2,at=seq(0,50,by=10))
box()

Xsq$residuals  # Pearson residuals - both is a little higher, neither is a little higher
Xsq$stdres 
