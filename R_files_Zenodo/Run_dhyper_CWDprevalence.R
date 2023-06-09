#Application 1 Estimating prevalence

#Run_dhyper_CWDprevalence.R

require(R2jags)
require(coda)

source("NordfjellaPrevalenceModel.R")
#model: "dhyper_3AgeClass_prev2016.bug"     
#model: "dhyper_3Ageclass_2016prior_prev2017.bug"
      
### PREPARING JAGS MODEL RUNS; 
# MCMC settings
niter <- 300000
nthin <- 3
nburn <- 200000
nchains <- 3

source("Modelling_dSe.R")
betaDist <- dSE_betadist  #Use one of 3 variants of dSE_betadist (Estimating_dSe.R)
#betaDist<-dSE_betadistExp24 #Use one of 3 variants of dSE_betadist (Estimating_dSe.R)

#Notice - if testing the variant dSE_betadist_longI3, slight changes are needed in the two models in NordfjellaPrevalenceModel.R. Specifically,
#Remove the line: dSeobexRLN1 ~ dbeta(Se.betaRLN1[1],Se.betaRLN1[2]) 
#Set: dSeAvg[1] <- 0

Se.betaRLN1<-c(betaDist["dRLNobexC",1],betaDist["dRLNobexC",2])
Se.betaRLN2<-c(betaDist["dRLNobexY",1],betaDist["dRLNobexY",2])
Se.betaRLN3<-c(betaDist["dRLNobex",1],betaDist["dRLNobex",2])
Se.betaO1<-c(betaDist["dObexC",1],betaDist["dObexC",2])
Se.betaO2<-c(betaDist["dObexY",1],betaDist["dObexY",2])
Se.betaO3<-c(betaDist["dObex",1],betaDist["dObex",2])

p0obex <-betaDist[c(4:6),"PrZeros"]
p0obexRLN <-betaDist[c(1:3),"PrZeros"]
p0obex.RLNpos <-betaDist[c(4:6),"Pr_Obex0.RLNpos"]

#########################
# Bundle data
#########################

#Nr of calves, yearlings and adults
antRein<-c(480,298,1329)
sd_Rein<- c(15.0,24.5,19.6)

my<-c(0,0,2)   #nr of cwd positive
mk<-c(81,23,187) #nr of animals tested
pRLN = 0.33  #nr of samples tested with RLN tissue included

#The prior distributions of prevalence ratios between age classes
#Alt1Y from  APPENDIX S5
betapar2<-c("2.65","3.50") #mode=0.33
xrange2<-0.8
x2min<-0

#Alt1C from  APPENDIX S5
betapar1<-c("1.40","4.60") #mode=0.05
xrange1<-0.5
x1min<-0

# Nordfjella sone 1
prev2016.data.n1 <- list(
                      My = my,   #nr of animals tested cwd-positive
                      Mk = mk,     #nr of animals tested
                      antM = antRein,    #Nr of animals per age class (estimated from pop model)
                      sd_antM = sd_Rein,  
                      N=2107,      #Pop size
                      sd_N=37.9,   #sd Pop size                      
                      pRLN = pRLN,  #nr of samples tested with RLN tissue included
                                   
                      #beta-distribution - prevalence ratio between age classes (1-calf, 2-yearling)                                                    
                      Se.betaPrev1 = c(betapar1[[1]][1],betapar1[[2]][1]),
                      Se.betaPrev2 = c(betapar2[[1]][1],betapar2[[2]][1]),
                      xrange2=xrange2,
                      xrange1=xrange1,
                      x2min=x2min,
                      x1min=x1min,
                      
                      #beta-distribution - dSe for detectable stages
                      Se.betaRLN1 = Se.betaRLN1[1:2],
                      Se.betaRLN2= Se.betaRLN2[1:2],
                      Se.betaRLN3 = Se.betaRLN3[1:2],
                      Se.betaO2= Se.betaO2[1:2],
                      Se.betaO3 = Se.betaO3[1:2],
                      
                      #Probability of infected animals being in 
                      #early infection stages  
                      p0obexRLN  =  p0obexRLN,         #prob of stage 0
                      p0obex  =  p0obex,               #prob of stage 0+1
                      pRLNpos.0obex = p0obex.RLNpos    #prob of stage 1
                      
                      )


# Initial values
initsN <- function() {list(TPall3 = runif(1,0.06,0.15),
          ba1=runif(1,0.04,0.1),
          ba2=runif(1,0.25,0.5),
          Ni=runif(1, prev2016.data.n1$N-200,prev2016.data.n1$N+200),
          antMi2=prev2016.data.n1$antM)
          }
          
prev_nf1_M1 <- jags(data=prev2016.data.n1, inits=initsN, 
                    parameters.to.save=c("AP","TP","TPall", 
                    "Mn","M0","M1","nr_Stage0","nr_Stage1"), 
                    model.file="dhyper_3AgeClass_prev2016.bug",n.chain=nchains, 
                    n.iter=niter, 
                    n.burnin=nburn, DIC=TRUE, working.directory=NULL, 
                    jags.seed = 123, refresh = niter/50, progress.bar = "text", digits=5)


nf1_M1.mcmc <- as.mcmc(prev_nf1_M1)
stat_nf1dSeM1<-cbind(summary(nf1_M1.mcmc)$statistics,summary(nf1_M1.mcmc)$quantiles)
est2016<-round(stat_nf1dSeM1,5)
est2016

########################################
#Update prevalence for 2017
#########################################

#prior from 2016
stat_nf1M1<-est2016

#The prior distributions of prevalence ratios between age classes
#Alt1C from  APPENDIX S5
betapar1<-c("1.40","4.60") #mode=0.05
xrange1<-0.5
x1min<-0

mean_prev2016<-round(c(stat_nf1M1["TPall[1]","Mean"],stat_nf1M1["TPall[2]","Mean"],stat_nf1M1["TPall[3]","Mean"]),5)
sd_prev2016<-round(c(stat_nf1M1["TPall[1]","SD"],stat_nf1M1["TPall[2]","SD"],stat_nf1M1["TPall[3]","SD"]),5)


antRein<-c(427,379,1345)
sdRein<- c(29.5,15.1,38.1)

my<-c(0,0,4)

mk<-c(118,82,351)

pRLN = 0.93  

# Nordfjella sone 1
                                       
prev2017.dat.n1 <- list(
                      My = my,   #nr of animals tested cwd-positive
                      Mk = mk,     #nr of animals tested 
                      antM = antRein,    #Nr of animals per age class (estimated from pop model)
                      sd_antM = sdRein,  
                      N = 2151,      #Pop size
                      sd_N = 105.1,  #sd Pop size
                      
                      pRLN = pRLN,  #nr of samples tested with RLN tissue included
                      
                      m_Mprev=mean_prev2016,
                      sd_Mprev=sd_prev2016,
                      
                      NAstage=75,  #Nr with unknown age class
                                  
                      #beta-distribution - prevalence ratio between calves and adults 
                      Se.betaPrev1 = c(betapar1[[1]][1],betapar1[[2]][1]),
                      xrange1=xrange1,
                      x1min=x1min,
                     
                      #beta-distribution - dSe for for detectable stages
                      Se.betaRLN1 = Se.betaRLN1[1:2],
                      Se.betaRLN2= Se.betaRLN2[1:2],
                      Se.betaRLN3 = Se.betaRLN3[1:2],
                      Se.betaO2= Se.betaO2[1:2],
                      Se.betaO3 = Se.betaO3[1:2],
                      
                      #Probability of infected animals being in early disease stages  
                      p0obexRLN  =  p0obexRLN,        #prob of stage 0
                      p0obex  =  p0obex,              #prob of stage 0+1
                      pRLNpos.0obex = p0obex.RLNpos   #prob of stage 1
                      
                      )


pdat<-prev2017.dat.n1

# Initial values

initsN17p <- function() {list(
          Ni=runif(1,pdat$N-200,pdat$N+200),
          antMi2=pdat$antM,
          TPall1=runif(3,pdat$m_Mprev-pdat$sd_Mprev,pdat$m_Mprev+pdat$sd_Mprev) 
          )         
}          
          


prev2017_nf1_M2 <- jags(data=prev2017.dat.n1, inits=initsN17p, 
                    parameters.to.save=c(
                    "TPall","mean_tpall","mean_tpall_YAd",
                    "TP","mean_tp",
                    "AP","mean_ap",
                    "Prev_stage1","Prev_obexPos", 
                    "Mn","M0","M1",
                    "nr_Stage1" ,"nr_obexPos"),                   
                    model.file="dhyper_3AgeClass_2016prior_prev2017.bug",
                    n.chain=nchains, n.iter=niter, 
                    n.burnin=nburn, DIC=TRUE, working.directory=NULL, 
                    jags.seed = 123, refresh = niter/50, progress.bar = "text", 
                    digits=8)


nf1_M2.mcmc <- as.mcmc(prev2017_nf1_M2)
stat_nf1M2<-cbind(summary(nf1_M2.mcmc)$statistics,summary(nf1_M2.mcmc)$quantiles)
est2017<-round(stat_nf1M2,5)
est2017


priorP<-matrix(NA,3,ncol(est2017),dimnames=list(c("priorTPall1","priorTPall2","priorTPall3"),colnames(est2017)))
priorP[1:3,1]<-mean_prev2016
priorP[1:3,2]<-sd_prev2016

est2017<-rbind(est2017,priorP)
est2017




