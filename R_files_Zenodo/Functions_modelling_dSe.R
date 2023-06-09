#Functions_modelling_dSe.R 
#Here is the functions required to run the two R-scripts Modelling_dSe.R 
#and Application2_PrDetectingDisease.R

#Draw n random values from pert-distribution defined by min, max and mode
rpert <- function( n, x.min, x.max, x.mode, lambda = 4 ){
    x.range <- x.max - x.min
    if( x.range == 0 ) return( rep( x.min, n ))
    mu <- ( x.min + x.max + lambda * x.mode ) / ( lambda + 2 )
    # special case if mu = mode
    if( mu == x.mode ){
        v <- ( lambda / 2 ) + 1
    }
    else {
        v <- (( mu - x.min ) * ( 2 * x.mode - x.min - x.max )) /
            (( x.mode - mu ) * ( x.max - x.min ))
    }
    w <- ( v * ( x.max - mu )) / ( mu - x.min )
    return ( rbeta( n, v, w ) * x.range + x.min )
}


################################################################################
#Infection curves (scores) for obex and RLN
################################################################################
#Make the curve for 4 yr
x <- seq(0,(48+24*2),0.25)  #48 = 2 yr (time-step unit 1 = 1/2 month)

#Obex score  ####
b = -5.97
m = 0.271 # slope
yobex <- exp((b + m*x)) / (1 + exp((b + m*x)))

#RLN score ######
bL= -4
mL = 0.5 # slope
yL <- 0.7*exp((bL + mL*x)) / (1 + exp((bL + mL*x)))


################################################################################
#Function to estimate Sensitivity by rescaling curves  
#See Table 1 for input values.
#The result is a table (dSeTab) with number of months since infection, running time step and the corresponding expected Se for samples including RLN, and for samples with low and high quality samples of obex
################################################################################

#Analytical test sensitivity
#testSe_RLNobex=0.95
#testSe_obex=0.925

#x.time = time step vector
xmax.2yr<-length(seq(0,48,0.25))

Make_dSeTable <- function(yL=yL,x1=x.time,yobex=yobex,Selow=0.25,HighQ_Selow=0.5,
     LowQ_Selow=0.25,testSe_RLNobex=0.95,testSe_obex=0.925){
     
     #Estimate dSe for RLN (including obex) by rescaling the infection curve 
     #yL<0.25 -> 0
     #Values of yL>= 0.25 are rescaled between Selow and testSev

     #Selow=0.25
     s_yL<-yL
     testSev=testSe_RLNobex
     
     #Scale according to value at 2 yr
     s_yL[yL>=0.25] <-  ((testSev-Selow)*(s_yL[yL>=0.25]-min(s_yL[yL>=0.25]))
      /(s_yL[xmax.2yr]-min(s_yL[yL>=0.25])))+Selow
     s_yL[yL<0.25]<-0


     #Estimate dSe for obex by rescaling the infection curve 
     #yobex<0.25 -> 0
     #Values of yobex >= 0.25 are rescaled between Selow and testSevO
     #Selow differs between low and high quality samples 

     #HighQ_Selow=0.5
     #LowQ_Selow=0.25
     testSevO=testSe_obex
     
     s_yOH=s_yOL=yobex
     s_yOH[yobex>=0.25] <-  ((testSevO-HighQ_Selow)*(s_yOH[yobex>=0.25]-min(s_yOH[yobex>=0.25]))/((s_yOH[xmax.2yr])-min(s_yOH[yobex>=0.25])))+HighQ_Selow
     s_yOH[yobex<0.25]<-0
     
     s_yOH[yobex>testSevO]<-testSevO
     
     s_yOL[yobex>=0.25] <-  ((testSevO-LowQ_Selow)*(s_yOL[yobex>=0.25]-min(s_yOL[yobex>=0.25]))/((s_yOL[xmax.2yr])-min(s_yOL[yobex>=0.25])))+LowQ_Selow
     s_yOL[yobex<0.25]<-0

     s_yOL[yobex>testSevO]<-testSevO

     #Running time step
     timestep <- seq(0,length(yL)-1,1)      #2 years = 48 months, 4 time step of each month
     
     #Sensitivity function for discrete time points
     yRLN=s_yL
     yObexLowq=s_yOL
     yObexHighq=s_yOH

     mpi<-x1/2   #Number of months since exposure (months post infection) 
      
     dSeTab<-cbind(mpi,timestep,yRLN,yObexLowq,yObexHighq)
     dSeTab
     }

     
###############################################################################
#Estimate Sensitivity by rescaling infection curves
###############################################################################


x <- seq(0,(48+24*2),0.25)  #4 yr (time-step unit 1 = 1/2 month)
x.time=x

dSeTab<-Make_dSeTable(yL=yL,yobex=yobex,x1=x.time)


#Remove time step=0
dSeTab1<-dSeTab[-1,]

################################################################################
#Function to add stochasticity to SE scores  (see Table 1 for input values)
################################################################################

IncludeSeVar<-function(Y=Y1){
      scoreY=Y
      var_low=1.25
      var_up=0.6
      
for(i in 1:length(Y)){
      scoreY=Y
      if(Y[i]>0){
      scoreY[i]<-rpert(1,x.min=(Y[i]^var_low),x.max=(Y[i]^var_up),x.mode=Y[i])
      }
      }
      scoreY[scoreY>1]<-1
      return(scoreY)
      }


###############################################################################
#Simulate distribution of dSe-values
#1) Assuming time since infection are randomly drawn from time line of 2 years# 
###############################################################################
#dSeTab1 = table with expected Se (estimated by rescaling infection curves) for testing RLN and high and low-quality obex samples  
#nsim = nr of simulations
#x1.max = max time step (time since infection) 

#rpert-function are utilized to random draw values from a distribution defined by min, max and mode
#Input data from Table 1
#PrLQ.min=0.02,
#PrLQ.max=0.60,
#PrLQ.mode=0.22,


simulate_dSe <- function(nsim=nsim,dSeT=dSeTab1,x1.max=x1.max,PrLQ.min=PrLQ.min,PrLQ.max=PrLQ.max,PrLQ.mode=PrLQ.mode){
     
     #Draw random time since infection 
     exptime <- sample(seq(1,x1.max,1),size=nsim,replace=T)

     score_RLNobex <- dSeT[exptime,"yRLN"]
     score_obexHQ <- dSeT[exptime,"yObexHighq"]
     score_obexLQ <- dSeT[exptime,"yObexLowq"]

     #Add stochasticity        
     SE_RLNobex<-IncludeSeVar(Y=score_RLNobex)
     SE_obexLQ<-IncludeSeVar(Y=score_obexLQ)
     SE_obexHQ<-IncludeSeVar(Y=score_obexHQ)

     #Draw random proportion of LQ obex sample from given pert-distribution
     PrLQv<-rpert(n=nsim,x.min=PrLQ.min,x.max=PrLQ.max,x.mode=PrLQ.mode)
    
     #Adjust for proportion of LQ obex samples
     SE_obex<-PrLQv* SE_obexLQ+(1-PrLQv)* SE_obexHQ

     dSE <- cbind(exptime,SE_RLNobex,SE_obex,SE_obexLQ,SE_obexHQ,PrLQv)
     colnames(dSE) = c("exptime","dSE_RLN","dSE_obex","dSE_obexLQ","dSE_obexHQ","PropLQ")

     dSE
     
     }


x1 <- seq(0,(48),0.25)
x1A<-length(x1)-1

x1 <- seq(0,(34),0.25)
x1Y<-length(x1)-1

x1 <- seq(0,(10),0.25)
x1C<-length(x1)-1


#Input data from Table 1
PrLQ.min =  0.02
PrLQ.max =  0.60
PrLQ.mode = 0.22

nsim=10000

SeAdults<-simulate_dSe(nsim=nsim,dSeT=dSeTab1,x1.max=x1A,PrLQ.min,PrLQ.max,PrLQ.mode)

SeYearlings<-simulate_dSe(nsim=nsim,dSeT=dSeTab1,x1.max=x1Y,PrLQ.min,PrLQ.max,PrLQ.mode)

SeCalves<-simulate_dSe(nsim=nsim,dSeT=dSeTab1,x1.max=x1C,PrLQ.min,PrLQ.max,PrLQ.mode)

################################################################################
#Simulate distribution of dSe-values 
#2) Assuming time since infection are randomly drawn from time line of 3 years
################################################################################

#Time step for drawing random time since infection, assuming lenght of disease progression is 3 years
#48 = 3 year = 36 months =   72 half months

#(48/36)*24
#32 = 2 year = 24 months =   48 half months

#(48/36)*17
#22.7 =  17 months =   34 half months
#22.5 ~ 22.7

#(48/36)*5
#6.7 =  5 months =   10 half months
#6.75 ~ 6.7 

#Adults. Assume time since infection for adults are max 48, i.e. 3 years back in time
#Yearlings: 22.5 = 1 yr and 5 months at time of hunting
#Calves: 6.75 = 5 months at time of hunting

x1 <- seq(0,(48),0.25)
x1A<-length(x1)-1

x1 <- seq(0,(22.5),0.25)
x1Y_3yrI<-length(x1)-1

x1 <- seq(0,(6.75),0.25)
x1C_3yrI<-length(x1)-1

nsim=10000

Se3yrI_Adults<-simulate_dSe(nsim=nsim,dSeT=dSeTab1,x1.max=x1A,PrLQ.min,PrLQ.max,PrLQ.mode)

Se3yrI_Yearlings<-simulate_dSe(nsim=nsim,dSeT=dSeTab1,x1.max=x1Y_3yrI,PrLQ.min,PrLQ.max,PrLQ.mode)

Se3yrI_Calves<-simulate_dSe(nsim=nsim,dSeT=dSeTab1,x1.max=x1C_3yrI,PrLQ.min,PrLQ.max,PrLQ.mode)


################################################################################
#Simulate distribution of dSe-values 
#3) Assuming an exponential time decay distribution of time since infection 
################################################################################

simulateExpDecay_dSe <- function(nsim=nsim,dSeT=dSeTab1,x1.max=x1.max,meanInfected,PrLQ.min=PrLQ.min,PrLQ.max=PrLQ.max,PrLQ.mode=PrLQ.mode){
     
     #Draw random time since infection 
     #Make sure to draw sufficiently many (n*100), to compensate for exptime>max(x1) being excluded
     #mean time since infection defined by meanInfected (half months) * 4 time steps
     exptime<-round(rexp(nsim*100,rate=1/(meanInfected*4))) 
     #exptime<-exptime+1   #Trenger vel ikke denne når vi fjerner timestep=0? 
     exptime1<-exptime[exptime<=x1.max]
     exptime1<-exptime1[exptime1>0]
     exptime1<-exptime1[1:nsim]  #keep n simulated values

     score_RLNobex <- dSeT[exptime1,"yRLN"]
     score_obexHQ <- dSeT[exptime1,"yObexHighq"]
     score_obexLQ <- dSeT[exptime1,"yObexLowq"]

     #Add stochasticity        
     SE_RLNobex<-IncludeSeVar(Y=score_RLNobex)
     SE_obexLQ<-IncludeSeVar(Y=score_obexLQ)
     SE_obexHQ<-IncludeSeVar(Y=score_obexHQ)

     #Draw random proportion of LQ obex sample from given pert-distribution
     PrLQv<-rpert(n=nsim,x.min=PrLQ.min,x.max=PrLQ.max,x.mode=PrLQ.mode)
 
     #Adjust for proportion of LQ obex samples
     SE_obex<-PrLQv* SE_obexLQ+(1-PrLQv)* SE_obexHQ

     dSE <- cbind(exptime1,SE_RLNobex,SE_obex,SE_obexLQ,SE_obexHQ,PrLQv)
     colnames(dSE) = c("exptime","dSE_RLN","dSE_obex","dSE_obexLQ","dSE_obexHQ","PropLQ")

     dSE
     }


x1 <- seq(0,(48+24*2),0.25)
x1A<-length(x1)-1

x1 <- seq(0,(34),0.25)
x1Y<-length(x1)-1

x1 <- seq(0,(10),0.25)
x1C<-length(x1)-1


#Input data from Table 1
PrLQ.min =  0.02
PrLQ.max =  0.60
PrLQ.mode = 0.22


nsim=10000

SeExpAdults<-simulateExpDecay_dSe(nsim=nsim,dSeT=dSeTab1,x1.max=x1A,meanInfected=24,PrLQ.min,PrLQ.max,PrLQ.mode)

SeExpYearlings<-simulateExpDecay_dSe(nsim=nsim,dSeT=dSeTab1,x1.max=x1Y,meanInfected=24,PrLQ.min,PrLQ.max,PrLQ.mode)

SeExpCalves<-simulateExpDecay_dSe(nsim=nsim,dSeT=dSeTab1,x1.max=x1C,meanInfected=24,PrLQ.min,PrLQ.max,PrLQ.mode)


##############################################################################


