#set working directory
#setwd("E:/Eggeman - migration/Eggeman_Data and Rscripts")
install.packages("R2jags")


setwd("/Users/mark.hebblewhite/Dropbox/Ya Ha Tinda/oikos_revision/Eackercode/BayesCalfCode")
ls()
#read in data file
calf_surv<-read.csv("calf_surv.csv",header=T)
head(calf_surv)

##############################################################################################################
#variables considered
#Migrant	status	(migrant	1	resident	=0)	Year	(Early=1	Late	=	0)

require(R2jags)
load.module("glm")

#manipulate data to match structure used for MARK (i.e., an offset for the survival period (summer vs. winter))
newdata<-cbind((c(calf_surv$summer, calf_surv$winter)),c(rep("summer",107),rep("winter",107)),rep(calf_surv$migcode,2),rep(calf_surv$period,2))
newdata2<-as.data.frame(newdata)
names(newdata2)<-c("surv","period","migcode","Year")
newdata2$period2<-ifelse(newdata2$period=="summer",0,1)
N<-length(newdata2$surv)

# create list of data
data<-list(Y=as.numeric(as.vector(newdata2$surv)),N=N,period2=as.numeric(as.vector(newdata2$period2))
,migcode=as.numeric(as.vector(newdata2$migcode)),Year=as.numeric(as.vector(newdata2$Year)))

# intial values for MCMC
inits<-function(){

		beta0guess = rnorm(1,0,2)
		beta1guess = rnorm(1,0,2)
		beta2guess = rnorm(1,0,2)
		beta3guess = rnorm(1,0,2)

		list(
		
		beta0=beta0guess,
		beta1=beta1guess,
		beta2=beta2guess,
		beta3=beta3guess

)
}

# parameters to track	
params<-c("beta0","beta1","beta2","beta3","Surv.mig.earlyS","Surv.mig.lateS","Surv.res.earlyS","Surv.res.lateS",
"Surv.mig.earlyW","Surv.mig.lateW","Surv.res.earlyW","Surv.res.lateW")

# run jags function for GIBBS sampling
rr.log2 <- jags(data, inits, params, "logisticJAGS.txt", 
        n.chains=2, n.iter=50000, n.burnin=30000, n.thin=2)

#view MCMC output
rr.log2



##################################################################################################################

