#This code illustrates how the analysis was run for the smooth coated otters

#Install required packages
#choose a mirror
chooseCRANmirror()
#Install
install.packages("survival", dependencies=T)
install.packages("combinat", dependencies=T)

#Load in the code for running NBDA available at https://lalandlab.st-andrews.ac.uk/freeware/

#Input the association matrix for smooth coated otters
amg1<-matrix(data=c(0.0,	0.48,	0.37,	0.11, 0.32, 0.36, 0.36, 0.48, 0.00,	0.38,	0.07,	0.42,	0.36,	0.35, 0.37, 0.38, 0.00, 0.12, 0.41, 0.35, 0.32, 0.11, 0.07, 0.12, 0.00, 0.12, 0.09, 0.08, 0.32, 0.42, 0.41, 0.12, 0.00, 0.48, 0.48, 0.36, 0.36, 0.35, 0.09, 0.48, 0.00, 0.73, 0.36, 0.35, 0.32, 0.08, 0.48, 0.73, 0.00), nrow=7)

#Input the order of acquisition for the 6 tasks
oag1t1<-c(6,7,2,3,5,1,4)
oag1t2<-c(7,6,5,3)
oag1t3<-c(3,7,5,6,4)
oag1t4<-c(7,6,5,3,4)
oag1t5<-c(5,4,1,6,3,7)
oag1t6<-c(5,3,7,4,6,1)

#Input the individual-level variables for each otter
sex<-c(1,0,1,1,1,0,0)
age<-c(7,7,2,2,2,1,1)
parent<-c(1,1,0,0,0,0,0)
asoc<-cbind(sex,age,parent)

#Create the oadata objects
ng1t1 <-oaData(assMatrix=amg1, asoc=asoc, orderAcq=oag1t1, groupid="1", taskid="1",id=c("Mum","Dad","Sis1","Sis2","Sis3","Bro1","Bro2"))
ng1t2 <-oaData(assMatrix=amg1, asoc=asoc, orderAcq=oag1t2, groupid="1", taskid="2",id=c("Mum","Dad","Sis1","Sis2","Sis3","Bro1","Bro2"))
ng1t3 <-oaData(assMatrix=amg1, asoc=asoc, orderAcq=oag1t3, groupid="1", taskid="3",id=c("Mum","Dad","Sis1","Sis2","Sis3","Bro1","Bro2"))
ng1t4 <-oaData(assMatrix=amg1, asoc=asoc, orderAcq=oag1t4, groupid="1", taskid="4",id=c("Mum","Dad","Sis1","Sis2","Sis3","Bro1","Bro2"))
ng1t5 <-oaData(assMatrix=amg1, asoc=asoc, orderAcq=oag1t5, groupid="1", taskid="5",id=c("Mum","Dad","Sis1","Sis2","Sis3","Bro1","Bro2"))
ng1t6 <-oaData(assMatrix=amg1, asoc=asoc, orderAcq=oag1t6, groupid="1", taskid="6",id=c("Mum","Dad","Sis1","Sis2","Sis3","Bro1","Bro2"))

#Set up a matrix determining what social learning models will be considered
#Here we consider models with social learning equal in all tasks and different across all tasks
sParamMatrix<-rbind(rep(1,6),1:6)

#Run the AIC table function which fits all the combinations of models to be considered.
#i.e. all combinations of individual level variables, asocial, additive and multiplicative models, and for both rows of the sParamMatrix given above
aicTable1<-aicTable(data=c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6"), asocialVar=1:3, task=F, group=F, sParamMatrix= sParamMatrix, aic="aicc", pure=F)
aicTable1<-aicTable1[order(as.numeric(aicTable1[,7])),]
aicTable1[,8]<-as.numeric(aicTable1[,7])-as.numeric(aicTable1[1,7])
aicTable1

#Refit the best model
bestModel<-multiCoxFit(c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6"), formula=~.+parent)
summary(bestModel)

#Do random effects make much difference?
bestModelRE<-multiCoxFit(c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6"), formula=~.+parent+frailty(id))
summary(bestModelRE)
bestModelRE@coef

#Fit a model with differences in social learning rate (s) among tasks
diffSLModel<-multiCoxFit(c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6"),sParam=1:6, formula=~.+parent)
summary(diffSLModel)

#To test for an overall difference in s among tasks- take 2* the difference in logLik for with Social Transmission 
2*(37.456-35.246)

#Then look this up in a Chi-sq null distribution
pchisq(4.42,5,lower.tail=F)

#So there was little evidence of a differnce in social learning bewteen tasks (LRT: Chi-sq = 4.42; d.f.= 5; p = 0.490)

#A function to plot the associations for Fig4
plotAssociations<-function(data,lty=1, symbol=NULL, xlab="Acquisition event", ylab="Total connection to informed individuals",title=NULL,plotID=T,offset=c(0.1,0),xlim=NULL, ylim=NULL, titlePos=c(0,0)){
if(class(data)=="oaData"){oadata<-data}
if(class(data)=="taData"){oadata<-data@oadata}

if(is.null(xlim)){xlim<-c(0,max(oadata@coxdata$time2)+0.5)}

	if(is.null(symbol)){
		plot(oadata@coxdata$time2,oadata@coxdata$stMetric,col=oadata@coxdata$status+1,xlab=xlab, ylab=ylab,main="",xlim=xlim, ylim=ylim);
		points(oadata@coxdata$time2[oadata@coxdata$status==1],oadata@coxdata$stMetric[oadata@coxdata$status==1],col=2);
		lines(oadata@coxdata$time2[oadata@coxdata$status==1],oadata@coxdata$stMetric[oadata@coxdata$status==1],col=2, lty=lty);		
	}else{
		plot(oadata@coxdata$time2,oadata@coxdata$stMetric,col=oadata@coxdata$status+1,pch=as.numeric(as.factor(oadata@coxdata[,6+symbol])),xlab=xlab, ylab=ylab, main="",xlim=xlim, ylim=ylim);
		points(oadata@coxdata$time2[oadata@coxdata$status==1],oadata@coxdata$stMetric[oadata@coxdata$status==1],col=2,pch=as.numeric(as.factor(oadata@coxdata[,6+symbol]))[oadata@coxdata$status==1]);
		lines(oadata@coxdata$time2[oadata@coxdata$status==1],oadata@coxdata$stMetric[oadata@coxdata$status==1],col=2, lty=lty);
	}
	if(plotID){
			text(oadata@coxdata$time2[oadata@coxdata$status==1]+offset[1],oadata@coxdata$stMetric[oadata@coxdata$status==1]+offset[2],col=2,labels=oadata@coxdata$id[oadata@coxdata$status==1]);		
	}
	text(x=titlePos[1],y=titlePos[2],labels=title)
}


par(mfrow=c(3,2), mar=c(4,4,0.5,0.5))
plotAssociations(ng1t1,symbol=3,title="a) Task 1",titlePos=c(1.35,1.9), offset=c(0.4,0),xlab="",ylab="",xlim=c(1,7.5),ylim=c(0,2))
plotAssociations(ng1t2,symbol=3,title="b) Task 2",titlePos=c(1.35,1.9), offset=c(0.4,0),xlab="",ylab="",xlim=c(1,7.5),ylim=c(0,2))
plotAssociations(ng1t3,symbol=3,title="c) Task 3",titlePos=c(1.35,1.9), offset=c(0.4,0),xlab="",xlim=c(1,7.5),ylim=c(0,2))
plotAssociations(ng1t4,symbol=3,title="d) Task 4",titlePos=c(1.35,1.9), offset=c(0.4,0),xlab="",ylab="",xlim=c(1,7.5),ylim=c(0,2))
plotAssociations(ng1t5,symbol=3,title="e) Task 5",titlePos=c(1.35,1.9), offset=c(0.4,0),ylab="",xlim=c(1,7.5),ylim=c(0,2))
plotAssociations(ng1t6,symbol=3,title="f) Task 6",titlePos=c(1.35,1.9), offset=c(0.4,0),ylab="",xlim=c(1,7.5),ylim=c(0,2))

#Now obtain confidence intervals using profile likelihood technique
oadata<-combineOaCoxData(c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6"))

sVals<-logLik<-seq(0,10,0.1)

for(i in 1:length(sVals)){
	logLik[i]<-multiCoxLikelihood(sVals[i],oadata, formula=~.+parent)	
}

plot(sVals, logLik, type="l")
abline(h= 37.456+1.92, lty=2)

sVals<-logLik<-seq(0.4,0.5,0.001)
#0.44 lower limit

#we can see that no matter how high we push s it is still within the 95% CI
multiCoxLikelihood(9999999,oadata, formula=~.+parent)
37.456+1.92
[1] 37.46922
> 37.456+1.92
[1] 39.376
> 


#Now we estimate the number of events that occured by social transmission (excliuding the innovator)

	object<-oadata
	object@mldata<-object@coxdata[object@coxdata$status==1,]
	
		number<-sum((37.428352 *object@mldata$stMetric)/(1+ 37.428352 *object@mldata$stMetric))
		total<-dim(object@mldata)[1]-6

number/total

#96.0% social transmission excluding innovator

#lower limit of 95%CI [0.4-Inf]

	object<-oadata
	object@mldata<-object@coxdata[object@coxdata$status==1,]
	
		number<-sum((0.44 *object@mldata$stMetric)/(1+ 0.44 *object@mldata$stMetric))
		total<-dim(object@mldata)[1]-6

number/total

#25.7% social transmission


#Read in AIC table to get SE's and MLE's for model averaging across models with equal s across tasks 
aicTable<-read.csv("AICsmoothInput.csv", header=T)

SE<-MLE<-matrix(NA,nrow=dim(aicTable)[1],ncol=4)

for(i in 1:dim(aicTable)[1]){

	ilv<-c((!is.na(aicTable[i,2]))*1,(!is.na(aicTable[i,3]))*2,(!is.na(aicTable[i,4]))*3)
	ilv<-ilv[ilv>0]

	if(aicTable[i,6]=="social"){
		if(aicTable[i,1]|is.na(aicTable[i,1]=="NA")){
			model<-optim(par=rep(0,length(ilv)+1),addLikelihood,hessian=T, data=c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6"), bounded=F, asocialVar=ilv)
			se<-sqrt(diag(solve(model$hessian)))
			MLE[i,1]<-model$par[1]
			SE[i,1]<-se[1]
			MLE[i,ilv+1]<-model$par[-1]
			SE[i,ilv+1]<-se[-1]
		}
	
	}else{
			model<-optim(par=rep(0,length(ilv)),nulladdLikelihood,hessian=T, data=c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6"), bounded=F, asocialVar=ilv)
			se<-sqrt(diag(solve(model$hessian)))
			MLE[i,ilv+1]<-model$par
			SE[i,ilv+1]<-se				
	}
}

#Fill in Multiplicative models

i<-1
	ilv<-c((!is.na(aicTable[i,2]))*1,(!is.na(aicTable[i,3]))*2,(!is.na(aicTable[i,4]))*3)
	ilv<-ilv[ilv>0]

model<-optim(par=0,multiCoxLikelihood,hessian=T, oadata=combineOaCoxData(c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6")), bounded=F, formula=~.+parent)
se<-sqrt(diag(solve(model$hessian)))
MLE[i,1]<-model$par[1]
SE[i,1]<-se[1]

model<-multiCoxFit(c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6"), bounded=F, formula=~.+parent)
se<-as.vector(model@coef[,3])
mle<-as.vector(model@coef[,1])
MLE[i,ilv+1]<-mle
SE[i,ilv+1]<-se	


i<-2
	ilv<-c((!is.na(aicTable[i,2]))*1,(!is.na(aicTable[i,3]))*2,(!is.na(aicTable[i,4]))*3)
	ilv<-ilv[ilv>0]

model<-optim(par=0,multiCoxLikelihood,hessian=T, oadata=combineOaCoxData(c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6")), bounded=F, formula=~.+sex+parent)
se<-sqrt(diag(solve(model$hessian)))
MLE[i,1]<-model$par[1]
SE[i,1]<-se[1]

model<-multiCoxFit(c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6"), bounded=F, formula=~.+sex)
se<-as.vector(model@coef[,3])
mle<-as.vector(model@coef[,1])
MLE[i,ilv+1]<-mle
SE[i,ilv+1]<-se	

i<-10
	ilv<-c((!is.na(aicTable[i,2]))*1,(!is.na(aicTable[i,3]))*2,(!is.na(aicTable[i,4]))*3)
	ilv<-ilv[ilv>0]

model<-optim(par=0,multiCoxLikelihood,hessian=T, oadata=combineOaCoxData(c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6")), bounded=F, formula=~.+sex)
se<-sqrt(diag(solve(model$hessian)))
MLE[i,1]<-model$par[1]
SE[i,1]<-se[1]

model<-multiCoxFit(c("ng1t1", "ng1t2", "ng1t3", "ng1t4", "ng1t5", "ng1t6"), bounded=F, formula=~.+sex)
se<-as.vector(model@coef[,3])
mle<-as.vector(model@coef[,1])
MLE[i,ilv+1]<-mle
SE[i,ilv+1]<-se	

write.csv(MLE,"MLEsmoothCut.csv")
write.csv(SE,"SEsmoothCut.csv")
