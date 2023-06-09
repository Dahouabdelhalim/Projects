#################################################################################################
# Code for Martinez et al Oikos. Version December 20 2020.
#################################################################################################
# Functions
# Function to compute the BIC
bic	<-	function(loglik,par,n){
	
	log(n)*par-2*(loglik)
	
}

# Error bar function for Figures
error.bar <- function(x,y,upper,lower=upper,length=0.1,type=c("X","Y"),...){if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	stop("vectors must be same length")
	if(type=="Y"){
		arrows(x,upper, x,lower, angle=90, code=3, length=length, ...)
		}
	else if(type=="X"){
		arrows(upper,y,lower,y, angle=90, code=3, length=length, ...)
	}
	}

# Load required packages
library(survival)

#################################################################################################
# PART 1
# Analyses of the response of tamarin monkeys to heterospecific,
# intraspecific and control playbacks. Here we use a right censored survival analysis and compared three 
# different statistical hypothesis.
# fit1 (H1 in manuscript)	=	No difference among stimulus
# fit2 (H3 in manuscript)	=	Differences between Control and Alarm stimulus
# fit3 (H2 in manuscript)	=	Differences between Control and alarm and betwenn intra 
#								and heterospecific response
#################################################################################################

# Read the data
trials	<-	read.delim("SAFU_response.txt")
not.flee.SAFU	<-	which(trials$Time2flee==30)
status.SAFU	<-	rep(1,length(trials$Time2flee))
status.SAFU[not.flee.SAFU]<-0


# Fit the models
fit1	<-	survreg(Surv(Time2flee,status.SAFU)~1,data=trials,dist="exponential")
fit2	<-	survreg(Surv(Time2flee,status.SAFU)~Treatment-1,data=trials,dist="exponential")
fit3	<-	survreg(Surv(Time2flee,status.SAFU)~Treatment1-1,data=trials,dist="exponential")

# Select best model according to BIC
fit.list.monkeys		<-	list(fit1,fit2,fit3)

fit.loglik.monkeys	<-	unlist(lapply(fit.list.monkeys,function(x)x$loglik[2]))
fit.par.monkeys		<-	unlist(lapply(fit.list.monkeys,function(x)x$df))

fit.bic.monkeys		<-	bic(fit.loglik.monkeys,fit.par.monkeys,nrow(trials))
fit.dbic.monkeys		<-	fit.bic.monkeys-min(fit.bic.monkeys)

best.mod.monkeys		<-	fit.list.monkeys[[which(fit.dbic.monkeys==0)]]

# Predict probability of response at the end of the trial
prob.resp.monkeys		<-	survfit(Surv(Time2flee,status.SAFU)~Treatment1-1,data=trials)
summary(prob.resp.monkeys)

# For parametric bootstrap
nboot				<-	2000
exp.rates.fit1		<-	exp(coef(fit1))
exp.rates.fit2		<-	exp(coef(fit2))
exp.rates.fit3		<-	exp(coef(fit3))
boot.mat.fit1		<- 	matrix(NA,nrow=nboot,ncol=1,dimnames=list(1:nboot,"Control"))
boot.mat.fit2		<- 	matrix(NA,nrow=nboot,ncol=3,dimnames=list(1:nboot,c("Control","ST","THSY")))
boot.mat.fit3		<- 	matrix(NA,nrow=nboot,ncol=2,dimnames=list(1:nboot,c("Alarm","Control")))
fit3.boot.presponse	<-	numeric(nboot)

for(i in 1:nboot){
	
	# for fit 1
	time2flee.fit1						<-	rexp(21,rate=1/exp.rates.fit1)
	time2flee.fit1[time2flee.fit1>30]	<-	30
	fit1.status							<-	rep(1,21)
	fit1.status[time2flee.fit1==30]		<-	0	
	ith.fit1								<-	survreg(Surv(time2flee.fit1,fit1.status)~1)
	boot.mat.fit1[i,]					<-	exp(coef(ith.fit1))
	
	# for fit2
	time2flee.fit2						<-	rexp(21,rate=1/rep(exp.rates.fit2,each=7))
	time2flee.fit2[time2flee.fit2>30]	<-	30
	treatment							<-	rep(c("Control","ST","THSY"),each=7)	
	fit2.status							<-	rep(1,21)
	fit2.status[time2flee.fit2==30]		<-	0	
	ith.fit2 							<-	survreg(Surv(time2flee.fit2,fit2.status)~treatment-1)
	boot.mat.fit2[i,]					<-	exp(coef(ith.fit2))

	# for fit3
	fit3.rates							<-	c(rep(exp.rates.fit3[1],14),rep(exp.rates.fit3[2],7))			
	time2flee.fit3						<-	rexp(21,rate=1/fit3.rates)
	time2flee.fit3[time2flee.fit3>30]	<-	30
	treatment							<-	rep(c("Alarm","Control"),times=c(14,7))	
	fit3.status							<-	rep(1,21)
	fit3.status[time2flee.fit3==30]		<-	0	
	ith.fit3 							<-	survreg(Surv(time2flee.fit3,fit3.status)~treatment-1)
	boot.mat.fit3[i,]					<-	exp(coef(ith.fit3))
	fit3.survfit							<-	survfit(Surv(time2flee.fit3,fit3.status)~treatment-1)$surv
	fit3.boot.presponse[i]				<-	fit3.survfit	[length(fit3.survfit)-1]
}

time2flee.confint.fit1	<-	apply(boot.mat.fit1,2,quantile,prob=c(0.025,0.975))
time2flee.confint.fit2	<-	apply(boot.mat.fit2,2,quantile,prob=c(0.025,0.975))
time2flee.confint.fit3	<-	apply(boot.mat.fit3,2,quantile,prob=c(0.025,0.975))
fit3.presponse.confint	<-	quantile(fit3.boot.presponse,prob=c(0.025,0.975))

# Goodness of fit
gsq.monkeys		<-	-2*(c(fit.loglik.monkeys[1],fit.loglik.monkeys[3])-fit.loglik.monkeys[2])
p.val.monkeys	<-	1-pchisq(gsq.monkeys,c(2,1))

#################################################################################################
# PART 2:
# Modeling the movement of tamarins as a response to the alarm calls.
# ST behavioral response
#################################################################################################
# Data analysis
distance	<-	trials$Delta.distance+.Machine$double.xmin # Add a very small number such that the rate is not inf

all.mean		<-	mean(trials$Delta.distance)
move.mean	<-	by(distance,trials$Treatment,mean)
move.alarm	<-	mean(distance[trials$Treatment!="Control"])


loglik.fullmod	<-	sum(dexp(distance[trials$Treatment=="Control"],rate=1/move.mean[1],log=TRUE)) +
					sum(dexp(distance[trials$Treatment=="ST"],rate=1/move.mean[2],log=TRUE)) +
					sum(dexp(distance[trials$Treatment=="THSY"],rate=1/move.mean[3],log=TRUE))

bic(loglik.fullmod,3,21)
					
loglik.mod		<-	sum(dexp(distance[trials$Treatment1=="Control"],rate=1/move.mean[1],log=TRUE)) +
					sum(dexp(distance[trials$Treatment1=="Alarm"],rate=1/move.alarm,log=TRUE))
bic(loglik.mod,2,21)
loglik.null		<-	sum(dexp(distance,rate=1/all.mean,log=TRUE))
bic(loglik.null,1,21)

bic.mov	<-	c(bic(loglik.null,1,21),bic(loglik.fullmod,3,21),bic(loglik.mod,1,21))
dbic.mov<-	bic.mov-min(bic.mov)


# Parametric bootstrap
nboot		<-	2000
rates	<-	move.mean
boot.mat.null	<-	matrix(NA,nrow=nboot,ncol=1)
colnames(boot.mat.null)	<-	c("Trials")
boot.mat.mod		<-	matrix(NA,nrow=nboot,ncol=2)	
colnames(boot.mat.mod)	<-	c("Control","Alarm")
boot.mat.fullmod	<-	matrix(NA,nrow=nboot,ncol=3)
colnames(boot.mat.fullmod)	<-	c("Control","ST","THSY")

for(i in 1:nboot){
	
	# Null mod
	mov.null				<-	rexp(21,1/all.mean)
	boot.mat.null[i,]	<-	mean(mov.null)
	
	# Mod
	mov.mod.rates		<-	rep(c(move.mean[1],move.alarm),times=c(7,14))
	mov.mod				<-	rexp(21,1/mov.mod.rates)
	boot.mat.mod[i,]	<-	c(mean(mov.mod[1:7]),mean(mov.mod[8:21]))
	
	# Full mod
	mov.full.rates		<-	rep(move.mean,each=7)
	mov.full				<-	rexp(21,1/mov.full.rates)
	boot.mat.fullmod[i,]<-	c(mean(mov.full[1:7]),mean(mov.full[8:14]),mean(mov.full[15:21]))
	
}

distance.confint.null	<-	apply(boot.mat.null,2,quantile,prob=c(0.025,0.975))
distance.confint.mod	<-	apply(boot.mat.mod,2,quantile,prob=c(0.025,0.975))
distance.confint.fullmod<-	apply(boot.mat.fullmod,2,quantile,prob=c(0.025,0.975))

# Goodness of fit
gsq.mov		<-	-2*(c(loglik.null,loglik.mod)-loglik.fullmod)
p.val.mov	<-	1-pchisq(gsq.mov,c(2,1))

#################################################################################################
# PART 3:
# Analyses of the response of antshrikes to heterospecific
# homospecific and control playbacks. Here we use a right censored survival analysis and compares five 
# different statistical models to test for differences in the responses. 
# Read and get the data ready
trials2	<-	read.delim("THSY_response.txt")
trials2$Time2flee[trials2$Time2flee==0]	<-	0.01
not.flee	<-	which(trials2$Time2flee==30)
status	<-	rep(1,length(trials2$Time2flee))
status[not.flee]<-0


# Fit the models
fit4					<-	survreg(Surv(Time2flee,status)~1,data=trials2,dist="exponential")
fit5					<-	survreg(Surv(Time2flee,status)~Trial-1,data=trials2,dist="exponential")
fit6					<-	survreg(Surv(Time2flee,status)~Trial1-1,data=trials2,dist="exponential")

fit.list.birds		<-	list(fit4,fit5,fit6)

fit.loglik.birds		<-	unlist(lapply(fit.list.birds,function(x)x$loglik[2]))
fit.par.birds		<-	unlist(lapply(fit.list.birds,function(x)x$df))

fit.bic.birds		<-	bic(fit.loglik.birds,fit.par.birds,24)
fit.dbic.birds		<-	fit.bic.birds-min(fit.bic.birds)

best.mod.birds	<-	fit.list.birds[[which(fit.dbic.birds==0)]]
prob.resp.birds	<-	survfit(Surv(Time2flee,status)~Trial1-1,data=trials2)
summary(prob.resp.birds)

# For parametric bootstrap
nboot				<-	2000
exp.rates.fit4		<-	exp(coef(fit4))
exp.rates.fit5		<-	exp(coef(fit5))
exp.rates.fit6		<-	exp(coef(fit6))
boot.mat.fit4		<- 	matrix(NA,nrow=nboot,ncol=1,dimnames=list(1:nboot,"Control"))
boot.mat.fit5		<- 	matrix(NA,nrow=nboot,ncol=3,dimnames=list(1:nboot,c("Control","ST","THSY")))
boot.vec.fit5		<-	numeric(nboot)
boot.mat.fit6		<- 	matrix(NA,nrow=nboot,ncol=2,dimnames=list(1:nboot,c("Alarm","Control")))
fit6.boot.presponse	<-	numeric(nboot)

for(i in 1:nboot){
	
	# for fit 4
	time2flee.fit4						<-	rexp(24,rate=1/exp.rates.fit4)
	time2flee.fit4[time2flee.fit4>30]	<-	30
	fit4.status							<-	rep(1,24)
	fit4.status[time2flee.fit4==30]		<-	0	
	ith.fit4								<-	survreg(Surv(time2flee.fit4,fit4.status)~1)
	boot.mat.fit4[i,]					<-	exp(coef(ith.fit4))
	
	# for fit 5
	time2flee.fit5						<-	rexp(24,rate=1/rep(exp.rates.fit5,each=8))
	time2flee.fit5[time2flee.fit5>30]	<-	30
	treatment							<-	rep(c("Control","ST","THSY"),each=8)	
	fit5.status							<-	rep(1,24)
	fit5.status[time2flee.fit5==30]		<-	0	
	ith.fit5 							<-	survreg(Surv(time2flee.fit5,fit5.status)~treatment-1)
	boot.mat.fit5[i,]					<-	exp(coef(ith.fit5))
	ith.prob.vec							<-	survfit(Surv(time2flee.fit5,fit5.status)~treatment-1)$surv
	boot.vec.fit5						
	
	# for fit 6
	fit6.rates							<-	c(rep(exp.rates.fit6[1],8),rep(exp.rates.fit6[2],16))			
	time2flee.fit6						<-	rexp(24,rate=1/fit6.rates)
	time2flee.fit6[time2flee.fit6>30]	<-	30
	treatment							<-	rep(c("Control","Alarm"),times=c(8,16))	
	fit6.status							<-	rep(1,24)
	fit6.status[time2flee.fit6==30]		<-	0	
	ith.fit6 							<-	survreg(Surv(time2flee.fit6,fit6.status)~treatment-1)
	boot.mat.fit6[i,]					<-	exp(coef(ith.fit6))
	fit6.survfit							<-	survfit(Surv(time2flee.fit6,fit6.status)~treatment-1)$surv
	fit6.boot.presponse[i]				<-	fit6.survfit	[length(fit6.survfit)-1]
	
			
}

time2flee.confint.fit4	<-	apply(boot.mat.fit4,2,quantile,prob=c(0.025,0.975))
time2flee.confint.fit5	<-	apply(boot.mat.fit5,2,quantile,prob=c(0.025,0.975))
time2flee.confint.fit6	<-	apply(boot.mat.fit6,2,quantile,prob=c(0.025,0.975))
fit6.presponse.confint	<-	quantile(fit6.boot.presponse,prob=c(0.025,0.975))

# Goodness of fit
gsq.birds		<-	-2*(c(fit.loglik.birds[1],fit.loglik.birds[3])-fit.loglik.birds[2])
p.val.birds		<-	1-pchisq(gsq.birds,c(2,1))


#################################################################################################
# PART 4: Figure
#################################################################################################
jpeg("Figure1.tiff",width=6.5,height=6.5,units="in"
	,res=300,type="cairo",family="times")
mat	<-	matrix(c(1,1,2,2,0,3,3,0),nrow=2,ncol=4,byrow=TRUE)
layout(mat)
par(mgp=c(1.5,0.2,0),mar=c(3,3,2,2),tcl=-0.2,oma=c(1,1,1,1))
# Thamnomanes response
plot(0,0,ylim=c(0,60),xlim=c(0.5,3.5),xlab="",ylab="Time until response (Sec)",xaxt="n",bty="l",cex.lab=1.2)
axis(1,at=1:3,labels=c("Control\\nn=8","Alarm\\nn=8","Alarm\\nn=8"),mgp=c(3,1.2,0),cex.axis=1)
points(rep(1,8)-0.075,trials2$Time2flee[trials2$Trial=="Control"],pch=21,col="gray90",bg="gray90",cex=3)
points(rep(2,3)-0.075,trials2$Time2flee[trials2$Trial=="SAFU"&trials2$Time2flee<30],pch=21,col="gray90",bg="gray90",cex=3*(1/8))
points(2-0.075,30,pch=21,col="gray90",bg="gray90",cex=3*(5/8))
points(rep(3,4)-0.075,trials2$Time2flee[trials2$Trial=="THSY"&trials2$Time2flee<30],pch=21,col="gray90",bg="gray90",cex=3*(1/8))
points(3-0.075,30,pch=21,col="gray90",bg="gray90",cex=3*(4/8))
points(2:3+0.075,exp.rates.fit5[2:3],pch=19,cex=1.5)
error.bar(x=(2:3)+0.075,y=exp.rates.fit5[2:3],upper=exp.rates.fit5[2:3],lower=time2flee.confint.fit5[1,2:3]
			,type="Y",length=0.01)
mtext("A.",adj=-0.05,line=0,cex=1)
mtext(c("Saddle-backed\\nTamarin","Bluish-slate\\nAntshrike"),at=2:3,line=-1,cex=0.8)
mtext("Screaming\\nPiha",at=1,line=-1,cex=0.8)

# Saguinus response
# Latency to response
plot(0,0,ylim=c(0,60),xlim=c(0.5,3.5),xlab="",ylab="Time until response (Sec)",xaxt="n",bty="l",cex.lab=1.2)
axis(1,at=1:3,labels=c("Control\\nn=7","Alarm\\nn=7","Alarm\\nn=7"),mgp=c(3,1.2,0),cex.axis=1)
points(rep(1,7)-0.075,trials$Time2flee[trials$Treatment=="Control"&trials$Time2flee==30],pch=21,col="gray90",bg="gray90",cex=3)
points(rep(2,2)-0.075,trials$Time2flee[trials$Treatment=="ST"&trials$Time2flee==30],pch=21,col="gray90",bg="gray90",cex=3*(2/7))
points(rep(2,3)-0.075,trials$Time2flee[trials$Treatment=="ST"&trials$Time2flee==1],pch=21,col="gray90",bg="gray90",cex=3*(3/7))
points(rep(2,2)-0.075,trials$Time2flee[trials$Treatment=="ST"&trials$Time2flee>1&trials$Time2flee<30],pch=21,col="gray90",bg="gray90",cex=3*(1/7))
points(rep(3,2)-0.075,trials$Time2flee[trials$Treatment=="THSY"&trials$Time2flee==30],pch=21,col="gray90",bg="gray90",cex=3*(2/7))
points(rep(3,2)-0.075,trials$Time2flee[trials$Treatment=="THSY"&trials$Time2flee==0.5],pch=21,col="gray90",bg="gray90",cex=3*(2/7))
points(rep(3,3)-0.075,trials$Time2flee[trials$Treatment=="THSY"&trials$Time2flee>0.5&trials$Time2flee<30],pch=21,col="gray90",bg="gray90",cex=3*(1/7))
points(2:3+0.075,exp.rates.fit2[2:3],pch=19,cex=1.5)
error.bar(x=(2:3)+0.075,y=exp.rates.fit2[2:3],upper=time2flee.confint.fit2[2,2:3],lower=time2flee.confint.fit2[1,2:3]
			,type="Y",length=0.01)
mtext("B.",adj=-0.05,line=0,cex=1)
mtext(c("Saddle-backed\\nTamarin","Bluish-slate\\nAntshrike"),at=2:3,line=-1,cex=0.8)
mtext("Screaming\\nPiha",at=1,line=-1,cex=0.8)

# Movement
plot(0,0,ylim=c(0,15),xlim=c(0.5,3.5),xlab="",ylab=expression(Delta~Distance(m)),xaxt="n",bty="l",cex.lab=1.2)
axis(1,at=1:3,labels=c("Control\\nn=7","Alarm\\nn=7","Alarm\\nn=7"),mgp=c(3,1.2,0),cex.axis=1)
points(rep(1,7)-0.1,trials$Delta.distance[trials$Treatment=="Control"],pch=21,col="gray90",bg="gray90",cex=3)
points(rep(2,2)-0.075,trials$Delta.distance[trials$Treatment=="ST"&trials$Delta.distance==0],pch=21,col="gray90",bg="gray90",cex=3*(2/7))
points(rep(2,5)-0.075,trials$Delta.distance[trials$Treatment=="ST"&trials$Delta.distance>0],pch=21,col="gray90",bg="gray90",cex=3*(1/7))
points(rep(3,2)-0.075,trials$Delta.distance[trials$Treatment=="THSY"&trials$Delta.distance==0],pch=21,col="gray90",bg="gray90",cex=3*(2/7))
points(rep(3,5)-0.075,trials$Delta.distance[trials$Treatment=="THSY"&trials$Delta.distance>0],pch=21,col="gray90",bg="gray90",cex=3*(1/7))
points(1:3+0.075,move.mean,pch=19,cex=1.5)
error.bar(x=(2:3)+0.075,y=move.mean[2:3],upper=distance.confint.fullmod[2,2:3],lower=distance.confint.fullmod[1,2:3]
			,type="Y",length=0.01,lwd=0.75)
mtext("C.",adj=-0.05,line=0,cex=1)
mtext(c("Saddle-backed\\nTamarin","Bluish-slate\\nAntshrike"),at=2:3,line=-1,cex=0.8)
mtext("Screaming\\nPiha",at=1,line=-1,cex=0.8)
dev.off()