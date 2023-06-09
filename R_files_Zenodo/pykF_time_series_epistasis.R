#pykF fitness competiitons -- fitness effects across backgrounds
library(nlme)
library(minpack.lm)
library(lmerTest)
library(afex)
library(VCA)
library(car)

#input data
folder = "~/Dropbox/R/pykF_background/"
d = read.table(paste(folder,"pykF_time_series_data.txt", sep=""), header=T, skip=4)

#estimate fitness (r)
m.ref = log(d$ref.final*100^d$competition.days/d$ref.initial)
m.test = log(d$test.final*100^d$competition.days/d$test.initial)
r = m.test/m.ref
data = cbind(d,r)

#account for marker cost in 20K A+2 5K competitions
marker.control = mean(data[data$population=="A+2" & data$time=="ARA",ncol(data)]) - 1
data[which(data$population=="A+2" & data$time=="5"),ncol(data)] <- data[which(data$population=="A+2" & data$time=="5"),ncol(data)] + marker.control
#and remove A+2 marker control data
data=data[-which(data$population=="A+2" & data$time=="ARA"),]
#A+1 5K as difference between pykF-Ev and pykF-Anc estimates vs. ancestor
pykF.ev = which(data$population=="A+1" & data$time=="EV")
pykF.anc = which(data$population=="A+1" & data$time=="ANC")
data[pykF.ev,ncol(data)] <- data[pykF.ev,ncol(data)] - (data[pykF.anc,ncol(data)] -1)
#add proper time variable
data[pykF.ev,"time"]<-"5"
#except, A+1 pykF-anc data has very few counts at beginning time point -- delete
#delete A+1 5k Anc rows
data=data[-c(pykF.anc,pykF.ev),]

#SEM function
SEM=function(d){
	return(sd(d)/sqrt(length(d)))
}

#95CI function
CI95=function(d){
	SE=sd(d)/sqrt(length(d))
	return(SE*qt(0.975,df=length(d)-1))
}

#error bar function
error.bars=function(t,offset){
	return(c(t[offset]-t[6+offset],t[offset]+t[6+offset]))	
}
	
#summary table of fitness values per genotype and time
populations = unique(data$population)
times = unique(data$time)
fit.summary=matrix(NA, nrow=length(times), ncol=(2*length(populations)))
for(j in 1:length(populations)){
	for(i in 1: length(times)){
		fit.summary[i,j] = mean(data[data$population==populations[j] & data$time==times[i],ncol(data)])
		fit.summary[i,j+length(populations)] = SEM(d=data[data$population==populations[j] & data$time==times[i],ncol(data)])
	}
}
rownames(fit.summary)=times
colnames(fit.summary)=c(paste(populations,".r",sep=""),paste(populations,".95CI",sep=""))

#remove ancestor data from summary
fit.summary.noWT=fit.summary[-7,]
fit.summary.noWT= fit.summary.noWT[,-c(7,14)]

#########################################
#ANCOVA
#remove ancestor competition
data.NoAnc = data[-c(which(data$population=="606"),(which(data$time=="anc"))) ,]
data.Anc = data[-c(which(data$population=="606")) ,]

A.times=matrix(NA, ncol=1, nrow=nrow(data.NoAnc))
for(i in 1:nrow(data.NoAnc)){
	A.times[i,1]=as.numeric(as.character(data.NoAnc[i,3]))
}

A.r = data.NoAnc[,ncol(data.NoAnc)]

A.population = data.NoAnc[,"population"]

A.data = data.frame(A.times,A.r,A.population)

#####################################################
#ANOVA

#global ANOVA testing for time and population*time interaction effect
m4 = aov(r ~ population + time + block + population:time, data = data.NoAnc)
m5 = aov(r ~ population + time + block + population:time, data = data.Anc)
summary(m4)
Anova(m4, type=2) #type II SS

summary(m5)
summary(Anova(m5, type=2)) #type II SS

##Same, but with population as a random effect
#Crossed random effects are handled in an easier way than in lme (+(1|a)+(1|b)).
m4.r.full = lmer(r ~ time + (1|population) + (1|population:time), data = data.Anc, REML=F)
summary(m4.r.full)
m4.r.noInt = lmer(r ~ time + (1|population), data = data.Anc, REML=F)
m4.r.noPop = lm(r ~ time, data = data.Anc)
m4.r.noTime = lmer(r ~ 1 + (1|population), data = data.Anc)

summary(m4.r.full)
summary(m4.r.noInt)
anova(m4.r.full,m4.r.noInt)
anova(m4.r.full,m4.r.noPop)
anova(m4.r.full,m4.r.noTime)

#using afex to give p-values
mixed(A.r ~ A.times + (1|A.population), data = A.data)



#fit ANOVA seperately for each time point
#test if there is variation over time
d2=data
d2[,3]=as.factor(d2[,3])
d2 = d2[-c(which(d2$population=="606")),]
times = unique(d2$time)

#population as a fixed effect
ANOVA.summary = NULL
for(i in 1:length(times)){
	m5 = summary(aov(r ~ population, data = d2[d2$time==times[i],]))
	p= m5[[1]][[5]][1]
	F.stat= m5[[1]][[4]][1]
	SS= m5[[1]][[2]][1]
	MS= m5[[1]][[3]][1]
	df = c(m5[[1]][[1]][1],m5[[1]][[1]][2])
	ANOVA.summary = rbind(ANOVA.summary,c(df=df,SS=SS,MS=MS,F.stat=F.stat,p=p))
}
rownames(ANOVA.summary) = paste("gen",times,sep="")

#population as a random effect
ANOVA.r.summary = NULL
for(i in 1:length(times)){
	#normalize data to mean
	mod.data = d2[d2$time==times[i],]
	mod.data["r"]<- mod.data["r"]/mean(mod.data[,ncol(mod.data)])
	model= lmer(r ~ 1+ (1|population) , data = mod.data)
	m5 = summary(model)
	ci = confint(model, method="profile")
	ci.2 = ci[c(1,4)]
	VCAinference(anovaMM(r ~ 1+ (population), d2[d2$time==times[i],]))
	ANOVA.r.summary = rbind(ANOVA.r.summary,c(vc=sqrt(m5$varcor[[1]][1])*100,lci=ci.2[1],uci=ci.2[2]))
}
rownames(ANOVA.r.summary) = paste("gen",times,sep="")


#####################################################
#plots and fit ANCOVA seperately for each population
#It's probably only useful to test for a trend and whether there is evidence for a peaked relationship
time = c(2,5,10,20,30,50)
par(mfrow=c(2,3))

for(i in 1:length(populations)){
y = fit.summary.noWT[1:6,paste(populations[i],".r", sep="")]
w = 1/fit.summary.noWT[1:6,paste(populations[i],".SEM", sep="")]

m.linear.single = aov(y ~ time , weights=w)
summary(m.linear.single)

m.quad.single = aov(y ~ poly(time,2,raw=T), weights = w)
summary(m.quad.single)

#Asymmetric peaked function
#m.gamma.single =  nlsLM(y ~ max*exp(-time/tau1)*(1-exp(-(time-toff)/tau2)), start=list(max=1.15,tau1=4,tau2=2,toff=2), weights=w, control = list(maxiter = 100))
#summary(m.gamma.single)

#print(AIC(m.linear.single,m.quad.single, m.gamma.single))
print(anova(m.linear.single, m.quad.single))

#plot single curve and model fits
par(mar=c(4,4,0.5,0.5))

plot.times=c(2,5,10,20,30,50,0)
plot(0,
	xlab="", ylab="",
	xaxt='n', yaxt='n',
	ylim=c(0.9,1.4), xlim=c(-5,55),
	col="white")
for(j in 1:nrow(fit.summary.noWT)){
		error = fit.summary.noWT[j,paste(populations[i],".SEM", sep="")]
		fit.estimate = fit.summary.noWT[j,paste(populations[i],".r", sep="")]
		lines(c(plot.times[j],plot.times[j]),c(fit.estimate-error, fit.estimate+error))	
	}

color=c("blue","purple","khaki","dark green","green","red","black")
#for(j in 1:(ncol(fit.summary.noWT)/2)){
		points(plot.times, fit.summary.noWT[,paste(populations[i],".r", sep="")], col=color[i], pch=16)	
		
#	}

axis(1,at=(plot.times),label=FALSE)
axis(2,at=c(1,1.15,1.3),label=FALSE)
text(plot.times, par("usr")[3]-0.045,plot.times, xpd=T, srt=0, cex=0.85)	
mtext(c(1.00,1.15,1.3),2, at=c(1,1.15,1.3), line=0.5, cex=0.85)
mtext("Generations (x1000)",1, line=2)
mtext("Relative fitness",2, line=2.2)
legend("topright",pch=16,bty='n',legend=populations[i], cex=0.75, col=color[i])

#plot model predictions
#curve(predict(m.quad.single, newdata = data.frame(time = x)), from = 5, to = 50, add = TRUE, lty=2, col="blue")

#curve(predict(m.gamma.single, newdata = data.frame(time = x)), from = 5, to = 50, add = TRUE, lty=2, col="dark green")

#curve(predict(m.linear.single, newdata = data.frame(time = x)), from = 5, to = 50, add = TRUE, lty=2, col="red")


}

########################################
#Variance in pykF effect explained by genotype at each time point - omega^2

omega = NULL
test.times=as.vector(unique(A.times))
for(i in 1:length(test.times)){
	m6 = summary(aov(A.r ~ A.population, data = A.data[A.data$A.times == test.times[i],]))	
	temp = (m6[[1]][[2]][1] - m6[[1]][[1]][1] * m6[[1]][[3]][2]) / (sum(m6[[1]][[2]]) + m6[[1]][[3]][2])
	temp2 = 100*round(temp,4)
	omega=c(omega,temp2)
}
message("% variance explained by genotype at each time point:")
message(paste("\\nTime",test.times,":",omega,sep=" "))

#and using SD to test overall variation
stdev = NULL
test.times=as.vector(unique(A.times))
for(i in 1:length(test.times)){
	m6 =  A.data[A.data$A.times == test.times[i],]
	#take means by population
	temp = tapply(m6$A.r, m6$A.population, FUN= mean)
	temp2 = sd(na.omit(as.vector(temp)))
	stdev=c(stdev,temp2)
}
message("standard deviation at each time point:")
message(paste("\\nTime",test.times,":",stdev,sep=" "))

#########################################
#plots
#effects in combined ancestor and 20K

#plot 1 -- complex plot with fitness & omega^2 
fig<-c("pykF_time_series_epistasis&vc.SEM.pdf")
pdf(paste(folder,fig,sep=""),onefile=FALSE,width=4,height=4)

par(mar=c(0,3.5,0.5,0.5))
par(fig=c(0,1,0.35,1))

jit=0.1
#remove ancestor data
fit.summary.noWT=fit.summary[-7,]
fit.summary.noWT= fit.summary.noWT[,-c(7,14)]

plot.times=c(2,5,10,20,30,50,0)
symbols = c(15,1, 17, 18 ,16 ,5)

jit=0.1
plot(0,
	xlab="", ylab="",
	xaxt='n', yaxt='n',
	ylim=c(0.9,1.4), xlim=c(-5,55),
	col="white")
for(j in 1:(ncol(fit.summary.noWT)/2)){
	for(i in 1:nrow(fit.summary.noWT)){	
		lines(c(plot.times[i]+(jit*j),plot.times[i]+(jit*j)), error.bars(t= fit.summary.noWT[i,], offset=j), col="dark grey")	
	}
}
color=c("blue","purple","khaki","dark green","green","red","black")
for(j in 1:(ncol(fit.summary.noWT)/2)){
		points(plot.times+(jit*j), fit.summary.noWT[,j], col=color[j], pch=symbols[j])	
	}

legend("topright",bty='n',legend=c("Ara-1", "Ara-6", "Ara+1", "Ara+2", "Ara+3", "Ara+5"), cex=0.75, col=color, pch=symbols)
axis(1,at=(plot.times),label=FALSE, tck=+0.03)
axis(2,at=c(1,1.15, 1.3),label=FALSE)
mtext(c(0,0.15,0.3),2, at=c(1,1.15, 1.3), line=0.5, cex=0.85)
mtext(expression(Effect~of~italic(pykF)~mutation~(s)),2, line=2)
#text(-4,1.36,"A", font=2)

par(mar=c(4,3.5,0.1,0.5))
par(fig=c(0,1,0,0.35), new=TRUE)
#convert any negative omega values to zero
omega[omega<0] <- 0
#population = fixed
#plot(plot.times[1:6],omega[1:6]/100, xlab="", ylab="", xlim=c(-5,55), ylim=c(-0.05,1), col="black", pch=16,xaxt="n", yaxt="n")

#population = random
plot(plot.times[1:7],ANOVA.r.summary[,1], xlab="", ylab="", xlim=c(-5,55), ylim=c(-1,15), col="black", pch=16,xaxt="n", yaxt="n")
for(i in 1:length(plot.times)){
	lines(c(plot.times[i],plot.times[i]),c(ANOVA.r.summary[i,2]*100,ANOVA.r.summary[i,3]*100))
}

axis(1,at=(plot.times),label=FALSE)
axis(2,at=c(0,10),label=FALSE)
text(plot.times, par("usr")[3]-6,plot.times, xpd=T, srt=0, cex=0.85)	
mtext(c(0,10),2, at=c(0,10), line=0.5, cex=0.85)
mtext(expression("Generations ("%*%"1000)"),1, line=1.8)
#mtext(expression(omega ^ 2),2, line=2)
mtext(expression(sigma %*%'100'),2, line=2)
#text(-4,12,"B", font=2)

dev.off()

#plot 2 -- simple plot with fitness only
fig<-c("pykF_time_series_epistasis.pdf")
pdf(paste(folder,fig,sep=""),onefile=FALSE,width=4,height=3.5)

par(mar=c(4,3.5,0.5,0.5))
jit=0.1
#remove ancestor data
fit.summary.noWT=fit.summary[-7,]
fit.summary.noWT= fit.summary.noWT[,-c(7,14)]

plot.times=c(2,5,10,20,30,50,0)
symbols = c(15,1, 17, 18 ,16 ,5)

jit=0.1
plot(0,
	xlab="", ylab="",
	xaxt='n', yaxt='n',
	ylim=c(0.9,1.4), xlim=c(-5,55),
	col="white")
for(j in 1:(ncol(fit.summary.noWT)/2)){
	for(i in 1:nrow(fit.summary.noWT)){	
		lines(c(plot.times[i]+(jit*j),plot.times[i]+(jit*j)),error.bars(t= fit.summary.noWT[i,], offset=j), col="dark grey")	
	}
}
color=c("blue","purple","khaki","dark green","green","red","black")
for(j in 1:(ncol(fit.summary.noWT)/2)){
		points(plot.times+(jit*j), fit.summary.noWT[,j], col=color[j], pch=symbols[j])	
	}

axis(1,at=(plot.times),label=FALSE)
axis(2,at=c(1,1.15,1.3),label=FALSE)
text(plot.times, par("usr")[3]-0.045,plot.times, xpd=T, srt=0, cex=0.85)	
mtext(c(0,0.15,0.3),2, at=c(1,1.15,1.3), line=0.5, cex=0.85)
mtext("Generations (x1000)",1, line=1.8)
mtext(expression(Effect~of~italic(pykF)~mutation),2, line=2)
legend("topright",bty='n',legend=c("Ara-1", "Ara-6", "Ara+1", "Ara+2", "Ara+3", "Ara+5"), cex=0.75, col=color, pch=symbols)

#plot model predictions
x=seq(2,55, length.out=100)
p=predict(m3,list(A.times=x))
pred2=data.frame(x,p)[order(x),]
#lines(pred2[,1]+1,pred2[,2],lty=2, col="black")


dev.off()

##try to determine groups of fitness effects
fit = fit.summary.noWT[,1:6]
#turn to distance matrix then Hclust?
clust = hclust(dist(t(fit), method="euclidean"))
plot(clust)

