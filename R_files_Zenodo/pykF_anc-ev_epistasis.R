#pykF fitness competiitons -- fitness effects across backgrounds
library(lme4)
library(pwr)
library(multcomp)
library(car)

#input data
folder = "~/Dropbox/R/pykF_background/"
d = read.table(paste(folder,"pykF_anc-20K_data.txt", sep=""), header=T, skip=3)

#estimate fitness (r)
m.ref = log(d$ref.final*100^d$competition.days/d$ref.initial)
m.test = log(d$test.final*100^d$competition.days/d$test.initial)
r = m.test/m.ref
data = cbind(d,r)

#account for marker cost in 20K A-2 competitions
marker.control = mean(data[data$reference.ID=="REL8594A",ncol(data)]) - 1
data[which(data$reference.ID=="TC1512"),ncol(data)] <- data[which(data$reference.ID=="TC1512"),ncol(data)] + marker.control
#and remove marker control data
data=data[-which(data$reference.ID=="REL8594A"),]

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
error.bars=function(t,bg){
	avg=ifelse(bg=="anc",c("anc.r"),c("20K.r"))
	err=ifelse(bg=="anc",c("anc.CI"),c("20K.CI"))
	return(c(t[avg]-t[err],t[avg]+t[err]))	
}
	
#summary table of fitness values per genotype and time
alleles = unique(data$test.genotype)
times = unique(data$time)
fit.summary=matrix(NA, nrow=length(alleles), ncol=4)
for(j in 1:length(times)){
	for(i in 1: length(alleles)){
		fit.summary[i,j] = mean(data[data$test.genotype==alleles[i] & data$time==times[j],ncol(data)])
		fit.summary[i,j+2] = CI95(d=data[data$test.genotype==alleles[i] & data$time==times[j],ncol(data)])
	}
}
rownames(fit.summary)=alleles
colnames(fit.summary)=c("anc.r","20K.r","anc.CI","20K.CI")



#summary table of fitness values per mutation in ancestor
mutations = unique(data[data$time=="0","test.ID"])
fit.summary.anc=matrix(NA, nrow=length(mutations), ncol=2)
	for(i in 1: length(mutations)){
		fit.summary.anc[i,1] = mean(data[data$test.ID==mutations[i] & data$time=="0" & data$date=="110812",ncol(data)])
		fit.summary.anc[i,2] = CI95(d=data[data$test.ID==mutations[i] & data$time=="0" & data$date=="110812",ncol(data)])
	}

#order fit.summary rows for plotting
allele.order = c("WT","P70T","P70Q","D127N","I264F","A301S-1","A301S-2","A301S-3","A301T","G381A","T462I","Del")
allele.match = match(allele.order,row.names(fit.summary))
fit.summary = (fit.summary[allele.match,])
rownames(fit.summary.anc)=c("A301S","A301T","D127N","Del","G381A","I264F","P70Q","P70T","T462I","Anc")
colnames(fit.summary.anc)=c("anc.r","anc.CI")

#order fit.summary rows for plotting
allele.order.anc = c("Anc","P70T","P70Q","D127N","I264F","A301S","A301T","G381A","T462I","Del")
allele.match = match(allele.order.anc,row.names(fit.summary.anc))
fit.summary.anc = (fit.summary.anc[allele.match,])

#########################################
#ANOVA
#remove ancestor competition
data.NoAnc = data[-which(data$test.genotype=="WT"),]
#ancestor background
m=aov(r ~ test.ID + date,data=data.NoAnc[data.NoAnc$time==0,])
summary(m)
Anova(m, type=2)
#remove deletion allele
m1=aov(r ~ test.ID + date,data=data.NoAnc[data.NoAnc$time==0 & data.NoAnc$test.genotype!="Del",])
summary(m1)
#20K background
m3=aov(r ~ test.genotype,data=data.NoAnc[data.NoAnc$time==20,])
summary(m3)
#20K backgrounds -- compare mutator v non-mutator
m11 = aov(r ~ mutator,data=data.NoAnc[data.NoAnc$time==20,])
summary(m11)
#20K backgrounds -- mutation nested under mutator v non-mutator 
m12 = aov(r ~ mutator + mutator:test.ID,data=data.NoAnc[data.NoAnc$time==20,])
summary(m12)
#20K background remove deletion allele
m4=aov(r ~ test.genotype,data=data.NoAnc[data.NoAnc$time==20 & data.NoAnc$test.genotype!="Del",])
summary(m4)
#Comparison of 3 populations with A301S mutation
#...20K
d = data.NoAnc[data.NoAnc$time==20,]
m7=aov(r ~ test.genotype,data=d[d$test.genotype=="A301S-1" | d$test.genotype=="A301S-2" | d$test.genotype=="A301S-3",])
summary(m7)
#...Ancestor
d = data.NoAnc[data.NoAnc$time==0 & data.NoAnc$date=="110812",]
m8=aov(r ~ test.genotype,data=d[d$test.genotype=="A301S-1" | d$test.genotype=="A301S-2" | d$test.genotype=="A301S-3",])
summary(m8)

#remove 2 of the 3 sets of A301S estimates from first block (i.e., leaving 4 estimates in each block, as for all other ancestor-mutation combinations)
d.0.noAnc = data[data$time=="0" & data$test.genotype!="WT",]
d.0.noAnc.8A301S = d.0.noAnc[-c(9:16),]
#one way
m9 = aov(r ~ test.ID + date,data=d.0.noAnc.8A301S[d.0.noAnc.8A301S$test.genotype!="Del",])
summary(m9)
Anova(m9, type = 2)
#nested under domain in ancestor -- 8 A301S mutations
m10 = aov(r ~ Domain + date + Domain:test.ID,data=d.0.noAnc.8A301S[d.0.noAnc.8A301S$test.genotype!="Del",])
summary(m10)
Anova(m10, type=2)

#test!!!!!
library(ez)
ezANOVA(data=data.NoAnc, dv = r, wid = ref.initial, between = .(time, test.genotype), type=3, detailed= T)
m11 = lm(r ~ time * test.genotype, data=data.NoAnc)
drop1(m11, .~. , test="F")
Anova(m11, type = 2)
anova(m11)



#Anc background mutation nested under domain omit deletion allele
#*Reported in manuscript
m6=aov(r ~ date + Domain + Domain:test.ID,data=data.NoAnc[data.NoAnc$time==0 & data.NoAnc$test.genotype!="Del",])
summary(m6)
Anova(m6, type=2)
#20K background mutation nested under domain omit deletion allele
#*Reported in manuscript
m5=aov(r ~ Domain/test.genotype,data=data.NoAnc[data.NoAnc$time==20 & data.NoAnc$test.genotype!="Del",])
summary(m5)
Anova(m5, type=2)
#Size of effects -- Anc
#omega^2: not accounting for domain
sum_m1 = summary(m1)[[1]]
SSm.m1 = sum_m1[["Sum Sq"]][1]
SSt.m1 = sum_m1[["Sum Sq"]][1]+sum_m1[["Sum Sq"]][2]+sum_m1[["Sum Sq"]][3]
DFm.m1 = sum_m1[["Df"]][1]
MSr.m1 = sum_m1[["Mean Sq"]][3] #includes data blocking term
W2.anc = (SSm.m1-DFm.m1*MSr.m1)/(SSt.m1+MSr.m1)

#mean pair-wise difference between alleles
#Anc
HSD.anc = TukeyHSD(aov(m1),"test.ID")
mean(abs(HSD.anc$test.ID[,1])) *100
#20K
HSD.ev = TukeyHSD(aov(m4),"test.genotype")
mean(abs(HSD.ev$test.genotype[,1])) *100

#Dunnetts test - comparing deletion effect to point mutations - TC1091 is deletion allele
summary(glht(m, linfct = mcp(test.ID = "Dunnett")))

#omega^2: domain
omega.domain.anc=(summary(m6)[[1]][[2]][1] - (summary(m6)[[1]][[1]][1]*summary(m6)[[1]][[3]][3]))/(sum(summary(m6)[[1]][[2]]) + summary(m6)[[1]][[3]][3])
#omega^2: background nested in domain
omega.background.anc=(summary(m6)[[1]][[2]][2] - (summary(m6)[[1]][[1]][2]*summary(m6)[[1]][[3]][3]))/(sum(summary(m6)[[1]][[2]]) + summary(m6)[[1]][[3]][3])

#Size of effects -- 20K
#omega^2: not accounting for domain
sum_m4 = summary(m4)[[1]]
SSm.m4 = sum_m4[["Sum Sq"]][1]
SSt.m4 = sum_m4[["Sum Sq"]][1]+sum_m1[["Sum Sq"]][2]
DFm.m4 = sum_m4[["Df"]][1]
MSr.m4 = sum_m4[["Mean Sq"]][2] 
W2.20K = (SSm.m4-DFm.m4*MSr.m4)/(SSt.m4+MSr.m4)
#mean pair-wise difference between alleles
HSD.20K = TukeyHSD(aov(m4),"test.genotype")
mean(abs(HSD.20K$test.genotype[,1]))

#20K
#omega^2: domain
omega.domain.20K=(summary(m5)[[1]][[2]][1] - (summary(m5)[[1]][[1]][1]*summary(m5)[[1]][[3]][3]))/(sum(summary(m5)[[1]][[2]]) + summary(m5)[[1]][[3]][3])
#omega^2: background nested in domain
omega.background.20K=(summary(m5)[[1]][[2]][2] - (summary(m5)[[1]][[1]][2]*summary(m5)[[1]][[3]][3]))/(sum(summary(m5)[[1]][[2]]) + summary(m5)[[1]][[3]][3])

#Anc
omega.domain.anc=(summary(m6)[[1]][[2]][1] - (summary(m6)[[1]][[1]][1]*summary(m6)[[1]][[3]][3]))/(sum(summary(m6)[[1]][[2]]) + summary(m6)[[1]][[3]][3])
#omega^2: background nested in domain
omega.background.anc=(summary(m6)[[1]][[2]][2] - (summary(m6)[[1]][[1]][2]*summary(m6)[[1]][[3]][3]))/(sum(summary(m6)[[1]][[2]]) + summary(m6)[[1]][[3]][3])


TukeyHSD(aov(m5),"Domain") #20K
TukeyHSD(aov(m6),"Domain") #Anc

#paired t-test, allele in anc and evolved background
#comparison only includes estimates collected in the same experimental block
ttest=vector(mode="numeric", length=length(allele.order))
for(i in 2:length(allele.order)){ #start at 2 to skip WT 
	temp = data[data$test.genotype==allele.order[i],]
	temp1 = t.test(temp[temp$time==0 & temp$date=="110812","r"], temp[temp$time==20,"r"], paired=F, var.equal=F)
	ttest[i] = temp1$p.value
}
#########################################
#mean mutation effects
point.mutations=mean(data.NoAnc[data.NoAnc$test.genotype != "Del" & data.NoAnc$reference.ID == "REL607",ncol(data.NoAnc)])
print(paste("Mean point mutation effect in Anc: ", round(point.mutations,3),sep=""))
deletion.mutation=mean(data[data$test.genotype == "Del" & data$reference.ID == "REL607",ncol(data)])
print(paste("Mean deletion mutation effect in Anc: ", round(deletion.mutation,3),sep=""))

#########################################
#plots
fit.summary.noWT=fit.summary[-1,]

#Two panel -- effects in ancestor and effects in 20K
fig<-c("pykF_anc-ev_epistasis-2panel-ANC-20K.pdf")
pdf(paste(folder,fig,sep=""),onefile=F,width=4,height=5)
par(oma=c(0.75,0,0,0))
par(mar=c(2.25,3,0.5,0.5))
par(mfrow=c(2,1))
#ancestor panel
fit.summary.noWT.anc=fit.summary.anc[-1,]
plot(fit.summary.noWT.anc[,1],
	xlab="", ylab="",
	xaxt='n', yaxt='n',
	ylim=c(0.95,1.33), xlim=c(1,9.5),
	col="white")
for(i in 1:nrow(fit.summary.noWT.anc)){	
	lines(c(i,i),error.bars(t= fit.summary.noWT.anc[i,], bg="anc"))	
}
points(fit.summary.noWT.anc[,1], pch=21, lwd=2, col=ifelse(rownames(fit.summary.noWT.anc)=="WT","grey", "red"), bg="white")
axis(1,at=(1:nrow(fit.summary.noWT.anc)),label=F)
axis(2,at=c(1,1.15,1.3),label=F)
x.labels = allele.order.anc[2:10]
x.labels[9]<-"Indel"
text((1:nrow(fit.summary.noWT.anc)), par("usr")[3]-0.03,x.labels, xpd=T, srt=45, cex=0.7, adj=c(1,1))	
mtext(c(0,0.15,0.3),2, at=c(1,1.15,1.3), line=0.65, cex=0.85)
text(1,1.315,"A")

#evolved panel
par(mar=c(3.25,3,0.25,0.5))
plot(fit.summary.noWT[,2],
	xlab="", ylab="",
	xaxt='n', yaxt='n',
	ylim=c(0.95,1.33), xlim=c(1,11.5),
	col="white")
for(i in 1:nrow(fit.summary.noWT)){	
	lines(c(i,i),error.bars(t= fit.summary.noWT[i,], bg="20K"))	
}
points(1:11, fit.summary.noWT[,2], pch=16, col="blue", cex=1.1)
axis(1,at=(1:nrow(fit.summary.noWT)),label=F)
axis(2,at=c(1,1.15,1.3),label=F)
#text((1:nrow(fit.summary.noWT)), par("usr")[3]-0.03,allele.order[2:12], xpd=T, srt=45, cex=0.65, adj=c(1,1))	
x.labels.ev = allele.order[c(2:5,9:12)]
x.labels.ev[8] <- "Indel" 
#text(c(1:4,8:11), par("usr")[3]-0.03,x.labels.ev, xpd=T, srt=45, cex=0.7, adj=c(1,1))	
par(lheight=.8)
#text(c(4.8,5.85,6.9), par("usr")[3]-0.03,c("A301S\\n (Ara-5)", "A301S\\n (Ara+1)", "A301S\\n (Ara+5)"), xpd=T, srt=45, cex=0.6, adj=c(1,1))	
text(seq(0.9,12, by=1.02), par("usr")[3]-0.03,
	c("P70T\\n (Ara+6)","P70Q\\n (Ara+3)","D127N\\n (Ara-3)","I264F\\n (Ara-2)","A301S\\n (Ara-5)", "A301S\\n (Ara+1)", "A301S\\n (Ara+5)","A301T\\n (Ara-6)","G381A\\n (Ara+2)","T462I\\n (Ara-4)","Indel\\n (Ara-1)"), 
	xpd=T, srt=45, cex=0.6, adj=c(1,1))	


mtext(c(0,0.15,0.3),2, at=c(1,1.15,1.3), line=0.65, cex=0.85)
mtext(expression(italic(pykF)~mutation),1, line=3)
mtext(expression(Effect~of~italic(pykF)~mutation~(s)),2, line=1.8, adj=-3.5)
#mtext("Relative fitness",2, line=2, adj=2.25)
#asterisks to indicate difference in fitness across backgrounds
for(i in 2:length(ttest)){
		text(i-1,0.96,ifelse(ttest[i]<=0.001,"***",
			ifelse(ttest[i]<=0.01,"**",
			ifelse(ttest[i]<=0.05,"*", ""))))
}
text(1,1.315,"B")

dev.off()

#Alternative combined panel -- Anc. A301S as line
fig<-c("pykF_anc-ev_epistasis-ANC_1A301S_20K-_anc.only.pdf")
pdf(paste(folder,fig,sep=""),onefile=F,width=4,height=3.5)
par(mar=c(4,3.5,0.5,0.5))
jit=0.1

fit.summary.noWT.anc=fit.summary.anc[-1,]
plot(fit.summary.noWT.anc[,1],
	xlab="", ylab="",
	xaxt='n', yaxt='n',
	ylim=c(0.95,1.33), xlim=c(1,11.5),
	col="white")
for(i in 1:nrow(fit.summary.noWT.anc)){	
	if(i<=5){lines(c(i,i),error.bars(t= fit.summary.noWT.anc[i,], bg="anc"))
		}else{
			lines(c(i+2,i+2),error.bars(t= fit.summary.noWT.anc[i,], bg="anc"))}	
}
lines(c(5,7),c(fit.summary.noWT.anc["A301S",1],fit.summary.noWT.anc["A301S",1]), col="red", lwd=2)
points(c(1:5,8:11),fit.summary.noWT.anc[,1], pch=21, lwd=2, col=ifelse(rownames(fit.summary.noWT)=="WT","grey", "red"), bg="white")

#text((1:nrow(fit.summary.noWT.anc)), par("usr")[3]-0.03,allele.order.anc[2:10], xpd=T, srt=45, cex=0.7, adj=c(1,1))	

#for(i in 1:nrow(fit.summary.noWT)){	
#	lines(c(i+jit,i+jit),error.bars(t= fit.summary.noWT[i,], bg="20K"))	
#}
#points(1:11+jit, fit.summary.noWT[,2], pch=16, col="blue", cex=1.1)
axis(1,at=(1:nrow(fit.summary.noWT)),label=F)
axis(2,at=c(1,1.15,1.3),label=F)
#text((1:nrow(fit.summary.noWT)), par("usr")[3]-0.03,allele.order[2:12], xpd=T, srt=45, cex=0.65, adj=c(1,1))	
text(c(1:4,8:11), par("usr")[3]-0.03,allele.order[c(2:5,9:12)], xpd=T, srt=45, cex=0.7, adj=c(1,1))	
par(lheight=.8)
text(c(4.8,5.85,6.9), par("usr")[3]-0.03,c("A301S\\n (Ara-5)", "A301S\\n (Ara+1)", "A301S\\n (Ara+5)"), xpd=T, srt=45, cex=0.6, adj=c(1,1))	

mtext(c(1.00,1.15,1.3),2, at=c(1,1.15,1.3), line=0.65, cex=0.85)
mtext(expression(italic(pykF)~mutation),1, line=3)
mtext("Relative fitness",2, line=2, adj=0.45)
#asterisks to indicate difference in fitness across backgrounds
#for(i in 2:length(ttest)){
#		text(i-1,0.96,ifelse(ttest[i]<=0.001,"***",
#			ifelse(ttest[i]<=0.01,"**",
#			ifelse(ttest[i]<=0.05,"*", ""))))
#}
dev.off()



#effects in combined ancestor and 20K
fig<-c("pykF_anc-ev_epistasis-ANC-20K.pdf")
pdf(paste(folder,fig,sep=""),onefile=F,width=4,height=3.5)

par(mar=c(5,4,0.5,0.5))
jit=0.1
plot(fit.summary.noWT[,1],
	xlab="", ylab="",
	xaxt='n', yaxt='n',
	ylim=c(0.95,1.35), xlim=c(1,11.5),
	col="white")
for(i in 1:nrow(fit.summary.noWT)){	
	lines(c(i,i),error.bars(t= fit.summary.noWT[i,], bg="anc"))	
}
for(i in 1:nrow(fit.summary.noWT)){	
	lines(c(i+jit,i+jit),error.bars(t= fit.summary.noWT[i,], bg="20K"))	
}
points(fit.summary.noWT[,1], pch=21, lwd=2, col=ifelse(rownames(fit.summary.noWT)=="WT","grey", "red"), bg="white")
points(1:11+jit, fit.summary.noWT[,2], pch=16, col="blue", cex=1.1)
axis(1,at=(1:nrow(fit.summary.noWT)),label=F)
axis(2,at=c(1,1.15,1.3),label=F)
text((1:nrow(fit.summary.noWT)), par("usr")[3]-0.03,allele.order[2:12], xpd=T, srt=45, cex=0.8, adj=c(1,0))	
mtext(c(1.00,1.15,1.3),2, at=c(1,1.15,1.3), line=0.65, cex=0.85)
mtext(expression(italic(pykF)~allele),1, line=3.5)
mtext("Relative fitness",2, line=2.2)
#asterisks to indicate difference in fitness across backgrounds
for(i in 2:length(ttest)){
		text(i-1,0.96,ifelse(ttest[i]<=0.001,"***",
			ifelse(ttest[i]<=0.01,"**",
			ifelse(ttest[i]<=0.05,"*", ""))))
}

dev.off()


#effects in ancestor
fig<-c("pykF_anc-ev_epistasis-ANC.pdf")
pdf(paste(folder,fig,sep=""),onefile=F,width=4,height=3.5)

par(oma=c(0,0,0,0))
par(mar=c(5,4,0.5,0.5))
plot(fit.summary[,1],
	xlab="", ylab="",
	xaxt='n', yaxt='n',
	ylim=c(1,1.17),
	col="white")
for(i in 1:nrow(fit.summary)){	
	lines(c(i,i),error.bars(t=fit.summary[i,], bg="anc"))	
}
points(fit.summary[,1], pch=1, col=ifelse(rownames(fit.summary)=="WT","grey", "red"), lwd=2)
axis(1,at=(1:nrow(fit.summary)),label=F)
axis(2,at=c(1,1.05,1.1),label=F)
text((1:nrow(fit.summary)), par("usr")[3]-0.011,allele.order, xpd=T, srt=45, cex=0.8, adj=c(1,0))	
mtext(c(1.00,1.05,1.1),2, at=c(1,1.05,1.1), line=0.65, cex=0.85)
mtext(expression(italic(pykF)~allele),1, line=3.5)
mtext("Relative fitness",2, line=2.2)

dev.off()

#effects in 20K
fig<-c("pykF_anc-ev_epistasis-20K.pdf")
pdf(paste(folder,fig,sep=""),onefile=F,width=4,height=3.5)

par(oma=c(0,0,0,0))
par(mar=c(5,4,0.5,0.5))
plot(fit.summary[,2],
	xlab="", ylab="",
	xaxt='n', yaxt='n',
	ylim=c(0.95,1.35), col="blue", cex=1.1,
	pch=16)
for(i in 1:nrow(fit.summary)){	
	lines(c(i,i),error.bars(t=fit.summary[i,], bg="20K"))	
}
axis(1,at=(1:nrow(fit.summary)),label=F)
axis(2,at=c(1,1.15,1.3),label=F)
text((1:nrow(fit.summary)), par("usr")[3]-0.03,allele.order, xpd=T, srt=45, cex=0.8, adj=c(1,0))	
mtext(c(1.00,1.15,1.30),2, at=c(1,1.15,1.3), line=0.65, cex=0.85)
mtext(expression(italic(pykF)~allele),1, line=3.5)
mtext("Relative fitness",2, line=2.2)

dev.off()

