###### This is the complete code for all talbes, figures, and analyses presented in the paper
	### "Experimental Tests of the Function & Flexibility of Song Consistency in a Wild Bird"
	### by Conor C Taff and Corey R Freeman-Gallant. Please see the accompanying
	### readme file for explanation of column headers.

###### Set the working directory to folder where data files are stored

	setwd("")
	
### Load packages to be used in analysis

	library(lme4)
	library(rethinking)
	library(bbmle)
	
###### Tell r whether to plot to tiff files or to the r window
	## if this value = 1, plots go to appropriately scaled tiffs

	plot.to.dev<-0
	
### MALE PLAYBACK DATA AND ANALYSIS

	### Read in the data files to be used in analyses
	
		male.d<-read.csv("MalePlayback.csv")
		
	### Create new columns that will be used in analyses below
	
		male.d$During<-ifelse(male.d$Period=="During",1,0)
		male.d$Post<-ifelse(male.d$Period=="Post",1,0)
		male.d$SRate<-ifelse(male.d$Rate=="Fast",1,0)
		male.d$SCons<-ifelse(male.d$Consistency=="High",1,0)
		for(i in 1:nrow(male.d)){male.d$PerIn5[i]<-
			male.d$In5[i]/sum(c(male.d$In5[i],male.d$X5to10[i],male.d$Over10[i]))*100}
		for(i in 1:nrow(male.d)){male.d$Per5to10[i]<-
			male.d$X5to10[i]/sum(c(male.d$In5[i],male.d$X5to10[i],male.d$Over10[i]))*100}
		for(i in 1:nrow(male.d)){male.d$PerOver10[i]<-
			male.d$Over10[i]/sum(c(male.d$In5[i],male.d$X5to10[i],male.d$Over10[i]))*100}
		male.d$lSongs<-log(male.d$Songs+1)
		male.d$lChips<-log(male.d$Chips+1)
		male.d$lRattles<-log(male.d$Rattles+1)
		male.d$lMoves<-log(male.d$Moves+1)
		male.d$Prox<-(male.d$PerIn5*0+male.d$Per5to10*1+male.d$PerOver10*2)/100
		male.d$cSongs<-(male.d$Songs-mean(male.d$Songs))/sd(male.d$Songs)
		male.d$cChips<-(male.d$Chips-mean(male.d$Chips))/sd(male.d$Chips)
		male.d$cRattles<-(male.d$Rattles-mean(male.d$Rattles))/sd(male.d$Rattles)
		male.d$cMoves<-(male.d$Moves-mean(male.d$Moves))/sd(male.d$Moves)
		male.d$cProx<-(male.d$Prox-mean(male.d$Prox))/sd(male.d$Prox)
			
	### Principal component analysis of male playback data based on log trasnformed songs, chips
		## rattles, moves from zone to zone, and (not log transformed) proximity score
	
		#pca.male<-male.d[,c(34:38)]
	
		#pc.male<-princomp(pca.male,scores=TRUE,cor=TRUE)
		#male.d$PC1<-pc.male$scores[,1]
		#male.d$PC2<-pc.male$scores[,2]
		
		pca.malec<-male.d[,c(39:43)]
		pc.male<-princomp(pca.malec,scores=TRUE,cor=TRUE)
		male.d$PC1<-pc.male$scores[,1]
		male.d$PC2<-pc.male$scores[,2]
		
	### Raw comparison of proximity distance without considering rate or consistency
		## Note that proximity score ranges from 0 (all time spent in 5 meters of speaker) to
		## 2 (all time spent greater than 10 meters from speaker).
	
		heights<-c(1.8973,1.087079,1.025556)
			
		barplot(heights,names.arg=c("Pre","During","Post"),ylim=c(0,2),
			ylab="Proximity Score",xlab="Playback Period")
		abline(h=0)
		
		bar<-0.03
		
		lines(rep(.7,2),c(1.8973-.059,1.8973+.059))
		lines(c(.7-bar,.7+bar),rep(1.8973-.059,2))
		lines(c(.7+bar,.7-bar),rep(1.8973+.059,2))
		
		lines(rep(1.9,2),c(1.087079-.0948,1.087079+.0948))
		lines(c(1.9-bar,1.9+bar),rep(1.087079-.0948,2))
		lines(c(1.9+bar,1.9-bar),rep(1.087079+.0948,2))
		
		lines(rep(3.1,2),c(1.02556-.1447,1.02556+.1447))
		lines(c(3.1-bar,3.1+bar),rep(1.02556-.1447,2))
		lines(c(3.1-bar,3.1+bar),rep(1.02556+.1447,2))
		
	### Fit glmm's with different response variables
	
		Response<-male.d$PC1  ## Change the response variable here to ramify through models below
			## options are: lSongs, lChips, lRattles, lMoves, Prox, PC1, PC2
		
		m<-lmer(Response~During+Post+SRate+SCons+(1|Male)+(1|TrialNumber)+(1|FocalStage),data=male.d)
		m.r.During<-lmer(Response~Post+SRate+SCons+(1|Male)+(1|TrialNumber)+(1|FocalStage),data=male.d)
		m.r.Post<-lmer(Response~During+SRate+SCons+(1|Male)+(1|TrialNumber)+(1|FocalStage),data=male.d)
		m.r.SRate<-lmer(Response~During+Post+SCons+(1|Male)+(1|TrialNumber)+(1|FocalStage),data=male.d)
		m.r.SCons<-lmer(Response~During+Post+SRate+(1|Male)+(1|TrialNumber)+(1|FocalStage),data=male.d)
				
	#Random effect
	
		m2<-lmer(Response~During+Post+SRate+SCons+(1|TrialNumber)+(1|FocalStage),data=male.d)
		anova(m,m2)
		
		summary(m)
		anova(m,m.r.During)
		anova(m,m.r.Post)
		anova(m,m.r.SRate)
		anova(m,m.r.SCons)
		
		post.m<-sample.naive.posterior(m,100000)
		
		## Point estimate for Pre-playback Response	
			pre.ci<-HPDI(post.m[,1]+post.m[,4]*mean(male.d$SRate)+post.m[,5]*mean(male.d$SCons))
			pre.mu<-mean(post.m[,1]+post.m[,4]*mean(male.d$SRate)+post.m[,5]*mean(male.d$SCons))			
		## Point estimate for During Playback Response		
			dur.ci<-HPDI(post.m[,1]+post.m[,2]+post.m[,4]*mean(male.d$SRate)+post.m[,5]*mean(male.d$SCons))
			dur.mu<-mean(post.m[,1]+post.m[,2]+post.m[,4]*mean(male.d$SRate)+post.m[,5]*mean(male.d$SCons))
			## Fast vs. slow holding consistency constant		
				## Fast
					dur.f.ci<-HPDI(post.m[,1]+post.m[,2]+post.m[,4]+post.m[,5]*mean(male.d$SCons))
					dur.f.mu<-mean(post.m[,1]+post.m[,2]+post.m[,4]+post.m[,5]*mean(male.d$SCons))				
				## Slow				
					dur.s.ci<-HPDI(post.m[,1]+post.m[,2]+post.m[,5]*mean(male.d$SCons))
					dur.s.mu<-mean(post.m[,1]+post.m[,2]+post.m[,5]*mean(male.d$SCons))					
			## Consistent vs. inconsistent holding rate constant			
				## Consistent					
					dur.h.ci<-HPDI(post.m[,1]+post.m[,2]+post.m[,4]*mean(male.d$SRate)+post.m[,5])
					dur.h.mu<-mean(post.m[,1]+post.m[,2]+post.m[,4]*mean(male.d$SRate)+post.m[,5])				
				## Inconsistent				
					dur.l.ci<-HPDI(post.m[,1]+post.m[,2]+post.m[,4]*mean(male.d$SRate))
					dur.l.mu<-mean(post.m[,1]+post.m[,2]+post.m[,4]*mean(male.d$SRate))			
		## Point estimate for Post Playback Response		
			post.ci<-HPDI(post.m[,1]+post.m[,3]+post.m[,4]*mean(male.d$SRate)+post.m[,5]*mean(male.d$SCons))
			post.mu<-mean(post.m[,1]+post.m[,3]+post.m[,4]*mean(male.d$SRate)+post.m[,5]*mean(male.d$SCons))		
			## Fast vs. slow holding consistency constant		
				## Fast		
					post.f.ci<-HPDI(post.m[,1]+post.m[,3]+post.m[,4]+post.m[,5]*mean(male.d$SCons))
					post.f.mu<-mean(post.m[,1]+post.m[,3]+post.m[,4]+post.m[,5]*mean(male.d$SCons))			
				## Slow			
					post.s.ci<-HPDI(post.m[,1]+post.m[,3]+post.m[,5]*mean(male.d$SCons))
					post.s.mu<-mean(post.m[,1]+post.m[,3]+post.m[,5]*mean(male.d$SCons))				
			## Consistent vs. inconsistent holding rate constant		
				## Consistent				
					post.h.ci<-HPDI(post.m[,1]+post.m[,3]+post.m[,4]*mean(male.d$SRate)+post.m[,5])
					post.h.mu<-mean(post.m[,1]+post.m[,3]+post.m[,4]*mean(male.d$SRate)+post.m[,5])				
				## Inconsistent			
					post.l.ci<-HPDI(post.m[,1]+post.m[,3]+post.m[,4]*mean(male.d$SRate))
					post.l.mu<-mean(post.m[,1]+post.m[,3]+post.m[,4]*mean(male.d$SRate))				
		## Build a table of these point estimates	
				m.table<-matrix(nrow=3,ncol=11)
				colnames(m.table)<-c("Pre","During","Post","Dur.Fast","Dur.Slow","Dur.Hi","Dur.Lo",
					"Post.Fast","Post.Slow","Post.Hi","Post.Lo")
				rownames(m.table)<-c("Mu","LowCI","HighCI")
				m.table[,1]<-c(pre.mu,pre.ci)
				m.table[,2]<-c(dur.mu,dur.ci)
				m.table[,3]<-c(post.mu,post.ci)
				m.table[,4]<-c(dur.f.mu,dur.f.ci)
				m.table[,5]<-c(dur.s.mu,dur.s.ci)
				m.table[,6]<-c(dur.h.mu,dur.h.ci)
				m.table[,7]<-c(dur.l.mu,dur.l.ci)
				m.table[,8]<-c(post.f.mu,post.f.ci)
				m.table[,9]<-c(post.s.mu,post.s.ci)
				m.table[,10]<-c(post.h.mu,post.h.ci)
				m.table[,11]<-c(post.l.mu,post.l.ci)				
				m.table	
				
			## Compare difference between pre and during playback
				post.m$Dif.pre.dur<-(post.m[,1]+post.m[,2]+post.m[,4]*mean(male.d$SRate)+
					post.m[,5]*mean(male.d$SCons))-
					(post.m[,1]+post.m[,4]*mean(male.d$SRate)+post.m[,5]*mean(male.d$SCons))
				j<-ifelse(mean(post.m$Dif.pre.dur)>0,1,0)
				ifelse(j==1,((100000-nrow(subset(post.m,post.m$Dif.pre.dur>0)))*2/100000),
					((100000-nrow(subset(post.m,post.m$Dif.pre.dur<0)))*2/100000))	
					
			## Compare difference between pre and post playback
				post.m$Dif.pre.post<-(post.m[,1]+post.m[,3]+post.m[,4]*mean(male.d$SRate)+
					post.m[,5]*mean(male.d$SCons))-
					(post.m[,1]+post.m[,4]*mean(male.d$SRate)+post.m[,5]*mean(male.d$SCons))
				j<-ifelse(mean(post.m$Dif.pre.post)>0,1,0)
				ifelse(j==1,((100000-nrow(subset(post.m,post.m$Dif.pre.post>0)))*2/100000),
					((100000-nrow(subset(post.m,post.m$Dif.pre.post<0)))*2/100000))	
					
			## Compare difference between during and post playback
				post.m$Dif.dur.post<-(post.m[,1]+post.m[,3]+post.m[,4]*mean(male.d$SRate)+
					post.m[,5]*mean(male.d$SCons))-
					(post.m[,1]+post.m[,2]+post.m[,4]*mean(male.d$SRate)+
					post.m[,5]*mean(male.d$SCons))
				j<-ifelse(mean(post.m$Dif.dur.post)>0,1,0)
				ifelse(j==1,((100000-nrow(subset(post.m,post.m$Dif.dur.post>0)))*2/100000),
					((100000-nrow(subset(post.m,post.m$Dif.dur.post<0)))*2/100000))	
					
	### T-test of SPCC before and after male experiments
	
		male.songs<-male.d[c(9,14,16,20,21,23,24,59,64,66,70:71,73:74),]
		
		t.test(subset(male.songs$Songs,male.songs$Period=="During"),
			subset(male.songs$Songs,male.songs$Period=="Pre"),paired=TRUE)
		
		t.test(subset(male.d$SPCC,male.d$Period=="Pre"),subset(male.d$SPCC,male.d$Period=="Post"),
			paired=TRUE)
			
	### barplot of PC1 and PC2 model predicted values for pre-during-post playback
	
		if(plot.to.dev==1){	
			pdf("MalePCFig.pdf",width=7.4,height=7.7)	
		}
	
		par(mfrow=c(1,2))
		
		## PC1
		
			heights<-c(0.86,-1.27,.18)
			low<-c(0.39,-1.73,-0.27)
			high<-c(1.32,-0.81,0.65)
				
			barplot(heights,names.arg=c("Pre","During","Post"),ylim=c(-2,3),
				ylab="PC1",xlab="Playback Period")
			abline(h=0)
			
			bar<-0.03
			
			lines(rep(.7,2),c(low[1],high[1]))
			lines(c(.7-bar,.7+bar),rep(low[1],2))
			lines(c(.7+bar,.7-bar),rep(high[1],2))
			
			lines(rep(1.9,2),c(low[2],high[2]))
			lines(c(1.9-bar,1.9+bar),rep(low[2],2))
			lines(c(1.9+bar,1.9-bar),rep(high[2],2))
			
			lines(rep(3.1,2),c(low[3],high[3]))
			lines(c(3.1-bar,3.1+bar),rep(low[3],2))
			lines(c(3.1-bar,3.1+bar),rep(high[3],2))
			
			lines(c(.7,1.9),rep(2,2))
			lines(c(.7,3.1),rep(2.3,2))
			lines(c(1.9,3.1),rep(2.6,2))
			
			text(1.3,2.1,"P < 0.001")
			text(1.9,2.4,"P = 0.002")
			text(2.5,2.7,"P < 0.001")
			
		## PC2
			
			heights<-c(.81,-.18,-.52)
			low<-c(.40,-.58,-.92)
			high<-c(1.21,0.23,-.12)
				
			barplot(heights,names.arg=c("Pre","During","Post"),ylim=c(-1,2),
				ylab="PC2",xlab="Playback Period",yaxp=c(-1,2,3))
			abline(h=0)
			
			bar<-0.03
			
			lines(rep(.7,2),c(low[1],high[1]))
			lines(c(.7-bar,.7+bar),rep(low[1],2))
			lines(c(.7+bar,.7-bar),rep(high[1],2))
			
			lines(rep(1.9,2),c(low[2],high[2]))
			lines(c(1.9-bar,1.9+bar),rep(low[2],2))
			lines(c(1.9+bar,1.9-bar),rep(high[2],2))
			
			lines(rep(3.1,2),c(low[3],high[3]))
			lines(c(3.1-bar,3.1+bar),rep(low[3],2))
			lines(c(3.1-bar,3.1+bar),rep(high[3],2))
			
			lines(c(.7,1.9),rep(1.4,2))
			lines(c(.7,3.1),rep(1.6,2))
			lines(c(1.9,3.1),rep(1.8,2))
			
			text(1.3,1.45,"P < 0.001")
			text(1.9,1.65,"P < 0.001")
			text(2.5,1.85,"ns")
			
			if(plot.to.dev==1){dev.off()}

	
### FEMALE PRESENTATION DATA AND ANALYSIS

	### Load the data
	
		fem.d<-read.csv("FemalePresent.csv")
		
	### Create new column headers before analysis
	
		fem.d$During<-ifelse(fem.d$Period=="During",1,0)
		fem.d$Post<-ifelse(fem.d$Period=="Post",1,0)
		for(i in 1:nrow(fem.d)){fem.d$PerIn5[i]<-round(
			fem.d$In5[i]/sum(c(fem.d$In5[i],fem.d$X5to10[i],fem.d$Over10[i]))*100)}
		for(i in 1:nrow(fem.d)){fem.d$Per5to10[i]<-round(
			fem.d$X5to10[i]/sum(c(fem.d$In5[i],fem.d$X5to10[i],fem.d$Over10[i]))*100)}
		for(i in 1:nrow(fem.d)){fem.d$PerOver10[i]<-round(
			fem.d$Over10[i]/sum(c(fem.d$In5[i],fem.d$X5to10[i],fem.d$Over10[i]))*100)}
		fem.d$Prox<-(fem.d$PerIn5*0+fem.d$Per5to10*1+fem.d$PerOver10*2)/100
		fem.d$lSongs<-log(fem.d$Songs+1)
		fem.d$lChips<-log(fem.d$Chips+1)
		fem.d$lRattles<-log(fem.d$Rattles+1)
		fem.d$lMoves<-log(fem.d$Moves+1)
		
	### Principal component analysis of behavioral results
		
		fem.d2<-fem.d[1:18,]
	
		pca.fem<-fem.d2[,c(21:25)]
	
		pc.fem<-princomp(pca.fem,scores=TRUE,cor=TRUE)
		fem.d2$PC1<-pc.fem$scores[,1]
		fem.d2$PC2<-pc.fem$scores[,2]
		
	### Simple t-tests to compare behaviors pre and during trial for female presentations
		
		t.test(subset(fem.d2$Songs,fem.d2$Period=="Pre"),subset(fem.d2$Songs,fem.d2$Period=="During"),
			paired=TRUE)
		t.test(subset(fem.d2$Prox,fem.d2$Period=="Pre"),subset(fem.d2$Prox,fem.d2$Period=="During"),
			paired=TRUE)	
		t.test(subset(fem.d2$Chips,fem.d2$Period=="Pre"),subset(fem.d2$Chips,fem.d2$Period=="During"),
			paired=TRUE)
		t.test(subset(fem.d2$Rattles,fem.d2$Period=="Pre"),subset(fem.d2$Rattles,fem.d2$Period=="During"),
			paired=TRUE)
		t.test(subset(fem.d2$Moves,fem.d2$Period=="Pre"),subset(fem.d2$Moves,fem.d2$Period=="During"),
			paired=TRUE)
			
		t.test(subset(fem.d2$PC1,fem.d2$Period=="Pre"),subset(fem.d2$PC1,fem.d2$Period=="During"),
			paired=TRUE)
		t.test(subset(fem.d2$PC2,fem.d2$Period=="Pre"),subset(fem.d2$PC2,fem.d2$Period=="During"),
			paired=TRUE)
		
	### Fit some models for SPCC
	
		m.Period<-lmer(SPCC~Post+During+(1|Male),data=fem.d)
		m.Post.r<-lmer(SPCC~During+(1|Male),data=fem.d)
		m.During.r<-lmer(SPCC~Post+(1|Male),data=fem.d)
		
		anova(m.Period,m.Post.r)
		anova(m.Period,m.During.r)
		
		post.p<-sample.naive.posterior(m.Period,100000)
		
		pre.ci<-HPDI(post.p[,1])
		pre.mu<-mean(post.p[,1])
		
		post.ci<-HPDI(post.p[,1]+post.p[,2])
		post.mu<-mean(post.p[,1]+post.p[,2])
		
		dur.ci<-HPDI(post.p[,1]+post.p[,3])
		dur.mu<-mean(post.p[,1]+post.p[,3])
		
		fem.table<-matrix(nrow=3,ncol=3)
		colnames(fem.table)<-c("Pre","During","Post")
		rownames(fem.table)<-c("Mu","LoCI","HiCI")
		
		fem.table[,1]<-c(pre.mu,pre.ci)
		fem.table[,2]<-c(dur.mu,dur.ci)
		fem.table[,3]<-c(post.mu,post.ci)
	
	### Barplot of SPCC by presentation period for female data
	
		if(plot.to.dev==1){
			pdf("SPCCFig.pdf",width=4.2,height=6.4)
		}
		
		barplot(fem.table[1,],names.arg=c("Pre","During","Post"),xlab="Presentation Period",
			ylim=c(0,1),ylab="Song Consistency (SPCC)",col="gray90")
		abline(h=0)
		
		bar<-0.03
		
		lines(rep(.7,2),c(fem.table[2,1],fem.table[3,1]))
		lines(c(.7-bar,.7+bar),rep(fem.table[2,1],2))
		lines(c(.7+bar,.7-bar),rep(fem.table[3,1],2))
		
		lines(rep(1.9,2),c(fem.table[2,2],fem.table[3,2]))
		lines(c(1.9-bar,1.9+bar),rep(fem.table[2,2],2))
		lines(c(1.9+bar,1.9-bar),rep(fem.table[3,2],2))
		
		lines(rep(3.1,2),c(fem.table[2,3],fem.table[3,3]))
		lines(c(3.1-bar,3.1+bar),rep(fem.table[2,3],2))
		lines(c(3.1-bar,3.1+bar),rep(fem.table[3,3],2))
		
		HPDI((post.p[,1]+post.p[,2])-(post.p[,1]+post.p[,3]))
		
		lines(c(0.6,3.2),rep(0.85,2))
		lines(c(1.8,3.2),rep(0.79,2))
		
		text(2.5,0.81,"* P < 0.01")
		text(1.9,0.87,"* P < 0.01")
		
		points(rep(.8,9),fem.d[1:9,15],pch=0,col="gray65")
		points(rep(1.9,9),fem.d[10:18,15],pch=0,col="gray65")
		points(rep(3.0,9),fem.d[19:27,15],pch=0,col="gray65")
		
		for(i in 1:9){
			lines(c(.8,1.9),c(fem.d[0+i,15],fem.d[9+i,15]),lty=2,col="gray60")
			lines(c(1.9,3.0),c(fem.d[9+i,15],fem.d[18+i,15]),lty=2,col="gray60")
		}
		
		if(plot.to.dev==1){dev.off()}
		
		
		
		
	