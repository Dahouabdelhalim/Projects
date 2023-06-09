###### This is the complete code for all talbes, figures, and analyses presented in the paper
	### "Fluctuations in neighbourhood fertility generate variable signaling effort"
	### by Conor C Taff, Gail L Patricell, and Corey R Freeman-Gallant. Please see the accompanying
	### readme file for explanation of column headers. Code was last run in R for Mac OS X
	### version 3.0.2 GUI 1.62 Snow Leopard Build.

###### Set the working directory to folder where data files are stored

	setwd("")

###### Load packages that will (or might) be used. Other dependent packages loaded automatically 
	## include Matrix v 1.1-0; MASS v 7.3-29; and lattice v 0.20-23

	library(lme4)  ## v 1.0-5
	library(bbmle)  ## v 1.0.15
	library(rethinking)  ## v 1.30
	library(coda)  ## v 0.16-1
	library(plotrix)  ## v 3.5-1
	library(glmmADMB)  ## v 0.8.0

###### Load data and parse into the sets that will be used

	d<-read.csv("ByMale.csv")    ## All Data
	d$JDate<-d$JDate-120  ## Recalculate date so that May 1 = 1
	
	dawnd<-subset(d,d$IncludeDawn=="Yes")   # Exclude days with no or partial dawn recording
	dawnd2<-subset(dawnd,dawnd$Imp==0)   # Exclude males with testosterone implates
	dawnd3<-subset(dawnd2,dawnd2$TotalDawn>1)    # Exclude days with 0 detected focal songs
	dawnd4<-subset(dawnd,dawnd$TotalDawn>1)  # Dataset including T implanted males
	
	dayd<-subset(d,d$FullDay=="Yes")   # Exclude days with no or partial day recording
	dayd2<-subset(dayd,dayd$Imp==0)    # Exclude males with testosterone implants
	dayd3<-subset(dayd2,dayd2$TotalDay>1)  # Exclude days with 0 detected focal songs
	dayd4<-subset(dayd,dayd$TotalDay>1) # Dataset including T implanted males
	
###### Tell r whether to plot to tiff files or to the r window
	## if this value = 1, plots go to appropriately scaled tiffs

	plot.to.dev<-1

###### Define a function to create AICc tables from glmmadmb fit models

		AICcTable<-function(x){
			OutName<-matrix(nrow=length(x),ncol=8)
			names<-c("LogL","K","AICc","dAIC","expAIC","Weight","N","Groups")
			colnames(OutName)<-names
			rownames(OutName)<-x
			likelihood<-vector(length=length(x))
			for(i in 1:length(x)){
				likelihood[i]<-as.numeric(get(x[i])$loglik[1])
			}
			parameters<-vector(length=length(x))
			for(i in 1:length(x)){
				parameters[i]<-as.numeric(get(x[i])$npar[1])
			}	
			n<-vector(length=length(x))
			for(i in 1:length(x)){
				n[i]<-as.numeric(get(x[i])$n[1])
			}
			groups<-vector(length=length(x))
			for(i in 1:length(x)){
				groups[i]<-as.numeric(get(x[i])$q[1])
			}
			OutName[,1]<-likelihood
			OutName[,2]<-parameters
			OutName[,3]<-(-2)*OutName[,1]+2*OutName[,2]+
				((2*OutName[,2]*(OutName[,2]+1))/
					(30-OutName[,2]-1))
			OutName[,4]<-round(OutName[,3]-min(OutName[,3]),2)
			OutName[,5]<-exp(-.5*OutName[,4])
			stand<-sum(OutName[,5])
			OutName[,7]<-n
			OutName[,8]<-groups
			OutName[,6]<-round(OutName[,5]/stand,3)
			temp<-order(OutName[,4])
			OutName<-OutName[temp,]	
		}
		
###### Table 1: Candidate model sets for dawn and daytime singing

	##### Dawn singing models
	
			dawndata<-dawnd3  # subset of data to use for model fits
			
		## 'Null' model with intercept only
			dawn.m0<-glmmadmb(TotalDawn~1+(1|Code1),family="nbinom",data=dawndata)
		## Abiotic factors only model with temperature and date
			dawn.m1<-glmmadmb(TotalDawn~NightLow+DawnHigh+JDate+(1|Code1),family="nbinom",data=dawndata)
		## Within pair social context only
			dawn.m2<-glmmadmb(TotalDawn~FF+NF+(1|Code1),family="nbinom",data=dawndata)
		## Both within and neighboring social context considered
			dawn.m3<-glmmadmb(TotalDawn~F400+FF+NF+(1|Code1),family="nbinom",data=dawndata)
		## Full model including all predictors
			dawn.m4<-glmmadmb(TotalDawn~F400+FF+NF+NightLow+DawnHigh+JDate+(1|Code1),family="nbinom",data=dawndata)
			
		Dawn.Table<-AICcTable(c("dawn.m0","dawn.m1","dawn.m2","dawn.m3","dawn.m4"))
		
	##### Daytime singing models
	
			daydata<-dayd3  # subset of data to use for model fits
			
		## 'Null' model with intercept only
			day.m0<-glmmadmb(TotalDay~1+(1|Code1),family="nbinom",data=daydata)
		## Abiotic factors only model with temperature and date
			day.m1<-glmmadmb(TotalDay~NightLow+DayHigh+JDate+(1|Code1),family="nbinom",data=daydata)
		## Within pair social context only
			day.m2<-glmmadmb(TotalDay~FF+NF+(1|Code1),family="nbinom",data=daydata)
		## Both within and neighboring social context considered
			day.m3<-glmmadmb(TotalDay~F400+FF+NF+(1|Code1),family="nbinom",data=daydata)
		## Full model including all predictors
			day.m4<-glmmadmb(TotalDay~F400+FF+NF+NightLow+DayHigh+JDate+(1|Code1),family="nbinom",data=daydata)
			
		Day.Table<-AICcTable(c("day.m0","day.m1","day.m2","day.m3","day.m4"))
		
####### Not reported in detail in the paper. This is the exact same analysis as in Table 1 above, except that males
	## that received Testosterone implants are also included. Each candidate model now includes a predictor variable
	## that specifies whether a male had a testosterone implant or not.		

	##### Dawn singing models
	
			dawndata2<-dawnd4  # subset of data to use for model fits
			
		## 'Null' model with intercept only
			dawn.m0.1<-glmmadmb(TotalDawn~1+(1|Code1),family="nbinom",data=dawndata2)
		## Abiotic factors only model with temperature and date
			dawn.m1.1<-glmmadmb(TotalDawn~NightLow+DawnHigh+JDate+Imp+(1|Code1),family="nbinom",data=dawndata2)
		## Within pair social context only
			dawn.m2.1<-glmmadmb(TotalDawn~FF+NF+Imp+(1|Code1),family="nbinom",data=dawndata2)
		## Both within and neighboring social context considered
			dawn.m3.1<-glmmadmb(TotalDawn~F400+FF+NF+Imp+(1|Code1),family="nbinom",data=dawndata2)
		## Full model including all predictors
			dawn.m4.1<-glmmadmb(TotalDawn~F400+FF+NF+NightLow+DawnHigh+JDate+Imp+(1|Code1),family="nbinom",data=dawndata2)
			
		Dawn.Table2<-AICcTable(c("dawn.m0.1","dawn.m1.1","dawn.m2.1","dawn.m3.1","dawn.m4.1"))
		
	##### Daytime singing models
	
			daydata2<-dayd4  # subset of data to use for model fits
			
		## 'Null' model with intercept only
			day.m0.1<-glmmadmb(TotalDay~1+(1|Code1),family="nbinom",data=daydata2)
		## Abiotic factors only model with temperature and date
			day.m1.1<-glmmadmb(TotalDay~NightLow+DayHigh+JDate+Imp+(1|Code1),family="nbinom",data=daydata2)
		## Within pair social context only
			day.m2.1<-glmmadmb(TotalDay~FF+NF+Imp+(1|Code1),family="nbinom",data=daydata2)
		## Both within and neighboring social context considered
			day.m3.1<-glmmadmb(TotalDay~F400+FF+NF+Imp+(1|Code1),family="nbinom",data=daydata2)
		## Full model including all predictors
			day.m4.1<-glmmadmb(TotalDay~F400+FF+NF+NightLow+DayHigh+JDate+Imp+(1|Code1),family="nbinom",data=daydata2)
			
		Day.Table2<-AICcTable(c("day.m0.1","day.m1.1","day.m2.1","day.m3.1","day.m4.1"))

###### Table 2: 

	### Fit model coefficients from best supported models above...	
	
		summary(dawn.m3)
		summary(day.m3)	
		
#### Figure 1. singing profiles by context

	#### NOTE: this figure uses three separate raw data files that have recordings in bins rather than as daily sums.
	#### These bins are just taken from raw observation data with all males pooled.
	#### NOTE: references to 'z' stage here and throughout are code for unmated stage
	### NOTE: times are given in seconds in this file midrise = seconds after sunrise, midpoint = seconds after midnight

		## Read data files
		
			dprofile<-read.csv("ByBin.csv")
			
		## Manipulate data to get ready for plotting
			
			dprofile$Hour<-dprofile$MidRise/60/60
			
			dprofile$z.upper<-(dprofile$Z.Mean+dprofile$Z.Std.Error)/2.5
			dprofile$z.lower<-(dprofile$Z.Mean-dprofile$Z.Std.Error)/2.5
			dprofile$f.upper<-(dprofile$F.Mean+dprofile$F.Std.Error)/2.5
			dprofile$f.lower<-(dprofile$F.Mean-dprofile$F.Std.Error)/2.5
			dprofile$nf.upper<-(dprofile$NF.Mean+dprofile$NF.Std.Error)/2.5
			dprofile$nf.lower<-(dprofile$NF.Mean-dprofile$NF.Std.Error)/2.5
			
			h<-dprofile$Hour
			
			col.1<-rgb(230,159,0,maxColorValue=255)
			col.2<-rgb(86,180,233,maxColorValue=255)
			col.3<-rgb(0,158,115,maxColorValue=255)
			
			## Plot loess smoothed lines for each context
			
				if(plot.to.dev==1){
					tiff("Figure1.tiff",width=10,height=4,units="in",compression="none",res=500)
				}
			
				plot(dprofile$Z.Mean/2.5~dprofile$Hour,type="n",ylim=c(0,5),xlab="Hours After Sunrise",ylab="Songs Rate (per minute)")
				mz<-loess(dprofile$Z.Mean/2.5~dprofile$Hour,se=TRUE,span=.2)
				zloe<-predict(mz,h,se=TRUE)
				lines(h,zloe$fit,col=col.2,lwd=2,lty=4)
				
				mnf<-loess(dprofile$NF.Mean/2.5~dprofile$Hour,se=TRUE,span=.2)
				nfloe<-predict(mnf,h,se=TRUE)
				lines(h,nfloe$fit,col=col.1,lwd=2,lty=5)
				
				mf<-loess(dprofile$F.Mean/2.5~dprofile$Hour,se=TRUE,span=.2)
				floe<-predict(mf,h,se=TRUE)
				lines(h,floe$fit,col=col.3,lwd=2)
			
			### Minus std error ofmean
			
				mzl<-loess(dprofile$z.lower~dprofile$Hour,se=TRUE,span=.2)
				zloel<-predict(mzl,h,se=TRUE)
				lines(h,zloel$fit,col=col.2,lwd=1,lty=3)
				
				mnfl<-loess(dprofile$nf.lower~dprofile$Hour,se=TRUE,span=.2)
				nfloel<-predict(mnfl,h,se=TRUE)
				lines(h,nfloel$fit,col=col.1,lwd=1,lty=3)
				
				mfl<-loess(dprofile$f.lower~dprofile$Hour,se=TRUE,span=.2)
				floel<-predict(mfl,h,se=TRUE)
				lines(h,floel$fit,col=col.3,lwd=1,lty=3)
			
			## Plus std error of mean
			
				mzu<-loess(dprofile$z.upper~dprofile$Hour,se=TRUE,span=.2)
				zloeu<-predict(mzu,h,se=TRUE)
				lines(h,zloeu$fit,col=col.2,lwd=1,lty=3)
				
				mnfu<-loess(dprofile$nf.upper~dprofile$Hour,se=TRUE,span=.2)
				nfloeu<-predict(mnfu,h,se=TRUE)
				lines(h,nfloeu$fit,col=col.1,lwd=1,lty=3)
				
				mfu<-loess(dprofile$f.upper~dprofile$Hour,se=TRUE,span=.2)
				floeu<-predict(mfu,h,se=TRUE)
				lines(h,floeu$fit,col=col.3,lwd=1,lty=3)
				
			## Plotting details
			
				legend(11.4,5,c("Unmated","Fertile Mate","Post-Fertile Mate","Recording Period"),
					col=c(col.2,col.3,col.1,"black"),lty=c(4,1,5,1),
					lwd=c(2,2,2,2))
					
				points(h,rep(0,length(h)),pch=15,cex=.5)	
				
				if(plot.to.dev==1){dev.off()}
		
#### Figure 2. Empirical observations of singing by social context

	## Dawn Chorus
	
		dawn.z.mean<-mean(subset(dawndata$TotalDawn,dawndata$Un==1))
		dawn.z.se<-sd(subset(dawndata$TotalDawn,dawndata$Un==1))/sqrt(length(subset(dawndata$TotalDawn,dawndata$Un==1)))
		
		dawn.f.mean<-mean(subset(dawndata$TotalDawn,dawndata$FF==1))
		dawn.f.se<-sd(subset(dawndata$TotalDawn,dawndata$FF==1))/sqrt(length(subset(dawndata$TotalDawn,dawndata$FF==1)))
		
		dawn.nf.mean<-mean(subset(dawndata$TotalDawn,dawndata$NF==1))
		dawn.nf.se<-sd(subset(dawndata$TotalDawn,dawndata$NF==1))/sqrt(length(subset(dawndata$TotalDawn,dawndata$NF==1)))
		
	## Daytime Song
	
		day.z.mean<-mean(subset(daydata$TotalDay,daydata$Un==1))
		day.z.se<-sd(subset(daydata$TotalDay,daydata$Un==1))/sqrt(length(subset(daydata$TotalDay,daydata$Un==1)))
		
		day.f.mean<-mean(subset(daydata$TotalDay,daydata$FF==1))
		day.f.se<-sd(subset(daydata$TotalDay,daydata$FF==1))/sqrt(length(subset(daydata$TotalDay,daydata$FF==1)))
		
		day.nf.mean<-mean(subset(daydata$TotalDay,daydata$NF==1))
		day.nf.se<-sd(subset(daydata$TotalDay,daydata$NF==1))/sqrt(length(subset(daydata$TotalDay,daydata$NF==1)))
		
	## Combine into one matrix
	
		social.est<-matrix(nrow=3,ncol=6)
		colnames(social.est)<-c("Dawn.Z","Dawn.F","Dawn.NF","Day.Z","Day.F","Day.NF")
		rownames(social.est)<-c("Mean","Mean-SE","Mean+SE")
		social.est[1,]<-c(dawn.z.mean,dawn.f.mean,dawn.nf.mean,day.z.mean,day.f.mean,day.nf.mean)
		social.est[2,]<-c(dawn.z.mean-dawn.z.se,dawn.f.mean-dawn.f.se,dawn.nf.mean-dawn.nf.se,
			day.z.mean-day.z.se,day.f.mean-day.f.se,day.nf.mean-day.nf.se)
		social.est[3,]<-c(dawn.z.mean+dawn.z.se,dawn.f.mean+dawn.f.se,dawn.nf.mean+dawn.nf.se,
			day.z.mean+day.z.se,day.f.mean+day.f.se,day.nf.mean+day.nf.se)
			
	### Plot the two song averages
	
		if(plot.to.dev==1){
			tiff("Figure2.tiff",width=8.6,height=4.7,units="in",compression="none",res=450)
		}
	
		par(mfrow=c(1,2))
		labels<-c("Unmated","Fertile","Post-Fertile")
		bar<-0.06   # Set width of error bar top and bottoms
		
		barplot(social.est[1,1:3],ylab="Songs Produced",ylim=c(0,200),names.arg=labels,xlab="Dawn Chorus")
			
			lines(c(0.7,0.7),social.est[2:3,1])
			lines(c(0.7-bar,0.7+bar),rep(social.est[3,1],2))
			lines(c(0.7-bar,0.7+bar),rep(social.est[2,1],2))
			
			lines(c(1.9,1.9),social.est[2:3,2])
			lines(c(1.9-bar,1.9+bar),rep(social.est[3,2],2))
			lines(c(1.9-bar,1.9+bar),rep(social.est[2,2],2))
			
			lines(c(3.1,3.1),social.est[2:3,3])
			lines(c(3.1-bar,3.1+bar),rep(social.est[3,3],2))
			lines(c(3.1-bar,3.1+bar),rep(social.est[2,3],2))
			
			corner.label(label="(a)")
			
		barplot(social.est[1,4:6],ylab="Songs Produced",ylim=c(0,500),names.arg=labels,xlab="Daytime")
	
			lines(c(0.7,0.7),social.est[2:3,4])
			lines(c(0.7-bar,0.7+bar),rep(social.est[3,4],2))
			lines(c(0.7-bar,0.7+bar),rep(social.est[2,4],2))
			
			lines(c(1.9,1.9),social.est[2:3,5])
			lines(c(1.9-bar,1.9+bar),rep(social.est[3,5],2))
			lines(c(1.9-bar,1.9+bar),rep(social.est[2,5],2))
			
			lines(c(3.1,3.1),social.est[2:3,6])
			lines(c(3.1-bar,3.1+bar),rep(social.est[3,6],2))
			lines(c(3.1-bar,3.1+bar),rep(social.est[2,6],2))
			
			corner.label(label="(b)")
			
			if(plot.to.dev==1){dev.off()}
			
#### Figure 3. Plot predictions from day.m3 acros neighborhood fertility

		if(plot.to.dev==1){
			tiff("Figure3.tiff",width=9.8,height=5.3,units="in",res=450,compression="none")
		}
		
		par(mfrow=c(1,2))
		
	## Plot the subset of empirical observations of dawn song during unmated stage only
		## Note these are just for illustration, lines are from fit model
	
		dawn.plot<-subset(dawndata,dawndata$Un==1)
		plot(dawn.plot$F400,dawn.plot$TotalDawn,col=col.alpha("slateblue",0.5),pch=16,
			xlab="Fertile Females in 400m",ylab="Dawn Chorus Songs") 
		
	## Calculate and plot MLE and +/- SE of parameter estimates from best model	
			
		r<-seq(from=-0.5,to=8,by=0.1)
		m.dawn<-sapply(r,function(z)exp(fixef(dawn.m3)[1]+fixef(dawn.m3)[2]*z))
		c.l.dawn<-sapply(r,function(z)exp(fixef(dawn.m3)[1]-stdEr(dawn.m3)[1]
			+(fixef(dawn.m3)[2]*z-stdEr(dawn.m3)[2])))
		c.h.dawn<-sapply(r,function(z)exp(fixef(dawn.m3)[1]+stdEr(dawn.m3)[1]
			+(fixef(dawn.m3)[2]*z+stdEr(dawn.m3)[2])))
		
		lines(r,m.dawn)	
		lines(r,c.l.dawn,lty=2)
		lines(r,c.h.dawn,lty=2)
		
	## Plot empirical means and stderrors for z observations 
	
		dawn.m.se<-matrix(nrow=4,ncol=7)
		for(i in 1:7){
			dawn.m.se[1,i]<-mean(subset(dawn.plot$TotalDawn,dawn.plot$F400==i-1))
		}
		serrors<-vector(length=7)
		for(i in 1:7){
			serrors[i]<-sd(subset(dawn.plot$TotalDawn,dawn.plot$F400==i-1))/
				sqrt(length(subset(dawn.plot$TotalDawn,dawn.plot$F400==i-1)))
		}
		dawn.m.se[2,]<-dawn.m.se[1,]+serrors
		dawn.m.se[3,]<-dawn.m.se[1,]-serrors
		for(i in 1:7){
			dawn.m.se[4,i]<-length(subset(dawn.plot$TotalDawn,dawn.plot$F400==i-1))
		}
		
		for(i in 1:6){
			lines(rep(i-1+0.1,2),c(dawn.m.se[2,i],dawn.m.se[3,i]))
			points(i-1+0.1,dawn.m.se[1,i],pch=16)
		}
		points(6+0.1,dawn.m.se[1,7],pch=16)		
		
		corner.label("(a)")

	## Plot the subset of empirical observations of full day song during unmated stage only
		## Note these are just for illustration, lines are from fit model
	
		day.plot<-subset(daydata,daydata$Un==1)
		plot(day.plot$F400,day.plot$TotalDay,col=col.alpha("slateblue",0.5),pch=16,
			xlab="Fertile Females in 400m",ylab="Daytime Songs") 
		
	## Calculate and plot MLE and +/- SE of parameter estimates from best model	
			
		r<-seq(from=-0.5,to=8,by=0.1)
		m.day<-sapply(r,function(z)exp(fixef(day.m3)[1]+fixef(day.m3)[2]*z))
		c.l.day<-sapply(r,function(z)exp(fixef(day.m3)[1]-stdEr(day.m3)[1]
			+(fixef(day.m3)[2]*z-stdEr(day.m3)[2])))
		c.h.day<-sapply(r,function(z)exp(fixef(day.m3)[1]+stdEr(day.m3)[1]
			+(fixef(day.m3)[2]*z+stdEr(day.m3)[2])))
		
		lines(r,m.day)	
		lines(r,c.l.day,lty=2)
		lines(r,c.h.day,lty=2)
		
	## Plot empirical means and stderrors for z observations 
	
		day.m.se<-matrix(nrow=4,ncol=7)
		for(i in 1:7){
			day.m.se[1,i]<-mean(subset(day.plot$TotalDay,day.plot$F400==i-1))
		}
		serrors<-vector(length=7)
		for(i in 1:7){
			serrors[i]<-sd(subset(day.plot$TotalDay,day.plot$F400==i-1))/
				sqrt(length(subset(day.plot$TotalDay,day.plot$F400==i-1)))
		}
		day.m.se[2,]<-day.m.se[1,]+serrors
		day.m.se[3,]<-day.m.se[1,]-serrors
		for(i in 1:7){
			day.m.se[4,i]<-length(subset(day.plot$TotalDay,day.plot$F400==i-1))
		}
		
		for(i in 1:6){
			lines(rep(i-1+0.1,2),c(day.m.se[2,i],day.m.se[3,i]))
			points(i-1+0.1,day.m.se[1,i],pch=16)
		}
		points(6+0.1,day.m.se[1,7],pch=16)	
		
		corner.label("(b)")	
		
	if(plot.to.dev==1){dev.off()}
	
					