
library(paleofire)

#set working directoryâ€”need to change for your own uses
setwd('/Users/patrickrafter/Documents/R/R/global14C/r_test_for_Zenmodo_and_Github')
# you have to set this to your working directory


#import data
dat<- read.csv('/Users/patrickrafter/Documents/R/R/global14C/r_test_for_Zenmodo_and_Github/rafter-2022-Global-14C-Compilation-FIN.csv', skip=1) #Data needs to be downloaded from PANGAEA here: https://doi.pangaea.de/10.1594/PANGAEA.946522
# this upload only begins reading on the 2nd row because the first row contains metadata about each column

dat<- subset(dat, neutral.density >= 27.5 & DELTA14Cage.atmos != 0 & threesigma1=='0' & sort=='0')
# only selects the data deeper than the 27.5 density surface, non-zero values, values that have been previously calculated to be >3 sigma the binned average, and data failing the preliminary curation (i.e., values with 14C/C higher than the contemporary atmosphere)

colnames(dat)
# dat<- subset(dat, neutral.density >= 27.5 & DELTA14Cage.atmos != -999)

intcal20<-read.table("/Users/patrickrafter/Documents/R/R/intcal/14C-intcal20.14c.txt")
colnames(intcal20)<- c("age","14C_age","error","D14C","sigma")
 

basins<- c('ALL','Atlantic','Pacific','SO')
depths<- c('ALL','deep','mid')
cutoffdepth<- 28
iteration.number<- 1000


pdf(paste('plot_1-MCMC-Age-Density-bootstrap-simpleplot.pdf',sep=''),width=7.2, height=9.6, useDingbats=FALSE, encoding="WinAnsi")
#dev.new(width=7.2, height=9.6)
par(mfrow=c(4,3))
par(mar=c(3,3,0.5,0.5),xpd=FALSE); par(ps = 9, cex = 1, cex.main = 1); par(mgp=c(1.5,0.5,0)); par(las=1); par(tck=-0.02)











for(k in 1:length(basins)){
	for(j in 1:length(depths)){
#k<- 1; j<-1
		basin<- basins[k]
		depth<- depths[j]

#select basin
#ALL, NA=0.11, NP=0.26, SoHemi=0.63
#basin<- 'NA'

#select depth
#all, intermediate=<2000, deep=>2000
#depth<- 'deep'


		#bin interval settings
		start<- 0
		end<- 27000
		interval<- 500

		int_s<- seq(start, end-interval, by=interval)
		int_e<- seq(start+interval, end, by=interval)
		int_age<- colMeans(rbind(int_s, int_e))
	
		#removes global outliers from basin bins
		if(j == 1 & k == 1) {dat<- dat} else {
			outliersALL<- read.csv('csv-14Cbin_outliersALL_ALLbootstrap.csv')
			dat<- subset(dat, !(DELTA14Cage.atmos %in% outliersALL$DELTA14Cage.atmos))
			}

		#subsets basin data 
		dat2<- if(basin == 'ALL') {dat} else
    		  #if(basin == 'NA') {subset(dat, latitude > -20)} else
    		  if(basin == 'Atlantic') {subset(dat, separate.the.basins == 1)} else
    		  if(basin == 'Pacific') {subset(dat, separate.the.basins == 2)} else
    		  #if(basin == 'Indian') {subset(dat, separate.the.basins == 3)} else
		  if(basin == 'SO') {subset(dat, separate.the.basins == 4)} 
		    		  
		dat2<- if(depth == 'ALL') {dat2} else
				if(depth == 'mid') {subset(dat2, neutral.density <= cutoffdepth )} else
				if(depth == 'deep') {subset(dat2, neutral.density > cutoffdepth)}
		
		
		plot(-999, -999, xlim=c(0,25000),ylim=rev(c(0,5000)), ylab='14C_vent_age', sep='',
xlab='calendar age (yr)') 
	axis(3,labels=FALSE,tick=TRUE, line=FALSE);axis(4,labels=FALSE,tick=TRUE, line=FALSE)
	grid (NULL,NULL, lty = 6, col = adjustcolor("cornsilk2", alpha=0.7))
	

		#bins data
		r<-matrix(ncol=6,nrow=length(int_s))
		
		s<- matrix(ncol=ncol(dat),nrow=0)
		s1<- matrix(ncol=ncol(dat),nrow=0)
		
		for (i in 1: nrow(r)){
			#i<- 1
			d<- subset(dat2, cal.age >= int_s[i] & cal.age < int_e[i])
			d<- subset(d,DELTA14Cage.atmos > -999 )
			#makes sure not to include the "missing data" values
			int_mean<- mean(d$DELTA14Cage.atmos, na.rm=TRUE)
			int_sd<- sd(d$DELTA14Cage.atmos, na.rm=TRUE)
			r[i,1]<- int_mean
			r[i,2]<- int_sd
			r[i,3]<- length(d$DELTA14Cage.atmos)
			points(d$cal.age, d$DELTA14Cage.atmos, pch=16, cex=0.7, col=adjustcolor('grey47', alpha=0.3))
			
			if(length(d$DELTA14Cage.atmos) > 1) { #this includes the data if there is only one data point in bin
				d2<- subset(d, DELTA14Cage.atmos > int_mean-3*int_sd & DELTA14Cage.atmos < int_mean+3*int_sd)
				int_mean3s<- mean(d2$DELTA14Cage.atmos, na.rm=TRUE)
				int_sd3s<- sd(d2$DELTA14Cage.atmos, na.rm=TRUE)
					d2.1<- subset(d2, DELTA14Cage.atmos > int_mean3s-3*int_sd3s & DELTA14Cage.atmos < int_mean3s+3*int_sd3s)
					int_mean3s2<- mean(d2.1$DELTA14Cage.atmos, na.rm=TRUE)
					int_sd3s2<- sd(d2.1$DELTA14Cage.atmos, na.rm=TRUE)
						d2.2<- subset(d2.1, DELTA14Cage.atmos > int_mean3s2-3*int_sd3s2 & DELTA14Cage.atmos < int_mean3s2+3*int_sd3s2)
						int_mean3s3<- mean(d2.2$DELTA14Cage.atmos, na.rm=TRUE)
						int_sd3s3<- sd(d2.2$DELTA14Cage.atmos, na.rm=TRUE)
							d2.3<- subset(d2.2, DELTA14Cage.atmos > int_mean3s3-3*int_sd3s3 & DELTA14Cage.atmos < int_mean3s3+3*int_sd3s3)
				
				
				
				r[i,4] <- mean(d2$DELTA14Cage.atmos, na.rm=TRUE)
				r[i,5] <- sd(d2$DELTA14Cage.atmos, na.rm=TRUE)
				r[i,6]	<- length(d2$DELTA14Cage.atmos)
				#points(d2$cal.age, d2$DELTA14Cage.atmos, pch=16, cex=0.7, col=adjustcolor('grey37', alpha=0.3))
				#d3<- subset(d, DELTA14Cage.atmos < int_mean-3*int_sd | DELTA14Cage.atmos > int_mean+3*int_sd)
				#points(d3$cal.age, d3$DELTA14Cage.atmos, pch=16, cex=0.8, col=adjustcolor('red3', alpha=0.6))
				#d3<- subset(d, DELTA14Cage.atmos < int_mean3s2-3*int_sd3s2 | DELTA14Cage.atmos > int_mean3s2+3*int_sd3s2)
				#points(d3$cal.age, d3$DELTA14Cage.atmos, pch=16, cex=0.8, col=adjustcolor('red3', alpha=0.6))
				
				d3.1<- subset(d, !(DELTA14Cage.atmos %in% d2$DELTA14Cage.atmos))
				
				points(d3.1$cal.age, d3.1$DELTA14Cage.atmos, pch=16, cex=1.1, col=adjustcolor('coral1', alpha=0.3))
				
				d3<- subset(d, !(DELTA14Cage.atmos %in% d2.3$DELTA14Cage.atmos))
				points(d3$cal.age, d3$DELTA14Cage.atmos, pch=16, cex=0.8, col=adjustcolor('red3', alpha=0.6))
				} else {
					d2<- d
					r[i,4] <- NA
					r[i,5] <- NA
					r[i,6]	<- length(d$DELTA14Cage.atmos)
				}
				s<- rbind(s, d2) 
				s1<- rbind(s1, d3.1) 
					
		}
		

		points(int_age, r[,1], pch=1,cex=1, col=adjustcolor('grey17',alpha=0.7), lwd=0.8)
		#bin sd
		for(i in 1:length(int_age)){
			lines(c(int_age[i],int_age[i]),c(r[i,4]+r[i,5],r[i,4]-r[i,5]), col=adjustcolor('grey57',alpha=0.9), lty=2)
		}
		#bin se
		for(i in 1:length(int_age)){
			lines(c(int_age[i],int_age[i]),c(r[i,4]+(r[i,5]/sqrt(r[i,6])),r[i,4]-(r[i,5]/sqrt(r[i,6]))), col=adjustcolor('grey17',alpha=0.9), lty=1)
		}
		#bin means excluding outliers
		points(int_age, r[,4], pch=16, col=adjustcolor('grey17',alpha=0.7), cex=1.1)
		
		
		
		
		
		
		
		
		#fit loess - fit to all data with mcmc resample of age and 14C errors + bootstrap 
		
		iterations<- iteration.number
		
		s2<-matrix(ncol=iterations,nrow=length(int_s))
		
		for(it in 1:iterations){ 
		#it<- 1
			dat3<- s
			#bootstrap data
			dat3<- dat3[sample(nrow(dat3), nrow(dat3),replace=TRUE), ]
			# dat3<- na.omit(dat3)
			dat3<- subset(dat3, DELTA14Cage.atmos > -999 )
			
			#mcmc age and 14C erors
			dat3$cal.age <- rnorm(length(dat3$cal.age), dat3$cal.age, dat3$cal.age.err)
			#plot(dat3$cal.age, dat3$cal.age.err)
			dat3$DELTA14Cage.atmos <- rnorm(length(dat3$DELTA14Cage.atmos), dat3$DELTA14Cage.atmos,dat3$DELTA14Cage.atmos.err)
			#plot(dat3$proxy.14C.age, dat3$proxy.14c.age.err)
			#dat3$DELTA14Cage.atmos<-(((exp(-dat3$proxy.14C.age/8033))/(exp(-(dat3$cal.age)/8266)))-1)*1000
			#plot(s$DELTA14Cage.atmos, dat3$DELTA14Cage.atmos)
			#plot(dat3$cal.age, dat3$DELTA14Cage.atmos)

	

		
			alpha <- bestLoess(loess(DELTA14Cage.atmos~cal.age, data=dat3))$minimum
	
			l<-loess(DELTA14Cage.atmos~cal.age, span=alpha, data=dat3)
			#l<-loess(D14Cmean~int_age, span=alpha, data=d2,weights=n)
			#p<-predict(l, newdat=data.frame(cal.age= int_age), se=FALSE)
			p<-predict(l, newdat=data.frame(cal.age= int_age), se=TRUE)
			if(it %in% round(runif((iterations+1)/10 , 1, iterations)) == TRUE){
				lines(int_age, p$fit, lty=1, col=adjustcolor('grey57',alpha=0.1), lwd=0.5)} else {NA}
		 	s2[,it] <- rnorm(length(p$fit), p$fit, p$se.fit)
		 	#s2[,it] <- p$fit
		}

		#loess quantiles 
		loess_upr95<- apply(s2, 1, quantile, probs=0.95, na.rm=TRUE) #95% range
		loess_lwr95<- apply(s2, 1, quantile, probs=0.05, na.rm=TRUE) #95% range
		loess_upr68<- apply(s2, 1, quantile, probs=0.68, na.rm=TRUE) #67% range
		loess_lwr68<- apply(s2, 1, quantile, probs=0.32, na.rm=TRUE) #67% range
		loess_fit<- apply(s2, 1, quantile, probs=0.5, na.rm=TRUE) #median value

		lines(int_age, loess_fit, lwd=1.3, col='black') 
		lines(int_age, loess_upr95, lty=1, col=adjustcolor('darkorange',alpha=0.6), lwd=0.6)
		lines(int_age, loess_lwr95, lty=1, col=adjustcolor('darkorange',alpha=0.6), lwd=0.6)
		lines(int_age, loess_upr68, lty=1, col=adjustcolor('firebrick2',alpha=0.6), lwd=0.6)
		lines(int_age, loess_lwr68, lty=1, col=adjustcolor('firebrick2',alpha=0.6), lwd=0.6)
		
#legends
		if(j == 1 & k == 1) {legend('topleft',pch=c(NA,16,16),lty=c(1,NA,NA),col=c('seagreen','grey57','coral1'),c('INTCAL20','all','3s outliers'), cex=0.8, bty='n')}
		
		if(j == 2 & k == 1) {legend('topleft',pch=c(16,1, NA,NA) ,lty=c(NA,NA,1,2),col=c('grey17','grey17','grey17','grey57','firebrick2'),c('bin - 3s outliers removed','bin - all','bin sd','bin se'), cex=0.8, bty='n')}
		
			if(j == 3 & k == 1) {legend('topleft',pch=c(NA,NA,NA),lty=c(1,1,1),col=c('black','firebrick2','darkorange'),c('LOESS 50th percentile','LOESS 68%', 'LOESS 95%'), cex=0.8, bty='n')}
		
		text(6000,-400,as.character(basin), cex=1.1, font=2); text(6000,-525,as.character(depth), cex=1.1)


		#export data
		out<-cbind(int_age, r, loess_fit, loess_upr95, loess_lwr95, loess_upr68, loess_lwr68)
		colnames(out)<- c('int_age','mean_all','sd_all','n_all','mean','sd','n','loess_fit', 'loess_upr95', 'loess_lwr95', 'loess_upr68', 'loess_lwr68')
# add standard error to output dataframe
		out <- as.data.frame(out)
		se <- out$sd/sqrt(out$n)
		out.new <- data.frame(out,loess_se = se)


#IF MAKING A LONG RUN, USE THESE TO OUTPUT THE RESULTS
		write.csv(out.new, paste('csv-DELTA14Cage-atmos_500bin_',as.character(basin),'_',as.character(depth),'bootstrap.csv',sep=''), row.names=FALSE)
		
		write.csv(s1, paste('csv-DELTA14Cage-atmos_500bin_',as.character(basin),'_',as.character(depth),'_outliersbootstrap.csv',sep=''), row.names=FALSE)
		
		write.csv(s1, paste('csv-14Cbin_outliers',as.character(basin),'_',as.character(depth),'bootstrap.csv',sep=''), row.names=FALSE)
		
		
		
	}
}

dev.off()







