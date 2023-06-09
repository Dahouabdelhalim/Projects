


#set working directory
setwd('/Users/patrickrafter/Documents/R/R/global14C/r_test_for_Zenmodo_and_Github')

#import data
dat<- read.csv('/Users/patrickrafter/Documents/R/R/global14C/r_test_for_Zenmodo_and_Github/rafter-2022-Global-14C-Compilation-FIN.csv', skip=1) #This was the file name we used, but you will have to re-name what you've downloaded from PANGAEA here: https://doi.pangaea.de/10.1594/
colnames(dat)

intcal20<-read.table("/Users/patrickrafter/Documents/R/R/intcal/14C-intcal20.14c.txt")
colnames(intcal20)<- c("age","14C_age","error","D14C","sigma") 

etc<- read.csv('/Users/patrickrafter/Documents/R/R/14C-etc/14c-global-age-density-etc.csv')
colnames(etc)



pdf(paste('plot_2-AGE-plotting-loess.pdf',sep=''), width=2.5, height=5, useDingbats=FALSE, encoding="WinAnsi") #make one big plot

# ONLY RUN SECOND AFTER global14C_bin.R

if (!require(dplyr)) install.packages('dplyr')
require(dplyr)
library(dplyr)




###################
#simple multipanel plot

#dev.new(width=7.2, height=4.8)
par(mfrow=c(2,1)) #use mfrow to define simple matrix of panels to fill - here 4 cols, 3 rows
par(mar=c(3,3,0.5,0.5),xpd=FALSE); par(ps = 9, cex = 0.7, cex.main = 1); par(mgp=c(1.5,0.5,0)); par(las=1); par(tck=0.02)



basins<- c('ALL','Atlantic','Pacific','SoSubtrop','SO')
depths<- c('ALL','deep','mid')

ymin<- -300
ymax<- 4000
xmin<- 0
xmax<- 25000
yshade<-ymax-500
text1<-ymin+450
text2<-ymin+400

cutoffdepth<- 2000

ydend<- 11700; ydstart<- 12900; hs1end1<- 14800; hs1start1<- 16200; LGM<-20500; hs2end<-23300; hs2start<-26500

#####
#plot
#import binned data
nad<- read.csv('csv-DELTA14Cage-atmos_500bin_Atlantic_deepbootstrap.csv')
nad<- filter(nad, int_age <25000)
nai<- read.csv('csv-DELTA14Cage-atmos_500bin_Atlantic_midbootstrap.csv')
nai<- filter(nai, int_age <25000)

npd<- read.csv('csv-DELTA14Cage-atmos_500bin_Pacific_deepbootstrap.csv')
npd<- filter(npd, int_age <25000)
npi<- read.csv('csv-DELTA14Cage-atmos_500bin_Pacific_midbootstrap.csv')
npi<- filter(npi, int_age <25000)

sod<- read.csv('csv-DELTA14Cage-atmos_500bin_SO_deepbootstrap.csv')
sod<- filter(sod, int_age <25000)
soi<- read.csv('csv-DELTA14Cage-atmos_500bin_SO_midbootstrap.csv')
soi<- filter(soi, int_age < 25000)



# cols		1		2 Atl		3 PO			4			5 SO		6			7 PO	inner	8 ATLinner	9 SO inner
cols<- c('black','goldenrod','#6e8b38',		'seagreen',	'#cdaa7d',	'tan4',		'#006400',	'#8b6508',	'#8b3e2f')











#FIRST ROW
#FIRST ROW FIRST PLOT

#All mid Depth Averages

plot(-999,-999,pch=1,col=adjustcolor("grey99", alpha=0.01),xlim=c(0,25000),ylim=rev(c(-125,3500)), ylab=expression(paste(proxyatmos^14,Cage, sep='')), xlab='calendar age (years BP)')

axis(3,labels=FALSE,tick=TRUE, line=FALSE);axis(4,labels=FALSE,tick=TRUE, line=FALSE)
#grid (NULL,NULL, lty = 6, col = adjustcolor("cornsilk2", alpha=0.7))
rect(ydend,ymin,ydstart,ymax,col=adjustcolor("dimgrey",alpha=0.2),border=FALSE)
rect(hs1end1,ymin, hs1start1,ymax,col=adjustcolor("dimgrey",alpha=0.2),border=FALSE)
rect(hs1end1,ymin, hs1start1 +1750,ymax,col=adjustcolor("dimgrey",alpha=0.2),border=FALSE)
rect(hs2end,ymin, hs2start,ymax,col=adjustcolor("dimgrey",alpha=0.2),border=FALSE)



#southern ocean
soi<- soi[is.na(soi$loess_fit)==FALSE,]
#soi<- filter(soi, soi$int_age < 21000)
polygon(x=c(soi$int_age,rev(soi$int_age)),y=c(soi$loess_upr95,rev(soi$loess_lwr95)),col=adjustcolor(cols[5],alpha=0.5),border=NA)
polygon(x=c(soi$int_age,rev(soi$int_age)),y=c(soi$loess_upr68,rev(soi$loess_lwr68)),col=adjustcolor(cols[5],alpha=0.7),border=NA)
lines(soi$int_age, soi$loess_fit, col=adjustcolor(cols[9], alpha=1), lwd=1.3, lty=c("11"))


#Pacific
npi<- npi[is.na(npi$loess_fit)==FALSE,]

polygon(x=c(npi$int_age,rev(npi$int_age)),y=c(npi$loess_upr95,rev(npi$loess_lwr95)),col = adjustcolor(cols[3],alpha=0.5),border=NA)
polygon(x=c(npi$int_age,rev(npi$int_age)),y=c(npi$loess_upr68,rev(npi$loess_lwr68)),col = adjustcolor(cols[3],alpha=0.7),border=NA)

lines(npi$int_age, npi$loess_fit, col=adjustcolor(cols[7], alpha=1), lwd=1.3, lty=c("11"))


#Atlantic
nai<- nai[is.na(nai$loess_fit)==FALSE,]
#nai<- filter(nai, int_age >2000)
#Atlantic 2 sigma envelope
polygon(x=c(nai$int_age,rev(nai$int_age)),y=c(nai$loess_upr95,rev(nai$loess_lwr95)),col = adjustcolor(cols[2],alpha=0.5),border=NA)
polygon(x=c(nai$int_age,rev(nai$int_age)),y=c(nai$loess_upr68,rev(nai$loess_lwr68)),col = adjustcolor(cols[2],alpha=0.7),border=NA)
lines(nai$int_age, nai$loess_fit, col=adjustcolor(cols[8], alpha=1), lwd=1.3, lty=c("11"))

#accoutrement
#modern error bars
lines(c(etc$age[2],etc$age[2]),c(etc$atlmid[2]+etc$atlmiderr[2],etc$atlmid[2]-etc$atlmiderr[2]), col=adjustcolor(cols[2],alpha=1))
lines(c(etc$age[2],etc$age[2]),c(etc$somid[2]+etc$somiderr[2],etc$somid[2]-etc$somiderr[2]), col=adjustcolor(cols[9],alpha=1))
lines(c(etc$age[2],etc$age[2]),c(etc$pacmid[2]+etc$pacmiderr[2],etc$pacmid[2]-etc$pacmiderr[2]), col=adjustcolor(cols[7],alpha=1))

points(etc$age[2],etc$atlmid[2], pch=0, cex=2, col=(cols[2]))
points(etc$age[2],etc$atlmid[2], pch=15, cex=1.8, col='white')

points(etc$age[2],etc$somid[2], pch=1, cex=2, col=(cols[9]))
points(etc$age[2],etc$somid[2], pch=16, cex=1.8, col='white')

points(etc$age[2],etc$pacmid[2], pch=5, cex=1.6, col=(cols[7]))
points(etc$age[2],etc$pacmid[2], pch=18, cex=2, col='white')



text(12000,text1, 'Pacific,Atlantic,&Southern 27.5-28 density 40S-60N')
#text(12000,text2, cex=0.7, '40S-60N')
text(mean(c(ydend, 	ydstart)), 			yshade, cex=0.7, 'YD')
text(mean(c(hs1end1,ydstart)), 			yshade, cex=0.7, 'BA')
text(LGM,								yshade, cex=0.7, 'LGM')
text(mean(c(hs1end1,hs1start1 +1750)),	yshade, cex=0.7, 'HS1')
text(mean(c(hs2end,	hs2start)),			yshade, cex=0.7, 'HS2')

{legend('bottomleft',lty=c(3,3,3),col=c(cols[2],cols[9],cols[7]),c('Atlantic','Southern','Pacific'),cex=0.9, bty='n')}




















#second row

#SECOND ROW FIRST PLOT
#All Deep Depth Averages

plot(-999,-999,pch=1,col=adjustcolor("grey99", alpha=0.01),xlim=c(0,25000),ylim=rev(c(-125,3500)), ylab=expression(paste(proxyatmos^14,Cage, sep='')),
xlab='calendar age (years BP)')

axis(3,labels=FALSE,tick=TRUE, line=FALSE);axis(4,labels=FALSE,tick=TRUE, line=FALSE)
#grid (NULL,NULL, lty = 6, col = adjustcolor("cornsilk2", alpha=0.7))
rect(ydend,ymin,ydstart,ymax,col=adjustcolor("dimgrey",alpha=0.2),border=FALSE)
rect(hs1end1,ymin, hs1start1,ymax,col=adjustcolor("dimgrey",alpha=0.2),border=FALSE)
rect(hs1end1,ymin, hs1start1 +1750,ymax,col=adjustcolor("dimgrey",alpha=0.2),border=FALSE)
rect(hs2end,ymin, hs2start,ymax,col=adjustcolor("dimgrey",alpha=0.2),border=FALSE)

#modern value lines
#lines(c(xmin,xmax),c(etc$atldeep[2],etc$atldeep[2]), col=adjustcolor('brown1',alpha=1))
#lines(c(xmin,xmax),c(etc$sodeep[2],etc$sodeep[2]), col=adjustcolor('salmon4',alpha=1))
#lines(c(xmin,xmax),c(etc$pacdeep[2],etc$pacdeep[2]), col=adjustcolor('dodgerblue4',alpha=1))


#SOUTHERN OCEAN
sod<- sod[is.na(sod$loess_fit)==FALSE,]
#sod<- filter(sod, int_age > 2000)
polygon(x=c(sod$int_age,rev(sod$int_age)),y=c(sod$loess_upr95,rev(sod$loess_lwr95)),col=adjustcolor(cols[5],alpha=0.5),border=NA)
polygon(x=c(sod$int_age,rev(sod$int_age)),y=c(sod$loess_upr68,rev(sod$loess_lwr68)),col = adjustcolor(cols[5],alpha=0.7),border=NA)
lines(sod$int_age, sod$loess_fit, col=adjustcolor(cols[9], alpha=1), lwd=1.3, lty=1)


#PACIFIC
npd<- npd[is.na(npd$loess_fit)==FALSE,]
polygon(x=c(npd$int_age,rev(npd$int_age)),y=c(npd$loess_upr95,rev(npd$loess_lwr95)),col=adjustcolor(cols[3],alpha=0.5),border=NA)
polygon(x=c(npd$int_age,rev(npd$int_age)),y=c(npd$loess_upr68,rev(npd$loess_lwr68)),col=adjustcolor(cols[3],alpha=0.7),border=NA)
lines(npd$int_age, npd$loess_fit, col=adjustcolor(cols[7], alpha=1), lwd=1.3, lty=1)



#ATLANTIC
nad<- nad[is.na(nad$loess_fit)==FALSE,]
polygon(x=c(nad$int_age,rev(nad$int_age)),y=c(nad$loess_upr95,rev(nad$loess_lwr95)),col=adjustcolor(cols[2],alpha=0.5),border=NA)
polygon(x=c(nad$int_age,rev(nad$int_age)),y=c(nad$loess_upr68,rev(nad$loess_lwr68)),col=adjustcolor(cols[2],alpha=0.7),border=NA)
#polygon(x=c(nad$int_age,rev(nad$int_age)),y=c(nad$loess_fit+nad$loess_se,rev(nad$loess_fit-nad$loess_se)),col = adjustcolor(cols[2],alpha=0.7),border=NA)
lines(nad$int_age, nad$loess_fit, col=adjustcolor(cols[8], alpha=1), lwd=1.3, lty=1)

#errors
lines(c(etc$age[2],etc$age[2]),c(etc$atldeep[2]+etc$atldeeperr[2],etc$atldeep[2]-etc$atldeeperr[2]), col=adjustcolor(cols[2],alpha=0.7))
lines(c(etc$age[2],etc$age[2]),c(etc$pacdeep[2]+etc$pacdeeperr[2],etc$pacdeep[2]-etc$pacdeeperr[2]), col=adjustcolor(cols[7],alpha=0.7))
lines(c(etc$age[2],etc$age[2]),c(etc$sodeep[2]+etc$sodeeperr[2],etc$sodeep[2]-etc$sodeeperr[2]), col=adjustcolor(cols[9],alpha=0.7))

points(etc$age[2],etc$atldeep[2], pch=15, cex=1.9, col=cols[2])
points(etc$age[2],etc$atldeep[2], pch=0, cex=2, col='black')

points(etc$age[2],etc$sodeep[2], pch=16, cex=1.9, col=cols[9])
points(etc$age[2],etc$sodeep[2], pch=21, cex=2, col='black')

points(etc$age[2],etc$pacdeep[2], pch=18, cex=2, col=cols[7])
points(etc$age[2],etc$pacdeep[2], pch=23, cex=1.8, col='black')



text(12000,text1, 'Pacific,Atlantic,&Southern density >28')
text(12000,text2, cex=0.7, '')
text(mean(c(ydend, 		ydstart)), 			yshade, cex=0.7, 'YD')
text(mean(c(hs1end1, 	ydstart)), 			yshade, cex=0.7, 'BA')
text(mean(c(hs1end1, 	hs1start1 +1750)),	yshade, cex=0.7, 'HS1')
text(LGM,						 			yshade, cex=0.7, 'LGM')
text(mean(c(hs2end, 	hs2start)),			yshade, cex=0.7, 'HS2')

{legend('bottomleft',lty=c(1,1,1),col=c(cols[2],cols[9],cols[7]),c('Atlantic','Southern','Pacific'),cex=0.9, bty='n')}










##########
dev.off() # this closes/saves the big plot device
#############






