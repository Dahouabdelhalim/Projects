#R code written by Will Gray and Patrick Rafter (prafter@uci.edu) for updating sediment core calendar age estimates based on planktic foraminifera 14C/C -- here calibrated to Intcal20 via Marine20 using the Bacon package in R


setwd('/Users/patrickrafter/Documents/R/R/global14C/r_test_for_Zenmodo_and_Github/rafter-2022-Global-14C-Compilation-FIN.csv') # set your working directory
#Data needs to be downloaded from PANGAEA here: https://doi.pangaea.de/10.1594/PANGAEA.946522
require(rbacon)
require(Hmisc)

dat_all<- read.csv('IN-prafter-2020-Fossil-D14C-Bacon-R-4.0.csv') #import your data
colnames(dat_all)

dat<- subset(dat_all, chro.meth == 3) #subset plank 14C chronologies
stationname<- unique(dat$station) #get list of station names

#pull out depth/age data and create file structure for Bacon
dir.create('Cores/') #create Cores folder first time running code
for(i in 1:length(stationname)){
	#i=1
	d<- subset(dat, station == stationname[i])
	labID<- as.character(paste(stationname[i],'_',seq(nrow(d)), sep=''))
	age<- d$planktic.14C.age
	error<- d$planktic.14C.age.err
	depth<- if(d$depth.dated.mid[1] > -999) {d$depth.dated.mid} else {rowMeans(cbind(d$depth.dated.top, d$depth.dated.bottom))} #mid depth calculator
	dR<- d$Delta.R
	dR_err<- d$Delta.R.err
	cc<-d$cc
	out<- data.frame(labID=labID, age=age, error=error, depth=depth, cc=cc, delta.R=dR,d.STD=dR_err)
	out<- out[(age > -999),] #remove depths with no plank 14C
	out<- out[order(depth),]
	out<- na.omit(out) #remove any NA rows
	dir.create(paste('Cores/Station',as.character(stationname[i]),sep='')) #create Station folder and write data to csv
	write.csv(out, paste('Cores/Station',as.character(stationname[i]),'/Station',as.character(stationname[i]),'.csv',sep=''),row.names=FALSE)
	}
	


#BACON those dudes (this will take awhile if you are running the full dataset)
 for(i in 1:length(stationname)){
	 #i=1
	 corename<- paste('Station',as.character(stationname[i]),sep='')
	 #first if statement skips cores with  bacon'd age models already done - comment out to run all
	if(length(list.files(paste('Cores/Station',as.character(stationname[i]),sep=''), pattern='_ages.txt')) < 1) {
		 d<- read.csv(paste('Cores/Station',as.character(stationname[i]),'/Station',as.character(stationname[i]),'.csv',sep=''))
		 d<- na.omit(d)
		 if(nrow(d) > 1){ #make sure more than one date 
			 core.length<- max(d$depth) - min(d$depth)
			 sections<- if (core.length < 25) {1} else  #sections determines n intervals
				 if(core.length >= 25 & core.length < 50){2} else 
				 if(core.length >= 50 & core.length < 250){5} else 
				 if(core.length >= 250 & core.length < 500){10} else 
				 if(core.length >= 500 & core.length < 1000){25} else 
				 if(core.length >= 1000){50}
			 approx14Cyrs<- (max(d$age) - min(d$age))
			 res.approx<- approx14Cyrs/core.length #gives approx sed rate as prior
			 Bacon(core=corename, thick=sections, ask=FALSE, suggest=FALSE, acc.mean=res.approx)
			 # Bacon(core=corename, cc=2, thick=sections, ask=FALSE, suggest=FALSE, acc.mean=res.approx)
			 }
		 } #comment out to run all
	 }


# #make list of cores with no Bacon age model (i.e. single dates)
 s<-vector(mode='numeric',length=0)
 for(i in 1:length(stationname)){
	 #i=1
	 corename<- paste('Station',as.character(stationname[i]),sep='')
	 if(length(list.files(paste('Cores/Station',as.character(stationname[i]),sep=''), pattern='_ages.txt')) < 1) {
		 s<- rbind(s,corename)
		 }
	 }
	
 write.csv(data.frame(site.name=s), 'cores_that_did_not_bacon.csv', row.names=FALSE )


# #output age model data to file with benthic 14C data
dat_all<- read.csv('IN-prafter-2020-Fossil-D14C-Bacon-R-4.0.csv') #necessary?
colnames(dat_all)
dat<-dat_all
stationname<- unique(dat$station) #get list of station names

dat$cal.age.m20.median<- NA
dat$cal.age.m20.min<- NA
dat$cal.age.m20.max<- NA

for(i in 1:length(stationname)){
	 #i=2
	 corename<- paste('Station',as.character(stationname[i]),sep='')
	 d<- subset(dat, station == stationname[i])
	 rownumbers<- as.numeric(rownames(d))
	 if(d$chro.meth[1] == 3){
		 if(length(list.files(paste('Cores/Station',as.character(stationname[i]),sep=''), pattern='_ages.txt')) >= 1) {
			 depth_out<- if(d$depth.dated.mid[1] > -999) {d$depth.dated.mid} else {rowMeans(cbind(d$depth.dated.top, d$depth.dated.bottom))} 
			 age_file<- list.files(paste('Cores/Station',as.character(stationname[i]),sep=''), pattern='_ages.txt')
			 ages<- read.table(paste('Cores/Station',as.character(stationname[i]),'/',age_file,sep=''), header=TRUE)
			 #ages<- rbind(c(0,-50,ages[1,3],mean(c(-50, ages[1,3])),mean(c(-50, ages[1,3]))), ages) 
			 if(length(ages$depth) > 1) {
				 #median<- approx(ages$depth, ages$median, xout=depth_out)$y
				 #min<- approx(ages$depth, ages$min, xout=depth_out)$y
				 #max<- approx(ages$depth, ages$min, xout=depth_out)$y
				 median<- approxExtrap(ages$depth, ages$median, xout=depth_out)$y
				 min<- approxExtrap(ages$depth, ages$min, xout=depth_out)$y
				 max<- approxExtrap(ages$depth, ages$max, xout=depth_out)$y
				 if(median[1] > -50) {
					 dat$cal.age.m20.median[rownumbers]<- median
					 dat$cal.age.m20.min[rownumbers]<- min
					 dat$cal.age.m20.max[rownumbers]<- max
					 } else if(median[1] < -50) {
						 ages<- rbind(c(0,-50,50,0,0), ages) #sets 0 depth to 0 if extrpolation leads to ages 							<-50 
						 #ages<- rbind(c(0,-50,ages[1,3],0,mean(c(-50, ages[1,3]))), ages) 
						 median<- approxExtrap(ages$depth, ages$median, xout=depth_out)$y
						 min<- approxExtrap(ages$depth, ages$min, xout=depth_out)$y
						 max<- approxExtrap(ages$depth, ages$max, xout=depth_out)$y 
						 dat$cal.age.m20.median[rownumbers]<- median
						 dat$cal.age.m20.min[rownumbers]<- min
						 dat$cal.age.m20.max[rownumbers]<- max
						 #plot(ages$median, ages$depth, xlim=c(-50, 40000)); points(median, depth_out, col='red')
						 }				
				 }
			 }
		 }
	 }
	
 write.csv(dat,'OUT-prafter-2020-Fossil-D14C-Bacon-R-4.0.csv', row.names=FALSE)







