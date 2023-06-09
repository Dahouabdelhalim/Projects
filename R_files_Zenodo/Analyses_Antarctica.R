# Analyses_Antarctica.R - Script to replicate all analyses from Molina et al. (2022) Heat tolerance of marine ectotherms in a warming Antarctica. Global Change Biology. 


  #### load data #### 

   library(openxlsx)
   library(tidyverse)
   library(ape)
   library(phytools)
   library(nlme)
   library(qpcR)
   library(emmeans)	
	
	
  # Opening dataset removing tunicate "Cnemidocarpa verrucosa" 
  # Removing single values for phyla Annelida, Nemertea, Cnidaria 
    
	  df <- read.xlsx("Dataset_Antarctic_thermal_tolerance.xlsx")[,-c(22,23)]
	  df <- df[-which(df$sp=="Cnemidocarpa verrucosa"),]
	  df <- df[-which(df$phylum=="Cnidaria" | df$phylum=="Annelida" | df$phylum=="Nemertea"),]
 	  df$sp.id <- as.factor(paste(df$sp,with(df, ave(seq_along(sp), sp, FUN=seq_along))))
 	  df$logt <- log10(df$t)
	
	 
  # Building phylogeny based on taxonomic data
  # Sorting df to match the order of the phylogeny
   
	  taxa <- ~bil.cni/prot.deut/phylum/class/order/family/genus/sp/sp.id
	  df[,1:9] <- lapply(df[,1:9], factor)
	  tr.df <- as.phylo.formula(taxa,data= df)
	  tr.df <- compute.brlen(tr.df)
	  df <- df[match(tr.df$tip.label[tr.df$edge[tr.df$edge[, 2] <= Ntip(tr.df), 2]],df$sp.id),]

	  
  # Merging data and including index for repeated measurements per species 
 	  

	  d1 <- gls(ctmax ~ logt + te + mort,data=df,correlation=corPagel(0,tr.df,fixed=TRUE,form=~paste(sp.id)),method="ML")	
	  d2 <- gls(ctmax ~ logt + te + mort,data=df,correlation=corPagel(1,tr.df,fixed=TRUE,form=~paste(sp.id)),method="ML")	
	  d3 <- gls(ctmax ~ logt + phylum + te + mort,data=df,correlation=corPagel(0,tr.df,fixed=TRUE,form=~paste(sp.id)),method="ML")	
	  d4 <- gls(ctmax ~ logt + phylum + te + mort,data=df,correlation=corPagel(1,tr.df,fixed=TRUE,form=~paste(sp.id)),method="ML")	
	  d5 <- gls(ctmax ~ logt * phylum + te + mort,data=df,correlation=corPagel(0,tr.df,fixed=TRUE,form=~paste(sp.id)),method="ML")	
	  d6 <- gls(ctmax ~ logt * phylum + te + mort,data=df,correlation=corPagel(1,tr.df,fixed=TRUE,form=~paste(sp.id)),method="ML")	

   	  aic <- AIC(d1,d2,d3,d4,d5,d6)
	  aic <- cbind(aic,akaike.weights(aic[,2]))
	  aic
  	
  # Obtaining z from mean adjusted estimates from model 6 (different z per phylum + phylognetic signal)
  	
	  z <- data.frame(emtrends(d6,pairwise~phylum,var="logt")$emtrends)
	  z


  # Recoding z per species to calculate tko = 1 min and 30 days	
  # Three steps
  # 1. Obtain z per species
  # 2. Calculate the standardized time assuming a static assay
  # 3. Interpolating Tko at 1 min and 30 days employing Tko and std time
  	
		df$z <- recode(df$phylum, "Arthropoda" = -z[1,2],
							      "Brachiopoda" = -z[2,2],
							      "Chordata" = -z[3,2], 
							      "Echinodermata" = -z[4,2],
							      "Mollusca"= -z[5,2])


   # Standard time adjusted
  
    	std.time.max <- function(z,t0,ctmax,ramp){ta <- seq(t0, ctmax,by=ramp); sum(1/10^((ctmax-ta[-length(ta)])/z))}	
		df$timestd <- df$t				
		for(i in which(df$ctmax != 'NA')){
    	if(df$te[i] == "dynamic"){df$timestd[i] <- std.time.max(df$z[i], df$t0[i], df$ctmax[i], df$ramp[i])}}


   # Tko for 1 min and 30 days
   
	    df$tko.1min <- df$ctmax - df$z*log10(1/df$timestd)
   		df$tko.30d <- df$ctmax - df$z*log10(30*1440/df$timestd)

	
   # Putting standardization to a test (assay effect in d6 is no longer significant as expected)
   
		df$logtstd <- log10(df$timestd)
		d7 <- gls(ctmax ~ logtstd * phylum + te + mort,data=df,correlation=corPagel(1,tr.df,fixed=TRUE,form=~paste(sp.id)),method="ML")	
		summary(d7)

	
	# Phylogenetic signal in tko.1min and tko.30d
		
		gls(tko.1min ~ phylum,data=df,correlation=corPagel(1,tr.df,form=~paste(sp.id)),method="ML")
		gls(tko.30d ~ phylum,data=df,correlation=corPagel(1,tr.df,form=~paste(sp.id)),method="ML")
		
		

# ----------------------------- FIG. 2 -----------------------------

	# New dataset collapsing values per species

		df.sp <- df[duplicated(df$sp)==FALSE,1:8]
		df.sp <- df.sp[order(df.sp$sp),]
		df.sp$nsp <- table(df$sp)
		df.sp$nsp <- paste(df.sp$sp," (",df.sp$nsp,")",sep="")
		df.sp$z <- tapply(df$z,df$sp,mean)
		df.sp$tko.1min <- tapply(df$tko.1min,df$sp,mean)
		df.sp$tko.30d <- tapply(df$tko.30d,df$sp,mean)	
		df.sp$col <- paste(recode(df.sp$phylum, "Arthropoda" = "darkred",
							      "Brachiopoda" = "grey62",
							      "Chordata" = "royalblue", 
							      "Echinodermata" = "mediumaquamarine",
							      "Mollusca"= "indianred1"))
		df$nsp <- df.sp$nsp[match(df$sp,df.sp$sp)]

		
	# Phylogeny at the species level
	
		taxa <- ~bil.cni/prot.deut/phylum/class/order/family/genus/sp
 		tr.df <- as.phylo.formula(taxa,data= df.sp)
	  	df.sp <- df.sp[match(tr.df$tip.label,df.sp$sp),] 		
 		tr.df$tip.label <- df.sp$nsp
 		tr.df <- ladderize(tr.df)
		order.sp <- tr.df$tip.label[tr.df$edge[tr.df$edge[, 2] <= Ntip(tr.df), 2]]			# order required for loop below
				
		quartz(,8.5,6.5)
		par(mar=c(4,0,1.5,0.1))  
		layout(matrix(c(1,2,3),1,3),widths=c(2,0.4,0.5))
		plot(tr.df,show.tip.label=T,cex=1.2, label.offset = 0.01) 


	# Range of tko between 1 min and 30 days including empirical data (loop)		
	
		par(xpd=NA)  
 		plot(df.sp$tko.30d,1:length(df.sp$sp),pch=21,bg=NA,col="white",
	       cex=1.5,xlim=c(-5,45),yaxt="n",ylab="", cex.lab=1.2,
   		   xlab=expression(paste("Standardized T"[ko]*" ("*degree*"C)")),bty="n")
	 	j <- tr.df$edge[tr.df$edge[, 2] <= Ntip(tr.df), 2]
	 	segments(df.sp$tko.30d[j], 1:length(df.sp$sp), df.sp$tko.1min[j],1:length(df.sp$sp), col=df.sp$col[j],lwd=6)

		arrows(46,36,5,36,length=0.05)
		text(0,36.7,"30 days",cex=1.1)
		text(48,36.7,"1 min",cex=1.1)	
		
			for(i in 1:length(order.sp)){
				ctmax <- df$ctmax[which(df$nsp==order.sp[i])]
				pch <- ifelse(df$te[which(df$nsp==order.sp[i])]=="dynamic",21,24)
				points(ctmax,rep(i,length(ctmax)),pch=pch,cex=1.5,bg=df.sp$col[j][i])}

	# Captions including calculated z from emmeans(d6) above
			
		plot(df.sp$tko.30d,1:length(df.sp$sp),pch=21,bg=NA,col="white",
	       cex=1.5,xlim=c(-5,5),yaxt="n",xlab="",ylab="",bty="n", axes=FALSE)
			text(-4,35,"Brachiopoda",cex=1.2, col="grey62",adj=0)
			text(-4,34,paste("z =",round(-z[2,2],2),"\\u00B1",round(z[2,3],2),"\\u00BAC"),cex=1.2,col="grey62",adj=0) 
			text(-4,31.4,"Mollusca",cex=1.2, col="indianred1",adj=0)
			text(-4,30.4,paste("z =",round(-z[5,2],2),"\\u00B1",round(z[5,3],2),"\\u00BAC"),cex=1.2,col="indianred1",adj=0) 
			text(-4,25.5,"Arthropoda",cex=1.2, col="darkred",adj=0)
			text(-4,24.5,paste("z =",round(-z[1,2],2),"\\u00B1",round(z[1,3],2),"\\u00BAC"),cex=1.2,col="darkred",adj=0) 
			text(-4,19.5,"Echinodermata",cex=1.2, col="mediumaquamarine",adj=0)
			text(-4,18.5,paste("z =",round(-z[4,2],2),"\\u00B1",round(z[4,3],2),"\\u00BAC"),cex=1.2,col="mediumaquamarine",adj=0) 
			text(-4,9.5,"Chordata",cex=1.2, col="royalblue1",adj=0)		
			text(-4,8.5,paste("z =",round(-z[3,2],2),"\\u00B1",round(z[3,3],2),"\\u00BAC"),cex=1.2,col="royalblue1",adj=0) 

			legend(-3,4.7,legend=c("Static"),pch=c(24),col = c("black"), pt.bg = c("black"), cex=1.2, bty = "n")
			legend(-3,3.7,legend=c("Dynamic"),pch=c(21),col = c("black"), pt.bg = c("black"), cex=1.2, bty = "n")





# ----------------------------- FIG. 3 -----------------------------

	# Comparison of CTmax and z between animals at low latitudes (Rezende et al 2014) and this study
  	
  	# Data from Rezende et al. 2014

	  	rez <- read.xlsx("Dataset_SM_RezendeEA2014.xlsx")
 	 	rez$col <- paste(recode(rez$class, "Insecta" = "darkred",
							      "Actinopterygii" = "royalblue",
							      "Bivalvia"= "indianred1"))

		  quartz(,10.8,5)
		  par(mfrow=c(1,2),mar=c(4,4.5,1,1))		
		  plot(NA,NA, xlim=c(0,10),ylim=c(0,65),ylab=expression(paste("CT"[max]*" ("*degree*"C)")),xlab=expression(paste("z ("*degree*"C)")),cex.lab=1.2,cex.axis=0.9,cex.main=2,xaxt="n",xaxs="i",yaxs="i")
	 	  axis(1,at=c(0,2,4,6,8,10),labels = c("0","2","4","6","8", "10"),font=1,las=1,cex.axis=0.9)	
		  points(rez$z, rez$CTmax,cex=2,bg="white", col=rez$col, pch = 22)
		  abline(a = 31, b = 3.2, col = "black", lty=5)
		  
 		  points(df.sp$z, df.sp$tko.1min,col="black",pch=22,bg=df.sp$col,cex=2)
 		  abline(lm(df.sp$tko.1min ~ df.sp$z),lty=2)
 		  
 		 
 	# Legend
 	  
 		  text(0.35,58,"Rezende et al. 2014",cex=0.9,adj=c(0,0))
 		  text(8,35.5,"This study",cex=0.9)
 		  legend(6.5,23,legend=c("Arthropoda","Mollusca", "Chordata","Echinodermata", "Brachiopoda"),
							col=c("darkred","indianred1","royalblue1", "mediumaquamarine","grey62"),
         					pch=15,cex=0.9, bty = "n", text.font = 1)

		plot(c(-1,-1),c(0,0),xlim=c(0,6),ylim=c(0,65),ylab=expression(paste("T"[ko]*" ("*degree*"C)")),
			xlab=expression(paste("log"[10]*"Time")),xaxs="i",yaxs="i",xaxt="n",cex.lab=1.2,cex.axis=0.9)
		for(i in 1:56){
			abline(as.numeric(rez$CTmax[i]),-as.numeric(rez$z[i]), col=rez$col[i],lwd=0.8)}
		polygon(c(0,0,6,6),c(0,80,80,0),col=rgb(1,1,1,0.8))
		
		for(i in 1:35){
		abline(df.sp$tko.1min[i],-df.sp$z[i],col=df.sp$col[i],lwd=1.2)}
		axis(1,at=log10(c(1,60,24*60,24*60*10,24*60*100)),labels=c("1 min","1 hour","1 day","10 days","100 days"),cex.axis=0.9)
	 	text(2.1,53,"Rezende et al. 2014",cex=0.9,adj=c(0,0))
 		text(1,6.5,"This study",cex=0.9)
		segments(log10(24*60*30),0,log10(24*60*30),59,lty=2)
		text(log10(24*60*30),62,"30 days",cex=0.8)





# ----------------------------- FIG. 4 -----------------------------
	
	
	# Open temperature data
	# Open 'Thermal landscape functions.R' described in Rezende et al. (2020) for simulations (last opened on 31/08/2022)
	
		setwd(paste(getwd(),"/Simulated mortality",sep=""))
		source("https://datadryad.org/stash/downloads/file_stream/392706")
	
	
	# Intertidal - Kuklinski & Balazy 2014
	# Dataset available at https://ars.els-cdn.com/content/image/1-s2.0-S1385110113001809-mmc1.xlsx (last opened on 31/08/2022)
	# Interpolate for 1-min resolution 
	
	itemp <- read.xlsx("https://ars.els-cdn.com/content/image/1-s2.0-S1385110113001809-mmc1.xlsx",startRow=4)
	itemp$date <- seq(as.POSIXct("2010-12-07 19:00"), as.POSIXct("2011-03-18 08:45"), by = "5 min")
		
		date <- seq(as.POSIXct("2010-12-07 19:00"), as.POSIXct("2011-03-18 08:49"), by = "1 min")
		time <- spline(itemp$A1,n=length(itemp$date)*5)$x	
		ta.a1 <- spline(itemp$A1,n=length(itemp$date)*5)$y
		ta.b1 <- spline(itemp$B1,n=length(itemp$date)*5)$y
		ta.c1 <- spline(itemp$C1,n=length(itemp$date)*5)$y		
		ta.a2 <- spline(itemp$A2,n=length(itemp$date)*5)$y
		ta.b2 <- spline(itemp$B2,n=length(itemp$date)*5)$y		
		ta.c2 <- spline(itemp$C2,n=length(itemp$date)*5)$y
		ta.a3 <- spline(itemp$A3,n=length(itemp$date)*5)$y
		ta.b3 <- spline(itemp$B3,n=length(itemp$date)*5)$y
		ta.c3 <- spline(itemp$C3,n=length(itemp$date)*5)$y		
		ta.a4 <- spline(itemp$A4,n=length(itemp$date)*5)$y
		ta.b4 <- spline(itemp$B4,n=length(itemp$date)*5)$y		
		ta.c4 <- spline(itemp$C4,n=length(itemp$date)*5)$y

		
	# Subtidal - Cardenas et al. 2018
	# Datasets available at "https://dfzljdn9uc3pi.cloudfront.net/2018/4289/1/PY1_raw_data.csv"
	#						"https://dfzljdn9uc3pi.cloudfront.net/2018/4289/1/PY2_raw_data.csv" (last opened on 31/08/2022)				
	# py1: 10 meters, py2: 20 meters
	# py1: every 1 h; py2: every 2h
	# Values selected from 7 December until March 18 (to agree with intertidal data)*
	# Years changed to 2010 and 2011 to standardize date across datasets 

		py1 <- read.csv("https://dfzljdn9uc3pi.cloudfront.net/2018/4289/1/PY1_raw_data.csv")
		py1$date <- as.POSIXct(py1$Datetime,format = "%Y-%m-%d-%H-%M")
		py1 <- py1[py1$date > "2016-12-07 18:00" & py1$date < "2017-03-18 08:45",]
		date.py1 <- seq(as.POSIXct("2010-12-07 19:00"), as.POSIXct("2011-02-17 19:59:00"), by = "1 min")
		time.py1 <- spline(py1$Temp,n=length(py1$Temp)*60)$x		
		ta.py1 <- spline(py1$Temp,n=length(py1$Temp)*60)$y
	
	
	# py2 remove first 60 min to synchronize everything

		py2 <- read.csv("https://dfzljdn9uc3pi.cloudfront.net/2018/4289/1/PY2_raw_data.csv")
		py2$date <- as.POSIXct(py2$Dt20m,format = "%Y-%m-%d-%H-%M")
		py2 <- py2[py2$date > "2016-12-07 18:00" & py2$date < "2017-03-18 08:45",]
		date.py2 <- seq(as.POSIXct("2010-12-07 18:00"), as.POSIXct("2011-01-24 09:59:00"), by = "1 min")[-c(1:60)]
		time.py2 <- spline(py2$tem20m,n=length(py2$tem20m)*120)$x[-c(1:60)]
		ta.py2 <- spline(py2$tem20m,n=length(py2$tem20m)*120)$y[-c(1:60)]

	# Merging ta datasets
	# Increasing size of ta vector for subtidal with 'NA' 
	# Removing first 300 data from 19:00 to 23:59 to start at 00:00

		ta.py1[(length(time.py1)+1):length(time)] <- NA
		ta.py2[(length(time.py2)+1):length(time)] <- NA		
		date <- date[-c(1:300)]
		ta <- cbind(ta.a1,ta.b1,ta.c1,ta.a2,ta.b2,ta.c2,ta.a3,ta.b3,ta.c3,ta.a4,ta.b4,ta.c4,ta.py1,ta.py2)[-c(1:300),]
		
		
		
	# ---- TDT curves -----		
 
	# Obtain average CTmax value for each taxon
	# phylum_i = 1, 2, 3, 4 and 5 for "Arthropoda","Brachiopoda","Chordata","Echinodermata","Mollusca", respectively
	
		ctmax <- tapply(df.sp$tko.1min,df.sp$phylum,mean)		
 	    ctmax <- data.frame("phylum"=c("Arthropoda","Brachiopoda","Chordata","Echinodermata","Mollusca"),"ctmax"= ctmax)


	# Running thermal tolerance landscape
		
		quartz()
		ind <- 1000
		for(i in 1:5){	
			static <- rep(c(round(ctmax[i,2],2),round(ctmax[i,2]+z[i,2],2), round(ctmax[i,2]+z[i,2]*2,2)),each=ind)
			time <- c(rnorm(ind,1,1/4), rnorm(ind,10,10/4), rnorm(ind,100,100/4))	
			assign(paste0("tl",i),tolerance.landscape(static, time))}
	
	
	
	# ---- Mortality in variable temperatures ----		
	
	# Running variable landscapes per phylum
	# phylum_i = 1, 2, 3, 4 and 5 for "Arthropoda","Brachiopoda","Chordata","Echinodermata","Mollusca", respectively
	# ta_j = ta.a1,ta.b1,ta.c1,ta.a2,ta.b2,ta.c2,ta.a3,ta.b3,ta.c3,ta.a4,ta.b4,ta.c4,ta.py1,ta.py2 (14 columns total)
		
		
	# # Intertidal temperatures with recovery (mortality simulations reset every 24 h) 
		# day <- rep(1:(length(ta[,10])/1440),each=1440)
		# for(i in 1:5){
			# for(j in 1:12){	
				# dyn <- matrix(,1,3)		
					# for(xx in 1:(length(ta[,j])/1440)){
			   				# d <- dynamic.landscape(ta[,j][xx==day], eval(parse(text=paste0("tl",i))))
			   				# dyn <- rbind(dyn,cbind(d$time,d$ta,d$alive))}
				# dyn <- dyn[-1,]
   	     		# dyn[which(is.na(dyn[,3])),3] <- 0
				# dyn <- dyn[seq(0,dim(dyn)[1],60),]
				# assign(paste0("dl",i,"_ta",j),na.omit(dyn))
				# write.table(eval(parse(text=paste0("dl",i,"_ta",j))),file=paste0("dl",i,"_ta",j,".txt"))}}
		
	# # Subtidal temperatures without recovery (mortality throughout the period without resetting)
	# # Simulated present, +0.5 and +1.0ºC.
		# for(i in 1:5){
			# for(j in 13:14){			
			   	# s <- dynamic.landscape(na.omit(ta[,j]), eval(parse(text=paste0("tl",i))))
			   	# s <- cbind(s$time,s$ta,s$alive)
			   	# s <- s[seq(0,dim(s)[1],60),]
			   	# assign(paste0("s",i,"_ta",j),s)
				# write.table(eval(parse(text=paste0("s",i,"_ta",j))),file=paste0("s",i,"_ta",j,"_0.txt"))}}

		# for(i in 1:5){
			# for(j in 13:14){			
			   	# s <- dynamic.landscape(na.omit(ta[,j])+1.0, eval(parse(text=paste0("tl",i))))
			   	# s <- cbind(s$time,s$ta,s$alive)
			   	# s <- s[seq(0,dim(s)[1],60),]
			   	# assign(paste0("s",i,"_ta",j),s)
				# write.table(eval(parse(text=paste0("s",i,"_ta",j))),file=paste0("s",i,"_ta",j,"_1.txt"))}}

		# for(i in 1:5){
			# for(j in 13:14){			
			   	# s <- dynamic.landscape(na.omit(ta[,j])+2.0, eval(parse(text=paste0("tl",i))))
			   	# s <- cbind(s$time,s$ta,s$alive)
			   	# s <- s[seq(0,dim(s)[1],60),]
			   	# assign(paste0("s",i,"_ta",j),s)
				# write.table(eval(parse(text=paste0("s",i,"_ta",j))),file=paste0("s",i,"_ta",j,"_2.txt"))}}


	# ---- Plotting results ----

	
	# Opening files saved above
	# Results in two objects: 'ta' with 12 temperatures data corresponding to logger and 'surv' with survival per phylum/logger	
	
		# Intertidal
	    # Obtaining output 'surv' with survival probability with 1 h resolution
	    
		phy <- c("art","bra","cho","ech","mol")
		logger <- c("a1","b1","c1","a2","b2","c2","a3","b3","c3","a4","b4","c4")
	
			ta <- TRUE; for(j in 1:12){
					ta <- cbind(ta,read.table(file=paste0("dl",1,"_ta",j,".txt"))[,2])}
			ta <- ta[,-1]; colnames(ta) <- logger 
			surv <- TRUE; for(i in 1:5){
					for(j in 1:12){
					surv <- cbind(surv,read.table(file=paste0("dl",i,"_ta",j,".txt"))[,3])}}
			surv <- surv[,-1]; colnames(surv) <- paste(rep(phy,each=length(logger)),rep(logger,length(phy)),sep=".")
			head(ta); head(surv)


		# Plot 1	

				surv <- surv[seq(0,dim(surv)[1],24),]			
				quartz(,11.5,7)		
				layout(matrix(c(1,2,4,5,3,3,6,6),4,2),widths=c(2,1),heights=c(1,0.5,1,0.5))
				par(mar=c(1,6,2,1),oma=c(5,0,2,0),cex.lab=1.4,cex.axis=1.1)
				plot(ta[,1],type="l",ylim=c(-5,22),col=NA,bty="n",xaxt="n",xlab="",ylab="Temperature (ºC)")
					for(i in 1:12){points(ta[,i],type="l",col="dark gray")}			
				par(xpd=NA)
					text(-24,26,"Intertidal",cex=1.8,adj=c(0,0))
					# text(-24,24,"Intermittent acute stress, 24 h recovery",cex=1.4,adj=c(0,0))
					text(24*80,15.5,"12 loggers\\n(2010 - 2011)")
				par(xpd=FALSE)
				
				plot(100 - surv[,1],type="l",ylim=c(0,100), xaxt="n",yaxt="n",yaxs="i",bty="l",ylab="Daily\\nMortality (%)")	
					for(i in 25:36){polygon(c(0:length(surv[,i]),0),c(0,100 - surv[,i],0),col="royalblue",border="royalblue")}
					for(i in 1:12){polygon(c(0:length(surv[,i]),0),c(0,100 - surv[,i],0),col="darkred",border="darkred")}
					for(i in 49:60){polygon(c(0:length(surv[,i]),0),c(0,100 - surv[,i],0),col="indianred1",border="indianred1")}
					for(i in 37:48){polygon(c(0:length(surv[,i]),0)-5,c(0,100 - surv[,i],0),col="brown",border="mediumaquamarine")}			
					for(i in 13:24){polygon(c(0:length(surv[,i]),0),c(0,100 - surv[,i],0),col="gray",border="gray")}
				axis(1,at=c(15,46,74),labels=NA) 
				axis(2,at=c(0,100))
				abline(0,0)
								
				par(mar=c(1,4.5,3,3))
				col <- c("darkred","gray","royalblue","mediumaquamarine","indianred1")	
				boxplot(apply((100 -surv),2,mean,na.rm=TRUE) ~ as.factor(rep(1:5,each=12)),
					log="y",ylim=10^(c(-8,2)),ylab="Mean Daily Mortality (%)",xaxt="n",col=col,border=col)
					axis(1,labels=NA)				
					abline(median(log10(unlist(na.omit(100 -surv)))),0)
					text(1:5,10^c(1.3,-1.5,1.8,-1.5,-1.1),c("a,b","b","a","b","b"),cex=1.2)
					text(4,10^-7,expression(paste(italic("F")["4,55"]*" = 3.49, "*italic("P ")*"= 0.013")),cex=1.3)

				# ANOVA comparing mean mortality rates			
					
					m  <- lm(log10(apply((100 -surv),2,mean,na.rm=TRUE)) ~ as.factor(rep(1:5,each=12)))
					summary(aov(m))				
					TukeyHSD(aov(m))

				# Number of days with mortality above 50%

					days <- TRUE
					for(i in 1:60){
					days[i] <- length(which(surv[,i]<50))}
					tapply(days,rep(1:5,each=12),sum)
					

		# Subtidal
	 		
		phy <- c("art","bra","cho","ech","mol")
		logger <- c("phy1","phy2")
		scenario <- c("0","1","2")	
			 
			 s <- TRUE
			 for(i in 1:5){
					s0 <- read.table(file=paste0("s",i,"_ta",13,"_0.txt"))
					s1 <- read.table(file=paste0("s",i,"_ta",13,"_1.txt"))
					s2 <- read.table(file=paste0("s",i,"_ta",13,"_2.txt"))
					sphy1 <- cbind(s0,s1,s2)
					
					s3 <- read.table(file=paste0("s",i,"_ta",14,"_0.txt"))
					s4 <- read.table(file=paste0("s",i,"_ta",14,"_1.txt"))
					s5 <- read.table(file=paste0("s",i,"_ta",14,"_2.txt"))
					sphy2 <- cbind(s3,s4,s5)
					names(sphy2) <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9")
					sphy2 <- rbind(sphy2,as.data.frame(matrix(NA,nrow(sphy1)-nrow(sphy2),9)))
					sphy12 <- cbind(sphy1,sphy2)	
					s <- cbind(s,sphy12)}
				
			# Whole dataset 	
			s <- s[,-1]	
			colnames(s) <- paste(rep(phy,each=18), rep(rep(logger,each=9),5), rep(rep(scenario,each=3),10), rep(c("time","ta","surv"),10),sep=".")
			
			dates <- seq(as.Date("2010-12-18"), as.Date("2011-03-27"), by="days")	
			ta.s <- cbind(s$art.phy1.0.ta,s$art.phy1.1.ta,s$art.phy1.2.ta,s$art.phy2.0.ta,s$art.phy2.1.ta,s$art.phy2.2.ta) 	
			ta.s <- ta.s[seq(0,dim(ta.s)[1],24),]
			surv.s <- s[,3*(1:30)]
			surv.s <- surv.s[seq(0,dim(surv.s)[1],24),]
				
				par(mar=c(1,6,2,1))
				plot(1:71,ta.s[,1],type="l",ylim=c(-5,23),col=NA,xlim=c(0,100),bty="n",xaxt="n",yaxt="n",xlab="",ylab="Temperature (ºC)")
						for(i in 1:6){points(1:71,ta.s[,i],type="l",col="dark gray",lwd=1.5)}
				par(xpd=NA)
					text(30,-2,"2 loggers\\n(2016 - 2017)")
					text(72,2.2*(1:3),c("current","+1ºC","+2ºC"),adj=c(0,0))
					text(-1,20,"Subtidal",cex=2,adj=c(0,0))
					# text(-1,18,"Persistent chronic stress, no recovery",cex=1.4,adj=c(0,0))
				par(xpd=FALSE)

				axis(2,at=c(-5,0,5,10,15))		
				col <- rep(c("darkred","gray","royalblue","mediumaquamarine","indianred1"),each=6)
				plot(100 - surv.s[,1],type="l",ylim=c(0,20),xlim=c(0,100), xlab="",xaxt="n",yaxt="n",yaxs="i",bty="l",ylab="Cumulative\\nMortality (%)")	
						for(i in c(12:7,24:19,30:25,6:1,18:13)){
							polygon(c(0,1:length(surv.s[,i]),length(surv.s[,i])),c(0,ifelse(is.na(100 - surv.s[,i]),0,100 - surv.s[,i]),0),
							col=col[i],border="black",lwd=0.5)}
				axis(1,at=c(15,46,74),labels=NA) 
				axis(2,at=c(0,20))
				par(xpd=NA)
					segments(71,0,71,65,lty=2)
					text(c(7,30,60,90),rep(-7,4),c("Dec","Jan","Feb","Mar"),cex=1.2)
				par(xpd=FALSE)

				par(mar=c(1,4.5,3,3))	
				col <- c("darkred","gray","royalblue","mediumaquamarine","indianred1")
				boxplot(apply((100 -surv.s),2,mean,na.rm=TRUE) ~ as.factor(rep(1:5,each=6)),
					ylim=10^(c(-8.5,2)),log="y",ylab="Mean Daily Mortality (%)",xaxt="n",col=col,border=col)
					axis(1,at=1:5,labels=c("Arth","Brach","Chord","Echin","Moll"),cex.axis=1.2)	
					abline(median(log10(unlist(na.omit(100 -surv.s)))),0)
					text(1:5,10^c(-5,1.9,-0.9,0.2,-0.1),c("a","b","c","c","c"),cex=1.2)
					text(4,10^-7,expression(paste(italic("F")["4,25"]*" = 107.2, "*italic("P ")*"< 0.001")),cex=1.3)
				
					
			# ANOVA comparing mean mortality rates
				
					m1  <- lm(log10(apply((100 -surv.s),2,mean,na.rm=TRUE)) ~ as.factor(rep(1:5,each=6)))
					summary(aov(m1))				
					TukeyHSD(aov(m1))


			# Correlation between tolerance to acute vs chronic stress 

					a <- boxplot(log10(apply((100 -surv),2,mean,na.rm=TRUE)) ~ as.factor(rep(1:5,each=12)),plot=FALSE)$stats[3,]
					b <- boxplot(log10(apply((100 -surv.s),2,mean,na.rm=TRUE)) ~ as.factor(rep(1:5,each=6)),plot=FALSE)$stats[3,]
					cor.test(a,b)

				setwd("..")
			
		
