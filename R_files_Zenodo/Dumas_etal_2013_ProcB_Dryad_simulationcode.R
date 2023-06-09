# This script provides the code for simulations featured in Dumas, Ross-Gillespie & Kuemmerli MS: "Switching between apparently redundant iron-uptake mechanisms benefits bacteria in changeable environments". This version was produced on 15.05.2013, using R 3.0.0 on Mac OS X 10.8.3

# Note: first download and install package 'MCMCglmm' and any dependencies. Our code uses the truncated normal distribution, for which this package provides some functions. 

######################## This function simulates competition between 'players' (typically, strains) where an environment switches over time. Switch rate can be set.
###### Simulation ###### Growth-rates per player per environment must be specified (optinally with error margins). There is no interaction between players.
######################## Pop size is constant. Output includes traces (optional) and comparison of final proportions and/or integrals over time.  

Simulation <- function(RatesTable, Players="all", Environments="default", TotalRuns=200, TotalTime=100, TracePlots=T, Cols=c("black", "red", "green3", "blue"), LetGRatesFlux=T, ResTime=list("Mean"=3, "CV"=0.5, "Min"=1, "Max"=7), ExtinctLimit=1e-10, Plots=c("Endpoints", "Integrals"))
	{
	Hmisc::requirePackage(MCMCglmm)	
	if (Players[1]=="all") {Players <- levels(RatesTable[,1])}	
	if (!all(Players %in% levels(RatesTable[,1]))) {stop("One or more values specified for 'Players' invalid")}
	if (Environments[[1]]=="default") {Environments <- levels(RatesTable[,2])}
	if (!all(levels(Environments) %in% levels(RatesTable[,2]))) {stop("One or more levels of 'Environments' vector invalid")}	
	if (LetGRatesFlux & any(is.null(RatesTable[,4:5]), is.na(RatesTable[,4:5]))) {LetGRatesFlux <- FALSE; warning("Rate CIs missing. Fixed rates used.")}	
	ResTimeMean <- ResTime[[1]]; ResTimeSD <- ResTime[[1]]*ResTime[[2]]; ResTimeMin <- ResTime[[3]]; ResTimeMax <- ResTime[[4]]
	SimResults <- list(Runs=list()); 
	for (Run in 1:TotalRuns)
		{
		Simul <- cbind(data.frame(V1=0, V2=NA), as.data.frame(matrix(data=rep(1/length(Players),length(Players)), ncol=length(Players)))); 
		names(Simul) <- c("Time","Environment", Players)
		if (TracePlots) 
			{plot(x=1:TotalTime, y=rep(-1,TotalTime), xlab="Time steps", ylab="Player proportions", ylim=c(0,1)); 
				legend(x=0, y=1, Players, text.col=Cols, pch=NA, bty="n");
				points(x=rep(Simul$Time, length(Players)), y=Simul[,3:ncol(Simul)], col=Cols, pch=21)}
		while (nrow(Simul)<TotalTime+1)
			{
			ResTime <- trunc(rtnorm(n=1, mean=ResTimeMean, sd=ResTimeSD, lower=ResTimeMin, upper=ResTimeMax))
			CurrEnvir <- sample(Environments, 1); Simul$Environment <- CurrEnvir	
			while (nrow(Simul)<TotalTime+1)
				{
				for (i in 0:ResTime)
					{
					NextLine <- cbind(data.frame(V1=nrow(Simul), V2=CurrEnvir), as.data.frame(matrix(data=rep(NA,length(Players)), ncol=length(Players)))); 
					names(NextLine) <- c("Time","Environment", Players)
					for (j in Players)
						{ 	if (LetGRatesFlux)
								{ NextLine[[j]] <- Simul[[j]][nrow(Simul)] * rtnorm(1, mean=RatesTable[(RatesTable[,1]==j & RatesTable[,2]==CurrEnvir),'Est'], 
									sd=diff(as.numeric(RatesTable[(RatesTable[[1]]==j & RatesTable[[2]]==CurrEnvir),4:5]))/2, lower=0)} 
							else
								{ NextLine[[j]] <- Simul[[j]][nrow(Simul)] * RatesTable[(RatesTable[,1]==j & RatesTable[,2]==CurrEnvir),'Est']} 	}	
					NextLine[,3:ncol(NextLine)] <- NextLine[,3:ncol(NextLine)] / sum(NextLine[,3:ncol(NextLine)])  
							 # i.e. Death affects all strains uniformly. All values are scaled so that sum = 1.
					NextLine[,3:ncol(NextLine)][NextLine[,3:ncol(NextLine)] < ExtinctLimit] <- 0		 
							# Here, we apply allee effects.
					NextLine[,3:ncol(NextLine)] <- NextLine[,3:ncol(NextLine)] / sum(NextLine[,3:ncol(NextLine)])   
							# i.e. We reapply population size standardisation, in case allee correction was made.
					if (TracePlots & NextLine$Time <= TotalTime+1) {points(x=rep(NextLine$Time, length(Players)), y=NextLine[,3:ncol(NextLine)], col=Cols, pch=21)}
					Simul <- rbind(Simul, NextLine)		
					}
				if (TracePlots) {abline(v=nrow(Simul)-1, lwd=0.5)}	
				ResTime <- trunc(rtnorm(n=1, mean=ResTimeMean, sd=ResTimeSD, lower=ResTimeMin, upper=ResTimeMax))
				CurrEnvir <- sample(Environments, 1)
				}
			}		
		SimResults$Runs[[Run]] <- Simul[1:(TotalTime+1),]	
		}
	Endpoints <- data.frame()
	for (i in 1:length(SimResults$Runs)) {
		NextLine <- cbind(data.frame("Run"=i), as.data.frame(matrix(data=rep(NA,length(Players)), ncol=length(Players))))
		names(NextLine) <- c("Run", Players)
		for (j in Players) {NextLine[,j] <- SimResults$Runs[[i]][nrow(SimResults$Runs[[i]]),j]}
		Endpoints <- rbind(Endpoints, NextLine) }
	Integrals <- data.frame()
	for (i in 1:length(SimResults$Runs)) {
		NextLine <- cbind(data.frame("Run"=i), as.data.frame(matrix(data=rep(NA,length(Players)), ncol=length(Players))))
		names(NextLine) <- c("Run", Players)
		for (j in Players) {NextLine[,j] <- sum(SimResults$Runs[[i]][2:nrow(SimResults$Runs[[i]]),j])}
		Integrals <- rbind(Integrals, NextLine)	}
	if ("Endpoints" %in% Plots)	
		{quartz(); boxplot(Endpoints[,2:ncol(Endpoints)], col=Cols, ylab=paste("Final proportion after", nrow(SimResults$Runs[[1]])-1, "timesteps (",nrow(Endpoints), "runs )"))}
	if ("Integrals" %in% Plots)		
		{quartz(); boxplot(Integrals[,2:ncol(Endpoints)], col=Cols, ylab=paste("Integral of proportions over", nrow(SimResults$Runs[[1]])-1, "timesteps (",nrow(Integrals), "runs )"))}
	SimResults[["Endpoints"]] <- Endpoints
	SimResults[["Integrals"]] <- Integrals
	return(SimResults)	
	}

## Import strain- and environment-specific growth rates from monoculture assays (see Figure 2)
RatesTable <- read.table(file.choose(), header=T)
getwd()	# In the following sections of code, several data objects will be saved during simulation runs. Here's where the files will be saved. Change if desired.

## Warning: there simulations can take many hours to complete (overnight?), depending on processor speed

## Generate data for Figure 3a 
SwitchRateResults <- list()
for (SwitchRate in as.character(c(3, 5, 10, 20, 30, 50, 100))) 
	{
	SimWithRateChange <- Simulation(RatesTable, Players="all", Environments=levels(RatesTable$Environment), TotalRuns=1000, TotalTime=100, TracePlots=F, Cols=c("darkgrey", "red", "green3"), LetGRatesFlux=T, ResTime=list("Mean"=as.numeric(SwitchRate), "CV"=0.5, "Min"=c(rep(1,6),100)[which(as.character(c(3, 5, 10, 20, 30, 50, 100))==SwitchRate)], "Max"=100), ExtinctLimit=1e-10, Plots=c("Endpoints")) 
	print(colMeans(SimWithRateChange$Endpoints[,2:4]))
	SwitchRateResults[[SwitchRate]] <- as.numeric(unlist(lapply(SimWithRateChange[[1]], function(x) {rle(x$Environment)[[1]][1:length(rle(x$Environment)[[1]])-1]})))
	EndpointsTemp <- SimWithRateChange$Endpoints; #IntegralsTemp <- SimWithRateChange$Integrals; 
	save(EndpointsTemp, file=paste(paste("EndpointsSwitch",SwitchRate,sep="_"),".RData", sep=""))
	save(SwitchRateResults, file="SwitchRateTrace.RData")
	}; rm(SwitchRate); rm(EndpointsTemp); rm(SimWithRateChange); rm(SwitchRateResults)

## Generate data for Figure 3b
EndpointsRunLength <- matrix(nrow=3) 
for (RunLength in c(10, 50, 100, 500, 1000, 5000))      
	{
	SimVaryingRunLengths <- Simulation(RatesTable, Players="all", Environments=levels(RatesTable$Environment), TotalRuns=1000, TotalTime=RunLength, TracePlots=F, Cols=c("darkgrey", "red", "green3"), LetGRatesFlux=T, ResTime=list("Mean"=3, "CV"=0.5, "Min"=1, "Max"=100), ExtinctLimit=1e-10, Plots=c("Endpoints")) 
	EndpointsRunLength <- cbind(EndpointsRunLength, matrix(colMeans(SimVaryingRunLengths$Endpoints[,2:4]), nrow=3))
	print(colMeans(SimVaryingRunLengths$Endpoints[,2:4]))
	EndpointsTemp <- SimVaryingRunLengths$Endpoints; 
	save(EndpointsTemp, file=paste(paste("EndpointsRunLength",RunLength,sep="_"),".RData", sep=""))
	}; rm(RunLength); rm(EndpointsRunLength); rm(SimVaryingRunLengths); rm(EndpointsTemp)

## Generate data for Figure 3c
#BiasResults <- list()
#IntegralsBias <- matrix(nrow=3)
#EndpointsBias <- matrix(nrow=3)
EnvCounts <- as.data.frame(matrix(0, ncol=12)); names(EnvCounts) <- c("Bias", levels(RatesTable$Environment))  
for (Bias in c(1, 3, 5, 10, 30, 50))  
	{
	Environments <- append(levels(RatesTable$Environment), c(rep("5",Bias-1), rep("6",(Bias/2)-1)))
	SimWithLowpHBias <- Simulation(RatesTable, Players="all", Environments=Environments, TotalRuns=1000, TotalTime=100, TracePlots=F, Cols=c("darkgrey", "red", "green3"), LetGRatesFlux=T, ResTime=list("Mean"=3, "CV"=0.5, "Min"=1, "Max"=100), ExtinctLimit=1e-10, Plots=c("Endpoints"))
	print(colMeans(SimWithLowpHBias$Endpoints[,2:4]))
	Step1 <- unlist(lapply(SimWithLowpHBias[[1]], function(x) {summary(as.factor(x$Environment[-1]))}))
	Step2 <- data.frame(Category=names(Step1), Count=Step1)
	Step3 <- as.data.frame(as.list(by(Step2$Count, Step2$Category, sum))); names(Step3) <- sub("X", "", names(Step3))
	Step4 <- vector(); for (i in levels(RatesTable$Environment)) {if (i %in% names(Step3)) {Step4 <- append(Step4, Step3[,i])} else {Step4 <- append(Step4, NA)}}
	EnvCounts <- rbind(EnvCounts, c(Bias, Step4))	
	EndpointsTemp <- SimWithLowpHBias$Endpoints; 
	save(EndpointsTemp, file=paste(paste("EndpointsBias",Bias,sep="_"),".RData", sep=""))
	save(EnvCounts, file="EnvCounts.RData")			
	}; rm(Bias); rm(Environments); rm(EnvCounts); rm(EndpointsTemp); rm(i); rm(Step1); rm(Step2); rm(Step3); rm(Step4); rm(SimWithLowpHBias)

### FIGURE 4

quartz(width=10, height=8); layout(matrix(ncol=3, nrow=5, data=c(1,1,1,1,1,1,1,1,1,2,3,3,2,3,3), byrow=T)); 

# Figure 4a - Constancy of environment
par(mai=c(1,1,0.5,0.5), family="Times"); Cols <- c("black", "darkgrey", "lightgrey") 
EndPointsSummary <- data.frame()
stripchart(c(1,1,1), ylab="", xlab="constancy of environment", vertical=T, pch=NA, xlim=c(2.6,110)/100, ylim=c(-0.04,1.04), log="x", cex.lab=1.5, cex.axis=1.5); title(ylab="final proportion", line=4, cex.lab=1.5); 
mtext(as.character(c(3,5,10,20,30,50,100)/100), at=c(3,5,10,20,30,50,100)/100, line=1, side=1)
for (BaseLineFreq in c(3,5,10,20,30,50,100)) 
	{load(paste("EndpointsSwitch_",BaseLineFreq,".RData",sep="")); 
	stripchart(EndpointsTemp[,2:4], at=BaseLineFreq*c(0.925,1,1.075)/100, col=Cols, vertical=T, jitter=0.035*BaseLineFreq/100, method="jitter", pch=21, cex=0.9, add=T); 
	EndPointsSummary <- rbind(EndPointsSummary, c(BaseLineFreq, colMeans(EndpointsTemp)[2:4]))}; names(EndPointsSummary) <- c("Freq", "PA01", "pchEF", "pvdD")
for (i in 1:3) 
	{lines(x=EndPointsSummary$Freq*c(0.925,1,1.075)[i]/100, y=EndPointsSummary[,1+i], col=Cols[i], lwd=3.5)
	points(x=EndPointsSummary$Freq*c(0.925,1,1.075)[i]/100, y=EndPointsSummary[,1+i], col="black", lwd=2, pch=c(21,22,23)[i], cex=3, bg="white")}
mtext("(a)", line=0, adj=-0.05*(par('usr')[2]-par('usr')[1]), font=3)
rm(i); rm(EndPointsSummary); rm(BaseLineFreq); rm(Cols); rm(EndpointsTemp)

# Figure 4b - Run length	
par(mai=c(1,1,0.1,0.5), family="Times"); Cols <- c("black", "darkgrey", "lightgrey") 
EndpointsRunLength <- matrix(ncol=1, nrow=3, dimnames=list(c("PA01", "pchEF", "pvdD"),NULL))
for (RunLength in c(10, 50, 100, 500, 1000, 5000)) 
	{load(paste("EndpointsRunLength_",RunLength,".RData",sep="")); EndpointsRunLength <- cbind(EndpointsRunLength,colMeans(EndpointsTemp)[2:4])}
mtext(c(10,50,100,500,1000,5000), at=barplot(EndpointsRunLength[,-1], beside=F, las=1, xlab="", col=Cols, cex.lab=1, font.lab=2, cex.axis=1.5, cex.lab=1, las=2, names.arg=rep("",6)), side=1, line=0.5, las=2, cex=1)
title(ylab="final proportion", line=4, cex.lab=1.5); title(xlab="duration of competition", line=4, cex.lab=1.5); 
mtext("(b)", line=0,  adj=-0.06*(par('usr')[2]-par('usr')[1]), font=3)
rm(EndpointsRunLength); rm(RunLength); rm(Cols); rm(EndpointsTemp)

# Figure 4c - Bias towards low pH environments
par(mai=c(1,0.8,0.1,0.5), family="Times"); Cols <- c("black", "darkgrey", "lightgrey") 
EndPointsSummary <- data.frame()
stripchart(x=c(1,1,1), ylab="", xlab="", vertical=T, pch=NA, xlim=c(0.95,55), ylim=c(-0.04,1.04), log="x", font.lab=1, cex.lab=1.5, cex.axis=1.5); title(ylab="final proportion", line=4, font.lab=1, cex.lab=1.5); title(xlab="bias towards low-pH environments", line=4, cex.lab=1.5)
mtext(as.character(c(1,3,5,10,30,50)), at=c(1,3,5,10,30,50), line=1, side=1)
for (Bias in c(1,3,5,10,30,50)) #c(2,3,4,5,10,15,20,30,40,50,100))
	{load(paste("EndpointsBias_",Bias,".RData",sep="")); 
	stripchart(EndpointsTemp[,2:4], at=Bias*c(0.925,1,1.075), col=Cols, vertical=T, jitter=0.035*Bias, method="jitter", pch=21, cex=0.1, add=T); 
	EndPointsSummary <- rbind(EndPointsSummary, c(Bias, colMeans(EndpointsTemp)[2:4]))}; names(EndPointsSummary) <- c("Freq", "PA01", "pchEF", "pvdD")
for (i in 1:3) 
	{lines(x=EndPointsSummary$Freq*c(0.925,1,1.075)[i], y=EndPointsSummary[,1+i], col=Cols[i], lwd=3.5, cex=1.5)
	points(x=EndPointsSummary$Freq*c(0.925,1,1.075)[i], y=EndPointsSummary[,1+i], col="black", lwd=2, pch=c(21,22,23)[i], cex=3, bg="white")}
mtext("(c)", line=0, adj=-0.075*(par('usr')[2]-par('usr')[1]), font=3)
rm(EndPointsSummary); rm(i); rm(Bias); rm(EndpointsTemp); rm(Cols)



