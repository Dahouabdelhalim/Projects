
#a simulation script to explore hows different temperature curve shape propoerties affect thermal windows
#this script is to analyze ant abundance data from Sagehen Experimental Forest - 2015
#it also demonstrates a theoretical result - warming negatively affects diurnal species more than nocturnal species
#it contains all final analyses for the project
 

#written by Marshall McMunn in 2015-21 - mmcmunn@gmail.com

#clear all objects
#rm(list=ls())

#setwd to where script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load libraries
                      library("gplots")
                      library("RColorBrewer")
                      library("reshape2")
                      library("viridis")
                      library("overlapping")
                      library("gridExtra")
                      library("lubridate")
                      library("tidyverse")
                      library("ggridges")
                      library("grid")
                      library("plyr")
                      library("e1071")
                      library("scales")
                      library("egg")

    
#load data
     d <- read.csv("TempDataRaw.csv")
     dSamp <- read.csv("TempData_SampleAverages.csv")
     locationD <- read.csv("locationData.csv")
     longCTdf <- read.csv("antData.csv")
     
#time objects to R format
     d$dateTime <- as.POSIXct(d$dateTime)
     dSamp$sampleMiddle<-as.POSIXct(dSamp$sampleMiddle, format = "%m/%d/%y %H:%M")
     
     
#create dataframe for fitting nls to sagehen temperature data
     sagehenD<- data.frame(time = d$hour+d$minute/60, tempC = d$tempC, locationID = d$locationID)
     
###############
####fit nls to get a and b estimates from data.#####
    #sunrise timeset to 6am, sunset 6pm for simulation. Fitting from field data June-Oct to extract
     #just for parameters a and b, overnight cooling coefficient and lag for temperature from light
     #given 12 hour light, 12 hour dark, this to get reasonable shape parameters before normalization
          nightLogic <- sagehenD$time<=6 | sagehenD$time>= 18
          dayLogic <- sagehenD$time>6 & sagehenD$time<18
          nightBinary <- as.numeric(nightLogic)
          dayBinary <- as.numeric(dayLogic)
          
          #djusts time 0 to sunrise
          sagehenD$quantTime[nightLogic] <- ifelse( sagehenD[nightLogic, "time"] > 17,sagehenD[nightLogic, "time"] -18, 
                                                    sagehenD[nightLogic, "time"]  +  6)
          sagehenD$quantTime[dayLogic] <- sagehenD[dayLogic , "time" ] - 6
          
          tempObs <- sagehenD$tempC
          quantTime <- sagehenD$quantTime
          
          #caclulates average hourly means, pulls the lowest and highest temperature hours of the day
          lowHrlyTemp <- min(with(sagehenD, tapply(tempC , floor(time) , mean, na.rm = TRUE)))
          highHrlyTemp <- max(with(sagehenD, tapply(tempC , floor(time) , mean, na.rm = TRUE)))
          
          #fit the Parton-Logan from Sagehen data. We only use a and b from this fitting
          #to approximate the sagehen overnight cooling coefficient and temperature lag coefficitients
          nls1 <-nls( tempObs ~ ((highHrlyTemp)*(sin((pi*quantTime)/(12 + 2*a))) + lowHrlyTemp) * dayBinary + 
                nightBinary * (lowHrlyTemp + ((((highHrlyTemp - lowHrlyTemp)*(sin((pi*12)/(12+ 2*a))) + lowHrlyTemp))-lowHrlyTemp)*exp(-b*quantTime / 12))
               , start = list(a = 0.7, b = 0.9))
          
          coef(nls1)["a"]
          coef(nls1)["b"]
########
          
          
##############          
#daily temperature simulation######
          #generates a vector of 10,000 temperature points in the shape of a 
          #parton logan function (eq 1 and 2 in manuscript)
          
    #from Parton and Logan 1981
    PL_func <- function(Y=12,Z=12,Tx = highHrlyTemp, Tn = lowHrlyTemp,
                        a = coef(nls1)["a"], b = coef(nls1)["b"], samples = 10000
                        ,smoothEnd = TRUE){
      dayTimes <- seq(0, Y, length.out = floor((Y/24)*samples))
      nightTimes <- seq(0, Z, length.out = ceiling((Z/24)*samples ))
      
      dayTemps <- (Tx - Tn)*(sin((pi*dayTimes)/(Y + 2*a))) + Tn
      Ts = (Tx - Tn)*(sin((pi*Y)/(Y + 2*a))) + Tn
      
      nightTemps <- Tn + (Ts-Tn)*exp(-b*nightTimes / Z)
      Temps <- c(dayTemps, nightTemps)
      
      #smooth the transition between nighttime cooling and sunrise start temperature.
      #It affects a few dozen of 10000 temperature observations by hundredths of degree C, but avoids a sudden temp change at sunrise
      #if sunrise temperature not reached exactly overnight cooling (happens with some parameter choices)
      ifelse(smoothEnd == TRUE & Temps<Temps[samples], Temps[samples], Temps)
    }

#plot the output to ensure the function is working
          simTemp<-PL_func()
    normSimTemp <- (simTemp - min(simTemp) ) / (max(simTemp) - min(simTemp))
    #look at two normalized days
    plot(c(normSimTemp,normSimTemp), type ="l")

#################
    
    
############
####BASE SIMULATION FUNCTION - FORAGING DURATION BY WARMING#######
    #used in base model, species thermal breadth, daily thermal range, dirunally asymetric warming

    #a function to calculate the difference in foraging duration
          #between two daily temperature scenarios
     deltaForage_PL <- function(input = normSimTemp, dailyMax , dailyMin , forageMax, 
                                                forageMin, dayToDayVar = 0 , tempChangeNight, 
                                                tempChangeDay, samples = 10000){
                          
                          x = seq(0,2*pi , length = samples) #define a sequence of index values
                          y = input #from PL function
                          
                          tempChange <- ifelse(y<=median(input) , 
                                               (tempChangeNight-tempChangeDay)*((median(input)-y)/median(input)) + tempChangeDay  ,
                                               tempChangeDay) #create a temperature change vector
                          #this assumes night time temperature change is describing the change in minimum
                          #the rest of the change during the night is proportional to value of the sin wave
                          #results in largest change at minimum, no diurnal asymetry at infleciton points
                          dailyRange <- (dailyMax - dailyMin) #calculate daily range in temp
                          y <- (y*(dailyRange)) + (dailyMin) #scale the sin wave to the new temp profile
                          y <- y + rnorm(n = samples , mean = 0, sd = dayToDayVar*dailyRange)
                          yMod <- y+tempChange #apply the climate change scenario
                          propMod <- sum(yMod < forageMax & yMod > forageMin) / samples #calculate niche overlap in new scneario
                          prop <- sum(y < forageMax & y > forageMin) / samples #calculate niche overlap in old scenario
                          
                          c(prop*24*60, propMod*24*60, (propMod - prop) * 24*60) #calculate the difference, and convert to minutes per day
                          
    }
    
    #test the funciton - output is vector of minutes c(activity time before warming, activity time after warm, change in activity time)
    deltaForage_PL(dailyMax = 40, dailyMin = 10, forageMax=20, 
                   forageMin=15, dayToDayVar = 0 , tempChangeNight=4, 
                   tempChangeDay=4)
    
############
    
############    
#MODIFIED FORAGING DURATION BY WARMNING FOR NORMALLY DISTRIBUTED ACTIVITY ####
#
#modify function such that foraging activity is normally distributed instead of uniform
    #
    deltaForage_PL_norm <- function(input = normSimTemp, dailyMax , dailyMin , forageMax, 
                               forageMin, dayToDayVar = 0 , tempChangeNight, 
                               tempChangeDay, samples = 10000){
      
      x = seq(0,2*pi , length = samples) #define a sequence of index values
      y = input #from PL function
      
      tempChange <- ifelse(y<=median(input) , 
                           (tempChangeNight-tempChangeDay)*((median(input)-y)/median(input)) + tempChangeDay  ,
                           tempChangeDay) #create a temperature change vector
      #this assumes night time temperature change is describing the change in minimum
      #the rest of the change during the night is proportional to value of the sin wave
      #results in largest change at minimum, no diurnal asymetry at infleciton points
      dailyRange <- (dailyMax - dailyMin) #calculate daily range in temp
      y <- (y*(dailyRange)) + (dailyMin) #scale the sin wave to the new temp profile
      y <- y + rnorm(n = samples , mean = 0, sd = dayToDayVar*dailyRange)
      yMod <- y+tempChange #apply the climate change scenario
      
      forageMod <- yMod < forageMax & yMod > forageMin #calculate niche overlap in new scneario
      forage <- y < forageMax & y > forageMin #calculate niche overlap in old scenario
      
      #This creates a vector of binomial probabilities in warming scenario
      #that are themselves normally distributed with SD temperature range/4 and mean middle of range
      p<-dnorm(yMod, mean = (forageMax+forageMin)/2, sd = (forageMax-forageMin)/4)
      #normalizes activity probability such that activity prob at optim temperature = 1
      p<-p*(1/max(p))
      forageMod<-rbinom(length(p),size=1,p)
      
      #This creates a vector of binomial probabilities in original scenario
      p<-dnorm(y, mean = (forageMax+forageMin)/2, sd = (forageMax-forageMin)/4)
      p<-p*(1/max(p))
      forage<-rbinom(length(p),size=1,p)
      
      propMod <- sum(forageMod) / samples #calculate niche overlap in new scneario
      prop <- sum(forage) / samples #calculate niche overlap in old scenario
      
      
      
      c(prop*24*60, propMod*24*60, (propMod - prop) * 24*60) #calculate the difference, and convert to minutes per day
      
    }
##############


############
#MODEL INPUT PLOTTING FUNCTION######



############
#######Figure 1 - explanatory figure#####
############
    #this is a function to produce the vector of temperatures used in a model as an object for plots (for explanatory purposes and troubleshooting)
    delta_PL_Vector <- function(input = normSimTemp, dailyMax , dailyMin, dayToDayVar = 0 , tempChangeNight, 
                                tempChangeDay, samples = 10000){
      y = normSimTemp #global object from PL function
      
      tempChange <- ifelse(y<=median(input) , 
                           (tempChangeNight-tempChangeDay)*((median(input)-y)/median(input)) + tempChangeDay  ,
                           tempChangeDay) #create a temperature change vector
      #this assumes night time temperature change is describing the change in minimum
      #the rest of the change during the night is proportional to value of the sin wave
      #results in largest change at minimum, no diurnal asymetry at infleciton points
      dailyRange <- (dailyMax - dailyMin) #calculate daily range in temp
      y <- (y*(dailyRange)) + (dailyMin) #scale the sin wave to the new temp profile
      y <- y + rnorm(n = samples , mean = 0, sd = dayToDayVar*dailyRange)
      yMod <- y+tempChange #apply the climate change scenario
      yMod}
    #testing function works
    plot(delta_PL_Vector(dailyMax = 40, dailyMin = 10, tempChangeNight = 0, tempChangeDay = 0))

    
    #the base parameters used for figure 1, and the base simulation from which sensitivities are demonstrated
    # strategy - loop over for a 10 degree C foraging range, 0 to 50 degree range center
    
    #the input vectors for minimum and maximum activity temperatures of species
    forageMinVector <- seq(from = -5, to = 45, length.out =  400)
    forageMaxVector <- seq(from = 5, to = 55, length.out = 400)
    
    #running the simulation with species 10 degrees wide activity ranges
    #a loop inside to take each output from the deltaForage_PL function 
    #(activity time, modified acitvity time, change in activity time)
    amountShift <- matrix(nrow = 400, ncol = 3)
    start <- Sys.time()
    for( i in 1:400){
      for(j in 1:3 ){
        amountShift[i,j]  <-  deltaForage_PL(dailyMax = 40, dailyMin = 10, forageMax=forageMaxVector[i], 
                                             forageMin=forageMinVector[i], dayToDayVar = 0 ,
                                             tempChangeNight=2, 
                                             tempChangeDay=2)[j]
      }
    }
    
    
    colnames(amountShift) <- c("forageDur" , "forageDurMod" , "diffForage")
    amountShift<-as.data.frame(amountShift)
    amountShift$forageTemp    <- (forageMinVector+forageMaxVector) / 2
    
    mAmountShift<-melt(amountShift)
    
    
    #######the base simulation model - Figure 1C
    p3 <- ggplot(amountShift, aes(x = forageTemp , y = diffForage)) + 
      geom_line(size = 1.2,) + geom_hline(yintercept = 0, linetype =2) + 
      labs(y = "change in activity duration \\nafter warming (minutes)", x = "mean activity temperature (°C)") +
      theme_bw() +
      annotate("rect" , xmin = -Inf , xmax = Inf, ymax = Inf, ymin = 0, fill = "lightblue", alpha = 0.35)+
      annotate("rect" , xmin = -Inf , xmax = Inf, ymax = 0, ymin = -Inf, fill = "pink", alpha = 0.35)+
      
      annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
               y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], color = "lightyellow", cex = 14, pch = 16) +
      annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
               y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], cex = 14, pch = 1) +
      annotate("text" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
               y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], label = "2", size = 7)+
      
      annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
               y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], color = "lightyellow", cex = 14, pch = 16) +
      annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
               y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], cex = 14, pch = 1) +
      annotate("text" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
               y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], label = "1", size = 7) +ggtitle("results \\nsensitivity to warming")+ 
      theme(text = element_text(family = "Arial", size = 24),legend.text=element_text(size=24))
    
    pdf("figures/base_simulation.pdf", width = 8, height = 8)
    p3
    dev.off()
      #these are two tempearture scenarios, one with 2 degrees warming
      control<-delta_PL_Vector(input = normSimTemp, dailyMax=40 , dailyMin=10, dayToDayVar = 0 , tempChangeNight = 0, 
                               tempChangeDay=0, samples = 10000)
      evenWarm <-delta_PL_Vector(input = normSimTemp, dailyMax=40 , dailyMin=10, dayToDayVar = 0 , tempChangeNight = 2, 
                                 tempChangeDay=2, samples = 10000)
      
      #dataframe for plotting
      demoData<-cbind(control = control, warm = evenWarm)
      demoData<-melt(demoData)
      demoData <- cbind(demoData, index = seq(0,24, length.out = 10000))

      #a plot with two temperatures, and number labels to link to simulation plot
      p1 <- ggplot(demoData, aes(x = index, y = value, group = Var2))+
        scale_color_manual(values=c("black", "red"),guide = guide_legend(reverse=TRUE),labels=c("current", "2°C warming")) + labs(color = "temperature \\nscenario",x = "hours after sunrise",y = "temperature (°C)") + theme_bw() + 
        #these are all hypothetical ant distributions for the plot
        geom_line(aes(colour = Var2 ), size = 1.5)+
        
        annotate("rect" , xmin = -Inf, xmax = Inf, ymin = 29, ymax = 39, fill = "lightyellow")+
        annotate("text" ,label = "species 1 \\nactivity", x = 20, y = 34, size = 7) +
        
        annotate("rect" , xmin = -Inf, xmax = Inf, ymin = 15, ymax = 25, fill = "lightyellow")+
        annotate("text" ,label = "species 2 \\nactivity", x = 20, y = 20, size = 7, family = "Arial")+
        geom_line(aes(colour = Var2 ),linetype = "dotdash", size = 0.9) +ggtitle("simulating thermally restricted activity")+ 
        theme(text = element_text(family = "Arial", size = 24),legend.text=element_text(size=24))
      
      
      #dataframes for producing cumulative daily activity time
      control<-delta_PL_Vector(input = normSimTemp, dailyMax=40 , dailyMin=10, dayToDayVar = 0 , tempChangeNight = 0, 
                               tempChangeDay=0, samples = 10000)
      evenWarm <-delta_PL_Vector(input = normSimTemp, dailyMax=40 , dailyMin=10, dayToDayVar = 0 , tempChangeNight = 2, 
                                 tempChangeDay=2, samples = 10000)
      
      
      #dataframe for plotting species cumulative foraging
      demoData<-cbind(control = control , evenWarm = evenWarm)
      demoData <- as.data.frame(demoData) 
      demoData<-melt(demoData)
      
      #these are true false vectors for active or not for example species in figure 1
      demoData$sp2<-demoData$value>15&demoData$value<25
      demoData$sp1<-demoData$value>28&demoData$value<38
      demoData<-melt(demoData, measure.vars = c("sp1" , "sp2"))
      colnames(demoData) <- c("scenario" , "temperature" , "species" , "active")
      demoData$timeHour <- rep( seq(from = 0 , to = 24 , length.out = 10000) , 4)
      
      demoData$interaction_spSc <- paste(demoData$species, demoData$scenario, sep = ".")
      
      demoData$cumSumA<-unlist(with(demoData , tapply(active , interaction_spSc, cumsum )))

#a plot with two temperatures, and number labels to link to simulation plot
      p2 <- ggplot(demoData, aes(x = timeHour, y = cumSumA *(1440/10000), 
                                 color = scenario, group = interaction(scenario , species))) + geom_line(size = 1.2) +
        facet_wrap(~species, ncol = 1, labeller = as_labeller(c("sp1" = "species 1" , "sp2" = "species 2")))+
                                  scale_color_manual(values=c("black", "red"),guide = guide_legend(reverse=TRUE),
                                                     labels=c("current", "2°C warming")) +
                     labs(color = "temperature \\nscenario",x = "hours after sunrise",
                        y = "cumulative activity (minutes)") + theme_bw() +
                  theme(strip.background =element_rect(fill="lightyellow"))+
                ggtitle("cumulative activity since sunrise")+ theme(text = element_text(family = "Arial", size = 24),legend.text=element_text(size=24))
      #necessary to bump the third panel down to lower right
      blank <- grid.rect(gp=gpar(col="white"))
      
      #figure 1
      png("figures/figure1_Explain_revision.png", width = 16, height = 12, res = 1200, units = "in")
      grid.arrange(p1,blank,p2,p3, ncol = 2)
      dev.off()
     
      
      
      
#senstivity analyses
      #basic strategy is to vary one parameter at a time and compare
      # to base model. all simulations are +2 degree C.
      
################

      
#sensitivity - normal distribution instead of uniform######
#same parameteres as base simulation, but being plugged into the normal distribution version of deltaForage_PL
      forageMinVector <- seq(from = -5, to = 45, length.out =  400)
      forageMaxVector <- seq(from = 5, to = 55, length.out = 400)
#a new object to catch loop output      
      amountShiftnorm <- matrix(nrow = 400, ncol = 3)
      #this simulation is stochastic, so we need to set a seed
      set.seed(27134)
      for( i in 1:400){
        for(j in 1:3 ){
          #here we are calling our modified simulation with the normally distributed binomial vector (see above)
          amountShiftnorm[i,j]  <-  deltaForage_PL_norm(dailyMax = 40, dailyMin = 10, forageMax=forageMaxVector[i], 
                                               forageMin=forageMinVector[i], dayToDayVar = 0 ,
                                               tempChangeNight=2, 
                                               tempChangeDay=2)[j]
        }
      }
      
      #rename columns of output
      colnames(amountShiftnorm) <- c("forageDur" , "forageDurMod" , "diffForage")
      amountShiftnorm<-as.data.frame(amountShiftnorm)
      amountShiftnorm$forageTemp <- (forageMinVector+forageMaxVector) / 2
      
      #create a column to allow rbind of normal and uniform together
      #this makes it easier to compare results on single panel ggplot
      amountShift$distType <- "uniform"
      amountShiftnorm$distType <- "normal"
      amountShiftnorm<-rbind(amountShift , amountShiftnorm)
      
      #drop the label column since this object "amountShift" is the base model output
      # and is reused in subsequent sensitivity analyses
      amountShift<-amountShift[,colnames(amountShift)!="distType"]
      
      png("figures/norm_2degree_sensitivity.png",height = 8 ,width = 8 , units = "in", res = 500)
      p2 <- ggplot(amountShiftnorm, aes(x = forageTemp , y = diffForage, group = distType)) + 
        geom_line(size = 1.2,aes(linetype = distType)) + geom_hline(yintercept = 0, linetype =2) + 
        labs(y = "change in daily activity \\nduration after warming (minutes)", x = "median activity temperature (°C)") +
        theme_bw() +
        annotate("rect" , xmin = -Inf , xmax = Inf, ymax = Inf, ymin = 0, fill = "lightblue", alpha = 0.35)+
        annotate("rect" , xmin = -Inf , xmax = Inf, ymax = 0, ymin = -Inf, fill = "pink", alpha = 0.35)+
        #add the example species to the base simulation to enable comparison
        annotate("point" , x = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>20)==1 , "forageTemp"],
                 y = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>20&amountShiftnorm$distType=="uniform")==1 , "diffForage"], color = "lightyellow", cex = 8, pch = 16) +
        annotate("point" , x = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>20)==1 , "forageTemp"],
                 y = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>20&amountShiftnorm$distType=="uniform")==1 , "diffForage"], cex = 8, pch = 1) +
        annotate("text" , x = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>20)==1 , "forageTemp"],
                 y = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>20&amountShiftnorm$distType=="uniform")==1 , "diffForage"], label = "2")+
        
        annotate("point" , x = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>33)==1 , "forageTemp"],
                 y = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>33&amountShiftnorm$distType=="uniform")==1 , "diffForage"], color = "lightyellow", cex = 8, pch = 16) +
        annotate("point" , x = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>33)==1 , "forageTemp"],
                 y = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>33&amountShiftnorm$distType=="uniform")==1 , "diffForage"], cex = 8, pch = 1) +
        annotate("text" , x = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>33)==1 , "forageTemp"],
                 y = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>33&amountShiftnorm$distType=="uniform")==1 , "diffForage"], label = "1") +
        ggtitle("normal vs. uniform")+
        labs(linetype = "distribution type")+ labs(tag = "B")+theme(text = element_text(size=28),legend.text=element_text(size=28))+
        guides(scale = "none", color=guide_legend(override.aes=list(size=2)))+theme(text = element_text(family="Arial"))
      
p2      
dev.off()
################

#sensitivity - diurnal asymmetry (3 and 1)######

#fig 2a
#look at the effect diurnally asymmetric warming. Still 2 degrees average, but 3 at night and 1 during day.
forageMinVector <- seq(from = -5, to = 45, length.out =  400)
forageMaxVector <- seq(from = 5, to = 55, length.out = 400)

amountShiftDiurnAsym <- matrix(nrow = 400, ncol = 3)

for( i in 1:400){
  for(j in 1:3 ){
    #this is how we add diurnal asymmetry, have more warming at night than day. with 12 hours night and day, it averages close to the same warming
    #this isn't exact, because there is a dampening function that is applied during the night
    #without this, you create a jagged change in slope between the day and the nighttime cooling, which is unrealistic and
    #alters the model output
    amountShiftDiurnAsym[i,j]  <-  deltaForage_PL(dailyMax = 40, dailyMin = 10, forageMax=forageMaxVector[i], 
                                              forageMin=forageMinVector[i], dayToDayVar = 0 ,
                                              tempChangeNight=3, 
                                              tempChangeDay=1)[j]
  }
}

#save sensitivity analysis
colnames(amountShiftDiurnAsym) <- c("forageDur" , "forageDurMod" , "diffForage")
amountShiftDiurnAsym<-as.data.frame(amountShiftDiurnAsym)
amountShiftDiurnAsym$forageTemp    <- (forageMinVector+forageMaxVector) / 2

#paste together with the base model
amountShiftDiurnAsym$scenario<-"warm nights"
amountShift$scenario<-"even warming"
amountShiftDiurnAsym<-rbind(amountShiftDiurnAsym , amountShift)
#drop this column so base model can be reused
amountShift<-amountShift[,colnames(amountShift)!="scenario"]


mAmountShift<-melt(amountShiftDiurnAsym)

#add the example species to the base simulation to enable comparison
png("figures/diurnAsym_1_3_degree_sensitivity.png",height = 8 ,width = 8, units = "in", res = 500 )
p3 <- ggplot(amountShiftDiurnAsym, aes(x = forageTemp , y = diffForage, group = scenario)) + 
  geom_line(size = 1.2,aes(linetype = scenario)) + geom_hline(yintercept = 0, linetype =2) + 
  labs(y = "change in activity duration \\nafter warming (minutes)", x = "median activity temperature (°C)") +
  theme_bw() +
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = Inf, ymin = 0, fill = "lightblue", alpha = 0.35)+
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = 0, ymin = -Inf, fill = "pink", alpha = 0.35)+
  
  annotate("point" , x = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>20&amountShiftDiurnAsym$scenario=="even warming")==1 , "forageTemp"],
           y = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>20&amountShiftDiurnAsym$scenario=="even warming")==1 , "diffForage"], color = "lightyellow", cex = 8, pch = 16) +
  annotate("point" , x = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>20&amountShiftDiurnAsym$scenario=="warm nights")==1 , "forageTemp"],
           y = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>20&amountShiftDiurnAsym$scenario=="even warming")==1 , "diffForage"], cex = 8, pch = 1) +
  annotate("text" , x = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>20&amountShiftDiurnAsym$scenario=="even warming")==1 , "forageTemp"],
           y = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>20&amountShiftDiurnAsym$scenario=="even warming")==1 , "diffForage"], label = "2")+
  
  annotate("point" , x = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>33&amountShiftDiurnAsym$scenario=="warm nights")==1 , "forageTemp"],
           y = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>33&amountShiftDiurnAsym$scenario=="even warming")==1 , "diffForage"], color = "lightyellow", cex = 8, pch = 16) +
  annotate("point" , x = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>33&amountShiftDiurnAsym$scenario=="even warming")==1 , "forageTemp"],
           y = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>33&amountShiftDiurnAsym$scenario=="even warming")==1 , "diffForage"], cex = 8, pch = 1) +
  annotate("text" , x = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>33&amountShiftDiurnAsym$scenario=="even warming")==1 , "forageTemp"],
           y = amountShiftDiurnAsym[cumsum(amountShiftDiurnAsym$forageTemp>33&amountShiftDiurnAsym$scenario=="even warming")==1 , "diffForage"], label = "1") +
  ggtitle("diurnal asymmetry")+ labs(tag = "C")+theme(text = element_text(size=28),legend.text=element_text(size=28))+
  guides(scale = "none", color=guide_legend(override.aes=list(size=2)))+theme(text = element_text(family="Arial"))

p3      
dev.off()

################

#sensitivity - daily variance#######
#compare low, high, and zero between day temperature variance
forageMinVector <- seq(from = -5, to = 45, length.out =  400)
forageMaxVector <- seq(from = 5, to = 55, length.out = 400)

#high daily variance scenario
amountShiftHighVar <- matrix(nrow = 400, ncol = 3)
set.seed(27134)
for( i in 1:400){
  for(j in 1:3 ){
    #day-to-day variation is just adding normally distributed noise to temperature observations.
    #the noise is added as SD = dayToDayVar * daily temperature range. So we first add noise with SD = 0.3*30 = 9 degree C
    amountShiftHighVar[i,j]  <-  deltaForage_PL(dailyMax = 40, dailyMin = 10, forageMax=forageMaxVector[i], 
                                         forageMin=forageMinVector[i], dayToDayVar = 0.3 ,
                                         tempChangeNight=3, 
                                         tempChangeDay=1)[j]
  }
}

#assemble sensitivity model outputs for plotting
colnames(amountShiftHighVar) <- c("forageDur" , "forageDurMod" , "diffForage")
amountShiftHighVar<-as.data.frame(amountShiftHighVar)
amountShiftHighVar$forageTemp    <- (forageMinVector+forageMaxVector) / 2


#low daily variance scenario
set.seed(27134)
amountShiftLowVar <- matrix(nrow = 400, ncol = 3)
for( i in 1:400){
  for(j in 1:3 ){
    #the noise is added for low variance is SD = dayToDayVar * daily temperature range. So we first add noise with SD = 0.1*30 = 3 degree C
    amountShiftLowVar[i,j]  <-  deltaForage_PL(dailyMax = 40, dailyMin = 10, forageMax=forageMaxVector[i], 
                                                forageMin=forageMinVector[i], dayToDayVar = 0.1 ,
                                                tempChangeNight=3, 
                                                tempChangeDay=1)[j]
  }
}


#assemble sensitivity model outputs for plotting
colnames(amountShiftLowVar) <- c("forageDur" , "forageDurMod" , "diffForage")
amountShiftLowVar<-as.data.frame(amountShiftLowVar)
amountShiftLowVar$forageTemp    <- (forageMinVector+forageMaxVector) / 2


#label the groups for plotting in ggplot
amountShift$variance<-0
amountShiftLowVar$variance<-0.1
amountShiftHighVar$variance<-0.3

#rbind the 3 simulations together
amountShiftVariance<-rbind(amountShift, amountShiftLowVar, amountShiftHighVar)
amountShift<-amountShift[,colnames(amountShift)!="variance"]


png("figures/dailVar_sensitivty.png",height = 8 ,width = 8 , units = "in", res = 500)
p1 <- ggplot(amountShiftVariance, aes(x = forageTemp , y = diffForage, group = as.factor(variance))) + 
  geom_line(size = 1.2,aes(color = as.factor(variance))) + geom_hline(yintercept = 0, linetype =2) + 
  labs(y = "change in activity duration \\nafter warming (minutes)", x = "median activity temperature (°C)") +
  theme_bw() +
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = Inf, ymin = 0, fill = "lightblue", alpha = 0.35)+
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = 0, ymin = -Inf, fill = "pink", alpha = 0.35)+scale_color_brewer(palette = "Dark2")+
 #label example species 
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], color = "lightyellow", cex = 8, pch = 16) +
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], cex = 8, pch = 1) +
  annotate("text" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], label = "2")+
  
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], color = "lightyellow", cex = 8, pch = 16) +
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], cex = 8, pch = 1) +
  annotate("text" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], label = "1") +ggtitle("day-to-day variation in temperature")+
  labs(color = "daily temperature\\nvariation \\nproportion \\nSD / daily range")+ labs(tag = "A")+theme(text = element_text(size=28),legend.text=element_text(size=28))+
  guides(scale = "none", color=guide_legend(override.aes=list(size=2)))+theme(text = element_text(family="Arial"))
p1      
dev.off()


################

#sensitivity - activity breadth######
#compare the effects of altering foraging breadth (base 10 degree, also run 5 and 20)

#narrow foraging - 5 degree activity windows
forageMinVector <- seq(from = -2.5, to = 47.5, length.out =  400)
forageMaxVector <- seq(from =2.5, to = 52.5, length.out = 400)

amountShiftNarrow <- matrix(nrow = 400, ncol = 3)
start <- Sys.time()

#here we are looping through species that are active at different temperatures, all with small activity ranges (5C)
for( i in 1:400){
  for(j in 1:3 ){
    amountShiftNarrow[i,j]  <-  deltaForage_PL(dailyMax = 40, dailyMin = 10, forageMax=forageMaxVector[i], 
                                         forageMin=forageMinVector[i], dayToDayVar = 0 ,
                                         tempChangeNight=2, 
                                         tempChangeDay=2)[j]
  }
}
colnames(amountShiftNarrow) <- c("forageDur" , "forageDurMod" , "diffForage")
amountShiftNarrow<-as.data.frame(amountShiftNarrow)
amountShiftNarrow$forageTemp    <- (forageMinVector+forageMaxVector) / 2

#broad foraging - 20 degree foraging window
forageMinVector <- seq(from = -10, to = 40, length.out =  400)
forageMaxVector <- seq(from = 10, to = 60, length.out = 400)

#here we are looping through species that are active at different temperatures, all with small activity ranges (20C)
amountShiftBroad <- matrix(nrow = 400, ncol = 3)
for( i in 1:400){
  for(j in 1:3 ){
    amountShiftBroad[i,j]  <-  deltaForage_PL(dailyMax = 40, dailyMin = 10, forageMax=forageMaxVector[i], 
                                               forageMin=forageMinVector[i], dayToDayVar = 0 ,
                                               tempChangeNight=2, 
                                               tempChangeDay=2)[j]
  }
}

colnames(amountShiftBroad) <- c("forageDur" , "forageDurMod" , "diffForage")
amountShiftBroad<-as.data.frame(amountShiftBroad)
amountShiftBroad$forageTemp    <- (forageMinVector+forageMaxVector) / 2

#these are labels for the different scenarios
amountShiftNarrow$foragingBreadth <- 5
amountShiftBroad$foragingBreadth <- 20
amountShift$foragingBreadth <- 10
amountShiftBreadth<-rbind(amountShift , amountShiftBroad, amountShiftNarrow)


png("figures/breadth_sensitivity.png",height = 8 ,width = 8 , units = "in", res = 500)
p4 <- ggplot(amountShiftBreadth, aes(x = forageTemp , y = diffForage , group = as.factor(foragingBreadth))) + 
  geom_line(size = 1.2,aes(color = as.factor(foragingBreadth))) + geom_hline(yintercept = 0, linetype =2) + 
  labs(y = "change in activity duration \\nafter warming (minutes)", x = "median activity temperature (°C)") +
  theme_bw() +
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = Inf, ymin = 0, fill = "lightblue", alpha = 0.35)+
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = 0, ymin = -Inf, fill = "pink", alpha = 0.35)+
    ggtitle("species activity breadth")+
  labs(color = "activity breadth\\n°C")+ scale_color_brewer(palette = "Dark2")+
  
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], color = "lightyellow", cex = 20, pch = 16) +
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], cex = 20, pch = 1) +
  annotate("text" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp",
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], label = "2", size = 20])+
  
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], color = "lightyellow", cex = 20, pch = 16) +
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], cex = 20, pch = 1) +
  annotate("text" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], label = "1", size = 20)+theme(text = element_text(size=28),legend.text=element_text(size=28))+
  guides(scale = "none", color=guide_legend(override.aes=list(size=2)))+theme(text = element_text(family="Arial"))
p4     
dev.off()

#drop the label column from the base simulation for reuse
amountShift<-amountShift[,colnames(amountShift)!="foragingBreadth"]


################

#sensitivity - daily range######
#compare the effects of altering daily max-min (range)

#narrow foraging - 5 degree activity window
forageMinVector <- seq(from = -5, to = 45, length.out =  400)
forageMaxVector <- seq(from = 5, to = 55, length.out = 400)

amountShiftNarrow <- matrix(nrow = 400, ncol = 3)
start <- Sys.time()

#here we are looping through scenarios with narrow daily temperature ranges
for( i in 1:400){
  for(j in 1:3 ){
    amountShiftNarrow[i,j]  <-  deltaForage_PL(dailyMax = 35, dailyMin = 15, forageMax=forageMaxVector[i], 
                                               forageMin=forageMinVector[i], dayToDayVar = 0 ,
                                               tempChangeNight=2, 
                                               tempChangeDay=2)[j]
  }
}
colnames(amountShiftNarrow) <- c("forageDur" , "forageDurMod" , "diffForage")
amountShiftNarrow<-as.data.frame(amountShiftNarrow)
amountShiftNarrow$forageTemp    <- (forageMinVector+forageMaxVector) / 2

#broad foraging - 20 degree foraging window
forageMinVector <- seq(from = -5, to = 45, length.out =  400)
forageMaxVector <- seq(from = 5, to = 55, length.out = 400)

#here we are looping through scenarios with broad daily temperature ranges
amountShiftBroad <- matrix(nrow = 400, ncol = 3)
for( i in 1:400){
  for(j in 1:3 ){
    amountShiftBroad[i,j]  <-  deltaForage_PL(dailyMax = 45, dailyMin = 5, forageMax=forageMaxVector[i], 
                                              forageMin=forageMinVector[i], dayToDayVar = 0 ,
                                              tempChangeNight=2, 
                                              tempChangeDay=2)[j]
  }
}

colnames(amountShiftBroad) <- c("forageDur" , "forageDurMod" , "diffForage")
amountShiftBroad<-as.data.frame(amountShiftBroad)
amountShiftBroad$forageTemp    <- (forageMinVector+forageMaxVector) / 2

#these are labels for the different scenarios
amountShiftNarrow$dailyRange <- 20
amountShiftBroad$dailyRange <- 40
amountShift$dailyRange <- 30
#amountShift<-amountShift[ , -which(colnames(amountShift)=="variance")]
amountShiftBreadth<-rbind(amountShift , amountShiftBroad, amountShiftNarrow)

png("figures/dailyRange_sensitivity.png",height = 8 ,width = 8 )
p5 <- ggplot(amountShiftBreadth, aes(x = forageTemp , y = diffForage , group = as.factor(dailyRange))) + 
  geom_line(size = 1.2,aes(color = as.factor(dailyRange))) + geom_hline(yintercept = 0, linetype =2) + 
  labs(y = "change in activity duration \\nafter warming (minutes)", x = "median activity temperature (°C)") +
  theme_bw() +
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = Inf, ymin = 0, fill = "lightblue", alpha = 0.35)+
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = 0, ymin = -Inf, fill = "pink", alpha = 0.35)+
  ggtitle("daily temperature range")+
  labs(color = "daily \\ntemperature \\nrange°C")+ labs(tag = "E")+scale_color_brewer(palette = "Dark2")+
  
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], color = "lightyellow", cex = 8, pch = 16) +
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], cex = 8, pch = 1) +
  annotate("text" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], label = "2")+
  
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], color = "lightyellow", cex = 8, pch = 16) +
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], cex = 8, pch = 1) +
  annotate("text" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], label = "1")+theme(text = element_text(size=28),legend.text=element_text(size=28))+
  guides(scale = "none", color=guide_legend(override.aes=list(size=2)))+theme(text = element_text(family="Arial"))
p5     
dev.off()
amountShift<-amountShift[ , colnames(amountShift)!="dailyRange"]

png("figures/figure2.png" ,  width = 22, height = 24, units = "in", res = 1200)
grid.draw(ggarrange(p1,p2,p3,p4,p5, ncol = 2))
dev.off()

###############

#ant data plots#####
#these species changed names
longCTdf[longCTdf$species=="Formica.fusca","species"]<-"Formica.subaenescens"
longCTdf[longCTdf$species=="Camponotus.laevigatus","species"]<-"Camponotus.laevissimus"

#only want to use species with enough counts to estimate activity temperature
speciesThermal<-names(table(longCTdf$species)[table(longCTdf$species)>20])


#normalize ant data for an overall average activity distribution
#vectors of means and SD by species

#spMeans<-with( longCTdf , tapply(meanHrlyTemp , species , mean))
#spSDs<-with( longCTdf , tapply(meanHrlyTemp , species , sd))

#mean(spSDs)
#for loop - each species temperature observations are grouped into an element within a list
listSppTemp <- list()
for(i in speciesThermal){
  listSppTemp[[i]] <- longCTdf[longCTdf$species ==i ,"meanHrlyTemp" ] 
}

#here we add environmental temperature observations as an element of the list
listSppTemp[["Environment"]] <- d[!is.na(d$tempC) , "tempC"]
#and a warming scenario as an element of the list
listSppTemp[["EnvironmentHot"]] <- d[!is.na(d$tempC) , "tempC"] + 2

#this is the empirical estimation of density overlap in secnarios
#it is automatically performed for all pairwise comparisons within the list
overlapActTemp <- overlap(listSppTemp , 1000)

#we want to retrieve just the output where species are being compared to environment or environment+2
#split up names
splitName <- strsplit(names(overlapActTemp$OV) , "-")
#first elements in each name
name1 <- sapply(splitName, "[", 1)
#second elements in each name
name2 <- sapply(splitName, "[", 2)

#confirm where in the object we will retrieve the overlap estimates from
overlapActTemp$DD$dominance
overlapActTemp$OV

#convert the messy list output into a dataframe
overlapDF <- data.frame(name1 = name1, name2 = name2, overlap = overlapActTemp$OV)

#we really just care about how each species overlaps with environment and environment+2
envOverlap<-overlapDF[overlapDF$name1=="Environment" | overlapDF$name1=="EnvironmentHot" | overlapDF$name2=="Environment" | overlapDF$name2=="EnvironmentHot" , ]
temp <- envOverlap[envOverlap$name1!="Environment" & envOverlap$name1!="EnvironmentHot" , ]
colnames(temp) <- c("name2" , "name1" , "overlap")
envOverlap <- envOverlap[-which(envOverlap$name1!="Environment" & envOverlap$name1!="EnvironmentHot") , ]

envOverlap <- rbind( temp , envOverlap)

#reorganize dataframe to enable calculation of the difference between regular scenario and hot scenario
envOverlapW <- envOverlap %>% spread(name1, overlap)

#caluclate difference between scenarios
envOverlapW$deltaHotEnv <- envOverlapW$EnvironmentHot - envOverlapW$Environment

#calculate mean forage temperatures for x axis of plot
meanForageTemp <- with(longCTdf , tapply(tempC, species, mean))

#add mean forage temperatures to dataframe for plotting
envOverlapW$meanForage <- meanForageTemp[match(envOverlapW$name2 , names(meanForageTemp) ) ]   

#subset just the species that are being looked at (20 or more individuals)
speciesShift <- subset(envOverlapW,name2 %in% speciesThermal)

#how many hours changed in foraging
speciesShift$EnvironmentHot / speciesShift$Environment
dev.off()
#convert change in foraging duration in minutes
envOverlapW$forageMins <- envOverlapW$deltaHotEnv*24*60
#linear model of relationship
lmOverlapMean<-with(envOverlapW, lm(forageMins~meanForage))
summary(lmOverlapMean)

#panel C of figure 3
p3<-ggplot(data = envOverlapW, aes(x = meanForage, y = forageMins, color = meanForage))+scale_colour_viridis(option = "plasma", limits=c(0, 60))+
              geom_point()+theme_bw()+ ggtitle("predicted change in species activity duration \\nunder 2°C simulated warming")+
  geom_abline(intercept = coef(lmOverlapMean)[1], slope = coef(lmOverlapMean)[2])+
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "mean foraging temperature(°C)", y = "change in duration \\nof daily activity (minutes)")+ theme(legend.position="none")+
theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=24),legend.text=element_text(size=24))+theme(text = element_text(family="Arial"))
p3

#ggridge plots
#by species - activity temperature distributions (panel B of figure 3)
p4<-ggplot(subset(longCTdf,species %in% speciesThermal), aes(x = meanHrlyTemp, y = reorder(species, meanHrlyTemp, FUN=median), fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1) +
  scale_x_continuous(expand = c(0.01, 0), limits = c(-5,65)) +theme_bw()+
  scale_y_discrete(expand = c(0.01, 0)) +ggtitle("activity temperature by species")+
  scale_fill_viridis(name = "temperature (°C)", option = "C") +
  labs( x = "temperature (°C)",y="") +
  theme_ridges(font_size = 24, grid = TRUE) + theme(axis.title.y = element_blank())+
  scale_y_discrete(labels = rev(c("Temnothorax rugatulus" , "Formica argentea" , "Formica lasiodes", "Formica subaenescens", "Formica accreta",
                                  "Formica obscuripes", "Formica sibylla" , "Manica invidia" , "Formica ravida" , "Formica sp. c.f. sibylla",
                                  "Temnothorax nevadensis","Camponotus modoc", "Tapinoma sessile" ,
                                  "Temnothorax nitens" ,"Myrmica tahoensis" )))  +coord_flip()+
  theme(axis.text.x = element_text(angle = 270))+theme(panel.border=element_rect(fill=NA,color="black", size=0.5, 
                                                                      linetype="solid"))+theme(text = element_text(24),legend.text=element_text(size=24))+theme(text = element_text(family="Arial"))+
  guides(shape = guide_legend(override.aes = list(size = 10)))
p4



#a quantitative time of day vector for plotting
d$timeofDay<-as.POSIXlt(d$sampleMiddle)$hour+as.POSIXlt(d$sampleMiddle)$min/60+as.POSIXlt(d$sampleMiddle)$sec/(24*60)
as.POSIXlt(d$sampleMiddle)
#ant present vector
d$antPresent<-d$sampleID%in%longCTdf$sampleID

#this resets zero on plot to 6am (around sunrise)
d$timeofDayPlot<-ifelse(d$timeofDay<6 , d$timeofDay+18, d$timeofDay-6)

#plot of ants present vs absent in hourly samples
p6<-ggplot(data = d, aes(x = timeofDayPlot, y = meanHrlyTemp, color = antPresent))+
geom_point(alpha = 0.03, size = 0.4)+theme_bw()+
  ylab("mean hourly temperature(°C)")+labs( x = "hours since sunrise")+
  scale_color_manual(name = "",values=c("red", "blue"),guide = guide_legend(reverse=TRUE),
                     labels=c("ants absent", "ants present"))+
  ggtitle("ant activity by time of day and temperature")+
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 0.2) ) )+ylim(c(-5,65))+theme(text = element_text(size=24),legend.text=element_text(size=24))+theme(text = element_text(family="Arial"))
p6

#plot site data  
d$timePlot<-as.POSIXct(d$time, format = "%H:%M:%S" )

#a function to calculate quantiles
qs <- c(0.05,0.25,0.50,0.75,0.95)


#temperature quantile plot
p7<-ggplot(data = d, aes(x = timeofDayPlot  , y = tempC))+
  geom_point(alpha = 0.03, size = 0.2)+theme_bw()+ 
  geom_quantile(quantiles=qs, formula=y ~ poly(x, 6), colour="blue") +
  xlab("time since sunrise")+ylab("temperature(°C)")+
  ggtitle("temperature quantiles through the day")+ ylim(c(-5,65))+ scale_colour_brewer(type = "div", palette = "Dark2")+theme(text = element_text(size=24),legend.text=element_text(size=24))+theme(text = element_text(family="Arial"))
p7


png(file = "figures/forageTempCombined.png", width = 20, height = 16, units = "in",res=1200)
grid.draw(ggarrange(p6,p4,p7,p3, ncol = 2))
dev.off()


dim(longCTdf)
unique(longCTdf$species)

sumMed<-with(longCTdf , tapply(meanHrlyTemp , species ,median, na.rm = TRUE ))
sumLength<-with(longCTdf , tapply(meanHrlyTemp , species ,length))
sumMean<-with(longCTdf , tapply(meanHrlyTemp , species ,mean, na.rm=TRUE))
sumVar<-with(longCTdf , tapply(meanHrlyTemp , species ,var, na.rm=TRUE))
sumSkew<-with(longCTdf , tapply(meanHrlyTemp , species ,skewness, na.rm=TRUE))
sumKurt<-with(longCTdf , tapply(meanHrlyTemp , species ,kurtosis, na.rm=TRUE))

sumTable<-cbind(sumLength,sumMean,sumMed,sumVar,sumSkew,sumKurt )
sumTable<-sumTable[  order(-sumLength) , ]
colnames(sumTable)<-c( "ants collected", "mean active temperature", "median active temperature",
                      "v ariance","skewness", "kurtosis")
write.csv(sumTable, "figures/summaryAntDistributions.csv")


################

#Supplemental Figure 2
#CT max####
CT<-longCTdf[match(unique(longCTdf$species) , longCTdf$species), c("species","CTmin", "CTmax","CTmaxN","CTminN") ]
medTemps<-data.frame(medTemps=with(longCTdf , tapply(meanHrlyTemp , species ,median, na.rm=TRUE )))
CT<-cbind(CT, medTemps)
CT$CTavg<-(CT$CTmin+CT$CTmax)/2
CT<-CT[CT$species%in%speciesThermal , ]
lmCTmin<-with(CT, lm(medTemps~CTmin))
summary(lmCTmin)

lmCTmax<-with(CT, lm(medTemps~CTmax))
summary(lmCTmax)
cor.test(CT$CTmin, CT$CTmax)
p10<-ggplot(data = CT, aes(x = CTmax, y = medTemps)) +geom_point() +geom_smooth(method = "lm") +
  ylab("median species \\nactivity temperature (ºC)") +xlab("CTmax (ºC)") + theme_bw()
pdf("figures/CTmax.pdf",height = 8 ,width = 8 )
p10
dev.off()

p11<-ggplot(data = CT, aes(x = CTmin, y = medTemps)) +geom_point() +geom_smooth(method = "lm") +
  ylab("median species \\nactivity temperature (ºC)") +xlab("CTmin (ºC)") + theme_bw()
pdf("figures/CTmin.pdf",height = 8 ,width = 8 )
p11
dev.off()

##############


#Supplemental Figure 3
#sensitivity - accumulation over season (using air temp, ignoring snow/rain)####

write.csv(dailySumm , "sagehen_dailySummaries.csv" )
#looks fine for a temperatre climate, but this parameterization (for a and b)
#does not work for polar sites (extreme night/day length)
plot(PL_func(Y=10, Z=14))
plot(PL_func(Y=10, Z=14))

#upsampling for instantaneous temp data at sagehen
#create a vector for each day of the year with sagehen day length and min/max data
out<-list()
for(i in 1:365){
  out[[i]]<-PL_func(Y=dailySumm$dayLength[i],Z=dailySumm$nightLength[i],Tx = dailySumm$dailyMax_tempC[i], Tn = dailySumm$dailyMin_tempC[i])
}

#seem like they all look good
plot(out[[10]])
plot(out[[241]])

#between two daily temperature scenarios
deltaForage_PL_sagehen <- function(input = out[[140]],forageMax,forageMin,
                                   tempChange=2, samples = 10000){
  y<-input
  yMod <- y+tempChange #apply the climate change scenario
  propMod <- sum(yMod < forageMax & yMod > forageMin) / samples #calculate niche overlap in new scneario
  prop <- sum(y < forageMax & y > forageMin) / samples #calculate niche overlap in old scenario
  c(prop*24*60, propMod*24*60, (propMod - prop) * 24*60) #calculate the difference, and convert to minutes per day
  
}

#new function seems to work
deltaForage_PL_sagehen(forageMax = 40, forageMin=-5)

#a loop inside to take each output from the deltaForage_PL function 
#(activity time, modified acitvity time, change in activity time)
forageMinVector <- seq(from = -35, to = 45, length.out =  800)
forageMaxVector <- seq(from = -25, to = 55, length.out = 800)

amountShiftSage <- array(dim = c(800, 3,365))
for(h in 1:365){
  for( i in 1:800){
    for(j in 1:3 ){
      amountShiftSage[i,j,h]  <-  deltaForage_PL_sagehen(input = out[[h]],forageMax=forageMaxVector[i], 
                                                         forageMin=forageMinVector[i], 
                                                         tempChange=2)[j]
    }
  }
}



yearAvg<-apply(amountShiftSage,2,rowMeans, na.rm=TRUE)
#summerAvg<-apply(amountShiftSage[,,172:264],2,rowMeans, na.rm=TRUE)
#fallAvg<-apply(amountShiftSage[,,265:354],2,rowMeans, na.rm=TRUE)
#winterAvg<-apply(amountShiftSage[,,c(355:365,1:79)],2,rowMeans, na.rm=TRUE)
#springAvg<-apply(amountShiftSage[,,c(80:171)],2,rowMeans, na.rm=TRUE)

yearAvg<-as.data.frame(yearAvg)
#summerAvg<-as.data.frame(summerAvg)
#fallAvg<-as.data.frame(fallAvg)
#winterAvg<-as.data.frame(winterAvg)
#springAvg<-as.data.frame(springAvg)

yearAvg$season<-"allYear"
#summerAvg$season<-"summer"
#fallAvg$season<-"fall"
#winterAvg$season<-"winter"
#springAvg$season<-"spring"

colnames(yearAvg)[3]<-"changeMinPerDay"
#colnames(summerAvg)[3]<-"changeMinPerDay"
#colnames(springAvg)[3]<-"changeMinPerDay"
#colnames(fallAvg)[3]<-"changeMinPerDay"
#colnames(winterAvg)[3]<-"changeMinPerDay"

#allYearAvg<-rbind(yearAvg , summerAvg, fallAvg, winterAvg, springAvg)
allYearAvg <- yearAvg
allYearAvg$forageTemp    <- (forageMinVector+forageMaxVector) / 2

plot(allYearAvg$changeMinPerDay)
p8 <- ggplot(allYearAvg, aes(x = forageTemp , y = changeMinPerDay)) + 
  geom_hline(yintercept = 0, linetype =2) +geom_line()+
  labs(y = "change in daily activity \\nduration after warming (minutes)", x = "median activity temperature (°C)") +
  theme_bw()+ggtitle("averaged effect across year \\nlocal air temperature")+
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = Inf, ymin = 0, fill = "lightblue", alpha = 0.35)+
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = 0, ymin = -Inf, fill = "pink", alpha = 0.35)
pdf("figures/seasonal_sensitivity.pdf",height = 6 ,width = 8 )
p8
dev.off()




#############


################
#Supplemental Figure 4
#sensitivity - day length#####


#15 light, 9 dark - 9 light, 15 dark, - 12 light, 12 dark
nightLength <- c(9, 12, 15)  
dayLength <- c(15, 12, 9)
#upsampling for instantaneous temp data at sagehen
#create a vector for each day of the year with sagehen day length and min/max data
out<-list()
for(i in 1:3){
  out[[i]]<-PL_func(Y=dayLength[i],Z=nightLength[i],Tx = 40, Tn = 10)
}

#seem like they all look good
plot(out[[3]])


#between two daily temperature scenarios
deltaForage_PL_dayLength <- function(input = out[[1]],forageMax,forageMin,
                                     tempChange=2, samples = 10000){
  y<-input
  yMod <- y+tempChange #apply the climate change scenario
  propMod <- sum(yMod < forageMax & yMod > forageMin) / samples #calculate niche overlap in new scneario
  prop <- sum(y < forageMax & y > forageMin) / samples #calculate niche overlap in old scenario
  c(prop*24*60, propMod*24*60, (propMod - prop) * 24*60) #calculate the difference, and convert to minutes per day
  
}

#new function seems to work
deltaForage_PL_dayLength(forageMax = 40, forageMin=-5)

#a loop inside to take each output from the deltaForage_PL function 
#(activity time, modified acitvity time, change in activity time)
forageMinVector <- seq(from = -5, to = 45, length.out =  400)
forageMaxVector <- seq(from = 5, to = 55, length.out = 400)

amountShiftDayLength <- array(dim = c(400, 3,365))
for(h in 1:3){
  for( i in 1:400){
    for(j in 1:3 ){
      amountShiftDayLength[i,j,h]  <-  deltaForage_PL_sagehen(input = out[[h]],forageMax=forageMaxVector[i], 
                                                              forageMin=forageMinVector[i], 
                                                              tempChange=2)[j]
    }
  }
}




day15night9<-as.data.frame(amountShiftDayLength[,,1])
day12night12<-as.data.frame(amountShiftDayLength[,,2])
day9night12<-as.data.frame(amountShiftDayLength[,,3])


day15night9$hoursDaylight<-"15"
day12night12$hoursDaylight<-"12"
day9night12$hoursDaylight<-"9"

colnames(day15night9)[3]<-"changeMinPerDay"
colnames(day12night12)[3]<-"changeMinPerDay"
colnames(day9night12)[3]<-"changeMinPerDay"

allDayLengths<-rbind(day15night9, day12night12, day9night12)

allDayLengths$forageTemp    <- (forageMinVector+forageMaxVector) / 2

plot(allDayLengths$changeMinPerDay)
p9 <- ggplot(allDayLengths, aes(x = forageTemp , y = changeMinPerDay, group = as.factor(hoursDaylight))) + 
  geom_hline(yintercept = 0, linetype =2) +geom_line(aes(color = as.factor(hoursDaylight)))+
  labs(y = "change in daily activity \\nduration after warming (minutes)", x = "median activity temperature (°C)") +
  theme_bw()+
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = Inf, ymin = 0, fill = "lightblue", alpha = 0.35)+
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = 0, ymin = -Inf, fill = "pink", alpha = 0.35)+
  #add the example species to the base simulation to enable comparison
  annotate("point" , x = allDayLengths[cumsum(allDayLengths$forageTemp>20)==1 , "forageTemp"],
           y = allDayLengths[cumsum(allDayLengths$forageTemp>20&allDayLengths$hoursDaylight=="12")==1 , "changeMinPerDay"], color = "lightyellow", cex = 8, pch = 16) +
  annotate("point" , x = amountShiftnorm[cumsum(allDayLengths$forageTemp>20)==1 , "forageTemp"],
           y = allDayLengths[cumsum(allDayLengths$forageTemp>20&allDayLengths$hoursDaylight=="12")==1 , "changeMinPerDay"], cex = 8, pch = 1) +
  annotate("text" , x = amountShiftnorm[cumsum(amountShiftnorm$forageTemp>20)==1 , "forageTemp"],
           y = allDayLengths[cumsum(allDayLengths$forageTemp>20&allDayLengths$hoursDaylight=="12")==1 , "changeMinPerDay"], label = "2")+
  
  annotate("point" , x = allDayLengths[cumsum(allDayLengths$forageTemp>33)==1 , "forageTemp"],
           y = allDayLengths[cumsum(allDayLengths$forageTemp>33&allDayLengths$hoursDaylight=="12")==1 , "changeMinPerDay"], color = "lightyellow", cex = 8, pch = 16) +
  annotate("point" , x = allDayLengths[cumsum(allDayLengths$forageTemp>33)==1 , "forageTemp"],
           y = allDayLengths[cumsum(allDayLengths$forageTemp>33&allDayLengths$hoursDaylight=="12")==1 , "changeMinPerDay"], cex = 8, pch = 1) +
  annotate("text" , x = allDayLengths[cumsum(amountShiftnorm$forageTemp>33)==1 , "forageTemp"],
           y = allDayLengths[cumsum(allDayLengths$forageTemp>33&allDayLengths$hoursDaylight=="12")==1 , "changeMinPerDay"], label = "1") +
  ggtitle("day length")+scale_color_brewer(palette = "Dark2")+
  labs(color = "hours of \\ndaylight")+
pdf("figures/dayLength_sensitivity.pdf",height = 6 ,width = 8 )
p9
dev.off()

################


#Supplemental Figure 5
#sensitivity - temperate vs tropical#####
#compare the effects of altering daily max-min (range)

#narrow foraging - 5 degree activity window
forageMinVector <- seq(from = -2.5, to = 47.5, length.out =  400)
forageMaxVector <- seq(from =2.5, to = 52.5, length.out = 400)

amountShiftTropic <- matrix(nrow = 400, ncol = 3)
start <- Sys.time()

#here we are looping through scenarios with narrow daily temperature ranges
for( i in 1:400){
  for(j in 1:3 ){
    amountShiftTropic[i,j]  <-  deltaForage_PL(dailyMax = 35, dailyMin = 15, forageMax=forageMaxVector[i], 
                                               forageMin=forageMinVector[i], dayToDayVar = 0 ,
                                               tempChangeNight=2, 
                                               tempChangeDay=2)[j]
  }
}
colnames(amountShiftTropic) <- c("forageDur" , "forageDurMod" , "diffForage")
amountShiftTropic<-as.data.frame(amountShiftTropic)
amountShiftTropic$forageTemp    <- (forageMinVector+forageMaxVector) / 2

#broad foraging - 20 degree foraging window
forageMinVector <- seq(from = -10, to = 40, length.out =  400)
forageMaxVector <- seq(from = 10, to = 60, length.out = 400)
amountShiftTemperate <- matrix(nrow = 400, ncol = 3)

#here we are looping through scenarios with broad daily temperature ranges
amountShiftTemperate <- matrix(nrow = 400, ncol = 3)
for( i in 1:400){
  for(j in 1:3 ){
    amountShiftTemperate[i,j]  <-  deltaForage_PL(dailyMax = 45, dailyMin = 5, forageMax=forageMaxVector[i], 
                                                  forageMin=forageMinVector[i], dayToDayVar = 0 ,
                                                  tempChangeNight=2, 
                                                  tempChangeDay=2)[j]
  }
}

colnames(amountShiftTemperate) <- c("forageDur" , "forageDurMod" , "diffForage")
amountShiftTemperate<-as.data.frame(amountShiftTemperate)
amountShiftTemperate$forageTemp    <- (forageMinVector+forageMaxVector) / 2

#these are labels for the different scenarios
amountShiftTropic$biome <- "tropic_5Cactive_20Cdayrange"
amountShiftTemperate$biome <- "temperate_20Cactive_40Cdayrange"
amountShift$biome <- "base_10Cactive_30Cdayrange"
#amountShift<-amountShift[ , -which(colnames(amountShift)=="variance")]
amountShiftBreadth<-rbind(amountShift , amountShiftTropic, amountShiftTemperate)

pdf("figures/biome_sensitivity.pdf",height = 5 ,width = 7 )
p6 <- ggplot(amountShiftBreadth, aes(x = forageTemp , y = diffForage , group = as.factor(biome))) + 
  geom_line(aes(color = as.factor(biome))) + geom_hline(yintercept = 0, linetype =2) + 
  labs(y = "change in activity duration \\nafter warming (minutes)", x = "median activity temperature (°C)") +
  theme_bw() +
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = Inf, ymin = 0, fill = "lightblue", alpha = 0.35)+
  annotate("rect" , xmin = -Inf , xmax = Inf, ymax = 0, ymin = -Inf, fill = "pink", alpha = 0.35)+
  ggtitle("daily temperature range")+
  labs(color = "biome")+ scale_color_brewer(palette = "Dark2")+
  
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], color = "lightyellow", cex = 8, pch = 16) +
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], cex = 8, pch = 1) +
  annotate("text" , x = amountShift[cumsum(amountShift$forageTemp>20)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>20)==1 , "diffForage"], label = "2")+
  
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], color = "lightyellow", cex = 8, pch = 16) +
  annotate("point" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], cex = 8, pch = 1) +
  annotate("text" , x = amountShift[cumsum(amountShift$forageTemp>33)==1 , "forageTemp"],
           y = amountShift[cumsum(amountShift$forageTemp>33)==1 , "diffForage"], label = "1")
p6     
dev.off()
amountShift<-amountShift[ , colnames(amountShift)!="dailyRange"]



#Figure 2 - all the senstivity analyses together
pdf("figures/all_sensitivity.pdf",height = 9 ,width = 18,onefile = FALSE )
#drop column for future reuse

grid.draw(ggarrange(p1, p2,p3,p4,p5, nrow = 2))
dev.off()

pdf("figures/breadth_sensitivity.pdf",height = 8 ,width = 8 )
p4 
#ant data
################

