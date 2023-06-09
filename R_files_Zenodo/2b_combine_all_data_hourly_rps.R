#NOT USED
#made with R version 4.0.3 (2020-10-10)
## modified on 14.03.2022 with modified ev charging cost
# produce energysourcesCO2",thetaval,"_rps.csv" for FIG8

########################################################################
#loads graph
setwd("./inputs")
filename <- 'timepoints.csv'
timelabeldf <- read.csv(filename, header=TRUE, sep=",")
filename <- 'loads.csv'
loadsdf <- read.csv(filename, header=TRUE, sep=",")
names(loadsdf)[names(loadsdf) == 'TIMEPOINT'] <- 'timepoint_id'
timelabeldf <- merge(timelabeldf, loadsdf, by=c("timepoint_id"))

filename <- 'timeseries.csv'
sampleweightdf <- read.csv(filename, header=TRUE, sep=",")
sampleweightdf$samplew <- round(sampleweightdf$ts_scale_to_period/sum(sampleweightdf$ts_scale_to_period),2)
names(sampleweightdf)[names(sampleweightdf) == 'TIMESERIES'] <- 'timeseries'

timelabeldf <- merge(timelabeldf, sampleweightdf, by=c("timeseries"))
timelabeldf$hour <- hour(timelabeldf$timestamp)
timelabeldf$month <- month(timelabeldf$timestamp)
timelabeldf$day <- day(timelabeldf$timestamp)
w <- timelabeldf[,c("timeseries","ts_scale_to_period","month","day","hour","samplew")]
print("load input data")
setwd('..') 
###############################################################################

directory <- c("outputs_theta01","outputs_theta05","outputs_theta2")
collist <- c("month"	,	"day"	,
             "hour"	,	"load_zone"	,
             "period"	,	"timepoint"	,
             "Biodiesel"	,	"Coal"	,
             "Diesel"	,	"LNG"	,
             "LSFO"	,	"Pellet.Biomass"	,
             "Battery"	,	"MSW"	,
             "SUN"	,	"WND"	,
             "curtail_Battery"	,	"curtail_MSW"	,
             "curtail_SUN"	,	"curtail_WND"	,
             "ZoneTotalCentralDispatch"	,	"ZoneTotalDistributedDispatch"	,
              "DRUnservedLoad"	,	"DispatchFuelCellMW"	,
              "FlexibleDemand"	,	"StorageNetCharge"	,
              "ChargeEVs"	,	"RunElectrolyzerMW"	,
              "LiquifyHydrogenMW"	,	"offered.price.energy"	,
             "final.price.energy", "final.q.energy"	,
              "final.q.energy.up"	,	"final.q.energy.down"	,
              "peak_day"	,	"base_load"	,
              "base_price"	,	"net.final.mc"	,
              "net.final.price"	,	"net.final.q"	,
              "timeseries"	,	"ts_scale_to_period"	,
              "samplew"	,	"theta"	,
              "category"	,	"ev"	,
              "load"	,	"cost"	,
              "pricing"	,	"scen"	,
              "renshare")

for (dirfold in directory) {
  dirfold0 <- paste("./outputshpc/",dirfold,"/",sep="")
  #setwd(dirfold)
  print(getwd())
  thetaval <- as.numeric(substr(dirfold, 14,16))
  thetaval=ifelse(thetaval!=2,thetaval/10,thetaval)
  print(paste("working with theta:",thetaval))
  
  evsq <- c("half", "full","2016")
  loadsq <- c("2045")
  categorysq <- c("rps")
  scenariosq <- seq(1,3)
  costsq <- c("current", "future")
  pricingsq <- c("flat","dynamic")
  scalem <- 1000000
  
  rpstargetsq <- seq(10,80,by=10)
  rpstargetsq<- c(rpstargetsq,seq(82,98,by=2))
  rpstargetsq <- paste("0",rpstargetsq,sep="") 
  rpstargetsq<- c("000",rpstargetsq,"100")
  
  print(rpstargetsq)
  
  if (file.exists(paste("data/energysourcesCO2",thetaval,"_rps.csv", sep = "")))
    recent <- read.csv(paste("data/energysourcesCO2",thetaval,"_rps.csv", sep = ""), header=TRUE, sep=",")
  
  for (cat in categorysq) {
    for (ev in evsq) {
      for(flatdyn in pricingsq) {
        for (load in loadsq) {
          for (cost in costsq) {
            for (scen in scenariosq){
              for (rpstarget in rpstargetsq){
                if (exists("recent")){
                  if (any(recent$renshare==rpstarget&recent$category==cat&recent$ev==ev&recent$pricing==flatdyn&recent$load==load&recent$cost==cost&recent$scen==scen)){
                    print(paste("Already run",thetaval,load,ev,cat,cost,flatdyn,scen))
                  } else {
                    if ((load=="2007" & ev!="half")){
                      print(paste("Skip",thetaval,load,ev,cat,cost,flatdyn,scen))
                    } else {
                      tryCatch({ 
                        print(paste(rpstarget,thetaval,load,ev,cat,cost,flatdyn,scen))
  
                        #gather final iteration in housely energy sources file
                        filename_e <- paste(dirfold0,'energy_sources_',cat,'_',rpstarget,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
                        dfen <- read.csv(filename_e, header=TRUE, sep=",")
                        colnames(dfen)[colnames(dfen) == "timepoint_label"] <- "timepoint"
                        dfen <- transform(dfen, hour=as.numeric(format(
                          strptime(dfen$timepoint, format="%Y-%m-%d %H:%M"), format="%H")))
                        dfen <- transform(dfen,month =as.numeric(format(
                          strptime(dfen$timepoint, format="%Y-%m-%d %H:%M"), format="%m")))
                        dfen <- transform(dfen,day =as.numeric(format(
                          strptime(dfen$timepoint, format="%Y-%m-%d %H:%M"), format="%d")))
                        dfen <- transform(dfen, net.final.mc=final.mc.energy+final.mc.energy.down-final.mc.energy.up)
                        dfen <- transform(dfen, net.final.price=offered.price.energy+offered.price.energy.down-offered.price.energy.up)
                        dfen <- transform(dfen, net.final.q=final.q.energy) 
                        dfenb <- merge(dfen, w,by=c("month","day","hour"), all.x=TRUE )

                        dfenb <- transform(dfenb, theta=thetaval)
                        dfenb <- transform(dfenb, category=cat)
                        dfenb <- transform(dfenb, ev=ev)
                        dfenb <- transform(dfenb, load=load)
                        dfenb <- transform(dfenb, cost=cost)
                        dfenb <- transform(dfenb, pricing=flatdyn)
                        dfenb <- transform(dfenb, scen=scen)
                        dfenb <- transform(dfenb, renshare=as.character(rpstarget))
                        dfenb <- dfenb[,collist]
                        ifelse(exists("recent"), 
                               recent <- rbind(recent, dfenb),
                               recent <- dfenb)
                        
                        write.csv(recent, file = paste("data/energysourcesCO2",thetaval,"_rps.csv", sep = ""), row.names = FALSE)
                        
                      }, error=function(e){cat("ERROR :",conditionMessage(e),thetaval,load,ev,cat,cost,flatdyn,scen, "\\n")})
                    }
                  }
                } else {
                  #similar as above
                  if ((load=="2007" & ev!="half")){
                    print(paste("Skip",thetaval,load,ev,cat,cost,flatdyn,scen))
                  } else {
                    tryCatch({ 
                      print(paste(rpstarget,thetaval,load,ev,cat,cost,flatdyn,scen))

                      #gather final iteration in housely energy sources file
                      filename_e <- paste(dirfold0,'energy_sources_',cat,'_',rpstarget,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
                      dfen <- read.csv(filename_e, header=TRUE, sep=",")
                      colnames(dfen)[colnames(dfen) == "timepoint_label"] <- "timepoint"
                      dfen <- transform(dfen, hour=as.numeric(format(
                        strptime(dfen$timepoint, format="%Y-%m-%d %H:%M"), format="%H")))
                      dfen <- transform(dfen,month =as.numeric(format(
                        strptime(dfen$timepoint, format="%Y-%m-%d %H:%M"), format="%m")))
                      dfen <- transform(dfen,day =as.numeric(format(
                        strptime(dfen$timepoint, format="%Y-%m-%d %H:%M"), format="%d")))
                      dfen <- transform(dfen, net.final.mc=final.mc.energy+final.mc.energy.down-final.mc.energy.up)
                      dfen <- transform(dfen, net.final.price=offered.price.energy+offered.price.energy.down-offered.price.energy.up)
                      dfen <- transform(dfen, net.final.q=final.q.energy) 
                      dfenb <- merge(dfen, w,by=c("month","day","hour"), all.x=TRUE )
                      
                      dfenb <- transform(dfenb, theta=thetaval)
                      dfenb <- transform(dfenb, category=cat)
                      dfenb <- transform(dfenb, ev=ev)
                      dfenb <- transform(dfenb, load=load)
                      dfenb <- transform(dfenb, cost=cost)
                      dfenb <- transform(dfenb, pricing=flatdyn)
                      dfenb <- transform(dfenb, scen=scen)
                      dfenb <- transform(dfenb, renshare=as.character(rpstarget))
                      dfenb <- dfenb[,collist]
                      ifelse(exists("recent"), 
                             recent <- rbind(recent, dfenb),
                             recent <- dfenb)
                      
                      write.csv(recent, file = paste("data/energysourcesCO2",thetaval,"_rps.csv", sep = ""), row.names = FALSE)
                      
                    }, error=function(e){cat("ERROR :",conditionMessage(e),thetaval,load,ev,cat,cost,flatdyn,scen, "\\n")})
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  rm(recent)
}  
