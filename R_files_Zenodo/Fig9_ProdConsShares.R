########################################################################
#Figure 9: Production and consumption shares by sources with 
            #pessimistic interhour demand flexibility.

########################################################################
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
########################################################################


#directory <- c("outputs_theta01","outputs_theta05","outputs_theta2")
directory <- c("outputs_theta01")
for (dirfold0 in directory) {
  dirfold <- paste("./outputshpc/",dirfold0,"/",sep="")
  #setwd(dirfold)
  print(getwd())
  thetaval <- as.numeric(substr(dirfold0, 14,16))
  thetaval=ifelse(thetaval!=2,thetaval/10,thetaval)
  print(paste("working with theta:",thetaval))
  
  # evsq <- c("half", "full","2016")
  # loadsq <- c("2007", "2045")
  # categorysq <- c("rps_100","fossil", "free")
  # scenariosq <- seq(1,3)
  # costsq <- c("current", "future")
  # pricingsq <- c("flat","dynamic")
  
  evsq <- c("half")
  loadsq <- c("2045")
  categorysq <- c("rps_100","fossil", "free")
  scenariosq <- c(3)
  costsq <- c("future")
  pricingsq <- c("flat","dynamic")
  
  scalem <- 1000000
  
  for (cat in categorysq) {
    for (ev in evsq) {
      for(flatdyn in pricingsq) {
        for (load in loadsq) {
          for (cost in costsq) {
            for (scen in scenariosq){
              tryCatch({ 
                if (load=="2007" & ev!="half"){
                  print("Skip")
                } else {
                  
                  print(paste(thetaval,load,ev,cat,cost,flatdyn,scen))
                  
                  #gather final iteration in housely energy sources file
                  filename_e <- paste(dirfold,'/energy_sources_',cat,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
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
                  dfen <- transform(dfen, batnet=Battery-StorageNetCharge)
                  dfen <- transform(dfen, DischargeBattery=ifelse(batnet>=0, batnet,0))
                  dfen <- transform(dfen, ChargeBattery=ifelse(batnet<0, batnet*-1,0))
                  dfen <- transform(dfen, evexp=dfen$ChargeEVs*dfen$net.final.price)
                  dfen <- merge(dfen, w,by=c("month","day","hour"), all.x=TRUE )
                  
                  cat <- ifelse(cat=="rps","rps_100",cat) 

                  names(dfen)[names(dfen) == 'LiquifyHydrogenMW'] <- 'LiquifyHidrogren'
                  names(dfen)[names(dfen) == 'RunElectrolyzerMW'] <- 'ProduceHydrogen'
                  colnames(dfen)[colnames(dfen) == "ZoneTotalCentralDispatch"] <- "LZ_NetDispatch"
                  names(dfen)[names(dfen) == 'DispatchFuelCellMW'] <- 'FuelCell'
                  dfprice <- dfen %>%
                    summarize(across(c(final.price.energy), 
                                     list(w.mean = ~weighted.mean(., w = ts_scale_to_period), 
                                          sd = ~sd(., na.rm = TRUE),
                                          quant25 = ~quantile(., probs = 0.25),
                                          quant50 = ~quantile(., probs = 0.5),
                                          quant75 = ~quantile(., probs = 0.75),
                                          quant90 = ~quantile(., probs = 0.9)), .groups = 'drop'))
                  dfprice$theta=thetaval
                  dfprice$category=cat
                  dfprice$ev=ev
                  dfprice$load=load
                  dfprice$cost=cost
                  dfprice$pricing=flatdyn
                  dfprice$scen=scen
                  
                  dfgen <- dfen %>%
                    summarise(across(c(Biodiesel, Coal,Diesel,LNG,LSFO,Motor_Diesel,
                                       Motor_Gasoline,Pellet.Biomass,Battery,MSW,SUN,WND,
                                       curtail_Battery,curtail_MSW,curtail_SUN,curtail_WND,LZ_NetDispatch,
                                       FlexibleDemand,ChargeBattery,LiquifyHidrogren,ProduceHydrogen,ChargeEVs,FuelCell), 
                                     list( ~ sum(. * ts_scale_to_period)), .groups = 'drop'))
                  colnames(dfgen) <- c("Biodiesel", "Coal","Diesel","LNG","LSFO","Motor_Diesel",
                                       "Motor_Gasoline","Pellet.Biomass","Battery","MSW","SUN","WND",
                                       "curtail_Battery","curtail_MSW","curtail_SUN","curtail_WND","LZ_NetDispatch",
                                       "FlexibleDemand","ChargeBattery","LiquifyHidrogren","ProduceHydrogen","ChargeEVs","FuelCell")
                  dfgen$theta=thetaval
                  dfgen$category=cat
                  dfgen$ev=ev
                  dfgen$load=load
                  dfgen$cost=cost
                  dfgen$pricing=flatdyn
                  dfgen$scen=scen
                  ifelse(exists("d"), 
                         d <- rbind(d, dfgen),
                         d <- dfgen)
                  ifelse(exists("dp"), 
                         dp <- rbind(dp, dfprice),
                         dp <- dfprice)
                   write.csv(d, file = paste("data/gen",thetaval,".csv", sep = ""), row.names = FALSE)
                   write.csv(dp, file = paste("data/price",thetaval,".csv", sep = ""), row.names = FALSE)
                  
                }}, error=function(e){cat("ERROR :",conditionMessage(e),thetaval,load,ev,cat,cost,flatdyn,scen, "\\n")})
            }
          }
        }
      }
    }
  }
}

####share of production

{
x <- d
x$Curtailment <- x$curtail_MSW+x$curtail_SUN+x$curtail_WND
for (c in c("Curtailment","Biodiesel", "Coal","Diesel","LNG","LSFO","Motor_Diesel",
                  "Motor_Gasoline","Pellet.Biomass","Battery","MSW","SUN","WND","FuelCell")) {
x[[c]] <-x[[c]]/(x$LZ_NetDispatch+x$FuelCell)
}
names(x)[names(x) == 'WND'] <- 'Wind'
names(x)[names(x) == 'SUN'] <- 'SolarPV'
names(x)[names(x) == 'MSW'] <- 'H.Power'


x <- x[,c("theta","category","cost","scen","pricing","load","ev","Curtailment","Biodiesel", "Coal","Diesel","LNG","LSFO","Motor_Diesel",
          "Motor_Gasoline","Pellet.Biomass","Battery","H.Power","SolarPV","Wind","FuelCell")]
x <- melt(x, id=c("theta","category","cost","scen","pricing","load","ev"))

x <- transform(x, category=ifelse(x$category=="free","Unconstrained",paste(x$category)))
x <- transform(x, category=ifelse(x$category=="rps_100","100% Clean",paste(x$category)))
x <- transform(x, category=ifelse(x$category=="fossil","Fossil",paste(x$category)))
x <- transform(x, pricing=ifelse(x$pricing=="dynamic","RTP",paste(x$pricing)))
x <- transform(x, pricing=ifelse(x$pricing=="flat","Flat Price",paste(x$pricing)))
x <- transform(x, categoryord=ifelse(x$category=="100% Clean",3, ifelse(x$category=="Fossil",1,2)))
x <- transform(x, category = reorder(category, categoryord))

ggplot(x[x$scen==3,], aes(x=pricing,y=value, fill = variable)) + 
  geom_bar(stat = "identity")+ theme_classic()+
  facet_grid(~category, scale = "free_y", space = "fixed") +
  xlab("") +
  ylab("Share of total production") +
  scale_fill_manual(name="Source",
                    values=c(Diesel="#999999", Coal="#000000", LSFO="#CCCCCC", LNG="#666666", H.Power="#CC9900", Biodiesel="#99CC00", 
                             Pellet.Biomass="#003300", Wind="#3399FF", SolarPV="#FFFF00", 
                             Battery="#990000",FuelCell="#DDA0DD",Curtailment="#FFFF99"))+
  scale_y_continuous(n.breaks = 6)
ggsave(paste("plot/Fig9a_Generators50EVFuturePessimistictheta01.pdf", sep = ""),height=4.1,width=8)


####share of demand
x <- d
x$totdemand <- x$ProduceHydrogen+x$ChargeEVs+x$LiquifyHidrogren+x$ChargeBattery+x$FlexibleDemand
for (c in c("ProduceHydrogen","ChargeEVs", 
            "LiquifyHidrogren", "ChargeBattery", "FlexibleDemand")) {
  x[[c]] <-x[[c]]/x$totdemand
}

x <- x[,c("theta","category","cost","scen","pricing","load","ev","ProduceHydrogen","ChargeEVs", 
          "LiquifyHidrogren", "ChargeBattery", "FlexibleDemand")]
x <- melt(x, id=c("theta","category","cost","scen","pricing","load","ev"))

x <- transform(x, category=ifelse(x$category=="free","Unconstrained",paste(x$category)))
x <- transform(x, category=ifelse(x$category=="rps_100","100% Clean",paste(x$category)))
x <- transform(x, category=ifelse(x$category=="fossil","Fossil",paste(x$category)))
x <- transform(x, pricing=ifelse(x$pricing=="dynamic","RTP",paste(x$pricing)))
x <- transform(x, pricing=ifelse(x$pricing=="flat","Flat Price",paste(x$pricing)))
x <- transform(x, categoryord=ifelse(x$category=="100% Clean",3, ifelse(x$category=="Fossil",1,1)))
x <- transform(x, category = reorder(category, categoryord))

ggplot(x[x$scen==3,], aes(x=pricing,y=value, fill = variable)) + 
  geom_bar(stat = "identity")+ theme_classic()+
  facet_grid(~category, scale = "free_y", space = "fixed") +
  xlab("") +
  ylab("Share of total consumption") +
  scale_fill_manual(name="Source",
                    values=c(FlexibleDemand="#CCCCCC",
                             ChargeBattery="#990000",
                             LiquifyHidrogren="#FF9933",
                             ChargeEVs="#666600",
                             ProduceHydrogen="#0033FF"))+
  scale_y_continuous(n.breaks = 6)
ggsave(paste("plot/Fig9b_Demand50EVFuturePessimistictheta01.pdf", sep = ""),height=4.1,width=8)
}
dev.off()
rm(d)
rm(dp)
