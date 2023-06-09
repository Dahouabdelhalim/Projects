#######################################################################
#Table 5: Main Results: Comparison of prices, quantities, and surplus with at and RTP pricing.
#Figure 4: Surplus gain from real time pricing under different policy, cost and demand flexibility scenarios.
#Figure 5: Cost of 100 percent renewable energy system under different policy, cost and demand flexibility scenarios.
# produce basesum*csv. for Figure 6
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
w <- timelabeldf[,c("month","day","hour","ts_scale_to_period")]
#w <- timelabeldf[which(timelabeldf$hour!=0),c("month","day","hour","samplew")]
setwd('..') 
###############################################################################
#baseline expenditure is around $1.3 billion
baselineexp <- 1300000000
directory <- c("./outputs_theta01")
for (dirfold in directory) {
  #setwd(dirfold)
  print(getwd())
  thetaval <- as.numeric(substr(dirfold, 16,18))
  thetaval=ifelse(thetaval!=2,thetaval/10,thetaval)
  
  recent <- read.csv(paste("data/energysourcesCO2",thetaval,".csv", sep = ""), header=TRUE, sep=",")
  recent <- transform(recent, days=days_in_month(as.Date(paste('2045-',recent$month,'-','1',sep=""),"%Y-%m-%d")))
  recent <- transform(recent, tots=-tots)
  recent <- merge(recent,w,by=c("month","day","hour"))
  
  ###cs by group###
  x<-recent
  idlist <- list(x$theta,x$ev,x$load,x$category, x$scen, x$cost, x$pricing, x$month, x$days, x$ts_scale_to_period.x,x$csmonth,x$evmonth)
  #take a mean because it is actually daily data but copied 24 times to match 24 hours dimension
  x0 <- aggregate(x[,c("evexp")], idlist, sum)
  colnames(x0) <- c("theta","ev","load","category","scen","cost","pricing","month", "days",
                    "samplew","csmonth","evmonth",
                    "evexp")
  x <- aggregate(x[,c("wtphighflex","wtpmidflex","wtpinflex",
                      "wtptot","pqhighflex","pqmidflex","pqinflex")], idlist, mean)
  colnames(x) <- c("theta","ev","load","category","scen","cost","pricing","month", "days",
                   "samplew","csmonth","evmonth",
                   "highflex","midflex","inflex","wtptot",
                   "highexp","midexp","infexp")
  bylist <- c("theta","ev","load","category","scen","cost","pricing","month", "days",
              "samplew","csmonth","evmonth")
  x <- merge(x0, x, by=c(bylist))
  
  x <- transform(x,highflex=(x$highflex*x$samplew))
  x <- transform(x,midflex=(x$midflex*x$samplew))
  x <- transform(x,inflex=(x$inflex*x$samplew))
  x <- transform(x,wtptot=(x$wtptot*x$samplew))
  
  x <- transform(x,highexp=(x$highexp*x$samplew))
  x <- transform(x,midexp=(x$midexp*x$samplew))
  x <- transform(x,infexp=(x$infexp*x$samplew))
  x <- transform(x,totexp=(x$highexp+x$midexp+x$infexp))
  x <- transform(x,evexp=(x$evexp*x$samplew))
  
  #x <- transform(x,csmonth=(x$csmonth*x$samplew))
  #x <- transform(x,evmonth=(x$evmonth*x$samplew))
  
  idlist <- list(x$theta,x$ev,x$load,x$category, x$scen, x$cost, x$pricing)
  #take a sum because we want to get the total
  x <- aggregate(x[,c("csmonth","evmonth","evexp","highflex","midflex","inflex",
                      "wtptot","highexp","midexp","infexp","totexp")], idlist, sum)
  colnames(x) <- c("theta","ev","load","category","scen","cost","pricing","csmonth","evmonth","evexp",
                   "cshighflex","csmidflex","csinflex","cstot","highexp","midexp","infexp","totexp")
  x <- transform(x, totcs=csmonth+evmonth)
  x1 <- x[,c("theta","ev","load","category","scen","cost","pricing","cshighflex",
             "csmidflex","csinflex","cstot","highexp","midexp","infexp","totexp","evexp")]
  
  #dbaseexp <- merge(dbaseexp, x1[which(x1$cost=="future" & x1$category=="fossil" & x1$pricing=="flat"),], by=c("scen"), all.x=TRUE, all.y=TRUE)
  dbaseexp <- x1[which(x1$cost=="future" & x1$category=="fossil" & x1$pricing=="flat"),]
  names(dbaseexp)[names(dbaseexp) == 'evexp'] <- 'evexp1'
  
  ########create latex summary result to get producer surplus##############      
  
  ########################################
  #REPORT TABLE# for basic case "half EV and 2045 projected load"
  ########################################
  # 1. Fossil - optimistic - current - flat
  # 2. Fossil - optimistic - current - variable
  # 3. Fossil - optimistic - future - flat
  # 4. Fossil - optimistic - future - variable
  #(a) % renewable; (b) average price; (c) price SD; (c) average load; \\Delta PS; \\Delta CS; \\Delta TS;
  
  dfen <- recent
  #we do not need to scale those two because they are already scaled properly internally in SWITCH
  #dfen <- transform(dfen,csmonth=(dfen$csmonth*dfen$ts_scale_to_period.x))
  #dfen <- transform(dfen,evmonth=(dfen$evmonth*dfen$ts_scale_to_period.x))
  
  idlist <- list(dfen$theta,dfen$ev,dfen$load,dfen$category, dfen$scen, dfen$cost, dfen$pricing, dfen$month)
  dfensumqmonth <- aggregate(dfen[,c("net.final.q")], idlist, mean)
  
  colnames(dfensumqmonth) <- c("theta","ev","load","category","scen","cost","pricing","month", "sumoq")
  pjoin <- dfen
  idlist <- list(pjoin$theta,pjoin$ev,pjoin$load,pjoin$category, pjoin$scen, pjoin$cost, pjoin$pricing)
  bylist <- c("Group.1","Group.2","Group.3","Group.4","Group.5","Group.6","Group.7")
  pjoin <- merge(aggregate(pjoin[,c("net.final.price")], idlist, mean),
                 aggregate(pjoin[,c("net.final.price")], idlist, sd), by=c(bylist))
  idlist <- list(dfensumqmonth$theta,dfensumqmonth$ev,dfensumqmonth$load,dfensumqmonth$category, dfensumqmonth$scen, dfensumqmonth$cost, dfensumqmonth$pricing)
  pjoin <- merge(pjoin,
                 aggregate(dfensumqmonth[,c("sumoq")], idlist, mean), 
                 by=c(bylist))
  colnames(pjoin) <- c("theta","ev","load","category","scen","cost","pricing","mean", "sd","meanq")
  idlist <- list(dfen$theta,dfen$ev,dfen$load,dfen$category, dfen$scen, dfen$cost, dfen$pricing, dfen$month, dfen$day)
  pjoinevcs <- aggregate(dfen[,c("csmonth","evmonth")], idlist, mean)
  colnames(pjoinevcs) <- c("theta","ev","load","category","scen","cost","pricing","month","day","csmonth","evmonth")
  idlist <- list(pjoinevcs$theta,pjoinevcs$ev,pjoinevcs$load,pjoinevcs$category, pjoinevcs$scen, pjoinevcs$cost, pjoinevcs$pricing)
  pjoinevcs1 <- aggregate(pjoinevcs[,c("csmonth","evmonth")], idlist, sum)
  colnames(pjoinevcs1) <- c("theta","ev","load","category","scen","cost","pricing","csmonth","evmonth")
  bylist <- c("theta","ev","load","category","scen","cost","pricing")
  pjoin <- merge(pjoin, pjoinevcs1,by=c(bylist))
  
  idlist <- list(dfen$theta,dfen$ev,dfen$load,dfen$category, dfen$scen, dfen$cost, dfen$pricing)
  pjoinerenshare <- aggregate(dfen[,c("renshare","tots","co2")], idlist, mean)
  colnames(pjoinerenshare) <- c("theta","ev","load","category","scen","cost","pricing","Renewable_Share","tots","co2")
  
  ###################################################################
  # % expenditure
  joinall <- merge(pjoinerenshare[,c("theta","ev","load","category","scen","cost","pricing","Renewable_Share","tots","co2")],
                   pjoin[,c("theta","ev","load","category","scen","cost","pricing","mean","meanq","sd","csmonth","evmonth")],
                   by=c("theta","ev","load","category","scen","cost","pricing"), all.y=TRUE)
  #joinall <- merge(joinall,dfj[,c("theta","ev","load","category","scen","cost","pricing","tots")], by=c("theta","ev","load","category","scen","cost","pricing"))
  joinall <- transform(joinall, evmonth=-evmonth)
  joinall <- transform(joinall, csmonth=-csmonth)
  joinall <- transform(joinall, dps=tots-evmonth-csmonth)
  joinall <- transform(joinall, csmonth1=evmonth+csmonth)
  joinall <- subset(joinall, select=-c(csmonth))
  joinall <- transform(joinall, csmonth=csmonth1)
  joinall <- subset(joinall, select=-c(csmonth1))
  joinall <- merge(joinall, x1, by=c("theta","ev","load","category","scen","cost","pricing"))
  #combinetable <- joinall
  joinall <- joinall[,colnames(joinall)!='evexp']
  
  #baseline is fossil future flat
  basesum <- joinall[which(joinall$cost=="future" & joinall$category=="fossil" & joinall$pricing=="flat"),]
  basesum <- merge(basesum[,c("theta","ev","load","scen","csmonth","evmonth","tots","dps",
                              "cshighflex","csmidflex","csinflex",
                              "highexp","midexp","infexp","totexp","co2")], 
                   dbaseexp[which(dbaseexp$cost=="future" & dbaseexp$category=="fossil" & dbaseexp$pricing=="flat"),
                            c("theta","ev","load","scen","evexp1")], by=c("theta","ev","load","scen"), all.x=TRUE, all.y=FALSE)
  colnames(basesum) <- c("theta","ev","load","scen","csmonthb","evmonthb","totsb","dpsb",
                         "cshighflexb","csmidflexb","csinflexb",
                         "highexpb","midexpb","infexpb","totexpb","co2","evexp")
  
  basesum$highexpb <- (basesum$highexpb/basesum$totexpb)*(baselineexp+basesum$evexp)
  basesum$midexpb <- (basesum$midexpb/basesum$totexpb)*(baselineexp+basesum$evexp)
  basesum$infexpb <- (basesum$infexpb/basesum$totexpb)*(baselineexp+basesum$evexp)
  basesum$totexpb <- baselineexp+basesum$evexp
  #setwd('..') 
  write.csv(basesum, file = paste("data/basesum",thetaval,".csv", sep = ""), row.names = FALSE)
  
  joinall <- merge(joinall, basesum, by=c("theta","ev","load","scen"), all.x=TRUE, all.y=TRUE)
  joinall <- transform(joinall, dcs = 1*round((csmonth-csmonthb)/totexpb*100,1))
  joinall <- transform(joinall, evcost = -1*round((evmonth-evmonthb)/evexp*100,1))
  joinall <- transform(joinall, dps = 1*round((dps-dpsb)/totexpb*100,1))
  joinall <- transform(joinall, tots =1*round((tots-totsb)/totexpb*100,1))
  joinall <- transform(joinall, high =1*round((cshighflex-cshighflexb)/highexpb*100,1))
  joinall <- transform(joinall, mid =1*round((csmidflex-csmidflexb)/midexpb*100,1))
  joinall <- transform(joinall, inf =1*round((csinflex-csinflexb)/infexpb*100,1))
  
  #joinall <- transform(joinall, tot3=ifelse(high!=0,inf+high+mid,0))
  #joinall <- transform(joinall, inf =ifelse(high!=0,(inf/tot3) * dcs,0))
  #joinall <- transform(joinall, high =ifelse(high!=0,(high/tot3) * dcs,0))
  #joinall <- transform(joinall, mid =ifelse(high!=0,(mid/tot3) * dcs,0))
  joinall <- transform(joinall, mean = round(mean,0))
  joinall <- transform(joinall, meanq = round(meanq,0))
  joinall <- transform(joinall, sd = round(sd,0))
  write.csv(joinall, file = paste("table_",thetaval,".csv", sep = ""), row.names = FALSE)
  joinall <- joinall[,c("theta","ev","load","category","scen","cost","pricing","Renewable_Share",
                        "mean","meanq","sd","dcs","evcost","dps","tots","high","mid","inf")]
  maintable <- joinall
  #value of dynamic for total surplus
  dynbase <- joinall[which(joinall$pricing=="flat"),c("theta","ev","load","category","scen","cost","tots")]
  colnames(dynbase) <- c("theta","ev","load","category","scen","cost","totsb")
  joinall <- merge(joinall, dynbase, by=c("theta","ev","load","category","scen","cost"), all.x=TRUE, all.y=FALSE)
  joinall <- transform(joinall, dyntot = round((tots-totsb),2))
  joinall <- joinall[,c("theta","ev","load","category","cost","scen","pricing",
                        "Renewable_Share","mean","meanq","sd","dcs","evcost","dps","tots","high","mid","inf","dyntot")]
  
  #value of dynamic for consumer
  dyncons <- maintable[which(maintable$pricing=="flat"),c("theta","ev","load","category","scen","cost","dcs")]
  colnames(dyncons) <- c("theta","ev","load","category","scen","cost","dcsb")
  dyncons <- merge(maintable, dyncons, by=c("theta","ev","load","category","scen","cost"), all.x=TRUE, all.y=FALSE)
  dyncons <- transform(dyncons, dyncons = round((dcs-dcsb),2))
  dyncons <- dyncons[,c("theta","ev","load","category","cost","scen","pricing",
                        "Renewable_Share","mean","meanq","sd","dcs","evcost","dps","tots","high","mid","inf","dyncons")]
  
  #value of dynamic for producer
  dynprod <- maintable[which(maintable$pricing=="flat"),c("theta","ev","load","category","scen","cost","dps")]
  colnames(dynprod) <- c("theta","ev","load","category","scen","cost","dpsb")
  dynprod <- merge(maintable, dynprod, by=c("theta","ev","load","category","scen","cost"), all.x=TRUE, all.y=FALSE)
  dynprod <- transform(dynprod, dynprod = round((dps-dpsb),2))
  dynprod <- dynprod[,c("theta","ev","load","category","cost","scen","pricing",
                        "Renewable_Share","mean","meanq","sd","dcs","evcost","dps","tots","high","mid","inf","dynprod")]
  
  graphdata <- joinall
  
  #write.csv(graphdata, file = paste("tables/maintable_theta=",thetaval,".csv", sep = ""), row.names = FALSE)
  
  forlatex <- graphdata
  #naming the figures
  num1 <- 4
  num2 <- 5
  source(file = "Fig4-5.R")
  source(file = "Tab5.R")
  source(file = "Fig7_DistPricesQ.R")
  #dev.off()
  
  ##############################################################################
  ##APPENDIX FigS1a Share of hours with Marginal Cost < 5 or <1 cents/KWh
  ##############################################################################
  pricebins <- transform(recent, mc50=ifelse(net.final.mc<50, 1,0))
  pricebins <- transform(pricebins, mc10=ifelse(net.final.mc<10, 1,0))
  idlist <- list(pricebins$category, pricebins$scen, pricebins$cost,pricebins$pricing,pricebins$load,pricebins$ev)
  pricebins <- aggregate(list(pricebins$mc50,pricebins$mc10), idlist, mean)
  colnames(pricebins) <- c("category","scen","cost","pricing","load","ev","sharemc50","sharemc10")
  pricebins <- transform(pricebins, category=ifelse(pricebins$category=="free","Unconstrained",paste(pricebins$category)))
  pricebins <- transform(pricebins, category=ifelse(pricebins$category=="rps_100","100% Clean",paste(pricebins$category)))
  pricebins <- transform(pricebins, category=ifelse(pricebins$category=="fossil","Fossil",paste(pricebins$category)))
  pricebins <- transform(pricebins, Demand_flexibility=ifelse(pricebins$scen==1,"1-Optimistic", ifelse(pricebins$scen==2,"2-Moderate","3-Pessimistic")))
  pricebins <- transform(pricebins, pricing=ifelse(pricebins$pricing=="dynamic","Real Time",paste(pricebins$pricing)))
  pricebins <- transform(pricebins, pricing=ifelse(pricebins$pricing=="flat","Flat",paste(pricebins$pricing)))
  pricebins <- transform(pricebins, pricing=paste(pricebins$pricing, "pricing"))
  pricebins <- transform(pricebins, cost=ifelse(pricebins$cost=="current","2016 Cost","2045 Cost"))
  pricebins <- transform(pricebins, categoryord=ifelse(pricebins$category=="100% Clean",2, ifelse(pricebins$category=="Fossil",1,3)))
  pricebins <- transform(pricebins, category = reorder(category, categoryord))
  
  ggplot(data = pricebins) + 
    geom_point(aes(category, sharemc10,color="olivedrab4"), size=6, position = "jitter", show.legend=F)+
    geom_point(aes(category, sharemc50,color="darkseagreen2"), size=4, position = "jitter")+
    facet_grid(cost~ pricing, scale = "free_y", space = "fixed") +
    scale_size_area() + 
    xlab("") +
    ylab("Share of hours with marginal cost") + theme_classic()+
    scale_color_manual("Price",values=c("olivedrab4","darkseagreen2"),labels=c("<1 cents/KWh","<5 cents/KWh"),
                       guide=guide_legend(override.aes = list(color=c("darkseagreen2","olivedrab4"))))+
    theme(legend.direction="horizontal",legend.position = "bottom", 
          legend.title = element_text(colour = "black", size = 15),
          legend.text = element_text(colour = "black", size = 15),
          axis.ticks=element_blank(), legend.key.width = unit(0.5,"cm"), 
          plot.title = element_text(hjust = -10),
          axis.text=element_text(size=13, angle = 0),
          axis.title=element_text(size=15),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15))+
    geom_hline(yintercept=0, col = "gray80")
  ggsave(paste("plot/FigS1a_MarginalCost_theta",thetaval,".pdf",sep=""),height=6.1,width=10)
}


################################################################################
##FOR APPENDIX
###############################################################################
#baseline expenditure is around $1.3 billion
baselineexp <- 1300000000
directory <- c("./outputs_theta05","./outputs_theta2")
for (dirfold in directory) {
  #setwd(dirfold)
  print(getwd())
  thetaval <- as.numeric(substr(dirfold, 16,18))
  thetaval=ifelse(thetaval!=2,thetaval/10,thetaval)
  
  recent <- read.csv(paste("data/energysourcesCO2",thetaval,".csv", sep = ""), header=TRUE, sep=",")
  recent <- transform(recent, days=days_in_month(as.Date(paste('2045-',recent$month,'-','1',sep=""),"%Y-%m-%d")))
  recent <- transform(recent, tots=-tots)
  recent <- merge(recent,w,by=c("month","day","hour"))
  
  ###cs by group###
  x<-recent
  idlist <- list(x$theta,x$ev,x$load,x$category, x$scen, x$cost, x$pricing, x$month, x$days, x$ts_scale_to_period.x,x$csmonth,x$evmonth)
  #take a mean because it is actually daily data but copied 24 times to match 24 hours dimension
  x0 <- aggregate(x[,c("evexp")], idlist, sum)
  colnames(x0) <- c("theta","ev","load","category","scen","cost","pricing","month", "days",
                    "samplew","csmonth","evmonth",
                    "evexp")
  x <- aggregate(x[,c("wtphighflex","wtpmidflex","wtpinflex",
                      "wtptot","pqhighflex","pqmidflex","pqinflex")], idlist, mean)
  colnames(x) <- c("theta","ev","load","category","scen","cost","pricing","month", "days",
                   "samplew","csmonth","evmonth",
                   "highflex","midflex","inflex","wtptot",
                   "highexp","midexp","infexp")
  bylist <- c("theta","ev","load","category","scen","cost","pricing","month", "days",
              "samplew","csmonth","evmonth")
  x <- merge(x0, x, by=c(bylist))
  
  x <- transform(x,highflex=(x$highflex*x$samplew))
  x <- transform(x,midflex=(x$midflex*x$samplew))
  x <- transform(x,inflex=(x$inflex*x$samplew))
  x <- transform(x,wtptot=(x$wtptot*x$samplew))
  
  x <- transform(x,highexp=(x$highexp*x$samplew))
  x <- transform(x,midexp=(x$midexp*x$samplew))
  x <- transform(x,infexp=(x$infexp*x$samplew))
  x <- transform(x,totexp=(x$highexp+x$midexp+x$infexp))
  x <- transform(x,evexp=(x$evexp*x$samplew))
  
  #x <- transform(x,csmonth=(x$csmonth*x$samplew))
  #x <- transform(x,evmonth=(x$evmonth*x$samplew))
  
  idlist <- list(x$theta,x$ev,x$load,x$category, x$scen, x$cost, x$pricing)
  #take a sum because we want to get the total
  x <- aggregate(x[,c("csmonth","evmonth","evexp","highflex","midflex","inflex",
                      "wtptot","highexp","midexp","infexp","totexp")], idlist, sum)
  colnames(x) <- c("theta","ev","load","category","scen","cost","pricing","csmonth","evmonth","evexp",
                   "cshighflex","csmidflex","csinflex","cstot","highexp","midexp","infexp","totexp")
  x <- transform(x, totcs=csmonth+evmonth)
  x1 <- x[,c("theta","ev","load","category","scen","cost","pricing","cshighflex",
             "csmidflex","csinflex","cstot","highexp","midexp","infexp","totexp","evexp")]
  
  #dbaseexp <- merge(dbaseexp, x1[which(x1$cost=="future" & x1$category=="fossil" & x1$pricing=="flat"),], by=c("scen"), all.x=TRUE, all.y=TRUE)
  dbaseexp <- x1[which(x1$cost=="future" & x1$category=="fossil" & x1$pricing=="flat"),]
  names(dbaseexp)[names(dbaseexp) == 'evexp'] <- 'evexp1'
  
  ########create latex summary result to get producer surplus##############      
  
  ########################################
  #REPORT TABLE# for basic case "half EV and 2045 projected load"
  ########################################
  # 1. Fossil - optimistic - current - flat
  # 2. Fossil - optimistic - current - variable
  # 3. Fossil - optimistic - future - flat
  # 4. Fossil - optimistic - future - variable
  #(a) % renewable; (b) average price; (c) price SD; (c) average load; \\Delta PS; \\Delta CS; \\Delta TS;
  
  dfen <- recent
  #we do not need to scale those two because they are already scaled properly internally in SWITCH
  #dfen <- transform(dfen,csmonth=(dfen$csmonth*dfen$ts_scale_to_period.x))
  #dfen <- transform(dfen,evmonth=(dfen$evmonth*dfen$ts_scale_to_period.x))
  
  idlist <- list(dfen$theta,dfen$ev,dfen$load,dfen$category, dfen$scen, dfen$cost, dfen$pricing, dfen$month)
  dfensumqmonth <- aggregate(dfen[,c("net.final.q")], idlist, mean)
  
  colnames(dfensumqmonth) <- c("theta","ev","load","category","scen","cost","pricing","month", "sumoq")
  pjoin <- dfen
  idlist <- list(pjoin$theta,pjoin$ev,pjoin$load,pjoin$category, pjoin$scen, pjoin$cost, pjoin$pricing)
  bylist <- c("Group.1","Group.2","Group.3","Group.4","Group.5","Group.6","Group.7")
  pjoin <- merge(aggregate(pjoin[,c("net.final.price")], idlist, mean),
                 aggregate(pjoin[,c("net.final.price")], idlist, sd), by=c(bylist))
  idlist <- list(dfensumqmonth$theta,dfensumqmonth$ev,dfensumqmonth$load,dfensumqmonth$category, dfensumqmonth$scen, dfensumqmonth$cost, dfensumqmonth$pricing)
  pjoin <- merge(pjoin,
                 aggregate(dfensumqmonth[,c("sumoq")], idlist, mean), 
                 by=c(bylist))
  colnames(pjoin) <- c("theta","ev","load","category","scen","cost","pricing","mean", "sd","meanq")
  idlist <- list(dfen$theta,dfen$ev,dfen$load,dfen$category, dfen$scen, dfen$cost, dfen$pricing, dfen$month, dfen$day)
  pjoinevcs <- aggregate(dfen[,c("csmonth","evmonth")], idlist, mean)
  colnames(pjoinevcs) <- c("theta","ev","load","category","scen","cost","pricing","month","day","csmonth","evmonth")
  idlist <- list(pjoinevcs$theta,pjoinevcs$ev,pjoinevcs$load,pjoinevcs$category, pjoinevcs$scen, pjoinevcs$cost, pjoinevcs$pricing)
  pjoinevcs1 <- aggregate(pjoinevcs[,c("csmonth","evmonth")], idlist, sum)
  colnames(pjoinevcs1) <- c("theta","ev","load","category","scen","cost","pricing","csmonth","evmonth")
  bylist <- c("theta","ev","load","category","scen","cost","pricing")
  pjoin <- merge(pjoin, pjoinevcs1,by=c(bylist))
  
  idlist <- list(dfen$theta,dfen$ev,dfen$load,dfen$category, dfen$scen, dfen$cost, dfen$pricing)
  pjoinerenshare <- aggregate(dfen[,c("renshare","tots","co2")], idlist, mean)
  colnames(pjoinerenshare) <- c("theta","ev","load","category","scen","cost","pricing","Renewable_Share","tots","co2")
  
  ###################################################################
  # % expenditure
  joinall <- merge(pjoinerenshare[,c("theta","ev","load","category","scen","cost","pricing","Renewable_Share","tots","co2")],
                   pjoin[,c("theta","ev","load","category","scen","cost","pricing","mean","meanq","sd","csmonth","evmonth")],
                   by=c("theta","ev","load","category","scen","cost","pricing"), all.y=TRUE)
  #joinall <- merge(joinall,dfj[,c("theta","ev","load","category","scen","cost","pricing","tots")], by=c("theta","ev","load","category","scen","cost","pricing"))
  joinall <- transform(joinall, evmonth=-evmonth)
  joinall <- transform(joinall, csmonth=-csmonth)
  joinall <- transform(joinall, dps=tots-evmonth-csmonth)
  joinall <- transform(joinall, csmonth1=evmonth+csmonth)
  joinall <- subset(joinall, select=-c(csmonth))
  joinall <- transform(joinall, csmonth=csmonth1)
  joinall <- subset(joinall, select=-c(csmonth1))
  joinall <- merge(joinall, x1, by=c("theta","ev","load","category","scen","cost","pricing"))
  #combinetable <- joinall
  joinall <- joinall[,colnames(joinall)!='evexp']
  
  #baseline is fossil future flat
  basesum <- joinall[which(joinall$cost=="future" & joinall$category=="fossil" & joinall$pricing=="flat"),]
  basesum <- merge(basesum[,c("theta","ev","load","scen","csmonth","evmonth","tots","dps",
                              "cshighflex","csmidflex","csinflex",
                              "highexp","midexp","infexp","totexp","co2")], 
                   dbaseexp[which(dbaseexp$cost=="future" & dbaseexp$category=="fossil" & dbaseexp$pricing=="flat"),
                            c("theta","ev","load","scen","evexp1")], by=c("theta","ev","load","scen"), all.x=TRUE, all.y=FALSE)
  colnames(basesum) <- c("theta","ev","load","scen","csmonthb","evmonthb","totsb","dpsb",
                         "cshighflexb","csmidflexb","csinflexb",
                         "highexpb","midexpb","infexpb","totexpb","co2","evexp")
  
  basesum$highexpb <- (basesum$highexpb/basesum$totexpb)*(baselineexp+basesum$evexp)
  basesum$midexpb <- (basesum$midexpb/basesum$totexpb)*(baselineexp+basesum$evexp)
  basesum$infexpb <- (basesum$infexpb/basesum$totexpb)*(baselineexp+basesum$evexp)
  basesum$totexpb <- baselineexp+basesum$evexp
  #setwd('..') 
  write.csv(basesum, file = paste("data/basesum",thetaval,".csv", sep = ""), row.names = FALSE)
  
  joinall <- merge(joinall, basesum, by=c("theta","ev","load","scen"), all.x=TRUE, all.y=TRUE)
  joinall <- transform(joinall, dcs = 1*round((csmonth-csmonthb)/totexpb*100,1))
  joinall <- transform(joinall, evcost = -1*round((evmonth-evmonthb)/evexp*100,1))
  joinall <- transform(joinall, dps = 1*round((dps-dpsb)/totexpb*100,1))
  joinall <- transform(joinall, tots =1*round((tots-totsb)/totexpb*100,1))
  joinall <- transform(joinall, high =1*round((cshighflex-cshighflexb)/highexpb*100,1))
  joinall <- transform(joinall, mid =1*round((csmidflex-csmidflexb)/midexpb*100,1))
  joinall <- transform(joinall, inf =1*round((csinflex-csinflexb)/infexpb*100,1))
  
  #joinall <- transform(joinall, tot3=ifelse(high!=0,inf+high+mid,0))
  #joinall <- transform(joinall, inf =ifelse(high!=0,(inf/tot3) * dcs,0))
  #joinall <- transform(joinall, high =ifelse(high!=0,(high/tot3) * dcs,0))
  #joinall <- transform(joinall, mid =ifelse(high!=0,(mid/tot3) * dcs,0))
  joinall <- transform(joinall, mean = round(mean,0))
  joinall <- transform(joinall, meanq = round(meanq,0))
  joinall <- transform(joinall, sd = round(sd,0))
  write.csv(joinall, file = paste("table_",thetaval,".csv", sep = ""), row.names = FALSE)
  joinall <- joinall[,c("theta","ev","load","category","scen","cost","pricing","Renewable_Share",
                        "mean","meanq","sd","dcs","evcost","dps","tots","high","mid","inf")]
  maintable <- joinall
  #value of dynamic for total surplus
  dynbase <- joinall[which(joinall$pricing=="flat"),c("theta","ev","load","category","scen","cost","tots")]
  colnames(dynbase) <- c("theta","ev","load","category","scen","cost","totsb")
  joinall <- merge(joinall, dynbase, by=c("theta","ev","load","category","scen","cost"), all.x=TRUE, all.y=FALSE)
  joinall <- transform(joinall, dyntot = round((tots-totsb),2))
  joinall <- joinall[,c("theta","ev","load","category","cost","scen","pricing",
                        "Renewable_Share","mean","meanq","sd","dcs","evcost","dps","tots","high","mid","inf","dyntot")]
  
  #value of dynamic for consumer
  dyncons <- maintable[which(maintable$pricing=="flat"),c("theta","ev","load","category","scen","cost","dcs")]
  colnames(dyncons) <- c("theta","ev","load","category","scen","cost","dcsb")
  dyncons <- merge(maintable, dyncons, by=c("theta","ev","load","category","scen","cost"), all.x=TRUE, all.y=FALSE)
  dyncons <- transform(dyncons, dyncons = round((dcs-dcsb),2))
  dyncons <- dyncons[,c("theta","ev","load","category","cost","scen","pricing",
                        "Renewable_Share","mean","meanq","sd","dcs","evcost","dps","tots","high","mid","inf","dyncons")]
  
  #value of dynamic for producer
  dynprod <- maintable[which(maintable$pricing=="flat"),c("theta","ev","load","category","scen","cost","dps")]
  colnames(dynprod) <- c("theta","ev","load","category","scen","cost","dpsb")
  dynprod <- merge(maintable, dynprod, by=c("theta","ev","load","category","scen","cost"), all.x=TRUE, all.y=FALSE)
  dynprod <- transform(dynprod, dynprod = round((dps-dpsb),2))
  dynprod <- dynprod[,c("theta","ev","load","category","cost","scen","pricing",
                        "Renewable_Share","mean","meanq","sd","dcs","evcost","dps","tots","high","mid","inf","dynprod")]
  
  graphdata <- joinall
  
  #write.csv(graphdata, file = paste("tables/maintable_theta=",thetaval,".csv", sep = ""), row.names = FALSE)
  
  forlatex <- graphdata
  #naming the figures
  num1 <- ifelse(thetaval!=2,"S3","S5")
  num2 <- ifelse(thetaval!=2,"S4","S6")
  source(file = "Fig4-5.R")
  source(file = "Tab5.R")
  
  ##############################################################################
  ##APPENDIX FigS1b Share of hours with Marginal Cost < 5 or < 1 cents/KWh
  ##############################################################################
  pricebins <- transform(recent, mc50=ifelse(net.final.mc<50, 1,0))
  pricebins <- transform(pricebins, mc10=ifelse(net.final.mc<10, 1,0))
  idlist <- list(pricebins$category, pricebins$scen, pricebins$cost,pricebins$pricing,pricebins$load,pricebins$ev)
  pricebins <- aggregate(list(pricebins$mc50,pricebins$mc10), idlist, mean)
  colnames(pricebins) <- c("category","scen","cost","pricing","load","ev","sharemc50","sharemc10")
  pricebins <- transform(pricebins, category=ifelse(pricebins$category=="free","Unconstrained",paste(pricebins$category)))
  pricebins <- transform(pricebins, category=ifelse(pricebins$category=="rps_100","100% Clean",paste(pricebins$category)))
  pricebins <- transform(pricebins, category=ifelse(pricebins$category=="fossil","Fossil",paste(pricebins$category)))
  pricebins <- transform(pricebins, Demand_flexibility=ifelse(pricebins$scen==1,"1-Optimistic", ifelse(pricebins$scen==2,"2-Moderate","3-Pessimistic")))
  pricebins <- transform(pricebins, pricing=ifelse(pricebins$pricing=="dynamic","Real Time",paste(pricebins$pricing)))
  pricebins <- transform(pricebins, pricing=ifelse(pricebins$pricing=="flat","Flat",paste(pricebins$pricing)))
  pricebins <- transform(pricebins, pricing=paste(pricebins$pricing, "pricing"))
  pricebins <- transform(pricebins, cost=ifelse(pricebins$cost=="current","2016 Cost","2045 Cost"))
  pricebins <- transform(pricebins, categoryord=ifelse(pricebins$category=="100% Clean",2, ifelse(pricebins$category=="Fossil",1,3)))
  pricebins <- transform(pricebins, category = reorder(category, categoryord))
  
  ggplot(data = pricebins) + 
    geom_point(aes(category, sharemc10,color="olivedrab4"), size=6, position = "jitter", show.legend=F)+
    geom_point(aes(category, sharemc50,color="darkseagreen2"), size=4, position = "jitter")+
    facet_grid(cost~ pricing, scale = "free_y", space = "fixed") +
    scale_size_area() + 
    xlab("") +
    ylab("Share of hours with marginal cost") + theme_classic()+
    scale_color_manual("Price",values=c("olivedrab4","darkseagreen2"),labels=c("<1 cents/KWh","<5 cents/KWh"),
                       guide=guide_legend(override.aes = list(color=c("darkseagreen2","olivedrab4"))))+
    theme(legend.direction="horizontal",legend.position = "bottom", 
          legend.title = element_text(colour = "black", size = 15),
          legend.text = element_text(colour = "black", size = 15),
          axis.ticks=element_blank(), legend.key.width = unit(0.5,"cm"), 
          plot.title = element_text(hjust = -10),
          axis.text=element_text(size=13, angle = 0),
          axis.title=element_text(size=15),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15))+
    geom_hline(yintercept=0, col = "gray80")
  ggsave(paste("plot/FigS1b_MarginalCost_theta",thetaval,".pdf",sep=""),height=6.1,width=10)
  
}
