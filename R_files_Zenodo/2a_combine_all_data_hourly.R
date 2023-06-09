#made with R version 4.0.3 (2020-10-10)
## combine all results from HPC
## need files: demand_summary_final, dual_costs, energy_sources
########################################################################
#loads time point weights
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
for (dirfold in directory) {
  dirfold0 <- paste("./outputshpc/",dirfold,"/",sep="")
  #setwd(dirfold)
  print(getwd())
  thetaval <- as.numeric(substr(dirfold, 14,16))
  thetaval=ifelse(thetaval!=2,thetaval/10,thetaval)
  print(paste("working with theta:",thetaval))
  
  evsq <- c("half", "full","2016")
  loadsq <- c("2007", "2045")
  categorysq <- c("rps_100","fossil", "free")
  scenariosq <- seq(1,3)
  costsq <- c("current", "future")
  pricingsq <- c("flat","dynamic")
  scalem <- 1000000
  
  if (file.exists(paste("data/energysourcesCO2",thetaval,".csv", sep = "")))
  recent <- read.csv(paste("data/energysourcesCO2",thetaval,".csv", sep = ""), header=TRUE, sep=",")

  for (cat in categorysq) {
    for (ev in evsq) {
      for(flatdyn in pricingsq) {
        for (load in loadsq) {
          for (cost in costsq) {
            for (scen in scenariosq){
              if (exists("recent")){
                if (any(recent$category==cat&recent$ev==ev&recent$pricing==flatdyn&recent$load==load&recent$cost==cost&recent$scen==scen)){
                  print(paste("Already run",thetaval,load,ev,cat,cost,flatdyn,scen))
                } else {
                  if ((load=="2007" & ev!="half")){
                    print(paste("Skip",thetaval,load,ev,cat,cost,flatdyn,scen))
                  } else {
                    tryCatch({ 
                      print(paste(thetaval,load,ev,cat,cost,flatdyn,scen))
                      cat <- ifelse(cat=="rps","rps_100",cat) 
                      
                      #gather summary. DONT FORGET THERE IS A SET OF SAME CODE BELOW SO YOU HAVE TO CHANGE IT TOO
                      filename <- paste(dirfold0,'demand_response_summary_final_',cat,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
                      dfw <- read.csv(filename, header=TRUE, sep=",")
                      co2 <- as.numeric(dfw$co2_per_year_2045)
                      
                      # ifelse(exists("sumtab"), 
                      #        sumtab <- rbind(sumtab, dfw),
                      #        sumtab <- dfw)
                      # dfw$theta=thetaval
                      # dfw$category=cat
                      # dfw$ev=ev
                      # dfw$load=load
                      # dfw$cost=cost
                      # dfw$pricing=flatdyn
                      # dfw$scen=scen
                      
                      #gather cs and ev charge cost by sample days
                      filename_d <- paste(dirfold0,'dual_costs_',cat,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
                      dfdual <- data.table(read.csv(filename_d, header=FALSE, sep=",",skip=1)[,1:6])
                      colnames(dfdual) <- c("constraint","constraint_second_key","direction","bound","dual","total_cost")
                      csmonth <- data.frame(dfdual[constraint %like% "DR_Convex_Bid_Weight"])[,c("constraint_second_key","total_cost")]
                      evmonth <- data.frame(dfdual[constraint %like% "EV"])[,c("constraint_second_key","total_cost")]
                      evmonth <- transform(evmonth,total_cost=abs(total_cost))
                      evmonth <- transform(evmonth, timeseries=as.numeric(substr(evmonth$constraint_second_key,1,8)))
                      csmonth <- transform(csmonth, timeseries=as.numeric(gsub("]", "", csmonth$constraint_second_key)))
                      names(csmonth)[names(csmonth) == 'total_cost'] <- 'csmonth'
                      names(evmonth)[names(evmonth) == 'total_cost'] <- 'evmonth'
                      idlist <- list(evmonth$timeseries)
                      evmonth <- aggregate(evmonth[,c("evmonth")], idlist, sum)
                      colnames(evmonth) <- c("timeseries","evmonth")
                      #write.csv(evmonth, file = paste("evmonth",thetaval,".csv", sep = ""), row.names = FALSE)
                      
                      csev <- merge(w,csmonth, by=c("timeseries"), all.x=TRUE)
                      csev <- csev[,c("timeseries","ts_scale_to_period","month","day","hour","samplew","csmonth")]
                      csev[is.na(csev)] <- 0
                      csev <- merge(csev, evmonth, by=c("timeseries"), all.x=TRUE)
                      csev <- csev[,c("timeseries","ts_scale_to_period","month","day","hour","samplew","csmonth","evmonth")]
                      csev[is.na(csev)] <- 0
                      
                      csev1 <- csev
                      csev1$theta=thetaval
                      csev1$category=cat
                      csev1$ev=ev
                      csev1$load=load
                      csev1$cost=cost
                      csev1$pricing=flatdyn
                      csev1$scen=scen
                      ifelse(exists("csev0"), 
                             csev0 <- rbind(csev0, csev1),
                             csev0 <- csev1)
                      #write.csv(csev0, file = paste("data/csev_",thetaval,".csv", sep = ""), row.names = FALSE)
                      
                      
                      #gather final iteration in housely energy sources file
                      filename_e <- paste(dirfold0,'energy_sources_',cat,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
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
                      dfenb <- dfen[,c("month","day","hour","net.final.mc","net.final.price","offered.price.energy","offered.price.energy.down", "offered.price.energy.up","final.price.energy","net.final.q","base_load","base_price","evexp")]
                      #define id as sampledays identification
                      dfenb <- setDT(dfenb)[, id := .GRP, by = c("month","day")]
                      
                      #Compute each sample day expenditure and CS customers
                      for (i in unique(dfenb$id)) {
                        datause0 <- dfenb[which(dfenb$id==i),]
                        month <- unique(datause0$month)
                        day <- unique(datause0$day)
                        p1 <- list(datause0$offered.price.energy, datause0$offered.price.energy.up, datause0$offered.price.energy.down)
                        p0 <- list(datause0$base_price, datause0$base_price*0, datause0$base_price*0)
                        q0 <- list(datause0$base_load, datause0$base_load*0, datause0$base_load*0)
                        #integralcs will give output: (list(wtp, cs, q.dot.p(q1, p1), q1))
                        #pqxxx indicate expenditure, from output list[[3]]
                        #wtpxxx indicate cs, from output list[[2]]
                        #the value is total by sample days, but since the table is in hourly, therefore we repeat it 24 times
                        ifelse(i==1, csen <- cbind(qtot=datause0$net.final.q,month=month, day=day, hour=seq(0,23),
                                                   pqhighflex=rep((integralcs(p1, p0, q0, 1, month, scen)[[3]]),24),
                                                   pqmidflex=rep((integralcs(p1, p0, q0, 2, month, scen)[[3]]),24),
                                                   pqinflex=rep((integralcs(p1, p0, q0, 3, month, scen)[[3]]),24),
                                                   wtphighflex=rep(integralcs(p1, p0, q0, 1, month, scen)[[2]],24),
                                                   wtpmidflex=rep(integralcs(p1, p0, q0, 2, month, scen)[[2]],24),
                                                   wtpinflex=rep(integralcs(p1, p0, q0, 3, month, scen)[[2]],24),
                                                   wtptot=rep(integralcs(p1, p0, q0, 0, month, scen)[[2]],24)),
                               csen <- rbind(csen, cbind(qtot=datause0$net.final.q,month=month, day=day, hour=seq(0,23),
                                                         pqhighflex=rep((integralcs(p1, p0, q0, 1, month, scen)[[3]]),24),
                                                         pqmidflex=rep((integralcs(p1, p0, q0, 2, month, scen)[[3]]),24),
                                                         pqinflex=rep((integralcs(p1, p0, q0, 3, month, scen)[[3]]),24),
                                                         wtphighflex=rep(integralcs(p1, p0, q0, 1, month, scen)[[2]],24),
                                                         wtpmidflex=rep(integralcs(p1, p0, q0, 2, month, scen)[[2]],24),
                                                         wtpinflex=rep(integralcs(p1, p0, q0, 3, month, scen)[[2]],24),
                                                         wtptot=rep(integralcs(p1, p0, q0, 0, month, scen)[[2]],24))))
                      }
                      
                      #print("monthly consumer surplus calculated")
                      dfenb <- merge(dfenb, csen, by=c("month","day","hour"))
                      dfenb <- transform(dfenb, theta=thetaval)
                      dfenb <- transform(dfenb, category=cat)
                      dfenb <- transform(dfenb, ev=ev)
                      dfenb <- transform(dfenb, load=load)
                      dfenb <- transform(dfenb, cost=cost)
                      dfenb <- transform(dfenb, pricing=flatdyn)
                      dfenb <- transform(dfenb, scen=scen)
                      dfenb <- transform(dfenb, renshare=round(dfw$renewable_share_2045*100,1))
                      dfenb <- transform(dfenb, co2=co2)
                      dfenb <- merge(dfenb,csev, by=c("month","day","hour"), all.x=TRUE)
                      dfenb <- transform(dfenb, tots=dfw$total_cost)
                      
                      ifelse(exists("recent"), 
                             recent <- rbind(recent, dfenb),
                             recent <- dfenb)
                      
                      #making plots
                      #sourceR(file = "production_plots.R")
                      #dfenb2 <- transform(dfenb2, days=days_in_month(as.Date(paste('2045-',dfenb2$month,'-','1',sep=""),"%Y-%m-%d")))
                      #print("hourly price and quantity table loaded")
                      write.csv(recent, file = paste("data/energysourcesCO2",thetaval,".csv", sep = ""), row.names = FALSE)
                      
                    }, error=function(e){cat("ERROR :",conditionMessage(e),thetaval,load,ev,cat,cost,flatdyn,scen, "\\n")})
                  }
                }
              } else {
                #similar as above
              if ((load=="2007" & ev!="half")){
                print(paste("Skip",thetaval,load,ev,cat,cost,flatdyn,scen))
                } else {
                tryCatch({ 
                  print(paste(thetaval,load,ev,cat,cost,flatdyn,scen))
		              cat <- ifelse(cat=="rps","rps_100",cat) 
		              
		              #gather summary
		              filename <- paste(dirfold0,'demand_response_summary_final_',cat,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
		              dfw <- read.csv(filename, header=TRUE, sep=",")
		              co2 <- as.numeric(dfw$co2_per_year_2045)
		              
		              # ifelse(exists("sumtab"), 
		              #        sumtab <- rbind(sumtab, dfw),
		              #        sumtab <- dfw)
		              # dfw$theta=thetaval
		              # dfw$category=cat
		              # dfw$ev=ev
		              # dfw$load=load
		              # dfw$cost=cost
		              # dfw$pricing=flatdyn
		              # dfw$scen=scen

		              #gather cs and ev charge cost by sample days
		              filename_d <- paste(dirfold0,'dual_costs_',cat,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
		              dfdual <- data.table(read.csv(filename_d, header=FALSE, sep=",",skip=1)[,1:6])
		              colnames(dfdual) <- c("constraint","constraint_second_key","direction","bound","dual","total_cost")
		              csmonth <- data.frame(dfdual[constraint %like% "DR_Convex_Bid_Weight"])[,c("constraint_second_key","total_cost")]
		              evmonth <- data.frame(dfdual[constraint %like% "EV"])[,c("constraint_second_key","total_cost")]
		              evmonth <- transform(evmonth,total_cost=abs(total_cost))
		              evmonth <- transform(evmonth, timeseries=as.numeric(substr(evmonth$constraint_second_key,1,8)))
		              csmonth <- transform(csmonth, timeseries=as.numeric(gsub("]", "", csmonth$constraint_second_key)))
		              names(csmonth)[names(csmonth) == 'total_cost'] <- 'csmonth'
		              names(evmonth)[names(evmonth) == 'total_cost'] <- 'evmonth'
		              idlist <- list(evmonth$timeseries)
		              evmonth <- aggregate(evmonth[,c("evmonth")], idlist, sum)
		              colnames(evmonth) <- c("timeseries","evmonth")
		              #write.csv(evmonth, file = paste("evmonth",thetaval,".csv", sep = ""), row.names = FALSE)
		              
		              csev <- merge(w,csmonth, by=c("timeseries"), all.x=TRUE)
		              csev <- csev[,c("timeseries","ts_scale_to_period","month","day","hour","samplew","csmonth")]
		              csev[is.na(csev)] <- 0
		              csev <- merge(csev, evmonth, by=c("timeseries"), all.x=TRUE)
		              csev <- csev[,c("timeseries","ts_scale_to_period","month","day","hour","samplew","csmonth","evmonth")]
		              csev[is.na(csev)] <- 0
		              
		              csev1 <- csev
		              csev1$theta=thetaval
		              csev1$category=cat
		              csev1$ev=ev
		              csev1$load=load
		              csev1$cost=cost
		              csev1$pricing=flatdyn
		              csev1$scen=scen
		              ifelse(exists("csev0"), 
		                     csev0 <- rbind(csev0, csev1),
		                     csev0 <- csev1)
		              #write.csv(csev0, file = paste("data/csev_",thetaval,".csv", sep = ""), row.names = FALSE)
		              
		              #gather final iteration in housely energy sources file
		              filename_e <- paste(dirfold0,'energy_sources_',cat,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
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
		              dfenb <- dfen[,c("month","day","hour","net.final.mc","net.final.price","offered.price.energy","offered.price.energy.down", "offered.price.energy.up","final.price.energy","net.final.q","base_load","base_price","evexp")]
		              #define id as sampledays identification
		              dfenb <- setDT(dfenb)[, id := .GRP, by = c("month","day")]
		              
		              #Compute each sample day expenditure and CS customers
		              for (i in unique(dfenb$id)) {
		                datause0 <- dfenb[which(dfenb$id==i),]
		                month <- unique(datause0$month)
		                day <- unique(datause0$day)
		                p1 <- list(datause0$offered.price.energy, datause0$offered.price.energy.up, datause0$offered.price.energy.down)
		                p0 <- list(datause0$base_price, datause0$base_price*0, datause0$base_price*0)
		                q0 <- list(datause0$base_load, datause0$base_load*0, datause0$base_load*0)
		                #integralcs function is loaded from flexshares.R
		                # integralcs <- function(p1, p0, q0, sigmath, month, scenario) {
		                #   # calculate the line integral of q dot p from p0 to p1, with the specified demand parameters
		                #   iqp <- integrate(
		                #     Vectorize(q.dot.p.for.t, vectorize.args=c('t')), 0, 1,
		                #     p0, p1, q0[[1]], p0[[1]], sigmath, month, scenario)$value
		                #   # get the quantities bid for energy and up and down reserves (assuming p0 and q0 are baselines)
		                #   q1 <- demandbytype(q0[[1]], p0[[1]], p1[[1]], p1[[2]], p1[[3]], sigmath, month, scenario)
		                #   # calculate wtp for this bid
		                #   wtp <- q.dot.p(q1, p1) - q.dot.p(q0, p0) - iqp
		                #   cs <- wtp - q.dot.p(q1, p1) #CS contains q.dot.p(q0, p0)
		                #   return(list(wtp, cs, q.dot.p(q1, p1), q1))
		                # }
		                #we will exctract 2 output: expenditure and CS customers
		                #pqxxx indicate expenditure, from output list[[3]]
		                #wtpxxx indicate cs, from output list[[2]]
		                #the value is total by sample days, but since the table is in hourly, therefore we repeat it 24 times
		                ifelse(i==1, csen <- cbind(qtot=datause0$net.final.q,month=month, day=day, hour=seq(0,23),
		                                           pqhighflex=rep((integralcs(p1, p0, q0, 1, month, scen)[[3]]),24),
		                                           pqmidflex=rep((integralcs(p1, p0, q0, 2, month, scen)[[3]]),24),
		                                           pqinflex=rep((integralcs(p1, p0, q0, 3, month, scen)[[3]]),24),
		                                           wtphighflex=rep(integralcs(p1, p0, q0, 1, month, scen)[[2]],24),
		                                           wtpmidflex=rep(integralcs(p1, p0, q0, 2, month, scen)[[2]],24),
		                                           wtpinflex=rep(integralcs(p1, p0, q0, 3, month, scen)[[2]],24),
		                                           wtptot=rep(integralcs(p1, p0, q0, 0, month, scen)[[2]],24)),
		                       csen <- rbind(csen, cbind(qtot=datause0$net.final.q,month=month, day=day, hour=seq(0,23),
		                                                 pqhighflex=rep((integralcs(p1, p0, q0, 1, month, scen)[[3]]),24),
		                                                 pqmidflex=rep((integralcs(p1, p0, q0, 2, month, scen)[[3]]),24),
		                                                 pqinflex=rep((integralcs(p1, p0, q0, 3, month, scen)[[3]]),24),
		                                                 wtphighflex=rep(integralcs(p1, p0, q0, 1, month, scen)[[2]],24),
		                                                 wtpmidflex=rep(integralcs(p1, p0, q0, 2, month, scen)[[2]],24),
		                                                 wtpinflex=rep(integralcs(p1, p0, q0, 3, month, scen)[[2]],24),
		                                                 wtptot=rep(integralcs(p1, p0, q0, 0, month, scen)[[2]],24))))
		              }
		              
		              #print("monthly consumer surplus calculated")
		              dfenb <- merge(dfenb, csen, by=c("month","day","hour"))
		              dfenb <- transform(dfenb, theta=thetaval)
		              dfenb <- transform(dfenb, category=cat)
		              dfenb <- transform(dfenb, ev=ev)
		              dfenb <- transform(dfenb, load=load)
		              dfenb <- transform(dfenb, cost=cost)
		              dfenb <- transform(dfenb, pricing=flatdyn)
		              dfenb <- transform(dfenb, scen=scen)
		              dfenb <- transform(dfenb, renshare=round(dfw$renewable_share_2045*100,1))
		              dfenb <- transform(dfenb, co2=co2)
		              dfenb <- merge(dfenb,csev, by=c("month","day","hour"), all.x=TRUE)
		              dfenb <- transform(dfenb, tots=dfw$total_cost)
		              
                  ifelse(exists("recent"), 
                         recent <- rbind(recent, dfenb),
                         recent <- dfenb)

                  #making plots
                  #sourceR(file = "production_plots.R")
                  #dfenb2 <- transform(dfenb2, days=days_in_month(as.Date(paste('2045-',dfenb2$month,'-','1',sep=""),"%Y-%m-%d")))
                  #print("hourly price and quantity table loaded")
                  write.csv(recent, file = paste("data/energysourcesCO2",thetaval,".csv", sep = ""), row.names = FALSE)
                  
                }, error=function(e){cat("ERROR :",conditionMessage(e),thetaval,load,ev,cat,cost,flatdyn,scen, "\\n")})
                }
              }
            }
          }
        }
      }
    }
  }
  #setwd('..') 
  rm(recent)
}  