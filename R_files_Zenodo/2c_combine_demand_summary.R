#made with R version 4.0.5 (2021-03-31)
## modified on 22 Nov 2022
# need demand_response_summary_final
# produce demandsummaryrps.csv for FIG8

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
  categorysq <- c("rps")
  scenariosq <- seq(1,3)
  costsq <- c("current", "future")
  pricingsq <- c("flat","dynamic")
  scalem <- 1000000
  
  rpstargetsq <- seq(10,80,by=10)
  rpstargetsq<- c(rpstargetsq,seq(82,98,by=2))
  rpstargetsq <- paste("0",rpstargetsq,sep="") 
  rpstargetsq<- c("000",rpstargetsq)
  rpstargetsq<- c(rpstargetsq,"100")
  
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
                    if (cat=="rps"){ 
                      for (rpstarget in rpstargetsq){
                        filename <- paste(dirfold0,'demand_response_summary_final_',cat,'_',rpstarget,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
                        #print(filename)
                        dfw <- read.csv(filename, header=TRUE, sep=",")
                        dfw$theta=thetaval
                        dfw$category=cat
                        dfw$ev=ev
                        dfw$load=load
                        dfw$cost=cost
                        dfw$pricing=flatdyn
                        dfw$scen=scen
                        dfw$rpstarget=rpstarget
                        ifelse(exists("sumtab"), sumtab <- rbind(sumtab, dfw), sumtab <- dfw)
                        }
                      } else {
                        filename <- paste(dirfold0,'demand_response_summary_final_',cat,'_',cost,'_cost_',flatdyn,'_',load,'_load_',ev,'_ev_scen',scen,'.csv',sep="")
                        #print(filename)
                        dfw <- read.csv(filename, header=TRUE, sep=",")
                        dfw$theta=thetaval
                        dfw$category=cat
                        dfw$ev=ev
                        dfw$load=load
                        dfw$cost=cost
                        dfw$pricing=flatdyn
                        dfw$scen=scen
                        dfw$rpstarget=""
                        ifelse(exists("sumtab"), 
                               sumtab <- rbind(sumtab, dfw),
                               sumtab <- dfw)
                      }
                  }
                }, error=function(e){cat("ERROR :",conditionMessage(e),thetaval,load,ev,cat,cost,flatdyn,scen, "\\n")})
            }
          }
        }
      }
    }
  }
  #setwd('..') 
}  
write.csv(sumtab, file = paste("data/demandsummaryrps.csv"), row.names = FALSE)
