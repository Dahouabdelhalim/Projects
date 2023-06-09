wormboot<-function(reps, mydata){
  # Expects a number of reps for the bootstrap
  # and a vector of worm CFU data for individuals
  capp<-length(mydata)
  batch5<-rep(0,reps)
  batch10<-rep(0,reps)
  batch20<-rep(0,reps)
  batch50<-rep(0, reps)
  for(i in 1:reps){
    idx5<-sample(1:capp,5,replace=TRUE)
    idx10<-sample(1:capp,10,replace=TRUE)
    idx20<-sample(1:capp,20,replace=TRUE)
    idx50<-sample(1:capp,50,replace=TRUE)
    batch5[i]<-mean(mydata[idx5])
    batch10[i]<-mean(mydata[idx10])
    batch20[i]<-mean(mydata[idx20])
    batch50[i]<-mean(mydata[idx50])
    }
  batch5log<-log10(batch5+1)
  batch10log<-log10(batch10+1)
  batch20log<-log10(batch20+1) 
  batch50log<-log10(batch50+1) 
  batch<-c(rep(1,times=capp), rep(5, times=reps), rep(10, times=reps), rep(20, times=reps), rep(50, times=reps))
  logCFU<-log10(mydata*runif(capp,down, up))
  logCount<-c(logCFU, batch5log, batch10log, batch20log, batch50log)
  Count<-10^logCount
  dataSet<-data.frame(batch, Count, logCount)
  return(dataSet)
}