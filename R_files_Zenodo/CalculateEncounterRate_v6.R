#Calculates the encounter rate in terms of porpoise positive seconds per second of deployment
# time, on a position by deployment by month basis

#Author: Len Thomas
# Update: 27/6/2013 - updated input files to use files produced on 14-05-2014
# Update: 15/7/2014 - increased robustness of read.effort.file so it can cope with MinutesON=0
# Update: 17/7/2014 - now calculates statistics by minute within day
# Update: 22/08/2014 - updated to use new dataset
# Update: 15/10/2014 - updated to use new dataset produced on 13-10-2014

effort.file.name<-"detections and environment - validated and cropped - 20141013.txt"
click.times.file.name<-"click details - validated and croppped - 20141013.txt"
effort.bymonth.file.name<-"effort.bymonth.bymin.txt"
click.seconds.bymonth.file.name<-"click.seconds.bymonth.bymin.txt"
n.bymonth.file.name<-"n.bymonth.bymin.txt"

#Number of minutes in a day
n.minutes<-1440

read.effort.file<-TRUE
if(read.effort.file){
  #Read in effort file, and work out number of seconds on by position, deployment, month
  #Write out record for each month
  con<-file(effort.file.name,"r")
  #read first line
  oneLine<-readLines(con,n=1)
  
  #Initalize
  old.deployment<-old.year<-old.month<-""
  onsecs<-numeric(n.minutes)
  cat("deployment\\tyear\\tmonth\\tminute\\teffort.secs\\n",file=effort.bymonth.file.name,append=F)
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    tok<-strsplit(oneLine,"\\t")
    #Only parse lines that contain the deployment number, separated by a dash
    if(substr(tok[[1]][1],5,5)=="-"){
  
      #collect new data
      new.deployment<-substr(tok[[1]][1],1,6)
      new.date.time<-strptime(tok[[1]][2],format="%d/%m/%Y %H:%M")
      new.year<-strftime(new.date.time,format="%Y")
      new.month<-strftime(new.date.time,format="%m")
      #if different from the old data, write out info for the old month
      if(!((old.deployment==new.deployment)&(old.year==new.year)&(old.month==new.month))){
       if(!(old.deployment=="")){
         for(i in 1:n.minutes){
           cat(paste(old.deployment,old.year,old.month,i,onsecs[i],sep="\\t"),"\\n",file=effort.bymonth.file.name,append=T)
         }
       }
       old.deployment<-new.deployment;old.year<-new.year;old.month<-new.month
       onsecs[1:n.minutes]<-0
      }
  
      #Add this effort, in seconds
      minute<-as.numeric(strftime(new.date.time,format="%H"))*60+as.numeric(strftime(new.date.time,format="%M"))+1
      #tok[[6]] is MinutesON; tok[[10]] is %TimeLost
      MinutesOn<-as.numeric(tok[[1]][6])
      #I do it this way because when MinutesOn is 0 the line has a different
      # format, and as.numeric(tok[[1]][10]) is NA - e.g., line 935321.
      if(MinutesOn>0){
        onsecs[minute]<-onsecs[minute]+MinutesOn*60*(1-as.numeric(tok[[1]][10])/100)
      }
    }
  }
  close(con)
  #Write out the last line
  if(!(old.deployment=="")){
    for(i in 1:n.minutes){
      cat(paste(old.deployment,old.year,old.month,i,onsecs[i],sep="\\t"),"\\n",file=effort.bymonth.file.name,append=T)
    }
  }
}

read.click.file<-TRUE
if(read.click.file){
  #Read in click times file, and work out number of porpoise positive seconds by position, deployment, month
  #Write out record for each month
  #Read in effort file, and work out number of seconds on by position, deployment, month
  #Write out record for each month
  con<-file(click.times.file.name,"r")
  #read first line
  oneLine<-readLines(con,n=1)
  #Initalize
  old.deployment<-old.year<-old.month<-old.date.time<-""
  old.second<--1
  clicksecs<-numeric(n.minutes)
  cat("deployment\\tyear\\tmonth\\tminute\\tclick.secs\\n",file=click.seconds.bymonth.file.name,append=F)
  
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    tok<-strsplit(oneLine,"\\t")
    #Only parse lines that contain the deployment number, separated by a dash
    if(substr(tok[[1]][1],5,5)=="-"){
      
      #collect new data
      new.deployment<-tok[[1]][1]
      new.date.time<-strptime(tok[[1]][2],format="%d/%m/%Y %H:%M")
      new.year<-strftime(new.date.time,format="%Y")
      new.month<-strftime(new.date.time,format="%m")
      new.second<-as.numeric(tok[[1]][4])%/%1E6
      #see if it's a new month or deployment
      is.new.record<-!((old.deployment==new.deployment)&(old.year==new.year)&(old.month==new.month))
      #it's not a new record if it's the first record in the dataset
      if(old.deployment=="") is.new.record<-F
      
      #if it's not a new record, but it is a new second, then add 1 to the porpoise
      # positive seconds
      if(!is.new.record){
        if(!((as.character(new.date.time)==old.date.time)&(new.second==old.second))){
          minute<-as.numeric(strftime(new.date.time,format="%H"))*60+as.numeric(strftime(new.date.time,format="%M"))+1
          clicksecs[minute]<-clicksecs[minute]+1
        }
      }
        
      #if it's a new record, write out the old one
      if(is.new.record){
        for(i in 1:n.minutes){
          cat(paste(old.deployment,old.year,old.month,i,clicksecs[i],sep="\\t"),"\\n",file=click.seconds.bymonth.file.name,append=T)
        }
      }
      #update the old records
      if(is.new.record|(old.deployment=="")){
        old.deployment<-new.deployment;old.year<-new.year;old.month<-new.month
        old.date.time<-new.date.time;old.second<-new.second
        minute<-as.numeric(strftime(new.date.time,format="%H"))*60+as.numeric(strftime(new.date.time,format="%M"))+1
        clicksecs[1:n.minutes]<-0
        clicksecs[minute]<-1
      }
    
    }
  }
  close(con)
  #Write out the last line
  if(!(old.deployment=="")){
    for(i in 1:n.minutes){
      cat(paste(old.deployment,old.year,old.month,i,clicksecs[i],sep="\\t"),"\\n",file=click.seconds.bymonth.file.name,append=T)
    }
  }
}

merge.files<-TRUE
run.checks<-FALSE #Note - this can take a *very* long time if you say TRUE
if(merge.files){
    effort.bymonth<-read.table(effort.bymonth.file.name,header=TRUE,sep="\\t",colClasses=c("character",rep("numeric",4)))
    if(run.checks){
      #Check you only have one unique deployment, year, month, minute
      tmp<-duplicated(effort.bymonth[,1:4])
      if(sum(tmp)>0) {
        #Will get duplicates if the effort file was not sorted by deployment, year and month
        cat("Warning - duplicated records in effort.bymonth\\n")
        tmp.dup<-effort.bymonth[tmp,]
        tmp.unique<-effort.bymonth[!tmp,]
        #Go through each duplicate, working out which unique value it corresponds to, and combine
        for(i in 1:dim(tmp.dup)[1]){
          if(i%%1000==0) cat(i,date(),"\\n")
          #work out which entry it is in the uniques
          ind<-which(tmp.unique[,1]==tmp.dup[i,1]&
            tmp.unique[,2]==tmp.dup[i,2]&
            tmp.unique[,3]==tmp.dup[i,3]&
            tmp.unique[,4]==tmp.dup[i,4])
          if(length(ind)!=1) stop("Didn't find a match\\n")
          #add in the effort.secs from the duplicate
          tmp.unique[ind,5]<-tmp.unique[ind,5]+tmp.dup[i,5]
        }
        #Now they should be unique, and the total effort.secs the same, so can replace effort.bymonth
        #Check efort.secs the same
        if(!(sum(tmp.unique[,5])==sum(effort.bymonth[,5]))) stop("Something went wrong when getting rid of duplicates!\\n")
        effort.bymonth<-tmp.unique
        rm(tmp.unique,tmp.dup,tmp)
      }
    }
    click.seconds.bymonth<-read.table(click.seconds.bymonth.file.name,header=T,sep="\\t",colClasses=c("character",rep("numeric",4)))
    if(run.checks){
      tmp<-duplicated(click.seconds.bymonth[,1:4])
      if(sum(tmp)>0){
        stop("Duplicate records in click.seconds.bymonth\\n")
      }
    }
    n.bymonth<-merge(effort.bymonth,click.seconds.bymonth,sort=F,all=T)
    #NA for click secs means none heard in that period
    n.bymonth$click.secs[is.na(n.bymonth$click.secs)]<-0
    #Order them - sort leaves them in a strange order
    n.bymonth<-n.bymonth[order(n.bymonth[,1],n.bymonth[,2],n.bymonth[,3],n.bymonth[,4]),]
    #Check for duplicates
    if(dim(n.bymonth)[1]!=dim(unique(n.bymonth[,1:4]))[1]) stop("Merged file rows not unique\\n")
    #Write out
    write.table(n.bymonth,file=n.bymonth.file.name,row.names=F,col.names=T,sep="\\t",quote=F)
}


