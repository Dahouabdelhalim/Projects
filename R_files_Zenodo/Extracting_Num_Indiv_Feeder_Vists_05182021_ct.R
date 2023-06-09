# Extracting all feeder visits, include date and hour for each visit
# HAD SOME MISTAKEs EARLIER!!!!
  # "under_10s" portion of loop was written wrong! And no ensuring next row was same date!
  # under_10s is now under 4s, which is how I did the captive extractions in 2014
  # Also, now get rid of false reads before calculating bouts, as sometime little shifts occur,
  # which aren't really the bird leaving or getting chased away, but used to make separate bouts...
  # also, fixed issue that 79111 wasn't in original database
  # Bird 2551-79108 wasn't in the database for some fucking reason until 111414, now it is...
  # Omitting Nature center feeder from all analysis

rm(list=ls())
setwd("C:/Users/court/Desktop/HOFI Contact Rates/Aviary RFID 2017")

# this gives a list of files whose size is > 10 bytes, which excludes all blanks, but no others
data.files<-list.files(pattern = "DATA", full.names=TRUE, ignore.case = F)
#data.files<-data.files[-grep("./NA", data.files)] ##something funky going on with this line
length(data.files)
data.files.info<-file.info(data.files)
data.files.no.blanks<-dimnames(data.files.info[data.files.info$size>10,])[[1]]
length(data.files.no.blanks)

data.files<-list.files(pattern = "DATA", full.names=TRUE, ignore.case = F)
length(data.files)
data.files.info<-file.info(data.files)
data.files.no.blanks<-dimnames(data.files.info[data.files.info$size>10,])[[1]]
length(data.files.no.blanks)

# this is the database of deployed PIT tags, NOTE:  now up through 020214!!!
# IF YOU REOPEN THIS FILE IN EXCEL, IT WILL MESS UP THE PIT TAG NUMBERS
# VERY IMPORTANT TO REALIZE THIS
# changed all PIT tags to have a 0 in front of them in excel
# also, fixed issue that 79111 wasn't in original database
##PIT_tags<-read.csv("C:/Users/court/Desktop/HOFI Contact Rates/Deployed PIT Tags_100717.csv")
PIT_tags<-read.csv("C:/Users/court/Desktop/HOFI Contact Rates/Aviary PIT List 2017--NEW PIT TAGS.csv")

unique(PIT_tags$PIT.Tag..)
uniquePITs<-as.character(na.omit(unique(PIT_tags$PIT.Tag..)))


# also, increase the max memory size to deal with this mess
memory.limit(size=2000)
memory.limit()
# a useful function for the loop below
position<-function(x,y){
    x[y]
 }


# Need some rule about what constitutes a continuous visit
# To start: No gaps of > 4s or other bird in between

# start with an example file


filex<-read.table(data.files.no.blanks[17],sep=" ",skip=1,header=F,col.names=c("PIT_Tag","Date","Time"))
    Reader<-rep(substr(data.files.no.blanks[17], start=3, stop=6),nrow(filex))
    Site<-rep(substr(data.files.no.blanks[17], start=3, stop=4),nrow(filex))
    Port<-rep(substr(data.files.no.blanks[17], start=5, stop=6),nrow(filex))
    hr<-as.numeric(lapply(strsplit(as.character(filex$Time),":"),position,1))
#says NAs induced fixed
#which(is.na(hr)==T)
filex[1:100,]
    mindec<-as.numeric(lapply(strsplit(as.character(filex$Time),":"),position,2))/60
    secdec<-as.numeric(lapply(strsplit(as.character(filex$Time),":"),position,3))/3600
    Dectime<-hr+mindec+secdec
  filex<-cbind(filex,Dectime,Reader,Site,Port)
  filex<-filex[as.character(filex$PIT_Tag) %in% uniquePITs==T | as.character(filex$PIT_Tag) %in% substr(uniquePITs, start=2, stop=10)==T,]


    samebird_as_above<-c(NA,as.character(filex$PIT_Tag)[2:length(filex$PIT_Tag)]==as.character(filex$PIT_Tag)[1:(length(filex$PIT_Tag)-1)])
    samedate_as_above<-c(NA,as.character(filex$Date)[2:length(filex$Date)]==as.character(filex$Date)[1:(length(filex$Date)-1)])
    under4s_from_above<-c(NA,filex$Dectime[1:(length(filex$Dectime)-1)]>(filex$Dectime[2:length(filex$Dectime)]-(4.00001/3600))) # changed from ealiers
    #under2s_from_above<-c(NA,filex$Dectime[1:(length(filex$Dectime)-1)]>(filex$Dectime[2:length(filex$Dectime)]-(2.00001/3600))) 
    same_bout_as_above<-samebird_as_above==TRUE & samedate_as_above==TRUE & under4s_from_above==TRUE

    bout_starts<-c(1,which(same_bout_as_above %in% FALSE))   # yields index locations of when bouts start, including first line always
    bout_ends<-c(which(same_bout_as_above %in% FALSE)-1,nrow(filex))     # yields index locations of when bouts end, including last line always

    Bout_duration_sec<-3600*(filex$Dectime[bout_ends]-filex$Dectime[bout_starts])   # yields duration of each bout, in seconds
    Bout_duration_sec[Bout_duration_sec==0]<-0.5                    # makes all bouts calculated at 0 = 0.5, = maximum down time of reader

    bouts_filex<-data.frame(PIT_Tag=as.character(filex[bout_starts,1]),Date=filex[bout_starts,2],Start_Time=filex[bout_starts,3],
      End_Time=filex[bout_ends,3],Reader=filex[bout_starts,5],Site=filex[bout_starts,6],Port=filex[bout_starts,7],Bout_duration_s=Bout_duration_sec)

    bouts_filex
    
    

# make that a loop
a<-Sys.time()

bouts<-data.frame(PIT_Tag=NULL,Date=NULL,Start_Time=NULL,
      End_Time=NULL,Reader=NULL,Site=NULL,Port=NULL,Bout_duration_s=NULL)

for (i in 1:length(data.files.no.blanks)) {

filex<-read.table(data.files.no.blanks[i],sep=" ",skip=1,header=F,col.names=c("PIT_Tag","Date","Time"))
    Reader<-rep(substr(data.files.no.blanks[i], start=3, stop=6),nrow(filex))
    Site<-rep(substr(data.files.no.blanks[i], start=3, stop=4),nrow(filex))
    Port<-rep(substr(data.files.no.blanks[i], start=5, stop=6),nrow(filex))
    hr<-as.numeric(lapply(strsplit(as.character(filex$Time),":"),position,1))
    mindec<-as.numeric(lapply(strsplit(as.character(filex$Time),":"),position,2))/60
    secdec<-as.numeric(lapply(strsplit(as.character(filex$Time),":"),position,3))/3600
    Dectime<-hr+mindec+secdec
  filex<-cbind(filex,Dectime,Reader,Site,Port)
  filex<-filex[as.character(filex$PIT_Tag) %in% uniquePITs==T | as.character(filex$PIT_Tag) %in% substr(uniquePITs, start=2, stop=10)==T,]
ifelse(nrow(filex)<2,print(paste("file #",i,data.files.no.blanks[i],"NO VISITS")) & next,print(paste("file #",i,data.files.no.blanks[i],"OK")))

    samebird_as_above<-c(NA,as.character(filex$PIT_Tag)[2:length(filex$PIT_Tag)]==as.character(filex$PIT_Tag)[1:(length(filex$PIT_Tag)-1)])
    samedate_as_above<-c(NA,as.character(filex$Date)[2:length(filex$Date)]==as.character(filex$Date)[1:(length(filex$Date)-1)])
    under4s_from_above<-c(NA,filex$Dectime[1:(length(filex$Dectime)-1)]>(filex$Dectime[2:length(filex$Dectime)]-(4.00001/3600)))  # changed from earliers
    #under2s_from_above<-c(NA,filex$Dectime[1:(length(filex$Dectime)-1)]>(filex$Dectime[2:length(filex$Dectime)]-(2.00001/3600)))  # changed from earliers
    same_bout_as_above<-samebird_as_above==TRUE & samedate_as_above==TRUE & under4s_from_above==TRUE

    bout_starts<-c(1,which(same_bout_as_above %in% FALSE))   # yields index locations of when bouts start, including first line always
    bout_ends<-c(which(same_bout_as_above %in% FALSE)-1,nrow(filex))     # yields index locations of when bouts end, including last line always

    Bout_duration_sec<-3600*(filex$Dectime[bout_ends]-filex$Dectime[bout_starts])   # yields duration of each bout, in seconds
    Bout_duration_sec[Bout_duration_sec==0]<-0.5                    # makes all bouts calculated at 0 = 0.5, = maximum down time of reader

    bouts_filex<-data.frame(PIT_Tag=as.character(filex[bout_starts,1]),Date=filex[bout_starts,2],Start_Time=filex[bout_starts,3],
      End_Time=filex[bout_ends,3],Reader=filex[bout_starts,5],Site=filex[bout_starts,6],Port=filex[bout_starts,7],Bout_duration_s=Bout_duration_sec)

    bouts<-rbind(bouts,bouts_filex)
#    warns<-warnings()
#    print(warns)
}

b<-Sys.time()
b-a
# 56.64 s

bouts[6630,]
bouts[is.na(bouts)==T,]

nrow(bouts[bouts$PIT_Tag=="0418253BE2",])
nrow(bouts[bouts$Site=="NA",])

unique(bouts$PIT_Tag)

#these lines just add 0 to the front of the PIT tag to standardize
#currently PIT Tags from previous experiment. 2017 PIT Tags already have leading 0s
bouts$PIT_Tag[bouts$PIT_Tag=="418257540"]<-"0418257540"
bouts$PIT_Tag[bouts$PIT_Tag=="418257462"]<-"0418257462"
bouts$PIT_Tag[bouts$PIT_Tag=="418252932"]<-"0418252932"
bouts$PIT_Tag[bouts$PIT_Tag=="418252737"]<-"0418252737"
bouts$PIT_Tag[bouts$PIT_Tag=="418257224"]<-"0418257224"
bouts$PIT_Tag[bouts$PIT_Tag=="418257508"]<-"0418257508"
bouts$PIT_Tag[bouts$PIT_Tag=="418257444"]<-"0418257444"
bouts$PIT_Tag[bouts$PIT_Tag=="418256854"]<-"0418256854"
bouts$PIT_Tag[bouts$PIT_Tag=="418253508"]<-"0418253508"
bouts$PIT_Tag[bouts$PIT_Tag=="418257046"]<-"0418257046"
bouts$PIT_Tag[bouts$PIT_Tag=="418256912"]<-"0418256912"
bouts$PIT_Tag[bouts$PIT_Tag=="418256271"]<-"0418256271"
bouts$PIT_Tag[bouts$PIT_Tag=="418255130"]<-"0418255130"
bouts$PIT_Tag[bouts$PIT_Tag=="418256522"]<-"0418256522"
bouts$PIT_Tag[bouts$PIT_Tag=="418257630"]<-"0418257630"
bouts$PIT_Tag[bouts$PIT_Tag=="418253595"]<-"0418253595"
bouts$PIT_Tag[bouts$PIT_Tag=="418256244"]<-"0418256244"
bouts$PIT_Tag[bouts$PIT_Tag=="4182573E7"]<-"04182573E7"

# changed all PIT tags to have a 0 in front of them

write.csv(bouts,"C:/Users/court/Desktop/HOFI Contact Rates/feeding_bouts_4sRule_2017.csv",row.names=F)
# double-check that all PIT tags to have a 0 in front of them in excel
# these bouts used the "NEW PIT TAGS" file update from 2018
#reload bouts data after fixing PIT tag numbers in Excel
bouts<-read.csv("C:/Users/Courtney/Courtney/Documents/House Finches VT/Aviary Contact Rates/2017/RFID/Index data_pre-inoc_TEST/feeding bouts_100617.csv")

unique(bouts$PIT_Tag)

library(doBy)

names(bouts)
sums <- summaryBy(Bout_duration_sec ~ PIT_Tag, FUN=mean, data=bouts)
sums



# DIDN'T TRY ANY OF THE CODE BELOW TODAY
# THAT'S OLD STUFF, BUT MAY STILL BE GOOD.










# from that can we get number of bouts for each bird, avg dur of bout, total time on feeder, # unique feeders
# for all of below: ONLY BIRDS THAT SHOWED UP AT THE FEEDERS, NO
# could add 0s for those bids, by appending a list of those birds that didn't show up...

num_bouts<-tapply(bouts$Bout_duration_s,as.character(bouts$PIT_Tag),length)
# for some reason, we need the as.character, or else it looks at misreads as well as real PIT_Tags...
hist(num_bouts)
length(num_bouts[num_bouts==1 | num_bouts==2])
# this looks right


total_time<-tapply(bouts$Bout_duration_s,as.character(bouts$PIT_Tag),sum)
sort(total_time)
write.csv(total_time,"C:/Users/Courtney/Courtney/Documents/House Finches VT/Aviary Contact Rates/2017/RFID/Index data_pre-inoc_TEST/total feeding time_100617.csv")


mean_bout_length<-tapply(bouts$Bout_duration_s,as.character(bouts$PIT_Tag),mean)
sort(mean_bout_length)

num.unique<-function(x){
  length(unique(x))
}

num_feeders_visited<-tapply(bouts$Site,as.character(bouts$PIT_Tag),num.unique)
hist(num_feeders_visited)

feeder_IDs_visited<-tapply(bouts$Site,as.character(bouts$PIT_Tag),I)
feeder_IDs_visited[dimnames(feeder_IDs_visited)[[1]]=="04179BAE9B" |
  dimnames(feeder_IDs_visited)[[1]]=="04179BBB8A" | dimnames(feeder_IDs_visited)[[1]]=="04179BAA5B"]
# birds that got sick only visited GR and VE

num_days_w_data<-tapply(bouts$Date,as.character(bouts$PIT_Tag),num.unique)
hist(num_days_w_data)

feeders_div_day<-num_feeders_visited/num_days_w_data

feeders_div_day[feeders_div_day==2]


layout(matrix(1:4,nrow=2,ncol=2,byrow=T))
par(xpd=T)
hist(num_bouts/num_days_w_data,col="grey",xlab="Number of feeding bouts/day",
  ylab="Number of birds",main="Number of feeding bouts/day",labels=F,breaks=100)
hist(total_time/num_days_w_data,col="grey",xlab="Total time on feeder/day",
  ylab="Number of birds",main="Total time on feeder/day",labels=F,breaks=100)
hist(mean_bout_length/num_days_w_data,col="grey",xlab="Mean feeding bout length/day",
  ylab="Number of birds",main="Mean feeding bout length/day",labels=F,breaks=100)
hist(num_feeders_visited/num_days_w_data,col="grey",xlab="Number of feeders visited/day",
  ylab="Number of birds",main="Number of feeders visited/day",labels=F,breaks=100)


time_rate<-as.array(total_time/num_days_w_data)
time_rate[dimnames(time_rate)[[1]]=="04179BAE9B" |
  dimnames(time_rate)[[1]]=="04179BBB8A" | dimnames(time_rate)[[1]]=="04179BAA5B"]
time_rate_pcntl<-as.array(rank(time_rate)/max(rank(time_rate)))
time_rate_pcntl[dimnames(time_rate_pcntl)[[1]]=="04179BAE9B" |
  dimnames(time_rate_pcntl)[[1]]=="04179BBB8A" | dimnames(time_rate_pcntl)[[1]]=="04179BAA5B"]

time_rate<-as.array(total_time/num_days_w_data)
sort(time_rate)

num_bouts_rate<-num_bouts/num_days_w_data
sort(num_bouts_rate)

num_bouts_rate_pcntl<-as.array(rank(num_bouts_rate)/max(rank(num_bouts_rate)))
num_bouts_rate_pcntl[dimnames(num_bouts_rate_pcntl)[[1]]=="04179BAE9B" |
  dimnames(num_bouts_rate_pcntl)[[1]]=="04179BBB8A" | dimnames(num_bouts_rate_pcntl)[[1]]=="04179BAA5B"]


feeder_IDs_visited[dimnames(feeder_IDs_visited)[[1]]=="04179BAE9B" |
  dimnames(feeder_IDs_visited)[[1]]=="04179BBB8A" | dimnames(feeder_IDs_visited)[[1]]=="04179BAA5B"]



mean_bout_length_rate<-mean_bout_length/num_days_w_data
mean_bout_length_rate_pcntl<-as.array(rank(mean_bout_length_rate)/max(rank(mean_bout_length_rate)))
mean_bout_length_rate_pcntl[dimnames(mean_bout_length_rate_pcntl)[[1]]=="04179BAE9B" |
  dimnames(mean_bout_length_rate_pcntl)[[1]]=="04179BBB8A" | dimnames(mean_bout_length_rate_pcntl)[[1]]=="04179BAA5B"]


num_feeders_visited_rate<-num_feeders_visited/num_days_w_data
num_feeders_visitedrate_pcntl<-as.array(rank(num_feeders_visited_rate)/max(rank(num_feeders_visited_rate)))
num_feeders_visitedrate_pcntl[dimnames(num_feeders_visitedrate_pcntl)[[1]]=="04179BAE9B" |
  dimnames(num_feeders_visitedrate_pcntl)[[1]]=="04179BBB8A" | dimnames(num_feeders_visitedrate_pcntl)[[1]]=="04179BAA5B"]


# don't know if this is good stats or not, but can you estimate the probability
# of drawing 3 values from one of those distributions and having all 3 be in the top X%?

# how to do 3 random nums without replacement?
output<-matrix(rep(0,300000),nrow=100000)

for (i in 1:100000){
output[i,]<-sample(1:118, 3, replace=F)
}
hist(output)
# so that's uniform, good

# actually doing the sampling

output<-rep(0,10000)

for (i in 1:10000){
  randnums<-sample(1:118,3,replace=F)
  pcntl<-(rank(num_bouts_rate)/max(rank(num_bouts_rate)))[randnums]
  all_top<-ifelse(sum(ifelse(pcntl>0.9,1,0))==3,1,0)
  output[i]<-all_top
}

length(output[output==1])/10000
# 6e-4
# does that mean that the randomly choosing 3 with pcntl > 0.9 is 0.0006?


# time on feeder:
output<-rep(0,10000)

for (i in 1:10000){
  randnums<-sample(1:118,3,replace=F)
  pcntl<-(rank(time_rate)/max(rank(time_rate)))[randnums]
  all_top<-ifelse(sum(ifelse(pcntl>0.85,1,0))==3,1,0)
  output[i]<-all_top
}

length(output[output==1])/10000
# .003
# does that mean probability of randomly choosing 3 with pcntl >0.85 is 0.003?


# bout length:
output<-rep(0,10000)

for (i in 1:10000){
  randnums<-sample(1:118,3,replace=F)
  pcntl<-(rank(mean_bout_length_rate)/max(rank(mean_bout_length_rate)))[randnums]
  all_top<-ifelse(sum(ifelse(pcntl>0.4,1,0))==3,1,0)
  output[i]<-all_top
}

length(output[output==1])/10000
# 0.215

# feeders visited:
output<-rep(0,10000)

for (i in 1:10000){
  randnums<-sample(1:118,3,replace=F)
  pcntl<-(rank(num_feeders_visited_rate)/max(rank(num_feeders_visited_rate)))[randnums]
  all_top<-ifelse(sum(ifelse(pcntl>0.3,1,0))==3,1,0)
  output[i]<-all_top
}

length(output[output==1])/10000
# 0.3317


