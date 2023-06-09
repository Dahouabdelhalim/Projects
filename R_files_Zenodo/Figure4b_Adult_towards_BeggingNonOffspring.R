#### Figure 4b
##### Adult->Begging Non-Offspring


##### Load the data
dataQ <- read.csv("dataQ.csv",stringsAsFactors=FALSE, sep=",")





##### Add obd prob and expected prob for this test
dataQ$Expected_Probability_to_feed_Begging_UnrelatedID_FedYES<-dataQ$Num_Begging_Unrelated/(dataQ$Num_Attending_IDUnknown+dataQ$Num_Attending_IDUnrelated)
dataQ$Probability_to_feed_Begging_UnrelatedID_FedYES<-dataQ$Num_Fed_Begging_UnrelatedIDs/(dataQ$Num_Fed_IDUnrelated+dataQ$Num_Fed_Unknown)



##### Prepare the dataframe
dataQ_BeggingALL<-dataQ[which(dataQ$Num_Begging_IDOffspring>0|dataQ$Num_Begging_Unrelated>0|dataQ$Num_NQ.begging >0),]


########## how selective do adults feed among Unrelated ID juveniles?
##### No NQ Duration
##### AMONG Unrelated juveniles with backpack
dataQ_BeggingUnrelated<-dataQ_BeggingALL[which(dataQ_BeggingALL$Num_Begging_Unrelated>0),]
dataQ_BeggingUnrelated_FedYES<-dataQ_BeggingUnrelated[which(dataQ_BeggingUnrelated$fed=="y"),]
dataQ_BeggingUnrelated_FedYES_NoNQDuration<-dataQ_BeggingUnrelated_FedYES[which(dataQ_BeggingUnrelated_FedYES$No_NQOffspring_Period=="y"),]
dataQ_BeggingUnrelated_FedYES_NoNQDuration<-dataQ_BeggingUnrelated_FedYES_NoNQDuration[which(dataQ_BeggingUnrelated_FedYES_NoNQDuration$Num_Attending_IDUnrelated>dataQ_BeggingUnrelated_FedYES_NoNQDuration$Num_Begging_Unrelated),]
dataQ_BeggingUnrelated_FedYES_NoNQDuration<-dataQ_BeggingUnrelated_FedYES_NoNQDuration[which(dataQ_BeggingUnrelated_FedYES_NoNQDuration$Num_Fed_IDUnrelated>0 | dataQ_BeggingUnrelated_FedYES_NoNQDuration$Num_Fed_Unknown>0),]

dataQ_BeggingUnrelated_FedYES_NoNQDuration<-dataQ_BeggingUnrelated_FedYES_NoNQDuration[which(dataQ_BeggingUnrelated_FedYES_NoNQDuration$Expected_Probability_to_feed_Begging_UnrelatedID_FedYES>0&dataQ_BeggingUnrelated_FedYES_NoNQDuration$Expected_Probability_to_feed_Begging_UnrelatedID_FedYES<1),]




##### Basic information
nrow(dataQ_BeggingUnrelated_FedYES_NoNQDuration)#41 #the number of feeding events fulfilling the criteria
length(unique(dataQ_BeggingUnrelated_FedYES_NoNQDuration$AdultID))#11 unique subject adults were observed across aviaries in total
n_juveniles<-as.numeric(unique(unlist(strsplit(dataQ_BeggingUnrelated_FedYES_NoNQDuration$IDchicks_AnyTime," "))))
n_juveniles<-n_juveniles[!is.na(n_juveniles)]
length(n_juveniles)#51 unique juveniles with backpack were observed across aviaries in total




##### [Plot] (these data points do NOT include events in which offspring was attend)
## data plotting
data <- dataQ_BeggingUnrelated_FedYES_NoNQDuration

(x_min<-min(data$Expected_Probability_to_feed_Begging_UnrelatedID_FedYES))
(x_max<-max(data$Expected_Probability_to_feed_Begging_UnrelatedID_FedYES))
newx<-seq(0.05,0.70,0.05)

library(lme4)
mod.obs<-glmer(Probability_to_feed_Begging_UnrelatedID_FedYES ~ Expected_Probability_to_feed_Begging_UnrelatedID_FedYES + (1|Aviary) + (1|AdultID), data=data, family=binomial)
newdf<-with(data, expand.grid(Expected_Probability_to_feed_Begging_UnrelatedID_FedYES =newx, Aviary=unique(Aviary),AdultID=unique(AdultID)))
newdf$y <- predict(mod.obs, newdf, type="response",re.form= ~(1|Aviary)+ (1|AdultID),allow.new.levels = T)
newy <- sapply(by(newdf$y,newdf$Expected_Probability_to_feed_Begging_UnrelatedID_FedYES,mean),identity)
newdf.obs <- newdf

BeggingUnrelated<-summary(mod.obs)

#bootstrapping
boots<-10000#try 10000
xys<-matrix(NA,nrow=boots,ncol=length(newx))
for(i in 1:boots){
  
  mod.rand <- NULL
  while(is.null(mod.rand)) {
    data.rand<-data[sample(1:nrow(data),nrow(data),replace=TRUE),]
    mod.rand = tryCatch({
      glmer(Probability_to_feed_Begging_UnrelatedID_FedYES ~ Expected_Probability_to_feed_Begging_UnrelatedID_FedYES + (1|Aviary) + (1|AdultID), data=data.rand, family=binomial)
    },error = function(e) {
      NULL
    })
  }
  
  newdf<-with(data.rand, expand.grid(Expected_Probability_to_feed_Begging_UnrelatedID_FedYES =newx, Aviary=unique(Aviary),AdultID=unique(AdultID)))
  newdf$y <- predict(mod.rand, newdf, type="response",re.form= ~(1|Aviary)+ (1|AdultID),allow.new.levels = T)
  xys[i,]<-sapply(by(newdf$y,newdf$Expected_Probability_to_feed_Begging_UnrelatedID_FedYES,mean),identity)
  
}
conf.ints<-apply(xys,2,quantile,c(0.025,0.975))

# summary of model
BeggingUnrelated<-summary(mod.obs)

pdf("ProbabilityToFeedBeggingIDUnrelated_AdultWithoutNQoffspringSituation.pdf")
par(mar=c(5, 5, 5, 5))
par(bty="l")
par(lwd = 2)

plot(NULL, xlim=c(0.03,0.72),ylim=c(0,1), xlab="P(begging non-offspring | attending)", ylab="P(begging non-offspring | fed)")
polygon(x=c(newx,rev(newx)), y=c(conf.ints[1,],rev(conf.ints[2,])), col="lightpink1",border="#FF0000", lwd=0.5)
points(data$Expected_Probability_to_feed_Begging_UnrelatedID_FedYES, data$Probability_to_feed_Begging_UnrelatedID_FedYES, pch=21, bg="darkgrey", colour="black", lwd=0.3, cex=0.7)


lines(newx,newy,col="red",lwd=2)
lines(newx,newx,col="black",lwd=2)

counts <- table(round(data$Expected_Probability_to_feed_Begging_UnrelatedID_FedYES*20)/20)
comb <- data.frame(newx=newx, newy=newy, counts=0)
comb$counts <- counts[match(as.character(comb$newx),names(counts))]
comb$counts[is.na(comb$counts)] <- 0
mean.selectivity <- weighted.mean(newy-newx,comb$counts)
CI.selectivity <- c(weighted.mean(conf.ints[1,]-newx,comb$counts), weighted.mean(conf.ints[2,]-newx,comb$counts))
comb$P <- NA
for (i in 1:nrow(comb)) {
  comb$P[i] <- sum((xys[,i]-newx[i])<0)/nrow(xys)
  if (comb$P[i] >= 0.975 | comb$P[i] <= 0.025) {
    text(x=comb$newx[i],y=comb$newy[i],"*",col="red", cex=1.5)
  }
}
P <- weighted.mean(comb$P,comb$counts) #p-value to report
dev.off()

BeggingUnrelated
(mean_selectivity_BeggingUnrelated<-mean(newy-newx))
(CI.selectivity_BeggingUnrelated<-colMeans(t(conf.ints)-newx))
(P_BeggingUnrelated<-P)


BeggingUnrelated_dataframe<-data.frame(EventID=dataQ_BeggingUnrelated_FedYES_NoNQDuration$EventID,Aviary=dataQ_BeggingUnrelated_FedYES_NoNQDuration$Aviary,AdultID=dataQ_BeggingUnrelated_FedYES_NoNQDuration$AdultID,Expected_Probability=dataQ_BeggingUnrelated_FedYES_NoNQDuration$Expected_Probability_to_feed_Begging_UnrelatedID_FedYES,Observed_Probability=dataQ_BeggingUnrelated_FedYES_NoNQDuration$Probability_to_feed_Begging_UnrelatedID_FedYES)
View(BeggingUnrelated_dataframe)




# identify the number of juveniles in dataQ_BeggingUnrelated_FedYES_NoNQDuration 
#Aviary 1
dataQ_BeggingUnrelated_FedYES_NoNQDuration_A1<-dataQ_BeggingUnrelated_FedYES_NoNQDuration[(dataQ_BeggingUnrelated_FedYES_NoNQDuration$Aviary==1),]
begging_unrelatedID_A1<-NULL
unrelatedID_A1<-NULL
for(i in 1:nrow(dataQ_BeggingUnrelated_FedYES_NoNQDuration_A1)){
  temp<-dataQ_BeggingUnrelated_FedYES_NoNQDuration_A1[i,]
  begging_unrelatedID_temp<-NULL
  begging_unrelatedID_temp<-unlist(strsplit(temp$Begging_UnrelatedID," "))
  begging_unrelatedID_A1<-c(begging_unrelatedID_A1,begging_unrelatedID_temp)
  unrelatedID_temp<-NULL
  unrelatedID_temp<-unlist(strsplit(temp$Attending_IDUnrelated," "))
  unrelatedID_A1<-c(unrelatedID_A1,unrelatedID_temp)
}
begging_unrelatedID_A1<-unique(begging_unrelatedID_A1)
unrelatedID_A1<-unique(unrelatedID_A1)
length(begging_unrelatedID_A1) #14
length(unrelatedID_A1)#22
#Aviary 2
dataQ_BeggingUnrelated_FedYES_NoNQDuration_A2<-dataQ_BeggingUnrelated_FedYES_NoNQDuration[(dataQ_BeggingUnrelated_FedYES_NoNQDuration$Aviary==2),]
begging_unrelatedID_A2<-NULL
unrelatedID_A2<-NULL
for(i in 1:nrow(dataQ_BeggingUnrelated_FedYES_NoNQDuration_A2)){
  temp<-dataQ_BeggingUnrelated_FedYES_NoNQDuration_A2[i,]
  begging_unrelatedID_temp<-NULL
  begging_unrelatedID_temp<-unlist(strsplit(temp$Begging_UnrelatedID," "))
  begging_unrelatedID_A2<-c(begging_unrelatedID_A2,begging_unrelatedID_temp)
  unrelatedID_temp<-NULL
  unrelatedID_temp<-unlist(strsplit(temp$Attending_IDUnrelated," "))
  unrelatedID_A2<-c(unrelatedID_A2,unrelatedID_temp)
}
begging_unrelatedID_A2<-unique(begging_unrelatedID_A2)
unrelatedID_A2<-unique(unrelatedID_A2)
length(begging_unrelatedID_A2) #13
length(unrelatedID_A2)#14
#Aviary 3
dataQ_BeggingUnrelated_FedYES_NoNQDuration_A3<-dataQ_BeggingUnrelated_FedYES_NoNQDuration[(dataQ_BeggingUnrelated_FedYES_NoNQDuration$Aviary==3),]
begging_unrelatedID_A3<-NULL
unrelatedID_A3<-NULL
for(i in 1:nrow(dataQ_BeggingUnrelated_FedYES_NoNQDuration_A3)){
  temp<-dataQ_BeggingUnrelated_FedYES_NoNQDuration_A3[i,]
  begging_unrelatedID_temp<-NULL
  begging_unrelatedID_temp<-unlist(strsplit(temp$Begging_UnrelatedID," "))
  begging_unrelatedID_A3<-c(begging_unrelatedID_A3,begging_unrelatedID_temp)
  unrelatedID_temp<-NULL
  unrelatedID_temp<-unlist(strsplit(temp$Attending_IDUnrelated," "))
  unrelatedID_A3<-c(unrelatedID_A3,unrelatedID_temp)
}
begging_unrelatedID_A3<-unique(begging_unrelatedID_A3)
unrelatedID_A3<-unique(unrelatedID_A3)
length(begging_unrelatedID_A3) #2
length(unrelatedID_A3)#6
#Aviary 4
dataQ_BeggingUnrelated_FedYES_NoNQDuration_A4<-dataQ_BeggingUnrelated_FedYES_NoNQDuration[(dataQ_BeggingUnrelated_FedYES_NoNQDuration$Aviary==4),]
begging_unrelatedID_A4<-NULL
unrelatedID_A4<-NULL
for(i in 1:nrow(dataQ_BeggingUnrelated_FedYES_NoNQDuration_A4)){
  temp<-dataQ_BeggingUnrelated_FedYES_NoNQDuration_A4[i,]
  begging_unrelatedID_temp<-NULL
  begging_unrelatedID_temp<-unlist(strsplit(temp$Begging_UnrelatedID," "))
  begging_unrelatedID_A4<-c(begging_unrelatedID_A4,begging_unrelatedID_temp)
  unrelatedID_temp<-NULL
  unrelatedID_temp<-unlist(strsplit(temp$Attending_IDUnrelated," "))
  unrelatedID_A4<-c(unrelatedID_A4,unrelatedID_temp)
}
begging_unrelatedID_A4<-unique(begging_unrelatedID_A4)
unrelatedID_A4<-unique(unrelatedID_A4)
length(begging_unrelatedID_A4) #11
length(unrelatedID_A4)#14
### total
begging_unrelatedID<-c(begging_unrelatedID_A1,begging_unrelatedID_A2,begging_unrelatedID_A3,begging_unrelatedID_A4)
unrelatedID<-c(unrelatedID_A1,unrelatedID_A2,unrelatedID_A3,unrelatedID_A4)
length(begging_unrelatedID)#40
length(unrelatedID)#56

