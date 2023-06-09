#### Model 1 ####
##### Figure 3a: Adult->Own Offspring


##### Load the data
# feeding data
dataQ <- read.csv("dataQ.csv",stringsAsFactors=FALSE, sep=",")





# Weird individuals # i.e. outliers who disproportinally fed non-offpsring
dataQ <- dataQ[which(!(dataQ$AdultID %in% c("1","3","14","144"))),]



dataQ$Expected_Probability_to_feed_IDOffspring_FedYES<-NA
for(i in 1:nrow(dataQ)){
  temp<-dataQ[i,]
  if(temp$fed=="y"){
    if(temp$No_NQOffspring_Period=="y"){
      temp$Expected_Probability_to_feed_IDOffspring_FedYES<-temp$Num_Attending_IDOffspring/(temp$Num_Attending_IDOffspring+temp$Num_Attending_IDUnknown+temp$Num_Attending_IDUnrelated)
      if(is.nan(temp$Expected_Probability_to_feed_IDOffspring_FedYES)){
        temp$Expected_Probability_to_feed_IDOffspring_FedYES<-0
      }
    }else if(temp$No_NQOffspring_Period=="n"){
      temp$Expected_Probability_to_feed_IDOffspring_FedYES<-NA
    }
  }
  dataQ[i,]<-temp
}



##### Prepare the dataframe
dataQ_IDOffspring<-dataQ[!is.na(dataQ$Expected_Probability_to_feed_IDOffspring_FedYES),]#Fed==YES
dataQ_IDOffspring<-dataQ_IDOffspring[which(dataQ_IDOffspring$Num_Attending_IDOffspring>0),]#at least 1 Offspring attending
dataQ_IDOffspring<-dataQ_IDOffspring[which(dataQ_IDOffspring$Num_Attending_IDUnrelated>0|dataQ_IDOffspring$Num_Attending_IDUnknown>0),]#at least 1 Unrelated attending (NQs are all unrelated because the adult does not have any NQ_Offspring at that time)



##### [Observed] Probability
dataQ_IDOffspring$Probability_to_feed_IDOffspring_FedYES<-dataQ_IDOffspring$Num_Fed_IDOffspring/(dataQ_IDOffspring$Num_Fed_IDOffspring+dataQ_IDOffspring$Num_Fed_IDUnrelated+dataQ_IDOffspring$Num_Fed_Unknown)

##### Basic information
nrow(dataQ_IDOffspring)#71 #the number of feeding events fulfilling the criteria
length(unique(dataQ_IDOffspring$AdultID))#24 #the number of unique subject adults were observed across aviaries in total
n_juveniles<-as.numeric(unique(unlist(strsplit(dataQ_IDOffspring$IDchicks_AnyTime," "))))
n_juveniles<-n_juveniles[!is.na(n_juveniles)]
length(n_juveniles)#67 unique juveniles with backpack were observed across aviaries in total



##### Generalized linear mixed model
### Random effects: Aviary and AdultID
library(lme4)

(x_min<-min(dataQ_IDOffspring$Expected_Probability_to_feed_IDOffspring_FedYES))
(x_max<-max(dataQ_IDOffspring$Expected_Probability_to_feed_IDOffspring_FedYES))


newx<-seq(0.05,0.75,0.05)

data <- dataQ_IDOffspring
mod.obs<-glmer(Probability_to_feed_IDOffspring_FedYES ~ Expected_Probability_to_feed_IDOffspring_FedYES + (1|Aviary) + (1|AdultID), data=data, family=binomial)
newdf<-with(data, expand.grid(Expected_Probability_to_feed_IDOffspring_FedYES =newx, Aviary=unique(Aviary),AdultID=unique(AdultID)))
newdf$y <- predict(mod.obs, newdf, type="response",re.form= ~(1|Aviary)+ (1|AdultID))

newy <- sapply(by(newdf$y,newdf$Expected_Probability_to_feed_IDOffspring_FedYES,mean),identity)
newdf.obs <- newdf


#bootstrapping
boots<-10000 #I will try 10000
xys<-matrix(NA,nrow=boots,ncol=length(newx))
for(i in 1:boots){
  mod.rand <- NULL
  while(is.null(mod.rand)) {
    data.rand<-data[sample(1:nrow(data),nrow(data),replace=TRUE),]
    mod.rand = tryCatch({
      glmer(Probability_to_feed_IDOffspring_FedYES ~ Expected_Probability_to_feed_IDOffspring_FedYES + (1|Aviary) + (1|AdultID), data=data.rand, family=binomial)
    },error = function(e) {
      NULL
    })
  }
  newdf<-with(data.rand, expand.grid(Expected_Probability_to_feed_IDOffspring_FedYES =newx, Aviary=unique(Aviary),AdultID=unique(AdultID)))
  newdf$y <- predict(mod.rand, newdf, type="response",re.form= ~(1|Aviary)+ (1|AdultID))
  xys[i,]<-sapply(by(newdf$y,newdf$Expected_Probability_to_feed_IDOffspring_FedYES,mean),identity)
}
conf.ints<-apply(xys,2,quantile,c(0.025,0.975))

# summary of model
OwnOffsp<-summary(mod.obs)






##### [Plot]
##### Black line: the probability expected by chance
##### Red line: the observed probability predicted by the model
##### Dot: each opportunity to selecte a juvenile to feed


par(mar=c(5, 5, 5, 5))
par(bty="l")
par(lwd = 2)

plot(NULL, xlim=c(0.03,0.78),ylim=c(0,1), xlab="P(offspring | attending)", ylab="P(offspring | fed)")
polygon(x=c(newx,rev(newx)), y=c(conf.ints[1,],rev(conf.ints[2,])), col="lightpink1",border="#FF0000", lwd=0.5)
points(data$Expected_Probability_to_feed_IDOffspring_FedYES, data$Probability_to_feed_IDOffspring_FedYES, pch=21, bg="darkgrey", colour="black", lwd=0.3, cex=0.7)
lines(newx,newy,col="red",lwd=2)
lines(newx,newx,col="black",lwd=2)

counts <- table(round(data$Expected_Probability_to_feed_IDOffspring_FedYES*20)/20)
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


OwnOffsp
(mean_selectivity_OwnOffsp<- mean(newy-newx))
(CI.selectivity_OwnOffsp<-colMeans(t(conf.ints)-newx))
(P_OwnOffsp<-P)






#### Model 2 ####
#### Supplemental Information Section 7
##### do juveniles show begging behaviour selectively to their own parents?
### create the dataframe, ($ChickID, $Offspring(Yes or No), $Number_Beg(#Event), $Num_quiet(#Event))
### in how many events did each juvenile show begging to their own parents and/or to non-parents?
### in how many events did each juvenile keep quiet to their own parents and/or to non-parents?

# feeding data
dataQ <- read.csv("dataQ.csv", stringsAsFactors=FALSE, sep=",")

dataQ_NoNQDuration<-dataQ[which(dataQ$No_NQOffspring_Period=="y"),]
All_Observed_IDJuveniles<-NULL
All_Observed_IDJuveniles_Aviary<-NULL
for(i in 1:nrow(dataQ_NoNQDuration)){
  temp<-dataQ_NoNQDuration[i,]
  Attending_IDOffspring<-unlist(strsplit(as.character(temp$Attending_IDOffspring), " "))
  All_Observed_IDJuveniles<-c(All_Observed_IDJuveniles,Attending_IDOffspring)
  All_Observed_IDJuveniles_Aviary<-c(All_Observed_IDJuveniles_Aviary,rep(temp$Aviary,length(Attending_IDOffspring)))
}

Selective_Begging<-as.data.frame(matrix(nrow = length(All_Observed_IDJuveniles_Aviary), ncol = 5))
colnames(Selective_Begging)<-c("Aviary","ChickID", "Offspring", "Number_Beg", "Num_quiet")
Selective_Begging$Aviary<-All_Observed_IDJuveniles_Aviary
Selective_Begging$ChickID<-All_Observed_IDJuveniles
Selective_Begging$ChickID[which(Selective_Begging$ChickID=="NA")]<-NA
Selective_Begging<-Selective_Begging[which(!is.na(Selective_Begging$ChickID)),]
Selective_Begging<-unique(Selective_Begging)
Selective_Begging<-rbind(Selective_Begging,Selective_Begging)
Selective_Begging$Offspring[1:(nrow(Selective_Begging)/2)]<-"y"
Selective_Begging$Offspring[(nrow(Selective_Begging)/2+1):(nrow(Selective_Begging))]<-"n"
Selective_Begging$Number_Beg<-0
Selective_Begging$Num_quiet<-0


for(i in 1:nrow(dataQ_NoNQDuration)){
  temp<-dataQ_NoNQDuration[i,]
  Aviary<-temp$Aviary
  Begging_IDOffspring<-unlist(strsplit(as.character(temp$Begging_IDOffspring), " "))
  Begging_IDUnrelated<-unlist(strsplit(as.character(temp$Begging_UnrelatedID), " "))
  Quiet_IDOffspring<-unlist(strsplit(as.character(temp$Quiet_IDOffspring), " "))
  Quiet_IDUnrelated<-unlist(strsplit(as.character(temp$Quiet_IDUnrelated), " "))
  
  if(length(Begging_IDOffspring)>0){
    if(!is.na(Begging_IDOffspring)){
      for(j in 1:length(Begging_IDOffspring)){
        Selective_Begging$Number_Beg[which(Selective_Begging$Aviary==Aviary&Selective_Begging$ChickID==Begging_IDOffspring[j]&Selective_Begging$Offspring=="y")]<-Selective_Begging$Number_Beg[which(Selective_Begging$Aviary==Aviary&Selective_Begging$ChickID==Begging_IDOffspring[j]&Selective_Begging$Offspring=="y")]+1
      }
    }
  }
  if(length(Begging_IDUnrelated)>0){
    if(!is.na(Begging_IDUnrelated)){
      for(j in 1:length(Begging_IDUnrelated)){
        Selective_Begging$Number_Beg[which(Selective_Begging$Aviary==Aviary&Selective_Begging$ChickID==Begging_IDUnrelated[j]&Selective_Begging$Offspring=="n")]<-Selective_Begging$Number_Beg[which(Selective_Begging$Aviary==Aviary&Selective_Begging$ChickID==Begging_IDUnrelated[j]&Selective_Begging$Offspring=="n")]+1
      }
    }
  }
  if(length(Quiet_IDOffspring)>0){
    if(!is.na(Quiet_IDOffspring)){
      for(j in 1:length(Quiet_IDOffspring)){
        Selective_Begging$Num_quiet[which(Selective_Begging$Aviary==Aviary&Selective_Begging$ChickID==Quiet_IDOffspring[j]&Selective_Begging$Offspring=="y")]<-Selective_Begging$Num_quiet[which(Selective_Begging$Aviary==Aviary&Selective_Begging$ChickID==Quiet_IDOffspring[j]&Selective_Begging$Offspring=="y")]+1
      }
    }
  }
  if(length(Quiet_IDUnrelated)>0){
    if(!is.na(Quiet_IDUnrelated)){
      for(j in 1:length(Quiet_IDUnrelated)){
        Selective_Begging$Num_quiet[which(Selective_Begging$Aviary==Aviary&Selective_Begging$ChickID==Quiet_IDUnrelated[j]&Selective_Begging$Offspring=="n")]<-Selective_Begging$Num_quiet[which(Selective_Begging$Aviary==Aviary&Selective_Begging$ChickID==Quiet_IDUnrelated[j]&Selective_Begging$Offspring=="n")]+1
      }
    }
  }
}


##### Basic information -- how many adults/juveniles? 
Aviary_Selective_Begging<-Selective_Begging[Selective_Begging$Aviary=="1",]
A1<-length(unique(Aviary_Selective_Begging$ChickID))#19 chicks in aviary 1
Aviary_Selective_Begging<-Selective_Begging[Selective_Begging$Aviary=="2",]
A2<-length(unique(Aviary_Selective_Begging$ChickID))#24 chicks in aviary 2
Aviary_Selective_Begging<-Selective_Begging[Selective_Begging$Aviary=="3",]
A3<-length(unique(Aviary_Selective_Begging$ChickID))#18 chicks in aviary 3
Aviary_Selective_Begging<-Selective_Begging[Selective_Begging$Aviary=="4",]
A4<-length(unique(Aviary_Selective_Begging$ChickID))#24 chicks in aviary 4
A1+A2+A3+A4#90 chicks across aviaries
length(dataQ_NoNQDuration$EventID)#1267 feeding events were included into this analysis #those events occurred in the period in which the subject adults did not have any offspring without backpack(s)
length(unique(dataQ_NoNQDuration$AdultID[which(dataQ_NoNQDuration$Aviary=="1")]))+length(unique(dataQ_NoNQDuration$AdultID[which(dataQ_NoNQDuration$Aviary=="2")]))+length(unique(dataQ_NoNQDuration$AdultID[which(dataQ_NoNQDuration$Aviary=="3")]))+length(unique(dataQ_NoNQDuration$AdultID[which(dataQ_NoNQDuration$Aviary=="4")]))#76 unique adults





#####ã€€Run glm model (offspringID as a random effect) 
library(lme4)
model_begging <- glmer(cbind(Number_Beg, Num_quiet) ~ Offspring + (1|ChickID) , data=Selective_Begging, family=binomial) #model2
OffspBegging<-summary(model_begging)

effect.size <- exp(fixef(model_begging)[2]) ## odds ratio - how much more likely is a chick to beg to a parent than a non-parent

#bootstrapping
boots<-10000#try 10000
effect.size.boot<-rep(NA,boots)
for(i in 1:boots){
  model_begging_boot<-NULL
  chickIDs <- unique(Selective_Begging$ChickID)
  while(is.null(model_begging_boot)) {
    chickIDs.boot <- sample(chickIDs,length(chickIDs),replace=TRUE)
    data.boot<-Selective_Begging[which(Selective_Begging$ChickID %in% chickIDs.boot),]
    model_begging_boot=tryCatch({
      glmer(cbind(Number_Beg, Num_quiet) ~ Offspring + (1|ChickID) , data=data.boot, family=binomial)
    },error = function(e) {
      NULL
    })
  }
  effect.size.boot[i]<-exp(fixef(model_begging_boot)[2])
}
(P <- sum(effect.size.boot<=1)/boots)
(conf.ints<-quantile(effect.size.boot,c(0.025,0.975)))

OffspBegging
(CI.selectivity_OffspBegging<-conf.ints)
(P_OffspBegging<-P)





Selective_Begging$Percentage_beg <- Selective_Begging$Number_Beg/(Selective_Begging$Num_quiet+Selective_Begging$Number_Beg)*100
#towards parents
to_Parents <- Selective_Begging$Percentage_beg[which(Selective_Begging$Offspring=="y")]
range(to_Parents)
unique(ave(to_Parents))
#towards non-parents
to_NonParents <- Selective_Begging$Percentage_beg[which(Selective_Begging$Offspring=="n")]
range(to_NonParents[!is.nan(to_NonParents)])
unique(ave(to_NonParents[!is.nan(to_NonParents)]))

