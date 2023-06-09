#### 
##### Figure 4a
##### Adult->Under-developed Non-offspring (i.e. non-offspring without backpack)


##### Load the data
# feeding data
dataQ <- read.csv("dataQ.csv",stringsAsFactors=FALSE, sep=",")




dataQ$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES<-NA
for(i in 1:nrow(dataQ)){
  temp<-dataQ[i,]
  
    if(temp$fed=="y"){
      if(temp$No_NQOffspring_Period=="y"){
        temp$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES<-temp$Num_Attending_IDUnrelated/(temp$Num_Attending_IDUnrelated+temp$Num_Attending_IDUnknown)
        if(is.nan(temp$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES)){
          temp$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES<-0
        }
      }else if(temp$No_NQOffspring_Period=="n"){
        temp$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES<-NA
      }
    }
  
  dataQ[i,]<-temp
}




##### Prepare the dataframe
dataQ_BackpackUnrelated<-dataQ[which(dataQ$fed=="y"),]#Fed==YES
dataQ_BackpackUnrelated<-dataQ[!is.na(dataQ$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES),]#Fed==YES
dataQ_BackpackUnrelated<-dataQ_BackpackUnrelated[which(dataQ_BackpackUnrelated$Num_Attending_IDUnknown>0&dataQ_BackpackUnrelated$Num_Attending_IDUnrelated>0),]#at least 1 unrelated attending (NQs are all unrelated because the adult does not have any NQ_Offspring at that time)
dataQ_BackpackUnrelated<-dataQ_BackpackUnrelated[which(dataQ_BackpackUnrelated$Num_Fed_IDUnrelated>0 | dataQ_BackpackUnrelated$Num_Fed_Unknown>0),]#at least 1 non-offspring was fed

##### [Observed] Probability
dataQ_BackpackUnrelated$Probability_to_feed_WithBackpackUnrelated_FedYES<-dataQ_BackpackUnrelated$Num_Fed_IDUnrelated/(dataQ_BackpackUnrelated$Num_Fed_IDUnrelated+dataQ_BackpackUnrelated$Num_Fed_Unknown)



##### Basic information
nrow(dataQ_BackpackUnrelated)#92 #the number of feeding events fulfilling the criteria
length(unique(dataQ_BackpackUnrelated$AdultID))#17 unique subject adults were observed across aviaries in total
n_juveniles<-as.numeric(unique(unlist(strsplit(dataQ_BackpackUnrelated$IDchicks_AnyTime," "))))
n_juveniles<-n_juveniles[!is.na(n_juveniles)]
length(n_juveniles)#57 unique juveniles with backpack were observed across aviaries in total


library(lme4)
## data plotting
data <- dataQ_BackpackUnrelated
(x_min<-min(data$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES))
(x_max<-max(data$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES))
newx<-seq(0.10,0.90,0.05)


mod.obs<-glmer(Probability_to_feed_WithBackpackUnrelated_FedYES ~ Expected_Probability_to_feed_WithBackpackUnrelated_FedYES + (1|Aviary) + (1|AdultID), data=data, family=binomial)

newdf<-with(data, expand.grid(Expected_Probability_to_feed_WithBackpackUnrelated_FedYES =newx, AdultID=unique(AdultID)))
newdf$y <- predict(mod.obs, newdf, type="response",re.form= ~(1|AdultID),allow.new.levels = T)
newy <- sapply(by(newdf$y,newdf$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES,mean),identity)
newdf.obs <- newdf



#bootstrapping
boots<-10000#try 100000
xys<-matrix(NA,nrow=boots,ncol=length(newx))
for(i in 1:boots){
  mod.rand <- NULL
  while(is.null(mod.rand)) {
    data.rand<-data[sample(1:nrow(data),nrow(data),replace=TRUE),]
    mod.rand = tryCatch({
      glmer(Probability_to_feed_WithBackpackUnrelated_FedYES ~ Expected_Probability_to_feed_WithBackpackUnrelated_FedYES + (1|AdultID), data=data.rand, family=binomial)
    },error = function(e) {
      NULL
    })
  }
  newdf<-with(data.rand, expand.grid(Expected_Probability_to_feed_WithBackpackUnrelated_FedYES =newx, AdultID=unique(AdultID)))

  newdf$y <- predict(mod.rand, newdf, type="response",re.form= ~(1|AdultID),allow.new.levels = T)
  xys[i,]<-sapply(by(newdf$y,newdf$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES,mean),identity)

}
conf.ints<-apply(xys,2,quantile,c(0.025,0.975))

# summary of model
DevNonOppsp_OfffspAttending<-summary(mod.obs)


pdf("Probability_to_feed_WithBackpackUnrelated_WithoutNQoffspringSituation_OffspringAttending.pdf")
par(mar=c(5, 5, 5, 5))
par(bty="l")
par(lwd = 2)

plot(NULL, xlim=c(0.08,0.92),ylim=c(0,1), xlab="P(non-offspring with backpack | attending)", ylab="P(non-offspring with backpack | fed)")
polygon(x=c(newx,rev(newx)), y=c(conf.ints[1,],rev(conf.ints[2,])), col="lightpink1",border="#FF0000", lwd=0.5)
points(data$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES, data$Probability_to_feed_WithBackpackUnrelated_FedYES, pch=21, bg="darkgrey", colour="black", lwd=0.3, cex=0.7)



lines(newx,newy,col="red",lwd=2)
lines(newx,newx,col="black",lwd=2)


counts <- table(round(data$Expected_Probability_to_feed_WithBackpackUnrelated_FedYES*20)/20)

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


DevNonOppsp_OfffspAttending
(mean_selectivity_DevNonOppsp_OfffspAttending<-mean(newy-newx))
(CI.selectivity_DevNonOppsp_OfffspAttending<-colMeans(t(conf.ints)-newx))
(P_DevNonOppsp_OfffspAttending<-P)
DevNonOppsp_OfffspAttending_dataframe<-data.frame(EventID=dataQ_BackpackUnrelated$EventID,Aviary=dataQ_BackpackUnrelated$Aviary,AdultID=dataQ_BackpackUnrelated$AdultID,Expected_Probability=dataQ_BackpackUnrelated$Expected_Probability_to_feed_GeneticOffspring_FedYES,Observed_Probability=dataQ_BackpackUnrelated$Probability_to_feed_GeneticOffspring_FedYES)
View(DevNonOppsp_OfffspAttending_dataframe)




#### model5 ####
### Do under-developed juveniles show begging behaviour to non-parents more than developed juveniles? 

### create the dataframe, ($EventID, $Under_developed(Yes(=without backpack) or No(=with backpackID)), $Number_Beg(#juvenile), $Number_quiet(#juvenile))
# feeding events in which there was no offspring without backpack (for each subject adult)
dataQ_NoNQDuration<-dataQ[which(dataQ$No_NQOffspring_Period=="y"),]
# feeding events in which there was at least one non-offspring with backpack or at least one non-offspring without backpack
dataQ_NoNQDuration_NONoffspring<-dataQ_NoNQDuration[which(dataQ_NoNQDuration$Num_Attending_IDUnrelated>0 | dataQ_NoNQDuration$Num_NQ.anytime>0),]

# new dataframe: ($EventID, $Under_developed(Yes(=without backpack) or No(=with backpackID)), $Number_Beg(#juvenile), $Number_quiet(#juvenile))
All_Events<-dataQ_NoNQDuration_NONoffspring$EventID #should be the same length as unique(dataQ_NoNQDuration_NONoffspring$EventID)
All_Aviary<-dataQ_NoNQDuration_NONoffspring$Aviary
Dev_Beg<-as.data.frame(matrix(nrow = length(All_Events), ncol = 5))
colnames(Dev_Beg)<-c("Aviary","EventID", "Under_developed", "Number_Beg", "Num_quiet")
Dev_Beg$EventID<-All_Events
Dev_Beg$Aviary<-All_Aviary
Dev_Beg<-rbind(Dev_Beg,Dev_Beg)
Dev_Beg$Under_developed[1:(nrow(Dev_Beg)/2)]<-"y"#without backpack (i.e. NQ)
Dev_Beg$Under_developed[(nrow(Dev_Beg)/2+1):(nrow(Dev_Beg))]<-"n" #with backpack (i.e. UnrelatedID)

for(i in 1:nrow(dataQ_NoNQDuration_NONoffspring)){
  temp<-dataQ_NoNQDuration_NONoffspring[i,]
  Num_Attending_NQ<-temp$Num_NQ.anytime
  Num_Attending_IDUnrelated<-temp$Num_Attending_IDUnrelated
  
  Num_begging_NQ<-temp$Num_NQ.begging
  Num_begging_UnrelatedID<-temp$Num_Begging_Unrelated
  Num_quiet_NQ<-Num_Attending_NQ-Num_begging_NQ
  Num_quiet_UnrelatedID<-Num_Attending_IDUnrelated-Num_begging_UnrelatedID
  
  Event<-temp$EventID
  Aviary<-temp$Aviary
  
  Dev_Beg$Number_Beg[which(Dev_Beg$Aviary==Aviary & Dev_Beg$EventID==Event & Dev_Beg$Under_developed=="y")]<-Num_begging_NQ
  Dev_Beg$Number_Beg[which(Dev_Beg$Aviary==Aviary & Dev_Beg$EventID==Event & Dev_Beg$Under_developed=="n")]<-Num_begging_UnrelatedID
  Dev_Beg$Num_quiet[which(Dev_Beg$Aviary==Aviary & Dev_Beg$EventID==Event & Dev_Beg$Under_developed=="y")]<-Num_quiet_NQ
  Dev_Beg$Num_quiet[which(Dev_Beg$Aviary==Aviary & Dev_Beg$EventID==Event & Dev_Beg$Under_developed=="n")]<-Num_quiet_UnrelatedID
}
Dev_Beg<-Dev_Beg[which(Dev_Beg$Num_quiet>=0),]#On purposly, excluding the events I did make mistakes during the observation as much as I can
Dev_Beg<-Dev_Beg[which(Dev_Beg$Number_Beg>=0),]#On purposly, excluding the events I did make mistakes during the observation as much as I can



##### Basic information -- how many adults/juveniles? 
Aviary_Dev_Beg<-Dev_Beg[Dev_Beg$Aviary=="1",]
(A1<-length(unique(Aviary_Dev_Beg$EventID)))#in aviary 1
Aviary_Dev_Beg<-Dev_Beg[Dev_Beg$Aviary=="2",]
(A2<-length(unique(Aviary_Dev_Beg$EventID)))#in aviary 2
Aviary_Dev_Beg<-Dev_Beg[Dev_Beg$Aviary=="3",]
(A3<-length(unique(Aviary_Dev_Beg$EventID)))#in aviary 3
Aviary_Dev_Beg<-Dev_Beg[Dev_Beg$Aviary=="4",]
(A4<-length(unique(Aviary_Dev_Beg$EventID)))#in aviary 4



#####ã€€Run glm model (offspringID as a random effect) 
model_dev_begging <- glmer(cbind(Number_Beg, Num_quiet) ~ Under_developed + (1|EventID) , data=Dev_Beg, family=binomial) #model5

effect.size <- exp(fixef(model_dev_begging)[2]) ## odds ratio - how much more likely is a chick to beg to a parent than a non-parent

DevNonOppspBegging_OfffspAttending<-summary(model_dev_begging)

#bootstrapping
boots<-10000#try 10000
effect.size.boot<-rep(NA,boots)
for(i in 1:boots){
  model_dev_begging_boot<-NULL
  events <- unique(Dev_Beg$EventID)
  while(is.null(model_dev_begging_boot)) {
    events.boot <- sample(events,length(events)*2,replace=TRUE) # keep this if using small data, assuming roughly half chicks have no backpack
    data.boot<-Dev_Beg[which(Dev_Beg$EventID %in% events.boot),]
    model_dev_begging_boot=tryCatch({
      glmer(cbind(Number_Beg, Num_quiet) ~ Under_developed + (1|EventID) , data=data.boot, family=binomial)
    },error = function(e) {
      NULL
    })
  }
  effect.size.boot[i]<-exp(fixef(model_dev_begging_boot)[2])
}


(model_dev_begging<-summary(model_dev_begging))
(CI.selectivity_model_dev_begging<-quantile(effect.size.boot,c(0.025,0.975)))
(P_model_dev_begging<-sum(effect.size.boot<=1)/boots)

