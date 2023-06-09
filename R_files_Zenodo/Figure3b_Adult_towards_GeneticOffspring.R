#### Model 3
##### Figure 3b #####
##### Adult->Genetic Offspring


##### Load the data
# feeding data
dataQ <- read.csv("dataQ.csv",stringsAsFactors=FALSE, sep=",")





dataQ$Expected_Probability_to_feed_GeneticOffspring_FedYES<-NA
for(i in 1:nrow(dataQ)){
  temp<-dataQ[i,]
  if(temp$fed=="y"){
    if(temp$No_NQGeneticOffspring_Period=="y"){
      temp$Expected_Probability_to_feed_GeneticOffspring_FedYES<-temp$Num_Attending_GeneticOffspring/(temp$Num_Attending_IDUnknown+temp$Num_Attending_IDUnrelated)
      
      if(is.nan(temp$Expected_Probability_to_feed_GeneticOffspring_FedYES)){
        temp$Expected_Probability_to_feed_GeneticOffspring_FedYES<-0
      }
    }else if(temp$No_NQGeneticOffspring_Period=="n"){
      temp$Expected_Probability_to_feed_GeneticOffspring_FedYES<-NA
    }
  }
  dataQ[i,]<-temp
}



##### Prepare the dataframe
dataQ_Genetic<-dataQ[!is.na(dataQ$Expected_Probability_to_feed_GeneticOffspring_FedYES),]#NoNQGeneticOffspring==YES
dataQ_Genetic<-dataQ_Genetic[which(dataQ_Genetic$fed=="y"),]
dataQ_Genetic<-dataQ_Genetic[which(dataQ_Genetic$No_NQGeneticOffspring_Period=="y"),]
dataQ_Genetic<-dataQ_Genetic[which(dataQ_Genetic$Num_Attending_GeneticOffspring>0),]#at least 1 Genetic Offspring attending
dataQ_Genetic<-dataQ_Genetic[which(((dataQ_Genetic$Num_Attending_IDUnrelated+dataQ_Genetic$Num_Attending_IDUnknown+dataQ_Genetic$Num_Attending_IDOffspring)-dataQ_Genetic$Num_Attending_GeneticOffspring)>0),]#at least 1 Unrelated attending (NQs are all unrelated because the adult does not have any NQ_Offspring at that time)
dataQ_Genetic<-dataQ_Genetic[which(dataQ_Genetic$Num_Fed_IDUnrelated>0 | dataQ_Genetic$Num_Fed_Unknown>0),]#at least 1 Unrelated was fed

##### [Observed] Probability
dataQ_Genetic$Probability_to_feed_GeneticOffspring_FedYES<-dataQ_Genetic$Num_Fed_GeneticOffspring/(dataQ_Genetic$Num_Fed_IDUnrelated+dataQ_Genetic$Num_Fed_Unknown)



##### Basic information -- how many events? how many unique subject adults?
nrow(dataQ_Genetic) #41 #the number of feeding events fulfilling the criteria
length(unique(dataQ_Genetic$AdultID))  #9 unique subject adults were observed across aviaries in total
AdultSex <- data.frame(dataQ_Genetic$AdultID,dataQ_Genetic$AdultsSex)
AdultSex <- table(AdultSex) # 4 females # 5 males
n_juveniles<-as.numeric(unique(unlist(strsplit(dataQ_Genetic$IDchicks_AnyTime," "))))
n_juveniles<-n_juveniles[!is.na(n_juveniles)]
length(n_juveniles)  #41 unique juveniles with backpack were observed across aviaries in total
table(dataQ_Genetic$Num_Fed_GeneticOffspring)# 0:22, 1:18, 2:1
0*22+1*18+2*1#20 events to genetic offspring


##### Generalized linear mixed model
### Random effects: Aviary and AdultID
library(lme4)

(x_min<-min(data$Expected_Probability_to_feed_GeneticOffspring_FedYES))
(x_max<-max(data$Expected_Probability_to_feed_GeneticOffspring_FedYES))

newx<-seq(0.05,1,0.05)

data <- dataQ_Genetic
mod.obs<-glmer(Probability_to_feed_GeneticOffspring_FedYES ~ Expected_Probability_to_feed_GeneticOffspring_FedYES + (1|Aviary) + (1|AdultID), data=data, family=binomial)
summary(mod.obs)
newdf<-with(data, expand.grid(Expected_Probability_to_feed_GeneticOffspring_FedYES =newx, Aviary=unique(Aviary),AdultID=unique(AdultID)))
newdf$y <- predict(mod.obs, newdf, type="response",re.form= ~(1|Aviary)+ (1|AdultID))

newy <- sapply(by(newdf$y,newdf$Expected_Probability_to_feed_GeneticOffspring_FedYES,mean),identity)
newdf.obs <- newdf


#bootstrapping
boots<-10000 #I will try 10000
xys<-matrix(NA,nrow=boots,ncol=length(newx))
for(i in 1:boots){
  mod.rand <- NULL
  while(is.null(mod.rand)) {
    data.rand<-data[sample(1:nrow(data),nrow(data),replace=TRUE),]
    mod.rand = tryCatch({
      glmer(Probability_to_feed_GeneticOffspring_FedYES ~ Expected_Probability_to_feed_GeneticOffspring_FedYES + (1|Aviary) + (1|AdultID), data=data.rand, family=binomial)
    },error = function(e) {
      NULL
    })
  }
  newdf<-with(data.rand, expand.grid(Expected_Probability_to_feed_GeneticOffspring_FedYES =newx, Aviary=unique(Aviary),AdultID=unique(AdultID)))
  newdf$y <- predict(mod.rand, newdf, type="response",re.form= ~(1|Aviary)+ (1|AdultID))
  xys[i,]<-sapply(by(newdf$y,newdf$Expected_Probability_to_feed_GeneticOffspring_FedYES,mean),identity)
}
conf.ints<-apply(xys,2,quantile,c(0.025,0.975))

# summary of model
GeneticOffsp<-summary(mod.obs)





##### [Plot]
##### Black line: the probability expected by chance
##### Red line: the observed probability predicted by the model
##### Dot: each opportunity to selecte a juvenile to feed
pdf("ProbabilityToFeedGeneticOffspring_AdultWithoutNQoffspringSituation_weirdomarked.pdf")
par(mar=c(5, 5, 5, 5))
par(bty="l")
par(lwd = 2)

plot(NULL, xlim=c(0.03,1.02),ylim=c(0,1), xlab="P(genetic offspring | attending)", ylab="P(genetic offspring | fed)")
polygon(x=c(newx,rev(newx)), y=c(conf.ints[1,],rev(conf.ints[2,])), col="lightpink1",border="#FF0000", lwd=0.5)
points(data$Expected_Probability_to_feed_GeneticOffspring_FedYES, data$Probability_to_feed_GeneticOffspring_FedYES, pch=21, bg="darkgrey", colour="black", lwd=0.3, cex=0.7)


lines(newx,newy,col="red",lwd=2)
lines(newx,newx,col="black",lwd=2)

counts <- table(round(data$Expected_Probability_to_feed_GeneticOffspring_FedYES*20)/20)
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


GeneticOffsp
(mean_selectivity_GeneticOffsp<-mean(newy-newx))
(CI.selectivity_GeneticOffsp<-colMeans(t(conf.ints)-newx))
(P_GeneticOffsp<-P)


GeneticOffsp_dataframe<-data.frame(EventID=dataQ_Genetic$EventID,Aviary=dataQ_Genetic$Aviary,AdultID=dataQ_Genetic$AdultID,Expected_Probability=dataQ_Genetic$Expected_Probability_to_feed_GeneticOffspring_FedYES,Observed_Probability=dataQ_Genetic$Probability_to_feed_GeneticOffspring_FedYES)
View(GeneticOffsp_dataframe)






#### Supplemental Information Section 9 ####
### do juveniles show begging behaviour selectively to their genetic parents? ###
### create the dataframe, ($ChickID, $Offspring(Yes or No), $Number_Beg(#Event), $Num_quiet(#Event))
### in how many events did each juvenile show begging to their genetic parents and/or to non-genetic parents?
### in how many events did each juvenile keep quiet to their genetic parents and/or to non-genetic parents?

# feeding data
dataQ <- read.csv("dataQ.csv",stringsAsFactors=FALSE, sep=",")

dataQ_NoNQDuration<-dataQ[which(dataQ$No_NQGeneticOffspring_Period=="y"),]
All_Observed_IDJuveniles<-NULL
All_Observed_IDJuveniles_Aviary<-NULL
for(i in 1:nrow(dataQ_NoNQDuration)){
  temp<-dataQ_NoNQDuration[i,]
  Attending_GeneticOffspring<-unlist(strsplit(as.character(temp$Attending_GeneticOffspring), " "))
  All_Observed_IDJuveniles<-c(All_Observed_IDJuveniles,Attending_GeneticOffspring)
  All_Observed_IDJuveniles_Aviary<-c(All_Observed_IDJuveniles_Aviary,rep(temp$Aviary,length(Attending_GeneticOffspring)))
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



dataQ_NoNQDuration$Begging_GeneticOffspring <- NA
dataQ_NoNQDuration$Begging_NonGeneticOffspring <- NA
dataQ_NoNQDuration$Quiet_GeneticOffspring <- NA
dataQ_NoNQDuration$Quiet_NonGeneticOffspring <- NA
for(i in 1:nrow(dataQ_NoNQDuration)){
  temp <- dataQ_NoNQDuration[i,]
  
  Attending_GeneticOffspring <- unlist(strsplit(temp$Attending_GeneticOffspring, " "))
  IDchicks_Begging <- unlist(strsplit(temp$IDchicks_Begging, " "))
  Begging_GeneticOffspring <- intersect(Attending_GeneticOffspring,IDchicks_Begging)
  Begging_NonGeneticOffspring <- setdiff(IDchicks_Begging,Attending_GeneticOffspring)
  Quiet_GeneticOffspring <- setdiff(Attending_GeneticOffspring,Begging_GeneticOffspring)
  Quiet_NonGeneticOffspring <- setdiff(Attending_GeneticOffspring,IDchicks_Begging)
  
  if(length(Begging_GeneticOffspring>0)){
    temp$Begging_GeneticOffspring <- paste0(Begging_GeneticOffspring, collapse = " ")
  }
  if(length(Begging_NonGeneticOffspring>0)){
    temp$Begging_NonGeneticOffspring <- paste0(Begging_NonGeneticOffspring, collapse = " ")
  }
  if(length(Quiet_GeneticOffspring>0)){
    temp$Quiet_GeneticOffspring <- paste0(Quiet_GeneticOffspring, collapse = " ")
  }
  if(length(Quiet_NonGeneticOffspring>0)){
    temp$Quiet_NonGeneticOffspring <- paste0(Quiet_NonGeneticOffspring, collapse = " ")
  }
  
  dataQ_NoNQDuration[i,] <- temp
  
}



for(i in 1:nrow(dataQ_NoNQDuration)){
  temp<-dataQ_NoNQDuration[i,]
  Aviary<-temp$Aviary
  Begging_IDOffspring<-unlist(strsplit(as.character(temp$Begging_GeneticOffspring), " "))
  Begging_IDUnrelated<-unlist(strsplit(as.character(temp$Begging_NonGeneticOffspring), " "))
  Quiet_IDOffspring<-unlist(strsplit(as.character(temp$Quiet_GeneticOffspring), " "))
  Quiet_IDUnrelated<-unlist(strsplit(as.character(temp$Quiet_NonGeneticOffspring), " "))
  
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
A1<-length(unique(Aviary_Selective_Begging$ChickID))#how many chicks in aviary 1
Aviary_Selective_Begging<-Selective_Begging[Selective_Begging$Aviary=="2",]
A2<-length(unique(Aviary_Selective_Begging$ChickID))#in aviary 2
Aviary_Selective_Begging<-Selective_Begging[Selective_Begging$Aviary=="3",]
A3<-length(unique(Aviary_Selective_Begging$ChickID))#in aviary 3
Aviary_Selective_Begging<-Selective_Begging[Selective_Begging$Aviary=="4",]
A4<-length(unique(Aviary_Selective_Begging$ChickID))#in aviary 4
A1+A2+A3+A4
length(dataQ_NoNQDuration$EventID)
length(unique(dataQ_NoNQDuration$AdultID[which(dataQ_NoNQDuration$Aviary=="1")]))+length(unique(dataQ_NoNQDuration$AdultID[which(dataQ_NoNQDuration$Aviary=="2")]))+length(unique(dataQ_NoNQDuration$AdultID[which(dataQ_NoNQDuration$Aviary=="3")]))+length(unique(dataQ_NoNQDuration$AdultID[which(dataQ_NoNQDuration$Aviary=="4")]))





#####ã€€Run glm model (offspringID as a random effect) 
library(lme4)
model_begging <- glmer(cbind(Number_Beg, Num_quiet) ~ Offspring + (1|ChickID) , data=Selective_Begging, family=binomial) #model
OffspBegging<-summary(model_begging)

effect.size <- exp(fixef(model_begging)[2]) ## odds ratio - how much more likely is a chick to beg to a genetic parent than to others

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
CI.selectivity_OffspBegging<-conf.ints
P_OffspBegging<-P

OffspBegging
CI.selectivity_OffspBegging
P_OffspBegging

