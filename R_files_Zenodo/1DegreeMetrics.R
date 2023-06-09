###Get averaged degree, actual degree, infected degree

options(stringsAsFactors = FALSE)

SurveyFile<-read.csv("allsightings.csv")
FocalFile<-read.csv("allfocals.csv")

calves<-unique(FocalFile$dolphin_id)
n<-length(calves)

#format dates
SurveyFile$date<-as.Date(SurveyFile$date)
FocalFile$birthdate<-as.Date(FocalFile$birthdate)

#calculate average degree by resampling sets of 5 surveys
indivs<-split(SurveyFile, SurveyFile$dolphin_id)

results<-list()

set.seed(2019) 

for (i in 1:n){
  
  focal<-calves[i]
  start<-FocalFile$birthdate[which(FocalFile$dolphin_id==focal)]
  end<-start+365
  surveys<-indivs[[focal]]
  trim_surveys<-surveys[(surveys$date>=start & surveys$date<=end),]
  observation_ids<-unique(trim_surveys$observation_id)
  #number of surveys
  total<-length(observation_ids)
  
  if(length(observation_ids)<5) {next} else{
    
    subsample_degree<-replicate(1000,{
      
      sobservation_ids<-sample(observation_ids, 5)
      associates<-SurveyFile$dolphin_id[which(SurveyFile$observation_id %in% sobservation_ids)]
      nassoc<-length(unique(associates[which(associates!=focal)]))
      return(nassoc)})
  }
  
  mean_degree<-mean(unlist(subsample_degree))
  
  #total and infected number of associates
  associates<-SurveyFile[,c("dolphin_id", "infection_status")][which(SurveyFile$observation_id %in% observation_ids),]
  assoc<-aggregate(infection_status~dolphin_id, data=associates, max)
  assoc<-assoc[assoc$dolphin_id!=focal,]

  ind_results<-data.frame(dolphin_id=focal, 
                          avg_degree=mean_degree,
                          true_degree=nrow(assoc),
                          infected_degree=sum(assoc$infection_status),
                          total_surveys=total)
  results[[i]]<-ind_results
}

modeldata<-as.data.frame(do.call("rbind", results))

modeldata$proportion_infected<-modeldata$infected_degree/modeldata$true_degree

#Merge with FocalFile

modeldata<-merge(FocalFile, modeldata, by="dolphin_id")

write.csv(modeldata, "modeldata.csv", row.names = FALSE)
