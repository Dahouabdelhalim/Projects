#Randomize data using curveball algorithm

options(stringsAsFactors = FALSE)

sf<-read.csv("allsightings.csv")
modeldata<-read.csv("modeldata.csv")
modeldatasexed<-modeldata[-which(is.na(modeldata$sex)),]
modeldatasexed$birthdate<-as.Date(modeldatasexed$birthdate)

#set up empty matrices

calves<-sort(modeldatasexed$dolphin_id)
n<-length(calves)

months<-unique(sf$twomonth_interval)
degree_matrix<-matrix(NA, nrow=n, ncol=1000)
infected_degree_matrix<-matrix(NA, nrow=n, ncol=1000)
rownames(degree_matrix)<-calves

tp<-aggregate(observation_id~twomonth_interval, data=sf, function(x) length(unique(x)))
mean(tp$observation_id) #169 sightings per period

source('3aRandomizingFunctions.R')

set.seed(2019)

for (q in 1:1000) {
  
  randsurveys<-make_random(sf)
  
  #now calculate degree and infected degree for all focal individuals
  
  indivs<-split(randsurveys, randsurveys$dolphin_id)
  
  results<-list()
  
  for (i in 1:n){

    focal<-calves[i]
    start<-modeldatasexed$birthdate[which(modeldatasexed$dolphin_id==focal)]
    end<-start+365
    surveys<-indivs[[focal]]
    trim_surveys<-surveys[(surveys$date>=start & surveys$date<=end),]
    observation_ids<-unique(trim_surveys$observation_id)
    
    #total and infected number of associates
    associates<-randsurveys[,c("dolphin_id", "infection_status")][which(randsurveys$observation_id %in% observation_ids),]
    assoc<-aggregate(infection_status~dolphin_id, data=associates, max)
    assoc<-assoc[assoc$dolphin_id!=focal,]
    
    ind_results<-data.frame(dolphin_id=focal, 
                            true_degree=nrow(assoc),
                            infected_degree=sum(assoc$infection_status))
    results[[i]]<-ind_results
    }
  
  fakemodeldata<-as.data.frame(do.call("rbind", results))
  
  degree_matrix[,q]<-fakemodeldata$true_degree
  infected_degree_matrix[,q]<-fakemodeldata$infected_degree
}

proportion_infected_matrix<-infected_degree_matrix/degree_matrix

write.csv(proportion_infected_matrix, "curveball_proportion_infected.csv")
