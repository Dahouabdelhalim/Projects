#This script splits the original database into the 39 trials included in the BeechCOSTe52

#set working directory 
setwd("/WorkingDirectory/")

#Read databases
Fsyl <- read.csv("Fsylvatica.csv")


#List available trials 
trials <- unique(Fsyl$Trial)
names(trials)= paste0("V",1:38)


#export as .csv all trials
for (i in 1:length(trials)) {

	print(trials[i])
	fs <- Fsyl[(Fsyl$Trial == trials[i]),]
	print(unique(fs$Trial))
	write.csv(fs, file = paste("trials/Fsylvatica_trial_",i,".csv",sep="" ), sep="\\t", col.names=TRUE)
}


	

