## LOAD EVERYTHING IN
setwd('/Volumes/data/BCM/EyeTracking_Movies/AQ/Analysis')

# load in trial by trial data for movie
movie_data_70 = read.csv('/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Analysis/Mouth_Eye_Movie_Movie.txt')

# load in trial by trial data for clear ALL (98)
clear_data_98 = read.csv('/Volumes/data/BCM/EyeTracking_Movies/AQ/Analysis/AllExpt_ByTrial.csv'); 

# load in AQ info for movie data
# (there are 105 here - corresponds to new num col above - but only 98 after
# taking out the duplicates & no shows to the eye tracking portion)
AQ_mov = read.csv('/Volumes/data/BCM/EyeTracking_Movies/AQ/AQ_Score_Movies.csv');
movie_data_70$AQ = AQ_mov 

library(lme4)
library(lmerTest)
library(car)
library(MuMIn)

# Syllable Task
lme_aq_bytrial <- lmer(PC_Low ~ AQ + (1|Talker/Syll) + (1|Code) + (0+Run|Code), data=clear_data_98)

# Sentence Task
for (i in 1:70){
  
  for (j in 1:length(movie_data_70$Sub)){
    
    if (movie_data_70[j,]$Sub == i) {
      movie_data_70[j,]$AQ <- AQ_mov[i,]$score
    }
  }
}

lme_aq_bytrial_mov <- lmer(P_Lower ~ AQ + (1|Movie) + (1|Sub) + (0+Run|Sub), data=movie_data_70)
