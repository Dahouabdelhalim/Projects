#install.packages('readxl')
#install.packages('gridExtra')
#install.pachakges('ggplot2')
#install.packages('RcppRoll')
require(readxl)
require(RcppRoll) # for efficient moving average

# Read meta-file
meta <- read.csv(file = "Meta_data/meta_trials.csv")
colnames(meta)[1] <- 'trial_ID'
meta$start_playback <- as.POSIXct(strptime(meta$start_playback, "%d-%m-%Y %H:%M:%S"))
meta <- subset(meta, meta$trial_ID != 8)
meta$interval<-as.numeric(as.character(meta$interval))

# Start times last (novel) exposure
last_e <- data.frame('Treatment' = c(rep("exposure", 6), "control"),
                     'Interval' = c(5, 10, 20, 40, 80, 160, NA),
                     'Start_last' = c(5400,5400,5400,5425,5410,5460,5400)) # in seconds

df.all <- NULL

# Select trial
for(trial in meta$trial_ID){ #
  print(paste("Trial: ", trial, " (out of ", nrow(meta), ")", sep = ""))
  file_name <- as.character(meta$file_name[meta$trial_ID == trial])
  
  # Determine which rows to read
  rng <- read.csv2(file = paste('Raw_data/', file_name, sep = ""), header = F)
  rng$nchar <- nchar(as.character(rng$V1))
  start <- min(which(rng$nchar == median(rng$nchar)))
  end <- max(which(rng$nchar == median(rng$nchar)))
  
  # READ DATA
  df <- read.table(file = paste('Raw_data/', file_name, sep = ""), sep = "", header = F, 
                   skip = start, nrows = end-start+1)
  rng <- NULL
  
  # Recognize date and time as such
  df$V1 <- substring(df$V1, 2)
  df$V2 <- substring(df$V2, 1, nchar(as.character(df$V2))-1)
  df$DateTime <- as.POSIXct(paste0(df$V1, df$V2, sep = " "), tz = "CET")
  
  # Print time range
  range(df$DateTime)
  
  # Rename mussel columns
  colnames(df)[3:10] <- c("Mussel_1", "Mussel_2", "Mussel_3", "Mussel_4", "Mussel_5", "Mussel_6",
                          "Mussel_7", "Mussel_8")
  
  # From wide to long
  df <- df[,c(3:10,12)]
  library(tidyr)
  df <- df %>% gather(Mussel, Open, -c(DateTime))
  df$Mussel <- as.factor(sub("Mussel_", "", df$Mussel))
  
  # Convert data to actual distance
  df$Open <- sqrt(1/(0.000005*df$Open+0.0005))
  
  # Absolute difference in distance
  for(i in unique(df$Mussel)){ # Check if mussel/sensor was valid
    if(meta[meta$file_name == file_name,(6+as.numeric(i))] == F){ # check if mussel was valid (from the meta file)
      df$Open[df$Mussel == i] <- NA
    }else{
      df$Open[df$Mussel == i] <- df$Open[df$Mussel == i]-min(df$Open[df$Mussel == i])
    }
  }
  
  # Add timing of the sound
  df$TrialTime <- (df$DateTime - meta$start_playback[meta$trial_ID == trial]) # minutes

  #remove mussels that open less than 1mm or if the mussel was not open for 25% of max. at 23 min
  for(i in unique(df$Mussel)){
    if(!is.na(mean(df$Open[df$Mussel==i]))){
      if(max(df$Open[df$Mussel==i])<0.5){
        df$Open[df$Mussel==i]<-NA 
      }
    }
  }
  
  # Add unique id for each mussel 
  df$Trial_mussel <- paste(trial, '_', df$Mussel, sep = "")
  df$Treatment <- meta$treatment[meta$trial_ID == trial]
  
  # Add meta-info
  df$Interval <- meta$interval[meta$trial_ID == trial]
  
  # Bind row to final df
  if(exists('df.all') && is.data.frame(get('df.all'))){
    df.all <- rbind(df.all, df)
  }else{
    df.all <- df
  }
} # end loop all trials

df.all$Treatment <- as.factor(df.all$Treatment)

df.all$TrialTime2 <- round(df.all$TrialTime/20)

library(dplyr)
df.all.ds <- df.all %>%
  group_by(round(TrialTime2), Trial_mussel) %>%
  summarise(mean = median(Open, na.rm = T), Treatment = unique(Treatment), Interval = unique(Interval))

group_means <- NULL
group_means <- df.all %>%
  group_by(round(TrialTime2), Treatment) %>%
  summarise(mean = mean(Open, na.rm = T),
            quantile_25 = quantile(Open, na.rm = T)[[2]],
            quantile_75 = quantile(Open, na.rm = T)[[4]])

colnames(group_means)[2] <- 'Treatment2'

require(ggplot2)
p2 <- ggplot(data = group_means, aes(x = `round(TrialTime2)`/3, y = mean, color = `Treatment2`))+
  geom_ribbon(aes(ymin = quantile_25, ymax = quantile_75, fill = Treatment2), alpha = .25, color = NA, show.legend = F)+
  geom_path(size = 1.1)+
  #stat_summary(geom="ribbon", fun.ymin="quantile_25", fun.ymax="quantile_75", aes(fill=Treatment2), alpha=0.3)+
  ylab("Open (mm)")+xlab("Time (min)")+
  scale_x_continuous(breaks=seq(0,90,10), limits = c(0,95))+
  scale_y_continuous(limits = c(0,4))+
  ggtitle('Mean valve gape over time')+
  scale_color_discrete(name = 'Treatment', labels = c(
    paste("Control (n =", length(unique(df.all$Trial_mussel[df.all$Treatment == 'controle' & !is.na(df.all$Open)])), ")", sep = ""), 
    paste("Exposure (n =", length(unique(df.all$Trial_mussel[df.all$Treatment == 'exposure' & !is.na(df.all$Open)])), ")", sep = "")))+
  theme_bw()
p2
ggsave('Figures/Fig_s1_HQ.jpg', p2, png(), units = "in", width = 8, height = 5, dpi = 600)
