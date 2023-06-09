#install.packages('readxl')
#install.packages('gridExtra')
#install.pachakges('ggplot2')
#install.packages('RcppRoll')
require(readxl)
require(RcppRoll) # for efficient moving average

vga <- NULL
vo <- NULL
reco <- NULL

# Read meta-file
meta <- read.csv(file = "Meta_data/meta_trials.csv")
colnames(meta)[1] <- 'trial_ID' # Sometimes R fails to read the first column name
meta$start_playback <- as.POSIXct(strptime(meta$start_playback, "%d-%m-%Y %H:%M:%S"))
meta <- subset(meta, meta$trial_ID != "8") # Trial 8 failed
meta$interval<-as.numeric(as.character(meta$interval))

# Start times last (novel) exposure
last_e <- data.frame('Treatment' = c(rep("exposure", 6), "control"),
                     'Interval' = c(5, 10, 20, 40, 80, 160, NA),
                     'Start_last' = c(5400,5400,5400,5425,5410,5460,5400)) # in seconds


# Loop to process all trials
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
  df$Open <- sqrt(1/(0.000005*df$Open+0.0005)) # = Calibration
  
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
  df$Sound <- NA
  for(i in unique(df$Mussel)){
    # Mark all but last exposure
    if(meta$treatment[meta$trial_ID == trial] == 'exposure'){
      for(j in seq((25*60), 
                   (last_e$Start_last[last_e$Interval == meta$interval[meta$trial_ID == trial] & last_e$Treatment == "exposure"] - 1),
                   (meta$interval[meta$trial_ID == trial]+5))){
        df$Sound[df$Mussel == i & df$TrialTime >= j][1] <- as.character(meta$first_tone[meta$trial_ID == trial])
      }
    }
    # Mark last exposure
    if(meta$treatment[meta$trial_ID == trial] == 'exposure'){
      if(meta$first_tone[meta$trial_ID == trial] == 200){
        last_freq <- 350
      }else{
        last_freq <- 200
      }
      df$Sound[df$Mussel == i & df$TrialTime > last_e$Start_last[last_e$Interval == meta$interval[meta$trial_ID == trial]][1]][1] <- as.character(last_freq)
    }else{ # control
      df$Sound[df$Mussel == i & df$TrialTime > last_e$Start_last[last_e$Treatment == 'control']][1] <- as.character(meta$first_tone[meta$trial_ID == trial])
    }
  }
  df$Sound[!is.na(df$Sound)] <- paste(df$Sound[!is.na(df$Sound)], 'Hz', sep = " ")
  df$Sound <- as.factor(df$Sound)
  df$SoundYN[!is.na(df$Sound)] <- 0
  
  
  ### Save 10 min per mussel to check for independence of the data
  require(xts)
  df.t.z <- NULL
  df.t <- df
  for(i in unique(df.t$Mussel)){
    if(max(df.t$Open[df.t$Mussel == i], na.rm = T) > 0.00){
    max.open <- max(df.t$Open[df.t$Mussel == i], na.rm = T)
    if(max.open < .5){
      df.t$Open[df.t$Mussel == i] <- NA
    }
    # switch to proportions
    df.t$Open[df.t$Mussel == i] <- df.t$Open[df.t$Mussel == i] / max.open
    
    # To regular time interval
    tenmin <- meta$start_playback[meta$trial_ID == trial] + seq((15*60), (25*60)-2, by = 2)
    observation <- xts(df.t$Open[df.t$Mussel == i], order.by = df.t$DateTime[df.t$Mussel == i])
    df.t.i <- as.data.frame(na.locf(merge(xts(,tenmin),observation))[tenmin])
    df.t.i <- data.frame('Mussel' = i, 'Open' = df.t.i$observation, 'TrialTime' = tenmin - meta$start_playback[meta$trial_ID == trial],
                         'trial' = trial)
    # merge individuals
    if(exists('df.t.z') && is.data.frame(get('df.t.z'))){
      df.t.z <- rbind(df.t.z, df.t.i)
    }else{
      df.t.z <- df.t.i
    }
    }
  }
  # merge different trials
  if(exists('df10min') && is.data.frame(get('df10min'))){
    df10min <- rbind(df10min, df.t.z)
  }else{
    df10min <- df.t.z
  }
  
  require(ggplot2)
  p.valve <- ggplot(data = df, aes(x = as.numeric(TrialTime)/60, y = Open))+
    geom_path()+
    geom_point(aes(x = as.numeric(TrialTime)/60, y = SoundYN, color = Sound))+
    facet_grid(Mussel ~ .)+
    ggtitle(paste("Trial", trial, "Valve gape", sep = " "))+
    ylab("Open (mm)")+xlab("Time (min)")+
    theme_bw()
  p.valve
  
  ######
  ## Absorption
  ######
  df2 <- read.csv("Raw_data/Absorption.csv")
  colnames(df2)[1] <- "Trial"
  
  # Select current trial only
  df2 <- subset(df2, df2$Trial == trial)
  
  # Make into wide format
  require(tidyr)
  df2 <- df2 %>% gather(Time, Absorption, -c(Trial, ID, Mussel))
  df2$Time <- sub("T_", "", df2$Time)
  df2$Time <- as.numeric(df2$Time)
  
  df2 <- subset(df2, df2$Absorption < 0.3) # = extreme outliers
  
  require(ggplot2)
  p.absorp <- ggplot(data = df2, aes(x = Time, y = Absorption, color = Mussel))+
    geom_point()+
    stat_summary(fun.y = mean, geom = "point", shape = 4, color = "red")+
    stat_summary(fun.y = mean, geom = "line", color = "red")+
    facet_grid(ID ~ .)+
    ggtitle("Filtration")+
    xlab("Time (min)")+
    theme_bw()
  p.absorp
  
  require(gridExtra)
  #print("Save plot")
  #overview <- grid.arrange(p.valve, p.absorp, nrow = 1, widths = c(0.7, 0.3))
  #ggsave(paste('Trial_', trial, '_new.png',sep = ""), overview, png(), units = "in", width = 8, height = 6.5)
  
  
  #######
  # Save valve gapes (0-60 & 60-95) to new df
  #######
  print("Get and save data")
  for(i in unique(df$Mussel)){
    
    # Get mean open
    m.open60 <- mean(df$Open[df$Mussel == i & df$TrialTime >= 0 & df$TrialTime < (60*60)], na.rm = T)
    m.open95 <- mean(df$Open[df$Mussel == i & df$TrialTime >= (60*60) & df$TrialTime < (95*60)], na.rm = T)
    
    # Get difference in absorption
    abs0 <- mean(df2$Absorption[df2$ID == i & df2$Time == 0], na.rm = T)
    abs60 <- mean(df2$Absorption[df2$ID == i & df2$Time == 60], na.rm = T) 
    abs95 <- mean(df2$Absorption[df2$ID == i & df2$Time == 95], na.rm = T)
    
    # maximum opening
    max.open <- max(df$Open[df$Mussel == i], na.rm = T)
    
    # Get sample rate
    s.rate <- length(df$Open[df$Mussel == i]) / ((as.numeric(range(df$DateTime[df$Mussel == i])[2]) - as.numeric(range(df$DateTime[df$Mussel == i])[1])) / 60)
    
    # Was there a mussel present? [it can be that the sensor failed, but that the mussel was present]
    if(is.na(meta[meta$file_name == file_name,(14+as.numeric(i))])){
      mussel.pres <- F
    }else{
      mussel.pres <- T
    }
    
    # Make temporary row for dataframe
    vga.temp <- data.frame('Trial' = trial, 'Mussel' = i, 'Treatment' = meta$treatment[meta$trial_ID == trial], 'Interval' = meta$interval[meta$trial_ID == trial], 'Start_tone' = meta$first_tone[meta$trial_ID == trial], 'Mussel_size' = as.numeric(meta[meta$trial_ID == trial, 14+as.numeric(i)]), 'Mean_gape_t60' = m.open60, 'Mean_gape_t95' = m.open95, 'Mean_abs_t0' = abs0, 'Mean_abs_t60' = abs60, 'Mean_abs_t95' = abs95, 'Max_open' = max.open, 'Sample_rate' = s.rate, 'Mussel_presence' = mussel.pres)
  
    # Bind row to actual df
    if(exists('vga') && is.data.frame(get('vga'))){
      vga <- rbind(vga, vga.temp)
    }else{
      vga <- vga.temp
    }
    
    #######
    # Do mussels close their valves in response to sound, new loop
    #######
    
    # first sound
    m.fb <- mean(df$Open[df$Mussel == i & df$TrialTime >= ((25*60)-20) & df$TrialTime < (25*60)], na.rm = T)
    m.fa <- mean(df$Open[df$Mussel == i & df$TrialTime >= (25*60) & df$TrialTime < ((25*60)+20)], na.rm = T)
    # second last sound
    if(meta$treatment[meta$trial_ID == trial] == "exposure"){
      interv <- meta$interval[meta$trial_ID == trial]
    }else if(meta$treatment[meta$trial_ID == trial] == "controle"){
      interv <- 20
    }
    m.slb <- mean(df$Open[df$Mussel == i & df$TrialTime >= (last_e$Start_last[last_e$Interval == interv][1]-interv-5-(20)) & df$TrialTime < (last_e$Start_last[last_e$Interval == interv][1]-interv-5)], na.rm = T)
    m.sla <- mean(df$Open[df$Mussel == i & df$TrialTime >= (last_e$Start_last[last_e$Interval == interv][1]-interv-5) & df$TrialTime < (last_e$Start_last[last_e$Interval == interv][1]-interv-5+(20))], na.rm = T)
    # last sound
    m.lb <- mean(df$Open[df$Mussel == i & df$TrialTime >= last_e$Start_last[last_e$Interval == interv][1]-(20) & df$TrialTime < last_e$Start_last[last_e$Interval == interv][1]], na.rm = T)
    m.la <- mean(df$Open[df$Mussel == i & df$TrialTime >= last_e$Start_last[last_e$Interval == interv][1] & df$TrialTime < last_e$Start_last[last_e$Interval == interv][1]+(20)], na.rm = T)
    
    # Make temporary row for dataframe
    vo.temp <- data.frame('Trial' = trial, 'Mussel' = i, 'Treatment' = meta$treatment[meta$trial_ID == trial], 'Interval' = meta$interval[meta$trial_ID == trial], 'Start_tone' = meta$first_tone[meta$trial_ID == trial], 'Mussel_size' = as.numeric(meta[meta$trial_ID == trial, 14+as.numeric(i)]), 'Mean_gape_bef_fst' = m.fb, 'Mean_gape_aft_fst' = m.fa, 'Mean_gape_bef_seclst' = m.slb, 'Mean_gape_aft_seclst' = m.sla, 'Mean_gape_bef_lst' = m.lb, 'Mean_gape_aft_lst' = m.la, 'Max_open' = max.open, 'Mussel_presence' = mussel.pres)
    
    # Bind row to actual df
    if(exists('vo') && is.data.frame(get('vo'))){
      vo <- rbind(vo, vo.temp)
    }else{
      vo <- vo.temp
    }
    
    #######
    # What is the recovery time to the sound?
    #######
    
    # make df with single individual
    df.sub <- subset(df, df$Mussel == i)
    
    bin <- 200 # bin size = sec/or no. of values
    BAv <- mean(df$Open[df$Mussel == i & df$TrialTime >= ((25*60)-bin) & df$TrialTime < (25*60)], na.rm = T)
    
    
    # Moving average
    TT <- roll_mean(df.sub$TrialTime[df.sub$TrialTime >= (25*60)], bin, fill = NA)
    MAv <- roll_mean(df.sub$Open[df.sub$TrialTime >= (25*60)], bin, fill = NA)
    MA <- data.frame('TrialTime' = TT, 'Valve_opening' = MAv)
    
    # Is there data of this individual?
    if(is.na(BAv)){
      Recov.min <- NA
    }else{
      # Does this individual recover at all?
      if(max(MA$Valve_opening, na.rm = T) > BAv){
        # At what point is the opening higher than during baseline?
        Recov.min <- ((MA$TrialTime[MA$Valve_opening > BAv & !is.na(MA$Valve_opening)][1]-min(MA$TrialTime, na.rm = T))/60)
      }else{
        Recov.min <- 70
      }
    }
      
    # Make temporary row for dataframe
    reco.temp <- data.frame('Trial' = trial, 'Mussel' = i, 'Treatment' = meta$treatment[meta$trial_ID == trial], 'Interval' = meta$interval[meta$trial_ID == trial], 'Start_tone' = meta$first_tone[meta$trial_ID == trial], 'Mussel_size' = as.numeric(meta[meta$trial_ID == trial, 14+as.numeric(i)]), 'Recovery_min' = Recov.min, 'Max_open' = max.open, 'Mean_gape_bef_fst' = m.fb, 'Mussel_presence' = mussel.pres)
    
    # Bind row to actual df
    if(exists('reco') && is.data.frame(get('reco'))){
      reco <- rbind(reco, reco.temp)
    }else{
      reco <- reco.temp
    }
    
  }
} # end loop all trials

write.csv(vga, file = "Processed_data/Valve_gape_Absorption.csv")
write.csv(vo, file = "Processed_data/Valve_opening_sound.csv")
write.csv(reco, file = "Processed_data/Valve_opening_recovery.csv")
write.csv(df10min, file = "Processed_data/Valve_gape_10min.csv")