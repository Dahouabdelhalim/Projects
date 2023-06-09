df <- read.csv(file = 'Processed_data/Valve_gape_10min.csv')[2:5]
df$Un_id <- paste(df$trial, df$Mussel, sep = '_')
# remove NA's
df <- na.omit(df)

df.ccf <- NULL
for(i in unique(df$Un_id)){
  print(i)
  trial <- df$trial[df$Un_id == i][1]
  mussel <- df$Mussel[df$Un_id == i][1]
  
  ## cross-correlate with other mussels from the same trial
  others <- unique(df$Mussel[df$trial == trial])
  others <- others[!others %in% mussel] # remove current mussel
  for(j in others){
    # did the combination already occur in other order?
    if(length(df.ccf$acf[df.ccf$Combination == paste(trial, '_', j, '-', trial, '_', mussel, sep = "")]) > 0){
      skip = T
    }else{
      skip = F
      df.ccf.t <- ccf(df$Open[df$Mussel == mussel & df$trial == trial], df$Open[df$Mussel == j & df$trial == trial], lag = 60, plot = F)
      df.ccf.t <- data.frame('acf' = df.ccf.t$acf, 'lag' = df.ccf.t$lag, 'Trial' = 'Same',
                             'Combination' = paste(i, '-', trial, '_', j, sep = ""))
      # Bind row to actual df
      if(exists('df.ccf') && is.data.frame(get('df.ccf'))){
        df.ccf <- rbind(df.ccf, df.ccf.t)
      }else{
        df.ccf <- df.ccf.t
      }
    }
  }
  
  ## With others mussels from a different trial
  othertrials <- unique(df$trial)
  othertrials <- othertrials[!othertrials %in% trial]
  # as many combinations as above
  for(j in 1:2){ # length(others)
    trial.r <- sample(othertrials, 1)
    mussel.r <- sample(unique(df$Mussel[df$trial == trial.r]), 1)
    df.ccf.t <- ccf(df$Open[df$Mussel == mussel & df$trial == trial], df$Open[df$Mussel == mussel.r & df$trial == trial.r], lag = 60, plot = F)
    df.ccf.t <- data.frame('acf' = df.ccf.t$acf, 'lag' = df.ccf.t$lag, 'Trial' = 'Different',
                           'Combination' = paste(i, '-', trial.r, '_', mussel.r, sep = ""))
    # Bind row to actual df
    if(exists('df.ccf') && is.data.frame(get('df.ccf'))){
      df.ccf <- rbind(df.ccf, df.ccf.t)
    }else{
      df.ccf <- df.ccf.t
    }
  }
}

# Get mean, median, 1st and 3th quantile
ccf.total <- data.frame('Lag' = rep(df.ccf.t$lag, 2), 'Trial' = c(rep('Same', length(df.ccf.t$lag)), rep('Different', length(df.ccf.t$lag))),
                        'Mean' = NA, 'Median' = NA, 'First_q' = NA, 'Third_q' = NA)
for(i in 1:nrow(ccf.total)){
  lag <- ccf.total$Lag[i]
  trial <- ccf.total$Trial[i]
  ccf.total$Mean[i] <- mean(abs(df.ccf$acf[df.ccf$lag == lag & df.ccf$Trial == trial]))
  ccf.total$Median[i] <- median(abs(df.ccf$acf[df.ccf$lag == lag & df.ccf$Trial == trial]))
  ccf.total$First_q[i] <- quantile(abs(df.ccf$acf[df.ccf$lag == lag & df.ccf$Trial == trial]))[2]
  ccf.total$Third_q[i] <- quantile(abs(df.ccf$acf[df.ccf$lag == lag & df.ccf$Trial == trial]))[4]
}

p1 <- ggplot(data = ccf.total, aes(x = Lag*2, y = Median, fill = Trial))+
  geom_line(aes(color = Trial), size = 1.1)+
  geom_ribbon(aes(ymin = First_q, ymax = Third_q), alpha = 0.2)+
  ylab('Pearson correlation coefficient')+
  xlab('Lag (s)')+
  ggtitle('Cross-correlations of valve gape behaviour')+
  ylim(0,.8)+xlim(-100, 100)+
  scale_color_manual(name="Individuals of", labels=c("another trial", "the same trial"), values = c('#f86def', '#6df876'))+
  scale_fill_manual(name="Individuals of", labels=c("another trial", "the same trial"), values = c('#f86def', '#6df876'))+
  theme_bw()
p1
ggsave('Figures/Fig_s3_HQ.jpg', p1, png(), units = "in", width = 8, height = 4.5, dpi = 600)


## numbers:
# Unique individuals
length(unique(df$Un_id))
# Combinaties per individual for 'same'
length(unique(df.ccf$Combination[df.ccf$Trial == 'Same']))
length(unique(df.ccf$Combination[df.ccf$Trial == 'Different']))
