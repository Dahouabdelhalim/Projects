#install.packages('devtools')
#require(devtools)
#install_git("https://gitlab.com/RTbecard/ezPSD.git", build_vignettes = F, force = T)
require(ezPSD)
#install.packages('ggplot2')
require(ggplot2)

# Read meta
meta <- read.csv2(file = "Meta_data/meta_recordings.csv")

# Load recording
files <- list.files(path = "Sound_recordings/", pattern = ".wav$")

# Recorder calibration
rec.calib <- data.frame('Sensitivity' = c(1, 4, 6, 10), 'Calib_constant' = c(25.36, 15.28, 5, -5.67))

PSD.df <- NULL
PSD.plot <- NULL
for(i in 1:length(meta$name)){
  print(paste("Read:", meta$name[i], "- rec sensitivity:", meta$sensitivity[i], sep = " "))
  
  # Read recording
  wav <- read.wav(file = paste('Sound_recordings/', as.character(meta$name[i]), sep = ""))
  
  # Enter calibration constants
  calib.hydrophone <- -158 # AS-1: -208; AS-1 + PA4: -158
  calib.recorder <- rec.calib$Calib_constant[rec.calib$Sensitivity == meta$sensitivity[i]]
  
  # Sample rate
  fs <- 44100

  ## Calibrate signal and apply bandpass filter
  data.calib <- ezPSD::calibrateSignal.Time(data = wav[[1]], calib = (calib.hydrophone - calib.recorder))

  # Where does the sound begin?
  st.sound <- min(which(data.calib[(0.2*fs):length(data.calib)] > 1*10^6.5))+0.2*fs
  print(paste('Start time = ', st.sound/fs, sep = ""))
  data.calib.200 <- data.calib[(st.sound+(fs*2.5)):(st.sound+(fs*5))] # get 5.5 sec of 200 hz sound
  data.calib.sil <- data.calib[(st.sound+(fs*12.5)):(st.sound+(fs*15))] # get 5.5 sec of silence
  data.calib.350 <- data.calib[(st.sound+(fs*22.5)):(st.sound+(fs*25))] # get 5.5 sec of 350 hz sound
  
  # Make PSD
  welch.200 <- ezWelch(data.calib.200, wl = (1024*6), olap = 0.5, fs = fs, windowType = 'Hann')
  welch.sil <- ezWelch(data.calib.sil, wl = (1024*6), olap = 0.5, fs = fs, windowType = 'Hann')
  welch.350 <- ezWelch(data.calib.350, wl = (1024*6), olap = 0.5, fs = fs, windowType = 'Hann')
  message('# of overlapping Window Segments: ',nrow(welch.350$WindowSegments))
  
  PSD.df.temp.200 <- data.frame('Frequency' = welch.200$Frequency, 'Power' = welch.200$PSD, 'File' = meta$name[i],
                            'Location' = meta$location[meta$name == meta$name[i]],
                            'Playback' = '200 Hz')
  PSD.df.temp.sil <- data.frame('Frequency' = welch.sil$Frequency, 'Power' = welch.sil$PSD, 'File' = meta$name[i],
                                'Location' = meta$location[meta$name == meta$name[i]],
                                'Playback' = 'Silence')
  PSD.df.temp.350 <- data.frame('Frequency' = welch.350$Frequency, 'Power' = welch.350$PSD, 'File' = meta$name[i],
                                'Location' = meta$location[meta$name == meta$name[i]],
                                'Playback' = '350 Hz')
  PSD.df.temp <- rbind(PSD.df.temp.200, PSD.df.temp.sil, PSD.df.temp.350)
  
  # Save
  if(exists('PSD.df') && is.data.frame(get('PSD.df'))){
    PSD.df <- rbind(PSD.df, PSD.df.temp)
  }else{
    PSD.df <- PSD.df.temp
  }
  
  # Determine sound level
  ref <- 1 # reference value of 1uPa for underwater measurements

  # Sound level based on PSD (unbandpassed)
  print(paste('200 Hz', 10*log10(sum(welch.200$Power[welch.200$Frequency >= 100 & welch.200$Frequency <= 600])/ref)), sep = " ")
  print(paste('350 Hz', 10*log10(sum(welch.350$Power[welch.350$Frequency >= 100 & welch.350$Frequency <= 600])/ref)), sep = " ")
  print(paste('Silence', 10*log10(sum(welch.sil$Power[welch.sil$Frequency >= 100 & welch.sil$Frequency <= 600])/ref)), sep = " ")
  print('--------------- next file ---------------')
}
PSD.df$Location <- as.factor(PSD.df$Location)
PSD.df$File_PB <- paste(PSD.df$File, PSD.df$Playback, sep = "_")

PSD.df$Playback <- as.factor(PSD.df$Playback)
PSD.df$Playback <- factor(PSD.df$Playback, levels = c("200 Hz", "350 Hz", "Silence"))

PSD.plot <- ggplot(data = PSD.df, aes(x = Frequency, y = 10*log10(Power/20), group = File_PB))+
  geom_line(aes(color = Playback, linetype = Location))+
  xlab("Frequency (Hz)")+ylab('SPL (dB re 1 \\u03bcPa/Hz)')+
  scale_x_continuous(breaks=seq(0,900,100), limits = c(0,900))+
  scale_y_continuous(limits = c(28,122))+
  ggtitle('Power spectral density (PSD)')+
  guides(linetype = FALSE)+
  theme_bw()+theme(legend.position="bottom")
PSD.plot
ggsave('Figures/Fig_2_HQ.jpg', PSD.plot, png(), units = "in", width = 4.3, height = 4.3, dpi = 600)
