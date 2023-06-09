#Calculating encounter rate from field and lab densities
library(lubridate)

data<-read.csv("...Field Demographics\\\\Host Density Foitzik MyData Combined.csv")

#Field: 
 #My data - nests per m2
  fdensity_md<-mean(data$hostspm2[data$DataSource=="MyData"])
  #0.2122066

  SD_fdensity_md<-sd(data$hostspm2[data$DataSource=="MyData"])
  
  length(data$hostspm2[data$DataSource=="MyData"])
  
  #95% CI intervals for my data - DENSITY
  UPCI_md<-fdensity_md+(1.96*(SD_fdensity_md)/sqrt(9))
  LOCI_md<-fdensity_md-(1.96*(SD_fdensity_md)/sqrt(9))
  
  #MAX/MIN 
  max_density_md<-max(data$hostspm2[data$DataSource=="MyData"])
  min_density_md<-min(data$hostspm2[data$DataSource=="MyData"])
  

 #Foitzik data
  fdensity_fd<-mean(data$hostspm2[data$DataSource=="Foitzik"])
  #0.4911258
  
  SD_fdensity_fd<-sd(data$hostspm2[data$DataSource=="Foitzik"])
  
  #MAX/MIN
  max_density_fd<-max(data$hostspm2[data$DataSource=="Foitzik"])
  min_density_fd<-min(data$hostspm2[data$DataSource=="Foitzik"])
  

#Lab
  #Density:
  ldensity<-17.78

  #Encounter Rate
  lencounter<-7476
  SD_lencounter<-6514.137867
  n=13
  
  UPCI_LER<-lencounter+1.96*( SD_lencounter/sqrt(n))
  LOCI_LER<-lencounter-1.96*( SD_lencounter/sqrt(n))
  
  
#Encounter Rate - MY DATA
  ERF_md<-lencounter*ldensity/fdensity_md
  
  seconds_to_period(ERF_md)
  #"7d 5H 59M 46.1992863359628S"
  
  #LAB ENCOUTNER 95% CI
    #UPCI
      UP_ERF_md1<-UPCI_LER*ldensity/fdensity_md
      seconds_to_period(UP_ERF_md1)
      #"10d 16H 24M 43.9339534103638S"
      
   #LOCI
      LO_ERF_md1<-LOCI_LER*ldensity/fdensity_md
      seconds_to_period(LO_ERF_md1)
      #"3d 19H 34M 48.4659311618889S"
  
  #DENSITY 95% CI
  #UPCI
    UP_ERF_md<-lencounter*ldensity/UPCI_md
    seconds_to_period(UP_ERF_md)
    #"5d 14H 32M 52.0435550636612S"
  #LOCI
    LO_ERF_md<-lencounter*ldensity/LOCI_md
    seconds_to_period(LO_ERF_md)
    #"10d 6H 10M 18.7059005649062S"
  
  #MAX
  MAX_ERF_md<-lencounter*ldensity/max_density_md
  seconds_to_period(MAX_ERF_md)
  #"4d 8H 23M 51.720280227717S"
  #MIN
  MIN_ERF_md<-lencounter*ldensity/min_density_md
  seconds_to_period(MIN_ERF_md)
  #"14d 11H 59M 32.404476222815S"
  
#Encounter Rate - Foitzik Data
  ERF_fd<-lencounter*ldensity/fdensity_fd
  
  seconds_to_period(ERF_fd)
  #"3d 3H 10M 50.1902764640399S"
  
  #MAX
  MAX_ERF_fd<-lencounter*ldensity/max_density_fd
  seconds_to_period(MAX_ERF_fd)
  #"1d 16H 37M 10.2310231023002S"
  #MIN
  MIN_ERF_fd<-lencounter*ldensity/min_density_fd
  seconds_to_period(MIN_ERF_fd)
  #"13d 20H 18M 29.5211963094771S"
  

#Encounter Rate - Combined
  fdensity_mean<-mean(data$hostspm2)
  SD_fdensity<-sd(data$hostspm2)
  ERF_all<-lencounter*ldensity/fdensity_mean

  
  seconds_to_period(ERF_all)
  #"4d 4H 59M 23.645689575118S"
  
  #MAX/MIN
  max_density_all<-max(data$hostspm2)
  min_density_all<-min(data$hostspm2)
  
  #MAX
  MAX_ERF_all<-lencounter*ldensity/max_density_all
  seconds_to_period(MAX_ERF_all)
  # 1d 16H 37M 10.2310231023002S
  #MIN
  MIN_ERF_all<-lencounter*ldensity/min_density_all
  seconds_to_period(MIN_ERF_all)
  #14d 11H 59M 32.404476222815S
  
  
  
  #LAB ENCOUTNER 95% CI
  #UPCI
  UP_ERF_all1<-UPCI_LER*ldensity/fdensity_mean
  seconds_to_period(UP_ERF_all1)
  #6d 4H 49M 31.3186154137366S"
  
  #LOCI
  LO_ERF_all1<-LOCI_LER*ldensity/fdensity_mean
  seconds_to_period(LO_ERF_all1)
  #"2d 5H 9M 15.9727637364995S"
  
  #Confidence Interval = mean +-1.96*std/sqrt n