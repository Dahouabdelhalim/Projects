### Code for processing of STAP raw files
#
# Author: Sebastian Duesing, edited by Christian Pilz, contact: pilz@tropos.de
#
# Information about the STAP are provided by the manufacturer Brechtel Inc. 
# Performance evaluation of the STAP are provided by Pilz et al. 2022 (https://doi.org/10.5194/amt-15-6889-2022)
# 
# STP correction to p = 1013.25 hPa and T = 273.15 K
# Azumi correction accounts for 25 % absorption enhancement by filter material
#
# Code requires manual plausibility check of data!!! 
###


recalcSTAPtoFile <-
  function(filename,
           stpcorrection = T,
           blackcorrection = T,
           averagingPeriod = 1,
           oldFile = F,
           azumi = F, 
           tau.use = T, 
           K2.use = T) {
    
    
    for (package in c('lubridate', 'tidyverse')) { #checks for packages installed
      if (!require(package, character.only=T, quietly=T)) {
        install.packages(package)
        library(package, character.only=T)
      }
    }
    
    
    absorptionTimePeriod <- function(period = 60,
                                     smp,
                                     ref,
                                     blk.smp,
                                     blk.ref,
                                     smp.p,
                                     smp.t,
                                     smp.flw,
                                     I0.ref,
                                     I0.smp,
                                     I0.blk.smp,
                                     I0.blk.ref,
                                     withBlack = TRUE,
                                     stp = TRUE,
                                     tau.use = TRUE,
                                     spot_area = 1.79451E-05,
                                     azumi = T,
                                     K2.use = T) {
      if (withBlack == T) {
        I0.smp <- I0.smp - I0.blk.smp
        I0.ref <- I0.ref - I0.blk.ref
        ref <- ref - blk.ref
        smp <- smp - blk.smp
      }
      i <- 1
      inv_mm <- rep(NA, period)
      while (i <= (length(smp) - period)) {

        smp_old <- mean(smp[(i):(i + (period - 1))], na.rm = T)
        ref_old <- mean(ref[(i):(i + (period - 1))], na.rm = T)
        
        smp_new <- mean(smp[(i + period):(i + (period * 2 - 1))], na.rm = T)
        ref_new <- mean(ref[(i + period):(i + (period * 2 - 1))], na.rm = T)
        
        p <- mean(smp.p[(i + period):(i + (period * 2 - 1))], na.rm = T)
        t <- mean(smp.t[(i + period):(i + (period * 2 - 1))], na.rm = T)
        if (stp == TRUE) {
          stp.corr <- p / 1013.25 * 273.15 / (t + 273.15)
        } else{
          stp.corr <- 1
        }
        flow_vol <-
          mean(smp.flw[(i + period):(i + (period * 2 - 1))], na.rm = T) * 0.001 /
          60 * period * stp.corr
        if (tau.use == T) {
          tau <- (smp_new / ref_new) / (I0.smp / I0.ref)
          
        } else{
          tau <- 1
        }
        if(azumi == T){
          azumi_corr = 1.25
        }else{
          azumi_corr = 1
        }
        if(K2.use == T){
          K2 = 1.22
        }else{
          K2 = 0.85
        }
        inv_mm <-
          c(
            inv_mm,
            1/azumi_corr * 1e6 * 0.85/K2 *1/ (1.0796 * tau + 0.71) * (spot_area / flow_vol) * log((smp_old / ref_old) / (smp_new / ref_new), base = exp(1))
            
          )
        
        i <- i + 1
      }
      if(period == 1){
        
      }else{
        inv_mm[(length(inv_mm)-(period-2)):length(inv_mm)]<-NA 
      }
      
      return(inv_mm)
    }
    
    tauTimePeriod <- function(period = 60,
                                     smp,
                                     ref,
                                     blk.smp,
                                     blk.ref,
                                     smp.p,
                                     smp.t,
                                     smp.flw,
                                     I0.ref,
                                     I0.smp,
                                     I0.blk.smp,
                                     I0.blk.ref,
                                     withBlack = TRUE,
                                     stp = TRUE,
                                     tau.use = TRUE,
                                     spot_area = 1.79451E-05,
                                     azumi = T,
                                     K2.use = T) {
      if (withBlack == T) {
        I0.smp <- I0.smp - I0.blk.smp
        I0.ref <- I0.ref - I0.blk.ref
        ref <- ref - blk.ref
        smp <- smp - blk.smp
      }
      i <- 1
      taulist <- rep(NA, period)
      while (i <= (length(smp) - period)) {
        
        smp_old <- mean(smp[(i):(i + (period - 1))], na.rm = T)
        ref_old <- mean(ref[(i):(i + (period - 1))], na.rm = T)
        
        smp_new <- mean(smp[(i + period):(i + (period * 2 - 1))], na.rm = T)
        ref_new <- mean(ref[(i + period):(i + (period * 2 - 1))], na.rm = T)
        
        p <- mean(smp.p[(i + period):(i + (period * 2 - 1))], na.rm = T)
        t <- mean(smp.t[(i + period):(i + (period * 2 - 1))], na.rm = T)
        if (stp == TRUE) {
          stp.corr <- p / 1013.25 * 273.15 / (t + 273.15)
        } else{
          stp.corr <- 1
        }
        flow_vol <-
          mean(smp.flw[(i + period):(i + (period * 2 - 1))], na.rm = T) * 0.001 /
          60 * period * stp.corr
        if (tau.use == T) {
          tau <- (smp_new / ref_new) / (I0.smp / I0.ref)
          
        } else{
          tau <- 1
        }
        if(azumi==T){
          azumi_corr = 1.25
        }else{
          azumi_corr = 1
        }
        if(K2.use == T){
          K2 = 1.22
        }else{
          K2 = 0.85
        }
        
        taulist <-
          c(
            taulist,tau)
        
        i <- i + 1
      }
      if(period == 1){
        
      }else{
        taulist[(length(taulist)-(period-2)):length(taulist)]<-NA 
      }
      
      return(taulist)
    }
    
    
    
    
    
    
    tmp <- getwd()
    setwd(dirname(filename))
    
    if (oldFile == T){
    ######## only for old files  
     initValues <-
       as.numeric(unlist(strsplit(scan(filename, what = "numeric", n = 50)[c(13:20)],"="))[c(2,4,6,8,10,12,14,16)]) #extract init values from file
     print(initValues)
     stap <-
      read.delim(filename,
                  skip = 25,
                  sep = "\\t",
                  header = T) # read STAP file
    
    
    } else {
      initValues <-
      as.numeric(scan(filename, what = "numeric", n = 50)[c(27, 29, 31, 33, 35, 37, 39, 41)]) #extract init values from file
      print(initValues)
    stap <-
      read.delim(filename,
                 skip = 29,
                 sep = "\\t",
                 header = T) # read STAP file
    
    
    }            
    
    if(nrow(stap)==0){
      print(paste0(basename(filename)," has no Entries! Skip."))
    }else if (nrow(stap)<(averagingPeriod*2)){
      print(paste0(basename(filename)," has not enough entries! Skip."))
    } else {
      stap$datum <-
        as.POSIXct(strptime(stap$X.YY.MM.DD, format = "%y/%m/%d", tz = "UTC"),
                   tz = "UTC")
      stap$time <-
        as.POSIXct(strptime(stap$HR.MN.SC, format = "%H:%M:%S", tz = "UTC"), tz =
                     "UTC")
      stap <- stap%>%mutate(`date` = datum + hour(time) * 3600 + minute(time) * 60 + second(time))
      
      stap <- stap %>%
        select(date, everything(),-X.YY.MM.DD, -HR.MN.SC, -datum ,-time)
      
      stap$red <-
        absorptionTimePeriod(
          period = averagingPeriod,
          smp = stap$red_smp,
          ref = stap$red_ref,
          blk.smp = stap$blk_smp,
          blk.ref = stap$blk_ref,
          smp.p = stap$smp_prs,
          smp.t = stap$smp_tmp,
          smp.flw = stap$smp_flw,
          I0.ref = initValues[2],
          I0.smp = initValues[2],
          I0.blk.smp = initValues[7],
          I0.blk.ref = initValues[8],
          withBlack = blackcorrection,
          stp = stpcorrection,
          tau.use = tau.use,
          azumi = azumi, 
          K2.use = K2.use
        )
      
      stap$green <-
        absorptionTimePeriod(
          period = averagingPeriod,
          smp = stap$grn_smp,
          ref = stap$grn_ref,
          blk.smp = stap$blk_smp,
          blk.ref = stap$blk_ref,
          smp.p = stap$smp_prs,
          smp.t = stap$smp_tmp,
          smp.flw = stap$smp_flw,
          I0.ref = initValues[4],
          I0.smp = initValues[3],
          I0.blk.smp = initValues[7],
          I0.blk.ref = initValues[8],
          withBlack = blackcorrection,
          stp = stpcorrection,
          tau.use = tau.use,
          azumi = azumi,
          K2.use = K2.use
        )
      
      stap$blue <-
        absorptionTimePeriod(
          period = averagingPeriod,
          smp = stap$blu_smp,
          ref = stap$blu_ref,
          blk.smp = stap$blk_smp,
          blk.ref = stap$blk_ref,
          smp.p = stap$smp_prs,
          smp.t = stap$smp_tmp,
          smp.flw = stap$smp_flw,
          I0.ref = initValues[6],
          I0.smp = initValues[5],
          I0.blk.smp = initValues[7],
          I0.blk.ref = initValues[8],
          withBlack = blackcorrection,
          stp = stpcorrection,
          tau.use = tau.use,
          azumi = azumi,
          K2.use = K2.use
        )
      
      stap$tau_red <-
        tauTimePeriod(
          period = averagingPeriod,
          smp = stap$red_smp,
          ref = stap$red_ref,
          blk.smp = stap$blk_smp,
          blk.ref = stap$blk_ref,
          smp.p = stap$smp_prs,
          smp.t = stap$smp_tmp,
          smp.flw = stap$smp_flw,
          I0.ref = initValues[2],
          I0.smp = initValues[2],
          I0.blk.smp = initValues[7],
          I0.blk.ref = initValues[8],
          withBlack = blackcorrection,
          stp = stpcorrection,
          tau.use = tau.use,
          azumi = azumi,
          K2.use = K2.use
        )
      
      stap$tau_green <-
        tauTimePeriod(
          period = averagingPeriod,
          smp = stap$grn_smp,
          ref = stap$grn_ref,
          blk.smp = stap$blk_smp,
          blk.ref = stap$blk_ref,
          smp.p = stap$smp_prs,
          smp.t = stap$smp_tmp,
          smp.flw = stap$smp_flw,
          I0.ref = initValues[4],
          I0.smp = initValues[3],
          I0.blk.smp = initValues[7],
          I0.blk.ref = initValues[8],
          withBlack = blackcorrection,
          stp = stpcorrection,
          tau.use = tau.use,
          azumi = azumi,
          K2.use = K2.use
        )
      
      stap$tau_blue <-
        tauTimePeriod(
          period = averagingPeriod,
          smp = stap$blu_smp,
          ref = stap$blu_ref,
          blk.smp = stap$blk_smp,
          blk.ref = stap$blk_ref,
          smp.p = stap$smp_prs,
          smp.t = stap$smp_tmp,
          smp.flw = stap$smp_flw,
          I0.ref = initValues[6],
          I0.smp = initValues[5],
          I0.blk.smp = initValues[7],
          I0.blk.ref = initValues[8],
          withBlack = blackcorrection,
          stp = stpcorrection,
          tau.use = tau.use,
          azumi = azumi,
          K2.use = K2.use
        )
      
      
      
      
      stap <- stap %>%
        select(`date`, red, green, blue ,tau_red,tau_green,tau_blue, everything())
      prfix <-
        sub(pattern = "(.*)\\\\..*$",
            replacement = "\\\\1",
            basename(filename))
      ##### 
      if (blackcorrection == T) { 
        prfix <- paste0(prfix, "_blk")
      }
      if (stpcorrection == T) {
        prfix <- paste0(prfix, "_stp")
      }
      if (azumi == T) {
        prfix <- paste0(prfix, "_azumi")
      }
      if (tau.use == T) {
        prfix <- paste0(prfix, "_tau")
      }
      if (K2.use == T) {
        prfix <- paste0(prfix, "_K2")
      }      
      write_delim(
        stap,
        paste0(prfix, "_", averagingPeriod, "s-avg_proc.dat"),
        append = F,
        col_names = T,
        delim = "\\t",
        na = ""
      )
      setwd(tmp)
    }
    
    return(stap)
  }
