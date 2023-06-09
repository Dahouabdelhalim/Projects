### Getting phenology estimates
#
# In this version, the function lets you estimate one of 3 phenophases (emergence, peak, senscence)
#
# Michael Stemkovski 2019-03-04
#
###

setwd("/home/michael/Documents/Grad School/Research Projects/bee_phenology")

### for pris
#setwd("~/bee_phenology")

#pop_data <- read.csv("Edited data/sp_time_series_2018_10_15_no_sings.csv") # doesn't include Osmia, Andrena, and Hylaeus
pop_data <- read.csv("Edited data/sp_time_series_2019_05_08_no_sings.csv")

#pop_data <- pop_data[which(pop_data$year == 2016 & pop_data$species == "Halictus virgatellus" & pop_data$site == "Hill"),]
#pop_data <- pop_data[which(pop_data$year == 2015 & pop_data$species == "Lasioglossum inconditum" & pop_data$site == "Hill"),]

#library(devtools)
#install_github("willpearse/phest")
#library(phest)
library(mgcv)
library(parallel)

# Functions -----------------------------------------------------------------

### new version that estimates the three phenophases separately
phenophase <- function(abundances, times, phase="onset", threshold=0.05, n_boot=1000, full_boot=FALSE, boot_iter=20){
  
  ### troubleshooting
  #print(paste(abundances,times,phase))
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  if(phase %!in% c("onset","peak","end")) stop("'phase' must be equal to 'onset', 'peak', or 'end'")
  if (max(abundances) == 0 | length(times) < 4){
    warning("Too few observations. Returning NA")
    return(list(phen_bounds = as.data.frame(matrix(NA,nrow=3,ncol=2)),
                phen_mean = as.data.frame(matrix(NA,nrow=1,ncol=2)),
                phen_sd = as.data.frame(matrix(NA,nrow=1,ncol=2)),
                smooth_par = NA, 
                max_pop = NA,
                max_pop_doy = NA,
                first_obs = NA,
                first_obs_doy = NA,
                last_obs = NA,
                last_obs_doy = NA,
                total_pop = NA))} 
  
  if(full_boot == TRUE) max_iter <- n_boot * boot_iter
  else max_iter <- n_boot
  
  time_vec <- min(times):max(times)
  
  # the first model fit, to be used to generate residuals to sample from
  gam_m <- gam(abundances ~ s(times,k=ifelse(length(times) < 5,4,5),bs="cr"))
  spar <- gam_m$sp
  gam_points <- predict.gam(gam_m, newdata=data.frame(times=times), type="response", se=F)
  residuals <- abundances - gam_points
  
  # record of phenology estimates from bootstrapping
  phen_rec <- matrix(data=NA,nrow=n_boot,ncol=2)
  colnames(phen_rec) <- c("estimate","value")
  boot <- 1
  
  for (i in 1:max_iter){
    
    ### troubleshooting
    #sprint(paste(i,"of",max_iter))
    
    if(boot > n_boot) break
    
    # 'artificial' dataset
    res_sample <- sample(residuals,replace=T)
    sim_data <- gam_points + res_sample 
    
    # new model fit
    gam_boot <- gam(sim_data ~ s(times,k=ifelse(length(times) < 5,4,5),bs="cr"))
    boot_fit <- predict.gam(gam_boot, newdata=data.frame(times=time_vec), type="response", se=F)
    
    # prepping to estimate phenophases
    maximum <- max(boot_fit)
    threshold_val <- threshold*maximum
    above_thresh <- which(boot_fit > threshold_val)
    
    # calculating phenophase and magnitude (value) estimates
    if(sum(above_thresh) > 0){
      if(phase == "onset"){
        estimate <- time_vec[min(above_thresh)]
        value <- threshold_val
        if (estimate <= min(times)){
          estimate <- NA
          value <- threshold_val
        } 
      } else if(phase == "end"){
        estimate <- time_vec[max(above_thresh)]
        value <- threshold_val
        if (estimate >= max(times)){
          estimate <- NA
          value <- threshold_val
        }
      } else if(phase == "peak"){
        estimate <- time_vec[which(boot_fit == maximum)]
        value <- maximum
        if (estimate <= min(times) | estimate >= max(times)){
          estimate <- NA
          value <- maximum
        } 
      } else{stop("phase must be equal to 'onset', 'peak', or 'end'")}
    } else {
      estimate <- NA
      value <- NA
    }

    ### troubleshooting
    #print(value)
    
    # checking whether we can record the output of this passbee phenology
    if(full_boot == FALSE | all(!is.na(c(estimate,value)))){
      phen_rec[boot,1] <- estimate # the doy estimate
      phen_rec[boot,2] <- value # the abundance estimate at the doy
      boot <- boot + 1
    }
    
  }
  
  # checking if bootstrapping ended too early, before the full phen_rec was filled
  if(full_boot & any(is.na(phen_rec))){
    warning("Reached iteration limit for full bootstrapping")
    return(list(phen_bounds = as.data.frame(matrix(NA,nrow=3,ncol=2)),
                phen_mean = as.data.frame(matrix(NA,nrow=1,ncol=2)),
                phen_sd = as.data.frame(matrix(NA,nrow=1,ncol=2)),
                smooth_par = NA, 
                max_pop = NA,
                max_pop_doy = NA,
                first_obs = NA,
                first_obs_doy = NA,
                last_obs = NA,
                last_obs_doy = NA,
                total_pop = NA))
  }
  
  # calculating estimate statistics
  phen_bounds <- apply(phen_rec,2,quantile,probs=c(0.025,0.5,0.975),na.rm=TRUE)
  phen_sd <- apply(phen_rec,2,sd,na.rm=TRUE)
  phen_mean <- apply(phen_rec,2,mean,na.rm=TRUE)
  
  # misc data metrics
  max_pop <- max(abundances)
  max_pop_doy <- 
    times[which(abundances == max_pop)][1]
  first_obs <- abundances[which(abundances > 0)][1]
  first_obs_doy <- times[which(abundances > 0)][1]
  last_obs <- tail(abundances[which(abundances > 0)],n=1)
  last_obs_doy <- tail(times[which(abundances > 0)],n=1)
  total_pop <- sum(abundances)
  num_points <- length(which(abundances > 0))
  
  output <- list(phen_bounds = phen_bounds,
                 phen_mean = phen_mean,
                 phen_sd = phen_sd,
                 smooth_par = spar, 
                 max_pop = max_pop,
                 max_pop_doy = max_pop_doy,
                 first_obs = first_obs,
                 first_obs_doy = first_obs_doy,
                 last_obs = last_obs,
                 last_obs_doy = last_obs_doy,
                 total_pop = total_pop) # for diagnistics: phen_rec = phen_rec
  
  return(output)
  
}

##### method exploration zone

# phenophase(pop_data$bowl_adj, pop_data$DOY, n_boot = 100, full_boot = T)
# #hist(a$phen_rec[,1],40)
# phenophase(pop_data$bowl_adj, pop_data$DOY, n_boot = 100, phase="peak")
# phenophase(pop_data$bowl_adj, pop_data$DOY, n_boot = 100, phase="end")
# phenophase(c(0,0,1,4,5,7,2,0,0,0), c(100,110,120,130,140,150,160,170,180,190), n_boot = 100)
# phenophase(c(0,0,1,4,5,7,2,0,0,0), c(100,110,120,130,140,150,160,170,180,190), phase="peak", n_boot = 100)
# phenophase(c(0,0,1,4,5,7,2,0,0,0), c(100,110,120,130,140,150,160,170,180,190), phase="end", n_boot = 100, full_boot = TRUE)$total_pop
# phenophase(c(0,0,1,4,5,7,2,0,0,0), c(100,110,120,130,140,150,160,170,180,190), phase="end", n_boot = 100, full_boot = TRUE, boot_iter = 1)
# 
# phenophase(c(0,0,1), c(100,110,120), n_boot = 100)


### simplified version that doesn't use bootstrapping - and can estimate all three phenophases at once
phenophase.bootless <- function(abundances, times, phase="all", threshold=0.05, prediction_int=0.01, spline_basis="cr", expand_time=FALSE, plot=FALSE, family="gaussian"){
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  if(phase %!in% c("onset","peak","end","all")) stop("'phase' must be equal to 'onset', 'peak', 'end', or 'all'")
  if(length(abundances) != length(times)) stop("The length of the abundance and time vectors must be equal.")
  if (max(abundances) == 0 | length(times) < 4){
    output <- rep(NA,9)
    names(output) <- c("onset", "onset_val", "onset_se", "peak", "peak_val", "peak_se", "end", "end_val", "end_se")
    warning("Too few observations. Must have at least 4 values in the time-series with a least one non-zero abundance value. Returning NA")
    return(as.list(output))
  }
  
  ### troubleshooting
  #print(paste(abundances,times,phase))
  
  # phase <- "all"
  # prediction_int <- 0.01
  # threshold <- 0.05
  # spline_basis <- "cr"
  # expand_time <- TRUE
  # plot <- TRUE
  # family="gaussian"
  
  #times <- pop_data$DOY
  
  #times <- seq(100,110,length.out = 6)
  
  #abundances <- pop_data$bowl_adj
  #abundances <- c(0,1,3,6,2,0)
  #abundances <- c(0,0,1,2,1,0)
  #abundances <- c(0,0,0,1,1,0)
  #abundances <- c(4,7,2,1,0,0)
  
  # gives a below zero curve
  # abundances <- c(0,1,2,1.5,0,0)
  # times <- 1:length(abundances)
  
  # gave error
  # abundances <- c(0.002382999, -0.002287320, 0.022773949, 0.041828991, 0.011120786, 0.006473611)  
  # times <- c(128, 142, 156, 170, 184, 198)
  
  ab_sd <- sd(abundances)
  time_sd <- sd(times)
  
  abundances <- scale(abundances, center = FALSE, scale = ab_sd)
  times <- scale(times, center = FALSE, scale = time_sd)
  if(expand_time){
    time_vec <- seq(min(times)-sd(times), max(times)+sd(times),by=prediction_int/time_sd)
  } else{
    time_vec <- seq(min(times),max(times),by=prediction_int/time_sd)
  }
  
  if(family == "gaussian"){
    gam_m <- gam(abundances ~ s(times,k=ifelse(length(times) < 5,4,5),bs=spline_basis))
  } else if(family == "nb"){
    gam_m <- gam(abundances ~ s(times,k=ifelse(length(times) < 5,4,5),bs=spline_basis), family=nb(theta = NULL, link = "log"))
  } else{
    gam_m <- gam(abundances ~ s(times,k=ifelse(length(times) < 5,4,5),bs=spline_basis))
  }
  
  gam_summary <- summary(gam_m)
  #se <- max(gam_summary$se)
  gam_fit <- as.numeric(predict.gam(gam_m, newdata=data.frame(times=time_vec), type="response", se=FALSE))
  # gam_fit[which(gam_fit<0)] <- 0 # rounding up negatives to zero because there's no such thing as negative abundance
  
  maximum <- max(gam_fit)
  threshold_val <- threshold*maximum
  above_thresh <- which(gam_fit > threshold_val)
  
  #plot(abundances ~ times, xlim=range(time_vec))
  #lines(gam_fit ~ time_vec, type="l")
  
  # calculating phenophase and magnitude (value) estimates
  if(sum(above_thresh) > 0){
    peak <- time_vec[which(gam_fit == maximum)]
    peak_pred <- predict.gam(gam_m, newdata=data.frame(times=peak), type="response", se=TRUE)
    peak_val <- peak_pred$fit
    peak_se <- peak_pred$se.fit
    
    time_elements_before_peak <- which(time_vec < peak)
    time_elements_after_peak <- which(time_vec > peak)
    
    onset <- time_vec[min(above_thresh)]
    if(any(gam_fit[which(time_vec == onset):which(time_vec == peak)] <= 0)){
      # test if any predicted values cross zero between max and onset
      onset <- max(time_vec[which(gam_fit[time_elements_before_peak] <= threshold_val)]) 
    } 
    onset_pred <- predict.gam(gam_m, newdata=data.frame(times=onset), type="response", se=TRUE)
    onset_val <- onset_pred$fit
    onset_se <- onset_pred$se.fit
    
    end <- time_vec[max(above_thresh)]
    if(any(gam_fit[which(time_vec == peak):which(time_vec == end)] <= 0)){
      # test if any predicted values cross zero between max and onset
      peak_element <- which(gam_fit == maximum)
      first_zero_after_peak <- min(length(time_elements_before_peak)+which(gam_fit[time_elements_after_peak] <= 0))
      last_threshold_element <- max(peak_element + which(gam_fit[peak_element:first_zero_after_peak] > threshold_val))
      end <- time_vec[last_threshold_element]
    }
    end_pred <- predict.gam(gam_m, newdata=data.frame(times=end), type="response", se=TRUE)
    end_val <- end_pred$fit
    end_se <- end_pred$se.fit
    
    if(onset <= min(time_vec)){onset <- NA; onset_val <- NA; onset_se <- NA} 
    if(peak <= min(time_vec) | peak >= max(time_vec)){peak <- NA; peak_val <- NA; peak_se <- NA} 
    if(end >= max(time_vec)){end <- NA; end_val <- NA; end_se <- NA} 
  } else{
    output <- rep(NA,9)
    names(output) <- c("onset", "onset_val", "onset_se", "peak", "peak_val", "peak_se", "end", "end_val", "end_se")
    return(as.list(output))
  }
  
  #paste("onset",onset,onset_val,onset_se,"peak",peak,peak_val,peak_se,"end",end,end_val,end_se)
  
  
  # back-transforming from scale
  output_scaled <- c(onset, peak, end, onset_se, peak_se, end_se, onset_val, peak_val, end_val)
  
  output <- rep(NA, 9)
  output[1:6] <- sapply(output_scaled[1:6], function(x) x*time_sd)
  output[7:9] <- sapply(output_scaled[7:9], function(x) x*ab_sd)
  names(output) <- c("onset", "peak", "end", "onset_se", "peak_se", "end_se", "onset_val", "peak_val", "end_val")
  
  if(plot){
    #par(mfrow=c(1,2))
    # gam_fit_se <- predict.gam(gam_m, newdata=data.frame(times=time_vec), type="response", se=TRUE)
    # plot(gam_fit_se$fit ~ time_vec)
    # plot(gam_fit_se$se.fit ~ time_vec)
    
    #plot.gam(gam_m, residuals=TRUE, all.terms = T)
    
    plot(times*time_sd, abundances*ab_sd, xlim=range(time_vec*time_sd))
    lines(time_vec*time_sd, gam_fit*ab_sd)
    abline(v=c(output[1], output[2], output[3]), col=c("blue", "green", "brown"))
    abline(v=c(output[1]-output[4], output[1]+output[4]),lty=2,col="blue")
    abline(v=c(output[2]-output[5], output[2]+output[5]),lty=2,col="green")
    abline(v=c(output[3]-output[6], output[3]+output[6]),lty=2,col="brown")
  }
  
  return(as.list(output))
  
}


# pop_data <- pop_data[which(pop_data$year == 2016 & pop_data$species == "Halictus virgatellus" & pop_data$site == "Hill"),]
# abs <- pop_data$bowl_adj
# times <- pop_data$DOY
# 
# a <- rnorm(1000)
# quantileCI(a, prob=0.025)
# quantileCI(a, prob=0.975)
# quantile(a,probs=c(0.025,0.5,0.975)) # what I'm already doing
# 
# b <- rep(c(-1,0,1), each=333)
# quantileCI(b, prob=0.025)
# quantileCI(b, prob=0.975)
# quantile(b,probs=c(0.025,0.5,0.975)) # what I'm already doing
# 
# abs <- c(0,2,5,10,2,3,0) # made up abundances
# times <- seq(100, 160, by=10) # made up times
# times_expanded <- rep(time, times=abs)
# 
# quantileCI(times_expanded, prob=0.05) # emergence?
# quantileCI(times_expanded, prob=0.5) # peak?
# quantileCI(times_expanded, prob=0.95) # senescence?
# 
# hist(times_expanded)
# plot(abs ~ times)
# abline()


#####

# Estimating phenologies, the new faster way -----------------------------------------------------------------

# wrapper to pass elements of a list to the phenophase function
phenophase.wrapper <- function(x, ...){
  args <- list(...)
  abundances <- x$bowl_adj
  times <- x$DOY
  # for boot phenophase estimator:
  #output <- phenophase(abundances, times, ...) # might have to be "args" instead of "..." (maybeeee???)
  
  # for bootless:
  output <- phenophase.bootless(abundances, times, ...)
  
  #output <- tryCatch(phenophase(abundances, times, ...), error=function(x) "i've gone wrong")
  #if(output=="i've gone wrong")
  #  return(list(error=x)) #which iteration gave me an error
  
  # this was causing an error, it looks like
  # write.table(paste(args$phase,paste(x$species,x$site,x$year,sep="_")[1],"mean:",output[[2]][1],"sd:",output[[3]][1],sep="_"),
  #             paste("Edited data/pheno_tracker/",paste(args$phase,paste(x$species,x$site,x$year,sep="_")[1],".txt",sep="_"),sep=""),
  #             row.names = FALSE, col.names = FALSE)
  print(paste(args[[1]],x$species,x$site,x$year)[1]) # to roughly track progress
  return(output)
}

unique_records <- unique(pop_data[,1:3])
u_rows <- 1:nrow(unique_records)


### this is very promising - work from this
pop_data_split <- split(pop_data, list(pop_data$site,pop_data$species,pop_data$year), drop=TRUE)
#pop_data_split$`Davids.Calliopsis teucrii.2009`

#pop_data_split <- pop_data_split[160:163] #160 seems to be where the problem is happening
#pop_data_split <- pop_data_split[161] #160 seems to be where the problem is happening
#Error in output[[2]] : subscript out of bounds
# In addition: Warning message:
#   In phenophase(abundances, times, ...) :
#   Show Traceback
# 
# Rerun with Debug
# Error in output[[2]] : subscript out of bounds

#length(pop_data_split)
#str(pop_data_split)




##### a test run
#a <- pop_data_split[1:5]
#a <- pop_data_split[49]
#output_em <- lapply(a, phenophase.wrapper, phase="onset", n_boot=100, full_boot=TRUE)
#output_peak <- lapply(a, phenophase.wrapper, phase="peak", n_boot=100, full_boot=TRUE)
#output_sen <- lapply(a, phenophase.wrapper, phase="end", n_boot=100, full_boot=TRUE)
#output_em <- mclapply(a, phenophase.wrapper, phase="onset", n_boot=100, full_boot=TRUE, mc.cores = 6)
#output_peak <- mclapply(a, phenophase.wrapper, phase="peak", n_boot=100, full_boot=TRUE, mc.cores = 6)
#output_sen <- mclapply(a, phenophase.wrapper, phase="end", n_boot=100, full_boot=TRUE, mc.cores = 6)

#### this is what I need to run ultimately
#output_em <- lapply(pop_data_split, phenophase.wrapper, phase="onset", n_boot=1000, full_boot=TRUE, boot_iter=20)
#output_peak <- lapply(pop_data_split, phenophase.wrapper, phase="peak", n_boot=1000, full_boot=TRUE, boot_iter=20)
#output_sen <- lapply(pop_data_split, phenophase.wrapper, phase="end", n_boot=1000, full_boot=TRUE, boot_iter=20)
#output_em_mc <- mclapply(pop_data_split, phenophase.wrapper, phase="onset", n_boot=1000, full_boot=TRUE, boot_iter=20, mc.cores = 4)
#output_peak <- mclapply(pop_data_split, phenophase.wrapper, phase="peak", n_boot=1000, full_boot=TRUE, boot_iter=20, mc.cores = 4)
#output_sen <- mclapply(pop_data_split, phenophase.wrapper, phase="end", n_boot=1000, full_boot=TRUE, boot_iter=20, mc.cores = 4)
output_em <- mclapply(pop_data_split, phenophase.wrapper, phase="onset", n_boot=1000, full_boot=TRUE, boot_iter=20, mc.cores = 12)
save.image("Edited data/em_done3.RData")
output_peak <- mclapply(pop_data_split, phenophase.wrapper, phase="peak", n_boot=1000, full_boot=TRUE, boot_iter=20, mc.cores = 12)
save.image("Edited data/peak_done3.RData")
output_sen <- mclapply(pop_data_split, phenophase.wrapper, phase="end", n_boot=1000, full_boot=TRUE, boot_iter=20, mc.cores = 12)
save.image("Edited data/sen_done3.RData")

# output <- sapply(a, phenophase.wrapper, n_boot=100, full_boot=TRUE)

# sum(is.na(output_em))
# sum(is.na(output_peak))
# sum(is.na(output_sen))
# 
# -which(is.na(output_sen))
# 
# as.numeric(which(is.na(output_em)))
# 
# length(which(sapply(output_em, inherits, what = "try-error"))) #only errors in emergence... weird 
# length(which(sapply(output_peak, inherits, what = "try-error")))
# length(which(sapply(output_sen, inherits, what = "try-error")))
# 
# first_error <- which(sapply(output_em, inherits, what = "try-error"))[1]
# errors <- which(sapply(output_em, inherits, what = "try-error"))
# output_em[first_error]
# output_peak[first_error]
# 
# output_em[errors] # looks like it was all "object 'value' not found"... I think I corrected it now... pleeeease work
# 
# problem_data <- pop_data_split[errors[1:4]]
# please_work <- mclapply(problem_data, phenophase.wrapper, phase="onset", n_boot=1000, full_boot=TRUE, boot_iter=20, mc.cores = 4)
# 
# mclapply(pop_data_split[160], phenophase.wrapper, phase="onset", n_boot=1000, full_boot=TRUE, boot_iter=20, mc.cores = 4)
# output_em <- output_em[-which(is.na(output_em))]
# output_em <- output_em[-which(is.na(output_em))]


##### bootless estimation

pop_data_split <- split(pop_data, list(pop_data$site,pop_data$species,pop_data$year), drop=TRUE)

dim(pop_data)
head(pop_data)

length(pop_data_split)

out_put_all <- lapply(pop_data_split, phenophase.wrapper, phase="all", plot=FALSE)

#####


data_summary_names <- c("year","site","species",
                        "em_doy","em_val","em_doy_sd","em_val_sd",
                        "peak_doy","peak_val","peak_doy_sd","peak_val_sd",
                        "sen_doy","sen_val","sen_doy_sd","sen_val_sd")

data_summary <- as.data.frame(matrix(nrow=length(output_em), ncol=length(data_summary_names)))
names(data_summary) <- data_summary_names

for(i in 1:length(output_em)){
  em <- output_em[i]
  peak <- output_peak[i]
  sen <- output_sen[i]

  info <- strsplit(names(em),"\\\\.")
  site <- info[[1]][1]
  species <- info[[1]][2]
  year <- info[[1]][3]

  data_summary_row <- rep(NA,length(data_summary_names))
  names(data_summary_row) <- data_summary_names

  #info
  data_summary_row["year"] <- year
  data_summary_row["site"] <- site
  data_summary_row["species"] <- species

  # emergence
  data_summary_row["em_doy"] <- em[[1]]$phen_mean[1]
  data_summary_row["em_val"] <- em[[1]]$phen_mean[2]
  data_summary_row["em_doy_sd"] <- em[[1]]$phen_sd[1]
  data_summary_row["em_val_sd"] <- em[[1]]$phen_sd[2]

  # peak
  data_summary_row["peak_doy"] <- peak[[1]]$phen_mean[1]
  data_summary_row["peak_val"] <- peak[[1]]$phen_mean[2]
  data_summary_row["peak_doy_sd"] <- peak[[1]]$phen_sd[1]
  data_summary_row["peak_val_sd"] <- peak[[1]]$phen_sd[2]

  # peak
  data_summary_row["sen_doy"] <- sen[[1]]$phen_mean[1]
  data_summary_row["sen_val"] <- sen[[1]]$phen_mean[2]
  data_summary_row["sen_doy_sd"] <- sen[[1]]$phen_sd[1]
  data_summary_row["sen_val_sd"] <- sen[[1]]$phen_sd[2]

  data_summary[i,] <- data_summary_row
}

head(data_summary)

write.csv(data_summary,"Edited data/data_summary_2019_05_08.csv",row.names = F)

save.image("Edited data/all_done_2019_05_08.RData")

#write.csv(data_summary,"Edited data/data_summary_2019_03_26_split_101_101.csv",row.names = F)


# rm nohup.out
# nohup Rscript phenophase_estimator.R &
# ls -l | wc -l





