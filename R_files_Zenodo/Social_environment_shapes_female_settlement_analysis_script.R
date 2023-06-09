# Script to perform resource selection function, model selection, and
# mother overlap Wilcoxon rank test in support of:
# "Social environment shapes female settlement decisions in a solitary carnivore"
# Hansen, J.E., Hertel, A.G., Frank, S.C., Kindberg, J., Zedrosser, A. 2021


# Load required packages --------------------------------------------------

library(plyr)
library(tidyverse)
library(car)
library(glmmTMB)
library(MuMIn)


# Import data -------------------------------------------------------------

social_df <- read.csv("social_environment.csv")

overlap_df <- read_csv("mother_overlap.csv")

# Scale continuous predictor variables ------------------------------------

social_df <- social_df %>% 
  mutate(across(c(relRatio, famIx, densDiff), scale))


# Check for multicollinearity ---------------------------------------------

cor(social_df[,c(4:7)], method = "spearman")

vif(lm(used ~ famIx + relRatio  + momPresent + densDiff, data = social_df))


# Fit base model and check interactions -----------------------------------

social_glmm <- glmmTMB(used ~ famIx + relRatio  +
                        momPresent + densDiff + (1|bearID),
                      family = binomial(),
                      data = social_df,
                      na.action = na.fail)
summary(social_glmm)



int_1 <- glmmTMB(used ~ famIx +  densDiff +
                   momPresent *  relRatio + (1|bearID),
                 family = binomial(),
                 data = social_df,
                 na.action = na.fail)
summary(int_1)


int_2 <- glmmTMB(used ~ relRatio  + densDiff +
                   momPresent *  famIx + (1|bearID),
                 family = binomial(),
                 data = social_df,
                 na.action = na.fail)
summary(int_2)


int_3 <- glmmTMB(used ~ famIx + relRatio +
                   momPresent * densDiff   + (1|bearID),
                 family = binomial(),
                 data = social_df,
                 na.action = na.fail)
summary(int_3)


int_4 <- glmmTMB(used ~  momPresent + densDiff +
                   famIx * relRatio + (1|bearID),
                 family = binomial(),
                 data = social_df,
                 na.action = na.fail)
summary(int_4)


int_5 <- glmmTMB(used ~  relRatio + momPresent +
                   famIx * densDiff + (1|bearID),
                 family = binomial(),
                 data = social_df,
                 na.action = na.fail)
summary(int_5)


int_6 <- glmmTMB(used ~  famIx + momPresent +
                   relRatio * densDiff + (1|bearID),
                 family = binomial(),
                 data = social_df,
                 na.action = na.fail)
summary(int_6)

AICc(social_glmm, int_1, int_2, int_3, int_4, int_5, int_6)


rm(int_1, int_2, int_3, int_4, int_5, int_6)


# Fit all subsets, model selection, and averaging -------------------------


soc_sub <- dredge(social_glmm, evaluate = TRUE, rank = "AICc") # fit/rank all subsets

soc_list <- get.models(soc_sub, subset = delta < 2) # select top model set

soc_avg <- model.avg(soc_list, fit = TRUE, revised.var = TRUE) # average over top models

summary(soc_avg)
confint(soc_avg)

r.squaredGLMM(social_glmm)

rm(soc_avg, soc_list, soc_sub, social_glmm)


# Permutations for variable importance ------------------------------------

# NB! This will take several hours to run!


# setup for permutation process

set.seed(42)
temp_mat <- matrix(nrow= 1000, ncol=17)
permutation_nb <- 1000 # permutation number
bearID <- unique(as.character(social_df$bearID))

for (i in 1:permutation_nb) {
  
  print(paste("Permutation Number", i))
  
  # permutation of y -------------------------------------
  
  sampleID <- sample(bearID, size = 51, replace = TRUE)
  sampleDF <- social_df[social_df$bearID %in% sampleID, ]
  
  # model averaging --------------------------------------
  
  reg <- glmmTMB(used ~ famIx + relRatio +
                   momPresent + densDiff + (1|bearID),
                 family = binomial(),
                 data = sampleDF,
                 na.action = na.fail)
  
  ms <- dredge(reg, evaluate = TRUE, rank = "AICc")
  confset <- ms[ms$delta <=2,]
  
  # control structure for model set < 2
  
  if(nrow(confset) > 1) {
    
    
    avgmod <- model.avg(confset)
    print(summary(avgmod))
    ci <- confint(avgmod)
    confint(avgmod)
    
    temp_mat[i,1] <- i
    temp_mat[i,2] <- summary(avgmod)$coefficients[8] # momPresent beta (cond. average)
    temp_mat[i,3:4] <-  ci[c(4,9)] # 95% CI momPresent
    temp_mat[i,5] <-  summary(avgmod)$sw["cond(momPresent)"] # sw momPresent
    temp_mat[i,6] <- summary(avgmod)$coefficients[6] # famIx beta
    temp_mat[i,7:8] <- ci[c(3,8)] # 95% CI famIx
    temp_mat[i,9] <- summary(avgmod)$sw["cond(famIx)"] # sw famIx
    temp_mat[i,10] <- summary(avgmod)$coefficients[10] # relRatio beta
    temp_mat[i,11:12] <- ci[c(5,10)] # 95% CI relRatio
    temp_mat[i,13] <- summary(avgmod)$sw["cond(relRatio)"] # sw relRatio
    temp_mat[i,14] <- summary(avgmod)$coefficients[4] # densDiff beta
    temp_mat[i,15:16] <- ci[c(2,7)] # 95% CI densDiff
    temp_mat[i,17] <- summary(avgmod)$sw["cond(densDiff)"] # sw densDiff
    
    
    ### computation progress (it could be very slow) -------
    
    cat(round(100*i/permutation_nb,2), "% \\n")
    flush.console()  
    
  } else {
    
    
    # only one model within delta 2 AICc, no averaging product
    
    
    temp_mat[i,1] <- i
    temp_mat[i,2:17] <- NA
    
    ### computation progress (it could be very slow) -------
    
    cat(round(100*i/permutation_nb,2), "% \\n")
    flush.console()  
    
    
  }
  
  
}


# change matrix to df and name columns

avg_df <- as.data.frame(temp_mat)

colnames(avg_df) <- c("Permutation", "momPresent_beta", "momPresent_ciLow",
                     "momPresent_ciHigh", "momPresent_sw", "famIx_beta", 
                     "famIx_ciLow", "famIx_ciHigh", "famIx_sw", 
                     "relRatio_beta", "relRatio_ciLow", "relRatio_ciHigh",
                     "relRatio_sw", "densDiff_beta", "densDiff_ciLow", "
                     densDiff_ciHigh", "densDiff_sw")

# get mean sw values

mean(avg_df[!(is.na(avg_df$momPresent_sw)),5]) # 0.9556886
mean(avg_df[!(is.na(avg_df$famIx_sw)),9]) #  0.9979466
mean(avg_df[!(is.na(avg_df$relRatio_sw)),13]) # 0.3448147
mean(avg_df[!(is.na(avg_df$densDiff_sw)),17]) # 0.9954417



# Analysis of home range overlap, alive versus dead mothers ---------------

alive_mother <- overlap_df %>% 
  filter(momStatus == "Alive",
         OGM > 0)
  
dead_mother <- overlap_df %>% 
  filter(momStatus == "Dead",
         OGM > 0)

wilcox.test(alive_mother$OGM, dead_mother$OGM, paired = FALSE)


# summary statistics

summary_overlap <- ddply(overlap_df, .(momStatus), summarise, 
                         mean=mean(OGM), 
                         median=median(OGM), 
                         sd = sd(OGM))
summary_overlap
