#Alrighty! Here are the models we are using for the 2020 data which 
#focuses on D. tinctorius presence over time! 

#First, let's load some packages
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(sjPlot)

pool_2020 <- read.csv("frogpool20_ECOEVO.csv", header = T, sep= ",") #CHLOE DATA
pool_2020$Year<- 2020
#then we also need to match up pool types. 
pool_2020[pool_2020$Type %in% c("Waka_Pool", "Buttress", "Tree_Hole"), "pool_type"] <- "alive"
pool_2020[pool_2020$Type %in% c("Palm", "Fallen_Tree"), "pool_type"] <- "dead"
pool_2020[pool_2020$Type %in% c("Andrius"), "pool_type"] <- "arboreal"


#############################
#### FP20 TIMESCALES MODELS ---- ####
#################################
#let's subset out arboreal data as well.
#also, let's go ahead and center pH
pool_2020_terr<- pool_2020 %>%
  tidyr::unite("pool_id", PoolID: TreeID, sep = "") %>%
  #watch out! 1 = death! so if tinc is more than 0, give 0!
  dplyr::filter(pool_type != "arboreal") %>%
  tidyr:: drop_na(pH) %>%
  mutate(pool_id = as.factor(pool_id),
         center_pH = pH - mean(pH),
         water_capacity = (4/3*pi* Depth * (Length/2) * (Width/2))/2, #estimate of semi-ellipsoid 
         s_area = (Length/2*Width/2)*pi, #divide length and width by 2 for radii
         sa_depth = s_area/Depth,
         tinc_y_n = ifelse(Dt_Tadpole_Num > 0, 1, 0),
         pool_type = as.factor(pool_type)) 

#Models----
#First, let's look at the effect of pool type + time + (rain?) on pH
#using only drop1

m1<- glmmTMB(pH ~ pool_type * Week * Dt_Tadpole_Num * water_capacity * sa_depth + (1| pool_id), 
             family = gaussian, data = pool_2020_terr)


drop1(m1, test = "Chi") #yeah

m2<- glmmTMB(pH ~ pool_type * Week * Dt_Tadpole_Num + water_capacity * sa_depth + (1| pool_id), 
             family = gaussian, data = pool_2020_terr)

drop1(m2, test = "Chi") #drop water capacity SA interaction

m3<- glmmTMB(pH ~ pool_type * Week * Dt_Tadpole_Num + water_capacity + sa_depth+ (1| pool_id), 
             family = gaussian, data = pool_2020_terr)

drop1(m3, test = "Chi") #drop water cap and sa:depth ratio

m4<- glmmTMB(pH ~ pool_type *Dt_Tadpole_Num * Week  + (1| pool_id), 
             family = gaussian, data = pool_2020_terr)

drop1(m4, test = "Chi")  #drop 3 way interaction

m5<- glmmTMB(pH ~ pool_type * Dt_Tadpole_Num + Week +  (1| pool_id), 
             family = gaussian, data = pool_2020_terr)

m5a<- glmmTMB(pH ~ pool_type + Dt_Tadpole_Num * Week +  (1| pool_id), 
              family = gaussian, data = pool_2020_terr)

m5b<- glmmTMB(pH ~ pool_type* Week + Dt_Tadpole_Num +  (1| pool_id), 
              family = gaussian, data = pool_2020_terr)

drop1(m5, test = "Chi")  #Drop last interaction
drop1(m5a, test = "Chi") 
drop1(m5b, test = "Chi") 

m6<- glmmTMB(pH ~ pool_type + Dt_Tadpole_Num + Week +  (1| pool_id), 
             family = gaussian, data = pool_2020_terr)

drop1(m6, test = "Chi")  #Drop tadpoles

m6a<- glmmTMB(pH ~ pool_type * Week  + (1| pool_id), 
              family = gaussian, data = pool_2020_terr) #let's just make sure that interaction isn't there
drop1(m6a, test = "Chi") #all good.


m7<- glmmTMB(pH ~ pool_type + Week  + (1| pool_id), 
             family = gaussian, data = pool_2020_terr)

drop1(m7, test = "Chi") #all good.

res <- simulateResiduals(m7, plot = T) #looks good
testOverdispersion(m7)
testZeroInflation(m7)

summary(m7)

tab_model(m7, 
          show.est = T, 
          string.est = "Odds Ratios")

tab_model(m7,
          show.r2 = F, 
          show.ngroups = F ,
          show.stat = "z",
          show.icc = F, 
          show.obs = F,
          string.ci = "CI",
          string.est = "Estimates",
          string.stat = "z",
          transform = NULL,
          pred.labels = c("(Intercept)", 
                          "Pool type [Dead]",
                          "Week"))
################################################################################################################
## TINC Presence over time----
###GLMMS tadpole detection over time
#this is only considering resampled pools, and thus excludes canopy pools-

m8<- glmmTMB(Dt_Tadpole_Num ~ pool_type * Week * pH * water_capacity* sa_depth +
               (1| pool_id), family = nbinom2, data = pool_2020_terr)

drop1(m8, test = "Chi") #nwaaaaay too overfit

#alright let's drop one interaction at a time here
#the drop1 is being used after every model here. I am just replacing the model names in 
#the function so I don't have to keep typing out the whole thing.
m9a<- glmmTMB(Dt_Tadpole_Num ~ pool_type * Week * pH + water_capacity * sa_depth +
                (1| pool_id), family = nbinom2, data = pool_2020_terr) #overfit

m9b<- glmmTMB(Dt_Tadpole_Num ~ pool_type * Week * pH + water_capacity + sa_depth +
                (1| pool_id), family = nbinom2, data = pool_2020_terr) #overfit

m9c<- glmmTMB(Dt_Tadpole_Num ~ pool_type * Week * pH  + water_capacity + #drop sa_depth
                (1| pool_id), family = nbinom2, data = pool_2020_terr)

m9d<- glmmTMB(Dt_Tadpole_Num ~ pool_type * Week * pH  + sa_depth + #drop water_capacity
                (1| pool_id), family = nbinom2, data = pool_2020_terr)

m9e<- glmmTMB(Dt_Tadpole_Num ~ pool_type + Week + pH + sa_depth+
               (1| pool_id), family = nbinom2, data = pool_2020_terr) 
#sa_depth and water_capacity are both predictors as proxies for how quickly pools would dry. 
#Where we would imagine that pools with a higher water capacity and a lower surface area to depth
#ratio would be less prone to dessication. Both of these predictors have been attempted to be fit
#in various ways here and just continuously show to not be informative in models.
#We also find the same thing when using a PCA on the 2019 data!

#let's continue on with the rest of the vars.

m9f<- glmmTMB(Dt_Tadpole_Num ~ pool_type * Week * pH  + 
                (1| pool_id), family = nbinom2, data = pool_2020_terr) #no 3way interactio

drop1(m9f, test = "Chi") #no week:pH interaction

m9g<- glmmTMB(Dt_Tadpole_Num ~ pool_type + Week * pH +
                (1| pool_id), family = nbinom2, data = pool_2020_terr) #no interaction betwn week and ph

m9h<- glmmTMB(Dt_Tadpole_Num ~ pool_type * Week + pH +
                (1| pool_id), family = nbinom2, data = pool_2020_terr)

drop1(m9g, test = "Chi") #no pool_type:week interaction

#Alright, so here we go, time is essential to our question and does not
#pull our residuals into space! Very nice. 
m10<- glmmTMB(Dt_Tadpole_Num ~ pool_type * pH + Week+
                 (1| pool_id), family = nbinom2, data = pool_2020_terr) 

drop1(m10, test = "Chi")

#no need to drop week, leads to model instability and also takes out 
#the whole point of the question (how is d. tinc tadpole numbers
#changing over time?)
step(m10, test = "Chi") #keep both.

res <- simulateResiduals(m10, plot = T) 

summary(m10)

tab_model(m10, 
          # transform = NULL, #CHECK OUT TRANSFORM EFFECT. 
          show.est = T, 
          string.est = "Estimate", 
          show.stat = "z", 
          string.stat = " z ", 
          show.icc = F,
          show.zeroinf = F, 
          show.obs = F, 
          show.r2 = F, 
          show.ngroups = F,
          dv.labels = "<I> D. tinctorius </I> tadpoles (count)", 
          pred.labels = c("(Intercept)", 
                          "Pool type [Dead]",
                          "pH",
                          "Week",
                          "Pool type [Dead]: pH"))

