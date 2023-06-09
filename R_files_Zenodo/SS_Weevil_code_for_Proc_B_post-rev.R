#START#####

#R code for the analysis in the manuscript "Social selection is density dependent but makes little contribution to total selection in New Zealand giraffe weevils"
#Data are provided as supplemental files
#Data are available at https://datadryad.org/stash/dataset/doi:10.5061/dryad.8vs71c3

#Load packages####
library(tidyverse)
library(lubridate)
library(car)
library(glmmTMB)
library(ggplot2)
library(grid)
library(gridExtra)

#Load and wrangle data####

setwd(....) #to be completed by user

weevil_raw_m = read.csv("weevil_obs_m.csv", header=T)#pre-xmas males
weevil_raw_f = read.csv("weevil_obs_f.csv", header=T)#pre-xmas females
weevil_raw_x_m = read.csv("weevil_obs_x_m.csv", header=T)#post-xmas males
weevil_raw_x_f = read.csv("weevil_obs_x_f.csv", header=T)#post-xmas females
measure_raw_x = read.csv("measurements_x.csv", header=T) #as measurements for post-xmas are in a separate file 


weevil_f = weevil_raw_f %>% mutate(sex = "F", season="pre", weevil=as.character(weevil), tree=as.character(tree)) %>%
  select(weevil, sex, date, time, season, tree, body_length,  num_mates = num_males, total_cop) 

weevil_m = weevil_raw_m %>% mutate(sex = "M", season="pre", weevil=as.character(weevil), tree=as.character(tree)) %>%
  select(weevil, sex, date, time, season, tree, body_length, num_mates = num_females, total_cop)   

measure_x = measure_raw_x %>% select(weevil, sex, body_length) %>%
  mutate(weevil = as.character(weevil))

weevil_x_f = weevil_raw_x_f %>% mutate(season="post") %>%
  left_join(measure_x, by = c("weevil")) %>%
  select(weevil,sex, date, time, season, tree, body_length, num_mates = num_males, total_cop)

weevil_x_m = weevil_raw_x_m  %>% mutate(season="post") %>%
  left_join(measure_x, by = c("weevil"))%>%
  select(weevil,sex, date, time, season, tree, body_length, num_mates = num_females, total_cop)

weevil_comb_both = rbind(weevil_f,weevil_m,weevil_x_f,weevil_x_m) 

for (i in 1:length(1:nrow(weevil_comb_both))) {
  weevil_comb_both$neigh_mean_blength[i] = mean(weevil_comb_both$body_length[weevil_comb_both$tree==weevil_comb_both$tree[i] & #same tree
                                                                       weevil_comb_both$date==weevil_comb_both$date[i] & #same day, with date don't need to also specify season
                                                                         weevil_comb_both$sex==weevil_comb_both$sex[i] & #same sex
                                                                        weevil_comb_both$weevil!=weevil_comb_both$weevil[i] ]) #NOT including focal individual
  
}

for (i in 1:length(1:nrow(weevil_comb_both))) {
  weevil_comb_both$sex_ratio[i] = nrow(weevil_comb_both[weevil_comb_both$tree==weevil_comb_both$tree[i] & #same tree
                                                        weevil_comb_both$date==weevil_comb_both$date[i] & #same day, with date don't need to also specify season
                                                        weevil_comb_both$sex=="M" , ]) / #counting males so higher SR = more male biased
                                                                                  #then divide by number of individuals
                                          nrow(weevil_comb_both[weevil_comb_both$tree==weevil_comb_both$tree[i] & #same tree
                                                                weevil_comb_both$date==weevil_comb_both$date[i]  #same day, with date don't need to also specify season
                                                                , ])
}

for (i in 1:length(1:nrow(weevil_comb_both))) {
  weevil_comb_both$density[i] = nrow(weevil_comb_both[
    weevil_comb_both$tree==weevil_comb_both$tree[i] & #same tree
    weevil_comb_both$date==weevil_comb_both$date[i] , ]) #& same day
}



#who was measured once?
weevil_once = weevil_comb_both %>%
  count(weevil) %>%
  filter(n == 1)

#who was only measured in the last week for each season?

weevil_last = weevil_comb_both %>% 
  group_by(season) %>%
  mutate(date = dmy(date),
         last_week = max(date) - 7) %>%
  group_by(weevil) %>%
  mutate(first_obs = min(date)) %>%
  filter(first_obs > last_week)


weevil_comb_complete = weevil_comb_both %>% filter(!is.na(body_length),
                                                   !is.na(neigh_mean_blength),
                                                   !is.infinite(sex_ratio),
                                                   !weevil %in% weevil_once$weevil,
                                                   !weevil %in% weevil_last$weevil) %>%
  mutate(body_length_z = scale(body_length),
         body_length_z2 = 0.5*body_length_z^2, #mean centre, then square, then half
         neigh_mean_blength_z = scale(neigh_mean_blength),
         neigh_mean_blength_z2 = 0.5*neigh_mean_blength_z^2,
         time2 = as.numeric(hms(time)),
         time2_z = scale(time2),
         time2_z2 = 0.5*time2_z^2)

#how many of each sex?
weevil_comb_complete %>%
  group_by(sex) %>%
  summarise(counts = n_distinct(weevil))

#not using an offset of mean fitness as in log-linear models using absolute fitness the regression coefficients on the predictor scale are directional selection gradients 
#see: https://www.biorxiv.org/content/10.1101/040618v1.full

#Is either sex under linear or quadratic social selection for body size? Model 1####


blength_mlselec_m1 = glmmTMB(num_mates ~  sex*(body_length_z + neigh_mean_blength_z +
                                                body_length_z2 + neigh_mean_blength_z2) +
                                                time2_z + time2_z2 +
                                                 (1|date) + (1|weevil) + (1|tree) ,
                              family="poisson",  
                             data = weevil_comb_complete)          
summary(blength_mlselec_m1)
Anova(blength_mlselec_m1, type="II")

hist(resid(blength_mlselec_m1)) #long tail but otherwise OK
plot(resid(blength_mlselec_m1), fitted(blength_mlselec_m1)) #looks like what you'd expect from a Pois model, bands of scores but no obvious directional pattern

#Does social selection interact with body size? Model 2####

blength_mlselec_m4 = glmmTMB(num_mates ~ sex * body_length_z * neigh_mean_blength_z +
                               time2_z + time2_z2 + 
                               (1|date) + (1|weevil) + (1|tree), 
                             family="poisson", 
                             data = weevil_comb_complete)            
summary(blength_mlselec_m4)
Anova(blength_mlselec_m4, type="II")

#Does social selection interact with male strategy? Model 3####

#Make new variable, guarding male (>40mm), sneaking male (<=40mm), or female
weevil_comb_complete$sex2 = ifelse(weevil_comb_complete$sex=="F", "female",
                                   ifelse(weevil_comb_complete$body_length>40, "Gmale", "Smale"))

blength_mlselec_m5 = glmmTMB(num_mates ~ sex2*(body_length_z + neigh_mean_blength_z) +
                               time2_z + time2_z2 + 
                               (1|date) + (1|weevil) + (1|tree), 
                             family="poisson", 
                             data = weevil_comb_complete)            
summary(blength_mlselec_m5)
Anova(blength_mlselec_m5, type="II")

#Does social selection interact with density? Model 4####

blength_mlselec_m6 = glmmTMB(num_mates ~ sex*scale(density)*
                               (body_length_z + neigh_mean_blength_z ) +
                               time2_z + time2_z2 +
                               (1|date) + (1|weevil) + (1|tree), 
                             family="poisson", 
                             data = weevil_comb_complete)            
summary(blength_mlselec_m6)
Anova(blength_mlselec_m6, type="II")

#Does social selection interact with sex ratio? Model 5####

blength_mlselec_m3 = glmmTMB(num_mates ~ sex*scale(sex_ratio)*
                               (body_length_z + neigh_mean_blength_z) +
                                time2_z + time2_z2 + 
                                 (1|date) + (1|weevil) + (1|tree), 
                             family="poisson",
                             data = weevil_comb_complete)         
summary(blength_mlselec_m3)
Anova(blength_mlselec_m3, type="II")


#Plot how social selection changes with density####

with(weevil_comb_complete, hist(density)) #natural break after 40, 

weevil_comb_complete$densityf = ifelse(weevil_comb_complete$density>40, "High", "Low")

weevil_comb_complete$cop_pred = predict(blength_mlselec_m6, type="response") #using model with density

num_mates_den_f = ggplot(weevil_comb_complete, 
                         aes(x = neigh_mean_blength, 
                             y = num_mates, 
                             col = densityf)) +
  geom_jitter(height=0.05) + 
  theme_classic() +
  scale_color_manual(values = c("dodgerblue", "black")) +
  xlab("Mean size of rivals (mm)") + 
  ylab ("Number of mating partners") +
  labs(col = "Density") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.key.size =  unit(0.5, "in")) +
  stat_smooth(inherit.aes=F, aes(x = neigh_mean_blength, 
                                 y = cop_pred, 
                                 group = densityf, 
                                col=densityf),
              show.legend = T,
              method = "glm", se = T,
              method.args = list(family = "poisson")) +
  facet_wrap(~sex, scale="free") 


#Estimating the interactant covariance####

with(weevil_comb_complete[weevil_comb_complete$sex=="F",], 
     cor.test(scale(body_length),scale(neigh_mean_blength)))
with(weevil_comb_complete[weevil_comb_complete$sex=="M",],
     cor.test(scale(body_length),scale(neigh_mean_blength))) 

#test whether assortment changes with density with a model

weevil_comb_complete$neigh_mean_blength_z = as.numeric(weevil_comb_complete$neigh_mean_blength_z)

assort_m1 = glmmTMB(neigh_mean_blength_z ~ sex * body_length_z * scale(density) + 
                               (1|date) + (1|weevil) + (1|tree), 
                             data = weevil_comb_complete)            
summary(assort_m1)
Anova(assort_m1, type="II")
#sig body length-density interaction, so relationship gets more positive as density increase. Starts negative, changes to be positive

#Plot how phenotypic assortment changes with density####

assort_den_f = ggplot(weevil_comb_complete, 
                         aes(x = body_length,
                             y = neigh_mean_blength, 
                             col = densityf)) +
  geom_jitter(height=0.05) + 
  theme_classic() +
  scale_color_manual(values = c("dodgerblue", "black")) +
  ylab("Mean size of rivals (mm)") + 
  xlab ("Body length (mm)") +
  labs(col = "Density") +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        legend.key.size =  unit(0.5, "in")) +
  stat_smooth(aes(group = densityf, 
                  col = densityf), 
              show.legend = T,
              method = "lm", se = T) +
  facet_wrap(~sex, scale="free") 


#So the negative SS at high densities may have limited impact if assortment is near zero and stays that way across densities

#What are the overall selection differentials?####

#Females:
(fixef(blength_mlselec_m1)[[1]][3] + 
   fixef(blength_mlselec_m1)[[1]][4]* with(weevil_comb_complete[weevil_comb_complete$sex=="F",], 
                                           cor.test(scale(body_length),scale(neigh_mean_blength))$est) )*
  var(weevil_comb_complete$body_length[weevil_comb_complete$sex == "F"])
# = 5.568

#Males:
(fixef(blength_mlselec_m1)[[1]][3] + fixef(blength_mlselec_m1)[[1]][9] +
    (fixef(blength_mlselec_m1)[[1]][4] + fixef(blength_mlselec_m1)[[1]][10])* 
    with(weevil_comb_complete[weevil_comb_complete$sex=="M",], 
         cor.test(scale(body_length),scale(neigh_mean_blength))$est) )*
  var(weevil_comb_complete$body_length[weevil_comb_complete$sex == "M"])
# = 41.177



#Plot how selection differentials change with density####


#range density (scaled) from min to max
density = seq(min(as.numeric(scale(weevil_comb_complete$density))),
              max(as.numeric(scale(weevil_comb_complete$density))),0.1)

#Selection differentials for females
ds_f =  fixef(blength_mlselec_m6)[[1]][4] + density * fixef(blength_mlselec_m6)[[1]][11]
ss_f = fixef(blength_mlselec_m6)[[1]][5] + density * fixef(blength_mlselec_m6)[[1]][12]
cij_f = fixef(assort_m1)[[1]][3] + density * fixef(assort_m1)[[1]][7]
var_f = var(weevil_comb_complete$body_length[weevil_comb_complete$sex == "F"])
S_f = var_f*(ds_f+ss_f*cij_f) #gradual decrease with density

#Selection differentials for males
ds_m = fixef(blength_mlselec_m6)[[1]][4] + fixef(blength_mlselec_m6)[[1]][9] +
  density * (fixef(blength_mlselec_m6)[[1]][11] + fixef(blength_mlselec_m6)[[1]][13])
ss_m = fixef(blength_mlselec_m6)[[1]][5] + fixef(blength_mlselec_m6)[[1]][10] + 
  density * (fixef(blength_mlselec_m6)[[1]][12] + fixef(blength_mlselec_m6)[[1]][14])
cij_m = fixef(assort_m1)[[1]][3] +  fixef(assort_m1)[[1]][5] +
  density * (fixef(assort_m1)[[1]][7] + fixef(assort_m1)[[1]][8])
var_m = var(weevil_comb_complete$body_length[weevil_comb_complete$sex == "M"])
S_m = var_m*(ds_m+ss_m*cij_m) #increase with density


real_density = density * sd(weevil_comb_complete$density) + mean(weevil_comb_complete$density)

#single panel
par(mfrow=c(1,3), mar = c(5,5,4,2))

plot(S_f ~ real_density, 
     ylim = c(0,60),type="l", 
     xlab = "",
     ylab = "Total selection differential",
     cex.lab = 1.5,
     frame = F, las=1, main="a.")
box(which = "plot", bty = "l")
points(S_m ~ real_density, lty = 2, type="l")

plot(ds_f*var_f ~ real_density, 
     ylim = c(0,60),type="l", 
     xlab = "Density (weevils per tree)",
     ylab = "Direct selection differential",
     cex.lab = 1.5,
     frame = F, las=1, main="b.")
box(which = "plot", bty = "l")
points(ds_m*var_m ~ real_density, lty = 2, type="l")

plot(ss_f*cij_f*var_f ~ real_density, 
     ylim = c(-1.5,60),type="l", 
     xlab = "",
     ylab = "Social selection differential",
     cex.lab = 1.5,
     frame = F, las=1, main="c.")
box(which = "plot", bty = "l")
points(ss_m*cij_m*var_m ~ real_density, lty = 2, type="l")
legend("topleft",
       legend=c("Females", "Males"),
       cex = 1.5,
       col="black",text.col="black",
       lty=c(1,2),
       #bty="n",
       inset = c(0.05,0.05))

#END####

