

# packages
library(ggplot2)
library(plyr)


# Lets try a two trait model
# for negative, random and positive assortment

# (A)  matching selection 
# (i.e. the most competitive trait in intraspecific competition is the most competitive trait when competing with heterospecifics.
#
#
# (1) with 2000 intraspecific competitors at a 50:50 ratio
#     with zero interspecific competitors
# (2) with 1500 intraspecific competitors at a 50:50 ratio
#     with 500 random interspecific competitors (i.e. 25% inter)
# (3) with 1500 intraspecific competitors at a 50:50 ratio
#     with 500 random intraspecific competitors removed (i.e. 25%)


# (B)  contrasting selection 
# (i.e. the most competitive trait in intraspecific competition is the least competitive trait when competing with heterospecifics.
#
#
# (1) with 1500 intraspecific competitors at a 50:50 ratio
#     with 500 random interspecific competitors (i.e. 25% inter)


# (C)
#
#
# plot.

########################
# (A1) All Intraspecific 
########################


# negative assortment
# set seed
set.seed(1986)
# results
neg_intra <- rep(NA, 1000)
for(j in 1:length(neg_intra)){
# create 1000 focal individuals, 50:50 with traits 0 or 1
focal <-  rep(0:1,each = 500)
# create 1000 competitor individuals,  50:50 with traits 0 or 1
comp <-  rep(1:0, each = 500)
# pair them in dyads with negative assortment
negative1 <-  data.frame(focal = focal, comp = comp)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
negative2 <- data.frame(focal = comp, comp = focal)
negative <- rbind(negative1, negative2)
# calculate fintess (W) based on traits
negative$W <- ifelse(negative$focal > negative$comp, 2,
                ifelse(negative$focal == negative$comp, 1, 0) )
# calculate relative fitness
w <- negative$W/mean(negative$W)
# mean standardised trait
z <- negative$focal/mean(negative$focal)
# mean standardised selection gradient
neg_intra[[j]] <- coef(lm(w ~ z))[[2]]
}


# random
# set seed
set.seed(1986)
# results
ran_intra <- rep(NA, 1000)
for(j in 1:length(ran_intra)){
# create 1000 focal individuals, randomly with traits 0 or 1
focal <-  sample(0:1, 1000, replace = TRUE)
# pair them with the a competitor, ensuring 50:50 traits 0 or 2
# so each each dyad is represented once
comp <- sample(c(rep(0, 1000-sum(focal == 0)), rep(1, 1000-sum(focal == 1)) ))
random1 <-  data.frame(focal = focal, comp = comp)
# calculate fintess (W) based on traits
random1$W <- ifelse(random1$focal > random1$comp, 2,
                ifelse(random1$focal == random1$comp, 1, 0) )
# create a data.frame with the other 1000 competitors as focals
random2 <- data.frame(focal = comp, comp = focal)
# and calculate fintess based on previous outcome
random2$W <- ifelse(random1$W == 2, 0, ifelse(random1$W == 0, 2, 1) )
# bind them together so we have each pair represented twice
# with each individual being the focal
random <- rbind(random1, random2)
# calculate relative fitness
w <- random$W/mean(random$W)
# mean standardised trait
z <- random$focal/mean(random$focal)
# mean standardised selection gradient
ran_intra[[j]] <- coef(lm(w ~ z))[[2]]
}


# positive
# set seed
set.seed(1986)
# results
pos_intra <- rep(NA, 1000)
for(j in 1:length(pos_intra)){
# create 1000 focal individuals, 50:50 with traits 0 or 1
focal <-  rep(0:1,each = 500)
# create 1000 competitor indidivuals,  50:50 with traits 0 or 1
comp <-  rep(0:1, each = 500)
# pair them in dyads with positive assortment
positive1 <-  data.frame(focal = focal, comp = comp)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
positive2 <- data.frame(focal = comp, comp = focal)
positive <- rbind(positive1, positive2)
# calculate fintess (W) based on traits
positive$W <- ifelse(positive$focal > positive$comp, 2,
                ifelse(positive$focal == positive$comp, 1, 0) )
# calculate relative fitness
w <- positive$W/mean(positive$W)
# mean standardised trait
z <- positive$focal/mean(positive$focal)
# mean standardised selection gradient
pos_intra[[j]] <- coef(lm(w ~ z))[[2]]
}



# have a look
par(mfrow=c(1,3))
plot(neg_intra, ylim = c(-2,2))
plot(ran_intra, ylim = c(-2,2))
plot(pos_intra, ylim = c(-2,2))



################################################################################
# (A2) With interspecific
################################################################################



# negative
# set seed
set.seed(1986)
# results
neg_inter <- rep(NA, 1000)
for(j in 1:length(neg_inter)){
# create 1000 focal individuals, 50:50 with traits 0 or 1
focal <-  rep(0:1,each = 500)
# create 1000 competitor individuals,  50:50 with traits 0 or 1
comp <-  rep(1:0, each = 500)
# pair them in dyads with negative assortment
negative1 <-  data.frame(focal = focal, comp = comp)
# sample 500 dyads where a focal or competitor will be swapped
dyads <- sample(1:1000, 500)
# create a vector to hold the which individual to swap
swaps <- rep(NA, 500)
for(i in 1:length(dyads)){
# 1 = column 1 which is focal
# 2 = column 2 which is competitor
swaps[i] <- sample( c(1, 2), 1 )
}
# now select new interspecific traits also from a population with 50:50 0 or 1
inter_traits <- sample(rep(0:1,each = 1000), 500)
# now change the traits of the swapped individuals
for(k in 1:length(dyads)){
negative1[, swaps[k] ][ dyads[k] ] <- inter_traits[k]
}
# note which individual was swapped
negative1$swapped_focal <- NA
negative1$swapped_focal[dyads] <- ifelse(swaps == 1, "focal_swapped", NA)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
negative2 <- data.frame(focal = negative1$comp, 
                        comp = negative1$focal
                        )
# note which individual was swapped
negative2$swapped_focal <- NA
negative2$swapped_focal[dyads] <- ifelse(swaps == 2, "focal_swapped", NA)                                   
# add together so each pair represented twice, each individual being the focal
negative <- rbind(negative1, negative2)
# calculate fintess (W) based on traits
negative$W <- ifelse(negative$focal > negative$comp, 2,
                ifelse(negative$focal == negative$comp, 1, 0) )
# subset only focal individuals that are focal species
negative <- negative[ is.na(negative$swapped_focal) ,]
# calculate relative fitness
w <- negative$W/mean(negative$W)
# mean standardised trait
z <- negative$focal/mean(negative$focal)
# mean standardised selection gradient
neg_inter[[j]] <- coef(lm(w ~ z))[[2]]
}



# random
# set seed
set.seed(1986)
# results
ran_inter <- rep(NA, 1000)
for(j in 1:length(ran_inter)){
# create 1000 focal individuals, randomly with traits 0 or 1
focal <-  sample(0:1, 1000, replace = TRUE)
# pair them with the a competitor, ensuring 50:50 traits 0 or 2
# so each each dyad is represented once
comp <- sample(c(rep(0, 1000-sum(focal == 0)), rep(1, 1000-sum(focal == 1)) ) )
random1 <-  data.frame(focal = focal, comp = comp)
# sample 500 dyads where a focal or competitor will be swapped
dyads <- sample(1:1000, 500)
# create a vector to hold the which individual to swap
swaps <- rep(NA, 500)
for(i in 1:length(dyads)){
# 1 = column 1 which is focal
# 2 = column 2 which is competitor
swaps[i] <- sample( c(1, 2), 1 )
}
# now select new interspecific traits also from a population with 50:50  0 or 1
inter_traits <- sample(rep(0:1,each = 1000), 500)
# now change the traits of the swapped individuals
for(k in 1:length(dyads)){
random1[, swaps[k] ][ dyads[k] ] <- inter_traits[k]
}
# note which individual was swapped
random1$swapped_focal <- NA
random1$swapped_focal[dyads] <- ifelse(swaps == 1, "focal_swapped", NA)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
random2 <- data.frame(focal = random1$comp, 
                        comp = random1$focal
                        )
# note which individual was swapped
random2$swapped_focal <- NA
random2$swapped_focal[dyads] <- ifelse(swaps == 2, "focal_swapped", NA)                                   
# add together so each pair represented twice, each individual being the focal
random <- rbind(random1, random2)
# calculate fintess (W) based on traits
random$W <- ifelse(random$focal > random$comp, 2,
                ifelse(random$focal == random$comp, 1, 0) )
# subset only focal individuals that are focal species
random <- random[ is.na(random$swapped_focal) ,]
# calculate relative fitness
w <- random$W/mean(random$W)
# mean standardised trait
z <- random$focal/mean(random$focal)
# mean standardised selection gradient
ran_inter[[j]] <- coef(lm(w ~ z))[[2]]
}



# positive
# set seed
set.seed(1986)
# results
pos_inter <- rep(NA, 1000)
for(j in 1:length(pos_inter)){
# create 1000 focal individuals, 50:50 with traits 0 or 1
focal <-  rep(0:1,each = 500)
# create 1000 competitor indidivuals,  50:50 with traits 0 or 1
comp <-  rep(0:1, each = 500)
# pair them in dyads with positive assortment
positive1 <-  data.frame(focal = focal, comp = comp)
# sample 500 dyads where a focal or competitor will be swapped
dyads <- sample(1:1000, 500)
# create a vector to hold the which individual to swap
swaps <- rep(NA, 500)
for(i in 1:length(dyads)){
# 1 = column 1 which is focal
# 2 = column 2 which is competitor
swaps[i] <- sample( c(1, 2), 1 )
}
# now select new interspecific traits also from a population with 50:50  0 or 1
inter_traits <- sample(rep(0:1,each = 1000), 500)
# now change the traits of the swapped individuals
for(k in 1:length(dyads)){
positive1[, swaps[k] ][ dyads[k] ] <- inter_traits[k]
}
# note which individual was swapped
positive1$swapped_focal <- NA
positive1$swapped_focal[dyads] <- ifelse(swaps == 1, "focal_swapped", NA)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
positive2 <- data.frame(focal = positive1$comp, 
                        comp = positive1$focal
                        )
# note which individual was swapped
positive2$swapped_focal <- NA
positive2$swapped_focal[dyads] <- ifelse(swaps == 2, "focal_swapped", NA)                                   
# add together so each pair represented twice, each individual being the focal
positive <- rbind(positive1, positive2)
# calculate fintess (W) based on traits
positive$W <- ifelse(positive$focal > positive$comp, 2,
                ifelse(positive$focal == positive$comp, 1, 0) )
# subset only focal individuals that are focal species
positive <- positive[ is.na(positive$swapped_focal) ,]
# calculate relative fitness
w <- positive$W/mean(positive$W)
# mean standardised trait
z <- positive$focal/mean(positive$focal)
# mean standardised selection gradient
pos_inter[[j]] <- coef(lm(w ~ z))[[2]]
}


# have a look
par(mfrow=c(1,3))
plot(neg_inter, ylim = c(-2,2))
plot(ran_inter, ylim = c(-2,2))
plot(pos_inter, ylim = c(-2,2))




############################
# (A3) Intraspecific control
#############################

# do everything as for inter
# i.e. pick random individuals
# mark them as interspecifics
# exclude them from analyses
# only thing that is different is that we do not change the traits of individuals
# this allows us to account for any differences in selection gradients that may occur from random removal of 500 individuals in interspecific


# negative
# set seed
set.seed(1986)
# results
neg_intra_control <- rep(NA, 1000)
for(j in 1:length(neg_intra_control)){
# create 1000 focal individuals, 50:50 with traits 0 or 1
focal <-  rep(0:1,each = 500)
# create 1000 competitor individuals,  50:50 with traits 0 or 1
comp <-  rep(1:0, each = 500)
# pair them in dyads with negative assortment
negative1 <-  data.frame(focal = focal, comp = comp)
# sample 500 dyads where a focal or competitor will be swapped
dyads <- sample(1:1000, 500)
# create a vector to hold the which individual to swap
swaps <- rep(NA, 500)
for(i in 1:length(dyads)){
# 1 = column 1 which is focal
# 2 = column 2 which is competitor
swaps[i] <- sample( c(1, 2), 1 )
}
# 
#
# 
#
#
#
# note which individual was swapped
negative1$swapped_focal <- NA
negative1$swapped_focal[dyads] <- ifelse(swaps == 1, "focal_swapped", NA)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
negative2 <- data.frame(focal = negative1$comp, 
                        comp = negative1$focal
                        )
# note which individual was swapped
negative2$swapped_focal <- NA
negative2$swapped_focal[dyads] <- ifelse(swaps == 2, "focal_swapped", NA)                                   
# add together so each pair represented twice, each individual being the focal
negative <- rbind(negative1, negative2)
# calculate fintess (W) based on traits
negative$W <- ifelse(negative$focal > negative$comp, 2,
                ifelse(negative$focal == negative$comp, 1, 0) )
# subset only focal individuals that are focal species
negative <- negative[ is.na(negative$swapped_focal) ,]
# calculate relative fitness
w <- negative$W/mean(negative$W)
# mean standardised trait
z <- negative$focal/mean(negative$focal)
# mean standardised selection gradient
neg_intra_control[[j]] <- coef(lm(w ~ z))[[2]]
}




# random
# set seed
set.seed(1986)
# results
ran_intra_control <- rep(NA, 1000)
for(j in 1:length(ran_intra_control)){
# create 1000 focal individuals, randomly with traits 0 or 1
focal <-  sample(0:1, 1000, replace = TRUE)
# pair them with the a competitor, ensuring 50:50 traits 0 or 2
# so each each dyad is represented once
comp <- sample(c(rep(0, 1000-sum(focal == 0)), rep(1, 1000-sum(focal == 1)) ) )
random1 <-  data.frame(focal = focal, comp = comp)
# sample 500 dyads where a focal or competitor will be swapped
dyads <- sample(1:1000, 500)
# create a vector to hold the which individual to swap
swaps <- rep(NA, 500)
for(i in 1:length(dyads)){
# 1 = column 1 which is focal
# 2 = column 2 which is competitor
swaps[i] <- sample( c(1, 2), 1 )
}
# 
#
# 
#
#
#
# note which individual was swapped
random1$swapped_focal <- NA
random1$swapped_focal[dyads] <- ifelse(swaps == 1, "focal_swapped", NA)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
random2 <- data.frame(focal = random1$comp, 
                        comp = random1$focal
                        )
# note which individual was swapped
random2$swapped_focal <- NA
random2$swapped_focal[dyads] <- ifelse(swaps == 2, "focal_swapped", NA)                                   
# add together so each pair represented twice, each individual being the focal
random <- rbind(random1, random2)
# calculate fintess (W) based on traits
random$W <- ifelse(random$focal > random$comp, 2,
                ifelse(random$focal == random$comp, 1, 0) )
# subset only focal individuals that are focal species
random <- random[ is.na(random$swapped_focal) ,]
# calculate relative fitness
w <- random$W/mean(random$W)
# mean standardised trait
z <- random$focal/mean(random$focal)
# mean standardised selection gradient
ran_intra_control[[j]] <- coef(lm(w ~ z))[[2]]
}




# positive
# set seed
set.seed(1986)
# results
pos_intra_control <- rep(NA, 1000)
for(j in 1:length(pos_intra_control)){
# create 1000 focal individuals, 50:50 with traits 0 or 1
focal <-  rep(0:1,each = 500)
# create 1000 competitor indidivuals,  50:50 with traits 0 or 1
comp <-  rep(0:1, each = 500)
# pair them in dyads with positive assortment
positive1 <-  data.frame(focal = focal, comp = comp)
# sample 500 dyads where a focal or competitor will be swapped
dyads <- sample(1:1000, 500)
# create a vector to hold the which individual to swap
swaps <- rep(NA, 500)
for(i in 1:length(dyads)){
# 1 = column 1 which is focal
# 2 = column 2 which is competitor
swaps[i] <- sample( c(1, 2), 1 )
}
# 
#
# 
#
#
#
# note which individual was swapped
positive1$swapped_focal <- NA
positive1$swapped_focal[dyads] <- ifelse(swaps == 1, "focal_swapped", NA)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
positive2 <- data.frame(focal = positive1$comp, 
                        comp = positive1$focal
                        )
# note which individual was swapped
positive2$swapped_focal <- NA
positive2$swapped_focal[dyads] <- ifelse(swaps == 2, "focal_swapped", NA)                                   
# add together so each pair represented twice, each individual being the focal
positive <- rbind(positive1, positive2)
# calculate fintess (W) based on traits
positive$W <- ifelse(positive$focal > positive$comp, 2,
                ifelse(positive$focal == positive$comp, 1, 0) )
# subset only focal individuals that are focal species
positive <- positive[ is.na(positive$swapped_focal) ,]
# calculate relative fitness
w <- positive$W/mean(positive$W)
# mean standardised trait
z <- positive$focal/mean(positive$focal)
# mean standardised selection gradient
pos_intra_control[[j]] <- coef(lm(w ~ z))[[2]]
}



# have a look
par(mfrow=c(1,3))
plot(neg_intra_control, ylim = c(-2,2))
plot(ran_intra_control, ylim = c(-2,2))
plot(pos_intra_control, ylim = c(-2,2))


######################################
# (B1) contrasting with interspecific 
######################################


# negative
# set seed
set.seed(1986)
# results
neg_inter_opp <- rep(NA, 1000)
for(j in 1:length(neg_inter_opp)){
# create 1000 focal individuals, 50:50 with traits 0 or 1
focal <-  rep(0:1,each = 500)
# create 1000 competitor individuals,  50:50 with traits 0 or 1
comp <-  rep(1:0, each = 500)
# pair them in dyads with negative assortment
negative1 <-  data.frame(focal = focal, comp = comp)
# sample 500 dyads where a focal or competitor will be swapped
dyads <- sample(1:1000, 500)
# create a vector to hold the which individual to swap
swaps <- rep(NA, 500)
for(i in 1:length(dyads)){
# 1 = column 1 which is focal
# 2 = column 2 which is competitor
swaps[i] <- sample( c(1, 2), 1 )
}
# now select new interspecific traits also from a population with 50:50  0 or 1
inter_traits <- sample(rep(0:1,each = 1000), 500)
# now change the traits of the swapped individuals
for(k in 1:length(dyads)){
negative1[, swaps[k] ][ dyads[k] ] <- inter_traits[k]
}
# note which individual was swapped
negative1$swapped_focal <- NA
negative1$swapped_focal[dyads] <- ifelse(swaps == 1, "focal_swapped", NA)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
negative2 <- data.frame(focal = negative1$comp, 
                        comp = negative1$focal
                        )
# note which individual was swapped
negative2$swapped_focal <- NA
negative2$swapped_focal[dyads] <- ifelse(swaps == 2, "focal_swapped", NA)                                   
#
#
#
# calculate fintess (W) based on traits for all dyads in one direction
# >>> for intraspecific competition
negative1$W <- ifelse(negative1$focal > negative1$comp, 2,
                ifelse(negative1$focal == negative1$comp, 1, 0) )
# >>> correct this for interspecific dyads
negative1$W[dyads] <- 
  ifelse(negative1$focal[dyads]  < negative1$comp[dyads], 2,
    ifelse(negative1$focal[dyads]  == negative1$comp[dyads] , 1, 0) )
# replicate this fitness for all dyads in other direction
# >>> for intraspecific competition
negative2$W <- ifelse(negative2$focal > negative2$comp, 2,
                ifelse(negative2$focal == negative2$comp, 1, 0) )
# >>> correct this for interspecific dyads
negative2$W[dyads] <- 
  ifelse(negative2$focal[dyads]  < negative2$comp[dyads], 2,
    ifelse(negative2$focal[dyads]  == negative2$comp[dyads] , 1, 0) )
# bins together
negative <- rbind(negative1, negative2)
#
#
# subset only focal individuals that are focal species
negative <- negative[ is.na(negative$swapped_focal) ,]
# calculate relative fitness
w <- negative$W/mean(negative$W)
# mean standardised trait
z <- negative$focal/mean(negative$focal)
# mean standardised selection gradient
neg_inter_opp[[j]] <- coef(lm(w ~ z))[[2]]
}


# random
# set seed
set.seed(1986)
# results
ran_inter_opp <- rep(NA, 1000)
for(j in 1:length(ran_inter_opp)){
# create 1000 focal individuals, randomly with traits 0 or 1
focal <-  sample(0:1, 1000, replace = TRUE)
# pair them with the a competitor, ensuring 50:50 traits 0 or 2
# so each each dyad is represented once
comp <- sample(c(rep(0, 1000-sum(focal == 0)), rep(1, 1000-sum(focal == 1)) ) )
random1 <-  data.frame(focal = focal, comp = comp)
# sample 500 dyads where a focal or competitor will be swapped
dyads <- sample(1:1000, 500)
# create a vector to hold the which individual to swap
swaps <- rep(NA, 500)
for(i in 1:length(dyads)){
# 1 = column 1 which is focal
# 2 = column 2 which is competitor
swaps[i] <- sample( c(1, 2), 1 )
}
# now select new interspecific traits also from a population with 50:50  0 or 1
inter_traits <- sample(rep(0:1,each = 1000), 500)
# now change the traits of the swapped individuals
for(k in 1:length(dyads)){
random1[, swaps[k] ][ dyads[k] ] <- inter_traits[k]
}
# note which individual was swapped
random1$swapped_focal <- NA
random1$swapped_focal[dyads] <- ifelse(swaps == 1, "focal_swapped", NA)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
random2 <- data.frame(focal = random1$comp, 
                        comp = random1$focal
                        )
# note which individual was swapped
random2$swapped_focal <- NA
random2$swapped_focal[dyads] <- ifelse(swaps == 2, "focal_swapped", NA)                                   
#
#
#
# calculate fintess (W) based on traits for all dyads in one direction
# >>> for intraspecific competition
random1$W <- ifelse(random1$focal > random1$comp, 2,
                ifelse(random1$focal == random1$comp, 1, 0) )
# >>> correct this for interspecific dyads
random1$W[dyads] <- 
  ifelse(random1$focal[dyads]  < random1$comp[dyads], 2,
    ifelse(random1$focal[dyads]  == random1$comp[dyads] , 1, 0) )
# replicate this fitness for all dyads in other direction
# >>> for intraspecific competition
random2$W <- ifelse(random2$focal > random2$comp, 2,
                ifelse(random2$focal == random2$comp, 1, 0) )
# >>> correct this for interspecific dyads
random2$W[dyads] <- 
  ifelse(random2$focal[dyads]  < random2$comp[dyads], 2,
    ifelse(random2$focal[dyads]  == random2$comp[dyads] , 1, 0) )
# bins together
random <- rbind(random1, random2)
#
#
# subset only focal individuals that are focal species
# subset only focal individuals that are focal species
random <- random[ is.na(random$swapped_focal) ,]
# calculate relative fitness
w <- random$W/mean(random$W)
# mean standardised trait
z <- random$focal/mean(random$focal)
# mean standardised selection gradient
ran_inter_opp[[j]] <- coef(lm(w ~ z))[[2]]
}



# positive
# set seed
set.seed(1986)
# results
pos_inter_opp <- rep(NA, 1000)
for(j in 1:length(pos_inter_opp)){
# create 1000 focal individuals, 50:50 with traits 0 or 1
focal <-  rep(0:1,each = 500)
# create 1000 competitor indidivuals,  50:50 with traits 0 or 1
comp <-  rep(0:1, each = 500)
# pair them in dyads with positive assortment
positive1 <-  data.frame(focal = focal, comp = comp)
# sample 500 dyads where a focal or competitor will be swapped
dyads <- sample(1:1000, 500)
# create a vector to hold the which individual to swap
swaps <- rep(NA, 500)
for(i in 1:length(dyads)){
# 1 = column 1 which is focal
# 2 = column 2 which is competitor
swaps[i] <- sample( c(1, 2), 1 )
}
# now select new interspecific traits also from a population with 50:50  0 or 1
inter_traits <- sample(rep(0:1,each = 1000), 500)
# now change the traits of the swapped individuals
for(k in 1:length(dyads)){
positive1[, swaps[k] ][ dyads[k] ] <- inter_traits[k]
}
# note which individual was swapped
positive1$swapped_focal <- NA
positive1$swapped_focal[dyads] <- ifelse(swaps == 1, "focal_swapped", NA)
# duplicate dyads so we have each pair represented twice
# with each individual being the focal
positive2 <- data.frame(focal = positive1$comp, 
                        comp = positive1$focal
                        )
# note which individual was swapped
positive2$swapped_focal <- NA
positive2$swapped_focal[dyads] <- ifelse(swaps == 2, "focal_swapped", NA)                                   
#
#
#
# calculate fintess (W) based on traits for all dyads in one direction
# >>> for intraspecific competition
positive1$W <- ifelse(positive1$focal > positive1$comp, 2,
                ifelse(positive1$focal == positive1$comp, 1, 0) )
# >>> correct this for interspecific dyads
positive1$W[dyads] <- 
  ifelse(positive1$focal[dyads]  < positive1$comp[dyads], 2,
    ifelse(positive1$focal[dyads]  == positive1$comp[dyads] , 1, 0) )
# replicate this fitness for all dyads in other direction
# >>> for intraspecific competition
positive2$W <- ifelse(positive2$focal > positive2$comp, 2,
                ifelse(positive2$focal == positive2$comp, 1, 0) )
# >>> correct this for interspecific dyads
positive2$W[dyads] <- 
  ifelse(positive2$focal[dyads]  < positive2$comp[dyads], 2,
    ifelse(positive2$focal[dyads]  == positive2$comp[dyads] , 1, 0) )
# bins together
positive <- rbind(positive1, positive2)
#
#
# subset only focal individuals that are focal species
positive <- positive[ is.na(positive$swapped_focal) ,]
# calculate relative fitness
w <- positive$W/mean(positive$W)
# mean standardised trait
z <- positive$focal/mean(positive$focal)
# mean standardised selection gradient
pos_inter_opp[[j]] <- coef(lm(w ~ z))[[2]]
}


# have a look
par(mfrow=c(1,3))
plot(neg_inter_opp, ylim = c(-2,2))
plot(ran_inter_opp, ylim = c(-2,2))
plot(pos_inter_opp, ylim = c(-2,2))



########################
# (C) Plot
########################

# quick look
par(mfrow=c(2,2))
# intra
plot(neg_intra, ylim = c(-1.5, 1.5), col = "purple")
points(ran_intra, ylim = c(-1.5, 1.5))
points(pos_intra, ylim = c(-1.5, 1.5), col = "green")
text(500,1.4,"intraspecfic only, matching")
# inter
plot(neg_inter, ylim = c(-1.5, 1.5), col = "purple")
points(ran_inter, ylim = c(-1.5, 1.5))
points(pos_inter, ylim = c(-1.5, 1.5), col = "green")
text(500,1.4,"interspecfic, matching")
# intra control
plot(neg_intra_control, ylim = c(-1.5, 1.5), col = "purple")
points(ran_intra_control, ylim = c(-1.5, 1.5))
points(pos_intra_control, ylim = c(-1.5, 1.5), col = "green")
text(500,1.4,"intraspecfic control, matching")
# inter opp
plot(neg_inter_opp, ylim = c(-1.5, 1.5), col = "purple")
points(ran_inter_opp, ylim = c(-1.5, 1.5))
points(pos_inter_opp, ylim = c(-1.5, 1.5), col = "green")
text(500,1.4,"interspecfic contrasting")


# organise data for combined plot
# data
# intraspecific
intra <- data.frame(B = c(neg_intra, ran_intra, pos_intra),
                    r = rep(c("Negative", "Random", "Positive"), each = 1000),
                    sim = rep("Intraspecfic")
                    )
# interspecific matching
inter <- data.frame(B = c(neg_inter, ran_inter, pos_inter),
                    r = rep(c("Negative", "Random", "Positive"), each = 1000),
                    sim = rep("Interspecfic_matching")
                    )
# intraspecific control
intra_con <- data.frame(B = c(neg_intra_control, ran_intra_control,
                              pos_intra_control),
                        r = rep(c("Negative","Random","Positive"), each = 1000),
                        sim = rep("Intraspecfic control")
                        )
# interspecific contrasting
inter_opp <- data.frame(B = c(neg_inter_opp, ran_inter_opp, pos_inter_opp),
                    r = rep(c("Negative", "Random", "Positive"), each = 1000),
                    sim = rep("Interspecfic_contrasting")
                    )        

# combine all data
alldata <- rbind(intra,inter,intra_con,inter_opp)

# look
ggplot(alldata, aes(x = r,  y = B)) +
  geom_point() +
  facet_wrap(~ sim)

# look
ggplot(alldata, aes(x = r,  y = B)) +
  geom_boxplot(aes(colour = sim), position  = "dodge")

# gets mean and 95% range for all simulations
alldata_means <- ddply(alldata, .(sim, r), summarise,
                       means = mean(B),
                       U = quantile(B, 0.975),
                       L = quantile(B, 0.025)
                       )

# look at quantiles
quantiles_gg <- ggplot(alldata_means, aes(x = r,  y = means, group = sim)) +
                  geom_errorbar(aes(ymin = L, ymax = U), 
                                width = 0.3, 
                                position = position_dodge(width = 1)) +
                  geom_point(aes(colour = sim), 
                             position = position_dodge(width = 1))
# plot
print(quantiles_gg)

# change order for plot
alldata_means$r <- factor(alldata_means$r, 
                     levels = c("Positive","Random","Negative"))
alldata_means$sim <- factor(alldata_means$sim, 
                     levels = c("Intraspecfic",
                                "Intraspecfic control",
                                "Interspecfic_matching",
                                "Interspecfic_contrasting"))
# look again                    
quantiles_gg <- ggplot(alldata_means, aes(x = r,  y = means, group = sim)) +
                  geom_errorbar(aes(ymin = L, ymax = U), 
                                width = 0.8, 
                                position = position_dodge(width = 1)) +
                  geom_point(aes(fill = sim), shape = 21, size = 3,
                             position = position_dodge(width = 1)) +
                  scale_fill_manual("", values = c("white","grey",
                                                   "black","blue")) +
                  ylab("Selection gradient") + 
                  xlab("Assortment")
# plot
print(quantiles_gg)




## plot only the control as the final point
inter_same <- data.frame(B = c(neg_inter, ran_inter, pos_inter),
                    r = rep(c("Negative", "Random", "Positive"), each = 1000),
                    sim = rep("Intra- & Interspecific (matching)")
                    )
intra_con <- data.frame(B = c(neg_intra_control, ran_intra_control,
                              pos_intra_control),
                        r = rep(c("Negative","Random","Positive"), each = 1000),
                        sim = rep("Intraspecific")
                        )
inter_opp <- data.frame(B = c(neg_inter_opp, ran_inter_opp, pos_inter_opp),
                    r = rep(c("Negative", "Random", "Positive"), each = 1000),
                    sim = rep("Intra- & Interspecific (contrasting)")
                    )
# add together
inter_n_con <- rbind(inter_same,intra_con,inter_opp)
# gets 95% range 
inter_n_con_means <- ddply(inter_n_con, .(sim, r), summarise,
                       means = mean(B),
                       U = quantile(B, 0.975),
                       L = quantile(B, 0.025)
                       )
# order levels
inter_n_con_means$r <- factor(inter_n_con_means$r, 
                     levels = c("Positive","Random","Negative"))
inter_n_con_means$sim <- factor(inter_n_con_means$sim, 
                     levels = c("Intraspecific", 
                     "Intra- & Interspecific (matching)",
                     "Intra- & Interspecific (contrasting)"))
# plot                    
# plot                    
quantiles_gg <- ggplot(inter_n_con_means, aes(x = r,  y = means, group = sim)) +
                  geom_hline(yintercept = 0, linetype = "dashed", 
                             colour  = "grey") +
                  geom_errorbar(aes(ymin = L, ymax = U,colour = sim), 
                                width = 0.2) +
                  geom_line(aes(colour = sim), size = 0.5) +
                  geom_point(aes(colour = sim, shape = sim, size = sim), 
                                 fill = "white") +
                  scale_colour_manual("", values = c("#1f78b4ff",
                                                     "#000000",
                                                     "#000000")) +
                  scale_shape_manual("", values = c(20,21,24)) +
                  scale_size_manual("", values = c(1.5,2,2)) +
                  ylab(expression(paste("Selection gradient (", beta, ")"))) + 
                  xlab("Assortment") + 
                  theme(legend.position = "bottom") +
                  guides(colour=guide_legend(nrow=3,byrow=TRUE)) +
                  scale_y_continuous(breaks = round(seq(-1.2, 1.2, by = 0.1),1))

# look
print(quantiles_gg)


## stop. :)

#################









