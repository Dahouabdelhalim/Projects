

#  Data and script to accompany the manuscript: 
#  "Winter associations predict social and extra-pair mating patterns in a wild songbird" '


#  Load required packages
sapply(c('lme4','sna','data.table','dplyr'),
       function(x) suppressMessages(require (x ,character.only =TRUE) ) )

#  Load data 
#list incluidng all datasets 
#(dataset excluding social pairs, dataset including social pairs, the generated spatial null model and the raw data collected at the bird feeders)
load("all_data.RData")




#   1. Do winter associations predict spatial breeding arrangement?
####--------------------------------------------------------------------------------------------------

#  Select data 
df <- all_data[[1]]

# Create matrices
pro <- df[,c("ID1","ID2","rank")]   # proximity data
spa <- df[,c("ID1","ID2","total_activity")]   # spatial overlap data
wa <- df[,c("ID1","ID2","rank_asso_avg")]    # winter associations data

prox <- with(pro, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- rank
  out
})
diag(prox) <- NA

winter_asso <- with(wa, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- rank_asso_avg
  out
})
diag(winter_asso) <- NA

spatial_overlap <- with(spa, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- total_activity
  out
})
diag(spatial_overlap) <- NA


#   Linear matrix regression model
x<-rgraph(221,2)
x[1,,] <- winter_asso
x[2,,] <- spatial_overlap

lm <- netlm(prox, x, nullhyp = "classical", test.statistic = "beta", mode="graph")
summary(lm)
obs <- lm$coefficients[2]  # coefficient of observed data


# Permutation test

# load spatial null model, see main text for details
spatial_null <- all_data[[3]]
ma_co <- rep(0,1000)

for (i in c(1:1000)){
  b <- spatial_null[[i]]
  b$ID1 <- as.factor(b$ID1)
  b$ID2 <- as.factor(b$ID2)
  bb <- b[order(ID1),]
  bb <- bb[order(ID2),]
  h <- merge(pro, bb, by.x=c("ID1","ID2"),by.y=c("ID1","ID2"), all.x=TRUE)
  h <- unique(h)
  hh <- h[,c("ID1","ID2","ranked_asso")]
  b_matrix <- with(hh, {
    out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                  dimnames=list(levels(ID1), levels(ID2)))
    out[cbind(ID1, ID2)] <- ranked_asso
    out
  })
  diag(b_matrix) <- NA
  x<-rgraph(221,2)
  x[1,,] <- b_matrix
  x[2,,] <- spatial_overlap
  lm_random <- netlm(prox, x, nullhyp = "classical", test.statistic = "beta", mode="graph")
  ma_co[i] <- lm_random[[1:2]]  
  
}

hist(ma_co,breaks=100)
abline(v=obs, col="red")
sum(abs(obs) < abs(ma_co))/1000  # p-value

####--------------------------------------------------------------------------------------------------



#   2. Do winter associations predict future social pairs?
####--------------------------------------------------------------------------------------------------

# Select data
df <- all_data[[2]]

# Social pair matrix
social_pairs <- with(df, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- social_pair
  out
})

# Winter association matrix
winter_asso <- with(df, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- rank_asso_avg
  out
})

#  Logistic network regression model
kk_m <- netlogit(social_pairs, winter_asso, nullhyp = "classical",test.statistic = "beta") # tests based on classical asymptotics
summary(kk_m)
obs <- kk_m[[1:2]]  # to get the coefficient


#  Permutation
data_to_swap <- df[order(rank)]  
rand_coef_winter <- rep(0,1000)

for (i in 1:1000) {
  
  data_random <- data_to_swap[, random:=sample(rank_asso_avg),by=rank]
  d_rand <- data_random[,c("ID1","ID2","random")]
  m_rand <- with(d_rand, {
    out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                  dimnames=list(levels(ID1), levels(ID2)))
    out[cbind(ID1, ID2)] <- random
    out
  })
  
  lm_r <- netlogit(social_pairs, m_rand, nullhyp = "classical",test.statistic = "beta")
  rand_coef_winter[i] <- lm_r[[1:2]]
  
}

hist(rand_coef_winter,breaks=100)
abline(v=obs, col="red")
sum(abs(obs) < abs(rand_coef_winter))/1000   # p-value


####--------------------------------------------------------------------------------------------------



#   3. Do winter associations predict extra-pair paternity?
####--------------------------------------------------------------------------------------------------

#  Select data 
df <- all_data[[1]]
# If analyses should be repeated selecting only the 1st an 2nd direct neighbours
#df <- df[rank<3,]

epp <- df[,c("ID1","ID2","epp")]
distance <- df[,c("ID1","ID2","rank")]
age <- df[,c("ID1","ID2","age")]
box_visit <- df[,c("ID1","ID2","box_visit")]
winter_asso <- df[,c("ID1","ID2","rank_asso_avg")]
diff_arrival <- df[,c("ID1","ID2","diff_arrival")]

# Create matrices
epp1 <- with(epp, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- epp
  out
})

distance <- with(distance, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- rank
  out
})

age <- with(age, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- age
  out
})

box_visit <- with(box_visit, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- box_visit
  out
})

winter_asso <- with(winter_asso, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- rank_asso_avg
  out
})

diff_arrival <- with(diff_arrival, {
  out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                dimnames=list(levels(ID1), levels(ID2)))
  out[cbind(ID1, ID2)] <- diff_arrival
  out
})

x<-rgraph(221,5)
x[1,,] <- distance
x[2,,] <- age
x[3,,] <- box_visit
x[4,,] <- winter_asso
x[5,,] <- diff_arrival

lm <- netlogit(epp1, x, nullhyp = "classical", test.statistic = "beta", mode="graph")
summary(lm)
obs_coef_winter <- lm$coefficients[5]


# Permutation

data_to_swap <- df[order(rank)]  
rand_coef_winter <- rep(0,1000)

for (i in 1:1000) {
  
  data_random <- data_to_swap[, random:=sample(rank_asso_avg),by=rank]
  d_rand <- data_random[,c("ID1","ID2","random")]
  winter <- with(d_rand, {
    out <- matrix(nrow=nlevels(ID1), ncol=nlevels(ID2),
                  dimnames=list(levels(ID1), levels(ID2)))
    out[cbind(ID1, ID2)] <- random
    out
  })
  
  x<-rgraph(221,5)
  x[1,,] <- distance
  x[2,,] <- age
  x[3,,] <- box_visit
  x[4,,] <- winter
  x[5,,] <- diff_arrival
  
  lm_r <- netlogit(epp1, x, nullhyp = "classical", test.statistic = "beta", mode="graph")
  rand_coef_winter[i] <- lm_r$coefficients[5]
  
}

hist(rand_coef_winter, breaks=100)
abline(v=obs_coef_winter, col="red")
sum(abs(obs_coef_winter) < abs(rand_coef_winter))/1000  # p-value

####--------------------------------------------------------------------------------------------------



#   4. Comparing the association strength between social partners, extra-pair partners and other close neighbours
####--------------------------------------------------------------------------------------------------

# Select data
df <- all_data[[1]]

# EP partners
e <- df[df$epp==1,]
e$group <- 2
e <- e[,c("group","rank_asso_avg","ID1","ID2")]
e$ID1 <- as.character(e$ID1)
e$ID2 <- as.character(e$ID2)
e <- e[ID1<=ID2,]
e <- unique(e)

# Social pairs
soc_pair <- all_data[[2]]
breed <- soc_pair[soc_pair$social_pair==1,]
s <- breed[,c("rank_asso_avg","ID1","ID2")]
s$ID1 <- as.character(s$ID1)
s$ID2 <- as.character(s$ID2)
s <- s[ID1<=ID2,]  
s$group <- 1

# Direct neighbours
r1 <- df[df$rank==1,]
r1 <- r1[r1$epp==0,]
r1$group <- 3
r1 <- r1[,c("group","rank_asso_avg","ID1","ID2")]
r1$ID1 <- as.character(r1$ID1)
r1$ID2 <- as.character(r1$ID2)
r1 <- r1[ID1<=ID2,]

# 2nd direct neighbours
r2 <- df[df$rank==2,]
r2 <- r2[r2$epp==0,]
r2$group <- 4
r2 <- r2[,c("group","rank_asso_avg","ID1","ID2")]
r2$ID1 <- as.character(r2$ID1)
r2$ID2 <- as.character(r2$ID2)
r2 <- r2[ID1<=ID2,]

df_groups <- rbind(s,e,r1,r2)


# Permutation

fin <- list()

for (i in 1:1000) {
  
  data_random<- df_groups[sample(nrow(df_groups),nrow(df_groups),replace = TRUE)]
  df_groups$random <- data_random$group
  df_groups[,mean_r:= mean(rank_asso_avg), by=random]
  mean <- df_groups[,.N,by=.(random,mean_r)]
  fin[[i]] <- mean
  
}

dat <- lapply(fin, as.data.frame) 
dat <- do.call(rbind, dat)

g4 <- dat[dat$random==4,]
g1 <- dat[dat$random==1,]
g3 <- dat[dat$random==3,]
g2 <- dat[dat$random==2,]
all_groups <- cbind(g1=g1,g2=g2,g3=g3,g4=g4)

# all 1000 random means
all_groups$g1_g2 <- all_groups$g1.mean_r - all_groups$g2.mean_r
all_groups$g1_g3 <- all_groups$g1.mean_r - all_groups$g3.mean_r
all_groups$g2_g3 <- all_groups$g2.mean_r - all_groups$g3.mean_r
all_groups$g1_g4 <- all_groups$g1.mean_r - all_groups$g4.mean_r
all_groups$g2_g4 <- all_groups$g2.mean_r - all_groups$g4.mean_r
all_groups$g3_g4 <- all_groups$g3.mean_r - all_groups$g4.mean_r

# observed difference between means
df_groups[,mean:= mean(rank_asso_avg), by=group]
mean <- df_groups[,.N,by=.(group,mean)]
mean<- as.data.frame(mean)
mean1_2 <- mean[1,2] - mean[2,2]
mean1_3 <- mean[1,2] - mean[3,2]
mean2_3 <- mean[2,2] - mean[3,2]
mean1_4 <- mean[1,2] - mean[4,2]
mean2_4 <- mean[2,2] - mean[4,2]
mean3_4 <- mean[3,2] - mean[4,2]

# P values of pairwise comparison of mean association strength
sum(abs(mean1_2) < abs(all_groups$g1_g2))/1000  
sum(abs(mean1_3) < abs(all_groups$g1_g3))/1000  
sum(abs(mean2_3) < abs(all_groups$g2_g3))/1000  
sum(abs(mean1_4) < abs(all_groups$g1_g4))/1000  
sum(abs(mean2_4) < abs(all_groups$g2_g4))/1000  
sum(abs(mean3_4) < abs(all_groups$g3_g4))/1000  

####--------------------------------------------------------------------------------------------------



#   5. Temporal changes in the social network
####--------------------------------------------------------------------------------------------------

#  Performed the same way as analyses in section 2. and 3. but including networks generated from a restricted time period

####--------------------------------------------------------------------------------------------------
