library("tidyverse")
library("gplots") # for plotCI
library("scales") # for transparency
library("lme4") # for linear mixed-effect models
library("glmmTMB") # for zero-inflation and hurdle models
library("DHARMa") # for model fits
library("emmeans") # for model comparisons
library("mapdata") # for map data 
library("geosphere") # for distances between points
par(mar=c(5, 4, 4, 2), oma=c(2, 2, 2, 2))
options(scipen = 999)

# define a function that runs multiple model fitting tests 
test_model_fit <- function(model) {
  print(summary(model))
  sim_res <- simulateResiduals(model)
  testZeroInflation(sim_res)
  testOutliers(sim_res)
  testDispersion(sim_res) # simulation test for over/under dispersion
  plot(sim_res)
  hist(resid(model))
}


#############################################
# --- Make maps of the populations used --- #
#############################################

# load population data
mdata <- read.csv(file="data/Ipo_RI_survey_pops.csv", stringsAsFactors = TRUE, na.strings = c("NA", ""))

par(mar=c(1, 1, 1, 1), oma=c(1, 1, 1, 1))
map("worldHires","usa", xlim=c(-105,-75),ylim=c(15,40), col="gray90", fill=TRUE)
map("worldHires","mexico", xlim=c(-105,-75),ylim=c(15,40), col="gray95", fill=TRUE, add=TRUE)
points(mdata$long, mdata$lat, pch=19, cex=1.5, col=as.numeric(mdata$type))
text(mdata$long, mdata$lat+1, mdata$name)

map("worldHires", "usa", xlim=c(-86,-75), ylim=c(33,38), col="gray90", fill=TRUE)
points(mdata$long, mdata$lat, pch=19, cex=1.5, col=as.numeric(mdata$type))
text(mdata$long, mdata$lat + 0.25, mdata$name)

par(mar=c(5, 4, 4, 2), oma=c(2, 2, 2, 2))

##############################
# --- Analyze cross data --- #
##############################

# load data 
cross_data <- read.csv(file="data/Ipo_RI_survey.csv", stringsAsFactors = TRUE, na.strings = c("NA", ""))
cross_data$r.cross.type <- ordered(cross_data$r.cross.type, levels = c("ACxAC","ACxSC","SCxAC","SCxSC","ACxAL","ACxSL","ALxAC","SLxAC","SCxAL","SCxSL","ALxSC","SLxSC","ALxAL","SLxAL","ALxSL","SLxSL"))

# filter out the population "Ocean"
cross_data <- cross_data %>% filter(mom.loc != "Ocean", dad.loc != "Ocean") %>% droplevels()

# look at seed weights
hist(c(cross_data$seed1, cross_data$seed1, cross_data$seed1, cross_data$seed1), breaks = 30, xlab = "Seed weight (mg)", main = NULL)
abline(v = 10, col = "blue", lty = 2)

# make means data
means_dat <- cross_data %>% group_by(focal.plant, mom, mom.pop, mom.type, dad, dad.pop, dad.type, r.sp, cross.type, r.cross.type, cross.loc, intraspecific) %>%
  summarise(mean = mean(n.seeds)) %>%
  as.data.frame %>% filter(., complete.cases(.))
means_dat$r.cross.type <- ordered(means_dat$r.cross.type, levels = c("ACxAC","ACxSC","SCxAC","SCxSC","ACxAL","ACxSL","ALxAC","SLxAC","SCxAL","SCxSL","ALxSC","SLxSC","ALxAL","SLxAL","ALxSL","SLxSL"))

# look at population averages (use data with Ocean)
list1 <- cross_data %>% group_by(cross.ind) %>% 
  summarise(cross.loc = first(cross.loc), cross.type = first(cross.type), prop = mean(min1seed)) %>%
  group_by(cross.loc, cross.type) %>% summarise(prop.mean = mean(prop))
list2 <- means_dat %>% group_by(cross.loc, cross.type) %>% summarise(mean = mean(mean))
pop_averages <- full_join(list1, list2, by = "cross.loc")


#######################################################################
### --- Test for differences between different types of crosses --- ###
#######################################################################

# look at data
stripchart(mean ~ r.cross.type, data=means_dat, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per pod", 
           pch = 16, cex = 1.5, lty = 1, col = alpha(c(rep("#AD3C8C", 4), rep("#4259A7", 8), rep("#38B14A", 4)), 0.6))
mosaicplot(table(cross_data$r.cross.type, cross_data$n.seeds), las=2, main=NULL, ylab = "Number of seeds per pod", xlab = "Cross type")
mosaicplot(table(cross_data$r.sp, cross_data$n.seeds), las=2, main=NULL)

# test models of fruit set
full_B_model1 <- glmer(min1seed ~ r.sp + (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = cross_data)
full_B_model2 <- glmer(min1seed ~ r.sp + (1 | mom.pop/mom/cross.ind), family = "binomial", data = cross_data)
full_B_model3 <- glmer(min1seed ~ (1 | mom.pop/mom/cross.ind), family = "binomial", data = cross_data)
anova(full_B_model1, full_B_model2, full_B_model3)
test_model_fit(full_B_model2)
lsmeans(full_B_model2, pairwise ~ r.sp)
contrast(lsmeans(full_B_model2, "r.sp"), list(CvL=c(1,0,0,-1), CvH=c(1,-1,-1,0), LvH=c(0,-1,-1,1), H1vH2=c(0,1,-1,0)))
#test_model_fit(full_B_model1)
#anova(full_B_model1, glmer(min1seed ~ (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = cross_data))
#lsmeans(full_B_model1, pairwise ~ r.sp)
#contrast(lsmeans(full_B_model1, "r.sp"), list(CvL=c(1,0,0,-1), CvH=c(1,-1,-1,0), LvH=c(0,-1,-1,1), H1vH2=c(0,1,-1,0)))

# test models of mean seed number
full_M_model1 <- glmmTMB(mean ~ r.sp + (1 | mom.pop/mom) + (1 | dad.pop/dad),  ziformula=~r.sp, data = means_dat)
full_M_model2 <- glmmTMB(mean ~ r.sp + (1 | mom.pop/mom),  ziformula=~r.sp, data = means_dat)
full_M_model3 <- glmmTMB(mean ~ (1 | mom.pop/mom),  ziformula=~r.sp, data = means_dat)
anova(full_M_model1, full_M_model2, full_M_model3)
test_model_fit(full_M_model2)
lsmeans(full_M_model2, pairwise ~ r.sp)
contrast(lsmeans(full_M_model2, "r.sp"), list(CvL=c(1,0,0,-1), CvH=c(1,-1,-1,0), LvH=c(0,-1,-1,1), H1vH2=c(0,1,-1,0)))
#test_model_fit(full_M_model1)
#anova(full_M_model1, glmmTMB(mean ~ (1 | mom.pop/mom) + (1 | dad.pop/dad),  ziformula=~r.sp, data = means_dat))
#lsmeans(full_M_model1, pairwise ~ r.sp)
#contrast(lsmeans(full_M_model1, "r.sp"), list(CvL=c(1,0,0,-1), CvH=c(1,-1,-1,0), LvH=c(0,-1,-1,1), H1vH2=c(0,1,-1,0)))


#####################################################################################
### --- Test for differences between different types of interspecific crosses --- ###
#####################################################################################

# subset data
inter_means <- means_dat %>% filter(intraspecific == FALSE) %>% droplevels()
inter_cross <- cross_data %>% filter(intraspecific == FALSE) %>% droplevels()

# look at data
mosaicplot(table(inter_cross$r.cross.type, inter_cross$n.seeds), las=2, main=NULL)
t <- table(inter_cross$cross.type, inter_cross$min1seed) %>% as.data.frame.matrix()
names(t) <- c("fail", "success")
t %>% mutate(sum = fail + success, percent = fail/sum)
par(mfrow=c(2,2))
par(mar=c(2,2,2,2))
pie(c(t[1,1],t[1,2])); pie(c(t[3,1],t[3,2])); pie(c(t[2,1],t[2,2])); pie(c(t[4,1],t[4,2]))
par(mfrow=c(1,1))

# test models of fruit set
inter_B_model1 <- glmer(min1seed ~ cross.type + (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = inter_cross)
inter_B_model2 <- glmer(min1seed ~ cross.type + (1 | mom.pop/mom/cross.ind), family = "binomial", data = inter_cross)
inter_B_model3 <- glmer(min1seed ~ (1 | mom.pop/mom/cross.ind), family = "binomial", data = inter_cross)
anova(inter_B_model1, inter_B_model2, inter_B_model3)
test_model_fit(inter_B_model2)
lsmeans(inter_B_model2, pairwise ~ cross.type)
contrast(lsmeans(inter_B_model2, "cross.type"), list(AvS=c(1,0,0,-1), ACvSC=c(1,1,-1,-1), ALvSL=c(1,-1,1,-1)))
#test_model_fit(inter_B_model1)
#anova(inter_B_model1, glmer(min1seed ~ (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = inter_cross))
#lsmeans(inter_B_model1, pairwise ~ cross.type)
#contrast(lsmeans(inter_B_model1, "cross.type"), list(AvS=c(1,0,0,-1), ACvSC=c(1,1,-1,-1), ALvSL=c(1,-1,1,-1)))

# test models of mean seed number
inter_M_model1 <- glmmTMB(mean ~ cross.type + (1 | mom.pop/mom) + (1 | dad.pop/dad), ziformula=~1, family = tweedie(link = "log"), data = inter_means)
inter_M_model2 <- glmmTMB(mean ~ cross.type + (1 | mom.pop/mom), ziformula=~1, family = tweedie(link = "log"), data = inter_means)
inter_M_model3 <- glmmTMB(mean ~ (1 | mom.pop/mom), ziformula=~1, family = tweedie(link = "log"), data = inter_means)
anova(inter_M_model1, inter_M_model2, inter_M_model3)
test_model_fit(inter_M_model2)
lsmeans(inter_M_model2, pairwise ~ cross.type)
contrast(lsmeans(inter_M_model2, "cross.type"), list(AvS=c(1,0,0,-1), ACvSC=c(1,1,-1,-1), ALvSL=c(1,-1,1,-1)))
#test_model_fit(inter_M_model1)
#anova(inter_M_model1, glmmTMB(mean ~ (1 | mom.pop/mom) + (1 | dad.pop/dad), ziformula=~1, family = tweedie(link = "log"), data = inter_means))
#lsmeans(inter_M_model1, pairwise ~ cross.type)
#contrast(lsmeans(inter_M_model1, "cross.type"), list(AvS=c(1,0,0,-1), ACvSC=c(1,1,-1,-1), ALvSL=c(1,-1,1,-1)))


##################################################################################################
### find and plot mean and 95% CI's for each "type" of interspecific cross using bootstrapping ###
##################################################################################################

# make new categories
inter_means <- inter_means %>% mutate(A = (mom.type == "AC" & dad.type == "AL") | ((mom.type == "AL" & dad.type == "AC")),
                                      S = (mom.type == "SC" & dad.type == "SL") | ((mom.type == "SL" & dad.type == "SC")),
                                      AC = mom.type == "AC" | dad.type == "AC",
                                      SC = mom.type == "SC" | dad.type == "SC",
                                      AL = mom.type == "AL" | dad.type == "AL",
                                      SL = mom.type == "SL" | dad.type == "SL")

# find mean and 95% bootstrap confidence interval for each cross type
cross_type_summary <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("cross.type", "mean", "LCI", "UCI", "N"))
cross_type_list <- c("A", "S", "AC", "SC", "AL", "SL")
set.seed(1)
for (i in 1:6){
  cross_type_summary[i,1] <- cross_type_list[[i]]
  # subset data
  temp <- inter_means %>% filter(inter_means[,13+i] == TRUE)
  # find actual mean
  cross_type_summary[i,2] <- mean(temp$mean)
  # find 95% bootstrap confidence interval using the percentile method
  bstrap <- c()
  for (j in 1:10000){ bstrap <- c(bstrap, mean(sample(temp$mean,length(temp$mean),replace=T))) }
  cross_type_summary[i,3] <- quantile(bstrap,.025)
  cross_type_summary[i,4] <- quantile(bstrap,.975)
  cross_type_summary[i,5] <- nrow(temp)
}

# plot means and 95% confidence intervals for all comparisons
plotCI(c(1.1,1.9, 3.1,3.9,5.1,5.9), cross_type_summary$mean, ui = cross_type_summary$UCI, li = cross_type_summary$LCI,
       err="y", pch = 1, gap = 0, lwd = 2, cex = 1.5, sfrac = 0, las = 1, ylim = c(-0.02, 0.26), xlim = c(0.5, 6.5), xaxt = "n",
       ylab="Mean number of seeds per pod", xlab="Interspecific cross types")
axis(1, at = c(1.1,1.9,3.1,3.9,5.1,5.9), labels = cross_type_list)
segments(c(1.1, 3.1, 5.1), rep(0.23, 3), c(1.9, 3.9, 5.9), rep(0.23, 3))
text(1.5, 0.225, labels = "*", pos = 3)
text(3.5, 0.225, labels = "***", pos = 3)
text(5.5, 0.225, labels = "NS", pos = 3, cex = 0.75)
text(c(0.6,1.1,1.9,3.1,3.9,5.1,5.9), rep(0, 7), labels = c("N =", cross_type_summary$N), pos = 1, cex = 0.75)


######################################################################################
### --- Test for differences between different types of cordatotriloba crosses --- ###
######################################################################################

# subset data
cor_means <- means_dat %>% filter((mom.type == "AC" | mom.type == "SC") & (dad.type == "AC" | dad.type == "SC")) %>% droplevels()
cor_cross <- cross_data %>% filter((mom.type == "AC" | mom.type == "SC") & (dad.type == "AC" | dad.type == "SC")) %>% droplevels()

# look at data
stripchart(mean ~ cross.type, data=cor_means, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per pod", pch = 16, cex = 1, lty = 1, col = alpha("blue", 0.6))
plot(table(cor_cross$cross.type, cor_cross$n.seeds))
t2 <- table(cor_cross$cross.type, cor_cross$min1seed) %>% as.data.frame.matrix()
names(t2) <- c("fail", "success")
t2 %>% mutate(sum = fail + success, percent = fail/sum)
par(mfrow=c(1,3))
par(mar=c(2,2,2,2))
pie(c(t2[1,1],t2[1,2])); pie(c(t2[2,1],t2[2,2])); pie(c(t2[3,1],t2[3,2]))
par(mfrow=c(1,1))

# test models of fruit set
cor_B_model1 <- glmer(min1seed ~ cross.type + (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = cor_cross)
cor_B_model2 <- glmer(min1seed ~ cross.type + (1 | mom.pop/mom/cross.ind), family = "binomial", data = cor_cross)
cor_B_model3 <- glmer(min1seed ~ (1 | mom.pop/mom/cross.ind), family = "binomial", data = cor_cross)
anova(cor_B_model1, cor_B_model2, cor_B_model3)
test_model_fit(cor_B_model2)
lsmeans(cor_B_model2, pairwise ~ cross.type)
#test_model_fit(cor_B_model1)
#anova(cor_B_model1, glmer(min1seed ~ (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = cor_cross))
#lsmeans(cor_B_model1, pairwise ~ cross.type)

# test models of mean seed number
cor_M_model1 <- lmer(mean ~ cross.type + (1 | mom.pop/mom) + (1 | dad.pop/dad), data = cor_means)
cor_M_model2 <- lmer(mean ~ cross.type + (1 | mom.pop/mom), data = cor_means)
cor_M_model3 <- lmer(mean ~ (1 | mom.pop/mom), data = cor_means)
anova(cor_M_model1, cor_M_model2, cor_M_model3)
test_model_fit(cor_M_model2)
lsmeans(cor_M_model2, pairwise ~ cross.type)
#test_model_fit(cor_M_model1)
#anova(cor_M_model1, lmer(mean ~ (1 | mom.pop/mom) + (1 | dad.pop/dad), data = cor_means))
#lsmeans(cor_M_model1, pairwise ~ cross.type)


################################################################################
### --- Test for differences between different types of lacunosa crosses --- ###
################################################################################

# subset data
lac_means <- means_dat %>% filter((mom.type == "AL" | mom.type == "SL") & (dad.type == "AL" | dad.type == "SL")) %>% droplevels()
lac_cross <- cross_data %>% filter((mom.type == "AL" | mom.type == "SL") & (dad.type == "AL" | dad.type == "SL")) %>% droplevels()

# look at data
stripchart(mean ~ cross.type, data=lac_means, las = 2, vertical = TRUE, method = "jitter", ylab = "Mean number of seeds produced per pod", pch = 16, cex = 1, lty = 1, col = alpha("blue", 0.6))
plot(table(lac_cross$cross.type, lac_cross$n.seeds))
t3 <- table(lac_cross$cross.type, lac_cross$min1seed) %>% as.data.frame.matrix()
names(t3) <- c("fail", "success")
t3 %>% mutate(sum = fail + success, percent = fail/sum)
par(mfrow=c(1,3))
par(mar=c(2,2,2,2))
pie(c(27,91)); pie(c(67,176)); pie(c(35,86))
pie(c(t3[1,1],t3[1,2])); pie(c(t3[2,1],t3[2,2])); pie(c(t3[3,1],t3[3,2]))
par(mfrow=c(1,1))

# test models of fruit set
lac_B_model1 <- glmer(min1seed ~ cross.type + (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = lac_cross)
lac_B_model2 <- glmer(min1seed ~ cross.type + (1 | mom.pop/mom/cross.ind), family = "binomial", data = lac_cross)
lac_B_model3 <- glmer(min1seed ~ (1 | mom.pop/mom/cross.ind), family = "binomial", data = lac_cross)
anova(lac_B_model1, lac_B_model2, lac_B_model3)
test_model_fit(lac_B_model2)
lsmeans(lac_B_model2, pairwise ~ cross.type)
#test_model_fit(lac_B_model1)
#anova(lac_B_model1, glmer(min1seed ~ (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = lac_cross))
#lsmeans(lac_B_model1, pairwise ~ cross.type)

# test models of mean seed number
lac_M_model1 <- lmer(mean ~ cross.type + (1 | mom.pop/mom) + (1 | dad.pop/dad), data = lac_means)
lac_M_model2 <- lmer(mean ~ cross.type + (1 | mom.pop/mom), data = lac_means)
lac_M_model3 <- lmer(mean ~ (1 | mom.pop/mom), data = lac_means)
anova(lac_M_model1, lac_M_model2, lac_M_model3)
test_model_fit(lac_M_model2)
lsmeans(lac_M_model2, pairwise ~ cross.type)
#test_model_fit(lac_M_model1)
#anova(lac_M_model1, lmer(mean ~ (1 | mom.pop/mom) + (1 | dad.pop/dad), data = lac_means))
#lsmeans(lac_M_model1, pairwise ~ cross.type)


#############################################################
### --- Figure for both sets of intraspecific crosses --- ###
#############################################################

par(mfrow=c(1,2))
estimates1 <- lsmeans(cor_M_model2, pairwise ~ cross.type)[[1]] %>% data.frame()
stripchart(mean ~ cross.type, data=cor_means, las = 2, vertical = TRUE, method = "jitter", 
           ylab = "Mean number of seeds produced per pod", ylim = c(0, 4), xlim = c(0.5, 3.5),
           pch = 16, cex = 1.5, lty = 1, col = alpha("#AD3C8C", 0.6))
plotCI(c(1.2, 2.2, 3.2), estimates1$lsmean, ui=estimates1$upper.CL, li=estimates1$lower.CL, add=TRUE, 
       err="y", pch = 1, gap = 0, lwd = 2, cex = 1.5, sfrac = 0, las = 1)
estimates2 <- lsmeans(lac_M_model2, pairwise ~ cross.type)[[1]] %>% data.frame()
stripchart(mean ~ cross.type, data=lac_means, las = 2, vertical = TRUE, method = "jitter", 
           ylab = "Mean number of seeds produced per pod", ylim = c(0, 4), xlim = c(0.5, 3.5),
           pch = 16, cex = 1.5, lty = 1, col = alpha("#38B14A", 0.6))
plotCI(c(1.2, 2.2, 3.2), estimates2$lsmean, ui=estimates2$upper.CL, li=estimates2$lower.CL, add=TRUE, 
       err="y", pch = 1, gap = 0, lwd = 2, cex = 1.5, sfrac = 0, las = 1)
par(mfrow=c(1,1))


####################################################################
### --- Test for an effect of geographic distance on crosses --- ###
####################################################################

# load population locations 
pop_data <- read.csv(file="data/Ipo_RI_survey_pops.csv", stringsAsFactors = TRUE, na.strings = c("NA", "")) %>%
  select(name, lat, long)

# add locations to means data
intra_means <- means_dat %>% filter(intraspecific == TRUE) %>%
  mutate(mom_pop = gsub('.{2}$', '', .$mom.pop)) %>%
  left_join(., pop_data, by = c("mom_pop" = "name")) %>%
  mutate(dad_pop = gsub('.{2}$', '', .$dad.pop)) %>%
  left_join(., pop_data, by = c("dad_pop" = "name")) %>% 
  mutate(pch = case_when(cross.type == "ACxxAC" | cross.type == "ALxxAL" ~ 15,
                         cross.type == "ACxxSC" | cross.type == "ALxxSL" ~ 16,
                         TRUE ~ 17)) %>% droplevels()

# add locations to cross data
intra_crosses <- cross_data %>% filter(intraspecific == TRUE) %>%
  mutate(mom_pop = gsub('.{2}$', '', .$mom.pop)) %>%
  left_join(., pop_data, by = c("mom_pop" = "name")) %>%
  mutate(dad_pop = gsub('.{2}$', '', .$dad.pop)) %>%
  left_join(., pop_data, by = c("dad_pop" = "name")) %>% 
  mutate(pch = case_when(cross.type == "ACxxAC" | cross.type == "ALxxAL" ~ 15,
                         cross.type == "ACxxSC" | cross.type == "ALxxSL" ~ 16,
                         TRUE ~ 17)) %>% droplevels()

# calculate distances
for (i in 1:nrow(intra_means)) {
  intra_means[i,"dist"] <- distm(c(intra_means[i,"long.x"], intra_means[i,"lat.x"]), c(intra_means[i,"long.y"], intra_means[i,"lat.y"]))/1000
}
for (i in 1:nrow(intra_crosses)) {
  intra_crosses[i,"dist"] <- distm(c(intra_crosses[i,"long.x"], intra_crosses[i,"lat.x"]), c(intra_crosses[i,"long.y"], intra_crosses[i,"lat.y"]))/1000
}

# test for effect on fruit set in I cordatotriloba
intra_cross_C <- intra_crosses %>% filter(r.sp == "CxC")
cor_dist_B_model1 <- glmer(min1seed ~ dist + (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = intra_cross_C)
cor_dist_B_model2 <- glmer(min1seed ~ dist + (1 | mom.pop/mom/cross.ind), family = "binomial", data = intra_cross_C)
cor_dist_B_model3 <- glmer(min1seed ~ (1 | mom.pop/mom/cross.ind), family = "binomial", data = intra_cross_C)
anova(cor_dist_B_model1, cor_dist_B_model2, cor_dist_B_model3)
test_model_fit(cor_dist_B_model2)
#anova(cor_dist_B_model1, glmer(min1seed ~ (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = intra_cross_C))
#test_model_fit(cor_dist_B_model1)

# test for effect on mean seeds in I cordatotriloba
intra_means_C <- intra_means %>% filter(r.sp == "CxC")
cor_dist_M_model1 <- lmer(mean ~ dist + (1 | mom.pop/mom) + (1 | dad.pop/dad), data = intra_means_C)
cor_dist_M_model2 <- lmer(mean ~ dist + (1 | mom.pop/mom), data = intra_means_C)
cor_dist_M_model3 <- lmer(mean ~ (1 | mom.pop/mom), data = intra_means_C)
anova(cor_dist_M_model1, cor_dist_M_model2, cor_dist_M_model3)
test_model_fit(cor_dist_M_model2)
#anova(cor_dist_M_model1, lmer(mean ~ (1 | mom.pop/mom) + (1 | dad.pop/dad), data = intra_means_C))
#test_model_fit(cor_dist_M_model1)

# test for effect on fruit set in I lacunosa
intra_cross_L <- intra_crosses %>% filter(r.sp == "LxL")
lac_dist_B_model1 <- glmer(min1seed ~ dist + (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = intra_cross_L)
lac_dist_B_model2 <- glmer(min1seed ~ dist + (1 | mom.pop/mom/cross.ind), family = "binomial", data = intra_cross_L)
lac_dist_B_model3 <- glmer(min1seed ~ (1 | mom.pop/mom/cross.ind), family = "binomial", data = intra_cross_L)
anova(lac_dist_B_model1, lac_dist_B_model2, lac_dist_B_model3)
test_model_fit(lac_dist_B_model2)
#anova(lac_dist_B_model1, glmer(min1seed ~ (1 | mom.pop/mom/cross.ind) + (1 | dad.pop/dad), family = "binomial", data = intra_cross_L))
#test_model_fit(lac_dist_B_model1)

# test for effect on mean seeeds in I lacunosa
intra_means_L <- intra_means %>% filter(r.sp == "LxL")
lac_dist_M_model1 <- lmer(mean ~ dist + (1 | mom/mom.pop) + (1 | dad.pop/dad), data = intra_means_L)
lac_dist_M_model2 <- lmer(mean ~ dist + (1 | mom/mom.pop), data = intra_means_L)
lac_dist_M_model3 <- lmer(mean ~ (1 | mom/mom.pop), data = intra_means_L)
anova(lac_dist_M_model1, lac_dist_M_model2, lac_dist_M_model3)
test_model_fit(lac_dist_M_model2)
#anova(lac_dist_M_model1, lmer(mean ~ (1 | mom/mom.pop) + (1 | dad.pop/dad), data = intra_means_L))
#test_model_fit(lac_dist_M_model1)

# plot relationahip between number of seeds and distance
par(mfrow=c(2,1))
plot(intra_means_C$mean ~ (intra_means_C$dist), xlim = c(0, 2500), ylim = c(0,4), xlab = "Distance (km)", ylab = "Mean number of seeds",
     pch = intra_means_C$pch, cex = 1.5, lty = 1, col = alpha("#AD3C8C", 0.6))
abline(summary(cor_dist_M_model2)[[10]][1,1], summary(cor_dist_M_model2)[[10]][2,1])
plot(intra_means_L$mean ~ intra_means_L$dist, xlim = c(0, 2500), ylim = c(0,4), xlab = "Distance (km)", ylab = "Mean number of seeds",
     pch = intra_means_C$pch, cex = 1.5, lty = 1, col = alpha("#38B14A", 0.6))
abline(summary(lac_dist_M_model2)[[10]][1,1], summary(lac_dist_M_model2)[[10]][2,1])
legend(2000, 4, legend = c("allo", "cross", "sym"), pch = 15:17, col = alpha("black", 0.6))
par(mfrow=c(1,1))


#########################################
### --- Test for F1 cross success --- ###
#########################################

# load data 
F1_data <- read.csv(file="data/Ipo_RI_F1_crosses.csv", stringsAsFactors = TRUE, na.strings = c("NA", "")) 
F1_data$type2 <- ordered(F1_data$type2, levels = c("CC", "BCC", "WF1", "BCL", "LL", "BS"))
mosaicplot(table(F1_data$type2, F1_data$n.seeds), las=2, main=NULL, ylab = "Number of seeds per pod", xlab = "Cross type")

# test models of fruit set
F1_B_model1 <- glmer(min1seed ~ type2 + (1 | mom) + (1 | dad), family = "binomial", data = F1_data)
F1_B_model2 <- glmer(min1seed ~ (1 | mom) + (1 | dad), family = "binomial", data = F1_data)
anova(F1_B_model1, F1_B_model2)
test_model_fit(F1_B_model1)
lsmeans(F1_B_model1, pairwise ~ type2)

# test models of mean seeds
F1_N_model1 <- glmer(n.seeds ~ type2 + (1 | mom) + (1 | dad), family = "poisson", data = F1_data)
F1_N_model2 <- glmer(n.seeds ~ (1 | mom) + (1 | dad), family = "poisson", data = F1_data)
anova(F1_N_model1, F1_N_model2)
test_model_fit(F1_B_model1)
lsmeans(F1_N_model1, pairwise ~ type2)

