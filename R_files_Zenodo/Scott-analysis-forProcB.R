library(effects)

#########################################
####### phenology & OSR graphs ##########
#########################################

# get OSR data for mating season (July-Sept) 2016

OSRdat = read.csv("Scott-data1-OSR-procB.csv", header = T)
head(OSRdat)
tail(OSRdat)
# date = day of mating season (year-month-day)
# males = number of mature males on each date   
# females = number of females considered receptive on each date 
# (we assumed females were recpetive during a 14-day window)
# x = index for graphing (accounting for calendar dates we did not sample)

# calculate OSR as ratio of males / females
OSRdat$ratiomf = OSRdat$males / OSRdat$females
head(OSRdat)
tail(OSRdat)
OSRdat$ratiomf[OSRdat$ratiomf == Inf] <- NA # replace Inf values with NA

# plot data (figure 1)

par(mfrow=c(3, 1), oma=c(4, 1, 1, 1), mar=c(1, 4, 0, 0))
par(las = 2, mgp = c(2.75, 0.75, 0), tck = -0.025, 
    cex.axis = 1.25, cex.lab = 1.25, cex = 1, lwd = 1.5)

# males
plot(OSRdat$x, OSRdat$males, pch = 21, bg = "white", type = "b",
             xlim = c(0, 73), ylim = c(-0.02*250, 250 + 0.02*250),
             xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
             ylab = "males", bty = "l", cex = 0.75)
axis(1, at = seq(1:72), labels = rep("", 72), tck = -0.025, 
     lwd.ticks = 1.5)
axis(2, at = seq(0, 250, by = 50), labels = seq(0, 250, by = 50), 
     mgp = c(3, .5, 0), tck = -0.025, lwd.ticks = 1.5)

# females
plot(OSRdat$x, OSRdat$females, pch = 21, bg = "grey50", type = "b",
             xlim = c(0, 73), ylim = c(-1, 40 + 1),
             xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
             ylab = "females", bty = "l", cex = 0.75)
axis(1, at = seq(1:72), labels = rep("", 72), lwd.ticks = 1.5)
axis(2, at = seq(0, 40, by = 10), labels = seq(0, 40, by = 10), 
     mgp = c(3, .5, 0), lwd.ticks = 1.5)

# ratio of males to females 
plot(OSRdat$x, OSRdat$ratiomf, pch = 21, bg = "black", type = "b",
             xlim = c(0, 73), ylim = c(-0.02*25, 25 + 0.02*25),
             xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
             ylab = "males/females", bty = "l", cex = 0.75)
axis(1, at = seq(1:72), labels = rep("", 72), lwd.ticks = 1.5)
axis(1, at = OSRdat$x, labels = OSRdat$date, lwd.ticks = 1.5, 
     mgp = c(3, .5, 0), cex.axis = 0.5)
axis(2, at = seq(0, 25, by = 5), labels = seq(0, 25, by = 5), 
     mgp = c(3, .5, 0), lwd.ticks = 1.5)


########################################
### calculating body condition index ###
########################################

# read data
size = read.csv("Scott-data2-malesize-procB.csv", header = T)
str(size)
head(size)

# id = spider id
# legs = mean tibia-patella length in mm
# mass = body mass in mg
# type = lab-reared ("lab") or field-collected ("field")

# plot log(mass) against leg length
par(mfrow=c(1, 1), oma=c(1, 1, 1, 1), mar=c(4, 4, 0, 0))
par(las = 1)
plot(log(mass) ~ legs, data = size)

# calculate condition index as residuals of log(mass) against size
cond = lm(log(mass) ~ legs, data = size)
summary(cond)

abline(cond, lty = 1, lwd = 1.5)

size$cond = residuals(cond)

########################################
##### field experiment 1: 2016 #########
########################################

# this is the analysis of the experimental data shown in Figs. 2a,c and 4a,c

gr16 = read.csv("Scott-data3-race2016-procB.csv", header = T)
head(gr16)
tail(gr16)

# id = spider id 
# run.minutes = time from male's release to recapture in minutes 
# legs = male leg length in mm 
# cage = code of cage at which male was captured
# distance = distance at which male was released (in m from females) 
# recap = whether male was recaptured (1 = yes; 0 = no) 
# speed = average speed of male (distance/run.minutes) in m/min


# distribution of males among cages
table(gr16$cage) # see data file 8 for distribution of cages in space

########################################################
### recapture rates for 2016 data 
########################################################

# number released by distance
table(gr16$distance)
sum(table(gr16$distance))

# number recaptured by distance
table(gr16$recap, gr16$distance)

# calculate proportion recaptured by distance
round(table(gr16$recap, gr16$distance)[2,] / table(gr16$distance), 2)

# total recapture rate
round(sum(table(gr16$recap, gr16$distance)[2,]) / sum(table(gr16$distance)), 2)

# GLM: recapture rate ~ distance and male leg length
recap16 = glm(recap ~ distance + legs, data = gr16, family = binomial)
summary(recap16) 

# graph the effects
plot(allEffects(recap16), type = "response")

#############################################
# plot recapture rate against distance (figure 4a)

# actual data for recapture graph
dots16 = round(table(gr16$recap, gr16$distance)[2,] / table(gr16$distance), 2)

par(mar = c(4,4,1,1), oma = c(0,0,0,0))
par(las = 1)

newdist = seq(5, 65, 5)
meanlegs16 = mean(gr16$legs, na.rm = T)
newlegs16 = rep(meanlegs16, length(newdist))
newDF16 = data.frame(distance = newdist, legs = newlegs16)

# calculate predicted values and approximate 95% CI
preds <- predict(recap16, newdata = newDF16, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit
fit2 <- recap16$family$linkinv(fit) #inverse link to get back on normal scale
upr2 <- recap16$family$linkinv(upr)
lwr2 <- recap16$family$linkinv(lwr)

plot(0, xlim = c(5, 65), ylim = c(0, 1),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
     bty = "n", xlab = "", ylab = "")
box(bty = "l", lwd = 1.5)

lines(newdist, fit2, col = "black", lwd = 2)
points(seq(10, 60, 10), dots16)

lines(newdist, upr2, col = "grey55", lty = 2, lwd = 1.5)
lines(newdist, lwr2, col = "grey55", lty = 2, lwd = 1.5)

axis(1, mgp = c(3, 0.75, 0), cex.axis = 1,  tck = -0.015, lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(10, 60, 10), at = seq(10, 60, 10))
axis(2, mgp=c(3, 0.75, 0), cex.axis = 1,  tck = -0.015, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(0, 1.0, 0.2), at = seq(0, 1.0, 0.2))
title(ylab = "recapture rate", line = 2.5, cex.lab = 1)
title(xlab = "distance from females (m)", line = 2, cex.lab = 1)
#############################################

########################################################
### analysis of speed for 2016 data 
########################################################

speed16a = lm(speed ~ distance + legs, data = gr16)
# check residuals
hist(speed16a$residuals) # a bit skewed
plot(speed16a$fitted.values, speed16a$residuals) # heteroskedasticity

speed16 = lm(log(speed) ~ distance + legs, data = gr16)
hist(speed16$residuals) # now normal
plot(speed16$fitted.values, speed16$residuals) # now homoskedastic
summary(speed16)

# graph the effects
plot(allEffects(speed16), type = "response")

######################################################
# plot speed against distance (figure 1c)

par(mfrow = c(1,1))
par(las=1)
par(mar = c(4,4,1,1), oma = c(0,0,0,0))
plot(speed ~ distance, data = gr16, pch = 1, cex = 1, 
     col = "black", xlim = c(5, 65), ylim = c(0, 1.5), lwd = 1.5,
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", 
     bty = "n", xlab = "", ylab = "")
box(bty = "l", lwd = 1.5)

newpreddist16 = predict(speed16, newdata = newDF16, interval = "confidence")
lines(newdist, exp(newpreddist16[,1]), col = "black", lwd = 2)
lines(newdist, exp(newpreddist16[,2]), col = "grey55", lty = 2, lwd = 1.5)
lines(newdist, exp(newpreddist16[,3]), col = "grey55", lty = 2, lwd = 1.5)
title(xlab = "distance from females (m)", line = 1.75, cex.lab = 1)
title(ylab = "speed (m/min)", line = 2.75, cex.lab = 1)
axis(1, mgp = c(3, 0.75, 0), cex.axis = 1,  tck = -0.015, lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(10, 60, 10), at = seq(10, 60, 10))
axis(2, mgp = c(3, 0.75, 0), cex.axis = 1,  tck=-0.015, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(0, 1.5, 0.25), at = seq(0, 1.5, 0.25))
######################################################

########################################
##### field experiment 2: 2017 #########
########################################

# this is the analysis of the experimental data shown in Figs. 2b,d and 4b,d

# read data
gr17 = read.csv("Scott-data4-race2017-procB.csv", header = T)
str(gr17)

# id = spider id 
# distance = distance at which male was released (in m from females) 
# cage = location of cage where captured (0 = farthest West; 40 = farthest East)
# run.minutes = time from male's release to recapture in minutes 
# mass = body mass (weighed live) in mg
# legs = leg length in mm 
# speed = average speed of male (distance/run.minutes) in m/min
# recap = whether male was recaptured (1 = yes; 0 = no) 


# add body condition index (calculated above using a larger population of males
# including both lab-reared and field caught males)

gr17 = merge(gr17, size[c("id", "cond")], by = "id", all = F)
head(gr17)

# plot condition index against leg length to check independence
plot(cond ~ legs, data = gr17)

# distribution of males among cages
table(gr17$cage) # see also data file 8 for cage distributions in space 

########################################################
### recapture rates for 2017 data 
########################################################

# number released by distance
table(gr17$distance)
sum(table(gr17$distance))

# number recaptured by distance
table(gr17$recap, gr17$distance)

# calculate proportion recaptured by distance
round(table(gr17$recap, gr17$distance)[2,] / table(gr17$distance), 2)

# total recapture rate
round(sum(table(gr17$recap, gr17$distance)[2,]) / sum(table(gr17$distance)), 2)


# GLM: recapture rate ~ distance + male leg length + condition index
recap17 = glm(recap ~ distance + legs + cond, data = gr17, family = binomial)
summary(recap17) 

# graph the effects
plot(allEffects(recap17), type = "response")

#############################################
# plot recapture rate against distance (figure 4b)


# actual data for recapture graph
dots17 = round(table(gr17$recap, gr17$distance)[2,] / table(gr17$distance), 2)

# calculate predicted values and approximate 95% CI
newdist = seq(5, 65, 1)
meanlegs17 = mean(gr17$legs)
newlegs17 = rep(meanlegs17, length(newdist))
meancond17 = mean(gr17$cond)
newcond17 = rep(meancond17, length(newdist))
newDF17 = data.frame(distance = newdist, legs = newlegs17, cond = newcond17)

preds <- predict(recap17, newdata = newDF17, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit
fit2 <- recap17$family$linkinv(fit) #inverse link to get back on normal scale
upr2 <- recap17$family$linkinv(upr)
lwr2 <- recap17$family$linkinv(lwr)

# plot
par(mar = c(4,4,1,1), oma = c(0,0,0,0))
par(las = 1)
plot(0, xlim = c(5, 65), ylim = c(0, 1),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
     bty = "n", xlab = "", ylab = "")
box(bty = "l", lwd = 1.5)

lines(newdist, fit2, col = "black", lwd = 1.5)
points(seq(10, 60, 10), dots17)

lines(newdist, upr2, col = "grey55", lty = 2, lwd = 1.5)
lines(newdist, lwr2, col = "grey55", lty = 2, lwd = 1.5)

axis(1, mgp = c(3, 0.75, 0), cex.axis = 1, tck = -0.015, lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(10, 60, 10), at = seq(10, 60, 10))
axis(2, mgp =c(3, 0.75, 0), cex.axis = 1,  tck = -0.015, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(0, 1.0, 0.2), at = seq(0, 1.0, 0.2))
title(ylab = "recapture rate", line = 2.75, cex.lab = 1)
title(xlab = "distance from females (m)", line = 1.75, cex.lab = 1)
####################################################

########################################################
### analysis of speed for 2017 data 
########################################################

# check that leg length and body condition index are not correlated
plot(cond ~ legs, data = gr17)
summary(lm(cond ~ legs, data = gr17)) # R^2 = 0; p = 0.88  
abline(lm(cond ~ legs, data = gr17)) # no relationship

speed17a = lm(speed ~ distance + legs + cond, data = gr17)
# check residuals
hist(speed17a$residuals) # a bit skewed
plot(speed17a$fitted.values, speed17a$residuals) # not homoskedastic
summary(speed17a)

speed17 = lm(log(speed) ~ distance + legs + cond, data = gr17)
hist(speed17$residuals) # now normal
plot(speed17$fitted.values, speed17$residuals) # a bit better
summary(speed17)

# graph the effects
plot(allEffects(speed17), type = "response")

######################################################
# plot speed against distance (figure 4d)

par(mfrow = c(1,1))
par(las=1)
par(mar = c(4,4,1,1), oma = c(0,0,0,0))

# plot actual data
plot(speed ~ distance, data = gr17, pch = 1, cex = 1, 
     col = "black", xlim = c(5, 65), ylim = c(0, 0.25), lwd = 1.5,
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", 
     bty = "n", xlab = "", ylab = "")
box(bty = "l", lwd = 1.5)

# calculate and plot predicted values and 95% confidence bands 
newpreddist17 = predict(speed17, newdata = newDF17, interval = "confidence")
lines(newdist, exp(newpreddist17[,1]), col = "black", lwd = 2)
lines(newdist, exp(newpreddist17[,2]), col = "grey55", lty = 2, lwd = 1.5)
lines(newdist, exp(newpreddist17[,3]), col = "grey55", lty = 2, lwd = 1.5)
title(xlab = "distance from females (m)", line = 1.5, cex.lab = 1)
title(ylab = "speed (m/min)", line = 2.5, cex.lab = 1)
axis(1, mgp=c(3, 0.5, 0), cex.axis=1,  tck=-0.015, lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(10, 60, 10), at = seq(10, 60, 10))
axis(2, mgp=c(3, 0.75, 0), cex.axis=1,  tck=-0.015, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(0, 0.25, 0.05), at = seq(0, 0.25, 0.05))
######################################################


###########################################
##### additional field experiment (2016b)
###########################################

# this is the analysis of the experimental data shown in supplementary 
# figures S3 and S4
# males in this experiment were released 25, 50, 75, and 100 m from females
# and there were only 4 females emitting pheromones in this experiment


# read data
gr16b = read.csv("Scott-data5-race2016b-procB.csv", header = T)
head(gr16b)

# id = spider id 
# run.minutes = time from male's release to recapture in minutes 
# legs = leg length of male in mm 
# cage = code of cage at which male was captured
# distance = distance at which male was released (in m from females) 
# recap = whether male was recaptured (1 = yes; 0 = no) 
# speed = average speed of male (distance/run.minutes) in m/min


########################################################
### recapture rates for 2016b (extra) data 
########################################################

# number released by distance
table(gr16b$distance)
sum(table(gr16b$distance))

# number recaptured by distance
table(gr16b$recap, gr16b$distance)

# calculate proportion recaptured by distance
round(table(gr16b$recap, gr16b$distance)[2,] / table(gr16b$distance), 2)

# total recapture rate
round(sum(table(gr16b$recap, gr16b$distance)[2,]) / sum(table(gr16b$distance)), 2)

# GLM: recapture rate ~ distance + male leg length 
recap16b = glm(recap ~ distance + legs, data = gr16b, family = binomial)
summary(recap16b)

# graph the effects
plot(allEffects(recap16b), type = "response")

#############################################
# plot recapture rate against distance (figure S4a)

# actual data for recapture graph
dotsgr16b = round(table(gr16b$recap, gr16b$distance)[2,] / table(gr16b$distance), 2)

newdist16b = seq(5, 105, 1)
meanlegs16b = mean(gr16b$legs)
newlegs16b = rep(meanlegs16b, length(newdist16b))
newDF16b = data.frame(distance = newdist16b, legs = newlegs16b)

preds <- predict(recap16b, newdata = newDF16b, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit
fit2 <- recap16b$family$linkinv(fit) #inverse link to get back on normal scale
upr2 <- recap16b$family$linkinv(upr)
lwr2 <- recap16b$family$linkinv(lwr)

plot(0, cex = 0, xlim = c(5, 105), ylim = c(0, 1),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", 
     bty = "n", xlab = "", ylab = "")
box(bty = "l", lwd = 1.5)

lines(newdist16b, fit2, col = "black", lwd = 1.5)
points(seq(25, 100, 25), dotsgr16b)

lines(newdist16b, upr2, col = "grey55", lty = 2, lwd = 1.5)
lines(newdist16b, lwr2, col = "grey55", lty = 2, lwd = 1.5)

axis(1, mgp = c(3, 0.75, 0), cex.axis = 1,  tck = -0.015, lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(10, 100, 10), at = seq(10, 100, 10))
axis(2, mgp=c(3, 0.75, 0), cex.axis = 1,  tck = -0.015, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(0, 1.0, 0.2), at = seq(0, 1.0, 0.2))
title(ylab = "recapture rate", line = 2.5, cex.lab = 1)
title(xlab = "distance from females (m)", line = 1.75, cex.lab = 1)

########################################################
### analysis of speed for 2016b (extra) data 
########################################################

speed16b = lm(speed ~ distance + legs, data = gr16b)
# check residuals
hist(speed16b$residuals) # ok
plot(speed16b$fitted.values, speed16b$residuals) # ok

summary(speed16b)

plot(allEffects(speed16b), type = "response")

####################################################
# plot recapture speed against distance (figure S4b)

# plot actual data
plot(speed ~ distance, data = gr16b, pch = 1, cex = 1, 
     col = "black", ylim = c(0, 0.75), xlim = c(5, 65), lwd = 1.5,
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", 
     bty = "n", xlab = "", ylab = "")
box(lwd = 1.5, bty = "l") 

newdist16b = seq(5, 65, 5)
meanlegs16b = mean(gr16b$legs)
newlegs16b = rep(meanlegs16b, length(newdist16b))
newDF16b = data.frame(distance = newdist16b, legs = newlegs16b)

newpreddist2 = predict(speed16b, newdata = newDF16b, interval = "confidence")
lines(newdist16b, newpreddist2[,1], col = "black", lwd = 1.5)
lines(newdist16b, newpreddist2[,2], col = "grey55", lty = 2, lwd = 1.5)
lines(newdist16b, newpreddist2[,3], col = "grey55", lty = 2, lwd = 1.5)
title(xlab = "distance from females (m)", line = 1.5, cex.lab = 1)
title(ylab = "speed (m/min)", line = 2.5, cex.lab = 1)
axis(1, mgp=c(3, 0.75, 0), cex.axis = 1,  tck = -0.015, lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(10, 60, 10), at = seq(10, 60, 10))
axis(2, mgp=c(3, 0.75, 0), cex.axis = 1,  tck = -0.015, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(0, 0.75, 0.25), at = seq(0, 0.75, 0.25))


##########################################
### analyses of laboratory experiments ###
##########################################

#######################################################
### x-race 1 (two widow males) with pheromone and wind 
#######################################################

# this is the analysis for the data shown in Fig. 3b,g,h

# read data
xrace1 = read.csv("Scott-data6-xrace-procB.csv", header = T)
head(xrace1)

# legs = male tibia-patella length (mm)
# mass = male mass (mg)
# id = male id
# rep = replicate number
# order = whether male was first (1) or second (2) to enter the maze
# timeaftX = time to travel from intersection of X to end of maze (min)
# speedaftX = speed to travel 1.4 m from intersection of X to end of maze (m/min)
# choice = arm of X chosen after intersection (left or right)
# agreement = whether male 2 chose same side as male 1 (yes) or not (no)


# treat order as a factor
xrace1$order = as.factor(xrace1$order)
str(xrace1)

# look at number of second males who followed silk of first male
table(xrace1$agreement)

# do binomial test (null: probability of choosing same side is 50%)
binom.test(x = 21, n = 22, alternative = "two.sided", conf.level = 0.95)

# generate id for males used in this experiment
xrace1$id = paste(xrace1$rep, xrace1$order, sep = ".")
xrace1$id

# merge with condition data calculated above using male id
xrace1 = merge(xrace1, size) # note size data for 2 males is missing

# plot condition index against leg length to check independence
par(mar = c(4,4,1,1))
plot(cond ~ legs, data = xrace1)
rl = lm(cond ~ legs, data = xrace1)
summary(rl) # R^2 = 0.34; P < 0.0001 
abline(rl$coefficients)

# analyze size and condition index separately because they are 
# not independent for this data set

# look at size and order (= silk availability) first
speedxl = lm(speedaftX ~ legs + order, data = xrace1)
plot(speedxl$residuals)
plot(speedxl$residuals, speedxl$fitted.values)
summary(speedxl) 
plot(allEffects(speedxl), type = "response")

# now look at condition index and order (= silk availability)
speedxc = lm(speedaftX ~ cond + order, data = xrace1)
plot(speedxc$residuals)
plot(speedxc$residuals, speedxc$fitted.values)
summary(speedxc) # r2 = .23
plot(allEffects(speedxc), type = "response")

# graph these results all together (data shown in fig 3g,h)
# divide data into 2 sets 
m1data = subset(xrace1, order == 1)
m2data = subset(xrace1, order == 2)

par(mfrow = c(1, 2), las = 1)
par(mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))
plot(speedaftX ~ legs, data = xrace1, lwd = 1.5, cex = 1.25,
     pch = 21, bg = c("white", "grey50")[as.numeric(order)],
     xlim = c(4, 7), ylim = c(0, 3),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
     bty = "n", xlab = "", ylab = "")
size1 = lm(speedaftX ~ legs, data = m1data)
abline(size1$coefficients, lwd = 1.5, col = "black", lty = 2)
size2 = lm(speedaftX ~ legs, data = m2data)
abline(size2$coefficients, lwd = 1.5, col = "black")

box(lwd = 1.5, bty = "l") 
axis(1, mgp = c(3, 0.8, 0), cex.axis = 1.25,  tck = -0.04, lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(4, 7, 1), at = seq(4, 7, 1))
axis(1, mgp = c(3, 0.8, 0), cex.axis = 1.25,  tck = -0.04, lwd = 1.5, lwd.ticks = 1.5,
     labels = rep("", 7), at = seq(4, 7, 0.5))
axis(2, mgp = c(3, 0.8, 0), cex.axis = 1.25,  tck = -0.04, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(0.5, 3, 1), at = seq(0.5, 3, 1))
axis(2, mgp = c(3, 0.8, 0), cex.axis = 1.25,  tck = -0.04, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = c("0", "1.0", "2.0", "3.0"), at = seq(0, 3, 1))
title(xlab = "tibia-patella length (mm)", line = 2.5, cex.lab = 1.25)
title(ylab = "speed (m/min)", line = 2.5, cex.lab = 1.25)
text(4.25, 2.75, labels = "g", cex = 1.5)

plot(speedaftX ~ cond, data = xrace1, lwd = 1.5, cex = 1.25,
     pch = 21, bg = c("white", "grey50")[as.numeric(order)],
     xlim = c(-1, 0), ylim = c(0, 3),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
     bty = "n", xlab = "", ylab = "")
cond1 = lm(speedaftX ~ cond, data = m1data)
abline(cond1$coefficients, lwd = 1.5, col = "black", lty = 2)
cond2 = lm(speedaftX ~ cond, data = m2data)
abline(cond2$coefficients, lwd = 1.5, col = "black")

box(lwd = 1.5, bty = "l") 
axis(1, mgp = c(3, 0.8, 0), cex.axis = 1.25,  tck = -0.04, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = seq(-.75, 1, 0.25), at = seq(-.75, 1, 0.25))
axis(1, mgp = c(3, 0.8, 0), cex.axis = 1.25,  tck = -0.04, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = "-1.0", at = -1)
axis(2, mgp = c(3, 0.8, 0), cex.axis = 1.25,  tck = -0.04, 
     lwd = 1.5, lwd.ticks = 1.5,
     labels = rep("", 7), at = seq(0, 3, 0.5))
title(xlab = "body condition index", line = 2.5, cex.lab = 1.25)
text(-0.9, 2.75, labels = "h", cex = 1.5)

#######################################################
### x-race 2 (two widow males) with wind only 
#######################################################

# data: n = 22; followed silk = 13

# do binomial test (null: probability of choosing same side is 50%)
binom.test(x = 13, n = 22, alternative = "two.sided", conf.level = 0.95)
# 0.59 chose silk; p = 0.52

#######################################################
### x-race 3 (two widow males) no wind and no pheromone 
#######################################################

# data: n = 22; followed silk = 18

# do binomial test (null: probability of choosing same side is 50%)
binom.test(x = 18, n = 22, alternative = "two.sided", conf.level = 0.95)
# 0.818 chose silk; p = 0.0043

#######################################################
### x-race 4 & 5 widow males and Steatoda grossa silk 
#######################################################

# race 4: widow males chose between Steatoda silk and control string
# data: n = 32; followed silk = 16

# do binomial test (null: probability of choosing same side is 50%)
binom.test(x = 16, n = 32, alternative = "two.sided", conf.level = 0.95)
# 0.5 chose silk; p = 1

# race 4: widow males chose between Steatoda silk and widow silk
# note we ran 16 replicates but in one the male did not sample both options
# data: n = 15; followed silk = 13 

# do binomial test (null: probability of choosing same side is 50%)
binom.test(x = 13, n = 15, alternative = "two.sided", conf.level = 0.95)
# 0.867 chose silk; p = 0.0074


########################################
####### male movements #################
########################################

# read data
mm = read.csv("Scott-data7-male-movements-procB.csv", header = T)
head(mm)

# date = day male was marked (year-month-day)
# location = location (log or rock ID) where male was first found
# type = male first maked as penultimat (= pmtomm) or mature (= mm)
# number = id code number
# colour = paint colour
# totalfemales = total number of female webs on which male was observed
# start = start location (log or rock ID)
# move1 = location (log or rock ID) after first movement
# move2 = location (log or rock ID) after second movement (if any)
# dist1 = distance traveled between start and move1 locations (in metres) 
# dist2 = distance traveled between move1 and move2 locations (in metres)
# totaldist = dist 1 + dist 2


# graph total distances moved (Supplemental figure S2)

par(mar= c(4, 4, 2, 2))
par(mfrow = c(1,1))
par(las = 1, lwd = 1.75)
hist(mm$totaldist, breaks = 10, xaxt = "n", yaxt = "n", lwd = 1.5,
     xaxs = "i", yaxs = "i", ylim = c(0, 18), 
     col = c(rep("grey", 3), rep("white", 11)),
     xlab = "", ylab = "", main = "")
axis(1, at = seq(0, 280, by = 20), labels = seq(0, 280, by = 20),
     mgp = c(3, 0.7, 0), cex.axis = 1.25, tck = -0.02, lwd = 1.75, lwd.ticks = 1.75)
axis(2, at = seq(0, 18, by = 2), labels = seq(0, 18, by = 2),
     mgp = c(3, 0.7, 0), cex.axis = 1.25, tck = -0.02, lwd = 1.75, lwd.ticks = 1.75)
title(xlab = "Distance (m)", line = 2.5, cex.lab = 1.25)
title(ylab = "Frequency", line = 2.5, cex.lab = 1.25)

# graph total distances moved for movements < 60 m (inset of fig S2)

par(mar= c(4, 4, 1, 1))
par(las = 1, lwd = 1.5)
hist(mm$totaldist, breaks = 40, lwd = 1.5,
     xlim = c(0, 60), ylim = c(0, 6), 
     xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", col = "grey",
     xlab = "", ylab = "", main = "")
axis(1, at = seq(0, 60, by = 5), labels = rep("", 13),
     mgp = c(3, 0.7, 0), cex.axis = 1.25, tck = -0.015, lwd = 1.5, lwd.ticks = 1.5)
axis(1, at = seq(0, 60, by = 10), labels = seq(0, 60, by = 10),
     mgp = c(3, 0.7, 0), cex.axis = 1.25, tck = -0.015, lwd = 1.5, lwd.ticks = 1.5)
axis(2, at = seq(0, 6, by = 1), labels = seq(0, 6, by = 1),
     mgp = c(3, 0.7, 0), cex.axis = 1.15, tck = -0.015, lwd = 1.5, lwd.ticks = 1.5)
title(xlab = "Distance (m)", line = 2, cex.lab = 1.25)
title(ylab = "Frequency", line = 2, cex.lab = 1.25)



