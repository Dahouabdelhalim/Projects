


################################################################
# create vectors
################################################################
data <- read.csv("X:\\\\Public\\\\Cheyenne\\\\survival_input.csv")
summary(data)



################################################################
# read mean clutch
################################################################
mu <- read.csv("X:\\\\Public\\\\Data_proofing_scripts\\\\Manuscripts\\\\Cheyenne\\\\mean_clutch_vols.csv", header = T)
mu <- mu[,2:4]
summary(mu)
mu$year <- mu$year - 1986
library(plyr)
data <- join(data, mu, type = 'left', by = c('year','nest'))


################################################################
# z-standardize volume by position
################################################################
data$meanVOL <- scale(data$muVOL)
summary(data$meanVOL)

data$resVOL <- data$VOL - data$muVOL
data$rVOL <- scale(data$resVOL)
summary(data$resVOL)
summary(data$rVOL)
summary(data$meanVOL)
# change all NA's to 0's for rVol and meanVOL
data$rVOL[is.na(data$rVOL)] <- 0
data$meanVOL[is.na(data$meanVOL)] <- 0
sd(data$Gmass)
mean(data$Gmass)
################################################################
# major effect
################################################################
data <- subset(data, !is.na(Gmass) & sex != 'U' & !is.na(age))

summary(data)

par(family = 'serif', mfrow = c(1,1))
plot(data$Gmass ~ data$age)

################################################################
# create brood level random effect
################################################################
data$brood <- paste(data$nest,'*',data$year,'*',data$parent1m, sep = "")
data$brood <- as.factor(data$brood); summary(data$brood)
data$brood <- as.numeric(data$brood)




################################################################
# fix 'among'
################################################################
# something
among <- data.frame(data$year, data$tagA)
among <- unique(among)
among[25,] <- c(15, 0)
among <- among[order(among$data.year),]
among <- among[!duplicated(among$data.year),]





################################################################
# peak at some things real quick
################################################################
summary(data$state)
summary(data$age)
summary(data$age)
summary(as.factor(data$pos))

summary(lm(Gmass ~ age, data = data))

################################################################
# response variable
################################################################
mass <- data$Gmass

################################################################
# covariates
################################################################
names(data)

# position in each laying sequence
pos <- data$pos; summary(pos)

# egg volume
# vol <- data$vol.all; summary(vol)
meanVOL <- data$meanVOL
rVOL <- data$rVOL

# indexed year
year <- data$year; summary(as.factor(year))

# state at tagging
S <- as.numeric(data$state); summary(S)
S[S==3] <- 1 # make dry and wet goslings both goslings
summary(S)

# brood random effect
brood <- data$brood
n.brood <- max(data$brood)

# z-standardized hatch date within each year
W <- data$tagW; summary(W)

# z-standardized mean hatch date of each year
A <- among[,2]; summary(A)

# gosling age at capture
age <- data$age

# gosling gender
sex <- rep(0, nrow(data))
for (i in 1:length(sex)){
  if(data$sex[i] == 'M'){sex[i] <- 1}
}






# install.packages('jagsUI')
library(jagsUI)
################################################################
# model
################################################################
sink("growth_model.jags")
cat("
    model {
    
    ###########################################################
    # annual temporal component
    ###########################################################
    alpha.mu ~ dnorm(0, 0.001)
    beta.mu.trend ~ dnorm(0, 0.001)
    beta.mu.phenology ~ dnorm(0, 0.001)
    tau.mass.t <- pow(sigma.mass.t, -2)
    sigma.mass.t ~ dunif(0,5)

    for (t in 1:n.years){
      mu.mass[t] <- alpha.mu + beta.mu.trend * t + beta.mu.phenology * A[t]
      eps.mass[t] ~ dnorm(mu.mass[t], tau.mass.t)
    }

    ###########################################################
    # model
    ###########################################################
    tau.mass.i <- pow(sigma.mass.i, -2)
    sigma.mass.i ~ dunif(0,200)

    beta.within ~ dnorm(0, 0.001)
    beta.sex ~ dnorm(0, 0.001)

    alpha <- 43.6

    for (j in 1:14){
      lay[j] ~ dnorm(0, 0.001)

    }

    lay[15] <- 0



    # state[1] ~ dnorm(0, 0.001)
    # state[2] <- 0


    sigma.brood ~ dunif(0,10)
    var.brood <- pow(sigma.brood, 2)
    tau.brood <- pow(sigma.brood, -2)

#################################################################################
# brood effect
#################################################################################
    for (j in 1:n.brood){
      b[j] ~ dnorm(0, tau.brood)
    }



    #################################################################################
    # individual, brood, and temporal effects on slope
    #################################################################################
    for (i in 1:n){
      mass[i] ~ dnorm(mu[i], tau.mass.i)
      slope[i] <- beta.within * W[i] + beta.sex * sex[i] + beta.meanVOL * meanVOL[i] + eps.mass[year[i]] + 
                  lay[pos[i]] + b[brood[i]] + beta.rVOL * rVOL[i] + beta.int * meanVOL[i] * rVOL[i]
      mu[i] <- alpha + slope[i] * age[i]
    }


    beta.meanVOL ~ dnorm(0, 0.0001)
    beta.rVOL ~  dnorm(0, 0.001)
    beta.int ~ dnorm(0, 0.001)



    # for (){
    #  pred.VOL <- slope[i] + beta.meanVOL * meanVOL[i] + beta.rVOL * rVOL[i] + beta.int * meanVOL * rVOL
    # }






    ###########################################################
    # model
    ###########################################################
    for (j in 1:14){
      pred.lay[j] <- 43.6 + (eps.mass[11] + lay[j]) * 30
    }

    for (k in 1:100){
      pred.vol[k] <- 43.6 + (eps.mass[11] + beta.meanVOL * egg.vol[k]) * 30
    }
    



    for (t in 1:n.years){
      pred.t[t] <- 43.6 + (alpha.mu + beta.mu.trend * t) * 30

      pred.annual[t] <- 43.6 + (eps.mass[t]) * 30
    }


    for (k in 1:100){
      pred.w[k] <- 43.6 + (eps.mass[11] + beta.within * hatch.date[k]) * 30
    }

    }
    
    ",fill = TRUE)
sink()

################################################################
# export and plot results
################################################################
#Bundle data
# dat <- list(m.sex = m.sex, m.vol = m.vol, m.W = m.W, m.Ai = m.Ai, mass = mass,
#             R = R, pos = pos, vol = vol, year = year, S = S, W = W, n.years = max(year), n = length(R), A = A, Ai = Ai, age = age, sex = sex)
egg.vol <- seq(min(meanVOL), max(meanVOL), length.out = 100)
hatch.date <- seq(min(W), max(W), length.out = 100)
dat <- list(mass = mass, A = A, n = length(sex), sex = sex, pos = pos,
            meanVOL = meanVOL, rVOL = rVOL, age = age, W = W, n.years = 21, year = year, S = S, 
            brood = brood, n.brood = n.brood, hatch.date = hatch.date, egg.vol = egg.vol)


hatch.date

# Initial values
inits <- function() list()

# Parameters to monitor
pars <- c('alpha.mu','beta.mu.trend','beta.mu.phenology',
          'beta.meanVOL','beta.rVOL','beta.int','beta.within','beta.sex','alpha.slope','sigma.mass.i','sigma.mass.t','pred.mass',
          'lay','state','eps.mass','mu.mass','pred.lay','pred.vol','pred.t','pred.w','pred.annual') #added beta.within for summary statistics


######################################################################################
# MCMC settings
######################################################################################
n.chains <- 2
n.thin <- 5
n.adapt <- 100
n.iter <- 10000
n.burnin <- 5000


######################################################################################
# compile model! 
######################################################################################
m <- jags(dat, inits, pars, "growth_model.jags", 
          n.chains = n.chains, n.adapt = n.adapt, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
summary(m)

tmp <- data.frame(summary(m))
tmp <- m
tmp <- (m$summary)
write.csv(tmp, "X:\\\\Public\\\\Cheyenne\\\\Results_m.csv")

str(m)
summary(as.factor(R))
hist(m$sims.list$beta.within)
m$mean$beta.within
m$q2.5$beta.within
m$q97.5$beta.within
m$sd$beta.within
m$f$beta.within


# m <- update(m, n.iter = 10000)

str(m)


par(mfrow = c(1,4), mar = c(5.1,5.1,2.1,2.1))
plot(m$mean$pred.lay[1:2], xlim = c(0.5,5.5), ylim = c(550,675), pch = 19, ylab = 'Mass (g)', cex.lab = 2)
arrows(seq(1,2), m$q2.5$pred.lay[1:2], seq(1,2), m$q97.5$pred.lay[1:2], lty = 2, length = 0)
plot(m$mean$pred.lay[3:5], xlim = c(0.5,5.5), ylim = c(550,675), pch = 19, ylab = 'Mass (g)', cex.lab = 2)
arrows(seq(1,3), m$q2.5$pred.lay[3:5], seq(1,3), m$q97.5$pred.lay[3:5], lty = 2, length = 0)
plot(m$mean$pred.lay[6:9], xlim = c(0.5,5.5), ylim = c(550,675), pch = 19, ylab = 'Mass (g)', cex.lab = 2)
arrows(seq(1,4), m$q2.5$pred.lay[6:9], seq(1,4), m$q97.5$pred.lay[6:9], lty = 2, length = 0)
plot(m$mean$pred.lay[10:14], xlim = c(0.5,5.5), ylim = c(550,675), pch = 19, ylab = 'Mass (g)', cex.lab = 2)
arrows(seq(1,5), m$q2.5$pred.lay[10:14], seq(1,5), m$q97.5$pred.lay[10:14], lty = 2, length = 0)












plot(m$mean$beta.vol[1:2], xlim = c(0.5,5.5), ylim = c(-0.5,3), pch = 19, ylab = expression(beta['volume']), cex.lab = 2)
arrows(seq(1,2), m$q2.5$beta.vol[1:2], seq(1,2), m$q97.5$beta.vol[1:2], lty = 2, length = 0)
plot(m$mean$beta.vol[3:5], xlim = c(0.5,5.5), ylim = c(-0.5,3), pch = 19, ylab = expression(beta['volume']), cex.lab = 2)
arrows(seq(1,3), m$q2.5$beta.vol[3:5], seq(1,3), m$q97.5$beta.vol[3:5], lty = 2, length = 0)
plot(m$mean$beta.vol[6:9], xlim = c(0.5,5.5), ylim = c(-0.5,3), pch = 19, ylab = expression(beta['volume']), cex.lab = 2)
arrows(seq(1,4), m$q2.5$beta.vol[6:9], seq(1,4), m$q97.5$beta.vol[6:9], lty = 2, length = 0)
plot(m$mean$beta.vol[10:14], xlim = c(0.5,5.5), ylim = c(-0.5,3), pch = 19, ylab = expression(beta['volume']), cex.lab = 2)
arrows(seq(1,5), m$q2.5$beta.vol[10:14], seq(1,5), m$q97.5$beta.vol[10:14], lty = 2, length = 0)






str(m)
str(m.phi)
par(family = 'sans', mar = c(5.1,5.1,2.1,2.1), mfrow = c(1,1))
plot(m$mean$pred.lay ~ m.phi$mean$pred.lay, pch = 19, ylim = c(575,675), xlim = c(0, 0.35), 
     ylab = expression(beta['growth']), xlab = expression(beta['phi']), cex.lab = 2)
arrows(m.phi$q2.5$pred.lay, m$mean$pred.lay, m.phi$q97.5$pred.lay, m$mean$pred.lay, length = 0, lty = 2)
arrows(m.phi$mean$pred.lay, m$q2.5$pred.lay, m.phi$mean$pred.lay, m$q97.5$pred.lay, length = 0, lty = 2)



plot(m$mean$eps.mass ~ m.phi$mean$eps.phi, pch = 19, ylim = c(12,28), xlim = c(-8, 0), 
     ylab = expression(eps['growth']), xlab = expression(eps['phi']), cex.lab = 2)
arrows(m.phi$q2.5$eps.phi, m$mean$eps.mass, m.phi$q97.5$eps.phi, m$mean$eps.mass, length = 0, lty = 2)
arrows(m.phi$mean$eps.phi, m$q2.5$eps.mass, m.phi$mean$eps.phi, m$q97.5$eps.mass, length = 0, lty = 2)





















par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(m$mean$pred.vol[,1] ~ seq(1,100), xaxt = 'n', pch = 19, xlim = c(1, 100), ylim = c(300,900), 
     xlab = "Egg Volume", ylab = 'Gosling Mass @ 30 days', type = 'n', cex.lab = 2)

lines(m$mean$pred.vol[,1], lwd = 2)
lines(m$mean$pred.vol[,2], lwd = 2)

lines(m$mean$pred.vol[,3], lwd = 2)
lines(m$mean$pred.vol[,4], lwd = 2)
lines(m$mean$pred.vol[,5], lwd = 2)

lines(m$mean$pred.vol[,6], lwd = 2)
lines(m$mean$pred.vol[,7], lwd = 2)
lines(m$mean$pred.vol[,8], lwd = 2)
lines(m$mean$pred.vol[,9], lwd = 2)

lines(m$mean$pred.vol[,10], lwd = 2)
lines(m$mean$pred.vol[,11], lwd = 2)
lines(m$mean$pred.vol[,12], lwd = 2)
lines(m$mean$pred.vol[,13], lwd = 2)
lines(m$mean$pred.vol[,14], lwd = 2)

lines(m$mean$pred.vol[,15], lwd = 2)



# lines(m$q2.5$pred.vol, lwd = 2, lty = 2)
# lines(m$q97.5$pred.vol, lwd = 2, lty = 2)
back.vol <- egg.vol * 5.5 + 84
simple.vol <- c(55,67.5,80,92.5,102.5)
axis(side = 1, at = seq(1,100,length.out = 5), labels = simple.vol)

print(back.vol)



par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(m$mean$pred.annual ~ seq(1987,2007), pch = 19, ylim = c(450,800),
     xlab = "Year", ylab = 'Female Gosling Mass (g) @ 30 days', cex.lab = 2)
arrows(seq(1987,2007), m$q2.5$pred.annual, seq(1987,2007), m$q97.5$pred.annual, lty = 2, length = 0)
lines(m$mean$pred.t ~ seq(1987,2007))
lines(m$q2.5$pred.t ~ seq(1987,2007), lty = 2)
lines(m$q97.5$pred.t ~ seq(1987,2007), lty = 2)




#####################################################################################
# Plots for Manuscript
####################################################################################


# gosling mass at 30 days ~ hatch date
hatch.date
par(mfrow = c(1,2), mar = c(5.1,5.1,2.1,2.1))
plot(plogis(m.phi$mean$pred.w+1.5) ~ seq(1,100), xaxt = 'n', pch = 19, ylim = c(0,0.4), 
     xlab = "Hatch Date", ylab = 'Survival', type = 'n', cex.lab = 2)
# points(mass ~ W)
lines(plogis(m.phi$mean$pred.w+1.5), lwd = 2)
lines(plogis(m.phi$q2.5$pred.w+1.5), lwd = 2, lty = 2)
lines(plogis(m.phi$q97.5$pred.w+1.5), lwd = 2, lty = 2)
hatch.back <- c(-3,-1.5,0,1.5,3)
axis(side = 1, at = seq(1,100,length.out = 5), labels = hatch.back)

plot(m$mean$pred.w ~ seq(1,100), xaxt = 'n', pch = 19, ylim = c(475,788), 
     xlab = "Hatch Date", ylab = 'Mass', type = 'n', cex.lab = 2)
# points(mass ~ W)
lines(m$mean$pred.w, lwd = 2)
lines(m$q2.5$pred.w, lwd = 2, lty = 2)
lines(m$q97.5$pred.w, lwd = 2, lty = 2)
hatch.back <- c(-3,-1.5,0,1.5,3)
axis(side = 1, at = seq(1,100,length.out = 5), labels = hatch.back)



# mean clutch volume
# mean clutch volume
egg.vol <- seq(min(meanVOL), max(meanVOL), length.out = 100)
real.vol <- seq(min(data$muVOL, na.rm = T), max(data$muVOL, na.rm = T), length.out = 5)

par(mfrow = c(1,2), mar = c(5.1,5.1,2.1,2.1))
plot(plogis(m.phi$mean$pred.vol+1.5) ~ seq(1,100), xaxt = 'n', pch = 19, ylim = c(0, 0.4), 
     xlab = expression(mu~Egg~Volume~(cm^3)), ylab = "Survival", type = 'n', cex.lab = 2)
lines(plogis(m.phi$mean$pred.vol+1.5), lwd = 2)
lines(plogis(m.phi$q2.5$pred.vol+1.5), lwd = 2, lty = 2)
lines(plogis(m.phi$q97.5$pred.vol+1.5), lwd = 2, lty = 2)
real.vol = round(real.vol)
axis(side = 1, at = seq(1,100,length.out = 5), labels = real.vol)


egg.vol <- seq(min(meanVOL), max(meanVOL), length.out = 100)
real.vol <- seq(min(data$muVOL, na.rm = T), max(data$muVOL, na.rm = T), length.out = 5)

#par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(m$mean$pred.vol ~ seq(1,100), xaxt = 'n', pch = 19, ylim = c(550,750), 
     xlab = expression(mu~Egg~Volume~(cm^3)), ylab = 'Mass', type = 'n', cex.lab = 2)
lines(m$mean$pred.vol, lwd = 2)
lines(m$q2.5$pred.vol, lwd = 2, lty = 2)
lines(m$q97.5$pred.vol, lwd = 2, lty = 2)
real.vol = round(real.vol)
axis(side = 1, at = seq(1,100,length.out = 5), labels = real.vol)




# year with trend
# year with trend

par(mfrow = c(1,2), mar = c(5.1,5.1,2.1,2.1))
plot(plogis(m.phi$mean$pred.annual+1.5) ~ seq(1987,2007), pch = 19, ylim = c(0, 0.4),
     xlab = "Year", ylab = "Survival", cex.lab = 2)
arrows(seq(1987,2007), plogis(m.phi$q2.5$pred.annual+1.5), seq(1987,2007), plogis(m.phi$q97.5$pred.annual+1.5), lty = 2, length = 0)
lines(plogis(m.phi$mean$pred.t+1.5) ~ seq(1987,2007))
lines(plogis(m.phi$q2.5$pred.t+1.5) ~ seq(1987,2007), lty = 2)
lines(plogis(m.phi$q97.5$pred.t+1.5) ~ seq(1987,2007), lty = 2)

#par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(m$mean$pred.annual ~ seq(1987,2007), pch = 19, ylim = c(450,800),
     xlab = "Year", ylab = 'Mass', cex.lab = 2)
arrows(seq(1987,2007), m$q2.5$pred.annual, seq(1987,2007), m$q97.5$pred.annual, lty = 2, length = 0)
lines(m$mean$pred.t ~ seq(1987,2007))
lines(m$q2.5$pred.t ~ seq(1987,2007), lty = 2)
lines(m$q97.5$pred.t ~ seq(1987,2007), lty = 2)


# lay position growth and survival

par(mfrow = c(2,4), mar = c(5.1,5.1,2.1,2.1))
plot(plogis(m.phi$mean$pred.lay[1:2]+1.5), xlim = c(0.5,5.5), ylim = c(0,0.3), pch = 19, ylab = "Survival", cex.lab = 2)
arrows(seq(1,2), plogis(m.phi$q2.5$pred.lay[1:2]+1.5), seq(1,2), plogis(m.phi$q97.5$pred.lay[1:2]+1.5), lty = 2, length = 0)
plot(plogis(m.phi$mean$pred.lay[3:5]+1.5), xlim = c(0.5,5.5), ylim = c(0,0.3), pch = 19, ylab = "Survival", cex.lab = 2)
arrows(seq(1,3), plogis(m.phi$q2.5$pred.lay[3:5]+1.5), seq(1,3), plogis(m.phi$q97.5$pred.lay[3:5]+1.5), lty = 2, length = 0)
plot(plogis(m.phi$mean$pred.lay[6:9]+1.5), xlim = c(0.5,5.5), ylim = c(0,0.3), pch = 19, ylab = "Survival", cex.lab = 2)
arrows(seq(1,4), plogis(m.phi$q2.5$pred.lay[6:9]+1.5), seq(1,4), plogis(m.phi$q97.5$pred.lay[6:9]+1.5), lty = 2, length = 0)
plot(plogis(m.phi$mean$pred.lay[10:14]+1.5), xlim = c(0.5,5.5), ylim = c(0,0.3), pch = 19, ylab = "Survival", cex.lab = 2)
arrows(seq(1,5), plogis(m.phi$q2.5$pred.lay[10:14]+1.5), seq(1,5), plogis(m.phi$q97.5$pred.lay[10:14]+1.5), lty = 2, length = 0)

plot(m$mean$pred.lay[1:2], xlim = c(0.5,5.5), ylim = c(560,685), pch = 19, ylab = 'Predicted Mass', cex.lab = 2)
arrows(seq(1,2), m$q2.5$pred.lay[1:2], seq(1,2), m$q97.5$pred.lay[1:2], lty = 2, length = 0)
plot(m$mean$pred.lay[3:5], xlim = c(0.5,5.5), ylim = c(560,685), pch = 19, ylab = 'Predicted Mass', cex.lab = 2)
arrows(seq(1,3), m$q2.5$pred.lay[3:5], seq(1,3), m$q97.5$pred.lay[3:5], lty = 2, length = 0)
plot(m$mean$pred.lay[6:9], xlim = c(0.5,5.5), ylim = c(560,685), pch = 19, ylab = 'Predicted Mass', cex.lab = 2)
arrows(seq(1,4), m$q2.5$pred.lay[6:9], seq(1,4), m$q97.5$pred.lay[6:9], lty = 2, length = 0)
plot(m$mean$pred.lay[10:14], xlim = c(0.5,5.5), ylim = c(560,685), pch = 19, ylab = 'Predicted Mass', cex.lab = 2)
arrows(seq(1,5), m$q2.5$pred.lay[10:14], seq(1,5), m$q97.5$pred.lay[10:14], lty = 2, length = 0)


# lay postion betwen phi and growth
windows()
shape <- c(17, 19, 17, 15, 19, 17, 15, 15, 19, 17, 15, 15, 15, 19)
color <- c('#D55E00', '#009E73', '#D55E00', '#0072B2', '#009E73', '#D55E00', '#0072B2', '#0072B2', '#009E73', '#D55E00', '#0072B2', '#0072B2', '#0072B2', '#009E73')
par(mfrow = c(1,1), mar = c(5.1, 5.1, 2.1, 2.1))
plot(plogis(m.phi$mean$pred.lay[1:14] + 1.5) ~ m$mean$pred.lay[1:14], pch= shape,
     xlab = "Mass", ylab = "Survival",
     xlim = c(575,680), ylim = c(0,0.3), cex.lab = 2, cex = 1.5, col= color, type = 'n')
# add legend
legend(575, 0.3, legend = c("First", "Other", "Last"), col =c('#D55E00','#0072B2','#009E73'), bty = "n",lty= NULL, cex = 1.2, text.col = "black", pch = c(17,15,19), pt.cex = 2, horiz = F, inset = c(0.01, 0.01))

arrows(m$mean$pred.lay[1:14],plogis(m.phi$q2.5$pred.lay[1:14]+1.5), m$mean$pred.lay[1:14], plogis(m.phi$q97.5$pred.lay[1:14]+1.5), lty = 2, length = 0, col= color)
arrows(m$q2.5$pred.lay[1:14],plogis(m.phi$mean$pred.lay[1:14] +1.5), m$q97.5$pred.lay[1:14],plogis(m.phi$mean$pred.lay[1:14]+1.5), lty = 2, length = 0, col = color)
points(plogis(m.phi$mean$pred.lay[1:14] + 1.5) ~ m$mean$pred.lay[1:14], 
       pch = shape, col = color, cex = 2)




#############################################################################
# separate
#############################################################################

windows()
par(mfrow = c(2,2), mar = c(5.1, 5.1, 2.1, 2.1))

# 2 egg clutch

plot(plogis(m.phi$mean$pred.lay[1:2] + 1.5) ~ m$mean$pred.lay[1:2], pch= shape,
     xlab = "Mass", ylab = "Survival",
     xlim = c(555,680), ylim = c(0,0.3), cex.lab = 1.5, cex = 1.5, col= color, type = 'n')
# add legend
legend(560, 0.3, legend = c("First", "Last"), col =c('#D55E00','#009E73'), bty = "n",lty= NULL, cex = 1.2, text.col = "black", pch = c(17,19), pt.cex = 2, horiz = F, inset = c(0.01, 0.01))

arrows(m$mean$pred.lay[1:2],plogis(m.phi$q2.5$pred.lay[1:2]+1.5), m$mean$pred.lay[1:2], plogis(m.phi$q97.5$pred.lay[1:2]+1.5), lty = 2, length = 0, col= color[1:2])
arrows(m$q2.5$pred.lay[1:2],plogis(m.phi$mean$pred.lay[1:2] +1.5), m$q97.5$pred.lay[1:2],plogis(m.phi$mean$pred.lay[1:2]+1.5), lty = 2, length = 0, col = color[1:2])
points(plogis(m.phi$mean$pred.lay[1:2] + 1.5) ~ m$mean$pred.lay[1:2], 
       pch = shape[1:2], col = color[1:2], cex = 2)


# 3 egg plot


plot(plogis(m.phi$mean$pred.lay[3:5] + 1.5) ~ m$mean$pred.lay[3:5], pch= shape,
     xlab = "Mass", ylab = "Survival",
     xlim = c(555,680), ylim = c(0,0.3), cex.lab = 1.5, cex = 1.5, col= color, type = 'n')
# add legend
legend(560, 0.3, legend = c("First", "Other", "Last"), col =c('#D55E00','#0072B2','#009E73'), bty = "n",lty= NULL, cex = 1.2, text.col = "black", pch = c(17,15,19), pt.cex = 2, horiz = F, inset = c(0.01, 0.01))

arrows(m$mean$pred.lay[3:5],plogis(m.phi$q2.5$pred.lay[3:5]+1.5), m$mean$pred.lay[3:5], plogis(m.phi$q97.5$pred.lay[3:5]+1.5), lty = 2, length = 0, col= color[3:5])
arrows(m$q2.5$pred.lay[3:5],plogis(m.phi$mean$pred.lay[3:5] +1.5), m$q97.5$pred.lay[3:5],plogis(m.phi$mean$pred.lay[3:5]+1.5), lty = 2, length = 0, col = color[3:5])
points(plogis(m.phi$mean$pred.lay[3:5] + 1.5) ~ m$mean$pred.lay[3:5], 
       pch = shape[3:5], col = color[3:5], cex = 2)


# 4 egg clutch

plot(plogis(m.phi$mean$pred.lay[6:9] + 1.5) ~ m$mean$pred.lay[6:9], pch= shape,
     xlab = "Mass", ylab = "Survival",
     xlim = c(555,680), ylim = c(0,0.3), cex.lab = 1.5, cex = 1.5, col= color, type = 'n')
# add legend
legend(560, 0.3, legend = c("First", "Other", "Last"), col =c('#D55E00','#0072B2','#009E73'), bty = "n",lty= NULL, cex = 1.2, text.col = "black", pch = c(17,15,19), pt.cex = 2, horiz = F, inset = c(0.01, 0.01))

arrows(m$mean$pred.lay[6:9],plogis(m.phi$q2.5$pred.lay[6:9]+1.5), m$mean$pred.lay[6:9], plogis(m.phi$q97.5$pred.lay[6:9]+1.5), lty = 2, length = 0, col= color[6:9])
arrows(m$q2.5$pred.lay[6:9],plogis(m.phi$mean$pred.lay[6:9] +1.5), m$q97.5$pred.lay[6:9],plogis(m.phi$mean$pred.lay[6:9]+1.5), lty = 2, length = 0, col = color[6:9])
points(plogis(m.phi$mean$pred.lay[6:9] + 1.5) ~ m$mean$pred.lay[6:9], 
       pch = shape[6:9], col = color[6:9], cex = 2)

# 5 egg clutch

plot(plogis(m.phi$mean$pred.lay[10:14] + 1.5) ~ m$mean$pred.lay[10:14], pch= shape,
     xlab = "Mass", ylab = "Survival",
     xlim = c(555,680), ylim = c(0,0.3), cex.lab = 1.5, cex = 1.5, col= color, type = 'n')
# add legend
legend(560, 0.3, legend = c("First", "Other", "Last"), col =c('#D55E00','#0072B2','#009E73'), bty = "n",lty= NULL, cex = 1.2, text.col = "black", pch = c(17,15,19), pt.cex = 2, horiz = F, inset = c(0.01, 0.01))

arrows(m$mean$pred.lay[10:14],plogis(m.phi$q2.5$pred.lay[10:14]+1.5), m$mean$pred.lay[10:14], plogis(m.phi$q97.5$pred.lay[10:14]+1.5), lty = 2, length = 0, col= color[10:14])
arrows(m$q2.5$pred.lay[10:14],plogis(m.phi$mean$pred.lay[10:14] +1.5), m$q97.5$pred.lay[10:14],plogis(m.phi$mean$pred.lay[10:14]+1.5), lty = 2, length = 0, col = color[10:14])
points(plogis(m.phi$mean$pred.lay[10:14] + 1.5) ~ m$mean$pred.lay[10:14], 
       pch = shape[10:14], col = color[10:14], cex = 2)




# add labels
name <- c("2,1","2,2", "3,1", "3,2", "3,3","4,1","4,2","4,3","4,4", "5,1","5,2","5,3","5,4","5,5")

library(calibrate)
textxy(m$mean$pred.lay[1:14], plogis(m.phi$mean$pred.lay[1:14] + 1.5), labs=name, cx = 0.5, dcol = "red", m = c(4, 4))


# year correlation between phi and growth
par(mfrow = c(1,1), mar = c(5.1, 5.1, 2.1, 2.1))
plot(plogis(m.phi$mean$pred.annual + 1.5) ~ m$mean$pred.annual, xlim = c(450,850), ylim = c(-0.01,0.41), pch= 19,
     xlab = expression(Mass[year]), ylab = expression(Survival[year]), cex.lab = 2)
arrows(m$mean$pred.annual,plogis(m.phi$q2.5$pred.annual +1.5), m$mean$pred.annual, plogis(m.phi$q97.5$pred.annual +1.5), lty = 2, length = 0)
arrows(m$q2.5$pred.annual,plogis(m.phi$mean$pred.annual +1.5), m$q97.5$pred.annual,plogis(m.phi$mean$pred.annual+1.5), lty = 2, length = 0)
abline(lm(plogis(m.phi$mean$pred.annual + 1.5) ~ m$mean$pred.annual), col = "blue")
#year.name <- c("1987", "1988", "1989", "1990", "1991", "1992", "1993", "1994", "1995", "1996", "1997",
               #"1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007")

#textxy(m$mean$pred.annual, m.phi$mean$pred.annual+ 1.5, labs=year.name, cx = 0.5, dcol = "red", m = c(4, 4))

# within year



# residual egg volume

egg.vol <- seq(min(meanVOL), max(meanVOL), length.out = 100)
real.vol <- seq(min(data$muVOL, na.rm = T), max(data$muVOL, na.rm = T), length.out = 5)

par(mfrow = c(1,2), mar = c(5.1,5.1,2.1,2.1))
plot(plogis(m.phi$mean$pred.vol+1.5) ~ seq(1,100), xaxt = 'n', pch = 19, ylim = c(0, 0.4), 
     xlab = expression(Mean~Clutch~Volume~(cm^3)), ylab = "Survival", type = 'n', cex.lab = 2)
lines(plogis(m.phi$mean$pred.vol+1.5), lwd = 2)
lines(plogis(m.phi$q2.5$pred.vol+1.5), lwd = 2, lty = 2)
lines(plogis(m.phi$q97.5$pred.vol+1.5), lwd = 2, lty = 2)
real.vol = round(real.vol)
axis(side = 1, at = seq(1,100,length.out = 5), labels = real.vol)


egg.vol <- seq(min(meanVOL), max(meanVOL), length.out = 100)
real.vol <- seq(min(data$muVOL, na.rm = T), max(data$muVOL, na.rm = T), length.out = 5)

#par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(m$mean$pred.vol ~ seq(1,100), xaxt = 'n', pch = 19, ylim = c(550,750), 
     xlab = expression(Mean~Clutch~Volume~(cm^3)), ylab = 'Mass', type = 'n', cex.lab = 2)
lines(m$mean$pred.vol, lwd = 2)
lines(m$q2.5$pred.vol, lwd = 2, lty = 2)
lines(m$q97.5$pred.vol, lwd = 2, lty = 2)
real.vol = round(real.vol)
axis(side = 1, at = seq(1,100,length.out = 5), labels = real.vol)





print(back.vol)


plot(mass ~ W, cex = 0.75, cex.lab = 2,
     xlab = 'Z-standardized Hatch Date', ylab = 'Gosling Mass')

