###############################################################
##### notebook for analyzing presence/absence of infections ###
#####        in virus diversity experiment (RPV & PAV)      ###
###############################################################

setwd("C:/Users/straussa/Documents/Research/Seabloom-Borer Lab/Virus Diversity Rstar/3x virus/manuscript/Ecology and Evolution/data package")

library(gplots) # for error bars
library(plyr) # for data manipulation
library(plot3D) # for 3D graphics


###############################################################
######### 1) read in data and trim missing values #############
###############################################################

data <- read.csv("vir.div.np.csv", na.strings=c("","NA"))
missing <- data[is.na(data$inf),] # which are missing?
data <- data[!is.na(data$inf),] # omit the NA's
# stupidly, need N and P to be capital
colnames(data)[1] <- "N"
colnames(data)[2] <- "P"
# convert levels of N and P into concentrations (uM)
data$N <- ifelse(data$N==1,7.5,ifelse(data$N==2,52.5,375))
data$P <- ifelse(data$P==1,1,ifelse(data$P==2,7,50))
unique(log(data$N))
unique(log(data$P))
# define N:P ratios
data$NP <- data$N/data$P
data$PN <- data$P/data$N
unique(data$NP)
unique(log(data$NP))
hist(data$NP)
hist(log(data$NP))

# infections as numeric
data$inf <- ifelse(data$inf=="y",1,0)
# for interpreting single first, then mixed:
data$div <- relevel(data$div, "single")

# for interpreting intercept as lowest level N & P:
data$logN <- log(data$N)-log(7.5)


###############################################################
#  2) summarize prevalence & s.e. for single and coinfections #
###############################################################

data.RPV <- data[data$vir=="RPV",]
data.PAV <- data[data$vir=="PAV",]
data.SGV <- data[data$vir=="SGV",]

# calculate prevalence for each group 
data$trID <- paste(data$N, data$P, data$vir, data$div) 
unique(data$trID)
treats <- unique(data$trID)
prev.sum <- ddply(data, .(trID), summarize, N=N[1], P=P[1], NP=NP[1], vir=vir[1], div=div[1],
                  prev=sum(inf)/length(inf),  
                  prev.se=sqrt((sum(inf)/length(inf)*(1-sum(inf)/length(inf)))/length(inf)))
prev.sum <- prev.sum[order(prev.sum$N,prev.sum$P),] # order

RPV.sin <- prev.sum[prev.sum$div=="single" & prev.sum$vir=="RPV",]
RPV.mix <- prev.sum[prev.sum$div=="mix" & prev.sum$vir=="RPV",]
PAV.sin <- prev.sum[prev.sum$div=="single" & prev.sum$vir=="PAV",]
PAV.mix <- prev.sum[prev.sum$div=="mix" & prev.sum$vir=="PAV",]
SGV.sin <- prev.sum[prev.sum$div=="single" & prev.sum$vir=="SGV",]
SGV.mix <- prev.sum[prev.sum$div=="mix" & prev.sum$vir=="SGV",]

mean(SGV.sin$prev)
mean(RPV.sin$prev)

# summary of coinfections:
data.mix <- data[data$div=="mix",]
data.mix$pID <- paste(data.mix$N, data.mix$P, data.mix$rep) # plant ID (for coinfections)
mix.plants <- unique(data.mix$pID)
coinfections <- data.frame(plant=mix.plants)
for(i in 1:length(mix.plants)){
  t.data <- data.mix[data.mix$pID==mix.plants[i],]
  coinfections$N[i] <- t.data$N[1]; coinfections$P[i] <- t.data$P[1]
  coinfections$NP[i] <- t.data$NP[1]
  coinfections$rep[i] <- t.data$rep[1]; coinfections$tot.mass[i] <- t.data$tot.mass[1]
  coinfections$co.any[i] <- ifelse(sum(t.data$inf)>1,1,0) # any kind of coinfection <1
  coinfections$co.RPV.PAV[i] <- ifelse(sum(t.data[t.data$vir!="SGV",]$inf)>1,1,0)
  coinfections$co.RPV.SGV[i] <- ifelse(sum(t.data[t.data$vir!="PAV",]$inf)>1,1,0)
  coinfections$co.PAV.SGV[i] <- ifelse(sum(t.data[t.data$vir!="RPV",]$inf)>1,1,0)
  coinfections$co.all[i] <- ifelse(sum(t.data$inf)>2,1,0) # all three viruses =3
  coinfections$richness[i] <- sum(t.data$inf)
  coinfections$atleast1[i] <- ifelse(sum(t.data$inf)>0,1,0)
}
coinfections$logN <- log(coinfections$N)-log(7.5)

sum(coinfections$co.any) # 60! much better than just 13 RPV + PAV coinfections
sum(coinfections$co.all) # 6
sum(coinfections$co.RPV.PAV) # 13
sum(coinfections$co.RPV.SGV) # 48
sum(coinfections$co.PAV.SGV) # 11

mean(coinfections$richness)
mean(coinfections[coinfections$richness>0,]$richness)

coinfections$trID <- paste(coinfections$N, coinfections$P)
coin.sum <- ddply(coinfections, .(trID), summarize, N=N[1], P=P[1], NP=NP[1],
                  prev.any=sum(co.any)/length(co.any), prev.all=sum(co.all)/length(co.all), 
                  prev.RPV.PAV=sum(co.RPV.PAV)/length(co.RPV.PAV), prev.RPV.SGV=sum(co.RPV.SGV)/length(co.RPV.SGV), 
                  prev.PAV.SGV=sum(co.PAV.SGV)/length(co.PAV.SGV),
                  prev.atleast1=sum(atleast1)/length(atleast1),
                  v.rich=mean(richness), v.rich.no0=mean(richness/atleast1, na.rm=T),
                  prev.any.se=sqrt((sum(co.any)/length(co.any)*(1-sum(co.any)/length(co.any)))/length(co.any)),
                  prev.all.se=sqrt((sum(co.all)/length(co.all)*(1-sum(co.all)/length(co.all)))/length(co.all)),
                  prev.RPV.PAV.se=sqrt((sum(co.RPV.PAV)/length(co.RPV.PAV)*(1-sum(co.RPV.PAV)/length(co.RPV.PAV)))/length(co.RPV.PAV)),
                  prev.RPV.SGV.se=sqrt((sum(co.RPV.SGV)/length(co.RPV.SGV)*(1-sum(co.RPV.SGV)/length(co.RPV.SGV)))/length(co.RPV.SGV)),
                  prev.PAV.SGV.se=sqrt((sum(co.PAV.SGV)/length(co.PAV.SGV)*(1-sum(co.PAV.SGV)/length(co.PAV.SGV)))/length(co.PAV.SGV)),
                  prev.atleast1.se=sqrt((sum(atleast1)/length(atleast1)*(1-sum(atleast1)/length(atleast1)))/length(atleast1)),
                  v.rich.se=sd(richness)/sqrt(length(richness)),
                  v.rich.no0.se=sd(richness/atleast1, na.rm=T)/sqrt(sum(atleast1)))


###############################################################
###################### 3) linear models  ######################
###############################################################

# RPV:
summary(glm(inf ~ (logN + log(P) + div)^2, data=data.RPV, family=binomial(link="logit"))) # plus some significant interactions. AIC 329
summary(glm(inf ~ log(NP) * div, data=data.RPV, family=binomial(link="logit")))
summary(glm(inf ~ log(NP) + div, data=data.RPV, family=binomial(link="logit"))) # both significant if removing interaction
summary(glm(inf ~ logN * log(P), data=data.RPV[data.RPV$div=="single",], family=binomial(link="logit"))) 
summary(glm(inf ~ logN + log(P), data=data.RPV[data.RPV$div=="single",], family=binomial(link="logit"))) 
summary(glm(inf ~ log(NP), data=data.RPV[data.RPV$div=="single",], family=binomial(link="logit"))) 
summary(glm(inf ~ logN * log(P), data=data.RPV[data.RPV$div=="mix",], family=binomial(link="logit"))) 
summary(glm(inf ~ logN + log(P), data=data.RPV[data.RPV$div=="mix",], family=binomial(link="logit"))) 
summary(glm(inf ~ log(NP), data=data.RPV[data.RPV$div=="mix",], family=binomial(link="logit"))) 
# artificially inflating sample size:
summary(glm(inf ~ logN * log(P), data=rbind(data.RPV[data.RPV$div=="single",],
                                            data.RPV[data.RPV$div=="single",]), family=binomial(link="logit"))) 
summary(glm(inf ~ logN + log(P), data=rbind(data.RPV[data.RPV$div=="single",],
                                            data.RPV[data.RPV$div=="single",]), family=binomial(link="logit"))) 
summary(glm(inf ~ log(NP), data=rbind(data.RPV[data.RPV$div=="single",],
                                      data.RPV[data.RPV$div=="single",]), family=binomial(link="logit"))) 

# SGV:
summary(glm(inf ~ (logN + log(P) + div)^2, data=data.SGV, family=binomial(link="logit"))) # plus some significant interactions. AIC 329
summary(glm(inf ~ log(NP) * div, data=data.SGV, family=binomial(link="logit")))
summary(glm(inf ~ logN * log(P), data=data.SGV[data.SGV$div=="single",], family=binomial(link="logit"))) 
summary(glm(inf ~ logN + log(P), data=data.SGV[data.SGV$div=="single",], family=binomial(link="logit"))) 
summary(glm(inf ~ log(NP), data=data.SGV[data.SGV$div=="single",], family=binomial(link="logit"))) 
summary(glm(inf ~ logN * log(P), data=data.SGV[data.SGV$div=="mix",], family=binomial(link="logit"))) 
summary(glm(inf ~ log(NP), data=data.SGV[data.SGV$div=="mix",], family=binomial(link="logit"))) 
# artificially inflating sample size:
summary(glm(inf ~ logN * log(P), data=rbind(data.SGV[data.SGV$div=="single",],
                                            data.SGV[data.SGV$div=="single",]), family=binomial(link="logit"))) 
summary(glm(inf ~ logN + log(P), data=rbind(data.SGV[data.SGV$div=="single",],
                                            data.SGV[data.SGV$div=="single",]), family=binomial(link="logit"))) 
summary(glm(inf ~ log(NP), data=rbind(data.SGV[data.SGV$div=="single",],
                                      data.SGV[data.SGV$div=="single",]), family=binomial(link="logit"))) 

# PAV:
summary(glm(inf ~ (logN + log(P) + div)^2, data=data.PAV, family=binomial(link="logit"))) # plus some significant interactions. AIC 329
summary(glm(inf ~ (logN + log(P) + div), data=data.PAV, family=binomial(link="logit"))) #N significant if remove interaction
summary(glm(inf ~ log(NP) * div, data=data.PAV, family=binomial(link="logit")))
summary(glm(inf ~ log(NP) + div, data=data.PAV, family=binomial(link="logit"))) # NS
summary(glm(inf ~ logN * log(P), data=data.PAV[data.PAV$div=="single",], family=binomial(link="logit"))) 
summary(glm(inf ~ logN + log(P), data=data.PAV[data.PAV$div=="single",], family=binomial(link="logit"))) 
summary(glm(inf ~ log(NP), data=data.PAV[data.PAV$div=="single",], family=binomial(link="logit"))) 
summary(glm(inf ~ logN * log(P), data=data.PAV[data.PAV$div=="mix",], family=binomial(link="logit"))) 
summary(glm(inf ~ log(NP), data=data.PAV[data.PAV$div=="mix",], family=binomial(link="logit"))) 
# artificially inflating sample size:
summary(glm(inf ~ logN * log(P), data=rbind(data.PAV[data.PAV$div=="single",],
                                            data.PAV[data.PAV$div=="single",]), family=binomial(link="logit"))) 
summary(glm(inf ~ logN + log(P), data=rbind(data.PAV[data.PAV$div=="single",],
                                            data.PAV[data.PAV$div=="single",]), family=binomial(link="logit"))) 
summary(glm(inf ~ log(NP), data=rbind(data.PAV[data.PAV$div=="single",],
                                      data.PAV[data.PAV$div=="single",]), family=binomial(link="logit"))) 

# more post-hoc for SGV:
summary(glm(inf ~ logN * div, data=data.SGV[data.SGV$P==50,], family=binomial(link="logit")))
summary(glm(inf ~ log(P) * div, data=data.SGV[data.SGV$N==7.5,], family=binomial(link="logit")))

# N x P effects on all single infections:
data$vir <- relevel(data$vir, "RPV")
# summary(glm(inf ~ (logN + log(P) + vir)^2, data=data[data$div=="single",], 
#             family=binomial(link="logit"))) # nothing significant. 
summary(glm(inf ~ (logN + log(P) + vir), data=data[data$div=="single",], 
            family=binomial(link="logit"))) # only sig is higher prev for RPV
summary(glm(inf ~ (logN + log(P)), data=data[data$div=="single",], 
            family=binomial(link="logit"))) # 
summary(glm(inf ~ log(NP) + vir, data=data[data$div=="single",], 
            family=binomial(link="logit"))) # 
summary(glm(inf ~ log(NP), data=data[data$div=="single",], 
            family=binomial(link="logit"))) # 


# N x P effects on viral richness:
summary(glm(atleast1 ~ logN * log(P), data=coinfections, family=binomial(link="logit"))) # interaction
summary(glm(richness ~ logN * log(P), data=coinfections, family = "poisson")) # interaction! close to .05
summary(glm(richness ~ logN * log(P), data=coinfections[coinfections$atleast1==1,], family = "poisson")) 
summary(glm(richness ~ logN + log(P), data=coinfections[coinfections$atleast1==1,], family = "poisson")) 

summary(glm(atleast1 ~ log(NP), data=coinfections, family=binomial(link="logit"))) # interaction
summary(glm(richness ~ log(NP), data=coinfections, family = "poisson")) # interaction! close to .05
summary(glm(richness ~ log(NP), data=coinfections[coinfections$atleast1==1,], family = "poisson")) 

# N x P effects on prevalence of specific coinfections:
summary(glm(co.any ~ logN * log(P), data=coinfections, family=binomial(link="logit"))) # interaction
summary(glm(co.any ~ log(NP), data=coinfections, family=binomial(link="logit"))) # trend
summary(glm(co.all ~ logN * log(P), data=coinfections, family=binomial(link="logit"))) # nothing
summary(glm(co.all ~ logN + log(P), data=coinfections, family=binomial(link="logit"))) # nothing
summary(glm(co.all ~ log(NP), data=coinfections, family=binomial(link="logit"))) # nothing
summary(glm(co.RPV.PAV ~ logN * log(P), data=coinfections, family=binomial(link="logit")))
summary(glm(co.RPV.PAV ~ logN + log(P), data=coinfections, family=binomial(link="logit")))
summary(glm(co.RPV.PAV ~ log(NP), data=coinfections, family=binomial(link="logit")))
summary(glm(co.RPV.SGV ~ logN * log(P), data=coinfections, family=binomial(link="logit"))) # interaction
summary(glm(co.RPV.SGV ~ log(NP), data=coinfections, family=binomial(link="logit"))) # yes
summary(glm(co.PAV.SGV ~ logN * log(P), data=coinfections, family=binomial(link="logit")))
summary(glm(co.PAV.SGV ~ logN + log(P), data=coinfections, family=binomial(link="logit")))
summary(glm(co.PAV.SGV ~ log(NP), data=coinfections, family=binomial(link="logit")))


###############################################################
###### 4) prediction planes for heat maps and 3D plots  #######
###############################################################


##########################
# For the heat maps AND 3D planes:
# (including pseudo-R2 tests)

# figure 1A
lm.RPV.sin <- glm(inf ~ log(N) * log(P), data=data.RPV[data.RPV$div=="single",], 
                  family=binomial(link="logit"))
summary(lm.RPV.sin) 
null.RPV.sin <- glm(inf ~1, data=data.RPV[data.RPV$div=="single",], 
                    family=binomial(link="logit"))
1-logLik(lm.RPV.sin)/logLik(null.RPV.sin) # R2 = 0.070

# figure 1B
lm.SGV.sin <- glm(inf ~ log(N) * log(P), data=data.SGV[data.SGV$div=="single",], 
                  family=binomial(link="logit")) 
summary(lm.SGV.sin) 
null.SGV.sin <- glm(inf ~1, data=data.SGV[data.SGV$div=="single",], 
                    family=binomial(link="logit"))
1-logLik(lm.SGV.sin)/logLik(null.SGV.sin) # 0.031

# figure 1C
lm.PAV.sin <- glm(inf ~ log(N) * log(P), data=data.PAV[data.PAV$div=="single",], 
                  family=binomial(link="logit")) 
summary(lm.PAV.sin) 
null.PAV.sin <- glm(inf ~1, data=data.PAV[data.PAV$div=="single",], 
                    family=binomial(link="logit"))
1-logLik(lm.PAV.sin)/logLik(null.PAV.sin) #0.021

# figure 1D
lm.any.sin <- glm(inf ~ log(N) * log(P), data=data[data$div=="single",], 
                  family=binomial(link="logit")) 
summary(lm.any.sin) 
null.any.sin <- glm(inf ~1, data=data[data$div=="single",], 
                    family=binomial(link="logit"))
1-logLik(lm.any.sin)/logLik(null.any.sin) # 0.007

# extension to test for diferences among viruses in single infections:
summary(glm(inf ~ log(N) * log(P) + vir, data=data[data$div=="single",], 
            family=binomial(link="logit"))) # RPV is higher prev
summary(glm(inf ~ log(N) * log(P) + vir, data=data[data$div=="single" & data$vir!="RPV",], 
            family=binomial(link="logit"))) # but no difference PAV vs. SGV

# figure 1E
lm.atleast1 <- glm(atleast1 ~ log(N) * log(P), data=coinfections, 
                   family=binomial(link="logit")) 
summary(lm.atleast1) 
null.atleast1 <- glm(atleast1 ~1, data=coinfections, 
                     family=binomial(link="logit"))
1-logLik(lm.atleast1)/logLik(null.atleast1) # 0.13


##########################
# More for 3D graphics:
# each virus in co-inoculation
lm.RPV.mix <- glm(inf ~ log(N) * log(P), data=data.RPV[data.RPV$div=="mix",], 
                  family=binomial(link="logit")) 
lm.PAV.mix <- glm(inf ~ log(N) * log(P), data=data.PAV[data.PAV$div=="mix",], 
                  family=binomial(link="logit")) 
lm.SGV.mix <- glm(inf ~ log(N) * log(P), data=data.SGV[data.SGV$div=="mix",], 
                  family=binomial(link="logit")) 

# planes for coinfection models:
lm.co.RPV.PAV <- glm(co.RPV.PAV ~ log(N) * log(P), data=coinfections, 
                     family=binomial(link="logit")) 
lm.co.RPV.SGV <- glm(co.RPV.SGV ~ log(N) * log(P), data=coinfections, 
                     family=binomial(link="logit")) 
lm.co.PAV.SGV <- glm(co.PAV.SGV ~ log(N) * log(P), data=coinfections, 
                     family=binomial(link="logit")) 
lm.co.all <- glm(co.all ~ log(N) * log(P), data=coinfections, 
                 family=binomial(link="logit")) 
lm.co.any <- glm(co.any ~ log(N) * log(P), data=coinfections, 
                 family=binomial(link="logit")) 
lm.rich <- glm(richness ~ log(N) * log(P), data=coinfections, 
               family="poisson") 
lm.rich.no0 <- glm(richness ~ log(N) * log(P), 
                   data=coinfections[coinfections$atleast1==1,], 
                   family="poisson") 

##########################
# for all planes:
grid.lines=25 # too many? hard to see some data points, but planes are clear
#grid.lines=10 # too few. harder to see planes, and some data points still hard to see anyway
x.pred <- exp(seq(log(min(RPV.sin$N)), log(max(RPV.sin$N)), length.out = grid.lines))
y.pred <- exp(seq(log(min(RPV.sin$P)), log(max(RPV.sin$P)), length.out = grid.lines))
xy <- expand.grid(N = x.pred, P = y.pred)

##########################
# predicted values for each virus separately:
z.pred.RPV.sin <- matrix(predict.glm(lm.RPV.sin, newdata = xy, type = "response"), 
                         nrow = grid.lines, ncol = grid.lines)
z.pred.RPV.mix <- matrix(predict.glm(lm.RPV.mix, newdata = xy, type = "response"), 
                         nrow = grid.lines, ncol = grid.lines)
z.pred.PAV.sin <- matrix(predict.glm(lm.PAV.sin, newdata = xy, type = "response"), 
                         nrow = grid.lines, ncol = grid.lines)
z.pred.PAV.mix <- matrix(predict.glm(lm.PAV.mix, newdata = xy, type = "response"), 
                         nrow = grid.lines, ncol = grid.lines)
z.pred.SGV.sin <- matrix(predict.glm(lm.SGV.sin, newdata = xy, type = "response"), 
                         nrow = grid.lines, ncol = grid.lines)
z.pred.SGV.mix <- matrix(predict.glm(lm.SGV.mix, newdata = xy, type = "response"), 
                         nrow = grid.lines, ncol = grid.lines)
z.pred.any.sin <- matrix(predict.glm(lm.any.sin, newdata = xy, type = "response"), 
                         nrow = grid.lines, ncol = grid.lines)

##########################
# predicted values for coinfections:
z.pred.co.RPV.PAV <- matrix(predict.glm(lm.co.RPV.PAV, newdata = xy, type = "response"), 
                            nrow = grid.lines, ncol = grid.lines)
z.pred.co.RPV.SGV <- matrix(predict.glm(lm.co.RPV.SGV, newdata = xy, type = "response"), 
                            nrow = grid.lines, ncol = grid.lines)
z.pred.co.PAV.SGV <- matrix(predict.glm(lm.co.PAV.SGV, newdata = xy, type = "response"), 
                            nrow = grid.lines, ncol = grid.lines)
z.pred.co.any <- matrix(predict.glm(lm.co.any, newdata = xy, type = "response"), 
                        nrow = grid.lines, ncol = grid.lines)
z.pred.co.all <- matrix(predict.glm(lm.co.all, newdata = xy, type = "response"), 
                        nrow = grid.lines, ncol = grid.lines)
z.pred.atleast1 <- matrix(predict.glm(lm.atleast1, newdata = xy, type = "response"), 
                          nrow = grid.lines, ncol = grid.lines)
z.pred.rich <- matrix(predict.glm(lm.rich, newdata = xy, type = "response"), 
                      nrow = grid.lines, ncol = grid.lines)
z.pred.rich.no0 <- matrix(predict.glm(lm.rich.no0, newdata = xy, type = "response"), 
                          nrow = grid.lines, ncol = grid.lines)


###########################################################
################## 5) HEAT MAPS  ##########################
###########################################################

#################################################
## FIG 1: risk maps each alone & all together ###

setwd("C:/Users/straussa/Documents/Research/Seabloom-Borer Lab/Virus Diversity Rstar/3x virus/figures")
png("heatmaps.onescale.png", width = 6, height = 5, res = 600, units='in')

m <- rbind(c(1,1,2,2,3,3),c(4,4,4,5,5,5),c(4,4,4,5,5,5))
layout(m)
par(mar=c(1,2,2,2), oma=c(0,0,0,1), las=1) # save 500w x 400t
#### RPV ####
scatter2D(c(2,log(xy[["N"]]),7), c(0,log(xy[["P"]]),5),
          colvar=c(0,z.pred.RPV.sin,1), # .49
          cex=2, pch=15, xaxt="n", yaxt="n",
          xlim=log(c(7.5,375)), ylim=log(c(1,50)), xlab="", ylab="",
          colkey=list(side=4, at=c(0.1,0.5,.9)),
          col = ramp.col(c("white", "red", "black")))
axis(side=1, at=log(c(7.5,52.5,375)), cex.axis=1, labels=F)
axis(side=2, at=log(c(1,7,50)), cex.axis=1, labels=F)
#### SGV ####
scatter2D(c(2,log(xy[["N"]]),7), c(0,log(xy[["P"]]),5),
          colvar=c(0,z.pred.SGV.sin,1), #0.1
          cex=2, pch=15, xaxt="n", yaxt="n",
          xlim=log(c(7.5,375)), ylim=log(c(1,50)), xlab="", ylab="",
          colkey=list(side=4, at=c(0.1,0.5,0.9)),
          col = ramp.col(c("white", "red", "black")))
axis(side=1, at=log(c(7.5,52.5,375)), cex.axis=1, labels=F)
axis(side=2, at=log(c(1,7,50)), cex.axis=1, labels=F)
#### PAV ####
scatter2D(c(2,log(xy[["N"]]),7), c(0,log(xy[["P"]]),5),
          colvar=c(0,z.pred.PAV.sin,1), #.25
          cex=2, pch=15, xaxt="n", yaxt="n",
          xlim=log(c(7.5,375)), ylim=log(c(1,50)), xlab="", ylab="",
          colkey=list(side=4, at=c(0.1,0.5,0.9)),
          col = ramp.col(c("white", "red", "black")))
axis(side=1, at=log(c(7.5,52.5,375)), cex.axis=1, labels=F)
axis(side=2, at=log(c(1,7,50)), cex.axis=1, labels=F)
#### any single infection ####
par(mar=c(4,2,4,2))
scatter2D(c(2,log(xy[["N"]]),7), c(0,log(xy[["P"]]),5),
          colvar=c(0,z.pred.any.sin,1), #.45?
          cex=2, pch=15, xaxt="n", yaxt="n",
          xlim=log(c(7.5,375)), ylim=log(c(1,50)), xlab="", ylab="",
          colkey=list(side=4, at=c(0.1,0.5,0.9)),
          col = ramp.col(c("white", "red", "black")))
axis(side=1, at=log(c(7.5,52.5,375)), cex.axis=1, labels=F)
axis(side=2, at=log(c(1,7,50)), cex.axis=1, labels=F)
#### ACTUAL MULTIPLE INFECTIONS ####
scatter2D(c(2,log(xy[["N"]]),7), c(0,log(xy[["P"]]),5),
          colvar=c(0,z.pred.atleast1,1), #.5
          cex=2, pch=15, xaxt="n", yaxt="n",
          xlim=log(c(7.5,375)), ylim=log(c(1,50)), xlab="", ylab="",
          colkey=list(side=4, at=c(0.1,0.5,.9)),
          col = ramp.col(c("white", "red", "black")))
axis(side=1, at=log(c(7.5,52.5,375)), cex.axis=1, labels=F)
axis(side=2, at=log(c(1,7,50)), cex.axis=1, labels=F)

dev.off()


###############################################################
####################### 6) 3-D plots  #########################
###############################################################

# nice tutorial here:
# http://www.sthda.com/english/wiki/impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization

##################################################
## FIG 2: planes for single and co-inoculations ##

setwd("C:/Users/straussa/Documents/Research/Seabloom-Borer Lab/Virus Diversity Rstar/figures")
png("3D.png", width = 6, height = 3, res = 600, units='in')
par(mfrow=c(1,3),mar=c(1,1,0.5,0.5), oma=c(0,0.5,0,0)) # save 600 x 300

# RPV
scatter3D(log(RPV.sin$N), log(RPV.sin$P), RPV.sin$prev, colvar=NULL,
          pch=21, bg="purple", cex=2, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0,1.2),
          theta=30, phi=16, lty=3, 
          xlab="", ylab="", zlab="", nticks=3, ticktype="detailed", 
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.RPV.sin,  
                      facets = NA, col="purple"))
scatter3D(log(RPV.mix$N), log(RPV.mix$P), RPV.mix$prev, colvar=NULL, add=T,
          pch=22, bg="orange", cex=2, type="h", lty=3,
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.RPV.mix,
                      facets = NA, col="orange"))

# SGV
scatter3D(log(SGV.sin$N), log(SGV.sin$P), SGV.sin$prev, colvar=NULL,
          pch=21, bg="purple", cex=2, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0,1.2),
          theta=30, phi=16, lty=3, # CI=SGV.sin.CI
          xlab="", ylab="", zlab="", nticks=3, ticktype="detailed", 
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.SGV.sin,  
                      facets = NA, col="purple"))
scatter3D(log(RPV.sin$N), log(RPV.sin$P), SGV.mix$prev, colvar=NULL, add=T,
          pch=22, bg="orange", cex=2, type="h", lty=3, 
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.SGV.mix,  
                      facets = NA, col="orange"))

# PAV
scatter3D(log(RPV.sin$N), log(RPV.sin$P), PAV.sin$prev, colvar=NULL,
          pch=21, bg="purple", cex=2, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0,0.5),
          theta=30, phi=16, lty=3, # CI=PAV.sin.CI
          xlab="", ylab="", zlab="", nticks=3, ticktype="detailed", 
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.PAV.sin,  
                      facets = NA, col="purple"))
scatter3D(log(RPV.sin$N), log(RPV.sin$P), PAV.mix$prev, colvar=NULL, add=T,
          pch=22, bg="orange", cex=2, type="h", lty=3, 
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.PAV.mix,  
                      facets = NA, col="orange"))

dev.off()

###########################################################
# alternative 2D visualizaion as interaction plot

setwd("C:/Users/straussa/Documents/Research/Seabloom-Borer Lab/Virus Diversity Rstar/3x virus/figures")
png("2D.png", width = 8, height = 6, res = 600, units='in')

par(mfrow=c(2,3),mar=c(4,2.5,0.5,0.5), oma=c(0,0.5,0,0)) # save 600 x 300

plot(log(RPV.sin$N)+log(RPV.sin$P)/10, RPV.sin$prev, xaxt="n", yaxt="n", ann=F,
     xlim=c(1.5, 7), ylim=c(-.05,1.1), cex=0)
plotCI(x=log(RPV.sin$N)+log(RPV.sin$P)/10, y=RPV.sin$prev,
       uiw=RPV.sin$prev.se, add=T, err='y', sfrac=0, gap=0)
plotCI(x=log(RPV.mix$N)+log(RPV.mix$P)/10, y=RPV.mix$prev,
       uiw=RPV.mix$prev.se, add=T, err='y', sfrac=0, gap=0)
points(log(RPV.sin[RPV.sin$P==1,]$N)+log(RPV.sin[RPV.sin$P==1,]$P)/10, RPV.sin[RPV.sin$P==1,]$prev, 
       bg="purple", cex=2, pch=22, type="o", lty=2) 
points(log(RPV.sin[RPV.sin$P==7,]$N)+log(RPV.sin[RPV.sin$P==7,]$P)/10, RPV.sin[RPV.sin$P==7,]$prev, 
       bg="purple", cex=2, pch=24, type="o", lty=2) 
points(log(RPV.sin[RPV.sin$P==50,]$N)+log(RPV.sin[RPV.sin$P==50,]$P)/10, RPV.sin[RPV.sin$P==50,]$prev, 
       bg="purple", cex=2, pch=25, type="o", lty=2) 
points(log(RPV.mix[RPV.mix$P==1,]$N)+log(RPV.mix[RPV.mix$P==1,]$P)/10, RPV.mix[RPV.mix$P==1,]$prev, 
       bg="orange", cex=2, pch=22, type="o", lty=2) 
points(log(RPV.mix[RPV.mix$P==7,]$N)+log(RPV.mix[RPV.mix$P==7,]$P)/10, RPV.mix[RPV.mix$P==7,]$prev, 
       bg="orange", cex=2, pch=24, type="o", lty=2) 
points(log(RPV.mix[RPV.mix$P==50,]$N)+log(RPV.mix[RPV.mix$P==50,]$P)/10, RPV.mix[RPV.mix$P==50,]$prev, 
       bg="orange", cex=2, pch=25, type="o", lty=2) 
axis(side=1, at=c(2,4,6), cex.axis=1.2)
axis(side=2, at=c(0,.5,1), cex.axis=1.2)

plot(log(SGV.sin$N)+log(SGV.sin$P)/10, SGV.sin$prev, xaxt="n", yaxt="n", ann=F,
     xlim=c(1.5, 7), ylim=c(-.05,1.1), cex=0)
plotCI(x=log(SGV.sin$N)+log(SGV.sin$P)/10, y=SGV.sin$prev,
       uiw=SGV.sin$prev.se, add=T, err='y', sfrac=0, gap=0)
plotCI(x=log(SGV.mix$N)+log(SGV.mix$P)/10, y=SGV.mix$prev,
       uiw=SGV.mix$prev.se, add=T, err='y', sfrac=0, gap=0)
points(log(SGV.sin[SGV.sin$P==1,]$N)+log(SGV.sin[SGV.sin$P==1,]$P)/10, SGV.sin[SGV.sin$P==1,]$prev, 
       bg="purple", cex=2, pch=22, type="o", lty=2) 
points(log(SGV.sin[SGV.sin$P==7,]$N)+log(SGV.sin[SGV.sin$P==7,]$P)/10, SGV.sin[SGV.sin$P==7,]$prev, 
       bg="purple", cex=2, pch=24, type="o", lty=2) 
points(log(SGV.sin[SGV.sin$P==50,]$N)+log(SGV.sin[SGV.sin$P==50,]$P)/10, SGV.sin[SGV.sin$P==50,]$prev, 
       bg="purple", cex=2, pch=25, type="o", lty=2) 
points(log(SGV.mix[SGV.mix$P==1,]$N)+log(SGV.mix[SGV.mix$P==1,]$P)/10, SGV.mix[SGV.mix$P==1,]$prev, 
       bg="orange", cex=2, pch=22, type="o", lty=2) 
points(log(SGV.mix[SGV.mix$P==7,]$N)+log(SGV.mix[SGV.mix$P==7,]$P)/10, SGV.mix[SGV.mix$P==7,]$prev, 
       bg="orange", cex=2, pch=24, type="o", lty=2) 
points(log(SGV.mix[SGV.mix$P==50,]$N)+log(SGV.mix[SGV.mix$P==50,]$P)/10, SGV.mix[SGV.mix$P==50,]$prev, 
       bg="orange", cex=2, pch=25, type="o", lty=2) 
axis(side=1, at=c(2,4,6), cex.axis=1.2)
axis(side=2, at=c(0,.5,1), cex.axis=1.2)

plot(log(PAV.sin$N)+log(PAV.sin$P)/10, PAV.sin$prev, xaxt="n", yaxt="n", ann=F,
     xlim=c(1.5, 7), ylim=c(-0.05,0.5), cex=0)
plotCI(x=log(PAV.sin$N)+log(PAV.sin$P)/10, y=PAV.sin$prev,
       uiw=PAV.sin$prev.se, add=T, err='y', sfrac=0, gap=0)
plotCI(x=log(PAV.mix$N)+log(PAV.mix$P)/10, y=PAV.mix$prev,
       uiw=PAV.mix$prev.se, add=T, err='y', sfrac=0, gap=0)
points(log(PAV.sin[PAV.sin$P==1,]$N)+log(PAV.sin[PAV.sin$P==1,]$P)/10, PAV.sin[PAV.sin$P==1,]$prev, 
       bg="purple", cex=2, pch=22, type="o", lty=2) 
points(log(PAV.sin[PAV.sin$P==7,]$N)+log(PAV.sin[PAV.sin$P==7,]$P)/10, PAV.sin[PAV.sin$P==7,]$prev, 
       bg="purple", cex=2, pch=24, type="o", lty=2) 
points(log(PAV.sin[PAV.sin$P==50,]$N)+log(PAV.sin[PAV.sin$P==50,]$P)/10, PAV.sin[PAV.sin$P==50,]$prev, 
       bg="purple", cex=2, pch=25, type="o", lty=2) 
points(log(PAV.mix[PAV.mix$P==1,]$N)+log(PAV.mix[PAV.mix$P==1,]$P)/10, PAV.mix[PAV.mix$P==1,]$prev, 
       bg="orange", cex=2, pch=22, type="o", lty=2) 
points(log(PAV.mix[PAV.mix$P==7,]$N)+log(PAV.mix[PAV.mix$P==7,]$P)/10, PAV.mix[PAV.mix$P==7,]$prev, 
       bg="orange", cex=2, pch=24, type="o", lty=2) 
points(log(PAV.mix[PAV.mix$P==50,]$N)+log(PAV.mix[PAV.mix$P==50,]$P)/10, PAV.mix[PAV.mix$P==50,]$prev, 
       bg="orange", cex=2, pch=25, type="o", lty=2) 
axis(side=1, at=c(2,4,6), cex.axis=1.2)
axis(side=2, at=c(0,.2,.4), cex.axis=1.2)

plot(log(RPV.sin$P)+log(RPV.sin$N)/20, RPV.sin$prev, xaxt="n", yaxt="n", ann=F,
     xlim=c(-.5, 5), ylim=c(-.05,1.1), cex=0)
plotCI(x=log(RPV.sin$P)+log(RPV.sin$N)/20, y=RPV.sin$prev,
       uiw=RPV.sin$prev.se, add=T, err='y', sfrac=0, gap=0)
plotCI(x=log(RPV.mix$P)+log(RPV.mix$N)/20, y=RPV.mix$prev,
       uiw=RPV.mix$prev.se, add=T, err='y', sfrac=0, gap=0)
points(log(RPV.sin[RPV.sin$N==7.5,]$P)+log(RPV.sin[RPV.sin$N==7.5,]$N)/20, RPV.sin[RPV.sin$N==7.5,]$prev, 
       bg="purple", cex=2, pch=22, type="o", lty=2) 
points(log(RPV.sin[RPV.sin$N==52.5,]$P)+log(RPV.sin[RPV.sin$N==52.5,]$N)/20, RPV.sin[RPV.sin$N==52.5,]$prev, 
       bg="purple", cex=2, pch=24, type="o", lty=2) 
points(log(RPV.sin[RPV.sin$N==375,]$P)+log(RPV.sin[RPV.sin$N==375,]$N)/20, RPV.sin[RPV.sin$N==375,]$prev, 
       bg="purple", cex=2, pch=25, type="o", lty=2) 
points(log(RPV.mix[RPV.mix$N==7.5,]$P)+log(RPV.mix[RPV.mix$N==7.5,]$N)/20, RPV.mix[RPV.mix$N==7.5,]$prev, 
       bg="orange", cex=2, pch=22, type="o", lty=2) 
points(log(RPV.mix[RPV.mix$N==52.5,]$P)+log(RPV.mix[RPV.mix$N==52.5,]$N)/20, RPV.mix[RPV.mix$N==52.5,]$prev, 
       bg="orange", cex=2, pch=24, type="o", lty=2) 
points(log(RPV.mix[RPV.mix$N==375,]$P)+log(RPV.mix[RPV.mix$N==375,]$N)/20, RPV.mix[RPV.mix$N==375,]$prev, 
       bg="orange", cex=2, pch=25, type="o", lty=2) 
axis(side=1, at=c(0,2,4), cex.axis=1.2)
axis(side=2, at=c(0,.5,1), cex.axis=1.2)

plot(log(SGV.sin$P)+log(SGV.sin$N)/20, SGV.sin$prev, xaxt="n", yaxt="n", ann=F,
     xlim=c(-.5, 5), ylim=c(-.05,1.1), cex=0)
plotCI(x=log(SGV.sin$P)+log(SGV.sin$N)/20, y=SGV.sin$prev,
       uiw=SGV.sin$prev.se, add=T, err='y', sfrac=0, gap=0)
plotCI(x=log(SGV.mix$P)+log(SGV.mix$N)/20, y=SGV.mix$prev,
       uiw=SGV.mix$prev.se, add=T, err='y', sfrac=0, gap=0)
points(log(SGV.sin[SGV.sin$N==7.5,]$P)+log(SGV.sin[SGV.sin$N==7.5,]$N)/20, SGV.sin[SGV.sin$N==7.5,]$prev, 
       bg="purple", cex=2, pch=22, type="o", lty=2) 
points(log(SGV.sin[SGV.sin$N==52.5,]$P)+log(SGV.sin[SGV.sin$N==52.5,]$N)/20, SGV.sin[SGV.sin$N==52.5,]$prev, 
       bg="purple", cex=2, pch=24, type="o", lty=2) 
points(log(SGV.sin[SGV.sin$N==375,]$P)+log(SGV.sin[SGV.sin$N==375,]$N)/20, SGV.sin[SGV.sin$N==375,]$prev, 
       bg="purple", cex=2, pch=25, type="o", lty=2) 
points(log(SGV.mix[SGV.mix$N==7.5,]$P)+log(SGV.mix[SGV.mix$N==7.5,]$N)/20, SGV.mix[SGV.mix$N==7.5,]$prev, 
       bg="orange", cex=2, pch=22, type="o", lty=2) 
points(log(SGV.mix[SGV.mix$N==52.5,]$P)+log(SGV.mix[SGV.mix$N==52.5,]$N)/20, SGV.mix[SGV.mix$N==52.5,]$prev, 
       bg="orange", cex=2, pch=24, type="o", lty=2) 
points(log(SGV.mix[SGV.mix$N==375,]$P)+log(SGV.mix[SGV.mix$N==375,]$N)/20, SGV.mix[SGV.mix$N==375,]$prev, 
       bg="orange", cex=2, pch=25, type="o", lty=2)
axis(side=1, at=c(0,2,4), cex.axis=1.2)
axis(side=2, at=c(0,.5,1), cex.axis=1.2)

plot(log(PAV.sin$P)+log(PAV.sin$N)/20, PAV.sin$prev, xaxt="n", yaxt="n", ann=F,
     xlim=c(-.5, 5), ylim=c(-.05,0.5), cex=0)
plotCI(x=log(PAV.sin$P)+log(PAV.sin$N)/20, y=PAV.sin$prev,
       uiw=PAV.sin$prev.se, add=T, err='y', sfrac=0, gap=0)
plotCI(x=log(PAV.mix$P)+log(PAV.mix$N)/20, y=PAV.mix$prev,
       uiw=PAV.mix$prev.se, add=T, err='y', sfrac=0, gap=0)
points(log(PAV.sin[PAV.sin$N==7.5,]$P)+log(PAV.sin[PAV.sin$N==7.5,]$N)/20, PAV.sin[PAV.sin$N==7.5,]$prev, 
       bg="purple", cex=2, pch=22, type="o", lty=2) 
points(log(PAV.sin[PAV.sin$N==52.5,]$P)+log(PAV.sin[PAV.sin$N==52.5,]$N)/20, PAV.sin[PAV.sin$N==52.5,]$prev, 
       bg="purple", cex=2, pch=24, type="o", lty=2) 
points(log(PAV.sin[PAV.sin$N==375,]$P)+log(PAV.sin[PAV.sin$N==375,]$N)/20, PAV.sin[PAV.sin$N==375,]$prev, 
       bg="purple", cex=2, pch=25, type="o", lty=2) 
points(log(PAV.mix[PAV.mix$N==7.5,]$P)+log(PAV.mix[PAV.mix$N==7.5,]$N)/20, PAV.mix[PAV.mix$N==7.5,]$prev, 
       bg="orange", cex=2, pch=22, type="o", lty=2) 
points(log(PAV.mix[PAV.mix$N==52.5,]$P)+log(PAV.mix[PAV.mix$N==52.5,]$N)/20, PAV.mix[PAV.mix$N==52.5,]$prev, 
       bg="orange", cex=2, pch=24, type="o", lty=2) 
points(log(PAV.mix[PAV.mix$N==375,]$P)+log(PAV.mix[PAV.mix$N==375,]$N)/20, PAV.mix[PAV.mix$N==375,]$prev, 
       bg="orange", cex=2, pch=25, type="o", lty=2)
axis(side=1, at=c(0,2,4), cex.axis=1.2)
axis(side=2, at=c(0,.2,.4), cex.axis=1.2)

dev.off()


############################################################
## FIG 3: prevalence and diversity in co-inoculated hosts ##

# confidence intervals:
rich.CI <- list(z = matrix(nrow = length(coin.sum$v.rich.se),
                           data = rep(coin.sum$v.rich.se,2)))
rich.no0.CI <- list(z = matrix(nrow = length(coin.sum$v.rich.no0.se),
                               data = rep(coin.sum$v.rich.no0.se,2)))
atleast1.CI <- list(z = matrix(nrow = length(coin.sum$prev.atleast1.se),
                               data = rep(coin.sum$prev.atleast1.se,2)))

setwd("C:/Users/straussa/Documents/Research/Seabloom-Borer Lab/Virus Diversity Rstar/figures")
png("3D.div.png", width = 8, height = 3, res = 600, units='in')
par(mfrow=c(1,3),mar=c(1,3.5,0.5,0.5), oma=c(0,0.5,0,0))

# at least one infection:
scatter3D(log(coin.sum$N), log(coin.sum$P), coin.sum$prev.atleast1, colvar=NULL,
          pch=23, bg="green3", cex=1, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0.5,1.05),
          theta=30, phi=25, lty=3, CI=atleast1.CI,
          xlab="", ylab="", zlab="", nticks=3, ticktype="detailed", 
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.atleast1,  
                      facets = NA, col="darkgreen"))

# diversity including the 0's:
scatter3D(log(coin.sum$N), log(coin.sum$P), coin.sum$v.rich, colvar=NULL,
          pch=23, bg="green3", cex=1, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0.5,2.1),
          theta=30, phi=25, lty=3, CI=rich.CI,
          xlab="", ylab="", zlab="", nticks=3, ticktype="detailed", 
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.rich,  
                      facets = NA, col="darkgreen"))

# diversity excluding the 0's:
scatter3D(log(coin.sum$N), log(coin.sum$P), coin.sum$v.rich.no0, colvar=NULL,
          pch=23, bg="green3", cex=1, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0.5,2.1),
          theta=30, phi=25, lty=3, CI=rich.no0.CI,
          xlab="", ylab="", zlab="", nticks=3, ticktype="detailed", 
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.rich.no0,  
                      facets = NA, col="darkgreen"))

dev.off()


################################################
## FIG SX planes for 5 types of coinfections ##

setwd("C:/Users/straussa/Documents/Research/Seabloom-Borer Lab/Virus Diversity Rstar/figures")
png("3D.co.png", width = 6, height = 4, res = 600, units='in')
par(mfrow=c(2,3),mar=c(1,1,3,0), oma=c(0,0,0,0)) #to make plots fill out space better
# save 600 x 400

scatter3D(log(coin.sum$N), log(coin.sum$P), coin.sum$prev.any, colvar=NULL,
          pch=23, bg="maroon4", cex=2, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0,1.1),
          theta=30, phi=16, lty=3, 
          xlab="", ylab="", zlab="", ticktype="detailed", nticks=3,
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.co.any,  
                      facets = NA, col="maroon4"))
scatter3D(log(coin.sum$N), log(coin.sum$P), coin.sum$prev.RPV.SGV, colvar=NULL,
          pch=23, bg="maroon4", cex=2, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0,1.1),
          theta=30, phi=16, lty=3, 
          xlab="", ylab="", zlab="", ticktype="detailed", nticks=3,
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.co.RPV.SGV,  
                      facets = NA, col="maroon4"))
scatter3D(log(coin.sum$N), log(coin.sum$P), coin.sum$prev.PAV.SGV, colvar=NULL,
          pch=23, bg="maroon4", cex=2, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0,0.22),
          theta=30, phi=16, lty=3, 
          xlab="", ylab="", zlab="", ticktype="detailed", nticks=3,
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.co.PAV.SGV,  
                      facets = NA, col="maroon4"))
scatter3D(log(coin.sum$N), log(coin.sum$P), coin.sum$prev.RPV.PAV, colvar=NULL,
          pch=23, bg="maroon4", cex=2, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0,0.22),
          theta=30, phi=16, lty=3, 
          xlab="", ylab="", zlab="", ticktype="detailed", nticks=3,
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.co.RPV.PAV,  
                      facets = NA, col="maroon4"))
scatter3D(log(coin.sum$N), log(coin.sum$P), coin.sum$prev.all, colvar=NULL,
          pch=23, bg="maroon4", cex=2, bty='b2', type="h",
          xlim=c(log(7.5)-0.5,log(375)+0.5), ylim=c(log(1)-0.5,log(50)+0.5), zlim=c(0,0.13),
          theta=30, phi=16, lty=3, 
          xlab="", ylab="", zlab="", ticktype="detailed", nticks=3,
          surf = list(x = log(x.pred), y = log(y.pred), z = z.pred.co.all,  
                      facets = NA, col="maroon4"))

dev.off()


######################################################################
################### 7) 2D PLOTS ALONG NP RATIO  ######################
######################################################################

setwd("C:/Users/straussa/Documents/Research/Seabloom-Borer Lab/Virus Diversity Rstar/figures")
png("2D-NP.png", width = 6, height = 6, res = 600, units='in')
par(mfrow=c(3,2),mar=c(0,0.5,1,1), oma=c(2.5,2,0,0)) 
# save 600 x 600

# prevalence of RPV over NP, with single (purple) and mixed (orange):
plot(log(prev.sum$NP), prev.sum$prev, col = "white", xlab = "", ylab = "", 
     ylim = c(0,1.2), xlim = c(-2.2,6.5), cex=2, pch=22, bg="white", xaxt="n", yaxt="n") 
axis(side=2, at=c(0,.5,1), cex.axis=1.2)
plotCI(x=log(RPV.sin$NP), y=RPV.sin$prev, uiw=RPV.sin$prev.se, add=T, err='y', sfrac=0, gap=0)
plotCI(x=log(RPV.mix$NP), y=RPV.mix$prev, uiw=RPV.mix$prev.se, add=T, err='y', sfrac=0, gap=0)
points(log(RPV.sin$NP), RPV.sin$prev, cex=2, pch=24, bg="purple")
points(log(RPV.mix$NP), RPV.mix$prev, cex=2, pch=25, bg="orange")

# prevalence of co any over NP
plot(log(coin.sum$NP), coin.sum$prev, col = "white", xlab = "", ylab = "", 
     ylim = c(0,1.2), xlim = c(-2.2,6.5), cex=2, pch=22, bg="white", xaxt="n", yaxt="n") 
#axis(side=2, at=c(0,.5,1), cex.axis=1.2)
plotCI(x=log(coin.sum$NP), y=coin.sum$prev.any, uiw=coin.sum$prev.any.se, 
       add=T, err='y', sfrac=0, gap=0)
points(log(coin.sum$NP), coin.sum$prev.any, cex=2, pch=23, bg="maroon4")

# prevalence of SGV over NP, with single (blue) and mixed (orange):
plot(log(prev.sum$NP), prev.sum$prev, col = "white", xlab = "", ylab = "", 
     ylim = c(0,1.2), xlim = c(-2.2,6.5), cex=2, pch=22, bg="white", xaxt="n", yaxt="n") 
axis(side=2, at=c(0,.5,1), cex.axis=1.2)
plotCI(x=log(SGV.sin$NP), y=SGV.sin$prev, uiw=SGV.sin$prev.se, add=T, err='y', sfrac=0, gap=0)
plotCI(x=log(SGV.mix$NP), y=SGV.mix$prev, uiw=SGV.mix$prev.se, add=T, err='y', sfrac=0, gap=0)
points(log(SGV.sin$NP), SGV.sin$prev, cex=2, pch=24, bg="purple")
points(log(SGV.mix$NP), SGV.mix$prev, cex=2, pch=25, bg="orange")

# prevalence of co SGV+RPV over NP
plot(log(coin.sum$NP), coin.sum$prev, col = "white", xlab = "", ylab = "", 
     ylim = c(0,1.2), xlim = c(-2.2,6.5), cex=2, pch=22, bg="white", xaxt="n", yaxt="n") 
#axis(side=2, at=c(0,.5,1), cex.axis=1.2)
plotCI(x=log(coin.sum$NP), y=coin.sum$prev.RPV.SGV, uiw=coin.sum$prev.RPV.SGV.se, 
       add=T, err='y', sfrac=0, gap=0)
points(log(coin.sum$NP), coin.sum$prev.RPV.SGV, cex=2, pch=23, bg="maroon4")

# prevalence of PAV over NP, with single (blue) and mixed (orange):
plot(log(prev.sum$NP), prev.sum$prev, col = "white", xlab = "", ylab = "", 
     ylim = c(0,0.55), xlim = c(-2.2,6.5), cex=2, pch=22, bg="white", xaxt="n", yaxt="n") 
axis(side=2, at=c(0,.2,0.4), cex.axis=1.2)
axis(side=1, at=c(-2,0,2,4,6), cex.axis=1.2)
plotCI(x=log(PAV.sin$NP), y=PAV.sin$prev, uiw=PAV.sin$prev.se, add=T, err='y', sfrac=0, gap=0)
plotCI(x=log(PAV.mix$NP), y=PAV.mix$prev, uiw=PAV.mix$prev.se, add=T, err='y', sfrac=0, gap=0)
points(log(PAV.sin$NP), PAV.sin$prev, cex=2, pch=24, bg="purple")
points(log(PAV.mix$NP), PAV.mix$prev, cex=2, pch=25, bg="orange")

# prevalence of co all over NP
plot(log(coin.sum$NP), coin.sum$prev, col = "white", xlab = "", ylab = "", 
     ylim = c(0,0.55), xlim = c(-2.2,6.5), cex=2, pch=22, bg="white", xaxt="n", yaxt="n") 
#axis(side=2, at=c(0,.1,.2), cex.axis=1.2)
axis(side=1, at=c(-2,0,2,4,6), cex.axis=1.2)
plotCI(x=log(coin.sum$NP), y=coin.sum$prev.all, uiw=coin.sum$prev.all.se, 
       add=T, err='y', sfrac=0, gap=0)
points(log(coin.sum$NP), coin.sum$prev.all, cex=2, pch=23, bg="maroon4")

dev.off()
