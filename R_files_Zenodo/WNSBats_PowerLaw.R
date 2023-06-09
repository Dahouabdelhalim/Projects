##Multi-scale model of WNS

##Two-part script: Part 1 fits distributions to the cave census data from NY, Turner et al 2011 
##Part 2 fits distribtuions to dispersal data from Nealy, found no reason to differ from exponential tail. Effectively uniform dispersal when ignoring philopatry

##Load library
library(poweRlaw)
library(fitdistrplus)

##---Part 1
##Data from Turner et al 2011, transcribed by A. Kramer
pop.data<-read.csv("NY lucifugus populations from Turner et al 2011.csv")
pops<-pop.data$PreWNS.count
#Change zeros to 1
pops<-ifelse(pops==0,1,pops)
saveRDS(pops,"NY_Ml_pops.rds")

#Fit alternative distributions to PreWNS.count
##Continuous power law object
cpl_pop=conpl$new(pops) 
est=estimate_xmin(cpl_pop) #estimate the cutoff
cpl_pop$setXmin(est)
bs_p <- bootstrap_p(cpl_pop, no_of_sims = 5000, threads = 8, seed = 5) #Could be powerlaw

#Continuous log normal object
cln_pop=conlnorm$new(pops)
est=estimate_xmin(cln_pop)
cln_pop$setXmin(est)
bs_ln <- bootstrap_p(cln_pop, no_of_sims = 5000, threads = 8, seed = 6) #Could be powerlaw

#Continuous exponential object
cex_pop=conexp$new(pops)
est=estimate_xmin(cex_pop)
cex_pop$setXmin(est)
bs_exp <- bootstrap_p(cex_pop, no_of_sims = 5000, threads = 8, seed = 6) #Could be lognormal, GOF prefers this over power law

##Visually examine the dispersal data and 3 potential distributions
plot(cpl_pop)
lines(cpl_pop,col="green")
lines(cln_pop,col="blue")
lines(cex_pop,col="red")
#Appears to be clear support for log-normal or power law, but tests don't rule out exponential, why?



##Do comparisons between the models as in package documentation, these likelihood ratio tests are not useful because we don't care about number of parameters, only want best fit
cln_pop$setXmin(cex_pop$getXmin()) #Set lognormal xmin to power law xmin
est=estimate_pars(cln_pop)
cln_pop$setPars(est)
comp_ln_pl=compare_distributions(cpl_pop,cln_pop)

cex_pop$setXmin(cex_pop$getXmin())
est=estimate_pars(cex_pop)
cex_pop$setPars(est)
comp_ex_pl=compare_distributions(cpl_pop,cex_pop)

comp_ex_ln=compare_distributions(cln_pop,cex_pop)
##Despite different fits to data there is not statistical significance



##Decided to try lognormal with "fitdistrplus" package
plotdist(pops,hist=TRUE,demp=TRUE)
descdist(pops, boot=1000) #This suggests beta distribution
lnfit <- fitdist(pops, "lnorm")
expfit <- fitdist(pops, "exp") #Can't fit
gammafit <- fitdist(pops, "gamma") #Can't fit
weifit <- fitdist(pops, "weibull") # Estimated but warnings

plot.legend <- c("lognormal", "Weibull")
denscomp(list(lnfit, weifit), legendtext = plot.legend)
qqcomp(list(lnfit, weifit), legendtext = plot.legend)
cdfcomp(list(lnfit, weifit), legendtext = plot.legend)
ppcomp(list(lnfit, weifit), legendtext = plot.legend)

##Goodness of fit statistics
gofstat(list(lnfit, weifit),fitnames = c("lnorm", "weibull")) #Slight reference for weibull in GOF

#Use Weibull with these parameters to generate population sizes
randCavePops<-function(n,distrib="spline"){
  if(distrib=="weibull") cave.pops<-rweibull(n, shape = 0.3155586, scale = 1508.5363225)
  if(distrib=="lnorm") cave.pops<-rlnorm(n,meanlog = 5.596804,sdlog=3.522228)
  if(distrib=="spline"){
    pops<-readRDS("NY_Ml_pops.rds")
    spline.est<-logspline(pops, lbound=1, ubound=200000)
    cave.pops<-rlogspline(n, spline.est)
  }
  return(cave.pops)
}



##Script for powerlaw distributed dispersal
##Don't have strong support for any of these, assume uniform distribution in the tail

##Load library
library(poweRlaw)

#Data from N..ly
data<-c(rep(25,15),rep(75,13),rep(225,6),rep(275,13),rep(325,24),rep(425,2),rep(525,3),rep(575,16))

##Continuous power law object
cpl_bat=conpl$new(data) 
est=estimate_xmin(cpl_bat) #estimate the cutoff
cpl_bat$setXmin(est)

#Continuous log normal object
cln_bat=conlnorm$new(data)
est=estimate_xmin(cln_bat)
cln_bat$setXmin(est)

#Continuous exponential object
cex_bat=conexp$new(data)
est=estimate_xmin(cex_bat)
cex_bat$setXmin(est)

##Visually examine the dispersal data and 3 potential distributions
plot(cpl_bat)
lines(cpl_bat,col="green")
lines(cln_bat,col="blue")
lines(cex_bat,col="red")

##Do comparisons between the models as in package documentation
cln_bat$setXmin(cpl_bat$getXmin()) #Set lognormal xmin to power law xmin
est=estimate_pars(cln_bat)
cln_bat$setPars(est)
comp_ln_pl=compare_distributions(cpl_bat,cln_bat)

cex_bat$setXmin(cpl_bat$getXmin())
est=estimate_pars(cex_bat)
cex_bat$setPars(est)
comp_ex_pl=compare_distributions(cpl_bat,cex_bat)

comp_ex_ln=compare_distributions(cpl_bat,cex_bat)

##No support for distribution other than exponential
##Export the object to use for random number generation in the simulation code
saveRDS(cex_bat,"ExponentialDispersal.Rds")

##Code for the simulation
exp_bat_disp<-readRDS("ExponentialDispersal.Rds")
##Use dist_rand to get random numbers from the distribution, this only generates the tail
##Second number is number of distances needed. Then get cave closest to that distance
distances<-dist_rand(exp_bat_disp,10)



