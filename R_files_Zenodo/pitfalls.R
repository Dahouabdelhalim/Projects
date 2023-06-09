###########################
#### Alex's approach:  #####
############################


setwd("C:/Users/straussa/Documents/Research/Coauthor Projects/Anita Porath-Krause/pitfalls/figures")
library(vegan) 
library(gplots)


#############################################
# 1) fit empirical species distribution     #
#     and simulate example communities with #
#     different degrees of skew             #
#############################################

# empirical data
setwd("C:/Users/straussa/Documents/Research/Coauthor Projects/Anita Porath-Krause/pitfalls/empirical data")
data <- read.csv("msb-diversity-4sites-otu-sums.csv")
hist(data$value)
emp_n_ind <- sum(data$value) # number of individuals
emp_n_sp <- nrow(data) # number of species

# https://rdrr.io/rforge/vegan/man/radfit.html
emp.fit <- radfit(x=data$value, family=poisson); emp.fit
emp.fit
# okay, so Mandelbrot fits best, and requires 3 params: c, gamma, beta
emp.fit$models$Mandelbrot$coefficients
emp_c <- emp.fit$models$Mandelbrot$coefficients[1]
emp_gamma <- emp.fit$models$Mandelbrot$coefficients[2]
emp_beta <- emp.fit$models$Mandelbrot$coefficients[3]

# the emperical distribution, as fit by the Mandelbrot
# found the equation for abundance of x on the radfit help page
plot(seq(1,emp_n_sp,1), emp_n_ind*emp_c*(seq(1,emp_n_sp,1)+emp_beta)^emp_gamma)
# on a log x scale
plot(log(seq(1,emp_n_sp,1)), emp_n_ind*emp_c*(seq(1,emp_n_sp,1)+emp_beta)^emp_gamma)

n_com <- 1 # level of replication. 
n_ind <- 1000000 # number of individuals in each communiy. setting to 10,000 for all. 
n_sp <-  30 #nrow(data) # number of species in the community. set to 10, 100, or 1000
n_sam <- 10 # number of individuals to sample. set to 100, 200, 500, 1000, 2500, 5000, 7500, 10000

# three different extremes of skew:
p_skew <- (seq(1,emp_n_sp,1)+emp_beta)^(emp_gamma*2)
p_mod <- (seq(1,emp_n_sp,1)+emp_beta)^(emp_gamma)
p_grad <- (seq(1,emp_n_sp,1)+emp_beta)^(emp_gamma/2)
# and three simulated communities:
ext_sam <- rmultinom(n=n_com, size=n_ind, prob=p_skew[seq(1, emp_n_sp, emp_n_sp/n_sp)])
mod_sam <- rmultinom(n=n_com, size=n_ind, prob=p_mod[seq(1, emp_n_sp, emp_n_sp/n_sp)])
grad_sam <- rmultinom(n=n_com, size=n_ind, prob=p_grad[seq(1, emp_n_sp, emp_n_sp/n_sp)])

png("eg.sp.dist.png", width = 6, height = 2, res = 600, units='in')

par(mfrow=c(1,3), mar=c(0,0,1,1), oma=c(2.5,2.5,0,0), las=1)

plot(log(ext_sam+1), ylim=c(0,15), type="h", xaxt="n", yaxt="n", ann=F)
axis(side=2, at=c(0,5,10,15))
axis(side=1, at=c(0,10,20,30))

plot(log(mod_sam+1), ylim=c(0,15), type="h", xaxt="n", yaxt="n", ann=F)
axis(side=1, at=c(0,10,20,30))

plot(log(grad_sam+1), ylim=c(0,15), type="h", xaxt="n", yaxt="n", ann=F)
axis(side=1, at=c(0,10,20,30))

dev.off()


###########################################
## 2) simulate communities that vary in   #
#     skew and richness, subsample, and   #
#     see when diversity indices converge #
###########################################

# number of individuals in each communiy. 
# Empirical data has n_ind=13M, so we'll try 1M, 13M, and 26M
n_ind <- 1000000
#n_ind <- 13000000
#n_ind <- 26000000

n_com <- 10 # number of communities to simulate. unit of replication for error bars.
# levels of replication depends on size of community (for RAM for my laptop.. 
# or rewrite everything to be parallel processed)
# 1M: n_com=100
# 13M: n_com=20
# 26M: n_com=10

abund.rar.num <- 1000 #number of individuals to sample for rarefactions


##### TO EXTEND ALL THE WAY TO FULLY SAMPLED COMMUNITIES:
# p.list <- data.frame(n_sam=rep(c(100,200,500,1000,2500,5000,7500,10000,20000,
#                                  50000,100000,200000,500000,1000000),9)) # for n_ind=1M 
# p.list <- data.frame(n_sam=rep(c(100,200,500,1000,2500,5000,7500,10000,20000,
#                                  50000,100000,200000,500000,1000000,
#                                  2000000,5000000,10000000,13000000),9)) # for n_ind=13M
# p.list <- data.frame(n_sam=rep(c(100,200,500,1000,2500,5000,7500,10000,20000,
#                                  50000,100000,200000,500000,1000000,
#                                  2000000,5000000,10000000,20000000,26000000),9)) # for n_ind=26M

##### TO EXTEND ONLY TO 10,000 indidivuals samples (for Fig3):
 p.list <- data.frame(n_sam=rep(c(100,200,500,1000,2500,5000,7500,10000),9))  

sam_levs <- length(unique(p.list$n_sam))
p.list$n_sp = c(rep(30,3*sam_levs), rep(300,3*sam_levs), rep(3000,3*sam_levs))
p.list$skew <- rep(c(rep(emp_gamma*2,sam_levs), rep(emp_gamma,sam_levs), rep(emp_gamma/2,sam_levs)),3)

for(i in 1:nrow(p.list)){
  print(i)
  
  #1) simulate full communities
  mcom <- rmultinom(n=n_com, size=n_ind,  
                    prob=(seq(1, emp_n_sp+0.5, emp_n_sp/p.list$n_sp[i])+emp_beta)^(p.list$skew[i]))
  # need to add the 0.5 because otherwise off by 1 for # sp (make up for starting 1 not 0)
  
  #2) calculate diversity indices for full community (to plot as asymtotes):
  p.list$rich_full[i] <- median(apply(mcom, 2, function(x) length(x[x>0])))
  p.list$simp_full[i] <- median(diversity(mcom, index="invsimpson", MARGIN=2))
  p.list$shan_full[i] <- median(diversity(mcom, index="shannon", MARGIN=2))
  
  #3) sub-sample 
  rownames(mcom) <- paste0("sp",1:nrow(mcom))
  mcom_long <- apply(mcom, 2, function(x) rep(rownames(mcom), x))
  mcom_samp <- apply(mcom_long, 2, function(x) sample(x, size=p.list$n_sam[i], replace=F))
  
  #4) calculate richness of samples
  richness <- apply(mcom_samp, 2, function(x) length(unique(x))) # richness 
  p.list$rich_med[i] <- median(richness)
  p.list$rich_2.5[i] <- quantile(richness, probs=seq(0.025,0.975,0.95))[1]
  p.list$rich_97.5[i] <- quantile(richness, probs=seq(0.025,0.975,0.95))[2]
  
  #5) calculate diversity indices of samples:
  mcom_samp_short <- apply(mcom_samp, 2, function(x) table(factor(x, levels = paste0("sp",1:p.list$n_sp[i]))))
  simp <- diversity(mcom_samp_short, index="invsimpson", MARGIN=2) # margin 2 says columns
  p.list$simp_med[i] <- median(simp)
  p.list$simp_2.5[i] <- quantile(simp, probs=seq(0.025,0.975,0.95))[1]
  p.list$simp_97.5[i] <- quantile(simp, probs=seq(0.025,0.975,0.95))[2]
  shan <- diversity(mcom_samp_short, index="shannon", MARGIN=2) # margin 2 says columns
  p.list$shan_med[i] <- median(shan)
  p.list$shan_2.5[i] <- quantile(shan, probs=seq(0.025,0.975,0.95))[1]
  p.list$shan_97.5[i] <- quantile(shan, probs=seq(0.025,0.975,0.95))[2]
  rar <- ifelse(p.list$n_sam[i] < abund.rar.num, NA, # only do the rarifaction if 
                rarefy(mcom_samp_short, abund.rar.num, MARGIN=2)) # margin 2 says columns
  p.list$rar_med[i] <- median(rar, na.rm=T)
  p.list$rar_2.5[i] <- quantile(rar, probs=seq(0.025,0.975,0.95), na.rm=T)[1]
  p.list$rar_97.5[i] <- quantile(rar, probs=seq(0.025,0.975,0.95), na.rm=T)[2]
  
  #6) calculate composition distance indices between sample and full:
  # not sure how to vectorize, so making my own loop:
  distances <- data.frame(com=seq(1,n_com,1))
  for(j in 1:n_com){
    spsite <- merge(mcom[,j], mcom_samp_short[,j], by=0, all=T)
    distances$bray[j] <- vegdist(wisconsin(t(spsite[,-1])), method="bray")
    distances$jac[j] <- vegdist(t(spsite[,-1]), method="jaccard", binary=TRUE)
  }
  p.list$bray_med[i] <- 1-median(distances$bray)
  p.list$bray_2.5[i] <- 1-quantile(distances$bray, probs=seq(0.025,0.975,0.95))[1]
  p.list$bray_97.5[i] <- 1-quantile(distances$bray, probs=seq(0.025,0.975,0.95))[2]
  p.list$jac_med[i] <- 1-median(distances$jac)
  p.list$jac_2.5[i] <- 1-quantile(distances$jac, probs=seq(0.025,0.975,0.95))[1]
  p.list$jac_97.5[i] <- 1-quantile(distances$jac, probs=seq(0.025,0.975,0.95))[2]
  
  }

p.list


#write.csv(p.list, file = "pitfalls.1M.100%.100reps.csv")
#write.csv(p.list, file = "pitfalls.13M.100%.20reps.csv")
#write.csv(p.list, file = "pitfalls.26M.100%.10reps.csv")

write.csv(p.list, file = "pitfalls.1M.to10K.100reps.csv")
write.csv(p.list, file = "pitfalls.13M.to10K.20reps.csv")
write.csv(p.list, file = "pitfalls.26M.to10K.10reps.csv")
