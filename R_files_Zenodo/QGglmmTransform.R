# load packages
library(QGglmm)
library(parallel)
library(MCMCglmm)
library(coda)
library(purrr)
library(repurrrsive)
library(ggplot2)
library(dplyr)
library(tictoc)

# load data
load("cebus.RData")

# load previously run models
load("Models1_fullModels.RData")
load("Models1_reducedModels.RData")

# load qgglmm objects
load("QGglmmTransform.RData")

################################################################################
################################################################################
################################################################################

### Calculate Intra-class correlation coefficients (posterior distributions)
### ICC.obs will be for each random effect from first chain of each model

# value for n.obs parameter in QGicc() and QGparams()
# this in needed for multinomial models (binomN.logit)
# no "right" answer for what to use, so will do with the average from the dataset
obs <- round(mean(cebus$n)) ### mean #obs per subject per month (i.e., 32)

# function for generating posterior distributions for ICC estimates
## multinomial2 models, 
## not accounting for fixed effects variance in denominator
genICC_M1 <- function(m){  #m=model, col=variance component
  mname <- deparse(substitute(m))  #store name of model
  rlist <- names(posterior.mode(m$VCV))  #store list of random effects
  (print(mname))
  tmp_df <- data.frame()
  for (r in rlist){
    (print(r))
    tmp <- pmap_dfr(list(mu       = m[["Sol"]][ , "(Intercept)"],
                         var.comp = m[["VCV"]][ , r],
                         var.p    = rowSums(m[["VCV"]])),
                    QGicc,
                    model = "binomN.logit",
                    n.obs = obs,
                    verbose = FALSE)
    tmp$component <- r
    tmp_df <- rbind(tmp_df,tmp)
  }
  tmp_df$model <- mname
  tmp_df$FE.var <- "No"
  tmp_df$model.type <- "binomN.logit"
  tmp_df
}

# use function on first chain of models
M_ICC1   <- genICC_M1(m09a[[1]])
M_ICC1.0 <- genICC_M1(m09a.0[[1]])
M_ICC1.1 <- genICC_M1(m09a.1[[1]])
M_ICC1.2 <- genICC_M1(m09a.2[[1]])
M_ICC1.3 <- genICC_M1(m09a.3[[1]])
M_ICC1.4 <- genICC_M1(m09a.4[[1]])

## multinomial2 models, accounting for fixed effects variance in denominator
### this will take a loooong time, because their is no closed form solution
## it would take approximately 9 days to integrate over the entire posterior!!
## and that for just ONE variance component. Our models have 10!
## so, we will run QGicc on a random sample of 1000 iterations instead of 4000
### and we will run these on multiple cores to speed things up
set.seed(2002)   #28 Days Later :p
ilst <- sort(sample.int(4000, 1000, replace=T))  #random sample index
mname <- deparse(substitute(m09a[[1]]))  #store name of model
rlist <- names(posterior.mode(m$VCV))    #store list of random effects
tic()
M_ICC2x <- mclapply(1:length(rlist), function(i) {
	pmap_dfr(list(predict  = predictF[ilst],
								var.comp = m[["VCV"]][ , rlist[i]][ilst],
								var.p    = rowSums(m[["VCV"]])[ilst]),
					 QGicc,
					 model = "binomN.logit",
					 n.obs = obs, 
					 verbose = FALSE)
}, mc.cores=5)
toc()  #379176.42 sec elapsed
M_ICC2 <- mapply(cbind, 
								 M_ICC2x, 
								 "component"=rlist, 
								 model=mname, 
								 SIMPLIFY=F)

################################################################################
################################################################################
################################################################################

# function for generating posterior distributions for ICC estimates
## poisson models, 
## not accounting for fixed effects variance in denominator
genICC_P1 <- function(m){  #m=model
  mname <- deparse(substitute(m))  #store name of model
  rlist <- names(posterior.mode(m$VCV))  #store list of random effects
  (print(mname))
  tmp_df <- data.frame()
  for (r in rlist){
    (print(r))
    tmp <- pmap_dfr(list(mu       = m[["Sol"]][ , "(Intercept)"],
                         var.comp = m[["VCV"]][ , r],
                         var.p    = rowSums(m[["VCV"]])),
                    QGicc,
                    model = "Poisson.log",
                    verbose = FALSE)
    tmp$component <- r
    tmp_df <- rbind(tmp_df,tmp)
  }
  tmp_df$model <- mname
  tmp_df$FE.var <- "No"
  tmp_df$model.type <- "Poisson.log"
  tmp_df
}

# use function on first chain of models
P_ICC1 <- genICC_P1(p09a[[1]])
#P_ICC1.0 <- genICC_P1(p09a.0[[1]])
P_ICC1.1 <- genICC_P1(p09a.1[[1]])
P_ICC1.2 <- genICC_P1(p09a.2[[1]])
P_ICC1.3 <- genICC_P1(p09a.3[[1]])
P_ICC1.4 <- genICC_P1(p09a.4[[1]])

# function for generating posterior distributions for ICC estimates
## poisson models, including fixed effects variance in denominator
genICC_P2 <- function(m){  #m=model, col=variance component
  mname <- deparse(substitute(m))  #store name of model
  print(mname)
  rlist <- colnames(m$VCV)  #store list of random effects
  flist <- colnames(m$Sol)  #store list of fixed effects
  flist <- flist[flist != "log(n)"]  #remove log(n)
  X  <- m[["X"]][,flist]
  predictF <- map(1:nrow(m[["Sol"]]),
  								~ as.vector(X %*% m[["Sol"]][.,flist]))
  tmp_df <- data.frame()
  for (r in rlist){
    (print(r))
    tmp <- pmap_dfr(list(predict  = predictF,
                         var.comp = m[["VCV"]][ , r],
                         var.p    = rowSums(m[["VCV"]])),
                    QGicc,
                    model = "Poisson.log",
                    verbose = FALSE)
    tmp$component <- r
    tmp_df <- rbind(tmp_df,tmp)
  }
  tmp_df$model <- mname
  tmp_df$FE.var <- "Yes"
  tmp_df$model.type <- "Poisson.log"
  tmp_df
}

# use function on first chain of models
P_ICC2 <- genICC_P2(p09a[[1]])
P_ICC2.0 <- genICC_P2(p09a.0[[1]])
P_ICC2.1 <- genICC_P2(p09a.1[[1]])
P_ICC2.2 <- genICC_P2(p09a.2[[1]])
P_ICC2.3 <- genICC_P2(p09a.3[[1]])
P_ICC2.4 <- genICC_P2(p09a.4[[1]])

###############################################################################
###############################################################################
################### values for Table 2 (h2) ############################
###############################################################################
###############################################################################

### heritability estimates

### NOTE: not including the variance from the fixed effects in the denominator
###############################################################################

obs <- 32
QGpost_M1 <- function(model,column){
  pdist <-
    pmap_dfr(list(mu    = model[["Sol"]][ , "(Intercept)"], 
                  var.a = model[["VCV"]][ , column],
                  var.p = rowSums(model[["VCV"]])),
    				 QGparams,
             model = "binomN.logit", 
             n.obs = obs, 
             verbose = FALSE)
  return(pdist)
}

M_h2_params1 <- QGpost_M1(m09a[[1]], "animal")
M_h2_params1.0 <- QGpost_M1(m09a.0[[1]], "animal")
M_h2_params1.1 <- QGpost_M1(m09a.1[[1]], "animal")
M_h2_params1.2 <- QGpost_M1(m09a.2[[1]], "animal")
M_h2_params1.3 <- QGpost_M1(m09a.3[[1]], "animal")
M_h2_params1.4 <- QGpost_M1(m09a.4[[1]], "animal")
# mean.obs	Phenotypic mean on the observed scale.
# var.obs	  Phenotypic variance on the observed scale.
# var.a.obs	Additive genetic variance on the observed scale.
# h2.obs	  Heritability on the observed scale.

posterior.mode(as.mcmc(M_h2_params1))
#   mean.obs    var.obs  var.a.obs     h2.obs 
# 15.9030805 27.1298392  3.1479818  0.1098704 
HPDinterval(as.mcmc(M_h2_params1))
#                lower      upper
# mean.obs  14.6967997 17.0654444
# var.obs   25.4454628 30.8002737
# var.a.obs  1.8704563  4.2162449
# h2.obs     0.0682225  0.1511407
# attr(,"Probability")
# [1] 0.95

QGpost_P1 <- function(model,column){
  pdist <-
    pmap_dfr(list(mu    = model[["Sol"]][ , "(Intercept)"], 
                  var.a = model[["VCV"]][ , column],
                  var.p = rowSums(model[["VCV"]])),
    				 QGparams,
             model = "Poisson.log", 
             verbose = FALSE)
  return(pdist)
}

P_h2_params1 <-   QGpost_P1(p09a[[1]],  "animal")
P_h2_params1.0 <- QGpost_P1(p09a.0[[1]],"animal")
P_h2_params1.1 <- QGpost_P1(p09a.1[[1]],"animal")
P_h2_params1.2 <- QGpost_P1(p09a.2[[1]],"animal")
P_h2_params1.3 <- QGpost_P1(p09a.3[[1]],"animal")
P_h2_params1.4 <- QGpost_P1(p09a.4[[1]],"animal")

posterior.mode(as.mcmc(P_h2_params1))
#   mean.obs    var.obs  var.a.obs     h2.obs 
# 0.92992220 1.25682310 0.02722049 0.02210571 
HPDinterval(as.mcmc(P_h2_params1))
#                lower      upper
# mean.obs  0.82840956 1.07924638
# var.obs   1.04567662 1.47703081
# var.a.obs 0.01753936 0.04068608
# h2.obs    0.01534259 0.03030106
# attr(,"Probability")
# [1] 0.95

### NOTE: including the variance from the fixed effects in the denominator ###
###############################################################################

obs <- round(mean(cebus$n))
QGpost_M2 <- function(model,column){
  X <- model[["X"]]
  predictF <- map(1:nrow(model[["Sol"]]),
                 ~ as.vector(X %*% model[["Sol"]][., ]))
  pdist <-
    pmap_dfr(list(predict = predictF,
                  var.a   = model[["VCV"]][ , column],
                  var.p   = rowSums(model[["VCV"]])),
    				 QGparams,
             model = "binomN.logit", 
             n.obs = obs, 
             verbose = FALSE)
  return(pdist)
}

# this will take several hours to run because there is no closed from solution
## for itegrating over the posterior
M_h2_FEvarYes_m09a <- QGpost_M2(m09a[[1]], "animal")
M_h2_FEvarYes_m09a.0 <- QGpost_M2(m09a.0[[1]], "animal")
M_h2_FEvarYes_m09a.1 <- QGpost_M2(m09a.1[[1]], "animal")
M_h2_FEvarYes_m09a.2 <- QGpost_M2(m09a.2[[1]], "animal")
M_h2_FEvarYes_m09a.3 <- QGpost_M2(m09a.3[[1]], "animal")
M_h2_FEvarYes_m09a.4 <- QGpost_M2(m09a.4[[1]], "animal")

posterior.mode(as.mcmc(M_h2_FEvarYes_m09a))
#    mean.obs     var.obs   var.a.obs      h2.obs 
# 16.93705245 32.76292011  3.02092027  0.07949638 
HPDinterval(as.mcmc(M_h2_FEvarYes_m09a))
#                 lower      upper
# mean.obs  15.67461661 17.7517351
# var.obs   31.06685266 36.1129277
# var.a.obs  1.82368764  4.0960864
# h2.obs     0.05522767  0.1216867
# attr(,"Probability")
# [1] 0.95

m <- p09a[[1]]
flist <- colnames(m$Sol)  #store list of fixed effects
flist <- flist[flist != "log(n)"]  #remove log(n)

QGpost_P2 <- function(model,column){
	X  <- m[["X"]][,flist]
	predictF <- map(1:nrow(m[["Sol"]]),
									~ as.vector(X %*% m[["Sol"]][.,flist]))
  pdist <-
    pmap_dfr(list(predict = predictF,
                  var.a   = model[["VCV"]][ , column],
                  var.p   = rowSums(model[["VCV"]])),
    				 QGparams,
             model = "Poisson.log", 
             verbose = FALSE)
  return(pdist)
}

P_h2_params2 <- QGpost_P2(p09a[[1]], "animal")
P_h2_params2.0 <- QGpost_P2(p09a.0[[1]], "animal")
P_h2_params2.1 <- QGpost_P2(p09a.1[[1]], "animal")
P_h2_params2.2 <- QGpost_P2(p09a.2[[1]], "animal")
P_h2_params2.3 <- QGpost_P2(p09a.3[[1]], "animal")
P_h2_params2.4 <- QGpost_P2(p09a.4[[1]], "animal")

posterior.mode(as.mcmc(P_h2_params2))
#   mean.obs    var.obs  var.a.obs     h2.obs 
# 1.02639523 1.47237493 0.03231422 0.02396625 
HPDinterval(as.mcmc(P_h2_params2))
#                lower      upper
# mean.obs  0.92528413 1.18497095
# var.obs   1.23809433 1.76023643
# var.a.obs 0.02183478 0.04891258
# h2.obs    0.01582036 0.03037728
# attr(,"Probability")
# [1] 0.95

### h2 on latent scale, w/o fixed effects variance in denominator ###
#####################################################################

posterior.mode(m09a[[1]]$VCV[,'animal']/rowSums(m09a[[1]]$VCV))  #0.1523473 
HPDinterval(m09a[[1]]$VCV[,'animal']/rowSums(m09a[[1]]$VCV))
#           lower     upper
# var1 0.09373562 0.2069232
# attr(,"Probability")
# [1] 0.95

posterior.mode(p09a[[1]]$VCV[,'animal']/rowSums(p09a[[1]]$VCV))
# 0.1131741 
HPDinterval(p09a[[1]]$VCV[,'animal']/rowSums(p09a[[1]]$VCV))
#           lower     upper
# var1 0.07608581 0.1493976
# attr(,"Probability")
# [1] 0.95

### h2 on latent scale, with fixed effects variance in denominator ###
#####################################################################

# functions taken from de Villemeruil 2021 tutorial
compute_varpred <- function(beta, design_matrix) {
	var(as.vector(design_matrix %*% beta))
}
## from matrix notation m=XB, where
## m: vector of predicted means (for each datapoint)
## X (design matrix): nrows=datapoints, ncolumns= n predictor variables, plus one
## beta (column): vector of parameters (one for each preditor variable)

# fixed effects variance
m <- m09a[[1]]
vf <- apply(m[["Sol"]], 1, compute_varpred, design_matrix = m[["X"]])
vr <- rowSums(m[["VCV"]])
vpFE <- vr + vf

h2_m09a_latent_FEvarYes <- m[["VCV"]][ , "animal"] / vpFE
posterior.mode(h2_m09a_latent_FEvarYes)  
#0.09573489 
HPDinterval(h2_m09a_latent_FEvarYes)
#           lower     upper
# var1 0.06829465 0.1535191
# attr(,"Probability")
# [1] 0.95

m <- p09a[[1]]
flist <- colnames(m$Sol)  #store list of fixed effects
flist <- flist[flist != "log(n)"]  #remove log(n)
X  <- m[["X"]][,flist]
vf <- apply(m[["Sol"]][,flist], 1, compute_varpred, design_matrix = X)
vr <- rowSums(m[["VCV"]])
vpFE <- vr + vf

h2_p09a_latent_FEvarYes <- m[["VCV"]][ , "animal"] / vpFE
posterior.mode(h2_p09a_latent_FEvarYes)
# 0.09411328  
HPDinterval(h2_p09a_latent_FEvarYes)
#           lower     upper
# var1 0.06462044 0.1249515
# attr(,"Probability")
# [1] 0.95

###############################################################################
###############################################################################
########## Supp Fig: ICC and h2 estimates by n.obs ############################
###############################################################################
###############################################################################

## plot repeatability and h2 point estimates, based on variable n.obs
## useful to see, since estimates are around the number of successes for a given number of observations

m  <- m09a[[1]]
mu <- mean(m[["Sol"]][ , "(Intercept)"])
vi <- mean(m[["VCV"]][ , "id"])
vm <- mean(m[["VCV"]][ , "Mother"])
va <- mean(m[["VCV"]][ , "animal"])
vp <- mean(rowSums(m[["VCV"]]))
QGicc(mu = mu, var.comp=vi+vm+va, var.p=vp, model = "binomN.logit", n.obs = obs)
#   mean.obs  var.obs var.comp.obs   icc.obs
# 1 15.49333 26.71796     4.459834 0.1669227
QGparams(mu=mu, var.a=va, var.p=vp, model = "binomN.logit", n.obs = obs)
#   mean.obs  var.obs var.a.obs    h2.obs
# 1 15.49333 26.71796  3.141351 0.1175745

n <- max(cebus$n)

rlst <- data.frame()
for (i in seq(n)){
  x <- QGicc(mu = mu, var.comp=vi+vm+va, var.p=vp, model = "binomN.logit", n.obs = i, verbose = F)
  rlst <- rbind(rlst,x)
}
rlst$n.obs <- as.numeric(rownames(rlst))
rlst$value <- "icc.obs"
rlst$var.a.obs <- NaN
names(rlst)[names(rlst) == "icc.obs"] <- "estimate"

hlst <- data.frame()
for (i in seq(n)){
  x <- QGparams(mu = mu, var.a = va, var.p=vp, model = "binomN.logit", n.obs = i, verbose = F)
  hlst <- rbind(hlst,x)
}
hlst$n.obs <- as.numeric(rownames(hlst))
hlst$value <- "h2.obs"
hlst$var.comp.obs <- NaN
names(hlst)[names(hlst) == "h2.obs"] <- "estimate"

rlst <- rbind(rlst,hlst)

sub <- subset(rlst,n.obs<351)
sub <- droplevels(sub)

p_obs <- ggplot(sub,  
                aes(x=n.obs, y=estimate, color=value)) + 
  geom_point() + 
  geom_line() + 
  theme_classic(base_size=14) + 
  ylim(c(0,0.25)) + 
  scale_x_continuous(breaks=seq(0, 350, 50),
  									 minor_breaks = seq(0, 350, 10)) +
  theme(panel.grid.major = element_line(color = 'grey',
                                        size = 0.25,
                                        linetype = 1),
  			panel.grid.minor.x = element_line(color = 'grey',
  																			size = 0.25,
  																			linetype = "longdash")) + 
  labs(title="Repeatability and heritability estimates by n.obs parameter") + 
	theme(legend.justification=c(1,0.1), legend.position=c(1,0.1))
(p_obs + scale_color_discrete(
	name = "", 
	labels = c(expression(italic(h)*{}^2), "ICC")))

ggsave("supp/SI_Fig4_h2_ICC_byOBS.pdf", height=7, width=7)
 
###############################################################################
###############################################################################
###############################################################################

### Repeatability estimates, on the data scale
### poisson models, including fixed effects variance in denominator

m <- p09a[[1]]
var.mu  <- m[["Sol"]][, "(Intercept)"]
var.tot <- rowSums(m[["VCV"]])
vRepS   <- m[["VCV"]][,'id'] + m[["VCV"]][,'Mother'] + m[["VCV"]][,'animal'] + 
	m[["VCV"]][,'id:Year'] + m[["VCV"]][,'Mother:Group:Year']
vRepL   <- m[["VCV"]][,'id'] + m[["VCV"]][,'Mother'] + m[["VCV"]][,'animal']

flist <- colnames(m$Sol)  #store list of fixed effects
flist <- flist[flist != "log(n)"]  #remove log(n)
# note this means this is an estimate for a sample of n=1
X  <- m[["X"]][,flist]
predictF <- map(1:nrow(m[["Sol"]]),
								~ as.vector(X %*% m[["Sol"]][.,flist]))

tmp <- pmap_dfr(list(predict  = predictF,
										 var.comp = vRepS,
										 var.p    = rowSums(m[["VCV"]])),
								QGicc,
								model = "Poisson.log",
								verbose = FALSE)
posterior.mode(as.mcmc(tmp$icc.obs))
# 0.04380545 
HPDinterval(as.mcmc(tmp$icc.obs))
# 					lower      upper
# var1 0.03691611 0.05212301
# attr(,"Probability")
# [1] 0.95

tmp <- pmap_dfr(list(predict  = predictF,
										 var.comp = vRepL,
										 var.p    = var.tot),
								QGicc,
								model = "Poisson.log",
								verbose = FALSE)
posterior.mode(as.mcmc(tmp$icc.obs))
# 0.0305924 
HPDinterval(as.mcmc(tmp$icc.obs))
# 					lower     upper
# var1 0.02443837 0.0380926
# attr(,"Probability")
# [1] 0.95

tmp <- pmap_dfr(
	list(
		mu       = var.mu,
		var.comp = vRepS,
		var.p    = var.tot
	),
	QGicc,
	model = "Poisson.log")
posterior.mode(as.mcmc(tmp$icc.obs))
# 0.04413263
HPDinterval(as.mcmc(tmp$icc.obs))
#           lower      upper
# var1 0.03639573 0.05219457
# attr(,"Probability")
# [1] 0.95

tmp <- pmap_dfr(
	list(
		mu       = var.mu,
		var.comp = vRepL,
		var.p    = var.tot
	),
	QGicc,
	model = "Poisson.log")
posterior.mode(as.mcmc(tmp$icc.obs))
# 0.02900323
HPDinterval(as.mcmc(tmp$icc.obs))
#           lower      upper
# var1 0.02370325 0.03760494
# attr(,"Probability")
# [1] 0.95

### Repeatability estimates, on the data scale
### multinomial2, including fixed effects variance in denominator

m <- m09a[[1]]
var.mu  <- m[["Sol"]][, "(Intercept)"]
var.tot <- rowSums(m[["VCV"]])
vRepS   <- m[["VCV"]][,'id'] + m[["VCV"]][,'Mother'] + m[["VCV"]][,'animal'] + 
	m[["VCV"]][,'id:Year'] + m[["VCV"]][,'Mother:Group:Year']
vRepL   <- m[["VCV"]][,'id'] + m[["VCV"]][,'Mother'] + m[["VCV"]][,'animal']

X  <- m[["X"]]
predict <- map(1:nrow(m[["Sol"]]), 
							 ~ as.vector(X %*% m[["Sol"]][., ]))

tmp <- pmap_dfr(list(predict  = predict,
										 var.comp = vRepS,
										 var.p    = var.tot),
								QGicc,
								model = "binomN.logit",
								n.obs = obs)
posterior.mode(as.mcmc(tmp$icc.obs))
HPDinterval(as.mcmc(tmp$icc.obs))

tmp <- pmap_dfr(list(predict  = predict,
										 var.comp = vRepL,
										 var.p    = var.tot),
								QGicc,
								model = "binomN.logit",
								n.obs = obs)
posterior.mode(as.mcmc(tmp$icc.obs))
HPDinterval(as.mcmc(tmp$icc.obs))

##############################################################################

## transformations to data scale for 'multinomial2' reduced models
## using a sample of 200 posterior draws to reduce time
## each one may take a day or longer to run
## if you run these, you may wish to save each as an .rds object before
### proceeding to run the next one in case your computer crashes

rlist <- names(posterior.mode(m09a[[1]]$VCV))  #store list of random effects

##### ##### ##### ##### #####

set.seed(2001)   #
ilst <- sort(sample.int(4000, 200, replace=T))  #random sample index
mname <- deparse(substitute(m09a.0[[1]]))  #store name of model
m <- m09a.0[[1]]
X <- m[["X"]]
predictF <- map(1:nrow(m[["Sol"]]),
								~ as.vector(X %*% m[["Sol"]][., ]))
tic()
M_ICC2x <- mclapply(1:length(rlist), function(i) { 
	pmap_dfr(list(predict  = predictF[ilst],
								var.comp = m[["VCV"]][ , rlist[i]][ilst],
								var.p    = rowSums(m[["VCV"]])[ilst]),
					 QGicc,
					 model = "binomN.logit",
					 n.obs = obs, 
					 verbose = FALSE)
}, mc.cores=5)
toc()  #88626.548 sec elapsed
M_ICC2.0 <- mapply(cbind, 
									 M_ICC2x, 
									 "component"=rlist, 
									 model=mname, 
									 SIMPLIFY=F)

##### ##### ##### ##### #####

set.seed(2002)   #28 Days Later :p
ilst <- sort(sample.int(4000, 200, replace=T))  #random sample index
print(ilst)
mname <- deparse(substitute(m09a.1[[1]]))  #store name of model
m <- m09a.1[[1]]
X <- m[["X"]]
predictF <- map(1:nrow(m[["Sol"]]),
								~ as.vector(X %*% m[["Sol"]][., ]))
tic()
M_ICC2x <- mclapply(1:length(rlist), function(i) {
	pmap_dfr(list(predict  = predictF[ilst],
								var.comp = m[["VCV"]][ , rlist[i]][ilst],
								var.p    = rowSums(m[["VCV"]])[ilst]),
					 QGicc,
					 model = "binomN.logit",
					 n.obs = obs, 
					 verbose = FALSE)
}, mc.cores=5)
toc()  #87106.52 sec elapsed
M_ICC2.1 <- mapply(cbind, 
									 M_ICC2x, 
									 "component"=rlist, 
									 model=mname, 
									 SIMPLIFY=F)

##### ##### ##### ##### #####

set.seed(15)   #
ilst <- sort(sample.int(4000, 200, replace=T))  #random sample index
mname <- deparse(substitute(m09a.2[[1]]))  #store name of model
m <- m09a.2[[1]]
X <- m[["X"]]
predictF <- map(1:nrow(m[["Sol"]]),
								~ as.vector(X %*% m[["Sol"]][., ]))
tic()
M_ICC2x <- mclapply(1:length(rlist), function(i) {
	pmap_dfr(list(predict  = predictF[ilst],
								var.comp = m[["VCV"]][ , rlist[i]][ilst],
								var.p    = rowSums(m[["VCV"]])[ilst]),
					 QGicc,
					 model = "binomN.logit",
					 n.obs = obs, 
					 verbose = FALSE)
}, mc.cores=5)
toc()  #81605.252 sec elapsed
M_ICC2.2 <- mapply(cbind, 
									 M_ICC2x, 
									 "component"=rlist, 
									 model=mname, 
									 SIMPLIFY=F)

##### ##### ##### ##### #####

set.seed(1983)   #
ilst <- sort(sample.int(4000, 200, replace=T))  #random sample index
mname <- deparse(substitute(m09a.3[[1]]))  #store name of model
m <- m09a.3[[1]]
X <- m[["X"]]
predictF <- map(1:nrow(m[["Sol"]]),
								~ as.vector(X %*% m[["Sol"]][., ]))
tic()
M_ICC2x <- mclapply(1:length(rlist), function(i) {
	pmap_dfr(list(predict  = predictF[ilst],
								var.comp = m[["VCV"]][ , rlist[i]][ilst],
								var.p    = rowSums(m[["VCV"]])[ilst]),
					 QGicc,
					 model = "binomN.logit",
					 n.obs = obs, 
					 verbose = FALSE)
}, mc.cores=5)
toc()  #86334.702 sec elapsed
M_ICC2.3 <- mapply(cbind, 
									 M_ICC2x, 
									 "component"=rlist, 
									 model=mname, 
									 SIMPLIFY=F)

##### ##### ##### ##### #####

set.seed(101)  
ilst <- sort(sample.int(4000, 200, replace=T))  #random sample index
mname <- deparse(substitute(m09a.4[[1]]))  #store name of model
m <- m09a.4[[1]]
X <- m[["X"]]
predictF <- map(1:nrow(m[["Sol"]]),
								~ as.vector(X %*% m[["Sol"]][., ]))
tic()
M_ICC2x <- mclapply(1:length(rlist), function(i) {
	pmap_dfr(list(predict  = predictF[ilst],
								var.comp = m[["VCV"]][ , rlist[i]][ilst],
								var.p    = rowSums(m[["VCV"]])[ilst]),
					 QGicc,
					 model = "binomN.logit",
					 n.obs = obs, 
					 verbose = FALSE)
}, mc.cores=5)
toc()  #84143.67 sec elapsed
M_ICC2.4 <- mapply(cbind, 
									 M_ICC2x, 
									 "component"=rlist, 
									 model=mname, 
									 SIMPLIFY=F)