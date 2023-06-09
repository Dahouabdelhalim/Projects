### RÃ¸d-Eriksen et al.: Predator interactions in alpine and Arctic tundra in relation to fluctuating prey
###
### Occupancy.R - Multi-species occupancy modelling following Rota et al. (2016)

library(unmarked)
library(AICcmodavg)


##### IMPORT DATA #####
#alpine data
df.alp <- read.csv2("data_alpine.csv", header = T, sep = ",", dec = ".", fileEncoding = "UTF-8")

#arctic data
df.arc <- read.csv2("data_arctic.csv", header = T, sep = ",", dec = ".", fileEncoding = "UTF-8")


##### DEFINE OCCUPANCY FRAMEWORK STRUCTURE #####
#define species matrices
spec.alp.rf <- as.matrix(df.alp[,c(114:167)]); colnames(spec.alp.rf) <- NULL #alpine - red fox
spec.alp.af <- as.matrix(df.alp[,c(168:221)]); colnames(spec.alp.af) <- NULL #alpine - arctic fox
spec.alp.wv <- as.matrix(df.alp[,c(222:275)]); colnames(spec.alp.wv) <- NULL #alpine - wolverine
spec.alp.ge <- as.matrix(df.alp[,c(276:329)]); colnames(spec.alp.ge) <- NULL #alpine - golden eagle
spec.arc.rf <- as.matrix(df.arc[,c(80:116)]); colnames(spec.arc.rf) <- NULL #arctic - red fox
spec.arc.af <- as.matrix(df.arc[,c(117:153)]); colnames(spec.arc.af) <- NULL #arctic - arctic fox
spec.arc.wv <- as.matrix(df.arc[,c(154:190)]); colnames(spec.arc.wv) <- NULL #arctic - wolverine
spec.arc.ge <- as.matrix(df.arc[,c(191:227)]); colnames(spec.arc.ge) <- NULL #arctic - golden eagle

#add matrices to list structure
#alpine
y.alp <- list()
y.alp[[1]] <- spec.alp.rf
y.alp[[2]] <- spec.alp.af
y.alp[[3]] <- spec.alp.wv
y.alp[[4]] <- spec.alp.ge
names(y.alp) <- c('redfox','arcfox','wolver','eagles')

#arctic
y.arc <- list()
y.arc[[1]] <- spec.arc.rf
y.arc[[2]] <- spec.arc.af
y.arc[[3]] <- spec.arc.wv
y.arc[[4]] <- spec.arc.ge
names(y.arc) <- c('redfox','arcfox','wolver','eagles')


#define detection covariate matrices
obs.alp.doy <- as.matrix(df.alp[,c(6:59)]); colnames(obs.alp.doy) <- NULL #day of year
obs.alp.carc <- as.matrix(df.alp[,c(60:113)]); colnames(obs.alp.carc) <- NULL #carcass status (absent/present)
obs.arc.doy <- as.matrix(df.arc[,c(6:42)]); colnames(obs.arc.doy) <- NULL #day of year
obs.arc.carc <- as.matrix(df.arc[,c(43:79)]); colnames(obs.arc.carc) <- NULL #carcass status (absent/present)

#add matrices to list structure
#alpine
det.covs.alp <- list()
det.covs.alp[[1]] <- obs.alp.doy #day of year
det.covs.alp[[2]] <- obs.alp.carc #carcass status (absent/present)
names(det.covs.alp) <- c('doy','carcass')

#arctic
det.covs.arc <- list()
det.covs.arc[[1]] <- obs.arc.doy #day of year
det.covs.arc[[2]] <- obs.arc.carc #carcass status (absent/present)
names(det.covs.arc) <- c('doy','carcass')


#define site covariate data frames
#center-scale numerical site covariates to maximize convergence
site.covs.alp <- data.frame(rods = scale(df.alp$rod_idx, center = T, scale = T), #alpine - rodent index
                            lalo = scale(df.alp$lalo, center = T, scale = T)) #alpine - long*lat interaction
site.covs.arc <- data.frame(rods = scale(df.arc$rod_idx, center = T, scale = T), #arctic - rodent index
                            lalo = scale(df.arc$lalo, center = T, scale = T)) #arctic - long*lat interaction


##### MULTI-SPECIES OCCUPANCY MODELLING #####
#create unmarked occuMulti data frames
udf.alp = unmarkedFrameOccuMulti(y = y.alp, siteCovs = site.covs.alp, obsCovs = det.covs.alp) #ALPINE region
udf.arc = unmarkedFrameOccuMulti(y = y.arc, siteCovs = site.covs.arc, obsCovs = det.covs.arc) #ARCTIC region

#inspect data frames
summary(udf.alp)
summary(udf.arc)

#define the ad hoc models
#0. detection
dd.doy.carc <- c('~doy+carcass','~doy+carcass','~doy+carcass','~doy+carcass')
#1. marginal without rodents
oc.lalo <- c('~lalo','~lalo','~lalo','~lalo','0','0','0','0','0','0','0','0','0','0','0')
#2. marginal with rodents
oc.lalo.rods <- c('~lalo+rods','~lalo+rods','~lalo+rods','~lalo+rods','0','0','0','0','0','0','0','0','0','0','0')
#3. pairwise conditional without rodents (intercept for pairwise interactions)
oc.lalo.cond.2.norods <- c('~lalo','~lalo','~lalo','~lalo','~1','~1','~1','~1','~1','~1','0','0','0','0','0')
#4. pairwise conditional with rodents
oc.lalo.cond.2.rods <- c('~lalo','~lalo','~lalo','~lalo','~rods','~rods','~rods','~rods','~rods','~rods','0','0','0','0','0')
#5. full conditional without rodents (intercept for all interactions)
oc.lalo.cond.4.norods <- c('~lalo','~lalo','~lalo','~lalo','~1','~1','~1','~1','~1','~1','~1','~1','~1','~1','~1')
#6. full conditional with rodents (intercept for pairwise interactions)
oc.lalo.cond.4.rods <- c('~lalo','~lalo','~lalo','~lalo','~1','~1','~1','~1','~1','~1','~rods','~rods','~rods','~rods','~rods')

#FIT MODELS - ALPINE
#1: marginal without rodents
alp.marg.norod <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo, data = udf.alp, method = "BFGS", control = list(maxit = 2000))

#2: marginal with rodents
alp.marg <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo.rods, data = udf.alp, method = "BFGS", control = list(maxit = 2000))

#3: pairwise conditional without rodents (intercept pairwise only)
alp.pair.norod <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo.cond.2.norods, data = udf.alp, method = "BFGS", control = list(maxit = 2000))

#4: pairwise conditional with rodents
alp.pair <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo.cond.2.rods, data = udf.alp, method = "BFGS", control = list(maxit = 2000))

#5: full conditional without rodents (intercept only)
alp.full.norod <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo.cond.4.norods, data = udf.alp, method = "BFGS", control = list(maxit = 2000))

#6: full conditional with rodents (intercept pairwise)
alp.full <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo.cond.4.rods, data = udf.alp, method = "BFGS", control = list(maxit = 2000))

#model list
mods.alp = list(marg.norod = alp.marg.norod, marg = alp.marg, pair.norod = alp.pair.norod,
                pair = alp.pair, full.norod = alp.full.norod, full = alp.full)
aictab(cand.set = mods.alp) #inspect the models

## MODEL VALIDATION ##
#cross-validation of model list
#results from cross-validation indicates whether the data is underfit or overfit to the model, i.e. how well the model performs
#calculates the Root-Mean-Square-Error (RMSE) and Mean Absolute Error (MAE)
mList.alp <- fitList(fits = mods.alp)
kfold.alp <- crossVal(mList.alp, method = "Kfold", folds = 5, parallel = TRUE, sort = "increasing", holdoutPct = .25)
kfold.alp #inspect results


#FIT MODELS - ARCTIC
#1: marginal without rodents
arc.marg.norod <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo, data = udf.arc, method = "BFGS", control = list(maxit = 2000))

#1: marginal with rodents
arc.marg <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo.rods, data = udf.arc, method = "BFGS", control = list(maxit = 2000))

#2: pairwise conditional without rodents
arc.pair.norod <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo.cond.2.norods, data = udf.arc, method = "BFGS", control = list(maxit = 2000))

#2: pairwise conditional with rodents
arc.pair <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo.cond.2.rods, data = udf.arc, method = "BFGS", control = list(maxit = 2000))

#3: full conditional without rodents
arc.full.norod <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo.cond.4.norods, data = udf.arc, method = "BFGS", control = list(maxit = 2000))

#3: full conditional with rodents
arc.full <- occuMulti(detformulas = dd.doy.carc, stateformulas = oc.lalo.cond.4.rods, data = udf.arc, method = "BFGS", control = list(maxit = 2000))

#model list
mods.arc = list(marg.norod = arc.marg.norod, marg = arc.marg, pair.norod = arc.pair.norod,
                pair = arc.pair, full.norod = arc.full.norod, full = arc.full)
aictab(cand.set = mods.arc) #inspect the models

## MODEL VALIDATION ##
#cross-validation of model list
#results from cross-validation indicates whether the data is underfit or overfit to the model, i.e. how well the model performs
#calculates the Root-Mean-Square-Error (RMSE) and Mean Absolute Error (MAE)
mList.arc <- fitList(fits = mods.arc)
kfold.arc <- crossVal(mList.arc, method = "Kfold", folds = 5, parallel = TRUE, sort = "increasing", holdoutPct = .25)
kfold.arc #inspect results


##### DESCRIPTIVE STATS (MARGINAL) #####
#ALPINE
#average detection probabilities
alp.det.pred.af <- predict(alp.marg.norod, type = "det", species = 'arcfox', appendData = TRUE, nsim = 5000)
alp.det.pred.af <- data.frame(Predicted = mean(alp.det.pred.af$Predicted, na.rm = T), SE = mean(alp.det.pred.af$SE, na.rm = T), lower = mean(alp.det.pred.af$lower, na.rm = T), upper = mean(alp.det.pred.af$upper, na.rm = T))
alp.det.pred.af$species <- "Arctic fox"
alp.det.pred.rf <- predict(alp.marg.norod, type = "det", species = 'redfox', appendData = TRUE, nsim = 5000)
alp.det.pred.rf <- data.frame(Predicted = mean(alp.det.pred.rf$Predicted, na.rm = T), SE = mean(alp.det.pred.rf$SE, na.rm = T), lower = mean(alp.det.pred.rf$lower, na.rm = T), upper = mean(alp.det.pred.rf$upper, na.rm = T))
alp.det.pred.rf$species <- "Red fox"
alp.det.pred.wo <- predict(alp.marg.norod, type = "det", species = 'wolver', appendData = TRUE, nsim = 5000)
alp.det.pred.wo <- data.frame(Predicted = mean(alp.det.pred.wo$Predicted, na.rm = T), SE = mean(alp.det.pred.wo$SE, na.rm = T), lower = mean(alp.det.pred.wo$lower, na.rm = T), upper = mean(alp.det.pred.wo$upper, na.rm = T))
alp.det.pred.wo$species <- "Wolverine"
alp.det.pred.ea <- predict(alp.marg.norod, type = "det", species = 'eagles', appendData = TRUE, nsim = 5000)
alp.det.pred.ea <- data.frame(Predicted = mean(alp.det.pred.ea$Predicted, na.rm = T), SE = mean(alp.det.pred.ea$SE, na.rm = T), lower = mean(alp.det.pred.ea$lower, na.rm = T), upper = mean(alp.det.pred.ea$upper, na.rm = T))
alp.det.pred.ea$species <- "Golden eagle"
(alp.det.pred.df <- rbind(alp.det.pred.af, alp.det.pred.rf, alp.det.pred.wo, alp.det.pred.ea))

#average marginal occupancy probabilities
alp.pred.af <- predict(alp.marg.norod, type = "state", species = 'arcfox', appendData = TRUE, nsim = 5000)
alp.pred.af.df <- data.frame(Predicted = mean(alp.pred.af$Predicted), SE = mean(alp.pred.af$SE), lower = mean(alp.pred.af$lower), upper = mean(alp.pred.af$upper))
alp.pred.af.df$species <- "Arctic fox"
alp.pred.rf <- predict(alp.marg.norod, type = "state", species = 'redfox', appendData = TRUE, nsim = 5000)
alp.pred.rf.df <- data.frame(Predicted = mean(alp.pred.rf$Predicted), SE = mean(alp.pred.rf$SE), lower = mean(alp.pred.rf$lower), upper = mean(alp.pred.rf$upper))
alp.pred.rf.df$species <- "Red fox"
alp.pred.wo <- predict(alp.marg.norod, type = "state", species = 'wolver', appendData = TRUE, nsim = 5000)
alp.pred.wo.df <- data.frame(Predicted = mean(alp.pred.wo$Predicted), SE = mean(alp.pred.wo$SE), lower = mean(alp.pred.wo$lower), upper = mean(alp.pred.wo$upper))
alp.pred.wo.df$species <- "Wolverine"
alp.pred.ea <- predict(alp.marg.norod, type = "state", species = 'eagles', appendData = TRUE, nsim = 5000)
alp.pred.ea.df <- data.frame(Predicted = mean(alp.pred.ea$Predicted), SE = mean(alp.pred.ea$SE), lower = mean(alp.pred.ea$lower), upper = mean(alp.pred.ea$upper))
alp.pred.ea.df$species <- "Golden eagle"
(alp.pred.df <- rbind(alp.pred.af.df, alp.pred.rf.df, alp.pred.wo.df, alp.pred.ea.df))

#ARCTIC
#average detection probabilities
arc.det.pred.af <- predict(arc.marg.norod, type = "det", species = 'arcfox', appendData = TRUE, nsim = 5000)
arc.det.pred.af <- data.frame(Predicted = mean(arc.det.pred.af$Predicted, na.rm = T), SE = mean(arc.det.pred.af$SE, na.rm = T), lower = mean(arc.det.pred.af$lower, na.rm = T), upper = mean(arc.det.pred.af$upper, na.rm = T))
arc.det.pred.af$species <- "Arctic fox"
arc.det.pred.rf <- predict(arc.marg.norod, type = "det", species = 'redfox', appendData = TRUE, nsim = 5000)
arc.det.pred.rf <- data.frame(Predicted = mean(arc.det.pred.rf$Predicted, na.rm = T), SE = mean(arc.det.pred.rf$SE, na.rm = T), lower = mean(arc.det.pred.rf$lower, na.rm = T), upper = mean(arc.det.pred.rf$upper, na.rm = T))
arc.det.pred.rf$species <- "Red fox"
arc.det.pred.wo <- predict(arc.marg.norod, type = "det", species = 'wolver', appendData = TRUE, nsim = 5000)
arc.det.pred.wo <- data.frame(Predicted = mean(arc.det.pred.wo$Predicted, na.rm = T), SE = mean(arc.det.pred.wo$SE, na.rm = T), lower = mean(arc.det.pred.wo$lower, na.rm = T), upper = mean(arc.det.pred.wo$upper, na.rm = T))
arc.det.pred.wo$species <- "Wolverine"
arc.det.pred.ea <- predict(arc.marg.norod, type = "det", species = 'eagles', appendData = TRUE, nsim = 5000)
arc.det.pred.ea <- data.frame(Predicted = mean(arc.det.pred.ea$Predicted, na.rm = T), SE = mean(arc.det.pred.ea$SE, na.rm = T), lower = mean(arc.det.pred.ea$lower, na.rm = T), upper = mean(arc.det.pred.ea$upper, na.rm = T))
arc.det.pred.ea$species <- "Golden eagle"
(arc.det.pred.df <- rbind(arc.det.pred.af, arc.det.pred.rf, arc.det.pred.wo, arc.det.pred.ea))

#average marginal occupancy probabilities
arc.pred.af <- predict(arc.marg.norod, type = "state", species = 'arcfox', appendData = TRUE, nsim = 5000)
arc.pred.af.df <- data.frame(Predicted = mean(arc.pred.af$Predicted), SE = mean(arc.pred.af$SE), lower = mean(arc.pred.af$lower), upper = mean(arc.pred.af$upper))
arc.pred.af.df$species <- "Arctic fox"
arc.pred.rf <- predict(arc.marg.norod, type = "state", species = 'redfox', appendData = TRUE, nsim = 5000)
arc.pred.rf.df <- data.frame(Predicted = mean(arc.pred.rf$Predicted), SE = mean(arc.pred.rf$SE), lower = mean(arc.pred.rf$lower), upper = mean(arc.pred.rf$upper))
arc.pred.rf.df$species <- "Red fox"
arc.pred.wo <- predict(arc.marg.norod, type = "state", species = 'wolver', appendData = TRUE, nsim = 5000)
arc.pred.wo.df <- data.frame(Predicted = mean(arc.pred.wo$Predicted), SE = mean(arc.pred.wo$SE), lower = mean(arc.pred.wo$lower), upper = mean(arc.pred.wo$upper))
arc.pred.wo.df$species <- "Wolverine"
arc.pred.ea <- predict(arc.marg.norod, type = "state", species = 'eagles', appendData = TRUE, nsim = 5000)
arc.pred.ea.df <- data.frame(Predicted = mean(arc.pred.ea$Predicted), SE = mean(arc.pred.ea$SE), lower = mean(arc.pred.ea$lower), upper = mean(arc.pred.ea$upper))
arc.pred.ea.df$species <- "Golden eagle"
(arc.pred.df <- rbind(arc.pred.af.df, arc.pred.rf.df, arc.pred.wo.df, arc.pred.ea.df))


# NB! We have excluded the code concerning plots and further predictions from this script.
# Please contact the corresponding author if these parts of the script are of interest.
