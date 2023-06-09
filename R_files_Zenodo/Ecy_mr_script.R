# load relevant packages
library(RMark)
library(tidyverse)

# Read mark-recapture data from  file
mark_data <- read_csv(file.path('dryad_mark_recapture_dataset.csv'),
                      col_types = 'fffc') %>%
  rename(ch = capture_history)

# Read file containing time intervals between samples
sample_lags <- read_csv(file.path('dryad_mark_recapture_lags.csv'),
                        col_types = 'ffDd')
############################################################
########### CAI Mark Recapture Models ######################
############################################################

# Process data for cai
caich <- process.data(mark_data[mark_data$site == 'CAI',], 
                      model = 'POPAN', 
                      time.intervals = sample_lags$lag30[sample_lags$site == 'CAI' & sample_lags$reach == 'C' & !is.na(sample_lags$lag30)],
                      groups = 'reach')

# Set-up design matrix
design <- make.design.data(caich)

# Detection probability formulae
p.int <- list(formula=~1,link='sin')
p.time <- list(formula=~-1+time, link='sin')
p.g <- list(formula=~-1+group, link='sin')
p.tpg <- list(formula=~-1+time:group, link='sin')

# Survival probability formula
Phi.sin <- list(formula=~-1+time:group,link='sin')

# Probability of entry formula
pent.timegroup <- list(formula=~time*group)

# Fit all models with 4 detection probabilities
caimod <- create.model.list('POPAN')
caiout <- mark.wrapper(caimod,data=caich, ddl=design)
caiout

# Fit top model from file (remove filename argument to fit from data)
topmod <- mark(caich, ddl=design, model="POPAN", 
               model.parameters=list('Phi'=Phi.sin,'p'=p.time,'pent'=pent.timegroup, 'N'=list(formula=~1)),
               filename='dryad_mark/CAI/mark005')

# Compute derived parameters, abundance and density
Nhat.cai <- popan.derived(caich, topmod)
cai.pop <- Nhat.cai$N
cai.pop$D <- cai.pop$D.se <- cai.pop$D.lcl <- cai.pop$D.ucl <- NA
cai.pop$length <- c(rep(35,18), rep(73, 18))
# Compute density with reach length
cai.pop[,c('D','D.lcl','D.ucl')] <- cai.pop[,c('N',"LCL","UCL")]/cai.pop$length
# Rescale standard error for density
cai.pop$D.se <- sqrt(cai.pop$se^2*cai.pop$length^(-2))

# Add some variables
cai.pop$date <- NA
cai.pop$date[cai.pop$reach == 'C'] <- sample_lags$date[sample_lags$site == "CAI" & sample_lags$reach == 'C']
cai.pop$date[cai.pop$reach=='I'] <- sample_lags$date[sample_lags$site=="CAI" & sample_lags$reach=='I']
cai.pop$date <- as.Date(cai.pop$date, origin='1970-01-01')
cai.pop$site  <- 'CAI'
cai.pop$Dlog <- log(cai.pop$D)

# Define differencing operator between control and introduction, ignore first and last samples
oper <- rbind(diag(c(0,rep(1,16),0),nrow=18,ncol=18),diag(c(0,rep(-1,16),0),nrow=18,ncol=18))
# Take the difference
dD.c <- t(oper) %*% cai.pop$Dlog
# Scale vcv matrix of abundance with length of each reach
Dvcv <- Nhat.cai$N.vcv * rbind(cbind(matrix(35^-2,nrow=18,ncol=18), matrix((35*73)^-1,nrow=18,ncol=18)),cbind(matrix((35*73)^-1,nrow=18,ncol=18), matrix(73^-2,nrow=18,ncol=18)))

# create a list of function to apply the delta method
dml <- list(~log(x1), ~log(x2), ~log(x3), ~log(x4),~log(x5), ~log(x6), ~log(x7), ~log(x8), ~log(x9), ~log(x10),~log(x11), ~log(x12), ~log(x13), ~log(x14), ~log(x15), ~log(x16), ~log(x17), ~log(x18), ~log(x19), ~log(x20), ~log(x21), ~log(x22), ~log(x23), ~log(x24), ~log(x25), ~log(x26), ~log(x27), ~log(x28), ~log(x29), ~log(x30), ~log(x31), ~log(x32), ~log(x33), ~log(x34), ~log(x35), ~log(x36))
# Apply delta method to find log scale vcv of density
Dlog.vcv <- msm::deltamethod(dml, mean = cai.pop$D, Dvcv, ses=F)

# Calculate density difference vcv
dD.c.vcv <- t(oper) %*% Dlog.vcv %*% oper
dD.se <- sqrt(diag(dD.c.vcv))
dD.c.vcv <- dD.c.vcv[2:17,2:17]
dD.cc <- data.frame(site='CAI', dD=dD.c, 'se'=dD.se)
dD.cc$date <- cai.pop$date[1:18]

############################################################
########### TAY Mark Recapture Models ######################
############################################################

# Repeat the same analysis for Taylor
# Process data
taych <- process.data(mark_data[mark_data$site == 'TAY',], 
                      model = 'POPAN', 
                      time.intervals = sample_lags$lag30[sample_lags$site == 'TAY' & sample_lags$reach == 'C' & !is.na(sample_lags$lag30)],
                      groups = 'reach')
# Set-up design matrix
design <- make.design.data(taych)

# Detection probability formulae
p.int <- list(formula=~1,link='sin')
p.time <- list(formula=~-1+time, link='sin')
p.g <- list(formula=~-1+group, link='sin')
p.tpg <- list(formula=~-1+time:group, link='sin')

# Survival probability formula
Phi.sin <- list(formula=~-1+time:group,link='sin')

# Probability of entry formula
pent.timegroup <- list(formula=~time*group)

# Fit all models with 4 detection probabilities
taymod <- create.model.list('POPAN')
tayout <- mark.wrapper(taymod, data=taych, ddl=design)
tayout

# Fit top model
topmod <- mark(taych, ddl=design, model="POPAN", 
               model.parameters=list('Phi'=Phi.sin,'p'=p.tpg,'pent'=pent.timegroup, 'N'=list(formula=~1)),
               filename = 'dryad_mark/TAY/mark037')

# Note sometimes mark fails to converge, producing singularities, inestimable SEs or the wrong number of estimated parameters
# This loop fits a model to data nmod times and extracts information about the common model fitting problems.
# This procedure can be used find problematic model fits and select a fit that has converged.

# nmod <- 50
# problems <- data.frame('mod' = 1:nmod, 'sing' = rep(NA,nmod), 'se0' = rep(NA,nmod), 'npar' = rep(NA,nmod))
# for (i in 1:nmod){
#   mod <- mark(taych, ddl=design, model='POPAN', 
#               model.parameters=list('Phi'=Phi.sin, 'p'=p.tpg, 'pent'=pent.timegroup, 'N'=list(formula=~1)), 
#               filename=paste0(ifelse(i<10,'mark00','mark0'),i))
#   problems$sing[i] <- length(mod$results$singular[!is.na(mod$results$singular)])
#   problems$se0[i] <- sum(mod$results$beta$se==0)
#   problems$npar[i] <- mod$results$npar.unadjusted
# }

# Compute derived parameters, abundance and density
Nhat.tay <- popan.derived(taych, topmod)
tay.pop <- Nhat.tay$N
tay.pop$D <- tay.pop$D.se <- tay.pop$D.lcl <- tay.pop$D.ucl <- NA
tay.pop$length <- c(rep(26,18),rep(50, 18))
# Compute density with reach length
tay.pop[,c('D','D.lcl','D.ucl')] <- tay.pop[,c('N',"LCL","UCL")]/tay.pop$length
# Rescale standard error for density
tay.pop$D.se <- sqrt(tay.pop$se^2*tay.pop$length^(-2))

# Add some variables
tay.pop$date <- NA
tay.pop$date[tay.pop$reach == 'C'] <- sample_lags$date[sample_lags$site == "TAY" & sample_lags$reach == 'C']
tay.pop$date[tay.pop$reach=='I'] <- sample_lags$date[sample_lags$site=="TAY" & sample_lags$reach=='I']
tay.pop$date <- as.Date(tay.pop$date, origin='1970-01-01')
tay.pop$site  <- 'TAY'
tay.pop$Dlog <- log(tay.pop$D)

# Define differencing operator between control and introduction, ignore first and last samples
oper <- rbind(diag(c(0,rep(1,16),0),nrow=18,ncol=18),diag(c(0,rep(-1,16),0),nrow=18,ncol=18))
# Take the difference
dD.t <- t(oper) %*% tay.pop$Dlog
# Scale vcv matrix of abundance with length of each reach
Dvcv <- Nhat.tay$N.vcv * rbind(cbind(matrix(26^-2,nrow=18,ncol=18), matrix((26*50)^-1,nrow=18,ncol=18)),cbind(matrix((26*50)^-1,nrow=18,ncol=18), matrix(50^-2,nrow=18,ncol=18)))

# create a list of function to apply the delta method
dml <- list(~log(x1), ~log(x2), ~log(x3), ~log(x4),~log(x5), ~log(x6), ~log(x7), ~log(x8), ~log(x9), ~log(x10),~log(x11), ~log(x12), ~log(x13), ~log(x14), ~log(x15), ~log(x16), ~log(x17), ~log(x18), ~log(x19), ~log(x20), ~log(x21), ~log(x22), ~log(x23), ~log(x24), ~log(x25), ~log(x26), ~log(x27), ~log(x28), ~log(x29), ~log(x30), ~log(x31), ~log(x32), ~log(x33), ~log(x34), ~log(x35), ~log(x36))
# Apply delta method to find log scale vcv of density
Dlog.vcv <- msm::deltamethod(dml, mean = tay.pop$D, Dvcv, ses=F)

# Calculate density difference vcv
dD.t.vcv <- t(oper) %*% Dlog.vcv %*% oper
dD.se <- sqrt(diag(dD.t.vcv))
dD.t.vcv <- dD.t.vcv[2:17,2:17]
dD.tt <- data.frame(site='TAY', dD=dD.t, 'se'=dD.se)
dD.tt$date <- tay.pop$date[1:18]

############################################################
########### LOL Mark Recapture Models ######################
############################################################

# Repeat the same analysis for Taylor
# Process data
lolch <- process.data(mark_data[mark_data$site == 'LOL',], 
                      model = 'POPAN', 
                      time.intervals = sample_lags$lag30[sample_lags$site == 'LOL' & sample_lags$reach == 'C' & !is.na(sample_lags$lag30)],
                      groups = 'reach')
# Set-up design matrix
design <- make.design.data(lolch)

# Detection probability formulae
p.int <- list(formula=~1,link='sin')
p.time <- list(formula=~-1+time, link='sin')
p.g <- list(formula=~-1+group, link='sin')
p.tpg <- list(formula=~-1+time:group, link='sin')

# Survival probability formula
Phi.sin <- list(formula=~-1+time:group,link='sin')

# Probability of entry formula
pent.timegroup <- list(formula=~time*group)

# Fit all models with 4 detection probabilities
lolmod <- create.model.list('POPAN')
lolout <- mark.wrapper(lolmod, data=lolch, ddl=design)
lolout

# Fit top model
topmod <- mark(lolch, ddl=design, model="POPAN", 
               model.parameters=list('Phi'=Phi.sin,'p'=p.tpg,'pent'=pent.timegroup, 'N'=list(formula=~1)),
               filename = 'dryad_mark/LOL/mark001')

# Note sometimes mark fails to converge, producing singularities, inestimable SEs or the wrong number of estimated parameters
# This loop fits a model to data nmod times and extracts information about the common model fitting problems.
# This procedure can be used find problematic model fits and select a fit that has converged.

# nmod <- 50
# problems <- data.frame('mod' = 1:nmod, 'sing' = rep(NA,nmod), 'se0' = rep(NA,nmod), 'npar' = rep(NA,nmod))
# for (i in 1:nmod){
#   mod <- mark(lolch, ddl=design, model='POPAN', 
#               model.parameters=list('Phi'=Phi.sin, 'p'=p.tpg, 'pent'=pent.timegroup, 'N'=list(formula=~1)), 
#               filename=paste0(ifelse(i<10,'mark00','mark0'),i))
#   problems$sing[i] <- length(mod$results$singular[!is.na(mod$results$singular)])
#   problems$se0[i] <- sum(mod$results$beta$se==0)
#   problems$npar[i] <- mod$results$npar.unadjusted
# }

# Compute derived parameters, abundance and density
Nhat.lol <- popan.derived(lolch, topmod)
lol.pop <- Nhat.lol$N
lol.pop$D <- lol.pop$D.se <- lol.pop$D.lcl <- lol.pop$D.ucl <- NA
lol.pop$length <- c(rep(62, 27),rep(72,27))
# Compute density with reach length
lol.pop[,c('D','D.lcl','D.ucl')] <- lol.pop[,c('N',"LCL","UCL")]/lol.pop$length
# Rescale standard error for density
lol.pop$D.se <- sqrt(lol.pop$se^2*lol.pop$length^(-2))

# Add some variables
lol.pop$date <- NA
lol.pop$date[lol.pop$reach == 'C'] <- sample_lags$date[sample_lags$site == "LOL" & sample_lags$reach == 'C']
lol.pop$date[lol.pop$reach=='I'] <- sample_lags$date[sample_lags$site=="LOL" & sample_lags$reach=='I']
lol.pop$date <- as.Date(lol.pop$date, origin='1970-01-01')
lol.pop$site  <- 'LOL'
lol.pop$Dlog <- log(lol.pop$D)

# Define differencing operator between control and introduction, ignore first and last samples
oper <- rbind(diag(c(0,rep(1,25),0),nrow=27,ncol=27),diag(c(0,rep(-1,25),0),nrow=27,ncol=27))
# Take the difference
dD.l <- t(oper) %*% lol.pop$Dlog
# Scale vcv matrix of abundance with length of each reach
Dvcv <- Nhat.lol$N.vcv * rbind(cbind(matrix(62^-2,nrow=27,ncol=27), matrix((62*72)^-1,nrow=27,ncol=27)),cbind(matrix((62*72)^-1,nrow=27,ncol=27), matrix(72^-2,nrow=27,ncol=27)))

# create a list of function to apply the delta method
dml <- list(~log(x1), ~log(x2), ~log(x3), ~log(x4),~log(x5), ~log(x6), ~log(x7), ~log(x8), ~log(x9), ~log(x10),~log(x11), ~log(x12), ~log(x13), ~log(x14), ~log(x15), ~log(x16), ~log(x17), ~log(x18), ~log(x19), ~log(x20), ~log(x21), ~log(x22), ~log(x23), ~log(x24), ~log(x25), ~log(x26), ~log(x27), ~log(x28), ~log(x29), ~log(x30), ~log(x31), ~log(x32), ~log(x33), ~log(x34), ~log(x35), ~log(x36), ~log(x37), ~log(x38), ~log(x39), ~log(x40), ~log(x41), ~log(x42), ~log(x43), ~log(x44), ~log(x45), ~log(x46), ~log(x47), ~log(x48), ~log(x49), ~log(x50), ~log(x51), ~log(x52), ~log(x53), ~log(x54))
# Apply delta method to find log scale vcv of density
Dlog.vcv <- msm::deltamethod(dml, mean = lol.pop$D, Dvcv, ses=F)

# Calculate density difference vcv
dD.l.vcv <- t(oper) %*% Dlog.vcv %*% oper
dD.se <- sqrt(diag(dD.l.vcv))
dD.l.vcv <- dD.l.vcv[2:26,2:26]
dD.ll <- data.frame(site='LOL', dD=dD.l, 'se'=dD.se)
dD.ll$date <- lol.pop$date[1:27]

############################################################
########### UPL Mark Recapture Models ######################
############################################################

# Repeat the same analysis for Taylor
# Process data
uplch <- process.data(mark_data[mark_data$site == 'UPL',], 
                      model = 'POPAN', 
                      time.intervals = sample_lags$lag30[sample_lags$site == 'UPL' & sample_lags$reach == 'C' & !is.na(sample_lags$lag30)],
                      groups = 'reach')
# Set-up design matrix
design <- make.design.data(uplch)

# Detection probability formulae
p.int <- list(formula=~1,link='sin')
p.time <- list(formula=~-1+time, link='sin')
p.g <- list(formula=~-1+group, link='sin')
p.tpg <- list(formula=~-1+time:group, link='sin')

# Survival probability formula
Phi.sin <- list(formula=~-1+time:group,link='sin')

# Probability of entry formula
pent.timegroup <- list(formula=~time*group)

# Fit all models with 4 detection probabilities
uplmod <- create.model.list('POPAN')
uplout <- mark.wrapper(uplmod, data=uplch, ddl=design)
uplout

# Fit top model
topmod <- mark(uplch, ddl=design, model="POPAN", 
               model.parameters=list('Phi'=Phi.sin,'p'=p.tpg,'pent'=pent.timegroup, 'N'=list(formula=~1)),
               filename = 'dryad_mark/UPL/mark045')

# Note sometimes mark fails to converge, producing singularities, inestimable SEs or the wrong number of estimated parameters
# This loop fits a model to data nmod times and extracts information about the common model fitting problems.
# This procedure can be used find problematic model fits and select a fit that has converged.

# nmod <- 50
# problems <- data.frame('mod' = 1:nmod, 'sing' = rep(NA,nmod), 'se0' = rep(NA,nmod), 'npar' = rep(NA,nmod))
# for (i in 1:nmod){
#   mod <- mark(uplch, ddl=design, model='POPAN', 
#               model.parameters=list('Phi'=Phi.sin, 'p'=p.tpg, 'pent'=pent.timegroup, 'N'=list(formula=~1)), 
#               filename=paste0(ifelse(i<10,'mark00','mark0'),i))
#   problems$sing[i] <- length(mod$results$singular[!is.na(mod$results$singular)])
#   problems$se0[i] <- sum(mod$results$beta$se==0)
#   problems$npar[i] <- mod$results$npar.unadjusted
# }

# Compute derived parameters, abundance and density
Nhat.upl <- popan.derived(uplch, topmod)
upl.pop <- Nhat.upl$N
upl.pop$D <- upl.pop$D.se <- upl.pop$D.lcl <- upl.pop$D.ucl <- NA
upl.pop$length <- c(rep(52,27), rep(74,27))
# Compute density with reach length
upl.pop[,c('D','D.lcl','D.ucl')] <- upl.pop[,c('N',"LCL","UCL")]/upl.pop$length
# Rescale standard error for density
upl.pop$D.se <- sqrt(upl.pop$se^2*upl.pop$length^(-2))

# Add some variables
upl.pop$date <- NA
upl.pop$date[upl.pop$reach == 'C'] <- sample_lags$date[sample_lags$site == "UPL" & sample_lags$reach == 'C']
upl.pop$date[upl.pop$reach=='I'] <- sample_lags$date[sample_lags$site=="UPL" & sample_lags$reach=='I']
upl.pop$date <- as.Date(upl.pop$date, origin='1970-01-01')
upl.pop$site  <- 'UPL'
upl.pop$Dlog <- log(upl.pop$D)

# Define differencing operator between control and introduction, ignore first and last samples
oper <- rbind(diag(c(0,rep(1,25),0),nrow=27,ncol=27),diag(c(0,rep(-1,25),0),nrow=27,ncol=27))
# Take the difference
dD.u <- t(oper) %*% upl.pop$Dlog
# Scale vcv matrix of abundance with length of each reach
Dvcv <- Nhat.upl$N.vcv * rbind(cbind(matrix(52^-2,nrow=27,ncol=27), matrix((52*74)^-1,nrow=27,ncol=27)),cbind(matrix((52*74)^-1,nrow=27,ncol=27), matrix(74^-2,nrow=27,ncol=27)))

# create a list of function to apply the delta method
dml <- list(~log(x1), ~log(x2), ~log(x3), ~log(x4),~log(x5), ~log(x6), ~log(x7), ~log(x8), ~log(x9), ~log(x10),~log(x11), ~log(x12), ~log(x13), ~log(x14), ~log(x15), ~log(x16), ~log(x17), ~log(x18), ~log(x19), ~log(x20), ~log(x21), ~log(x22), ~log(x23), ~log(x24), ~log(x25), ~log(x26), ~log(x27), ~log(x28), ~log(x29), ~log(x30), ~log(x31), ~log(x32), ~log(x33), ~log(x34), ~log(x35), ~log(x36), ~log(x37), ~log(x38), ~log(x39), ~log(x40), ~log(x41), ~log(x42), ~log(x43), ~log(x44), ~log(x45), ~log(x46), ~log(x47), ~log(x48), ~log(x49), ~log(x50), ~log(x51), ~log(x52), ~log(x53), ~log(x54))
# Apply delta method to find log scale vcv of density
Dlog.vcv <- msm::deltamethod(dml, mean = upl.pop$D, Dvcv, ses=F)

# Calculate density difference vcv
dD.u.vcv <- t(oper) %*% Dlog.vcv %*% oper
dD.se <- sqrt(diag(dD.u.vcv))
dD.u.vcv <- dD.u.vcv[2:26,2:26]
dD.uu <- data.frame(site='UPL', dD=dD.u, 'se'=dD.se)
dD.uu$date <- upl.pop$date[1:27]

######################################################
### Assemble results #################################
######################################################

dD <- rbind(dD.cc, dD.tt, dD.ll, dD.uu)
# Remove first and last observations from each stream
dD <- dD[dD$dD!=0,]

# Add additional variables for plotting
dD$facet <- 1
dD$facet[dD$site %in% c('LOL','UPL')] <- 2
dD$canopy <- 'closed'
dD$canopy[dD$site %in% c('TAY','UPL')] <- 'open'
dD$canopy <- factor(dD$canopy)
dD$pair <- factor(dD$facet)
dD$gupin <- NA
dD$gupin[dD$site%in% c('CAI','TAY')] <- as.Date('2009-03-18')
dD$gupin[dD$site %in% c('LOL','UPL')] <- as.Date('2008-03-28')
dD$guptime <- (as.numeric(dD$date)-dD$gupin)/30
dD$tsi <- ifelse(dD$guptime<0,0,dD$guptime)

facet_names <- c('2'='LOL & UPL','1' = 'CAI & TAY')

# Define matrices of zeros to construct vcv of all streams
mat0ct1 <- matrix(0, nrow=16,ncol=16)
mat0ct2 <- matrix(0, nrow=25,ncol=16)
mat0lu1 <- matrix(0, nrow=25, ncol=25)
mat0lu2 <- matrix(0, nrow=16, ncol=25)
dD.vcv <- cbind(rbind(dD.c.vcv,mat0ct1, mat0ct2, mat0ct2), 
                rbind(mat0ct1,dD.t.vcv,mat0ct2,mat0ct2), 
                rbind(mat0lu2,mat0lu2,dD.l.vcv,mat0lu1), 
                rbind(mat0lu2,mat0lu2, mat0lu1, dD.u.vcv))


# Fit variance components model with reml
vcomp.final = var.components.reml(theta=dD$dD, design=model.matrix(~pair*canopy*tsi,dD), vcv=dD.vcv)

# Compute the linear predictor and sd on log scale
dD$linpred = model.matrix(~pair*canopy*tsi, dD) %*% as.matrix(vcomp.final$beta$Estimate,ncol=1)
predvcv = model.matrix(~pair*canopy*tsi,dD)%*%as.matrix(vcomp.final$vcv.beta)%*%t(model.matrix(~pair*canopy*tsi,dD))
dD$predsd = sqrt(diag(predvcv))

# Transform to density scale
dD$resp = exp(-dD$linpred)*(1+vcomp.final$sigmasq/2)
dD$respse = sqrt(exp(-2*dD$linpred)*(1+vcomp.final$sigmasq/2)^2)*dD$predsd

# Make Figure 2
tiff('Fig2.tif', height = 5600, width = 8000, units = 'px', res = 800, compression = 'lzw')
ggplot(dD, aes(guptime, exp(-dD), color=site, shape=site)) +
  geom_point(size=2, alpha = 0.25) + 
  geom_errorbar(aes(ymin = exp(-dD-se), ymax = exp(-dD+se)), width=1, alpha = 0.25) + 
  geom_vline(aes(xintercept=0),linetype='dashed',color='red')  + 
  geom_line(aes(y = resp, color = site, linetype = site))+ 
  geom_ribbon(aes(ymin = exp(-linpred+qnorm(.025)*predsd)*(1+vcomp.final$sigmasq/2),
                  ymax = exp(-linpred+qnorm(.975)*predsd)*(1+vcomp.final$sigmasq/2), fill=site), alpha=.2, color=NA) + 
  xlab('Months since introduction') + 
  ylab(expression(~frac('Introduction','Control')~ 'Density ( '~frac('fish',m)~')')) +  
  scale_color_manual('',values=c('CAI'='darkgoldenrod1','TAY'='darkmagenta','LOL'='darkgoldenrod1','UPL'='darkmagenta')) +  
  scale_fill_manual('',values=c('CAI'='darkgoldenrod1','TAY'='darkmagenta','LOL'='darkgoldenrod1','UPL'='darkmagenta')) + 
  scale_shape_manual('', values=c('CAI'=16,'TAY'=1,'LOL'=17,'UPL'=2)) + 
  scale_linetype_manual('', values = c('CAI' = 'dashed', 'TAY' = 'solid', 'LOL' = 'dashed', 'UPL' = 'solid')) +
  facet_wrap(~facet,ncol=1,nrow=2, labeller= as_labeller(facet_names), scales='free_y') + 
  theme_minimal() +
  theme(panel.spacing=unit(2,'lines'), 
        strip.text=element_text(face='bold',size=16), 
        legend.text=element_text(size=13), 
        axis.text=element_text(size=13), 
        axis.title=element_text(size=20))
dev.off()

# Compute slope of Intro/Control density for each stream, overall and contrasts
stream.slope <- matrix(vcomp.final$beta$Estimate, nrow=1) %*% cbind(c(rep(0,3), 1, rep(0, 4)), c(rep(0,3),1,0,1,0,0), c(rep(0,3),1,0,0,1,0),c(rep(0,3),1,0,1,1,1))
slope.vcv <- t(cbind(c(rep(0,3), 1, rep(0, 4)), c(rep(0,3),1,0,1,0,0), c(rep(0,3),1,0,0,1,0),c(rep(0,3),1,0,1,1,1))) %*% as.matrix(vcomp.final$vcv.beta) %*% cbind(c(rep(0,3), 1, rep(0, 4)), c(rep(0,3),1,0,1,0,0), c(rep(0,3),1,0,0,1,0),c(rep(0,3),1,0,1,1,1))

# Contrast and marginal operator matrix
oper.slope <- cbind('LOL-CAI'=c(1,-1,0,0), 'LOL-UPL'=c(1,0,-1,0), 'LOL-TAY'=c(1,0,0,-1), 'CAI-UPL'=c(0,1,-1,0),'CAI-TAY'=c(0,1,0,-1),'UPL-TAY'=c(0,0,1,-1),'canopy'=c(1,1,-1,-1),'pair'=c(1,-1,1,-1),'pairxcanopy'=c(1,-1,-1,1), 'marg'=c(.25,.25,.25,.25), 'LOL'=c(1,0,0,0), 'CAI'=c(0,1,0,0), 'UPL'=c(0,0,1,0), 'TAY'=c(0,0,0,1))
dslope <- stream.slope %*% oper.slope
dslope.vcv <- t(oper.slope) %*% slope.vcv %*% oper.slope
data.frame(est = -c(dslope), std.dev = sqrt(diag(dslope.vcv)), chi2 = c(dslope)^2/diag(dslope.vcv), pval = pchisq(c(dslope)^2/diag(dslope.vcv), df=1, lower.tail=F))

