# Load packages
library(tidyverse)
library(parallel)
library(nlme)
library(emmeans)

# Set number of cores for parallel processing
nodes <- detectCores() - 1

# Read data and rename columns for brevity
rec.sum <- read.csv(file.path('dryad_size_distribution_dataset.csv')) %>%
  rename('pair' = 'stream_pair', 
         'canopy' = 'canopy_treatment', 
         'hab' = 'habitat', 
         'tsi.m' = 'time_since_introduction', 
         'ncap' = 'no_captured',
         'mean_tl' = 'mean_total_length') %>%
  mutate(tsitr = ifelse(tsi.m < 0, 0, tsi.m),
         season = factor(season),
         pair = factor(pair),
         rhab = factor(paste(reach, hab)),
         srhab = factor(paste(site, reach, hab)))


### Begin model fitting
# Fit base model with all interactions
ms0 <- gls(mean_tl~tsitr*hab*reach*pair*canopy, data = rec.sum)
# Update model include autocorrelation structure within each time series
ms0vc <- update(ms0, correlation=corCAR1(form=~tsi.m|srhab))
# Likelihood ratio test
anova(ms0,ms0vc)
# Work with autocorrelation model
ms0 <- ms0vc

# Examine residuals of model
res <- residuals(ms0,'pearson')
fit <- fitted(ms0)
plot(fit,res)
plot(rec.sum$tsitr,res)
boxplot(res~rec.sum$reach)
boxplot(res~rec.sum$hab)
boxplot(res~rec.sum$season)
boxplot(res~rec.sum$pair)
boxplot(res~rec.sum$canopy)
boxplot(res~rec.sum$reach*rec.sum$hab*rec.sum$season*rec.sum$site)
# Residual correlations across reaches and habitats of each steram
error <- cbind(rec.sum, res)
errorw <- spread(error[,c('srhab','tsi.m','res')], key=srhab, res)
round(cor(errorw[,-1],use='pairwise.complete.obs'),2)

# Set-up all possible combinations of variance models
varcombos <- c(combn(c('reach','season','hab','pair','canopy'),1,FUN=paste,simplify=T,collapse='*'), # Single term models
               combn(c('reach','season','hab','pair','canopy'),2,FUN=paste,simplify=T,collapse='*'), # 2-way interactions
               combn(c('reach','season','hab','pair','canopy'),3,FUN=paste,simplify=T,collapse='*'), # 3-way interactions
               combn(c('reach','season','hab','pair','canopy'),4,FUN=paste,simplify=T,collapse='*'), # 4-way interactions
               combn(c('reach','season','hab','pair','canopy'),5,FUN=paste,simplify=T,collapse='*')) # 5-way interactions

# Set-up empty lists to hold the variance formulae
varIDlist <- list()
varexplist <- list()
varfitlist <- list()
varcomblist <- list()

# Populate variance formulae
varexplist[[1]]  <- varExp(form=~tsitr)
varfitlist[[1]] <- varExp(form=~fitted(.))

for (i in 1:length(varcombos)){
  varIDlist[[i]] <- varIdent(form=as.formula(paste0('~1|',varcombos[i])))
  varexplist[[i+1]] <- varExp(form=as.formula(paste0('~tsitr|',varcombos[i])))
  varfitlist[[i+1]] <- varExp(form=as.formula(paste0('~fitted(.)|',varcombos[i])))
  varcomblist[[i]] <- varComb(varIDlist[[i]],varexplist[[1]])
  varcomblist[[i+length(varcombos)]] <- varComb(varIDlist[[i]], varfitlist[[1]])
}

# Concatenate into single list
varstruc <- do.call(c,list(varIDlist,varexplist, varfitlist, varcomblist))

# Fit variance structures
# Set-up parallel cluster
cl <- makeCluster(nodes,type='FORK')
# Export model and variance formulae to cores
clusterExport(cl,list('ms0','varstruc'))

# Fit model with the variance formulae in parallel
ms0list <- mclapply(1:length(varstruc), 
                    FUN= function(x){
                      tryCatch(update(ms0,weights=varstruc[[x]]), error=function(err) NA)
                    }, mc.cores=nodes)
# Shut down cluster
stopCluster(cl)
# Remove list elements that are not model objects with class 'gls' (from models that failed to fit)
ms0list <- ms0list[which(lapply(ms0list,class)=='gls')]
# Find the top models in the list
bbmle::AICtab(ms0list, base=T, weights=T)
# Compare to base model
bbmle::AICtab(ms0,ms0list[[24]])
# Examine top model
summary(ms0list[[24]])

# Residual inspection
res <- residuals(ms0list[[24]],type='pearson')
fit <- fitted(ms0list[[24]], type='pearson')
plot(fit,res)
plot(rec.sum$tsitr,res)
boxplot(res~rec.sum$reach)
boxplot(res~rec.sum$hab)
boxplot(res~rec.sum$season)
boxplot(res~rec.sum$pair)
boxplot(res~rec.sum$canopy)
boxplot(res~rec.sum$search)
ggplot(data= data.frame(res), aes(x=res))+geom_histogram()
ggplot(data= data.frame(res))+geom_qq(aes(sample=res))

# Residual correlations across reaches and habitats within each stream
error <- cbind(rec.sum, res)
error$srhab <- paste(error$site, error$reach, error$hab)
errorw <- spread(error[,c('srhab','tsi.m','res')], key=srhab, res)
apply(errorw[,-1], 2, acf,lag.max=12,na.action=na.pass,ylim=c(-2,2))
round(cor(errorw[,-1],use='pairwise.complete.obs'),2)

# Coeffecients for Table S4
summary(ref_grid(ms0list[[24]], at=list(tsitr=0))) %>% 
  arrange(pair, canopy, reach) %>% 
  mutate(prediction = round(prediction,1), SE=round(SE, 2))

summary(emtrends(ms0list[[24]], ~pair*canopy*reach*hab, var='tsitr', mode = 'df.error')) %>% 
  arrange(pair, canopy, reach) %>% 
  mutate(tsitr.trend = round(tsitr.trend,2), SE=round(SE, 3))

# Contrast analysis of temporal trend in mean body size
emtrends(ms0list[[24]], pairwise~reach*hab, var='tsitr', mode = 'df.error')
emtrends(ms0list[[24]], pairwise~reach*hab|pair*canopy, var='tsitr', mode = 'df.error')

emtrends(ms0list[[24]], ~reach*hab, var='tsitr', mode = 'df.error',
         contr = list('hab' = c(1,1,-1,-1), 
                    'reach' = c(1,-1,1,-1), 
                    'reachxhab' =c(1,-1,-1,1)))
emtrends(ms0list[[24]], ~reach*hab*canopy, var='tsitr', mode = 'df.error',
         contr = list('habxcanopy'=c(1,1,-1,-1,-1,-1, 1, 1), 
                    'reachxcanopy' = c(1,-1,1,-1, -1,1,-1,1),
                    'habxreachxcanopy'=c(1,-1,-1,1,-1,1,1,-1)))
emtrends(ms0list[[24]], ~reach*hab*pair, var='tsitr', mode = 'df.error',
         contr = list('habxpair'=c(1,1,-1,-1,-1,-1, 1, 1), 
                    'reachxpair' = c(1,-1,1,-1, -1,1,-1,1),
                    'habxreachxpair'=c(1,-1,-1,1,-1,1,1,-1)))

emtrends(ms0list[[24]], ~reach*hab*pair*canopy, var='tsitr', mode = 'df.error',
         contr = list('habxpairxcanopy'=c(1,1,-1,-1,-1,-1, 1, 1,-1,-1, 1, 1, 1, 1, -1,-1), 
                      'reachxpairxcanopy' = c(1,-1,1,-1, -1,1,-1,1, -1,1,-1,1,1,-1,1,-1), 
                      'habxreachxpairxcanopy' =c(1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1,1)))


# Variance estimates for groups
# Including seasonal differences in variance
# group.vars <- coef(ms0list[[24]]$modelStruct$varStruct, unconstrained=F) %>%
#   enframe() %>%
#   separate(name, into=c('season','pair','canopy'), sep='\\\\.|\\\\*') 

# Variance estimates for groups averaging over seasons
# group.vars.aseasonal <- coef(ms0list[[24]]$modelStruct$varStruct, unconstrained=F) %>%
#   enframe() %>%
#   separate(name, into=c('season','pair','canopy'), sep='\\\\.|\\\\*') %>% 
#   group_by(pair, canopy) %>% 
#   summarise(A = mean(value))


# Predictions for figure
# Set-up grid for predictions
predframe <- expand.grid('tsitr'=seq(0,max(rec.sum$tsitr),.5), 
                        'hab'=factor(c('P','R')),
                        'reach'=factor(c('C','I')),
                        'pair'=factor(c('1','2')), 
                        'canopy' = factor(c('intact','thinned')), 
                        'season'=mean(as.numeric(unique(rec.sum[,c('tsi.m','season')])[,2]))-1) # Use mean seasonality for predictions

# Define design matrix
Xmat <- model.matrix(formula(ms0list[[24]])[-2], predframe)
# Extract model coeffecients
betamat <- coef(ms0list[[24]])
# Find predictions
predmat <- Xmat %*% betamat
# Store in predframe
predframe$pred <- predmat
# Find vcov of predicitons
predvcv <- Xmat %*% vcov(ms0list[[24]]) %*% t(Xmat)
predframe$beta.sd <- sqrt(diag(predvcv))
predframe$low <- predframe$pred + qnorm(0.025)*predframe$beta.sd
predframe$high <- predframe$pred + qnorm(0.975)*predframe$beta.sd

# Project estimates at time 0 to pre-guppy samples and add to predictions
predframe <- predframe %>% 
  filter(tsitr==0) %>% 
  mutate(tsitr = ifelse(pair==1, -11, -2.5)) %>% 
  bind_rows(predframe)

# Add site variable
predframe$site = NA
predframe$site[predframe$pair=='2'&predframe$canopy=='intact'] <- 'CAI'
predframe$site[predframe$pair=='2'&predframe$canopy=='thinned'] <- 'TAY'
predframe$site[predframe$pair=='1'&predframe$canopy=='intact'] <- 'LOL'
predframe$site[predframe$pair=='1'&predframe$canopy=='thinned'] <- 'UPL'

# Add reach-habitat variable
predframe$rhab <- paste(predframe$reach, predframe$hab)

# Clip to end point for second stream pair
predframe <- predframe[!(predframe$site %in% c('CAI','TAY') & predframe$tsitr>37),]

# Construct figure
tiff('Fig6.tif', height = 5600, width = 8000, units ='px', res=800, compression = 'lzw')
ggplot(data=predframe, aes(x=tsitr,y=pred,color=rhab,shape=rhab)) + 
  geom_point(aes(x=tsi.m,y=mean_tl),data=rec.sum, alpha = 0.25) +
  geom_line(aes(linetype = rhab)) +
  geom_ribbon(aes(ymin = low, ymax=high, fill=rhab), color = NA, alpha=.1) +
  scale_shape_manual('', values=c('C P'=16,'C R'=2,'I P'=16,'I R'=2), 
                     labels=c('control pool','control riffle','introduction pool','introduction riffle')) + 
  scale_linetype_manual('',values=c('C P'='solid','C R'='dashed','I P'='solid','I R'='dashed'), 
                        labels=c('control pool','control riffle','introduction pool','introduction riffle')) +
  scale_fill_manual('',values=c('C P'='blue','C R'='blue','I P'='red','I R'='red'), 
                    labels=c('control pool','control riffle','introduction pool','introduction riffle')) +
  scale_color_manual('', values=c('C P'='blue','C R'='blue','I P'='red','I R'='red'), 
                     labels=c('control pool','control riffle','introduction pool','introduction riffle')) +
  theme_minimal() +
  theme(panel.spacing=unit(2,'lines'), 
        strip.text=element_text(face='bold',size=16), 
        legend.text=element_text(size=13), 
        axis.text=element_text(size=13), 
        axis.title=element_text(size=20))+
  xlab('Months since introduction') +
  ylab('Mean TL (mm)') +
  facet_wrap(~site, scales='free') + 
  scale_x_continuous(breaks =c(-10, 0, 15, 30, 45)) + 
  scale_y_continuous(breaks = c(40, 50,60))
dev.off()

# Percent change over the course of the study
# Extract the rows of the design matrix for start and end points
Xchange <- matrix(Xmat[Xmat[,'tsitr']==0 | Xmat[,'tsitr'] == ifelse(Xmat[,'pair2']==1, 36, 48)], ncol=ncol(Xmat), 
                 dimnames=list(NULL, colnames(Xmat)))
# Find predictions at those points
pchange <- Xchange %*% betamat
# Find vcov at those points
pchange.vcv <- Xchange %*% vcov(ms0list[[24]]) %*% t(Xchange)

# Define difference operator for end - start points
oper.diff <- diag(-1,nrow=nrow(Xchange),ncol=nrow(Xchange))[,seq(1,nrow(Xchange),2)] + diag(1,nrow=nrow(Xchange),ncol=nrow(Xchange))[,seq(2,nrow(Xchange),2)]
# Compute difference
diff <- t(oper.diff) %*% pchange
# Compute variance of teh difference
diff.vcv <- t(oper.diff) %*% pchange.vcv %*% oper.diff
# Divide by start value
pdiff <- diff/pchange[seq(1,nrow(Xchange),2)]
# Adjust vcv to reflect division
pdiff.vcv <- diff.vcv * (1/pchange[seq(1,nrow(Xchange),2)]) %*% t(1/pchange[seq(1,nrow(Xchange),2)])
# Find se of % difference
pdiff.se <- sqrt(diag(pdiff.vcv))
# Put everything in a data.frame
perc <- data.frame(Xchange[seq(1,nrow(Xchange),2),c('habR','reachI','pair2','canopythinned')], pdiff=pdiff, se=pdiff.se)
perc$habR <- ifelse(perc$habR==1,'R','P')
perc$reachI <- ifelse(perc$reachI==1,'I','C')
perc$pair2 <- ifelse(perc$pair2==1,2,1)
perc$canopythinned <- ifelse(perc$canopythinned==1,'thinned','intact')
# Define averaging operator matrix
oper.avg <- cbind(rep(c(1/4,0,0,0),4), rep(c(0,1/4,0,0),4), rep(c(0,0,1/4,0),4), rep(c(0,0,0,1/4),4))
# Compute average % change
avg.perc <- t(oper.avg) %*% pdiff
avg.perc.vcv <- t(oper.avg) %*% pdiff.vcv %*% oper.avg
ans <- data.frame(hab=rep(c('P','R'),2), reach=rep(c('C','I'),each=2), pchange=avg.perc,se = sqrt(diag(avg.perc.vcv)))

