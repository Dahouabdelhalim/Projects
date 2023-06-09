# Load packages
library(tidyverse)
library(parallel)
library(nlme)
library(emmeans)
library(chron)
# Set number of cores for parallel processing
nodes <- detectCores() - 1

# Read data and rename columns for brevity
newmark.sum <- read.csv(file.path('dryad_recruit_dataset.csv')) %>%
  rename('pair' = 'stream_pair', 
         'canopy' = 'canopy_treatment', 
         'hab' = 'habitat', 
         'tsi.m' = 'time_since_introduction', 
         'search' = 'search_period', 
         'len' = 'stream_length', 
         'N' = 'recruit_n', 
         'adult' = 'adult_n', 
         'D' = 'recruit_density', 
         'lnD' = 'log_recruit_density') %>%
  mutate(tsitr = ifelse(tsi.m < 0, 0, tsi.m),
         pair = factor(pair),
         pcrh = factor(paste(pair, canopy, reach, hab)),
         search = factor(search),
         season = factor(season))

### Begin model fitting
# Fit base model with all interactions
mr0 <- gls(lnD~tsitr*reach*hab*pair*canopy+adult, data=newmark.sum)
# Update model include autocorrelation structure within each time series
mr0ac <- update(mr0, correlation = corCAR1(form=~tsi.m|pcrh))
# Likelihood ratio test
anova(mr0, mr0ac)
# Work with autocorrelation model
mr0 <- mr0ac

# Examine residuals of model
res <- residuals(mr0,type='pearson')
fit <- fitted(mr0, type='pearson')
plot(fit,res)
plot(newmark.sum$tsitr,res)
plot(newmark.sum$adult,res, col=as.numeric(newmark.sum$season))
boxplot(res~newmark.sum$reach)
boxplot(res~newmark.sum$hab)
boxplot(res~newmark.sum$season)
boxplot(res~newmark.sum$pair)
boxplot(res~newmark.sum$canopy)
boxplot(res~newmark.sum$search)
boxplot(res~newmark.sum$reach*newmark.sum$hab*newmark.sum$season*newmark.sum$site)
ggplot(data= data.frame(res), aes(x=res))+geom_histogram()
ggplot(data= data.frame(res))+geom_qq(aes(sample=res))

# correlations across reaches and habitats within each stream
error <- cbind(newmark.sum, res)
error$srhab <- paste(error$site, error$reach, error$hab)
errorw <- spread(error[,c('srhab','tsi.m','res')], key=srhab, res)
apply(errorw[,-1], 2, acf,lag.max=12,na.action=na.pass,ylim=c(-2,2))
round(cor(errorw[,-1],use='pairwise.complete.obs'),2)

# Set-up all possible combinations of variance models
varcombos <- c(combn(c('reach','season','hab','pair','canopy','search'),1,FUN=paste,simplify=T,collapse='*'), # Single term models
              combn(c('reach','season','hab','pair','canopy','search'),2,FUN=paste,simplify=T,collapse='*'), # 2-way interactions
              combn(c('reach','season','hab','pair','canopy','search'),3,FUN=paste,simplify=T,collapse='*'), # 3-way interactions
              combn(c('reach','season','hab','pair','canopy','search'),4,FUN=paste,simplify=T,collapse='*'), # 4-way interactions
              combn(c('reach','season','hab','pair','canopy','search'),5,FUN=paste,simplify=T,collapse='*'), # 5-way interactions
              combn(c('reach','season','hab','pair','canopy','search'),6,FUN=paste,simplify=T,collapse='*')) # 6-way interaction

# Set-up empty lists to hold the variance formulae
varIDlist <- list()
varexplist <- list()
varadlist <- list()
varfitlist <- list()
varcomblist <- list()

# Populate variance formulae
varexplist[[1]]  <- varExp(form=~tsitr)
varadlist[[1]] <- varExp(form=~adult)
varfitlist[[1]] <- varExp(form=~fitted(.))

for (i in 1:length(varcombos)){
  varIDlist[[i]] <- varIdent(form=as.formula(paste0('~1|',varcombos[i])))
  varexplist[[i+1]] <- varExp(form=as.formula(paste0('~tsitr|',varcombos[i])))
  varadlist[[i+1]] <- varExp(form=as.formula(paste0('~adult|',varcombos[i])))
  varfitlist[[i+1]] <- varExp(form=as.formula(paste0('~fitted(.)|',varcombos[i])))
  varcomblist[[i]] <- varComb(varIDlist[[i]],varexplist[[1]])
  varcomblist[[i+length(varcombos)]] <- varComb(varIDlist[[i]], varfitlist[[1]])
  varcomblist[[i+2*length(varcombos)]] <- varComb(varIDlist[[i]], varadlist[[1]])
}

# Concatenate into single list
varstruc <- do.call(c,list(varIDlist,varexplist,varadlist, varfitlist, varcomblist))
# Fit variance structures
# Set-up parallel cluster
cl <- makeCluster(nodes,type='FORK')
# Export model object and variance formulae to each core
clusterExport(cl,list('mr0','varstruc'))

# Fit model with the variance formulae in parallel
mr0list <- mclapply(1:length(varstruc), 
                    FUN=function(x){
                      tryCatch(update(mr0,weights=varstruc[[x]]), error=function(err) NA)},
                    mc.cores=nodes)
# Shut down cluster
stopCluster(cl)
# Remove list elements that are not model objects with class 'gls' (from models that failed to fit)
mr0list <- mr0list[which(lapply(mr0list,class)=='gls')]
# Find the top models in the list
bbmle::AICtab(mr0list, base=T, weights=T)
# Compare to base model
bbmle::AICtab(mr0list[[304]], mr0, base = T, weights = T)
# Refit final model for convenience
mr0final <- update(mr0a, weights=varComb(varIdent(form = ~1|reach*season*hab*canopy*search), varExp(form=~fitted(.))))

# Examine top model
summary(mr0final)

# Residual diagnostics
res <- residuals(mr0final,type='pearson')
fit <- fitted(mr0final, type='pearson')
plot(fit,res)
plot(newmark.sum$tsitr,res)
plot(newmark.sum$adult,res, col=as.numeric(newmark.sum$season))
boxplot(res~newmark.sum$reach)
boxplot(res~newmark.sum$hab)
boxplot(res~newmark.sum$season)
boxplot(res~newmark.sum$pair)
boxplot(res~newmark.sum$canopy)
boxplot(res~newmark.sum$search)
boxplot(res~newmark.sum$reach*newmark.sum$hab*newmark.sum$season*newmark.sum$site)
ggplot(data= data.frame(res), aes(x=res))+geom_histogram()
ggplot(data= data.frame(res))+geom_qq(aes(sample=res))
# correlations across reaches and habitats within each stream
error <- cbind(newmark.sum, res)
error$srhab <- paste(error$site, error$reach, error$hab)
errorw <- spread(error[,c('srhab','tsi.m','res')], key=srhab, res)
apply(errorw[,-1], 2, acf,lag.max=12,na.action=na.pass,ylim=c(-2,2))
round(cor(errorw[,-1],use='pairwise.complete.obs'),2)

# Coeffecients for Table S2
# Intercept
summary(ref_grid(mr0final, at=list(tsitr=0))) %>% 
  arrange(pair, canopy, reach) %>% 
  mutate(prediction = round(prediction,2), SE=round(SE, 3))
# Slope with respect to time
summary(emtrends(mr0list[[304]], ~pair*canopy*reach*hab, var='tsitr', mode='df.error')) %>% 
  arrange(pair, canopy, reach) %>% 
  mutate(tsitr.trend = round(tsitr.trend,3), SE=round(SE, 4))

# Contrast analysis of temporal trend
emtrends(mr0final, pairwise~reach*hab, var='tsitr', mode = 'df.error', side = '=', adjust='none')
emtrends(mr0final, ~reach*hab, var='tsitr',
         contr=list('hab' = c(1,1,-1,-1), 
                    'reach' = c(1,-1,1,-1), 
                    'reachxhab' = c(1,-1,-1,1)), 
         mode='df.error', side = '>')
emtrends(mr0final, ~reach*hab*canopy, var='tsitr',
         contr=list('habxcanopy' = c(1,1,-1,-1,-1,-1, 1, 1), 
                    'reachxcanopy' = c(1,-1,1,-1, -1,1,-1,1), 
                    'habxreachxcanopy'=c(1,-1,-1,1,-1,1,1,-1)), 
         mode='df.error')
emtrends(mr0final, ~reach*hab*pair, var='tsitr',
         contr=list('habxpair' = c(1,1,-1,-1,-1,-1, 1, 1), 
                    'reachxpair' = c(1,-1,1,-1, -1,1,-1,1),
                    'habxreachxpair'=c(1,-1,-1,1,-1,1,1,-1)), 
         mode='df.error')
emtrends(mr0final, ~reach*hab*pair*canopy, var='tsitr',
         contr=list('habxpairxcanopy' = c(1,1,-1,-1,-1,-1, 1, 1,-1,-1, 1, 1, 1, 1, -1,-1), 
                    'reachxpairxcanopy' = c(1,-1,1,-1, -1,1,-1,1, -1,1,-1,1,1,-1,1,-1), 
                    'habxreachxpairxcanopy' =c(1,-1,-1,1,-1,1,1,-1,-1,1,1,-1,1,-1,-1,1)), 
         mode='df.error')

# Calculate estimated % change over course of study
# Set-up grid of points across groups to make predictions
predframe.recdelta <- expand.grid('tsitr' = seq(0,max(newmark.sum$tsitr),.5), 
                                  'hab' = factor(c('P','R')),
                                  'reach' = factor(c('C','I')),
                                  'pair' = factor(c('1','2')), 
                                  'canopy' = factor(c('intact','thinned')), 
                                  'season'=mean(as.numeric(unique(newmark.sum[,c('tsi.m','season')])[,2]))-1, # Use the average seasonality over the study period
                                  'search' = factor(c(0,1)), 
                                  adult=mean(unique(newmark.sum[,c('tsi.m','adult')])[,2])) # use the average number of adult captures over the study period
# Take only the introduction and end point of each study site
predframe.recdelta <- predframe.recdelta[(predframe.recdelta$tsitr == 0 | 
                                            predframe.recdelta$tsitr == ifelse(predframe.recdelta$pair == 1, 48, 36)) & 
                                           predframe.recdelta$search == ifelse(predframe.recdelta$pair == 1 & predframe.recdelta$tsitr == 0, 0, 1),]
# Order the remaining rows
predframe.recdelta <- predframe.recdelta[order(predframe.recdelta$pair, predframe.recdelta$canopy, predframe.recdelta$reach, predframe.recdelta$hab),]
# Construct model design matrix
Xmat.cond <- model.matrix(formula(mr0final)[-2], predframe.recdelta)
# Extract model coefficients
beta.cond <- coef(mr0final)
# Find predictions
pred.cond <- Xmat.cond %*% beta.cond
# compute variance covariance of predictions
pred.vcv <- Xmat.cond %*% vcov(mr0final) %*% t(Xmat.cond)

# Set-up matrix operator to take difference of start and end point for each study group
oper.diff <- diag(-1, nrow=nrow(Xmat.cond), ncol=nrow(Xmat.cond))[,seq(1,nrow(Xmat.cond),2)] + 
  diag(1, nrow=nrow(Xmat.cond),ncol=nrow(Xmat.cond))[,seq(2,nrow(Xmat.cond),2)]
# Compute the difference and variance-covariance of the predictions
diff <- t(oper.diff) %*% pred.cond
diff.vcv <- t(oper.diff) %*% pred.vcv %*% oper.diff

# Combine the design data with the differences and standard errors
ratio <- cbind(predframe.recdelta[seq(1,nrow(predframe.recdelta),2),],diff=diff, se=sqrt(diag(diff.vcv)))

# Define a matrix operator to average across the different streams
oper.avg <- cbind(rep(c(1/4,0,0,0),4), rep(c(0,1/4,0,0),4), rep(c(0,0,1/4,0),4), rep(c(0,0,0,1/4),4))
# Find the average percent change and variance-covariance in each reach-habitat
avg.perc <- t(oper.avg) %*% diff
avg.perc.vcv <- t(oper.avg) %*% diff.vcv %*% oper.avg
# Table with percent change and 95% confidence bounds
ans <- data.frame(hab=rep(c('P','R'),2), 
                  reach=rep(c('C','I'),each=2), 
                  pchange=exp(avg.perc)-1,
                  lcl = exp(avg.perc+qnorm(.025)*sqrt(diag(avg.perc.vcv)))-1, 
                  ucl = exp(avg.perc+qnorm(.975)*sqrt(diag(avg.perc.vcv)))-1) %>% 
  mutate(pchange = round(100*pchange,1), 
         lcl = round(100*lcl,2), 
         ucl = round(100*ucl,2))

# Variance estimates for groups
# Including seasonal differences in variance
# group.vars <- coef(mr0final$modelStruct$varStruct, unconstrained=F) %>% 
#   enframe() %>% 
#   separate(name, into=c('prefix','reach','season','hab','canopy','search'), sep='\\\\.|\\\\*') %>% 
#   pivot_wider(id_cols= c('reach','season','hab','canopy','search'), names_from = prefix) %>% 
#   fill(B, .direction='up') %>% 
#   filter(reach != 'expon') 

# Variance estimates for groups averaging over seasons
group.vars.aseasonal <- coef(mr0final$modelStruct$varStruct, unconstrained=F) %>% 
  enframe() %>% 
  separate(name, into=c('prefix','reach','season','hab','canopy','search'), sep='\\\\.|\\\\*') %>% 
  pivot_wider(id_cols= c('reach','season','hab','canopy','search'), names_from = prefix) %>% 
  fill(B, .direction='up') %>% 
  filter(reach != 'expon') %>% 
  group_by(reach, hab, canopy, search) %>% 
  summarise(A = mean(A), 
            B = mean(B))

# Predictions for figure
# Set-up grid for predictions
predframe.rec <- expand.grid('tsitr' = seq(0,max(newmark.sum$tsitr),.1), 
                            'hab' = factor(c('P','R')),
                            'reach' = factor(c('C','I')),
                            'pair' = factor(c('1','2')), 
                            'canopy' = factor(c('intact','thinned')), 
                            adult = mean(unique(newmark.sum[,c('tsi.m','adult')])[,2])) # Use average adult capture numbers

# Define design matrix
Xmat <- model.matrix(formula(mr0final)[-2], predframe.rec)
# Extract model coefficients
betamat <- coef(mr0final)
# Make predictions
predframe.rec$pred <- c(Xmat %*% betamat)
# Find vcov of predicitons
predvcv <- Xmat %*% vcov(mr0final) %*% t(Xmat)
predframe.rec$beta.var <- diag(predvcv)

# Project Estimates at time 0 to pre-guppy samples and add to predictions
# Then add previously calculated variance covariates and calculate expectations, SDs and CIs
predframe.rec <- predframe.rec %>% 
  filter(tsitr==0) %>% 
  dplyr::select(-tsitr) %>% 
  right_join(expand.grid('tsitr' = seq(-11,-0.5,0.1), 
                         'hab'=factor(c('P','R')),
                         'reach'=factor(c('C','I')),
                         'pair'=factor(c('1','2')), 
                         'canopy' = factor(c('intact','thinned')), 
                         adult=mean(unique(newmark.sum[,c('tsi.m','adult')])[,2]))) %>% 
  filter(!(pair==2 & tsitr < -2.5)) %>% 
  bind_rows(predframe.rec) %>% 
  mutate(origin.m = 3, 
         origin.d = ifelse(pair == '1', 28, 18), 
         origin.y = ifelse(pair == '1', 2008, 2009), 
         origin.t = dates(paste(origin.m, origin.d, origin.y, sep='/')), 
         date = dates(tsitr*30) + as.numeric(origin.t), 
         search = ifelse(years(date) >= 2009, '1', '0')) %>% 
  left_join(group.vars.aseasonal) %>% 
  mutate(A = ifelse(is.na(A),1,A), 
         sigma = summary(mr0final)$sigma, 
         A.sigma2 = (A*sigma)^2) %>%
  fill(B, .direction = 'downup') %>% 
  mutate(sigma2 = A.sigma2*exp(2*B*pred)) %>% 
  mutate(expect = exp(pred)*(1+sigma2/2)-1, 
         expect.sd = sqrt(beta.var), 
         low = exp(pred + qnorm(.025)*expect.sd)*(1+sigma2/2)-1, 
         high = exp(pred + qnorm(.975)*expect.sd)*(1+sigma2/2)-1) 


# Add site to predframe.rec
predframe.rec$site <- NA
predframe.rec$site[predframe.rec$pair=='2'&predframe.rec$canopy=='intact'] <- 'CAI'
predframe.rec$site[predframe.rec$pair=='2'&predframe.rec$canopy=='thinned'] <- 'TAY'
predframe.rec$site[predframe.rec$pair=='1'&predframe.rec$canopy=='intact'] <- 'LOL'
predframe.rec$site[predframe.rec$pair=='1'&predframe.rec$canopy=='thinned'] <- 'UPL'
# Clip to end point for second stream pair
predframe.rec <- predframe.rec[!(predframe.rec$site %in% c('CAI','TAY') & predframe.rec$tsitr>37),]

# Create a reach-habitat variable
predframe.rec$rhab <- paste(predframe.rec$reach, predframe.rec$hab)
newmark.sum$rhab <- paste(newmark.sum$reach, newmark.sum$hab)

# Construct figure
tiff('Fig4.tif', height = 5600, width = 8000, units = 'px', res= 800, compression ='lzw')
ggplot(data=predframe.rec, aes(x=tsitr,y=expect,color=rhab,shape=rhab)) + 
  geom_point(aes(x=tsi.m,y=D),data=newmark.sum, alpha = 0.25) +
  geom_line(aes(linetype = rhab)) +
  geom_ribbon(aes(ymin = low, ymax=high, fill=rhab), alpha=.1, color=NA) +
  scale_shape_manual('',values=c('C P'=16,'C R'=2,'I P'=16,'I R'=2), 
                     labels=c('control pool','control riffle','introduction pool','introduction riffle')) + 
  scale_linetype_manual('', values=c('C P'='solid','C R'='dashed','I P'='solid','I R'='dashed'), 
                        labels=c('control pool','control riffle','introduction pool','introduction riffle')) +
  scale_fill_manual('',values=c('C P'='blue','C R'='blue','I P'='red','I R'='red'), 
                    labels=c('control pool','control riffle','introduction pool','introduction riffle')) +
  scale_color_manual('',values=c('C P'='blue','C R'='blue','I P'='red','I R'='red'), 
                     labels=c('control pool','control riffle','introduction pool','introduction riffle')) +
  theme_minimal() +
  theme(panel.spacing=unit(2,'lines'), 
        strip.text=element_text(face='bold',size=16), 
        legend.text=element_text(size=13), 
        axis.text=element_text(size=13), 
        axis.title=element_text(size=20)) +
  xlab('Months since introduction') +
  ylab(expression('Captured Killifish Recruit Density ('~frac(recruits,m)~')')) +
  facet_wrap(~site, scales='free') + 
  scale_x_continuous(breaks = c(0,15,30,45)) +
  scale_y_continuous(breaks = seq(0,1.5,0.5))
dev.off()


