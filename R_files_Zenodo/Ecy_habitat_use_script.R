# Load packages
library(tidyverse)
library(parallel)
library(nlme)
library(emmeans)
# Set number of cores for parallel processing
nodes <- detectCores() - 1

# Read data  and wrangle into analysis format
kgsum <- read.csv(file.path('dryad_habitat_use_dataset.csv')) %>% 
  rename('pair' = 'stream_pair', 
         'canopy' = 'canopy_treatment', 
         'hab' = 'habitat', 
         'sp' = 'species', 
         'tsi.m' = 'time_since_introduction', 
         'len' = 'stream_length', 
         'ncap' = 'no_captured', 
         'capdens' = 'capture_density', 
         'log.cd' = 'log_capture_density') %>%
  pivot_wider(id_cols = c(site:reach, sp:tsi.m), 
              names_from = hab, values_from = log.cd, values_fill = 0) %>%
  mutate(log.ratio = P - R,
         tsitr = ifelse(tsi.m < 0, 0, tsi.m),
         pair = factor(pair),
         legend = factor(paste(sp, reach)),
         ssr = factor(paste(sp, site, reach))) %>%
  select(-c(P:R))

### Begin model fitting
# Fit base model with all interactions
mh0 <- gls(log.ratio~legend*pair*canopy*tsitr, data=kgsum)
# Update model include autocorrelation structure within each time series
mh0ac <- update(mh0, correlation=corCAR1(form=~tsi.m|ssr))
# Likelihood ratio test  
anova(mh0, mh0ac)
# Work with autocorrelation model
mh0 <- mh0ac
## Examine residuals of model
res = residuals(mh0,'pearson')
fit = fitted(mh0)
plot(fit,res)
plot(kgsum$tsi.m, res, col=as.numeric(kgsum$legend))
boxplot(res~kgsum$reach)
boxplot(res~kgsum$sp)
boxplot(res~kgsum$legend)
boxplot(res~kgsum$site)
boxplot(res~kgsum$ssr)
ggplot(data= data.frame(res), aes(x=res)) +
  geom_histogram() +
  facet_wrap(~kgsum$legend)
ggplot(data= data.frame(res)) +
  geom_qq(aes(sample=res))

# Set-up all possible combinations of variance models
varcombos <- c(combn(c('reach','sp','pair','canopy'),1,FUN=paste,simplify=T,collapse='*'), # Single term models
              combn(c('reach','sp','pair','canopy'),2,FUN=paste,simplify=T,collapse='*'), # 2-way interactions
              combn(c('reach','sp','pair','canopy'),3,FUN=paste,simplify=T,collapse='*'), # 3-way interactions
              combn(c('reach','sp','pair','canopy'),4,FUN=paste,simplify=T,collapse='*')) # 4-way interactions
varcombos <- gsub(pattern='reach\\\\*sp',replacement='legend',varcombos) # Replace reach*sp with legend
# Set-up empty lists to hold the variance formulae
varIDlist <- list()
varexplist <- list()
varfitlist <- list()
varcomblist <- list()

# Populate variance formulae
varexplist[[1]] <- varExp(form=~tsitr)
varfitlist[[1]] <- varExp(form=~fitted(.))

for (i in 1:length(varcombos)){
  varIDlist[[i]] <- varIdent(form=as.formula(paste0('~1|',varcombos[i])))
  varexplist[[i+1]] <- varExp(form=as.formula(paste0('~tsitr|',varcombos[i])))
  varfitlist[[i+1]] <- varExp(form=as.formula(paste0('~fitted(.)|',varcombos[i])))
  varcomblist[[i]] <- varComb(varIDlist[[i]],varexplist[[1]])
  varcomblist[[i+length(varcombos)]] <- varComb(varIDlist[[i]], varfitlist[[1]])
}
# Turn variance formulae into one list
varstruc <- do.call(c,list(varIDlist,varexplist, varfitlist, varcomblist))

# Set-up parallel cluster
cl <- makeCluster(nodes,type='FORK')
# Export model and variance structures to each core
clusterExport(cl,list('mh0','varstruc'))

# Fit model with the variance formulae in parallel
mh0list <- mclapply(1:length(varstruc), 
                 FUN=function(x){
                   tryCatch(update(mh0,weights=varstruc[[x]]), error=function(err) NA)},
                 mc.cores=nodes)
# Shut down cluster
stopCluster(cl)
# Remove list elements that are not model objects with class 'gls' (from models that failed to fit)
mh0list <- mh0list[which(lapply(mh0list,class)=='gls')]
# Find the top models in the list
bbmle::AICtab(mh0list, base = T, weights=T)
# Compare to base model
bbmle::AICtab(mh0list[[77]], mh0, base=T, weights=T)

# Examine top model
summary(mh0list[[77]])

# Residual inspection
res <- residuals(mh0list[[77]],'pearson')
fit <- fitted(mh0list[[77]])
plot(fit,res)
plot(kgsum$tsi.m, res, col=as.numeric(kgsum$legend))
boxplot(res~kgsum$reach)
boxplot(res~kgsum$sp)
boxplot(res~kgsum$legend)
boxplot(res~kgsum$site)
boxplot(res~kgsum$ssr)
ggplot(data= data.frame(res), aes(x=res)) +
  geom_histogram() + 
  facet_wrap(~kgsum$legend)
ggplot(data= data.frame(res)) + 
  geom_qq(aes(sample=res))

#Extract coefficients for variance structure of top model
# Column 'A' will contain the group-wise variance estimates
# Column 'B' will contain the exponential effect of the fitted values on variance
group.vars <- coef(mh0list[[77]]$modelStruct$varStruct, unconstrained=F) %>% 
  enframe() %>% 
  separate(name, into = c('prefix','sp','reach','pair','canopy'), sep = ' |\\\\.|\\\\*') %>% 
  pivot_wider(id_cols = c(sp, reach, pair, canopy), names_from=prefix) %>% 
  fill(B, .direction='up') %>% 
  filter(sp != 'expon')

# Find the fitted values for each sp, reach, pair and canopy combination
# at the mean tsitr
fits.comp <- kgsum %>% 
  group_by(sp, reach) %>% 
  modelr::data_grid(pair, canopy) %>% 
  mutate(legend = paste(sp, reach), 
         tsitr = mean(kgsum$tsitr)) %>% 
  predict(object = mh0list[[77]], newdata = .) %>% 
  data.frame(fit = .)

# Find the residual standard deviation associated with each legend
# at the mean tsitr
sigma.mid <- kgsum %>% 
  group_by(sp, reach) %>% 
  modelr::data_grid(pair, canopy) %>% 
  mutate(legend = paste(sp, reach), 
         tsitr = mean(kgsum$tsitr)) %>% 
  left_join(group.vars) %>% 
  mutate(A = ifelse(is.na(A),1,A), 
         sigma = summary(mh0list[[77]])$sigma, 
         A.sigma2 = (A*sigma)^2) %>% 
  fill(B) %>% 
  bind_cols(fits.comp) %>% 
  mutate(sigma2 = A.sigma2*exp(2*B*fit)) %>%
  summarise(sigma2 = mean(sigma2)) %>%
  mutate(sigma = sqrt(sigma2))

# Coeffecients for Table S1
summary(ref_grid(mh0list[[77]], at=list(tsitr = 0))) %>% 
  arrange(pair, canopy, legend) %>% 
  mutate(prediction = round(prediction,2), 
         SE=round(SE, 3))

summary(emtrends(mh0list[[77]], ~pair*canopy*legend, var='tsitr', mode = 'df.error')) %>% 
  arrange(pair, canopy, legend) %>% 
  mutate(tsitr.trend = round(tsitr.trend,3), 
         SE=round(SE, 4))

# Contrast analysis
# Differences between legend at mean tsi.m
emmeans(update(ref_grid(mh0list[[77]]), tran = 'log'), 
        pairwise~legend, type='response', bias.adj =T, sigma = sigma.mid$sigma, side = '>', adjust='none') 
# or on the log scale
contrast(emmeans(mh0list[[77]], ~legend*tsitr), interaction=c('pairwise'), by='tsitr')
# For each stream...
contrast(emmeans(mh0list[[77]], ~legend*pair*canopy), interaction=c('pairwise','pairwise'), by=c('pair','canopy'))

# Comparison of slopes
emtrends(mh0list[[77]], pairwise~legend, var='tsitr', mode = 'df.error')

# Test of slope across all legend combinations
summary(emtrends(mh0list[[77]], ~1, var='tsitr', mode = 'df.error')) %>% mutate(t.ratio = tsitr.trend/SE, p.value = pt(t.ratio, df, lower.tail=T))
# For each stream...
contrast(emtrends(mh0list[[77]], ~legend*pair*canopy, var='tsitr', mode = 'df.error'), interaction=c('pairwise'), by=c('pair','canopy'))


#### Compute predictions for figure
# Create grid of points across all the treatment groups 
predframe <- expand.grid('tsitr' = seq(0,max(kgsum$tsitr)+1,.5),
                        'legend' = factor(c('Guppy I','Killifish C','Killifish I')),
                        'pair' = factor(c('1','2')), 
                        'canopy' = factor(c('intact','thinned')))

# Find design matrix for the predictions
Xmat <- model.matrix(formula(mh0list[[77]])[-2], predframe)
# Beta estimates
betamat <- coef(mh0list[[77]])
# Find predictions and store
predmat <- Xmat %*% (betamat)
predframe$pred <- c(predmat)
# Find variance of predictions and store
predvcv <- Xmat %*% vcov(mh0list[[77]]) %*% t(Xmat)
predframe$beta.var <- diag(predvcv)

# Use previously computed variances for each group and find adjust variance accordingly to find expectation and confidence interval
predframe <- predframe %>% 
  separate(legend, into=c('sp','reach'), remove = F) %>% 
  left_join(group.vars) %>% 
  mutate(A = ifelse(is.na(A),1,A), 
         sigma = summary(mh0list[[77]])$sigma, 
         A.sigma2 = (A*sigma)^2) %>% 
  fill(B) %>% 
  mutate(sigma2 = A.sigma2*exp(2*B*pred), 
         expect = exp(pred)*(1+sigma2/2), 
         expect.sd = sqrt(beta.var), 
         low = exp(pred + qnorm(.025)*sqrt(beta.var))*(1+sigma2/2), 
         high = exp(pred + qnorm(.975)*sqrt(beta.var))*(1+sigma2/2))

# Project estimates at time 0 to pre-guppy samples and add to predictions
predframe <- predframe %>% 
  filter(tsitr==0) %>% 
  mutate(tsitr = ifelse(pair==1, -11, -2.5)) %>% 
  bind_rows(predframe)

# Add site variable to predictions
predframe$site <- NA
predframe$site[predframe$pair=='2'&predframe$canopy=='intact'] <- 'CAI'
predframe$site[predframe$pair=='2'&predframe$canopy=='thinned'] <- 'TAY'
predframe$site[predframe$pair=='1'&predframe$canopy=='intact'] <- 'LOL'
predframe$site[predframe$pair=='1'&predframe$canopy=='thinned'] <- 'UPL'

# Adjust the end point based upon the site
predframe <- predframe[!(predframe$site %in% c('CAI','TAY') & predframe$tsitr>37),]

# Remove guppy points prior to introduction
predframe <- predframe[!(predframe$sp == 'Guppy' & predframe$tsitr < 0),]

# Calculate pool:riffle ratios from data
raw.ratios <- read.csv(file.path('dryad_habitat_use_dataset.csv')) %>% 
  rename('pair' = 'stream_pair', 
         'canopy' = 'canopy_treatment', 
         'hab' = 'habitat', 
         'sp' = 'species', 
         'tsi.m' = 'time_since_introduction', 
         'len' = 'stream_length', 
         'ncap' = 'no_captured', 
         'capdens' = 'capture_density', 
         'log.cd' = 'log_capture_density') %>%
  pivot_wider(id_cols = c(site:reach, sp:tsi.m), 
              names_from = hab, values_from = capdens) %>% 
  mutate(ratio = P/R, 
         legend = factor(paste(sp, reach))) 



tiff('Fig3.tif', height = 5600, width = 8000, units = 'px', res = 800, compression = 'lzw')
ggplot(data=predframe, aes(x=tsitr, y=expect, color=legend, shape=legend)) + 
  geom_point(aes(x=tsi.m,y=ratio),data=raw.ratios, alpha = 0.25) + 
  geom_line(aes(linetype = legend)) + 
  geom_ribbon(aes(ymin = low, ymax=high, fill=legend), color = NA, alpha=.1, show.legend=c('fill'=T)) + 
  geom_hline(aes(yintercept=1), color='black',alpha=0.5) +
  scale_shape_manual('',values=c('Guppy I'=17,'Killifish C'=15,'Killifish I'=16), 
                     labels=c('guppy introduction','killifish control','killifish introduction')) + 
  scale_linetype_manual('',values=c('Killifish C'='solid','Killifish I'='dashed', 'Guppy I'='dotted'), 
                        labels=c('guppy introduction','killifish control','killifish introduction')) + 
  scale_discrete_manual(aesthetics=c('color','fill'), name='', 
                        values=c('Killifish C'='blue','Killifish I'='red', 'Guppy I'= 'peru'), 
                        labels=c('guppy introduction','killifish control','killifish introduction')) + 
  theme_minimal() +
  theme(panel.spacing=unit(2,'lines'), 
        strip.text=element_text(face='bold',size=16), 
        legend.text=element_text(size=13), 
        axis.text=element_text(size=13), 
        axis.title=element_text(size=20)) +
  xlab('Months since introduction') +
  ylab(expression(~frac(Pool,Riffle)~' Density of Captured Fish'))+
  facet_wrap(~site, scales='free_x') + 
  scale_y_continuous(trans='log10', breaks =c(0.1,1,10,100,1000), labels = c('0.1','1.0','10','100','1000')) + 
  scale_x_continuous(breaks = c(0,15,30,45))
dev.off()



