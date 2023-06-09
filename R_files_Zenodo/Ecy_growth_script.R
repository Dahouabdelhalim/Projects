# Load packages
library(tidyverse)
library(parallel)
library(lme4)
library(gganimate)

# Read data and rename columns for brevity
gro <- read.csv(file.path('dryad_growth_dataset.csv')) %>%
  rename('pair' = 'stream_pair', 
         'canopy' = 'canopy_treatment', 
         'hab' = 'habitat', 
         'tsi.c' = 'time_since_introduction_centered', 
         'tl' = 'total_length_initial', 
         'lagtl' = 'total_length_final', 
         'dt' = 'time_interval', 
         'rate' = 'growth_rate', 
         'lrate' = 'log_growth_rate') %>%
  mutate(pair = factor(pair),
         fish_id = factor(fish_id),
         rhab = factor(paste(reach, hab)),
         tsi.m = tsi.c + 18.32083,
         tlc = tl - mean(tl))


# Begin model fitting

# First, use a random intercept for fish_id
# This is for demonstration purposes. There's insufficient variation data to 
# fit the full random intercept model that also includes a temporal trend

# Define auxiliary functions that break growth into two phases
# juvenile
b.js <- function(x,m) ifelse(x<m,m-x,0)
# adult
b.as <- function(x,m) ifelse(x<m,0,x-m)

# Define model as a function of the breakpoint, m
gmod.re <- function(m){
  mod <- lmer(lrate~b.js(tlc,m)*site*reach*hab + b.as(tlc,m)*site*reach*hab + (1|fish_id), data=gro)
  -logLik(mod)
}
# Find the optimal breakpoint
growth.opt <- optim(par=0,lower=-25,upper=25, fn = gmod.re, hessian=T, method='Brent')

# Errors indicative that there's insufficient information 
# Store the estimate
m <- growth.opt$par

# Fit model using the breakpoint estimate
mod.re <- lmer(lrate~b.js(tlc,m)*site*reach*hab + b.as(tlc,m)*site*reach*hab + (1|fish_id), data=gro)
summary(mod.re)
# Note the extremely small variance with respect to individual
# Revert to a simple lm approach and fit a model with a temporal trend

# Define model as a function of the breakpoint, m
gmods <- function(m){
  mod <- lm(lrate~b.js(tlc,m)*site*reach*hab*tsi.c + b.as(tlc,m)*site*reach*hab*tsi.c, data=gro)
  -logLik(mod)
}
# Find the optimal breakpoint
growth.opt <- optim(par=0,lower=-25,upper=25, fn = gmods, hessian=T, method='Brent')
# Store the estimate
m <- growth.opt$par

# Estimates of uncertainty in breakpoint based upon 
# Hessian of optimization
m.sd <- sqrt(solve(growth.opt$hessian))
# and profile confidence intervals
root <- function(m,value){gmods(m)-value}
m.ci <- c(uniroot(root, lower=-30, upper=m, value=growth.opt$value+qchisq(.95,1))$root, # lower 95% confidence bound
         uniroot(root, lower=m, upper=30, value=growth.opt$value+qchisq(.95,1))$root)+mean(gro$tl) # upper 95% confidence bound

# Fit model using the breakpoint estimate
mod <- lm(lrate~b.js(tlc,m)*site*reach*hab*tsi.c+b.as(tlc,m)*site*reach*hab*tsi.c, data=gro)
# Have a look at the estimates
summary(mod)


# Animation of predictions over time
# Create a grid of points to make predictions
newdat <- expand.grid('site'=factor(c('CAI','LOL','TAY','UPL')),
                      'reach'=factor(c('C','I')),
                      'hab'=factor(c('P','R')),
                      'tl'=seq(10,100,.1), 
                      tsi.m=seq(0,36,1)) %>% 
  mutate('rhab'=factor(paste(reach,hab)), 
         tlc=(tl-mean(gro$tl)), 
         tsi.c = tsi.m-18.32083) %>% 
  filter(!(site == 'CAI' & reach == 'I' & hab == 'R' & tlc < m))
# Make predictions
pred <- predict(mod, newdata=newdat, interval='confidence',level=0.95)
# Adjust predictions for log transform
pred <- exp(pred)*(1+as.numeric(sigma(mod))^2/2)-3.5
colnames(pred) <- c('pred','lower','upper')
# Add predictions and confidence intervals to point grid
newdat <- cbind(newdat,pred)

# Set-up animation
animate_time  <-  newdat %>% 
  ggplot(aes(x=tl,y=pred,color=rhab,shape=rhab,linetype=rhab, group=rhab)) + 
  geom_line(aes(linetype=rhab)) + 
  scale_shape_manual('',values=c('C P'=16,'C R'=2,'I P'=16,'I R'=2), 
                     labels=c('control pool','control riffle','introduction pool','introduction riffle')) + 
  scale_linetype_manual('',values=c('C P'='solid','C R'='dashed','I P'='solid','I R'='dashed'), 
                        labels=c('control pool','control riffle','introduction pool','introduction riffle')) + 
  scale_fill_manual('',values=c('C P'='blue','C R'='blue','I P'='red','I R'='red'), 
                    labels=c('control pool','control riffle','introduction pool','introduction riffle')) + 
  scale_color_manual('',values=c('C P'='blue','C R'='blue','I P'='red','I R'='red'), 
                     labels=c('control pool','control riffle','introduction pool','introduction riffle')) + 
  theme_minimal() + 
  theme(panel.spacing=unit(2,'lines'),
        strip.text=element_text(face='bold',size=16), 
        legend.text=element_text(size=16), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=20)) + 
  xlab('Total Length (mm)') + 
  ylab(expression('Monthly Growth Rate ('~frac('mm','mo.')~')')) +
  facet_wrap(~site, scales='free_y') + 
  geom_point(data = transform(gro, tsi.m=NULL), aes(x=tl,y=rate, color=rhab, shape=rhab), alpha=.1) + 
  transition_time(tsi.m) + 
  labs(title = 'Time since introduction: {frame_time} months') 
# To create a gif...
# animate(animate_time, nframes = 42, fps = 2.5, width = 720, height = 480, res = 72, renderer = gifski_renderer(loop=T), end_pause = 5, dev='png', type='quartz')
# To create a movie file
animate(animate_time, nframes = 42, fps = 2.5, width = 720, height = 480, res = 72, renderer = av_renderer(), end_pause = 5, dev='png', type='quartz')
# Save animation in mpeg format
anim_save('KilliGrowth_lineanim.mpeg')

# Figure 5
# Pattern of growth rates at 25 mm, 40 mm and 55 mm TL
# Set-up grid for predictions
newdat = expand.grid('site'=factor(c('CAI','LOL','TAY','UPL')),
                     'reach'=factor(c('C','I')),
                     'hab'=factor(c('P','R')),
                     'tl'=c(25,40, 55), 
                     tsi=seq(4,36,.5)) %>% 
  mutate('rhab'=factor(paste(reach,hab)), 
         tlc=(tl-mean(gro$tl)), tsi.c = tsi-18.32083) %>% 
  filter(!(site == 'CAI' & reach == 'I' & hab == 'R' & tlc < m))
# Make predictions on grid
pred <- predict(mod, newdata=newdat, interval='confidence',level=0.95)
# Adjust predictions for log transform
pred <- exp(pred)*(1+as.numeric(sigma(mod))^2/2)-3.5
colnames(pred) <- c('pred','lower','upper')
# Add predictions to grid
newdat <- cbind(newdat,pred)

# Modify dataset for plotting data points in the representative size panels
gro.graph  <- gro %>% 
  mutate(tsi=NULL, 
         tl.class = ifelse(tl <= 32.5, 25, 
                           ifelse(tl > 32.5 & tl <= 47.5, 40,
                                  ifelse(tl > 47.5, 55, NA))), 
         tl.ptsize = 1/(1+abs(tl-tl.class)), 
         tl.lab = paste('TL =', tl.class,'mm')) %>% 
  filter(!is.na(tl.class))

tiff('Fig5.tif', height = 6400, width = 8000, units = 'px', res =800, compression = 'lzw')
newdat %>% 
  filter(!(site %in% c('LOL','UPL') & tsi<10), !(site %in% c('CAI','TAY') & tsi>24)) %>% 
  mutate(tl.lab = paste('TL =', tl, 'mm')) %>% 
  ggplot(aes(x=tsi,y=pred,color=rhab,shape=rhab,linetype=rhab, group=rhab)) +
  geom_line(aes(linetype=rhab)) +
  geom_ribbon(aes(ymin = lower, ymax=upper, fill=rhab), color = NA, linetype='solid', alpha=.1) + 
  scale_shape_manual('',values=c('C P'=16,'C R'=2,'I P'=16,'I R'=2), 
                     labels=c('control pool','control riffle','introduction pool','introduction riffle')) + 
  scale_linetype_manual('',values=c('C P'='solid','C R'='dashed','I P'='solid','I R'='dashed'), 
                        labels=c('control pool','control riffle','introduction pool','introduction riffle')) +
  scale_fill_manual('',values=c('C P'='blue','C R'='blue','I P'='red','I R'='red'), 
                    labels=c('control pool','control riffle','introduction pool','introduction riffle')) +
  scale_color_manual('',values=c('C P'='blue','C R'='blue','I P'='red','I R'='red'), 
                     labels=c('control pool','control riffle','introduction pool','introduction riffle')) +
  theme_minimal() +
  theme(panel.spacing=unit(2,'lines'), 
        strip.text=element_text(face='bold',size=10), 
        legend.text=element_text(size=13), 
        axis.text=element_text(size=13), 
        axis.title=element_text(size=20)) + 
  xlab('Time since introduction (months)') + 
  ylab(expression('Monthly Growth Rate ('~frac('mm','mo.')~')')) +
  facet_grid(site~tl.lab, scales='free_y', switch=NULL) + 
  geom_point(data=gro.graph, aes(x=tsi.m, y=rate, color=rhab, shape=rhab, size=tl.ptsize, alpha=tl.ptsize), 
             show.legend=c('color'=T,'shape'=T,'size'=F, 'alpha'=F)) +
  scale_alpha(range=c(.1, .3)) +
  scale_size(range=c(.5,3)) + 
  scale_y_continuous(breaks = seq(-5,15,5))
dev.off()

# Find coeffecients for Table S3
# Combination of factors at which to generate predictions to get coeffiencients
datmm.time <- expand.grid('site'=factor(c('CAI','LOL','TAY','UPL')),
                          'reach'=factor(c('C','I')),
                          'hab'=factor(c('P','R')),
                          'tlc'=0, # Mean size
                          'tsi.c'=c(0,1)) %>% # Time 0,1
  arrange(site,reach,hab)
# Extract beta estimates
beta <- coef(mod)
# Find design matrix
Xmat <- model.matrix(formula(mod)[-2], datmm.time)
# Remove intercept from every other row so that only time trend is present
Xmat[seq(2,nrow(Xmat),2),] <- Xmat[seq(2,nrow(Xmat),2),]-Xmat[seq(1,nrow(Xmat),2),]
# Compute predicted values
predinttime <- Xmat %*% (beta)
# compute vcov matrix
predinttime.vcv <- (Xmat) %*% vcov(mod) %*% t(Xmat)
# Extract the intercept predictions and vcov matrix
predint <- predinttime[seq(1,length(predinttime),2)]
predint.vcv <- predinttime.vcv[seq(1,length(predinttime),2),seq(1,length(predinttime),2)]
# Extract the time predicitons and vcov matrix
predtime <- predinttime[seq(2,length(predinttime),2)]
predtime.vcv <- predinttime.vcv[seq(2,length(predinttime),2),seq(2,length(predinttime),2)]

# Now assess size:time interactions
datmm.sizetime <- expand.grid('site'=factor(c('CAI','LOL','TAY','UPL')),
                              'reach'=factor(c('C','I')),
                              'hab'=factor(c('P','R')),
                              'tlc'=c(-1,0,1)+m, # +/- 1 unit from the breakpoint estimate
                              'tsi.c'=c(0,1)) %>% 
  arrange(site,reach,hab,tsi.c)
# Define model matrix on the new points
Xmat <- model.matrix(formula(mod)[-2], datmm.sizetime)
# Adjust model matrix to remove the estimate at the breakpoint from the points on either side
Xmat[seq(1,nrow(Xmat),6),] <- Xmat[seq(1,nrow(Xmat),6),]-Xmat[seq(2,nrow(Xmat),6),]
Xmat[seq(3,nrow(Xmat),6),] <- Xmat[seq(3,nrow(Xmat),6),]-Xmat[seq(2,nrow(Xmat),6),]
Xmat[seq(4,nrow(Xmat),6),] <- Xmat[seq(4,nrow(Xmat),6),]-Xmat[seq(5,nrow(Xmat),6),]
Xmat[seq(6,nrow(Xmat),6),] <- Xmat[seq(6,nrow(Xmat),6),]-Xmat[seq(5,nrow(Xmat),6),]

# Compute predictions and vcov matrix
predcoef <- Xmat %*% (beta)
predcoef.vcv <- (Xmat) %*% vcov(mod) %*% t(Xmat)

# Extract predictions and vcov matrices for the breakpoint, juvenile and adult phases at time 0 and 1
vcv.bpt0 <-  predcoef.vcv[seq(2,length(predcoef),6),seq(2,length(predcoef),6)]
vcv.bjst0 <- predcoef.vcv[seq(1,length(predcoef),6),seq(1,length(predcoef),6)]
vcv.bast0 <- predcoef.vcv[seq(3,length(predcoef),6),seq(3,length(predcoef),6)]
bhat.bpt0 <- predcoef[seq(2,length(predcoef),6)]
bhat.jst0 <- predcoef[seq(1,length(predcoef),6)]
bhat.ast0 <- predcoef[seq(3,length(predcoef),6)]
vcv.bpt1 <- predcoef.vcv[seq(5,length(predcoef),6),seq(5,length(predcoef),6)]
vcv.bjst1 <- predcoef.vcv[seq(4,length(predcoef),6),seq(4,length(predcoef),6)]
vcv.bast1 <- predcoef.vcv[seq(6,length(predcoef),6),seq(6,length(predcoef),6)]
bhat.bpt1 <- predcoef[seq(5,length(predcoef),6)]
bhat.jst1 <- predcoef[seq(4,length(predcoef),6)]
bhat.ast1 <- predcoef[seq(6,length(predcoef),6)]

# Define an operator to find the size:time interactions
sizetime.oper <- cbind(diag(-1,nrow=nrow(Xmat),ncol=nrow(Xmat))[,seq(1,nrow(Xmat),6)] + 
                        diag(-1,nrow=nrow(Xmat),ncol=nrow(Xmat))[,seq(4,nrow(Xmat),6)],
                      diag(-1,nrow=nrow(Xmat),ncol=nrow(Xmat))[,seq(3,nrow(Xmat),6)] + 
                        diag(-1,nrow=nrow(Xmat),ncol=nrow(Xmat))[,seq(6,nrow(Xmat),6)])

# Find the size:time interactions and the vcov matrix
bhat.st <- t(sizetime.oper) %*% (predcoef)
vcv.st <- t(sizetime.oper) %*% predcoef.vcv %*% sizetime.oper
bhat.stj <- c(bhat.st)[1:16]
vcv.stj <- vcv.st[1:16, 1:16]
bhat.sta <- c(bhat.st)[17:32]
vcv.sta <- vcv.st[17:32, 17:32]
 # Put all parameter estimates in a table
all.par <- bind_rows(
  expand.grid('parameter' = c('y-intercept', 'time'), 
              stage='-', '
              stream'=factor(c('CAI','LOL','TAY','UPL')),
              'reach'=factor(c('C','I')),
              'hab'=factor(c('P','R'))), 
  expand.grid('parameter' = c('size', 'time:size'), 
              stage= c('Juvenile','Adult'),  
              'stream'=factor(c('CAI','LOL','TAY','UPL')),
              'reach'=factor(c('C','I')),
              'hab'=factor(c('P','R')))
  ) %>% 
  mutate(parameter=factor(parameter, levels=c('breakpoint','y-intercept','time','size','time:size')), 
         stage=factor(stage, levels=c('-','Juvenile','Adult'))) %>% 
  arrange(parameter, stage, stream, reach, hab,stage)  %>% 
  mutate(bhat = c(predint,predtime, bhat.jst0, bhat.ast0, bhat.stj, bhat.sta), 
         v = c(diag(predint.vcv), diag(predtime.vcv), diag(vcv.bjst0), diag(vcv.bast0), diag(vcv.stj), diag(vcv.sta))) %>% 
  mutate(se = sqrt(v)) %>% 
  bind_rows(
    data.frame(parameter='breakpoint',
               stage='-',
               stream='-',
               reach='-',
               hab='-', 
               bhat=m+mean(gro$tl), 
               se= m.sd)) %>% 
  mutate(bhat = round(bhat,3), 
         se = round(se,4), 
         parameter=factor(parameter, levels=c('breakpoint','y-intercept','time','size','time:size')), 
         stage=factor(stage, levels=c('-','Juvenile','Adult'))) %>% 
  arrange(parameter, stage, stream, reach, hab)

# contrasts for 25, 40, 55 mm individuals @ mean time
# Define a grid for contrasts
datmm.size <- expand.grid('site'=factor(c('CAI','LOL','TAY','UPL')),
                         'reach'=factor(c('C','I')),
                         'hab'=factor(c('P','R')),
                         'tl'=c(25,40,55), 
                         'tsi.c'=c(0,1)) %>% 
  mutate(tlc = tl-mean(gro$tl)) %>% 
  arrange(site,reach,hab,tlc)
# Extract coefficients
beta  <- coef(mod)
# Define design matrix
Xmat <- model.matrix(formula(mod)[-2], datmm.size)
# Remove intercept from every other row so that only time trend is present
Xmat[seq(2,nrow(Xmat),2),] <- Xmat[seq(2,nrow(Xmat),2),]-Xmat[seq(1,nrow(Xmat),2),]
# Calculate predicted values
predintsize <- Xmat %*% (beta)
# and variance covariance matrix
predintsize.vcv <- (Xmat) %*% vcov(mod) %*% t(Xmat)
# Extract the intercepts and vcov
predint.s <- predintsize[seq(1,length(predintsize),2)]
predint.s.vcv <- predintsize.vcv[seq(1,length(predintsize),2),seq(1,length(predintsize),2)]
# Extract the temporal trend and vcov
predtime.s <- predintsize[seq(2,length(predintsize),2)]
predtime.s.vcv <- predintsize.vcv[seq(2,length(predintsize),2),seq(2,length(predintsize),2)]

# Breakup intercepts by size class: sm = 25 mm, md = 40 mm, lg = 55 mm
predint.sm <- predint.s[seq(1,length(predint.s),3)]
predint.sm.vcv <- predint.s.vcv[seq(1,length(predint.s),3),seq(1,length(predint.s),3)]
predint.md <- predint.s[seq(2,length(predint.s),3)]
predint.md.vcv <- predint.s.vcv[seq(2,length(predint.s),3),seq(2,length(predint.s),3)]
predint.lg <- predint.s[seq(3,length(predint.s),3)]
predint.lg.vcv <- predint.s.vcv[seq(3,length(predint.s),3),seq(3,length(predint.s),3)]
# Repeat for temporal trends
predtime.sm <- predtime.s[seq(1,length(predtime.s),3)]
predtime.sm.vcv <- predtime.s.vcv[seq(1,length(predtime.s),3),seq(1,length(predtime.s),3)]
predtime.md <- predtime.s[seq(2,length(predtime.s),3)]
predtime.md.vcv <- predtime.s.vcv[seq(2,length(predtime.s),3),seq(2,length(predtime.s),3)]
predtime.lg <- predtime.s[seq(3,length(predtime.s),3)]
predtime.lg.vcv <- predtime.s.vcv[seq(3,length(predtime.s),3),seq(3,length(predtime.s),3)]

# Define contrast operator matrices
# Differences between reaches or reaches x habitats within each stream
oper.win <- cbind('CAIr'=c(1,1,-1,-1,rep(0,12)), 
                 'CAIrh'=c(1,-1,-1,1,rep(0,12)), 
                 'LOLr'=c(rep(0,4),1,1,-1,-1,rep(0,8)), 
                 'LOLrh'=c(rep(0,4),1,-1,-1,1,rep(0,8)), 
                 'TAYr'=c(rep(0,8),1,1,-1,-1,rep(0,4)), 
                 'TAYrh'=c(rep(0,8),1,-1,-1,1,rep(0,4)), 
                 'UPLr'=c(rep(0,12),1,1,-1,-1), 
                 'UPLrh'=c(rep(0,12),1,-1,-1,1))
# Differences between reach or reaches x habitats across streams
oper.avg <- cbind('reach'=rep(c(1,1,-1,-1),4),
                 'hab'=rep(c(1,-1,1,-1),4),
                 'rhab'=rep(c(1,-1,-1,1),4), 
                 'reach_noC'=c(rep(0,4),rep(c(1,1,-1,-1),3)), # Ignore caigual reach with insufficient data
                 'hab_noC'=c(rep(0,4),rep(c(1,-1,1,-1),3)), # Ignore caigual reach with insufficient data
                 'rhab_noC'=c(rep(0,4), rep(c(1,-1,-1,1),3))) # Ignore caigual reach with insufficient data

# Differences between reaches-habitat combinations
oper.rhdiff <- cbind('CP-CR'=rep(c(1,-1,0,0),4), 
                    'CP-IP'=rep(c(1,0,-1,0),4),
                    'CP-IR'=rep(c(1,0,0,-1)),
                    'CR-IP'=rep(c(0,1,-1,0),4),
                    'CR-IR'=rep(c(0,1,0,-1),4),
                    'IP-IR'=rep(c(0,0,1,-1),4))
# Excluding problematic Caigual
oper.rhdiff_noC <- cbind('CP-CR_noC'=c(rep(0,4),rep(c(1,-1,0,0),3)), 
                        'CP-IP_noC'=c(rep(0,4),rep(c(1,0,-1,0),3)),
                        'CP-IR_noC'=c(rep(0,4),rep(c(1,0,0,-1),3)),
                        'CR-IP_noC'=c(rep(0,4),rep(c(0,1,-1,0),3)),
                        'CR-IR_noC'=c(rep(0,4),rep(c(0,1,0,-1),3)),
                        'IP-IR_noC'=c(rep(0,4),rep(c(0,0,1,-1),3)))

# Combine contrast operators into on big matrix
oper <- cbind(oper.avg, oper.rhdiff, oper.rhdiff_noC, oper.win)

# Contrasts at the time series midpoint for the three sizes
contrint.sm <- t(oper) %*% predint.sm
contrint.sm.vcv <- t(oper) %*% predint.sm.vcv %*% oper
contrint.md <- t(oper) %*% predint.md
contrint.md.vcv <- t(oper) %*% predint.md.vcv %*% oper
contrint.lg <- t(oper) %*% predint.lg
contrint.lg.vcv <- t(oper) %*% predint.lg.vcv %*% oper

# Contrasts of the temporal trend for the three sizes
contrtime.sm <- t(oper) %*% predtime.sm
contrtime.sm.vcv <- t(oper) %*% predtime.sm.vcv %*% oper
contrtime.md <- t(oper) %*% predtime.md
contrtime.md.vcv <- t(oper) %*% predtime.md.vcv %*% oper
contrtime.lg <- t(oper) %*% predtime.lg
contrtime.lg.vcv <- t(oper) %*% predtime.lg.vcv %*% oper

# Assemble table of contrasts with significance
tab.out <- data.frame(term = rep(c('intercept','time'), each=3*ncol(oper)),
                      size = rep(rep(unique(datmm.size$tl), each = ncol(oper)),2), 
                      'contrast' = rep(colnames(oper),6), 
                      est = c(contrint.sm, contrint.md, contrint.lg, contrtime.sm, contrtime.md, contrtime.lg), 
                      se = c(sqrt(diag(contrint.sm.vcv)), sqrt(diag(contrint.md.vcv)), sqrt(diag(contrint.lg.vcv)), 
                             sqrt(diag(contrtime.sm.vcv)), sqrt(diag(contrtime.md.vcv)), sqrt(diag(contrtime.lg.vcv)))) %>% 
  mutate(t = est/se, ddf = df.residual(mod), p = 2*(1-pt(abs(t),df=ddf)))

# Marginalize the estimates across the streams
oper.rh  <- cbind('C P'=rep(c(1,0,0,0),4)/4,
                  'C R'=rep(c(0,1,0,0),4)/4, 
                  'I P'=rep(c(0,0,1,0),4)/4,
                  'I R'=rep(c(0,0,0,1),4)/4)
# Exclude problematic reach-habitat from Caigual
oper.rh_noC  <- cbind('C P_noC'=c(rep(0,4),rep(c(1,0,0,0),3)/3),
                      'C R_noC'=c(rep(0,4),rep(c(0,1,0,0),3)/3), 
                      'I P_noC'=c(rep(0,4),rep(c(0,0,1,0),3)/3),
                      'I R_noC'=c(rep(0,4),rep(c(0,0,0,1),3)/3))

# Compute marginal intercepts and uncertainty for three sizes
int.marg.sm <- t(predint.sm) %*% oper.rh_noC
int.marg.sm.vcv <- t(oper.rh_noC) %*% predint.sm.vcv %*% oper.rh_noC
int.marg.md <- t(predint.md) %*% oper.rh
int.marg.md.vcv <- t(oper.rh) %*% predint.md.vcv %*% oper.rh
int.marg.lg <- t(predint.lg) %*% oper.rh
int.marg.lg.vcv <- t(oper.rh) %*% predint.lg.vcv %*% oper.rh

# Compute marginal temporal trends and uncertainty for three sizes
time.marg.sm <- t(predtime.sm) %*% oper.rh_noC
time.marg.sm.vcv <- t(oper.rh_noC) %*% predtime.sm.vcv %*% oper.rh_noC
time.marg.md <- t(predtime.md) %*% oper.rh
time.marg.md.vcv <- t(oper.rh) %*% predtime.md.vcv %*% oper.rh
time.marg.lg <- t(predtime.lg) %*% oper.rh
time.marg.lg.vcv <- t(oper.rh) %*% predtime.lg.vcv %*% oper.rh

# Assemble marginalized estimates and uncertainties
all.marg.int <- c(int.marg.sm, int.marg.md, int.marg.lg)
all.marg.time <- c(time.marg.sm, time.marg.md, time.marg.lg)
all.marg.int.se <- c(sqrt(diag(int.marg.sm.vcv)), sqrt(diag(int.marg.md.vcv)), sqrt(diag(int.marg.lg.vcv)))
all.marg.time.se <- c(sqrt(diag(time.marg.sm.vcv)), sqrt(diag(time.marg.md.vcv)), sqrt(diag(time.marg.lg.vcv)))

data.frame(term = rep(c('intercept','time'), each = 12), 
           size = rep(rep(c(25,40,55), each=4),2), 
           reach = rep(rep(c('C','I'),each=2), 6), 
           hab = rep(rep(c('P','R'),2), 6), 
           est = c(exp(all.marg.int)*(1+as.numeric(sigma(mod))^2/2)-3.5, # Adjust estimate for log transformation
                   100*(exp(all.marg.time)-1)), # % change over time
           lower=c(exp(all.marg.int+qt(.025,df=df.residual(mod))*all.marg.int.se)*(1+as.numeric(sigma(mod))^2/2)-3.5, 
                   100*(exp(all.marg.time+qt(.025,df=df.residual(mod))*all.marg.time.se)-1)), 
           upper=c(exp(all.marg.int+qt(.975,df=df.residual(mod))*all.marg.int.se)*(1+as.numeric(sigma(mod))^2/2)-3.5, 
                   100*(exp(all.marg.time+qt(.975,df=df.residual(mod))*all.marg.time.se)-1))) %>% 
  mutate(est = round(est,2), 
         lower = round(lower, 3), 
         upper = round(upper, 3))

