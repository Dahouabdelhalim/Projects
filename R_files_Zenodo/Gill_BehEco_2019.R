setwd('X:/Dropbox (Personal)/crypticfolder_encrypted/Doktorarbeit/R/JD_workspace/')
setwd('C:/Users/gill/Dropbox (Personal)/JaapLisa/Review3/R')

require(ggplot2)
require(doBy)
require(lme4)
require(arm)

## Fig. 1, Model 1, Table 1: Resident females spent most time alone inside their nest-box post-fertile ---------------------------------------------
fertile <- read.csv('1_fem_fertile.csv')
fertile$date <- as.Date(fertile$date)
dursumm <- read.csv('2_2017_cam_all_summary_timeinnest.csv')
dursumm$date <- as.Date(dursumm$date)
dursumm$dur <- as.numeric(as.character(dursumm$dur))
dursumm$Female <- fertile$female_henderson[match(paste(dursumm$cam, dursumm$date), paste(fertile$cam, fertile$date))]
dursumm$Female <- factor(dursumm$Female, levels =c('pre-fertile', 'fertile', 'post-fertile'))
dursumm$tabsimplestage <- factor(dursumm$tabsimplestage,levels=c('pre-laying','laying','incubating'))
dursumm$dur_int <- as.integer(as.character(dursumm$dur))
dursumm$hours <- dursumm$dur/60
dursumm1 <- dursumm[!(dursumm$date %in% unique(dursumm$date)[1] & dursumm$cam %in% 'cam09') ,] #remove NAs (missing observations)


plotmod <- function(x)
{ par(mfrow=c(2,2))
  plot(fitted(mod), resid(mod))
  abline(h=0, lty=2)
  qqnorm(resid(mod), main='normal qq-plot, residuals')
  qqline(resid(mod))
  scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))
  qqnorm(unlist(ranef(mod)), main='normal qq-plot, random effects')
  qqline(unlist(ranef(mod)))
  par(mfrow=c(1,1))}

#Model 1:
mod <- lmer(log(dur+1) ~ whos_in * Female + (1|Nest.box), data = dursumm1)
modnull<- lmer(log(dur+1) ~ whos_in + Female + (1|Nest.box), data = dursumm1)
modbasic <- lmer(log(dur+1) ~ 1 + (1|Nest.box) , data = dursumm1)
anova(mod, modnull, modbasic)
# refitting model(s) with ML (instead of REML)
# Data: dursumm1
# Models:
#   modbasic: log(dur + 1) ~ 1 + (1 | Nest.box)
# modnull: log(dur + 1) ~ whos_in + Female + (1 | Nest.box)
# mod: log(dur + 1) ~ whos_in * Female + (1 | Nest.box)
# Df    AIC    BIC   logLik deviance  Chisq Chi Df Pr(>Chisq)    
# modbasic  3 2110.2 2123.0 -1052.12   2104.2                             
# modnull   7 1890.5 1920.3  -938.23   1876.5 227.77      4  < 2.2e-16 ***
#   mod      11 1448.5 1495.3  -713.22   1426.5 450.01      4  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

summary(mod)
plot(resid(mod))
plotmod(mod)
fixef(mod)
ranef(mod)
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
apply(bsim@fixef,2,quantile,prob=c(0.025,0.5,0.975))

#credible intervals:
newdat <- expand.grid(Female=levels(dursumm$Female), whos_in =levels(dursumm$whos_in))
Xmat <- model.matrix(~whos_in*Female, data=newdat)
dim(Xmat)
predmat <- matrix(ncol=nsim, nrow=nrow (newdat))
newdat$pred <- Xmat%*% fixef(mod)
for(i in 1:nsim) predmat[,i] <-(Xmat%*%bsim@fixef[i,])
newdat$lower <- apply(predmat, 1, quantile, prob=0.025)
newdat$upper <- apply(predmat, 1, quantile, prob=0.975)

#Table 1:
newdat

#Fig. 1:
ggplot()+
  geom_jitter(data= dursumm1, aes(x=Female, y=dur+1, shape =Nest.box), height=0, width =0.15, size =1.8, colour='gray60')+
  facet_wrap(~whos_in)+
  scale_shape_manual(name  ="Nest-box", breaks=levels(dursumm$Nest.box),labels=
                       levels(dursumm$Nest.box), values =c(0:4,6,8,9))+
  theme_bw()+  #theme(axis.text=element_text(size=14,colour='black'),axis.title=element_text(size=16,face="bold"))+
  geom_point(data=newdat, aes(x= Female, y= exp(pred),colour= whos_in), size = 3)+
  geom_segment(data=newdat, aes(x = Female, y = exp(newdat$lower), xend = Female, 
                                yend = exp(newdat$upper), colour = whos_in), size =1.5)+
  scale_y_log10(name=('Time spent in nest [minutes]'))

# ## R2 following Nakagawa and Schielzeth 2013:
RS <- function(x){
  Fixed <- fixef(mod)[2] * getME(mod,"X")[, 2] + fixef(mod)[3] * getME(mod,"X")[, 3] +
    fixef(mod)[4] * getME(mod,"X")[, 4] + fixef(mod)[5] * getME(mod,"X")[, 5]+
    fixef(mod)[6] * getME(mod,"X")[, 6] + fixef(mod)[7] * getME(mod,"X")[, 7]+
    fixef(mod)[8] * getME(mod,"X")[, 8] + fixef(mod)[9] * getME(mod,"X")[, 9]
  VarF <- var(Fixed)
 paste("marginal: ", VarF/(VarF + VarCorr(mod)$Nest.box[1] + attr(VarCorr(mod), "sc")^2),"; ",
 "conditional: ",(VarF + VarCorr(mod)$Nest.box[1])/(VarF + VarCorr(mod)$Nest.box[1] +
                                                    (attr(VarCorr(mod), "sc")^2)), sep = "")
 }
 RS(mod) # =  "marginal: 0.72193172482612; conditional: 0.721931724826121"


 
## Fig. 2, Model 2, Table 2: Within pairs, full copulations were increasingly replaced by copulation attempts with progressing breeding ---------------------------------------------
 fertile <- read.csv('1_fem_fertile.csv')
 fertile$date <- as.Date(fertile$date)
 fertile$nr_EPevents <- as.numeric(as.character(fertile$nr_EPevents))
 fertile$nr_IPevents <- as.numeric(as.character(fertile$nr_IPevents))
 
 cops <- read.csv('3_2017_copulations_datetime.csv')
 cops$date <- as.Date(cops$date) 
 cops$female_henderson <- fertile$female_henderson[(match(paste(cops$date, cops$cam),
                                                          paste(fertile$date,fertile$cam)))]
 copsummary <- data.frame(table(cops$Nest.box, cops$IEP,  cops$att.cop, cops$female_henderson))
 names(copsummary) <- c("Nest.box","IEP","att.cop","Female","nr_events")
 copsummary$Female <- factor(copsummary$Female, levels =c('pre-fertile', 'fertile', 'post-fertile'))
 IP <-   copsummary[copsummary$IEP =='IP_IP',]
 
 #proportion of copulations out of all events:
 IPatt <- IP[IP$att.cop =='attempts',]
 IPcop <- IP[IP$att.cop =='copulations',]
 IPatt$nr_atts <- IPatt$nr_events
 IPatt$nr_cops <- IPcop$nr_events[match(paste(IPatt$Nest.box, IPatt$Female),paste(IPcop$Nest.box, IPcop$Female))]
 IPatt$nr_events <- IPatt$nr_atts+ IPatt$nr_cops
 IPatt$p <- IPatt$nr_cops/(IPatt$nr_atts+IPatt$nr_cops)
 
 plotmod <- function(x)
 { par(mfrow=c(2,2))
   plot(fitted(mod), resid(mod))
   abline(h=0, lty=2)
   qqnorm(resid(mod), main='normal qq-plot, residuals')
   qqline(resid(mod))
   scatter.smooth(fitted(mod), sqrt(abs(resid(mod))))
   qqnorm(unlist(ranef(mod)), main='normal qq-plot, random effects')
   qqline(unlist(ranef(mod)))
   par(mfrow=c(1,1))}
 
 
 ###binomial:
 mod <- glmer(cbind(IPatt$nr_cops, IPatt$nr_atts) ~ Female + (1|Nest.box) , data= IPatt, family =binomial)
 modnull <- glmer(cbind(IPatt$nr_cops, IPatt$nr_atts) ~1+ (1|Nest.box) , data= IPatt, family =binomial)
 anova(modnull, mod, test ='Chisq')
 
 summary(mod)
 anova(mod)
 plot(density(resid(mod)))
 plotmod(mod)
 
 # ## R2 following Nakagawa and Schielzeth 2013:
 RS <- function(x){
   Fixed <- fixef(mod)[2] * getME(mod,"X")[, 2] + fixef(mod)[3] * getME(mod,"X")[, 3]
   VarF <- var(Fixed)
   paste("marginal: ", VarF/(VarF + VarCorr(mod)$Nest.box[1] + attr(VarCorr(mod), "sc")^2),"; ",
         "conditional: ",(VarF + VarCorr(mod)$Nest.box[1])/(VarF + VarCorr(mod)$Nest.box[1] +
                                                              (attr(VarCorr(mod), "sc")^2)), sep = "")
 }
 RS(mod) #  "marginal: 0.233897876050841; conditional: 0.734229141847923"
 
 
 nsim <- 2000
 bsim <- sim(mod, n.sim=nsim)
 apply(bsim@fixef,2,quantile,prob=c(0.025,0.975))
 # (Intercept) Femalefertile Femalepost-fertile
 # 2.5%   -0.4916125   -1.41822482          -2.943497
 # 97.5%   1.9231908    0.07762834          -1.399099
 
 
 newdat <- expand.grid(Female=levels(IPatt$Female))
 Xmat <- model.matrix(~Female, data=newdat)
 dim(Xmat)
 predmat <- matrix(ncol=nsim, nrow=nrow (newdat))
 
 for(i in 1:nsim) predmat[,i] <- plogis(Xmat %*% bsim@fixef [i,])
 newdat$lower <- apply(predmat, 1, quantile, prob=0.025)
 newdat$upper <- (apply(predmat, 1, quantile, prob=0.975))
 newdat$pred <- plogis(Xmat %*% fixef(mod))
 
 #Table 2:
 newdat 
 # Female      lower     upper      pred
 # pre-fertile 0.38850976 0.8840796 0.6765385
 #     fertile 0.28249401 0.7378330 0.5113340
 #post-fertile 0.07647847 0.3901926 0.1881760
 
 #Fig. 2:
 ggplot()+
   geom_point(data= IPatt, aes(x=Female, y=p, shape =Nest.box,colour = Nest.box) ,size =3)+
   scale_shape_manual(name  ="Nest-box", breaks=levels(IPatt$Nest.box),labels=
                        levels(IPatt$Nest.box), values =c(0:4,6,8,9))+
   geom_line(data=IPatt, aes(x=Female, y=p, group = Nest.box, colour =Nest.box))+
   
   theme_bw()+  scale_x_discrete(name=('Female'))+scale_y_continuous(name='IP cop/(IP cop + IP att)')+
   geom_point(aes(x= c(1:3)+.1, y= c((newdat$pred))), size = 3)+
   geom_segment(aes(x = c(1:3)+.1, y = c((newdat$lower)), xend = c(1:3)+.1,
                    yend = c((newdat$upper))), size =1.2)
 
 
 
 

## Fig. 3, correlations: Negative association between a pair's chick fledging success and the occurrences of extra-pair sexual behaviour ---------------------------------------------
nests2017 <- read.csv('4_workfile_2017_nests.csv')
fertile <- read.csv('1_fem_fertile.csv')
fertile$date <- as.Date(fertile$date)
fertile$nr_EPevents <- as.numeric(as.character(fertile$nr_EPevents))
fertile$nr_IPevents <- as.numeric(as.character(fertile$nr_IPevents))

nestscams <- nests2017[nests2017$Nest.box %in% c(levels(nests2017$Nest.box)),]
nestscams$p <- nestscams$max_fledged/nestscams$max_chicks
nestscams$survived <- nestscams$num_fledged
nestscams$died <- nestscams$max_chicks - nestscams$survived


#Check: number of EP events in relation to clutch size: 
cor(nestscams$nr_EPevents, nestscams$max_eggs) #0.02099755
cor.test(nestscams$nr_EPevents, nestscams$max_eggs)
# Pearson's product-moment correlation
# 
# data:  nestscams$nr_EPevents and nestscams$max_eggs
# t = 0.051445, df = 6, p-value = 0.9606
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.6939432  0.7150897
# sample estimates:
# cor 
# 0.02099755 

#Fig.  (not shown in paper)
ggplot()+ geom_point(data = nestscams ,
                     aes(x=nr_EPevents, y= max_eggs, shape =Nest.box),size=3, stroke =1.1)+#, width =0, height = 0.4)+
  scale_shape_manual(name  ="Nest-box", breaks=levels(nestscams$Nest.box),labels=
                       levels(nestscams$Nest.box), values =c(0:4,6,8,9))+
  scale_colour_manual(values =c('black','blue'))+
  theme_bw()+
  scale_y_continuous(name=('Clutch size'))+ scale_x_continuous(name=('EP events'))+
geom_abline(intercept= 4.843373, slope = 0.006487)



# Correlation: number of EP events in relation to chick survival probability (per nest)
cor(nestscams$nr_EPevents , nestscams$p) # -0.7693077
cor.test(nestscams$nr_EPevents , nestscams$p)
# Pearson's product-moment correlation
# 
# data:  nestscams$nr_EPevents and nestscams$p
# t = -2.9496, df = 6, p-value = 0.02563
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.9558205 -0.1411578
# sample estimates:
# cor 
# -0.7693077 

# Fig. 3:
ggplot()+ geom_point(data = nestscams ,
                     aes(x=nr_EPevents, y= p, shape =Nest.box),size=3, stroke =1.1)+#, width =0, height = 0.4)+
  scale_shape_manual(name  ="Nest-box", breaks=levels(nestscams$Nest.box),labels=
                       levels(nestscams$Nest.box), values =c(0:4,6,8,9))+
  theme_bw(base_size = 18)+
  scale_y_continuous(name=('Proportion of chicks fledged'))+ scale_x_continuous(name=('EP events'))+
  geom_abline(intercept= 0.41807, slope = -0.03832) +
  annotate("text", x=9, y=0.65, label= "Pearson's r: -0.7693077", size = 6)
  


## Fig.S1: Nest-box occupancy and extra-pair (EP) sexual behaviour over breeding stages ---------------------------------------------

#time spent:
dur <- read.csv('5_workfile_all_durations.csv'); 
dur$date <- as.Date(dur$date, format ='%Y-%m-%d');
alone <- read.csv('6_2017_fem_alone_summary.csv')
alone$date <- as.Date(alone$date, format ='%Y-%m-%d')
alone$pc_alone <- (alone$time_fem_alone*60)/alone$timeinnest

#breeding stages:
tabsimplestage <- read.csv('7_workfile_breedingstages.csv')
tabsimplestage$date <- as.Date(tabsimplestage$date, format ='%Y-%m-%d')
tabsimplestage <- tabsimplestage[tabsimplestage$sum >0,]
tabsimplestage$simplestage <- factor(tabsimplestage$simplestage,levels=c('pre-laying','laying','incubating'))

#chick information:
nn1 <- read.csv('8_summary_nestchecks_filmedcams.csv')
nn1$date_hatched <- as.Date(nn1$date_hatched);
nn1$date_laid <- as.Date(nn1$date_laid)

#copulations:
cops <- read.csv('9_2017_copulations_datetime.csv'); 

cops$date <- as.Date(cops$date, format ='%Y-%m-%d')
cops <- cops[cops$IEP %in% c('EP_EP','EP_IP', 'IP_IP'),];cops$IEP <- factor(cops$IEP, levels =c('IP_IP','EP_IP', 'EP_EP'))
cops <- cops[cops$att.cop %in% c('attempts','copulations'),]

#incoming EPCs:
cops_EP_in <- cops[cops$IEP %in% 'EP_IP',]
cops_EP_in$simplestage <- tabsimplestage[tabsimplestage$sum >0,]$simplestage[match(paste(cops_EP_in$date, cops_EP_in$camera),
                                                                                   paste(tabsimplestage[tabsimplestage$sum >0,]$date, tabsimplestage[tabsimplestage$sum >0,]$cam))]
cops_EP_in$pc_spent <- dur$pc_spent[match(paste(cops_EP_in$date, cops_EP_in$cam), paste(dur$date, dur$cam))]

#outgoing EPCs:
residentmales <- read.csv('10_names_residentmales.csv')
cops_EP_out <- cops_EP_in[cops_EP_in$ring_male %in% residentmales$ring_male,]
cops_EP_out$male_from_cam <- residentmales$cam[match(cops_EP_out$ring_male, residentmales$ring_male)]
cops_EP_out$cam <- cops_EP_out$male_from_cam
cops_EP_out$simplestage <- tabsimplestage[tabsimplestage$sum >0,]$simplestage[match(paste(cops_EP_out$date, cops_EP_out$camera),
                                                                                    paste(tabsimplestage[tabsimplestage$sum >0,]$date, tabsimplestage[tabsimplestage$sum >0,]$cam))]

# Fig. S1:
mypalette = c("lightgrey","#FFFF99", "#F0E442", "#999999")

ggplot()+ ggtitle("Residents' time spent in nest, breeding stages, \\n hatching success and EP occurrences")+ 
  theme_bw() + facet_wrap(~Nest.box,ncol=2)+ 
  # breeding stages:
  geom_bar(data= tabsimplestage,aes(x=date, y=sum,fill=simplestage), stat='identity')+
  scale_fill_manual(values = mypalette , name="Breeding stage")+
  scale_y_continuous(name='Proportion of day spent in nest')+
  
  # time in nest for males (turquoise) and females (red):
  geom_line(data= dur[dur$IP=='IP' ,], aes(x=date, y=pc_spent,group = ring, colour=sex),size=1.1)  +  
  
  # eggs hatched:
  geom_jitter(data=nn1[nn1$hatched %in% c('yes'),], aes(y= 0.051, x=date_laid, group=egg.ID), size = 1.5, width = 0.05, height = 0, 
              shape =16, stroke =1.5)+
  geom_jitter(data=nn1[nn1$hatched %in% c('no'),], aes(y= 0.04, x=date_laid, group=egg.ID), size = 2, width = 0.05, height = 0, 
              shape =1)+
  
  # non-resident males intruding in each camera:
  geom_jitter(data= cops_EP_in, aes(x=date, y=0.85, group = ring_male), colour='black', size=1.8, shape = 25,width = 0.0, height = 0.1, 
              stroke =1.25)+
  #resident males being filmed by other cameras:  
  geom_jitter(data= cops_EP_out, aes(x=date, y=0.5, group= ring_male), colour='blue', size=1.8, shape = 24,width = 0.0, height = 0.1,
              stroke =1.25)




## Figs. S2 and S3:  ---------------------------------------------
fertile <- read.csv('1_fem_fertile.csv')
fertile$date <- as.Date(fertile$date)
fertile$nr_EPevents <- as.numeric(as.character(fertile$nr_EPevents))
fertile$nr_IPevents <- as.numeric(as.character(fertile$nr_IPevents))

cops <- read.csv('9_2017_copulations_datetime.csv')
cops$date <- as.Date(cops$date, format ='%Y-%m-%d')
cops <- cops[cops$IEP %in% c('EP_EP','EP_IP', 'IP_IP'),];
cops$IEP <- factor(cops$IEP, levels =c('IP_IP','EP_IP', 'EP_EP'))
cops <- cops[cops$att.cop %in% c('attempts','copulations'),]
cops$Female <- fertile$female_henderson[(match(paste(cops$date, cops$Nest.box),
                                                   paste(fertile$date,fertile$Nest.box)))]
copsummary <- data.frame(table(cops$Nest.box, cops$IEP,  cops$att.cop, cops$Female))
names(copsummary) <- c("Nest.box","IEP","att.cop","Female","nr_events")
copsummary$Female <- factor(copsummary$Female, levels =c('pre-fertile', 'fertile', 'post-fertile'))

EP <-   copsummary[copsummary$IEP =='EP_IP',]
IP <-   copsummary[copsummary$IEP =='IP_IP',]

tabsimplestage <- read.csv('7_workfile_breedingstages.csv')
tabsimplestage$date <- as.Date(tabsimplestage$date, format ='%Y-%m-%d')
tabsimplestage <- tabsimplestage[tabsimplestage$sum >0,]
tabsimplestage$simplestage <- factor(tabsimplestage$simplestage,levels=c('pre-laying','laying','incubating'))


dur <- read.csv('5_workfile_all_durations.csv');
dur$date <- as.Date(dur$date, format ='%Y-%m-%d');
#exact times:
 cops$timeMM_dec <- cops$timeMM/60
 cops$time_dec <- as.numeric(paste(cops$timeHH+cops$timeMM_dec))
#first and last entries:
min <- read.csv('11_entries_exits_summary.csv')
min$date <- as.Date(min$date)

#Fig.S2: Residents' (F, M) 
mypalette = c("lightgrey","#FFFF99", "#F0E442", "#999999")

ggplot()+ 
  geom_segment(data= tabsimplestage[tabsimplestage$sum >0,], 
               aes(x=date, xend=date, y=5, yend=sum*21,colour=simplestage), size =4.5)+
  theme_bw() + facet_wrap(~Nest.box,ncol=2)+ 
  scale_color_manual(values = mypalette , name="Breeding stage")+
  scale_x_date(name ='Date', date_breaks = '4 days', date_labels = "%d/%m")+
  scale_y_continuous(name='Time of day')+
  ggtitle("Residents' first and last entries, breeding stages, time of day and \\n intra- and extra-pair copulation attempts and full copulations")+ 

#  resident female's first and last entry per day (red lines):
  geom_line(data =min[min$sex =='f',], aes (y=earliest, x=date),  colour ='#F8766D', size =1)  +
  geom_line(data =min[min$sex =='f',], aes (y=latestexit, x=date),  colour ='#F8766D', size =1)  +
#  resident male's first and last entry per day (turquoise lines):  
  geom_line(data =min[min$sex =='m',], aes (y=latestexit, x=date),  colour ='#00BFC4', size =1)  +
  geom_line(data =min[min$sex =='m',], aes (y=earliest, x=date),  colour ='#00BFC4', size =1)  +
# IP copulation attempts (black triangles; i.e. resident male with resident female):
  geom_jitter(data =cops[cops$att.cop%in% 'attempts'& cops$IEP%in% 'IP_IP',], aes(y=time_dec, x=date), shape =2, size =1.5, height =0, width =0.2)+
# IP copulations (black dots; i.e. resident male with resident female):
  geom_jitter(data =cops[cops$att.cop%in% 'copulations'& cops$IEP%in% 'IP_IP',], aes(y=time_dec, x=date), shape =16, size =2, height =0, width =0.2)+
# EP copulation attempts (red triangles, i.e. non-resident male with resident female):  
  geom_jitter(data =cops[cops$att.cop%in% 'attempts'& cops$IEP%in% 'EP_IP',], aes(y=time_dec, x=date), shape =17,colour='red', size =2, height =0, width =0.2)+
# EP copulations (red dots, i.e. non-resident male with resident female):  
  geom_jitter(data =cops[cops$att.cop%in% 'copulations'& cops$IEP%in% 'EP_IP',], aes(y=time_dec, x=date), shape =16, colour='red', size =2, height =0, width =0.2)+
# non-resident IP copulation attempts (blue triangles, i.e. non-resident male with non-resident female):  
  geom_jitter(data =cops[cops$att.cop%in% 'attempts'& cops$IEP%in% 'EP_EP',], aes(y=time_dec, x=date), shape =17,colour='blue', size =2, height =0, width =0.2)+
# non-resident IP copulations (blue dots, i.e. non-resident male with non-resident female):  
  geom_jitter(data =cops[cops$att.cop%in% 'copulations'& cops$IEP%in% 'EP_EP',], aes(y=time_dec, x=date), shape =16, colour='blue',size =2, height =0, width =0.2)



## Fig. S3: Boxplots of IP and EP copulation attempts and full copulations at different breeding stages (summary)
## intra-pair:
ggplot()+
  geom_boxplot(data= IP, aes(x=Female, y=nr_events))+
  geom_point(data= IP, aes(x=Female, y=nr_events,colour=Nest.box, shape =Nest.box),size=2)+
  geom_line(data= IP, aes(x=Female, y=nr_events, colour=Nest.box, group=Nest.box))+
  scale_shape_manual(name  ="Nest-box", breaks=levels(IP$Nest.box),labels=
                       levels(IP$Nest.box), values =c(0:4,6,8,9))+
  facet_wrap(~att.cop)+
  theme_bw()+  scale_x_discrete(name=('Female'))+  
  scale_y_continuous(name=('IP frequency'))
## extra-pair:
ggplot()+
  geom_boxplot(data= EP, aes(x=Female, y=nr_events))+
  geom_point(data= EP, aes(x=Female, y=nr_events,colour=Nest.box, shape =Nest.box),size=2)+
  geom_line(data= EP, aes(x=Female, y=nr_events, colour=Nest.box, group=Nest.box))+
  scale_shape_manual(name  ="Nest-box", breaks=levels(EP$Nest.box),labels=
                       levels(EP$Nest.box), values =c(0:4,6,8,9))+
  facet_wrap(~att.cop)+
  theme_bw()+  scale_x_discrete(name=('Female'))+  
  scale_y_continuous(name=('EP frequency'), breaks = c(0,3,6,9))



## Additional information -----

#time between IP male exit, EP enter and IP male return:
innest <- read.csv('12_allcams_entries.csv'); innest$date <- as.Date(innest$date, format = '%Y-%m-%d')
IP <- innest[innest$IP =='IP',]
IP$cam_male <- paste(IP$cam, "m", sep='_')
IP$cam_female <- paste(IP$cam, "f", sep='_')
innest$resident_male <- IP$cam_male[match(innest$cam, IP$cam)]
innest$resident_female <- IP$cam_female[match(innest$cam, IP$cam)]
innest$date_cam_ring <- paste(innest$date,innest$cam,innest$ring)
innest$timein <- strptime(innest$in.,  format = '%H:%M')
innest$timeout <- strptime(innest$out,  format = '%H:%M')



EP <- innest[innest$IP =='EP',]
EP$sexual <- NA
EP$sexual[EP$intrusion_sexual %in% c('attempt', 'cop','cop, fight', 'copulation')] <- 'yes'
EP$sexual[EP$intrusion_sexual %in% c('fight_both', 'fight_fem','no')] <- 'no'
EP$sexual <- factor(EP$sexual, levels =c('yes','no'))

EP$EPafterIP_hr <- substr(EP$numEpin_after_Ipout,  1,2)
EP$EPafterIP_min <- substr(EP$numEpin_after_Ipout,  4,5)
EP$EPafterIP_sec <- substr(EP$numEpin_after_Ipout,  7,8)
EP$after_hr <- as.numeric(EP$EPafterIP_hr)
EP$after_min <- as.numeric(EP$EPafterIP_min)/60
EP$after_min[EP$EP$EPafterIP_min %in% c('00') ]  <- 0
EP$after_sec <- as.numeric(EP$EPafterIP_sec)
EP$after_minsec <- EP$after_min + EP$after_sec
EP$after <- EP$after_hr + EP$after_min + EP$after_sec
mean(EP$after) # 0.6641026 hours = 39.84615 mins TOTAL
sd(EP$after)*60 # 90.4514 min
median(EP$after)*60 # 16 min
mean(EP[EP$after <= 8,]$after) # 0.4672653 hours = 28.03592 EXCLUDING NIGHTS
sd(EP[EP$after <= 8,]$after)*60 # 33.37693 min
median(EP[EP$after <= 8,]$after)*60 # 15 min

EP$EPbeforeIP_hr <- substr(EP$NumEpout_before_IPin,  1,2)
EP$EPbeforeIP_min <- substr(EP$NumEpout_before_IPin,  4,5)
EP$EPbeforeIP_sec <- substr(EP$NumEpout_before_IPin,  7,8)
EP$before_hr <- as.numeric(EP$EPbeforeIP_hr)
EP$before_min <- as.numeric(EP$EPbeforeIP_min)/60
EP$before_min[EP$EP$EPbeforeIP_min %in% c('00') ]  <- 0
EP$before_sec <- as.numeric(EP$EPbeforeIP_sec)
EP$before_minsec <- EP$before_min + EP$before_sec
EP$before <- EP$before_hr + EP$before_min + EP$before_sec
mean(EP$before) # 0.346401 hours = 20.78406 mins TOTAL
sd(EP$before)*60 # 60.69652 min
median(EP$before)*60 # 6 min
mean(EP[EP$before <= 8,]$before) # 0.2669779 hours = 16.01867 minutes EXCLUDING NIGHTS
sd(EP[EP$before <= 8,]$before)*60 #  24.84523 min
median(EP[EP$before <= 8,]$before)*60 #5 min


## intrusions:
intr <- read.csv('13_2017_intrusions.csv')
intr$date <- as.Date(intr$date, format ='%Y-%m-%d')
table(intr$IP, intr$female_innest) # 1099 events (excluding kestrel and unknown cases), 45 of which when female was inside
intr <- intr[!intr$ring %in% 'kestrel',]
intr <- intr[!intr$ring %in% 'unknown',]






