#Van der Burg et al.
#Genomic architecture of a genetically assimilated seasonal color pattern
#Supplemental script
#hormone analysis

horm<-read.table("aaz3017_vanderburg_etal_Ecdysone_HPLC_results_RED_PLAS.txt", header=T)
head(horm)
attach(horm)


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

hormsum<-ddply(horm,c("region.temp","point"),summarise, 
                N = sum(!is.na(X20_E_norm)),
                mean = mean(X20_E_norm, na.rm=TRUE), 
                sd = sd(X20_E_norm,na.rm=TRUE),
                se = sd / sqrt(N))
hormdf<-data.frame(horm)
attach(hormdf)


t.test(x=horm[horm$BATCH=="one"&horm$Line=="NC",2],
       y=horm[horm$BATCH=="two"&horm$Line=="NC",2])




t.test(x=horm[horm$region.temp=="nc.w"&horm$point=="a.three",2],
       y=horm[horm$region.temp=="red.w"&horm$point=="a.three",2])

t.test(x=horm[horm$region.temp=="nc.w"&horm$point=="b.five",2],
       y=horm[horm$region.temp=="red.w"&horm$point=="b.five",2])

t.test(x=horm[horm$region.temp=="nc.w"&horm$point=="c.seven",2],
       y=horm[horm$region.temp=="red.w"&horm$point=="c.seven",2])

t.test(x=horm[horm$region.temp=="nc.w"&horm$point=="d.nine",2],
       y=horm[horm$region.temp=="red.w"&horm$point=="d.nine",2])



hormdfc<-summarySE(hormdf, measurevar="X20.OH.Ecdysone", groupvars=c("region.temp","point"), na.rm=TRUE)
attach(hormdfc)

horm_red<-subset(hormdfc,region.temp=="nc.w"|region.temp=="red.w"|region.temp=="nc.c")
pd <- position_dodge(.1)

ggplot(horm_red, aes(x=point, y=X20.OH.Ecdysone, shape=region.temp, colour=region.temp, ymax = max(270))) + 
  geom_point(position=pd, size =3 ) +
  geom_line(position=pd, aes(group=region.temp)) +
  geom_errorbar(aes(ymin=X20.OH.Ecdysone-se, ymax=X20.OH.Ecdysone + se), width=.1, position=pd) +  theme_bw() +
  
  theme_bw() + 
  theme(legend.justification=c(1,0), legend.position=c(.93,.65),
        legend.title=element_text(size=12),
        legend.text=(element_text(size=12))) +
  scale_color_manual(name = "", 
                     labels = c("cold","warm","red"), 
                     values = c("#397438",'#72499d',"#cb2026")) +
  scale_shape_manual(name = "",
                     labels = c("cold","warm","red"), 
                     values = c(16, 15,17)) +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.line.x = element_line(size = .5, color = 'black', linetype = 1)) +
  theme(axis.line.y = element_line(size = .5, color = 'black', linetype = 1)) +
  theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size=12), axis.text.y = element_text(size=12)) +
  xlab("relative development time") + ylab("20-OH Ecdysone Titer (pg/uL)") +
  scale_x_discrete(breaks=c("a.three", "b.five", "c.seven",'d.nine'), 
                   labels=c("12%", "24%", "36%","48%")) +
  scale_y_continuous(breaks=seq(0,270,30))
