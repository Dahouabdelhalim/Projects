## load required packages
require(plyr);require(dplyr);require(lme4); require(lmerTest); require(MuMIn); require(ggplot2)


## load in data
load(url("https://github.com/GeladaResearchProject/Tinsley_Johnson_2018_Behav_Ecol/blob/master/seasonality_in_geladas_data.RData?raw=true"))



# There are 5 R objects:
# gc_fem_to: Glucocorticoid metabolites data, with sample values ("GC" - ng/g) for each sample, along with corresponding individual ID code, unit, date of sample collection, age of individual at time of collection, reproductive state ("repro": C (cycling), L (lactating), P (pregnant)), weather variables (mean max and min temp over the 30 days prior to sample collection; total rain over the 90 days prior to sample collection), and whether or not a takeover occurred within the 30 days prior to sample collection ("to_yn"). 

# swell_rate: Number of females that returned to cycling (i.e., first post-partum swelling) for each month ("num_swells") and the number of females in the study population in each month ("num_fem"). We counted the number of females the returned to cycling following takeovers separately (within 90 days, "to_yn"). Each month has the following corresponding weather variables: mean max and min temp over the 30 days prior to sample collection; total rain over the 90 days prior to sample collection.

# concep_rate: The same as "swell_rate" except for the number of conceptions per month.

# birth_rate: The same as "swell_rate" except for the number of births per month, and the takeover window reflects the 90 days around the corresponding conception. 

# death_analysis: For all births to individuals mothers (ID) observed before 2014, the year of the birth, whether the birth occurred within 9 months of a takeover ("DOB.to_yn"); whether the infant died before turning 2 years of age ("INF.DEATH"), whether the infant death was likely due to infanticide ("DOD.to_yn"), and whether the birth occurred during the 3-month ecological peak (Aug-Sep-Oct: "peak3"). 


#############################################################################
###                 GCM Model Selection                                   ###
#############################################################################

gc0<-lmer(log(GC)~
            (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc1<-lmer(log(GC)~
            poly(AGE,2)*repro +
            (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc2<-lmer(log(GC)~
            scale(maxt_1)+
            poly(AGE,2)*repro +
            (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc3<-lmer(log(GC)~
            scale(mint_1)+
            poly(AGE,2)*repro +
            (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc4<-lmer(log(GC)~
            scale(rain)+
            poly(AGE,2)*repro +
            (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc5<-lmer(log(GC)~
            scale(maxt_1)+scale(mint_1)+
            poly(AGE,2)*repro +
            (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc6<-lmer(log(GC)~
            scale(maxt_1)*scale(mint_1)+
            poly(AGE,2)*repro +
            (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc7<-lmer(log(GC)~
            scale(maxt_1)+scale(rain)+
            poly(AGE,2)*repro +
            (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc8<-lmer(log(GC)~
            scale(maxt_1)*scale(rain)+
            poly(AGE,2)*repro +
            (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc9<-lmer(log(GC)~
            scale(mint_1)+scale(rain)+
            poly(AGE,2)*repro +
            (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc10<-lmer(log(GC)~
             scale(mint_1)*scale(rain)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc11<-lmer(log(GC)~
             scale(maxt_1)+scale(mint_1)+scale(rain)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc12<-lmer(log(GC)~
             scale(maxt_1)*scale(mint_1)+
             scale(rain)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc13<-lmer(log(GC)~
             scale(maxt_1)*scale(rain)+
             scale(mint_1)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc14<-lmer(log(GC)~
             scale(mint_1)*scale(rain)+
             scale(maxt_1)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc15<-lmer(log(GC)~
             to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc16<-lmer(log(GC)~
             scale(maxt_1)+to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc17<-lmer(log(GC)~
             scale(maxt_1)*to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc18<-lmer(log(GC)~
             scale(mint_1)+to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc19<-lmer(log(GC)~
             scale(mint_1)*to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc20<-lmer(log(GC)~
             scale(rain)+to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc21<-lmer(log(GC)~
             scale(rain)*to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc22<-lmer(log(GC)~
             scale(maxt_1)+scale(mint_1)+to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc23<-lmer(log(GC)~
             scale(maxt_1)*scale(mint_1)+
             to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc24<-lmer(log(GC)~
             scale(maxt_1)*to_yn+
             scale(mint_1)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc25<-lmer(log(GC)~
             scale(mint_1)*to_yn+
             scale(maxt_1)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc26<-lmer(log(GC)~
             scale(maxt_1)+scale(rain)+to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc27<-lmer(log(GC)~
             scale(maxt_1)*scale(rain)+
             to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc28<-lmer(log(GC)~
             scale(maxt_1)*to_yn+
             scale(rain)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc29<-lmer(log(GC)~
             scale(rain)*to_yn+
             scale(maxt_1)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc30<-lmer(log(GC)~
             scale(mint_1)+scale(rain)+to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc31<-lmer(log(GC)~
             scale(mint_1)*scale(rain)+
             to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc32<-lmer(log(GC)~
             scale(mint_1)*to_yn+
             scale(rain)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc33<-lmer(log(GC)~
             scale(rain)*to_yn+
             scale(mint_1)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc34<-lmer(log(GC)~
             scale(maxt_1)+scale(mint_1)+scale(rain)+to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc35<-lmer(log(GC)~
             scale(maxt_1)*scale(mint_1)+
             scale(rain)+to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc36<-lmer(log(GC)~
             scale(maxt_1)*scale(rain)+
             scale(mint_1)+to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc37<-lmer(log(GC)~
             scale(maxt_1)*to_yn+
             scale(mint_1)+scale(rain)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc38<-lmer(log(GC)~
             scale(mint_1)*scale(rain)+
             scale(maxt_1)+to_yn+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc39<-lmer(log(GC)~
             scale(mint_1)*to_yn+
             scale(maxt_1)+scale(rain)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

gc40<-lmer(log(GC)~
             scale(rain)*to_yn+
             scale(maxt_1)+scale(mint_1)+
             poly(AGE,2)*repro +
             (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to)

## GC MODEL COMPARISONS ----
gc_out.put<-model.sel(gc0,gc1,gc2,gc3,gc4,gc5,gc6,gc7,gc8,gc9,gc10,gc11,gc12,gc13,gc14,gc15,gc16,gc17,gc18,gc19,gc20,gc21,gc22,gc23,gc24,gc25,gc26,gc27,gc28,gc29,gc30,gc31,gc32,gc33,gc34,gc35,gc36,gc37,gc38,gc39,gc40, rank="AIC")

## Table 2: To get the degrees of freedom, log-likelihood, delta-AIC and weight for the top models (where delta-AIC is less than 6):
gc_out.put

## Table 3/Table S1: To get importance of each predictor, and the number of models containing them:
importance(gc_out.put)

## Table 3/Table S1, Fig. 3: To get estimate, adjusted standard error, z-value and p-value for each predictor, see "Model-averaged coefficients: (full average)":
gc_MA.ests<-model.avg(gc_out.put, revised.var = TRUE) 
summary(gc_MA.ests)

## Table 3/Table S1, Fig. 3: To get confidence intervals for all predictors:
confint(gc_MA.ests, full=T)



#############################################################################
###                 Fig. 2bcd: GCM Figures                                ###
#############################################################################

## GCMs: RESIDUALS FOR GRAPHS ----
# Control for effects of reproductive state, age, and takeover (see effect of temp)
gc_fem_to$gc_resid_notemp<-resid(lmer(log(GC)~
                                 poly(AGE,2)*repro +
                                 to_yn +
                                 (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to))

# Control for effects of reproductive state, age, and temperature (see effect of takeover)
gc_fem_to$gc_resid_noto<-resid(lmer(log(GC)~
                                      poly(AGE,2)*repro+
                                      scale(mint_1)+scale(maxt_1)+
                                      (1|ID)+(1|UNIT)+(1|Year),data=gc_fem_to))


gc_graphs_new <- gc_fem_to %>%
  select(Month:GC,gc_resid_notemp,gc_resid_noto)


## Code to make function: "summarySE"----
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
#  library(plyr)
  
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
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Calculate 95% confidence intervals for GCM graphs ----
gc_graphs_summary<-summarySE(gc_graphs_new, measurevar="gc_resid_notemp", groupvars=c("Month")) ## gives you 95% confidence intervals

## FIG 2B ----
ggplot(gc_graphs_summary, aes(x=as.factor(Month),y=gc_resid_notemp, group=1)) +
  geom_rect(data = gc_graphs_summary, aes(xmin = Month-0.5, xmax = Month+0.5, ymin = -Inf, ymax = Inf, fill = month), fill=c("#abd9e9","white","#fdae61","#fdae61","#fdae61","#fdae61","white","#abd9e9","#abd9e9","white","white","#abd9e9"), alpha=0.35) +
  scale_x_discrete(labels=c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  geom_text(aes(x = 4.5, y = 0.14, label = "HOT")) +
  geom_text(aes(x = 8.5, y = -0.1, label = "COLD")) +
  geom_text(aes(x = 12, y = -0.1, label = "COLD")) +
  geom_text(aes(x = 1, y = -0.1, label = "COLD")) +
  theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0)) +
  labs(x = "Months (2006-2014)", y = "Residual GCMs", title = "") +
  geom_errorbar(aes(ymin=gc_resid_notemp-ci, ymax=gc_resid_notemp+ci), width=.1) +
  geom_line() +
  geom_point(shape=20, size=4) 

## PREP FOR FIG 2C and FIG 2D ----
gc_to<-summarySE(gc_fem_to, measurevar=c("gc_resid_noto"), groupvars=c("to_yn")) ## gives you 95% confidence intervals

gc_to[,1][gc_to[,1] == "N"] <- "NO\\nTAKEOVER"
gc_to[,1][gc_to[,1] == "Y"] <- "TAKEOVER"

gc_fem_to$temp<-gc_fem_to$Month
gc_fem_to[,15][gc_fem_to[,15] == "1"] <- "COLD"
gc_fem_to[,15][gc_fem_to[,15] == "2"] <- "MED"
gc_fem_to[,15][gc_fem_to[,15] == "3"] <- "HOT"
gc_fem_to[,15][gc_fem_to[,15] == "4"] <- "HOT"
gc_fem_to[,15][gc_fem_to[,15] == "5"] <- "HOT"
gc_fem_to[,15][gc_fem_to[,15] == "6"] <- "HOT"
gc_fem_to[,15][gc_fem_to[,15] == "7"] <- "MED"
gc_fem_to[,15][gc_fem_to[,15] == "8"] <- "COLD"
gc_fem_to[,15][gc_fem_to[,15] == "9"] <- "COLD"
gc_fem_to[,15][gc_fem_to[,15] == "10"] <- "MED"
gc_fem_to[,15][gc_fem_to[,15] == "11"] <- "MED"
gc_fem_to[,15][gc_fem_to[,15] == "12"] <- "COLD"

gc_temp<-summarySE(gc_fem_to[gc_fem_to$temp!="MED",], measurevar=c("gc_resid_notemp"), groupvars=c("temp")) ## gives you 95% confidence intervals

gc_temp$temp=as.factor(gc_temp$temp)
gc_temp$temp = factor(gc_temp$temp, levels(gc_temp$temp)[c(2,1)])

d_to<-data.frame(x=gc_to$to_yn, y=gc_to$gc_resid_noto, ci=gc_to$ci, se=gc_to$se)
d_temp<-data.frame(x=gc_temp$temp, y=gc_temp$gc_resid_notemp, ci=gc_temp$ci, se=gc_temp$se)

d_to$panel <- "c. Takeover"
d_temp$panel <- "d. Temperature"


d_gc <- rbind(d_to, d_temp)

d_gc$x=as.factor(d_gc$x)


## FIG 2C and FIG 2D ----
ggplot(data = d_gc, mapping = aes(x = as.factor(x), y = y, fill=x)) +
  facet_wrap(~panel, ncol=2, scales="free_x") + 
  geom_hline(yintercept = 0, colour = "#636363", linetype = 3) +
  xlab("") +
  ylab("Residual GCMs") +
  geom_boxplot(aes(ymin = y-ci, ymax = y+ci, upper = y+se, middle = y, lower = y-se, group = panel, fill=x, shape=panel, width=0.5), stat="identity", colour = "#636363", show.legend=FALSE) +
  scale_fill_manual(values = c("grey","#ffffff","#fdae61","#abd9e9")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), plot.title = element_text(hjust = 0), strip.text.x=element_text(face = "bold",size = "16",hjust = 0)) 

#############################################################################
###                 Cycling Model Selection                               ###
#############################################################################

cy0<-glmer(cbind(num_swells,num_fem-num_swells)~
             (1|mth_yr),family="binomial",data=swell_rate)

cy1<-glmer(cbind(num_swells,num_fem-num_swells)~
             scale(avg_maxt)+
             (1|mth_yr),family="binomial",data=swell_rate)

cy2<-glmer(cbind(num_swells,num_fem-num_swells)~
             scale(avg_mint)+
             (1|mth_yr),family="binomial",data=swell_rate)

cy3<-glmer(cbind(num_swells,num_fem-num_swells)~
             scale(avg_rain)+
             (1|mth_yr),family="binomial",data=swell_rate)

cy4<-glmer(cbind(num_swells,num_fem-num_swells)~
             scale(avg_maxt)+scale(avg_mint)+
             (1|mth_yr),family="binomial",data=swell_rate)

cy5<-glmer(cbind(num_swells,num_fem-num_swells)~
             scale(avg_maxt)*scale(avg_mint)+
             (1|mth_yr),family="binomial",data=swell_rate)

cy6<-glmer(cbind(num_swells,num_fem-num_swells)~
             scale(avg_maxt)+scale(avg_rain)+
             (1|mth_yr),family="binomial",data=swell_rate)

cy7<-glmer(cbind(num_swells,num_fem-num_swells)~
             scale(avg_maxt)*scale(avg_rain)+
             (1|mth_yr),family="binomial",data=swell_rate)

cy8<-glmer(cbind(num_swells,num_fem-num_swells)~
             scale(avg_mint)+scale(avg_rain)+
             (1|mth_yr),family="binomial",data=swell_rate)

cy9<-glmer(cbind(num_swells,num_fem-num_swells)~
             scale(avg_mint)*scale(avg_rain)+
             (1|mth_yr),family="binomial",data=swell_rate)

cy10<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)+scale(avg_mint)+scale(avg_rain)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy11<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)*scale(avg_mint)+
              scale(avg_rain)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy12<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)*scale(avg_rain)+
              scale(avg_mint)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy13<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_mint)*scale(avg_rain)+
              scale(avg_maxt)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy14<-glmer(cbind(num_swells,num_fem-num_swells)~
              to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy15<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)+to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy16<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)*to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy17<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_mint)+to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy18<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_mint)*to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy19<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_rain)+to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy20<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_rain)*to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy21<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)+scale(avg_mint)+to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy22<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)*scale(avg_mint)+
              to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy23<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)*to_yn+
              scale(avg_mint)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy24<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_mint)*to_yn+
              scale(avg_maxt)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy25<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)+scale(avg_rain)+to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy26<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)*scale(avg_rain)+
              to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy27<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)*to_yn+
              scale(avg_rain)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy28<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_rain)*to_yn+
              scale(avg_maxt)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy29<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_mint)+scale(avg_rain)+to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy30<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_mint)*scale(avg_rain)+
              to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy31<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_mint)*to_yn+
              scale(avg_rain)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy32<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_rain)*to_yn+
              scale(avg_mint)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy33<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)+scale(avg_mint)+scale(avg_rain)+to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy34<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)*scale(avg_mint)+
              scale(avg_rain)+to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy35<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)*scale(avg_rain)+
              scale(avg_mint)+to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy36<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_maxt)*to_yn+
              scale(avg_mint)+scale(avg_rain)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy37<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_mint)*scale(avg_rain)+
              scale(avg_maxt)+to_yn+
              (1|mth_yr),family="binomial",data=swell_rate)

cy38<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_mint)*to_yn+
              scale(avg_maxt)+scale(avg_rain)+
              (1|mth_yr),family="binomial",data=swell_rate)

cy39<-glmer(cbind(num_swells,num_fem-num_swells)~
              scale(avg_rain)*to_yn+
              scale(avg_maxt)+scale(avg_mint)+
              (1|mth_yr),family="binomial",data=swell_rate)

## CYCLING MODEL COMPARISONS ----
cy_out.put<-model.sel(cy0,cy1,cy2,cy3,cy4,cy5,cy6,cy7,cy8,cy9,cy10,cy11,cy12,cy13,cy14,cy15,cy16,cy17,cy18,cy19,cy20,cy21,cy22,cy23,cy24,cy25,cy26,cy27,cy28,cy29,cy30,cy31,cy32,cy33,cy34,cy35,cy36,cy37,cy38,cy39)

## Table 2: To get the degrees of freedom, log-likelihood, delta-AICc and weight for the top models (where delta-AICc is less than 6):
cy_out.put

## Table 3/Table S1: To get importance of each predictor, and the number of models containing them:
importance(cy_out.put)

## Table 3/Table S1, Fig. 3: To get estimate, adjusted standard error, z-value and p-value for each predictor, see "Model-averaged coefficients: (full average)":
cy_MA.ests<-model.avg(cy_out.put, revised.var = TRUE) 
summary(cy_MA.ests)

## Table 3/Table S1, Fig. 3: To get confidence intervals for all predictors:
confint(cy_MA.ests, full=T)



#############################################################################
###                 Conceptions Model Selection                           ###
#############################################################################

c0<-glmer(cbind(num_concep,num_fem-num_concep)~
            (1|mth_yr),family="binomial",data=concep_rate)

c1<-glmer(cbind(num_concep,num_fem-num_concep)~
            scale(avg_maxt)+
            (1|mth_yr),family="binomial",data=concep_rate)

c2<-glmer(cbind(num_concep,num_fem-num_concep)~
            scale(avg_mint)+
            (1|mth_yr),family="binomial",data=concep_rate)

c3<-glmer(cbind(num_concep,num_fem-num_concep)~
            scale(avg_rain)+
            (1|mth_yr),family="binomial",data=concep_rate)

c4<-glmer(cbind(num_concep,num_fem-num_concep)~
            scale(avg_maxt)+scale(avg_mint)+
            (1|mth_yr),family="binomial",data=concep_rate)

c5<-glmer(cbind(num_concep,num_fem-num_concep)~
            scale(avg_maxt)*scale(avg_mint)+
            (1|mth_yr),family="binomial",data=concep_rate)

c6<-glmer(cbind(num_concep,num_fem-num_concep)~
            scale(avg_maxt)+scale(avg_rain)+
            (1|mth_yr),family="binomial",data=concep_rate)

c7<-glmer(cbind(num_concep,num_fem-num_concep)~
            scale(avg_maxt)*scale(avg_rain)+
            (1|mth_yr),family="binomial",data=concep_rate)

c8<-glmer(cbind(num_concep,num_fem-num_concep)~
            scale(avg_mint)+scale(avg_rain)+
            (1|mth_yr),family="binomial",data=concep_rate)

c9<-glmer(cbind(num_concep,num_fem-num_concep)~
            scale(avg_mint)*scale(avg_rain)+
            (1|mth_yr),family="binomial",data=concep_rate)

c10<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)+scale(avg_mint)+scale(avg_rain)+
             (1|mth_yr),family="binomial",data=concep_rate)

c11<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)*scale(avg_mint)+
             scale(avg_rain)+
             (1|mth_yr),family="binomial",data=concep_rate)

c12<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)*scale(avg_rain)+
             scale(avg_mint)+
             (1|mth_yr),family="binomial",data=concep_rate)

c13<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_mint)*scale(avg_rain)+
             scale(avg_maxt)+
             (1|mth_yr),family="binomial",data=concep_rate)

c14<-glmer(cbind(num_concep,num_fem-num_concep)~
             to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c15<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)+to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c16<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)*to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c17<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_mint)+to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c18<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_mint)*to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c19<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_rain)+to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c20<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_rain)*to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c21<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)+scale(avg_mint)+to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c22<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)*scale(avg_mint)+
             to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c23<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)*to_yn+
             scale(avg_mint)+
             (1|mth_yr),family="binomial",data=concep_rate)

c24<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_mint)*to_yn+
             scale(avg_maxt)+
             (1|mth_yr),family="binomial",data=concep_rate)

c25<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)+scale(avg_rain)+to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c26<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)*scale(avg_rain)+
             to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c27<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)*to_yn+
             scale(avg_rain)+
             (1|mth_yr),family="binomial",data=concep_rate)

c28<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_rain)*to_yn+
             scale(avg_maxt)+
             (1|mth_yr),family="binomial",data=concep_rate)

c29<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_mint)+scale(avg_rain)+to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c30<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_mint)*scale(avg_rain)+
             to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c31<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_mint)*to_yn+
             scale(avg_rain)+
             (1|mth_yr),family="binomial",data=concep_rate)

c32<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_rain)*to_yn+
             scale(avg_mint)+
             (1|mth_yr),family="binomial",data=concep_rate)

c33<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)+scale(avg_mint)+scale(avg_rain)+to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c34<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)*scale(avg_mint)+
             scale(avg_rain)+to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c35<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)*scale(avg_rain)+
             scale(avg_mint)+to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c36<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_maxt)*to_yn+
             scale(avg_mint)+scale(avg_rain)+
             (1|mth_yr),family="binomial",data=concep_rate)

c37<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_mint)*scale(avg_rain)+
             scale(avg_maxt)+to_yn+
             (1|mth_yr),family="binomial",data=concep_rate)

c38<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_mint)*to_yn+
             scale(avg_maxt)+scale(avg_rain)+
             (1|mth_yr),family="binomial",data=concep_rate)

c39<-glmer(cbind(num_concep,num_fem-num_concep)~
             scale(avg_rain)*to_yn+
             scale(avg_maxt)+scale(avg_mint)+
             (1|mth_yr),family="binomial",data=concep_rate)

## CONCEPTION MODEL COMPARISONS ----
c_out.put<-model.sel(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,c35,c36,c37,c38,c39)

## Table 2: To get the degrees of freedom, log-likelihood, delta-AICc and weight for the top models (where delta-AICc is less than 6):
c_out.put

## Table 3/Table S1: To get importance of each predictor, and the number of models containing them:
importance(c_out.put)

## Table 3/Table S1, Fig. 3: To get estimate, adjusted standard error, z-value and p-value for each predictor, see "Model-averaged coefficients: (full average)":
c_MA.ests<-model.avg(c_out.put, revised.var = TRUE) 
summary(c_MA.ests)

## Table 3/Table S1, Fig. 3: To get confidence intervals for all predictors:
confint(c_MA.ests, full=T)


#############################################################################
###                 Births Model Selection                                ###
#############################################################################

b0<-glmer(cbind(num_births,num_fem-num_births)~
            (1|mth_yr),family="binomial",data=birth_rate)

b1<-glmer(cbind(num_births,num_fem-num_births)~
            scale(avg_maxt)+
            (1|mth_yr),family="binomial",data=birth_rate)

b2<-glmer(cbind(num_births,num_fem-num_births)~
            scale(avg_mint)+
            (1|mth_yr),family="binomial",data=birth_rate)

b3<-glmer(cbind(num_births,num_fem-num_births)~
            scale(avg_rain)+
            (1|mth_yr),family="binomial",data=birth_rate)

b4<-glmer(cbind(num_births,num_fem-num_births)~
            scale(avg_maxt)+scale(avg_mint)+
            (1|mth_yr),family="binomial",data=birth_rate)

b5<-glmer(cbind(num_births,num_fem-num_births)~
            scale(avg_maxt)*scale(avg_mint)+
            (1|mth_yr),family="binomial",data=birth_rate)

b6<-glmer(cbind(num_births,num_fem-num_births)~
            scale(avg_maxt)+scale(avg_rain)+
            (1|mth_yr),family="binomial",data=birth_rate)

b7<-glmer(cbind(num_births,num_fem-num_births)~
            scale(avg_maxt)*scale(avg_rain)+
            (1|mth_yr),family="binomial",data=birth_rate)

b8<-glmer(cbind(num_births,num_fem-num_births)~
            scale(avg_mint)+scale(avg_rain)+
            (1|mth_yr),family="binomial",data=birth_rate)

b9<-glmer(cbind(num_births,num_fem-num_births)~
            scale(avg_mint)*scale(avg_rain)+
            (1|mth_yr),family="binomial",data=birth_rate)

b10<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)+scale(avg_mint)+scale(avg_rain)+
             (1|mth_yr),family="binomial",data=birth_rate)

b11<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)*scale(avg_mint)+
             scale(avg_rain)+
             (1|mth_yr),family="binomial",data=birth_rate)

b12<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)*scale(avg_rain)+
             scale(avg_mint)+
             (1|mth_yr),family="binomial",data=birth_rate)

b13<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_mint)*scale(avg_rain)+
             scale(avg_maxt)+
             (1|mth_yr),family="binomial",data=birth_rate)

b14<-glmer(cbind(num_births,num_fem-num_births)~
             to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b15<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)+to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b16<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)*to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b17<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_mint)+to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b18<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_mint)*to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b19<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_rain)+to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b20<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_rain)*to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b21<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)+scale(avg_mint)+to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b22<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)*scale(avg_mint)+
             to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b23<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)*to_yn+
             scale(avg_mint)+
             (1|mth_yr),family="binomial",data=birth_rate)

b24<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_mint)*to_yn+
             scale(avg_maxt)+
             (1|mth_yr),family="binomial",data=birth_rate)

b25<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)+scale(avg_rain)+to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b26<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)*scale(avg_rain)+
             to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b27<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)*to_yn+
             scale(avg_rain)+
             (1|mth_yr),family="binomial",data=birth_rate)

b28<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_rain)*to_yn+
             scale(avg_maxt)+
             (1|mth_yr),family="binomial",data=birth_rate)

b29<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_mint)+scale(avg_rain)+to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b30<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_mint)*scale(avg_rain)+
             to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b31<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_mint)*to_yn+
             scale(avg_rain)+
             (1|mth_yr),family="binomial",data=birth_rate)

b32<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_rain)*to_yn+
             scale(avg_mint)+
             (1|mth_yr),family="binomial",data=birth_rate)

b33<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)+scale(avg_mint)+scale(avg_rain)+to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b34<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)*scale(avg_mint)+
             scale(avg_rain)+to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b35<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)*scale(avg_rain)+
             scale(avg_mint)+to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b36<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_maxt)*to_yn+
             scale(avg_mint)+scale(avg_rain)+
             (1|mth_yr),family="binomial",data=birth_rate)

b37<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_mint)*scale(avg_rain)+
             scale(avg_maxt)+to_yn+
             (1|mth_yr),family="binomial",data=birth_rate)

b38<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_mint)*to_yn+
             scale(avg_maxt)+scale(avg_rain)+
             (1|mth_yr),family="binomial",data=birth_rate)

b39<-glmer(cbind(num_births,num_fem-num_births)~
             scale(avg_rain)*to_yn+
             scale(avg_maxt)+scale(avg_mint)+
             (1|mth_yr),family="binomial",data=birth_rate)

## BIRTH MODEL COMPARISONS ----
b_out.put<-model.sel(b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22,b23,b24,b25,b26,b27,b28,b29,b30,b31,b32,b33,b34,b35,b36,b37,b38,b39)

## Table 2: To get the degrees of freedom, log-likelihood, delta-AICc and weight for the top models (where delta-AICc is less than 6):
b_out.put

## Table 3/Table S1: To get importance of each predictor, and the number of models containing them:
importance(b_out.put)

## Table 3/Table S1, Fig. 3: To get estimate, adjusted standard error, z-value and p-value for each predictor, see "Model-averaged coefficients: (full average)":
b_MA.ests<-model.avg(b_out.put, revised.var = TRUE) 
summary(b_MA.ests)

## Table 3/Table S1, Fig. 3: To get confidence intervals for all predictors:
confint(b_MA.ests, full=T)


#############################################################################
###                 Fig. 4ab: Reproductive events                         ###
#############################################################################
## Swell prep ----
swell_rate$to_label <- revalue(swell_rate$to_yn, c("N"="No Takeover", "Y"="Takeover"))

swell_graphs_mth<-ddply(swell_rate,.(month,to_label),summarize,n=sum(num_swells),avg_rate=mean(swell_rate),sd=sd(swell_rate),se=sd(swell_rate)/sqrt(length(swell_rate)))

swell_graphs_mth$to_label=as.factor(swell_graphs_mth$to_label)


## Concep prep ----
concep_rate$to_label <- revalue(concep_rate$to_yn, c("N"="No Takeover", "Y"="Takeover"))

concep_graphs_mth<-ddply(concep_rate,.(month,to_label),summarize,n=sum(num_concep),avg_rate=mean(concep_rate),sd=sd(concep_rate),se=sd(concep_rate)/sqrt(length(concep_rate)))

concep_graphs_mth$to_label=as.factor(concep_graphs_mth$to_label)


## Birth prep ----
birth_rate$to_label <- revalue(birth_rate$to_yn, c("N"="No Takeover", "Y"="Takeover"))

birth_graphs_mth<-ddply(birth_rate,.(month,to_label),summarize,n=sum(num_births),avg_rate=mean(birth_rate),sd=sd(birth_rate),se=sd(birth_rate)/sqrt(length(birth_rate)))

birth_graphs_mth$to_label=as.factor(birth_graphs_mth$to_label)


## Create datasets for graphs ----
d1a<-data.frame(x=as.factor(swell_graphs_mth[swell_graphs_mth$to_label=="No Takeover",]$month), y=swell_graphs_mth[swell_graphs_mth$to_label=="No Takeover",]$avg_rate, num=swell_graphs_mth[swell_graphs_mth$to_label=="No Takeover",]$n, z="Return to cycling", se=swell_graphs_mth[swell_graphs_mth$to_label=="No Takeover",]$se)
d1b<-data.frame(x=as.factor(concep_graphs_mth[concep_graphs_mth$to_label=="No Takeover",]$month), y=concep_graphs_mth[concep_graphs_mth$to_label=="No Takeover",]$avg_rate, num=concep_graphs_mth[concep_graphs_mth$to_label=="No Takeover",]$n, z="Conceptions", se=concep_graphs_mth[concep_graphs_mth$to_label=="No Takeover",]$se)
d1c<-data.frame(x=as.factor(birth_graphs_mth[birth_graphs_mth$to_label=="No Takeover",]$month), y=birth_graphs_mth[birth_graphs_mth$to_label=="No Takeover",]$avg_rate, num=birth_graphs_mth[birth_graphs_mth$to_label=="No Takeover",]$n, z="Births", se=birth_graphs_mth[birth_graphs_mth$to_label=="No Takeover",]$se)

d2a<-data.frame(x=as.factor(swell_graphs_mth[swell_graphs_mth$to_label=="Takeover",]$month), y=swell_graphs_mth[swell_graphs_mth$to_label=="Takeover",]$avg_rate, num=swell_graphs_mth[swell_graphs_mth$to_label=="Takeover",]$n, z="Return to cycling", se=swell_graphs_mth[swell_graphs_mth$to_label=="Takeover",]$se)
d2b<-data.frame(x=as.factor(concep_graphs_mth[concep_graphs_mth$to_label=="Takeover",]$month), y=concep_graphs_mth[concep_graphs_mth$to_label=="Takeover",]$avg_rate, num=concep_graphs_mth[concep_graphs_mth$to_label=="Takeover",]$n, z="Conceptions", se=concep_graphs_mth[concep_graphs_mth$to_label=="Takeover",]$se)
d2c<-data.frame(x=as.factor(birth_graphs_mth[birth_graphs_mth$to_label=="Takeover",]$month), y=birth_graphs_mth[birth_graphs_mth$to_label=="Takeover",]$avg_rate, num=birth_graphs_mth[birth_graphs_mth$to_label=="Takeover",]$n, z="Births", se=birth_graphs_mth[birth_graphs_mth$to_label=="Takeover",]$se)

d1a$panel <- "No Takeover"
d1b$panel <- "No Takeover"
d1c$panel <- "No Takeover"
d2a$panel <- "Takeover"
d2b$panel <- "Takeover"
d2c$panel <- "Takeover"

F3_ecol <- rbind(d1a, d1b, d1c)
F3_soc <- rbind(d2a, d2b, d2c)

F3 <- rbind(F3_ecol, F3_soc)

pd <- position_dodge(0.1)


## FIGURE 4: ----
# 4a: No takeover
ggplot(data = F3, mapping = aes(x = as.factor(x))) + 
  facet_wrap(~z, ncol=3) + 
  scale_x_discrete(labels=c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  geom_bar(data=F3[F3$panel=="No Takeover",], stat="identity", aes(y=num, group=z, fill=z), colour="#636363", show.legend=FALSE) +
  geom_errorbar(data=F3[F3$panel=="No Takeover",], aes(ymin=(y-se)/0.002, ymax=(y+se)/0.002, group=z), width=.1) +
  geom_line(data=F3[F3$panel=="No Takeover",],aes(y=y/0.002, group=z)) +
  geom_point(data=F3[F3$panel=="No Takeover",],aes(y=y/0.002, group=z), shape=15, size = 3, show.legend=FALSE) +
  scale_y_continuous(sec.axis=sec_axis(~.*0.002,name="Mean rate\\n(n per female)")) +
  scale_fill_manual(values=c("#f7f7f7","#f7f7f7","#f7f7f7"), guide=FALSE) +
  labs(x = "", y = "Total observed", title = "a. No takeover") +
  theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0)) 

# 4b: Takeover
ggplot(data = F3, mapping = aes(x = as.factor(x))) + 
  facet_wrap(~z, ncol=3) + 
  scale_x_discrete(labels=c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  geom_vline(xintercept=3, linetype=2) +
  geom_text(aes(x = 5.5, y = 17.75, label = "Takeover Peak")) +
  geom_bar(data=F3[F3$panel=="Takeover",], stat="identity", aes(y=num, group=z, fill=z), colour="#636363", show.legend=FALSE) +
  geom_errorbar(data=F3[F3$panel=="Takeover",], aes(ymin=(y-se)/0.02, ymax=(y+se)/0.02, group=z), width=.1) +
  geom_line(data=F3[F3$panel=="Takeover",],aes(y=y/0.02, group=z)) +
  geom_point(data=F3[F3$panel=="Takeover",],aes(y=y/0.02, group=z), shape=15, size = 3, show.legend=FALSE) +
  scale_y_continuous(sec.axis=sec_axis(~.*0.02,name="Mean rate\\n(n per female)")) +
  scale_fill_manual(values=c("#f7f7f7","#f7f7f7","#f7f7f7"), guide=FALSE) +
  labs(x = "Months (2006-2014)", y = "Total observed", title = "b. Takeover") +
  theme(legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0)) 


#############################################################################
###                 Infant Death Models                                   ###
#############################################################################
## All births
summary(glmer(as.factor(INF.DEATH) ~
                as.factor(peak3) +
                (1|DOB.yr) + (1|ID),family="binomial",data=death_analysis))

## No infanticides
summary(glmer(as.factor(INF.DEATH) ~
                as.factor(peak3) +
                (1|DOB.yr) + (1|ID),family="binomial",data=death_analysis[death_analysis$DOD.to_yn!="Y",]))

## No takeover pre-birth ('ecological peak') AND no infanticides
summary(glmer(as.factor(INF.DEATH) ~
                as.factor(peak3) +
                (1|DOB.yr) + (1|ID),family="binomial",data=death_analysis[death_analysis$DOB.to_yn=="N"&death_analysis$DOD.to_yn!="Y",]))

