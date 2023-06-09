#Script for 
#Timbres of tolerance: Co-feeding tolerance in chimpanzees and bonobos is better explained by group-specific variation than species differences 
#(van Leeuwen et al., PRSB 2023)

rm(list=ls())
setwd("~/Desktop/Co-feeding tolerance MS_PRSB/Analyses")
xdata=read.table("van_Leeuwen_et_al._dataset_co-feeding_tolerance.txt", header=T, sep="\\t")

#############
str(xdata)
source("diagnostic_fcns.r")
source("helpers.r")
nrow(xdata)
unique(xdata$group)
xdata$group=as.factor(xdata$group)
library(ggplot2)
library(ggthemes)

xdata$ids_not=xdata$group_size-xdata$ids_nr
xdata$ids_prop=xdata$ids_nr/xdata$group_size
sum(is.na(xdata$ids_prop))
xdata=subset(xdata, !is.na(xdata$ids_prop))

#preparing variables (z-transforming covariates, dummy coding factors for inclusiong in random effect structure)
xdata$z.scan=as.vector(scale(xdata$scan))
xdata$gr.session=as.factor(paste(xdata$group, xdata$session, sep="_"))
range(table(xdata$session))
xdata$group=as.factor(xdata$group)
xdata$species=ifelse(xdata$group=="BBG1" | xdata$group=="BBG2" | xdata$group=="3" | xdata$group=="4" | xdata$group=="1" | xdata$group=="2" | xdata$group=="A-Chimps" | xdata$group=="Antwerp" | xdata$group=="chimpB", "chimpanzee", "bonobo")
levels(xdata$group)=c("Chimp_C1","Chimp_C2","Chimp_C3", "Chimp_C4", "ChimpA_LPZ", "Chimp_Antwerp", "Chimp_BB1", "Chimp_BB2", "Bonobo_LPZ", "ChimpB_LPZ", "Bonobo_F1", "Bonobo_F2", "Bonobo_Lola1","Bonobo_Lola2","Bonobo_Lola3","Bonobo_PLD")
table(xdata$species, xdata$group)
xdata$species=as.factor(xdata$species)
xdata$setting=as.factor(xdata$setting)

#MODELING COFEEDING TOLERANCE

#Roger Mundry's answer (see acknowledgements paper): not with beta error bc beta models continuous proportions while we are dealing with discrete binomials! but check under/overdispersion

#########

#MODEL 1 -- SPECIES EFFECTS

library(lme4)
contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))

full=glmer(cbind(ids_nr, ids_not)~z.scan*species*setting+(1+z.scan||gr.session)+(1+z.scan||group), data=xdata, family=binomial, control=contr)

#cbind(fixef(full),confint(full))

# Computing profile confidence intervals ...
                                                      # 2.5 %    97.5 %
# .sig01                              -0.51400835  0.24569381 0.3588517
# .sig02                              -0.22308254  0.06035447 0.1681963
# .sig03                               0.34705997  0.41000323 0.9092227
# .sig04                              -0.78164415  0.24309605 0.5611446
# (Intercept)                         -0.06449846 -1.24219026 0.2165706
# z.scan                              -0.23175661 -0.66818818 0.2235118
# specieschimpanzee                   -0.84204061 -0.62065076 1.3080156
# settingzoo                          -0.52550418 -2.01496672 0.1854922
# z.scan:specieschimpanzee            -0.51400835 -0.65678587 0.5208260
# z.scan:settingzoo                   -0.22308254 -0.84862018 0.3640892
# specieschimpanzee:settingzoo         0.34705997 -2.13647438 0.4682869
# z.scan:specieschimpanzee:settingzoo -0.78164415 -1.32747176 0.2833701

null=glmer(cbind(ids_nr, ids_not)~1+(1+z.scan||gr.session)+(1+z.scan||group), data=xdata, family=binomial, control=contr)

anova(full, null, test="Chisq")
summary(full)

# Data: xdata
# Models:
# null: cbind(ids_nr, ids_not) ~ 1 + (1 + z.scan || gr.session) + (1 + z.scan || group)
# full: cbind(ids_nr, ids_not) ~ z.scan * species * setting + (1 + z.scan || gr.session) + (1 + z.scan || group)
     # npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# null    5 4269.7 4295.1 -2129.9   4259.7                         
# full   12 4250.4 4311.2 -2113.2   4226.4 33.322  7  2.306e-05 ***

singles=drop1(full, test="Chisq")
# Model:
# cbind(ids_nr, ids_not) ~ z.scan * species * setting + (1 + z.scan || 
    # gr.session) + (1 + z.scan || group)
                       # npar    AIC    LRT Pr(Chi)
# <none>                      4250.4               
# z.scan:species:setting    1 4250.2 1.7584  0.1848

red=glmer(cbind(ids_nr, ids_not)~z.scan*species*setting-
z.scan:species:setting
+(1+z.scan||gr.session)+(1+z.scan||group), data=xdata, family=binomial, control=contr)

singles=drop1(red, test="Chisq")
#sum(is.na(xdata))

# z.scan:species     1 4250.8 2.6330 0.10467  
# z.scan:setting     1 4254.0 5.8109 0.01593 *
# species:setting    1 4249.8 1.6103 0.20445 

#strategies differ across settings, logical consequence of plot VS swing

#test for SETTING EFFECT

red2=glmer(cbind(ids_nr, ids_not)~z.scan+species+setting+(1+z.scan||gr.session)+(1+z.scan||group), data=xdata, family=binomial, control=contr)

summary(red2)

cbind(fixef(red2),confint(red2))
#									2.5 %     97.5 %
# (Intercept)       -0.26271284 -0.92349135  0.3954923
# z.scan            -0.55779106 -0.83679158 -0.2889004
# specieschimpanzee -0.09399814 -0.79721508  0.6139545
# settingzoo        -1.26904367 -1.98166742 -0.5729780

singles=drop1(red2, test="Chisq")

# Model:
# cbind(ids_nr, ids_not) ~ z.scan + species + setting + (1 + z.scan || 
    # gr.session) + (1 + z.scan || group)
        # npar    AIC     LRT   Pr(Chi)    
# <none>       4253.3                      
# z.scan     1 4263.6 12.3125 0.0004499 ***
# species    1 4251.4  0.0776 0.7805118    
# setting    1 4261.5 10.2075 0.0013987 ** 

summary(red2)

# Random effects:
 # Groups       Name        Variance Std.Dev.
 # gr.session   (Intercept) 0.08852  0.2975  
 # gr.session.1 z.scan      0.01318  0.1148  
 # group        (Intercept) 0.42122  0.6490  
 # group.1      z.scan      0.25786  0.5078  
# Number of obs: 1175, groups:  gr.session, 147; group, 16

# Fixed effects:
                  # Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -0.2627     0.3140  -0.837 0.402729    
# z.scan             -0.5578     0.1307  -4.267 1.98e-05 ***
# specieschimpanzee  -0.0940     0.3365  -0.279 0.779963    
# settingzoo         -1.2690     0.3357  -3.780 0.000157 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#plot MAIN effect of SETTING

library(ggplot2)
xdata=droplevels(xdata)

str(xdata)
head(xdata)
xdata$proportion=xdata$ids_nr/xdata$group_size

ggplot(xdata, aes(x=setting, y= proportion, color=setting))+
	geom_point()+
	geom_boxplot()+
	scale_y_continuous(breaks = seq(0, 0.85, 0.2), limits=c(0,0.85))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	#xlab("Time point (15s scans)")+
 	#scale_x_continuous(breaks = seq(1, 8, 1))+
 	#geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T, aes(group=T))+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)))
 	
ggsave("Setting main effect.pdf", width=5, height=5)

tapply(xdata$ids_prop, xdata$setting, mean)
tapply(xdata$ids_prop, xdata$setting, sd)

#plot interaction SETTING*SCAN

library(ggplot2)
xdata=droplevels(xdata)

str(xdata)
head(xdata)
xdata$proportion=xdata$ids_nr/xdata$group_size

#SCAN*SETTING

ggplot(xdata, aes(x=scan, y= proportion, group=scan, color=setting))+
	#geom_point()+
	geom_boxplot()+
	scale_y_continuous(breaks = seq(0, 0.85, 0.2), limits=c(0,0.85))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Time point (15s scans)")+
 	scale_x_continuous(breaks = seq(1, 8, 1))+
 	geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T, aes(group=T))+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)))+
 	facet_wrap(~ setting)
 	
ggsave("Scan setting interaction_SM.pdf", width=7, height=7)

#continue with backward omissions

red3=glmer(cbind(ids_nr, ids_not)~z.scan*setting+species+(1+z.scan||gr.session)+(1+z.scan||group), data=xdata, family=binomial, control=contr)

singles=drop1(red3, test="Chisq")

# species           1 4248.4 0.0773  0.78099  
# z.scan:setting    1 4253.3 5.0029  0.02531 *

#no species differences

ggplot(xdata, aes(x=species, y= proportion, color=species))+
	geom_point()+
	geom_boxplot()+
	scale_y_continuous(breaks = seq(0, 1, 0.2), limits=c(0,1))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	#xlab("Time point (15s scans)")+
 	#scale_x_continuous(breaks = seq(1, 8, 1))+
 	#geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T, aes(group=T))+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)))
 	
ggsave("No species differences.pdf", width=5, height=5)

summary(red3)

# Random effects:
 # Groups       Name        Variance Std.Dev.
 # gr.session   (Intercept) 0.08854  0.2976  
 # gr.session.1 z.scan      0.01322  0.1150  
 # group        (Intercept) 0.42143  0.6492  
 # group.1      z.scan      0.18345  0.4283  
# Number of obs: 1175, groups:  gr.session, 147; group, 16

# Fixed effects:
                  # Estimate Std. Error z value Pr(>|z|)    
# (Intercept)       -0.26274    0.31402  -0.837 0.402766    
# z.scan            -0.26018    0.16383  -1.588 0.112248    
# settingzoo        -1.27091    0.33579  -3.785 0.000154 ***
# specieschimpanzee -0.09375    0.33649  -0.279 0.780532    
# z.scan:settingzoo -0.53820    0.22326  -2.411 0.015925 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Correlation of Fixed Effects:
            # (Intr) z.scan sttngz spcsch
# z.scan       0.000                     
# settingzoo  -0.592  0.000              
# specschmpnz -0.611  0.000  0.010       
# z.scn:sttng -0.001 -0.734  0.019  0.001


tapply(xdata$ids_prop, xdata$species, mean)
tapply(xdata$ids_prop, xdata$species, sd)

> tapply(xdata$ids_prop, xdata$species, mean)
    bonobo chimpanzee 
 0.3381437  0.3384791 
> tapply(xdata$ids_prop, xdata$species, sd)
    bonobo chimpanzee 
 0.1666205  0.2255431 

#check for overdispersion
#function from Ben Bolker (also checked with function from Roger Mundry)
overdisp_fun <- function(model) {
    rdf <- df.residual(model)
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(full)

       # chisq        ratio          rdf            p 
 # 701.8624537    0.6034931 1163.0000000    1.0000000 
 
#Ben bolker: "mild underdispersion is sometimes ignored, since it tends in general to lead to conservative rather than anti-conservative results"


#MODEL 2 -- GROUP EFFECTS

full=glmer(cbind(ids_nr, ids_not)~z.scan*group+(1+z.scan||gr.session)+(1|setting), data=xdata, family=binomial, control=contr)
null=glmer(cbind(ids_nr, ids_not)~1+(1+z.scan||gr.session)+(1|setting), data=xdata, family=binomial, control=contr)

summary(full)

#cbind(fixef(full), confint(full))

anova(full, null, test="Chisq")

     # npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# null    4 4602.7 4623.0 -2297.4   4594.7                         
# full   35 4165.8 4343.2 -2047.9   4095.8 498.94 31  < 2.2e-16 ***

singles=drop1(full, test="Chisq")

#z.scan:group   15 4362.2 226.42 < 2.2e-16 ***

red=glmer(cbind(ids_nr, ids_not)~z.scan+group+(1+z.scan||gr.session)+(1|setting), data=xdata, family=binomial, control=contr)

drop1(red, test="Chisq")

# cbind(ids_nr, ids_not) ~ z.scan + group + (1 + z.scan || gr.session) + 
    # (1 | setting)
       # npar    AIC     LRT   Pr(Chi)    
# <none>      4362.2                      
# z.scan    1 4458.9  98.657 < 2.2e-16 ***
# group    15 4506.8 174.592 < 2.2e-16 ***

summary(red)

# Random effects:
 # Groups       Name        Variance Std.Dev.
 # gr.session   (Intercept) 0.07851  0.2802  
 # gr.session.1 z.scan      0.16007  0.4001  
 # setting      (Intercept) 0.00000  0.0000  
# Number of obs: 1175, groups:  gr.session, 147; setting, 2

# Fixed effects:
                   # Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -0.63982    0.09357  -6.838 8.04e-12 ***
# z.scan             -0.44433    0.03975 -11.178  < 2e-16 ***
# groupChimp_C2       0.95787    0.13056   7.336 2.19e-13 ***
# groupChimp_C3       0.20798    0.16000   1.300  0.19363    
# groupChimp_C4       0.70702    0.13844   5.107 3.27e-07 ***
# groupChimpA_LPZ    -0.86251    0.16224  -5.316 1.06e-07 ***
# groupChimp_Antwerp -1.31219    0.18424  -7.122 1.06e-12 ***
# groupChimp_BB1     -0.19033    0.15863  -1.200  0.23021    
# groupChimp_BB2     -1.20776    0.17785  -6.791 1.12e-11 ***
# groupBonobo_LPZ    -2.14707    0.25180  -8.527  < 2e-16 ***
# groupChimpB_LPZ    -1.41858    0.21471  -6.607 3.92e-11 ***
# groupBonobo_F1     -0.32689    0.21418  -1.526  0.12695    
# groupBonobo_F2      0.14208    0.17334   0.820  0.41238    
# groupBonobo_Lola1  -0.31097    0.13346  -2.330  0.01980 *  
# groupBonobo_Lola2   0.30607    0.13449   2.276  0.02286 *  
# groupBonobo_Lola3   0.37029    0.13348   2.774  0.00554 ** 
# groupBonobo_PLD    -0.21089    0.16652  -1.266  0.20534    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#GROUP COMPARISONS

xdata$proportion=xdata$ids_nr/xdata$group_size
xdata=droplevels(xdata)
tapply(xdata$proportion, xdata$group, mean)

     # Chimp_C1      Chimp_C2      Chimp_C3      Chimp_C4    ChimpA_LPZ Chimp_Antwerp 
    # 0.3532197     0.5762987     0.4000000     0.5147569     0.2279412     0.1500000 
    # Chimp_BB1     Chimp_BB2    Bonobo_LPZ    ChimpB_LPZ     Bonobo_F1     Bonobo_F2 
    # 0.3161376     0.1770833     0.0787037     0.1458333     0.2833333     0.3828125 
 # Bonobo_Lola1  Bonobo_Lola2  Bonobo_Lola3    Bonobo_PLD 
    # 0.2832341     0.4210754     0.4381868     0.3177083  
tapply(xdata$proportion, xdata$group, sd)
     # Chimp_C1      Chimp_C2      Chimp_C3      Chimp_C4    ChimpA_LPZ Chimp_Antwerp 
   # 0.13620073    0.12152105    0.16996732    0.20050914    0.21183213    0.13333333 
    # Chimp_BB1     Chimp_BB2    Bonobo_LPZ    ChimpB_LPZ     Bonobo_F1     Bonobo_F2 
   # 0.12081593    0.20493256    0.12957333    0.18425693    0.13182912    0.14450453 
 # Bonobo_Lola1  Bonobo_Lola2  Bonobo_Lola3    Bonobo_PLD 
   # 0.06990686    0.08725916    0.15194978    0.19282362 
  
library(emmeans)
emmeans(red, pairwise ~ group, adjust="tukey")

 # Chimp_C1 sanctuary - Chimp_C2 sanctuary          -0.9579 0.131 Inf  -7.336  <.0001
 # Chimp_C1 sanctuary - Chimp_C3 sanctuary          -0.2080 0.160 Inf  -1.300  0.9960
 # Chimp_C1 sanctuary - Chimp_C4 sanctuary          -0.7070 0.138 Inf  -5.107  <.0001
 # Chimp_C1 sanctuary - Bonobo_Lola1 sanctuary       0.3110 0.133 Inf   2.330  0.6010
 # Chimp_C1 sanctuary - Bonobo_Lola2 sanctuary      -0.3061 0.134 Inf  -2.276  0.6420
 # Chimp_C1 sanctuary - Bonobo_Lola3 sanctuary      -0.3703 0.133 Inf  -2.774  0.2840
 # Chimp_C1 sanctuary - ChimpA_LPZ zoo               0.8625 0.162 Inf   5.316  <.0001
 # Chimp_C1 sanctuary - Chimp_Antwerp zoo            1.3122 0.184 Inf   7.122  <.0001
 # Chimp_C1 sanctuary - Chimp_BB1 zoo                0.1903 0.159 Inf   1.200  0.9983
 # Chimp_C1 sanctuary - Chimp_BB2 zoo                1.2078 0.178 Inf   6.791  <.0001
 # Chimp_C1 sanctuary - ChimpB_LPZ zoo               1.4186 0.215 Inf   6.607  <.0001
 # Chimp_C1 sanctuary - Bonobo_LPZ zoo               2.1471 0.252 Inf   8.527  <.0001
 # Chimp_C1 sanctuary - Bonobo_F1 zoo                0.3269 0.214 Inf   1.526  0.9797
 # Chimp_C1 sanctuary - Bonobo_F2 zoo               -0.1421 0.173 Inf  -0.820  1.0000
 # Chimp_C1 sanctuary - Bonobo_PLD zoo               0.2109 0.167 Inf   1.266  0.9970
 # Chimp_C2 sanctuary - Chimp_C3 sanctuary           0.7499 0.159 Inf   4.729  0.0003
 # Chimp_C2 sanctuary - Chimp_C4 sanctuary           0.2508 0.137 Inf   1.834  0.9042
 # Chimp_C2 sanctuary - Bonobo_Lola1 sanctuary       1.2688 0.132 Inf   9.633  <.0001
 # Chimp_C2 sanctuary - Bonobo_Lola2 sanctuary       0.6518 0.133 Inf   4.910  0.0001
 # Chimp_C2 sanctuary - Bonobo_Lola3 sanctuary       0.5876 0.132 Inf   4.460  0.0009
 # Chimp_C2 sanctuary - ChimpA_LPZ zoo               1.8204 0.161 Inf  11.311  <.0001
 # Chimp_C2 sanctuary - Chimp_Antwerp zoo            2.2701 0.183 Inf  12.399  <.0001
 # Chimp_C2 sanctuary - Chimp_BB1 zoo                1.1482 0.157 Inf   7.304  <.0001
 # Chimp_C2 sanctuary - Chimp_BB2 zoo                2.1656 0.177 Inf  12.253  <.0001
 # Chimp_C2 sanctuary - ChimpB_LPZ zoo               2.3764 0.214 Inf  11.117  <.0001
 # Chimp_C2 sanctuary - Bonobo_LPZ zoo               3.1049 0.251 Inf  12.369  <.0001
 # Chimp_C2 sanctuary - Bonobo_F1 zoo                1.2848 0.213 Inf   6.028  <.0001
 # Chimp_C2 sanctuary - Bonobo_F2 zoo                0.8158 0.172 Inf   4.743  0.0002
 # Chimp_C2 sanctuary - Bonobo_PLD zoo               1.1688 0.165 Inf   7.077  <.0001
 # Chimp_C3 sanctuary - Chimp_C4 sanctuary          -0.4990 0.165 Inf  -3.023  0.1593
 # Chimp_C3 sanctuary - Bonobo_Lola1 sanctuary       0.5190 0.161 Inf   3.224  0.0922
 # Chimp_C3 sanctuary - Bonobo_Lola2 sanctuary      -0.0981 0.162 Inf  -0.606  1.0000
 # Chimp_C3 sanctuary - Bonobo_Lola3 sanctuary      -0.1623 0.161 Inf  -1.008  0.9998
 # Chimp_C3 sanctuary - ChimpA_LPZ zoo               1.0705 0.186 Inf   5.769  <.0001
 # Chimp_C3 sanctuary - Chimp_Antwerp zoo            1.5202 0.205 Inf   7.413  <.0001
 # Chimp_C3 sanctuary - Chimp_BB1 zoo                0.3983 0.182 Inf   2.184  0.7087
 # Chimp_C3 sanctuary - Chimp_BB2 zoo                1.4157 0.199 Inf   7.101  <.0001
 # Chimp_C3 sanctuary - ChimpB_LPZ zoo               1.6266 0.233 Inf   6.986  <.0001
 # Chimp_C3 sanctuary - Bonobo_LPZ zoo               2.3551 0.267 Inf   8.806  <.0001
 # Chimp_C3 sanctuary - Bonobo_F1 zoo                0.5349 0.232 Inf   2.302  0.6220
 # Chimp_C3 sanctuary - Bonobo_F2 zoo                0.0659 0.195 Inf   0.337  1.0000
 # Chimp_C3 sanctuary - Bonobo_PLD zoo               0.4189 0.189 Inf   2.213  0.6881
 # Chimp_C4 sanctuary - Bonobo_Lola1 sanctuary       1.0180 0.140 Inf   7.296  <.0001
 # Chimp_C4 sanctuary - Bonobo_Lola2 sanctuary       0.4010 0.140 Inf   2.854  0.2387
 # Chimp_C4 sanctuary - Bonobo_Lola3 sanctuary       0.3367 0.140 Inf   2.413  0.5376
 # Chimp_C4 sanctuary - ChimpA_LPZ zoo               1.5695 0.167 Inf   9.375  <.0001
 # Chimp_C4 sanctuary - Chimp_Antwerp zoo            2.0192 0.189 Inf  10.695  <.0001
 # Chimp_C4 sanctuary - Chimp_BB1 zoo                0.8974 0.164 Inf   5.478  <.0001
 # Chimp_C4 sanctuary - Chimp_BB2 zoo                1.9148 0.183 Inf  10.482  <.0001
 # Chimp_C4 sanctuary - ChimpB_LPZ zoo               2.1256 0.219 Inf   9.719  <.0001
 # Chimp_C4 sanctuary - Bonobo_LPZ zoo               2.8541 0.255 Inf  11.182  <.0001
 # Chimp_C4 sanctuary - Bonobo_F1 zoo                1.0339 0.218 Inf   4.742  0.0002
 # Chimp_C4 sanctuary - Bonobo_F2 zoo                0.5649 0.178 Inf   3.173  0.1067
 # Chimp_C4 sanctuary - Bonobo_PLD zoo               0.9179 0.171 Inf   5.354  <.0001
 # Bonobo_Lola1 sanctuary - Bonobo_Lola2 sanctuary  -0.6170 0.136 Inf  -4.550  0.0006
 # Bonobo_Lola1 sanctuary - Bonobo_Lola3 sanctuary  -0.6813 0.135 Inf  -5.061  <.0001
 # Bonobo_Lola1 sanctuary - ChimpA_LPZ zoo           0.5515 0.163 Inf   3.378  0.0582
 # Bonobo_Lola1 sanctuary - Chimp_Antwerp zoo        1.0012 0.185 Inf   5.408  <.0001
 # Bonobo_Lola1 sanctuary - Chimp_BB1 zoo           -0.1206 0.160 Inf  -0.756  1.0000
 # Bonobo_Lola1 sanctuary - Chimp_BB2 zoo            0.8968 0.179 Inf   5.014  0.0001
 # Bonobo_Lola1 sanctuary - ChimpB_LPZ zoo           1.1076 0.216 Inf   5.139  <.0001
 # Bonobo_Lola1 sanctuary - Bonobo_LPZ zoo           1.8361 0.253 Inf   7.271  <.0001
 # Bonobo_Lola1 sanctuary - Bonobo_F1 zoo            0.0159 0.215 Inf   0.074  1.0000
 # Bonobo_Lola1 sanctuary - Bonobo_F2 zoo           -0.4531 0.174 Inf  -2.600  0.3980
 # Bonobo_Lola1 sanctuary - Bonobo_PLD zoo          -0.1001 0.167 Inf  -0.598  1.0000
 # Bonobo_Lola2 sanctuary - Bonobo_Lola3 sanctuary  -0.0642 0.136 Inf  -0.473  1.0000
 # Bonobo_Lola2 sanctuary - ChimpA_LPZ zoo           1.1686 0.164 Inf   7.119  <.0001
 # Bonobo_Lola2 sanctuary - Chimp_Antwerp zoo        1.6183 0.186 Inf   8.705  <.0001
 # Bonobo_Lola2 sanctuary - Chimp_BB1 zoo            0.4964 0.160 Inf   3.093  0.1325
 # Bonobo_Lola2 sanctuary - Chimp_BB2 zoo            1.5138 0.180 Inf   8.426  <.0001
 # Bonobo_Lola2 sanctuary - ChimpB_LPZ zoo           1.7246 0.216 Inf   7.977  <.0001
 # Bonobo_Lola2 sanctuary - Bonobo_LPZ zoo           2.4531 0.253 Inf   9.693  <.0001
 # Bonobo_Lola2 sanctuary - Bonobo_F1 zoo            0.6330 0.216 Inf   2.937  0.1970
 # Bonobo_Lola2 sanctuary - Bonobo_F2 zoo            0.1640 0.175 Inf   0.937  0.9999
 # Bonobo_Lola2 sanctuary - Bonobo_PLD zoo           0.5170 0.168 Inf   3.072  0.1402
 # Bonobo_Lola3 sanctuary - ChimpA_LPZ zoo           1.2328 0.163 Inf   7.551  <.0001
 # Bonobo_Lola3 sanctuary - Chimp_Antwerp zoo        1.6825 0.185 Inf   9.089  <.0001
 # Bonobo_Lola3 sanctuary - Chimp_BB1 zoo            0.5606 0.160 Inf   3.512  0.0378
 # Bonobo_Lola3 sanctuary - Chimp_BB2 zoo            1.5781 0.179 Inf   8.825  <.0001
 # Bonobo_Lola3 sanctuary - ChimpB_LPZ zoo           1.7889 0.215 Inf   8.302  <.0001
 # Bonobo_Lola3 sanctuary - Bonobo_LPZ zoo           2.5174 0.252 Inf   9.971  <.0001
 # Bonobo_Lola3 sanctuary - Bonobo_F1 zoo            0.6972 0.215 Inf   3.244  0.0871
 # Bonobo_Lola3 sanctuary - Bonobo_F2 zoo            0.2282 0.174 Inf   1.310  0.9956
 # Bonobo_Lola3 sanctuary - Bonobo_PLD zoo           0.5812 0.167 Inf   3.471  0.0433
 # ChimpA_LPZ zoo - Chimp_Antwerp zoo                0.4497 0.206 Inf   2.187  0.7063
 # ChimpA_LPZ zoo - Chimp_BB1 zoo                   -0.6722 0.184 Inf  -3.647  0.0238
 # ChimpA_LPZ zoo - Chimp_BB2 zoo                    0.3452 0.199 Inf   1.733  0.9382
 # ChimpA_LPZ zoo - ChimpB_LPZ zoo                   0.5561 0.232 Inf   2.395  0.5515
 # ChimpA_LPZ zoo - Bonobo_LPZ zoo                   1.2846 0.267 Inf   4.818  0.0002
 # ChimpA_LPZ zoo - Bonobo_F1 zoo                   -0.5356 0.234 Inf  -2.290  0.6311
 # ChimpA_LPZ zoo - Bonobo_F2 zoo                   -1.0046 0.197 Inf  -5.095  <.0001
 # ChimpA_LPZ zoo - Bonobo_PLD zoo                  -0.6516 0.191 Inf  -3.414  0.0519
 # Chimp_Antwerp zoo - Chimp_BB1 zoo                -1.1219 0.204 Inf  -5.502  <.0001
 # Chimp_Antwerp zoo - Chimp_BB2 zoo                -0.1044 0.218 Inf  -0.480  1.0000
 # Chimp_Antwerp zoo - ChimpB_LPZ zoo                0.1064 0.248 Inf   0.429  1.0000
 # Chimp_Antwerp zoo - Bonobo_LPZ zoo                0.8349 0.281 Inf   2.975  0.1793
 # Chimp_Antwerp zoo - Bonobo_F1 zoo                -0.9853 0.250 Inf  -3.947  0.0078
 # Chimp_Antwerp zoo - Bonobo_F2 zoo                -1.4543 0.216 Inf  -6.745  <.0001
 # Chimp_Antwerp zoo - Bonobo_PLD zoo               -1.1013 0.210 Inf  -5.247  <.0001
 # Chimp_BB1 zoo - Chimp_BB2 zoo                     1.0174 0.198 Inf   5.135  <.0001
 # Chimp_BB1 zoo - ChimpB_LPZ zoo                    1.2282 0.232 Inf   5.299  <.0001
 # Chimp_BB1 zoo - Bonobo_LPZ zoo                    1.9567 0.266 Inf   7.343  <.0001
 # Chimp_BB1 zoo - Bonobo_F1 zoo                     0.1366 0.231 Inf   0.590  1.0000
 # Chimp_BB1 zoo - Bonobo_F2 zoo                    -0.3324 0.194 Inf  -1.712  0.9442
 # Chimp_BB1 zoo - Bonobo_PLD zoo                    0.0206 0.188 Inf   0.109  1.0000
 # Chimp_BB2 zoo - ChimpB_LPZ zoo                    0.2108 0.242 Inf   0.870  1.0000
 # Chimp_BB2 zoo - Bonobo_LPZ zoo                    0.9393 0.275 Inf   3.410  0.0526
 # Chimp_BB2 zoo - Bonobo_F1 zoo                    -0.8809 0.245 Inf  -3.596  0.0285
 # Chimp_BB2 zoo - Bonobo_F2 zoo                    -1.3499 0.210 Inf  -6.421  <.0001
 # Chimp_BB2 zoo - Bonobo_PLD zoo                   -0.9969 0.204 Inf  -4.882  0.0001
 # ChimpB_LPZ zoo - Bonobo_LPZ zoo                   0.7285 0.300 Inf   2.431  0.5239
 # ChimpB_LPZ zoo - Bonobo_F1 zoo                   -1.0917 0.273 Inf  -4.001  0.0064
 # ChimpB_LPZ zoo - Bonobo_F2 zoo                   -1.5607 0.242 Inf  -6.444  <.0001
 # ChimpB_LPZ zoo - Bonobo_PLD zoo                  -1.2077 0.237 Inf  -5.099  <.0001
 # Bonobo_LPZ zoo - Bonobo_F1 zoo                   -1.8202 0.303 Inf  -6.009  <.0001
 # Bonobo_LPZ zoo - Bonobo_F2 zoo                   -2.2892 0.276 Inf  -8.306  <.0001
 # Bonobo_LPZ zoo - Bonobo_PLD zoo                  -1.9362 0.271 Inf  -7.148  <.0001
 # Bonobo_F1 zoo - Bonobo_F2 zoo                    -0.4690 0.242 Inf  -1.940  0.8572
 # Bonobo_F1 zoo - Bonobo_PLD zoo                   -0.1160 0.237 Inf  -0.490  1.0000
 # Bonobo_F2 zoo - Bonobo_PLD zoo                    0.3530 0.201 Inf   1.759  0.9306

# Results are given on the log odds ratio (not the response) scale. 
# P value adjustment: tukey method for comparing a family of 16 estimates 

#split emmeans comparison between sanctuary-housed and zoo-housed groups because of the effect of "setting" (usage of different experimental assays)

#PLOT

library(ggplot2)
xdata=droplevels(xdata)

#all groups together

levels(xdata$group)
xdata$proportion=xdata$ids_nr/xdata$group_size

xdata$group=factor(xdata$group, levels=c("Chimp_C1","Chimp_C2","Chimp_C3", "Chimp_C4", "ChimpA_LPZ", "Chimp_Antwerp", "Chimp_BB1", "Chimp_BB2", "ChimpB_LPZ", "Bonobo_LPZ", "Bonobo_F1", "Bonobo_F2", "Bonobo_Lola1","Bonobo_Lola2","Bonobo_Lola3","Bonobo_PLD"))

ggplot(xdata, aes(x=group, y= proportion, color=species))+
	geom_point()+
	geom_boxplot(width=0.5)+
	scale_y_continuous(breaks = seq(0, 1, 0.2), limits=c(0,1))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Group")+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 14, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)), axis.text.x=element_text(angle=90,hjust=1))+
 	facet_wrap(~setting, ncol = 1, strip.position = "right", scales = "free_x")

ggsave("Co-feeding totals per group_facet wrap.pdf", width=9, height=9)

#slopes per species

ggplot(xdata, aes(x=scan, y= proportion, group=scan, color=species))+
	#geom_point()+
	geom_boxplot()+
	scale_y_continuous(breaks = seq(0, 0.85, 0.2), limits=c(0,0.85))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Time point (15s scans)")+
 	scale_x_continuous(breaks = seq(1, 8, 1))+
 	geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T, aes(group=T))+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)))+
 	facet_wrap(~ species)

ggsave("Cofeeding per species_slopes.pdf", width=7, height=7)

#slopes (original figures)

ggplot(xdata, aes(x=scan, y= proportion, group=scan, color=species))+
	geom_point()+
	geom_boxplot()+
	stat_sum(aes(size=factor(..n..)), geom="point", show.legend=FALSE)+
	scale_size_discrete(range = c(0.05, 0.9))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Time point (15s scans)")+
 	scale_x_continuous(breaks = seq(1, 8, 1))+
 	scale_y_continuous(breaks = seq(0, 0.8, 0.2))+
 	geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x + I(x^2), 	se=T, aes(group=T))+
 	scale_colour_pander() +
 	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)))+
 	facet_wrap(~ group)

ggsave("Cofeeding per group_slopes.pdf", width=7, height=7)


#MODEL 3 -- GROUP-LEVEL INDICES EFFECTS ON COFEEDING TOLERANCE

rm(list=ls())
setwd("~/Desktop/Co-feeding tolerance MS_PRSB/Analyses")
ydata=read.table("van_Leeuwen_et_al._dataset_co-feeding_tolerance_group indices.txt", header=T, sep="\\t")
str(ydata)

ydata$z.group.size=as.vector(scale(ydata$group_size))
ydata$z.male.ratio=as.vector(scale(ydata$male.ratio))
ydata$z.mean.age=as.vector(scale(ydata$mean.age))
ydata$z.modularity =as.vector(scale(ydata$modularity))
ydata$z.degree=as.vector(scale(ydata$degree))
ydata$z.prop.kin.dyads =as.vector(scale(ydata$prop.kin.dyads))

#omit one observation

sum(is.na(ydata$mean.prop))
nrow(ydata)
ydata=subset(ydata, !is.na(ydata$mean.prop))

sort(unique(ydata$group))
ydata$setting=ifelse(ydata$group=="Bonobo_F1" |ydata$group=="Bonobo_F2" |ydata$group=="Bonobo_LPZ" |ydata$group=="Bonobo_PLD" |ydata$group=="Chimp_Antwerp" |ydata$group=="Chimp_BB1" |ydata$group=="Chimp_BB2" |ydata$group=="ChimpA_LPZ" |ydata$group=="ChimpB_LPZ", c("zoo"), c("sanctuary"))

library(lme4)
contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))

str(ydata)
ydata=subset(ydata, !is.na(ydata$z.degree))

full=glmer(cbind(ids_nr, ids_not)~z.scan*setting+species*(z.group.size+z.modularity+z.degree+ z.prop.kin.dyads+ z.mean.age+ z.male.ratio)+(1+z.scan||gr.session)+(1|group), data=ydata, family=binomial, control=contr)

null=glmer(cbind(ids_nr, ids_not)~z.scan*setting+(1+z.scan||gr.session)+(1|group), data=ydata, family=binomial, control=contr)

anova(null, full, test="Chisq")

     # npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# null    7 3965.7 4000.5 -1975.8   3951.7                        
# full   20 3957.2 4056.7 -1958.6   3917.2 34.465 13   0.001022 **

confint(full)

# (Intercept)                        -5.85429132 -1.226168290
# z.scan                             -0.30201678 -0.120836823
# settingzoo                         -1.89834053  2.040222154
# specieschimpanzee                   1.38084848  4.998597696
# z.group.size                       -1.56504843  0.215749763
# z.modularity                              -Inf -1.259149046
# z.degree                           -0.69061005  0.353549938
# z.prop.kin.dyads                   -1.80924572  0.003241254
# z.mean.age                         -6.88587870 -1.058784312
# z.male.ratio                       -4.69465805 -1.258926881
# z.scan:settingzoo                  -0.67578470 -0.385450267
# specieschimpanzee:z.group.size     -0.07876191  1.817313189
# specieschimpanzee:z.modularity      1.45169350  4.850312009
# specieschimpanzee:z.degree         -0.78047898  0.524474625
# specieschimpanzee:z.prop.kin.dyads -0.11243785  2.426390440
# specieschimpanzee:z.mean.age        0.95717565  4.337203837
# specieschimpanzee:z.male.ratio      1.40703328  4.834346206


#from Ben Bolker
overdisp_fun <- function(model) {
    rdf <- df.residual(model)
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(full)#0.5932728

drop1(full, test="Chisq")

# z.scan:setting              1 4000.4 45.201 1.778e-11 ***
# species:z.group.size        1 3958.4  3.158  0.075543 .  
# species:z.modularity        1 3962.7  7.468  0.006282 ** 
# species:z.degree            1 3955.4  0.152  0.696662    
# species:z.prop.kin.dyads    1 3957.7  2.534  0.111398    
# species:z.mean.age          1 3964.3  9.073  0.002595 ** 
# species:z.male.ratio        1 3965.7 10.535  0.001171 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(full)

#MODULARITY

str(ydata)
range(ydata$modularity)
ydata$proportion=ydata$ids_nr/ydata$group_size

ggplot(ydata, aes(x= modularity, y= proportion, group=species, color=species))+
	geom_point()+
	geom_count(show.legend = F)+
	#geom_boxplot()+
	scale_y_continuous(breaks = seq(0, 0.9, 0.2), limits=c(0,0.9))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Modularity")+
 	#scale_x_continuous(breaks = seq(0, 0.85, 0.2), limits=c(0,0.85))+
 	geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T)+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)), axis.title.x = element_text(margin = margin(t=10, r=10, b=20, l=15)))
 	
ggsave("Modularity species interaction_sig.pdf", width=7, height=7)


#MEAN AGE
str(ydata)

ggplot(ydata, aes(x= mean.age, y= proportion, group=species, color=species))+
	geom_point()+
	geom_count(show.legend = F)+
	#geom_boxplot()+
	scale_y_continuous(breaks = seq(0, 0.9, 0.2), limits=c(0,0.9))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Age (group average)")+
 	scale_x_continuous(breaks = seq(10, 35, 5), limits=c(10,35))+
 	geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T)+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)), axis.title.x = element_text(margin = margin(t=10, r=10, b=20, l=15)))
 	
ggsave("Mean age species interaction_sig_better scaled x-axis.pdf", width=7, height=7)


#MALE RATIO

ggplot(ydata, aes(x= male.ratio, y= proportion, group=species, color=species))+
	geom_point()+
	geom_count(show.legend = F)+
	#geom_boxplot()+
	scale_y_continuous(breaks = seq(0, 0.9, 0.2), limits=c(0,0.9))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Sex ratio (proportion males)")+
 	#scale_x_continuous(breaks = seq(0, 0.85, 0.2), limits=c(0,0.85))+
 	geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T)+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)), axis.title.x = element_text(margin = margin(t=10, r=10, b=20, l=15)))
 	
ggsave("Male ratio species interaction_sig.pdf", width=7, height=7)


#main effects

red=glmer(cbind(ids_nr, ids_not)~z.scan*setting+z.group.size+z.degree+z.prop.kin.dyads+species*(z.modularity+ z.mean.age+ z.male.ratio)+(1+z.scan||gr.session)+(1|group), data=ydata, family=binomial, control=contr)

drop1(red, test="Chisq")

summary(red)

# z.group.size            1 3957.0  3.107  0.077939 .  
# z.degree                1 3957.5  3.559  0.059228 .  
# z.prop.kin.dyads        1 3956.8  2.910  0.088023 .  
# z.scan:setting          1 3998.9 44.963 2.008e-11 ***
# species:z.modularity    1 3968.4 14.489  0.000141 ***
# species:z.mean.age      1 3965.4 11.543  0.000680 ***
# species:z.male.ratio    1 3970.3 16.377 5.191e-05 ***

# Number of obs: 1071, groups:  gr.session, 134; group, 15

# Fixed effects:
                               # Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                    -2.41008    0.55247  -4.362 1.29e-05 ***
# z.scan                         -0.21069    0.04554  -4.627 3.71e-06 ***
# settingzoo                     -0.70689    0.40647  -1.739 0.082021 .  
# z.group.size                    0.19672    0.09919   1.983 0.047328 *  
# z.degree                       -0.23731    0.11936  -1.988 0.046788 *  
# z.prop.kin.dyads               -0.28254    0.16071  -1.758 0.078734 .  
# specieschimpanzee               2.03417    0.46411   4.383 1.17e-05 ***
# z.modularity                   -2.20905    0.39450  -5.600 2.15e-08 ***
# z.mean.age                     -1.90049    0.53514  -3.551 0.000383 ***
# z.male.ratio                   -1.95527    0.43074  -4.539 5.64e-06 ***
# z.scan:settingzoo              -0.52654    0.07306  -7.207 5.71e-13 ***
# specieschimpanzee:z.modularity  2.05143    0.44155   4.646 3.38e-06 ***
# specieschimpanzee:z.mean.age    1.35119    0.35058   3.854 0.000116 ***
# specieschimpanzee:z.male.ratio  2.01516    0.40300   5.000 5.72e-07 ***

confint(red)

# (Intercept)                    -3.58591038 -1.28610352
# z.scan                         -0.30185921 -0.12095573
# settingzoo                     -1.54855502  0.15298714
# z.group.size                   -0.02383749  0.39464580
# z.degree                       -0.57482534  0.01247654
# z.prop.kin.dyads               -0.62742720  0.05568436
# specieschimpanzee               1.09049357  3.01992665
# z.modularity                   -3.04391148 -1.39966910
# z.mean.age                     -3.03441744 -0.81956677
# z.male.ratio                   -2.86643324 -1.08817843
# z.scan:settingzoo              -0.67320867 -0.38365279
# specieschimpanzee:z.modularity  1.15187358  3.00398014
# specieschimpanzee:z.mean.age    0.64218341  2.09551067
# specieschimpanzee:z.male.ratio  1.20121624  2.87097038

#Main effect TRENDS

#GROUP SIZE

str(ydata)

ggplot(ydata, aes(x= group_size, y= proportion))+
	geom_point()+
	geom_count(show.legend = F)+
	#geom_boxplot()+
	scale_y_continuous(breaks = seq(0, 0.9, 0.2), limits=c(0,0.9))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Group size")+
 	#scale_x_continuous(breaks = seq(0, 0.85, 0.2), limits=c(0,0.85))+
 	geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T)+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)), axis.title.x = element_text(margin = margin(t=10, r=10, b=20, l=15)))
 	
ggsave("Group size main effect n=15.pdf", width=7, height=7)

#without largest group (outlier) - group 2 CWOT

zdata=subset(ydata, ydata$group!="Chimp_C2")
red=glmer(cbind(ids_nr, ids_not)~z.scan*setting+z.group.size+z.degree+z.prop.kin.dyads+species*(z.modularity+ z.mean.age+ z.male.ratio)+(1+z.scan||gr.session)+(1|group), data=zdata, family=binomial, control=contr)

drop1(red, test="Chisq")

#z.group.size            1 3437.4  0.093 0.7602819 

str(ydata)

ggplot(zdata, aes(x= group_size, y= proportion))+
	geom_point()+
	geom_count(show.legend = F)+
	#geom_boxplot()+
	scale_y_continuous(breaks = seq(0, 0.9, 0.2), limits=c(0,0.9))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Group size")+
 	#scale_x_continuous(breaks = seq(0, 0.85, 0.2), limits=c(0,0.85))+
 	geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T)+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)), axis.title.x = element_text(margin = margin(t=10, r=10, b=20, l=15)))
 	
ggsave("Group size main effect n=14.pdf", width=7, height=7)

#SOCIAL BONDING

str(ydata)

ggplot(ydata, aes(x= degree, y= proportion))+
	geom_point()+
	geom_count(show.legend = F)+
	#geom_boxplot()+
	scale_y_continuous(breaks = seq(0, 0.9, 0.2), limits=c(0,0.9))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Bond strength (group average)")+
 	#scale_x_continuous(breaks = seq(0, 0.85, 0.2), limits=c(0,0.85))+
 	geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T)+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 16, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)), axis.title.x = element_text(margin = margin(t=10, r=10, b=20, l=15)))
 	
ggsave("Social bonding main effect.pdf", width=7, height=7)

#RELATEDNESS

library(ggplot2)
library(ggthemes)
ydata=droplevels(ydata)

str(ydata)
head(ydata)
ydata$proportion=ydata$ids_nr/ydata$group_size
ydata$prop.kin.dyads_prop=ydata$prop.kin.dyads/100
range(ydata$prop.kin.dyads_prop)

ydata$preds = predict(red, ydata, type = "response")
library(ggeffects)
library(tidyverse)

str(ydata)

#Simpson's paradox: prop dyads actually positive corr but negative estimate due to other variables
#plot marginal effect

df=ggpredict(red, terms="z.prop.kin.dyads")

ggplot(df, aes(x, predicted))+
	geom_point()+
	geom_count(show.legend = F)+
	#geom_boxplot()+
	geom_smooth(method="glm", method.args = list(family = "binomial"), formula = y ~ x, 	se=T, aes(group=T))+
	geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1)+
	scale_y_continuous(breaks = seq(0, 0.25, 0.05), limits=c(0,0.25))+
	#theme(axis.text.x=element_text(angle=0,hjust=1))+
 	ylab("Co-feeding tolerance (in proportion)")+
 	xlab("Maternal kinship (proportion of dyads)")+
 	#scale_x_continuous(breaks = seq(0, 1, 0.1))+
 	scale_colour_pander() +
	scale_fill_pander() +
 	theme_pander(base_size = 15, base_family = "Times")+
 	theme(axis.title.y.left = element_text(margin = margin(0, 10, 20, 15)), axis.title.x = element_text(margin = margin(t=10, r=10, b=20, l=15)))
 	#facet_wrap(~scan)

ggsave("Kinship main effect_marginal effect.pdf", width=7, height=7)


#check for random effect of GROUP

full=glmer(cbind(ids_nr, ids_not)~z.scan*setting+species*(z.group.size+z.modularity+z.degree+ z.prop.kin.dyads+ z.mean.age+ z.male.ratio)+(1+z.scan||gr.session)+(1|group), data=ydata, family=binomial, control=contr)

group=glmer(cbind(ids_nr, ids_not)~z.scan*setting+species*(z.group.size+z.modularity+z.degree+ z.prop.kin.dyads+ z.mean.age+ z.male.ratio)+(1+z.scan||gr.session), data=ydata, family=binomial, control=contr)

anova(full, group, test="Chisq")

      # npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# group   19 3959.9 4054.5 -1961.0   3921.9                       
# full    20 3957.2 4056.7 -1958.6   3917.2 4.6962  1    0.03023 *

#STABILITY

#model stability - jack-knifing

#interaction model

res=glm(cbind(ids_nr, ids_not)~z.scan*setting+species*(z.group.size+z.modularity+z.degree+ z.prop.kin.dyads+ z.mean.age+ z.male.ratio), data=ydata, family=binomial)

summary(res)

subs=levels(as.factor(ydata$group))

all.coeffs=matrix(NA, nrow=length(subs), ncol=length(coefficients(res)))
colnames(all.coeffs)=names(coefficients(res))
for(i in 1:length(subs)){
	ires=glm(cbind(ids_nr, ids_not)~z.scan*setting+species*(z.group.size+z.modularity+z.degree+ z.prop.kin.dyads+ z.mean.age+ z.male.ratio), data=subset(ydata, group!=subs[i]), family=binomial)
	ires=coefficients(ires)
	all.coeffs[i, names(ires)]=ires
}
round(cbind(coefficients(res), t(apply(all.coeffs, 2, range, na.rm=T))), 3)
m.stab=cbind(coefficients(res), t(apply(all.coeffs, 2, range, na.rm=T)))

                                     # [,1]    [,2]   [,3]
# (Intercept)                        -2.663 -50.934  2.096
# z.scan                             -0.174  -0.199 -0.122
# settingzoo                         -0.708  -4.710 40.050
# specieschimpanzee                   2.319  -0.217 40.938
# z.group.size                       -0.606  -0.607 43.864
# z.modularity                       -2.170 -38.204  0.595
# z.degree                           -0.088 -10.840  3.094
# z.prop.kin.dyads                   -0.612 -18.699  5.735
# z.mean.age                         -2.348 -44.511  1.812
# z.male.ratio                       -2.178 -31.164  0.686
# z.scan:settingzoo                  -0.511  -0.627 -0.405
# specieschimpanzee:z.group.size      0.876 -43.593  1.400
# specieschimpanzee:z.modularity      1.938  -2.033 37.972
# specieschimpanzee:z.degree         -0.107 -12.877 10.645
# specieschimpanzee:z.prop.kin.dyads  0.452  -5.896  8.508
# specieschimpanzee:z.mean.age        1.980  -2.094  6.411
# specieschimpanzee:z.male.ratio      2.253  -2.320 18.708

#make GGPLOT

xx=as.data.frame(m.stab)
str(xx)
#xx=read.table("Glmm_modelstab.txt", header=T, sep="\\t")
xx[1:13,]
#xx=xx[1:4,]
colnames(xx)=c("orig", "min", "max")
row.names(xx)

library(ggplot2)
library(plyr)

segment_data = data.frame(
    x = xx$min,
    xend = xx$max, 
    y = unname(row.names(xx)),
    yend = unname(row.names(xx)))

library(ggplot2)
library(scales)

ggplot(xx, aes(x=orig, y= row.names(xx)))+
	geom_point()+
	#geom_point(position = position_jitter(w = 0.3, h = 0), colour="darkgrey", alpha=1/2)+
	#geom_boxplot()+
	theme_bw()+
	theme(text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=10))+
	geom_segment(data=segment_data, aes(x = x, y = y, xend = xend, yend = yend, col="red"), size=1, show.legend = F)+
	geom_vline(xintercept=0, linetype="dashed", color = "blue", size=0.3)+
	#scale_y_continuous(breaks = seq(-5, 15, by = 5), limits=c(-5,15))+
	#scale_x_continuous(breaks = seq(-1.5, 1, by = 0.25), limits=c(-1.5,1))+
 	ylab("")+
 	xlab("Estimate range")
 	
ggsave("Model Stability Group Omission Model 3 interactions n=15_Oct2022.pdf", width = 6, height = 8)

#main effects

res=glm(cbind(ids_nr, ids_not)~z.scan*setting+z.group.size+z.degree+z.prop.kin.dyads+species*(z.modularity+ z.mean.age+ z.male.ratio), data=ydata, family=binomial)

summary(res)

subs=levels(as.factor(ydata$group))

all.coeffs=matrix(NA, nrow=length(subs), ncol=length(coefficients(res)))
colnames(all.coeffs)=names(coefficients(res))
for(i in 1:length(subs)){
	ires=glm(cbind(ids_nr, ids_not)~z.scan*setting+z.group.size+z.degree+z.prop.kin.dyads+species*(z.modularity+ z.mean.age+ z.male.ratio), data=subset(ydata, group!=subs[i]), family=binomial)
	ires=coefficients(ires)
	all.coeffs[i, names(ires)]=ires
}
round(cbind(coefficients(res), t(apply(all.coeffs, 2, range, na.rm=T))), 3)
m.stab=cbind(coefficients(res), t(apply(all.coeffs, 2, range, na.rm=T)))

                                 # [,1]   [,2]   [,3]
# (Intercept)                    -2.103 -3.628 -1.542
# z.scan                         -0.174 -0.200 -0.122
# settingzoo                     -0.771 -1.282  0.186
# z.group.size                    0.249  0.084  0.330
# z.degree                       -0.215 -0.651  0.046
# z.prop.kin.dyads               -0.217 -0.609 -0.015
# specieschimpanzee               1.799  1.295  3.473
# z.modularity                   -1.984 -3.485 -1.682
# z.mean.age                     -1.573 -3.113 -1.057
# z.male.ratio                   -1.667 -2.919 -1.257
# z.scan:settingzoo              -0.510 -0.628 -0.397
# specieschimpanzee:z.modularity  1.738  1.392  3.260
# specieschimpanzee:z.mean.age    1.174  0.706  1.943
# specieschimpanzee:z.male.ratio  1.708  1.321  3.378


#make GGPLOT

xx=as.data.frame(m.stab)
str(xx)
#xx=read.table("Glmm_modelstab.txt", header=T, sep="\\t")
xx[1:13,]
#xx=xx[1:4,]
colnames(xx)=c("orig", "min", "max")
row.names(xx)

library(ggplot2)
library(plyr)

segment_data = data.frame(
    x = xx$min,
    xend = xx$max, 
    y = unname(row.names(xx)),
    yend = unname(row.names(xx)))

library(ggplot2)
library(scales)

ggplot(xx, aes(x=orig, y= row.names(xx)))+
	geom_point()+
	#geom_point(position = position_jitter(w = 0.3, h = 0), colour="darkgrey", alpha=1/2)+
	#geom_boxplot()+
	theme_bw()+
	theme(text=element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title=element_text(size=10))+
	geom_segment(data=segment_data, aes(x = x, y = y, xend = xend, yend = yend, col="red"), size=1, show.legend = F)+
	geom_vline(xintercept=0, linetype="dashed", color = "blue", size=0.3)+
	#scale_y_continuous(breaks = seq(-5, 15, by = 5), limits=c(-5,15))+
	#scale_x_continuous(breaks = seq(-1.5, 1, by = 0.25), limits=c(-1.5,1))+
 	ylab("")+
 	xlab("Estimate range")
 	
ggsave("Model Stability Group Omission Model 3 main effects n=15_Oct2022_new.pdf", width = 6, height = 8)
