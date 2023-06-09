

################################################################################
################                                                ################
################    MODELS AND FIGURES FOR FISCHER ET. AL.      ###############
################   "FEMALE FITNESS COSTS LINKED TO RESOURCE     ###############
################    DEFENCE AND RELATEDNESS OF COMPETITORS"     ###############
################                                                ################
################################################################################


############
## LOAD THE REQUIRED LIBRARIES

library(lme4)
library(ggplot2)
library(emmeans)
library(plyr) 
library(nortest) 
library(car) 
library(afex)



###############################################################################
##                                                                           ##
##                    IMPORTING AND PREPARING THE DATAFILES                  ## 
##                                                                           ##
###############################################################################


############
## IMPORTING THE DATASET RELATED TO ALL ANALYSIS INCLUDING FITNESS COSTS

maf.d1<-read.csv("Fischer_FitnessCosts_Data.csv", header=T)
str(maf.d1)
maf.d1$nr.comm.pups <- as.numeric(maf.d1$nr.comm.pups)

############
## IMPORTING THE DATASET RELATED TO ALL ANALYSIS INCLUDING THE SCENT MARKING 
## ASSAY 

ccm.c.ur.re<-read.csv("Fischer_ScentMarking_Data.csv", header=T)
str(ccm.c.ur.re)

############
## IMPORTING THE DATASET RELATED TO ALL ANALYSIS INCLUDING OFFSPRING MASS

WWP.W<-read.csv("Fischer_WeaningMass_Data.csv", header=T)

############
## IMPORTING THE DATASET RELATED TO THE ANALYSIS INCLUDING SLEEPING NEST SITE
## ASSOCIATIONS

SN <- read.csv("Fischer_Sleep_Data.csv", header=T)
str(SN)

############
## IMPORTING THE DATASET RELATED TO FIGURE S1 "ACTIVITY OF FEMALES DURING PEAK
## MATERNAL INVESTMENT
lall <- read.csv("Fischer_Activity_Data.csv", header=T)

############
## CLEAN THE DATA STRUCTURE

# SUBSET DATASET "maf.d1" TO SHOW DATA FOR COMMUNAL NESTS ONLY
maf.d1.c<-subset(maf.d1,comm.nest>0)
maf.d1.c<-droplevels(maf.d1.c)
str(maf.d1)

# SUBSET DATASET "maf.d1" TO SHOW ONLY DATA FOR FEMALES THAT WEANED OFFSPRING
maf.d1.b<-subset(maf.d1, nr.pups >= 1)      
str(maf.d1.b)

# CHANGE THE CONTRAST FOR THE FACTOR "stimulus_comp" ON THE DATASET "ccm.new" 
# TO CONDUCT THE ORTHOGONAL CONTRASTS

# MAKE A COPY OF THE DATASET "ccm.c.ur.re" CALLED "ccm.new"
ccm.new <- ccm.c.ur.re
str(ccm.new)

# MAKE "stimulus_comp" TO A FACTOR AND CHECK THE CONTRAST SETTINGS
ccm.new$stimulus_comp<-as.factor(ccm.new$stimulus_comp)
contrasts(ccm.new$stimulus_comp)

# CREATE THE MATRIX FOR THE NEW CONTRAST SETTING
newcontrast<-cbind(c(-1, -1, -1, 1, 1, 1), c(-2, 1, 1, 0, 0, 0), 
                    c(0, -1, 1, 0, 0, 0), c(0, 0, 0, -2, 1, 1), 
                    c(0, 0, 0, 0, -1, 1))

# CREATE THE COLUMN NAMES
colnames(newcontrast)<-c("HC->LC", "HC/control->HC/stimulus", 
                         "HC/RE->HC/UR", "LC/control->LC/stimulus", 
                         "LC/RE->LC/UR")

# CREATE THE ROW NAMES
rownames(newcontrast)<-c("HC_contr", "HC_RE", "HC_UR", "LC_contr",
                         "LC_RE", "LC_UR")

# CHANGE THE CONTRAST SETTINGS OF the factor "stimulus_comp" TO ORTHOGONAL 
# CONTRASTS ON THE DATASET "ccm.new"
contrasts(ccm.new$stimulus_comp) = newcontrast

# CREATE AVERAGE WEANING MASSES FOR EACH LITTER AND EACH FEMALE FROM THE
# DATAFILE "WWP.W"
wwp.ALL.w<-ddply(WWP.W, c("Competition", "Relatedness", "sister_ID", 
                        "Age_start", "Sex", "total", "dam_weight_w", "block"), 
                 summarise, 
                 weight = mean(Weight, na.rm = TRUE))  

# CHANGE "treat2" TO A FACTOR ON THE DATASET "SN"
SN$treat2 <- as.factor(SN$treat2)


###############################################################################
##                                                                           ##
##                          MODEL FOR TABLE 1                                ## 
##                     "NUMBER OF WEANED OFFSPRING"                          ##
##                                                                           ##
###############################################################################


############
## MODEL FOR TABLE 1 "NUMBER OF WEANED OFFSPRING"

Pups.m1.b<-lmer(log(nr.pups) ~  comp + rel
                + comm.nest +  (1|Nest) + (1|block), data=maf.d1.b)

# THE FIT OF THIS MODEL IS SINGULAR BECAUSE THE RANDOM FACTOR BLOCK IDENTITY
# ACCOUNTS FOR 0 VARIANCE. WE STILL THINK IT IS IMPORTANT TO INCLUDE THIS 
# RANDOM EFFECT BECAUSE IT STILL ACCOUNTS FOR THE MINOR DIFFERENCES BETWEEN 
# THE BLOCKS

# OBTAIN THE ESTIMATES
summary(Pups.m1.b)

# OBTAIN THE P-VALUES
drop1(Pups.m1.b, test="F")

# REMOVED FACTORS
# +  Weight_W
# + Age_start 

# CHECK THE MODEL
par(mfrow=c(2,2))
plot(resid(Pups.m1.b)~fitted(Pups.m1.b))
abline(h=0, lty=2)
hist(resid(Pups.m1.b))
qqnorm(resid(Pups.m1.b))
qqline(resid(Pups.m1.b))
resid.Pups.m1.b<-resid(Pups.m1.b)
shapiro.test(resid.Pups.m1.b)
lillie.test(resid.Pups.m1.b)


###############################################################################
##                                                                           ##
##                        MODEL AND FIGURE FOR TABLE 2A                      ## 
##                "ACTIVITY OF FEMALES DURING THE LIGHT PHASE"               ##
##                                                                           ##
###############################################################################


############
## CODE FOR FIGURE 2A "ACTIVITY OF FEMALES DURING THE LIGHT PHASE"

# REMOVE NAs FROM THE DATASET FOR THE PLOT
maf.d1_1 <- maf.d1[!is.na(maf.d1$activityL),]

# THE BOXPLOT
AcL.B1<-ggplot(maf.d1_1, aes(x=comp, y=activityL, 
                     fill=rel)) + 
                   geom_boxplot(position=position_dodge(0.9), fatten=3)
AcL.B1 + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill="white"),
              panel.border = element_rect(colour="white", fill=NA, size=1)) + 
              theme(axis.line = element_line(colour="black"))+ 
              xlab("\\nNest site availability") + 
              ylab("Activity during the light phase\\n")+ 
              labs(fill = "Relatedness of\\ncompetitors")+
              scale_fill_discrete(name ="Relatedness of\\ncompetitors", 
              labels = c("Related", "Unrelated")) + 
              scale_x_discrete(labels = c("Single nest site",
              "Multiple nest sites")) + 
               theme(axis.text.x = element_text(size = 12),
                     axis.text.y = element_text(size = 12),  
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14)) + 
              theme(legend.title=element_text(size=14), 
                    legend.text=element_text(size=12)) +
              theme(legend.key=element_blank())

############
## MODEL FOR TABLE 2A "ACTIVITY OF FEMALES DURING THE LIGHT PHASE"

AcL.m1 <- lmer(log(activityL) ~ rel + comp + log(time.in.nest.L) + (1|Trial), 
              data=maf.d1)

# OBTAIN THE ESTIMATES
summary(AcL.m1)

# OBTAIN THE P-VALUES
drop1(AcL.m1, test="F")

# CHECK THE MODEL
par(mfrow=c(2,2))
plot(resid(AcL.m1)~fitted(AcL.m1))
abline(h=0, lty=2)
hist(resid(AcL.m1))
qqnorm(resid(AcL.m1))
qqline(resid(AcL.m1))
resid.AcL.m1<-resid(AcL.m1)
shapiro.test(resid.AcL.m1)
lillie.test(resid.AcL.m1)


###############################################################################
##                                                                           ##
##                       MODEL AND FIGURE FOR TABLE 2B                       ## 
##                             "AREA SCENT MARKED"                           ##
##                                                                           ##
###############################################################################


############
## CODE FOR FIGURE 3 "AREA SCENT MARKED"

###### 
## MAKE A BOXPLOT

b.SM1<-ggplot(ccm.c.ur.re, aes(x=stimulus3, y=sum_area_marked_range, 
                           fill=comp)) + 
                       geom_boxplot(position=position_dodge(0.9), fatten=3)

b.SM1 +  theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill="white"),
                  panel.border = element_rect(colour="white", fill=NA, size=1)) + 
                  theme(axis.line = element_line(colour="black"))+ 
                  xlab("\\nStimulus presentation") + 
                  ylab("Area scent marked (cm^2)\\n")+ 
                  labs(fill = "Nest site availability\\n")+
                  scale_fill_manual(name = "Nest site availability", 
                        labels = c("Single nest site", "Multiple nest sites"),
                        values = c("#E6DBD1", "#878785")) + 
                  scale_x_discrete(labels = c("Control",
                  "Related urine", "Unrelated urine" )) + 
                 theme(axis.text.x = element_text(size = 12),
                        axis.text.y = element_text(size = 12),  
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_text(size = 14)) + 
                  theme(legend.title=element_text(size=14), 
                        legend.text=element_text(size=12)) +
                  theme(legend.key=element_blank())


############
## MODEL FOR TABLE 2B "AREA SCENT MARKED"

m1.SUM.a<-lmer(log(sum_area_marked_range+0.01) ~ stimulus3 * re * comp + 
               + (1|ID) + (1|Trial), data=ccm.c.ur.re)

# OBTAIN THE ESTIMATES
summary(m1.SUM.a)

# REMOVE THREE-WAY INTERACTION
m2.SUM.a<-update(m1.SUM.a,~. -stimulus3:re:comp)

# OBTAIN ESTIMATES
summary(m2.SUM.a)

# REMOVE TWO-WAY INTERACTION RELATEDNESS OF COMPETITORS * NEST SITE 
# AVAILABILITY
m3.SUM.a<-update(m2.SUM.a,~. -re:comp)

# OBTAIN ESTIMATES
summary(m3.SUM.a)

# REMOVE TWO-WAY INTERACTION STIMULUS * RELATEDNESS OF COMPETITORS
m4.SUM.a<-update(m3.SUM.a,~. -stimulus3:re)

# OBTAIN ESTIMATES
summary(m4.SUM.a)

# OBTAIN P-VALUES FOR THE TWO-WAY INTERACTION STIMULUS * NEST SITE 
# AVAILABILITY
drop1(m4.SUM.a, test="F")

### GET THE P-VALUES FOR ALL OTHER FACTOR WITH "mixed()"
m.SUM.a<-mixed(log(sum_area_marked_range+0.01) ~ stimulus3 + re + comp + 
               comp:stimulus3 + (1|ID) + (1|Trial), data=ccm.c.ur.re)

# DISPLAY THE P-VALUES
anova(m.SUM.a)

# CHECK THE MODEL
par(mfrow=c(2,2))
plot(resid(m1.SUM.a)~fitted(m1.SUM.a))
abline(h=0, lty=2)
hist(resid(m1.SUM.a))
qqnorm(resid(m1.SUM.a))
qqline(resid(m1.SUM.a))
resid.m1.SUM.a<-resid(m1.SUM.a)
shapiro.test(resid.m1.SUM.a)
lillie.test(resid.m1.SUM.a) 

# OBTAIN MARGINAL MEANS FOR TWO-WAY INTERACTION 
# STIMULUS * NEST SITE AVAILABILITY PRESENTED IN FIGURE 3 (NEED TO BACK
# TRANSFORM THE VALUES WITH exp())
lsmeans(m4.SUM.a, pairwise~stimulus3:comp, adjust="tukey")

############
## PERFORM ORTHOGONAL CONTRASTS TO INVESTIGATE THE SIGNIFICANT TWO-WAY 
## INTERACTION STIMULUS * NEST SITE AVAILABILITY

m.SUM.a.n<-lmer(log(sum_area_marked_range+0.01) ~ 
               stimulus_comp + (1|ID) + (1|Trial), data=ccm.new)

# OBTAIN ESTIMATES AND P-VALUES PRESENTED IN TABLE S7
summary(m.SUM.a.n)


###############################################################################
##                                                                           ##
##                       MODEL AND FIGURE FOR TABLE 3A                       ## 
##                        "TOTAL COMMUNAL LITTER SIZE"                       ##
##                                                                           ##
###############################################################################


############
## CODE FOR FIGURE 4 "TOTAL LITTER SIZES OF COMMUNAL NESTS"

# REMOVE NAs FROM THE DATASET FOR THE PLOT
maf.d1.c_1 <- maf.d1.c[!is.na(maf.d1.c$nr.comm.pups),]

# THE BOXPLOT
TL.b<-ggplot(maf.d1.c_1, aes(x=comp, y=nr.comm.pups, fill=rel)) + 
                   geom_boxplot(position=position_dodge(0.9), fatten=3) 

TL.b + theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white"),
             panel.border = element_rect(colour="white", fill=NA, size=1)) + 
             theme(axis.line = element_line(colour="black"))+ 
             xlab("\\nNest site availability") + 
             ylab("Total litter sizes of communal nests\\n")+ 
             labs(fill = "Relatedness of\\ncompetitors")+
             scale_fill_discrete(name ="Relatedness of\\ncompetitors", 
             labels = c("Related", "Unrelated")) + 
             scale_x_discrete(labels = c("Single nest site",
             "Multiple nest sites")) + 
             coord_cartesian(ylim=c(8,16)) +
             theme(axis.text.x = element_text(size = 12),
                   axis.text.y = element_text(size = 12),  
                   axis.title.x = element_text(size = 14),
                   axis.title.y = element_text(size = 14)) + 
             theme(legend.title=element_text(size=14), 
                   legend.text=element_text(size=12)) +
              theme(legend.key=element_blank())

############
## MODEL FOR TABLE 3A "TOTAL COMMUNAL LITTER SIZE"

NWOCN.1<-lm(log(nr.comm.pups) ~ rel + comp, data=maf.d1.c)

# OBTAIN ESTIMATES
summary(NWOCN.1)

# OBTAIN P-VALUES
drop1(NWOCN.1, test="F")

# CHECK THE MODEL
par(mfrow=c(2,2))
plot(resid(NWOCN.1)~fitted(NWOCN.1))
abline(h=0, lty=2)
hist(resid(NWOCN.1))
qqnorm(resid(NWOCN.1))
qqline(resid(NWOCN.1))
resid.NWOCN.1<-resid(NWOCN.1)
shapiro.test(resid.NWOCN.1)
lillie.test(resid.NWOCN.1)


###############################################################################
##                                                                           ##
##                            MODEL FOR TABLE 3B                             ## 
##                            "REPRODUCTIVE SKEW"                            ##
##                                                                           ##
###############################################################################


############
## MODEL FOR TABLE 3B "REPRODUCTIVE SKEW"

RS.CN1<-lm(log(diff.pups+1) ~ comp + rel + nr.comm.pups, 
             data=maf.d1.c)

# OBTAIN ESTIMATES
summary(RS.CN1)

# OBTAIN P-VALUES
drop1(RS.CN1, test="F")

# CHECK THE MODEL
par(mfrow=c(2,2))
plot(resid(RS.CN1)~fitted(RS.CN1))
abline(h=0, lty=2)
hist(resid(RS.CN1))
qqnorm(resid(RS.CN1))
qqline(resid(RS.CN1))
resid.RS.CN1<-resid(RS.CN1)
shapiro.test(resid.RS.CN1)
lillie.test(resid.RS.CN1)


###############################################################################
##                                                                           ##
##                       MODEL AND FIGURE FOR TABLE 4                        ## 
##                             "WEANING MASSES"                              ##
##                                                                           ##
###############################################################################


############
## CODE FOR FIGURE 2B "AVERAGE WEANING MASSES"

w.b<-ggplot(wwp.ALL.w, aes(x=Competition, y=weight, fill=Relatedness)) + 
                    geom_boxplot(position=position_dodge(0.9), fatten=3)

w.b + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill="white"),
              panel.border = element_rect(colour="white", fill=NA, size=1)) + 
              theme(axis.line = element_line(colour="black"))+ 
              xlab("\\nNest site availability") + 
              ylab("Average weaning mass (in g)\\n")+ 
              labs(fill = "Relatedness of\\ncompetitors")+
              scale_fill_discrete(name ="Relatedness of\\ncompetitors", 
              labels = c("Related", "Unrelated")) + 
              scale_x_discrete(labels = c("Single nest site",
              "Multiple nest sites")) + 
               theme(axis.text.x = element_text(size = 12),
                     axis.text.y = element_text(size = 12),  
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14)) + 
              theme(legend.title=element_text(size=14), 
                    legend.text=element_text(size=12))+
              theme(legend.key=element_blank())


############
## MODEL FOR TABLE 4 "WEANING MASSES"

wwp.m1.c<-lmer(log(weight) ~ Competition + Relatedness + Age_start + Sex +
             + total 
             + (1|sister_ID) + (1|block), data=wwp.ALL.w)

# THE FIT OF THIS MODEL IS SINGULAR BECAUSE THE RANDOM FACTOR BLOCK IDENTITY
# ACCOUNTS FOR 0 VARIANCE. WE STILL THINK IT IS IMPORTANT TO INCLUDE THIS 
# RANDOM EFFECT BECAUSE IT STILL ACCOUNTS FOR THE MINOR DIFFERENCES BETWEEN 
# THE BLOCKS

# OBTAIN THE ESTIMATES
summary(wwp.m1.c)

# OBTAIN THE P-VALUES
drop1(wwp.m1.c, test="F")

# REMOVED FACTORS
# + dam_weight_w

# CHECK THE MODEL
par(mfrow=c(2,2))
plot(resid(wwp.m1.c)~fitted(wwp.m1.c))
abline(h=0, lty=2)
hist(resid(wwp.m1.c))
qqnorm(resid(wwp.m1.c))
qqline(resid(wwp.m1.c))
resid.wwp.m1.c<-resid(wwp.m1.c)
shapiro.test(resid.wwp.m1.c)
lillie.test(resid.wwp.m1.c)


###############################################################################
##                                                                           ##
##                            MODEL FOR TABLE S4                             ## 
##               "WHETHER A FEMALE WEANED A LITTER (YES,NO)"                 ##
##                                                                           ##
###############################################################################


############
## MODEL FOR TABLE S4 "WHETHER A FEMALE WEANED A LITTER (YES,NO)"

LTB.m1<-glmer(pregnant ~ rel + comp + (1|Nest) + (1|block), 
              family=binomial, data=maf.d1)

# THE FIT OF THIS MODEL IS SINGULAR BECAUSE THE RANDOM FACTOR BLOCK IDENTITY
# ACCOUNTS FOR 0 VARIANCE. WE STILL THINK IT IS IMPORTANT TO INCLUDE THIS 
# RANDOM EFFECT BECAUSE IT STILL ACCOUNTS FOR THE MINOR DIFFERENCES BETWEEN 
# THE BLOCKS

# OBTAIN ESTIMATES
summary(LTB.m1)

# OBTAIN P-VALUES
drop1(LTB.m1, test="Chisq")


###############################################################################
##                                                                           ##
##                            MODEL FOR TABLE S5                             ## 
##              "WHETHER A COMMUNAL NEST WAS FORMED (YES,NO)"                ##
##                                                                           ##
###############################################################################


############
## MODEL FOR TABLE S5 "WHETHER A COMMUNAL NEST WAS FORMED (YES,NO)"

LTFC.m1<-glmer(comm.nest2 ~ rel + comp + (1|block), 
               family=binomial, data=maf.d1)

# OBTAIN THE ESTIMATES
summary(LTFC.m1)

#OBTAIN THE P-VALUES
drop1(LTFC.m1, test="Chisq")


###############################################################################
##                                                                           ##
##                            MODEL FOR TABLE S6A                            ## 
##                "ACTIVITY OF FEMALES DURING THE DARK PHASE"                ##
##                                                                           ##
###############################################################################


######
### MODEL FOR TABLE S6A "ACTIVITY OF FEMALES DURING THE DARK PHASE"

Ac.m1 <- lmer(log(activityD) ~ rel + comp + log(time.in.nest.D) + (1|Trial), 
              data=maf.d1)

# OBTAIN THE ESTIMATES
summary(Ac.m1)

# OBTAIN THE P-VALUES
drop1(Ac.m1, test="F")

# CHECK THE MODEL
par(mfrow=c(2,2))
plot(resid(Ac.m1)~fitted(Ac.m1))
abline(h=0, lty=2)
hist(resid(Ac.m1))
qqnorm(resid(Ac.m1))
qqline(resid(Ac.m1))
resid.Ac.m1<-resid(Ac.m1)
shapiro.test(resid.Ac.m1)
lillie.test(resid.Ac.m1)


###############################################################################
##                                                                           ##
##                            MODEL FOR TABLE S6B                            ## 
##                          "NUMBER OF SCENT MARKS"                          ##
##                                                                           ##
###############################################################################

######
### MODEL FOR TABLE S6B "NUMBER OF SCENT MARKS"

m1.SUM.n<-lmer(log(sum_number_scent_marks_range+1) ~ stimulus3 * re * comp + 
               + (1|ID) + (1|Trial), data=ccm.c.ur.re)

# OBTAIN ESTIMATES
summary(m1.SUM.n)

# REMOVE THREE-WAY INTERACTION STUMULUS * RELATEDNESS OF COMPETITORS * 
# NEST SITE AVAILABILITY
m2.SUM.n<-update(m1.SUM.n,~. -stimulus3:re:comp)

# OBTAIN ESTIMATES
summary(m2.SUM.n)

# REMOVE TWO-WAY INTERACTION RELATEDNESS OF COMPETITORS * NEST SITE 
# AVAILABILITY
m3.SUM.n<-update(m2.SUM.n,~. -re:comp)

# OBTAIN ESTIMATES
summary(m3.SUM.n)

# REMOVE TWO-WAY INTERACTION STUMULUS * RELATEDNESS OF COMPETITORS
m4.SUM.n<-update(m3.SUM.n,~. -stimulus3:re)

# OBTAIN ESTIMATES
summary(m4.SUM.n)

# CHECK WHETHER THE TWO-WAY INTERACTION STIMULUS * NEST SITE AVAILABILITY
# IS SIGNIFICANT
drop1(m4.SUM.n, test="F")

# REMOVE THE TWO-WAY INTERACTION STUMULUS * NEST SITE AVAILABILITY
m5.SUM.n<-update(m4.SUM.n,~. -stimulus3:comp)

# OBTAIN THE ESTIMATES
summary(m5.SUM.n)

# OBTAIN P-VALUES
drop1(m5.SUM.n, test="F")

# CHECK THE MODEL
par(mfrow=c(2,2))
plot(resid(m1.SUM.n)~fitted(m1.SUM.n))
abline(h=0, lty=2)
hist(resid(m1.SUM.n))
qqnorm(resid(m1.SUM.n))
qqline(resid(m1.SUM.n))
resid.m1.SUM.n<-resid(m1.SUM.n)
shapiro.test(resid.m1.SUM.n)
lillie.test(resid.m1.SUM.n) 


###############################################################################
##                                                                           ##
##                      MODEL AND FIGURE FOR TABLE S8                        ## 
## "SLEEPING NEST SITE ASSOCIATION BETWEEN OLDER AND YOUNGER SISTER PAIRS"   ##
##                                                                           ##
###############################################################################


############
## CODE FOR FIGURE S2 "ASSOCIATION INDEX BETWEEN OLDER AND YOUNGER SISTER 
## PAIR"

## CREATE THE RATIO OF HOW OFTEN OLDER AND YOUNGER SISTERS PAIRS WERE FOUND IN
## THE SAME NEST BOX DIVIDED BY THE NUMBER OF OBSERVATIONS
SN$ratio <- SN$found / (SN$found + SN$not.found)
 

# THE BOXPLOT
AI.SN.box<-ggplot(SN, aes(x=treat2, y=ratio, fill=treat2)) + 
                    geom_boxplot(position=position_dodge(0.9), fatten=3)

AI.SN.box + theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill="white"),
              panel.border = element_rect(colour="white", fill=NA, size=1)) + 
              theme(axis.line = element_line(colour="black"))+ 
              xlab("\\nRelatedness of competitors") + 
              ylab("Association index between older and younger sister pairs\\n")+ 
              scale_fill_manual(values = c("#F8766D", "#00BFC4")) + 
              scale_x_discrete(labels = c("Related", "Unrelated")) +
              theme(axis.text.x = element_text(size = 12),
                     axis.text.y = element_text(size = 12),  
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14)) + 
              theme(legend.title=element_text(size=14), 
                    legend.text=element_text(size=12)) +
              theme(legend.position="none")
             
######
### MODEL FOR TABLE S8 "SLEEPING NEST SITE ASSOCIATION BETWEEN OLDER AND 
### YOUNGER SISTER PAIRS"

# CREATING THE DEPENDENT VARIABLE WITH SUCCESSESS (NUMBER OF TIMES OLDER AND 
# YOUNGER SISTER PAIRS WERE FOUND TOGETHER AND FAILURES (NUMBER OF TIME OLDER
# AND YOUNGER SISTER PAIRS WERE NOT FOUND TOGETHER
y1 <- cbind(SN$found, SN$not.found)

# THE BINOMIAL MODEL TO ANALYSE PROPORTIONS
m4 <- glm(y1 ~ treat2, family=binomial, data=SN)

# OBTAIN THE ESTIMATES
summary(m4)

# OBTAIN THE P-VALUES
drop1(m4, test="Chisq")


###############################################################################
##                                                                           ##
##                      MODEL AND FIGURE FOR TABLE S9                        ## 
##   "TOTAL TIME SPENT IN NEST OF COMMUNAL BREEDING FEMALES FROM PND 0-14"   ##
##                                                                           ##
###############################################################################


############
## CODE FOR FIGURE S3 "TOTAL TIME SPENT IN NEST OF COMMUNAL BREEDING FEMALES"  

## CALCULATE THE TOTAL TIME SPENT IN NEST OF COOMUNAL BREEDING FEMALES ONLY
maf.d1.c$total.time.s <- maf.d1.c$time.in.nest.D + maf.d1.c$time.in.nest.L

# CONVERSION INTO HOURS
maf.d1.c$total.time.h <- maf.d1.c$total.time.s / 3600

# REMOVE NAs FROM THE DATASET FOR THE PLOT
maf.d1.c_2 <- maf.d1.c[!is.na(maf.d1.c$total.time.h),]

# THE BOXPLOT
TT.c<-ggplot(maf.d1.c_2, aes(x=rel, y=total.time.h, fill=rel)) + 
             geom_boxplot(position=position_dodge(0.9), fatten=3) 
TT.c + theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_rect(fill="white"),
             panel.border = element_rect(colour="white", fill=NA, size=1)) + 
             theme(axis.line = element_line(colour="black"))+ 
             theme(legend.position="none") +            
             xlab("\\nRelatedness of competitors") + 
             scale_y_continuous("Total time spent in nest box from PND 0-14\\nof communal breeding females (in hours)\\n", labels=scales::comma) +             
             scale_x_discrete(labels = c("Related","Unrelated")) + 
             theme(axis.text.x = element_text(size = 12),
                   axis.text.y = element_text(size = 12),  
                   axis.title.x = element_text(size = 14),
                   axis.title.y = element_text(size = 14)) + 
             theme(legend.title=element_text(size=14), 
                   legend.text=element_text(size=12)) +
             theme(legend.key=element_blank())

######
### MODEL FOR TABLE S9 "TOTAL TIME SPENT IN NEST OF COMMUNAL BREEDING FEMALES" 

time.c<-lmer(log(total.time.h) ~ rel  +(1|Trial), data=maf.d1.c)

# OBTAIN THE ESTIMATES
summary(time.c)

# OBTAIN THE P-VALUES
drop1(time.c,test="F")

# CHECK THE MODEL
par(mfrow=c(2,2))
plot(resid(time.c)~fitted(time.c))
abline(h=0, lty=2)
hist(resid(time.c))
qqnorm(resid(time.c))
qqline(resid(time.c))
resid.time.c<-resid(time.c)
shapiro.test(resid.time.c)
lillie.test(resid.time.c)


###############################################################################
##                                                                           ##
##                         CODE FOR FIGURE S1                                ## 
##          "ACTIVITY OF FEMALES DURING PEAK MATERNAL INVESTMENT"            ##
##                                                                           ##
###############################################################################

############
## CODE FOR FIGURE S1 

# CHANGE PND TO A FACTOR
lall$PND <- as.factor(lall$PND)

# THE BOXPLOT
AC <- ggplot(lall, aes(x=PND, y=activity, 
            fill=phase)) + 
                  geom_boxplot()
# NAME THE LEGEND
legend_title <- "Time of day"

AC + theme_classic()+
              xlab("\\nPost natal day") + 
              ylab("Activity of females\\n") + 
              coord_cartesian(ylim=c(0,2500)) +
              scale_fill_manual(legend_title, values=c("#808080", "seashell1"))+ 
              theme(axis.text.x = element_text(size = 12),
                     axis.text.y = element_text(size = 12),  
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14)) + 
              theme(legend.title=element_text(size=14), 
                    legend.text=element_text(size=12)) + 
               geom_vline(xintercept = 1.5, linetype=2) +
               geom_vline(xintercept = 2.5, linetype=2) +
               geom_vline(xintercept = 3.5, linetype=2) +
               geom_vline(xintercept = 4.5, linetype=2) +
               geom_vline(xintercept = 5.5, linetype=2) +
               geom_vline(xintercept = 6.5, linetype=2) +
               geom_vline(xintercept = 7.5, linetype=2) +
               geom_vline(xintercept = 8.5, linetype=2) +
               geom_vline(xintercept = 9.5, linetype=2) +
               geom_vline(xintercept = 10.5, linetype=2) +
               geom_vline(xintercept = 11.5, linetype=2) +
               geom_vline(xintercept = 12.5, linetype=2) +
               geom_vline(xintercept = 13.5, linetype=2) +
               geom_vline(xintercept = 14.5, linetype=2)


###############################################################################


