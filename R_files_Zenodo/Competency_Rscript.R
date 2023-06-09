library(ggplot2)
library(gridExtra)
library(plyr)
library(multcomp)
library(MASS)
library(lme4)
library(nlme)
library(phyloseq)
library(ellipse)
library(vegan)
#SummarySE and lm_eqn functions at end of script

#########################################################
##Dry weights larvae
#########################################################

larvalweight_Competency<-read.table("larvalweight_Competency.csv", header=T, sep=",", strip.white=T)
larvalweight_Competency$cross<-as.factor(larvalweight_Competency$cross)

#Reproductive crosses
(summary_weights<-summarySE(larvalweight_Competency, measurevar="weight", groupvars=c("cross"), na.rm=TRUE))

fig1_a_mean<-ggplot()+ geom_jitter(aes(cross, weight), data = larvalweight_Competency, colour = I("grey"),position = position_jitter(width = 0.05))+
  geom_point(data=summary_weights,aes(x=cross,ymin=weight, ymax=weight,y=weight,group=cross))+
  labs(x="", y="larval dry weight (micrograms)")+theme_classic()+theme(axis.text.x = element_blank())+theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+
  geom_pointrange(aes(x=summary_weights$cross, y=summary_weights$weight,ymin=summary_weights$weight-summary_weights$se, ymax=summary_weights$weight+summary_weights$se))


#Parental identity
(summary_weights_dam<-summarySE(larvalweight_Competency, measurevar="weight", groupvars=c("dam"), na.rm=TRUE))

fig2_a_mean<-ggplot()+ geom_jitter(aes(dam, weight), data = larvalweight_Competency, colour = I("grey"),position = position_jitter(width = 0.05))+
  geom_point(data=summary_weights_dam,aes(x=dam,ymin=weight, ymax=weight,y=weight,group=dam))+
  labs(x="", y="larval dry weight (micrograms)")+theme_classic()+theme(axis.text.x = element_blank())+theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+
  geom_pointrange(aes(x=summary_weights_dam$dam, y=summary_weights_dam$weight,ymin=summary_weights_dam$weight-summary_weights_dam$se, ymax=summary_weights_dam$weight+summary_weights_dam$se))

#lm's
#Pure vs hybrids
(summary_weights_purity<-summarySE(larvalweight_Competency, measurevar="weight", groupvars=c("purity"), na.rm=TRUE))
weights.pure.lm<-lm(weight~purity, data=larvalweight_Competency)
par(mfrow=c(2,3))
plot(weights.pure.lm, ask=F, which=1:6)
summary(weights.pure.lm) #p = 0.774
anova(weights.pure.lm) #provides Df and F value

#Population crosses
weights.cross.lm<-lm(weight~cross, data=larvalweight_Competency)
par(mfrow=c(2,3))
plot(weights.cross.lm, ask=F, which=1:6)
summary(weights.cross.lm) #base level is OxO 
summary(glht(weights.cross.lm, linfct=mcp(cross="Tukey")))

#Dams
weights.dam.lm<-lm(weight~dam, data=larvalweight_Competency)
plot(weights.dam.lm, ask=F, which=1:6)
summary(weights.dam.lm)
summary(glht(weights.dam.lm, linfct=mcp(dam="Tukey")))
 
#########################################################
##Larval culture survival
#########################################################

#3 culture reps per family. If all 3 culture reps alive by start of settlement = alive. Binary: all cultures alive or dead.
#All cultures for Family 25 (W7 x O5) and Family 29 (W7 x W11) died. So 23/25 families survived.

#Pure vs. hybrid mixtures. 
Purity_larvsurv <-
  matrix(c(1, 1, 12, 11),
         nrow = 2,
         dimnames =
           list(c("Pure", "Mix"),
                c("Dead", "Alive")))
fisher.test(Purity_larvsurv, alternative="two.sided", conf.int=TRUE, simulate.p.value=TRUE, B=2000)

#Population crosses
Cross_larvsurv <-
  matrix(c(0, 0, 1, 1, 6, 6, 5, 6),
         nrow = 4,
         dimnames =
           list(c("OxO", "OxW", "WxO", "WxW"),
                c("Dead", "Alive")))
fisher.test(Cross_larvsurv, hybrid=TRUE, alternative="two.sided", conf.int=TRUE, simulate.p.value=TRUE, B=2000)  

#Dams
dam_larvsurv <-
  matrix(c(3, 4, 3, 3, 2, 5, 3, 0,0,0,0,0,0,0,0,2),
         nrow = 8,
         dimnames =
           list(c("O3", "O5", "W11", "W5", "O6", "W10", "O4", "W7"),
                c("Alive", "Dead")))
fisher.test(dam_larvsurv, alternative="two.sided", conf.int=TRUE, simulate.p.value=TRUE, B=2000)

#Sires
sire_larvsurv <-
  matrix(c(3, 3, 4, 3, 4, 3, 1, 2,0,0,0,0,1,0,1,0),
         nrow = 8,
         dimnames =
           list(c("O4", "O3", "O6", "W10", "W11", "W5", "O5", "W7"),
                c("Alive", "Dead")))
fisher.test(sire_larvsurv, alternative="two.sided", conf.int=TRUE, simulate.p.value=TRUE, B=2000)


#########################################################
#Settlement 
#########################################################
#Removed family 25 and 29 becuase they didn't survive larval stage.NAs
#23 total families at this point. Three families didn't have any settlers by juvenile outplant stage. They kept swimming around. F21 (W10 x O5), F26 (W10 x O4), F30 (W5 x W7)
#20 families had settlers by 19 d.p.f
Settlment_Competency<-read.table("Settlement_Competency.csv", header=T, sep=",", strip.white=T)

##Pure vs. hybrid mixtures
#Aggregate settlement reps because reps are random effect in glmm
aggdata_fam<-ddply(Settlment_Competency, .(family), summarize, sum_sett=sum(settlementAbund))
aggdata_fam$purity<-c("pure","mix","pure","pure","pure","mix","mix","pure","pure","mix","mix","pure","pure","pure","mix","mix","mix","pure","mix","pure","mix","pure","mix")
(aggdata_fam)
(summarySE_aggdata_fam_purity<-summarySE(aggdata_fam, measurevar="sum_sett", groupvars=("purity")))
aggdata_fam$purity<-as.factor(aggdata_fam$purity)
kruskal.test(sum_sett~purity, data=aggdata_fam)

##Population crosses
aggdata_fam$cross<-c("pureo","mixo","pureo","pureo","pureo","mixo","mixw","purew","purew","mixo","mixw","purew","purew","purew","mixw","mixw","mixw","purew","mixo","pureo","mixo","pureo","mixo")
aggdata_fam$cross<-as.factor(aggdata_fam$cross)
aggdata_fam$cross <- factor(aggdata_fam$cross, levels = c("pureo", "mixo", "mixw", "purew")) #make sure to have the same order of factors on X axis for Figures 1 and 2
levels(aggdata_fam$cross)
(summarySE_aggdata_fam_cross<-summarySE(aggdata_fam, measurevar="sum_sett", groupvars=("cross")))
kruskal.test(sum_sett~cross, data=aggdata_fam) #Signif. at 0.059. check post-hoc
settlment.cross.glm_nb_aggdata_fam<-glm.nb(sum_sett~cross, data=aggdata_fam)
summary(settlment.cross.glm_nb_aggdata_fam)
summary(glht(settlment.cross.glm_nb_aggdata_fam, linfct=mcp(cross="Tukey")))

fig1_b_mean<-ggplot()+ geom_jitter(aes(cross, sum_sett), data = aggdata_fam, colour = I("grey"),position = position_jitter(width = 0.05))+
  geom_point(data=summarySE_aggdata_fam_cross,aes(x=cross,ymin=sum_sett, ymax=sum_sett,y=sum_sett,group=cross))+labs(x="", y="Number of settled juveniles")+
  theme_classic()+theme(axis.text.x = element_blank())+ylim(0,300)+theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+
  geom_pointrange(aes(x=summarySE_aggdata_fam_cross$cross, y=summarySE_aggdata_fam_cross$sum_sett,ymin=summarySE_aggdata_fam_cross$sum_sett-summarySE_aggdata_fam_cross$se, ymax=summarySE_aggdata_fam_cross$sum_sett+summarySE_aggdata_fam_cross$se))


##Dam, sire, and dam x sire effects
#Poisson, observation random effect to deal with overdispersion, dam, sire, dam x sire interaction as random effects
set.pois.6_M<-glmer(settlementAbund ~ (1|dam) + (1|sire) + (1|dam:sire) + (1|obs), family=poisson, data=Settlment_Competency)
X.var<- var(as.vector(lme4::fixef(set.pois.6_M) %*% t(set.pois.6_M@pp$X)))
## Extract the variance components for the random effects (not including the residuals)
Z.var <- sum(
  sapply(
    VarCorr(set.pois.6_M)[!sapply(unique(unlist(strsplit(names(ranef(set.pois.6_M)),":|/"))), function(l)
      length(unique(set.pois.6_M@frame[,l])) == nrow(set.pois.6_M@frame))],
    function(Sigma) {
      X <- model.matrix(set.pois.6_M)
      Z <- X[,rownames(Sigma)]
      sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X) } ) )
## Extract the variance componts for the residuals
R.var <- attr(lme4::VarCorr(set.pois.6_M), "sc")^2
## The marginal R2 (proportion of variance due to fixed effects)
(R2.marginal <- X.var/(X.var+Z.var+R.var))
# 0 % of variance due to fixed effects: no fixed effects
## The proportion of variance due to random effects
(R2.random <- Z.var/(X.var+Z.var+R.var))
# 61.9% % (mom and dad and interaciton and obs)
## The proportion of variance due to residuals
VarCorr(set.pois.6_M) 
#obs=0.91586= 28% = 0.91586/(0.91586+ 0.32216+1.23403 +0.77704)*100
#dam:sire=0.32216=9.91539%
#sire= 1.23403 = 37.9807 %
#dam=0.77704= 23.9156 %
#total=3.24909
#THIS ALL MAKES UP 61% of the total variation. So dam is 23.9% of 61% = 14.579
(R2.resid <- R.var/(X.var+Z.var+R.var))
# 0.38% just due to residuals
## The conditional R2 (proportion of variance due to fixed and random effects)
(R2.conditional <- (X.var+Z.var)/(X.var+Z.var+R.var)) #.61%

#with Dam O6 as reference level
Settlment_Competency$dam <- relevel(Settlment_Competency$dam, ref="O6")
set.pois.6<-glmer(settlementAbund ~ dam + (1|sire) + (1|obs), family=poisson, data=Settlment_Competency)
summary(set.pois.6)
summary(glht(set.pois.6), linfct=mcp(dam="Tukey"))

#with sire W10 as reference
Settlment_Competency$sire <- relevel(Settlment_Competency$sire, ref="W10")
set.pois.6.sire<-glmer(settlementAbund ~ sire + (1|dam) + (1|obs), family=poisson, data=Settlment_Competency)
summary(set.pois.6.sire) #some convergence warnings due to low sample sizes in some comparisons (only one family with sire O5 left at this stage)
summary(glht(set.pois.6.sire), linfct=mcp(sire="Tukey"))

#Dam figure averages
aggdata_fam$dam<-c("O4","O3","O5","O5","O5","O5","W11","W11","W11","O4","W10","W10","W10","W10","W10","W5","W5","W5","O4","O6","O6","O3","O3")
aggdata_fam$dam <- factor(aggdata_fam$dam, levels = c("O3","O4","O5","O6","W10","W11","W5")) #make sure to have the same order of factors on X axis for Figure 2
summarySE_settle_dam<-summarySE(aggdata_fam, measurevar="sum_sett", groupvars=("dam"))
fig2_b_mean<-ggplot()+ geom_jitter(aes(dam, sum_sett), data = aggdata_fam, colour = I("grey"),position = position_jitter(width = 0.05))+
  geom_point(data=summarySE_settle_dam,aes(x=dam,ymin=sum_sett, ymax=sum_sett,y=sum_sett,group=dam))+
  labs(x="", y="Number of settled recruits")+theme_classic()+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+
  geom_pointrange(aes(x=summarySE_settle_dam$dam, y=summarySE_settle_dam$sum_sett,ymin=summarySE_settle_dam$sum_sett-summarySE_settle_dam$se, ymax=summarySE_settle_dam$sum_sett+summarySE_settle_dam$se))

#sire figure averages
aggdata_fam$sire<-c("O6","W10","O4","O6","O3","W11","O6","W10","W5","W11","O5","W11","W7","W5","O4","O6","O3","W7","W5","O3","W10","O4","W11")
aggdata_fam$sire <- factor(aggdata_fam$sire, levels = c("O3","O4","O5","O6","W10","W11","W5", "W7")) #make sure to have the same order of factors on X axis for Figure 2
summarySE_settle_sire<-summarySE(aggdata_fam, measurevar="sum_sett", groupvars=("sire"))
fig2_c_mean<-ggplot()+ geom_jitter(aes(sire, sum_sett), data = aggdata_fam, colour = I("grey"),position = position_jitter(width = 0.05)) +
  geom_point(data=summarySE_settle_sire,aes(x=sire,ymin=sum_sett, ymax=sum_sett,y=sum_sett,group=sire))+
  labs(x="", y="Number of settled recruits")+theme_classic()+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+
  geom_pointrange(aes(x=summarySE_settle_sire$sire, y=summarySE_settle_sire$sum_sett,ymin=summarySE_settle_sire$sum_sett-summarySE_settle_sire$se, ymax=summarySE_settle_sire$sum_sett+summarySE_settle_sire$se))
fig2_c_mean

#########################################################
#Survival/mortality in field
#########################################################
#Removed Families 26 (W10 x O4), F30 (W5 x W7) and F21 (W10 x O5) because they didn't have any settlers. 
Mortality_Competency<-read.table("Mortality_Competency.csv", header=T, sep=",", strip.white=T)

#Purity vs. hybrid effects
kruskal.test(Percentmortality~purity, data=Mortality_Competency) #p =0.73
(summarySE_mort_purity<-summarySE(Mortality_Competency, measurevar="Percentmortality", groupvars=("purity")))

#Cross
kruskal.test(Percentmortality~cross, data=Mortality_Competency) #p =0.47
(summarySE_mort_cross<-summarySE(Mortality_Competency, measurevar="Percentmortality", groupvars=("cross")))

fig1_c_mean<-ggplot()+ geom_jitter(aes(cross, Percentmortality), data = Mortality_Competency, colour = I("grey"),position = position_jitter(width = 0.05)) +
  geom_point(data=summarySE_mort_cross,aes(x=cross,ymin=Percentmortality, ymax=Percentmortality,y=Percentmortality,group=cross))+labs(x="Population cross", y="Percent mortality of juveniles")+theme_classic()+
  geom_pointrange(aes(x=summarySE_mort_cross$cross, y=summarySE_mort_cross$Percentmortality,ymin=summarySE_mort_cross$Percentmortality-summarySE_mort_cross$se, ymax=summarySE_mort_cross$Percentmortality+summarySE_mort_cross$se))+
  scale_x_discrete(breaks=c("OO","OW","WO","WW"), labels=c("OO", "OW", "WO", "WW"))+theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
fig1_c_mean

Figure1.composite<-grid.arrange(fig1_a_mean,fig1_b_mean,fig1_c_mean, ncol=1)

#Dams
(summarySE_mort_dams<-summarySE(Mortality_Competency, measurevar="Percentmortality", groupvars=("dam")))
fig2_d_mean<-ggplot()+ geom_jitter(aes(dam, Percentmortality), data = Mortality_Competency, colour = I("grey"),position = position_jitter(width = 0.05)) +
  geom_point(data=summarySE_mort_dams,aes(x=dam,ymin=Percentmortality, ymax=Percentmortality,y=Percentmortality,group=dam))+labs(x="Parental Identity", y="Percent mortality of juveniles")+theme_classic()+
  geom_pointrange(aes(x=summarySE_mort_dams$dam, y=summarySE_mort_dams$Percentmortality,ymin=summarySE_mort_dams$Percentmortality-summarySE_mort_dams$se, ymax=summarySE_mort_dams$Percentmortality+summarySE_mort_dams$se))+
  scale_x_discrete(breaks=c("O3","O4","O5","O6", "W10", "W11", "W5", "W7"), labels=c("O3","O4","O5","O6", "W10", "W11", "W5", "W7"))+theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
fig2_d_mean

Figure2.composite<-grid.arrange(fig2_a_mean,fig2_b_mean, fig2_c_mean, fig2_d_mean, ncol=1) #fig2_c_mean has extra parent at end of x-axis (W7)
#best to inport fig2a, b and d into illustrator together than seperatly add in fig2_c_mean because it has 8 columns intead of 7

#Familial and parental effects on juvenile survivorhsip
a<-glm.nb(abundance~settlement, data=Mortality_Competency) # estimate with only fixed effects
theta <- a$theta
nb8.glmm.gHQ.M<- glmer(abundance ~ settlement +(1|dam) +(1|sire)+(1|dam:sire), data=Mortality_Competency, family=negative.binomial(theta))
summary(nb8.glmm.gHQ.M)
#negative binomial for overdispersion.

X.var<- var(as.vector(lme4::fixef(nb8.glmm.gHQ.M) %*% t(nb8.glmm.gHQ.M@pp$X)))
## Extract the variance components for the random effects (not including the residuals)
Z.var <- sum(
  sapply(
    VarCorr(nb8.glmm.gHQ.M)[!sapply(unique(unlist(strsplit(names(ranef(nb8.glmm.gHQ.M)),":|/"))), function(l)
      length(unique(nb8.glmm.gHQ.M@frame[,l])) == nrow(nb8.glmm.gHQ.M@frame))],
    function(Sigma) {
      X <- model.matrix(nb8.glmm.gHQ.M)
      Z <- X[,rownames(Sigma)]
      sum(diag(Z %*% Sigma %*% t(Z)))/nrow(X) } ) )
## Extract the variance componts for the residuals
R.var <- attr(lme4::VarCorr(nb8.glmm.gHQ.M), "sc")^2
## The marginal R2 (proportion of variance due to fixed effects)
(R2.marginal <- X.var/(X.var+Z.var+R.var))
#0.53 % of variance due to fixed effects, so just settlement here
## The proportion of variance due to random effects
(R2.random <- Z.var/(X.var+Z.var+R.var))
# 0.1278 % to sire, dam and interaction
## The proportion of variance due to residuals
(R2.resid <- R.var/(X.var+Z.var+R.var))
# 0.34 due to residuals
#VarCorr is random and residual.

VarCorr(nb8.glmm.gHQ.M) #VarCorr gives random effects and residuals
#dam*sire=0=0%
#sire=0 = 0%
#dam=0.43496==~37.956% of random effects (4.85%of total)
#residual=0.71098 ==~ 62.043% of random effects (7.929% of total)
#total=(0.71098+0.43496= 1.14594) 

## The conditional R2 (proportion of variance due to fixed and random effects)
(R2.conditional <- (X.var+Z.var)/(X.var+Z.var+R.var))
#0.6584823, 65.9 = 66% of variance due to fixed and random effects, fixed settlment is 53%, random is 13%, resid is 34%

#dam:sire variance=0
#sire variance=0
#fixed effect of Settlement signif.
#drop terms and refit

#set up deviation coding for dam identity
deviation<-with(Mortality_Competency, C(dam, sum, 7))
nb8.glm.devcode<-glm.nb(abundance ~ settlement*deviation, data = Mortality_Competency) #iteration level reached. 
summary(nb8.glm.devcode) 
#After dam included as fixed effect, settlement is no longer signficant alone, so drop

nb8.glm.zeros.devcode_nosett<-glm.nb(abundance ~ settlement:deviation+deviation, data = Mortality_Competency) #iteration level reached. 
summary(nb8.glm.zeros.devcode_nosett)
levels(deviation) #[1] "O3"  "O4"  "O5"  "O6"  "W10" "W11" "W5"

#Plotting model predictions from final model.
predframe<-with(Mortality_Competency, rbind(data.frame(expand.grid(deviation=levels(deviation),settlement=seq(min(settlement),max(settlement),length=100)))))
predframe$pred.response <- predict(nb8.glm.zeros.devcode_nosett,re.form=NA, newdata=predframe, type="response")
minmaxvals <- range(nb8.glm.zeros.devcode_nosett$settlement)
Figure3<-ggplot(predframe,aes(settlement,predframe$pred.response,color=deviation))+geom_line(size=1, aes(linetype=deviation))+
  scale_y_log10()+theme_classic()+scale_linetype_manual(values=c("solid","dashed","dotted","dotdash", "solid","dashed","dotted"))+
  scale_color_manual(values=c("black","black","black","black","grey60","grey60","grey60"))+
  labs(x="Number of settled larvae", y="Predicted log abundance of surviving juveniles")+
  scale_x_continuous(breaks=seq(0,500,50))+theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))

#################################
#Relationship between larval weight and juvenile survivorship in the field?
Weight_surv_corr_Competency<-read.table("Weight_surv_corr_Competency.csv", header=T, sep=",", strip.white=T)


#run lm_eqn function first at bottom of script

#Excluding larvae of familie that didn't settle (F26, F30, F21)
Weight_surv_corr_Competency_onlysurv<-Weight_surv_corr_Competency[1:39,]
tail(Weight_surv_corr_Competency_onlysurv) #confirm those larval from above listed families are gone
Supplementary.Figure1<-ggplot(Weight_surv_corr_Competency_onlysurv, aes(x=weight, y=Percmortality))+geom_point(size=4, aes(shape=dam))+theme_classic()+
  geom_smooth(method = "lm", se=FALSE, color="black")+ylab("Juvenile mortality (%)")+
  xlab("Larval dry weight (micrograms)")+geom_text(aes(x = 15, y = 110, label = lm_eqn(lm(Percmortality ~ weight, Weight_surv_corr_Competency_onlysurv))), parse = TRUE)+
  theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+scale_shape_manual(values=c(15,16,17,18,8,9,3))
Supplementary.Figure1 #R2=0.201

#Including larvae of families that didn't settle (F26, F30, F21)
Supplementary.Figure2<-ggplot(Weight_surv_corr_Competency, aes(x=weight, y=Percmortality))+geom_point(size=4, aes(color=dam))+
  theme_classic()+geom_smooth(method = "lm", se=FALSE, color="black")+ylab("Juvenile mortality (%)")+xlab("Larval dry weight (micrograms)")+
  geom_text(aes(x = 15, y = 110, label = lm_eqn(lm(Percmortality ~ weight, Weight_surv_corr_Competency))), parse = TRUE)+theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))
Supplementary.Figure2 #R2=0.328


#########################################################
#Symbiodinium community
#########################################################

###################
#All raw data uploaded to SRA
Bestworst_rawreads<-read.csv("Bestworst_rawreads.csv", row.names=1)
Bestworst_rawreads_phy<-otu_table(Bestworst_rawreads, taxa_are_rows=TRUE)

bestwortstsurv_metadata_Competency<-read.csv("bestwortstsurv_metadata_Competency.csv",row.names=1)
bestwortstsurv_metadata_Competency_phy<-sample_data(bestwortstsurv_metadata_Competency)

bestwortstsurv_taxonomy_Competency<-read.csv("bestwortstsurv_taxonomy_Competency.csv",row.names=1)
bestwortstsurv_taxonomy_Competency_phy<-as.matrix(bestwortstsurv_taxonomy_Competency)
bestwortstsurv_taxonomy_Competency_phy<-tax_table(bestwortstsurv_taxonomy_Competency_phy)

bestworst_merge_raw<-merge_phyloseq(Bestworst_rawreads_phy,bestwortstsurv_metadata_Competency_phy,bestwortstsurv_taxonomy_Competency_phy)
bestworst_merge_raw

famextract_subset_deseq<-phyloseq_to_deseq2(bestworst_merge_raw, ~family)
gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(famextract_subset_deseq), 1, gm_mean)
diagdds = estimateSizeFactors(famextract_subset_deseq, geoMeans = geoMeans)
famextract_subset_diagdds2 = DESeq(diagdds, fitType="local")

######now to specific comparisions
###F8 versus 3 best families
#f8xf1
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F8", "F1"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f8xf12
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F8", "F12"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f8xf4
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F8", "F4"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

###F14 versus 3 best families
#f14xf1
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F14", "F1"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f14xf4 
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F14", "F4"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f14xf12
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F14", "F12"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

###F17 versus 3 best families
#f17xf1- no signficant differences
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F17", "F1"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
head(sigtab2)

#f17xf4 
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F17", "F4"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f17xf12=no diff
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F17", "F12"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f18xf1
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F18", "F1"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f18xf4
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F18", "F4"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f18xf12
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F18", "F12"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f28xf1
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F28", "F1"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f28xf4
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F28", "F4"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f28xf12
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F28", "F12"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2


#healthy comparisons
#f1xf4
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F1", "F4"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

#f1xf12=no diff
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F1", "F12"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2


#f4xf12
reslevels<-results(famextract_subset_diagdds2, contrast=c("family", "F4", "F12"))
reslev<-as.data.frame(reslevels)
resadj<-reslevels[order(reslevels$padj, na.last=NA),]
resadj2<-as.data.frame(resadj)
alpha=0.05
sigtab2<-resadj2[(resadj2$padj < alpha),]#
sigtab2

######################
######################
#output values from above Deseq2 tests imported into excel and averaged to make summary figure below
logchanges<-read.csv("foldchanges_Competency.csv") 
ggplot(logchanges, aes(x=OTUs, y=average.fold.change.low.to.high))+geom_bar(stat="identity", color="black", aes(fill=clade))+geom_errorbar(aes(ymin=logchanges$average.fold.change.low.to.high-logchanges$se, ymax=logchanges$average.fold.change.low.to.high+logchanges$se),colour="black", width=.2)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_fill_manual(values=c("black", "grey21", "grey82","grey45", "grey100"))+ylim(-9, 15) #+ scale_fill_grey()
#colors, order and text labels changed in photoshop

######################
#Chao1 scores
surv.scores<-estimate_richness(bestworst_merge_raw, measures = c("Observed", "Chao1"))
dim(surv.scores)#55x5
(mean_healthy<-mean(surv.scores[c(1:25,32:42,55), 2])) #12.75946 Chao1 
(mean_healthy_se<-mean(surv.scores[c(1:25,32:42,55), 3])) #2.790429 se  
(mean_not_healthy<-mean(surv.scores[c(26:31,43:54), 2])) #16.43 chao1
(mean_not_healthy_se<-mean(surv.scores[c(26:31,43:54), 3])) #1.71 se

######################
#NMDS plot
#start with Deseq2 shiftlog normalized matrix of counts. This normalization was done with full dataset of all samples but concatenated the normalized counts down to just samples of interest
bestwortstsurv_norm_geom1_shiftlog_Competency<-read.csv("bestwortstsurv_norm_geom1_shiftlog_Competency.csv", row.names=1)
bestwortstsurv_norm_geom1_shiftlog_Competency_phy<-otu_table(bestwortstsurv_norm_geom1_shiftlog_Competency, taxa_are_rows=TRUE)

bestworst_merge_phyloseq<-merge_phyloseq(bestwortstsurv_norm_geom1_shiftlog_Competency_phy,bestwortstsurv_metadata_Competency_phy,bestwortstsurv_taxonomy_Competency_phy)

ord2<-plot_ordination(bestworst_merge_phyloseq, ordinate(bestworst_merge_phyloseq, "NMDS", "bray"), color= "family", shape="surv")
orddata2 <- ord2$data 

#manually click through each layer for proper functionality
df_ell2 <- data.frame()
for(g in levels(orddata2$family)){
  df_ell2 <- rbind(df_ell2, cbind(as.data.frame(with(orddata2[orddata2$family==g,], ellipse(cor(NMDS1,NMDS2), 
                                                                                            scale=c(sd(NMDS1),sd(NMDS2)), 
                                                                                            centre=c(mean(NMDS1),mean(NMDS2))))),family=g))
}
head(df_ell2)

Figure4.NMDS.color<-ggplot(data = orddata2, aes(NMDS1, NMDS2))+geom_point(aes(color=family), size=3)+scale_color_manual(values=c("black","darkblue","#ff6f4b","#da40ae","#fc0040","#d94f34","lightblue","#e19028"))+ 
  geom_path(data=df_ell2, aes(x=x, y=y, color=family), size=1)+theme_classic()+geom_vline(xintercept = 0, size=1)+geom_hline(yintercept = 0, size=1)+theme(axis.line.x = element_line(color = "black"), axis.line.y = element_line(color = "black"))+theme(legend.key.size=unit(1.5,"cm"))#+theme(panel.grid.major=element_line(color="black")) #+scale_linetype_manual(values=c("solid","dashed","solid","dashed","dotted","twodash", "dotted","dotdash"))+
Figure4.NMDS.color


####Adonis Permutational mutivariate analysis of variance
ord<-plot_ordination(onlyplastic, ordinate(onlyplastic, "NMDS", "bray"), color="surv", shape="family")
ord #this gives stress =0.1086899. stress type 1 week ties.
#Murray: stress value of the first run is very low (as is the stress of any one of the candidate configurations\\

df<-as(sample_data(bestworst_merge_phyloseq),"data.frame")
d <- phyloseq::distance(bestworst_merge_phyloseq, method="bray")
adonis(d ~ surv+family, df)
#R2= .178=17.8% of variance explained by survivorship.p=0.001
#More variance explained by family. R2=0.21. also signficant. p=0.001


#############################################################################################################
#SummarySE function
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

#################################################
lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}
