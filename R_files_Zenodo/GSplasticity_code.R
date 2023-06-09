# analyses and code by Maggie R. Wagner
# maggie.r.wagner@gmail.com

####### Clear workspace ########
rm(list=ls())
####### Color palettes #######
GLSPalette<-c('purple','forest green','grey','orange')
sitePalette<- c("darkcyan", "#000000", "#E69F00", "#56B4E9")
sspPalette<-c('red','black')

####### Load packages #######
library(maps)
library(mapdata)
library(reshape2)
library(dplyr)
library(vegan)
library(ggplot2)
library(lme4)
library(lmerTest)
library(lsmeans)
library(car)
library(ggmap)
library(grid)
library(scales)
library(tidyr)
library(stringr)
library(ggjoy)
library(multcomp)
library(cowplot)

####### in case MASS is also loaded, redefine select() to point to dplyr #######
select<-dplyr::select


####### Function: PVEranef() #######
PVEranef<-function(model) {
  RanEffects<-names(ranef(model))
  df<-data.frame('Term'=character(),'PVE'=numeric())
  for (i in 1:length(RanEffects)) {
    term<-RanEffects[i]
    df<-rbind(df,data.frame('Term'=term,
                            'PVE'=subset(as.data.frame(VarCorr(model)),grp==term)$vcov/sum(as.data.frame(VarCorr(model))$vcov)))
  }
  return(rand(model)$rand.table %>% as.data.frame() %>% 
           merge(.,df,by.x='row.names',by.y='Term'))
}

####### Function: r2.LMM() => estimate R^2 for a model as correlation between fitted and observed values #######
r2.LMM <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

####### Function: simplify.results() #######
simplify.results<-function(lmm) {
  # extract the dependent variable name for this model:
  response<-colnames(model.frame(lmm))[1]
  # convert ANOVA output to data frame:
  fixedresults<-anova(lmm) %>% as.data.frame() %>% 
    mutate(Term=row.names(.),Response=response) %>% # store the Term names and dependent variable name as columns
    select(-ends_with('Sq')) %>% # remove Sums of Squares / Mean Squares columns
    rename(F.or.ChiSq=F.value, p.value.raw=`Pr(>F)`) # rename columns for merging with random-effects results
  # convert random effects tests to dataframe: 
  randresults<-rand(lmm)$rand.table %>% as.data.frame() %>% 
    mutate(Term=row.names(.),Response=response) %>% # store the Term names and dependent variable name as columns
    mutate(DenDF=Inf) %>% 
    rename(F.or.ChiSq=Chi.sq, NumDF=Chi.DF,p.value.raw=p.value) # rename columns for merging with fixed-effects results
  return(rbind(fixedresults,randresults)) # return the combined results for this model
}
####### Function: getLSmeans() #######
getLSmeans<-function(lmm) {
  # extract the dependent variable name for this model:
  response<-colnames(model.frame(lmm))[1]
  # Get LS means for Site:
  lsm<-lsmeans::lsmeans(lmm,~Site) %>% summary() %>% as.data.frame() %>%
    mutate(Genotype='all',DevStage='all',Term='Site',Response=response)
  # Add LS means for Genotype:
  lsmeans::lsmeans(lmm,~Genotype) %>% summary() %>% as.data.frame() %>%
    mutate(Site='all',DevStage='all',Term='Genotype',Response=response) %>%
    rbind(., lsm)->lsm
  # Add LS means for Genotype*Site:
  lsmeans::lsmeans(lmm,~Genotype:Site) %>% summary() %>% as.data.frame() %>%
    mutate(DevStage='all',Term='Genotype:Site',Response=response) %>%
    rbind(., lsm)->lsm
  # Add LS means for Developmental Stage:
  lsmeans::lsmeans(lmm,~DevStage) %>% summary() %>% as.data.frame() %>%
    mutate(Genotype='all',Site='all',Term='DevStage',Response=response) %>%
    rbind(., lsm)->lsm
  ### For sqrt.totalGS, convert estimates back to Total [GS] (real units)
  if (response=='sqrt.totalAGS') {
    lsm$lsmean<-lsm$lsmean^2
    lsm$lower.CL<-lsm$lower.CL^2
    lsm$upper.CL<-lsm$upper.CL^2
    lsm$SE<-lsm$SE^2
    lsm$Response='totalAGS' # note that the values have been changed (remove the square root)
  }
  return(lsm)
}
####### Function: getLSmeans.sqrt() #######
## For getting LS means of sqrt-transformed response variables
## define function to correctly square negative estimates:
negsq<-function(number) {
  if (number<0) { numsq <- -(number^2)}
  else {numsq <- number^2 }
  return(numsq)
}

getLSmeans.sqrt<-function(lmm) {
  # extract the dependent variable name for this model:
  response<-colnames(model.frame(lmm))[1]
  # Get LS means for Site:
  lsm<-lsmeans::lsmeans(lmm,~Site) %>% summary() %>% as.data.frame() %>%
    mutate(Genotype='all',DevStage='all',Term='Site',Response=response)
  # Add LS means for Genotype:
  lsmeans::lsmeans(lmm,~Genotype) %>% summary() %>% as.data.frame() %>%
    mutate(Site='all',DevStage='all',Term='Genotype',Response=response) %>%
    rbind(., lsm)->lsm
  # Add LS means for Genotype*Site:
  lsmeans::lsmeans(lmm,~Genotype:Site) %>% summary() %>% as.data.frame() %>%
    mutate(DevStage='all',Term='Genotype:Site',Response=response) %>%
    rbind(., lsm)->lsm
  # Add LS means for Developmental Stage:
  lsmeans::lsmeans(lmm,~DevStage) %>% summary() %>% as.data.frame() %>%
    mutate(Genotype='all',Site='all',Term='DevStage',Response=response) %>%
    rbind(., lsm)->lsm
  ### Convert estimates back to micromoles/mg (real units)
  lsm$upper.SE<-sapply((lsm$lsmean+lsm$SE),negsq) # square after adding S.E.
  lsm$lower.SE<-sapply((lsm$lsmean-lsm$SE),negsq) # square after subtracting S.E.
  lsm$lsmean<-sapply((lsm$lsmean),negsq)
  lsm$lower.CL<-sapply((lsm$lower.CL),negsq)
  lsm$upper.CL<-sapply((lsm$upper.CL),negsq)
  lsm$SE<-sapply((lsm$SE),negsq)
  return(lsm)
}
####### Make subdirectories #######
system('mkdir figures')
system('mkdir tables')
system('mkdir intermediate_data')

####### Session Info #######
sessionInfo()
"R version 3.3.2 (2016-10-31)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X Yosemite 10.10.5

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] cowplot_0.9.1      multcomp_1.4-8     TH.data_1.0-8      MASS_7.3-47        survival_2.41-3    mvtnorm_1.0-6     
[7] ggjoy_0.4.0        ggridges_0.4.1     stringr_1.2.0      tidyr_0.7.2        scales_0.5.0.9000  ggmap_2.6.1       
[13] car_2.1-6          lsmeans_2.27-61    lmerTest_2.0-36    lme4_1.1-14        Matrix_1.2-12      ggplot2_2.2.1.9000
[19] vegan_2.4-5        lattice_0.20-35    permute_0.9-4      dplyr_0.7.4        reshape2_1.4.3     mapdata_2.2-6     
[25] maps_3.2.0        

loaded via a namespace (and not attached):
[1] splines_3.3.2       Formula_1.2-2       assertthat_0.2.0    sp_1.2-5            latticeExtra_0.6-28
[6] backports_1.1.1     quantreg_5.34       glue_1.2.0          digest_0.6.13       RColorBrewer_1.1-2 
[11] checkmate_1.8.5     minqa_1.2.4         colorspace_1.3-2    sandwich_2.4-0      htmltools_0.3.6    
[16] plyr_1.8.4          pkgconfig_2.0.1     SparseM_1.77        purrr_0.2.4         xtable_1.8-2       
[21] jpeg_0.1-8          MatrixModels_0.4-1  htmlTable_1.11.0    tibble_1.3.4        mgcv_1.8-22        
[26] nnet_7.3-12         lazyeval_0.2.1      proto_1.0.0         pbkrtest_0.4-7      magrittr_1.5       
[31] estimability_1.2    nlme_3.1-131        foreign_0.8-69      tools_3.3.2         data.table_1.10.4-3
[36] geosphere_1.5-7     RgoogleMaps_1.4.1   munsell_0.4.3       cluster_2.0.6       bindrcpp_0.2       
[41] rlang_0.1.6.9003    nloptr_1.0.4        rstudioapi_0.7      rjson_0.2.15        htmlwidgets_0.9    
[46] base64enc_0.1-3     gtable_0.2.0        codetools_0.2-15    R6_2.2.2            gridExtra_2.3      
[51] zoo_1.8-0           knitr_1.17          bindr_0.1           Hmisc_4.0-3         stringi_1.1.6      
[56] parallel_3.3.2      Rcpp_0.12.14        mapproj_1.2-5       png_0.1-7           rpart_4.1-11       
[61] acepack_1.4.1       coda_0.19-1        "
####### Load Raw Data #######
allcensus<-as.data.frame(read.table('raw_data/all_census_clean.txt',sep='\\t',header=TRUE)) # cleaned data from all plants, all censuses 
allpeaks<-as.data.frame(read.table('raw_data/allpeaks_Field2014.txt',sep='\\t',header=TRUE)) # chromatograph peaks
plants<-as.data.frame(read.table('raw_data/plants_wFitness.txt',sep='\\t',header=TRUE))  # planting design / genotype info for each plant, with some summary phenotype data from entire summer
HPLC<-as.data.frame(read.table('raw_data/Field2014_HPLC.txt',sep='\\t',header=TRUE)) # guide to linking HPLC sample IDs to unique Plant IDs
ssp<-as.data.frame(read.table('raw_data/ssp_key.txt',sep='\\t',header=TRUE)) # key linking genotypes to subspecies
blks<-as.data.frame(read.table('raw_data/block_data.txt',sep='\\t',header=TRUE)) # environmental data on each block
merge(plants,ssp,by='Genotype') -> plants # add subspecies information

####### Process HPLC data #######
# consolidate chromatograph peak data into a data frame:
AGS<-reshape(allpeaks[,-1],idvar="HPLC_ID",v.names='Area',timevar='GS',direction='wide') # reshape data frame from melted to wide format (and toss out 'row' column)
AGS[is.na(AGS)]<-0 # replace NA with 0 for peak areas
AGS<-merge(AGS,HPLC,by='HPLC_ID',suffixes=c('','.y')) # merge with HPLC metadata
select(plants,Plant_ID,Genotype,Ssp,Block) %>% # merge with selected plant data
  merge(.,AGS,by='Plant_ID',suffixes=c('','.y')) -> AGS 

####### Standardize AGS data and calculate emergent properties of AGS profiles #######
## Standardize raw peak areas by sinigrin & sample weight ## 
## Relative response factors taken from Clarke et al. 2010, and C. Olson-Manning (personal communication)
amtSinigrin<-0.05 # amount of sinigrin in each sample (micromol)
# units of GS concentrations: micromol per mg dry weight

AGS<-filter(AGS,Area.sinigrin>0) %>% # Kick out any data points with no Sinigrin peak
  # standardize raw peak areas by sample weight and internal standard to get absolute concentrations (micromol/mg):
  mutate(conc1ME=(amtSinigrin*Area.1ME)/(Area.sinigrin*Weight*1), # relative response factor = 1
         conc1MP=(amtSinigrin*Area.1MP)/(Area.sinigrin*Weight*1), # relative response factor = 1
         conc6MSOH=(amtSinigrin*Area.6MSOH)/(Area.sinigrin*Weight*1), # relative response factor = 1
         conc2OH1ME=(amtSinigrin*Area.2OH1ME)/(Area.sinigrin*Weight*1.32)) %>% # relative response factor = 1.32
  # calculate emergent properties Total [AGS] and BC-ratio from absolute concentrations :
  mutate(totalAGS=conc1ME+conc1MP+conc2OH1ME+conc6MSOH, # sum of all aliphatic GS
         sqrt.totalAGS=sqrt(totalAGS), # square-root transform to improve homoscedasticity in downstream ANOVA
         # calculate BC-ratio based on peak areas (the Weight factor cancels out):
         # Doing it this way because BC-ratio can still be calculated for the 2 samples with missing Weight data:
         BCrat=((Area.2OH1ME/1.32)+Area.1ME+Area.1MP)/((Area.2OH1ME/1.32)+Area.1ME+Area.1MP+Area.6MSOH)) # divide by 1.32 for only compound with RRF != 1

### Calculate branched-chain diversity ###
# Because this is within-sample diversity, all peaks in each sample standardized by same weight & internal standard
# Therefore, can use just peak areas + relative response factors to calculate diversity
# Do it this way so that we can get BC-diversity even for samples with missing Weight data

AGS<-mutate(AGS,RRF.1ME=Area.1ME/1,RRF.1MP=Area.1MP/1,RRF.2OH1ME=Area.2OH1ME/1.32) # standardize peak areas by RRF only
## Calculate BC-diversity: 
AGS$BCdiv<-diversity(select(AGS,RRF.1ME,RRF.1MP,RRF.2OH1ME),index='shannon') # use Shannon index
## Remove temporary RRF-normalized columns:
AGS<-select(AGS,-RRF.1ME,-RRF.1MP,-RRF.2OH1ME)

####### Combine AGS data with census data #######
## Goal: get the height and developmental stage of each individual plant at the time the leaf sample was taken for HPLC

# Change format of tissue collection dates to match dates of censuses:
levels(AGS$Harvested)<-gsub('-Jun-14','jun2014',levels(AGS$Harvested))
AGS$Harvested<-factor(as.Date(AGS$Harvested,"%d%b%Y"))

AGS<-mutate(AGS,obsID=paste0(Plant_ID,'.',Harvested)) # make new column: 'obsID' to link to unique observation ID from census data (one plant, one date)

# Merge:
AGS<-merge(AGS,select(allcensus,obsID,H,Stg),by='obsID',all.x=TRUE) %>% rename(DevStage=Stg) # get height and developmental stage from that date

# consolidate all reproducing plants to same developmental stage for analysis
AGS$DevStage<-factor(plyr::mapvalues(AGS$DevStage,from=c('FO','FS1','FS2','SO'),to=c('FF','FF','FF','FF')))
AGS$DevStage<-ordered(AGS$DevStage, levels = c('R', 'B', 'FF')) # order the developmental stages

####### Save clean AGS data #######
save(AGS,file='intermediate_data/AGS_data_clean.RData')

####### Make cartoon of 25 genotypes randomly arrayed in a 10x5 plot #######
dat <- as.data.frame(matrix(ncol=2, nrow=50))
colnames(dat)<-c('x','y')
dat$x<-c(1:10) # set x coordinates
arrange(dat,x)->dat
dat$y<-c(1:5) # set y coordinates

pdf(file='figures/Figure1b.pdf',width=10,height=5)
plot(dat$y~dat$x,type='n')
points(x=dat$x,y=dat$y,pch=sample(c(1:25,1:25),replace=FALSE), # use 25 different shapes
       cex=2.5,bg='grey',lwd=3)
dev.off()
rm(dat) # clean up
####### Supp. Figure 4: block environmental data #######
## Summarize vegetation density and herbivory intensity for each block:
mutate(plants,RosDamPC=(RosDLN/RosLN)*RosAD/100) %>% # First, calculate rosette damage percentage
  select(Block,Veg,RosDamPC) %>% group_by(Block) %>% # Calculate mean vegetation cover, rosette herbivory of each block
  summarize(BlkVeg=mean(Veg,na.rm=TRUE),BlkHerbivory=mean(RosDamPC,na.rm=TRUE)) %>% 
  ungroup %>% as.data.frame() -> veg
# Plot everything:
pdf(file='figures/Supp_Fig4.pdf',width=9,height=7)
merge(veg,blks) %>%
  rename(Veg=BlkVeg,PlantDiv=Div,Herbivory=BlkHerbivory) %>%
  gather(key='variable',value='value',pH,conductivity,NO3_N,P,K,Ca,Mg,S,Na,Fe,Zn,Mn,Cu,PlantDiv,Veg,Herbivory) %>%
  ggplot(.,aes(x=1,color=Site,y=value))+
  facet_wrap(~variable,scales='free')+
  geom_boxplot(varwidth=TRUE)+
  scale_color_manual(values=sitePalette)+
  guides(color=guide_legend(ncol=4))+
  theme_classic()+
  theme(legend.position='top',legend.text=element_text(size=16))+
  theme(legend.title=element_text(face='bold',size=20))+
  theme(axis.title=element_blank(),axis.text.x=element_blank())+
  theme(axis.text.y = element_text(size=14))+
  theme(strip.background = element_rect(fill='grey90',color='grey90'))+
  theme(strip.text=element_text(face='bold',size=18))+
  theme(legend.background = element_rect(fill='grey90'))
dev.off()

####### How complete is the replication of genotypes among sites? #######
table(AGS$Site,AGS$Genotype)
# Only one cell with 0 observations: genotype "Bay Horse Saddle" (BHS) at Jam 
### Need to remove BHS before ANOVA, however, can keep BHS for reaction
### norm analysis (because there are still enough observations in blocks at
### other sites to fit a reaction norm) and phenotypic selection analysis.
### Remove BHS after making copy of df with BHS:
AGS.withBHS<-AGS
AGS.noBHS<-filter(AGS,Genotype!='BHS') %>% mutate(Genotype=factor(Genotype))
table(AGS.noBHS$Site,AGS.noBHS$Genotype) # no more empty cells
####### Figure 2b: Glucosinolate profile bar chart #######
pdf(file="figures/Figure2b.pdf",width=11,height=5)
arrange(AGS.withBHS,BCrat) %>%
  within(.,Plant_ID<-as.ordered(factor(Plant_ID,levels=Plant_ID[order(BCrat)]))) %>%
  melt(.,measure.vars=c('conc2OH1ME','conc1MP','conc1ME','conc6MSOH'),
       variable.name="GStype",value.name="conc",na.rm=TRUE) %>%
  within(.,GStype<-as.ordered(factor(GStype, 
                                     levels=c('conc6MSOH','conc1ME','conc1MP','conc2OH1ME')))) %>% 
  ggplot(.,aes(x=Plant_ID,y=conc/totalAGS,fill=GStype))+
  facet_wrap(~Ssp,scales='free')+
  geom_bar(stat='identity',width=1)+
  scale_fill_manual(values=c('orange','grey', "forest green","purple"),
                    labels=c('6MSOH','1ME', '1MP  ','2OH1ME  '),
                    name='')+
  scale_y_continuous(labels = scales::percent)+ylab("")+
  theme_classic()+
  theme(axis.text.x = element_blank(),axis.ticks = element_blank(),
        axis.title.x=element_text(size=24,face='bold'),
        axis.title.y=element_text(size=24,face='bold'),
        axis.text.y=element_text(size=20),
        legend.text=element_text(size=22),legend.position='top')
dev.off()
####### Figure 2c: Glucosinolate profile bubble chart #######
pdf(file="figures/Figure2c.pdf",width=11,height=5)
ggplot(AGS.withBHS,aes(x=BCrat,y=BCdiv,color=Ssp,size=totalAGS))+
  geom_point(alpha=0.4)+
  scale_color_manual(values=c('red','black'),guide=FALSE)+
  scale_size_continuous(name="Total\\n[AGS]",guide=FALSE)+
  xlab("BC-ratio")+ylab("BC-diversity")+
  theme_classic()+
  theme(axis.title.y=element_text(size=24,face='bold'))+
  theme(axis.title.x=element_text(size=24,face='bold'))+
  theme(axis.text.x=element_text(size=20))+
  theme(axis.text.y=element_text(size=20))+
  theme(legend.title=element_text(size=22,face='bold'))+
  theme(legend.text=element_text(size=20))
dev.off()

####### ------ PARTITIONING VARIANCE IN GLUCOSINOLATE PROFILES ------ #######
## Note: Block is implicitly nested within Site (by unique names for all levels of Block factor)
####### REML linear mixed model for Total [AGS] #######
GxE.Site.totalAGS<-lmer(sqrt.totalAGS~Genotype*Site+DevStage+H+(1|Block)+(1|Block:Genotype)+(1|Batch),data=AGS.noBHS)
qqnorm(resid(GxE.Site.totalAGS));qqline(resid(GxE.Site.totalAGS))
plot(resid(GxE.Site.totalAGS)~fitted(GxE.Site.totalAGS)) 
hist(resid(GxE.Site.totalAGS),50) # looks good

r2.LMM(GxE.Site.totalAGS) # 0.490

anova(GxE.Site.totalAGS) # Test fixed effects

PVEranef(GxE.Site.totalAGS) # Test random effects and % variance explained:

####### Do we get the same results excluding the total [AGS] outlier SAD12? ####
# YES
GxE.Site.totalAGS.noSAD<-lmer(sqrt.totalAGS~Genotype*Site+DevStage+H+(1|Block)+(1|Block:Genotype)+(1|Batch),data=filter(AGS.noBHS,Genotype!='SAD'))
anova(GxE.Site.totalAGS.noSAD)
PVEranef(GxE.Site.totalAGS.noSAD) # Test random effects and % variance explained:

####### REML linear mixed model for BC-ratio #######
GxE.Site.BCrat<-lmer(BCrat~Genotype*Site+DevStage+H+(1|Block)+(1|Block:Genotype)+(1|Batch),data=AGS.noBHS)
qqnorm(resid(GxE.Site.BCrat));qqline(resid(GxE.Site.BCrat)) # 
plot(resid(GxE.Site.BCrat)~fitted(GxE.Site.BCrat)) # tails but no obvious heteroscedasticity
hist(resid(GxE.Site.BCrat),50)

r2.LMM(GxE.Site.BCrat) # 0.989

anova(GxE.Site.BCrat) # Test fixed effects:

PVEranef(GxE.Site.BCrat) # Test random effects and % variance explained:

####### REML linear mixed model for BC-diversity #######
GxE.Site.BCdiv<-lmer(BCdiv~Genotype*Site+DevStage+H+(1|Block)+(1|Block:Genotype)+(1|Batch),data=AGS.noBHS)
qqnorm(resid(GxE.Site.BCdiv));qqline(resid(GxE.Site.BCdiv))
plot(resid(GxE.Site.BCdiv)~fitted(GxE.Site.BCdiv)) 
hist(resid(GxE.Site.BCdiv),50)
r2.LMM(GxE.Site.BCdiv) # 0.89
anova(GxE.Site.BCdiv) # Test fixed effects:
PVEranef(GxE.Site.BCdiv) # Test random effects and % variance explained:

####### Combine results of all REML models and adjust p-values for 3 comparisons: #######
rbind(simplify.results(GxE.Site.totalAGS),
      simplify.results(GxE.Site.BCrat),
      simplify.results(GxE.Site.BCdiv)) -> GSvar.results
## Correct p-values for each term:
group_by(GSvar.results,Term) %>%
  mutate(p.value.adj=p.adjust(p.value.raw,'holm')) %>%
  ungroup %>% as.data.frame -> GSvar.results

## Save results:
save(GSvar.results,file='tables/stats_GSvar_models.RData')

## Save as table:
mutate(GSvar.results,DenDF=round(DenDF),F.or.ChiSq=round(F.or.ChiSq,digits=2)) %>% 
  select(-p.value.raw) %>% 
  mutate(p.value.adj=signif(p.value.adj,digits=3)) %>% 
  write.table(file='tables/stats_GSvar_models.txt',sep='\\t',row.names=FALSE,col.names=TRUE)
####### Extract LS means from REML linear mixed models#######
# Combine all LS means into one data frame:
## Start with BC-ratio and BC-diversity:
rbind(getLSmeans(GxE.Site.BCrat),getLSmeans(GxE.Site.BCdiv)) %>%
  mutate(upper.SE=lsmean+SE,lower.SE=lsmean-SE) -> LSM.all

## Back-transform Total [AGS] estimates to be on same scale as raw measurements and add to dataframe:
getLSmeans(GxE.Site.totalAGS) %>%
  mutate(upper.SE=sapply((lsmean+SE),negsq), # see source file for negsq()
         lower.SE=sapply((lsmean-SE),negsq), 
         lsmean=sapply(lsmean,negsq),
         lower.CL=sapply(lower.CL,negsq),
         upper.CL=sapply(upper.CL,negsq),
         SE=sapply(SE,negsq),
         Response='totalAGS') %>% ### Rename Total [AGS] to indicate values have been back-transformed 
rbind(.,LSM.all) -> LSM.all

### Merge with subspecies data and re-order by subspecies
full_join(LSM.all,ssp,by='Genotype') %>% 
  arrange(Ssp) %>%
  mutate(Genotype=factor(Genotype,Genotype))  -> LSM.all

### Save LS means
save(LSM.all,file='intermediate_data/LSmeans_all_GxE.RData')

####### Figure 3a: Main-effect LS means #######
pdf(file='figures/Figure3a.pdf',width=9,height=5)
filter(LSM.all,Term=='Genotype') %>% 
  mutate(Response=plyr::mapvalues(Response,from=c('BCrat','BCdiv','totalAGS'),
                                  to=c(' BC-ratio','BC-diversity',' Total [AGS]'))) %>%
  ggplot(.,aes(x=Genotype,y=lsmean,color=Ssp,alpha=Ssp))+
  facet_wrap(~Response,ncol=3,scales='free')+
  geom_point(size=2,position=position_dodge(width=1))+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0.0,size=.5,position=position_dodge(width=1))+
  theme_classic()+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_blank())+
  theme(axis.title.y=element_blank(),axis.text.y=element_text(size=22))+
  theme(legend.title= element_text(size=30),legend.text=element_text(size=26))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(2.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))+
  theme(strip.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))+
  geom_point(data=filter(LSM.all,Term=='Site') %>% 
               mutate(Response=plyr::mapvalues(Response,from=c('BCrat','totalAGS','BCdiv'),
                                               to=c(' BC-ratio',' Total [AGS]','BC-diversity')))%>%
               mutate(Genotype=plyr::mapvalues(Site,from=c('Ald','Jam','Mah','Sil'),to=c('EGM','MIL','BHM','THT'))),
             aes(x=Genotype,y=lsmean,color=Site,alpha=Site),size=4,shape=17)+
  geom_errorbar(data=filter(LSM.all,Term=='Site') %>% 
                  mutate(Response=plyr::mapvalues(Response,from=c('BCrat','totalAGS','BCdiv'),
                                                  to=c(' BC-ratio',' Total [AGS]','BC-diversity')))%>%
                  mutate(Genotype=plyr::mapvalues(Site,from=c('Ald','Jam','Mah','Sil'),to=c('EGM','MIL','BHM','THT'))),
                aes(x=Genotype,ymin=lower.CL,ymax=upper.CL,color=Site),width=0.0,size=1.5,position=position_dodge(width=1))+
  scale_alpha_manual(values=c(1,0.3,1,1,1,0.3),guide=FALSE)+
  scale_color_manual(values=c(sitePalette[1],"red",sitePalette[2:4],"black"),guide=FALSE)
dev.off()
####### Figure 3b: Genotype*Site LS means #######
pdf(file='figures/Figure3b.pdf',width=9,height=5)
filter(LSM.all,Term=='Genotype:Site') %>% 
  mutate(Response=plyr::mapvalues(Response,from=c('BCrat','BCdiv','totalAGS'),
                                  to=c(' BC-ratio','BC-diversity',' Total [AGS]'))) %>%
  ggplot(.,aes(x=Genotype,y=lsmean,color=Site))+
  facet_wrap(~Response,ncol=1,scales='free')+
  geom_point(position=position_dodge(width=.9))+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0,position=position_dodge(width=.9))+
  geom_vline(xintercept=0.5+1:24,linetype='dotted',color='grey')+
  scale_color_manual(values=sitePalette,guide=FALSE)+
  theme_classic()+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_blank())+
  theme(axis.title.y=element_blank(),axis.text.y=element_text(size=14))
dev.off()
####### Supp. Figure 7: Developmental Stage #######
pdf(file='figures/Supp_Fig7.pdf',width=9,height=7)
filter(LSM.all,Term=='DevStage') %>% 
  mutate(DevStage=ifelse(DevStage=='R',' rosette', # rename stages for nicer graph
                         ifelse(DevStage=='B','bolting',
                                'fruiting /\\nflowering'))) %>%
  mutate(Response=ifelse(Response=='BCrat',' BC-ratio', # rename traits
                         ifelse(Response=='BCdiv','BC-diversity',' Total [AGS]'))) %>%
  ggplot(.,aes(x=DevStage,y=lsmean))+
  facet_wrap(~Response,ncol=3,scales='free')+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0,size=1.5)+
  ylab("Least-squares mean (95% CI)")+
  xlab("Developmental stage")+
  theme_classic()+
  theme(strip.background=element_rect(fill='grey90'),strip.text=element_text(size=18))+
  theme(axis.title=element_text(size=18),axis.text.x=element_text(size=12,angle=30,vjust=0.8,hjust=0.8))
dev.off()

####### Check for genetic correlation between traits using Genotype LS means #######
## Extract all Genotype LS means
filter(LSM.all,Term=='Genotype') %>% 
  select(Genotype,lsmean,Response,Ssp) %>%
  spread(key=Response,value=lsmean) -> genoLSM

# subset genotypes: only those with branched-chain functionality ("gain of function allele only"):
filter(genoLSM,!Genotype%in%c('MAH','PAR','SAD')) -> genoLSM.GOFonly

## Spearman rank genetic correlations: all genotypes included
cor.test(~BCdiv+BCrat,data=genoLSM,method='spearman') 
cor.test(~totalAGS+BCrat,data=genoLSM,method='spearman') 
cor.test(~totalAGS+BCdiv,data=genoLSM,method='spearman') 

p.adjust(c(0.02672,0.02283,0.00249),method='holm') # adjust P-values

## Spearman rank genetic correlations: excluding genotypes with no branched-chain functionality
cor.test(~BCdiv+BCrat,data=genoLSM.GOFonly,method='spearman') 
cor.test(~totalAGS+BCrat,data=genoLSM.GOFonly,method='spearman') 
cor.test(~totalAGS+BCdiv,data=genoLSM.GOFonly,method='spearman') 

p.adjust(c(0.4185,0.05792,0.004684),method='holm') # correct p-values

####### Extract block BLUPs from REML linear mixed models #######
### BLUPs for Blocks (mean for all genotypes):
# Start with BC-ratio
ranef(GxE.Site.BCrat)$Block %>% as.data.frame %>% 
  mutate(Block=row.names(.)) %>% 
  rename(BLUP.BCrat=`(Intercept)`) -> blockBLUPs

# Now do Total AGS
## Keep sqrt-transformed values for now- will un-transform after adding BLUPs to LS means
ranef(GxE.Site.totalAGS)$Block %>% as.data.frame %>% 
  mutate(Block=row.names(.)) %>% 
  rename(BLUP.sqrt.totalAGS=`(Intercept)`) %>%
  merge(.,blockBLUPs,by='Block',suffixes=c('','.y')) -> blockBLUPs

# Now BC diversity
ranef(GxE.Site.BCdiv)$Block %>% as.data.frame %>% 
  mutate(Block=row.names(.)) %>% 
  rename(BLUP.BCdiv=`(Intercept)`) %>%
  merge(.,blockBLUPs,by='Block',suffixes=c('','.y')) -> blockBLUPs

## Add Site info for each block
blockBLUPs<-merge(blockBLUPs,select(blks,Block,Site),by='Block')

## Combine BLUPs with site LS means
filter(LSM.all,Term=='Site') %>% # start with LS means
  unique() %>% select(Site,lsmean,Response) %>%
  spread(key=Response,value=lsmean) %>% 
  merge(.,blockBLUPs,by='Site',suffixes=c('','.y'),all=T) %>%  # merge with BLUPs
  ## add BLUPs to LS means 
  mutate(est.BCrat=BCrat+BLUP.BCrat,
         est.totalAGS=(sqrt(totalAGS)+BLUP.sqrt.totalAGS)^2, # temporarily re-transform (sqrt transform) the Total [AGS] LS mean to put it on same scale as BLUP, then square the sum to put on original scale
         est.BCdiv=BCdiv+BLUP.BCdiv) %>%
  ## Get rid of extra columns:
  select(Site,Block,starts_with('est'))-> blockBLUPs

## Save BLUPs
save(blockBLUPs,file='intermediate_data/block_BLUPs.RData')

####### Extract Genotype*Block BLUPs from REML linear mixed models #######
### BLUPs for Geno*Block:
# Start with BC-ratio
ranef(GxE.Site.BCrat)$`Block:Genotype` %>% as.data.frame %>% 
  mutate(BlockGeno=row.names(.)) %>% 
  rename(BLUP.BCrat=`(Intercept)`) -> genoblockBLUPs

# Now do Total AGS
ranef(GxE.Site.totalAGS)$`Block:Genotype` %>% as.data.frame %>% 
  mutate(BlockGeno=row.names(.)) %>% 
  rename(BLUP.sqrt.totalAGS=`(Intercept)`) %>%
  merge(.,genoblockBLUPs,by='BlockGeno',suffixes=c('','.y')) -> genoblockBLUPs

# Now BC diversity
ranef(GxE.Site.BCdiv)$`Block:Genotype` %>% as.data.frame %>% 
  mutate(BlockGeno=row.names(.)) %>% 
  rename(BLUP.BCdiv=`(Intercept)`) %>%
  merge(.,genoblockBLUPs,by='BlockGeno',suffixes=c('','.y')) -> genoblockBLUPs

## Add Site info for each block
mutate(genoblockBLUPs,Block=substr(BlockGeno,start=1,stop=4)) %>% # extract Block number
  mutate(Genotype=substr(BlockGeno,start=6,stop=8)) %>% # extract Genotype 
  merge(.,select(blks,Block,Site),by='Block') -> genoblockBLUPs 

## Merge with subspecies info:
merge(genoblockBLUPs,ssp,by='Genotype') -> genoblockBLUPs

## Back-transform Total [AGS] BLUPs
mutate(genoblockBLUPs,BLUP.totalAGS=sapply(BLUP.sqrt.totalAGS,negsq)) -> genoblockBLUPs

## Get rid of extra columns:
select(genoblockBLUPs,Site,Block,Genotype,Ssp,starts_with('BLUP.'),-BLUP.sqrt.totalAGS)-> genoblockBLUPs

## Save genotype-block BLUPs:
save(genoblockBLUPs,file='intermediate_data/genoblock_BLUPs.RData')

####### Supp. Figure 7a: Block BLUPs #######
pdf(file='figures/Supp_Fig7a.pdf',width=11,height=5)
ggplot(blockBLUPs,aes(x=est.totalAGS,y=est.BCrat,color=Site))+
  geom_point(size=2,alpha=0.7)+
  scale_color_manual(values=sitePalette)+
  theme_classic()+
  xlab("Total [AGS] (block BLUPs)")+
  ylab("BC-ratio (block BLUPs)")+
  theme(axis.title.y=element_text(size=21),axis.title.x=element_text(size=21))+
  theme(axis.text.y=element_text(size=15),axis.text.x=element_text(size=15))+
  theme(legend.background=element_rect(fill='grey90'))+
  theme(legend.title=element_text(size=21),legend.text=element_text(size=15))
dev.off()
####### Supp. Figure 7b: Geno*Block BLUPs #######
# order by subspecies:
arrange(genoblockBLUPs,Ssp) %>%
  mutate(Genotype=factor(Genotype,Genotype)) -> genoblockBLUPs 

pdf(file='figures/Supp_Fig7b.pdf',width=11,height=7)
melt(genoblockBLUPs,measure.vars=c('BLUP.BCrat','BLUP.totalAGS','BLUP.BCdiv')) %>%
  mutate(variable=plyr::mapvalues(variable,from=c('BLUP.BCdiv','BLUP.BCrat','BLUP.totalAGS'),
                                  to=c('BC-\\ndiversity',' BC-\\nratio',' Total\\n[AGS]'))) %>%
  ggplot(.,aes(x=Genotype,y=value,color=Ssp))+
  facet_grid(variable~Site,scales='free')+
  geom_boxplot()+
  scale_color_manual(values=sspPalette,guide=FALSE)+
  ylab("Block deviation from mean trait value")+
  theme_classic()+
  theme(axis.title.y=element_text(size=21),axis.title.x=element_text(size=21))+
  theme(axis.text.y=element_text(size=15))+
  theme(axis.text.x=element_blank())+
  theme(panel.spacing=unit(1.5,'lines'))+
  theme(strip.background = element_rect(fill='grey90'),strip.text=element_text(size=21))
dev.off()

####### ------ PHENOTYPIC SELECTION ANALYSIS - LINEAR ------ #######
####### Summarize fitness #######
# What proportion of plants survived the winter?
table(plants$Site,plants$Survival)
"    N   Y
Ald 457 643
Jam 445 355
Mah 353 747
Sil 290 710     "
# Ald = 58.5% survival
# Jam = 44.4%
# Mah = 67.9%
# Sil = 71.0%

# What proportion of surviving plants reproduced?
filter(plants,Survival=='Y') %>% select(Survival,FrProd,Site) %>%
  mutate(FrProd=factor(ifelse(FrProd==0,'no_rep','reproduced'))) %>%
  table()
# Ald = 233/(410+233) = 36.2%
# Jam = 78/(277+78) = 22.0%
# Mah = 429/(318+429) = 57.4%
# Sil = 424/(286+424) = 59.7%

####### Prepare data for phenotypic selection analysis #######
## Combine glucosinolate data with fitness data from censuses:
## Exclude ambiguous fitness records (a few census records were ambiguous and couldn't be reconciled with previous census data- exclude these)
fitAGS <- select(plants,Plant_ID,Ambiguous,FrProd,Survival) %>% 
  merge(.,AGS.withBHS,by='Plant_ID') %>% # We can include BHS genotype for this- empty cells won't be a problem
  filter(Ambiguous=='no') # kick out 12 ambiguous fitness records

dim(fitAGS) # total of 1494 data points for phenotypic selection analysis

## Calculate individuals' relative fecundity within each site
group_by(fitAGS,Site) %>% mutate(relFrProd=FrProd/mean(FrProd,na.rm=TRUE)) -> fitAGS
### Mean relative fecundity = 1 in each site

## Standardize traits (mean=0,sd=1) within each site
group_by(fitAGS,Site) %>% 
  mutate(std.BCrat=(BCrat-mean(BCrat,na.rm=TRUE))/sd(BCrat,na.rm=TRUE),
         std.totalAGS=(totalAGS-mean(totalAGS,na.rm=TRUE))/sd(totalAGS,na.rm=TRUE),
         std.sqrt.totalAGS=(sqrt.totalAGS-mean(sqrt.totalAGS,na.rm=TRUE))/sd(sqrt.totalAGS,na.rm=TRUE),
         std.BCdiv=(BCdiv-mean(BCdiv,na.rm=TRUE))/sd(BCdiv,na.rm=TRUE)) %>%
  ungroup %>% as.data.frame -> fitAGS

####### Phenotypic selection analysis-- Estimate linear selection DIFFERENTIALS at each site ####### 
for(i in 1:4) {
  site<-c('Ald','Jam','Mah','Sil')[i]
  site.subset<-filter(fitAGS,Site==site) # look within one site at a time
  R.model<-lmer(relFrProd~std.BCrat+(1|Block),data=site.subset) # fit BC-ratio selection differential at site i
  D.model<-lmer(relFrProd~std.BCdiv+(1|Block),data=site.subset) # fit BC-diversity selection differential at site i
  T.model<-lmer(relFrProd~std.totalAGS+(1|Block),data=site.subset) # fit Total [AGS] selection differential at site i
  rbind(cbind(anova(R.model),data.frame('Beta'=fixef(R.model)[['std.BCrat']],'Int'=fixef(R.model)[['(Intercept)']])),
        cbind(anova(D.model),data.frame('Beta'=fixef(D.model)[['std.BCdiv']],'Int'=fixef(D.model)[['(Intercept)']])),
        cbind(anova(T.model),data.frame('Beta'=fixef(T.model)[['std.totalAGS']],'Int'=fixef(T.model)[['(Intercept)']]))) -> results
  results$Term<-row.names(results)
  results$Site<-site
  if(i==1) { LinPhenoSelDiffs<-results }
  else {LinPhenoSelDiffs<-rbind(LinPhenoSelDiffs,results)}
}
rm(results) # clean up
####### Permutation tests: site-specific selection DIFFERENTIALS #######
## Strategy: permute fecundity within blocks. regress separately on each trait. store coefficients & F values from each permutation
set.seed(777) #for reproducible random sampling
perm.data<-select(fitAGS,Site,Block,relFrProd,std.BCrat,std.totalAGS,std.BCdiv) # copy data
# Execute the loop below then go get a cup of coffee, this will take a few minutes
for (i in 1:1000) {
  group_by(perm.data,Block) %>% mutate(perm.BCrat=sample(std.BCrat),
                                       perm.totalAGS=sample(std.totalAGS),
                                       perm.BCdiv=sample(std.BCdiv)) %>% # permute trait values within each block
    ungroup %>% as.data.frame -> perm.data
  if (i==1) { # initialize data frame to hold results
    perm.results<-data.frame('Site'=character(),'Trait'=character(),'Permutation'=integer(),'Beta'=numeric(),'F'=numeric())
  }
  for(j in 1:4) { # divide data up by Site
    site<-c('Ald','Jam','Mah','Sil')[j]
    site.subset<-filter(perm.data,Site==site)
    # calculate BC-ratio selection differential at site j using permuted data
    perm.R<-(lmer(relFrProd~perm.BCrat+(1|Block),data=site.subset))
    perm.results<-rbind(perm.results,
                        data.frame('Site'=site,'Trait'='BCrat','Permutation'=i,'Beta'=fixef(perm.R)[2],'F'=anova(perm.R)$F.value))
    # calculate Total [AGS] selection differential at site j using permuted data
    perm.T<-(lmer(relFrProd~perm.totalAGS+(1|Block),data=site.subset))
    perm.results<-rbind(perm.results,
                        data.frame('Site'=site,'Trait'='totalAGS','Permutation'=i,'Beta'=fixef(perm.T)[2],'F'=anova(perm.T)$F.value))
    # calculate BC-diversity selection differential at site j using permuted data
    perm.D<-(lmer(relFrProd~perm.BCdiv+(1|Block),data=site.subset))
    perm.results<-rbind(perm.results,
                        data.frame('Site'=site,'Trait'='BCdiv','Permutation'=i,'Beta'=fixef(perm.D)[2],'F'=anova(perm.D)$F.value))
    # clean up:
    rm(perm.R,perm.T,perm.D,site.subset,site)
  }
  print(paste0(' * Permutation number ',i,' complete *'))
}
rm(D.model,R.model,T.model) # clean up
save(perm.results,file=paste0('intermediate_data/sitewise_selection_differentials_permutations_',date(),'.RData'))

## Compare to observed values
LinPhenoSelDiffs$P_perm=7 # initialize permutation test p-value with impossible value
for (j in 1:4){
  site<-c('Ald','Jam','Mah','Sil')[j]
  for (k in 1:3) {
    trait<-c('BCrat','totalAGS','BCdiv')[k]
    obsF<-filter(LinPhenoSelDiffs,Site==site,Term==paste0('std.',trait))$F.value
    greater_perm_F<- filter(perm.results,Site==site,Trait==trait,`F`>=obsF) # only keep permutations with greater F value than observed F value
    greater_perm_F<- dim(greater_perm_F)[1] # get number of permutations
    LinPhenoSelDiffs<-mutate(LinPhenoSelDiffs,P_perm=ifelse(Site==site & Term==paste0('std.',trait),greater_perm_F/1000,P_perm))
  }
}
rm(i,j,k,obsF) # clean up

## Correct for multiple comparisons
group_by(LinPhenoSelDiffs,Term) %>% mutate(P_perm_corrected=p.adjust(P_perm,method='holm')) %>%
  ungroup %>% as.data.frame() -> LinPhenoSelDiffs

filter(LinPhenoSelDiffs,Term=='std.totalAGS') # negative directional selection at all 4 sites
filter(LinPhenoSelDiffs,Term=='std.BCrat') # no directional selection
filter(LinPhenoSelDiffs,Term=='std.BCdiv') # negative directional selection at Ald, Sil; positive selection at Jam
####### Print results of phenotypic selection analysis (differentials) to table #######
select(LinPhenoSelDiffs,Term,Site,Beta,NumDF,DenDF,F.value,P_perm_corrected) %>% 
  mutate(Term=plyr::mapvalues(Term,from=c('std.BCrat','std.BCdiv','std.totalAGS'),to=c(' BC-ratio','BC-diversity',' Total [AGS]'))) %>%
  arrange(Term) %>% 
  mutate(DenDF=round(DenDF),Beta=round(Beta,3),F.value=round(F.value,2)) %>%
  write.table(.,file='tables/seldiffs.txt',sep='\\t',row.names=FALSE,col.names=TRUE)

####### Does selection differ among sites? LOCALLY STANDARDIZED phenotypic selection analysis #######
## All sites: linear selection DIFFERENTIALS - using LOCALLY-STANDARDIZED trait values, fitness
VarLinPhenoSelDiff<-rbind(
  anova(lmer(relFrProd~Site*std.BCrat+(1|Block),data=fitAGS)) %>% mutate(Term=row.names(.),Response='BCrat'),
  anova(lmer(relFrProd~Site*std.BCdiv+(1|Block),data=fitAGS)) %>% mutate(Term=row.names(.),Response='BCdiv'),
  anova(lmer(relFrProd~Site*std.totalAGS+(1|Block),data=fitAGS)) %>% mutate(Term=row.names(.),Response='totalAGS')
) %>% rename(P=`Pr(>F)`)
filter(VarLinPhenoSelDiff,P<0.05,Term!='(Intercept)')
"     Sum Sq   Mean Sq NumDF    DenDF  F.value            P           Term Response
1  79.87266  26.62422     3 1436.853 10.16084 1.268407e-06 Site:std.BCdiv    BCdiv
2 104.93904 104.93904     1 1483.745 41.40364 1.666804e-10   std.totalAGS totalAGS    "
## confirms results from site-specific analyses: 
#### overall selection on totalAGS, plus BCdiv selection varies by site

## Confirm with permutation test:
set.seed(777) #for reproducible random sampling
VarLinPhenoSel.perm.data<-select(fitAGS,Site,Block,relFrProd,std.BCrat,std.totalAGS,std.BCdiv) # copy data
for (i in 1:1000) {
  group_by(VarLinPhenoSel.perm.data,Block) %>% mutate(perm.BCrat=sample(std.BCrat),
                                                      perm.totalAGS=sample(std.totalAGS),
                                                      perm.BCdiv=sample(std.BCdiv)) %>% # permute trait values within each block
    ungroup %>% as.data.frame -> VarLinPhenoSel.perm.data
  if (i==1) { # initialize data frame to hold results
    VarLinPhenoSel.perm.results<-data.frame('Trait'=character(),'Permutation'=integer(),'F.value'=numeric())
  }
  R.model<-anova(lmer(relFrProd~Site*perm.BCrat+(1|Block),data=VarLinPhenoSel.perm.data)) %>% filter(grepl('^Site:perm.',row.names(.)))
  D.model<-anova(lmer(relFrProd~Site*perm.BCdiv+(1|Block),data=VarLinPhenoSel.perm.data)) %>% filter(grepl('^Site:perm.',row.names(.)))
  T.model<-anova(lmer(relFrProd~Site*perm.totalAGS+(1|Block),data=VarLinPhenoSel.perm.data)) %>% filter(grepl('^Site:perm.',row.names(.)))
  VarLinPhenoSel.perm.results<-rbind(VarLinPhenoSel.perm.results,
                                     data.frame('Trait'='BCrat','Permutation'=i,'F.value'=R.model$F.value),
                                     data.frame('Trait'='BCdiv','Permutation'=i,'F.value'=D.model$F.value),
                                     data.frame('Trait'='totalAGS','Permutation'=i,'F.value'=T.model$F.value))
  rm(R.model,D.model,T.model)
  print(paste0(' * Permutation number ',i,' complete *'))
}
rm(i, VarLinPhenoSel.perm.data) # clean up
save(VarLinPhenoSel.perm.results,file=paste0('intermediate_data/variable_selection_differentials_permutations_locally_standardized_TraitValues_',date(),'.RData'))

filter(VarLinPhenoSelDiff,grepl('^Site:',Term)) # view observed F values
# Now ask: How many permutations gave an F value at least as high as the observed F value?
filter(VarLinPhenoSel.perm.results,Trait=='BCrat',F.value>=2.094429) %>% dim # P=0.147
filter(VarLinPhenoSel.perm.results,Trait=='BCdiv',F.value>=10.160843) %>% dim # P<0.001
filter(VarLinPhenoSel.perm.results,Trait=='totalAGS',F.value>=1.482731) %>% dim # P=0.264

select(VarLinPhenoSelDiff,-ends_with('Sq')) %>%
  mutate(DenDF=round(DenDF),F.value=round(F.value,2),P=round(P,3)) %>% 
write.table(.,file='tables/variableselection_locally_standardized.txt',sep='\\t',row.names=FALSE,col.names=TRUE)
####### For Total [AGS] at Silver Creek: does excluding one high-AGS outlier affect conclusion? #######
linsel.totalAGS.Sil.noOutlier <- lmer(relFrProd~std.totalAGS+(1|Block),data=filter(fitAGS,Site=='Sil',std.totalAGS<8))
anova(linsel.totalAGS.Sil.noOutlier)
"Analysis of Variance Table of type III  with  Satterthwaite 
approximation for degrees of freedom
             Sum Sq Mean Sq NumDF  DenDF F.value    Pr(>F)    
std.totalAGS  33.84   33.84     1 355.89  33.354 1.676e-08 ***  "
summary(linsel.totalAGS.Sil.noOutlier)
"              Estimate Std. Error        df t value Pr(>|t|)    
(Intercept)    0.98729    0.07744  20.30000  12.750 3.81e-11 ***
std.totalAGS  -0.38558    0.06676 355.90000  -5.775 1.68e-08 *** "

####### ------ PHENOTYPIC SELECTION ANALYSIS - STABILIZING/DISRUPTIVE ------ #######
####### Phenotypic QUADRATIC selection analysis within each site ####### 
for(i in 1:4) {
  site<-c('Ald','Jam','Mah','Sil')[i]
  site.subset<-filter(fitAGS,Site==site) # look within one site at a time
  R.model<-lmer(relFrProd~std.BCrat+I(std.BCrat^2)+(1|Block),data=site.subset) # fit BC-ratio selection differential at site i
  D.model<-lmer(relFrProd~std.BCdiv+I(std.BCdiv^2)+(1|Block),data=site.subset) # fit BC-diversity selection differential at site i
  T.model<-lmer(relFrProd~std.totalAGS+I(std.totalAGS^2)+(1|Block),data=site.subset) # fit Total [AGS] selection differential at site i
  rbind(cbind(anova(R.model),data.frame('a'=fixef(R.model)[['I(std.BCrat^2)']],'b'=fixef(R.model)[['std.BCrat']],'c'=fixef(R.model)[['(Intercept)']])),
        cbind(anova(D.model),data.frame('a'=fixef(D.model)[['I(std.BCdiv^2)']],'b'=fixef(D.model)[['std.BCdiv']],'c'=fixef(D.model)[['(Intercept)']])),
        cbind(anova(T.model),data.frame('a'=fixef(T.model)[['I(std.totalAGS^2)']],'b'=fixef(T.model)[['std.totalAGS']],'c'=fixef(T.model)[['(Intercept)']]))) -> results
  results$Term<-row.names(results)
  results$Site<-site
  if(i==1) { QuadPhenoSelDiffs<-results }
  else {QuadPhenoSelDiffs<-rbind(QuadPhenoSelDiffs,results)}
}
QuadPhenoSelDiffs<-rename(QuadPhenoSelDiffs,P=`Pr(>F)`) %>% 
  filter(str_detect(Term,'I')) %>% # the statistical test of linear term is meaningless- discard and only keep quadratic term
  mutate(Term=ifelse(Term=='I(std.BCrat^2)','std.BCrat',
                     ifelse(Term=='I(std.totalAGS^2)','std.totalAGS','std.BCdiv'))) # rename terms for simplicity
####### Permutation tests: Quadratic selection at each site #######
## Strategy: permute fecundity within blocks. regress separately on each trait. store coefficients & F values from each permutation
set.seed(777) #for reproducible random sampling
perm.data<-select(fitAGS,Site,Block,relFrProd,std.BCrat,std.totalAGS,std.BCdiv) # copy data
for (i in 1:1000) {
  group_by(perm.data,Block) %>% mutate(perm.BCrat=sample(std.BCrat),
                                       perm.totalAGS=sample(std.totalAGS),
                                       perm.BCdiv=sample(std.BCdiv)) %>% # permute trait values within each block
    ungroup %>% as.data.frame -> perm.data
  if (i==1) { # initialize data frame to hold results
    perm.results<-data.frame('Site'=character(),'Trait'=character(),'Permutation'=integer(),'gamma'=numeric(),'F'=numeric()) # 'a' is the quadratic coefficient
  }
  for(j in 1:4) { # divide data up by Site
    site<-c('Ald','Jam','Mah','Sil')[j]
    site.subset<-filter(perm.data,Site==site)
    # calculate BC-ratio quadratic selection at site j using permuted data
    perm.R<-(lmer(relFrProd~perm.BCrat+I(perm.BCrat^2)+(1|Block),data=site.subset))
    perm.results<-rbind(perm.results,
                        data.frame('Site'=site,'Trait'='BCrat','Permutation'=i,'gamma'=fixef(perm.R)[3],'F'=anova(perm.R)$F.value[2]))
    # calculate Total [AGS] quadratic selection at site j using permuted data
    perm.T<-(lmer(relFrProd~perm.totalAGS+I(perm.totalAGS^2)+(1|Block),data=site.subset))
    perm.results<-rbind(perm.results,
                        data.frame('Site'=site,'Trait'='totalAGS','Permutation'=i,'gamma'=fixef(perm.T)[3],'F'=anova(perm.T)$F.value[2]))
    # calculate BC-diversity quadratic selection at site j using permuted data
    perm.D<-(lmer(relFrProd~perm.BCdiv+I(perm.BCdiv^2)+(1|Block),data=site.subset))
    perm.results<-rbind(perm.results,
                        data.frame('Site'=site,'Trait'='BCdiv','Permutation'=i,'gamma'=fixef(perm.D)[3],'F'=anova(perm.D)$F.value[2]))
    # clean up:
    rm(perm.R,perm.T,perm.D,site.subset,site)
  }
  print(paste0(' * Permutation number ',i,' complete *'))
}
rm(i,j) # clean up
save(perm.results,file=paste0('intermediate_data/sitewise_quadratic_selection_permutations_',date(),'.RData'))

## Compare to observed values
QuadPhenoSelDiffs$P_perm=7 # initialize permutation test p-value with impossible value
for (j in 1:4){
  site<-c('Ald','Jam','Mah','Sil')[j]
  for (k in 1:3) {
    trait<-c('BCrat','totalAGS','BCdiv')[k]
    obsF<-filter(QuadPhenoSelDiffs,Site==site,Term==paste0('std.',trait))$F.value
    greater_perm_F<- filter(perm.results,Site==site,Trait==trait,`F`>=obsF) # only keep permutations with greater F value than observed F value
    greater_perm_F<- dim(greater_perm_F)[1] # get number of permutations
    QuadPhenoSelDiffs<-mutate(QuadPhenoSelDiffs,P_perm=ifelse(Site==site & Term==paste0('std.',trait),greater_perm_F/1000,P_perm))
  }
} rm(j,k,site) # clean up

## Correct for multiple comparisons
group_by(QuadPhenoSelDiffs,Term) %>% mutate(P_perm_corrected=p.adjust(P_perm,method='holm')) %>%
  ungroup %>% as.data.frame() -> QuadPhenoSelDiffs

filter(QuadPhenoSelDiffs,Term=='std.totalAGS') # + quadratic (disruptive) selection at Ald, Sil
filter(QuadPhenoSelDiffs,Term=='std.BCrat') # + quadratic (disruptive) selection at Sil
filter(QuadPhenoSelDiffs,Term=='std.BCdiv') # no significant quadratic selection

####### Write Quadratic Selection results to table #######
write.table(QuadPhenoSelDiffs,'tables/quad_sel_results.txt',sep='\\t',col.names=TRUE,row.names=FALSE)

####### Confirm internal fitness minimum: Total [AGS] at Alder - * #######
### Get fecundity residuals after controlling for block variation
resid.Ald<-resid(lmer(relFrProd~(1|Block),data=filter(fitAGS,Site=='Ald'))) %>%
  as.data.frame(col.names=c('resid')) %>% 
  cbind(.,filter(fitAGS,Site=='Ald') %>% select(std.totalAGS)) 
colnames(resid.Ald)<-c('RelFrProdResid','std.totalAGS')

MOS.totalAGS.Ald <- MOStest(x=resid.Ald$std.totalAGS,y=resid.Ald$RelFrProdResid)
plot(MOS.totalAGS.Ald,which=1) # check fit
plot(MOS.totalAGS.Ald,which=3) # check residuals- looks OK
  
rm(resid.Ald,MOS.totalAGS.Ald) # clean up data objects used in this analysis
####### Confirm internal fitness minimum: Total [AGS] at Silver - *  #######
resid.Sil<-resid(lmer(relFrProd~(1|Block),data=filter(fitAGS,Site=='Sil'))) %>%
  as.data.frame(col.names=c('resid')) %>% 
  cbind(.,filter(fitAGS,Site=='Sil') %>% select(std.totalAGS)) 
colnames(resid.Sil)<-c('RelFrProdResid','std.totalAGS')

MOS.totalAGS.Sil <- MOStest(x=resid.Sil$std.totalAGS,y=resid.Sil$RelFrProdResid)
plot(MOS.totalAGS.Sil,which=1) # check fit
plot(MOS.totalAGS.Sil,which=3) # check residuals-  meh
MOS.totalAGS.Sil
rm(resid.Sil,MOS.totalAGS.Sil) # clean up data objects used in this analysis

### NOTE: from the plot of the fit, looks like one high-[AGS] outlier might be driving result- check what happens when it is removed:
####### Confirm internal fitness minimum: Total [AGS] at Silver - Excluding one high-[AGS] outlier * #######
####### Get residuals after controlling for block variation #######
resid.Sil.noOutlier<-resid(lmer(relFrProd~(1|Block),data=filter(fitAGS,Site=='Sil',std.totalAGS<8))) %>%
  as.data.frame(col.names=c('resid')) %>% 
  cbind(.,filter(fitAGS,Site=='Sil',std.totalAGS<8) %>% select(std.totalAGS)) 
colnames(resid.Sil.noOutlier)<-c('RelFrProdResid','std.totalAGS') # rename fitness residuals column for clarity

MOS.totalAGS.Sil.noOutlier <- MOStest(x=resid.Sil.noOutlier$std.totalAGS,y=resid.Sil.noOutlier$RelFrProdResid)
plot(MOS.totalAGS.Sil.noOutlier,which=1) # check fit- looks a lot better
plot(MOS.totalAGS.Sil.noOutlier,which=3) # check residuals-  look way better
MOS.totalAGS.Sil.noOutlier
rm(resid.Sil.noOutlier,MOS.totalAGS.Sil.noOutlier) # clean up data objects used in this analysis
####### Confirm internal fitness minimum: BC-ratio at Silver - * #######
### Get residuals after controlling for block variation
resid.Sil<-resid(lmer(relFrProd~(1|Block),data=filter(fitAGS,Site=='Sil'))) %>%
  as.data.frame(col.names=c('resid')) %>% 
  cbind(.,filter(fitAGS,Site=='Sil') %>% select(std.BCrat))
colnames(resid.Sil)<-c('RelFrProdResid','std.BCrat')

MOS.BCrat.Sil <- MOStest(x=resid.Sil$std.BCrat,y=resid.Sil$RelFrProdResid,control = list(maxit = 50))
plot(MOS.BCrat.Sil,which=1) # check fit
plot(MOS.BCrat.Sil,which=3) # check residuals - look OK
MOS.BCrat.Sil
rm(resid.Sil,MOS.BCrat.Sil) # clean up data objects used in this analysis

####### ------ TESTING FOR ADAPTIVE PLASTICITY AMONG SITES ------ #######
####### Reshape LS mean trait values for graphing fecundity changes due to plasticity: #######
# Total [AGS]:
genoLSmeans.totalAGS<-filter(LSM.all,Term=='Genotype',Response=='totalAGS') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% # get LS mean trait values for all genotypes (averaged across sites)
  # merge with Genotype*Site LS means (i.e., site-specific LS means for each genotype):
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Ald',Response=='totalAGS') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Ald=lsmean,lower.CL.Ald=lower.CL,upper.CL.Ald=upper.CL),by='Genotype') %>%
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Jam',Response=='totalAGS') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Jam=lsmean,lower.CL.Jam=lower.CL,upper.CL.Jam=upper.CL),by='Genotype') %>%
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Mah',Response=='totalAGS') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Mah=lsmean,lower.CL.Mah=lower.CL,upper.CL.Mah=upper.CL),by='Genotype') %>%
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Sil',Response=='totalAGS') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Sil=lsmean,lower.CL.Sil=lower.CL,upper.CL.Sil=upper.CL),by='Genotype') %>%
  mutate(Trait='totalAGS')

# BC-ratio:
genoLSmeans.BCrat<-filter(LSM.all,Term=='Genotype',Response=='BCrat') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% # get LS mean trait values for all genotypes (averaged across sites)
  # merge with Genotype*Site LS means (i.e., site-specific LS means for each genotype):
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Ald',Response=='BCrat') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Ald=lsmean,lower.CL.Ald=lower.CL,upper.CL.Ald=upper.CL),by='Genotype') %>%
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Jam',Response=='BCrat') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Jam=lsmean,lower.CL.Jam=lower.CL,upper.CL.Jam=upper.CL),by='Genotype') %>%
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Mah',Response=='BCrat') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Mah=lsmean,lower.CL.Mah=lower.CL,upper.CL.Mah=upper.CL),by='Genotype') %>%
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Sil',Response=='BCrat') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Sil=lsmean,lower.CL.Sil=lower.CL,upper.CL.Sil=upper.CL),by='Genotype') %>%
  mutate(Trait='BCrat')

# BC-diversity:
genoLSmeans.BCdiv<-filter(LSM.all,Term=='Genotype',Response=='BCdiv') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% # get LS mean trait values for all genotypes (averaged across sites)
  # merge with Genotype*Site LS means (i.e., site-specific LS means for each genotype):
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Ald',Response=='BCdiv') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Ald=lsmean,lower.CL.Ald=lower.CL,upper.CL.Ald=upper.CL),by='Genotype') %>%
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Jam',Response=='BCdiv') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Jam=lsmean,lower.CL.Jam=lower.CL,upper.CL.Jam=upper.CL),by='Genotype') %>%
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Mah',Response=='BCdiv') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Mah=lsmean,lower.CL.Mah=lower.CL,upper.CL.Mah=upper.CL),by='Genotype') %>%
  merge(.,filter(LSM.all,Term=='Genotype:Site', Site=='Sil',Response=='BCdiv') %>% select(Genotype,lsmean,lower.CL,upper.CL) %>% rename(lsmean.Sil=lsmean,lower.CL.Sil=lower.CL,upper.CL.Sil=upper.CL),by='Genotype') %>%
  mutate(Trait='BCdiv')

# combine all 3 traits:
genoLSmeans <- rbind(genoLSmeans.BCrat,genoLSmeans.totalAGS,genoLSmeans.BCdiv) %>%
  merge(.,ssp,by='Genotype') # add subspecies info
rm(genoLSmeans.BCrat,genoLSmeans.totalAGS,genoLSmeans.BCdiv) # clean up
####### Function for z-transformation of a given trait at a given site #######
ztransform <- function(traitvalue,site,trait) {
  traitvalue.z = (traitvalue - mean(filter(fitAGS,Site==site)[[trait]],na.rm=TRUE))/sd(filter(fitAGS,Site==site)[[trait]],na.rm=TRUE)
  return(traitvalue.z) }

####### Get fitness residuals at each site after controlling for variation due to blocks #######
resid.Ald<-resid(lmer(relFrProd~(1|Block),data=filter(fitAGS,Site=='Ald'))) %>%
  as.data.frame(col.names=c('resid')) %>% cbind(.,filter(fitAGS,Site=='Ald') %>% select(std.BCrat,std.totalAGS,std.BCdiv))
colnames(resid.Ald)<-c('RelFrProdResid','std.BCrat','std.totalAGS','std.BCdiv')

resid.Jam<-resid(lmer(relFrProd~(1|Block),data=filter(fitAGS,Site=='Jam'))) %>%
  as.data.frame(col.names=c('resid')) %>% cbind(.,filter(fitAGS,Site=='Jam') %>% select(std.BCrat,std.totalAGS,std.BCdiv))
colnames(resid.Jam)<-c('RelFrProdResid','std.BCrat','std.totalAGS','std.BCdiv')

resid.Mah<-resid(lmer(relFrProd~(1|Block),data=filter(fitAGS,Site=='Mah'))) %>%
  as.data.frame(col.names=c('resid')) %>% cbind(.,filter(fitAGS,Site=='Mah') %>% select(std.BCrat,std.totalAGS,std.BCdiv))
colnames(resid.Mah)<-c('RelFrProdResid','std.BCrat','std.totalAGS','std.BCdiv')

resid.Sil<-resid(lmer(relFrProd~(1|Block),data=filter(fitAGS,Site=='Sil'))) %>%
  as.data.frame(col.names=c('resid')) %>% cbind(.,filter(fitAGS,Site=='Sil') %>% select(std.BCrat,std.totalAGS,std.BCdiv))
colnames(resid.Sil)<-c('RelFrProdResid','std.BCrat','std.totalAGS','std.BCdiv')

# Make alternative residuals data frame for the case of Total [AGS] at Silver Creek, which has an outlier that skews results:
resid.Sil.noOutlier<-resid(lmer(relFrProd~(1|Block),data=filter(fitAGS,Site=='Sil',std.totalAGS<8))) %>%
  as.data.frame(col.names=c('resid')) %>% cbind(.,filter(fitAGS,Site=='Sil',std.totalAGS<8) %>% select(std.BCrat,std.totalAGS,std.BCdiv))
colnames(resid.Sil.noOutlier)<-c('RelFrProdResid','std.BCrat','std.totalAGS','std.BCdiv')

####### LINEAR ONLY: Function to graph fitness consequences of plasticity among sites for a specific trait at a specific site #######
lin.plast.fit <- function(trait, site,save=FALSE) {
  stdTrait <- paste0('std.',trait) # for calling up columns from data frames
  siteName <- ifelse(site=='Ald','Alder Creek',
                     ifelse(site=='Jam','Jackass Meadow',
                            ifelse(site=='Mah','Mahogany Valley','Silver Creek')))
  print(siteName) # debugging line
  traitName <- ifelse(trait=='BCrat','BC-ratio',
                      ifelse(trait=='BCdiv','BC-diversity','Total [AGS]')) # formal trait name for plot labels
  print(traitName) # debugging line
  ### STEP 1: Calculate selection differential for this trait at this site
  # Express relative fecundity as a function of trait values (either linear or quadratic)
  # Regress relative fecundity RESIDUALS (after controlling for block effects- see above) onto trait values:
  ## Designate name for this selection model:
  modelname<-paste0('linsel.',trait,'.',site,'.resid') # use a different prefix if we are fitting a quadratic model
  ## Load relative-fecundity residuals (after controlling for variation among blocks) for this site- these data frames were created previously (see above):
  if (trait=='totalAGS' & site=='Sil') { residFec <- resid.Sil.noOutlier } # for this one special case, use a data frame of fecundity residuals with one high-[AGS] outlier removed
  else { residFec<-get(paste0('resid.',site)) } # for all other trait-site combinations, use normal residuals dataframe
  ## Specify selection model formula:
  formula<-paste0('residFec$RelFrProdResid~residFec$std.',trait) # specify an additional quadratic term for quadratic selection models
  assign(modelname,lm(as.formula(formula))) # create the selection model object
  
  ## Save model coefficients for easy reference:
  coef.int<-coef(get(modelname))[['(Intercept)']] # save intercept
  coef.lin<-coef(get(modelname))[2] # save linear coefficient
  LinSelFunction <- function(x) coef.int + coef.lin*x
  
  ### STEP 2: Plot to depict AVERAGE effect of plasticity on relative fecundity (average for all genotypes)
  # Use the selection function (fitted above) to predict relative fitness for the BASELINE trait value (LS mean across all sites), which represents a scenario WITHOUT plasticity
  # AND for the site-specific mean trait value, which represents a scenario WITH plasticity 
  # Note this is the average for all genotypes.
  # Z-transform these LS mean trait values using data specific to this site:
  global.LSmean.traitValue <- ztransform(mean(filter(LSM.all,Term=='Site',Response==trait)$lsmean),site,trait) # use mean of Site LS means as baseline trait value
  siteSpecific.LSmean.traitValue <- ztransform(filter(LSM.all,Term=='Site',Response==trait,Site==site)$lsmean,site,trait) # compare to site-specific LS mean trait value
  # Plot selection function with baseline and site-specific trait values:
  selFunctionPlot <- stat_function(fun=LinSelFunction,color='red') 
  # Make ggplot object to summarize effect of plasticity on fitness for this trait, this site:
  selplot<-ggplot(residFec,aes(x=residFec[[stdTrait]],y=RelFrProdResid))+
    geom_point(alpha=0.25) +
    selFunctionPlot + # will fit either linear or quadratic regression depending on specified selection 'type' - see above
    geom_vline(xintercept=siteSpecific.LSmean.traitValue,color='#000000')+ # plot LS mean trait value for this site (scenario with plasticity)
    geom_vline(xintercept=global.LSmean.traitValue,color='grey60',linetype='dashed')+ # plot LS mean trait value averaged across sites (scenario without plasticity)
    theme_classic()+
    labs(x=paste0(traitName,' (z-transformed)'), y=paste0('Relative fecundity at ',siteName))
  
  ### STEP 3: Genotype-specific plasticity and selection
  # now summarize for each genotype: what is the estimated change in relative fitness for each genotype due to its plasticity? ###
  # Z-transform LS means for each genotype: overall genotype LS mean (i.e., 'non plastic genotype' scenario) and genotype*site-specific LS mean
  colnames(genoLSmeans)<-str_replace(colnames(genoLSmeans),paste0('.',site),'.site') # rename the needed columns (for this field site) to have the prefix ".site", so they're easier to call up
  
  genoPlastFit <- filter(genoLSmeans,Trait==trait) %>% # Start with data frame of all genotype LS means (both experiment-wide, and specific to each site- from Geno x Site interaction term) 
    select(Genotype,Ssp,lsmean,lower.CL,upper.CL,ends_with('.site')) # remove unneeded columns
  
  genoPlastFit <- mutate(genoPlastFit, Site=site,
                         z.lsmean.global=ztransform(lsmean,site,trait), # apply z-transformations to the global LS mean (baseline) using the observed values at this site
                         z.lower.CL.global=ztransform(lower.CL,site,trait), # lower end of 95% CI for global LS mean trait value
                         z.upper.CL.global=ztransform(upper.CL,site,trait), # upper end of CI
                         z.lsmean.site=ztransform(lsmean.site,site,trait), # Do the same for site-specific LS mean trait values
                         z.lower.CL.site=ztransform(lower.CL.site,site,trait), # lower end of 95% CI for site-specific LS mean trait value
                         z.upper.CL.site=ztransform(upper.CL.site,site,trait)) # upper end of CI
  genoPlastFit <- select(genoPlastFit,Genotype,Site,Ssp,starts_with('z')) %>%  # only keep z-transformed values
    mutate(Trait=trait,
           # plug in z-transformed experiment-wide LS mean to estimate relative fitness of a hypothetical "non plastic" version of this genotype:
           w.noPlast = LinSelFunction(z.lsmean.global), 
           # plug in z-transformed site-specific LS mean to estimate the genotype's true relative fitness in this site:
           w.withPlast = LinSelFunction(z.lsmean.site), 
           # the trait estimate has a 95% confidence interval; incorporate this uncertainty:
           # Depending on whether selection is negative or positive, the upper fitness estimate may come from either the upper or lower trait value estimate (and same for the lower bound of the fitness CI):
           w.withPlast.lower = LinSelFunction(z.lower.CL.site),
           w.withPlast.upper = LinSelFunction(z.upper.CL.site), 
           delta_w = w.withPlast-w.noPlast, # calculate estimated change in relative fitness * due to plasticity * : subtract the genotype's expected fitness if there were no plasticity from its true site-specific fitness to quantify the change in RELATIVE FITNESS due to plasticity
           delta_w_lower = w.withPlast.lower - w.noPlast, # 95% CI lower limit for the estimated change in relative fitness * due to plasticity * 
           delta_w_upper = w.withPlast.upper - w.noPlast) # 95% CI upper limit for the estimated change in relative fitness * due to plasticity * 
  # Designate each delta_w observation as significantly nonzero or not
  genoPlastFit<-mutate(genoPlastFit,PlastType=case_when(
    genoPlastFit$delta_w_lower<0 & genoPlastFit$delta_w_upper>0 ~ 'neither', # if 95% CI includes 0
    genoPlastFit$delta_w_lower>0 & genoPlastFit$delta_w_upper<0 ~ 'neither', # if 95% CI includes 0
    genoPlastFit$delta_w_lower>0 & genoPlastFit$delta_w_upper>0 ~ 'adaptive', # if entire 95% CI >0 --> adaptive plasticity
    genoPlastFit$delta_w_lower<0 & genoPlastFit$delta_w_upper<0 ~ 'nonadaptive')) # if entire 95% CI <0 --> nonadaptive plasticity
  genoPlastFit<-mutate(genoPlastFit,Ssp=ifelse(Ssp=='E','East\\nsubspecies', 'West\\nsubspecies')) # rename subspecies for labeling
  # plot it:
  genoPlastFitPlot<- ggplot(genoPlastFit,aes(x=Genotype,y=delta_w,color=PlastType)) +
    facet_wrap(~Ssp,scales='free_x')+
    geom_hline(yintercept=0,color='black')+
    geom_point(aes(size=PlastType))+
    geom_errorbar(aes(ymin=delta_w_lower,ymax=delta_w_upper),width=0)+
    scale_size_manual(values=c('adaptive'=2.5,'nonadaptive'=2.5,'neither'=1.5))+
    scale_color_manual(values=c('adaptive' = 'orchid','nonadaptive'='blue','neither'='grey70'))+
    theme_classic()+
    ylab('Change in local relative fecundity\\ndue to plasticity among sites (95% CI)')+
    theme(axis.text.x = element_text(angle=65,vjust=0.7,hjust=0.6))+
    theme(axis.title.x = element_blank(), legend.position='none')+
    theme(strip.background=element_blank(),strip.text=element_blank())
  # specify what to name the file:
  filename <- paste0('figures/linear_plastSelPlot_',trait,'_',site,'.pdf')
  if (save==TRUE) { ggsave(filename=filename,
                           plot_grid(selplot,genoPlastFitPlot),
                           width=9,height=3)}
  else {print(selplot); print(genoPlastFitPlot)}
  return(genoPlastFit) # return the dataframe 
}
####### QUADRATIC ONLY: Function to graph fitness consequences of plasticity among sites for a specific trait at a specific site #######
quad.plast.fit <- function(trait, site,save=FALSE) {
  stdTrait <- paste0('std.',trait) # for calling up columns from data frames
  siteName <- ifelse(site=='Ald','Alder Creek',
                     ifelse(site=='Jam','Jackass Meadow',
                            ifelse(site=='Mah','Mahogany Valley','Silver Creek')))
  print(siteName) # debugging line
  traitName <- ifelse(trait=='BCrat','BC-ratio',
                      ifelse(trait=='BCdiv','BC-diversity','Total [AGS]')) # formal trait name for plot labels
  print(traitName) # debugging line
  ### STEP 1: Calculate selection differential for this trait at this site
  # Express relative fecundity as a function of trait values (either linear or quadratic)
  # Regress relative fecundity RESIDUALS (after controlling for block effects- see above) onto trait values:
  ## Designate name for this selection model:
  modelname<-paste0('quadsel.',trait,'.',site,'.resid') # use a different prefix if we are fitting a quadratic model
  ## Load relative-fecundity residuals (after controlling for variation among blocks) for this site- these data frames were created previously (see above):
  if (trait=='totalAGS' & site=='Sil') { residFec <- resid.Sil.noOutlier } # for this one special case, use a data frame of fecundity residuals with one high-[AGS] outlier removed
  else { residFec<-get(paste0('resid.',site)) } # for all other trait-site combinations, use normal residuals dataframe
  ## Specify selection model formula:
  formula<-paste0('residFec$RelFrProdResid~residFec$std.',trait,' + I(residFec$std.',trait,'^2)') # specify an additional quadratic term for quadratic selection models
  assign(modelname,lm(as.formula(formula))) # create the selection model object
  
  ## Save model coefficients for easy reference:
  coef.int<-coef(get(modelname))[['(Intercept)']] # save intercept
  coef.lin<-coef(get(modelname))[2] # save linear coefficient
  coef.quad <- 2*coef(get(modelname))[3] # DOUBLE the quadratic coefficient (Stinchcombe et al. 2008 - Evolution) and save it
  z.vertex <- -coef.lin/(2*coef.quad) # location of min/max- from quadratic equation (-b/(2a))
  QuadSelFunction <- function(x) coef.int + coef.lin*x + coef.quad*(x^2)
  
  ### STEP 2: Plot to depict AVERAGE effect of plasticity on relative fecundity (average for all genotypes)
  # Use the selection function (fitted above) to predict relative fitness for the BASELINE trait value (LS mean across all sites), which represents a scenario WITHOUT plasticity
  # AND for the site-specific mean trait value, which represents a scenario WITH plasticity 
  # Note this is the average for all genotypes.
  # Z-transform these LS mean trait values using data specific to this site:
  global.LSmean.traitValue <- ztransform(mean(filter(LSM.all,Term=='Site',Response==trait)$lsmean),site,trait) # use mean of Site LS means as baseline trait value
  siteSpecific.LSmean.traitValue <- ztransform(filter(LSM.all,Term=='Site',Response==trait,Site==site)$lsmean,site,trait) # compare to site-specific LS mean trait value
  # Plot selection function with baseline and site-specific trait values:
  selFunctionPlot <- stat_function(fun=QuadSelFunction,color='red') # for quadratic selection- have to do it this way to incorporate doubled gamma
  # Make ggplot object to summarize effect of plasticity on fitness for this trait, this site:
  selplot<-ggplot(residFec,aes(x=residFec[[stdTrait]],y=RelFrProdResid))+
    geom_point(alpha=0.25) +
    selFunctionPlot + # will fit either linear or quadratic regression depending on specified selection 'type' - see above
    geom_vline(xintercept=siteSpecific.LSmean.traitValue,color='#000000')+ # plot LS mean trait value for this site (scenario with plasticity)
    geom_vline(xintercept=global.LSmean.traitValue,color='grey60',linetype='dashed')+ # plot LS mean trait value averaged across sites (scenario without plasticity)
    theme_classic()+
    labs(x=paste0(traitName,' (z-transformed)'), y=paste0('Relative fecundity at ',siteName))
  
  ### STEP 3: Genotype-specific plasticity and selection
  # now summarize for each genotype: what is the estimated change in relative fitness for each genotype due to its plasticity? ###
  # Z-transform LS means for each genotype: overall genotype LS mean (i.e., 'non plastic genotype' scenario) and genotype*site-specific LS mean
  colnames(genoLSmeans)<-str_replace(colnames(genoLSmeans),paste0('.',site),'.site') # rename the needed columns (for this field site) to have the prefix ".site", so they're easier to call up
  
  genoPlastFit <- filter(genoLSmeans,Trait==trait) %>% # Start with data frame of all genotype LS means (both experiment-wide, and specific to each site- from Geno x Site interaction term) 
    select(Genotype,Ssp,lsmean,lower.CL,upper.CL,ends_with('.site')) # remove unneeded columns
  
  genoPlastFit <- mutate(genoPlastFit, Site=site,
                         z.lsmean.global=ztransform(lsmean,site,trait), # apply z-transformations to the global LS mean (baseline) using the observed values at this site
                         z.lower.CL.global=ztransform(lower.CL,site,trait), # lower end of 95% CI for global LS mean trait value
                         z.upper.CL.global=ztransform(upper.CL,site,trait), # upper end of CI
                         z.lsmean.site=ztransform(lsmean.site,site,trait), # Do the same for site-specific LS mean trait values
                         z.lower.CL.site=ztransform(lower.CL.site,site,trait), # lower end of 95% CI for site-specific LS mean trait value
                         z.upper.CL.site=ztransform(upper.CL.site,site,trait), # upper end of CI
                         vertex.in.CI = ifelse(z.vertex %in% c(z.lower.CL.site,z.upper.CL.site),'yes','no')) # assess whether fitness max/min falls in trait value 95% CI
  genoPlastFit <- select(genoPlastFit,Genotype,Site,Ssp,starts_with('z'),vertex.in.CI) %>%  # only keep z-transformed values
    mutate(Trait=trait,
           # plug in z-transformed experiment-wide LS mean to estimate relative fitness of a hypothetical "non plastic" version of this genotype:
           w.noPlast = QuadSelFunction(z.lsmean.global), 
           # plug in z-transformed site-specific LS mean to estimate the genotype's true relative fitness in this site:
           w.withPlast = QuadSelFunction(z.lsmean.site),
           # the trait estimate has a 95% confidence interval; incorporate this uncertainty by evaluating fitness at CI boundaries:
           w.withPlast.lowerZ = QuadSelFunction(z.lower.CL.site), 
           w.withPlast.upperZ = QuadSelFunction(z.upper.CL.site))
  # We now have relative fecundity estimates for the site-specific LS mean trait value and also the boundaries of the 95% CI for trait value at this site.
  # Next step: Subtract the global mean trait value to translate these values into DELTA fitness (fitness change due to plasticity):
  genoPlastFit <- mutate(genoPlastFit, delta_w = w.withPlast - w.noPlast,
                         delta_w_upper = w.withPlast.upperZ - w.noPlast,
                         delta_w_lower = w.withPlast.lowerZ - w.noPlast)
  # Now: determine error bars (needed because for some genotypes, delta_w_lower > delta_w_upper -- the 'upper/lower' refer to limits of trait value CI, NOT fitness CI)
  # Also, for quadratic selection: if the site-specific mean trait value falls near the inflection point, both extremes of the trait value 95% CI may change fitness in the same direction. 
  # Account for this ^ by defining the 95% CI of fecundity as the minimum and maximum values of the fitness function evaluated over the 95% CI of trait values
  genoPlastFit<-mutate(genoPlastFit,
                       errorBarUp=ifelse(delta_w_upper>delta_w_lower,delta_w_upper,delta_w_lower),
                       errorBarDown=ifelse(delta_w_upper>delta_w_lower,delta_w_lower,delta_w_upper))
  # Finally, account for special cases in which the whole fitness function's maximum or minimum falls within the trait value 95% CI:
  if (coef.quad > 0) { # for disruptive selection, update fitness minimum
    genoPlastFit <- mutate(genoPlastFit, errorBarDown = ifelse(delta_w < delta_w_upper & delta_w < delta_w_lower,
                                                               QuadSelFunction(z.vertex) - w.noPlast, errorBarDown))
  } 
  else if (coef.quad < 0) { # for stabilizing selection, update fitness maximum
    genoPlastFit <- mutate(genoPlastFit, errorBarUp = ifelse(delta_w > delta_w_upper & delta_w > delta_w_lower,
                                                             QuadSelFunction(z.vertex) - w.noPlast, errorBarUp))
  }
  # Finally, characterize each genotype's plastic response at this site (relative to "non-plastic" hypothetical genotype) as "adaptive", "nonadaptive", or "neither":
  genoPlastFit<-mutate(genoPlastFit,PlastType=case_when(
    genoPlastFit$errorBarDown<0 & genoPlastFit$errorBarUp>0 ~ 'neither', # if delta_fecundity 95% CI includes 0
    genoPlastFit$errorBarDown>0 & genoPlastFit$errorBarUp>0 ~ 'adaptive', # if entire delta_fecundity 95% CI >0 --> adaptive plasticity
    genoPlastFit$errorBarDown<0 & genoPlastFit$errorBarUp<0 ~ 'nonadaptive')) # if entire delta_fecundity 95% CI <0 --> nonadaptive plasticity
  # Now, plot it:
  genoPlastFitPlot<- ggplot(genoPlastFit,aes(x=Genotype,y=delta_w,color=PlastType)) +
    facet_wrap(~Ssp,scales='free_x')+
    geom_hline(yintercept=0,color='black')+
    geom_point(aes(size=PlastType))+
    geom_errorbar(aes(ymin=errorBarDown,ymax=errorBarUp),width=0)+
    scale_size_manual(values=c('adaptive'=2.5,'nonadaptive'=2.5,'neither'=1.5))+
    scale_color_manual(values=c('adaptive' = 'orchid','nonadaptive'='blue','neither'='grey70'))+
    theme_classic()+
    ylab('Change in local relative fecundity\\ndue to plasticity among sites (95% CI)')+
    theme(axis.text.x = element_text(angle=65,vjust=0.7,hjust=0.6))+
    theme(axis.title.x = element_blank(), legend.position='none')+
    theme(strip.background=element_blank(),strip.text=element_blank())
  # specify what to name the file:
  filename <- paste0('figures/quadratic_plastSelPlot_',trait,'_',site,'.pdf')
  if (save==TRUE) { ggsave(filename=filename,
                           plot_grid(selplot,genoPlastFitPlot),
                           width=9,height=3)}
  else {print(selplot); print(genoPlastFitPlot)}
  return(genoPlastFit) # return the dataframe 
}

####### Analyze plasticity & fitness for all cases of selection, significant or not  #######
plastfit.linear<-
  rbind(lin.plast.fit('totalAGS','Ald',save=TRUE), # 
        lin.plast.fit('totalAGS','Jam',save=TRUE), # 
        lin.plast.fit('totalAGS','Mah',save=TRUE), # 
        lin.plast.fit('totalAGS','Sil',save=TRUE), # 
        lin.plast.fit('BCdiv','Ald',save=TRUE), # 
        lin.plast.fit('BCdiv','Jam',save=TRUE), # 
        lin.plast.fit('BCdiv','Mah',save=TRUE), # 
        lin.plast.fit('BCdiv','Sil',save=TRUE), # 
        lin.plast.fit('BCrat','Ald',save=TRUE), # 
        lin.plast.fit('BCrat','Jam',save=TRUE), #
        lin.plast.fit('BCrat','Mah',save=TRUE), # 
        lin.plast.fit('BCrat','Sil',save=TRUE)) #
# Quadratic:
plastfit.quadratic<-
  rbind(quad.plast.fit('totalAGS','Ald',save=TRUE), # 
        quad.plast.fit('totalAGS','Jam',save=TRUE), #
        quad.plast.fit('totalAGS','Mah',save=TRUE), # 
        quad.plast.fit('totalAGS','Sil',save=TRUE), # 
        quad.plast.fit('BCdiv','Ald',save=TRUE), # 
        quad.plast.fit('BCdiv','Jam',save=TRUE), # 
        quad.plast.fit('BCdiv','Mah',save=TRUE), # 
        quad.plast.fit('BCdiv','Sil',save=TRUE), #
        quad.plast.fit('BCrat','Ald',save=TRUE), # 
        quad.plast.fit('BCrat','Jam',save=TRUE), #
        quad.plast.fit('BCrat','Mah',save=TRUE), # 
        quad.plast.fit('BCrat','Sil',save=TRUE)) #
dev.off()

####### *** Results: linear selection #######
table(plastfit.linear$Ssp, plastfit.linear$PlastType)
"                    adaptive neither nonadaptive
East\\nsubspecies       42      87          15
West\\nsubspecies       15     123           6 "

## EAST
# Adaptive = 42/144 = 29.2%
# Nonadaptive = 15/144 = 10.4%
# Neutral = 87/144 = 60.4%

## WEST
# Adaptive = 15/144 = 10.4%
# Nonadaptive = 6/144 = 4.2%
# Neutral = 123/144 = 85.4%

## TOTAL
# Adaptive = 57/288 = 19.8%
# Nonadaptive = 21/288 = 7.3%
# Neutral = 210/288 = 72.9%

### Hypothesis test: overall, are adaptive and nonadaptive plasticity equally likely?
# Probability of observing 57 of 78 cases of non-neutral plasticity in adaptive direction due to random chance (prob=0.5):
binom.test(57, 78, p = 0.5,alternative = 'two.sided',conf.level = 0.95)


### Chi-squared test: Do subspecies differ in relative frequency of adaptive vs. nonadaptive plasticity?
PlastFitCounts.linear<-as.table(rbind(c(42,15),c(15,6))) # enter counts of adaptive vs. nonadaptive plasticity for each subspecies
dimnames(PlastFitCounts.linear)<-list(Ssp=c('E','W'),Type=c('Adaptive','Nonadaptive'))
chisq.test(PlastFitCounts.linear,simulate.p.value=TRUE,B=10000)
"X-squared = 0.039686, df = NA, p-value = 1" 

### Chi-squared test:
### Do subspecies differ in how often plasticity is neutral vs. non-neutral? 
PlastFitCounts.linear2<-as.table(rbind(c(87),c(123)))
dimnames(PlastFitCounts.linear2)<-list(Ssp=c('E','W'),Type=c('Yes'))
chisq.test(PlastFitCounts.linear2,simulate.p.value=TRUE,B=10000)

### Barplot: summarize by genotype:
pdf(file='figures/Figure4b.pdf',width=9,height=3)
plastfit.linear$PlastType<-ordered(plastfit.linear$PlastType, levels = c('adaptive', 'nonadaptive', 'neither')) # order the categories
ggplot(plastfit.linear,aes(x=Genotype,fill=PlastType))+
  facet_wrap(~Ssp,scales='free_x')+
  geom_bar(stat='count')+
  scale_fill_manual(values=c('orchid','blue','grey70'))+
  scale_y_continuous(breaks=c(0,5,10,15,20))+
  theme_classic()+
  theme(axis.text.x=element_text(size=15,angle=75,vjust=0.6))+
  theme(axis.text.y=element_text(size=15))+
  theme(axis.title.y=element_blank(),axis.title.x=element_text(size=20))+
  theme(legend.title=element_blank())+
  theme(strip.background=element_rect(fill='grey50',color='grey50'))
dev.off()

####### *** Results: quadratic selection #######
table(plastfit.quadratic$Ssp, plastfit.quadratic$PlastType)
"    adaptive neither nonadaptive
  E       47      88           9
  W       16     123           5 "

## EAST
# Adaptive = 47/144 = 32.6%
# Nonadaptive = 9/144 = 6.3%
# Neutral = 88/144 = 61.1%

## WEST
# Adaptive = 16/144 = 11.1%
# Nonadaptive = 5/144 = 3.5%
# Neutral = 123/144 = 85.4%

## TOTAL
# Adaptive = 63/288 = 21.9%
# Nonadaptive = 14/288 = 4.9%
# Neutral = 211/288 = 73.3%

### Hypothesis test: overall, are adaptive and nonadaptive plasticity equally likely?
# Probability of observing 62 of 77 cases of non-neutral plasticity in adaptive direction due to random chance (prob=0.5):
binom.test(63, 77, p = 0.5,alternative = 'two.sided',conf.level = 0.95)

### Chi-squared test: Do subspecies differ in relative frequency of adaptive vs. nonadaptive plasticity?
PlastFitCounts.quadratic<-as.table(rbind(c(47,9),c(16,5))) # enter counts of adaptive vs. nonadaptive plasticity for each subspecies
dimnames(PlastFitCounts.quadratic)<-list(Ssp=c('E','W'),Type=c('Adaptive','Nonadaptive'))
chisq.test(PlastFitCounts.quadratic,simulate.p.value=TRUE,B=10000)

### Chi-squared test: Do subspecies differ in how often plasticity is neutral vs. non-neutral? 
PlastFitCounts.quadratic2<-as.table(rbind(c(88),c(123)))
dimnames(PlastFitCounts.quadratic2)<-list(Ssp=c('E','W'),Type=c('Yes'))
chisq.test(PlastFitCounts.quadratic2,simulate.p.value=TRUE,B=10000)

### Barplot: summarize by genotype:
pdf(file='figures/Supp_Fig8b.pdf',width=9,height=3)
plastfit.quadratic$PlastType<-ordered(plastfit.quadratic$PlastType, levels = c('adaptive', 'nonadaptive', 'neither')) # order the categories
ggplot(plastfit.quadratic,aes(x=Genotype,fill=PlastType))+
  facet_wrap(~Ssp,scales='free_x')+
  geom_bar(stat='count')+
  scale_fill_manual(values=c('orchid','blue','grey70'))+
  scale_y_continuous(breaks=c(0,5,10,15,20))+
  theme_classic()+
  theme(axis.text.x=element_text(size=15,angle=75,vjust=0.6))+
  theme(axis.text.y=element_text(size=15))+
  theme(axis.title.y=element_blank(),axis.title.x=element_text(size=20))+
  theme(legend.title=element_blank())+
  theme(strip.background=element_rect(fill='grey50',color='grey50'))
dev.off()

####### Selection GRADIENTS- simultaneous selection on all traits #######
####### Phenotypic selection analysis-- GRADIENTS -- Estimate linear selection gradients at each site ####### 
# fit same model with all 3 traits at each site and save ANOVA results:
LinPhenoSelGradStats<-rbind(anova(lmer(relFrProd~std.BCrat+std.totalAGS+std.BCdiv+(1|Block),data=subset(fitAGS,Site=='Ald')))[1:6] %>% mutate(Trait=str_replace(row.names(.),'std.',''),Site='Ald'),
                            anova(lmer(relFrProd~std.BCrat+std.totalAGS+std.BCdiv+(1|Block),data=subset(fitAGS,Site=='Jam')))[1:6] %>% mutate(Trait=str_replace(row.names(.),'std.',''),Site='Jam'),
                            anova(lmer(relFrProd~std.BCrat+std.totalAGS+std.BCdiv+(1|Block),data=subset(fitAGS,Site=='Mah')))[1:6] %>% mutate(Trait=str_replace(row.names(.),'std.',''),Site='Mah'),
                            anova(lmer(relFrProd~std.BCrat+std.totalAGS+std.BCdiv+(1|Block),data=subset(fitAGS,Site=='Sil')))[1:6] %>% mutate(Trait=str_replace(row.names(.),'std.',''),Site='Sil'))
# correct P-values
group_by(LinPhenoSelGradStats,Trait) %>% mutate(P_corrected=p.adjust(`Pr(>F)`,method='holm')) %>%
  ungroup %>% as.data.frame() %>% filter(P_corrected<0.05)
# Results: congruent with selection differentials, except no hit for Total [AGS] at Jam, extra hit for BC-ratio at Mah
# Check directions of selection gradients:
lmer(relFrProd~std.BCrat+std.totalAGS+std.BCdiv+(1|Block),data=subset(fitAGS,Site=='Ald')) %>% summary()
# Alder Creek: negative selection on Total [AGS], negative selection on BC-diversity
lmer(relFrProd~std.BCrat+std.totalAGS+std.BCdiv+(1|Block),data=subset(fitAGS,Site=='Jam')) %>% summary()
# Jackass Meadow: positive selection on BC-diversity
lmer(relFrProd~std.BCrat+std.totalAGS+std.BCdiv+(1|Block),data=subset(fitAGS,Site=='Mah')) %>% summary()
# Mahogany Valley: negative selection on Total [AGS], negative selection on BC-ratio
lmer(relFrProd~std.BCrat+std.totalAGS+std.BCdiv+(1|Block),data=subset(fitAGS,Site=='Sil')) %>% summary()
# Silver Creek: negative selection on Total [AGS], negative selection on BC-diversity
# Conclusion: no disagreement about direction of selection on any trait or any site
####### Permutation tests: site-specific selection GRADIENTS #######
## Strategy: permute trait values within blocks. regress relative fecundity onto trait values (permute only 1 trait at a time). store coefficients & F values from each permutation
set.seed(777) #for reproducible random sampling
perm.data.grad<-select(fitAGS,Site,Block,relFrProd,std.BCrat,std.totalAGS,std.BCdiv) # copy data
for (i in 1:1000) {
  group_by(perm.data.grad,Block) %>% mutate(perm.BCrat=sample(std.BCrat),
                                            perm.totalAGS=sample(std.totalAGS),
                                            perm.BCdiv=sample(std.BCdiv)) %>% # permute trait values within each block
    ungroup %>% as.data.frame -> perm.data.grad
  if (i==1) { # initialize data frame to hold results
    perm.results.grad<-data.frame('Site'=character(),'Trait'=character(),'Permutation'=integer(),'Beta'=numeric(),'F'=numeric())
  }
  for(j in 1:4) { # divide data up by Site
    site<-c('Ald','Jam','Mah','Sil')[j]
    site.subset<-filter(perm.data.grad,Site==site)
    # calculate BC-ratio selection differential at site j using permuted values for BCrat only
    perm.R<-(lmer(relFrProd~perm.BCrat+std.totalAGS+std.BCdiv+(1|Block),data=site.subset))
    perm.results.grad<-rbind(perm.results.grad,
                             data.frame('Site'=site,'Trait'='BCrat','Permutation'=i,'Beta'=fixef(perm.R)['perm.BCrat'],'F'=anova(perm.R)['perm.BCrat',]$F.value))
    # calculate Total [AGS] selection differential at site j using permuted values for total GS only
    perm.T<-(lmer(relFrProd~std.BCrat+perm.totalAGS+std.BCdiv+(1|Block),data=site.subset))
    perm.results.grad<-rbind(perm.results.grad,
                             data.frame('Site'=site,'Trait'='totalAGS','Permutation'=i,'Beta'=fixef(perm.T)['perm.totalAGS'],'F'=anova(perm.T)['perm.totalAGS',]$F.value))
    # calculate BC-diversity selection differential at site j using permuted values for BCdiv only
    perm.D<-(lmer(relFrProd~std.BCrat+std.totalAGS+perm.BCdiv+(1|Block),data=site.subset))
    perm.results.grad<-rbind(perm.results.grad,
                             data.frame('Site'=site,'Trait'='BCdiv','Permutation'=i,'Beta'=fixef(perm.D)['perm.BCdiv'],'F'=anova(perm.D)['perm.BCdiv',]$F.value))
    # clean up:
    rm(perm.R,perm.T,perm.D,site.subset,site)
  }
  print(paste0(' * Permutation number ',i,' complete *'))
}
# Save the permutations:
save(perm.results.grad,file=paste0('intermediate_data/sitewise_selection_gradients_permutations_',date(),'.RData'))

## Compare to observed values
LinPhenoSelGradStats$P_perm=7 # initialize permutation test p-value with impossible value
for (j in 1:4){
  site<-c('Ald','Jam','Mah','Sil')[j]
  for (k in 1:3) {
    trait <- c('BCrat','totalAGS','BCdiv')[k]
    obsF <- filter(LinPhenoSelGradStats,Site==site,Trait==trait)$F.value # get the actual observed F value from non-permuted model
    greater_perm_F <- filter(perm.results.grad,Site==site,Trait==trait,`F`>=obsF) # only keep permutations with F value equal to/greater than observed F value
    greater_perm_F <- dim(greater_perm_F)[1] # get number of permutations with F value equal to/higher than observed F value
    LinPhenoSelGradStats<-mutate(LinPhenoSelGradStats,P_perm=ifelse(Site==site & Trait==trait,greater_perm_F/1000,P_perm))
  }
}
## Correct for multiple comparisons
group_by(LinPhenoSelGradStats,Trait) %>% mutate(P_perm_corrected=p.adjust(P_perm,method='holm')) %>%
  ungroup %>% as.data.frame() -> LinPhenoSelGradStats

filter(LinPhenoSelGradStats,Trait=='totalAGS') # selection at Ald, Mah, Sil
filter(LinPhenoSelGradStats,Trait=='BCrat') # selection at Mah
filter(LinPhenoSelGradStats,Trait=='BCdiv') # selection at Ald, Sil (P=0.08 at Jam)
####### Print results of phenotypic selection gradients to table #######
select(LinPhenoSelGradStats,Trait,Site,NumDF,DenDF,F.value,P_perm_corrected) %>% 
  mutate(Trait=plyr::mapvalues(Trait,from=c('BCrat','BCdiv','totalAGS'),to=c(' BC-ratio','BC-diversity',' Total [AGS]'))) %>%
  arrange(Site) %>% 
  mutate(DenDF=round(DenDF),F.value=round(F.value,2)) %>%
  write.table(.,file='tables/selectionGradients.txt',sep='\\t',row.names=FALSE,col.names=TRUE)

####### ------ REACTION NORM ANALYSIS: using pooled trait means ------ #######
####### Get EIs (Pooled Mean Trait Values): Calculate Block grand means for each trait #######
## Calculate EIs as pooled mean trait values:
group_by(AGS.withBHS,Block) %>% summarize(blkBCrat=mean(BCrat,na.rm=TRUE),
                                          blktotalAGS=mean(totalAGS,na.rm=TRUE),
                                          blkBCdiv=mean(BCdiv,na.rm=TRUE)) %>%
  ungroup %>% as.data.frame %>%
  merge(.,blks,by='Block') -> blkAGSmeans
#######  Calculate line mean trait values for each genotype in each block #######
group_by(AGS.withBHS,Genotype,Block) %>% summarize(BCrat=mean(BCrat,na.rm=TRUE),
                                                   totalAGS=mean(totalAGS,na.rm=TRUE),
                                                   BCdiv=mean(BCdiv,na.rm=TRUE)) %>%
  ungroup %>% as.data.frame %>% 
  merge(.,blkAGSmeans,by='Block',suffixes=c('','.y')) %>% 
  mutate(Genotype=factor(Genotype)) %>% 
  merge(., ssp,by='Genotype') %>% 
  rename(EI.BCrat=blkBCrat,EI.totalAGS=blktotalAGS,EI.BCdiv=blkBCdiv) -> genoblockMeans

save(genoblockMeans,file='intermediate_data/genoblock_Means.RData')
#######  Calculate linear reaction norm for each genotype, each trait - separately at each site #######
for(site in c('Ald','Jam','Mah','Sil')) { # one site at a time
  for(geno in levels(genoblockMeans$Genotype)) { # one genotype at a time
    if (!(site=='Jam' & geno=='BHS')) {  # skip the combination of genotype BHS at site Jam (no surviving individuals)
    print(paste0('Calculating reaction norms for genotype ',geno,' at site ',geno)) # testing line
    geno.subset<-filter(genoblockMeans,Genotype==geno,Site==site) # subset the data accordingly
    R.model<-lm(BCrat~EI.BCrat,data=geno.subset) # fit BC-ratio reaction norm for this genotype at this site
    D.model<-lm(BCdiv~EI.BCdiv,data=geno.subset) # fit BC-diversity reaction norm for this genotype at this site
    T.model<-lm(totalAGS~EI.totalAGS,data=geno.subset) # fit Total [AGS] reaction norm for this genotype at this site
    rbind(cbind(as.data.frame(t(summary(R.model)$coefficients[2,])),data.frame('RN.slope'=R.model$coefficients[['EI.BCrat']],'RN.intercept'=R.model$coefficients[['(Intercept)']],'Trait'='BCrat')),
          cbind(as.data.frame(t(summary(D.model)$coefficients[2,])),data.frame('RN.slope'=D.model$coefficients[['EI.BCdiv']],'RN.intercept'=D.model$coefficients[['(Intercept)']],'Trait'='BCdiv')),
          cbind(as.data.frame(t(summary(T.model)$coefficients[2,])),data.frame('RN.slope'=T.model$coefficients[['EI.totalAGS']],'RN.intercept'=T.model$coefficients[['(Intercept)']],'Trait'='totalAGS'))) -> RN.results
    RN.results$Genotype<-geno # keep track of the genotype
    RN.results$Site<-site # keep track of which site
    if(geno=='BHM' & site=='Ald') { RNs<-RN.results } # for first iteration, create new dataframe
    else {RNs<-rbind(RNs,RN.results)} # for subsequent iterations, add on to existing dataframe
    }  
  }
}

## merge results with subspecies info ##
RNs<-merge(RNs,ssp,by='Genotype')

## save reaction norms 
save(RNs,file='intermediate_data/RNs_Mean.RData')
#######  Calculate RN heights (i.e., trait value at 'average' block) within each site #######
# Calculate the average EI for each trait at each site
siteEImeans <- select(blkAGSmeans,blkBCrat,blkBCdiv,blktotalAGS,Site) %>% 
  group_by(Site) %>%
  summarize(meanEI.BCrat=mean(blkBCrat,na.rm=TRUE), # for each trait at each site, calculate average EI
            meanEI.BCdiv=mean(blkBCdiv,na.rm=TRUE),
            meanEI.totalAGS=mean(blktotalAGS,na.rm=TRUE)) %>%
  ungroup() %>% as.data.frame() %>% 
  gather(key='Trait',value='meanEI',starts_with('meanEI')) %>% # melt data frame into long format
  mutate(Trait=str_replace(Trait,'meanEI.','')) %>% # simplify trait names
  mutate(TraitSite=paste0(Trait,Site)) %>% # for merging purposes, make new column with trait-site combination (should be unique)
  select(-Site,-Trait) # get rid of unneeded columns

# Merge mean EIs into reaction norm dataframe and calculate heights:
RNs <- mutate(RNs,TraitSite=paste0(Trait,Site)) %>% # for merging purposes, make new column with trait-site combination
  merge(.,siteEImeans,by='TraitSite') %>% # merge meanEI values into RNs data frame
  mutate(RNheight=RN.intercept+RN.slope*meanEI) # evaluate each genotype's reaction norm at the mean EI for that trait at that site

## save reaction norms 
save(RNs,file='intermediate_data/RNs_Mean.RData')

#######  Plot reaction norms for each trait in each site #######
## Melt data into long format:
genoblockMeans.melted <- gather(genoblockMeans,key='Trait',value='genoblock.Mean',BCrat,totalAGS,BCdiv) %>%
  mutate(EI.Mean=ifelse(Trait=='BCrat',EI.BCrat, # match up EIs with the right traits
                        ifelse(Trait=='totalAGS',EI.totalAGS,EI.BCdiv))) %>%
  select(-EI.BCrat,-EI.totalAGS,-EI.BCdiv)  # remove unneeded columns

# Make the plot:
pdf(file='figures/Figure5.pdf',width=11,height=7)
filter(genoblockMeans.melted,!(Trait=='totalAGS' & genoblock.Mean > 0.3)) %>% # for plotting clarity, exclude one high-[AGS] outlier (genotype SAD in block s.02)
  ggplot(.,aes(x=EI.Mean,y=genoblock.Mean))+
  facet_wrap(Site~Trait,scales='free',ncol=3)+
  geom_point(aes(color=Ssp),alpha=0)+
  scale_color_manual(values=sspPalette,guide=FALSE)+
  geom_abline(data=RNs,mapping=aes(intercept=RN.intercept,slope=RN.slope,color=Ssp))+
  geom_abline(intercept=0,slope=1,color='blue',linetype='dashed',size=1)+ # by definition, the mean reaction norm has a slope of 1
  geom_vline(data=RNs,mapping=aes(xintercept = meanEI),linetype='dotted',size=1,color='purple')+
  xlab('EI (block pooled mean)')+ylab('Genotype-block pooled mean')+
  theme_classic()+
  theme(strip.background = element_blank(),strip.text=element_blank())+
  theme(axis.text.y=element_text(size=13),axis.text.x=element_text(size=13))+
  theme(panel.spacing = unit(1.5, 'lines'))
dev.off()
#######  Test for heterogeneity of reaction norm slopes for each trait (in each site) #######
for (site in c('Ald','Jam','Mah','Sil')) { # look at one site at a time
  genoblockMeans.subset <- filter(genoblockMeans, Site==site) # subset by site
  # perform ANOVA for each trait:
  BCrat.test <- lm(BCrat~Genotype*EI.BCrat,data=genoblockMeans.subset) %>% Anova(type='III') %>% as.data.frame() %>% mutate(Term=row.names(.),Trait='BCrat',Site=site) 
  BCdiv.test <- lm(BCdiv~Genotype*EI.BCdiv,data=genoblockMeans.subset) %>% Anova(type='III') %>% as.data.frame() %>% mutate(Term=row.names(.),Trait='BCdiv',Site=site) 
  totalAGS.test <- lm(totalAGS~Genotype*EI.totalAGS,data=genoblockMeans.subset) %>% Anova(type='III') %>% as.data.frame() %>% mutate(Term=row.names(.),Trait='totalAGS',Site=site)
  # Record stats for the 3 tests at this site in a data frame:
  site.results <- rbind(BCrat.test,BCdiv.test,totalAGS.test)
  if (site=='Ald') { # for first iteration, initialize brand new data frame:
    RN.het.results <- site.results }
  else { rbind(RN.het.results,site.results) -> RN.het.results } # otherwise add to existing data frame
}

# Adjust P-values to correct for multiple comparisons:
group_by(RN.het.results,Trait,Term) %>%
  mutate(P.adj=p.adjust(`Pr(>F)`,method='holm')) %>%
  ungroup() %>% as.data.frame() -> RN.het.results 

# View results:
RN.het.results %>% filter(str_detect(Term,'Genotype:EI.'),P.adj<0.05) # focus on Genotype*EI interactions to describe genetic variation for reaction norm slopes
"      Sum Sq Df  F value       Pr(>F)                 Term    Trait Site        P.adj
1 0.03690419 23 2.669063 0.0001576956 Genotype:EI.totalAGS totalAGS  Jam 0.0006307824
2 0.52978981 24 2.040350 0.0030818048    Genotype:EI.BCrat    BCrat  Mah 0.0123272192
3 0.03182203 24 1.905593 0.0069416469 Genotype:EI.totalAGS totalAGS  Mah 0.0208249406
4 0.43985694 24 1.954229 0.0055443743    Genotype:EI.BCrat    BCrat  Sil 0.0166331228"
# Conclusion: genetic variation for reaction norm slopes for:
### BC-ratio at Mah, Sil
### total [AGS] at Jam, Mah

####### ------ TESTING FOR SELECTION ON REACTION NORM SLOPES ------ #######
####### Check initial plant diameters (At time of transplanting) #######
filter(plants, init_diam==0,Survival=='N') %>% summary() # 141 plants were dead on arrival to the field; exclude these from selection analysis
# update DOA plants- make a separate Survival category for plants that were DOA:
mutate(plants,Survival=factor(ifelse(init_diam==0&Survival=='N','DOA',
                                     ifelse(Survival=='N','N','Y')))) -> plants

####### Calculate relative fitness for each genotype #######
## We already removed ambiguous fitness records from the dataset (above)
## Note whether each plant survived or died/went missing (ignore DOA plants):
filter(plants,Survival=='N') %>% dim() # 1404 plants out of 3859 died over the winter (3859 = 4000 - 141 DOA)

# Make data frame for use in genotypic selection analysis:
genofit<-filter(plants,Survival!='DOA') %>% #exclude DOA plants
  mutate(FrProd=ifelse(Survival=='Y',FrProd,NA)) %>% # don't include non-survivors in fecundity estimate
  mutate(Survival=ifelse(Survival=='Y',1,0)) %>% # recode survival as 1 or 0, for easier calculation of % survival
  group_by(Genotype,Site) %>% # calculate relative fitness of each genotype separately in each site 
  # get pct survival and mean fecundity for each genotype at each site:
  # (and also genotype's mean diameter at time of transplant- as a covariate)
  summarize(PctSurvival=mean(Survival,na.rm=TRUE),meanFrProd=mean(FrProd,na.rm=TRUE),meanDiam=mean(init_diam,na.rm=TRUE)) %>% 
  mutate(W=ifelse(PctSurvival>0,PctSurvival*meanFrProd,0)) %>% # multiply % survival by mean fecundity for each genotype in each site; also specify W=0 for the genotype with Fecundity = NA (i.e., BHS at Jam)
  ungroup %>% as.data.frame %>% 
  group_by(Site) %>% mutate(relW=W/mean(W,na.rm=TRUE),  # calculate relative fitness of each genotype at each site
                            relFrProd=meanFrProd/mean(meanFrProd,na.rm=TRUE), # calculate relative fecundity of each genotype at each site
                            relSurvival=PctSurvival/mean(PctSurvival,na.rm=TRUE)) %>% # calculate relative survival rate of each genotype at each site
  ungroup %>% as.data.frame %>%
  mutate(GenoSite=paste0(Genotype,Site)) # make new variable for merging purposes

####### Combine reaction norm data with fitness data #######
select(RNs,TraitSite,Site,Genotype,Ssp,RN.slope,RN.intercept,Trait,RNheight) %>%
  mutate(GenoSite=paste0(Genotype,Site)) %>% # make new variable for merging purposes
  merge(.,genofit,by='GenoSite',suffixes=c('','.y')) -> genofit 

####### Genotypic selection analysis on reaction norms ####### 
for(site in c('Ald','Jam','Mah','Sil')) {
  for (trait in c('BCrat','BCdiv','totalAGS')) {
    site.subset<-filter(genofit,Site==site,Trait==trait) # look within one site at a time, one trait at a time
    selmodel<-lm(relW~meanDiam+RNheight+RN.slope,data=site.subset)
    as.data.frame(summary(selmodel)$coefficients) %>% # store results
      mutate(Term=rownames(.),Trait=trait,Site=site) %>%
      filter(Term!='(Intercept)') -> results 
    if(site=='Ald' & trait=='BCrat') { RNsel.Fit<-results } # for first iteration, initialize new data frame
    else {RNsel.Fit<-rbind(RNsel.Fit,results)} # for other iterations, add on to existing dataframe
  }
}

# Correct P-values for multiple comparisons (testing for selection in 4 different sites):
rename(RNsel.Fit,P=`Pr(>|t|)`) %>% group_by(Term,Trait) %>%
  mutate(P_corrected=p.adjust(P,method='holm')) %>%
  ungroup %>% as.data.frame -> RNsel.Fit
write.table(RNsel.Fit,file='tables/RNselection_Fit.txt',sep='\\t',col.names=TRUE,row.names=FALSE)

####### Supp. Figure 9: Compare trait & environmental variances among sites #######
# Is within-site (among-block) heterogeneity different between Eastern habitats and Western habitats?
# Merge block-level vegetation/herbivory data with all block-level glucosinolate data and environmental data together in long format:
merge(veg,blkAGSmeans,by='Block') %>% 
  rename(PlantDiv=Div,Veg=BlkVeg,Herbivory=BlkHerbivory) %>%
  mutate(HabitatType=factor(ifelse(Site%in%c('Jam','Mah'),'E','W'))) %>% # Group sites into habitat types based on the subspecies endemic to that area
  gather(key='variable',value='value',pH,conductivity,NO3_N,P,K,Ca,Mg,S,Na,Fe,Zn,Mn,Cu,PlantDiv,Veg,Herbivory,blkBCrat,blktotalAGS,blkBCdiv) %>%
  mutate(variable=factor(variable)) -> blkdata.long

# Levene's test of residuals (after accounting for differences among sites)
ltests<-data.frame('F.value'=numeric(),'P'=numeric(),'Variable'=character())
res.blk.env<-data.frame('residuals'=numeric(),'Site'=factor(),'HabitatType'=factor(),'Variable'=factor())
for(var in levels(blkdata.long$variable)){
  blkdata.subset<-filter(blkdata.long,variable==var) # subset the data - look at one variable at a time
  model<-lm(value~Site,data=blkdata.subset) # fit a model- partition variance among sites. residuals = variation among blocks
  res<-as.data.frame(residuals(model)) %>% # store residuals
    cbind(.,select(blkdata.subset,Site,HabitatType)) %>% # add Site and Habitat Type data to residuals
    rename(residuals=`residuals(model)`) %>%
    mutate(Variable=var)
  res.blk.env<-rbind(res.blk.env,res)
  ltest<-leveneTest(residuals(model) ~ blkdata.subset$HabitatType) %>% # test for equal variances
    rename(P=`Pr(>F)`,F.value=`F value`) %>% mutate(Variable=var) %>% filter(!is.na(P)) %>% select(-Df)
  ltests<-rbind(ltests,ltest)
}
mutate(ltests,P.corrected=p.adjust(P,method='holm')) -> ltests
filter(ltests,P.corrected<0.05)

# Plot: 
pdf(file='figures/Supp_Fig9.pdf')
filter(res.blk.env) %>%
  ggplot(.,aes(x=residuals,y=HabitatType,color=HabitatType,fill=HabitatType,height=..density..))+
  geom_joy(stat='density',rel_min_height=0.01,scale=3)+
  facet_wrap(~Variable,scales='free_x')+
  scale_color_manual(values=c('red','black'))+
  scale_fill_manual(values=c('red','black'))+
  theme_classic()+
  theme(axis.title.y=element_blank())+
  theme(legend.position=c(0.9,0.1),legend.background=element_rect(fill='grey90'))+
  theme(legend.text=element_text(face='bold'),legend.title=element_text(face='bold'))+
  theme(strip.background=element_rect(fill='grey90'))+
  theme(strip.text=element_text(face='bold'))+theme(axis.text.y=element_blank())+
  theme(axis.text.x=element_text(angle=35,vjust=0.75))+
  theme(panel.spacing=unit(1,'lines'))
dev.off()

####### ------ Supplementary Analysis: Individual AGS compounds rather than emergent properties of GS profiles ------ #######
## Note: Block is implicitly nested within Site (by unique names for all levels of Block factor)
####### REML linear mixed model for conc1ME #######
GxE.Site.conc1ME<-lmer(sqrt(conc1ME)~Genotype*Site+DevStage+H+(1|Block)+(1|Block:Genotype)+(1|Batch),data=AGS.noBHS)
qqnorm(resid(GxE.Site.conc1ME));qqline(resid(GxE.Site.conc1ME))
plot(resid(GxE.Site.conc1ME)~fitted(GxE.Site.conc1ME)) # sqrt transformation improves residuals

r2.LMM(GxE.Site.conc1ME) # 0.65

anova(GxE.Site.conc1ME) # Test fixed effects:
PVEranef(GxE.Site.conc1ME) # Test random effects and % variance explained:

####### REML linear mixed model for conc1MP #######
GxE.Site.conc1MP<-lmer(sqrt(conc1MP)~Genotype*Site+DevStage+H+(1|Block)+(1|Block:Genotype)+(1|Batch),data=AGS.noBHS)
qqnorm(resid(GxE.Site.conc1MP));qqline(resid(GxE.Site.conc1MP))
plot(resid(GxE.Site.conc1MP)~fitted(GxE.Site.conc1MP))  # OK
hist(resid(GxE.Site.conc1MP),50) # looks OK

r2.LMM(GxE.Site.conc1MP) # 0.7003975

anova(GxE.Site.conc1MP) # Test fixed effects
PVEranef(GxE.Site.conc1MP) # Test random effects and % variance explained:

####### REML linear mixed model for conc6MSOH #######
GxE.Site.conc6MSOH<-lmer(sqrt(conc6MSOH)~Genotype*Site+DevStage+H+(1|Block)+(1|Block:Genotype)+(1|Batch),data=AGS.noBHS)
qqnorm(resid(GxE.Site.conc6MSOH));qqline(resid(GxE.Site.conc6MSOH))
plot(resid(GxE.Site.conc6MSOH)~fitted(GxE.Site.conc6MSOH))  
r2.LMM(GxE.Site.conc6MSOH) # 0.7193242

anova(GxE.Site.conc6MSOH) # Test fixed effects:
PVEranef(GxE.Site.conc6MSOH) # Test random effects and % variance explained:

####### REML linear mixed model for conc2OH1ME #######
GxE.Site.conc2OH1ME<-lmer(sqrt(conc2OH1ME)~Genotype*Site+DevStage+H+(1|Block)+(1|Block:Genotype)+(1|Batch),data=AGS.noBHS)
qqnorm(resid(GxE.Site.conc2OH1ME));qqline(resid(GxE.Site.conc2OH1ME))
plot(resid(GxE.Site.conc2OH1ME)~fitted(GxE.Site.conc2OH1ME))  
hist(resid(GxE.Site.conc2OH1ME),50) 

r2.LMM(GxE.Site.conc2OH1ME) # 0.7663687

anova(GxE.Site.conc2OH1ME) # Test fixed effects:
PVEranef(GxE.Site.conc2OH1ME) # Test random effects and % variance explained:
####### Combine results of all REML models and adjust p-values for 4 comparisons: #######
rbind(simplify.results(GxE.Site.conc1ME),
      simplify.results(GxE.Site.conc1MP),
      simplify.results(GxE.Site.conc6MSOH),
      simplify.results(GxE.Site.conc2OH1ME)) -> GSvar.results.conc
## Correct p-values for each term:
group_by(GSvar.results.conc,Term) %>%
  mutate(p.value.adj=p.adjust(p.value.raw,'holm')) %>%
  ungroup %>% as.data.frame -> GSvar.results.conc

## Save results:
save(GSvar.results.conc,file='tables/conc_stats_GSvar_models.RData')

## Save as table:
mutate(GSvar.results.conc,DenDF=round(DenDF),F.or.ChiSq=round(F.or.ChiSq,digits=2)) %>% 
  select(-p.value.raw) %>% 
  mutate(p.value.adj=signif(p.value.adj,digits=3)) %>% 
  write.table(file='tables/conc_stats_GSvar_models.txt',sep='\\t',row.names=FALSE,col.names=TRUE)

####### Extract LS means from REML linear mixed models (simultaneously back-transform) #######
# Combine all LS means into one data frame:
getLSmeans.sqrt(GxE.Site.conc1ME) %>% # the getLSmeans.sqrt() function simultaneously back-transforms
  rbind(.,getLSmeans.sqrt(GxE.Site.conc1MP)) %>%
  rbind(.,getLSmeans.sqrt(GxE.Site.conc6MSOH)) %>%
  rbind(.,getLSmeans.sqrt(GxE.Site.conc2OH1ME)) -> LSM.all.conc
####### Check for genetic correlation between compounds using Genotype LS means #######
## Extract all Genotype LS means
filter(LSM.all.conc,Term=='Genotype') %>% 
  select(Genotype,lsmean,Response) %>%
  spread(key=Response,value=lsmean) %>%
  rename(conc1ME=`sqrt(conc1ME)`,conc1MP=`sqrt(conc1MP)`,
         conc2OH1ME=`sqrt(conc2OH1ME)`,conc6MSOH=`sqrt(conc6MSOH)`) -> genoLSM.conc

# subset genotypes: only those with branched-chain functionality ("gain of function allele only"):
filter(genoLSM.conc,!Genotype%in%c('MAH','PAR','SAD')) -> genoLSM.GOFonly.conc

## Spearman rank correlations: all genotypes included
cor.test(~conc1MP+conc1ME,data=genoLSM.conc,method='spearman') #
cor.test(~conc6MSOH+conc1ME,data=genoLSM.conc,method='spearman') # 
cor.test(~conc1MP+conc6MSOH,data=genoLSM.conc,method='spearman') # 

cor.test(~conc1ME+conc2OH1ME,data=genoLSM.conc,method='spearman') # 
cor.test(~conc1MP+conc2OH1ME,data=genoLSM.conc,method='spearman') # 
cor.test(~conc6MSOH+conc2OH1ME,data=genoLSM.conc,method='spearman') # 

p.adjust(c(0.0002938,0.0004382,2.404e-6,0.02672,0.00101,0.003507),method='holm') # adjust p-values

## Pearson's product moment genetic correlations: excluding genotypes with no branched-chain functionality
cor.test(~conc6MSOH+conc1ME,data=genoLSM.GOFonly.conc,method='spearman') # p = 0.01315 rho= -0.54
cor.test(~conc1MP+conc1ME,data=genoLSM.GOFonly.conc,method='spearman') # p=0.006844 rho= 0.58
cor.test(~conc1MP+conc6MSOH,data=genoLSM.GOFonly.conc,method='spearman') # p=2.93e-6  rho= -0.82

cor.test(~conc1ME+conc2OH1ME,data=genoLSM.GOFonly.conc,method='spearman') # p=0.002827 rho= -0.63
cor.test(~conc1MP+conc2OH1ME,data=genoLSM.GOFonly.conc,method='spearman') # p=2.9e-6 rho= -0.82
cor.test(~conc6MSOH+conc2OH1ME,data=genoLSM.GOFonly.conc,method='spearman') # p=3.2e-5 rho= 0.79

p.adjust(c(0.01807,0.0117,2.93e-6,0.00249,2.9e-6,3.2e-5),method='holm')

####### Prepare LS means for plotting: #######
# First, rename response variables:
mutate(LSM.all.conc,Response=ifelse(Response=='sqrt(conc1ME)','conc1ME',
                                    ifelse(Response=='sqrt(conc1MP)','conc1MP',
                                           ifelse(Response=='sqrt(conc2OH1ME)','conc2OH1ME','conc6MSOH')))) %>%
### Merge with subspecies data and re-order by subspecies
full_join(.,ssp,by='Genotype') %>% 
  arrange(Ssp) %>%
  mutate(Genotype=factor(Genotype,Genotype))  -> LSM.all.conc

####### Supp. Figure 5a: Main-effect LS means #######
pdf(file='figures/Supp_Fig5a.pdf',width=9,height=5)
filter(LSM.all.conc,Term=='Genotype') %>% 
  ggplot(.,aes(x=Genotype,y=lsmean,color=Ssp,alpha=Ssp))+
  facet_wrap(~Response,ncol=2,scales='free')+
  geom_point(size=2,position=position_dodge(width=1))+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0.0,size=.5,position=position_dodge(width=1))+
  theme_classic()+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_blank())+
  theme(axis.title.y=element_blank(),axis.text.y=element_text(size=22))+
  theme(legend.title= element_text(size=30),legend.text=element_text(size=26))+
  theme(legend.key.height=unit(1.5,"lines"),legend.key.width=unit(2.5,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))+
  theme(strip.background = element_rect(fill="gray90",color="gray90"), strip.text=element_text(size=24,face="bold"))+
  geom_point(data=filter(LSM.all.conc,Term=='Site') %>% 
               mutate(Genotype=plyr::mapvalues(Site,from=c('Ald','Jam','Mah','Sil'),to=c('EGM','MIL','BHM','THT'))),
             aes(x=Genotype,y=lsmean,color=Site,alpha=Site),size=4,shape=17)+
  geom_errorbar(data=filter(LSM.all.conc,Term=='Site') %>% 
                mutate(Genotype=plyr::mapvalues(Site,from=c('Ald','Jam','Mah','Sil'),to=c('EGM','MIL','BHM','THT'))),
                aes(x=Genotype,ymin=lower.CL,ymax=upper.CL,color=Site),width=0.0,size=1.5,position=position_dodge(width=1))+
  scale_alpha_manual(values=c(1,0.3,1,1,1,0.3),guide=FALSE)+
  scale_color_manual(values=c(sitePalette[1],"red",sitePalette[2:4],"black"),guide=FALSE)
dev.off()

####### Supp. Figure 5b: Genotype*Site LS means #######
pdf(file='figures/Supp_Fig5b.pdf',width=9,height=5)
filter(LSM.all.conc,Term=='Genotype:Site') %>% 
  ggplot(.,aes(x=Genotype,y=lsmean,color=Site))+
  facet_wrap(~Response,ncol=2,scales='free')+
  geom_point(position=position_dodge(width=.9))+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0,position=position_dodge(width=.9))+
  geom_vline(xintercept=0.5+1:24,linetype='dotted',color='grey')+
  scale_color_manual(values=sitePalette,guide=FALSE)+
  ylab("micromoles/mg")+
  theme_classic()+
  theme(axis.title.x=element_text(size=24,face="bold"),axis.text.x=element_blank())+
  theme(axis.title.y=element_text(size=24,face="bold"),axis.text.y=element_text(size=14))+
  theme(strip.background = element_blank(),strip.text=element_blank())
dev.off()
