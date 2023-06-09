# SCRIPT USED FOR THE ESTIMATION OF HISTORICAL EXTRA-PAIR PATERNITY RATES
# IN A WESTERN HUMAN POPULATION IN FUNCTION OF DIFFERENT SOCIO-DEMOGRAPHIC FACTORS
# SEE PAPER "A historical-genetic reconstruction of human extra-pair paternity"
# by Maarten H.D. Larmuseau, Pieter van den Berg, Sofie Claerhout, Francesc Calafell, 
# Alessio Boattini, Leen Gruyters, Michiel Vandenbosch, Kelly Nivelle, Ronny Decorte, 
# Tom Wenseleers (Current Biology, 2019)
# all analyses were run using R version 3.5.3


#### 0. LOAD REQUIRED PACKAGES ####

library(maptools)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(scales)
library(afex)
library(mgcv)
library(broom)
library(car)
library(psych)
library(tidyr)
library(emmeans)
library(export)
# adjust working directory as appropriate to point to where Dryad data files are stored
setwd("~/Temp")


#### 1. LOAD SPATIAL DATA FOR SPATIAL PLOTS LATER ON & SET GRAPHICS OPTIONS ####

# load provincial regions
admin_bel_l0<-raster::getData('GADM',country='BEL',level=0) # country outlines
admin_bel_l1<-raster::getData('GADM',country='BEL',level=1)
admin_bel_l2<-raster::getData('GADM',country='BEL',level=2) # provinces
admin_bel_l2<-admin_bel_l2[,"NAME_2"]
names(admin_bel_l2)="NAME"
admin_bel_l2$NAME=c("Brussel","Antwerpen","Limburg","Oost-Vlaanderen","Vlaams Brabant","West-Vlaanderen","Waals Brabant","Henegauwen","Luik","Luxemburg","Namen")

admin_nl_l0<-raster::getData('GADM',country='NL',level=0)
admin_nl_l1<-raster::getData('GADM',country='NL',level=1) # provinces
admin_nl_l1<-admin_nl_l1[,"NAME_1"]
names(admin_nl_l1)="NAME"
admin_nl_l1_lakes <- subset(admin_nl_l1, admin_nl_l1$NAME  %in% c("Zeeuwse meren", "IJsselmeer", "Flevoland"))
admin_nl_l1 <- subset(admin_nl_l1, !admin_nl_l1$NAME  %in% c("Zeeuwse meren", "IJsselmeer","Flevoland"))
# we remove Zeeuwse meren en Ijsselmeer as they are water and Flevoland as it didn't exist yet
admin_nl_l2<-raster::getData('GADM',country='NL',level=2)

admin_fra_l0<-raster::getData('GADM',country='FRA',level=0)
admin_fra_l1<-raster::getData('GADM',country='FRA',level=1) # provinces
admin_fra_l2<-raster::getData('GADM',country='FRA',level=2)
admin_fra_l3<-raster::getData('GADM',country='FRA',level=3)
admin_fra_l4<-raster::getData('GADM',country='FRA',level=4)

admin_lux_l0<-raster::getData('GADM',country='LUX',level=0)
admin_lux_l1<-raster::getData('GADM',country='LUX',level=1) # provinces

admin_de_l0<-raster::getData('GADM',country='DE',level=0)
admin_de_l1<-raster::getData('GADM',country='DE',level=1) # provinces

# define which provinces to include
brabant = c("Antwerpen","Vlaams Brabant","Brussel","Noord-Brabant","Waals Brabant")
holland = c("Noord-Holland","Zuid-Holland","Utrecht", "Gelderland") # I added Utrecht & Gelderland
limburg = c("Limburg")
vlaanderen = c("Oost-Vlaanderen","West-Vlaanderen","Zeeland","Nord-Pas-de-Calais")
allprovinces = c(brabant, holland, limburg, vlaanderen)

# combined SpatialPolygonsDataFrame at province level for Belgium & Netherlands
admin_nl_l1B <- admin_nl_l1
admin_nl_l1B$NAME <- c("Drenthe", "Friesland", "Gelderland", "Groningen", "Nederlands Limburg",
                       "Noord-Brabant", "Noord-Holland", "Overijssel", "Utrecht", "Zeeland",
                       "Zuid-Holland")
admin_bel_l2B <- admin_bel_l2
admin_bel_l2B$NAME <- c("Brussel", "Antwerpen", "Limburg", "Oost-Vlaanderen", "Vlaams Brabant",
                        "West-Vlaanderen", "Waals Brabant", "Henegauwen", "Luik", "Luxemburg",
                        "Namen")
admin_prov_all <- spRbind(spChFIDs(admin_bel_l2B, make.unique(admin_bel_l2B$NAME)),
                          spChFIDs(admin_nl_l1B, make.unique(admin_nl_l1B$NAME)))
admin_prov_all$NAME # all Belgian & Dutch provinces
admin_prov_sel <- subset(admin_prov_all, admin_prov_all$NAME  %in% 
                           c("Antwerpen","Vlaams Brabant","Brussel","Waals Brabant","Noord-Brabant",
                             "Noord-Holland","Zuid-Holland","Utrecht","Gelderland","Limburg",
                             "Nederlands Limburg","Oost-Vlaanderen","West-Vlaanderen","Zeeland")) 
# provinces with EPP data
# we combine polygons of Vlaams Brabant, Brussel & Waals Brabant in new region "Brabant"
# since we have density estimates over time only for Brabant as a whole
grps <- c("Brabant","Antwerpen","Limburg",
          "Oost-Vlaanderen","Brabant","West-Vlaanderen",
          "Brabant","Gelderland","Nederlands Limburg",
          "Noord-Brabant","Noord-Holland","Utrecht",
          "Zeeland","Zuid-Holland")
admin_prov_sel <- unionSpatialPolygons(admin_prov_sel, grps)
mydata <- data.frame(NAME=row.names(admin_prov_sel))
rownames(mydata) <- row.names(admin_prov_sel)
admin_prov_sel <- SpatialPolygonsDataFrame(admin_prov_sel, mydata)

# lists of cities included in agglomeration of Antwerp or Brussels
antwagglo=c("Antwerpen","Antwerpen-Berchem","Berchem","Berendrecht","Borgerhout","Deurne","Ekeren","Hoboken","Merksem","Wilrijk") # cities included in agglomeration Antwerpen
bruagglo=c("Brussel","Sint-Jans-Molenbeek","Sint-Joost-ten-Node","Schaarbeek",
           "Elsene","Etterbeek","Sint-Gillis","Anderlecht","Koekelberg",            
           "Laken","Sint-Agatha-Berchem","Ganshoren","Jette","Evere",
           "Vorst", "Ukkel", "Watermaal-Bosvoorde", "Oudergem",              
           "Sint-Pieters-Woluwe", "Sint-Lambrechts-Woluwe") # cities included in agglomeration Brussel
# main cities for which we have most EPP data (n>=20) 
main_cities_urban = c("Amsterdam","Antwerpen","Ardooie","Beveren","Boom","Breda","Brugge","Brussel","Den Haag","Eindhoven","Gent","Izegem","Kortrijk","Leuven","Maastricht","Menen","Meulebeke","Oostende","Roeselare","Rotterdam","Ruiselede","Tielt","Tilburg","Torhout","Turnhout","Wingene") # we consider Boom as urban even though it wasn't in Infra-Clio database, based on pop size in 1850 (>5000) & density
main_cities_rural = c("Aartrijke","Beveren-Roeselare","Eernegem","Essen","Geel","Geluwe","Gits","Gullegem","Herentals","Herenthout","Ichtegem","Iddergem","Idegem","Kalmthout","Langdorp","Lichtervelde","Liedekerke","Maldegem", "Meer","Mortsel","Oisterwijk","Oudenburg","Pittem","Rijkevorsel","Waarschoot","Wuustwezel","Zedelgem","Zichem","Zwevezele")

# set plot options/colours for ggmaps
al <<- 1 
aloutl <<- 1
alprov <<- 0.7 # alpha coloured province fill, higher vals=darker col/less transparent
provcols <<- c("darkorchid","darkgreen","darkorange","turquoise4")
provfill <<- sapply(provcols, function(x) muted(x,c=70,l=60)) 
provoutl <<- sapply(provcols, function(x) muted(x,c=70,l=40)) 
allprovoutl <<- adjustcolor(c("grey75"),alpha.f=al)
countryoutl <<- adjustcolor(c("grey75"),alpha.f=al)
countryfill <<- muted("antiquewhite",l=95,c=1)
seafill <<- muted("dodgerblue",l=40,c=80) 
ruralfill <<- muted("green4",l=40,c=80) # for rural areas
legcols <<- c(rgb(0,0,255,maxColorValue = 255), # for legend, blue to red gradient
              rgb(0,124,244,maxColorValue = 255),
              rgb(219,0,189, maxColorValue = 255),
              rgb(230,0,133, maxColorValue = 255),
              rgb(255,0,0, maxColorValue = 255))  
legcols <<- c("blue","darkorchid","red") # for legend
textcol <<- "white"
longlims <<- c(2.2,7.2)
latlims <<- c(50.3,53.5)
dfboxsea <- data.frame( # to draw the sea
  xmin = longlims[[1]], xmax = longlims[[2]], ymin = latlims[[1]], ymax=latlims[[2]] )
lw <<- 1.2

# define a ggplot2 theme for our plots
mytheme = theme_bw() + theme(axis.line = element_line(colour = "black", size=0.5),
                             panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(),
                             panel.background = element_rect(colour = "black", size=0.5, 
                                                             fill=rgb(red=237,green=240,blue=240,maxColorValue=255)), # muted("steelblue1", l = 94, c = 10)),
                             panel.border = element_rect(colour = "black", size=0.5),
                             plot.background = element_rect(fill = "transparent"),
                             axis.text=element_text(size=16, colour = "black"),
                             axis.title=element_text(size=22, colour = "black"),  # face="bold"
                             legend.key=element_blank(), legend.background=element_blank(),
                             legend.justification = c(0, 1), legend.position = c(0.1, 0.98),
                             legend.text=element_text(size=12), legend.title=element_text(size=14),
                             legend.spacing.y = unit(-0.09, "cm"),
                             axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), colour = "black"),
                             axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), colour = "black"),
                             axis.text.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0), colour = "black"),
                             axis.text.x = element_text(margin = margin(t = 7, r = 0, b = 0, l = 0), 
                                                        angle = 0, hjust = 0.5,
                                                        colour = "black"),
                             axis.ticks.length = unit(.25, "cm"),
                             axis.ticks = element_line(colour = "black", size = 0.5) )
colrs <<- c(rgb(red = 220, green = 0, blue = 0, maxColorValue = 255),
            adjustcolor("royalblue2",alpha.f=1),
            adjustcolor("green3",alpha.f=1))



#### 2. DEFINE FUNCTION TO CALCULATE CENSORED BINOMIAL GLM (OR GAM) TO ESTIMATE THE EPP RATE ####
# Based on iteratively readjusting the observed nr of mismatches using an Empirical Bayes logic
# as in Lasmuseau et al. (2013) Proc B paper, but now allowing for arbitrarily complex binomial GLM 
# or GAM model, and using fractionate Bernouilli events to encode uncertainty in where EPP events occurred.
# niter=3 seems enough for convergence, but we use niter=10 to be on safe side

# The input data are expected to also contain columns maxn, NPERPSEUDOPAIR and PSEUDOPAIR, which
# contain respectively the theoretical maximum number of EPP events within a given genealogical pseudopair (maxn),
# the total number of meiosis that separate the two sampled individuals of a given pseudopair (NPERPSEUDOPAIR),
# and pseudopair ID (PSEUDOPAIR).

censbinomfit = function(form=as.formula(cbind(EPP,1) ~ 1), 
                        data, niter=10, type="glm") {
  require(afex)
  require(mgcv)
  dataupd <<- data
  set.seed(1)
  if (type=="glm") {fitorig=suppressWarnings(glm(formula=form,data=dataupd,family=binomial))} else {fitorig=suppressWarnings(gam(formula=form,data=dataupd,family=binomial))}
  if (type=="glm") {compl=complete.cases(fitorig$data)} else {compl=complete.cases(fitorig$model)} # take into account possible missing values
  mismatchidxs = which(data[compl,]$EPP>0) # ancestors in pairs with mismatch
  nmismatches = length(mismatchidxs) 
  for (iter in 1:(niter-1)) {
    epp = dataupd$EPP[compl][mismatchidxs]
    maxn = data[compl,]$maxn
    nsperpair = data[compl,]$NPERPSEUDOPAIR
    pseudopair = data[compl,]$PSEUDOPAIR 
    if (type=="glm") {fit1=suppressWarnings(glm(formula=form,data=dataupd,family=binomial))} else {fit1=suppressWarnings(gam(formula=form,data=dataupd,family=binomial))}
    preds=predict(fit1, type="response") # predicted EPP rates
    # now do Empirical Bayes update of EPP given estimated p's for all cases where there was a mismatch,
    # to allow for posibility of there having been more than 1 mismatch within a given pair
    for (j in 1:nmismatches) {
      pspair = pseudopair[[mismatchidxs[[j]]]]
      pest = plogis(mean(qlogis(preds[pspair==pseudopair]))) # we calculate avg pest over PSEUDOPAIR (calculated on logit scale & then backtransformed)
      n = maxn[[mismatchidxs[[j]]]] # this is the maximum possible nr of mismatches 
      epp[j] = sum( (dbinom(1:n,n,pest)/(1-dbinom(0,n,pest))) * (1:n) ) / nsperpair[[mismatchidxs[[j]]]] # updated epp to take into account nonsampling
    }
    dataupd$EPP[mismatchidxs] = epp # updated to take into account nonsampling
  }
  dataupd$EPP[dataupd$EPP>1]=1
  if (type=="glm") {fit1=suppressWarnings(glm(formula=form,data=dataupd,family=binomial))} else {fit1=suppressWarnings(gam(formula=form,data=dataupd,family=binomial))}
  preds=plogis(predict(fit1)) # predicted EPP rates
  dataupd$pred=preds
  return(list(fit=fitorig,fitupd=fit1,dataupd=dataupd)) # fitorig=original GLM fit, fitupd=fit with data that are corrected for nondetection
}



#### 3. TEST WHICH MODELS GIVE BEST AIC & COMPARE AGAINST BEST MODEL ####

data_dryad = read.csv("data_dryad.csv")
colnames(data_dryad)
# VARIABLES INCLUDED IN DATASET:
# PAIR = genealogical pair of presumed patrilineally related men that were
#        Y-chromosome genotyped
# PSEUDOPAIR = genealogical pseudopair; genealogical pairs were split up in
#             pseudopairs if additional genetic evidence was available that 
#             allowed us to infer that EPP events were absent in given 
#             parts of the genealogy
# SURNAME = surname, anonymyzed to conform to privacy regulations
# MISMATCH = if there was a Y chromosome mismatch within a given pseudopair or not
# NPERPAIR = total number of meioses that separate a given genealogical pair
# NPERPSEUDOPAIR = total number of meioses that separate a given genealogical pseudopair
# EPP = the observed probability for an EPP event to have happened, 
#       encoded as a fractionate Bernouilli event (=MISMATCH/NPERPSEUDOPAIR)
# maxn = theoretical maximum number of EPP events that could have occurred within
#        a given pseudopair
# YEAR = year of birth, rounded off to nearest decade to conform to privacy regulations
# SOCIAL_CLASS = socioeconomic class based on occupation of father 
#               (Farmers, Low incompe, Medium or high income or missing)
# DENSITY = density (inhabitants/km^2) in the cite of birth at birth
# URBAN_RURAL = RURAL or URBAN, i.e. whether the city of birth could be considered 
#               to lie in an urban region or not; is in the Infra Clio database, 
#               cities that had more than 5,000 inhabitants in 1850 were considered 
#               urban
# COUNTRY = country of birth ("Belgie", "Nederland" or "outside_studyarea")
# PROVINCE = province where city of birth was located ("Antwerpen", "Brabant", 
#            "Gelderland", "Limburg", "Nederlands Limburg", "Noord-Brabant",     
#            "Noord-Holland", "Oost-Vlaanderen", "Utrecht", "West-Vlaanderen",
#            "Zeeland", "Zuid-Holland", "outside_studyarea")  
# SAMPLING_CAMPAIGN = which of the two sampling campaigns the data were derived from

# NOTE: YEAR OF BIRTH WAS ROUNDED TO NEAREST DECADE, SURNAME WAS ANONYMIZED AND
# AGE_FATHER & AGE_MOTHER & ASSOCIATED MODELS WERE REMOVED TO CONFORM TO PRIVACY REGULATIONS
# DUE TO THE ROUNDING OFF OF YEAR OF BIRTH THE OUTPUT OF THE ANALYSES BELOW MAY
# DIFFER VERY SLIGHTLY FROM THE ONES IN THE PAPER; THE ANALYSIS OUTPUTS GIVEN BELOW
# AS COMMENTS WERE OBTAINED USING THE ORIGINAL DATA THAT WAS NOT ROUNDED OFF

fit_0 = censbinomfit(form=as.formula( 
  cbind(EPP,1) ~
    SOCIAL_CLASS
), data=data_dryad, type="glm", niter=10)
fit_1 = censbinomfit(form=as.formula( # BEST MODEL
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    log2(DENSITY)
), data=data_dryad, type="glm", niter=10)
fit_2 = censbinomfit(form=as.formula( 
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    DENSITY 
), data=data_dryad, type="glm", niter=10)
fit_3 = censbinomfit(form=as.formula(
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    URBAN_RURAL 
), data=data_dryad, type="glm", niter=10)
fit_4 = censbinomfit(form=as.formula( 
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    YEAR 
), data=data_dryad, type="glm", niter=10)
fit_7 = censbinomfit(form=as.formula( 
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    COUNTRY 
), data=data_dryad, type="glm", niter=10)
fit_8 = censbinomfit(form=as.formula( 
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    PROVINCE 
), data=data_dryad, type="glm", niter=10)
fit_9 = censbinomfit(form=as.formula( 
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    SAMPLING_CAMPAIGN 
), data=data_dryad, type="glm", niter=10)
fit_10 = censbinomfit(form=as.formula( 
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    log2(DENSITY) +
    YEAR 
), data=data_dryad, type="glm", niter=10)
fit_13 = censbinomfit(form=as.formula( 
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    log2(DENSITY) +
    COUNTRY
), data=data_dryad, type="glm", niter=10)
fit_14 = censbinomfit(form=as.formula( 
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    log2(DENSITY) +
    PROVINCE
), data=data_dryad, type="glm", niter=10)
fit_15 = censbinomfit(form=as.formula( 
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    log2(DENSITY) +
    SAMPLING_CAMPAIGN
), data=data_dryad, type="glm", niter=10)

# Let's now calculate deltaAICs, evidence ratios and Akaike weights for this list of models 
# cf. https://cran.r-project.org/web/packages/ESTER/vignettes/ESTER.html
nmodels=12
AICdf=AIC(fit_0$fitupd,fit_1$fitupd,fit_2$fitupd,fit_3$fitupd,fit_4$fitupd,
          fit_7$fitupd,fit_8$fitupd,fit_9$fitupd,fit_10$fitupd,
          fit_13$fitupd,fit_14$fitupd,fit_15$fitupd)
AICdf$deltaAIC=sapply(1:nmodels, function (i) AICdf[i,"AIC"]-min(AICdf[,"AIC"])) # delta AICs compared to best model
AICdf$evidence_ratio = exp(0.5*AICdf$deltaAIC) # evidence ratios
AICdf$AICweights = sapply(1:nmodels, function (i) exp(-0.5*AICdf$deltaAIC[i])/sum(exp(-0.5*AICdf$deltaAIC))) # Akaike weights
AICdf$model = c("SOCIOECONOMIC_CLASS",
                "SOCIOECONOMIC_CLASS+log2(DENSITY)",
                "SOCIOECONOMIC_CLASS+DENSITY",
                "SOCIOECONOMIC_CLASS+RURALvsURBAN",
                "SOCIOECONOMIC_CLASS+YEAR_OF_BIRTH",
                "SOCIOECONOMIC_CLASS+COUNTRY",
                "SOCIOECONOMIC_CLASS+PROVINCE",
                "SOCIOECONOMIC_CLASS+SAMPLING_CAMPAIGN",
                "SOCIOECONOMIC_CLASS+log2(DENSITY)+YEAR_OF_BIRTH",
                "SOCIOECONOMIC_CLASS+log2(DENSITY)+COUNTRY",
                "SOCIOECONOMIC_CLASS+log2(DENSITY)+PROVINCE",
                "SOCIOECONOMIC_CLASS+log2(DENSITY)+SAMPLING_CAMPAIGN")
AICdf_ordered = AICdf[order(AICdf$AICweights,decreasing=T),]

# also add p values based on Type II Anovas
fits = list(fit_0$fitupd, fit_1$fitupd, fit_2$fitupd, fit_3$fitupd, fit_4$fitupd, 
            fit_7$fitupd, fit_8$fitupd, fit_9$fitupd, fit_10$fitupd, 
            fit_13$fitupd, fit_14$fitupd, fit_15$fitupd)
fits_Anovas = do.call(rbind,lapply(1:length(fits), function(fitnr) data.frame(fitnr=fitnr-1, data.frame(tidy(Anova(fits[[fitnr]],type="II"))))))
fits_Anovas$term[fits_Anovas$term=="log2(DENSITY)"]="DENSITY"
fits_Anovas$term[fits_Anovas$term=="YEAR"]="YEAR_OF_BIRTH"
AICdf_ordered=cbind(AICdf_ordered, spread(fits_Anovas[,-c(3,4)], term, p.value)[order(AICdf$AICweights,decreasing=T),])
AICdf_ordered[is.na(AICdf_ordered)]=""
write.csv(AICdf_ordered, "Table AIC values.csv")

# check for collinearity among predictors for best model using generalized variance inflation factors
vif(fit_1$fitupd) # generalized inflation factors are close to 1 so there is no strong collinearity
#                   GVIF Df GVIF^(1/(2*Df))
# SOCIAL_CLASS  1.034657  3        1.005694
# log2(DENSITY) 1.034657  1        1.017181


table(data_dryad$SOCIAL_CLASS,(data_dryad$YEAR>=1800&data_dryad$YEAR<=2020))
#                         FALSE TRUE
# Farmers                   146  994
# Low income                 53  871
# Medium or high income     131  709
# missing                  1666 2248
table(data_dryad$SOCIAL_CLASS,(data_dryad$YEAR>=1800&data_dryad$YEAR<=2020))/rowSums(table(data_dryad$SOCIAL_CLASS,(data_dryad$YEAR>=1800&data_dryad$YEAR<=2220)))
#                              FALSE       TRUE
# Farmers                 0.12807018 0.87192982
# Low income              0.05735931 0.94264069
# Medium or high income   0.15595238 0.84404762
# missing                 0.42565151 0.57434849
table(data_dryad$SOCIAL_CLASS,data_dryad$URBAN_RURAL)
#                         RURAL URBAN
# Farmers                   878   262
# Low income                528   396
# Medium or high income     482   358
# missing                  2491  1423
table(data_dryad$SOCIAL_CLASS,data_dryad$URBAN_RURAL)/rowSums(table(data_dryad$SOCIAL_CLASS,data_dryad$URBAN_RURAL))
#                             RURAL     URBAN
# Farmers                 0.7701754 0.2298246
# Low income              0.5714286 0.4285714
# Medium or high income   0.5738095 0.4261905
# missing                 0.6364333 0.3635667



#### 4. FINAL BEST MODEL BASED ON AIC: MODEL WITH SOCIAL CLASS & DENSITY ####

fit_1 = censbinomfit(form=as.formula(
  cbind(EPP,1) ~
    SOCIAL_CLASS +
    log2(DENSITY)
), data=data_dryad, type="glm", niter=10)

max(fit_1$dataupd$EPP-data_dryad$EPP) # max amount of empirical Bayes adjustment=0.02049313

mean(fit_1$dataupd$EPP[data_dryad$EPP!=0]/data_dryad$EPP[data_dryad$EPP!=0]) # mean amount of empirical Bayes adjustment=factor of 1.11
range(fit_1$dataupd$EPP[data_dryad$EPP!=0]/data_dryad$EPP[data_dryad$EPP!=0]) # 1-1.23

mean((fit_1$dataupd$EPP-data_dryad$EPP)[data_dryad$EPP>0]) # avg amount of empirical Bayes adjustment=0.007511701
mean(((fit_1$dataupd$EPP-data_dryad$EPP)/data_dryad$EPP)[data_dryad$EPP>0]) # avg amount of empirical Bayes adjustment=11.3%
AIC(fit_1$fitupd) # 225.5675
range(predict(fit_1$fitupd,type="response")) # 0.003998653 0.059443548
# predicted range for fathers, medium/high income and low income:
range(predict(fit_1$fitupd,type="response")[data_dryad$SOCIAL_CLASS=="Farmers"]) # 0.005083179 0.016942769
range(predict(fit_1$fitupd,type="response")[data_dryad$SOCIAL_CLASS=="Medium or high income"]) # 0.003998653 0.014250723
range(predict(fit_1$fitupd,type="response")[data_dryad$SOCIAL_CLASS=="Low income"]) # 0.01833567 0.05944355
emmeans(fit_1$fitupd,~1,type="response") # average predicted EPP rate
# 1             prob      SE df  asymp.LCL asymp.UCL
#      overall 0.016 0.00231 Inf    0.0121    0.0213

# Model coefficient table:
summary(fit_1$fitupd) 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)Farmers                  -5.45453    0.48250 -11.305  < 2e-16 ***
#   SOCIAL_CLASSLow income               1.29632    0.35781   3.623 0.000291 ***
#   SOCIAL_CLASSMedium or high income   -0.17782    0.48620  -0.366 0.714563    
#   SOCIAL_CLASSmissing                  0.21410    0.34552   0.620 0.535490    
#   log2(DENSITY)                        0.10508    0.04966   2.116 0.034348 *  

# Model Anova table:
Anova(fit_1$fitupd,type="III")
# Analysis of Deviance Table (Type III tests)
# 
#  Response: cbind(EPP, 1)
#                 LR Chisq Df Pr(>Chisq)    
#   SOCIAL_CLASS   27.2872  3  5.125e-06 ***
#   log2(DENSITY)   4.4017  1     0.0359 *  

# Planned comparisons for comparisons among different social classes: 
contrast(emmeans(fit_1$fitupd,~SOCIAL_CLASS),method="pairwise", adjust="none") # planned contrasts
# contrast                           estimate    SE  df z.ratio p.value
# Farmers - Low income                 -1.296 0.358 Inf -3.623  0.0003 
# Farmers - Medium or high income       0.178 0.486 Inf  0.366  0.7146 
# Farmers - missing                    -0.214 0.346 Inf -0.620  0.5355 
# Low income - Medium or high income    1.474 0.409 Inf  3.604  0.0003 
# Low income - missing                  1.082 0.224 Inf  4.828  <.0001 
# Medium or high income - missing      -0.392 0.396 Inf -0.989  0.3228 


# PLOT MODEL PREDICTIONS
df=data_dryad
preds=predict(fit_1$fitupd, type="response")
df$preds=preds  
df=rbind(df[df$SOCIAL_CLASS=="missing",],
         df[df$SOCIAL_CLASS=="Medium or high income",],
         df[df$SOCIAL_CLASS=="Farmers",],
         df[df$SOCIAL_CLASS=="Low income",]
)
df$SOCIAL_CLASS=factor(df$SOCIAL_CLASS,levels=c("Low income","Medium or high income","Farmers","missing"))


# exploratory plot: violin plots of distribution of density in function of social class
datatoplot=data_dryad
datatoplot=datatoplot[datatoplot$SOCIAL_CLASS!="missing",]
datatoplot$SOCIAL_CLASS=droplevels(datatoplot$SOCIAL_CLASS)
datatoplot$SOCIAL_CLASS=factor(datatoplot$SOCIAL_CLASS, levels=c("Farmers","Medium or high income","Low income"), 
                               labels=c("Farmer","Med-High","Low"))
qplot(x=SOCIAL_CLASS, y=log2(DENSITY), data=datatoplot, fill=SOCIAL_CLASS, geom="violin") + mytheme + xlab("Social class") + ylab("Log2(density)") +  scale_fill_manual(name="SOCIAL CLASS",labels=c("Farmers","Medium or high income","Low income"),             values=c(colrs[[3]],colrs[[2]],colrs[[1]],adjustcolor("grey",alpha.f=1))) + guides(fill=FALSE)
graph2png(file="Violin plot_density in function of social class.png",width=9,height=6)
graph2ppt(file="Violin plot_density in function of social class.pptx",width=9,height=6)

# plot with model predictions
plot1 = qplot(data=df[df$YEAR>1380,],x=YEAR,y=100*preds,
              xlab=expression("Year of birth"),ylab="Extra-pair paternity frequency (%)",
              shape=I(21),stroke=I(0.5),
              fill=SOCIAL_CLASS,
              size=log10(DENSITY)) + mytheme +
  scale_fill_manual(name="SOCIAL CLASS",
                    labels=c("Low income","Medium or high income","Farmers","Unknown"),
                    values=c(colrs[[1]],
                             colrs[[2]],
                             colrs[[3]],
                             adjustcolor("grey",alpha.f=1))) +
  scale_radius(name=expression ("DENSITY (per "*km^2*")"),
               breaks=c(5,4,3,2,1),
               labels=c("100,000","10,000","1,000","100","10"),range=c(0.5,6),trans="identity") + 
  coord_cartesian(ylim=c(0,6.2),xlim=c(1380,2020),expand=FALSE) +  
  guides(size = guide_legend(override.aes = list(alpha = 1, colour=I("black"), fill=I("black"), linetype = 0)),
         fill = guide_legend(override.aes = list(linetype = 0, size=4)))


plot1
graph2png(file="final binom GLM model_model predictions in function of CLASS DENSITY by YEAR.png",width=11,height=6)
graph2ppt(file="final binom GLM model_model predictions in function of CLASS DENSITY by YEAR.pptx",width=11,height=6)


# PLOT PARTIAL EFFECTS OF DIFFERENT FACTORS BASED ON EXPECTED MARGINAL MEANS
# partial effect of SOCIOECONOMIC CLASS (with 95% confidence bands, ie one-sided 95% confidence intervals)
dfemmeansCLASS=data.frame(summary(emmeans(fit_1$fitupd,~SOCIAL_CLASS),
                                  weights="proportional", type="response", level=0.9))
dfemmeansCLASS 
# SOCIAL_CLASS                     prob          SE  df   asymp.LCL  asymp.UCL
# 1                 Farmers 0.011554383 0.003724880 Inf 0.006130654 0.02167180
# 2              Low income 0.040983144 0.007140470 Inf 0.029062192 0.05750432
# 3   Medium or high income 0.009690321 0.003576135 Inf 0.004691610 0.01990841
# 4                 missing 0.014273604 0.002078099 Inf 0.010724407 0.01897486
dfemmeansCLASS$SOCIAL_CLASS=factor(dfemmeansCLASS$SOCIAL_CLASS,
                                   levels=c("missing","Farmers","Medium or high income","Low income"),
                                   labels=c("Unknown","Farmer","Medium or high income","Low income"))
dfemmeansCLASS=dfemmeansCLASS[dfemmeansCLASS$SOCIAL_CLASS!="Unknown",]
dfemmeansCLASS$SOCIAL_CLASS=droplevels(dfemmeansCLASS$SOCIAL_CLASS)
plot2 = qplot(data=dfemmeansCLASS, x=SOCIAL_CLASS,
              y=100*prob,
              ylab="Extra-pair paternity frequency (%)",
              xlab="Social class",fill=SOCIAL_CLASS,geom="col",width=0.75) +
  geom_linerange(aes(x=as.numeric(SOCIAL_CLASS),ymin=100*asymp.LCL,ymax=100*asymp.UCL),lwd=0.5) +
  scale_y_continuous(breaks=c(0:8), expand=c(0,0)) +
  scale_fill_manual(name="SOCIAL CLASS",
                    values=c(
                      colrs[[3]],
                      colrs[[2]],
                      colrs[[1]])) +
  scale_x_discrete(labels= c("Farmer","Med-High","Low")) +
  coord_cartesian(ylim=c(0,6)) +
  mytheme 
plot2
graph2png(file="final binom GLM model_partial effect SOCIOECONOMIC CLASS.png",width=11,height=6)
graph2ppt(file="final binom GLM model_partial effect SOCIOECONOMIC CLASS.pptx",width=11,height=6)

# partial effect DENSITY (with 95% confidence bands, ie one-sided 95% confidence intervals)
dfemmeansDENS=data.frame(summary(emmeans(fit_1$fitupd,~DENSITY,
                                         at=list(DENSITY=10^seq(0,4.03,length.out=1000)),
                                         weights="proportional", type="response"), level=0.9)) # we show 95% confidence bounds here = one-sided 95% confidence limits
plot3 = qplot(data=dfemmeansDENS,x=sqrt(DENSITY),y=100*prob,ylab="Extra-pair paternity frequency (%)",
              xlab="Density (per km2)",col=I("red3"),geom="line",lwd=I(1), group=factor(df)) +
  geom_ribbon(aes(ymin=100*asymp.LCL,ymax=100*asymp.UCL),alpha=I(0.25),fill=I("red3"),colour=NA,lwd=0) + 
  geom_point(aes(colour=factor(df))) +
  scale_colour_manual(values=c(adjustcolor("white",alpha.f=0)), name=NULL, labels=NULL) +
  # scale_y_continuous(breaks=c(0:4)) +
  scale_x_continuous(breaks=sqrt(c(10,100,1000,10000)), labels=c("10","100","1,000","10,000")) +
  coord_cartesian(ylim=c(0,3.4), expand=FALSE) +  # xlim=c(log10(100), log10(15000)
  mytheme 
plot3
graph2png(file="final binom GLM model_partial effect DENSITY.png",width=11,height=6)
graph2ppt(file="final binom GLM model_partial effect DENSITY.pptx",width=11,height=6)

min(data_dryad$DENSITY) # 2.129112
max(data_dryad$DENSITY) # 12571.38
data.frame(summary(emmeans(fit_1$fitupd,~DENSITY,
                           at=list(DENSITY=c(min(data_dryad$DENSITY), max(data_dryad$DENSITY))),
                           weights="proportional", type="response"), level=0.9))
# predicted EPP for lowest & highest density: 0.6% & 2.3%
# DENSITY        prob          SE  df  asymp.LCL  asymp.UCL
# 1     2.129112 0.006284988 0.002390212 Inf 0.00335886 0.01173026
# 2 12571.376155 0.023046634 0.006251001 Inf 0.01472199 0.03590696



# final multipanel figure Fig. 2
p1 = ggplotGrob(plot1)
p2 = ggplotGrob(plot2+theme(axis.title.y=element_blank(), legend.position="none"))
p3 = ggplotGrob(plot3+theme(axis.title.y=element_blank(), legend.position="none"))
grid.arrange(p1, arrangeGrob(p2, p3, ncol=1, heights=c(1,1)), 
             ncol=2, widths=c(2,1))
graph2png(file="final binom GLM model_multipanel figure.png",dpi=1200,width=12,height=6)
graph2ppt(file="final binom GLM model_multipanel figure.pptx",width=12,height=6)



# PLOT PREDICTIONS OF MODEL PATIOTEMPORALLY, USING HISTORICAL DENSITY ESTIMATES FOR CITIES & PROVINCES (RURAL DENSITIES)
# we plot predictions for each year between 1750 and 1950 individually and then combine all the frames in a video

years = seq(1750, 1950, by=1) # years to make predictions for

# make predictions for farmers in rural areas in provinces
# first we load our historical estimates of the rural and averall average densities 
# in each province through time (period 1750-1950) [for details and data sources see paper]: 
predsprovstoplot = read.csv("data_dryad_densities_provinces.csv")
predsprovstoplot$DENSITY = predsprovstoplot$DENSITY_RURAL_PROV # use rural densities as densities
predsprovstoplot$URBAN_RURAL = "RURAL"
predsprovstoplot$SOCIAL_CLASS = "Farmers"
preds = predict(fit_1$fitupd, newdata=predsprovstoplot, type="response") # predicted EPP prob for rural areas per province for Farmers
colnames(predsprovstoplot)[2] = "id"
predsprovstoplot$EPP_PRED = preds
predsprovstoplot = predsprovstoplot[,c("YEAR","id","EPP_PRED")]
provdata = fortify(admin_prov_sel)
provdata = merge(predsprovstoplot, provdata, by="id")
provdata$YEAR = factor(provdata$YEAR)
range(predsprovstoplot$EPP_PRED) # 0.00716992 0.01005535
mean(predsprovstoplot$EPP_PRED) # 0.008557556=0.85% = avg prediction for farmers in rural areas

# make predictions for lower socioeconomic classes in larger cities (urbancities)
# first we load our historical estimates of the densities and number of inhabitants 
# in larger urbanized cities through time (period 1750-1950) [for details and data sources see paper]: 
predscities = read.csv("data_dryad_densities_cities.csv")
predscities$URBAN_RURAL = "URBAN"
predscities$SOCIAL_CLASS = "Low income"
preds = predict(fit_1$fitupd, newdata=predscities, type="response") # predicted EPP prob per city
predscities$EPP_PRED_LOW = preds
predscities$YEAR = factor(predscities$YEAR)
range(predscities$EPP_PRED_LOW) # 0.02003370 0.05944355
max(predscities$EPP_PRED_LOW)/min(predsprovstoplot$EPP_PRED) # 8.15 fold diff between low class in city vs farmers on countryside
predscities$EPP_PRED = predscities$EPP_PRED_LOW

# make predictions for higher socioeconomic classes in larger cities (urbancities)
predscities_high = predscities
predscities_high$SOCIAL_CLASS = "Medium or high income"
predscities_high$EPP_PRED_HIGH = predict(fit_1$fitupd, newdata=predscities_high, type="response") # predicted EPP prob per city
predscities$EPP_PRED_HIGH = predscities_high$EPP_PRED_HIGH
range(predscities$EPP_PRED_HIGH) # 0.004659185 0.014264882

# translate Dutch city names to English where possible & add some space here and there to prevent overlapping labels
predscities$CITY=as.character(predscities$CITY)
predscities$CITY[predscities$CITY=="Antwerpen"] = "Antwerp"
predscities$CITY[predscities$CITY=="Oostende"] = "Ostend"
predscities$CITY[predscities$CITY=="Gent"] = "Ghent"
predscities$CITY[predscities$CITY=="Brussel"] = "Brussels"
predscities$CITY[predscities$CITY=="Den Haag"] = "The Hague"
predscities$CITY[predscities$CITY=="Brugge"] = "     Bruges"
predscities$CITY[predscities$CITY=="Leuven"] = "              Leuven"
predscities$CITY = as.factor(predscities$CITY)
# cities on map to label (cities with high EPPs and sufficient data) :
citiestolabel = c("Amsterdam","Antwerp","Boom","Breda","     Bruges", 
                  "Brussels","The Hague","Ghent","Kortrijk", 
                  "              Leuven","Maastricht","Ostend","Rotterdam","Tilburg") 

# sort by NRINHABITANTS to have largest cities being plotted on top
predscities = predscities[order(predscities$NRINHABITANTS),]

# add symbol for legend
for (i in 1:length(years)) {
  predscities = rbind(predscities, NA)
  r = nrow(predscities)
  predscities$YEAR[r] = years[i]
  predscities$LONGITUDE[r] = 2.46
  predscities$LATITUDE[r] = 52.14
  predscities$EPP_PRED_LOW[r] = 0.055
  predscities$EPP_PRED_HIGH[r] = 0.0095
  predscities$NRINHABITANTS[r] = 300000
}

# convert NRINHABITANTS into two new variables specifying the absolute point sizes
predscities$SIZECITIES = 0.32*(predscities$NRINHABITANTS)^(1/4) # for possible use with scale_size_identity()
predscities$SIZECITIES_INNER = predscities$SIZECITIES/3

maxepp <<- 0.075 # maximum of EPP scale (across all years plotted)

# now make and export plots per year
# the individual frames were combined into a video using Quicktime Pro

for (year in years) {
  
  print(year)
  
  provdata_subs <<- provdata[provdata$YEAR==year,]
  predscities_subs <<- predscities[predscities$YEAR==year,]
  
  print(qplot(data=provdata_subs, geom="blank") + 
          geom_rect(data=dfboxsea,aes(group=NULL,x=NULL,y=NULL,
                                      xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                    fill=I(seafill),colour=NA)+ # sea
          geom_rect(data=data.frame(xmin=2.3,xmax=3.9,ymin=51.89,ymax=53.44),aes(group=NULL,x=NULL,y=NULL, # white legend box
                                                                                 xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                    fill=I("white"),colour=NA)+ 
          geom_rect(data=data.frame(xmin=2.465-0.075,xmax=2.465+0.075,
                                    ymin=52.29-0.05-0.3,ymax=52.29+0.05-0.3),aes(group=NULL,x=NULL,y=NULL,  # for rural legend
                                                                                 xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                    fill=ruralfill,colour=NA) +
          geom_text(data=data.frame(LONGITUDE=c(2.61,2.61,2.61), LATITUDE=c(52.15,52.09,51.98), 
                                    CITY=c("urban (low class)", "urban (med-high class)", "rural (farmers)"), 
                                    colour=c(I("red"),I("green4"),I("green4"))), 
                    aes(x=LONGITUDE, y=LATITUDE, label=CITY), colour=c(I("red"),ruralfill,ruralfill),
                    size=I(3), hjust=I(0), vjust=I(0),
                    show.legend=FALSE) + # urban vs rural legend
          geom_polygon(data=admin_nl_l0,aes(x=long,y=lat,group=group), # COUNTRIES/LAND, BACKGROUND
                       fill=countryfill,colour=NA,lwd=lw)+ # NETHERLANDS
          geom_polygon(data=admin_bel_l0,aes(x=long,y=lat,group=group),
                       fill=countryfill,colour=NA,lwd=lw)+ # BELGIUM
          geom_polygon(data=admin_de_l0,aes(x=long,y=lat,group=group),
                       fill=countryfill,colour=NA,lwd=lw)+ # GERMANY
          geom_polygon(data=admin_lux_l0,aes(x=long,y=lat,group=group),
                       fill=countryfill,colour=NA,lwd=lw)+ # LUXEMBURG
          geom_polygon(data=admin_fra_l0,aes(x=long,y=lat,group=group),
                       fill=countryfill,colour=NA,lwd=lw)+ # FRANCE
          geom_polygon(data=admin_nl_l1_lakes,aes(x=long,y=lat,group=group),
                       fill=I(seafill),colour=I(seafill))+ # NETHERLANDS WATER, IJSSELMEER & FLEVOLAND
          geom_polygon(data=admin_bel_l0,aes(x=long,y=lat,group=group), # COUNTRIES THICK OUTLINES
                       fill=NA, colour=countryoutl,lwd=lw)+ # BELGIUM
          geom_polygon(data=admin_de_l0,aes(x=long,y=lat,group=group),
                       fill=NA, colour=countryoutl,lwd=lw)+ # GERMANY
          geom_polygon(data=admin_lux_l0,aes(x=long,y=lat,group=group),
                       fill=NA, colour=countryoutl,lwd=lw)+ # LUXEMBURG
          geom_polygon(data=admin_fra_l0,aes(x=long,y=lat,group=group),
                       fill=NA, colour=countryoutl,lwd=lw)+ # FRANCE
          geom_polygon(data=provdata_subs, aes(x=long, y=lat, group=group, fill=EPP_PRED), colour = "grey")+ # PLOT RURAL PREDICTIONS
          geom_point(data=predscities_subs, aes(x=LONGITUDE, y=LATITUDE, fill=EPP_PRED_LOW, 
                                                size=SIZECITIES, stroke=0), shape=I(21), colour=adjustcolor("black",alpha.f=0)) +
          # PLOT PREDICTIONS FOR URBAN LOW CLASS
          geom_point(data=predscities_subs, aes(x=LONGITUDE, y=LATITUDE, fill=EPP_PRED_HIGH,   
                                                size=SIZECITIES_INNER, stroke=0), shape=I(21), colour=adjustcolor("black",alpha.f=0)) +
          # PLOT PREDICTIONS FOR URBAN MED-HIGH CLASS
          geom_text(data=predscities_subs[(predscities_subs$CITY %in% citiestolabel),], 
                    aes(x=LONGITUDE, y=LATITUDE, label=CITY, size=SIZECITIES/2.2), 
                    colour=I("black"), hjust=I(0.5), vjust=I(-1.3), show.legend=FALSE) + # city labels 
          scale_size_identity(name="nr. of inhabitants",
                              breaks = c(0.32*(1000)^(1/4),
                                         0.32*(10000)^(1/4), 
                                         0.32*(100000)^(1/4), 
                                         0.32*(500000)^(1/4)),
                              labels = c("1,000", "10,000",  
                                         "100,000", "500,000"), limits=c(1,500000)) +
          scale_fill_gradientn(name="EPP frequency (%)",
                               colours = # legcols, 
                                 colorRampPalette(c(ruralfill, "yellow3", "red", "magenta3"),
                                                  bias=0.95, # put lower value for more green
                                                  space="Lab",interpolate="linear")(100), 
                               breaks=c(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10),
                               labels=c("0","1","2","3","4","5","6","7","8","9","10"),
                               limits=c(0,maxepp)) +
          scale_colour_gradientn(name="EPP frequency (%)",
                                 colours = # legcols,
                                   colorRampPalette(c(ruralfill, "yellow3", "red", "magenta3"),
                                                    bias=0.95, # put lower value for more green
                                                    space="Lab",interpolate="linear")(100),
                                 breaks=c(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10),
                                 labels=c("0","1","2","3","4","5","6","7","8","9","10"),
                                 limits=c(0,maxepp)) +
          xlab("Longitude")+ylab("Latitude")+
          coord_quickmap(xlim=longlims,ylim=latlims,expand=FALSE)+
          theme(panel.background = element_rect(fill="grey90"), plot.title = element_text(hjust = 0.5, size=12)) +
          theme(legend.position = "right") + theme(legend.direction = "vertical") + 
          guides(size = guide_legend(ncol = 1, title.position="top", 
                                     title.hjust = 0, keywidth=1, keyheight=1,
                                     override.aes = list(alpha = 1, colour=NA, fill=I("black")), order=1),
                 fill = guide_colourbar(title.position="top", 
                                        title.hjust = 0, barwidth=7, barheight=0.8, direction = "horizontal", order=2),
                 colour = FALSE) +
          theme(axis.title=element_text(size=18), axis.text=element_text(size=12)) + 
          theme(legend.key=element_blank(), legend.background=element_blank()) +
          theme(legend.justification = c(0, 1), legend.position = c(0.02, 0.98)) + 
          theme(legend.text=element_text(size=8), 
                legend.title=element_text(size=10)) +
          theme(legend.spacing.y = unit(0.09, "cm")) + 
          theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text=element_text(colour = "black"), 
            axis.title=element_text(colour = "black"),  
            axis.title.y = element_text(colour = "black"),
            axis.title.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.text.x = element_text(colour = "black") 
          ) + ggtitle(paste0("Estimated human extra-pair paternity in ",as.character(year)))) 
  
  
  graph2png(file=paste0("video_EPP_",as.character(year),".png"), dpi=1200, width=6.2, height=6.2, font="Arial")
  
  dev.off()
  gc()
  
}



#### 5. COMPARISON WITH MODEL FIT AT THE GENEALOGICAL PAIR LEVEL ####

# model fit in function of most common social class per genealogical pair, geometric mean density, 
# avg year of birth, most common country per pair, most common province per pair, sampling campaign

data_pair = data.frame(PAIR = unique(data_dryad$PAIR))
data_pair$SOCIAL_CLASS = factor(sapply(data_pair$PAIR, function (pair) names(which.max(table(data_dryad[(data_dryad$PAIR==pair)&(data_dryad$SOCIAL_CLASS!="missing"),"SOCIAL_CLASS"])))), levels=c("Farmers","Medium or high income","Low income")) # most common social class per pair
data_pair$GMEANDENSITY = sapply(data_pair$PAIR, function (pair) geometric.mean(data_dryad[data_dryad$PAIR==pair,"DENSITY"])) # geometric mean density
data_pair$YEAR = sapply(data_pair$PAIR, function (pair) mean(data_dryad[data_dryad$PAIR==pair,"YEAR"])) # mean year of birth
data_pair$COUNTRY = factor(sapply(data_pair$PAIR, function (pair) names(which.max(table(data_dryad[(data_dryad$PAIR==pair)&(data_dryad$COUNTRY!="outside_studyarea"),"COUNTRY"]))))) # most common country
data_pair$PROVINCE_MAP = factor(sapply(data_pair$PAIR, function (pair) names(which.max(table(data_dryad[(data_dryad$PAIR==pair)&(data_dryad$PROVINCE_MAP!="outside_studyarea"),"PROVINCE_MAP"]))))) # most common province
data_pair$SAMPLING_CAMPAIGN = factor(sapply(data_pair$PAIR, function (pair) names(which.max(table(data_dryad[data_dryad$PAIR==pair,"SAMPLING_CAMPAIGN"])))))
data_pair$MISMATCH = sapply(data_pair$PAIR, function (pair) max(data_dryad[data_dryad$PAIR==pair,"MISMATCH"])) # note that extra info on where EPP happened is not taken into account here
data_pair$NPERPAIR = sapply(data_pair$PAIR, function (pair) mean(data_dryad[data_dryad$PAIR==pair,"NPERPAIR"]))

fitpair = glm(cbind(MISMATCH,NPERPAIR-1) ~ SOCIAL_CLASS + GMEANDENSITY, family=binomial, data=data_pair)
summary(fitpair)
# Coefficients:
#                                         Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                           -4.8599965  0.1985727 -24.475  < 2e-16 ***
# SOCIAL_CLASSMedium or high income     -0.0949720  0.3738703  -0.254  0.79948    
# SOCIAL_CLASSLow income                 1.2863763  0.2392197   5.377 7.56e-08 ***
# GMEANDENSITY                           0.0002192  0.0000820   2.673  0.00752 ** 

# note that there was no overdispersion in this model:
fitpair_quasi = glm(cbind(MISMATCH,NPERPAIR-1) ~ SOCIAL_CLASS + GMEANDENSITY, family=quasibinomial, data=data_pair)
summary(fitpair_quasi) # no overdispersion - dispersion parameter = 0.81, ie close to 1
# Coefficients:
#                                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         -4.860e+00  1.790e-01 -27.148  < 2e-16 ***
# SOCIAL_CLASSMedium or high income   -9.497e-02  3.371e-01  -0.282  0.77824    
# SOCIAL_CLASSLow income               1.286e+00  2.157e-01   5.965 4.59e-09 ***
# GMEANDENSITY                         2.192e-04  7.393e-05   2.965  0.00317 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for quasibinomial family taken to be 0.8127681)



# PLOT MODEL PREDICTIONS
df = data_pair
preds = predict(fitpair, type="response")
df$preds = preds  
df = rbind(df[df$SOCIAL_CLASS=="Medium or high income",],
           df[df$SOCIAL_CLASS=="Farmers",],
           df[df$SOCIAL_CLASS=="Low income",]
)
df$SOCIAL_CLASS=factor(df$SOCIAL_CLASS, levels=c("Low income","Medium or high income","Farmers"))

plot1 = qplot(data=df[df$YEAR>1600,],x=YEAR,y=100*preds,
              xlab=expression("Average year of birth"),ylab="Extra-pair paternity frequency (%)",
              shape=I(21),stroke=I(0.5),
              fill=SOCIAL_CLASS,
              size=log10(GMEANDENSITY)) + mytheme +
  scale_fill_manual(name="PREDOMINANT SOCIAL CLASS",
                    labels=c("Low income","Medium or high income","Farmers"),
                    values=c(colrs[[1]],
                             colrs[[2]],
                             colrs[[3]])) +
  scale_radius(name=expression ("GEOM. MEAN DENSITY (per "*km^2*")"),
               breaks=c(4,3,2,1),
               labels=c("10,000","1,000","100","10"),range=c(0.5,6),trans="sqrt") + 
  coord_cartesian(ylim=c(0,6.2),xlim=c(1600,2020),expand=FALSE) +  
  guides(size = guide_legend(override.aes = list(alpha = 1, colour=I("black"), fill=I("black"), linetype = 0)),
         fill = guide_legend(override.aes = list(linetype = 0, size=5)))


plot1
graph2png(file="GLM per pair_model predictions in function of CLASS DENSITY by YEAR.png",width=11,height=6)
graph2ppt(file="GLM per pair_model predictions in function of CLASS DENSITY by YEAR.pptx",width=11,height=6)


# PLOT PARTIAL EFFECTS OF DIFFERENT FACTORS BASED ON EXPECTED MARGINAL MEANS

# partial effect PREDOMINANT SOCIOECONOMIC CLASS (with 95% confidence bands, ie one-sided 95% confidence intervals)
dfemmeansCLASS=data.frame(summary(emmeans(fitpair,~SOCIAL_CLASS),
                                  weights="proportional", type="response", level=0.9))
dfemmeansCLASS 
#              SOCIAL_CLASS        prob          SE  df   asymp.LCL  asymp.UCL
# 1               Farmers   0.008465831 0.001655515 Inf 0.006134844 0.01167209
# 2 Medium or high income   0.007704723 0.002427362 Inf 0.004584772 0.01292025
# 3            Low income   0.029978467 0.003882650 Inf 0.024210622 0.03706822
dfemmeansCLASS$SOCIAL_CLASS=factor(dfemmeansCLASS$SOCIAL_CLASS,
                                   levels=c("Farmers","Medium or high income","Low income"),
                                   labels=c("Farmer","Medium or high income","Low income"))
plot2 = qplot(data=dfemmeansCLASS, x=SOCIAL_CLASS,
              y=100*prob,
              ylab="Extra-pair paternity frequency (%)",
              xlab="Social class",fill=SOCIAL_CLASS,geom="col",width=0.75) +
  geom_linerange(aes(x=as.numeric(SOCIAL_CLASS),ymin=100*asymp.LCL,ymax=100*asymp.UCL),lwd=0.5) +
  scale_y_continuous(breaks=c(0:8), expand=c(0,0)) +
  scale_fill_manual(name="SOCIAL CLASS",
                    values=c(
                      colrs[[3]],
                      colrs[[2]],
                      colrs[[1]])) +
  scale_x_discrete(labels= c("Farmer","Med-High","Low")) +
  coord_cartesian(ylim=c(0,4)) +
  mytheme 
plot2
graph2png(file="GLM per pair_partial effect SOCIOECONOMIC CLASS.png", width=11, height=6)
graph2ppt(file="GLM per pair_partial effect SOCIOECONOMIC CLASS.pptx", width=11, height=6)

# partial effect of GEOMETRIC MEAN DENSITY (with 95% confidence bands, ie one-sided 95% confidence intervals)
dfemmeansDENS=data.frame(summary(emmeans(fitpair,~GMEANDENSITY,
                                         at=list(GMEANDENSITY=10^seq(0,4.03,length.out=1000)), # 3.975
                                         weights="proportional", type="response"), level=0.9)) # we show 95% confidence bounds here = one-sided 95% confidence limits
plot3 = qplot(data=dfemmeansDENS,x=sqrt(GMEANDENSITY),y=100*prob,ylab="Extra-pair paternity frequency (%)",
              xlab="Density (per km2)",col=I("red3"),geom="line",lwd=I(1), group=factor(df)) +
  geom_ribbon(aes(ymin=100*asymp.LCL,ymax=100*asymp.UCL),alpha=I(0.25),fill=I("red3"),colour=NA,lwd=0) + 
  geom_point(aes(colour=factor(df))) +
  scale_colour_manual(values=c(adjustcolor("white",alpha.f=0)), name=NULL, labels=NULL) +
  scale_x_continuous(breaks=sqrt(c(100,1000,10000)), labels=c("100","1,000","10,000")) +
  coord_cartesian(ylim=c(0,5), expand=FALSE) +  
  mytheme 
plot3
graph2png(file="GLM per pair_partial effect DENSITY.png",width=11,height=6)
graph2ppt(file="GLM per pair_partial effect DENSITY.pptx",width=11,height=6)

min(data_pair$GMEANDENSITY) # 45.60468
log10(min(data_pair$GMEANDENSITY)) # 1.66
max(data_pair$GMEANDENSITY) # 9448.159
log10(max(data_pair$GMEANDENSITY)) # 3.975
data.frame(summary(emmeans(fitpair, ~GMEANDENSITY,
                           at=list(GMEANDENSITY=c(min(data_pair$GMEANDENSITY), max(data_pair$GMEANDENSITY))),
                           weights="proportional", type="response"), level=0.9))
# predicted EPP for lowest & highest per-pair geom mean density: 1.1% & 8.0%
#   GMEANDENSITY       prob          SE  df   asymp.LCL  asymp.UCL
# 1     45.60468 0.01100318 0.001411428 Inf 0.008908095 0.01358426
# 2   9448.15932 0.08034280 0.055194381 Inf 0.024930715 0.22987997



# multipanel figure
p1 = ggplotGrob(plot1)
p2 = ggplotGrob(plot2+theme(axis.title.y=element_blank(), legend.position="none"))
p3 = ggplotGrob(plot3+theme(axis.title.y=element_blank(), legend.position="none"))
grid.arrange(p1, arrangeGrob(p2, p3, ncol=1, heights=c(1,1)), 
             ncol=2, widths=c(2,1))
graph2png(file="GLM per pair_multipanel figure.png",dpi=1200,width=12,height=6)
graph2ppt(file="GLM per pair_multipanel figure.pptx",width=12,height=6)
