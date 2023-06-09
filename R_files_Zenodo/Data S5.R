rm(list=ls())

library(divDyn)
library(nlme)
library(chronosphere)
library(scales) 
library(metafor)
#library(glmulti) #There is a bug clash between R 4.2 and installing this, which crashes the workspace
#use earlier version of R or a patched version
library(rcompanion)
library(visreg)
library(readxl)
library(plotrix)
library(car)

#Specify working directory
setwd("C:/Users/")

#Functions ----

gen.diff.corr<-function(x,z){
  f1<-cor.test(x[-1],x[-length(x)]) #First order correlation
  new.f1<-x[-1]-f1$estimate*x[-length(x)] #Generalised differences temp
  f2<-cor.test(z[-1],z[-length(z)]) #First order correlation
  new.f2<-z[-1]-f2$estimate*z[-length(z)] #Generalised differences temp
  out<-cor.test(new.f1,new.f2,method="spearman")
  return(out)
}

#Function from http://stackoverflow.com/questions/743812/calculating-moving-average-in-r
ma.pres <- function(x,n=3){filter(x,rep(1/n,n), sides=2)}  

#Function to make probabilistic extinctions
make.probdat<-function(dat){
  
  subsdat2<-NA
  probdat<-slc.dat[1:2,]
  
  for (i in 43:94){ #43 from Permian to Pleistocene
    subsdat<-subset(dat,dat$slc==i) 
    
    #remove the line below for range-through (RT), as used in the paper. Add it for without RT
    #subsdat<-subset(subsdat,subsdat$occurrence==T)
    
    #This uses the sampling completeness (age and taxon) as the probability of true extinction.   REMOVE next four lines and end two to remove this step
    for (j in 1:length(taxnams)){
      
      subsdat$ext[subsdat$taxon==taxnamsall[j,1]]<-subsdat$ext[subsdat$taxon==taxnamsall[j,1]] * ((sam3t[i+1] + as.numeric(taxnamsall[j,2])+0.138762)/2)}  #+0.138762 is the difference between the means of tax and time
    
    subsdat3<-subsdat #Store this as just the extinctions from the current bin, else the below line will keep accumulating the taxa
    subsdat<-rbind(subsdat,subsdat2) #add subsdat2, from previous bin, to subsdat
    probdat<-rbind(probdat,subsdat) #store new extinction probability data for plotting etc
    
    #Add for sampling completeness
    subsdat2<-subsdat3[subsdat3$ext>0,] #The remainder of the probablility of LAD (via sampling completeness) is then forced into the next bin (unless a pulsed event)
    subsdat2$ext<- 1-subsdat2$ext
    subsdat2$slc<-subsdat2$slc+1 #update the slc for these remainder extinctions
    
    print(i)
  } 
  
  probdat[1:2,]<-NA
  
  return(probdat)
}

#Big function to run logistic regressions for all stages for whichever data (e.g. warm-water split from cold)
allregress<-function(dat,prob,therm.all=FALSE){
  #Vectors for McFadden's pseudo-R squared
  therm.sig.all<-matrix(NA,nrow=2,ncol=4)
  modloglik<-numeric(length=95)
  
  #Update the following with random effects
  nulldf<-numeric(length=95);b<-Sys.time()
  subsdat2<-NA
  
  for (i in 43:94){ #43 from Permian to Pleistocene
    subsdat<-subset(dat,dat$slc==i) 
    
    #remove the line below for RT. Add it for without RT
    #subsdat<-subset(subsdat,subsdat$occurrence==T)
    
    if(prob==TRUE){
      #This uses the sampling completeness (age and taxon) as the probability of true extinction.   REMOVE next four lines and end two to remove this step
      for (j in 1:length(taxnams)){
        
        subsdat$ext[subsdat$taxon==taxnamsall[j,1]]<-subsdat$ext[subsdat$taxon==taxnamsall[j,1]] * ((sam3t[i+1] + as.numeric(taxnamsall[j,2])+0.138762)/2)}  #+0.138762 is the difference between the means of tax and time
      
      subsdat3<-subsdat #Store this as just the extinctions from the current bin, else the below line will keep accumulating the taxa
      subsdat<-rbind(subsdat,subsdat2)} #add subsdat2, from previous bin, to subsdat
    
    if(length(table(subsdat$ext[is.na(subsdat$thermN)==F]))<=1) print(c(i,"too few occurrences"))
    
    if(length(table(subsdat$ext[is.na(subsdat$thermN)==F]))>1){
      
      nullmod <- glm(ext~1, family="binomial",data=subsdat[is.na(subsdat$ext)==F,])
      
      if(therm.all==FALSE){
        m5<-glm(ext~thermN,family="binomial",data=subsdat[!is.na(subsdat$ext)& !is.na(subsdat$thermN),])
      }
      if(therm.all==TRUE){
        #Polynomial of 3 is probably the best compromise (given the smaller age-scale sample sizes)
        m5<-glm(ext~poly(thermN,3),family="binomial",data=subsdat[!is.na(subsdat$ext)& !is.na(subsdat$thermN),])
      }
      
      modloglik[i]<-1-logLik(m5)/logLik(nullmod)
      
      tax.sig<-cbind(summary(m5)$coefficients[,c(1,2,4)],i) 
      therm.sig.all<-rbind(therm.sig.all,tax.sig)}#} 
    
    if(prob==TRUE){
      #Add for sampling completeness
      subsdat2<-subsdat3[subsdat3$ext>0,] #The remainder of the probablility of LAD (via sampling completeness) is then pushed into the next bin
      subsdat2$ext<- 1-subsdat2$ext 
      subsdat2$slc<-subsdat2$slc+1 #for completeness, update the slc for these remainder extinctions
    }
    
    print(i)
  } ;Sys.time()-b
  
  therm.sig.all<-therm.sig.all[-1:-2,]
  
  #To add stage names
  st.names.therm<-numeric()
  for (i in 43:94){
    st.names.therm[therm.sig.all[,4]==i]<-paste(stages[i,"stage"])}
  
  therm.sig.all<-cbind(therm.sig.all,st.names.therm)
  modloglik<-cbind(modloglik,stages[,4])
  
  #add the environ vars. glsdatT stages 18:91
  therm.sig.all<-cbind(row.names(therm.sig.all),therm.sig.all)
  therm.sig.all2<-as.data.frame(therm.sig.all)
  therm.sig.all2<-cbind(therm.sig.all2,
                        glsdatT$temperature[as.numeric(therm.sig.all[,5])-17],
                        glsdatT$tropshel[as.numeric(therm.sig.all[,5])-17],
                        glsdatT$tempshel[as.numeric(therm.sig.all[,5])-17],
                        glsdatT$troptempshel[as.numeric(therm.sig.all[,5])-17],
                        glsdatT$Occur.med[as.numeric(therm.sig.all[,5])-17],
                        glsdatT$slc.latran[as.numeric(therm.sig.all[,5])-17],
                        glsdatT$nrec[as.numeric(therm.sig.all[,5])-17],
                        glsdatT$Peterssili[as.numeric(therm.sig.all[,5])-17], #note that the carb is only really before the Triassic
                        glsdatT$Lat_max[as.numeric(therm.sig.all[,5])-17],
                        glsdatT$Lat_min[as.numeric(therm.sig.all[,5])-17])
  row.names(therm.sig.all2)<-paste(1:length(therm.sig.all2[,1]))
  names(therm.sig.all2)<-c("TermName","Estimate","SE","P","Slc","Stage","temperature","tropshel","tempshel","troptempshel","Occur.med","latran","nrec","Peterssili","Lat_max","Lat_min")
  return(therm.sig.all2)
}

#function for model selection in meta-analysis
rma.glmulti <- function(formula, data, ...)   rma(formula, sei=as.numeric(moddat$SE)^2, data=moddat, method="REML", ...)

#Inverse logit function and plotting from https://www.polyu.edu.hk/cbs/sjpolit/logisticregression.html
inverselogit <- function(x){ exp(x) / (1+exp(x) ) }



#Initial handling of PBDB data download ----

spdat <- read.csv(file="Data S1. pbdb_data.Carbonif_Holocene.csv", header=TRUE)

# load corresponding stages
data(stages) #from divDyn

##Stratigraphic binning ----
#This following the standard steps in the 'ddPhanero' vignette accompanying the divDyn package
#https://doi.org/10.5281/zenodo.2545983

data(keys)
# the 'stg' entries (lookup)
stgMin <- categorize(spdat[ ,"early_interval"], keys$stgInt)
stgMax <- categorize(spdat[ ,"late_interval"], keys$stgInt)
# convert to numeric
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)
# empty container
spdat$slc <- rep(NA, nrow(spdat))
# select entries, where
stgCondition <- c(
  # the early and late interval fields indicate the same stg
  which(stgMax==stgMin),
  # or the late_intervar field is empty
  which(stgMax==-1))
# in these entries, use the stg indicated by the early_interval
spdat$slc[stgCondition] <- stgMin[stgCondition]

sum(is.na(spdat$slc)==F)/length(spdat$slc)

#Filling in NAs with given binning column
for (i in stages$stg){
  spdat$slc[is.na(spdat$slc) & spdat$time_bins==stages$stage[i]]<-i }
sum(is.na(spdat$slc)==F)/length(spdat$slc) #proportion saved

# Exclude occurrences with bins unassigned
spdat <- subset(spdat, is.na(slc)==F)

# filter out terrestrial and freshwater environments. The following environments are expected to contain terrestrial taxa.
omitEnv <- c(
  "\\"floodplain\\"", "alluvial fan", "cave", "\\"channel\\"", "channel lag" ,
  "coarse channel fill", "crater lake", "crevasse splay", "dry floodplain",
  "delta plain", "dune", "eolian indet.", "fine channel fill", "fissure fill",
  "fluvial indet.", "fluvial-lacustrine indet.", "fluvial-deltaic indet.",
  "glacial", "interdune", "karst indet.", "lacustrine - large",
  "lacustrine - small", "lacustrine delta front", "lacustrine delta plain",
  "lacustrine deltaic indet.", "lacustrine indet.",
  "lacustrine interdistributary bay", "lacustrine prodelta", "levee", "loess",
  "mire/swamp", "pond", "sinkhole", "spring", "tar", "terrestrial indet.",
  "wet floodplain")
# omit the occurrences
spdat <- spdat[!spdat$environment%in%omitEnv, ]

#############Sorting out taxonomy ----

#Remove traces (keep all kinds of body fossils)
tracelist<-unique(spdat$pres_mode[grep("trace",spdat$pres_mode)])
spdat<-spdat[(spdat$pres_mode %in% unique(tracelist[-grep("body",tracelist)]))==FALSE,]

#marine invert clades and traits binned as in Reddin et al. 2021, DOI: 10.1111/gcb.15434
spdat$class[spdat$class%in%c("Sclerospongiae","Hexactinellida")]<-paste("Demospongea") 
spdat<-spdat[!spdat$genus%in%c(""," "),]#remove unlabelled genera
spdat<-spdat[!spdat$order%in%c("Unionida"),] #Remove near completely freshwater Unionids. Checked and several of these genera are freshwater only

ph.n <- c("Brachiopoda","Porifera","Echinodermata","Bryozoa")
cl.n <- c("Bivalvia", "Gastropoda","Anthozoa","Trilobita","Scaphopoda","Ostracoda","Xiphosura","Polychaeta","Polyplacophora","Rostroconchia","Tergomya","Malacostraca")
spdat <- subset(spdat, class %in% cl.n | phylum %in% ph.n)

#remove nektonic and planktonic
spdat<- spdat[!spdat$life_habit%in%c("nektobenthic", "nektonic", "planktonic" ) | spdat$class=="Bivalvia",] #remove nektonic or planktonic. keep the scallops
spdat<- spdat[!spdat$motility%in%c("passively mobile", "passively mobile, epibiont") ,] #removes the borers of driftwood and a couple of odd groups

#These artificial names or groupings ensure the smaller classes are not completely skipped by the regressions by being given NA here. Check original NAs by only running the first line below
spdat$taxon<-spdat$class
spdat$taxon[spdat$phylum=="Porifera"]<-paste("Porifera")
spdat$taxon[spdat$phylum=="Bryozoa"]<-paste("Bryozoa")
spdat$taxon[spdat$class%in%c("Edrioasteroidea","Ophiocistioidea","Stylophora")]<-paste("earlyEchino") 
spdat$taxon[spdat$class=="NO_CLASS_SPECIFIED" & spdat$phylum=="Echinodermata"]<-paste("earlyEchino") 
spdat$taxon[spdat$class=="Paterinata"]<-paste("Lingulata")
spdat$taxon[spdat$class=="Tergomya"]<-paste("Polyplacophora")
spdat$taxon[spdat$class=="Xiphosura"]<-paste("Trilobita") #similar sampling completeness expected
spdat$taxon[spdat$class%in%c("Craniata","Chileata")]<-paste("smallBrachio")
spdat$taxon[spdat$class=="NO_CLASS_SPECIFIED" & spdat$phylum=="Brachiopoda"]<-paste("smallBrachio") 

#### Get subgenus names and raise to genera
new.gen <- as.character(spdat$genus)
bsp <- strsplit(as.character(new.gen), " ")
sg <- sapply(bsp, "[", 1)
sg2 <- sapply(bsp, "[", 2)
sg2 <- gsub("\\\\(|\\\\)", "", sg2)
sg[!is.na(sg2)] <- sg2[!is.na(sg2)]
spdat$gen <- sg

#Save the clade sampling completeness
divstats<-divDyn(spdat,"clgen","slc")
sam3t<-divstats$samp3t
taxnams<-unique(spdat$taxon)
taxnamsall<-matrix(NA,ncol=2,nrow = length(taxnams))
for (i in 1:length(taxnams)){
  if(length(unique(spdat[spdat[,"taxon"]==taxnams[i],"slc"]))>1){ #does this taxon have occurrences in more than 2 stages? Otherwise leave 'NA'
    taxsam3t<-divDyn(spdat[spdat$taxon==taxnams[i],],"gen","slc")$samp3t
    taxnamsall[i,]<-c(taxnams[i],median(taxsam3t[!taxsam3t%in%c(1,0)],na.rm=T))# remove 1s, 0s and NAs for averaging
  }else{taxnamsall[i,1]<-taxnams[i]}}
taxnamsall[order(as.numeric(taxnamsall[,2])),]
#Sampling completeness == 0.5 are 50:50 whether their occurrence gets picked up at all
taxnamsall[as.numeric(taxnamsall[,2])<0.55 & is.na(taxnamsall[,2])==F,1] #But results stay the same when the worst sampled groups are removed
#"Porifera"       "Crinoidea"      "Rostroconchia"  "smallBrachio"   "Ophiuroidea"    "Asteroidea"     "Holothuroidea"  "Polyplacophora"

spdat$clgen <- paste(spdat$class, spdat$gen) #Makes clgen for next step

#Clearing up memory of unnecessary objects
rm(list=c("stgMin","stgMax","stgCondition","sg","sg2","new.gen","bsp"))

#Removing singletons (note that this should be done after sampling completeness is calculated)
x <- table(factor(spdat$clgen))
spdat <- subset(spdat, !spdat$clgen %in% names(x[x==1]))

#correcting a bad datum (not fixed in the PBDB yet)
spdat<-spdat[!(spdat$clgen%in% c("Anthozoa Paleofavosites") &spdat$slc>51),] #Triassic tabulate

# Delete duplicate genus occurrences within collections
t1 <- subset(spdat, select=c(collection_no, clgen))
spdat <- spdat[!duplicated(t1),]

#Merging the Hettangian and then splitting resultant Sinemurian
#First, bin data to the nearest 1myr
spdat$mean_ma<-(spdat$max_ma+spdat$min_ma)/2
spdat$trunc_ma<-trunc(spdat$mean_ma)
spdat$slc[spdat$slc==59]<-60
spdat$slc[spdat$slc==60 & spdat$trunc_ma<193.5]<-60.5
spdat$slc[spdat$slc==60]<-59
spdat$slc[spdat$slc==60.5]<-60
spdat<-spdat[,!names(spdat)%in%c("trunc_ma","mean_ma")]

data(keys)
# the environmental and lithological information can be categorized

# bathymetry
spdat$depenv <- categorize(spdat$environment,keys$bath)
spdat$depenv[spdat$environment%in%c("deep subtidal indet.","deep subtidal ramp","deep subtidal shelf","coastal indet.")] <- "mid" #Between fair weather and storm wave base 
# grain size
spdat$gra <- categorize(spdat$lithology1,keys$grain)

# reef or not?
spdat$reef <- categorize(spdat$environment, keys$reef)
spdat$reef[spdat$lith=="clastic" & spdat$environment=="marine indet."] <- "non-reef"

#Primary lithology
spdat$lith<-NA
spdat$lith[spdat$lithology1 %in% c("limestone" ,"carbonate" ,"packstone","wackestone","dolomite","lime mudstone",
                                     "grainstone","reef rocks" ,"framestone","bafflestone","bindstone","floatstone","rudstone" )]<-"Carbonate"
spdat$lith[spdat$lithology1 %in% c("shale","siltstone" ,"claystone","mudstone","conglomerate", "siliciclastic","phyllite",
                                     "sandstone","slate", "schist"  ,"quartzite")]<-"Clastic"
spdat$lith[spdat$lithology1 %in% c("marl","mixed carbonate-siliciclastic" )]<-"Mixed"



#Loading in the final data ---- 

#Read in the finished occurrence data
spdat1<-read.csv("Data S2. Reddin et al. final occurrence data.csv")

#min number of unique locations? the clgen that have only been sampled in one temperature cell through time may be useless for this analysis
tRans<-tapply(spdat1$SST,spdat1$clgen, range)
tRans<-unlist(lapply(tRans,diff))
spdat1<-spdat1[spdat1$clgen %in% names(tRans[tRans>0]),]

#Removing singletons (note that these should already be removed but the above steps may have recreated some)
x <- table(factor(spdat1$clgen)); length(x[x==1])
spdat1 <- subset(spdat1, !spdat1$clgen %in% names(x[x==1]))
rm(x)

#Read in the finished environmental data
glsdatT<-read.csv("Data S3. Finished env data.csv")

#Read in the additional environmental data
Extra.Env.dat<-read.csv("Data S4. Additional env data.csv")
Extra.Env.dat<-Extra.Env.dat[,-1]

## Hyperthermals and warm events ----

#stages including other warming events listed by Scotese et al. (2021) 
Hyperther.Scot<-c(51, 53, 56, 58, 62, 65, 67, 72, 74, 75, 76, 78, 83, 84, 85, 86, 90) 
hyperther<-c(51,56,58,61,74,76,83) #Bin numbers of the hyperthermal onsets: PT, Carn, TJ, PLiToa, Apt-Alb, Cen-Tur, PETM
excluhyperther<-c(52,53,62) #Bin numbers for exclusion where the hyperthermal continues


#Calculate extinction timings ----

#decide thermal mean per taxon per age or through time?
tOpts<-tapply(spdat1$SST,spdat1$clgen, median)
tRans<-tapply(spdat1$SST,spdat1$clgen, range)
tRans<-unlist(lapply(tRans,diff))
tQuRans<-tapply(spdat1$SST,spdat1$clgen, quantile,probs=c(0.25,0.75))
tQuRans<-unlist(lapply(tQuRans,diff))
names(tQuRans)<-gsub(".75%","",names(tQuRans),fixed=T)
tQuRans2<-tapply(spdat1$SST,spdat1$clgen, quantile,probs=c(0.05,0.95))
tMins<-unlist(lapply(tQuRans2,min))
tMaxs<-unlist(lapply(tQuRans2,max))
rm("tQuRans2")

#Makes range-through occurrence time series with trait entries for each taxon
slc.dat<-modeltab(spdat1, "clgen", "slc", taxvars = c("taxon","phylum","class","genus","order","composition","motility","diet","life_habit"), rt = TRUE, singletons = TRUE)

#remove pre-Carboniferous
spdat1<-spdat1[spdat1$slc>35,]
slc.dat<-slc.dat[slc.dat$slc>35,] 

#note that the RTs will obviously have no slc-based therm niche estimates. modeltab() copies their value from the first occurrence.
#But better idea is to give them their all-time optimum below

#Fill in the T optima, interquartile ranges, upper and lower tolerances (based on all records, assuming static thermal niche). Takes 4 minutes for me
b<-Sys.time()
slc.dat$thermN<-NA; slc.dat$thermRan<-NA; slc.dat$thermMin<-NA; slc.dat$thermMax<-NA
for (i in 1:length(slc.dat[,1])){
  slc.dat$thermN[slc.dat$clgen==names(tOpts[i])]<-tOpts[i]
  slc.dat$thermRan[slc.dat$clgen==names(tQuRans[i])]<-tQuRans[i] #using interquantile range rather than 'tRans'
  slc.dat$thermMin[slc.dat$clgen==names(tMins[i])]<-tMins[i] #0.1 and 0.9 quantile
  slc.dat$thermMax[slc.dat$clgen==names(tMaxs[i])]<-tMaxs[i]} 
Sys.time() - b

#Again remove those with the least information. should be none
slc.dat<-slc.dat[slc.dat$clgen %in% names(tRans[tRans>0]),]

#Make the extinctions probablistic (non-integer)
probdat<-make.probdat(slc.dat)


#Analysis. Checking which polynomial is optimal ----

#Non-hyperthermal data pooled
m5<-glm(ext~thermN,family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN) & !probdat$slc %in%c(hyperther,excluhyperther),])  #}
m10<-glm(ext~poly(thermN,2),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN) & !probdat$slc %in%c(hyperther,excluhyperther),])  #}
m15<-glm(ext~poly(thermN,3),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& !probdat$slc %in%c(hyperther,excluhyperther),])  #}
m20<-glm(ext~poly(thermN,4),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& !probdat$slc %in%c(hyperther,excluhyperther),])  #}
m25<-glm(ext~poly(thermN,5),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& !probdat$slc %in%c(hyperther,excluhyperther),])  #}
m30<-glm(ext~poly(thermN,6),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& !probdat$slc %in%c(hyperther,excluhyperther),])  #}
m35<-glm(ext~poly(thermN,7),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& !probdat$slc %in%c(hyperther,excluhyperther),])  #}
m40<-glm(ext~poly(thermN,8),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& !probdat$slc %in%c(hyperther,excluhyperther),])  #}
BIC(m5,m10,m15,m20,m25,m30,m35,m40) #15
newdat = data.frame(thermN = seq(min(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN),]$thermN), max(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN),]$thermN), length.out = 100))
newdat$nonhyppred = predict(m15, newdata = newdat,type="response")

#Hyperthermal data pooled
m5<-glm(ext~thermN,family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
m10<-glm(ext~poly(thermN,2),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN) & probdat$slc %in%hyperther,])  #}
m15<-glm(ext~poly(thermN,3),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
m20<-glm(ext~poly(thermN,4),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
m25<-glm(ext~poly(thermN,5),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
m30<-glm(ext~poly(thermN,6),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
m35<-glm(ext~poly(thermN,7),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
m40<-glm(ext~poly(thermN,8),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
m45<-glm(ext~poly(thermN,9),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
m50<-glm(ext~poly(thermN,10),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
BIC(m5,m10,m15,m20,m25,m30,m35,m40,m45,m50) #10
newdat$hyppred = predict(m10, newdata = newdat,type="response")

#To plot these for checking
plot(ext ~ thermN, data = slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermN)& slc.dat$slc %in%hyperther,],xlab="Genus thermal optimum, degrees Celsius",ylab="Extinction odds",type="n",main="(H) Multiple stages pooled, ",ylim=c(0,1))
with(newdat, lines(x = thermN, y = nonhyppred,lwd=2,col="blue"))
with(newdat, lines(x = thermN, y = hyppred,lwd=2,col="red"))

quarts<-quantile(tapply(spdat1$SST,spdat1$clgen, median),probs=c(0.1,0.9))
abline(v=quarts,lty="dashed")


###Fig. 1 ----
graphics.off()
png("Fig.1.Final.png",4,4,units = "in",res=300)

m15<-glm(ext~poly(thermN,3),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& !probdat$slc %in%c(hyperther,excluhyperther),])  #}
newdat = data.frame(thermN = seq(min(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN),]$thermN), max(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN),]$thermN), length.out = 100))
newdat$nonhyppred = predict(m15, newdata = newdat,type="response")
x<-predict(m15, newdata = newdat,se.fit=T)
newdat$nonhyppred2 = inverselogit(x[[1]]) #check this agrees with the response
newdat$nonhypu95 = inverselogit(x[[1]] + 1.96 * x[[2]])
newdat$nonhypl95 = inverselogit(x[[1]] - 1.96 * x[[2]])

m10<-glm(ext~poly(thermN,2),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN) & probdat$slc %in%hyperther,])  #}
#m20<-glm(ext~poly(thermN,3),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
#m30<-glm(ext~poly(thermN,4),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,])  #}
BIC(m10,m20,m30)

newdat$hyppred = predict(m10, newdata = newdat,type="response")
x<-predict(m10, newdata = newdat,se.fit=T)
newdat$hyppred2 = inverselogit(x[[1]]) #check this agrees with the response
newdat$hypu95 = inverselogit(x[[1]] + 1.96 * x[[2]])
newdat$hypl95 = inverselogit(x[[1]] - 1.96 * x[[2]])

length(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,1])
plot(ext ~ thermN, data = probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,],
     ylim=c(0,1.02),col=alpha("red",0.2),ylab="Extinction odds",xlab="Genus thermal optimum, degrees Celsius",pch=20)
#title(main=paste("Hyperthermal data pooled, n = ",paste(length(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,"clgen"]))),cex.main=0.9)
polygon(c(rev(newdat$thermN), newdat$thermN), c(rev(newdat$hypu95), newdat$hypl95), col = alpha("red",0.3), border = NA)
with(newdat, lines(x = thermN, y = nonhyppred,lwd=2,col="blue"))
with(newdat, lines(x = thermN, y = hyppred,lwd=2,col="red"))

quarts<-quantile(tapply(spdat1$SST,spdat1$clgen, median),probs=c(0.1,0.9))
abline(v=quarts,lty="dashed")

x<-hist(probdat$thermN,plot =F) #,main="(B)",xlab="Thermal optimum, degrees C" 
x$counts<-x$counts/20000 # 100000
plot(x,add=T,col=alpha("grey",0.5))
#abline(v=quarts,lty="dashed")
with(newdat, lines(x = thermN, y = nonhyppred,lwd=2,col="blue"))
with(newdat, lines(x = thermN, y = hyppred,lwd=2,col="red"))

dev.off()

##Zooming in to check check curve intersections on the last plot
plot(ext ~ thermN, data = probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,],xlab="",ylab="",ylim=c(0.08,0.16),col="red",type="n",xlim=c(7.5,16))
plot(ext ~ thermN, data = probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc %in%hyperther,],xlab="",ylab="",ylim=c(0.05,0.1),col="red",type="n",xlim=c(19,24))
with(newdat, lines(x = thermN, y = nonhyppred,lwd=2,col="blue"))
with(newdat, lines(x = thermN, y = hyppred,lwd=2,col="red"))

with(newdat, lines(x = thermN, y = nonhypu95,col="blue"))
with(newdat, lines(x = thermN, y = hypu95,col="red"))
with(newdat, lines(x = thermN, y = nonhypl95,col="blue"))
with(newdat, lines(x = thermN, y = hypl95,col="red"))

#Use locator() to check curve intersections
locator()


##Fig. 2. Hyperthermal log plots from GLMs ----
png("AllHyperthermalResponsesUpdate.png",height=9, width = 5,units = "in",res=300)

#Divide the screen in 1 line and 2 columns
par(  mfrow=c(4,2),   oma = c(2, 2, 0, 0)) 

#Make the margin around each graph a bit smaller
par(mar=c(2,2,2,2))

samran<-tapply(spdat1$SST[spdat1$slc%in%hyperther],spdat1$slc[spdat1$slc%in%hyperther],quantile,probs=c(0.10,0.90))
j=1 #;layout(matrix(1:8,4,2,byrow=T))
for (i in hyperther){
  m5<-glm(ext~poly(thermN,3),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc==i,])  #}
  newdat = data.frame(thermN = seq(min(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc==i,]$thermN), max(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc==i,]$thermN), length.out = 100))
  newdat$pred = predict(m5, newdata = newdat,type="response")
  x<-predict(m5, newdata = newdat,se.fit=T)
  newdat$pred2 = inverselogit(x[[1]]) #check this agrees with the response
  newdat$u95 = inverselogit(x[[1]] + 1.96 * x[[2]])
  newdat$l95 = inverselogit(x[[1]] - 1.96 * x[[2]])
  plot(ext ~ thermN, data = probdat[!is.na(probdat$ext)& !is.na(probdat$thermN) & probdat$slc==i,],xlab="",ylab="",col="red",cex=0.8,ylim=c(0,1.02))
  title(main=paste(c("(a)","(b)","(c)","(d)","(e)","(f)","(g)")[j]," ",stages[i,"stage"],", n =",length(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc==i,1])),cex.main=0.9)
  polygon(c(rev(newdat$thermN), newdat$thermN), c(rev(newdat$u95), newdat$l95), col = alpha("red",0.3), border = NA)
  with(newdat, lines(x = thermN, y = pred,lwd=2,col="red"))
  print(nagelkerke(m5))
  abline(v=samran[[j]],lty="dashed")
  j<-j+1
  
}

mtext("Extinction odds", outer = TRUE, cex = 1, side = 2 )
mtext("Genus thermal optimum, degrees Celsius", outer = TRUE,line=0.5, cex = 1, side = 1 )

dev.off()

##Fig. S3 ----

png("allstg.png",width=20,height=13,units="in",res=300)
par(mar=c(2.5,3,2,1))
layout(matrix(1:54,6,9,byrow=T))
for (i in 43:94){
  m5<-glm(ext~poly(thermN,3),family="binomial",data=probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc==i,])  #}
  newdat = data.frame(thermN = seq(min(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc==i,]$thermN), max(probdat[!is.na(probdat$ext)& !is.na(probdat$thermN)& probdat$slc==i,]$thermN), length.out = 100))
  newdat$pred = predict(m5, newdata = newdat,type="response")
  x<-predict(m5, newdata = newdat,se.fit=T)
  newdat$pred2 = inverselogit(x[[1]]) #check this agrees with the response
  newdat$u95 = inverselogit(x[[1]] + 1.96 * x[[2]])
  newdat$l95 = inverselogit(x[[1]] - 1.96 * x[[2]])
  
  plot(ext ~ thermN, data = probdat[!is.na(probdat$ext)& !is.na(probdat$thermN) & probdat$slc==i,],xlab="",ylab="",ylim=c(0,1),axes=F)
  polygon(c(rev(newdat$thermN), newdat$thermN), c(rev(newdat$u95), newdat$l95), col = alpha("grey",0.3), border = NA)
  axis(1,at=seq(0,40,10),cex.axis=1.5); axis(2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=1.5,las=1); title(main=stages$stage[i],cex.main=2) #short
  with(newdat, lines(x = thermN, y = pred,lwd=2))
  if(i %in% Hyperther.Scot[!Hyperther.Scot%in%hyperther]){
    with(newdat, lines(x = thermN, y = pred,lwd=2,col="#ca2c92"))
    polygon(c(rev(newdat$thermN), newdat$thermN), c(rev(newdat$u95), newdat$l95), col = alpha("pink",0.3), border = NA)
  }
  if(i %in% hyperther){
    with(newdat, lines(x = thermN, y = pred,lwd=2,col="red"))
    polygon(c(rev(newdat$thermN), newdat$thermN), c(rev(newdat$u95), newdat$l95), col = alpha("red",0.3), border = NA)
  }
}
dev.off()


#Analysis. Meta-regression split for hot and cold slopes ----

#split hot and cold based on a mathematically defined threshold
quantile(tapply(spdat1$SST,spdat1$clgen, median),probs=c(0.1,0.9))
quarts<-c(15,30)
#At the above cut-offs, linear relationships are optimal
slc.cold<-slc.dat[slc.dat$thermN<=quarts[1],] #these are mutually exclusive, 19.5 marks the 'cold' as the lower quartile
hist(slc.cold$thermN); coldoffset<-16  
hist(log(coldoffset-slc.cold$thermN)) 
slc.cold$thermN<-log(coldoffset-slc.cold$thermN) #note this transformation reverses the sign

slc.hot<-slc.dat[slc.dat$thermN>=quarts[2],]
hist(slc.hot$thermN) ; hotoffset<-29 
hist(log(slc.hot$thermN-hotoffset))
slc.hot$thermN<-log(slc.hot$thermN-hotoffset)

##Fig. S1. Histogram figure ----
layout(matrix(1:4,2,2))
hist(slc.dat[slc.dat$thermN<=quarts[1],]$thermN,main="(A) Cold genera",xlab="Temperature, degrees C")#;log(21-slc.cold$thermN)
hist(slc.cold$thermN,main="(C) Cold genera (transformed)",xlab="Temperature, degrees C",axes=F,xlim=log(coldoffset-c(-2,15)));axis(2)
axis(1,at=log(coldoffset-c(-2,5,10,13,15)),labels=c(-2,5,10,13,15))
hist(slc.dat[slc.dat$thermN>=quarts[2],]$thermN,main="(B) Warm genera",xlab="Temperature, degrees C")#;log(21-slc.cold$thermN)
hist(slc.hot$thermN,main="(D) Warm genera (transformed)",xlab="Temperature, degrees C",axes=F);axis(2)
axis(1,at=log(c(30,31,32,34)-hotoffset),labels=c(30,31,32,34))


cold.prob<-make.probdat(slc.cold)
hot.prob<-make.probdat(slc.hot)
all.prob<-make.probdat(slc.dat)

#Removing the poorly sampled clades?
#slc.cold2<-slc.cold[!slc.cold$taxon %in% taxnamsall[as.numeric(taxnamsall[,2])<0.55 & is.na(taxnamsall[,2])==F,1],]
#slc.hot2<-slc.hot[!slc.hot$taxon %in% taxnamsall[as.numeric(taxnamsall[,2])<0.55 & is.na(taxnamsall[,2])==F,1],]

#Running all the logistic regressions, one per stage, for cold-water and warm-water organisms separately
cold.sig.all<-allregress(slc.cold,TRUE)
hot.sig.all<-allregress(slc.hot,TRUE)

##Meta-analysis (over the many logistic regressions) and tabulating the meta-coefficients

#First average and tabulate warm-water meta-regression intercepts
#Non-Hyperthermal temperature at hot edge
moddat<-hot.sig.all[hot.sig.all$TermName=="(Intercept)"& !hot.sig.all$Slc%in%c(hyperther,excluhyperther,Hyperther.Scot),]
modback<-rma(yi=as.numeric(moddat$Estimate),
             sei=as.numeric(moddat$SE),
             method="REML")
#Hyperthermal temperature at hot edge
moddat<-hot.sig.all[hot.sig.all$TermName=="(Intercept)"&hot.sig.all$Slc%in%hyperther,]
modhyp<-rma(yi=as.numeric(moddat$Estimate),
            sei=as.numeric(moddat$SE),
            method="REML")
modback.pred<-predict(modback, transf=transf.ilogit, digits=3)
modhyp.pred<-predict(modhyp, transf=transf.ilogit, digits=3)
#proper interpretation of the back-transformed value (i.e., 0.036) is that it reflects the median infection risk (although it would often be interpreted as the average risk). 
#https://www.metafor-project.org/doku.php/analyses:stijnen2010
warmInterincr2<-rbind(c(modback.pred[[1]],modback.pred[[3]],modback.pred[[4]]),
                      c(modhyp.pred[[1]],modhyp.pred[[3]],modhyp.pred[[4]])
                      ,c(modhyp.pred[[1]]-modback.pred[[1]],modhyp.pred[[3]]-modback.pred[[3]],modhyp.pred[[4]]-modback.pred[[4]]))
warmInterincr<-rbind(c(modback$beta[1], modback$ci.lb[1], modback$ci.ub[1]),
                     c(modhyp$beta[1],modhyp$ci.lb[1], modhyp$ci.ub[1]),
                     c(modhyp$beta[1]-modback$beta[1],modhyp$ci.lb[1]-modback$beta[1], modhyp$ci.ub[1]-modback$beta[1]))
colnames(warmInterincr)<-c("beta","ci.lb","ci.ub")
colnames(warmInterincr2)<-c("beta","ci.lb","ci.ub")


#Average and tabulate warm-water meta-regression slopes
moddat<-hot.sig.all[hot.sig.all$TermName=="thermN"& !hot.sig.all$Slc%in%c(hyperther,excluhyperther,Hyperther.Scot),]
modback<-rma(yi=as.numeric(moddat$Estimate),
             sei=as.numeric(moddat$SE),
             method="REML")
summary(modback)

#Hyperthermal temperature at hot edge
moddat<-hot.sig.all[hot.sig.all$TermName=="thermN"&hot.sig.all$Slc%in%hyperther,]
modhyp<-rma(yi=as.numeric(moddat$Estimate),
            sei=as.numeric(moddat$SE),
            method="REML")
summary(modhyp)

#extinction risk increase for warm niches overlaps zero
modback.pred<-predict(modback, transf=transf.ilogit, digits=3)
modhyp.pred<-predict(modhyp, transf=transf.ilogit, digits=3)
warmincr2<-rbind(c(modback.pred[[1]],modback.pred[[3]],modback.pred[[4]]),
                 c(modhyp.pred[[1]],modhyp.pred[[3]],modhyp.pred[[4]])
                 ,c(modhyp.pred[[1]]-modback.pred[[1]],modhyp.pred[[3]]-modback.pred[[3]],modhyp.pred[[4]]-modback.pred[[4]]))
warmincr<-rbind(c(modback$beta[1], modback$ci.lb[1], modback$ci.ub[1]),
                c(modhyp$beta[1],modhyp$ci.lb[1], modhyp$ci.ub[1]),
                c(modhyp$beta[1]-modback$beta[1],modhyp$ci.lb[1]-modback$beta[1], modhyp$ci.ub[1]-modback$beta[1]))
colnames(warmincr)<-c("beta","ci.lb","ci.ub")
colnames(warmincr2)<-c("beta","ci.lb","ci.ub")


#Average and tabulate cold-water meta-regression intercepts
#General temperature at cold edge
moddat<-cold.sig.all[cold.sig.all$TermName=="(Intercept)"& !cold.sig.all$Slc%in%c(hyperther,excluhyperther,Hyperther.Scot),]
modback<-rma(yi=as.numeric(moddat$Estimate),
             sei=as.numeric(moddat$SE),
             method="REML")
#Hyperthermal temperature at cold edge
moddat<-cold.sig.all[cold.sig.all$TermName=="(Intercept)"&cold.sig.all$Slc%in%hyperther,]
modhyp<-rma(yi=as.numeric(moddat$Estimate),
            sei=as.numeric(moddat$SE),
            method="REML")
modback.pred<-predict(modback, transf=transf.ilogit, digits=3)
modhyp.pred<-predict(modhyp, transf=transf.ilogit, digits=3)
coldInterincr2<-rbind(c(modback.pred[[1]],modback.pred[[3]],modback.pred[[4]]),
                      c(modhyp.pred[[1]],modhyp.pred[[3]],modhyp.pred[[4]])
                      ,c(modhyp.pred[[1]]-modback.pred[[1]],modhyp.pred[[3]]-modback.pred[[3]],modhyp.pred[[4]]-modback.pred[[4]]))
coldInterincr<-rbind(c(modback$beta[1], modback$ci.lb[1], modback$ci.ub[1]),
                     c(modhyp$beta[1],modhyp$ci.lb[1], modhyp$ci.ub[1]),
                     c(modhyp$beta[1]-modback$beta[1], modhyp$ci.lb[1]-modback$beta[1], modhyp$ci.ub[1]-modback$beta[1]))
colnames(coldInterincr)<-c("beta","ci.lb","ci.ub")
colnames(coldInterincr2)<-c("beta","ci.lb","ci.ub")

#Average and tabulate cold-water meta-regression slopes
moddat<-cold.sig.all[cold.sig.all$TermName=="thermN"& !cold.sig.all$Slc%in%c(hyperther,excluhyperther,Hyperther.Scot),]
modback<-rma(yi=as.numeric(moddat$Estimate),
             sei=as.numeric(moddat$SE),
             method="REML")
summary(modback)

#Hyperthermal temperature at cold edge
moddat<-cold.sig.all[cold.sig.all$TermName=="thermN"&cold.sig.all$Slc%in%hyperther,]
modhyp<-rma(yi=as.numeric(moddat$Estimate),
            sei=as.numeric(moddat$SE),
            method="REML")
summary(modhyp)

#hyperthermal extinction risk increase for cold niches does not overlap zero
modback.pred<-predict(modback, transf=transf.ilogit, digits=3)
modhyp.pred<-predict(modhyp, transf=transf.ilogit, digits=3)
coldincr2<-rbind(c(modback.pred[[1]],modback.pred[[3]],modback.pred[[4]]),
                 c(modhyp.pred[[1]],modhyp.pred[[3]],modhyp.pred[[4]])
                 ,c(modhyp.pred[[1]]-modback.pred[[1]],modhyp.pred[[3]]-modback.pred[[3]],modhyp.pred[[4]]-modback.pred[[4]]))
coldincr<-rbind(c(modback$beta[1], modback$ci.lb[1], modback$ci.ub[1]),
                c(modhyp$beta[1],modhyp$ci.lb[1], modhyp$ci.ub[1]),
                c(modhyp$beta[1]-modback$beta[1], modhyp$ci.lb[1]-modback$beta[1], modhyp$ci.ub[1]-modback$beta[1]))
colnames(coldincr)<-c("beta","ci.lb","ci.ub")
colnames(coldincr2)<-c("beta","ci.lb","ci.ub")

warmInterincr
warmincr
coldInterincr
coldincr
#Intercepts negative because they are the weighted means of the slope intercepts

##Fig. 3 ----

png("Fig.3.Combined.png",6,6,units = "in",res=300) #6,6
layout(matrix(c(1,1,2,2,1,1,2,2,3,4,5,6),3,4,byrow=T))
par(mai=c(0.6,0.6,0.7,0))
dat<-slc.cold[!is.na(slc.cold$ext)& !is.na(slc.cold$thermN) & slc.cold$slc %in%hyperther,]
probdat<-cold.prob[!is.na(cold.prob$ext)& !is.na(cold.prob$thermN) & cold.prob$slc %in%hyperther,]
# Get the unique gpas
gpas <- sort( unique( probdat$thermN ) )

plot(ext ~ I(coldoffset-exp(thermN)), data = probdat,xlab="",ylab="Extinction risk, odds",main="(a) Cold-water genera",col="red",cex=0.8,ylim=c(0,1.02),las=1)

# Plot the model prediction for each GPA, with a line through them
polygon(c(rev(coldoffset-exp(gpas)), coldoffset-exp(gpas)), c(rev(inverselogit( coldInterincr[2,3] + gpas*coldincr[2,3] )), inverselogit( coldInterincr[2,2] + gpas*coldincr[2,2] )), col = alpha("red",0.3), border = NA)
lines( coldoffset-exp(gpas), inverselogit( coldInterincr[2,1] + gpas*coldincr[2,1] ), col="red", lwd=2 )

polygon(c(rev(coldoffset-exp(gpas)), coldoffset-exp(gpas)), c(rev(inverselogit( coldInterincr[1,3] + gpas*coldincr[1,3] )), inverselogit( coldInterincr[1,2] + gpas*coldincr[1,2] )), col = alpha("blue",0.3), border = NA)
lines( coldoffset-exp(gpas), inverselogit( coldInterincr[1,1] + gpas*coldincr[1,1] ), col="blue", lwd=2 )

par(mai=c(0.6,0.2,0.7,0.4))
##The hotwater plot
dat<-slc.hot[!is.na(slc.hot$ext)& !is.na(slc.hot$thermN) & slc.hot$slc %in%hyperther,]
probdat<-hot.prob[!is.na(hot.prob$ext)& !is.na(hot.prob$thermN) & hot.prob$slc %in%hyperther,]
# Get the unique gpas
gpas <- sort( unique( probdat$thermN ) )
plot(ext ~ I(exp(thermN)+hotoffset), data = probdat,xlab="Genus thermal optimum, degrees Celsius",ylab="",main="(b) Warm-water genera",axes=F,col="red",cex=0.8,ylim=c(0,1.02))
axis(1);axis(2,labels = F);box()

# Plot the model prediction for each GPA, with a line through them
polygon(c(rev(exp(gpas)+hotoffset), exp(gpas)+hotoffset), c(rev(inverselogit( warmInterincr[2,3] + gpas*warmincr[2,3] )), inverselogit( warmInterincr[2,2] + gpas*warmincr[2,2] )), col = alpha("red",0.3), border = NA)
lines( exp(gpas)+hotoffset, inverselogit( warmInterincr[2,1] + gpas*warmincr[2,1] ), col="red", lwd=2 )

polygon(c(rev(exp(gpas)+hotoffset), exp(gpas)+hotoffset), c(rev(inverselogit( warmInterincr[1,3] + gpas*warmincr[1,3] )), inverselogit( warmInterincr[1,2] + gpas*warmincr[1,2] )), col = alpha("blue",0.3), border = NA)
lines( exp(gpas)+hotoffset, inverselogit( warmInterincr[1,1] + gpas*warmincr[1,1] ), col="blue", lwd=2 )

#########four auxilliary plots

#Warm niche first two
par(mai=c(1,0.8,0.8,0.4)+0.02)
par(mar=c(4,4,2.5,0)+0.1)
#Cold niche, first two
plot(c(1.5,2,2.5),coldInterincr[,1],ylim=c(-5,0.8),axes=F,ylab="Extinction risk, log-odds",xlab="",main="(c) Intercept",xlim=c(1,3))
axis(2,las=2);abline(h=0)
arrows(c(1.5,2,2.5),coldInterincr[,2],c(1.5,2,2.5),coldInterincr[,3],length=0,lwd=2,col=c("blue","red","black"))
points(c(1.5,2,2.5),coldInterincr[,1],ylim=c(-0.1,0.2),pch=16,col=c("blue","red","black"))
box()

plot(c(1.5,2,2.5),coldincr[,1],ylim=c(0,1.3),axes=F,ylab="",xlab="",main="(d) Slope",xlim=c(1,3))
axis(2,las=2);abline(h=0)
arrows(c(1.5,2,2.5),coldincr[,2],c(1.5,2,2.5),coldincr[,3],length=0,lwd=2,col=c("blue","red","black"))
points(c(1.5,2,2.5),coldincr[,1],ylim=c(-0.1,0.2),pch=16,col=c("blue","red","black"))
box()

plot(c(1.5,2,2.5),warmInterincr[,1],ylim=c(-2.5,2),axes=F,ylab="",xlab="",main="(e) Intercept",xlim=c(1,3)) #ylim=c(-0.1,0.2)
axis(2,las=2);abline(h=0)
arrows(c(1.5,2,2.5),warmInterincr[,2],c(1.5,2,2.5),warmInterincr[,3],length=0,lwd=2,col=c("blue","red","black"))
points(c(1.5,2,2.5),warmInterincr[,1],ylim=c(-0.1,0.2),pch=16,col=c("blue","red","black"))
box()

plot(c(1.5,2,2.5),warmincr[,1],ylim=c(-0.5,0.75),axes=F,ylab="",xlab="",main="(f) Slope",xlim=c(1,3)) #ylim=c(-0.1,0.2)
axis(2,las=2);abline(h=0)
arrows(c(1.5,2,2.5),warmincr[,2],c(1.5,2,2.5),warmincr[,3],length=0,lwd=2,col=c("blue","red","black"))
points(c(1.5,2,2.5),warmincr[,1],ylim=c(-0.1,0.2),pch=16,col=c("blue","red","black"))
box()

dev.off()

##Deeper meta-analysis exploration ----

#Hyperthermal temperature at cold edge (these are the same as above)
moddat<-cold.sig.all[cold.sig.all$TermName=="thermN"&cold.sig.all$Slc%in%hyperther,]
coldslo<-rma(yi=as.numeric(moddat$Estimate),
             sei=as.numeric(moddat$SE),
             method="REML")

#Hyperthermal temperature at cold edge
moddat<-cold.sig.all[cold.sig.all$TermName=="(Intercept)"&cold.sig.all$Slc%in%hyperther,]
coldint<-rma(yi=as.numeric(moddat$Estimate),
             sei=as.numeric(moddat$SE),
             method="REML")

#Hyperthermal temperature at hot edge
moddat<-hot.sig.all[hot.sig.all$TermName=="thermN"&hot.sig.all$Slc%in%hyperther,]
hotslo<-rma(yi=as.numeric(moddat$Estimate),
            sei=as.numeric(moddat$SE),
            method="REML")

#Hyperthermal temperature at hot edge
moddat<-hot.sig.all[hot.sig.all$TermName=="(Intercept)"&hot.sig.all$Slc%in%hyperther,]
hotint<-rma(yi=as.numeric(moddat$Estimate),
            sei=as.numeric(moddat$SE),
            method="REML")

##repeatedly fit the specified model, leaving out one observation/study at a time.
#a useful tool to assess the influence of a single study on the meta-analysis results and for identifying potential outliers.
#important for the hyperthermals due to low sample size 
loo.coldslo<-rbind(as.numeric(c(NA,coldslo$beta,coldslo$se,coldslo$zval,coldslo$pval,coldslo$ci.lb,coldslo$ci.ub,coldslo$QE,coldslo$QEp,coldslo$tau2,coldslo$I2,coldslo$H2)),data.frame(stages$short[hyperther],leave1out(coldslo,digits =2))) 
loo.coldint<-rbind(as.numeric(c(NA,coldint$beta,coldint$se,coldint$zval,coldint$pval,coldint$ci.lb,coldint$ci.ub,coldint$QE,coldint$QEp,coldint$tau2,coldint$I2,coldint$H2)),data.frame(stages$short[hyperther],leave1out(coldint,digits =2))) 
loo.hotslo<-rbind(as.numeric(c(NA,hotslo$beta,hotslo$se,hotslo$zval,hotslo$pval,hotslo$ci.lb,hotslo$ci.ub,hotslo$QE,hotslo$QEp,hotslo$tau2,hotslo$I2,hotslo$H2)),data.frame(stages$short[hyperther],leave1out(hotslo,digits =2))) 
loo.hotint<-rbind(as.numeric(c(NA,hotint$beta,hotint$se,hotint$zval,hotint$pval,hotint$ci.lb,hotint$ci.ub,hotint$QE,hotint$QEp,hotint$tau2,hotint$I2,hotint$H2)),data.frame(stages$short[hyperther],leave1out(hotint,digits =2))) 
#write.csv(rbind(loo.coldslo,loo.coldint,loo.hotslo,loo.hotint),file="LOO.all.csv")

summary(coldslo)
summary(coldint)
summary(hotslo)
summary(hotint)
#-REVC tau squared, random effects error (between studies)
#-Cochran's Q is measure of heterogeneity (i.e. tests whether REVC = 0) i.e. whether to include moderators or not
#coldslo Q 0.1520, 40.29% 1.67. = slope unlikely affected by much through time
#coldint Q < .0001, 73.93% 3.84
#hotslo Q 0.0015, 68.30%, 3.15
#hotint Q < .0001, 84.37%, 6.40
#-I2 is the percentage of total variation due to true differences rather than/relative to chance/sampling
#I^2 (total heterogeneity / total variability):   40.29%
# A value of H2 = 1 indicates perfect homogene- ity among the studies

## Fig. 4. Forest plots ----
layout(matrix(1:4,2,2,byrow=T))
forest(coldslo,slab=cold.sig.all[cold.sig.all$TermName=="thermN"&cold.sig.all$Slc%in%hyperther,"Stage"],main="(a) Cold-water slope")
abline(v=0.1932854)
forest(coldint,slab=cold.sig.all[cold.sig.all$TermName=="thermN"&cold.sig.all$Slc%in%hyperther,"Stage"],main="(b) Cold-water intercept")
abline(v=-2.4985704)
forest(hotslo,slab=hot.sig.all[hot.sig.all$TermName=="thermN"&hot.sig.all$Slc%in%hyperther,"Stage"],main="(c) Hot-water slope")
abline(v=-0.06949686)
forest(hotint,slab=hot.sig.all[hot.sig.all$TermName=="thermN"&hot.sig.all$Slc%in%hyperther,"Stage"],main="(d) Hot-water intercept")
abline(v=-2.057657)


#Drivers: variable distributions ----

intercold<-cold.sig.all[cold.sig.all$TermName=="(Intercept)",c("Estimate","SE")]
interhot<-hot.sig.all[hot.sig.all$TermName=="(Intercept)",c("Estimate","SE")]
Ncold<-cold.sig.all[cold.sig.all$TermName=="thermN",c("Estimate","SE")]
Nhot<-hot.sig.all[hot.sig.all$TermName=="thermN",]
hotcolddat<-cbind(as.numeric(intercold[,1]),as.numeric(intercold[,2]),as.numeric(interhot[,1]),as.numeric(interhot[,2]),as.numeric(Ncold[,1]),as.numeric(Ncold[,2]),as.numeric(Nhot[,2]),as.numeric(Nhot[,3]),Nhot[,c(5:11,13:16)])
colnames(hotcolddat)<-c("intercoldB","intercoldSE","interhotB","interhotSE","nichecoldB","nichecoldSE","nichehotB","nichehotSE",colnames(hotcolddat)[9:19])

#get the right slc data and add log-transformed reef occurrence proportions
reefcoll<-spdat1[!duplicated(spdat1$collection_no),c("reef","slc")] #check collections
proptry<-tapply(reefcoll$reef,reefcoll$slc,function(g) sum(g=="reef",na.rm=T)/sum(g%in%c("reef","non-reef"),na.rm=T) ) #p-value = 0.006 after hyperthermal
hotcolddat$reefprop<-logit(proptry[row.names(proptry)%in%hotcolddat$Slc],adjust=0.003787879) #uses the lowest non-zero proportion

##Fig. 6a. Lat max ----
layout(matrix(1:3,1,3,byrow = T))
boxplot(cold.sig.all[cold.sig.all$TermName=="thermN" &!cold.sig.all$Slc%in%c(hyperther,excluhyperther),"Lat_max"],main="(a)",  #used to be 4-cold.sig.all[...
        cold.sig.all[cold.sig.all$TermName=="thermN" &cold.sig.all$Slc%in%hyperther,"Lat_max"],border=c("blue","red"),ylab="Latitudinal change, degrees",axes=F)
axis(1,at=1:2,labels=c("Non-hyper","Hyper")); box(); axis(2)
#axis(2,at=4-c(2.5,3,3.5),labels=round(100-exp(c(2.5,3,3.5)),digits = 1))   #log(100-samplvars[5:78,2])
points(rep(2,length(hyperther)),cold.sig.all[cold.sig.all$TermName=="thermN" &cold.sig.all$Slc%in%hyperther,"Lat_max"],pch=16,col="red")
text(rep(2,length(hyperther)),cold.sig.all[cold.sig.all$TermName=="thermN" &cold.sig.all$Slc%in%hyperther,"Lat_max"],labels=stages$short[hyperther],col="red",adj=c(-0.3,0.5))

##Fig. 6b. temperature ----
boxplot(cold.sig.all[cold.sig.all$TermName=="thermN" &!cold.sig.all$Slc%in%c(hyperther+1,hyperther,excluhyperther),"temperature"],main="(b)",
        cold.sig.all[cold.sig.all$TermName=="thermN" &cold.sig.all$Slc%in%(hyperther+1),"temperature"],border=c("blue","red"),ylab="Temperature, degrees C",axes=F)
axis(1,at=1:2,labels=c("Non-hyper","Hyper+1")); box()
axis(2,at=c(2.8,3,3.2,3.4,3.6),labels=round(exp(c(2.8,3,3.2,3.4,3.6)),digits = 1))   #log(100-samplvars[5:78,2])
points(rep(2,length(hyperther)),cold.sig.all[cold.sig.all$TermName=="thermN" &cold.sig.all$Slc%in%(hyperther+1),"temperature"],pch=16,col="red")
text(rep(2,length(hyperther)),cold.sig.all[cold.sig.all$TermName=="thermN" &cold.sig.all$Slc%in%(hyperther+1),"temperature"],labels=paste(stages$short[hyperther],"+1"),col="red",adj=c(-0.3,0.5))

##Fig. 6c. reefs proportion of all ----
proptry<-tapply(reefcoll$reef,reefcoll$slc,function(g) sum(g=="reef",na.rm=T)/sum(g%in%c("reef","non-reef"),na.rm=T) ) #p-value = 0.006 after hyperthermal
boxplot(proptry[-c(hyperther-35+1,hyperther-35,excluhyperther-35)],main="(c)",
        proptry[hyperther-35+1],border=c("blue","red"),ylab="Reef collections, proportion of total")
axis(1,at=1:2,labels=c("Non-hyper","Hyper+1")); box()
points(rep(2,length(hyperther)),proptry[hyperther-35+1],pch=16,col="red")
text(rep(2,length(hyperther)),proptry[hyperther-35+1],labels=paste(stages$short[hyperther],"+1"),col="red",adj=c(-0.3,0.5))
#Note the carnian is likely an outlier because it is followed by the 19.5 myr long Norian. However, carbonate


##Fig. S6. Pleistocene vs non-hyperthermal shelf area ----
png("PleistoceneShelfDistribution.png",width = 7,height = 7,units = "in",res=300)
layout(matrix(1:4,2,2,byrow=T))
boxplot(Extra.Env.dat[3]/1000000,ylab="Tropical shelf area, million km^2",main="(A)"); points(rep(1,length(hyperther)),Extra.Env.dat[3]/1000000,pch=16,col="red") 
points(1,Extra.Env.dat[95,3]/1000000,pch=8,cex=2,lwd=2) #tropshelfarea
text(rep(1,length(hyperther)),Extra.Env.dat[hyperther,3]/1000000,labels=stages$short[hyperther],col="red",adj=c(-0.3,0.5))

boxplot(Extra.Env.dat[4]/1000000,ylab="Temperate and polar shelf area, million km^2",main="(B)"); points(rep(1,length(hyperther)),Extra.Env.dat[4]/1000000,pch=16,col="red")
points(1,Extra.Env.dat[95,4]/1000000,pch=8,cex=2,lwd=2) #tempshelfarea 
text(rep(1,length(hyperther)),Extra.Env.dat[hyperther,4]/1000000,labels=stages$short[hyperther],col="red",adj=c(-0.3,0.5))

boxplot(log(Extra.Env.dat[4]/Extra.Env.dat[3]),ylab="Temperate:tropical shelf area, log ratio",main="(C)"); points(rep(1,length(hyperther)),log(Extra.Env.dat[hyperther,4]/Extra.Env.dat[hyperther,3]),pch=16,col="red") 
points(1,log(Extra.Env.dat[95,4]/Extra.Env.dat[95,3]),pch=8,cex=2,lwd=2) #troptempshelfarea. lower values are more tropical than temperate area
text(rep(1,length(hyperther)),log(Extra.Env.dat[hyperther,4]/Extra.Env.dat[hyperther,3]),labels=stages$short[hyperther],col="red",adj=c(-0.3,0.5))

boxplot(exp(hotcolddat$temperature),ylab="Global SST median, degrees Celsius",main="(D)")
points(rep(1,length(hyperther)),exp(hotcolddat$temperature[hotcolddat$Slc%in%hyperther]),pch=16,col="red") 
text(rep(1,length(hyperther)),exp(hotcolddat$temperature[hotcolddat$Slc%in%hyperther]),labels=stages$short[hyperther],col="red",adj=c(-0.3,0.5))

dev.off()

##Fig. S7. ----

png("Proportion of flooded shelf area.png",5,8,res=300,units="in")
layout(1:2)
tsplot(stages[43:91,],boxes="sys", shading="series",ylim=range(Extra.Env.dat[43:91,1]),ylab="Proportion",labels.args=list(cex=0.7))
lines(stages[43:91,"mid"],Extra.Env.dat[43:91,1],lwd=2); title(main="(a) Cold-water")
abline(h=Extra.Env.dat[95,1],lty="dashed");abline(v=stages[hyperther,"mid"],col="red")
median(Extra.Env.dat[43:91,1]);Extra.Env.dat[95,1] #<15

tsplot(stages[43:91,],boxes="sys", shading="series",ylim=range(Extra.Env.dat[43:91,2]),ylab="Proportion",labels.args=list(cex=0.7))
lines(stages[43:91,"mid"],Extra.Env.dat[43:91,2],lwd=2); title(main="(b) Warm-water")
abline(h=Extra.Env.dat[95,2],lty="dashed");abline(v=stages[hyperther,"mid"],col="red")
median(Extra.Env.dat[43:91,2]);Extra.Env.dat[95,2] #>30
dev.off()


#Meta-analysis model selection with or without temperature, for each coefficient ----

rma.glmulti <- function(formula, data, ...)   rma(formula, sei=nichecoldSE, data=hotcolddat, method="REML", ...)
res <- glmulti(nichecoldB ~ temperature+tempshel+Occur.med+nrec, data=hotcolddat,level=1, fitfunction=rma.glmulti, crit="bic")
res <- glmulti(nichecoldB ~ tempshel+Occur.med+nrec+Peterssili+Lat_max, data=hotcolddat,level=1, fitfunction=rma.glmulti, crit="bic")
summary(res@objects[[1]])

rma.glmulti <- function(formula, data, ...)   rma(formula, sei=nichehotSE, data=hotcolddat, method="REML", ...)
res <- glmulti(nichehotB ~  temperature+tropshel+troptempshel+Occur.med+nrec, data=hotcolddat,level=1, fitfunction=rma.glmulti, crit="bic")
res <- glmulti(nichehotB ~  tropshel+troptempshel+Occur.med+nrec+Peterssili+Lat_min+I(c(reefprop[-1],NA)), data=hotcolddat,level=1, fitfunction=rma.glmulti, crit="bic")
summary(res@objects[[1]])

#cold and hot intercepts
rma.glmulti <- function(formula, data, ...)   rma(formula, sei=intercoldSE, data=hotcolddat, method="REML", ...)
res <- glmulti(intercoldB ~ temperature+tempshel+Occur.med+nrec, data=hotcolddat,level=1, fitfunction=rma.glmulti, crit="bic") #crit="aicc"
res <- glmulti(intercoldB ~ tempshel+Occur.med+nrec+Peterssili+Lat_max, data=hotcolddat,level=1, fitfunction=rma.glmulti, crit="bic") #crit="aicc"
summary(res@objects[[1]])

rma.glmulti <- function(formula, data, ...)   rma(formula, sei=interhotSE, data=hotcolddat, method="REML", ...)
res <- glmulti(interhotB ~ temperature+tropshel+troptempshel+Occur.med+nrec, data=hotcolddat,level=1, fitfunction=rma.glmulti, crit="bic")
res <- glmulti(interhotB ~ tropshel+troptempshel+Occur.med+nrec+Peterssili+Lat_min+I(c(reefprop[-1],NA)), data=hotcolddat,level=1, fitfunction=rma.glmulti, crit="bic")
summary(res@objects[[1]])

print(res) #Overview of model selection for whichever object above
summary(res@objects[[1]]) #Detail on best model, but...

tmp <- weightable(res)
tmp[1:10,1] #to check which model has which variable
summary(res@objects[[9]]) 
#Go through each independent variable and see whether it is present in the top 10 optimal models:
length(grep("temperature",tmp[1:10,1]))/length(tmp[1:10,1]) 
length(grep("Lat_max",tmp[1:10,1]))/length(tmp[1:10,1]) 
length(grep("Lat_min",tmp[1:10,1]))/length(tmp[1:10,1]) 
length(grep("Peterssili",tmp[1:10,1]))/length(tmp[1:10,1]) 
length(grep("reefprop",tmp[1:10,1]))/length(tmp[1:10,1]) 

length(grep("tempshel",tmp[1:10,1]))/length(tmp[1:10,1]) 
length(grep("tropshel",tmp[1:10,1]))/length(tmp[1:10,1]) 
length(grep("troptempshel",tmp[1:10,1]))/length(tmp[1:10,1]) 

plot(res@objects[[5]],type="p")


##Time series with gls ----
hotcolddat$PeterssiliNEXT<- c(hotcolddat$Peterssili[2:49],-0.065833551,NA,NA,NA)
hotcolddat$reefpropNEXT<- c(hotcolddat$reefprop[-1],NA)

#without temp
mod.gls <- gls(nichecoldB ~ tempshel+Occur.med+nrec+Peterssili+Lat_max, data=hotcolddat[-50:-52,],na.action=na.omit,
               correlation=NULL, method="ML")
mod.gls <- gls(nichehotB ~ tropshel+troptempshel+Occur.med+nrec+Peterssili+Lat_min+reefpropNEXT, data=hotcolddat[-49:-52,],na.action=na.omit,
               correlation=corARMA(p=1), method="ML")
mod.gls <- gls(intercoldB ~ tempshel+Occur.med+nrec+Peterssili+Lat_max, data=hotcolddat[-49:-52,],na.action=na.omit,
               correlation=corARMA(p=1), method="ML")
mod.gls <- gls(interhotB ~ tropshel+troptempshel+Occur.med+nrec+Peterssili+Lat_min+reefpropNEXT, data=hotcolddat[-49:-52,],na.action=na.omit,
               correlation=corARMA(p=1), method="ML")

#with temp
mod.gls <- gls(nichecoldB ~ temperature+tempshel+Occur.med+nrec+PeterssiliNEXT+Lat_max, data=hotcolddat[-50:-52,],na.action=na.omit,
               correlation=NULL, method="ML")
mod.gls <- gls(nichehotB ~ temperature+tropshel+troptempshel+Occur.med+nrec+PeterssiliNEXT+Lat_min+reefpropNEXT, data=hotcolddat[-49:-52,],na.action=na.omit,
               correlation=corARMA(p=1), method="ML")
mod.gls <- gls(intercoldB ~ temperature+tempshel+Occur.med+nrec+PeterssiliNEXT+Lat_max, data=hotcolddat[-49:-52,],na.action=na.omit,
               correlation=corARMA(p=1), method="ML")
mod.gls <- gls(interhotB ~ temperature+tropshel+troptempshel+Occur.med+nrec+PeterssiliNEXT+Lat_min+reefpropNEXT, data=hotcolddat[-49:-52,],na.action=na.omit,
               correlation=corARMA(p=1), method="ML")

#Automate the dropping of terms
for(i in 1:length(mod.gls$coefficients)){
  dropping<-drop1(mod.gls)
  todrop<-row.names(dropping)[dropping$AIC==min(dropping$AIC)]
  print(todrop)
  if(!todrop=="<none>"){
    mod.gls <- update(mod.gls, paste(". ~ . -",todrop))
  }}
summary(mod.gls)
nagelkerke(mod.gls)

#For deciding the starting autoregression term
mod.gls.2 <- update(mod.gls, correlation=corARMA(p=2))
mod.gls.1 <- update(mod.gls, correlation=corARMA(p=1))
mod.gls.0 <- update(mod.gls, correlation=NULL)
anova(mod.gls,mod.gls.2, mod.gls.1,mod.gls.0)

##R2 and splitting out effects ----
#Obtaining a pseudo R2. pseudo R-squared measures are relative measures among similar models indicating how well the model explains the data
#Low R2 means little chance of predictive utility (too much variation in data)
nullmod <- gls(nichecoldB~1,data=hotcolddat,na.action=na.omit,correlation=NULL, method="ML")
nullmod <- gls(nichehotB~1,data=hotcolddat,na.action=na.omit,correlation=corARMA(p=1), method="ML")
nullmod <- gls(intercoldB~1,data=hotcolddat,na.action=na.omit,correlation=corARMA(p=1), method="ML")
nullmod <- gls(interhotB~1,data=hotcolddat,na.action=na.omit,correlation=corARMA(p=1), method="ML")
1-(as.numeric(logLik(mod.gls)/logLik(nullmod))) #McFadden

R2 <- cor(hotcolddat$nichecoldB[-50:-52],predict(mod.gls))^2
R2

R2.1 <- 1 - with(hotcolddat[-50:-52,], (sum((nichecoldB-predict(mod.gls))^2)/sum((nichecoldB-mean(nichecoldB,na.rm=T))^2,na.rm=T)))
R2.1

nagelkerke(mod.gls)


##Checking effect of continental shelf area alone i.e. not change in  continental shelf area ----
#switch to unscaled for prediction
hotcolddat$raster_area.temp<-Extra.Env.dat[as.numeric(hotcolddat$Slc),4]
hotcolddat$raster_area.trop<-Extra.Env.dat[as.numeric(hotcolddat$Slc),3]

#regressive terms are all optimal
Sub <- hotcolddat[-50:-52,]
Sub$raster_area.temp<-scale(Sub$raster_area.temp)
Sub$raster_area.trop<-scale(Sub$raster_area.trop)
Sub$temperature<-scale(Sub$temperature)
Sub$Occur.med<-scale(Sub$Occur.med)
Sub$nrec<-scale(Sub$nrec)

#In short, the change in it has little effect, but 
#To visualise the (always positive) interaction, use visreg()
Sub2<-Sub
names(Sub2)[names(Sub2)=="nichecoldSE"]<-paste("vi"); Sub2$vi<-Sub2$vi^2
rma.glmulti <- function(formula, data, ...)   rma.mv(formula, vi, data=data,  method="ML", ...)
res <- glmulti(nichecoldB ~ temperature+raster_area.temp*Occur.med+nrec, data=Sub2,level=1, fitfunction=rma.glmulti, crit="bic")
summary(res@objects[[1]])

Sub2<-Sub
names(Sub2)[names(Sub2)=="nichehotSE"]<-paste("vi"); Sub2$vi<-Sub2$vi^2
rma.glmulti <- function(formula, data, ...)   rma.mv(formula, vi, data=data,  method="ML", ...)
res <- glmulti(nichehotB ~ temperature+raster_area.trop*Occur.med+nrec, data=Sub2,level=1, fitfunction=rma.glmulti, crit="bic")
summary(res@objects[[1]])
res <- rma(nichehotB ~ temperature, data=Sub,method="REML",sei=intercoldSE)
summary(res)

#cold and hot intercepts
Sub2<-Sub
names(Sub2)[names(Sub2)=="intercoldSE"]<-paste("vi"); Sub2$vi<-Sub2$vi^2
rma.glmulti <- function(formula, data, ...)   rma.mv(formula, vi, data=data,  method="ML", ...)
res <- glmulti(intercoldB ~ temperature+raster_area.temp*Occur.med+nrec, data=Sub2,level=1, fitfunction=rma.glmulti, crit="bic")
summary(res@objects[[1]])
res <- rma(intercoldB ~ raster_area.temp*Occur.med+nrec, data=Sub,method="REML",sei=intercoldSE)
summary(res)

Sub2<-Sub
names(Sub2)[names(Sub2)=="interhotSE"]<-paste("vi"); Sub2$vi<-Sub2$vi^2
rma.glmulti <- function(formula, data, ...)   rma.mv(formula, vi, data=data,  method="ML", ...)
res <- glmulti(interhotB ~ temperature+raster_area.trop*Occur.med+nrec, data=Sub2,level=1, fitfunction=rma.glmulti, crit="bic")
summary(res@objects[[1]])
res <- rma(interhotB ~  temperature+Occur.med, data=Sub,sei=interhotSE,method="REML")
summary(res)

nullmod <- gls(nichecoldB~1,data=hotcolddat,na.action=na.omit,correlation=NULL, method="ML")
nullmod <- gls(nichehotB~1,data=hotcolddat,na.action=na.omit,correlation=corARMA(p=1), method="ML")
nullmod <- gls(intercoldB~1,data=hotcolddat,na.action=na.omit,correlation=corARMA(p=1), method="ML")
nullmod <- gls(interhotB~1,data=hotcolddat,na.action=na.omit,correlation=corARMA(p=1), method="ML")
1-(as.numeric(logLik(mod.gls)/logLik(nullmod))) #McFadden

####the interaction using gls (autocorrelation term already optimised)
mod.gls <- gls(nichecoldB ~ temperature+raster_area.temp*Occur.med+nrec, data=Sub,na.action=na.omit,
               correlation=NULL, method="ML")
mod.gls <- gls(nichehotB ~ temperature+raster_area.trop*Occur.med+nrec, data=Sub,na.action=na.omit,
               correlation=corARMA(p=1), method="ML")
mod.gls <- gls(intercoldB ~ temperature+raster_area.temp*Occur.med+nrec, data=Sub,na.action=na.omit,
               correlation=corARMA(p=1), method="ML")
mod.gls <- gls(interhotB ~ temperature+raster_area.trop*Occur.med+nrec, data=Sub,na.action=na.omit,
               correlation=corARMA(p=1), method="ML")
summary(mod.gls)
visreg(mod.gls,points=list(cex=1, pch=1))
visreg(mod.gls,"raster_area.temp",by="Occur.med") #By default, cross-sections are taken at the 10th, 50th, and 90th quantiles
visreg(mod.gls,"raster_area.trop",by="Occur.med") #By default, cross-sections are taken at the 10th, 50th, and 90th quantiles
nagelkerke(mod.gls)

#Model selection for gls
for(i in 1:length(mod.gls$coefficients)){
  dropping<-drop1(mod.gls)
  todrop<-row.names(dropping)[dropping$AIC==min(dropping$AIC)]
  print(todrop)
  if(!todrop=="<none>"){
    mod.gls <- update(mod.gls, paste(". ~ . -",todrop))
  }}


#Alternative analysis: Hyperthermals as a grouping variable ----
#Try having hyperthermal vs not as a categorical. Ideally this would be continuous
hypercat<-numeric(length(hotcolddat$Slc))
hypercat[hotcolddat$Slc%in%hyperther]<-1 
hypercat[hotcolddat$Slc%in%c(excluhyperther,Hyperther.Scot[!Hyperther.Scot%in%hyperther])]<-NA 
hotcolddat2<-cbind(hotcolddat,hypercat)

mod1<-rma(yi=intercoldB,mods = ~ hypercat,sei=intercoldSE,data=hotcolddat2[-50:-52,],method="REML")
mod0<-rma(yi=intercoldB,sei=intercoldSE,data=hotcolddat2[is.na(hotcolddat2$hypercat)==F & is.na(hotcolddat2$temperature)==F,],method="REML")
anova(mod1,mod0) #full preferred but not by BIC
summary(mod1)

mod1<-rma(yi=interhotB,mods = ~ hypercat,sei=interhotSE,data=hotcolddat2[-50:-52,],method="REML")
mod0<-rma(yi=interhotB,sei=interhotSE,data=hotcolddat2[is.na(hotcolddat2$hypercat)==F & is.na(hotcolddat2$temperature)==F,],method="REML")
anova(mod1,mod0) #full preferred
summary(mod1)

mod1<-rma(yi=nichecoldB,mods = ~ hypercat,sei=nichecoldSE,data=hotcolddat2[-50:-52,],method="REML")
mod0<-rma(yi=nichecoldB,sei=nichecoldSE,data=hotcolddat2[is.na(hotcolddat2$hypercat)==F & is.na(hotcolddat2$temperature)==F,],method="REML")
anova(mod1,mod0) #full preferred
summary(mod1) #the R2 is the effect of the moderator alone

mod1<-rma(yi=nichehotB,mods = ~ hypercat,sei=nichehotSE,data=hotcolddat2[-50:-52,],method="REML")
mod0<-rma(yi=nichehotB,sei=nichehotSE,data=hotcolddat2[is.na(hotcolddat2$hypercat)==F & is.na(hotcolddat2$temperature)==F,],method="REML")
anova(mod1,mod0) #reduced preferred


#Fig S4a. high tolerance plot ----

newdat = data.frame(thermMax = seq(min(slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax),]$thermMax), max(slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax),]$thermMax), length.out = 100))

#png("ThermSelectOverallCurves.png",width = 5, height = 9, units = "in", res = 300)
#layout(matrix(c(1,1,1,2),4,1))

plot(ext ~ thermMax, data = slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax),],xlab="Genus lower thermal tolerance, degrees Celsius",ylab="Extinction odds",main="(B)")

m10<-glm(ext~poly(thermMax,2),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax) & !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m15<-glm(ext~poly(thermMax,3),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m20<-glm(ext~poly(thermMax,4),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m25<-glm(ext~poly(thermMax,5),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m30<-glm(ext~poly(thermMax,6),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m35<-glm(ext~poly(thermMax,7),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m40<-glm(ext~poly(thermMax,8),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
BIC(m10,m15,m20,m25,m30,m35,m40) #15
newdat$pred = predict(m15, newdata = newdat,type="response")
with(newdat, lines(x = thermMax, y = pred,lwd=2,col="blue"))

m5<-glm(ext~thermMax,family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& slc.dat$slc %in%hyperther,])  #}
m10<-glm(ext~poly(thermMax,2),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax) & slc.dat$slc %in%hyperther,])  #}
m15<-glm(ext~poly(thermMax,3),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& slc.dat$slc %in%hyperther,])  #}
m20<-glm(ext~poly(thermMax,4),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& slc.dat$slc %in%hyperther,])  #}
m25<-glm(ext~poly(thermMax,5),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& slc.dat$slc %in%hyperther,])  #}
m30<-glm(ext~poly(thermMax,6),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& slc.dat$slc %in%hyperther,])  #}
m35<-glm(ext~poly(thermMax,7),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& slc.dat$slc %in%hyperther,])  #}
m40<-glm(ext~poly(thermMax,8),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMax)& slc.dat$slc %in%hyperther,])  #}
BIC(m10,m15,m20,m25,m30,m35,m40) #15
newdat$pred = predict(m15, newdata = newdat,type="response")
with(newdat, lines(x = thermMax, y = pred,lwd=2,col="red"))

##Fig. S4b. low tolerance plot ----

newdat = data.frame(thermMin = seq(min(slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin),]$thermMin), max(slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin),]$thermMin), length.out = 100))

plot(ext ~ thermMin, data = slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin),],xlab="Genus upper thermal tolerance, degrees Celsius",ylab="Extinction odds",main="(A)")

m10<-glm(ext~poly(thermMin,2),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin) & !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m15<-glm(ext~poly(thermMin,3),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m20<-glm(ext~poly(thermMin,4),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m25<-glm(ext~poly(thermMin,5),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m30<-glm(ext~poly(thermMin,6),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m35<-glm(ext~poly(thermMin,7),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
m40<-glm(ext~poly(thermMin,8),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& !slc.dat$slc %in%c(hyperther,excluhyperther),])  #}
#BIC(m10,m15,m20) #15
BIC(m10,m15,m20,m25,m30,m35,m40) #15
newdat$pred = predict(m20, newdata = newdat,type="response")
with(newdat, lines(x = thermMin, y = pred,lwd=2,col="blue"))

m5<-glm(ext~thermMin,family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& slc.dat$slc %in%hyperther,])  #}
m10<-glm(ext~poly(thermMin,2),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin) & slc.dat$slc %in%hyperther,])  #}
m15<-glm(ext~poly(thermMin,3),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& slc.dat$slc %in%hyperther,])  #}
m20<-glm(ext~poly(thermMin,4),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& slc.dat$slc %in%hyperther,])  #}
m25<-glm(ext~poly(thermMin,5),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& slc.dat$slc %in%hyperther,])  #}
m30<-glm(ext~poly(thermMin,6),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& slc.dat$slc %in%hyperther,])  #}
m35<-glm(ext~poly(thermMin,7),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& slc.dat$slc %in%hyperther,])  #}
m40<-glm(ext~poly(thermMin,8),family="binomial",data=slc.dat[!is.na(slc.dat$ext)& !is.na(slc.dat$thermMin)& slc.dat$slc %in%hyperther,])  #}
BIC(m10,m15,m20,m25,m30,m35,m40) #15
newdat$pred = predict(m15, newdata = newdat,type="response")
with(newdat, lines(x = thermMin, y = pred,lwd=2,col="red"))
