# CODE FOR PAQUETTE & HARGREAVES, ECOLOGY LETTERS 2021 ###################.
# Biotic interactions contribute more often to species’ warm vs. cool range limits #################.

#Cleaning, analysis, and figures
library(lattice)
library(lme4) #for mixed models
library(car) #for Anova tables
library(emmeans) #for lsmeans
library(visreg)

citation(package='lme4')
packageVersion('lme4')
citation(package='emmeans')
packageVersion('emmeans')
citation(package='visreg')
packageVersion('visreg')

setwd("~/Dropbox/A shared folders/2019 Alexandra Paquette/") 

#DATA EXPLANATION #########################
#Overview: data come from Cahill et al 2014 doi:10.1111/jbi.12231 + literature search up to end of 2019, which expanded Cahill data to include cool range limits, studies published after Cahill literature search, studies missed by Cahill literature search. See detailed notes on data compilation and extraction in Supplementary Material for paper
#Files: Cool data are in 1 csv, Warm data are in 2nd csv (merged in code) 
#Structure: 1 row per study x species x range limit (elevation or latitude, sometimes separated by continent/ocean) x factor assessed

#Note: Cahill includes some 'resurvey' type studies, eg Feely et al, Lenoir et al. Cahill interpreted a shift in warm limit as supporting the driver Temperature, even though there were no data about regional climate change at that time (did it actually get hotter?). We decided not to include resurveys as they've already been reviewed elsewhere

#Notes re changes & corrections to data during project:
# checked all studies where Genus=='Multiple' from Cahill, extracted data for single spp when possible, added data to cool data base when appropriate, and added notes to exclude some of Cahill's original lines of data or studies (detailed below and in excel):
# --Ferreira (cool) - took out precip and 'other' lines, as study only related these to diversity patterns, not distribution patterns 
# --Loidi - changed Temp supported to N from Cahill as only supported at cool limit, and added study to cool RL data
# --Morin et al 2007 - AH separated by species and added to cool data
# --Normand et al 2010 - EXCLUDE: can't untangle species and no factor supported for all species so cant classify as 'Y or N'
# --Bassler et al 2010 - study looks at entire range of species (doesnt distinguish warm vs cool RLs). Cahill included 1 summary line for all species for each variable tested in model & had support=N for everything other than temperature - while Temp was most important, other factors were important too, so those 'N' designations are incorrect. Use only the 18 species in Table 3 for whom individual GLMs were run. Of these only some species had their warm RL in survey, and only some had cool RL in survey. Still a bit odd as study doesnt distinguish between range limits but allow since was included in original Cahill analysis.
#--Attore - exclude. study doesnt break down range by range edge and factors are both supported (significant correlation with range) and not supported (species way underfills its range even on small island)
#--Bowman - excluded some unsupported data from warm RL, didnt add to cool as data not species-specific enough
#--Cadena - split results by study area and added to cool
#--Ebert precip result excluded (couldnt find support in paper)
#--Fields 2012 - exclude. native mussel does have low heat tolerance but this is never explicitly linked/compared to its warm RL 
#--Guisan 1998 - added competition to lower limit (no data on cold limit)
#--Hersteinsson 1992 - changed one 'Biotic Other' result (dont know which driver that was supposed to be and none signficant) 
#--Iszkulo 2009 - took out latitudinal data as no clear latitudinal range limit, and took out one 'Abiotic Other factor' as only precipitation supported (no cool limit data)
#--Jacobsen 1991 - added salinity 
#--Kristiansen 2004 exclude study - dont test temps low enough to be relevant to northern limit and dont relate lethal thermal limit to actual southern limit (just to eastern warm limit)
#--Li etal 2007 - exclude warm data as heat limit never explicitly related to location of warm RL.  
#--Marino - altered concl re precip & added to cool
#--Middleton - added Dispersal & Competition & cool data, changed study type from 'model' to "observation'
#--Mooney 1965 - added Temp (no cool RL data)
#--Morueta - changed conclusion re topography
#--Para-Olea 2005 - combined temp & preip -> Climate and added to cool data
#--Portner 2001 - exclude 'other' as cant figure out what would be
#--Swenson 2008 - added other oak sp, changed conclusion re precip & Other, added to cool data
#--Tannerfeld 2002 - exclude 'Other' cant figure out what would be in paper
#--Tolley 2010 - exclude elevation (too hard to decipher & related to lat) and other Abiotic (not sure what that was)
#--Wang 2993 - changed conclusion for warm RL & added to cool data
#--Xu 2001 - EXCLUDE study only maps range in China, but northern RL is in russia and south RL is in korea. 


#DATA RL DRIVERS/FACTORS (skip to 'Export/call in allRL' line634 if dont need to redo data) #################

#Data cool range limits drivers ------------------------------------------
cool.orig <- read.csv("Data/Paquette&Hargreaves data cool limits Dryad.csv", 
                      skipNul = TRUE, fileEncoding="latin1", 
                      header = TRUE, stringsAsFactors =T) 
cool.orig <- cool.orig[cool.orig$Genus!='',]; cool.orig <- droplevels(cool.orig) #drop empty rows
names(cool.orig) 
#check 'Authors' column - sometimes has problems when called in and has to be manually renamed
names(cool.orig)[names(cool.orig) == "ï..Authors"] <- "Authors"
names(cool.orig)
summary(cool.orig)
cool.orig[duplicated(cool.orig),] #should be none 

#check have unique dois for each study
nlevels(cool.orig$Title) #250 unique titles # DIFF ##############
nlevels(cool.orig$doi) #250 unique dois (should match # of titles) # DIFF ################
#if needed, find cases where same title has >1 doi etc
# test <- unique(cool.orig[, c('Authors', 'Year.pub', 'Title', 'doi')]); 
# test <- droplevels(test)
# dim(test)
# summary(test)

#__explanation of variables (full explanation in excel file) --------------------------
names(cool.orig)
#Authors / Year.pub / doi / Title / Abstract = information about the publication data are extracted from
#Genus = genus of study organism. Listed as 'Multiple' if data presented across genera & cant be separated
#Species = species epithet. 'sp.' means species is unique but unnamed, 'spp.' means data are merged across >1 species
#Number of species = mostly 1, but if Genus=Multiple or Species=spp., denotes how many species included in data
#General.env: F=freshwater, M=marine, T=terrestrial (changed to Environment in the final data)
#High.taxon = species classified into Plant/Algae; Animal; Fungus; Protist (Note Protist category not consistent, some algae technically protists)
#Taxon = more specific classification (e.g. bird) - altered for Cahill data as level inconsistent. 'Mixed' if data point could not be separated even by broad taxon (eg if study surveyed all plants, including herbs, trees, ferns, mosses etc)
#Taxon.description.if.mixed = more details if Taxon=Mixed
#Exotic = Y if the species was not native where it was studied, N=native (lots of NAs as Cahill didn't collect this)
#RL.Lat.Elev = is range limit latitudinal (cool=polar, warm = equatorward) or elevational
#RL = cool or warm
#Method.category = study method, using categories from Cahill et al. 'Field experiment', 'Lab experiment' (these were almost always thermal physiology studies that then compared a heat/cold threshold to mean/min/max temperatures at a given RL), 'Observation' or 'Model'
#Summary.method = brief summary of study method
#Driver.examined = specific potential range-limiting ecological factor. eg temperature, competition, precipitation. used 'Climate' when effects of Temp and Precip cant be disentangled from each other or other climate factors in results (often the case with SDMs). If study looked at >1 type of factor within a category these are rolled into a single data point for that driver. eg if tested Mean annual temp + max_temp + min temp we have 1 line for Temperature (support would be Y (supported) if any of the three were supported)
#Driver.notes = more details on driver, especially if Driver = Other
#Biotic.Abiotic = was ecological driver examined  biotic or abiotic (changed to Driver.type in the final dataset). D = dispersal (separated from 'Biotic' as it is an intrinsic feature of the focal species, rather than a feature of its environment)
#Driver.supported = Y if supported by the study, N if not. Cahill et al seem to list as Y if even partially supported (eg for 30% of species in multi-species study) - in those cases with mixed support across species, we separated results by species or excluded
#Results.notes = notes re results
#Data.location = where in study did we find data for the conclusion in 'Driver.supported'
#Data.quality = our attempt to note quality of data. did not do for all studies so dont use in analyses
#Country..state = country and state/province where data come from / range limit
#Continent = continent of range limit
#Search = how found the paper
#person.did.initial.cull = person who deemed whether study had useable data
#person.did.extraction = person who extracted data from study
#person.who double checked = person who verified the data - done for subset of studies/data points
#Study.inCahill2014 = was STUDY (not data point) in Cahill2014. Options are Y=yes, N.butrelevant (ie contains data on warm RLs, but either Cahill missed it or was published after their search), or N.notwarm (study doesnt contain any data on warm RLs)
#keep? = column used to exclude data if later deemed it no good for analysis

#additional variables in warm data only:
#NB: got rid of Cahill's 'Study class' as didnt find it a good measure of study quality and was hard to apply consistently
#Data.inCahill2014 = was this ROW of data in Cahill etal 2014?  Y if row was in Cahill unmodified, 'N' otherwise (eg study not in Cahill, or conclusions altered from Cahill, eg if disagree with their results re whether driver supported or not, or if obtained data for single species where Cahill only had grouped data, etc)
#Data.modified.from.Cahill. = 'NA' if study/data not in Cahill. If Data were in Cahill, answers can be: 'N - not checked' if we have not reviewed Cahill's conclusions; 'N' if we have reviewed Cahill conclusions and agree; 'exclude xyz' if data row was in Cahill but we dont want to include in ours as we were able to separate by species, disagree w Cahill conclusion, can't figure out where they got data (esp when Driver=Other) etc; Y if line is new and modified from Cahill

#__data cleaning cool RLs------------------------------
cool <- cool.orig #make new data base for cleaning 

#get rid of busy columns not included in the analysis
cool$Abstract <- cool$Search <- cool$Location.detailed <- NULL
cool$person.did.initial.cull <- cool$Data.location <- cool$Summary.method <- NULL
summary(cool)

#create column for species name
table(cool$Taxa.description.if.mixed[cool$Taxon=='Mixed']) #only 2 papers with multiple genera that cant be separated
cool$Sp.name <- as.factor(ifelse(cool$Taxon=='Mixed', as.character(cool$Taxa.description.if.mixed), 
                          ifelse(cool$Taxon!='Mixed' & cool$Genus=='Multiple', paste(cool$Taxon, cool$Species, sep='.'), 
                                 paste(cool$Genus, cool$Species, sep='.'))))
head(cool)
test <- cool[cool$Genus=='Multiple',]; test <- droplevels(test) #how much of data is grouped across genera?
test #not much (9 lines from 4 studies)

#create better column for general environment
levels(cool$General.env) 
cool$Environment <- as.factor(ifelse(cool$General.env=='T', 'terrestrial',
                          ifelse(cool$General.env=='M', 'marine', 'freshwater')))
summary(cool)
table(cool$Environment, cool$General.env) #make sure Gen env and Environment line up properly

#make sure taxa are assigned to Environment & Higher.taxon that makes sense
table(cool$Taxon, cool$Environment)
cool[cool$Taxon=='Mollusc' & cool$Environment=='terrestrial',] # Deroceras.panormitanum = land slug so correct
cool[cool$Taxon=='Reptile',] #freshwater reptile = turtles, marine reptile = sea snake so correct
#check taxa classificaiton 
table(cool$Taxon, cool$High.taxon)

#create better variable name for factor type (biotic abiotic) and add new category for dispersal limitation
table(cool$Driver.examined, cool$Biotic.Abiotic)
cool$Driver.type <- as.factor(ifelse(cool$Driver.examined=='Dispersal', 'dispersal',
                                 ifelse(cool$Biotic.Abiotic=='B', 'biotic', 'abiotic')))
table(cool$Biotic.Abiotic, cool$Driver.type) #check
table(cool$Driver.examined, cool$Driver.type) 

#create quantitative column for Driver supported (makes binomial analyses on this response less confusing)
cool$Driver.supported01 <- ifelse(cool$Driver.supported=='Y', 1, 0)
table(cool$Driver.supported, cool$Driver.supported01)

#check no blanks in doi (so can use as random effect for study)
cool[cool$doi=='',] #make sure no blanks
# cool.orig[cool.orig$doi=='',] #check if any exist
# cool$xdoi <- cool$doi
# cool$doi <- as.factor(ifelse(cool$Authors=='' & cool$Year.pub==1, 'X', as.character(cool$xdoi)))
# levels(cool$doi) #at least one blank
# cool$xdoi <- NULL

#check remaining columns
levels(cool$Continent) #check spelling
cool[cool$Continent=='',] #should be none
summary(cool)
cool[is.na(cool$Driver.supported),c('Title', 'Genus', 'Species', "Driver.examined", 'Driver.notes', 'Driver.supported', 'Result.notes')] #one abiotic NA = line of data where can't see map properly to extract answer
cool <- cool[!is.na(cool$Driver.supported),] #get rid of NAs
cool$x <- NULL
summary(cool)

#create final data frame for analyses
names(cool)
cool[duplicated(cool),] #should be none 
levels(cool$keep.)
coolRLs <- cool[!grepl('exclude',cool$keep.), #get rid of rows of data we decided are invalid
                c('Authors', 'Year.pub', 'Title', 'doi', 
                   'Genus', 'Species', 'Sp.name', 'High.taxon', 'Taxon', 
                   'RL', 'RL.Lat.Elev', 'Environment', 'Exotic', 'Country..state.', 'Continent',
                   'Driver.examined','Driver.notes','Driver.type','Driver.supported', 'Driver.supported01',
                   'Method.category',
                   'Study.inCahill2014','person.did.extraction','person.who.double.checked')]
coolRLs <- droplevels(coolRLs); summary(coolRLs)
colnames(coolRLs) <- c('Authors', 'Year.pub', 'Title', 'doi', 
                   'Genus', 'Species', 'Sp.name', 'High.taxon', 'Taxon', 
                   'RL', 'RL.Lat.Elev', 'Environment', 'Exotic', 'Country.state', 'Continent',
                   'Driver.examined','Driver.notes','Driver.type','Driver.supported', 'Driver.supported01',
                   'Method',
                   'Study.inCahill2014', 'data.extractor','double.checker')
summary(coolRLs)

dim(coolRLs) #1183 rows  
coolRLs[duplicated(coolRLs),] #should be none 
coolRLs <- droplevels(coolRLs) #should make no diff but just in case
nlevels(coolRLs$doi) #248 studies 
nlevels(coolRLs$Sp.name) #467 taxa (includes some multiple species)
nlevels(coolRLs$Genus) #335 (includes genus 'Multiple') #DIFF ##############

#do have any studies or taxa where only study dispersal?
test <- coolRLs[coolRLs$Driver.type!='dispersal',]; test <- droplevels(test)
nlevels(coolRLs$doi); nlevels(test$doi) #1 study - transplant experiment that didnt ever test or even speculate about ecological causes that can be defined as abiotic or biotic. eg say 'not enough microsites'
test0 <- aggregate(coolRLs[c('Driver.examined')], by=coolRLs[c('Sp.name','doi','Title','Authors','Driver.type','RL')], function(x) length(unique(x)))
head(test0)
test0$dispersal <- ifelse(test0$Driver.type=='dispersal',1,0) 
test0$abiotic <- ifelse(test0$Driver.type=='abiotic',1,0) 
test0$biotic <- ifelse(test0$Driver.type=='biotic',1,0) 
summary(test0)
test1 <- aggregate(test0[c('dispersal','abiotic','biotic')], by=test0[c('Sp.name','Title','Authors','doi','RL')], function(x) sum(x))
test1[test1$dispersal>0 & test1$abiotic==0 & test1$biotic==0,]
#were some tests where they never identified a factor that limited the range. 
#Kimbal et al 
#PEarson 1931 some sppRL

#Data warm range limits drivers ------------------------------------

#data are of 3 types: 
#1) Cahill et al 2014's Appendix, reformatted to fit APs more easily analysable format (1 row per driver), with some corrections made (eg Taxon made consistent); of these some results have been changed eg if disagreed with Cahill conclusions re driver support (did NOT check all studies), or if Cahill included data grouped by species but was possible to get  data for individual species 
#2) data that should have been in Cahill original study but that were missed; 
#3) data published since Cahill literature search
warm.orig <- read.csv("Data/Paquette&Hargreaves data warm limits Dryad.csv", skipNul=T, stringsAsFactors=T) 
summary(warm.orig)
warm.orig <- Filter(function(x)!all(is.na(x)), warm.orig) #if get a bunch of NA columns
names(warm.orig)
warm.orig[duplicated(warm.orig),] #should be none 

#check have unique dois for each study
nlevels(warm.orig$Title) #224 unique titles
nlevels(warm.orig$doi) #224 unique dois
nlevels(warm.orig$Authors) #220 unique Author lists
#double check cases where have repeated dois or Autho lists for title
test <- unique(warm.orig[, c('Authors', 'Year.pub', 'Title', 'doi')]); test <- droplevels(test)
dim(test) 
summary(test) #4 cases where same author list has >Title/doi, but no cases where title or doi repeated so must be independent papers


#__data cleaning warm limits --------------------------- ---
warm <- warm.orig #make new data base for cleaning 

#get rid of some of the long columns to make data base easier to read
warm$Abstract <- warm$Data.location <- NULL
warm$Summary.method <- warm$person.did.initial.cull <- NULL 
summary(warm)

#create column for species name
table(warm$Taxa.description.if.mixed[warm$Taxon=='Mixed'])
warm$Sp.name <- as.factor(ifelse(warm$Taxon=='Mixed', as.character(warm$Taxa.description.if.mixed), 
                          ifelse(warm$Taxon!='Mixed' & warm$Genus=='Multiple', paste(warm$Taxon, warm$Species, sep='.'), 
                                 paste(warm$Genus, warm$Species, sep='.'))))
head(warm)
test <- warm[warm$Genus=='Multiple',]; test <- droplevels(test); test
table(test$Taxon, test$Sp.name)

#create better column for general environment
warm$Environment <- as.factor(ifelse(warm$General.env=='T', 'terrestrial',
                          ifelse(warm$General.env=='M', 'marine', 'freshwater')))
summary(warm)
table(warm$Environment, warm$General.env) #make sure Gen env and Environment line up properly
warm$General.env <- NULL

#make sure taxa are assigned to Environment & Higher.taxon that makes sense
#NB: AH added 'evergreenAngiosperm' as wasnt sure if Cahill's Evergreen was being used to refer to gymnosperms or non-deciduous trees
table(warm$Taxon, warm$Environment) 
warm[warm$Taxon=='Mixed',] #these are multi-species papers
#check taxa classification - using different categories here compared to cool limits (which had fungi instead of protist)
table(warm$Taxon, warm$High.taxon) #keep as is for now - can resolve later once combined with cool limits data

#create better variable name for factor type (biotic abiotic) and add new category for dispersal limitation
table(warm$Driver.examined, warm$Biotic.Abiotic)
warm$Driver.type <- as.factor(ifelse(warm$Driver.examined=='Dispersal', 'dispersal',
                                 ifelse(warm$Biotic.Abiotic=='B', 'biotic', 'abiotic')))
table(warm$Biotic.Abiotic, warm$Driver.type) #check
table(warm$Driver.examined, warm$Driver.type) 

#create quantitative column for Driver supported, rename 'Driver.supp', get rid of NAs/unks
names(warm)
summary(warm$Driver.supported)
warm[warm$Driver.supported=='unk' | is.na(warm$Driver.supported),c(14,15,17:21,28)]
warm <- warm[warm$Driver.supported!='unk' & !is.na(warm$Driver.supported),]; warm <- droplevels(warm)
dim(warm); dim(warm.orig) #gets rid of ~8 rows
warm$Driver.supported01 <- ifelse(warm$Driver.supported=='Y', 1, 0)
table(warm$Driver.supported, warm$Driver.supported01)

#check no blanks in doi (so can use as random effect for study)
warm[warm$doi=='',] #check for  blanks

#check remaining columns
levels(warm$Continent)
levels(warm$Study.inCahill2014) #should all be Y, N, or N.butrelevant  
summary(warm) #missing lots of rows for Exotic (>400) as Cahill didnt record that

#create final data frame for analyses
names(warm)
summary(warm) #of 1175 rows of data, only 247 have not been doubled checked 
warmRLs <- warm[,c('Authors', 'Year.pub', 'Title', 'doi', 
                   'Genus', 'Species', 'Sp.name', 'High.taxon', 'Taxon', 
                   'RL', 'RL.Lat.Elev', 'Environment', 'Exotic', 'Country..state.', 'Continent',
                   'Driver.examined','Driver.notes','Driver.type','Driver.supported', 'Driver.supported01',
                   'Method.category', 
                   'Study.inCahill2014','Data.inCahill2014', 'Data.modified.from.Cahill.',
                   'person.did.extraction', 'person.who.double.checked')]
summary(warmRLs)
colnames(warmRLs) <- c('Authors', 'Year.pub', 'Title', 'doi', 
                   'Genus', 'Species', 'Sp.name', 'High.taxon', 'Taxon', 
                   'RL', 'RL.Lat.Elev', 'Environment', 'Invasive', 'Country.state', 'Continent',
                   'Driver.examined','Driver.notes','Driver.type','Driver.supported', 'Driver.supported01',
                   'Method', 
                   'Study.inCahill2014','Data.inCahill2014', 'Data.modified.from.Cahill',
                   'data.extractor', 'double.checker')
summary(warmRLs)

#double check dont have any studies or taxa where only study dispersal
test <- warmRLs[coolRLs$Driver.type!='dispersal',]; test <- droplevels(test)
nlevels(warmRLs$doi); nlevels(test$doi) #lose some studies - there were some transplant experiments that didnt ever test or even speculate about ecological causes that can be defined as abiotic or biotic. eg say 'not enough microsites'
table(test$Driver.examined, test$Driver.type)
test0 <- aggregate(warmRLs[c('Driver.examined')], by=warmRLs[c('Sp.name','doi','Title','Driver.type')], function(x) length(unique(x)))
head(test0)
test0$dispersal <- ifelse(test0$Driver.type=='dispersal',1,0) 
test0$abiotic <- ifelse(test0$Driver.type=='abiotic',1,0) 
test0$biotic <- ifelse(test0$Driver.type=='biotic',1,0) 
summary(test0)
test1 <- aggregate(test0[c('dispersal','abiotic','biotic')], by=test0[c('Sp.name','Title','doi')], function(x) sum(x))
test1[test1$dispersal>0 & test1$abiotic==0 & test1$biotic==0,] #just 1 study x species combos (then why does it look like we lose 10 on L398?)
#were some tests where they never identified a factor that limited the range. 
#Gauthier et al Ecotype differentiation 


#Merge cool & warm limit driver data---------------------------
#create columns to match warm data (just so merge successfully)
coolRLs$Data.modified.from.Cahill <- as.factor(NA)
coolRLs$Data.inCahill2014 <- as.factor('N')

names(coolRLs)
cool.formerge <- coolRLs[,c('Authors', 'Year.pub', 'Title', 'doi', 
                            'Genus', 'Species', 'Sp.name', 'High.taxon', 'Taxon', 
                            'Environment', 'RL', 'RL.Lat.Elev', 'Continent', 'Country.state',   
                            'Driver.examined','Driver.notes','Driver.type','Driver.supported', 'Driver.supported01',
                            'Method', 'Study.inCahill2014', 'Data.inCahill2014','Data.modified.from.Cahill',
                            'data.extractor','double.checker')]
head(cool.formerge)
summary(cool.formerge)
cool.formerge[duplicated(cool.formerge),] #should be none
levels(cool.formerge$Authors)

names(warmRLs)
warm.formerge <- warmRLs[,c('Authors', 'Year.pub', 'Title', 'doi', 
                            'Genus', 'Species', 'Sp.name', 'High.taxon', 'Taxon', 
                            'Environment', 'RL.Lat.Elev', 'Continent', 'Country.state', 'RL',  
                            'Driver.examined','Driver.notes','Driver.type','Driver.supported', 'Driver.supported01',
                            'Method', 'Study.inCahill2014', 'Data.inCahill2014','Data.modified.from.Cahill',
                            'data.extractor','double.checker')]
warm.formerge[duplicated(warm.formerge),] #should be none
#merge
allRLdata <- rbind(cool.formerge, warm.formerge)
summary(allRLdata)
allRLdata[allRLdata$High.taxon=='',] #should be none

#check doi vs Title numbers again - possible that used diff spellings of doi in each database (eg include 'doi' or not)
test <- aggregate(allRLdata[c('doi')], by=allRLdata[c('Title','Authors')], function(x) length(unique(x)))
summary(test) #should only be 1 doi per Title x Author contribution (ie no title taht appears >1 and doi max=1)
# test[test$doi > 1,]

#taxonomic details - make sure consistent in assigning 'Taxon' categories (inherited some from Cahill)
table(allRLdata$High.taxon, allRLdata$RL)
#check out the Protist
allRLdata[allRLdata$High.taxon=='Protist',] #protist = protozoan
table(allRLdata$Taxon, allRLdata$RL)
#fix Taxon (note some ferns might still be included as herbaceaous but dont think will do analysis by taxon so fine)
#fixed weird categories in Cahill too: ie cactus = succulent, lizard = reptile, grass = herbaceaous
allRLdata$Taxon <- as.factor(
                    ifelse(allRLdata$Taxon=='EvergreenAngiosperm', 'Evergreen Angiosperm', as.character(allRLdata$Taxon)))
table(allRLdata$Taxon, allRLdata$RL)
#what are the Plant Taxa?
test <- allRLdata[allRLdata$High.taxon=='Plant/Algae',]; test <- droplevels(test)
table(test$Taxon, test$RL)  

#__check mixed taxa studies where couldn't separate results even by Genus -------------------- ---
mix <- allRLdata[allRLdata$Genus=='Multiple',]; mix <- droplevels(mix)
table(mix$RL) #most studies with taxa all mixed together are from warm database
mix[,c('Authors','doi','Genus','Species','Sp.name','High.taxon','Taxon','RL','Method','Driver.examined','Data.inCahill2014')]
levels(mix$doi) #19 studies
levels(mix$Authors)
# #--Bassler etal 2010: EXCLUDE multi-spp data = iffy & was able to extract more specific data for subset of individual species (see notes at top)
# #--Bischoff-Basmann: EXCLUDE multi-spp data (Cahill), USE ungrouped data extracted by species
# #--Bruelheidi 2003 = EXCLUDE multi-spp data (Cahill), USE ungrouped data extracted by species
# #--Butterwick etal - EXCLUDE STUDY. thermal physiology only, only speculation re range in discussion, not enough detail to unpack by species
# #--Calosi etal 2010 - cant unpack by species so add to cold limit data
# #--Canham & Thomas 2010 #EXCLUDE: trees in US. present data on climate range but cant find any data on geographic range 
# #--Cunningham, HR; Rissler, LJ; Buckley, LB; Urban, MC
# #--Feely 2012: EXCLUDE = resurvey study, which have already been synthesized elsewhere.  ---
# #--Gignac etal 2000 - sphagnum species quantified and analysed together - cant pull apart but should be included in cool data
# #--Hargrove & Rotenberry 2011 - EXCLUDE multi-spp data (Cahill) AH separated by species and added to cool data
# #--Hsieh etal 2009 - EXCLUDE multispecies data - extracted by species (important as importance really varies)
# #--Kimura 2004 - EXCLUDE measure thermal tolerance and correlate to location of range edge but dont relate absolute value of tolerance to actual climate at range edge (just species form farther north more cold tolerant)
# #--Larsen 2012 - EXCLUDE multi-spp data (not all species had both RLs in summary so # spp = wrong, and not all warm limits responded the same way so universal Temp = supported also not the case). USE ungrouped data instead (added to cool data too)
# #--Lenoir etal 2010: EXCLUDE - resurvey data like Feely. 
# #Loehel 1998 #boreal tree but include evergreen & deciduous. AH adjusted data now = good
# #loidi etal 2010 = biogeography heathlands - cant break down by species
# #Monahan 2009 = EXCLUDE grouped data as only presented at whole-range level (cant separate warm vs cool limits). USE data for the 2 species shown in Fig 1
# #Morin et al 2007 - should be in cool RL data too... AH separated by species and added to cool data now = good
# #Normand etal 2009 - European plants. but only have warm range limit data
# #Tittensor etal 2009

#check 'Sp.name' for studies with multiple species
test <- allRLdata[allRLdata$Taxon=='Mixed',]; test <- droplevels(test)
summary(test) # all plant studies
summary(allRLdata)

#__get rid of Cahill data rows where couldnt verify or disagreed with conclusion  ---------
dim(allRLdata) #2358 rows
allRLdata[grepl('exclude', allRLdata$Data.modified.from.Cahill),]
allRL <- allRLdata[!grepl('exclude', allRLdata$Data.modified.from.Cahill),] #get rid of Cahill data that cant verify
dim(allRLdata); dim(allRL) #gets rid of >200 lines 
dim(allRL) #2125 lines 
allRL$Title <- NULL #dont need this anymore
summary(allRL)
table(allRL$Data.modified.from.Cahill, allRL$RL)
allRL$Data.modified.from.Cahill <- NULL
allRL <- droplevels(allRL)
#check for duplication errors
allRL[duplicated(allRL),] #should be none

#look again at studies with multiple species combined
mix <- allRL[allRL$Genus=='Multiple',]; mix <- droplevels(mix)
dim(mix) #22 lines of data
nlevels(mix$doi) #from 5 studies with data mixed across genera, 
summary(mix) #most = warm RL but not as lopsided since ungrouped lots of the Cahill data & added relevant studies to cool data set

#check Drivers
summary(allRL)
table(allRL$Driver.examined, allRL$Driver.type) #Climate for when Temp & Precip effects could not be separated
test <- allRL[allRL$Driver.examined=='Other' & allRL$Driver.notes=='',]; test <- droplevels(test)
table(test$RL, test$Driver.type) #no more unexplained 'Other' - ALH filled them all Nov/Dec 2020 
#create unified 'Other' description? (not sure we'll ever use it)
test <- allRL[allRL$Driver.examined=='Other',]; test <- droplevels(test)
table(test$Driver.notes, test$RL)
summary(allRL[,c('Driver.examined','Driver.notes')])
allRL$Driver.descrip <- as.factor(
        ifelse(allRL$Driver.examined=='Other' & grepl('abiotic habitat', allRL$Driver.notes, ignore.case=T), 'abiotic habitat',
        ifelse(allRL$Driver.examined=='Other' & grepl('fire', allRL$Driver.notes, ignore.case=T), 'fire',
        ifelse(allRL$Driver.examined=='Other' & grepl('salinity', allRL$Driver.notes, ignore.case=T), 'salinity',
        ifelse(allRL$Driver.examined=='Other' & grepl('oxygen', allRL$Driver.notes, ignore.case=T), 'oxygen',
        ifelse(allRL$Driver.examined=='Other' & grepl('topography', allRL$Driver.notes, ignore.case=T), 'topography',
        ifelse(allRL$Driver.examined=='Other' & grepl('nutrient', allRL$Driver.notes, ignore.case=T), 'nutrients',
        ifelse(allRL$Driver.examined=='Other' & grepl('cloud cover', allRL$Driver.notes, ignore.case=T), 'cloud cover',
        ifelse(allRL$Driver.examined=='Other' & grepl('snow', allRL$Driver.notes, ignore.case=T), 'snow',
        ifelse(allRL$Driver.notes=='Day length' | grepl('photoperiod', allRL$Driver.notes, ignore.case=T), 'photoperiod',
        ifelse(allRL$Driver.notes=='PAR' | grepl('light', allRL$Driver.notes, ignore.case=T) | grepl('Solar', allRL$Driver.notes), 'light',
        ifelse(allRL$Driver.notes=='pollination' | grepl('mutual', allRL$Driver.notes) | grepl('facilitat', allRL$Driver.notes, ignore.case=T), 'lack mutualists',
        ifelse(allRL$Driver.examined=='Other', as.character(allRL$Driver.notes), '')))))))))))))
table(allRL$Driver.descrip, allRL$RL)
#check data on mutualists
test <- allRL[allRL$Driver.descrip=='lack mutualists',]; test <- droplevels(test)
summary(test)
dim(test) #21 rows data total, 12 cool, 9 warm 
nlevels(test$doi) #14 studies - not bad
nlevels(test$Sp.name) #15 species
test
allRL[grepl('facilitat', allRL$Driver.notes, ignore.case =T),]

#code for checking whether certain studies included  
# allRL[grepl('Pauw', allRL$Authors),] #not there (in this case b/c there was not an actual RL in the study)
# allRL[grepl('Afkhami', allRL$Authors),] #not there (in this case b/c range limit was longitudinal)

#make simpler Driver column that also clusters abiotic & biotic
table(allRL$Driver.examined)
allRL$Driver <- as.factor(ifelse(allRL$Driver.examined=='Biogenic habitat', 'B.Habitat',
                          ifelse(allRL$Driver.examined=='Climate', 'A.Climate',
                          ifelse(allRL$Driver.examined=='Competition', 'B.Comp',
                          ifelse(allRL$Driver.examined=='Host/Food availability', 'B.HostFood',
                          ifelse(allRL$Driver.examined=='Parasitism/Disease', 'B.ParasitDis',
                          ifelse(allRL$Driver.examined=='Precipitation', 'A.Precip',
                          ifelse(allRL$Driver.examined=='Predation/Herbivory', 'B.PredHerb',
                          ifelse(allRL$Driver.examined=='Soil', 'A.Soil',
                          ifelse(allRL$Driver.examined=='Temperature', 'A.Temp', 
                          ifelse(allRL$Driver.type=='abiotic' & allRL$Driver.examined=='Other', 'A.Other', 
                          ifelse(allRL$Driver.type=='biotic' & allRL$Driver.examined=='Other', 'B.Other',
                                 'Dispersal'))))))))))))
table(allRL$Driver.examined, allRL$Driver)

#Make simplified Method columns that lists only the most powerful method used (field exp > lab exp > obs > model) (in some cases observational studies were more convincing than lab studies but only use this to identify field experiments so this is fine)
summary(allRL$Method)
allRL$Method.all <- allRL$Method
allRL$Method <- as.factor(ifelse(grepl('field', allRL$Method.all, ignore.case=T), 'Field experiment', 
                          ifelse(grepl('lab', allRL$Method.all, ignore.case=T), 'Lab experiment', 
                          ifelse(allRL$Method.all=='Model + Observation', 'Observation', 
                          as.character(allRL$Method.all)))))
table(allRL$Method, allRL$Method.all)


#__create random factors ----------------------- 
levels(allRL$Country.state)
levels(allRL$RL.Lat.Elev)
allRL$RL.Lat.Elevshort <- as.factor(ifelse(allRL$RL.Lat.Elev=='Elevational', 'Elev', 'Lat'))
levels(allRL$Continent)
allRL$Continentshort <- as.factor(ifelse(grepl(';', allRL$Continent), 'Multiple', as.character(allRL$Continent)))
table(allRL$Continent, allRL$Continentshort)
levels(allRL$Country.state)
#random factor for species-RL-location (use to count number of species-range limits)
allRL$sppRLloc <- as.factor(paste(allRL$Sp.name, allRL$Continentshort, allRL$RL.Lat.Elevshort, allRL$RL, sep='.'))
test <- aggregate(allRL[c('doi')], by=allRL[c('sppRLloc')], function(x) length(unique(x)))
hist(test$doi) #most range limits just studied once but some studied up to 5x
#random factor for RL-location (use when create support by factor to not merge across continents)
allRL$RLloc <- as.factor(paste(allRL$Continentshort, allRL$RL.Lat.Elevshort, allRL$RL, sep='.'))
test <- aggregate(allRL[c('doi')], by=allRL[c('RLloc')], function(x) length(unique(x)))
test
hist(test$doi)
allRL$RL.Lat.Elevshort <- allRL$Continentshort <- NULL
summary(allRL)
dim(allRL) #2125

#Export / call in allRL data on range limit drivers ------------------
#write.csv(allRL, 'Data/allRL processed.csv')
allRL <- read.csv("Data/allRL processed.csv", skipNul=T, stringsAsFactors=T)
allRL$X <- NULL
summary(allRL)
dim(allRL) #2125
nlevels(allRL$doi) #340 # DIFF #################

#__make df tallying support by factor type ------
#the 'was driver supported' analysis is good in that it only tests a factor as much as it was tested. BUT, if a RL was entirely driven by temperature and a study tested Temp, precip, and soil, the result would be <50% support for abiotic factors.  So that's not quite right either. So aggregate: if study looked at abiotic factor(s) and at least 1 abiotic factor was supported, Abiotic=Y, but if none supported Abiotic=N. Same for studies that looked at biotic factors. 
#adds up support for each study x sp x RL combo. if no abiotic driver tested, wont be a line. if at least 1 abiotic driver tested and none were supported, support=0
summary(allRL)
names(allRL)
#also want to keep track of Methods if want to subset only spp including experimental data
allRL$dataFexp <- ifelse(allRL$Method=='Field experiment', 1, 0)
#code to do in 1 step wouldnt work - do by brute force instead
allRLsum1 <- aggregate(allRL[c('Driver.supported01','dataFexp')], 
                      by=allRL[c('doi','Sp.name','Genus','Species','Continent','Environment','RL.Lat.Elev','RL','RLloc','Driver.type')], 
                      function(x) sum(x))
head(allRLsum1) #added columns show up 
colnames(allRLsum1)[colnames(allRLsum1)=='Driver.supported01'] <- 'sumDriver.support'; head(allRLsum1) 
dim(allRLsum1) #1447
#now get n drivers tested
n <- aggregate(allRL[c('Driver.supported01')], 
                      by=allRL[c('doi','Sp.name','Genus','Species','Continent','Environment','RL.Lat.Elev','RL','RLloc','Driver.type')], 
                      function(x) length(x))
head(n)
dim(n) #1447
colnames(n)[colnames(n)=='Driver.supported01'] <- 'nDrivers.tested'; head(n)
allRLsum <- merge(allRLsum1,n); head(allRLsum)
dim(allRLsum) 
#check out how many drivers tested per sp x RL x study
hist(allRLsum$nDrivers.tested) #max of 5 factors tested at given sp x RL per study (good)
bwplot(nDrivers.tested ~ Driver.type, data=allRLsum)
#make new variable to assess support by driver type: if Driver.supported >0, ==yes
allRLsum$Drivertype.supported01 <- ifelse(allRLsum$sumDriver.support>0, 1, 0)
summary(allRLsum)
bwplot(Driver.supported01 ~ Driver.type | RL, data=allRL)
bwplot(Drivertype.supported01 ~ Driver.type | RL, data=allRLsum) #support for abiotic goes way up
#make new variable to denote whether data includes field experiments if Driver.supported >0, ==yes
allRLsum$inclFexp <- as.factor(ifelse(allRLsum$dataFexp>0, 'Y', 'N'))
summary(allRLsum)
allRLsum$dataFexp <- NULL
allRL$dataFexp <- NULL
summary(allRL)
nlevels(allRL$doi) #340 # DIFF #################

#__make var denoting species RLs w data on both biotic & abiotic drivers --------- --
#species for which both were tested provide strongest tests of their relative importance
#should be only studies that assess both abiotic and biotic drivers, or only species for which have data on both? Think the 2nd - expect that most studies will look at both at same time, but if had species where one group of authors looked at only abio factors and another looked at only biotic, would still want to count it in these analyses. 
summary(allRLsum) #should already be only 1 line per sp x RL x driver type    
allRLsum$dataBio <- ifelse(allRLsum$Driver.type=='biotic', 1, 0) 
allRLsum$dataAbio <- ifelse(allRLsum$Driver.type=='abiotic', 1, 0)
    
#aggregate across Driver types (biotic & abiotic - not tracking dispersal in this case)
bothDrivers <- aggregate(allRLsum[c('dataBio','dataAbio')], 
                        by=allRLsum[c('Sp.name','Continent','RL.Lat.Elev','RL')], function(x) sum(x))
    head(bothDrivers)
    summary(bothDrivers) #max=4 different studies looked at abiotic Drivers for that species x RL
    hist(bothDrivers$dataAbio)
    bothDrivers[bothDrivers$dataAbio>=3,]
allRLsum[allRLsum$Sp.name=='Abies.lasiocarpa',] #extra rows come from multiple studies per RL
#create variables denoting whether a study looked at both abiotic AND biotic drivers
bothDrivers$spRL.dataAbioandBio <- as.factor(ifelse(bothDrivers$dataAbio>0 & bothDrivers$dataBio>0, 'Y', 'N'))
    summary(bothDrivers) 
    table(bothDrivers$spRL.dataAbioandBio, bothDrivers$RL) #wow - LOTS (>130 per RL) of species x RL combinations have data on both biotic and abiotic drivers
    test <- bothDrivers[bothDrivers$spRL.dataAbioandBio=='Y',]; test <- droplevels(test)
    summary(test)    
    nlevels(test$Sp.name) #>200 species! 
    bothDrivers$dataBio <- bothDrivers$dataAbio <- NULL #get rid of these for merging

#merge back onto main (allRL) data  
allRLformerge <- allRL
dim(allRLformerge) #2125 
dim(bothDrivers) #889
#merge
allRL <- merge(allRLformerge, bothDrivers)
dim(allRL) #2125 x 1 col extra
summary(allRL)
nlevels(allRL$doi) #340 # DIFF ##############

#merge back onto main (allRLsum) data
allRLsumformerge <- allRLsum
dim(allRLsumformerge) 
#merge
names(allRLsum)
names(bothDrivers)
allRLsum <- merge(allRLsumformerge, bothDrivers)
dim(allRLsum) 
summary(allRLsum)

#__make var denoting which species x location have data on both RL types (cool & warm) -------- ---
#should be only studies that assess both RLs, or only species for which have data on both RLs? Think the second. eg authors might present data on cool RL in one paper and warm RL in another paper
summary(allRLsum) #should already be only 1 line per sp x RL x driver type    
allRLsum$datacool <- ifelse(allRLsum$RL=='cool', 1, 0) 
allRLsum$datawarm <- ifelse(allRLsum$RL=='warm', 1, 0)

#aggregate across RL (cool & warm)
#expect to get sums of 1 (1 driver type at 1 spRL) to 6 (all 3 driver types (bio, abio, dispersal) at both RLs). but can get more if multiple studies done on same RL
bothRL <- aggregate(allRLsum[c('datacool','datawarm')], 
                    by=allRLsum[c('Sp.name','Continent','RL.Lat.Elev')], function(x) sum(x))
head(bothRL)
summary(bothRL) #getting some high values for data cool & data warm (>4)
hist(bothRL$datacool)
bothRL[bothRL$datacool>4,]
#check species with more rows than expected
allRLsum[allRLsum$Sp.name=='Abies.amabilis',] #extra rows come from dispersal
allRLsum[allRLsum$Sp.name=='Abies.lasiocarpa',] #extra rows come from multiple studies per RL
allRLsum[allRLsum$Sp.name=='Tsuga.heterophylla',] #extra rows come from multiple studies per RL + disp
#create variables denoting whether a study looked at both RL
bothRL$bothRL <- as.factor(ifelse(bothRL$datacool>0 & bothRL$datawarm>0, 'Y', 'N'))
summary(bothRL)
bothRL$datacool <- bothRL$datawarm <- NULL #remove so doesnt mess up merging

#merge back onto main (allRL) data
allRLformerge <- allRL
dim(allRLformerge) #2125
dim(bothRL) #702  
#merge
allRL <- merge(allRLformerge, bothRL)
dim(allRL) #2125 
summary(allRL)
table(allRL$spRL.dataAbioandBio, allRL$bothRL) #have 607 rows of data where both abio and bio factors were tested at both warm and cool RL!

#merge back onto main (allRLsum) data
allRLsumformerge <- allRLsum
dim(allRLsumformerge) #1447
#merge
names(allRLsum)
names(bothRL)
allRLsum <- merge(allRLsumformerge, bothRL)
dim(allRLsum) #should be same rows
summary(allRLsum)
table(allRLsum$spRL.dataAbioandBio, allRLsum$bothRL) #have >340 rows of data where both abio and bio factors were tested at both warm and cool RL!


#__get rid of dispersal---------------------------
#will help make sure dont get sample sizes wrong
allRLwD <- allRL
allRL <- allRLwD[allRLwD$Driver.type!='dispersal',]; allRL <- droplevels(allRL)
dim(allRLwD); dim(allRL) #gets rid of ~200 lines. final size = 1941
nlevels(allRLwD$doi); nlevels(allRL$doi) #lose 1 study, as expected

allRLsumwD <- allRLsum
allRLsum <- allRLsumwD[allRLsumwD$Driver.type!='dispersal',]; allRLsum <- droplevels(allRLsum)
dim(allRLsumwD); dim(allRLsum) #gets rid of ~200 lines


# DATA LOCATIONS ########################################
#lat & long taken from the range limit cited by the authors. If there were multiple study sites we used the highest latitude/elevation site where the species occurred for cool limits, and the lowest for warm. For the occasional global model or study over several countries or an ocean, we either tried to match coordinates to a map of the range they provided, or if it was truly global I did not put coordinates

#__data cool limits locations ---------------------
location.cool.orig <- read.csv("Data/Paquette&Hargreaves locations cool limits.csv", skipNul=T, stringsAsFactors=T)
#get rid of blank lines
location.cool.orig <- location.cool.orig[location.cool.orig$Genus!='',]; 
location.cool.orig <- droplevels(location.cool.orig) 
summary(location.cool.orig) #Taxon=mostly blank - had to add at end to deal with species name where Genus = multiple
location.cool.orig[duplicated(location.cool.orig),] #should be none

#keep only columns that dont duplicate column in driver data above 
names(location.cool.orig)
loc.cool <- location.cool.orig[,c('Authors', 'Year.pub', 'doi', 'Genus', 'Species', 'Taxon',
                                       'RL.Lat.Elev', 'Continent', 'Latitude', 'Longitude')]
#create species name as in lit search data
summary(loc.cool[,c('Genus','Taxon')])
loc.cool[loc.cool$Taxon=='Mixed',]
loc.cool[loc.cool$Genus=='Multiple',]
loc.cool$Sp.name <- as.factor(ifelse(loc.cool$doi=='10.1046/j.1365-2699.1998.2540735.x', 'southern trees', 
                              ifelse(loc.cool$doi=='10.1111/j.1654-1103.2010.01204.x', 'Heathland plants',    
                              ifelse(loc.cool$Taxon!='Mixed' & loc.cool$Genus=='Multiple', paste(loc.cool$Taxon, loc.cool$Species, sep='.'), 
                                     paste(loc.cool$Genus, loc.cool$Species, sep='.')))))
levels(loc.cool$Sp.name)
summary(loc.cool)
#check no blanks for DOI
loc.cool[loc.cool$doi=='',]

#check that latitudes make sense for cool limits
summary(loc.cool$Latitude)
hist(x=as.numeric(loc.cool$Latitude)) 
#have some cool limits that are really close to the equator - many are Elevational limits, which is fine. Alex checked the latitudinal ones:
loc.cool[loc.cool$RL.Lat.Elev=='Latitudinal' & loc.cool$Latitude > -10 & loc.cool$Latitude < 10, c('Authors','Year.pub','Genus','Species','RL.Lat.Elev','Latitude')]
#Asseh 3 spp: northern RLs in Cote D'Ivoire so these are correct
#Monahan Sporophila americana - really is a super equatorial species
#Freeman 3 spp: West Africa, northern limits = cool edge, so limits really are equatorial


#__data warm limits locations ---------------------
location.warm.orig <- read.csv("Data/Paquette&Hargreaves locations cool limits.csv", skipNul = TRUE, stringsAsFactors =T)
location.warm.orig <- Filter(function(x)!all(is.na(x)), location.warm.orig) #gets rid of any weird columns at end
summary(location.warm.orig)
#check cases with missing locations
test <- location.warm.orig[is.na(location.warm.orig$Latitude),] 
test <- droplevels(test); summary(test) #0 rows with missing data. from mixed-taxa model
nlevels(test$doi) #should be from 0 studies
location.warm.orig[location.warm.orig$Exclude=='Y',]
#check for duplicated rows
location.warm.orig[duplicated(location.warm.orig),] #should be none

#keep only columns that dont duplicate literature search data
names(location.warm.orig)
loc.warm <- location.warm.orig[
                              #loc.warm.orig$Exclude!='Y', #if want to get rid of Cahill excluded studies right away
                              c('Authors', 'Year.pub', 'doi', 'Genus', 'Species', 'Taxon', 
                                'RL.Lat.Elev', 'Continent', 'Latitude', 'Longitude')]
dim(location.warm.orig); dim(loc.warm) #gets rid of 14 rows if get rid of Cahill excluded studies
loc.warm <- droplevels(loc.warm)

#create species name as in lit search data
summary(loc.warm[,c('Genus','Taxon')])
loc.warm[loc.warm$Taxon=='Mixed',] #some studies have already been excluded from RL data. dont need name fixes for:
    #Bruelheide 2003
    #Feely 2012
    #Lenoir 2010
    #Normand
loc.cool[loc.warm$Genus=='Multiple',]
loc.warm$Sp.name <- as.factor(ifelse(loc.warm$doi=='10.1890/10-0312.1', 'temperate trees', #Canham
                              ifelse(loc.warm$doi=='10.1046/j.1365-2699.1998.2540735.x', 'boreal trees', #loehle
                              ifelse(loc.warm$doi=='10.1111/j.1654-1103.2010.01204.x', 'Heathland plants', #Loidi
                              ifelse(loc.warm$doi=='10.1111/j.1654-1103.2010.01201.x', 'European plants', #Lenoir
                              ifelse(loc.warm$Taxon!='Mixed' & loc.warm$Genus=='Multiple', paste(loc.warm$Taxon, loc.warm$Species, sep='.'), 
                                     paste(loc.warm$Genus, loc.warm$Species, sep='.')))))))
loc.warm[loc.warm$Genus=='Multiple',] 
#check no blanks for DOI
loc.warm[loc.warm$doi=='',]

head(loc.warm)

#check that latitudes make sense for warm limits
hist(loc.warm$Latitude)
#have some warm limits that are above 60 degrees latitude - fine if they are elevational limits but check the latitudinal ones
loc.warm[loc.warm$RL.Lat.Elev=='Latitudinal' & (loc.warm$Latitude > 60 | loc.warm$Latitude < -60), c('Authors','Year.pub','Genus','Species','RL.Lat.Elev','Latitude')]
#Bischoff-Basmann 1996 – just very polar (highest distribution in this paper is below South Georgia Island around -55 degrees)
#Van Dijk 1999 – genuinely just circles Antarctica


#__merge cool & warm location data-----------------------------

#create columns to identify where data come from
loc.cool$RL <- as.factor('cool')
loc.warm$RL <- as.factor('warm')

names(loc.cool)
names(loc.warm)

locations <- rbind(loc.cool, loc.warm) 
summary(locations)
locations$Year.pub <- as.factor(locations$Year.pub); levels(locations$Year.pub)
table(locations$RL.Lat.Elev, locations$RL) #very even split!
locations[locations$RL.Lat.Elev=='',] 
locations[duplicated(locations),] #there are 0 rows duplicated


#__merge locations & drivers data --------------------

summary(locations)
loc.formerge <- locations[!is.na(locations$Latitude),c('doi', 'Sp.name', 'RL', 'RL.Lat.Elev', 'Continent','Latitude', 'Longitude')]
loc.formerge <- droplevels(loc.formerge)
loc.formerge[duplicated(loc.formerge),] #should be none
dim(loc.formerge) #1026 rows # DIFF ###################
#create simpler version of allRL - makes easier to check merging errors and dont need all these columns
summary(allRL)
allRLformerge <- allRL
allRLformerge$sppRLloc <- allRLformerge$RLloc <- allRLformerge$Method.all <- allRLformerge$spRL.dataAbioandBio <- NULL
allRLformerge$Driver.supported <- allRLformerge$Driver <- allRLformerge$double.checker <-NULL
allRLformerge$Study.inCahill2014 <- allRLformerge$bothRL <-  NULL
summary(allRLformerge)
dim(allRLformerge) #1931
#compare number of studies in each
nlevels(loc.formerge$doi) #353 - more bc has some studies we later exclude # DIFF ##############
nlevels(allRLformerge$doi) #339 # DIFF ##################3
#merge
RLloc <- merge(allRLformerge, loc.formerge, all.x=T, by=c('doi', 'Sp.name', 'RL', 'RL.Lat.Elev', 'Continent'))
dim(allRLformerge); #1931 
dim(RLloc) #should add 2 columns & 0 rows. if not get ready for the nightmare of figuring our the duplicated rows... 
RLloc[duplicated(RLloc),] #should be none
# #must be a better way to find which rows are getting duplicated but I dont know what it is....
# summary(allRLformerge[,c('RL','High.taxon','Environment','Continent','Method')]) #4 extra rows in Plant cool land RLs, 2 North America, 2 Europe
# summary(RLloc[,c('RL','High.taxon','Environment','Continent','Method')])
# test <- RLloc[RLloc$RL=='cool' & RLloc$Taxon=='Herbaceous' & RLloc$Environment=='terrestrial' & RLloc$Method=='Model' &
#                 (RLloc$Continent=='North America' | RLloc$Continent=='Europe'),]; test <- droplevels(test);
# testall <- allRL[allRL$RL=='cool' & allRL$Taxon=='Herbaceous' & allRL$Environment=='terrestrial' & allRL$Method=='Model' &
#                 (allRL$Continent=='North America' | allRL$Continent=='Europe'),]; testall <- droplevels(testall);
# summary(testall[,c('Genus')]); summary(test[,c('Genus')]) #still problems with Ambrosia!
# #are there any NAs for Latitude? (YES - 13)
summary(loc.formerge$Latitude); summary(RLloc$Latitude)  #no NAs for Latitude in location data, so shouldnt be in RLloc
testNA <- RLloc[is.na(RLloc$Latitude),]; testNA <- droplevels(testNA)
testNA

RLloc$absLat <- abs(RLloc$Latitude)
RLloc$Latzone <- as.factor(ifelse(RLloc$absLat<23.5, 'tropical', ifelse(RLloc$absLat>66.5, 'polar', 'temperate')))
summary(RLloc)
xyplot(absLat ~ Latitude, data=RLloc)

#check doi numbers again - possible that used diff spellings of doi in each database (eg include 'doi' or not)
test <- aggregate(RLloc[c('doi')], by=RLloc[c('Authors','Year.pub')], function(x) length(unique(x)))
summary(test)
test[test$doi>1,]
RLloc[RLloc$Authors=='Allen, CD; Breshears DD',] #yup, 2 papers in same year
RLloc[RLloc$Authors=='Crozier, LG',] #yup, 2 papers in same year

#check for missing locations
summary(loc.formerge) #not missing any latitudes
summary(RLloc) #so should not be any NAs for latitude
# #if there are NAs check them with code below
# test <- RLloc[is.na(RLloc$Latitude),]; test <- droplevels(test); test
# test[test$doi=='10.1046/j.1469-8137.1998.00232.x',]
# loc.formerge[loc.formerge$doi=='10.1046/j.1469-8137.1998.00232.x' & loc.formerge$RL=='warm',]
# allRLformerge[allRLformerge$doi=='10.1046/j.1469-8137.1998.00232.x',]

# allRLloc <- RLloc
# RLloc <- allRLloc[!is.na(allRLloc$Latitude),]; RLloc <- droplevels(RLloc)
# summary(RLloc)
# summary(allRLloc)
# write.csv(allRLloc, 'Data/data check missing latitudes 21 03 25.csv')


#__merge by factor.type data ----------------- ---

nlevels(allRLsum$doi) #338
nlevels(loc.formerge$doi) #352 # DIF #################
RLlocsum <- merge(allRLsum, loc.formerge, all.x=T, all.y=F, by=c('doi', 'Sp.name', 'RL', 'RL.Lat.Elev','Continent'))
dim(allRLsum); #1263 x 20 # DIFF ##################3
dim(RLlocsum) #adds 2 columns - shouldnt add any rows
RLlocsum[duplicated(RLlocsum),] #0 rows duplicated
summary(RLlocsum)

#get rid of NA lines if any
RLlocsum <- RLlocsum[!is.na(RLlocsum$Latitude),]; RLlocsum <- droplevels(RLlocsum)
summary(RLlocsum)
#absoluate latitude 
RLlocsum$absLat <- abs(RLlocsum$Latitude)
RLlocsum$Latzone <- as.factor(ifelse(RLlocsum$absLat<23.5, 'tropical', ifelse(RLlocsum$absLat>66.5, 'polar', 'temperate')))
summary(RLlocsum)
#which are polar studies? (40 polar data points)
RLlocsum[RLlocsum$Latzone=='polar',]


# DATA SUMMARIZING ################################################################################

summary(allRL)
table(allRL$Driver.examined, allRL$Driver.type)
levels(allRL$doi) #make sure no 'doix AND doiy' ones (from Cahill)
allRL[grepl(' & ', allRL$doi),] 
hist(allRLsum$sumDriver.support)

#__Abstract------------------
#how many taxon-range limits?
summary(allRL)
nlevels(allRL$sppRLloc) #885 taxon range limits
#get % of range limits involving interactions from ANALYSES > 1b

#__Results text------------------
#examples of studies that assess and find support for multiple factors
summary(allRLsum)
test <- allRLsum[allRLsum$nDrivers.tested>1 & allRLsum$sumDriver.support>2,]
test <- droplevels(test); summary(test)
test

#1st paragraph Results (get fraction of RL interactions contribute to in Analyses > 1b below)
names(allRL) 
nlevels(allRL$doi) # n studies (338)
dim(allRL); n <- 1931 # #data points
nlevels(allRL$Sp.name) # n taxa (including some species groups) (654)
test <- allRL[allRL$Genus!='Multiple' & allRL$Species!='spp.',]; test <- droplevels(test); nlevels(test$Sp.name) # unique species or ssp 630
654-630 #24 groups of species
test <- allRL[allRL$Genus=='Multiple' | allRL$Species=='Multiple' | allRL$Species=='spp.',]; test <- droplevels(test); 
        levels(test$Sp.name) #examples of taxa grouped across species
test <- aggregate(allRL[c('Sp.name')], by=allRL[c('High.taxon')], function(x) length(unique(x))); test # # unique plants/algae (359), animals (292)
test <- aggregate(allRL[c('Sp.name')], by=allRL[c('High.taxon','Taxon')], function(x) length(unique(x))); test
#vascular land plants
p <- test[test$High.taxon=='Plant/Algae' & test$Taxon!='Algae' & test$Taxon!='Moss' & test$Taxon!='Seagrass',]; p; sum(p$Sp.name) #land vasc plants: 314
        314/359 # % of plants/algae that are terrestrial vascular plants
a <- test[test$High.taxon=='Animal' & (test$Taxon=='Amphibian' | test$Taxon=='Bird' | test$Taxon=='Fish' | test$Taxon=='Mammal' | test$Taxon=='Reptile'),]; a
      nvert = sum(a$Sp.name); nvert #vertebrates 157
      nvert/292 # % of animals that are vertebrates
table(allRL$Driver.type); 523/(n) # %biotic drivers
#get % of range limits biotic factors contribute to in ANALYSES > 1b
# % data from land, marine, freshwater
table(allRL$Environment); 1502/n; 368/n; 61/n 
test <- aggregate(allRL[c('doi')], by=allRL[c('Method')], function(x) length(unique(x))); test #unique dois per Method
        122/(nlevels(allRL$doi)) # % of studies that included field experiments
# % data from elevational vs latitudinal limits
table(allRL$RL.Lat.Elev); 904/n #lat range limits = 47%      
#number of factors assessed per study type (stated but no numbers given)
test <- aggregate(allRL[c('Driver.examined')], by=allRL[c('Method','doi')], length); head(test) #driver assessed per study for each study type
        bwplot(Driver.examined ~ Method, data=test) #some observational study assessed >200 drivers x range limits x taxa
        test[test$Driver.examined>100,] #3 studies
        allRL[allRL$doi=='10.1007/s12224-010-9059-4',c(1:10)]
        bwplot(Driver.examined ~ Method, data=test, ylim=c(0,50)) #Bassler etal >8 spp, assesed cool & warm limits for some & up to 7 factors per RL
        mean <- aggregate(test[c('Driver.examined')], by=test[c('Method')], mean); mean #driver assessed per study for each study type
table(allRL$Method); 545/n; #factor assessments from field exp


#__Table 1------------------
        
#Table 1: n studies
nlevels(allRL$doi) #338 distinct papers
test <- aggregate(allRL[c('doi')], by=allRL[c('RL')], function(x) length(unique(x))); test #this gives #unique dois per RL type

#Table 1: n taxa (counts rows aggregated across genera or spp)?
nlevels(allRL$Sp.name) #654 'taxa' groups
test <- aggregate(allRL[c('Sp.name')], by=allRL[c('RL')], function(x) length(unique(x))); test #this gives #unique taxa per RL type

#Table 1: n range limits (species x RL (cool/warm lat/elev) x continent) combinations?
nlevels(allRL$sppRLloc) 
test <- aggregate(allRL[c('sppRLloc')], by=allRL[c('RL')], function(x) length(unique(x))) 
test

#Table 1: n range limit x driver combinations (ie how much data for 'support by factor' analyses)?
table(allRL$Driver.type); 
nA=1408; nB=523; nB/(nB+nA) #28% of data is on biotic factors
table(allRL$Driver.type, allRL$RL)

# #Table 1: n range limit x driver combinations (ie how much data for 'support by factor type' analyses)?
# names(allRLsum)
# table(allRLsum$Driver.type); 
# nA=851; nB=412; nB/(nB+nA) #33 biotic%
# table(allRLsum$Driver.type, allRLsum$RL)

#__Table S1------------------
table(allRL$Environment, allRL$Method, allRL$Driver.type)
#final column (%biotic across methods)
table(allRL$Environment, allRL$Driver.type)
5/(56+5) #freshwater
75/(293+75) #marine
443/(1059+443) #terrestrial
#final row (%biotic across environments)
table(allRL$Method, allRL$Driver.type)
230/(315+230) #field exp
4/(199+4) #lab exp
118/(519+118) #models
171/(375+171) #observational data

#__Discussion -----------------
#how many studies look at mutualisms (paragraph 5 of Discussion)?
table(allRL$Driver.descrip)
test <- aggregate(allRL[c('doi')], by=allRL[c('Driver.descrip')], function(x) length(unique(x))); test #lack mutualists = 14 studies
nlevels(allRL$doi) #338 studies

#how do people study competition?
test <- allRL[allRL$Driver.examined=='Competition',]; test <- droplevels(test)
summary(test)
levels(test$Driver.notes)
test[grepl('several', test$Driver.notes),]



#__SI & other (not in paper)------------------

#SI 1.2 - When multiple methods contribute
summary(allRL$Method.all) #how many data points are conclusions derived from >1 method?
2+4+15+2+3+4+1+1+7+1+5

#SI 1.3 Data checking: how many of data points we collected were double checked?
test <- allRL[allRL$double.checker=='',]; dim(test)
notchecked <- 776
dim(allRL); n <- 1931 #data points (1 per factor)
(n-notchecked)/n #60% of data = double checked
#of data we collected
test <- allRL[allRL$data.extractor!='Cahill et al',]; test <- droplevels(test); 
dim(test); n <- 1584
notchecked <- test[test$double.checker=='',]; dim(notchecked)
notchecked <- 655
(n-notchecked)/n #58.6% of data we collected was double checked. so same frequency
# of Cahill data
test <- allRL[allRL$data.extractor=='Cahill et al',]; test <- droplevels(test); dim(test); n <- 347
summary(test$double.checker); notchecked <- 121
(n-notchecked)/n #65% of Cahill was double checked. so same frequency
# since we emphasize results from field experiments, how many of them were double checked?
levels(allRL$Method)
table(allRL$double.checker, allRL$Method) #36 rows not checked (all collected by ALH)
table(allRL$Method); nF=545 #of 546 rows
(nF-36)/nF #94% of field data checked 

#for map - how many cases were missing longitude? 
summary(loc.cool) #1 missing long (none missing lat)
summary(loc.warm) #0 missing long (none missing lat)
locations[is.na(locations$Longitude),]
summary(RLlocsum) #already excluded 1 study where could not get latitude or longitude
#Map text: how many range limits?
test <- aggregate(allRL[c('sppRLloc')], by=allRL[c('RL','Driver.type')], function(x) length(unique(x))) 
test

#additional stuff not in paper or SI
#are any species looked at in more than 1 study? There's one row per driver so this is a bit hard to figure out
test <- aggregate(allRL[c('doi')], by=allRL[c('Sp.name')], function(x) length(unique(x))) #this gives #unique dois per species
hist(test$doi) #by far the majority only studied in 1 paper, max = 5 papers
test[test$doi==5,] #check species w high numbers to make sure no doi problems
allRL[allRL$Sp.name=='Fagus.sylvatica',] #really is 5 studies...
allRL[allRL$Sp.name=='Semibalanus.balanoides',] #really is 5 studies...
test[test$doi==4,] #check species w high numbers to make sure no doi problems
allRL[allRL$Sp.name=='Picea.glauca',] #really is 5 studies...
allRL[allRL$Sp.name=='Populus.tremuloides',] #really are 4 studies (warm & cool data published in separate studies by same authors in same year)
allRL[allRL$Sp.name=='Tsuga.heterophylla',] #really are 4 studies
test[test$doi==3,] #check species w high numbers to make sure no doi problems

#How many species per study? (even more in some of the 'Multiple' Genus studies)
test <- aggregate(allRL[c('Sp.name')], by=allRL[c('doi')], function(x) length(unique(x))) 
hist(test$Sp.name) #biggest study has 30 species



# ANALYSES ################################################################################
#will want to assess overdispersion in binomial models. Ben Bolker gives code here: https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion
overdisp_Bolker <- function(model) {
    rdf <- df.residual(model)
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
    
#Q1) bio vs abio at cool vs warm RL (=full model) ######################

#1a) bothRL by driver: RESULTS Driver x RL = signif  -------------------
#do we have good enough data coverage?
table(allRL$Driver.type, allRL$RL) #yup. dispersal data already removed

#get sample size
n <- allRL[allRL$Driver.type!='dispersal',]; n <- droplevels(n); dim(n); 
#random factors for sp & study boundary=singular
coolvwarm.bydriver <- glmer(Driver.supported01 ~ Driver.type*RL + (1|Sp.name) + (1|doi),
                        family=binomial,
                        control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                        data=allRL) #boundary=singular
    #test significance of factor.type * RL interaction
    coolvwarm.bydrivernoX <- glmer(Driver.supported01 ~ Driver.type + RL + (1|Sp.name) + (1|doi),
                             family=binomial,
                             control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                             data=allRL) #singular
    anova(coolvwarm.bydriver, coolvwarm.bydrivernoX, test='Chisq') #interaction signif
    # # #check 'boundary singular' warning - goes away if get rid of species random effect, fixed effect results stay same 
    # coolvwarm.bydriver.nosp <- glmer(Driver.supported01 ~ Driver.type*RL + (1|doi),
    #                     family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
    #                     data=allRL) #boundary=singular
    # coolvwarm.bydrivernoX.nosp <- glmer(Driver.supported01 ~ Driver.type+RL + (1|doi),
    #                     family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
    #                     data=allRL) #boundary=singular
    # anova(coolvwarm.bydriver.nosp, coolvwarm.bydrivernoX.nosp, test='Chisq') #fixed effect results exactly the same
        visreg(coolvwarm.bydriver, xvar='Driver.type', type='conditional', ylab='support by driver', ylim=c(-3,3)) 
        visreg(coolvwarm.bydriver, xvar='Driver.type', by='RL', type='conditional', ylab='support by driver', ylim=c(-3,3), layout=c(2,1)) 
    #Prediction 1 (cool vs warm for given driver type)
    lsmeans(coolvwarm.bydriver, pairwise ~ RL | Driver.type, transform='response')
    #biowarm = .62, bio cool = .39. so not quite twice as often across study types
    (.619-.394)/.619 #bio supported 36% more often at warm vs cool limits
    #Prediction 2 (biotic vs abiotic at given RL)
    lsmeans(coolvwarm.bydriver, pairwise ~ Driver.type | RL, transform='response')
    #at cool: abiotic > biotic, abiotic > .5, biotic < .5 
    #at warm: abiotic = biotic, abiotic > .5, biotic = .5
    visreg(coolvwarm.bydriver, xvar='Driver.type', by='RL', type='conditional', ylab='support by driver', 
           ylim=c(-3,3), layout=c(2,1)) 
#double check lsmean results using bootstrapped CI (did but ignore. author of lsmeans says CI can be misleading in GLMMS as some random effects cancel out, so lsmeans more reliable https://stats.stackexchange.com/questions/181034/post-hoc-test-of-interaction-factor-in-binomial-glmm-with-proportions)

    
#_sensitivity tests ------------------------------------
#__i) only field experiments  ------
#(few lab tests of biotic factors at cool RLs):
table(allRL$RL, allRL$Method) #looks good at first but...
table(allRL$Driver.type, allRL$Method, allRL$RL) #very few tests of biotic drivers in Lab experiments. 

#get sample size
n <- allRL[allRL$Method=='Field experiment',]; n <- droplevels(n); dim(n); 
levels(n$doi) #122 studies (for Fig 2 caption)
coolvwarm.fexp <- glmer(Driver.supported01 ~ Driver.type*RL + (1|Sp.name) + (1|doi), 
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                    data=allRL[allRL$Driver.type!='dispersal' & allRL$Method=='Field experiment',]) #singular
    #reduced model to test significance
    coolvwarm.fexpnoX <- glmer(Driver.supported01 ~ Driver.type + RL + (1|Sp.name) + (1|doi),
                         family=binomial,
                    data=allRL[allRL$Driver.type!='dispersal' & allRL$Method=='Field experiment',]) #singular
    anova(coolvwarm.fexp, coolvwarm.fexpnoX, test='Chisq') #interaction SIGNIFICANT
    visreg(coolvwarm.fexp, xvar='Driver.type', by='RL', type='conditional', layout=c(2,1)) #at warm RLs support for abiotic drops & support for biotic increases 
    #Prediction 1 (cool vs warm for given driver type)
    lsmeans(coolvwarm.fexp, pairwise ~ RL | Driver.type, transform='response')
    #get increase in support for bio at warm vs cool limits
    (.753-.320)/.753 #bio supported 57.5% more often at warm vs cool limits
    #Prediction 2 (biotic vs abiotic at given RL)
    lsmeans(coolvwarm.fexp, pairwise ~ Driver.type | RL, transform='response')
    
     
#__ii) excl mixed-species data -----
#how many would this exclude?
test <- allRL[allRL$Species=='Multiple' | allRL$Species=='spp.',]; summary(test) #75 multi-species data points
#get sample size
n <- allRL[allRL$Species!='Multiple' & allRL$Species!='spp.',]; n <- droplevels(n); dim(n); 
#levels(n$Species)
coolvwarm.nomultsp <- glmer(Driver.supported01 ~ Driver.type*RL + (1|Sp.name) + (1|doi), 
                    family=binomial,
                    control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                    data=allRL[allRL$Species!='Multiple' & allRL$Species!='spp.',]) #singular
    #reduced model to test significance of interaction
    coolvwarm.nomultspnoX <- glmer(Driver.supported01 ~ Driver.type + RL + (1|Sp.name) + (1|doi),
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRL[allRL$Species!='Multiple' & allRL$Species!='spp.',]) #singular
    anova(coolvwarm.nomultsp, coolvwarm.nomultspnoX, test='Chisq') #interaction signif
    #Prediction 1 (cool vs warm for given driver type)
    lsmeans(coolvwarm.nomultsp, pairwise ~ RL | Driver.type, transform='response')
    #Prediction 2 (biotic vs abiotic at given RL)
    lsmeans(coolvwarm.nomultsp, pairwise ~ Driver.type | RL, transform='response')
        
    
#__iii) only sppRL w both driver types studied -----
#get sample size
n <- allRL[allRL$Driver.type!='dispersal' & allRL$spRL.dataAbioandBio=='Y',]; n <- droplevels(n); 
dim(n); 
test <- aggregate(n[c('doi')], by=n[c('Driver.type','Method')], function(x) length(unique(x))); test

coolvwarm.bothAB <- glmer(Driver.supported01 ~ Driver.type*RL + (1|Sp.name) + (1|doi), #no warning
                        family=binomial,
                        data=allRL[allRL$spRL.dataAbioandBio=='Y',]) 
coolvwarm.bothABnoX <- glmer(Driver.supported01 ~ Driver.type + RL + (1|Sp.name) + (1|doi), #no warning
                        family=binomial,
                        data=allRL[allRL$spRL.dataAbioandBio=='Y',]) 
    anova(coolvwarm.bothAB, coolvwarm.bothABnoX, test='Chisq') #signif
    visreg(coolvwarm.bothAB, xvar='Driver.type', by='RL', type='conditional', layout=c(2,1)) #more extreme
    #Prediction 1 (cool vs warm for given driver type)
    lsmeans(coolvwarm.bothAB, pairwise ~ RL | Driver.type, transform='response')
    #Prediction 2 (biotic vs abiotic at given RL)
    lsmeans(coolvwarm.bothAB, pairwise ~ Driver.type | RL, transform='response')

    
#__iv) only spp x region x RL type (Elev or Lat) where cool & warm RLs studied -----
summary(allRL)
table(allRL$Method, allRL$bothRL) #note can have diff # of dois per range limits as doing by RL, not by study
#get sample size - models and lab experiments are mostly testing abiotic drivers
n <- allRL[allRL$bothRL=='Y',]; n <- droplevels(n); dim(n); 
test <- aggregate(n[c('doi')], by=n[c('Driver.type','Method')], function(x) length(unique(x))); test

coolvwarm.bothRL <- glmer(Driver.supported01 ~ Driver.type*RL + (1|Sp.name) + (1|doi), 
                        family=binomial,
                        data=allRL[ allRL$bothRL=='Y',]) 
coolvwarm.bothRLnoX <- glmer(Driver.supported01 ~ Driver.type + RL + (1|Sp.name) + (1|doi), #singular
                        family=binomial,
                        data=allRL[allRL$bothRL=='Y',]) 
    anova(coolvwarm.bothRL, coolvwarm.bothRLnoX, test='Chisq') #interaction significant
    visreg(coolvwarm.bothRL, xvar='Driver.type', by='RL', type='conditional', layout=c(2,1)) #more extreme at cool
    #Prediction 1 (cool vs warm for given driver type)
    lsmeans(coolvwarm.bothRL, pairwise ~ RL | Driver.type, transform='response')
    visreg(coolvwarm.bothRL, xvar='RL', by='Driver.type', type='conditional', layout=c(2,1)) #more extreme at cool
    #Prediction 2 (biotic vs abiotic at given RL)
    lsmeans(coolvwarm.bothRL, pairwise ~ Driver.type | RL, transform='response')
        
    
# #__abio & bio cool & warm, include field exp data (close as we can get to ideal (IMO) data but dataset getting pretty small)  --------- --
# #get sample size - models and lab experiments are mostly testing abiotic drivers
# n <- allRL[allRL$bothRL=='Y' & allRL$spRL.dataAbioandBio=='Y' & allRL$Method=='Field experiment',]; n <- droplevels(n); 
# dim(n); #150 rows of data
# levels(n$doi) #24 studies
# 
# coolvwarm.typebestf <- glmer(Driver.supported01 ~ Driver.type*RL + (1|Sp.name) + (1|doi), #no warning
#                                 family=binomial,
#                                 data=allRL, subset=bothRL=='Y' & spRL.dataAbioandBio=='Y' & Method=='Field experiment') #singular
# coolvwarm.typebestf.noX <- glmer(Driver.supported01 ~ Driver.type + RL + (1|Sp.name) + (1|doi), #no warning
#                                 family=binomial,
#                                 data=allRL, subset=bothRL=='Y' & spRL.dataAbioandBio=='Y' & Method=='Field experiment') #singular
#       anova(coolvwarm.typebestf, coolvwarm.typebestf.noX, test='Chisq') #int signif
#       #Prediction 1 (cool vs warm for given driver type)
#       lsmeans(coolvwarm.typebestf, pairwise ~ RL | Driver.type, transform='response')
#       #both predictions supported
#       #Prediction 2 (biotic vs abiotic at given RL)
#       lsmeans(coolvwarm.typebestf, pairwise ~ Driver.type | RL, transform='response')


#1b) both RLs by driver.type -------------------------
bwplot(Drivertype.supported01 ~ Driver.type | RL, data=allRLsum)
summary(allRLsum)

#get sample size
n <- allRLsum[allRLsum$Driver.type!='dispersal',]; n <- droplevels(n); dim(n);
#with covariate for how many factors were tested at that range limit
coolvwarm.bydrivertype.n <- glmer(Drivertype.supported01 ~ Driver.type*RL + nDrivers.tested + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRLsum) #singular
    overdisp_Bolker(coolvwarm.bydrivertype.n) #not overdispersed
    coolvwarm.bydrivertype.n.noX <- glmer(Drivertype.supported01 ~ Driver.type+RL + nDrivers.tested + (1|Sp.name) + (1|doi), 
                                        family=binomial, #singular (max grad warning - new March 29)
                                        data=allRLsum)
      anova(coolvwarm.bydrivertype.n, coolvwarm.bydrivertype.n.noX, test='Chisq') #interaction even more signif w n driver
      visreg(coolvwarm.bydrivertype.n, xvar='RL', by='Driver.type', ylab='support by driver type', main='w n covariate', 
             type='conditional', #scale='response', #get error message with scale = response 
             layout=c(2,1)) #interaction more imp at warm RLs
      Anova(coolvwarm.bydrivertype.n, type='III') #n drivers tested is super signif, not surprisingly
      #Prediction 1 (cool vs warm for given driver type)
      lsmeans(coolvwarm.bydrivertype.n, pairwise ~ RL | Driver.type, transform='response') 
      #Prediction 2 (biotic vs abiotic at given RL)
      lsmeans(coolvwarm.bydrivertype.n, pairwise ~ Driver.type | RL, transform='response') 
      #significance of covariate
      coolvwarm.bydrivertype <- glmer(Drivertype.supported01 ~ Driver.type*RL + (1|Sp.name) + (1|doi), #singular
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRLsum)
      anova(coolvwarm.bydrivertype.n, coolvwarm.bydrivertype, test='Chisq') #significance of n covariate
      #how many range limits overall do biotic interactions contribute to? 64.6% (SI 4.1)
      lsmeans(coolvwarm.bydrivertype.n, pairwise ~ Driver.type, transform='response') 

      
#_sensitivity tests ------------------------------------
#__i) only in incl field experiments -----------------
#combine data across studies so cant only use field experiments, but can require that at least some of data come from field experiments
summary(allRLsum)
#get sample size
n <- allRLsum[allRLsum$inclFexp=='Y',]; n <- droplevels(n); dim(n);

coolvwarm.dtypefexp <- glmer(Drivertype.supported01 ~ Driver.type*RL + nDrivers.tested + (1|Sp.name) + (1|doi), 
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                data=allRLsum[allRLsum$inclFexp=='Y',]) #singular
      overdisp_Bolker(coolvwarm.dtypefexp) #not overdispersed
coolvwarm.dtypefexp.noX <- glmer(Drivertype.supported01 ~ Driver.type + RL + nDrivers.tested + (1|Sp.name) + (1|doi),   
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                data=allRLsum[allRLsum$inclFexp=='Y',]) #singular
      anova(coolvwarm.dtypefexp, coolvwarm.dtypefexp.noX, test='Chisq') #int signif
      #Prediction 1 (cool vs warm for given driver type)
      lsmeans(coolvwarm.dtypefexp, pairwise ~ RL | Driver.type, transform='response')
      #Prediction 2 (biotic vs abiotic at given RL)
      lsmeans(coolvwarm.dtypefexp, pairwise ~ Driver.type | RL, transform='response') 
      visreg(coolvwarm.dtypefexp, xvar='Driver.type', by='RL', type='conditional', layout=c(2,1)) 
      visreg(coolvwarm.dtypefexp, xvar='RL', by='Driver.type', type='conditional', layout=c(2,1)) 
      #how many range limits overall do biotic interactions contribute to?
      lsmeans(coolvwarm.dtypefexp, pairwise ~ Driver.type, transform='response')


#__ii) excl mixed-species data ----
allRLsum[allRLsum$Species=='Multiple' | allRLsum$Species=='spp.',]  
#get sample size
n <- allRLsum[allRLsum$Species!='Multiple' & allRLsum$Species!='spp.',]; n <- droplevels(n); 
dim(n); 

coolvwarm.dtypenomult <- glmer(Drivertype.supported01 ~ Driver.type*RL + nDrivers.tested + (1|Sp.name) + (1|doi), 
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                data=allRLsum[allRLsum$Species!='Multiple' & allRLsum$Species!='spp.',]) #singular
    overdisp_Bolker(coolvwarm.dtypenomult) #not overdispersed
    coolvwarm.dtypenomult.noX <- glmer(Drivertype.supported01 ~ Driver.type+RL + nDrivers.tested + (1|Sp.name) + (1|doi), 
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                data=allRLsum[allRLsum$Species!='Multiple' & allRLsum$Species!='spp.',]) #singular
      anova(coolvwarm.dtypenomult, coolvwarm.dtypenomult.noX, test='Chisq') 
      visreg(coolvwarm.dtypenomult, xvar='RL', by='Driver.type', type='conditional', layout=c(2,1)) 
      #Prediction 1 (cool vs warm for given driver type)
      lsmeans(coolvwarm.dtypenomult, pairwise ~ RL | Driver.type, transform='response') 
      #Prediction 2 (biotic vs abiotic at given RL)
      lsmeans(coolvwarm.dtypenomult, pairwise ~ Driver.type | RL, transform='response') 

      
#__iii) only spp x RL where both biotic & abiotic studied -------------
names(allRLsum)
#get sample size
n <- allRLsum[allRLsum$Driver.type!='dispersal' & allRLsum$spRL.dataAbioandBio=='Y',]; n <- droplevels(n); 
dim(n); 
coolvwarm.dtypebothAB <- glmer(Drivertype.supported01 ~ Driver.type*RL + nDrivers.tested + (1|Sp.name) + (1|doi), 
                                family=binomial,
                                data=allRLsum[allRLsum$spRL.dataAbioandBio=='Y',]) #no warnings
      overdisp_Bolker(coolvwarm.dtypebothAB) #not overdispersed
coolvwarm.dtypebothAB.noX <- glmer(Drivertype.supported01 ~ Driver.type + RL + nDrivers.tested + (1|Sp.name) + (1|doi), 
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                data=allRLsum[allRLsum$spRL.dataAbioandBio=='Y',]) #singular
      anova(coolvwarm.dtypebothAB, coolvwarm.dtypebothAB.noX, test='Chisq') #int signif
      visreg(coolvwarm.dtypebothAB, xvar='Driver.type', by='RL', type='conditional', layout=c(2,1)); abline(h=0) 
      #Prediction 1 (cool vs warm for given driver type)
      lsmeans(coolvwarm.dtypebothAB, pairwise ~ RL | Driver.type, transform='response')
      #Prediction 2 (biotic vs abiotic at given RL)
      lsmeans(coolvwarm.dtypebothAB, pairwise ~ Driver.type | RL, transform='response') 

      
#__iv) only spp x region x RL type (Elev or Lat) where cool & warm RLs studied -----
#get sample size
n <- allRLsum[allRLsum$bothRL=='Y',]; n <- droplevels(n); dim(n); 

coolvwarm.dtypebothRL <- glmer(Drivertype.supported01 ~ Driver.type*RL + nDrivers.tested + (1|Sp.name) + (1|doi), 
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                data=allRLsum[allRLsum$bothRL=='Y',]) #singular
      overdisp_Bolker(coolvwarm.dtypebothRL) #not overdispersed
coolvwarm.dtypebothRL.noX <- glmer(Drivertype.supported01 ~ Driver.type + RL + nDrivers.tested + (1|Sp.name) + (1|doi), 
                                family=binomial,
                                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                                data=allRLsum[allRLsum$bothRL=='Y',]) #singular
      anova(coolvwarm.dtypebothRL, coolvwarm.dtypebothRL.noX, test='Chisq') #int signif
      visreg(coolvwarm.dtypebothRL, xvar='Driver.type', by='RL', type='conditional', layout=c(2,1)) 
      #Prediction 1 (cool vs warm for given driver type)
      lsmeans(coolvwarm.dtypebothRL, pairwise ~ RL | Driver.type, transform='response')
      visreg(coolvwarm.dtypebothRL, xvar='RL', by='Driver.type', type='conditional', layout=c(2,1)) 
      #Prediction 2 (biotic vs abiotic at given RL)
      lsmeans(coolvwarm.dtypebothRL, pairwise ~ Driver.type | RL, transform='response') 
    
      
# #__abio & bio  cool & warm, include field exp data (close as we can get to ideal (IMO) data but dataset getting pretty small)----- --
# #get sample size
# n <- allRLsum[allRLsum$bothRL=='Y' & allRLsum$spRL.dataAbioandBio=='Y' & allRLsum$inclFexp=='Y',]; n <- droplevels(n); dim(n)
# nlevels(n$Sp.name) #species
# nlevels(n$doi) #only 21 studies, and they dont necessarily each look at bio & abio factors at warm & cool RLs
# 
# coolvwarm.dtypebestf <- glmer(Drivertype.supported01 ~ Driver.type*RL + (1|Sp.name) + (1|doi), #no warning
#                                 family=binomial,
#                                 data=allRLsum, subset=bothRL=='Y' & spRL.dataAbioandBio=='Y' & inclFexp=='Y') #singular
#       overdisp_Bolker(coolvwarm.dtypebestf) #not overdispersed
#       coolvwarm.dtypebestf.noX <- glmer(Drivertype.supported01 ~ Driver.type + RL + (1|Sp.name) + (1|doi), #no warning
#                                     family=binomial,
#                                     data=allRLsum, subset=bothRL=='Y' & spRL.dataAbioandBio=='Y' & inclFexp=='Y') #singular
#       anova(coolvwarm.dtypebestf, coolvwarm.dtypebestf.noX, test='Chisq') #int signif
#       #Prediction 1 (cool vs warm for given driver type)
#       lsmeans(coolvwarm.dtypebestf, pairwise ~ RL | Driver.type, transform='response')
#       #Prediction 2 (biotic vs abiotic at given RL)
#       lsmeans(coolvwarm.dtypebestf, pairwise ~ Driver.type | RL, transform='response')
#       visreg(coolvwarm.dtypebestf, xvar='Driver.type', by='RL', type='conditional', layout=c(2,1)) 
#       visreg(coolvwarm.dtypebestf, xvar='RL', by='Driver.type', type='conditional', layout=c(2,1)) 
      
      
#_Table S2: Lat vs Elev & by Env -----------------------------------------------------
#__Lat vs elev (all data) --------------------------
#__support by factor --------------------------------- ---
table(allRL$Driver.type, allRL$RL.Lat.Elev, allRL$RL) #least data is for latitude x biotic 
dim(allRL) #n
coolvwarm.LE <- glmer(Driver.supported01 ~ Driver.type*RL*RL.Lat.Elev + (1|Sp.name) + (1|doi), 
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          family=binomial, data=allRL) 
coolvwarm.LE.no3X <- glmer(Driver.supported01 ~ Driver.type*RL + Driver.type*RL.Lat.Elev + RL*RL.Lat.Elev + (1|Sp.name) + (1|doi), 
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          family=binomial, data=allRL) 
      anova(coolvwarm.LE, coolvwarm.LE.no3X, test='Chisq') #3 way interaction is NS
      Anova(coolvwarm.LE.no3X, type='III')
#__support by factor type --------------------------------- ---
dim(allRLsum) #n
coolvwarm.LE <- glmer(Drivertype.supported01 ~ Driver.type*RL*RL.Lat.Elev + nDrivers.tested +(1|Sp.name) +(1|doi), 
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          family=binomial, data=allRLsum) 
coolvwarm.LE.no3X <- glmer(Drivertype.supported01 ~ Driver.type*RL + Driver.type*RL.Lat.Elev + RL*RL.Lat.Elev + nDrivers.tested +(1|Sp.name) +(1|doi), 
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          family=binomial, data=allRLsum) 
      anova(coolvwarm.LE, coolvwarm.LE.no3X, test='Chisq') #3 way interaction is NS
      Anova(coolvwarm.LE.no3X, type='III')

#__Lat vs elev (field exp) --------------------------------
#__support by factor --------------------------------- ---
table(allRL$Driver.type, allRL$RL.Lat.Elev, allRL$RL) #least data is for latitude x biotic 
test <- allRL[allRL$Method=='Field experiment',]; test <- droplevels(test)
table(test$Driver.type, test$RL.Lat.Elev, test$RL) #least data is for latitude x warm (30 & 23 points)
#get sample size
n <- allRL[allRL$Method=='Field experiment',]; dim(n)
coolvwarm.LEf <- glmer(Driver.supported01 ~ Driver.type*RL*RL.Lat.Elev + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRL[allRL$Method=='Field experiment',]) 
coolvwarm.LEf.no3X <- glmer(Driver.supported01 ~ Driver.type*RL + Driver.type*RL.Lat.Elev + RL*RL.Lat.Elev + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRL[allRL$Method=='Field experiment',]) 
      anova(coolvwarm.LEf, coolvwarm.LEf.no3X, test='Chisq') #3 way interaction is NS
      Anova(coolvwarm.LEf.no3X, type='III')
#__support by factor type --------------------------------- ---
n <- allRLsum[allRLsum$inclFexp=='Y',]; dim(n)
coolvwarm.LEf <- glmer(Drivertype.supported01 ~ Driver.type*RL*RL.Lat.Elev + nDrivers.tested +(1|Sp.name) +(1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRLsum[allRLsum$inclFexp=='Y',]) 
coolvwarm.LEf.no3X <- glmer(Drivertype.supported01 ~ Driver.type*RL + Driver.type*RL.Lat.Elev + RL*RL.Lat.Elev + nDrivers.tested +(1|Sp.name) +(1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRLsum[allRLsum$inclFexp=='Y',]) 
      anova(coolvwarm.LEf, coolvwarm.LEf.no3X, test='Chisq') #3 way interaction is NS
      Anova(coolvwarm.LEf.no3X, type='III')
      

#__Land v marine (all studies) --------------------------
#__support by factor --------------------------------- ---
table(allRL$Driver.type, allRL$Environment, allRL$RL) #have little data for freshwater 
#get sample size - as excluded freshwater
n <- allRL[allRL$Environment!='freshwater',]; n <- droplevels(n); dim(n)
coolvwarm.env <- glmer(Driver.supported01 ~ Driver.type*RL*Environment + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRL[allRL$Environment!='freshwater',]) 
coolvwarm.env.no3X <- glmer(Driver.supported01 ~ Driver.type*RL + Driver.type*Environment + RL*Environment + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRL[allRL$Environment!='freshwater',]) 
      anova(coolvwarm.env, coolvwarm.env.no3X, test='Chisq') #3 way interaction is NS
      Anova(coolvwarm.env.no3X, type='III')
#__support by factor type --------------------------------- ---
table(allRLsum$Driver.type, allRLsum$Environment, allRLsum$RL) #have little data for freshwater 
#get sample size - as excluded freshwater
n <- allRLsum[allRLsum$Environment!='freshwater',]; n <- droplevels(n); dim(n)
coolvwarmtype.env <- glmer(Drivertype.supported01 ~ Driver.type*RL*Environment + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRLsum[allRLsum$Environment!='freshwater',]) 
coolvwarmtype.env.no3X <- glmer(Drivertype.supported01 ~ Driver.type*RL + Driver.type*Environment + RL*Environment + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRLsum[allRLsum$Environment!='freshwater',]) 
      anova(coolvwarmtype.env, coolvwarmtype.env.no3X, test='Chisq') #3 way interaction is NS
      Anova(coolvwarm.env.no3X, type='III')

#__Land v marine (field exp - decided to exclude as not many marine field exp) --------------------------
test <- allRL[allRL$Method=='Field experiment',]; test <- droplevels(test)
table(test$Driver.type, test$Environment, test$RL) #pretty low for marine 

      
## Q2) Which Drivers? (comp, pred, temp etc) ###############################

#look at categories            
table(allRL$Driver.examined, allRL$RL) #Other isnt meaningful to compare to, so exclude these 
table(allRL$Driver.examined, allRL$Driver.type)
#what is category 'other'?
test <- allRL[allRL$Driver.examined=='Other',]; test <- droplevels(test)
table(test$Driver.descrip, test$Driver.type)

bwplot(Driver.supported01 ~ RL | Driver.examined, data=allRL, subset=Driver.examined!='Other') 
bwplot(Driver.supported01 ~ RL | Driver.examined, data=allRL) #'Other' isnt meanginful in this context so exclude

#random factors for sp & doi 
coolvwarm.driver <- glmer(Driver.supported01 ~ Driver.examined*RL + (1|Sp.name) + (1|doi), #no warnings!
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRL, subset=Driver.examined!='Other') 
    #reduced model to test significance
    coolvwarm.drivernoX <- glmer(Driver.supported01 ~ Driver.examined + RL + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRL, subset=Driver.examined!='Other') 
    anova(coolvwarm.driver, coolvwarm.drivernoX, test='Chisq') #interaction signif
    lsmeans(coolvwarm.driver, pairwise ~ RL | Driver.examined, transform='response')
    #more support for one RL vs the other?
      #Biogenic habitat: NS
      #Climate: NS
      #Competition: NS (p=0.0826)
      #Host/Food: NS
      #Parasitism/Disease: warm > cool 
      #Precip: NS
      #Pred/Herbiv: warm > cool
      #Soil: NS
      #Temp: cool > warm
    lsmeans(coolvwarm.driver, pairwise ~ Driver.examined | RL, transform='response') 
    #at coolRL: Temp a  Climate ab  BioHab ab  Soil ab  Precip b  Food/Host abc  Comp b  PredH c  Pathogen c
        #NB: Temp almost > Soil (P = 0.061)
        # only Temp supported > 50% of time
    #at warmRL: no signif contrasts (pred/herbiv almost bigger than soil 0.0519)
        # BioHab, Predation/Herb, Temp supported > 50% of time 
    visreg(coolvwarm.driver, xvar='RL', by='Driver.examined', type='conditional')
    visreg(coolvwarm.driver, xvar='Driver.examined', by='RL', type='conditional')
    #get sample size
    n <- allRL[allRL$Driver.examined!='Other',]; n <- droplevels(n); dim(n)
    table(n$Driver.examined, n$RL) #will need for Fig 3

#__which drivers Field experiments -------------------------------
summary(allRL)
Fdat <- allRL[allRL$Method=='Field experiment',]; Fdat <- droplevels(Fdat)
table(Fdat$Driver.examined, Fdat$RL)
#cant do Host/Food (too few) or biogenic habitat (none)
#get sample size
n <- allRL[allRL$Driver.examined!='Other' & allRL$Driver.examined!='Host/Food availability',]; n <- droplevels(n); dim(n)
table(n$Driver.examined, n$RL)
coolvwarm.driverF <- glmer(Driver.supported01 ~ Driver.examined*RL + (1|Sp.name) + (1|doi), #singular
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=Fdat, subset=Driver.examined!='Other' & Driver.examined!='Host/Food availability')
    #reduced model to test significance
    coolvwarm.driverFnoX <- glmer(Driver.supported01 ~ Driver.examined + RL + (1|Sp.name) + (1|doi),
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=Fdat, subset=Driver.examined!='Other' & Driver.examined!='Host/Food availability')
    anova(coolvwarm.driverF, coolvwarm.driverFnoX, test='Chisq') #interaction signif
    visreg(coolvwarm.driverF, xvar='RL', by='Driver.examined', type='conditional')
    lsmeans(coolvwarm.driverF, pairwise ~ RL | Driver.examined, transform='response')
    #more support for one RL vs the other?
      #Climate: cool > warm
      #Competition: warm > cool
      #Parasitism/Disease: warm > cool
      #Precip: NS
      #Pred/Herbiv: warm > cool
      #Soil: NS
      #Temp: NS 
    lsmeans(coolvwarm.driverF, pairwise ~ Driver.examined | RL, transform='response')
    #F at cool:   Temp a Clim a Soil ab Precip ab Comp ab  Pred/Herb b Paras b
    #F at warmRL: Pred a Comp a  Precip a  Paras ab Temp ab   Soil b  Clim b
    visreg(coolvwarm.driverF, xvar='RL', by='Driver.examined', type='conditional')
    visreg(coolvwarm.driverF, xvar='Driver.examined', by='RL', type='conditional')
    

    
## Q3) Effect of Latitude ############################################################
#_byDriver x absLatitude--------------------------------
#use data with locations merged on
summary(RLloc) #already got rid of NAs for absLat
bwplot(Driver.supported01 ~ Latzone | Driver.type*RL, data=RLloc)
bwplot(sumDriver.support ~ Latzone | Driver.type*RL, data=RLlocsum)

#get sample size
dim(RLloc) #n 1931
#can use 'lsmeans' to compare slopes in interactions but need different call from emmeans package (etrends)
#https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html#covariates
coolvwarm.lat <- glmer(Driver.supported01 ~ Driver.type*RL*absLat + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=RLloc) #singular
coolvwarm.lat.no3X <- glmer(Driver.supported01 ~ Driver.type*RL + Driver.type*absLat + RL*absLat + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=RLloc) #singular
      anova(coolvwarm.lat, coolvwarm.lat.no3X, test='Chisq') #3way interaction signif (P=0.0292)
      #conditional plots all terms must be specified. ie need to plot for either warm or cool
      visreg(coolvwarm.lat, xvar='absLat', by='Driver.type', cond=list(RL='cool'), type='conditional', main='coolRLs')
      visreg(coolvwarm.lat, xvar='absLat', by='Driver.type', cond=list(RL='warm'), type='conditional', main='warm RLs')
      #look at significance
      emtrends(coolvwarm.lat, pairwise ~ RL | Driver.type, var='absLat')
      #for biotic factors, lat trends at cool and warm limits are NOT sig dif (P=0.0586) (results sensitive)
      #for abiotic factors, lat trends at cool and warm limits are NOT sig dif (P>0.3)
      emtrends(coolvwarm.lat, pairwise ~ Driver.type | RL, var='absLat')
      #at cool limits, the lat trend for biotic and abiotic factors is same (P=0.055) (results sensitive)
      #at warm limits, the lat trend for biotic and abiotic factors is same (P>0.3)
      #All four individual trends are NS (CI overlap 0)

#trend seems to be real (see below), so what is causing them?  Are their biases in the data that might be giving false signals?
#Data exploration: plot data not collapsed across hemispheres (ie Lat, not absolute latitude) 
summary(RLloc)
xyplot(Latitude ~ RL, data=RLloc)
#data have very similar distribution between the two types of range limits, so unlikely that bias in more/less southern hemisphere data causing the difference
#histogram(Latitude | RL, data=RLloc) #this one should work according to lattice help menu but doesnt
hist(loc.cool$Latitude, xlim=c(-76,76))
hist(loc.warm$Latitude, xlim=c(-76,76))
#split again by biotic/abiotic factors
xyplot(Latitude ~ RL | Driver.type, data=RLloc)
summary(RLloc)
#which driver examined most often?
table(RLloc$Driver.examined, RLloc$Driver.type)
#look at distribution of climate tests
summary(RLloc)
hist(RLloc$Latitude[RLloc$Driver.examined=='Climate' | RLloc$Driver.examined=='Temperature'])
hist(RLloc$Latitude[RLloc$Driver.examined=='Climate' | RLloc$Driver.examined=='Temperature'])
xyplot(Latitude ~ as.factor(RL), data=RLloc[RLloc$Driver.type=='Climate' | RLloc$Driver.type=='Temperature',])

#check effect of elevation?
#cant as dont have the data. also for all but field experiments would be pretty hard to get an 'average' elevation

      
#__SI.5 absLat field exp only - since those are the results we have most confidence in ------
summary(RLloc)
#get sample size
n <- RLloc[RLloc$Method=='Field experiment',]; n <- droplevels(n); dim(n) #545
coolvwarm.latF <- glmer(Driver.supported01 ~ Driver.type*RL*absLat + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=RLloc[RLloc$Method=='Field experiment',]) #singular
coolvwarm.latF.no3X <- glmer(Driver.supported01 ~ Driver.type*RL + Driver.type*absLat + RL*absLat + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=RLloc[RLloc$Method=='Field experiment',]) #does not like model: max grad 0.4
      anova(coolvwarm.latF, coolvwarm.latF.no3X, test='Chisq') #3way interaction signif
      emtrends(coolvwarm.latF, pairwise ~ RL | Driver.type, var='absLat')
      #for biotic factors, lat trends at cool vs warm limits NS
      #for abiotic factors, lat trends at cool vs warm limits are sig dif
      emtrends(coolvwarm.latF, pairwise ~ Driver.type | RL, var='absLat')
      #at cool limits, lat trend for biotic and abiotic factors is NOT diff (P=0.058)
      #at warm limits, lat trend for biotic and abiotic factors is diff (P=0.0057)
      #none of the individual trends is signif (CI overlap zero for all four)
      #conditional plots all terms must be specified. ie need to plot for either warm or cool
      visreg(coolvwarm.latF, xvar='absLat', by='RL', cond=list(Driver.type='biotic'), type='conditional', overlay=T, scale='response', ylim=c(0,1),
             ylab='Biotic support fieldexp') #similar to all studies   
      visreg(coolvwarm.latF, xvar='absLat', by='RL', cond=list(Driver.type='abiotic'), type='conditional', overlay=T, scale='response', ylim=c(0,1),
             ylab='Abiotic support fieldexp') #now both decline, and decline for warm limits much steeper   

      
#__SI.5 absLat Land - expect stronger effect on land than marine  ------
levels(RLloc$Environment)
#get sample size
n <- RLloc[RLloc$Environment=='terrestrial',]; n <- droplevels(n); dim(n) 
coolvwarm.latLand <- glmer(Driver.supported01 ~ Driver.type*RL*absLat + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=RLloc[RLloc$Environment=='terrestrial',]) #singular
coolvwarm.latLand.no3X <- glmer(Driver.supported01 ~ Driver.type*RL + Driver.type*absLat + RL*absLat + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=RLloc[RLloc$Environment=='terrestrial',]) #singular
      anova(coolvwarm.latLand, coolvwarm.latLand.no3X, test='Chisq') #3way interaction NS (changed March 26)
      emtrends(coolvwarm.latLand, ~ RL | Driver.type, var='absLat')
      #none of the four individual trends in significant
      Anova(coolvwarm.latLand.no3X, type='III')
      #conditional plots all terms must be specified. ie need to plot for either warm or cool
      visreg(coolvwarm.latLand, xvar='absLat', by='RL', cond=list(Driver.type='biotic'), type='conditional', overlay=T, scale='response', 
             ylim=c(0,1), ylab='Biotic support Land') #v similar to all studies   
      visreg(coolvwarm.latLand, xvar='absLat', by='RL', cond=list(Driver.type='abiotic'), type='conditional', overlay=T, scale='response', 
             ylim=c(0,1), ylab='Abiotic support Land') #v similar to all studies   
      visreg(coolvwarm.latLand, xvar='absLat', by='RL', cond=list(Driver.type='biotic'), type='conditional', overlay=T,
             ylab='Biotic support Land') #v similar to all studies   
      visreg(coolvwarm.latLand, xvar='absLat', by='RL', cond=list(Driver.type='abiotic'), type='conditional', overlay=T,
             ylab='Abiotic support Land') #v similar to all studies 
#check latitude x factor type interaction (original prediction)
coolvwarm.latLand.no2LatX <- glmer(Driver.supported01 ~ Driver.type*RL + absLat + (1|Sp.name) + (1|doi),
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=RLloc[RLloc$Environment=='terrestrial',]) #singular
      anova(coolvwarm.latLand.no2LatX, coolvwarm.latLand.no3X, test='Chisq') #NS
      Anova(coolvwarm.latLand.no2LatX, type='III')

      
#__SI.5 absLat Marine - expect weaker effect  ------
#get sample size
n <- RLloc[RLloc$Environment=='marine',]; n <- droplevels(n); dim(n) #368
summary(n$absLat) #min lat = 8.39
coolvwarm.latMarine <- glmer(Driver.supported01 ~ Driver.type*RL*absLat + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=RLloc[RLloc$Environment=='marine',]) #singular
coolvwarm.latMarine.no3X <- glmer(Driver.supported01 ~ Driver.type*RL + Driver.type*absLat + RL*absLat + (1|Sp.name) + (1|doi), 
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=RLloc[RLloc$Environment=='marine',]) #no warnings
      anova(coolvwarm.latMarine, coolvwarm.latMarine.no3X, test='Chisq') #3way interaction NOT signif
      Anova(coolvwarm.latMarine.no3X, type='III')
      emtrends(coolvwarm.latMarine, ~ RL | Driver.type, var='absLat')
      #all 4 trends NS
# #test importance of Driver.type x absLatitude (this is what we predicted - ie biotic factors become genrally more important as Lat decreases and abiotic factors become generally more important as Lat increases)
# coolvwarm.latMarine.noLatX <- glmer(Driver.supported01 ~ Driver.type*RL + absLat + (1|Sp.name) + (1|doi), 
#                           family=binomial,
#                           control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
#                           data=RLloc[RLloc$Environment=='marine',]) #no warnings
#       anova(coolvwarm.latMarine.noLatX, coolvwarm.latMarine.no3X, test='Chisq') #2way interactions NOT signif
#       Anova(coolvwarm.latMarine.noLatX, type='III') #now nothing but Driver.type is significant
    
            
            
#FIGURES ######################################################################

#Fig.1 Predictions ------------------------------
#_Fig1 grouped by driver (USE) ------- -
a=0.6; b=0.92 #alpha = opacity, so lower #s = more transparent
colC <- 'dodgerblue3'; colW <- 'red3'
colCbio <- colCabio <- adjustcolor('dodgerblue', alpha=a)
colWbio <- colWabio <- adjustcolor('red3', alpha=a)
coll = "dimgrey" #abline colour


quartz(height=6, width=9) #for Anna
layout(matrix(c(1,2,3,4,5,6,7,8), nrow=2, ncol=4, byrow=T), widths=c(2,2,.5,2), heights=c(1,1)); 
layout.show(8)

xbio=2; xabio=6; #centre of warm vs cool data (where tick marks go)
sep=.75 #separation btwn biotic & abiotic within RL
xbioC=xbio-sep; xbioW=xbio+sep; xabioC=xabio-sep; xabioW=xabio+sep
xmax=7.9; xlim=c(.6,xmax)
lwd=2 #prediction lines
RW=sep-.2 #half the width of each bar
cxPred=.7 #Prediction labels
cxCW = .8 #cool/warm labels (xaxis)
x=c(1, 100); y=c(-2, 100) #dummy data to create plot
marPred=c(2,2,3,1); marData=c(2,3,3,0) 

#par(oma=c(3,3,1,2), mar=marPred, bty='l')
par(oma=c(3,3,1.5,2), mar=marPred, bty='l')
#A) Topleft = Prediction 1: biotic more important at warm limits than cool limits
plot(x=x, y=y, type="n", xlim=xlim, xaxt='n', yaxt='n', xlab=NA, ylab=NA)
    axis(side=1, at=c(xbioC, xbioW, xabioC, xabioW), labels = NA) 
    adjAB=.12
    mtext(side=1, at=c(xbioC, xabioC)-adjAB, line=1, text=c('cool', 'cool'), cex=cxCW, col=colC)
    mtext(side=1, at=c(xbioW, xabioW)+adjAB, line=1, text=c('warm', 'warm'), cex=cxCW, col=colW)
#Cool biotic
rect(xleft=xbioC-RW, xright=xbioC+RW, ybottom=0, ytop=40, col=colCbio)
#Warm biotic
rect(xleft=xbioW-RW, xright=xbioW+RW, ybottom=0, ytop=70, col = colWbio)
    cxPan=.95; lPanlet= -1.2; lPan= -1.2; lPan=.4 #1st lPan = inside plot margins, 2nd = above plot
    mtext(text=c('Biotic','Abiotic'), side=3, line=lPan, at=c(xbio, xabio), cex=cxPan) #within panels
    lPl=80;  lines(x=c(xbioC,xbioW), y=c(lPl,lPl), lw=lwd, col='grey48') #prediction line
    lPtxt=-3.1; mtext(text='Prediction 1', side=3, line=lPtxt, at=2.1, cex=cxPred)
    mtext(text ="A", side=3, line=lPanlet, cex=cxPan, adj=-.1)

#B) Topmiddle = Hypothesis 2: biotic more important than abiotic at warm limits
plot(x=x, y=y, type="n", las=1, xlim=xlim, xaxt='n', yaxt='n', xlab=NA, ylab=NA)
    axis(side=1, at=c(xbioC, xbioW, xabioC, xabioW), labels =NA) 
    mtext(side=1, at=c(xbioC, xabioC)-adjAB, line=1, text=c('cool', 'cool'), cex=cxCW, col=colC)
    mtext(side=1, at=c(xbioW, xabioW)+adjAB, line=1, text=c('warm', 'warm'), cex=cxCW, col=colW)
#Warm abiotic
rect(xleft=xabioW-RW, xright=xabioW+RW, ybottom=0, ytop=40, col = colWabio)
#Warm biotic
rect(xleft=xbioW-RW, xright=xbioW+RW, ybottom=0, ytop=70, col=colWbio)
    mtext(text=c('Biotic','Abiotic'), side=3, line=lPan, at=c(xbio, xabio), cex=cxPan)
    lines(x=c(xbioW,xabioW), y=c(lPl,lPl), lw=lwd, col=colW)
    mtext(text='Prediction 2', side=3, line=lPtxt, at=4.8, cex=cxPred)
    mtext(text="B", side=3, line=lPanlet, cex=cxPan, adj=-.1)

#empty space
plot(x=x, y=y, type = "n", bty="n", xlim=xlim, xaxt='n', yaxt='n', xlab=NA, ylab=NA)

#E) Top right = Hargreaves study (elev + latitude combined)
par(mar=marData)
plot(x=x, y=y, type = "n", las=1, xlim=xlim, xaxt='n', yaxt='n', xlab=NA, ylab=NA)
      ylab=c('0','0.2','0.4','0.6','0.8','1.0'); axis(side=2, at=seq(0,100,by=20), labels=ylab, las=1)
      abline(h=50, col=coll)
      axis(side=1, at=c(xbioC, xbioW, xabioC, xabioW), labels=NA)
      mtext(side=1, at=c(xbioC, xabioC)-adjAB, line=1, text=c('cool', 'cool'), cex=cxCW, col=colC)
      mtext(side=1, at=c(xbioW, xabioW)+adjAB, line=1, text=c('warm', 'warm'), cex=cxCW, col=colW)
#Cool biotic: 18% of 38 elev = 7 + 0 lat = 7, denom = 38+8=46, %=7/46 = 15%
rect(xleft=xbioC-RW, xright=xbioC+RW, ybottom=0, ytop=15, col=colCbio, border="black")
#Warm biotic: 55% of 20 elev = 11 + 1 of 1 lat = 12/21= 57%
rect(xleft=xbioW-RW, xright=xbioW+RW, ybottom=0, ytop=57, col=colWbio)
#mtext(text='Hargreaves et al. 2014', side=3, line= lPanlet, at=4, cex=cxPan)
#mtext(text='Tested & supported Prediction 1', side=3, line=lconcl, at=4.2, cex=cxPred)
mtext(text=c('Biotic','Abiotic'), side=3, line=lPan, at=c(xbio, xabio), cex=cxPan)
lstudy=lPtxt+.8; mtext(text='Hargreaves et al. 2014', side=3, line=lstudy, at=4, cex=cxPred)
lconcl=lstudy-1.1; mtext(text='tested & supported Prediction 1', side=3, line=lconcl, at=4.2, cex=cxPred)
mtext(text="E", side=3, line=lPanlet, cex=cxPan, adj=-.21)
#lylab=2.2; mtext("% limits involving biotic factors", side=2, line=lylab, adj=.2, cex=cxAB) 
lylab=2.3; mtext("Prop limits involving biotic factors", at=15, side=2, line=lylab, adj=.2, cex=cxAB) 

#C) Bottomleft = Prediction 3: 
par(mar=marPred)
plot(x=x, y=y, type = "n", las=1, xlim=xlim, xaxt='n', yaxt='n', xlab=NA, ylab=NA)
    #axis(side=2, at=50, labels='50', tick=T, las=2)
    axis(side=2, at=50, labels='0.5', tick=T, las=2)
    axis(side=1, at=c(xbioC, xbioW, xabioC, xabioW), labels =NA) 
    mtext(side=1, at=c(xbioC, xabioC)-.1, line=1, text=c('cool', 'cool'), cex=cxAB, col=colC)
    mtext(side=1, at=c(xbioW, xabioW)+.1, line=1, text=c('warm', 'warm'), cex=cxAB, col=colW)
    lines(x=c(0,xmax-.1), y=c(50, 50), col=coll)
#Warm biotic
rect(xleft=xbioW-RW, xright=xbioW+RW, ybottom=0, ytop=70, col=colWbio)
lPl=xbioW-RW-.4; lines(x=c(lPl,lPl), y=c(50, 70), lw=lwd, col=colW)
#text(x=lPl-.4, y=50, labels='Prediction 3', cex=cxPred+.4, srt=90, adj=0) #above abline
text(x=lPl-.4, y=50, labels='Prediction 3', cex=cxPred+.4, srt=90, adj=.2) #centered over difference
mtext(text="C", side=3, line=lPanlet, cex=cxPan, adj=-.12)

#D) Bottom middle = all predictions combined
plot(x=x, y=y, type = "n", las=1, xlim=xlim, xaxt='n', yaxt='n', xlab=NA, ylab=NA)
    axis(side=2, at=50, labels='0.5', tick=T, las=2)
    axis(side=1, at=c(xbioC, xbioW, xabioC, xabioW), labels =NA) 
    mtext(side=1, at=c(xbioC, xabioC)-.1, line=1, text=c('cool', 'cool'), cex=cxAB, col=colC)
    mtext(side=1, at=c(xbioW, xabioW)+.1, line=1, text=c('warm', 'warm'), cex=cxAB, col=colW)
    lines(x=c(0,7.4), y=c(50, 50), col=coll) #reference line
#Cool abiotic
rect(xleft=xabioC-RW, xright=xabioC+RW, ybottom=0, ytop=70, col=colCabio)
#Cool biotic
rect(xleft=xbioC-RW, xright=xbioC+RW, ybottom=0, ytop=40, col=colCbio)
#Warm abiotic
rect(xleft=xabioW-RW, xright=xabioW+RW, ybottom=0, ytop=40, col=colWabio)
#Warm biotic
rect(xleft=xbioW-RW, xright=xbioW+RW, ybottom=0, ytop=70, col=colWbio)
mtext(text='All predictions supported', side=3, line=lPtxt, at=4.1, cex=cxPred)
mtext(text="D", side=3, line=lPanlet, cex=cxPan, adj=-.12)

#empty space
plot(x=x, y=y, type="n", bty="n", xlim=xlim, xaxt='n', yaxt='n', xlab=NA, ylab=NA)

#F) Bottomright = Cahill results
par(mar=marData)
plot(x=x, y=y, type = "n", las=1, xaxt='n', yaxt='n', xlim=xlim, xlab=NA, ylab=NA)
    axis(side=2, at=seq(0,100,by=20), labels=c('0','0.2','0.4','0.6','0.8','1.0'), las=1)
    axis(side=1, at=c(xbioC, xbioW, xabioC, xabioW), labels=NA)
    mtext(side=1, at=c(xbioC, xabioC)-.1, line=1, text=c('cool', 'cool'), cex=cxAB)
    mtext(side=1, at=c(xbioW, xabioW)+.1, line=1, text=c('warm', 'warm'), cex=cxAB)
    abline(h=50, col=coll)
#Warm abiotic
rect(xleft=xabioW-RW, xright=xabioW+RW, ybottom=0, ytop=79.2, col=colWabio)
#Warm biotic
rect(xleft=xbioW-RW, xright=xbioW+RW, ybottom=0, ytop=59, col=colWbio)
#add cool & warm labels to top
#mtext(text='Cahill et al. 2014', side=3, line= lPanlet, at=4, cex=cxPan)
#mtext(text = '?', side = 3, line = -13, at = xbio, cex = 2.5)
#mtext(text='Tested & refuted Prediction 2', side=3, line=lPtxt, at=4.2, cex=cxPred)
mtext(text='Cahill et al. 2014', side=3, line=lstudy, at=4, cex=cxPred)
mtext(text='tested & refuted Prediction 2', side=3, line=lconcl, at=4.2, cex=cxPred)
mtext(text="F", side=3, line=lPanlet, cex=cxPan, adj=-.21)
mtext("Proportion support by factor", side=2, line=lylab, adj=0.5, cex=cxAB) 

#add overall x and y labels 
cxtit <- 1; lbot=3.4; ltop=23.5
mtext("Type of range limit", side=1, line=lbot, adj=.5, cex=cxtit) #bigger adj -> left
mtext("Type of range limit", side=1, line=lbot, adj=-5.5, cex=cxtit)
#mtext("% Support for biotic or abiotic factors", side=2, line=48, adj=-1, cex=cxtit)#bigger line #s -> left
mtext("% Support for biotic or abiotic factors limiting cool vs. warm range limits", side=2, line=48, adj=0, cex=cxtit) #bigger line #s -> left
mtext("Empirical results", side=3, line=ltop, adj=.4, cex=cxtit+.2)
mtext("Predictions", side=3, line=ltop, adj=-3.5, cex=cxtit+.2) #bigger adj #s move further left


#Fig.2 MainResults RL x factor type, all vs field ---------------------------------------------
#__data -----------------------
#_panel 1 - byDriver by factor (Bio v Abio) (define which model want to use for plotting) -- 
  byDriv <- coolvwarm.bydriver; byDriv
  byDrivlsm <- data.frame(summary(lsmeans(byDriv, ~ Driver.type | RL, type='response'))); 
  byDrivlsm 
  lsmeans(byDriv, pairwise ~ RL | Driver.type, type='response') #def says contrast is signif
  summary(allRL)
  visreg(byDriv, xvar='RL', by="Driver.type", partial=T, layout=c(2,1))
  visreg(byDriv, xvar='RL', by="Driver.type",  partial=T, scale='response', layout=c(2,1))
  vFitbyD <- visreg(byDriv, xvar='RL', by="Driver.type",  partial=T, scale='response', layout=c(2,1))
  vFitbyD$fit$ABnum <- factor(ifelse(vFitbyD$fit$Driver.type=='biotic', 1, 2)); vFitbyD$fit #hack to get numeric
  adj=.15
  vFitbyD$fit$order <- ifelse(vFitbyD$fit$Driver.type=='biotic' & vFitbyD$fit$RL=='cool', 1-adj,
                       ifelse(vFitbyD$fit$Driver.type=='biotic' & vFitbyD$fit$RL=='warm', 2-adj,
                       ifelse(vFitbyD$fit$Driver.type=='abiotic' & vFitbyD$fit$RL=='cool', 3+adj, 4+adj)))
  head(vFitbyD$fit) #fitted means
  head(vFitbyD$res) #residuals
  vFitbyD$res$RLnum <- ifelse(vFitbyD$res$RL=='cool',1,2) #hack to get numeric xaxis so can jitter
  vFitbyD$res$Drivern <- ifelse(vFitbyD$res$Driver.type=='biotic',1,2) 
  vFitbyD$res$order <- ifelse(vFitbyD$res$Driver.type=='biotic' & vFitbyD$res$RL=='cool', 1-adj,
                       ifelse(vFitbyD$res$Driver.type=='biotic' & vFitbyD$res$RL=='warm', 2-adj,
                       ifelse(vFitbyD$res$Driver.type=='abiotic' & vFitbyD$res$RL=='cool', 3+adj, 4+adj)))

#_data panel 2 - by driver (field experiments) -------------- --- 
  byDrivF <- coolvwarm.fexp; byDrivF
  byDrivFlsm <- data.frame(summary(lsmeans(byDrivF, ~ Driver.type | RL, type='response'))); 
  byDrivFlsm 
  visreg(byDrivF, xvar='RL', by="Driver.type", partial=T, layout=c(2,1))
  visreg(byDrivF, xvar='RL', by="Driver.type",  partial=T, scale='response', layout=c(2,1))
  vFitbyDF <- visreg(byDrivF, xvar='RL', by="Driver.type",  partial=T, scale='response', layout=c(2,1))
  vFitbyDF$fit$ABnum <- factor(ifelse(vFitbyDF$fit$Driver.type=='biotic', 1, 2)); vFitbyDF$fit #hack to get numeric
  vFitbyDF$fit$order <- ifelse(vFitbyDF$fit$Driver.type=='biotic' & vFitbyDF$fit$RL=='cool', 1-adj,
                       ifelse(vFitbyDF$fit$Driver.type=='biotic' & vFitbyDF$fit$RL=='warm', 2-adj,
                       ifelse(vFitbyDF$fit$Driver.type=='abiotic' & vFitbyDF$fit$RL=='cool', 3+adj, 4+adj)))
  head(vFitbyDF$fit) #fitted means
  head(vFitbyDF$res) #residuals
  vFitbyDF$res$RLnum <- ifelse(vFitbyDF$res$RL=='cool',1,2) #hack to get numeric xaxis so can jitter
  vFitbyDF$res$Drivern <- ifelse(vFitbyDF$res$Driver.type=='biotic',1,2) 
  vFitbyDF$res$order <- ifelse(vFitbyDF$res$Driver.type=='biotic' & vFitbyDF$res$RL=='cool', 1-adj,
                       ifelse(vFitbyDF$res$Driver.type=='biotic' & vFitbyDF$res$RL=='warm', 2-adj,
                       ifelse(vFitbyDF$res$Driver.type=='abiotic' & vFitbyDF$res$RL=='cool', 3+adj, 4+adj)))
#__plot by driver all vs field -----------------
quartz(height=4, width=8)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow = T)); #layout.show(2)
    cxp=.5; 
    #xl=1.8; yl=2.6
    cxaxt=1.1; cxsig=.9 #size axis tick labels & model significance
    colCool='dodgerblue'; colWarm='red3'; aCI=.6; ap=.4 #lower values of alpha = more transparent
    colCoolCI <- adjustcolor(colCool, alpha.f=aCI); colWarmCI <- adjustcolor(colWarm, alpha.f=aCI) 
    colCoolp <- adjustcolor(colCool, alpha.f=ap);   colWarmp <- adjustcolor(colWarm, alpha.f=ap) 
    xct=1; r=.4 #defines width of rectangle
    xlim=c(0.2,4.6); x=c(1-adj,2-adj,3+adj,4+adj); xtxt=c('cool','warm','cool','warm')
    ymin=0; ylim=c(ymin,1);
    xtl=.6; lwd=1.1

par(oma=c(2,3.5,1.5,0), mar=c(1.3,1,1.5,.5), bty='l') 
#_A) all studies ------------ ---
plot(vFitbyD$fit$visregFit ~ vFitbyD$fit$order, xlim=xlim, ylim=ylim, las=1, pch='—', cex=2.5, xaxt='n') 
    #axis(side=2, at=seq(0,1,by=.2), labels=NA) #makes more cluttered
    axis(side=1, at=x, labels=NA)
    mtext(side=1, line=xtl, at=x, text=xtxt, cex=cxaxt, col=c(colC,colW,colC,colW)) #colours defined for fig1
    abline(h=.5)
    #mtext(at=xat, side=3, line=, text=c('Biotic factors','Abiotic factors'), cex=1.1)
    # biotic cool
    resBC <- vFitbyD$res$visregRes[vFitbyD$res$RL=='cool' & vFitbyD$res$Driver.type=='biotic']; head(resBC)
    points(jitter(x=rep(1-adj,length=length(resBC)),factor=20), 
           y=resBC, cex=cxp, pch=16, col=colCoolp)
    byDrivlsm; 
    rect(xleft=1-adj-r, xright=1-adj+r, ytop=byDrivlsm$asymp.UCL[2], ybottom=byDrivlsm$asymp.LCL[2],col=colCoolCI, border=NA)
    # biotic warm
    resBW <- vFitbyD$res$visregRes[vFitbyD$res$RL=='warm' & vFitbyD$res$Driver.type=='biotic']; head(resBW)
    points(jitter(x=rep(2-adj, length=length(resBW)), factor=9), 
                  y=resBW, cex=cxp, pch=16, col=colWarmp)
    rect(xleft=2-adj-r, xright=2-adj+r, ytop=byDrivlsm$asymp.UCL[4], ybottom=byDrivlsm$asymp.LCL[4], col=colWarmCI, border=NA)
    # abiotic cool
    resAC <- vFitbyD$res$visregRes[vFitbyD$res$RL=='cool' & vFitbyD$res$Driver.type=='abiotic']; head(resAC)
    points(jitter(x=rep(3+adj, length=length(resAC)), factor=5),
                  y=resAC, cex=cxp, pch=16, col=colCoolp)
    rect(xleft=3+adj-r, xright=3+adj+r, ytop=byDrivlsm$asymp.UCL[1], ybottom=byDrivlsm$asymp.LCL[1], col=colCoolCI, border=NA)
    #abiotic warm
    resAW <- vFitbyD$res$visregRes[vFitbyD$res$RL=='warm' & vFitbyD$res$Driver.type=='abiotic']; head(resAW)
    points(jitter(x=rep(4+adj, length=length(resAW)), factor=4),
                  y=resAW, cex=cxp, pch=16, col=colWarmp)
    rect(xleft=4+adj-r, xright=4+adj+r, ytop=byDrivlsm$asymp.UCL[3], ybottom=byDrivlsm$asymp.LCL[3], col=colWarmCI, border=NA)
    par(new=T); plot(vFitbyD$fit$visregFit ~ vFitbyD$fit$order, xlim=xlim, ylim=ylim, 
                       xaxt='n', yaxt='n', pch='—', cex=2.5)
    #significance
      #text(x=.1, y=ymin, labels='factor type x limit type P < 0.001', pos=4, cex=cxsig)
      #text(x=xat, y=c(rep(.98,2)), labels=c('*',''), cex=2.4)
      lP1=.96; lP2=.07;
      lines(x=x[1:2], y=c(lP1,lP1), lw=lwd, col='grey48') #prediction 1 lines: C < W for biotic
      lines(x=x[3:4], y=c(lP1,lP1), lw=lwd, col='grey48') #prediction 1 lines: C > W for abiotic (lsm says P=0.02 even though CI overlap)
      lines(x=c(x[1],x[3]), y=c(lP2,lP2), lw=lwd, col=colC) #prediction 2 lines: C < W for biotic
      lPlet=1.6; 
      mtext(at=.8, side=3, line=lPlet, text='A) All studies', cex=1.2) #left justified
      #mtext(at=2.5, side=3, line=lPlet, text='A) All studies', cex=1.2) #centered
      xat=c(1.5-adj,3.5+adj+adj); cxaxl=1.1
      mtext(at=xat, side=3, line=, text=c('Biotic','Abiotic'), cex=1.1)
      mtext(at=.5, side=2, line=3, text='Probability factor contributes to range limit', cex=cxaxl)
      #mtext(at=2.3, side=1, line=2, text='Type of range limit', cex=cxaxl)

#_B) field experiments ------ --      
plot(vFitbyDF$fit$visregFit ~ vFitbyDF$fit$order, xlim=xlim, ylim=ylim, las=1, pch='—', cex=2.5, xaxt='n', yaxt='n') 
    axis(side=1, at=x, labels=NA)
    mtext(side=1, line=xtl, at=x, text=xtxt, cex=cxaxt, col=c(colC,colW,colC,colW))
    abline(h=.5)
    # biotic cool
    resBCF <- vFitbyDF$res$visregRes[vFitbyDF$res$RL=='cool' & vFitbyDF$res$Driver.type=='biotic']; head(resBCF)
    points(jitter(x=rep(1-adj,length=length(resBCF)),factor=20),
                  y=resBCF, cex=cxp, pch=16, col=colCoolp)
    rect(xleft=1-adj-r, xright=1-adj+r, ytop=byDrivFlsm$asymp.UCL[2], ybottom=byDrivFlsm$asymp.LCL[2], 
         col=colCoolCI, border=NA)
    # biotic warm
    resBWF <- vFitbyDF$res$visregRes[vFitbyDF$res$RL=='warm' & vFitbyDF$res$Driver.type=='biotic']; head(resBWF)
    points(jitter(x=rep(2-adj, length=length(resBWF)), factor=9),
                  y=resBWF, cex=cxp, pch=16, col=colWarmp)
    rect(xleft=2-adj-r, xright=2-adj+r, ytop=byDrivFlsm$asymp.UCL[4], ybottom=byDrivFlsm$asymp.LCL[4], col=colWarmCI, border=NA)
    # abiotic cool
    resACF <- vFitbyDF$res$visregRes[vFitbyDF$res$RL=='cool' & vFitbyDF$res$Driver.type=='abiotic']; head(resACF)
    points(jitter(x=rep(3+adj, length=length(resACF)), factor=5),
                  y=resACF, cex=cxp, pch=16, col=colCoolp)
    rect(xleft=3+adj-r, xright=3+adj+r, ytop=byDrivFlsm$asymp.UCL[1], ybottom =byDrivFlsm$asymp.LCL[1], col=colCoolCI, border=NA)
    #abiotic warm
    resAWF <- vFitbyDF$res$visregRes[vFitbyDF$res$RL=='warm' & vFitbyDF$res$Driver.type=='abiotic']; head(resACF)
    points(jitter(x=rep(4+adj, length=length(resAWF)), factor=4),
                  y=resAWF, cex=cxp, pch=16, col=colWarmp)
    rect(xleft=4+adj-r, xright=4+adj+r, ytop=byDrivFlsm$asymp.UCL[3], ybottom =byDrivFlsm$asymp.LCL[3], col=colWarmCI, border=NA)
    par(new=T); plot(vFitbyDF$fit$visregFit ~ vFitbyDF$fit$order, xlim=xlim, ylim=ylim, 
                       xaxt='n', yaxt='n', pch='—', cex=2.5)
    #significance
    #text(x=.1, y=ymin, labels='factor type x limit type P < 0.001', pos=4, cex=cxsig)
    #text(x=xat, y=c(rep(.98,2)), labels=c('*','*'), cex=2.4) #by stars (abandonned)
      lines(x=x[1:2], y=c(lP1,lP1), lw=lwd, col='grey48') # cool < warm for biotic
      #lines(x=x[3:4], y=c(lP1,lP1), lw=lwd, col='grey48') # cool > warm for abiotic (no longer signif March 31)
      lines(x=c(x[1],x[3]), y=c(lP2,lP2), lw=lwd, col=colC) #prediction 2 lines: @cool abio > bio
      lines(x=c(x[2],x[4]), y=c(lP2,lP2)+.04, lw=lwd, col=colW) #prediction 2 lines: @warm abio < bio
    #labels  
      mtext(at=1.3, side=3, line=lPlet, text='B) Field experiments', cex=1.2) #left justified
      #mtext(at=2.5, side=3, line=lPlet, text='B) Field experiments', cex=1.2) #centered
      mtext(at=xat, side=3, line=, text=c('Biotic','Abiotic'), cex=1.1)
      #mtext(at=xat, side=3, line=, text=c('Biotic factors','Abiotic factors'), cex=1.1)
      #mtext(at=2.3, side=1, line=2, text='Type of range limit', cex=cxaxl)
      mtext(at=0, side=1, line=2, text='Type of range limit', cex=cxaxl) #one central label



#Fig.3 indiv drivers ---------------------------------------------

#__data Fig 3---------
summary(allRL)
levels(allRL$Driver.examined)
test <- allRL[allRL$Driver=='A.Other',]; test <- droplevels(test)
summary(test[,c('Driver','Driver.descrip')]) #bouyancy and upwellings all from one study
table(allRL$Driver, allRL$Driver.examined)

#redo model with Driver -> shorter code and easier to interpret with renamed variable
coolvwarm.driver <- glmer(Driver.supported01 ~ Driver*RL + (1|Sp.name) + (1|doi), #no warnings!
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=allRL, subset=Driver.examined!='Other') 

  driv <- coolvwarm.driver
  Drivlsm <- data.frame(summary(lsmeans(driv, ~ RL | Driver, type='response'))); 
  Drivlsm 
  # visreg(driv, xvar='RL', by="Driver",  partial=T, scale='response')
  # visreg(driv, xvar='Driver', by="RL",  partial=T, scale='response')
  adj4=0.18
  vDriv <- visreg(driv, xvar='RL', by="Driver",  partial=T, scale='response')
  head(vDriv$fit)
  BAgap=.5
  lsmeans(driv, ~ Driver | RL, type='response')  
  #if want to order things by range limit and then factors in order of importance:
  RLgap=.5
  vDriv$fit$orderbyRL <- ifelse(vDriv$fit$Driver=='A.Temp' & vDriv$fit$RL=='cool', 1,
                         ifelse(vDriv$fit$Driver=='A.Climate' & vDriv$fit$RL=='cool', 1.5,
                         ifelse(vDriv$fit$Driver=='B.Habitat' & vDriv$fit$RL=='cool', 2,        
                         ifelse(vDriv$fit$Driver=='A.Soil' & vDriv$fit$RL=='cool', 2.5,
                         ifelse(vDriv$fit$Driver=='A.Precip' & vDriv$fit$RL=='cool', 3,
                         ifelse(vDriv$fit$Driver=='B.HostFood' & vDriv$fit$RL=='cool', 3.5,       
                         ifelse(vDriv$fit$Driver=='B.Comp' & vDriv$fit$RL=='cool', 4,       
                         ifelse(vDriv$fit$Driver=='B.PredHerb' & vDriv$fit$RL=='cool', 4.5,       
                         ifelse(vDriv$fit$Driver=='B.ParasitDis' & vDriv$fit$RL=='cool', 5,
                         #warm
                         ifelse(vDriv$fit$Driver=='B.PredHerb' & vDriv$fit$RL=='warm', 5.5+RLgap,
                         ifelse(vDriv$fit$Driver=='B.Habitat' & vDriv$fit$RL=='warm', 6+RLgap,
                         ifelse(vDriv$fit$Driver=='A.Temp' & vDriv$fit$RL=='warm', 6.5+RLgap,
                         ifelse(vDriv$fit$Driver=='A.Precip' & vDriv$fit$RL=='warm', 7+RLgap,       
                         ifelse(vDriv$fit$Driver=='B.Comp' & vDriv$fit$RL=='warm', 7.5+RLgap,
                         ifelse(vDriv$fit$Driver=='B.HostFood' & vDriv$fit$RL=='warm', 8+RLgap,        
                         ifelse(vDriv$fit$Driver=='A.Climate' & vDriv$fit$RL=='warm', 8.5+RLgap,        
                         ifelse(vDriv$fit$Driver=='B.ParasitDis' & vDriv$fit$RL=='warm', 9+RLgap,        
                         ifelse(vDriv$fit$Driver=='A.Soil' & vDriv$fit$RL=='warm', 9.5+RLgap, 10))))))))))))))))))
  
#subset residual and CI data by factor first so plot code less crazy
#Abiotic factors
    resTemp <- vDriv$res[vDriv$res$Driver=='A.Temp',]; head(resTemp)
    lsmTemp <- Drivlsm[Drivlsm$Driver=='A.Temp',]; lsmTemp
    resClim <- vDriv$res[vDriv$res$Driver=='A.Climate',]; head(resClim)
    lsmClim <- Drivlsm[Drivlsm$Driver=='A.Climate',]; lsmClim
    resPrecip <- vDriv$res[vDriv$res$Driver=='A.Precip',]; head(resPrecip)
    lsmPrecip <- Drivlsm[Drivlsm$Driver=='A.Precip',]; lsmPrecip
    resSoil <- vDriv$res[vDriv$res$Driver=='A.Soil',]; head(resSoil)
    lsmSoil <- Drivlsm[Drivlsm$Driver=='A.Soil',]; lsmSoil
#Biotic factors
    resHab <- vDriv$res[vDriv$res$Driver=='B.Habitat',]; head(resHab) # biogenic habitat
    lsmHab <- Drivlsm[Drivlsm$Driver=='B.Habitat',]; #lsmHab
    resHost <- vDriv$res[vDriv$res$Driver=='B.HostFood',]; head(resHost)
    lsmHost <- Drivlsm[Drivlsm$Driver=='B.HostFood',]; lsmHost
    resDis <- vDriv$res[vDriv$res$Driver=='B.ParasitDis',]; head(resDis)
    lsmDis <- Drivlsm[Drivlsm$Driver=='B.ParasitDis',]; lsmDis
    resComp <- vDriv$res[vDriv$res$Driver=='B.Comp',]; head(resComp)
    lsmComp <- Drivlsm[Drivlsm$Driver=='B.Comp',]; lsmComp
    resPred <- vDriv$res[vDriv$res$Driver=='B.PredHerb',]; head(resComp)
    lsmPred <- Drivlsm[Drivlsm$Driver=='B.PredHerb',]; lsmComp

#__plot main results ---------
quartz(height=4.5, width=8) #if want labels rotated
    cxp=.5; xtck=c(seq(0,1))#if want to plot points have to extend x axis
    abw=1.5; xl=1.8; yl=2.6; xtl=.6
    cxaxt=.85; cxsig=.9 #size axis tick labels
    xct=1; r=.18 #defines width of rectangle
    xlim=c(.8,9.8); x=c(1,2,3,4,5)
    ymin=0; ylim=c(ymin,1) #for transformed scale
#original blue / red colours    
    aCI=.6; ap=.4 #lower values of alpha = more transparent
    colCool='dodgerblue'; colWarm='red3' #if using red blue colour scheme from rest of paper
    colCoolCI <- adjustcolor(colCool, alpha.f=aCI); colWarmCI <- adjustcolor(colWarm, alpha.f=aCI) 
    colCoolp <- adjustcolor(colCool, alpha.f=ap);   colWarmp <- adjustcolor(colWarm, alpha.f=ap) 
#biotic / abiotic colours    
    colA='grey'; colB='chartreuse3'; 
    colACI <- adjustcolor(colA, alpha.f=aCI); colBCI <- adjustcolor(colB, alpha.f=aCI) 
    colAp <- adjustcolor(colA, alpha.f=ap);   colBp <- adjustcolor(colB, alpha.f=ap) 

par(oma=c(1.9,3.5,1.4,0), mar=c(3.1,1,0,.5), bty='l') 
plot(vDriv$fit$visregFit ~ vDriv$fit$orderbyRL, xlim=xlim, ylim=ylim, las=1, pch='—', cex=1.3, xaxt='n', yaxt='n') 
    axis (side=1, at=c(seq(1,5,by=.5),seq(5.5,9.5,by=.5)+RLgap), labels=NA, tck=-0.02)
    #cool range limits
    at=seq(1,5,by=.5); cxxl=.7
    text(x=at, y=-.12, srt=35, adj=1, cex=cxxl, xpd=T,#rotated driver labels
       labels=c('Temp','Climate','Bio habitat','Soil','Moisture','Food/Hosts','Competition','Pred/Herbiv','Pathogens'))
    mtext (side=1, at=at, line=.1, cex=cxxl, text=c('A','A','B','A','A','B','B','B','B'), xpd=T)
    #warm range limits
    at=seq(5.5,9.5,by=.5)+RLgap; 
    text(x=at, y=-.12, srt=35, adj=1, cex=cxxl, xpd=T,
        labels=c('Pred/Herbiv','Bio habitat','Temp','Moisture','Competition','Food/Hosts','Climate','Pathogens','Soil'))
    mtext (side=1, at=at, line=.1, cex=cxxl, text=c('B','B','A','A','B','B','A','B','A'), xpd=T)
    abline(h=.5)
    axis(side=2, at=seq(0,1,by=.2), las=1)
    mtext(side=2, line=3, text='Probability factor contributes to range limit', cex=cxaxl)
#Cool limits first
    x=1; points(jitter(x=rep(x, length=length(resTemp$visregRes[resTemp$RL=='cool'])), factor=7),
                       y=resTemp$visregRes[resTemp$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmTemp$asymp.UCL[1], ybottom=lsmTemp$asymp.LCL[1], col=colCoolCI, border=NA)
    x=1.5; points(jitter(x=rep(x, length=length(resClim$visregRes[resClim$RL=='cool'])), factor=3),
                  y=resClim$visregRes[resClim$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmClim$asymp.UCL[1], ybottom=lsmClim$asymp.LCL[1], col=colCoolCI, border=NA)
    x=2; points(jitter(x=rep(x,length=length(resHab$visregRes[resHab$RL=='cool'])),factor=3),
                  y=resHab$visregRes[resHab$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmHab$asymp.UCL[1], ybottom=lsmHab$asymp.LCL[1], col=colCoolCI, border=NA)
    x=2.5; points(jitter(x=rep(x, length=length(resSoil$visregRes[resSoil$RL=='cool'])), factor=2),
                  y=resSoil$visregRes[resSoil$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmSoil$asymp.UCL[1], ybottom=lsmSoil$asymp.LCL[1], col=colCoolCI, border=NA)
     x=3; points(jitter(x=rep(x, length=length(resPrecip$visregRes[resPrecip$RL=='cool'])), factor=2),
                  y=resPrecip$visregRes[resPrecip$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmPrecip$asymp.UCL[1], ybottom=lsmPrecip$asymp.LCL[1], col=colCoolCI, border=NA)
    x=3.5; points(jitter(x=rep(x, length=length(resHost$visregRes[resHost$RL=='cool'])),factor=2),
                        y=resHost$visregRes[resHost$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmHost$asymp.UCL[1], ybottom=lsmHost$asymp.LCL[1], col=colCoolCI, border=NA)
    x=4; points(jitter(x=rep(x,length=length(resComp$visregRes[resComp$RL=='cool'])),factor=2),
                  y=resComp$visregRes[resComp$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmComp$asymp.UCL[1], ybottom=lsmComp$asymp.LCL[1], col=colCoolCI, border=NA)
    x=4.5; points(jitter(x=rep(x,length=length(resPred$visregRes[resPred$RL=='cool'])),factor=1.5),
                  y=resPred$visregRes[resPred$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmPred$asymp.UCL[1], ybottom=lsmPred$asymp.LCL[1], col=colCoolCI, border=NA)
    x=5; points(jitter(x=rep(x,length=length(resDis$visregRes[resDis$RL=='cool'])),factor=1.5),
                  y=resDis$visregRes[resDis$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmDis$asymp.UCL[1], ybottom=lsmDis$asymp.LCL[1], col=colCoolCI, border=NA)
#Warm range limits: 'Biohab Pred/H Temp Comp Moisture Food/Hosts Climate Soil Pathogens
    x=5.5+RLgap; points(jitter(x=rep(x,length=length(resPred$visregRes[resPred$RL=='warm'])),factor=.9), #Pred/Herbiv
                  y=resPred$visregRes[resPred$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmPred$asymp.UCL[2], ybottom=lsmPred$asymp.LCL[2], col=colWarmCI, border=NA)
    x=6+RLgap; points(jitter(x=rep(x,length=length(resHab$visregRes[resHab$RL=='warm'])),factor=1.3),
                  y=resHab$visregRes[resHab$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmHab$asymp.UCL[2], ybottom=lsmHab$asymp.LCL[2], col=colWarmCI, border=NA)
    x=6.5+RLgap; points(jitter(x=rep(x,length=length(resTemp$visregRes[resTemp$RL=='warm'])),factor=1),  #Temp
                  y=resTemp$visregRes[resTemp$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmTemp$asymp.UCL[2], ybottom=lsmTemp$asymp.LCL[2], col=colWarmCI, border=NA)
    x=7+RLgap; points(jitter(x=rep(x,length=length(resPrecip$visregRes[resPrecip$RL=='warm'])),factor=.9),
                  y=resPrecip$visregRes[resPrecip$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmPrecip$asymp.UCL[2], ybottom=lsmPrecip$asymp.LCL[2], col=colWarmCI, border=NA)
    #Competition
    x=7.5+RLgap; points(jitter(x=rep(x,length=length(resComp$visregRes[resComp$RL=='warm'])),factor=1),
                  y=resComp$visregRes[resComp$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmComp$asymp.UCL[2], ybottom=lsmComp$asymp.LCL[2], col=colWarmCI, border=NA)
    #Food/hosts
    x=8+RLgap; points(jitter(x=rep(x,length=length(resHost$visregRes[resHost$RL=='warm'])),factor=1),
                  y=resHost$visregRes[resHost$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmHost$asymp.UCL[2], ybottom=lsmHost$asymp.LCL[2], col=colWarmCI, border=NA)
    #Climate
    x=8.5+RLgap; points(jitter(x=rep(x,length=length(resClim$visregRes[resClim$RL=='warm'])),factor=.8),
                  y=resClim$visregRes[resClim$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmClim$asymp.UCL[2], ybottom=lsmClim$asymp.LCL[2], col=colWarmCI, border=NA)
    x=9.+RLgap; points(jitter(x=rep(x,length=length(resDis$visregRes[resDis$RL=='warm'])),factor=.7),
                  y=resDis$visregRes[resDis$RL=='warm'],cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmDis$asymp.UCL[2], ybottom=lsmDis$asymp.LCL[2], col=colWarmCI, border=NA)
    #soil
    x=9.5+RLgap; points(jitter(x=rep(x,length=length(resSoil$visregRes[resSoil$RL=='warm'])),factor=.8),
                  y=resSoil$visregRes[resSoil$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmSoil$asymp.UCL[2], ybottom=lsmSoil$asymp.LCL[2], col=colWarmCI, border=NA)
        
    mtext(at=c(1.8,6.8), side=3, line=.1, text=c('Cool range limits','Warm range limits'), cex=1.1)
    #text(x=.4, y=ymin, labels='factor x limit type P < 0.001', pos=4, cex=cxsig)
    #at coolRL: Temp a  Climate ab  BioHab ab  Soil ab  Precip b  Food/Host abc  Comp b  PredH c  Pathogen c
    text(x=seq(1,5,by=.5), y=rep(.99,5), labels=c('a','ab','ab','ab','b','abc','b','c','c'), cex=.7, col='dodgerblue3') #cool
    #text(x=seq(5.5,9.55,by=.5)+RLgap, y=rep(.1,5), labels=c('a','a','a','a','a','a','a','a','a'), cex=.7, col=colWarm)
    mtext(at=5, side=1, line=3.5, text='Potentially range-limiting factor assessed', cex=1.1)
    n <- allRL[allRL$Driver.examined!='Other',]; n <- droplevels(n); table(n$Driver, n$RL)
    text(x=seq(1,5,by=.5), y=rep(0,5), labels=c('379','67','34','56','134','41','106','39','17'), cex=.7, col='dodgerblue3') #cool
    text(x=seq(6,10,by=.5), y=rep(0,5), labels=c('29','36','310','152','99','26','47','14','49'), cex=.7, col=colWarm) 

        
#__plot again for field experiments -------------------    
coolvwarm.driverF <- glmer(Driver.supported01 ~ Driver*RL + (1|Sp.name) + (1|doi), #singular
                          family=binomial,
                          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)),
                          data=Fdat, subset=Driver.examined!='Other' & Driver.examined!='Host/Food availability') 
  driv <- coolvwarm.driverF
  Drivlsm <- data.frame(summary(lsmeans(driv, ~ RL | Driver, type='response'))); 
  Drivlsm 
  adj4=0.18
  vDriv <- visreg(driv, xvar='RL', by="Driver",  partial=T, scale='response')
  BAgap=.5;  RLgap=.5
  lsmeans(driv, pairwise ~ Driver | RL, type='response')  
  #if want to order things by range limit and then factors in order of importance:
    #F at cool:   Temp a Clim ab > Precip abc Soil abc Comp > bc Pred/Herb c Paras c
  vDriv$fit$orderbyRL <- ifelse(vDriv$fit$Driver=='A.Temp' & vDriv$fit$RL=='cool', 1,
                         ifelse(vDriv$fit$Driver=='A.Climate' & vDriv$fit$RL=='cool', 1.5,
                         ifelse(vDriv$fit$Driver=='A.Soil' & vDriv$fit$RL=='cool', 2,        
                         ifelse(vDriv$fit$Driver=='A.Precip' & vDriv$fit$RL=='cool', 2.5,
                         ifelse(vDriv$fit$Driver=='B.Comp' & vDriv$fit$RL=='cool', 3,
                         ifelse(vDriv$fit$Driver=='B.PredHerb' & vDriv$fit$RL=='cool', 3.5,       
                         ifelse(vDriv$fit$Driver=='B.ParasitDis' & vDriv$fit$RL=='cool', 4,       
    #F at warmRL: Comp a  Pred a Precip a  Paras ab Temp ab   Soil b  Clim b
                         ifelse(vDriv$fit$Driver=='B.PredHerb' & vDriv$fit$RL=='warm', 5.5+RLgap,
                         ifelse(vDriv$fit$Driver=='B.Comp' & vDriv$fit$RL=='warm', 6+RLgap,
                         ifelse(vDriv$fit$Driver=='A.Precip' & vDriv$fit$RL=='warm', 6.5+RLgap,
                         ifelse(vDriv$fit$Driver=='B.ParasitDis' & vDriv$fit$RL=='warm', 7+RLgap,
                         ifelse(vDriv$fit$Driver=='A.Temp' & vDriv$fit$RL=='warm', 7.5+RLgap,       
                         ifelse(vDriv$fit$Driver=='A.Climate' & vDriv$fit$RL=='warm', 8+RLgap,        
                         ifelse(vDriv$fit$Driver=='A.Soil' & vDriv$fit$RL=='warm', 8.5+RLgap, NA))))))))))))))
#subset residual and CI data by factor first so plot code less crazy
#Abiotic factors
    resTemp <- vDriv$res[vDriv$res$Driver=='A.Temp',]; head(resTemp)
    lsmTemp <- Drivlsm[Drivlsm$Driver=='A.Temp',]; lsmTemp
    resClim <- vDriv$res[vDriv$res$Driver=='A.Climate',]; head(resClim)
    lsmClim <- Drivlsm[Drivlsm$Driver=='A.Climate',]; lsmClim
    resPrecip <- vDriv$res[vDriv$res$Driver=='A.Precip',]; head(resPrecip)
    lsmPrecip <- Drivlsm[Drivlsm$Driver=='A.Precip',]; lsmPrecip
    resSoil <- vDriv$res[vDriv$res$Driver=='A.Soil',]; head(resSoil)
    lsmSoil <- Drivlsm[Drivlsm$Driver=='A.Soil',]; lsmSoil
#Biotic factors
    resDis <- vDriv$res[vDriv$res$Driver=='B.ParasitDis',]; head(resDis)
    lsmDis <- Drivlsm[Drivlsm$Driver=='B.ParasitDis',]; lsmDis
    resComp <- vDriv$res[vDriv$res$Driver=='B.Comp',]; head(resComp)
    lsmComp <- Drivlsm[Drivlsm$Driver=='B.Comp',]; lsmComp
    resPred <- vDriv$res[vDriv$res$Driver=='B.PredHerb',]; head(resComp)
    lsmPred <- Drivlsm[Drivlsm$Driver=='B.PredHerb',]; lsmComp
        
quartz(height=4.5, width=8) #if want labels rotated
    cxp=.5; xtck=c(seq(0,1))
    abw=1.5; xl=1.8; yl=2.6; xtl=.6
    cxaxt=.85; cxsig=.9; cxaxl=1.1 #size axis tick labels
    xct=1; r=.18 #defines width of rectangle
    xlim=c(.8,9.8); x=c(1,2,3,4,5)
    ymin=0; ylim=c(ymin,1) #for transformed scale
    #colours are set in code for Fig 3
 
par(oma=c(1.9,3.5,1.4,0), mar=c(3.1,1,0,.5), bty='l') 
plot(vDriv$fit$visregFit ~ vDriv$fit$orderbyRL, xlim=xlim, ylim=ylim, las=1, pch='—', cex=1.3, xaxt='n', yaxt='n')
    xC=seq(1,4,by=.5); xW=seq(5.5,8.5,by=.5)+RLgap #define ticks for cool & warm RL
    axis(side=1, at=c(xC,xW), labels=NA, tck=-0.02)
    #cool range limits
    cxxl=.7
    text(x=xC, y=-.12, srt=35, adj=1, cex=cxxl, xpd=T,#rotated driver labels
       labels=c('Temp','Climate','Soil','Moisture','Competition','Pred/Herbiv','Pathogens'))
    mtext (side=1, at=xC, line=.1, cex=cxxl, text=c('A','A','A','A','B','B','B'), xpd=T)
    #warm range limits
    text(x=xW, y=-.12, srt=35, adj=1, cex=cxxl, xpd=T,
        labels=c('Pred/Herbiv','Competition','Moisture','Pathogens','Temp','Climate','Soil'))
    mtext (side=1, at=xW, line=.1, cex=cxxl, text=c('B','B','A','B','A','A','A'), xpd=T)
    abline(h=.5)
    axis(side=2, at=seq(0,1,by=.2), las=1)
    mtext(side=2, line=3, text='Probability factor contributes to range limit', cex=cxaxl)
#Cool limits first
    x=1; points(jitter(x=rep(x, length=length(resTemp$visregRes[resTemp$RL=='cool'])), factor=7),
                       y=resTemp$visregRes[resTemp$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmTemp$asymp.UCL[1], ybottom=lsmTemp$asymp.LCL[1], col=colCoolCI, border=NA)
    x=1.5; points(jitter(x=rep(x, length=length(resClim$visregRes[resClim$RL=='cool'])), factor=3),
                  y=resClim$visregRes[resClim$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmClim$asymp.UCL[1], ybottom=lsmClim$asymp.LCL[1], col=colCoolCI, border=NA)
    x=2.; points(jitter(x=rep(x, length=length(resSoil$visregRes[resSoil$RL=='cool'])), factor=1.5),
                  y=resSoil$visregRes[resSoil$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmSoil$asymp.UCL[1], ybottom=lsmSoil$asymp.LCL[1], col=colCoolCI, border=NA)
    x=2.5; points(jitter(x=rep(x, length=length(resPrecip$visregRes[resPrecip$RL=='cool'])), factor=2),
                  y=resPrecip$visregRes[resPrecip$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmPrecip$asymp.UCL[1], ybottom=lsmPrecip$asymp.LCL[1], col=colCoolCI, border=NA)
    x=3; points(jitter(x=rep(x,length=length(resComp$visregRes[resComp$RL=='cool'])),factor=2),
                  y=resComp$visregRes[resComp$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmComp$asymp.UCL[1], ybottom=lsmComp$asymp.LCL[1], col=colCoolCI, border=NA)
    x=3.5; points(jitter(x=rep(x,length=length(resPred$visregRes[resPred$RL=='cool'])),factor=1.5),
                  y=resPred$visregRes[resPred$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmPred$asymp.UCL[1], ybottom=lsmPred$asymp.LCL[1], col=colCoolCI, border=NA)
    x=4; points(jitter(x=rep(x,length=length(resDis$visregRes[resDis$RL=='cool'])),factor=1.5),
                  y=resDis$visregRes[resDis$RL=='cool'], cex=cxp, pch=16, col=colCoolp)
        rect(xleft=x-r, xright=x+r, ytop=lsmDis$asymp.UCL[1], ybottom=lsmDis$asymp.LCL[1], col=colCoolCI, border=NA)
#Warm range limits
    #Pred/Herbiv
    x=5.5+RLgap; points(jitter(x=rep(x,length=length(resPred$visregRes[resPred$RL=='warm'])),factor=1),
                  y=resPred$visregRes[resPred$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmPred$asymp.UCL[2], ybottom=lsmPred$asymp.LCL[2], col=colWarmCI, border=NA)
    #Competition
    x=6+RLgap; points(jitter(x=rep(x,length=length(resComp$visregRes[resComp$RL=='warm'])),factor=1),
                  y=resComp$visregRes[resComp$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmComp$asymp.UCL[2], ybottom=lsmComp$asymp.LCL[2], col=colWarmCI, border=NA)
    x=6.5+RLgap; points(jitter(x=rep(x,length=length(resPrecip$visregRes[resPrecip$RL=='warm'])),factor=1),
                  y=resPrecip$visregRes[resPrecip$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmPrecip$asymp.UCL[2], ybottom=lsmPrecip$asymp.LCL[2], col=colWarmCI, border=NA)
    x=7+RLgap; points(jitter(x=rep(x,length=length(resDis$visregRes[resDis$RL=='warm'])),factor=.7),
                  y=resDis$visregRes[resDis$RL=='warm'],cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmDis$asymp.UCL[2], ybottom=lsmDis$asymp.LCL[2], col=colWarmCI, border=NA)
    #Temp
    x=7.5+RLgap; points(jitter(x=rep(x,length=length(resTemp$visregRes[resTemp$RL=='warm'])),factor=1),
                  y=resTemp$visregRes[resTemp$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmTemp$asymp.UCL[2], ybottom=lsmTemp$asymp.LCL[2], col=colWarmCI, border=NA)
    #Climate
    x=8.+RLgap; points(jitter(x=rep(x,length=length(resClim$visregRes[resClim$RL=='warm'])),factor=.8),
                  y=resClim$visregRes[resClim$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmClim$asymp.UCL[2], ybottom=lsmClim$asymp.LCL[2], col=colWarmCI, border=NA)
    #soil
    x=8.5+RLgap; points(jitter(x=rep(x,length=length(resSoil$visregRes[resSoil$RL=='warm'])),factor=.8),
                  y=resSoil$visregRes[resSoil$RL=='warm'], cex=cxp, pch=16, col=colWarmp)
        rect(xleft=x-r, xright=x+r, ytop=lsmSoil$asymp.UCL[2], ybottom=lsmSoil$asymp.LCL[2], col=colWarmCI, border=NA)
    mtext(at=c(1.8,6.8), side=3, line=.1, text=c('Cool range limits','Warm range limits'), cex=1.1)
    #significance letters
    #F at cool:   Temp a Clim a Soil ab Precip ab Comp ab  Pred/Herb b Paras b
    #F at warmRL: Pred a Comp a  Precip a  Paras ab Temp ab   Soil b  Clim b
    text(x=xC, y=rep(1,5), labels=c('a','a','ab','ab','ab','b','b'), cex=.7, col='dodgerblue3')     
    text(x=xW, y=rep(1,5), labels=c('a','a','a','ab','ab','b','b'), cex=.7, col=colWarm)     
    mtext(at=5, side=1, line=3.5, text='Potentially range-limiting factor assessed in field experiments', cex=1.1)
    #sample sizes
    n <- Fdat[Fdat$Driver.examined!='Other' & Fdat$Driver.examined!='Host/Food availability',]; 
    n <- droplevels(n); nC <- table(n$Driver[n$RL=='cool']); nC #gave up on automating further
    text(x=xC, y=rep(-.015,5), labels=c('117','20','19','33','77','37','11'), cex=.7, col='dodgerblue3')     
    nW <- table(n$Driver[n$RL=='warm']); nW
    text(x=xW, y=rep(-.015,5), labels=c('23','53','24','4','52','8','7'), cex=.7, col=colWarm)     
        

            
#Fig.4 absLatitude interaction -------------------------      
#need way to get confidence intervals around a line
#in SciAdv paper (2018) was able to get upper and lower CI from a GLMM binomial model, using visreg.  But if try that now (even on same model with same code) the columns for upper and lower CI are NA.  Think visreg stopped doing this

coolvwarm.lat #winning binomidal GLMM
#get data without CIs   
vLatBbt <- visreg(coolvwarm.lat, xvar='absLat', by='RL', cond=list(Driver.type='biotic'), type='conditional', scale='response', partial=T)    
vLatAbt <- visreg(coolvwarm.lat, xvar='absLat', by='RL', cond=list(Driver.type='abiotic'), type='conditional', scale='response', partial=T)    

#__Ben Bolker wrote a function for getting confidence intervals 'easyPredCI' -------
#https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html
##try workign through 'Gopher Totoise' example, bit over half way down the page. ends with plot with confidence intervals
#Grouse ticks examples right after looks promising too (poisson GLMM) but dont plot results except in first ggplot example
#BB: Computing confidence intervals on the predicted values is relatively easy if we’re willing to completely ignore the random effects, and the uncertainty of the random effects. Here is a generic function that extracts the relevant bits from a fitted model and returns the confidence intervals for predictions:
easyPredCI <- function(model, newdata=NULL, alpha=0.05) {
    ## baseline prediction, on the linear predictor (logit) scale:
    pred0 <- predict(model, re.form=NA, newdata=newdata)
    ## fixed-effects model matrix for new data
    X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
    beta <- fixef(model) ## fixed-effects coefficients
    V <- vcov(model)     ## variance-covariance matrix of beta
    pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
    ## inverse-link function
    linkinv <- family(model)$linkinv
    ## construct 95% Normal CIs on the link scale and
    ##  transform back to the response (probability) scale:
    crit <- -qnorm(alpha/2)
    linkinv(cbind(conf.low=pred0-crit*pred.se,
                  conf.high=pred0+crit*pred.se))
}

#make reference grid of data that want to predict across  
summary(RLloc[,c('Driver.type','RL','absLat')]) #absLat only goes to 76
pframe <- expand.grid(absLat=0:76, Driver.type=c('abiotic','biotic'), RL=c('cool','warm')); summary(pframe) #creates grid for estimating/plotting
#use predict to get predict lines
predline <- predict(coolvwarm.lat, #GLMM
                  newdata=pframe, #predict over data fram defined in line above
                  re.form=NA, #model is predicts across any level of random effect (stats.stackexchange.com/questions/452257/glmm-why-are-my-confidence-intervals-so-wide-when-my-p-values-are-so-small)
                  type="response")
#use function to get confidence intervals  
predCI <- easyPredCI(coolvwarm.lat, newdata=pframe) 
#put all together in one data frame
pred <- as.data.frame(cbind(pframe,predline,predCI)); head(pred)

#try simulation code - see if CI make more sense for significant trend
#result is biggish bootMer object (takes >1h to run)
set.seed(101)
simCI <- bootMer(coolvwarm.lat, #(expect to take 1.5 hours for 400!, >4.5h for 1000)
              FUN=function(x)
              predict(x, re.form=NA, newdata=pframe,
              type="response"),
#              nsim=400)
              nsim=1000)
head(simCI)
#extract the runs corresponding to 95% confidence intervals
bootCI <- t(apply(simCI$t, 2, quantile, c(0.025,0.975), na.rm=TRUE))
colnames(bootCI) <- c('conf.low.boot','conf.high.boot')
#add to data frame. since used same grid to predict just just stick these columns on
pred <- as.data.frame(cbind(pframe,predline,predCI,bootCI)); head(pred)

summary(pred); dim(pred)
write.csv(pred,'~/Q3 confidence intervals from 1000 simulations 21 06 19.csv')


#_All data biotic vs abioitc panels--------------------
quartz(height=4, width=8)
layout(matrix(c(1,2), nrow=1, ncol=2, byrow = T)); 
    cxp=.6; cxl=2; 
    cxaxt=.85; #size axis tick labels
    xlim=c(0,80); ymin=-.01; ylim=c(ymin,1.01)
    cxaxl=1.1
    colCool='dodgerblue'; colWarm='red3'; 
    colCoolCI <- adjustcolor(colCool, alpha.f=.3); 
    colWarmCI <- adjustcolor(colWarm, alpha.f=.3) 
    colgreyCI <- adjustcolor('grey', alpha.f=.3); 
    
par(oma=c(2,3.5,1,0), mar=c(1.5,1,1,.5), bty='l') 
# #what it should look like
# visreg(coolvwarm.lat, xvar='absLat', by='RL', cond=list(Driver.type='biotic'), type='conditional', scale='response', 
#        ylim=ylim, partial=T, overlay=T, main='Biotic'); abline(h=.5)
# visreg(coolvwarm.lat, xvar='absLat', by='RL', cond=list(Driver.type='abiotic'), type='conditional', scale='response',
#        ylim=ylim, yaxt='n', partial=T, overlay=T, main='Abio'); abline(h=.5)

#A) Biotic drivers ---------------- ---
fit = vLatBbt$fit; head(fit) #vLatBbt -> visregFit IS backtransformed
fitCool <- fit[fit$RL=='cool',]; fitWarm <- fit[fit$RL=='warm',]
res = vLatBbt$res; head(res) #need the backtransformed ones
resCool <- res[res$RL=='cool',]; resWarm <- res[res$RL=='warm',]
#plot
head(fitCool)
plot(fitCool$visregFit ~ fitCool$absLat, type='l', xlim=xlim, ylim=ylim, xlab='', col=colCool, las=1, lwd=cxl)
  abline(h=.5)
#cool limits confidence intervals from Predict
  head(pred)
  CI <- pred[pred$Driver.type=='biotic' & pred$RL=='cool',]; head(CI)
  # #CI ignoring random effects (smooth but bit wide - signif trend looks NS)
  # polygon(y=c(CI$conf.low, rev(CI$conf.high)), x=c(CI$absLat, rev(CI$absLat)), col=colCoolCI, border=NA)
  # #CI with random effects incl (bootstrapped from 1000 simulations)
  polygon(y=c(CI$conf.low.boot, rev(CI$conf.high.boot)), x=c(CI$absLat, rev(CI$absLat)), col=colCoolCI, border=NA)
  #lines(CI$predline ~ CI$absLat, col=colCool, lwd=cxl, lty=1) #check that line from predict = line from visreg
  points(y=resCool$visregRes, x=resCool$absLat, cex=cxp, pch=16, col=colCoolCI)
#warm limits
head(fitWarm)
lines(fitWarm$visregFit[fitWarm$absLat<64] ~ fitWarm$absLat[fitWarm$absLat<64], col=colWarm, lwd=cxl)
lines(fitWarm$visregFit[fitWarm$absLat>64] ~ fitWarm$absLat[fitWarm$absLat>64], col=colWarm, lwd=cxl, lty=3)
  CI <- pred[pred$Driver.type=='biotic' & pred$RL=='warm',]; head(CI)
  polygon(y=c(CI$conf.low.boot, rev(CI$conf.high.boot)), x=c(CI$absLat, rev(CI$absLat)), col=colWarmCI, border=NA)
  points(y=resWarm$visregRes, x=resWarm$absLat, cex=cxp, pch=16, col=colWarmCI)
  mtext(side=3, text='A) Biotic factors', at=15, cex=cxaxl, line=.5)
  mtext(side=2, line=3, text='Probability factor contributes to range limit', cex=cxaxl)
  
#B) Abiotic drivers
fit = vLatAbt$fit; head(fit) #vLatBbt -> visregFit IS backtransformed
fitCool <- fit[fit$RL=='cool',]; fitWarm <- fit[fit$RL=='warm',]
res = vLatAbt$res; head(res) #need the backtransformed ones
resCool <- res[res$RL=='cool',]; resWarm <- res[res$RL=='warm',]
plot(fitCool$visregFit ~ fitCool$absLat, type='l', xlim=xlim, ylim=ylim, xlab='', yaxt='n', col=colCool, lwd=cxl)#, lty=2)
  abline(h=.5)
  axis(side=2, labels=NA)
  #cool limits confidence intervals
  CI <- pred[pred$Driver.type=='abiotic' & pred$RL=='cool',]; head(CI)
  polygon(y=c(CI$conf.low.boot, rev(CI$conf.high.boot)), x=c(CI$absLat, rev(CI$absLat)), col=colCoolCI, border=NA)
  #lines(CI$predline ~ CI$absLat, col=colCool, lwd=cxl, lty=1) #check that line from predict = line from visreg
  points(y=resCool$visregRes, x=resCool$absLat, cex=cxp, pch=16, col=colCoolCI)
#warm limits
lines(fitWarm$visregFit[fitWarm$absLat<64] ~ fitWarm$absLat[fitWarm$absLat<64], col=colWarm, lwd=cxl)
lines(fitWarm$visregFit[fitWarm$absLat>64] ~ fitWarm$absLat[fitWarm$absLat>64], col=colWarm, lwd=cxl, lty=3)
  CI <- pred[pred$Driver.type=='abiotic' & pred$RL=='warm',]; head(CI)
  polygon(y=c(CI$conf.low.boot, rev(CI$conf.high.boot)), x=c(CI$absLat, rev(CI$absLat)), col=colWarmCI, border=NA)
  points(y=resWarm$visregRes, x=resWarm$absLat, cex=cxp, pch=16, col=colWarmCI)
  mtext(side=3, text='B) Abiotic factors', at=15, cex=cxaxl, line=.5)
  mtext(side=1, line=2.5, at=-1, text='Absolute Latitude (decimal degrees)', cex=cxaxl)

  
# Fig. S1 maps --------------------------

#(Alex did maps in separate script)
#for map legends: 
#number of cool vs warm range limits where biotic vs abiotic factors assessed
table(allRLsum$RL, allRLsum$Driver.type)
  
# Fig. S2 latitudinal effects field exp, land, marine --------------------------

#get visreg for predict lines and residuals
#Field exp only     
coolvwarm.latF 
vLatBbtF <- visreg(coolvwarm.latF, xvar='absLat', by='RL', cond=list(Driver.type='biotic'), type='conditional', scale='response', partial=T)    
vLatAbtF <- visreg(coolvwarm.latF, xvar='absLat', by='RL', cond=list(Driver.type='abiotic'), type='conditional', scale='response', partial=T)    

#LAnd studies only     
coolvwarm.latLand 
vLatBbtL <- visreg(coolvwarm.latLand, xvar='absLat', by='RL', cond=list(Driver.type='biotic'), type='conditional', scale='response', partial=T)    
vLatAbtL <- visreg(coolvwarm.latLand, xvar='absLat', by='RL', cond=list(Driver.type='abiotic'), type='conditional', scale='response', partial=T)    

#Marine studies only     
coolvwarm.latMarine
vLatBbtM <- visreg(coolvwarm.latMarine, xvar='absLat', by='RL', cond=list(Driver.type='biotic'), type='conditional', scale='response', partial=T)    
vLatAbtM <- visreg(coolvwarm.latMarine, xvar='absLat', by='RL', cond=list(Driver.type='abiotic'), type='conditional', scale='response', partial=T)  
#what its supposed to look like (remember colours are reversed)
visreg(coolvwarm.latMarine, xvar='absLat', by='RL', cond=list(Driver.type='biotic'), type='conditional', overlay=T, scale='response', ylim=c(0,1),
             ylab='Biotic support Marine') #opposite of  all study results   
visreg(coolvwarm.latMarine, xvar='absLat', by='RL', cond=list(Driver.type='abiotic'), type='conditional', overlay=T, scale='response', ylim=c(0,1),
             ylab='Abiotic support Marine') #now both decline   


#make reference grid of data that want to predict across  
summary(RLloc[,c('Driver.type','RL','absLat')]) #absLat only goes to 76
pframe <- expand.grid(absLat=0:76, Driver.type=c('abiotic','biotic'), RL=c('cool','warm')); summary(pframe) #creates grid for estimating/plotting
#Field exp
predlineF <- predict(coolvwarm.latF, #GLMM
                  newdata=pframe, re.form=NA, type="response")
predCIF <- easyPredCI(coolvwarm.latF, newdata=pframe); colnames(predCIF) <- c('conf.lowF','conf.highF') #use function to get confidence intervals 
predSI <- as.data.frame(cbind(pframe,predlineF,predCIF)); head(predSI) #put all together in one data frame
#Land
predlineL <- predict(coolvwarm.latLand, #GLMM
                  newdata=pframe, re.form=NA, type="response")
predCIL <- easyPredCI(coolvwarm.latLand, newdata=pframe); colnames(predCIL) <- c('conf.lowL','conf.highL') 
predSI <- as.data.frame(cbind(predSI,predlineL,predCIL)); head(predSI) #add to data frame
#Marine
predlineM <- predict(coolvwarm.latMarine, #GLMM
                  newdata=pframe, re.form=NA, type="response")
predCIM <- easyPredCI(coolvwarm.latMarine, newdata=pframe); colnames(predCIM) <- c('conf.lowM','conf.highM') 
predSI <- as.data.frame(cbind(predSI,predlineM,predCIM)); head(predSI) #add to data frame

#_SI plotting biotic vs abioitc panels, 1 row per data subset--------------------
quartz(height=9, width=7)
layout(matrix(c(1,2,3,4,5,6), nrow=3, ncol=2, byrow = T)); 
    cxp=.6; cxl=2; 
    cxaxt=.85; #size axis tick labels
    xlim=c(0,80); ymin=-.01; ylim=c(ymin,1.01)
    cxaxl=1
    colCool='dodgerblue'; colWarm='red3'; 
    colCoolCI <- adjustcolor(colCool, alpha.f=.3); 
    colWarmCI <- adjustcolor(colWarm, alpha.f=.3) 
    
par(oma=c(3,3.5,2,0), mar=c(2,1,2,.5), bty='l') 

#_A left) Field Biotic drivers -------------------
fit = vLatBbtF$fit; head(fit); fitCool <- fit[fit$RL=='cool',]; fitWarm <- fit[fit$RL=='warm',]
res = vLatBbtF$res; head(res); resCool <- res[res$RL=='cool',]; resWarm <- res[res$RL=='warm',]
#plot
head(fitCool)
plot(fitCool$visregFit ~ fitCool$absLat, type='l', xlim=xlim, ylim=ylim, xlab='', ylab='', col=colCool, las=1, lwd=cxl)
  abline(h=.5)
#cool limits confidence intervals from Predict
  head(predSI)
  CI.Bcool <- predSI[predSI$Driver.type=='biotic' & predSI$RL=='cool' & predSI$absLat<73,]; head(CI.Bcool)
  polygon(y=c(CI.Bcool$conf.lowF, rev(CI.Bcool$conf.highF)), x=c(CI.Bcool$absLat, rev(CI.Bcool$absLat)), col=colCoolCI, border=NA)
  points(y=resCool$visregRes, x=resCool$absLat, cex=cxp, pch=16, col=colCoolCI)
#warm limits
head(fitWarm)
lines(fitWarm$visregFit[fitWarm$absLat<60] ~ fitWarm$absLat[fitWarm$absLat<60], col=colWarm, lwd=cxl)
lines(fitWarm$visregFit[fitWarm$absLat>60] ~ fitWarm$absLat[fitWarm$absLat>60], col=colWarm, lwd=cxl, lty=3)
  CI.Bwarm <- predSI[predSI$Driver.type=='biotic' & predSI$RL=='warm' & predSI$absLat<73,]; head(CI.Bwarm)
  polygon(y=c(CI.Bwarm$conf.lowF, rev(CI.Bwarm$conf.highF)), x=c(CI.Bwarm$absLat, rev(CI.Bwarm$absLat)), col=colWarmCI, border=NA)
  #polygon(y=c(CI$conf.low.boot, rev(CI$conf.high.boot)), x=c(CI$absLat, rev(CI$absLat)), col=colWarmCI, border=NA)
  points(y=resWarm$visregRes, x=resWarm$absLat, cex=cxp, pch=16, col=colWarmCI)
  lp=.5; mtext(side=3, text='A) Field experiments', at=10, cex=cxaxl, line=lp)
  mtext(side=3, text='Biotic factors', at=40, cex=cxaxl, line=2)

#_A right) Field Abiotic drivers
fit = vLatAbtF$fit; head(fit); fitCool <- fit[fit$RL=='cool',]; fitWarm <- fit[fit$RL=='warm',]
res = vLatAbtF$res; head(res); resCool <- res[res$RL=='cool',]; resWarm <- res[res$RL=='warm',]
plot(fitCool$visregFit ~ fitCool$absLat, type='l', xlim=xlim, ylim=ylim, xlab='', yaxt='n', col=colCool, lwd=cxl)#, lty=2)
  abline(h=.5)
  axis(side=2, labels=NA)
  #cool limits confidence intervals
  CI.Acool <- predSI[predSI$Driver.type=='abiotic' & predSI$RL=='cool' & predSI$absLat<73,]; head(CI.Acool)
  polygon(y=c(CI.Acool$conf.lowF, rev(CI.Acool$conf.highF)), x=c(CI.Acool$absLat, rev(CI.Acool$absLat)), col=colCoolCI, border=NA)
  points(y=resCool$visregRes, x=resCool$absLat, cex=cxp, pch=16, col=colCoolCI)
#warm limits
lines(fitWarm$visregFit[fitWarm$absLat<60] ~ fitWarm$absLat[fitWarm$absLat<60], col=colWarm, lwd=cxl)
lines(fitWarm$visregFit[fitWarm$absLat>60] ~ fitWarm$absLat[fitWarm$absLat>60], col=colWarm, lwd=cxl, lty=3)
  CI.Awarm <- predSI[predSI$Driver.type=='abiotic' & predSI$RL=='warm' & predSI$absLat<73,]; head(CI.Awarm)
  polygon(y=c(CI.Awarm$conf.lowF, rev(CI.Awarm$conf.highF)), x=c(CI.Awarm$absLat, rev(CI.Awarm$absLat)), col=colWarmCI, border=NA)
  points(y=resWarm$visregRes, x=resWarm$absLat, cex=cxp, pch=16, col=colWarmCI)
  mtext(side=3, text='Abiotic factors', at=40, cex=cxaxl, line=2)

#_B left) Land Biotic drivers -------------------
fit = vLatBbtL$fit; head(fit); fitCool <- fit[fit$RL=='cool',]; fitWarm <- fit[fit$RL=='warm',]
res = vLatBbtL$res; head(res); resCool <- res[res$RL=='cool',]; resWarm <- res[res$RL=='warm',]
#plot
plot(fitCool$visregFit ~ fitCool$absLat, type='l', xlim=xlim, ylim=ylim, xlab='', ylab='', col=colCool, las=1, lwd=cxl)
  abline(h=.5)
#cool limits confidence intervals from Predict
  # #CI ignoring random effects (smooth but bit wide - signif trend looks NS)
  polygon(y=c(CI.Bcool$conf.lowL, rev(CI.Bcool$conf.highL)), x=c(CI.Bcool$absLat, rev(CI.Bcool$absLat)), col=colCoolCI, border=NA)
  #polygon(y=c(CI$conf.low.boot, rev(CI$conf.high.boot)), x=c(CI$absLat, rev(CI$absLat)), col=colCoolCI, border=NA)
  points(y=resCool$visregRes, x=resCool$absLat, cex=cxp, pch=16, col=colCoolCI)
#warm limits
lines(fitWarm$visregFit[fitWarm$absLat<64] ~ fitWarm$absLat[fitWarm$absLat<64], col=colWarm, lwd=cxl)
lines(fitWarm$visregFit[fitWarm$absLat>64] ~ fitWarm$absLat[fitWarm$absLat>64], col=colWarm, lwd=cxl, lty=3)
  polygon(y=c(CI.Bwarm$conf.lowL, rev(CI.Bwarm$conf.highL)), x=c(CI.Bwarm$absLat, rev(CI.Bwarm$absLat)), col=colWarmCI, border=NA)
  #polygon(y=c(CI$conf.low.boot, rev(CI$conf.high.boot)), x=c(CI$absLat, rev(CI$absLat)), col=colWarmCI, border=NA)
  points(y=resWarm$visregRes, x=resWarm$absLat, cex=cxp, pch=16, col=colWarmCI)
  lp=.5; mtext(side=3, text='B) Terrestrial environments', at=15, cex=cxaxl, line=lp)

#_B right) Land Abiotic drivers
fit = vLatAbtL$fit; head(fit); fitCool <- fit[fit$RL=='cool',]; fitWarm <- fit[fit$RL=='warm',]
res = vLatAbtL$res; head(res); resCool <- res[res$RL=='cool',]; resWarm <- res[res$RL=='warm',]
plot(fitCool$visregFit ~ fitCool$absLat, type='l', xlim=xlim, ylim=ylim, xlab='', yaxt='n', col=colCool, lwd=cxl)#, lty=2)
  abline(h=.5)
  axis(side=2, labels=NA)
  #cool limits confidence intervals
  polygon(y=c(CI.Acool$conf.lowL, rev(CI.Acool$conf.highL)), x=c(CI.Acool$absLat, rev(CI.Acool$absLat)), col=colCoolCI, border=NA)
  points(y=resCool$visregRes, x=resCool$absLat, cex=cxp, pch=16, col=colCoolCI)
#warm limits
lines(fitWarm$visregFit[fitWarm$absLat<63] ~ fitWarm$absLat[fitWarm$absLat<63], col=colWarm, lwd=cxl)
lines(fitWarm$visregFit[fitWarm$absLat>63] ~ fitWarm$absLat[fitWarm$absLat>63], col=colWarm, lwd=cxl, lty=3)
  polygon(y=c(CI.Awarm$conf.lowL, rev(CI.Awarm$conf.highL)), x=c(CI.Awarm$absLat, rev(CI.Awarm$absLat)), col=colWarmCI, border=NA)
  points(y=resWarm$visregRes, x=resWarm$absLat, cex=cxp, pch=16, col=colWarmCI)
  
#_C left) Marine Biotic drivers -------------------
fit = vLatBbtM$fit; head(fit); fitCool <- fit[fit$RL=='cool',]; fitWarm <- fit[fit$RL=='warm',]
res = vLatBbtM$res; head(res); resCool <- res[res$RL=='cool',]; resWarm <- res[res$RL=='warm',]
#plot
plot(fitCool$visregFit[fitCool$absLat>25] ~ fitCool$absLat[fitCool$absLat>25], type='l', 
     xlim=xlim, ylim=ylim, xlab='', ylab='', col=colCool, las=1, lwd=cxl)
     lines(fitCool$visregFit[fitCool$absLat<25] ~ fitWarm$absLat[fitWarm$absLat<25], col=colCool, lwd=cxl, lty=3)
  abline(h=.5)
#cool limits confidence intervals from Predict
  CI.Bcool <- predSI[predSI$Driver.type=='biotic' & predSI$RL=='cool'  & predSI$absLat>7 & predSI$absLat<78,]
  polygon(y=c(CI.Bcool$conf.lowM, rev(CI.Bcool$conf.highM)), x=c(CI.Bcool$absLat, rev(CI.Bcool$absLat)), col=colCoolCI, border=NA)
  points(y=resCool$visregRes, x=resCool$absLat, cex=cxp, pch=16, col=colCoolCI)
#warm limits 
lines(fitWarm$visregFit[fitWarm$absLat<64] ~ fitWarm$absLat[fitWarm$absLat<64], col=colWarm, lwd=cxl)
lines(fitWarm$visregFit[fitWarm$absLat>64] ~ fitWarm$absLat[fitWarm$absLat>64], col=colWarm, lwd=cxl, lty=3)
  CI.Bwarm <- predSI[predSI$Driver.type=='biotic' & predSI$RL=='warm'  & predSI$absLat>7 & predSI$absLat<78,]
  polygon(y=c(CI.Bwarm$conf.lowM, rev(CI.Bwarm$conf.highM)), x=c(CI.Bwarm$absLat, rev(CI.Bwarm$absLat)), col=colWarmCI, border=NA)
  points(y=resWarm$visregRes, x=resWarm$absLat, cex=cxp, pch=16, col=colWarmCI)
  lp=.5; mtext(side=3, text='C) Marine environments', at=15, cex=cxaxl, line=lp)

#C right) Marine Abiotic drivers
fit = vLatAbtM$fit; head(fit); fitCool <- fit[fit$RL=='cool',]; fitWarm <- fit[fit$RL=='warm',]
res = vLatAbtM$res; head(res); resCool <- res[res$RL=='cool',]; resWarm <- res[res$RL=='warm',]
plot(fitCool$visregFit[fitCool$absLat>26] ~ fitCool$absLat[fitCool$absLat>26], type='l', xlim=xlim, ylim=ylim, xlab='', yaxt='n', col=colCool, lwd=cxl)#, lty=2)
     lines(fitCool$visregFit[fitCool$absLat<26] ~ fitWarm$absLat[fitWarm$absLat<26], col=colCool, lwd=cxl, lty=3)
  abline(h=.5)
  axis(side=2, labels=NA)
  #cool limits confidence intervals
  CI.Acool <- predSI[predSI$Driver.type=='abiotic' & predSI$RL=='cool'  & predSI$absLat>7 & predSI$absLat<79,]
  polygon(y=c(CI.Acool$conf.lowM, rev(CI.Acool$conf.highM)), x=c(CI.Acool$absLat, rev(CI.Acool$absLat)), col=colCoolCI, border=NA)
  points(y=resCool$visregRes, x=resCool$absLat, cex=cxp, pch=16, col=colCoolCI)
#warm limits
lines(fitWarm$visregFit[fitWarm$absLat<64] ~ fitWarm$absLat[fitWarm$absLat<64], col=colWarm, lwd=cxl)
lines(fitWarm$visregFit[fitWarm$absLat>64] ~ fitWarm$absLat[fitWarm$absLat>64], col=colWarm, lwd=cxl, lty=3)
  CI.Awarm <- predSI[predSI$Driver.type=='abiotic' & predSI$RL=='warm'  & predSI$absLat>7 & predSI$absLat<78,]
  polygon(y=c(CI.Awarm$conf.lowM, rev(CI.Awarm$conf.highM)), x=c(CI.Awarm$absLat, rev(CI.Awarm$absLat)), col=colWarmCI, border=NA)
  points(y=resWarm$visregRes, x=resWarm$absLat, cex=cxp, pch=16, col=colWarmCI)

  mtext(side=1, line=2.5, at=-1, text='Absolute Latitude (decimal degrees)', cex=cxaxl)
