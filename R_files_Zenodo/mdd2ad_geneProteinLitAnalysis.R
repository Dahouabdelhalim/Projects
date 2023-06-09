###
###
###
###
#########################
### GENE/PROTEIN/ENZYME ANALYSIS
###
###

library(RPostgreSQL)
library(graph)
library(igraph)
library(tidyverse)
library(ggplot2)
library(stringi)
library(compiler)
library(gsubfn)
setwd('/home/smalec/Projects/semiramis')
#########
# DB setup
#########
drv <- dbDriver('PostgreSQL')
con <-
  dbConnect(
    drv,
    dbname = "causalehr",
    host = "localhost",
    port = 5432,
    user = "smalec",
    password = "mandarin"
  )

#########
# manually defined exclusions
#########
manuallyDefinedExclusions <- c('C0450442', 'C0243077', 'C0021083', 'C0013216', 'C0009429', 'C0199176', 'C0085104',
                               'C0017296', 'C0597357', 'C1519033', 'C1280903', 'C0003289', 'C0279025', 'C0279025',
                               'C0040808', 'C0009244', 'C0175723', 'C0011581', 'C0026868', 'C0683465', 'C1519595',
                               'C0042890', 'C0150321', 'C0581602', 'C1874451', 'C0012854', 'C0033972', 'C0279493',
                               'C0014935', 'C0086132', 'C0344315', 'C0016452', 'C1137094', 'C2718059', 'C0013806',
                               'C0183185', 'C0678908', 'C2827774', 'C0031765', 'C0003827', 'C0034991', 'C0581601',
                               'C0282402', 'C0150593', 'C1254359', 'C3858690', 'C1096024', 'C0020431', 'C1293131',
                               'C1964029', 'C0204727', 'C0005893', 'C0185027', 'C0439775', 'C0418981', 'C0152060',
                               'C0183210', 'C0870230', 'C0677850', 'C0376626', 'C0677850', 'C0949766', 'C0150133',
                               'C1515119', 'C0079613', 'C0455142', 'C0585941', 'C0279494', 'C0920425', 'C0033968', 
                               'C0935576', 'C1527374', 'C1269683', 'C0034618', 'C0376547', 'C0020431', 'C0185117',
                               'C0204514', 'C0439857', 'C0680038', 'C0846672', 'C0015618', 'C0182537', 'C0183115', 
                               'C0220908', 'C0279494', 'C0935576', 'C1268930', 'C0013218', 'C0034764', 'C0183210', 
                               'C0232164', 'C0281481', 'C1139730', 'C1511300', 'C1512629', 'C2986605', 'C3177188',
                               'C0030567', 'C0497327', 'C0474169', 'C0011269', 'C0599917', 'C0599917')

length(manuallyDefinedExclusions)

#########
# Covariates that are measured for subjects in sample from electronic health database
#########
getGoodCovars <- function(n) {
  fNames <- list.files(path = "data/r3dataset/data")
  i = 0
  goodData <- c()
  for (f in fNames) {
    if (sum(read.csv(file = paste("data/r3dataset/data/", f, sep = ""), header = FALSE)) >= n) 
    { #print(f) 
      i = i + 1 
      #print(i)
      goodData <- c(goodData, f) 
    }
  }
  #print(goodData)
  goodData <- gsubfn(goodData, pattern = "\\\\.csv", replacement = "")
  return(goodData)
}

#predicates <- " ('CAUSES', 'PREVENTS', 'PREDISPOSES') " 
# NOTE: 'TREATS' -- excluding because of the ambiguous causal meaning of the verb 'to treat'
predicates <- " ('CAUSES', 'PREDISPOSES', 'PREVENTS', 'TREATS') " # 'AFFECTS', 'DISRUPTS', 'AUGMENTS', 'TREATS', 'STIMULATES', 'INHIBITS') "

#predicates <- " ('CAUSES', 'PREDISPOSES', 'PREVENTS', 'TREATS', 'AFFECTS', 'DISRUPTS', 'AUGMENTS', 'TREATS', 'STIMULATES', 'INHIBITS') "
#########
# get confounders
#########
getConfounders <- function(predicates) {
  sqlTxt <- paste(" WITH ZX AS (SELECT cp.subject_cui AS covar, COUNT(*) AS scnt, cp.subject_name --, string_agg(cp.subject_semtype) as subject_semtype
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY subject_cui, subject_name),
                        ZY AS (SELECT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name --, string_agg(cp.subject_semtype) as subject_semtype
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ('C0002395', 'C0494463', 'C0276496')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.subject_cui, cp.subject_name)
                        SELECT DISTINCT ZX.covar, SUM(zx.scnt) as theta1, SUM(zy.scnt) AS theta2, lower(zx.subject_name) as covarname --, string_agg(zx.subject_semtype) 
                                       FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar 
                                       GROUP BY ZX.covar, lower(zx.subject_name)
                                       ORDER BY theta2 desc, theta1 desc;", sep = "")
  print(sqlTxt)
  confounders <- dbGetQuery(con, sqlTxt)
  write.table(file = "covariates/rawConfounders_for_depression2AD.txt", x = confounders, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(confounders)
}


confounders <- getConfounders(predicates)
print(confounders)
print(sort(confounders[,4]))

#########
#
#########
#getConfounders(predicates)

getColliders <- function(predicates) {
  colliders <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.object_cui AS covar, COUNT(*) AS scnt, cp.object_name, cp.object_semtype
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY object_cui, object_name, object_semtype),
                        ZY AS (SELECT cp.object_cui AS covar , COUNT(*) AS scnt, cp.object_name, cp.object_semtype
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ('C0002395', 'C0494463', 'C0276496')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.object_cui, cp.object_name, cp.object_semtype)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, zy.scnt AS theta2, lower(zx.object_name) as covarname, zx.object_semtype FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta2 desc, theta1 desc;", sep = ""))
  write.table(file = "covariates/rawColliders_for_depression2AD.txt", x = colliders, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(colliders)
}

#getColliders(predicates)

colliders <- getColliders(predicates)
print(colliders)
print(sort(colliders[,4]))

getMediators <- function(predicates) {
  mediators <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.object_cui AS covar, COUNT(*) AS scnt, cp.object_name as covarname
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY object_cui, object_name),
                        ZY AS (SELECT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name as covarname
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ('C0002395', 'C0494463', 'C0276496')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.subject_cui, cp.subject_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, zy.scnt AS theta2, lower(zx.covarname) as covarname FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta2 desc, theta1 desc;", sep = ""))
  write.table(file = "covariates/rawMediators_for_depression2AD.txt", x = mediators, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(mediators)
}

mediators <- getMediators(predicates)
print(mediators)
print(sort(mediators[,4]))


getMediators <- function(predicates) {
  mediators <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.object_cui AS covar, COUNT(*) AS scnt, cp.object_name as covarname
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(c.pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY object_cui, object_name),
                        ZY AS (SELECT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name as covarname
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ('C0002395', 'C0494463', 'C0276496') OR cp.object_name LIKE '%dement%' 
                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.subject_cui, cp.subject_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, zy.scnt AS theta2, lower(zx.covarname) as covarname FROM ZX AS zx, ZY AS zy WHERE zx.covar = zy.covar ORDER BY theta2 desc, theta1 desc;", sep = ""))
  write.table(file = "covariates/rawMediators_for_depression2AD.txt", x = mediators, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(mediators)
}

mediators <- getMediators(predicates)
print(mediators)
print(sort(mediators[,4]))

getPrecisionVars <- function(predicates) {
  precisionVars <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name as covarname
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.object_cui IN ('C0002395', 'C0494463', 'C0276496')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.subject_cui, cp.subject_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, lower(zx.covarname) as covarname FROM ZX AS zx ORDER BY theta1 desc;", sep = ""))
  write.table(file = "covariates/rawPrecisionVars_for_AD.txt", x = precisionVars, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(precisionVars)
}

precisionVars <- getPrecisionVars(predicates)

getOutcomeEffects <- function(predicates) {
  outcomeEffects <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.object_cui AS covar , COUNT(*) AS scnt, cp.object_name as covarname
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ('C0002395', 'C0494463', 'C0276496')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.object_cui, cp.object_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, lower(zx.covarname) as covarname FROM ZX AS zx ORDER BY theta1 desc;", sep = ""))
  #write.table(file = "covariates/rawPrecisionVars_for_AD.txt", x = precisionVars, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(outcomeEffects)
}


outcomeEffects <- getOutcomeEffects(predicates)


getExposureEffects <- function(predicates) {
  exposureEffects <- dbGetQuery(con, paste(" WITH ZX AS (SELECT cp.object_cui AS covar , COUNT(*) AS scnt, cp.object_name as covarname
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.object_cui, cp.object_name)
                        SELECT DISTINCT ZX.covar, zx.scnt as theta1, lower(zx.covarname) as covarname FROM ZX AS zx ORDER BY theta1 desc;", sep = ""))
  #write.table(file = "covariates/rawPrecisionVars_for_AD.txt", x = precisionVars, row.names = FALSE, col.names = FALSE, quote = FALSE)  
  return(exposureEffects)
}


exposureEffects <- getExposureEffects(predicates)

confounders[,4]
colliders[,4]
mediators[,4]
precisionVars[,3]


subjects <- dbGetQuery(con, paste(" SELECT DISTINCT cp.subject_cui AS covar , COUNT(*) AS scnt, cp.subject_name as covarname
                        FROM causalPredications cp, citations c
                        WHERE cp.pmid = c.pmid
                        AND CAST(pyear AS NUMERIC) >= 2010
                        AND cp.subject_cui IN ('C0011581', 'C1269683', 'C0871546', 'C0041696', 'C0588008', 'C0588006', 'C0221480', 'C0154409', 'C0024517')
                        AND cp.predicate IN ", predicates, "
                        GROUP BY cp.subject_cui, cp.subject_name
                        ORDER BY scnt DESC;
                        ", sep = ""))


print(confounders)
#confoundersCooked <- setdiff(setdiff(setdiff(confounders[,4], colliders[,4]), mediators[,4]), precisionVars[,3])
confoundersCooked <- setdiff(setdiff(confounders[,4], colliders[,4]), mediators[,4])
print(sort(confoundersCooked))

intersect(outcomeEffects[,3], confoundersCooked)
intersect(outcomeEffects[,3], precisionVarsCooked)

nrow(confounders)
length(confoundersCooked)

intersect(confounders[,4], mediators[,4])
#[1] "diabetes"                         "obesity"                          "hypertensive disease"             "diabetes mellitus"               
#[5] "hypercholesterolemia"             "metabolic syndrome"               "parkinson disease"                "cerebrovascular accident"        
#[9] "cardiovascular diseases"          "sleep disturbances"               "cerebral atrophy"                 "alzheimer's disease"             
#[13] "ischemic stroke"                  "heart failure"                    "rheumatoid arthritis"             "malnutrition"                    
#[17] "migraine disorders"               "virus diseases"                   "disability nos"                   "psoriasis"                       
#[21] "convulsions"                      "osteoporosis"                     "stress disorders, post-traumatic" "blind vision"                    
intersect(intersect(confounders[,4], mediators[,4]), colliders[,4])
#[1] "diabetes"                 "obesity"                  "hypertensive disease"     "diabetes mellitus"        "parkinson disease"        "cerebrovascular accident"
#[7] "cardiovascular diseases"  "sleep disturbances"       "cerebral atrophy"         "alzheimer's disease"      "ischemic stroke"          "rheumatoid arthritis"    
#[13] "malnutrition"             "disability nos"           "osteoporosis"            
#> 
intersect(intersect(confounders[,4], mediators[,4]), colliders[,4])
#[1] "diabetes"                 "obesity"                  "hypertensive disease"     "diabetes mellitus"        "parkinson disease"        "cerebrovascular accident"
#[7] "cardiovascular diseases"  "sleep disturbances"       "cerebral atrophy"         "alzheimer's disease"      "ischemic stroke"          "rheumatoid arthritis"    
#[13] "malnutrition"             "disability nos"           "osteoporosis"            
#> 
#  > 
#  > 


sort(setdiff(setdiff(confounders[,4], mediators[,4]), colliders[,4]))




setdiff(setdiff(mediators[,4], confounders[,4]), colliders[,4])
character(0)
setdiff(setdiff(mediators[,4], confounders[,4]), colliders[,4])
#Error in as.vector(y) : object 'confoundersrs' not found
setdiff(setdiff(mediators[,4], confounders[,4]), colliders[,4])
#[1] "septicemia"              "immune system diseases"  "wernicke encephalopathy"
#> 
#  > 
setdiff(setdiff(colliders[,4], confounders[,4]), mediators[,4])
#[1] "malignant neoplasms" "epilepsy"            "psychotic disorders" "suicide attempt"     "feeling suicidal"    "at risk for suicide"
#> 
#  > 
#  > 

#print confounder mediators
setdiff(intersect(confounders[,4], mediators[,4]), colliders[,4])
#[1] "hypercholesterolemia"             "metabolic syndrome"               "heart failure"                    "migraine disorders"              
#[5] "virus diseases"                   "psoriasis"                        "convulsions"                      "stress disorders, post-traumatic"
#[9] "blind vision" 


# print confounder colliders
setdiff(intersect(confounders[,4], colliders[,4]), mediators[,4])  

# print collider mediators
setdiff(intersect(colliders[,4], mediators[,4]), confounders[,4]) 




print(colliders)
collidersCooked <- setdiff(setdiff(colliders[,4], confounders[,4]), mediators[,4])
print(collidersCooked)

print(mediators)
mediatorsCooked <- setdiff(setdiff(mediators[,4], colliders[,4]), confounders[,4])
print(mediatorsCooked)

print(precisionVars)
precisionVarsCooked <- setdiff(setdiff(setdiff(precisionVars[,3], colliders[,4]), mediators[,4]), confounders[,4])
print(precisionVarsCooked)



#########

# NEXT: 
## run only with KB containing subset of clinical human AD studies >= 2010 (SemRep)
##
## Rerun SQL hitting only ScopedSemRep with same queries
##
##                                       # in both SemRep and KG               # in SemRep, not in KG              # not in SemRep, in KG           # not in Semrep or in KG
## confounders [in SemRepScoped    ]      
##             [Not in SemRepScoped]
##
## mediators   []
##
## colliders   []
##
## 
##
##




#####
# LEGACY




cuiLookup <- function(cuis) {
  covnames <- c()
  for (cui in cuis) {
    meaning <- dbGetQuery(con, paste("select distinct lower(covname) as covname from causalconcepts cc where cc.cui = '", cui, "'; ", sep = ""))
    #print(meaning[,1])
    covnames <- c(covnames, meaning[,1])
  }
  sort(unlist(covnames))
}
#cuiLookup(confounders)

getCovariates <- function(criterion, squelch) {
  
  fs <- list.files(path = "data/r3dataset/data/", pattern = ".csv")
  fs <- stringr::str_remove_all(string = fs, pattern = ".csv")
  
  ## run queries
  print("getting confounders")
  confounders <- getConfounders(predicates)
  colliders <- getColliders(predicates)
  mediators <- getMediators(predicates)
  precisionVars <- getPrecisionVars(predicates)
  
  if (criterion == "rawConfounders") {
    print("################################################")
    print(paste("#### Assembling covariates for CRITERION_", criterion))
    print("################################################")
    confounderNames <- confounders[,4]
    confounders <- confounders[,1]
    print("# of raw confounders")
    print(length(confounders))
    confounders <- intersect(confounders, fs)
    print("# of raw measured confounders")
    print(length(confounders))
    confounders <- setdiff(confounders, manuallyDefinedExclusions)
    print("# of confounders minus manually defined exclusions and synonyms")
    print(length(confounders))
    print("## the confounders ##")
    print(cuiLookup(confounders))
    print("############")
    covariates <- confounders
    
    print("# of raw covariates")
    print(length(covariates))
    
    print("## the covariates ##")
    print(cuiLookup(covariates))
    print("###################")
  } else if (criterion == "rawConfoundersPlusPrecisionVars") {
    #   criterion = "booka"
    print("################################################")
    print(paste("#### Assembling covariates for CRITERION_", criterion))
    print("################################################")
    confounders <- confounders[,1]
    print("# confounders raw")
    print(length(confounders))
    confounders <- intersect(confounders, fs)
    print("# measured variables")
    print(length(confounders))
    confounders <- setdiff(confounders, manuallyDefinedExclusions)
    print("# of confounders minus manually defined exclusions and synonyms")
    print(length(confounders))
    print("## the confounders ##")
    print(cuiLookup(confounders))
    print("############")
    
    precisionVars <- precisionVars[,1]
    print("# of raw precision variables")
    print(length(precisionVars))
    precisionVars <- intersect(precisionVars, fs)
    print("# of measured precision variables")
    print(length(precisionVars))
    precisionVars <- setdiff(precisionVars, manuallyDefinedExclusions)
    print("# of precision variables minus manually defined exclusions and synonyms")
    print(length(precisionVars))
    print("## the precision variables ##")
    print(cuiLookup(precisionVars))
    print("################")
    
    covariates = c(confounders, precisionVars)
    print("# of raw covariates")
    print(length(covariates))
    
    print("## the covariates ##")
    print(cuiLookup(covariates))
    print("###################")
    
  } else if (criterion == "cookedConfounders") {
    
    #    criterion = "booka"
    print("################################################")
    print(paste("#### Assembling covariates for CRITERION_", criterion))
    print("################################################")
    confounders <- confounders[,1]
    print("# confounders raw")
    print(length(confounders))
    confounders <- intersect(confounders, fs)
    print("# measured variables")
    print(length(confounders))
    confounders <- setdiff(confounders, manuallyDefinedExclusions)
    print("# of confounders minus manually defined exclusions and synonyms")
    print(length(confounders))
    print("## the confounders ##")
    print(cuiLookup(confounders))
    print("############")
    
    colliders <- colliders[,1]
    print("# of colliders raw")
    print(length(colliders))
    colliders <- intersect(colliders, fs)
    print("# of measured colliders")
    print(length(colliders))
    colliders <- setdiff(colliders, manuallyDefinedExclusions)
    print("# of colliders minus manually defined exclusions and synonyms")
    print(length(colliders))
    print("## the raw measured colliders ##")
    print(cuiLookup(colliders))
    print("##############")
    
    mediators <- mediators[,1]
    print("# of mediators raw")
    print(length(mediators))
    mediators <- intersect(mediators, fs)
    print("# of measured mediators")
    print(length(mediators))
    mediators <- setdiff(mediators, manuallyDefinedExclusions)
    print("# of mediators minus manually defined exclusions and synonyms")
    print(length(mediators))
    print("## the raw measured mediators ##")
    print(cuiLookup(mediators))
    print("############################")
    print("## ill-behaved covariates ##")
    print("############################")
    print("## the confounder/mediators according to SemMedDB ##")
    confounderMediators <- intersect(confounders, mediators)
    confounderMediators <- setdiff(confounderMediators, colliders)
    cuiLookup(confounderMediators)
    print("###############################")
    print("## the collider/confounders ##")
    print("##############################")
    colliderConfounders <- intersect(confounders, colliders)
    colliderConfounders <- setdiff(colliderConfounders, mediators)
    cuiLookup(colliderConfounders)
    print("###############################")
    print("## the collider/confounder/mediators (or 'chimeras') ##")
    print("##############################")
    chimeras <- intersect(confounderMediators, colliderConfounders)
    print(cuiLookup(chimeras))
    print("##############################")
    print("## well-behaved covariates ##")
    print("##############################")
    print("## the pure confounders ##")
    print("xxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    confounders <- setdiff(confounders, colliders)
    confounders <- setdiff(confounders, mediators)
    covariates <- confounders
    print(cuiLookup(confounders))
    
    print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    print("## the colliders ##")
    pureColliders <- setdiff(colliders, confounders)
    pureColliders <- setdiff(pureColliders, mediators)
    print(cuiLookup(pureColliders))
    print("## the mediators ##")
    pureMediators <- setdiff(mediators, confounders)
    pureMediators <- setdiff(mediators, confounders)
    print(cuiLookup(pureMediators))
    print("################################")
    print("# of raw covariates")
    print(length(covariates))
    
    print("## the covariates ##")
    print(cuiLookup(covariates))
    print("###################")
  } else if (criterion == "cookedConfoundersPlusPrecisionVars") {
    
    print("################################################")
    print(paste("#### Assembling covariates for CRITERION_", criterion))
    print("################################################")
    confounders <- confounders[,1]
    print("# confounders raw")
    print(length(confounders))
    confounders <- intersect(confounders, fs)
    print("# measured variables")
    print(length(confounders))
    confounders <- setdiff(confounders, manuallyDefinedExclusions)
    print("# of confounders minus manually defined exclusions and synonyms")
    print(length(confounders))
    print("## the confounders ##")
    print(cuiLookup(confounders))
    print("############")
    
    colliders <- colliders[,1]
    print("# of colliders raw")
    print(length(colliders))
    colliders <- intersect(colliders, fs)
    print("# of measured colliders")
    print(length(colliders))
    colliders <- setdiff(colliders, manuallyDefinedExclusions)
    print("# of colliders minus manually defined exclusions and synonyms")
    print(length(colliders))
    print("## the raw measured colliders ##")
    print(cuiLookup(colliders))
    print("##############")
    
    mediators <- mediators[,1]
    print("# of mediators raw")
    print(length(mediators))
    mediators <- intersect(mediators, fs)
    print("# of measured mediators")
    print(length(mediators))
    mediators <- setdiff(mediators, manuallyDefinedExclusions)
    print("# of mediators minus manually defined exclusions and synonyms")
    print(length(mediators))
    print("## the raw measured mediators ##")
    print(cuiLookup(mediators))
    print("############################")
    print("## ill-behaved covariates ##")
    print("############################")
    print("## the confounder/mediators according to SemMedDB ##")
    confounderMediators <- intersect(confounders, mediators)
    confounderMediators <- setdiff(confounderMediators, colliders)
    cuiLookup(confounderMediators)
    print("###############################")
    print("## the collider/confounders ##")
    print("##############################")
    colliderConfounders <- intersect(confounders, colliders)
    colliderConfounders <- setdiff(colliderConfounders, mediators)
    cuiLookup(colliderConfounders)
    print("###############################")
    print("## the collider/confounder/mediators (or 'chimeras') ##")
    print("##############################")
    chimeras <- intersect(confounderMediators, colliderConfounders)
    print(cuiLookup(chimeras))
    print("##############################")
    print("## well-behaved covariates ##")
    print("##############################")
    print("## the pure confounders ##")
    print("xxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    confounders <- setdiff(confounders, colliders)
    confounders <- setdiff(confounders, mediators)
    covariates <- confounders
    print(cuiLookup(confounders))
    
    print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
    print("## the colliders ##")
    pureColliders <- setdiff(colliders, confounders)
    pureColliders <- setdiff(pureColliders, mediators)
    print(cuiLookup(pureColliders))
    print("## the mediators ##")
    pureMediators <- setdiff(mediators, confounders)
    pureMediators <- setdiff(mediators, confounders)
    print(cuiLookup(pureMediators))
    print("################################")
    
    precisionVars <- precisionVars[,1]
    print("# of raw precision variables")
    print(length(precisionVars))
    precisionVars <- intersect(precisionVars, fs)
    print("# of measured precision variables")
    print(length(precisionVars))
    precisionVars <- setdiff(precisionVars, manuallyDefinedExclusions)
    print("# of precision variables minus manually defined exclusions and synonyms")
    print(length(precisionVars))
    print("## the precision variables ##")
    print(cuiLookup(precisionVars))
    print("################")
    
    confounders <- setdiff(confounders, colliders)
    confounders <- setdiff(confounders, mediators)
    precisionVars <- setdiff(precisionVars, colliders)
    precisionVars <- setdiff(precisionVars, mediators)
    
    covariates = c(confounders, precisionVars)
    print(length(covariates))
    
    
    print("# of raw covariates")
    print(length(covariates))
    
    print("## the covariates ##")
    print(cuiLookup(covariates))
    print("###################")
  }
  
  #print("MINUS variables not in data")
  print(length(covariates))
  # squelch <- 100
  goodCovars <- getGoodCovars(squelch)
  covariates <- intersect(covariates, goodCovars)
  print(paste("MINUS variables < ", squelch, " in data"))
  print(length(covariates))
  print("#####################")
  print("%%%%%% ADJ SET %%%%")
  print("###################")
  print(cuiLookup(covariates))
  print("###################")
  return(covariates)
}

gcc <- cmpfun(getCovariates)

rawConfoundersPlusPrecisionVars <- gcc("rawConfoundersPlusPrecisionVars", 100000)
rawConfoundersPlusPrecisionVars 
