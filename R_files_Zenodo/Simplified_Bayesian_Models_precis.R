###################################################################################################
#Bayesian hierarchical models for simplified models (ie. without microhabitat) to calculate trap variance

### This code requires rstan and rethinking. If using a Mac you may need to install XCode
###################################################################################################

######## Load packages
{library(rstan)
  library(rethinking) 
  library(devtools)
  library(lme4)
  library(dplyr)
  library(tidyr)}

######## Load all data
{morph <- read.csv("01_univMorph.csv")# morphological information
  trap <- read.csv("07.ecoloTrap.csv") # environmental information (substrate, flow, depth, etc.)
  parasite <- read.csv("02.sexParas.csv") # sex of fish
  PCNM_vectors <- read.csv("resultsPCNMvectors_allPCNMs.csv") # PCNMs to account for spatial distance
}

# Correct duplicates in trap data
trap <- trap[order(trap$fishID.univ),]
trap <- trap[-c(seq(1,39, by = 2)),]

##############################################################################################################################################################

##################################################################
### Size corrections for morphological traits, includes sex
##################################################################

#define the function to size-correct morphometric data
f.sizecorrect.withsex <- function(standard.length, vector.of.columns, dataframe.to.use, watershed, sex) {
  #Calculate overall mean standard length (Ls)
  Ls <- mean(standard.length, na.rm = TRUE) 
  #Call individual standard length
  L0 <- standard.length
  #Calculate common, within-group slope, b, for each trait. 
  b.vector.lmm <- vector()
  for (i in vector.of.columns) {
    abcd <- (dataframe.to.use[i])
    b.model <- lmer((abcd[,])~(standard.length)*sex + (1|watershed)) 
    ####lmer describes the fixed and random effects of the linear mixed model- watershed is random effect because of the bar notation | 
    b <- coef(summary(b.model))[2,1]
    b.vector.lmm <- c(b.vector.lmm, b)
  }
  # size correct
  xx <- dataframe.to.use  
  columnnames <- colnames(xx)
  j=1
  for (i in vector.of.columns) {
    M0 <- xx[,i] #grab the appropriate column of data
    Ms = M0 * ((Ls/L0)^b.vector.lmm[j]) #size correction formula
    j=j+1
    columnnames <- c(columnnames, paste(colnames(xx[i]), "sc", sep = "."))
    xx <- cbind(xx, Ms)
  }
  colnames(xx) <- columnnames # Rename the columns in the temporary dataframe xx
  return(xx) #Output a new dataframe with the name provided in "outputfilename"
}

### average across left and right pelvic spine length
morph$Pelvic.Spine.Length.Avg<-(morph$Left.Side.Pelvic.Spine.Length.mm+morph$Right.Side.Pelvic.Spine.Length.mm)/2

#size correct traits and make new data frame "morph.sc" where ".sc" is a size-corrected trait
morph.sc <- f.sizecorrect.withsex(morph$Standard.Length.mm, c(7:21, 27:36, 41,43,45, 51:55), morph, morph$watershed, parasite$Sex) #save traits we want to test
names(morph.sc) #row 55 is average pelvic spine length



##############################################################################################################################################################


##################################################################
### PCA of substrate and vegetation structure, for later use
##################################################################

#####Note, here we are only doing the simplified models, but this creates full data frame with habitat characteristics
# so that you could also do the microhabitat models with this data frame and so dimensions are consistent

#### make new columns for each letter (type of substrate) as T/F, use for PCA (A through J)
substrate.matrix <- data.frame(A =  grepl("A", trap$substrate),  B =  grepl("B", trap$substrate), C =  grepl("C", trap$substrate),
                               D =  grepl("D", trap$substrate),  E =  grepl("E", trap$substrate),  F =  grepl("F", trap$substrate),
                               G =  grepl("G", trap$substrate),  H =  grepl("H", trap$substrate),  I =  grepl("I", trap$substrate), 
                               J =  grepl("J", trap$substrate))
trap$substratePCA1 <- prcomp(substrate.matrix, scale = T)$x[,1] # First PCA axis
trap$substratePCA2 <- prcomp(substrate.matrix, scale = T)$x[,2] # Second PCA axis

#### do the same for vegetation matrix
vegetation.matrix <- data.frame(A =  grepl("A", trap$veg_structure),  D =  grepl("D", trap$veg_structure), E =  grepl("E", trap$veg_structure),
                                F =  grepl("F", trap$veg_structure),  Bh =  grepl("Bh", trap$veg_structure),  Cg =  grepl("Cg", trap$veg_structure),
                                Bg =  grepl("Bg", trap$veg_structure),  Blp =  grepl("Blp", trap$veg_structure),  Clp =  grepl("Clp", trap$veg_structure), 
                                Cros =  grepl("Cros", trap$veg_structure), Bspp =  grepl("Bspp", trap$veg_structure), Bsc =  grepl("Bsc", trap$veg_structure), 
                                Calg =  grepl("Calg", trap$veg_structure), Ch =  grepl("Ch", trap$veg_structure))
trap$vegetationPCA1 <- prcomp(vegetation.matrix, scale = T)$x[,1] # First PCA axis
trap$vegetationPCA2 <- prcomp(vegetation.matrix, scale = T)$x[,2] # Second PCA axis

names(trap) # in the trap data set you should now have PC 1 and 2 for substrate and for vegetation

####################################################################################################################################

##################################################################
### merge morphological and sex data (by fish) & habitat data (by trap)
##################################################################

### merge morph and trap data
dim(trap) #should be 2137
dim(morph) # should be 1463
morphtrap <- merge(morph.sc, trap, by = "fishID.univ", all.x = T)
dim(morphtrap)# should be 1463
names(morphtrap)


### make a new column in morphtrap and then in PCNM_vectors so I can merge the two data frames
head(morphtrap)
morphtrap$watershed_habitat_trap <- paste(morphtrap$watershed.x, morphtrap$habitat, morphtrap$trapID)
PCNM_vectors$watershed_habitat_trap <- paste(PCNM_vectors$watershed, PCNM_vectors$habitat, PCNM_vectors$Trap)
morphtrap<- merge (morphtrap, PCNM_vectors[,c(5,6,24)], by = "watershed_habitat_trap", all.x=T) #only saving PCNM 1 and 2

### merge in sex data from parasite dataframe
morphtrapsex <- merge(morphtrap, parasite[,c(1,4)], by = "fishID.univ", all.x = T) #adds in sex column

### Keep only the data where sex is known (some fish are questionable: M?, F?, NA) and the  convert to binary factor (0 or 1)
morphtrapsex <- subset(morphtrapsex, Sex == "F" | Sex == "M") #keep only rows where Sex column contains M or F value, some fish were not identified by sex
morphtrapsex$Sex<-factor(morphtrapsex$Sex) 
levels(morphtrapsex$Sex) <- c(1,0) #makes those levels into binary 1 or 0 (Female = 1)
morphtrapsex$Sex <- as.character(morphtrapsex$Sex)
morphtrapsex$Sex <- as.numeric(morphtrapsex$Sex) 

#choose the morph traits I will analyze
columns.to.analyze <- c(45, 47, 49, 57:82, 85:89) #includes average pelvic spine length
morph.sc.toanalyze <- morphtrapsex[,columns.to.analyze]



#############################################################################################################################################################################################

#ANALYSIS OF MORPHOLOGICAL TRAITS (distribution=NORMAL)

### There are 2 models for normally distributed, size-corrected traits. One for lake and one for stream.

### SEPARATE the lake and stream data sets 

#####STREAM
morph_stream  <- morph.sc.toanalyze[morphtrapsex$habitat == "stream",-c(1,2,3), ] #remove the count data in first three columns
trap_stream <- morphtrapsex[morphtrapsex$habitat == "stream",]
uniquetrap_stream <- paste(trap_stream$watershed.x, trap_stream$habitat, trap_stream$trapID)
uniquetrap_stream <- as.numeric(factor(uniquetrap_stream)) 
watershed_n_stream <- as.numeric(factor(trap_stream$watershed.x))
stream_metadata <- data.frame(uniquetrap_stream, watershed_n_stream, Sex = trap_stream$Sex)

results.storage.stream <- as.data.frame(matrix(nrow = length(columns.to.analyze), ncol = 26))


#####LAKE
morph_lake  <- morph.sc.toanalyze[morphtrapsex$habitat == "lake",-c(1,2,3)] #remove the count data in first three columns
trap_lake <- morphtrapsex[morphtrapsex$habitat == "lake",] 
uniquetrap_lake <- paste(trap_lake$watershed.x, trap_lake$habitat, trap_lake$trapID)
uniquetrap_lake <- as.numeric(factor(uniquetrap_lake))
watershed_n_lake <- as.numeric(factor(trap_lake$watershed.x))
lake_metadata <- data.frame(uniquetrap_lake, watershed_n_lake, Sex = trap_lake$Sex)

results.storage.lake <- as.data.frame(matrix(nrow = length(columns.to.analyze), ncol = 26))

for (i in 1:ncol(morph_stream)){
  
  trait <- colnames(morph_stream)[i] #trait name is the column name
  print(trait)
  
  modelstream <- alist( 
    Ystream ~ dnorm(mu, sigma),
    mu <- grandmean + a_sex*Sex + a_trap[uniquetrap_stream] + a_watershed[watershed_n_stream], 
    grandmean ~ dnorm(0,10),
    a_sex ~ dnorm(0,10),
    a_trap[uniquetrap_stream] ~ dnorm(0, sigma_trap),
    a_watershed[watershed_n_stream] ~ dnorm(0, sigma_watershed), 
    sigma_trap ~ dcauchy(0,1),
    sigma_watershed ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  )
  
  Ystream<-morph_stream[, i]
  stream.data.to.analyze <- na.omit(cbind(Ystream, stream_metadata))
  head(stream.data.to.analyze)
  model.out.stream <- map2stan(modelstream, data = stream.data.to.analyze, chains = 3, iter = 4000)
  results.stream <- precis(model.out.stream, digits = 5)
  results_i.stream <- c(trait, results.stream@output$Mean, results.stream@output$lower, results.stream@output$upper, results.stream@output$n_eff, results.stream@output$Rhat)
  results.storage.stream[i,] <- results_i.stream 
}


colnames(results.storage.stream) <- c("Trait", "Mean_grandmean", "Mean_a_sex", "Mean_sigma_trap", "Mean_sigma_watershed", "Mean_sigma", 
                                      "lower_grandmean",  "lower_a_sex", "lower_sigma_trap", "lower_sigma_watershed", "lower_sigma", 
                                      "upper_grandmean",  "upper_a_sex", "upper_sigma_trap", "upper_sigma_watershed", "upper_sigma",
                                      "n_eff_grandmean", "n_eff_a_sex", "n_eff_sigma_trap", "n_eff_sigma_watershed", "n_eff_sigma",
                                      "Rhat_grandmean", "Rhat_a_sex", "Rhat_sigma_trap", "Rhat_sigma_watershed", "Rhat_sigma")



##############################################################################
#LAKE MODEL
##############################################################################

for (i in 1:ncol(morph_lake)){
  
  trait <- colnames(morph_lake)[i] #trait name is the column name
  print(trait)
  
  modellake<- alist( 
    Ylake ~ dnorm(mu, sigma),
    mu <- grandmean +  a_sex*Sex + a_trap[uniquetrap_lake] + a_watershed[watershed_n_lake], 
    grandmean ~ dnorm(0,10), 
    a_sex ~ dnorm(0,10),
    a_trap[uniquetrap_lake] ~ dnorm(0, sigma_trap),
    a_watershed[watershed_n_lake] ~ dnorm(0, sigma_watershed), 
    sigma_trap ~ dcauchy(0,1),
    sigma_watershed ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  )
  
  Ylake<-morph_lake[, i]
  lake.data.to.analyze <- na.omit(cbind(Ylake, lake_metadata))
  head(lake.data.to.analyze)
  model.out.lake <- map2stan(modellake, data = lake.data.to.analyze, chains=3, iter=4000)
  results.lake <- precis(model.out.lake, digits = 5)
  results_i.lake <- c(trait, results.lake@output$Mean, results.lake@output$lower, results.lake@output$upper, results.lake@output$n_eff, results.lake@output$Rhat)
  results.storage.lake[i,] <- results_i.lake
} 


colnames(results.storage.lake) <- c("Trait", "Mean_grandmean", "Mean_a_sex", "Mean_sigma_trap", "Mean_sigma_watershed", "Mean_sigma", 
                                    "lower_grandmean",  "lower_a_sex", "lower_sigma_trap", "lower_sigma_watershed", "lower_sigma", 
                                    "upper_grandmean",  "upper_a_sex", "upper_sigma_trap", "upper_sigma_watershed", "upper_sigma",
                                    "n_eff_grandmean",  "n_eff_a_sex", "n_eff_sigma_trap", "n_eff_sigma_watershed", "n_eff_sigma",
                                    "Rhat_grandmean",   "Rhat_a_sex", "Rhat_sigma_trap", "Rhat_sigma_watershed", "Rhat_sigma")


#############################################################################################################################################################################################




#############################################################################################################################################################################################

#ANALYSIS OF MORPHOLOGICAL TRAITS (distribution=POISSON)

#############################################################################################################################################################################################

#DEFINE THE DATA SETS FOR THE COUNT DATA (plate number left and right and gill raker number)


#####STREAM
morph_stream_count  <- morph.sc.toanalyze[morphtrapsex$habitat == "stream", c(1,2,3)] #only keep plate count and gill raker count data
trap_stream <- morphtrapsex[morphtrapsex$habitat == "stream",]
uniquetrap_stream <- paste(trap_stream$watershed.x, trap_stream$habitat, trap_stream$trapID)
uniquetrap_stream <- as.numeric(factor(uniquetrap_stream)) 
watershed_n_stream <- as.numeric(factor(trap_stream$watershed.x))
stream_metadata_count <- data.frame(uniquetrap_stream, watershed_n_stream, Sex = trap_stream$Sex)

results.storage.stream.count <- as.data.frame(matrix(nrow = length(columns.to.analyze), ncol = 21)) 


#####LAKE
morph_lake_count  <- morph.sc.toanalyze[morphtrapsex$habitat == "lake", c(1,2,3)] #only want the first three columns 
trap_lake <- morphtrapsex[morphtrapsex$habitat == "lake",] 
uniquetrap_lake <- paste(trap_lake$watershed.x, trap_lake$habitat, trap_lake$trapID)
uniquetrap_lake <- as.numeric(factor(uniquetrap_lake))
watershed_n_lake <- as.numeric(factor(trap_lake$watershed.x))
lake_metadata_count <- data.frame(uniquetrap_lake, watershed_n_lake, Sex = trap_lake$Sex)

results.storage.lake.count <- as.data.frame(matrix(nrow = length(columns.to.analyze), ncol = 21)) 



##############################################################################
#STREAM COUNT MODEL
##############################################################################

for (i in 1:ncol(morph_stream_count)){
  trait <- colnames(morph_stream_count)[i] #trait name is the column name
  print(trait)
  
  modelstream_count <- alist( 
    Ystream ~ dpois(lambda), #using Poisson distribution for count data, no mu or sigma, only lambda. Log makes lambda linear
    log(lambda) <- grandmean + a_sex*Sex + a_trap[uniquetrap_stream] + a_watershed[watershed_n_stream], 
    grandmean ~ dnorm(0,10), 
    a_sex ~ dnorm(0,10),
    a_trap[uniquetrap_stream] ~ dnorm(0, sigma_trap),# mu_trap can go in for the 0
    a_watershed[watershed_n_stream] ~ dnorm(0, sigma_watershed), 
    sigma_trap ~ dcauchy(0,1),
    sigma_watershed ~ dcauchy(0,1)
  )
  
  Ystream<-morph_stream_count[, i]
  stream.data.to.analyze <- na.omit(cbind(Ystream, stream_metadata_count))
  head(stream.data.to.analyze)
  
  
  ### Add in starting values
  startingvals<- c(grandmean = mean (log(stream.data.to.analyze$Ystream)),
                   a_sex = summary(lm(log(stream.data.to.analyze$Ystream)~stream.data.to.analyze$Sex))
                   $coef[2,1],
                   sigma_trap = 0.1,
                   sigma_watershed = 0.1)
  
  model.out.stream.count <- map2stan(modelstream_count, data = stream.data.to.analyze, chains = 3, iter = 4000, start= startingvals) 
  results.stream.count <- precis(model.out.stream.count, digits = 5)
  results_i.stream.count <- c(trait, results.stream.count@output$Mean, results.stream.count@output$lower, results.stream.count@output$upper, results.stream.count@output$n_eff, results.stream.count@output$Rhat)
  results.storage.stream.count[i,] <- results_i.stream.count }


colnames(results.storage.stream.count) <- c("Trait", "Mean_grandmean", "Mean_a_sex", "Mean_sigma_trap", "Mean_sigma_watershed", 
                                            "lower_grandmean",  "lower_a_sex", "lower_sigma_trap", "lower_sigma_watershed", 
                                            "upper_grandmean",   "upper_a_sex", "upper_sigma_trap", "upper_sigma_watershed",
                                            "n_eff_grandmean", "n_eff_a_sex", "n_eff_sigma_trap", "n_eff_sigma_watershed",
                                            "Rhat_grandmean",  "Rhat_a_sex", "Rhat_sigma_trap", "Rhat_sigma_watershed")


##############################################################################
#LAKE COUNT MODEL
##############################################################################

for (i in 1:ncol(morph_lake_count)){
  trait <- colnames(morph_lake_count)[i] #trait name is the column name
  print(trait)
  
  modellake_count <- alist( 
    Ylake ~ dpois(lambda),
    log(lambda) <- grandmean + a_sex*Sex +  a_trap[uniquetrap_lake] + a_watershed[watershed_n_lake], 
    grandmean ~ dnorm(0,10), 
    a_sex ~ dnorm(0,10),
    a_trap[uniquetrap_lake] ~ dnorm(0, sigma_trap),# mu_trap can go in for the 0
    a_watershed[watershed_n_lake] ~ dnorm(0, sigma_watershed), 
    sigma_trap ~ dcauchy(0,1),
    sigma_watershed ~ dcauchy(0,1)
  )
  
  Ylake<-morph_lake_count[, i]
  lake.data.to.analyze <- na.omit(cbind(Ylake, lake_metadata_count))
  head(lake.data.to.analyze)
  
  
  ### Add in starting values
  startingvals<- c(grandmean = mean (log(Ylake), na.rm =T),
                   a_sex = summary(lm(log(Ylake)~lake_metadata_count$Sex))
                   $coef[2,1],
                   sigma_trap = 0.1,
                   sigma_watershed = 0.1)
  
  model.out.lake.count <- map2stan(modellake_count, data = lake.data.to.analyze, chains=3, iter=4000, start=startingvals)
  results.lake.count <- precis(model.out.lake.count, digits = 5)
  results_i.lake.count <- c(trait, results.lake.count@output$Mean, results.lake.count@output$lower, results.lake.count@output$upper, results.lake.count@output$n_eff, results.lake.count@output$Rhat)
  results.storage.lake.count[i,] <- results_i.lake.count
} 



colnames(results.storage.lake.count) <- c("Trait", "Mean_grandmean",  "Mean_a_sex", "Mean_sigma_trap", "Mean_sigma_watershed", 
                                          "lower_grandmean",  "lower_a_sex", "lower_sigma_trap", "lower_sigma_watershed", 
                                          "upper_grandmean", "upper_a_sex", "upper_sigma_trap", "upper_sigma_watershed",
                                          "n_eff_grandmean",  "n_eff_a_sex",  "n_eff_sigma_trap", "n_eff_sigma_watershed",
                                          "Rhat_grandmean", "Rhat_a_sex","Rhat_sigma_trap", "Rhat_sigma_watershed")


