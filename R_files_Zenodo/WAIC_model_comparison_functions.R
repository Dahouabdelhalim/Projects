####################################################################################################################################
# Does sigma trap term improve the model? Using WAIC comparison of models with and without sigma trap included
# we did not include the habitat characteristics in our model comparisons
####################################################################################################################################


##################################################################
### load packages

{ library(rstan)
  library(rethinking) 
  library(devtools)
  library(lme4)
  library(dplyr)
  library(tidyr)}

### load raw data
{ morph <- read.csv("01_univMorph.csv") #contains morph data
  trap <- read.csv("07.ecoloTrap.csv") #This csv contains the environmental information (substrate, flow, depth, etc.)
  parasite <- read.csv("02.sexParas.csv") #contains sex
  PCNM_vectors <- read.csv("resultsPCNMvectors_allPCNMs.csv") #This csv contains the PCNMs to account for spatial distance
}

# Correct some duplicates in trap data
trap <- trap[order(trap$fishID.univ),]
trap <- trap[-c(seq(1,39, by = 2)),]

##############################################################################################################################################################

##################################################################
#### Size corrections for morphological traits 
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

###change trait to average pelvic spine length
morph$Pelvic.Spine.Length.Avg<-(morph$Left.Side.Pelvic.Spine.Length.mm+morph$Right.Side.Pelvic.Spine.Length.mm)/2

morph.sc <- f.sizecorrect.withsex(morph$Standard.Length.mm, c(7:21, 27:36, 41,43,45, 51:55), morph, morph$watershed, parasite$Sex)
names(morph.sc) #row 55 is average pelvic spine length


##############################################################################################################################################################

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

### Keep only the data where sex is known (some fish are questionable: M?, F?, NA) and convert to binary 0 or 1
morphtrapsex <- subset(morphtrapsex, Sex == "F" | Sex == "M") #keep only rows where Sex column contains M or F value
morphtrapsex$Sex<-factor(morphtrapsex$Sex) #eliminates the extra levels so I now only have 2, M or F
levels(morphtrapsex$Sex) <- c(1,0) #makes those levels into binary 1 or 0 (F=1)
morphtrapsex$Sex<-as.character(morphtrapsex$Sex)
morphtrapsex$Sex<-as.numeric(morphtrapsex$Sex)

#choose the morph traits I will analyze
columns.to.analyze <- c(45, 47, 49, 57:82, 85:89) #includes average pelvic spine length
morph.sc.toanalyze <- morphtrapsex[,columns.to.analyze]

#########################################################################################################
#Older MODEL WITHOUT THE HABITAT CHARACTERISTICS
#########################################################################################################

##########################
#separate lake and stream morph and trap data so you can run two separate mdodels on each habitat type

#####STREAM
morph_stream  <- morph.sc.toanalyze[morphtrapsex$habitat == "stream",]
trap_stream <- morphtrapsex[morphtrapsex$habitat == "stream",]
uniquetrap_stream <- paste(trap_stream$watershed.x, trap_stream$habitat, trap_stream$trapID)
uniquetrap_stream <- as.numeric(factor(uniquetrap_stream)) 
watershed_n_stream <- as.numeric(factor(trap_stream$watershed.x))
stream_metadata <- data.frame(uniquetrap_stream, watershed_n_stream, Sex = trap_stream$Sex)


##########################
#can use these functions to calculate WAIC comparisons by trait. 

#model with trap variance term for count data (so trait_number = 1, 2 or 3)
model_withtrap_count_stream <- function(trait_number) {
  trait <- colnames(morph_stream)[trait_number] 
  print(trait)
  
  modelstream <- alist( 
    Ystream ~ dpois(lambda),
    log(lambda) <- grandmean + a_trap[uniquetrap_stream] + a_watershed[watershed_n_stream] + a_sex*Sex, 
    grandmean ~ dnorm(0,10), 
    a_sex ~ dnorm(0,10),
    a_trap[uniquetrap_stream] ~ dnorm(0, sigma_trap),
    a_watershed[watershed_n_stream] ~ dnorm(0, sigma_watershed), 
    sigma_trap ~ dcauchy(0,1),
    sigma_watershed ~ dcauchy(0,1)
  )
  
  Ystream<-morph_stream[, trait_number]
  stream.data.to.analyze <- data.frame(cbind(Ystream, stream_metadata))
  stream.data.to.analyze <-na.omit(stream.data.to.analyze)
  
  startingvals<- c(grandmean = mean(log(stream.data.to.analyze$Ystream)),
                   a_sex = summary(lm(log(stream.data.to.analyze$Ystream)~stream.data.to.analyze$Sex))
                   $coef[2,1],
                   sigma_watershed = 0.1,
                   sigma_trap = 0.1)
  
 
  
  output <- map2stan(modelstream, data = stream.data.to.analyze, start=startingvals)
  return(output)
}


#model without trap variance term for count data (so trait_number = 1, 2 or 3)
model_withouttrap_count_stream <- function(trait_number) {
trait <- colnames(morph_stream)[trait_number] 
print(trait)

modelstream <- alist( 
  Ystream ~ dpois(lambda),
  log(lambda) <- grandmean  + a_watershed[watershed_n_stream] + a_sex*Sex, 
  grandmean ~ dnorm(0,10),
  a_watershed[watershed_n_stream] ~ dnorm(0, sigma_watershed), 
  a_sex ~ dnorm(0,10),
  sigma_watershed ~ dcauchy(0,1)
)

Ystream<-morph_stream[, trait_number]
stream.data.to.analyze <- data.frame(cbind(Ystream, stream_metadata))
stream.data.to.analyze <-na.omit(stream.data.to.analyze)

startingvals<- c(grandmean = mean(log(stream.data.to.analyze$Ystream)),
                 a_sex = summary(lm(log(stream.data.to.analyze$Ystream)~stream.data.to.analyze$Sex))
                 $coef[2,1],
                 sigma_watershed = 0.1,
                 sigma_trap = 0.1)


output<- map2stan(modelstream, data = stream.data.to.analyze, start=startingvals)
return(output)
}


#Example for trait 1:
model.out1<- model_withtrap_count_stream(1)
model.out1WAIC<- model_withouttrap_count_stream(1)

#we did this for traits 1, 2, 3

################################################################################################################

##### For normally distributed traits... 

model_withtrap_stream <- function(trait_number) {
  trait <- colnames(morph_stream)[trait_number] 
  print(trait)
  modelstream <- alist( 
    Ystream ~ dnorm(mu, sigma),
    mu <- grandmean + a_trap[uniquetrap_stream] + a_watershed[watershed_n_stream] + a_sex*Sex, 
    grandmean ~ dnorm(0,10), 
    a_sex ~ dnorm(0,10),
    a_trap[uniquetrap_stream] ~ dnorm(0, sigma_trap),
    a_watershed[watershed_n_stream] ~ dnorm(0, sigma_watershed), 
    sigma_trap ~ dcauchy(0,1),
    sigma_watershed ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  )
  
  Ystream<-morph_stream[, trait_number]
  stream.data.to.analyze <- data.frame(cbind(Ystream, stream_metadata))
  stream.data.to.analyze <-na.omit(stream.data.to.analyze)
  output <- map2stan(modelstream, data = stream.data.to.analyze)
  return(output)
}

  
model_withouttrap_stream <- function(trait_number) {
  trait <- colnames(morph_stream)[trait_number] 
  print(trait)
  modelstream <- alist( 
    Ystream ~ dnorm(mu, sigma),
    mu <- grandmean + a_watershed[watershed_n_stream] + a_sex*Sex, 
    grandmean ~ dnorm(0,10), 
    a_sex ~ dnorm(0,10),
    a_watershed[watershed_n_stream] ~ dnorm(0, sigma_watershed), 
    sigma_watershed ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  )
  
  Ystream<-morph_stream[, trait_number]
  stream.data.to.analyze <- data.frame(cbind(Ystream, stream_metadata))
  stream.data.to.analyze <-na.omit(stream.data.to.analyze)
  output <- map2stan(modelstream, data = stream.data.to.analyze)
  return(output)
}
  

### do this for all normal traits (4-34)
  model.out4<- model_withtrap_stream(4)
  model.out4WAIC<- model_withouttrap_stream(4)
  
  
 
### Get WAIC values for streams
  
  results.storageWAIC_stream <- as.data.frame(matrix(nrow = 68, ncol = 26)) 
  
  
  streamWAIC<- compare(model.out1, model.out1WAIC, model.out2, model.out2WAIC, model.out3, model.out3WAIC, model.out4, model.out4WAIC, model.out5, model.out5WAIC, model.out6, model.out6WAIC, model.out7,
          model.out7WAIC, model.out8, model.out8WAIC, model.out9, model.out9WAIC,model.out10, model.out10WAIC,model.out11, model.out11WAIC,model.out12, model.out12WAIC ,model.out13, model.out13WAIC ,model.out14, model.out14WAIC,
          model.out15, model.out15WAIC,model.out16, model.out16WAIC,model.out17, model.out17WAIC,model.out18, model.out18WAIC,model.out19, model.out19WAIC, model.out20, model.out20WAIC,model.out21, model.out21WAIC,model.out22, model.out22WAIC,
          model.out23, model.out23WAIC,model.out24, model.out24WAIC,model.out25, model.out25WAIC,model.out26, model.out26WAIC,model.out27, model.out27WAIC,model.out28, model.out28WAIC,model.out29, model.out29WAIC,model.out30, model.out30WAIC,
          model.out31, model.out31WAIC,model.out32, model.out32WAIC,model.out33, model.out33WAIC,model.out34, model.out34WAIC)
  
  results.storageWAIC_stream<- streamWAIC@output
  plot(streamWAIC, SE=TRUE, dSE=TRUE)

  
##################################################################################################################################

  
  
  
##################################################################################################################################

  ###FOR LAKE TRAITS

##########################
  
###LAKE
morph_lake  <- morph.sc.toanalyze[morphtrapsex$habitat == "lake",]
trap_lake <- morphtrapsex[morphtrapsex$habitat == "lake",]
uniquetrap_lake <- paste(trap_lake$watershed.x, trap_lake$habitat, trap_lake$trapID)
uniquetrap_lake <- as.numeric(factor(uniquetrap_lake))
watershed_n_lake <- as.numeric(factor(trap_lake$watershed.x))
lake_metadata <- data.frame(uniquetrap_lake, watershed_n_lake, Sex = trap_lake$Sex)

##########################

##########################
#can use these functions to calculate WAIC comparisons by trait. 

#model with trap variance term for count data (so trait_number = 1, 2 or 3)
model_withtrap_count_lake <- function(trait_number) {
  trait <- colnames(morph_lake)[trait_number] 
  print(trait)
  
  modellake <- alist( 
    Ylake ~ dpois(lambda),
    log(lambda) <- grandmean + a_trap[uniquetrap_lake] + a_watershed[watershed_n_lake] + a_sex*Sex, 
    grandmean ~ dnorm(0,10), 
    a_sex ~ dnorm(0,10),
    a_trap[uniquetrap_lake] ~ dnorm(0, sigma_trap),
    a_watershed[watershed_n_lake] ~ dnorm(0, sigma_watershed), 
    sigma_trap ~ dcauchy(0,1),
    sigma_watershed ~ dcauchy(0,1)
  )
  
  Ylake<-morph_lake[, trait_number]
  lake.data.to.analyze <- data.frame(cbind(Ylake, lake_metadata))
  lake.data.to.analyze <-na.omit(lake.data.to.analyze)
  
  startingvals<- c(grandmean = mean(log(lake.data.to.analyze$Ylake)),
                   a_sex = summary(lm(log(lake.data.to.analyze$Ylake)~lake.data.to.analyze$Sex))
                   $coef[2,1],
                   sigma_watershed = 0.1,
                   sigma_trap = 0.1)
  
  output <- map2stan(modellake, data = lake.data.to.analyze, start=startingvals)
  return(output)
}


#model without trap variance term for count data (so trait_number = 1, 2 or 3)
model_withouttrap_count_lake <- function(trait_number) {
  trait <- colnames(morph_lake)[trait_number] 
  print(trait)
  
  modellake <- alist( 
    Ylake ~ dpois(lambda),
    log(lambda) <- grandmean  + a_watershed[watershed_n_lake] + a_sex*Sex, 
    grandmean ~ dnorm(0,10),
    a_watershed[watershed_n_lake] ~ dnorm(0, sigma_watershed), 
    a_sex ~ dnorm(0,10),
    sigma_watershed ~ dcauchy(0,1)
  )
  
  Ylake<-morph_lake[, trait_number]
  lake.data.to.analyze <- data.frame(cbind(Ylake, lake_metadata))
  lake.data.to.analyze <-na.omit(lake.data.to.analyze)
  
  startingvals<- c(grandmean = mean(log(lake.data.to.analyze$Ylake)),
                   a_sex = summary(lm(log(lake.data.to.analyze$Ylake)~lake.data.to.analyze$Sex))
                   $coef[2,1],
                   sigma_watershed = 0.1,
                   sigma_trap = 0.1)
  
 
  output<- map2stan(modellake, data = lake.data.to.analyze, start=startingvals)
  return(output)
}


#Example for trait 1:
model.out1<- model_withtrap_count_lake(1)
model.out1WAIC<- model_withouttrap_count_lake(1)
#we did this for traits 1, 2, 3

################################################################################################################

##### For normally distributed traits... 

model_withtrap_lake <- function(trait_number) {
  trait <- colnames(morph_lake)[trait_number] 
  print(trait)
  modellake <- alist( 
    Ylake ~ dnorm(mu, sigma),
    mu <- grandmean + a_trap[uniquetrap_lake] + a_watershed[watershed_n_lake] + a_sex*Sex, 
    grandmean ~ dnorm(0,10), 
    a_sex ~ dnorm(0,10),
    a_trap[uniquetrap_lake] ~ dnorm(0, sigma_trap),
    a_watershed[watershed_n_lake] ~ dnorm(0, sigma_watershed), 
    sigma_trap ~ dcauchy(0,1),
    sigma_watershed ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  )
  
  Ylake<-morph_lake[, trait_number]
  lake.data.to.analyze <- data.frame(cbind(Ylake, lake_metadata))
  lake.data.to.analyze <-na.omit(lake.data.to.analyze)
  output <- map2stan(modellake, data = lake.data.to.analyze)
  return(output)
}


model_withouttrap_lake <- function(trait_number) {
  trait <- colnames(morph_lake)[trait_number] 
  print(trait)
  modellake <- alist( 
    Ylake ~ dnorm(mu, sigma),
    mu <- grandmean + a_watershed[watershed_n_lake] + a_sex*Sex, 
    grandmean ~ dnorm(0,10), 
    a_sex ~ dnorm(0,10),
    a_watershed[watershed_n_lake] ~ dnorm(0, sigma_watershed), 
    sigma_watershed ~ dcauchy(0,1),
    sigma ~ dcauchy(0,1)
  )
  
  Ylake<-morph_lake[, trait_number]
  lake.data.to.analyze <- data.frame(cbind(Ylake, lake_metadata))
  lake.data.to.analyze <-na.omit(lake.data.to.analyze)
  output <- map2stan(modellake, data = lake.data.to.analyze)
  return(output)
}


### do this for all normal traits (4-34)
model.out4<- model_withtrap_lake(4)
model.out4WAIC<- model_withouttrap_lake(4)

########################################################################################################################

### get WAIC comparisons

results.storageWAIC_lake <- as.data.frame(matrix(nrow = 68, ncol = 26))


lakeWAIC<- compare(model.out1, model.out1WAIC, model.out2, model.out2WAIC, model.out3, model.out3WAIC, model.out4, model.out4WAIC, model.out5, model.out5WAIC, model.out6, model.out6WAIC, model.out7,
                     model.out7WAIC, model.out8, model.out8WAIC, model.out9, model.out9WAIC,model.out10, model.out10WAIC,model.out11, model.out11WAIC,model.out12, model.out12WAIC ,model.out13, model.out13WAIC ,model.out14, model.out14WAIC,
                     model.out15, model.out15WAIC,model.out16, model.out16WAIC,model.out17, model.out17WAIC,model.out18, model.out18WAIC,model.out19, model.out19WAIC, model.out20, model.out20WAIC,model.out21, model.out21WAIC,model.out22, model.out22WAIC,
                     model.out23, model.out23WAIC,model.out24, model.out24WAIC,model.out25, model.out25WAIC,model.out26, model.out26WAIC,model.out27, model.out27WAIC,model.out28, model.out28WAIC,model.out29, model.out29WAIC,model.out30, model.out30WAIC,
                     model.out31, model.out31WAIC,model.out32, model.out32WAIC,model.out33, model.out33WAIC,model.out34, model.out34WAIC)

results.storageWAIC_lake<- lakeWAIC@output
plot(lakeWAIC, SE=TRUE, dSE=TRUE)

