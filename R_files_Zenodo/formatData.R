require(reshape2)
merged <- read.delim("mergedUTM.tsv", stringsAsFactors=F)

# Make distance matrices for each population, then melt them into long form
# data frames

# Create a new matrix for the UTM coordinates and name the rows of the matrix
# by plant
preBurk1 <- as.matrix(cbind(merged$Northing[merged$population=="Burk1"], merged$Easting[merged$population=="Burk1"]))
rownames(preBurk1) <- merged$name[merged$population=="Burk1"]  

# Make a matrix of all pairwise distances between plants
distBurk1 <- as.matrix(dist(preBurk1))

# Use the melt function to convert the distance matrix to a long-form data frame
meltBurk1 <- melt(distBurk1)
burk1Out <- meltBurk1

# Create new variables for the sex of the focal and partner plants
burk1Out$focalSex <- NA
burk1Out$partnerSex <- NA
colnames(burk1Out) <- c("focal", "partner", "distance", "focalSex", "partnerSex")
burk1Out$focal <- as.character(burk1Out$focal)
burk1Out$partner <- as.character(burk1Out$partner)

# Iterate through the long-form data frame and add in the sex of focal and 
# partners
for (i in 1:nrow(burk1Out)) {
  fName <- burk1Out$focal[i]
  pName <- burk1Out$partner[i]
  burk1Out$focalSex[i] <- unique(merged$sex[merged$name==fName])
  burk1Out$partnerSex[i] <- unique(merged$sex[merged$name==pName])
}

# Create a new matrix for the UTM coordinates and name the rows of the matrix
# by plant (these are the same steps, but repeated for the Burk2 pop)
preBurk2 <- as.matrix(cbind(merged$Northing[merged$population=="Burk2"], merged$Easting[merged$population=="Burk2"]))
rownames(preBurk2) <- merged$name[merged$population=="Burk2"]  

# Make a matrix of all pairwise distances between plants
distBurk2 <- as.matrix(dist(preBurk2))

# Use the melt function to convert the distance matrix to a long-form data frame
meltBurk2 <- melt(distBurk2)
burk2Out <- meltBurk2

# Create new variables for the sex of the focal and partner plants
burk2Out$focalSex <- NA
burk2Out$partnerSex <- NA
colnames(burk2Out) <- c("focal", "partner", "distance", "focalSex", "partnerSex")
burk2Out$focal <- as.character(burk2Out$focal)
burk2Out$partner <- as.character(burk2Out$partner)

# Iterate through the long-form data frame and add in the sex of focal and 
# partners
for (i in 1:nrow(burk2Out)) {
  fName <- as.character(burk2Out$focal[i]) 
  pName <- as.character(burk2Out$partner[i])
  burk2Out$focalSex[i] <- unique(merged$sex[merged$name==fName])
  burk2Out$partnerSex[i] <- unique(merged$sex[merged$name==pName])
}

# combine all of the long-form data frames into one long data frame
longestFrame <- rbind(burk1Out, burk2Out)

# make a new data frame that has rows for all individuals in all populations
localSex <- as.character(merged$name)
localSex <- as.data.frame(localSex)
colnames(localSex) <- "name"

# initilize variables for local ratios and pop sizes (laborious)
localSex$sex <- NA
localSex$population <- NA
localSex$sexRatio <- NA
localSex$popSize <- NA
localSex$r025 <- NA
localSex$r05 <- NA
localSex$r10 <- NA
localSex$r15 <- NA
localSex$r20 <- NA
localSex$r25 <- NA
localSex$r30 <- NA
localSex$r35 <- NA
localSex$r40 <- NA
localSex$r45 <- NA
localSex$r50 <- NA
localSex$r55 <- NA
localSex$r60 <- NA

localSex$d025 <- NA
localSex$d05 <- NA
localSex$d10 <- NA
localSex$d15 <- NA
localSex$d20 <- NA
localSex$d25 <- NA
localSex$d30 <- NA
localSex$d35 <- NA
localSex$d40 <- NA
localSex$d45 <- NA
localSex$d50 <- NA
localSex$d55 <- NA
localSex$d60 <- NA


# look up sex and population name  from the source data frame
for (i in 1:nrow(localSex)) {
  name <- as.character(localSex$name[i])
  localSex$sex[localSex$name == name] <- as.character(merged$sex[merged$name == name])
  localSex$population[localSex$name == name] <- as.character(merged$population[merged$name == name])
}

# convert sex, and population into factors
localSex$sex <- factor(localSex$sex, levels=c("hermaphrodite", "female"))
localSex$population <- factor(localSex$population, levels=c("Burk1", "Burk2"))

# calculate sex ratio and pop size for each population
for (i in 1:length(levels(localSex$population))) {
  nHerm = 0
  nFem = 0
  pop = levels(localSex$population)[i]
  subByPop <- localSex[localSex$population==pop,]
  popSize = nrow(subByPop)
  for (j in 1:nrow(subByPop)) {
    if (subByPop$sex[j] == "hermaphrodite") {
      nHerm <- nHerm + 1
    } 
    
    if (subByPop$sex[j] == "female") {
      nFem <- nFem + 1
    }
  }
  localSex$sexRatio[localSex$population==pop] <- (nHerm / (nHerm + nFem))
  localSex$popSize[localSex$population==pop] <- popSize
}

# iterate through the list, subset all of the distances for that individual
# calculate the local sex ratios and densities at the different radii


for (i in 1:nrow(localSex)) {
  # A vector to specify the radius lengths 
  distVector = c(0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6)
  # A vector to store the individual sex ratios across spatial scales
  ratVector = rep(NA, length(distVector))
  # A vector to store individual sex ratios weighted by fruit number
  #fRatVector = rep(NA, length(distVector))
  # A vector for local population sizes
  denVector = rep(NA, length(distVector))
  # mean fitness
  #meanWVector = rep(NA, length(distVector))
  # relative fitness
  #relWVector = rep(NA, length(distVector))
  # The name of the individual this iteration
  name = localSex$name[i]
  # Subset the long-form data frame by individual name
  subByName <- longestFrame[longestFrame$focal==name,]
  # For each individual, iterate through all of the different spatial
  # scales, and calculate local sex ratios and population sizes
  for (j in 1:length(distVector)) {
    nHerm = 0
    nFem = 0
    fHerm = 0
    fFem = 0
    nfHerm = 0
    nfFem = 0
    subSub <- subByName[subByName$distance > 0 & subByName$distance < distVector[j],]
    for (k in 1:nrow(subSub)) {
      if (nrow(subSub) == 0) {
        nHerm <- NA
        pName <- NA
        nFem <- NA
      } else {
        if (subSub$partnerSex[k] == "hermaphrodite") {
          nHerm <- nHerm + 1
          pName <- subSub$partner[k]
          
        } 
        
        
        if (subSub$partnerSex[k] == "female") {
          nFem <- nFem + 1
          pName <- subSub$partner[k]
        }
      }
    }
    ratVector[j] <- (nHerm / (nHerm + nFem))
    denVector[j] <- (nHerm + nFem)
    
    # Transfer the local pop sizes and ratios to the
    # working wide-form data frame
    localSex$r025[localSex$name==name] <- ratVector[1]
    localSex$r05[localSex$name==name] <- ratVector[2]
    localSex$r10[localSex$name==name] <- ratVector[3]
    localSex$r15[localSex$name==name] <- ratVector[4]
    localSex$r20[localSex$name==name] <- ratVector[5]
    localSex$r25[localSex$name==name] <- ratVector[6]
    localSex$r30[localSex$name==name] <- ratVector[7]
    localSex$r35[localSex$name==name] <- ratVector[8]
    localSex$r40[localSex$name==name] <- ratVector[9]
    localSex$r45[localSex$name==name] <- ratVector[10]
    localSex$r50[localSex$name==name] <- ratVector[11]
    localSex$r55[localSex$name==name] <- ratVector[12]
    localSex$r60[localSex$name==name] <- ratVector[13]
    
    localSex$d025[localSex$name==name] <- denVector[1]
    localSex$d05[localSex$name==name] <- denVector[2]
    localSex$d10[localSex$name==name] <- denVector[3]
    localSex$d15[localSex$name==name] <- denVector[4]
    localSex$d20[localSex$name==name] <- denVector[5]
    localSex$d25[localSex$name==name] <- denVector[6]
    localSex$d30[localSex$name==name] <- denVector[7]
    localSex$d35[localSex$name==name] <- denVector[8]
    localSex$d40[localSex$name==name] <- denVector[9]
    localSex$d45[localSex$name==name] <- denVector[10]
    localSex$d50[localSex$name==name] <- denVector[11]
    localSex$d55[localSex$name==name] <- denVector[12]
    localSex$d60[localSex$name==name] <- denVector[13]
  }
}

# append the UTM coordinates to the localSex dataframe, for plotting later on
temp <- merged[, c(1, 5, 6)]
temp1 <- merge(localSex, temp, by="name", all.x=T)
localSex <- temp1
localSex$sex <- as.character(localSex$sex)
# for each individual in localSex, make a new row in "frame" that contains
# its name, pop, spatial scale, local ratio and pop size information
frame <- NA
for (i in 1:nrow(localSex)) {
  name <- as.character(localSex$name[i])
  pop <- as.character(localSex$population[i])
  frame <- rbind(frame, c(name, pop, 0.5, localSex$sex[i], localSex$r05[i], localSex$d05[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 1.0, localSex$sex[i], localSex$r10[i], localSex$d10[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 1.5, localSex$sex[i], localSex$r15[i], localSex$d15[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 2.0, localSex$sex[i], localSex$r20[i], localSex$d20[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 2.5, localSex$sex[i], localSex$r25[i], localSex$d25[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 3.0, localSex$sex[i], localSex$r30[i], localSex$d30[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 3.5, localSex$sex[i], localSex$r35[i], localSex$d35[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 4.0, localSex$sex[i], localSex$r40[i], localSex$d40[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 4.5, localSex$sex[i], localSex$r45[i], localSex$d45[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 5.0, localSex$sex[i], localSex$r50[i], localSex$d50[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 5.5, localSex$sex[i], localSex$r55[i], localSex$d55[i], localSex$Northing[i], localSex$Easting[i]))
  frame <- rbind(frame, c(name, pop, 6.0, localSex$sex[i], localSex$r60[i], localSex$d60[i], localSex$Northing[i], localSex$Easting[i]))
}
# remove rowname and convert to a data frame with correct variable types

rownames(frame) <- NULL
colnames(frame) <- c("name", "population", "level", "sex", "sexRatio", "density", "Northing", "Easting")
frame <- as.data.frame(frame, stringsAsFactors=FALSE)
frame <- frame[!is.na(frame$name),]
frame$population <- as.factor(frame$population)
frame$level <- as.factor(frame$level)
frame$sex <- as.factor(frame$sex)
frame$sexRatio <- as.numeric(frame$sexRatio)
frame$density <- as.numeric(frame$density)
frame$Northing <- as.numeric(frame$Northing)
frame$Easting <- as.numeric(frame$Easting)

varAnalysis <- frame[!is.na(frame$name),]
varAnalysis$numLevel <- as.numeric(as.character(varAnalysis$level))

# This code drops individuals that don't have sex ratios at all neighborhood sizes
names <- varAnalysis$name[varAnalysis$level==0.5 & !is.na(varAnalysis$sexRatio)]
tempFrame <- data.frame(rbind(rep(NA, 9)))
colnames(tempFrame) <- colnames(varAnalysis)
for (i in 1:length(names)) {
  tempFrame <- rbind(tempFrame, varAnalysis[varAnalysis$name==names[i],])
}
tempFrame <- tempFrame[!is.na(tempFrame$name),]
varAnalysis <- tempFrame

varAnalysis$level <- factor(varAnalysis$level)
varAnalysis$sex <- factor(varAnalysis$sex)
varAnalysis$population <- factor(varAnalysis$population)
varAnalysis$name <- factor(varAnalysis$name)
varAnalysis$sqrtSexRatio <- sqrt(varAnalysis$sexRatio)
varAnalysis$transSexRatio <- asin(varAnalysis$sqrtSexRatio)

# Create a new variable to encode sex as 0,1
varAnalysis$sexEncode <- NA
for (i in 1:nrow(varAnalysis)) {
  if (varAnalysis$sex[i] == "hermaphrodite") {
    varAnalysis$sexEncode[i] <- 1
  } else if (varAnalysis$sex[i] == "female") {
    varAnalysis$sexEncode[i] <- 0
  }
}

# A new variable to store the population sex ratio
varAnalysis$popSexRatio <- NA
for (i in 1:length(levels(varAnalysis$population))) {
  varAnalysis$popSexRatio[varAnalysis$population==levels(varAnalysis$population)[i]] <- mean(varAnalysis$sexEncode[varAnalysis$population==levels(varAnalysis$population)[i]])
}
varAnalysis$sexDev <- varAnalysis$sexEncode - varAnalysis$popSexRatio
varAnalysis$socDev <- varAnalysis$sexRatio - varAnalysis$popSexRatio

###################################################################
# Helper functions for calculating standardized selection gradients
# written by Malcolm Augat
# http://people.virginia.edu/~ma5ke/files/software/seln.r
###################################################################

# For each of the given values, subtracts the mean and divides by the standard
# deviation.
standardize = function(x) {
  ndim = max(length(dim(x)), 1)
  if (ndim == 1) return((x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
  margin = 2:ndim
  m = array(data=rep(apply(x, margin, mean, na.rm=TRUE), each=dim(x)[1]),
            dim=dim(x), dimnames=dimnames(x))
  s = array(data=rep(apply(x, margin, sd, na.rm=TRUE), each=dim(x)[1]),
            dim=dim(x), dimnames=dimnames(x))
  return((x - m) / s)
}

# Divides each of the given values by the mean.
relativize = function(x) {
  ndim = max(length(dim(x)), 1)
  if (ndim == 1) return(x / mean(x, na.rm=TRUE))
  margin = 2:ndim
  m = array(data=rep(apply(x, margin, mean, na.rm=TRUE), each=dim(x)[1]),
            dim=dim(x), dimnames=dimnames(x))
  return(x / m)
}