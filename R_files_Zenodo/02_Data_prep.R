### SUBSET AND CLEAN PBDB OCCURRENCE DATA### 

##load saved data files
load("pbdbAnimals.Rdata") #first and last occurrences for taxa in PBDB
data.sep <- read.delim(file="genera_range_through.txt", sep="\\t", header=TRUE) #Sepkoski-based range data with sizes

#create file with genus names and log lengths
myvars <- c("taxon_name", "max_length", "phylum", "class", "order", "family", "calc_max_vol") #extract size and higher taxonomy from range through dataset
y <- data.sep[myvars]
y$logsize <- log10(y$max_length) #transform size to logsize
y$logvol <- log10(y$calc_max_vol) #transform volume to logvolume
y <- subset(y, max_length>0) #remove genera lacking size values
y <- subset(y, calc_max_vol>0) #remove genera lacking size values
names(y) <- c("genus", "max_length", "phylum", "class", "order", "family", "calc_max_vol", "logsize", "logvol")

genera.x <- subset(genera, lad_age<=485.4 & fad_age>1)

##reconcile the higher taxonomy between the two datasets, prioritizing the higher taxonomy from Heim et al 2015
gen1 <- merge(genera.x, y, by="genus")
gen1$phylum.x <- as.character(gen1$phylum.x)
gen1$phylum.y <- as.character(gen1$phylum.y)
gen1$class.x <- as.character(gen1$class.x)
gen1$class.y <- as.character(gen1$class.y)
gen1$order.x <- as.character(gen1$order.x)
gen1$order.y <- as.character(gen1$order.y)
gen1$family.x <- as.character(gen1$family.x)
gen1$family.y <- as.character(gen1$family.y)

for (i in 1:length(gen1$genus)) {
  gen1$phylum[i] <- gen1$phylum.y[i] #apply higher taxonomy from Heim database
  if (gen1$phylum[i]=="") {gen1$phylum[i]==gen1$phylum.x[i]} #if higher taxonomy missing from Heim, apply PBDB values
  
  gen1$class[i] <- gen1$class.y[i] #apply higher taxonomy from Heim database
  if (gen1$class[i]=="") {gen1$class[i]==gen1$class.x[i]} #if higher taxonomy missing from Heim, apply PBDB values
  
  gen1$order[i] <- gen1$order.y[i] #apply higher taxonomy from Heim database
  if (gen1$order[i]=="") {gen1$order[i]==gen1$order.x[i]} #if higher taxonomy missing from Heim, apply PBDB values
  
  gen1$family[i] <- gen1$family.y[i] #apply higher taxonomy from Heim database
  if (gen1$family[i]=="") {gen1$family[i]==gen1$family.x[i]} #if higher taxonomy missing from Heim, apply PBDB values
  
  print(i)
}

##combine bony fish into a single informal class
gen1$class <- as.character(gen1$class)
gen1$class[gen1$class=="Osteichthyes"] <- 'bony fish'
gen1$class[gen1$class=="Actinopteri"] <- 'bony fish'
gen1$class[gen1$class=="Actinopterygii"] <- 'bony fish'

pbdb$class <- as.character(pbdb$class)
pbdb$class[pbdb$class=="Osteichthyes"] <- 'bony fish'
pbdb$class[pbdb$class=="Actinopteri"] <- 'bony fish'
pbdb$class[pbdb$class=="Actinopterygii"] <- 'bony fish'


##further cleaning of the resulting dataset
myvars <- c("genus", "fad_age", "lad_age", "logsize", "logvol", "phylum", "class", "order", "family")  #get reduced list of variable to use for analysis
gen2 <- unique(gen1[,myvars])  #subset dataset to desired variables, unique entries
gen5 <- subset(gen2, class=="Bivalvia" | class=="Cephalopoda" | class=="Gastropoda" |
class=="bony fish" | class=="Ostracoda" | class=="Trilobita" |
class=="Strophomenata" | class=="Rhynchonellata" | class=="Crinoidea" | class=="Echinoidea")

##subset paleobiology database occurrences to intervals of interest
pbdb.1 <- subset(pbdb, max_ma>=1 & min_ma<=485.4)
pbdb.1$genus <- pbdb.1$accepted_name
pbdb.1 <- subset(pbdb.1, select = -c(phylum, order, family, accepted_name))

##standardize the max_ma and min_ma for named stages
library(plyr)
yy.1 <- ddply(pbdb.1, c("early_interval"), function(df)c(max(df$max_ma), min(df$min_ma)))
names(yy.1) <- c("early_interval", "max_ma", "min_ma")
pbdb.2 <- subset(pbdb.1, select = -c(max_ma, min_ma))
pbdb.3 <- merge(pbdb.2, yy.1, by="early_interval")

##find first and last occurrences by genus in cleaned timescale
firstlast <- ddply(pbdb.3, c("genus", "class"), function(df)c(max(df$max_ma), min(df$max_ma)))
names(firstlast) <- c("genus", "class", "fad", "lad")
pbdb.4 <- merge(pbdb.3, firstlast, by=c("genus", "class"))

##merge size and higher taxonomy information with pbdb occurrences
gen6 <- merge(pbdb.4, gen5, by=c("genus", "class"))

##clean up first and last occurrence dates
for (i in 1:length(gen6$genus)) { 
  if (gen6$fad_age[i]>gen6$fad[i]) {gen6$fad[i]<-gen6$fad_age[i]} #if higher taxonomy missing from Heim, apply PBDB values
  if (gen6$lad_age[i]<1) {gen6$lad[i]<-1} #if higher taxonomy missing from Heim, apply PBDB values
  print(i)
}

write.csv(gen6, file="pbdb.cleaned.classes.csv")

gen7 <- unique(subset(gen6, select=c(genus, phylum, class, order, family, fad, lad, logsize, logvol)))
write.csv(gen7, file="genus.sizes.ranges.cleaned.classes.csv")

