###CREATE CAPTURE HISTORIES FOR EACH GENERA AND ASSIGN MASS EXTINCTION AND RECOVERY INTERVALS###

#read in occurrences and sizes
genera.1 <- read.csv("pbdb.cleaned.classes.csv")  
size.data <- read.csv("genus.sizes.ranges.cleaned.classes.csv")

#delimit the time intervals of interest
start.age <- c(485.4) #start time, in millions of years ago, of the study interval
end.age <- c(1) #end time, in millions of years ago, of the study interval

nd <- list()

#subset Ordovician through Pleistocene occurrences
genera <- subset(genera.1, max_ma<=start.age & min_ma>=end.age)
x <- table(unlist(genera$genus, genera$class, genera$logvol), genera$max_ma) #table genus occurrences by interval
x <- x[, match(colnames(x), rev(colnames(x)))]  # reverse the column order so time is going forward
x[x > 0] <- 1 # set counts > 1 equal to 1; presence-absence matrix

#matrix to hold binary extinction values, with 0 being background and 1 being mass extinction
mass.extinction <- matrix(rep(rep(0, ncol(x)), each=nrow(x)),ncol=ncol(x)) 
colnames(mass.extinction) <- colnames(x) 
mass.extinction[,"443.4"] <- 1  #end-Ordovician mass extinction
mass.extinction[,"372.2"] <- 1  #Late Devonian mass extinction
mass.extinction[,"252.17"] <- 1 #end-Permian mass extinction
mass.extinction[,"201.3"] <- 1 #end-Triassic mass extinction
colnames(mass.extinction) <- paste("mass.ext",1:ncol(x),sep="")

#matrix to hold binary origination values, with 0 being background and 1 being post-mass extinction
recovery.mass <- matrix(rep(rep(0, ncol(x)), each=nrow(x)),ncol=ncol(x)) 
colnames(recovery.mass) <- colnames(x) 
recovery.mass[,"440.8"] <- 1  #end-Ordovician recovery
recovery.mass[,"358.9"] <- 1  #Late Devonian recovery
recovery.mass[,"251.2"] <- 1 #end-Permian recovery
recovery.mass[,"199.3"] <- 1 #end-Triassic recovery
colnames(recovery.mass) <- paste("rec.mass",1:ncol(x),sep="")

#create a column called 'ch' with encounter history as a string O's and 1's'.
ch <- apply(x, 1, paste, collapse="") 
encounterHistory <- data.frame(genus=rownames(x), ch, nSight=rowSums(x), stringsAsFactors=FALSE)
encounterHistory <- cbind(encounterHistory, mass.extinction, recovery.mass)
encounterHistory <- encounterHistory[encounterHistory$nSight > 0,] #remove any taxon with no occurrences

#create a list of genera by class
myvars <- c("genus", "class", "phylum", "order", "logvol") 
class.list <- unique(genera[myvars]) #create a list of genera by class

#bring class assignments back into encounterHistory
encounterHistory <- merge(encounterHistory, class.list, by="genus") 

# Make sure to install and load the "dplyr" package before running this code!
library(plyr)
class.median <- ddply(encounterHistory, c("class"), function(df)c(median(df$logvol), mean(df$logvol), length(df$logvol))) #find median and mean size by class
names(class.median) <- c("class", "median", "mean", "n") #name output dataset
class.median <- class.median[class.median$n>1,] #restrict data to classes with enough genera
encounterHistory <- merge(encounterHistory, class.median, by="class") #merge median size information back into encounterHistory
encounterHistory$small <- 0 #create binary variable for larger or smaller than median
encounterHistory$small[encounterHistory$logvol>encounterHistory$median] <- 1 #assign new value to genera bigger than median

myvars <- c("ch", "class", "phylum", "logvol", colnames(mass.extinction), colnames(recovery.mass))
nd4 <- encounterHistory[myvars] #subset occs data frame to relevant columns

#All genera with logvol data
save(nd4, file="Logvol_encounter_history_First4.Rdata") 