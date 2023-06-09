## This script  merges two phytoplankton datasets from Chen & Laws L&O (2017; with some additions) and Kremer et al. L&O (2017) ####

#Read in data:

#Kremer et al. (2017) data
Kdat <- read.csv('Kremer2017.csv')

#Chen and Laws (2017) data with some additions
Cdat <- read.csv('Chen2017.csv')
Cdat <- Cdat[!is.na(Cdat$ID),]
Cdat <- Cdat[!is.na(Cdat$Genus),]
Cdat <- Cdat[!is.na(Cdat$Temp),]
Cdat <- Cdat[!is.na(Cdat$Growth),]

Cdat$Genus <- as.character(Cdat$Genus)
Kdat$genus <- as.character(Kdat$genus)

Cdat$Group                                   <- NA
Cdat$Group[Cdat$Phylum == 'Bacillariophyta'] <- 'Diatom'
Cdat$Group[Cdat$Phylum == 'Dinoflagellata']  <- 'Dino'
Cdat$Group[Cdat$Phylum == 'Cyanobacteria']   <- 'Cyan'
Cdat$Group[Cdat$Phylum == 'Chlorophyta']     <- 'Green'
Cdat$Group[Cdat$Phylum == 'Haptophyta']      <- 'Hapto'
Cdat$Group  <- as.factor(Cdat$Group)
Cdat        <- Cdat[!is.na(Cdat$Group),]

#Create new dataframe
newdat    <- matrix(NA, nr = nrow(Kdat) + nrow(Cdat), nc = 8)
newdat    <- as.data.frame(newdat)
colnames(newdat) <- c('ID','Group','Genus','Species','Habitat','Temperature','Growth','Volume')

for (i in 1:nrow(Cdat)){
    newdat[i, 'ID']          <- paste0('Chen',Cdat[i, 'ID']) 
    newdat[i, 'Group']       <- as.character(Cdat[i,'Group'])
    newdat[i, 'Genus']       <- as.character(Cdat[i,'Genus'])
    newdat[i, 'Species']     <- as.character(Cdat[i,'Species'])
    newdat[i, 'Habitat']     <- as.character(Cdat[i,'Habitat'])
    newdat[i, 'Temperature'] <- Cdat[i,'Temp']
    newdat[i, 'Growth']      <- Cdat[i,'Growth']
    newdat[i, 'Volume']      <- Cdat[i,'Volume']
    newdat[i, 'Lon']         <- Cdat[i,'Lon']
    newdat[i, 'Lat']         <- Cdat[i,'Lat']
}

#Remove NA values in the Kremer dataset
Kdat      <- Kdat[!is.na(Kdat$genus),]
Kdat      <- Kdat[!is.na(Kdat$temperature),]
Kdat      <- Kdat[!is.na(Kdat$r),]

######## Find duplicates between the two datasets
Kdat$frep <- FALSE
Kdat$Kind <- NA
for (i in 1:nrow(Kdat)){
   for (j in 1:nrow(Cdat)){
       if (round(Cdat[j, 'Temp'],0)   == round(Kdat[i, 'temperature'],0) && 
           round(Cdat[j, 'Growth'],2) == round(Kdat[i, 'r'], 2)          &&
                 Cdat[j, 'Genus']     ==       Kdat[i,'genus']){
                  Kdat[i, 'frep'] <- TRUE  #Found repetition
                  Kdat[i, 'Kind'] <- j     #The index for Chen data
                  break
           }
   }
}
######

IDs <- unique(Kdat$isolate.code)
k   <- nrow(newdat[!is.na(newdat$ID),])
for (i in 1:length(IDs)){
    cff  <- Kdat[Kdat$isolate.code == IDs[i], ]
    Ncff <- nrow(cff)
    cff[is.na(cff$frep),'frep'] <- FALSE
    
    #Merge the Kremer dataset into Chen dataset and keep the format same
    if (!any(cff$frep)){
        newdat[(k+1):(k+Ncff), 'ID']          <- paste0('Kremer',cff[, 'isolate.code']) 
        newdat[(k+1):(k+Ncff), 'Group']       <- as.character(cff[,'group'])
        newdat[(k+1):(k+Ncff), 'Genus']       <- as.character(cff[,'genus'])
        newdat[(k+1):(k+Ncff), 'Species']     <- paste(as.character(cff[,'genus']), as.character(cff[,'species']))
        newdat[(k+1):(k+Ncff), 'Habitat']     <- as.character(cff[,'environment'])
        newdat[(k+1):(k+Ncff), 'Temperature'] <- cff[,'temperature']
        newdat[(k+1):(k+Ncff), 'Growth']      <- cff[,'r']
        newdat[(k+1):(k+Ncff), 'Volume']      <- 10**cff[,'bv']
    }
    k <- k + Ncff
}
newdat$ID      <- as.factor(newdat$ID)
newdat         <- newdat[!is.na(newdat$ID),]
newdat$Group   <- as.factor(newdat$Group)
newdat$Habitat <- as.factor(newdat$Habitat)

###Harmonize the names of habitat and functional group
newdat$Habitat[newdat$Habitat == 'freshwater']   <- 'Freshwater'
newdat$Habitat[newdat$Habitat == 'marine']       <- 'Marine'
newdat$Group[newdat$Group   == 'Cyanobacteria']  <- 'Cyan'
newdat$Group[newdat$Group   == 'Diatoms']        <- 'Diatom'
newdat$Group[newdat$Group   == 'Greens']         <- 'Green'
newdat$Group[newdat$Group   == 'Dinoflagellates']<- 'Dino'

save(newdat, file = 'Merged_PHY.Rdata')
