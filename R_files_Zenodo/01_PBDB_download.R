### READ IN ALL OCCURRENCES FROM PALEOBIOLOGY DATABASE AT GENUS LEVEL ###

# function to eliminate occurrences that are not constrained to geological age
# modified from velociraptr package 
constrainAges <- function (Data, Timescale) 
{
  Data[, "early_interval"] <- as.character(Data[, "early_interval"])
  Data[, "late_interval"] <- as.character(Data[, "late_interval"])
  for (i in 1:nrow(Timescale)) {
    EarlyPos <- which(Data[, "max_ma"] > Timescale[i, "min_ma"] & Data[, "max_ma"] <= Timescale[i, "max_ma"])
    Data[EarlyPos, "early_interval"] <- as.character(Timescale[i, "interval_name"])
    LatePos <- which(Data[, "min_ma"] >= Timescale[i, "min_ma"] &  Data[, "min_ma"] < Timescale[i, "max_ma"])
    Data[LatePos, "late_interval"] <- as.character(Timescale[i,"interval_name"])
  }
  Data <- Data[Data[, "early_interval"] == Data[, "late_interval"], ]
  return(Data)
}

pbdb <- data.frame()
maxAge <- read.delim("https://paleobiodb.org/data1.2/intervals/single.tsv?name=Cambrian")
geoIntervals <- read.delim(paste("https://paleobiodb.org/data1.2/intervals/list.tsv?scale=1&scale_level=5&max_ma=",maxAge$max_ma[1], sep=""))

# subgenera are elevated to genus
for(i in 1:nrow(geoIntervals)) {
  pbdbUri <- URLencode(paste("https://paleobiodb.org/data1.2/occs/list.tsv?all_records&base_name=Animalia&interval=",geoIntervals$interval_name[i],"&show=paleoloc,classext&idreso=lump_gensub&limit=all", sep=""))
  temp <- read.delim(pbdbUri)
  pbdb <- rbind(pbdb, temp)
  print(as.character(geoIntervals$interval_name[i])) ## just to keep track in the console
}
pbdb <- unique(pbdb) # eliminate duplicates from poorly resolved genera
pbdb <- subset(pbdb, !is.na(genus)) # make sure there are not taxa without genus names
pbdb <- constrainAges(pbdb, geoIntervals) # eliminate occurrences not constrained to a single age
pbdb$genus_id <- pbdb$genus_no
pbdb$genus_id[!is.na(pbdb$subgenus_no)] <- pbdb$subgenus_no[!is.na(pbdb$subgenus_no)] # use subgenus number for unique genus id

# create data frame of just genera and their fads and lads
genera <- unique(pbdb[,match(c('genus_id','phylum','class','order', 'family', 'genus'),colnames(pbdb))]) # unique genera
fads <- tapply(pbdb$max_ma, pbdb$genus_id, max) # calculate FADs
lads <- tapply(pbdb$min_ma, pbdb$genus_id, min) # calculate LADs
genera$fad_age <- fads[match(genera$genus_id, names(fads))]
genera$lad_age <- lads[match(genera$genus_id, names(lads))]

save(pbdb, genera, geoIntervals, file="pbdbAnimals.Rdata")

#read in data files from Paleosize Database#
data.pbdb <- pbdb #PBDB data

write.csv(data.pbdb, file="pbdb.data.csv")
