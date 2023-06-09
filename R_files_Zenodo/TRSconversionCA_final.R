####
# Title: TRS Batch Conversion
# Author: Katie Pearson
# Date: February 12, 2021
# Purpose: This script cleans the verbatimCoordinates field into a standardized TRS value
# Then it compares the TRS value (with county) to a reference file and imports the
# centroid of the section.
# This version of the script is customized for California TRS data.

### Several critical exceptions apply ###
# The script is not able to address the following circumstances:
# - any information below the rank of section (it scrubs all quarter, half, etc. data before converting to coordinates)
# - when section is listed before the township and range
# - when multiple sections are listed
# - when TRS and UTM data are concatenated into a single field
# - when TRS data don't exist in the reference file (the reference file we used largely excludes sections that are not the standard size)
# - when no section is provided

#load packages
library(DataCombine)
library(gdata)
library(stringr)

#load your dataset that needs coordinates
dat <- read.csv("PATH")

#in my dataset, I first opened the file in Excel and used the function
#=IF(ISERROR(SEARCH("T",[CELL])),"","pickme")
#to flag all the rows that could potentially contain TRS data
#I then sorted the dataset to look at all the cells with "pickme" in that column.
#From there, I manually flagged all data that had TRS data first
#(not at the end of a cell also containing UTM values, due to the difficulty
#of parsing digits from digits). These I flagged as "x" in a TRS column.
#This line will limit my dataset to just those potentially "easy" targets.
fdat <- subset(dat,dat$TRS=="x")

#make a new column to hold the scrubbed data
fdat$trs_only <- fdat$verbatimCoordinates

#make everything uppercase to avoid case redundancy
fdat$trs_only <- toupper(fdat$trs_only)

#scrub excess information from field and standardize
#to township (T), range (R), and section (X)
tofrom <- matrix(ncol = 2, nrow = 76)
tofrom <- as.data.frame(tofrom)
colnames(tofrom) <- c("from","to")
tofrom$from[1] <- "\\\\."
tofrom$to[1] <- ""
tofrom$from[2] <- "SECTION "
tofrom$from[3] <- "SEC "
tofrom$from[4] <- "SECT "
tofrom$from[5] <- "SEC"
tofrom$to[2:5] <- "X"
tofrom$from[6] <- " "
tofrom$from[7] <- "M-CA"
tofrom$from[8] <- "MERIDIAN:MD"
tofrom$from[9] <- ","
tofrom$from[10] <- ";"
tofrom$from[11] <- "NW1/4"
tofrom$from[12] <- "NE1/4"
tofrom$from[13] <- "SW1/4"
tofrom$from[14] <- "SE1/4"
tofrom$from[15] <- "N1/2"
tofrom$from[16] <- "S1/2"
tofrom$from[17] <- "OF"
tofrom$to[6:17] <- ""
tofrom$from[18] <- "TWP"
tofrom$to[18] <- "T"
tofrom$from[19] <- "MT.DIABLOMERIDIAN"
tofrom$from[20] <- "NW/4"
tofrom$from[21] <- "NE/4"
tofrom$from[22] <- "SW/4"
tofrom$from[23] <- "SE/4"
tofrom$from[24] <- "N/2"
tofrom$from[25] <- "S/2"
tofrom$from[26] <- "TRS:"
tofrom$from[27] <- "NW1/16"
tofrom$from[28] <- "NE1/16"
tofrom$from[29] <- "SW1/16"
tofrom$from[30] <- "SE1/16"
tofrom$from[31] <- "NENE"
tofrom$from[32] <- "NENW"
tofrom$from[33] <- "NWNE"
tofrom$from[34] <- "NWNW"
tofrom$from[35] <- "SESE"
tofrom$from[36] <- "SWSW"
tofrom$from[37] <- "SESW"
tofrom$from[38] <- "SWSE"
tofrom$from[39] <- "NESE"
tofrom$from[40] <- "NWSW"
tofrom$from[41] <- "NESW"
tofrom$from[42] <- "NWSE"
tofrom$from[43] <- "SENE"
tofrom$from[44] <- "SWNW"
tofrom$from[45] <- "SENW"
tofrom$from[46] <- "SWNE"
tofrom$from[47] <- "W1/2"
tofrom$from[48] <- "W/2"
tofrom$from[49] <- "E1/2"
tofrom$from[50] <- "E/2"
tofrom$from[51] <- "W1/4"
tofrom$from[52] <- "N1/4"
tofrom$from[53] <- "E1/4"
tofrom$from[54] <- "S1/4"
tofrom$to[19:54] <- ""
tofrom$from[55] <- "R "
tofrom$to[55] <- "R"
tofrom$from[56] <- "T "
tofrom$to[56] <- "T"
tofrom$from[57] <- "S "
tofrom$to[57] <- "S"
tofrom$from[58] <- " N"
tofrom$to[58] <- "N"
tofrom$from[59] <- " S "
tofrom$to[59] <- "S "
tofrom$from[60] <- "SANBERNARDINOMERIDIAN "
tofrom$to[60] <- ""
tofrom$from[61] <- "NR"
tofrom$to[61] <- "N"
tofrom$from[62] <- "SR"
tofrom$to[62] <- "S"
tofrom$from[63] <- "WS"
tofrom$to[63] <- "W"
tofrom$from[64] <- "ES"
tofrom$to[64] <- "E"
tofrom$from[65] <- "NW"
tofrom$from[66] <- "SW"
tofrom$from[67] <- "NE"
tofrom$from[68] <- "SE"
tofrom$from[69] <- "/16"
tofrom$from[70] <- "/4"
tofrom$from[71] <- "/2"
tofrom$from[72] <- "AND"
tofrom$from[73] <- "EEDGE"
tofrom$from[74] <- "SEDGE"
tofrom$from[75] <- "NEDGE"
tofrom$from[76] <- "WEDGE"
tofrom$to[65:76] <- ""
fdat <- FindReplace(fdat, "trs_only", tofrom, from = "from", to = "to", exact = FALSE)

# if the range is listed before the township, switch the order
for(i in 1:dim(fdat)[1]){
  if(startsWith(fdat$trs_only[i],"R")){
    twp <- str_match(fdat$trs_only[i],"T(.*?)[NS]")
    twp <- twp[,1]
    rng <- str_match(fdat$trs_only[i],"R(.*?)[WE]")
    rng <- rng[,1]
    sec <- sub(".*X","",fdat$trs_only[i])
    fdat$trs_only[i] <- paste(twp,rng,sec, sep = "")
  }
}

# remove township (T), range (R), and section (X) delimiters
# also remove leading zeros in front of section and range
tofromfinal <- matrix(ncol = 2, nrow = 7)
tofromfinal <- as.data.frame(tofromfinal)
colnames(tofromfinal) <- c("from","to")
tofromfinal$from[1] <- "T"
tofromfinal$from[2] <- "R"
tofromfinal$from[3] <- "X"
tofromfinal$to[1:3] <- ""
tofromfinal$from[4] <- "N0"
tofromfinal$to[4] <- "N"
tofromfinal$from[5] <- "S0"
tofromfinal$to[5] <- "S"
tofromfinal$from[6] <- "E0"
tofromfinal$to[6] <- "E"
tofromfinal$from[7] <- "W0"
tofromfinal$to[7] <- "W"
fdat <- FindReplace(fdat, "trs_only", tofromfinal, from = "from", to = "to", exact = FALSE)

#remove some other last-minute issues (California-specific)
tofrommer <- matrix(ncol = 2, nrow = 3)
tofrommer <- as.data.frame(tofrommer)
colnames(tofrommer) <- c("from","to")
tofrommer$from[1] <- "MDIABLOMEIDIAN"
tofrommer$from[2] <- "SANBENADINOMEIDIAN"
tofrommer$from[3] <- "HUMBOLDMEIDIAN"
tofrommer$to[1:3] <- ""
fdat <- FindReplace(fdat, "trs_only", tofrommer, from = "from", to = "to", exact = FALSE)

#remove leading zeros
for(i in 1:dim(fdat)[1]){
  if(startsWith(fdat$trs_only[i],"0")){
    fdat$trs_only[i] <- sub(".", "", fdat$trs_only[i])
  }
}

###Part 2: Match TRS values to lat/longs

#load TRS data
#these must include a latitude, longitude, TRS pattern (e.g., 1N13E12), and county field
centroids <- read.csv("PATH")

#make TRS + county columns
centroids$trs_county <- as.character(paste(centroids$pattern,centroids$county, sep = " "))
centroids$trs_county <- toupper(centroids$trs_county)
fdat$trs_county <- as.character(paste(fdat$trs_only,fdat$county, sep = " "))
fdat$trs_county <- toupper(fdat$trs_county)

#remove instances of "county" and "co." from the county field
tofromco <- matrix(ncol = 2, nrow = 2)
tofromco <- as.data.frame(tofromco)
colnames(tofromco) <- c("from","to")
tofromco$from[1] <- " COUNTY"
tofromco$from[2] <- " CO\\\\."
tofromco$to[1:2] <- ""
fdat <- FindReplace(fdat, "trs_county", tofromco, from = "from", to = "to", exact = FALSE)

#add columns for the interpreted (int) data
fdat$decimalLatitude_int <- NA
fdat$decimalLongitude_int <- NA
fdat$coordinateUncertaintyInMeters_int <- NA
fdat$georeferencedBy_int <- NA
fdat$georeferenceRemarks_int <- NA
fdat$georeferenceProtocol_int <- NA
fdat$georeferenceSources_int <- NA
fdat$flag <- NA

#find the coordinate data!
#you will want to fill in the appropriate information for fields listed below
#specifically, georeferencedBy and georeferenceRemarks
for(i in 1:dim(fdat)[1]){
  #find the cell in the centroids file that matches the TRS + county pattern
  m <- match(fdat$trs_county[i],centroids$trs_county)
  #if no match, flag and move on
  if(is.na(m)){
    fdat$flag[i] <- "could not interpret"
    next
  } else {
    #if there is a match, assign the following information to the fields
    fdat$decimalLatitude_int[i] <- as.character(centroids$latitude[m])
    fdat$decimalLongitude_int[i] <- as.character(centroids$longitude[m])
    fdat$coordinateUncertaintyInMeters_int[i] <- "969"
    fdat$georeferencedBy_int[i] <- "kdpearso (via TRS conversion script)"
    fdat$georeferenceRemarks_int[i] <- "auto-generated from TRS by California Phenology Network; February 2021"
    fdat$georeferenceSources_int[i] <- "TRS from label"
    fdat$georeferenceProtocol_int[i] <- "DOI:10.5281/zenodo.4507032"
    #optional progress report
    print(i)
  }
}

#get a count of how many could not be interpreted
not_interpreted <- subset(fdat,fdat$flag=="could not interpret")
dim(not_interpreted)[1]
#calculate a percentage that COULD be interpreted
interpreted <- dim(fdat)[1] - dim(not_interpreted)[1]
filtered <- dim(fdat)[1]
unfiltered <- dim(dat)[1]
print(paste(round((interpreted/filtered)*100, digits = 1),"% successfully interpreted!",sep=""))
print(paste("This is ",round((interpreted/unfiltered)*100, digits = 1),"% of your total starting dataset.",sep=""))

#save this for posterity
write.csv(fdat, "C:/Users/Katie/Desktop/CA_interpreted_coords_v12.csv")

fat <- subset(fdat,is.na(fdat$flag))
thistab <- table(as.character(fat$collid))
thistab <- as.data.frame(thistab)
