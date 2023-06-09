#this code is meant to identify the species not previously identified in the parks where there were less than 50 (or other specified number) in the park
setwd("C:/Users/amkat/Desktop/VSFS/compare_iNat_NPS")

library(stringr) #use str_extract from this package

#table3 <- read.csv("./ACAD/ACAD_table3.csv", stringsAsFactors = FALSE)

parkname <- "ZION"
table3 <- read.csv(paste("./", parkname, "/", parkname, "_table3.csv", sep=""), stringsAsFactors = FALSE)
table3 <- table3[table3$Occurences.in.BISON != 'No data available',]

park_states <- read.csv("./20200928_migration/nps_boundary.csv", stringsAsFactors = FALSE)

BISON_in_state <- rep(NA, length(table3$Scientific_name))
indivdparkstate <- rep(NA, length(table3$Scientific_name))
speciesundercutoff <- data.frame(table3$Scientific_name, table3$Invasive.Record, BISON_in_state, 
                                 indivdparkstate, table3$Occurences.in.BISON, stringsAsFactors = FALSE)

cutoff <- 50 #put in the cutoff amount that you want for entries in the state

for (i in 1:length(table3$Scientific_name)){
  numberinBISON <- unlist(strsplit(table3$Occurences.in.BISON[i], split=";")) #splits up the BISON occurrences from table 3
  numberinBISON <- gsub("[[:space:]]", "", numberinBISON) #remove white spaces
  parkstate <- park_states$STATE[park_states$UNIT_CODE == parkname] #get park state from the database
  
  if (length(grep(parkstate, numberinBISON)) > 0){
  statedata <- as.integer(str_extract(numberinBISON[pmatch(parkstate, numberinBISON)], "[0-9]+$")) 
   if (statedata <= cutoff) {
      speciesundercutoff$BISON_in_state[i] <- statedata
      speciesundercutoff$indivdparkstate[i] <- parkstate
   }
  }
  
  if (length(grep(parkstate, numberinBISON)) <= 0){
    speciesundercutoff$BISON_in_state[i] <- 'No records'
    speciesundercutoff$indivdparkstate[i] <- parkstate
  }
}

speciesundercutoff <- na.omit(speciesundercutoff[speciesundercutoff$BISON_in_state != "NA",])
speciesundercutoff

write.csv(speciesundercutoff, paste("./20200928_migration/", parkname, "_migration_under_50.csv", sep=""))

rm(BISON_in_state, indivdparkstate, speciesundercutoff)
