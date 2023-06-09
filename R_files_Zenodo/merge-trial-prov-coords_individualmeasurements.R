#This script merges the individual phenotypic data with the geographical coordinates of the corresponding trials and provenances. 

#set working directory 
setwd("/WorkingDirectory/")


#Read databases
Fsyl <- read.csv("Fsylvatica.csv")
Trials <- read.csv("Trial_coords.csv")
Prov <- read.csv("Prov_coords.csv")

#Rename coordinates 
colnames(Trials)[14] <- "lonTrial"
colnames(Trials)[15] <- "latTrial"
colnames(Prov)[13] <- "lonProv"
colnames(Prov)[14] <- "latProv"

#Add tree age for each measure
Fsyl$H95cmAge <- (1995 - Fsyl$YearPlantation)
Fsyl$H96cmAge <- (1996 - Fsyl$YearPlantation)
Fsyl$H97cmAge <- (1997 - Fsyl$YearPlantation)
Fsyl$H98cmAge <- (1998 - Fsyl$YearPlantation)
Fsyl$H99cmAge <- (1999 - Fsyl$YearPlantation)
Fsyl$H00cmAge <- (2000 - Fsyl$YearPlantation)
Fsyl$H01cmAge <- (2001 - Fsyl$YearPlantation)
Fsyl$H02cmAge <- (2002 - Fsyl$YearPlantation)
Fsyl$H03cmAge <- (2003 - Fsyl$YearPlantation)
Fsyl$H04cmAge <- (2004 - Fsyl$YearPlantation)
Fsyl$H05cmAge <- (2005 - Fsyl$YearPlantation)
Fsyl$H06cmAge <- (2006 - Fsyl$YearPlantation)
Fsyl$H07cmAge <- (2007 - Fsyl$YearPlantation)
Fsyl$H08cmAge <- (2008 - Fsyl$YearPlantation)


#Merge geographical coordinates of trials and provenances with phenotypic data
new = merge(Fsyl,Trials, "Trial")
data = merge(new, Prov, by="ID_ProvCode")


#save your data with trial and provenance coordinates: 
write.csv(data, file = "Fsylvatica_coord.csv", sep="\\t", col.names=TRUE)


