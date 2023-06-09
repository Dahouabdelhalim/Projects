install.packages("maptools")
install.packages("sp")
install.packages("rgeos")
install.packages("adehabitatHR")

require(maptools)

setwd("C:/Workspace/hunt0176/Workspace/R Working Directory/Ss_HomeRange_Adehabitat")

fn <- "C:/Workspace/hunt0176/Workspace/R Working Directory/Ss_HomeRange_Adehabitat/Ss_HR_Adehabitat_setup_UTM.shp"

shp.pts <- readShapePoints(fn, verbose=T)
#MUST MATCH YOUR NUMBER OF RECORDS

xys <- coordinates(shp.pts)
names(shp.pts)

ids <- shp.pts$ID
idsp <- data.frame(ids) 
coordinates(idsp) <- xys
class(idsp)
head(idsp)

#Give a UTM reference system:
proj4string(idsp) <- CRS("+proj=utm +zone=50S +ellps=WGS84") #NB. NWC region is UTM50S
head(idsp)
idsp #will show you the list

library(adehabitatHR)

clu2 <- clusthr(idsp)
class(clu2)
clu2

length(clu2)  #NB. This number has to correspond with the total number of individuals (i.e. 50)

plot(clu2)

kud <- kernelUD(idsp[,1], 
         h = "href", 
         grid = 200, 
         same4all = TRUE, 
         hlim = c(0.1, 1.5), 
         kern = c("bivnorm"), 
         extent = 0.5) 
kud
image(kud)
#THIS WILL JUST PLOT THE HOME RANGES ON A GRID; CHECK THAT MAKES SENSE WITH YOUR INDIVIDUAL'S DISTRIBUTION

#You may need to alter your h-value, so you need to determine what 'href' values for each individual are
#You can query this specifically, where 1,2,21,52 etc are your individual IDs 
#You can do this for all or some IDs, so you get an idea of what h values href is using

kud[[1]]@h
kud[[2]]@h
kud[[21]]@h
kud[[52]]@h

#You then need to convert the object estUD (i.e. kud) in homeranges to shp to open and view in ArcMap
#The purpose of doing this is so you can see all individual homeranges and whether or not
#they are too expanded or overlap land

homerange<- getverticeshr(kud)
writePolyShape (homerange, "utilizationhref")

#If homeranges of IDs seem suitable (in ArcMap), then proveed with calculation of homerange overlap (below)
##############################################################################
#If homeranges of IDs not suitable (in ArcMap), you will need to estimate the UD with a different value of h:
#E.g. h=400

kudh400 <- kernelUD(idsp[,1], 
                h = 400, 
                grid = 200, 
                same4all = TRUE, 
                hlim = c(0.1, 1.5), 
                kern = c("bivnorm"), 
                extent = 0.5) 
homerange<- getverticeshr(kudh400)
writePolyShape (homerange, "utilizationh400")

#Re-check in ArcMap again and if needed, re-do above steps with different h-value until it is suitable.
##############################################################################

#Estimating kernel home range overlap:

kov<-kerneloverlaphr(kud, meth="UDOI",percent=95, conditional=TRUE) #THIS WILL CALCULATE THE OVERLAP

write.table(kov, file = "UDOI.csv", sep = ",", col.names = NA,qmethod = "double")
#THIS IS THE RESULTING FILE and WILL BE STORED IN THE DIRECTORY WHERE R FILES ARE STORED;; 
#then you have to open it with Excel and check the matrix, check that the names and the rest 
#matches exactly as the ID names in your input file in SOCPRoG. 
#For example - input file for R the IDs were 1,2,3,4 etc, but in SOCROG are S001, S002, S003 etc
# So in the resulting matrix you will have to add those letters/prefix numbers to each ID, use the concatenate function in Excel if needed.  

# Next steps:
#1. Upload the Excel csv file to SOCPROG in the GAI module.
#2. Remember to add the short name at the bottom of the upload screen or it will make problems (e.g. UDOI or HR, similar)
#3. Check that the matrix looks fine in the SOCPROG/DOS screen
#4. Add gregariousnes & temporal overlap predictor variables (if haven't already)
#5. Do MRQAP tests to check if variables are significant - if yes, click the GAI button
#6. If variables not significant, remove, re-do MRQAP
#7. Once all variables significant, click GAI button, you will enter next screen, then you can do the rest of your analyses
