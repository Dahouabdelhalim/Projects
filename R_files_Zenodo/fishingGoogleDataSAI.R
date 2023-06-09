suppressMessages(library(robis))
suppressMessages(library(sf))
suppressMessages(library(sqldf))
suppressMessages(library(cowplot))
suppressMessages(library(ggplot2))
suppressMessages(library(DT))
library(rgeos) #modified to eliminate wicket package




if (years!="ALL"){
  cat("selected years:",years,"\\n")
  tmp1<-dataVessel[(dataVessel$years == years) ,]
}else{
  tmp1<-dataVessel
}

if (seasons!="ALL"){
  cat("selected months:",seasons,"\\n")
  tmp1<-tmp1[(tmp1$months == seasons),]
}else{
  
}

write.csv(tmp1,"filtered_data.csv",fileEncoding = "UTF-8",row.names = F,quote = FALSE)
inputFilteredFile<-"filtered_data.csv"
dataTable<<-read.csv(inputFilteredFile,header=T,sep=',')

dataVesselAgg<- aggregate(fahs~Xcentroid+Ycentroid,dataTable,sum)

cat("identification of time limits \\n")
time<-dataTable[,4]

timeMin<-paste0(min(time)-occurrenceTimeRange,"-01-","01")
timeMax<-paste0(max(time)+occurrenceTimeRange,"-12-","31")

cat("analysis is performed with data from",timeMin,"to",timeMax,"\\n")

cat("filter on the bounding box \\n")

#ProjectionExtensionCoord<-wkt_coords(projectionExtension)

ProjectionExtensionCoordT <-readWKT(projectionExtension)
ProjectionExtensionCoord <- as.data.frame(coordinates(ProjectionExtensionCoordT@polygons[[1]]@Polygons[[1]]))

xmin<-min(ProjectionExtensionCoord$x)
xmax<-max(ProjectionExtensionCoord$x)
ymin<-min(ProjectionExtensionCoord$y)
ymax<-max(ProjectionExtensionCoord$y)

cat("filtering...\\n")

# filter the datasets to select only the xcentroid and ycentroid in the selected polygon
dataTableSubset<-dataVesselAgg[dataVesselAgg$Xcentroid>=xmin & dataVesselAgg$Xcentroid<=xmax & 
                                 dataVesselAgg$Ycentroid>=ymin & dataVesselAgg$Ycentroid<=ymax,]

if(dim(dataTableSubset)[1]==0){
  cat("!!! The bounding box is too narrow. Retrieved 0 vessel data. Please use a larger bounding box.")
  stopexecution
}		


geometries<-c()
for (coord in 1:nrow(dataTableSubset)) {
  
  x0<-dataTableSubset$Xcentroid[coord]
  y0<-dataTableSubset$Ycentroid[coord]
  
  x1<-x0-(res/2)
  y1<-y0-(res/2)
  
  x2<-x0-(res/2)
  y2<-y0+(res/2)
  
  x3<-x0+(res/2)
  y3<-y0+(res/2)
  
  x4<-x0+(res/2)
  y4<-y0-(res/2)
  
  geometriesPoly<-paste0("POLYGON((",x1," ",y1,",",x2," ",y2,",",x3," ",y3,",",x4," ",y4,",",x1," ",y1,"))")
  geometries<-c(geometries,geometriesPoly)
}

dataTableSubsetGeom<-cbind(dataTableSubset,geometries)
cat("selecting centroids as function of fishing hours \\n")
dataTableSubsetGeom$norm<-dataTableSubsetGeom$fahs/sum(dataTableSubsetGeom$fahs)

dataTableSubsetGeom<-dataTableSubsetGeom[which(dataTableSubsetGeom$norm>0),]

cat("selecting observations as specified by",sensitivity,"\\n")
dataTableSubsetGeom<-dataTableSubsetGeom[order(dataTableSubsetGeom$norm,decreasing = TRUE),]

if(sensitivity==100){
  dataTableSubsetGeomShort<-dataTableSubsetGeom[1:nrow(dataTableSubsetGeom),] 
}else if (sensitivity==80){
  p<-nrow(dataTableSubsetGeom)*80
  p<-round(p*0.01)
  dataTableSubsetGeomShort<-dataTableSubsetGeom[1:p,] 
} else {
  p<-nrow(dataTableSubsetGeom)*50
  p<-round(p*0.01)
  dataTableSubsetGeomShort<-dataTableSubsetGeom[1:p,] 
}

geometricmean<-exp(mean(log(dataTableSubsetGeomShort$fahs)))
lognormalsd<-(sd(log(dataTableSubsetGeomShort$fahs)))
lowerconf<-exp(mean(log(dataTableSubsetGeomShort$fahs))-lognormalsd)
upperconf<-exp(mean(log(dataTableSubsetGeomShort$fahs))+lognormalsd)#+1.92 lognormalsd

cat("selecting centroids greater than upperconf \\n")
dataTableHigh<-dataTableSubsetGeomShort[which(dataTableSubsetGeomShort$fahs>=upperconf),]

# all other fishing zones
dataTableAll<-dataTableSubsetGeomShort
rm(dataTableSubsetGeomShort)
rm(dataTableSubsetGeom)


## in case there are not high-fishing zones
if(dim(dataTableHigh)[1]==0){
  dataTableHigh<-dataTableAll[which(dataTableAll$norm==max(dataTableAll$norm)),]
  
}else{
  dataTableHigh<-dataTableHigh
}

cat("Number of fishing areas:", nrow(dataTableAll),"\\n")
cat("Number of most exploited fishing areas:", nrow(dataTableHigh),"\\n")


# FAO data
FAOlist<-read.delim(inputFAOList,header=T,sep=",",encoding="UTF-8")
FAOlist<-FAOlist[,c("TAXOCODE","Scientific_name","geometry")]
names(FAOlist)[names(FAOlist) == "Scientific_name"] <- "scientificName"
names(FAOlist)[names(FAOlist) == "geometry"] <- "geometries"

all_fao_obis_species<-NA
for (r in 1:nrow(dataTableHigh)){
  
  geometry<-dataTableHigh$geometries[r]
  fishingHours<-dataTableHigh$fahs[r]
  fishingHoursNorm<-dataTableHigh$norm[r]
  
  cat("retrieve list of species at high fishing areas from OBIS in #",r,"of",nrow(dataTableHigh),"\\n")
  speciesOBIS<-occurrence(geometry=geometry,fields = c("scientificName","class"), startdate =timeMin ,enddate =timeMax )
  
  if(dim(speciesOBIS)[1]==0){
    cat("occurrences retrieved for",r, "of",nrow(dataTableHigh),"are:",nrow(speciesOBIS),"\\n")
    next
  }
  
  speciesOBIS<-speciesOBIS[which(grepl("^[a-z]+ [a-z]+$",tolower(speciesOBIS$scientificName))),]
  names(speciesOBIS)[names(speciesOBIS) == "class"] <- "sp_class"
  
  # for the NA in sp_class
  speciesOBIS$sp_class<-ifelse(is.na(speciesOBIS$sp_class),speciesOBIS$scientificName,speciesOBIS$sp_class)
  
  # cat("Aggregate the Species and Class \\n")
  ## added 04/03
  cat("Aggregate the Species and Class \\n")
  if(dim(speciesOBIS)[1]==0){
    next
  }
  ####
  speciesOBIS$Nocc<-as.numeric("1")
  speciesOBISSum<-aggregate(Nocc~scientificName+sp_class,speciesOBIS,sum)
  
  
  cat("check if among the OBIS species list there are species belonging to the ASFIS List of Species for Fishery Statistics Purposes of FAO \\n")
  speciesOBISSum<-speciesOBISSum[speciesOBISSum$scientificName %in% FAOlist$scientificName,]
  
  if(dim(speciesOBISSum)[1]==0){
    cat("species belonging to FAO list:",nrow(speciesOBISSum),"\\n")
    next
  }else{
    
    cat("dataset OBIS after FAO:",nrow(speciesOBISSum),"\\n")
    speciesOBISSum$geometries<-geometry
    speciesOBISSum$fahs<-fishingHours
    speciesOBISSum$norm<-fishingHoursNorm
    
  }
  
  if (is.na(all_fao_obis_species)[1]){
    all_fao_obis_species<-speciesOBISSum
  }else{
    all_fao_obis_species<-rbind(all_fao_obis_species,speciesOBISSum)
  }
  cat("Number of species belonging to FAO list found at fishing areas:",nrow(all_fao_obis_species)," \\n")
  
}  

if(is.na(all_fao_obis_species) || dim(all_fao_obis_species)[1]==0){
  cat("!!! NO species was reported in OBIS for the selected area and time period.")
  stopexecution
}
cat("...intersect FAO polygons \\n")
#filter FAO SPECIES
FAOlist<-FAOlist[FAOlist$scientificName %in% all_fao_obis_species$scientificName,]
FAOlistSub<-FAOlist[FAOlist$geometries !="",]

if(dim(FAOlistSub)[1]==0){
  cat("!!! Among retrieved OBIS species, there are no FAO fish stocks. Your analysis cannot continue in this area. Please try selecting a larger area.")
  stopexecution
}
FAOlistSub$geometries<-as.character(FAOlistSub$geometries)
FAOlistSub_sf<-st_as_sf(FAOlistSub,wkt="geometries")
#FAOlistSub_sf_valid = (FAOlistSub_sf)#st_make_valid(FAOlistSub_sf)
FAOlistSub_sf_valid = FAOlistSub_sf
FAOlistSub_sf_valid %>% st_cast()

allspecies<-unique(all_fao_obis_species$scientificName)
all_fao_obis_species$is_stock<-F
isthereonestock<-FALSE
for (species in allspecies){
  
  fao_sn_geom<-FAOlistSub_sf_valid[FAOlistSub_sf_valid$scientificName==species,"geometries"][1]
  
  if (dim(fao_sn_geom)[1]==0){
    #cat("Species",species,"is not in FAO stock list\\n")
    all_fao_obis_species[all_fao_obis_species$scientificName==species,"is_stock"]<-F
    next
  }
  
  fao_sn_geomtxt<-st_as_text(fao_sn_geom$geometries)
  fao_sn_geomtxtsc<-st_as_sfc(fao_sn_geomtxt)#st_make_valid(st_as_sfc(fao_sn_geomtxt))
  
  species_geometries<-all_fao_obis_species[all_fao_obis_species$scientificName==species,"geometries"]
  isStock<-FALSE
  
  #faounion <- st_make_valid(fao_sn_geom)#st_make_valid(st_union(fao_sn_geom))#st_make_valid(st_union(fao_sn_geom))#st_make_valid(st_union(fao_sn_geom))
  
  for (sp_geom in species_geometries){
    #intersect
    sp_geom_sf<-st_as_sfc(sp_geom)
    #intersection<-st_intersection(fao_sn_geomtxtsc,sp_geom_sf)
    
    intersection<-st_as_text(st_intersection(sp_geom_sf,st_as_sfc(st_as_text(fao_sn_geomtxtsc %>% st_cast("POLYGON")))))
    
    if (length(intersection)>0)
    {
      isthereonestock<-TRUE
      isStock<-TRUE
      break
    }
  }
  
  cat("Species",species,"has intersection",isStock,"\\n")
  
  #if (species=="Coryphaenoides rupestris")
   # stooop
  
  all_fao_obis_species[all_fao_obis_species$scientificName==species,"is_stock"]<-isStock
}

if(!isthereonestock){
  cat("!!! Among retrieved OBIS species in the selected area, there are no official FAO fish stocks. Your analysis cannot continue in this area. Perhaps the bounding box is too small. Please try selecting a larger area.")
  stopexecution
}

cat("Searching stocks in non-high fishing zones..\\n")

all_fao_obis_species$high_fishing<-TRUE

allStocks<-unique(all_fao_obis_species[all_fao_obis_species$is_stock==T,"scientificName"])
all_fao_obis_species2<-all_fao_obis_species
all_fao_obis_species2$geometries<-as.character(all_fao_obis_species2$geometries)
allpoly<-st_as_text(st_union(st_as_sf(dataTableAll,wkt="geometries")))
highpoly<-st_as_text(st_union(st_as_sf(dataTableHigh,wkt="geometries")))
onePoly<-st_difference(st_as_sfc(allpoly),st_as_sfc(highpoly))
onePoly<-st_as_text(onePoly)

geometry<-onePoly
fishingHours<-sum(dataTableAll$fahs)-sum(dataTableHigh$fahs)
fishingHoursNorm<-NA

for (stock in allStocks){
  errorwithocc<<-F
  occurrenceOBIS<-NA
  tryCatch({
    occurrenceOBIS<-occurrence(scientificname = stock, geometry=geometry,fields = c("scientificName","class"),startdate=timeMin ,enddate= timeMax) 
  }, warning = function(w) {
    errorwithocc<<-T
    cat("WARNING when calling OBIS\\n")
  }, error = function(e) {
    errorwithocc<<-T
    cat("ERROR when calling OBIS\\n")
  }, finally = {
    
  })
  
  nOccurrences = 0
  className<-""
  if (!errorwithocc && !is.na(occurrenceOBIS) && dim(occurrenceOBIS)[1]>0){
    nOccurrences=dim(occurrenceOBIS)[1]
    className<-occurrenceOBIS$class[1]
  }
  
  row<-c(stock,className,nOccurrences,geometry,fishingHours,fishingHoursNorm,T,F)
  all_fao_obis_species2<-rbind(all_fao_obis_species2,row)
}

cat("Searching stocks in non-high fishing zones complete\\n")
cat("Filtering the species occurrences \\n")
all_fao_obis_species2$fahs<-as.numeric(all_fao_obis_species2$fahs)
all_fao_obis_species2$norm<-as.numeric(all_fao_obis_species2$norm)

all_fao_obis_species2$norm<-ifelse(is.na(all_fao_obis_species2$norm),0,all_fao_obis_species2$norm)

all_fao_obis_species2$Nocc<-as.numeric(all_fao_obis_species2$Nocc)

all_fao_obis_species2<-all_fao_obis_species2[which(all_fao_obis_species2$Nocc>minimumNumberOfRecords),]

if(dim(all_fao_obis_species2)[1]==0){
  cat("!!! The number of species occurrences in non-high fishing zones is less than the choosen 'minimum number of records'. Try to run the 
      analysis setting a lower minimum of records or increasing the sensitivity")
  stopexecution
}

numberofNonHighPolygons<-length(dataTableAll$fahs)-length(dataTableHigh$fahs)
numberofHighPolygons<-length(dataTableHigh$fahs)
all_fao_obis_species2$targetScore<-0
all_fao_obis_species2$targetFishery<-"Not impacted"

speciesAggregationForScoring<-aggregate(list(Nocc =all_fao_obis_species2$Nocc,fahs=all_fao_obis_species2$fahs), 
                                        by=list(scientificName=all_fao_obis_species2$scientificName,high_fishing=all_fao_obis_species2$high_fishing),sum)

speciesforscore<-unique(speciesAggregationForScoring$scientificName)
all_fao_obis_species3<-all_fao_obis_species2


for (species in speciesforscore){
  highfish<-speciesAggregationForScoring[which(speciesAggregationForScoring$scientificName==species & speciesAggregationForScoring$high_fishing==T),]
  nonhighfish<-speciesAggregationForScoring[which(speciesAggregationForScoring$scientificName==species & speciesAggregationForScoring$high_fishing==F),]
  
  numberofHighPolygons<-sqldf(paste0("select count(*) as npol from all_fao_obis_species2 where scientificName=='",species,"' AND high_fishing=='TRUE'"),drv="SQLite")
  
  numberofHighPolygons<-numberofHighPolygons$npol
  
  numberofTotalPolygons<-sqldf(paste0("select count(*) as npol, sum(fahs) as fahs from all_fao_obis_species2 where scientificName=='",species,"'"),drv="SQLite")
  
  if (dim(nonhighfish)[1]>0 && dim(highfish)[1]>0)
    targetScore<-((highfish$Nocc*highfish$fahs)/numberofHighPolygons)-((nonhighfish$Nocc*nonhighfish$fahs)/numberofNonHighPolygons)
  
  else if (dim(nonhighfish)[1]>0)
    targetScore<--((nonhighfish$Nocc*nonhighfish$fahs)/numberofNonHighPolygons)
  else
    targetScore<-((highfish$Nocc*highfish$fahs)/numberofHighPolygons)
  
  if (length(targetScore)==0)
    targetScore=0
  
  all_fao_obis_species3[all_fao_obis_species3$scientificName==species,"targetScore"]<-targetScore
  isStock<-all_fao_obis_species3[all_fao_obis_species3$scientificName==species,"is_stock"][1]
  if (isStock){
    targetFishery<-ifelse(targetScore>0,"Trawling","Non-Trawling")
  }else{
    targetFishery<-ifelse(targetScore>0,"By-Catch","None")
  }
  all_fao_obis_species3[all_fao_obis_species3$scientificName==species,"targetFishery"]<-targetFishery
  
  fishingpressure<-numberofTotalPolygons$fahs/numberofTotalPolygons$npol
  all_fao_obis_species3[all_fao_obis_species3$scientificName==species,"fishingpressurescore"]<-fishingpressure
  
}

targetByCatch<-all_fao_obis_species3$targetScore[which(all_fao_obis_species3$is_stock=='FALSE')]#speciesAggregationForTargetting[speciesAggregationForTargetting$targetFishery=='By-Catch',"targetScore"]

geometricmeanBC<-exp(mean(log(targetByCatch)))
lognormalsdBC<-(sd(log(targetByCatch)))
lowerconfBC<-exp(mean(log(targetByCatch))-lognormalsdBC)
upperconfBC<-exp(mean(log(targetByCatch))+lognormalsdBC)
all_fao_obis_species3$fishingpressure<-""
all_fao_obis_species3[which(all_fao_obis_species3$is_stock=='FALSE' & all_fao_obis_species3$targetScore>=upperconfBC),"fishingpressure"]<-"Impacted"

fishing_pressure_scores<-all_fao_obis_species3$fishingpressurescore[which(all_fao_obis_species3$is_stock=='TRUE')]
geometricmeanFP<-exp(mean(log(fishing_pressure_scores)))
lognormalsdFP<-(sd(log(fishing_pressure_scores)))
lowerconfFP<-exp(mean(log(fishing_pressure_scores))-lognormalsdFP)
upperconfFP<-exp(mean(log(fishing_pressure_scores))+lognormalsdFP)

all_fao_obis_species3[which(all_fao_obis_species3$is_stock=='TRUE'),"fishingpressure"]<-"Medium"
all_fao_obis_species3[which(all_fao_obis_species3$is_stock=='TRUE' & all_fao_obis_species3$fishingpressurescore>=upperconfFP),"fishingpressure"]<-"High"
all_fao_obis_species3[which(all_fao_obis_species3$is_stock=='TRUE' & all_fao_obis_species3$fishingpressurescore<=lowerconfFP),"fishingpressure"]<-"Low"


cat("Checking the conservation status of species \\n")

#IUCN<-checklist(scientificname=speciesforscore,redlist=TRUE)
errorinIUCN<<-F
tryCatch({
  IUCN<-checklist(scientificname=speciesforscore,redlist=TRUE)
}, warning = function(w) {
  errorinIUCN<<-T
  cat("WARNING None of the species is included in the 'IUCN Red List species'\\n")
}, error = function(e) {
  errorinIUCN<<-T
  cat("WARNING None of the species is included in the 'IUCN Red List species'\\n")
}, finally = {
  
})

if(errorinIUCN){
  all_fao_obis_species3$Threatened<- "FALSE"
}else{
  cat("Found several threatened species.\\n")
  all_fao_obis_species3$Threatened<- IUCN$scientificName[match(all_fao_obis_species3$scientificName,IUCN$scientificName)]
  all_fao_obis_species3$Threatened<-ifelse(is.na(all_fao_obis_species3$Threatened),"FALSE","TRUE")
}

cat("Computation finished \\n")

##############################################################################
##############################################################################
cat("Preparing for exporting plots and datasets \\n")
##############################################################################
##############################################################################

theme_map <- function(...) {
  theme_minimal() +
    theme(
      axis.title.x = element_text(color ="black",size = 12),#,15),
      axis.title.y = element_text(color ="black",size = 12),#15),
      axis.text=element_text(size=12),
      panel.grid.major = element_line(colour = "grey"),#element_blank(),
      panel.grid.minor = element_line(colour = "grey"),
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.title = element_text(size=12),#15),
      legend.text = element_text(size=11),#12),
      legend.background =element_blank(),
      panel.border = element_blank()
      #plot.caption=element_text(face = "italic",size=12)
    )
}


cat("Map of the most exploited fishing areas considering fishing effort \\n")

xmin<-max(min(dataTableSubset$Xcentroid)-10,-180)
xmax<-min(max(dataTableSubset$Xcentroid)+10,180)
ymin<-max(min(dataTableSubset$Ycentroid)-10,-90)
ymax<-min(max(dataTableSubset$Ycentroid)+10,90)


dataTableAll_sf<-st_as_sf(dataTableAll,wkt="geometries")

for (i in 1:nrow(dataTableAll_sf)){
  
  if(dataTableAll_sf$fahs[i]<=lowerconf){
    dataTableAll_sf$class[i]<-"Low Fishing Effort"
  } else if (dataTableAll_sf$fahs[i]>lowerconf && dataTableAll_sf$fahs[i]<=upperconf){
    dataTableAll_sf$class[i]<-"Medium Fishing Effort"
  } else {
    dataTableAll_sf$class[i]<-"High Fishing Effort"
    
  }
  
}

expl_loc<- ggplot(dataTableAll_sf,aes(Xcentroid,Ycentroid))+ 
  borders(fill="black",colour="black")+
  geom_sf(aes(fill = class))+
  coord_sf(xlim = c(xmin,xmax),ylim=c(ymin,ymax))+
  scale_fill_manual(values=c("Low Fishing Effort"="#ffffd4",
                             "Medium Fishing Effort"="#fe9929",
                             "High Fishing Effort"= "#993404"),
                    name = "")+
  labs(x ="Longitude", y = "Latitude",caption = paste("High, medium and low fishing effort areas"))+
  theme_map()+
  theme(panel.grid.major = element_line(colour = "grey"),plot.caption=element_text(face = "italic",size=10))


cat("Plots of species occurrences, richness and fishing hours across target fishery \\n")

all_fao_obis_species3$fahs<-as.numeric(all_fao_obis_species3$fahs)
all_fao_obis_species3$Nocc<-as.numeric(all_fao_obis_species3$Nocc)
all_fao_obis_species3$norm<-as.numeric(all_fao_obis_species3$norm)

SppFA<-all_fao_obis_species3
SppFA$fishingpressure<-ifelse(SppFA$fishingpressure=="","Other",SppFA$fishingpressure)
SppFA<-aggregate(list(fahs =SppFA$fahs,Nocc=SppFA$Nocc), 
                 by=list(scientificName=SppFA$scientificName,
                         fishingpressure=SppFA$fishingpressure, Threatened=SppFA$Threatened),sum)

SppFA$fishingpressure<-ifelse(SppFA$fishingpressure=="High","HPrs.Spp.",
                              ifelse( SppFA$fishingpressure=="Medium","MPrs.Spp.",
                                      ifelse( SppFA$fishingpressure=="Low","LPrs.Spp.",
                                              ifelse( SppFA$fishingpressure=="Other","Non-stock Spp.",
                                                      ifelse( SppFA$fishingpressure=="Impacted","Impacted Spp.",SppFA$fishingpressure)))))




HA_occ<-ggplot(SppFA, aes(x=fishingpressure,y=Nocc, fill=Threatened))+ 
  geom_bar(stat="identity",width = 0.3,position = position_dodge())+
  scale_fill_manual(values = c("#999999", "#f03b20"),labels=c("Not Threatened","Threatened"))+
  #scale_x_discrete(labels=c("HPrs.Spp.", "MPrs.Spp.", "Non-stock Spp.", "Bycatch Spp.")) +
  theme_map()+
  theme(axis.line = element_line(colour = "black"),panel.grid.minor =element_blank(),
        panel.grid.major=element_blank(),
        axis.text.x = element_text(), legend.title = element_blank(),
        axis.title.x=element_blank(),legend.position = "none")+
  labs(y = "Number of occurrences")


HA_fahs<-ggplot(SppFA, aes(x=fishingpressure,y=fahs, fill=Threatened))+ 
  geom_bar(stat="identity",width = 0.3,position = position_dodge())+
  scale_fill_manual(values = c("#999999", "#f03b20"),labels=c("Not Threatened","Threatened"))+
  # scale_x_discrete(labels=c("HPrs.Spp.", "MPrs.Spp.", "Non-stock Spp.", "Bycatch Spp.")) +
  theme_map()+
  theme(axis.line = element_line(colour = "black"),panel.grid.minor =element_blank(),
        panel.grid.major=element_blank(),
        axis.text.x = element_text(),legend.title = element_blank(),
        axis.title.x=element_blank(),legend.position = "none")+
  labs(y = "Fishing hours")

SppFA$id<-as.numeric("1")
SppFA2<-aggregate(list(id =SppFA$id), 
                  by=list(Threatened=SppFA$Threatened, fishingpressure=SppFA$fishingpressure),sum)


HA_rich<-ggplot(SppFA2, aes(x=fishingpressure,y=id, fill=Threatened))+ 
  geom_bar(stat="identity",width = 0.3,position = position_dodge())+
  scale_fill_manual(values = c("#999999", "#f03b20"), labels=c("Not Threatened","Threatened"))+
  #scale_x_discrete(labels=c("HPrs.Spp.", "MPrs.Spp.", "Non-stock Spp.", "Bycatch Spp.")) +
  theme_map()+
  theme(axis.line = element_line(colour = "black"),panel.grid.minor =element_blank(),
        panel.grid.major=element_blank(),
        axis.text.x = element_text(),legend.title = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "bottom", plot.caption=element_text(hjust = 0,face = "italic"))+
  labs(y = "Number of species",caption = "Fishing pressures on IUCN threatened and not threatened species.Abbreviations indicate: High Pressure Species (HPrs.Spp.), Medium Pressure Species (MPrs.Spp.),
       Low Pressure Species (LPrs.Spp), Non-stock Species (Non-stock Spp.), and Possible-Impacted Species (Impacted Spp.)")

BarPlots<-plot_grid(HA_fahs,HA_occ, HA_rich,labels = "AUTO",ncol = 1)


# round only numeric columns
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}


# aggregate by occurrences and fishing hours
summaryTable<-aggregate(list(fahs =all_fao_obis_species3$fahs,Nocc=all_fao_obis_species3$Nocc), 
                        by=list(scientificName=all_fao_obis_species3$scientificName,
                                className=all_fao_obis_species3$sp_class,
                                targetFishery=all_fao_obis_species3$targetFishery, 
                                Stock=all_fao_obis_species3$is_stock,
                                FishingPressures=all_fao_obis_species3$fishingpressure, 
                                Threatened=all_fao_obis_species3$Threatened,
                                fishingpressurescore=all_fao_obis_species3$fishingpressurescore),sum)

summaryTable<-round_df(summaryTable, digits = 3)

summaryTable <- summaryTable[order(summaryTable$Nocc,decreasing = TRUE),]

summaryTable$Threatened<-ifelse(summaryTable$Threatened=="TRUE","Y","N")
summaryTable$Stock<-ifelse(summaryTable$Stock=="TRUE","Y","N")

summaryTable<-summaryTable[order(summaryTable$fishingpressurescore,decreasing = T),]
summary_html <- DT::datatable(summaryTable[,c(2,1,9,8,4,6,5)],rownames = FALSE,
                              colnames = c('Class','Scientific name', 'N. of Occurrences', 'Fishing hours', 'Is a Stock', 'Is Threatened', 'Estimated Fishing Impact'),
                              caption = htmltools::tags$caption(
                                style = 'caption-side: bottom; text-align: left;',
                                htmltools::em('Summary table of fishing activity and impacts per species')),
                              options = list(
                                initComplete = JS(
                                  "function(settings, json) {",
                                  "$('body').css({'font-family': 'Calibri'});",
                                  "}"
                                ),
                                columnDefs = list(list(className = 'dt-center', targets = "_all")),
                                pageLength = 20,lengthChange = FALSE# for the show entries
                              )
)







