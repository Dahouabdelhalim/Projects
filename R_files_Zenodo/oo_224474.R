library(jsonlite)
library(tidyverse)

#get a list of barcodes by which to retrieve data from the API from the JPEG files in the dataset
setwd("IMG LOCATION")
files=list.files(pattern="*.jpg",recursive=T,include.dirs=F)
splitnames=strsplit(files,pattern="/")
filenames=sapply(seq(1,length(files)),function(x) splitnames[[x]][length(splitnames[[x]])])
filenames=gsub(".jpg","",filenames)

#Request the data from the API. You'll need an access token, which you can get through the API (https://api.laji.fi) 
#under APIUser->POST
rawdata=list(1)
for (i in 2:201) {
  rawdata[[i]]=fromJSON(
    paste0("https://api.laji.fi/v0/warehouse/query/single?documentId=http%3A%2F%2Fid.luomus.fi%2F",
           filenames[i-1],
           "&access_token=FILL_IN_ACCESS_TOKEN"))
}

#Convert the massive data list to a data frame
#given the different number of variables for each specimen, we're just going to full join them all
tdata=unlist(rawdata[[2]])
td=data.frame(t(tdata))
colnames(td)=names(tdata)
fulltd=td
for (i in 3:201) {
  tdata=unlist(rawdata[[i]])
  td=data.frame(t(tdata))
  colnames(td)=names(tdata)
  fulltd=full_join(fulltd,td)
}

#Filter out interesting properties that are named
lesstd=select(fulltd,
              "document.collectionId",
              "document.documentId",
              "document.gatherings.team",
              "document.gatherings.units.taxonVerbatim",
              "document.gatherings.units.author",
              "document.gatherings.units.det",
              "document.gatherings.units.recordBasis",
              "document.gatherings.units.unitId",
              "document.gatherings.units.typeSpecimen",
              "document.gatherings.units.media.licenseId",
              "document.gatherings.units.media.copyrightOwner",
              "document.gatherings.eventDate.begin",
              "document.gatherings.eventDate.end",
              "document.gatherings.conversions.year",
              "document.gatherings.locality",
              "document.gatherings.country",
              "document.createdDate",
              "document.modifiedDate",
              "document.linkings.editors.fullName",
              "document.notes",
              "document.gatherings.conversions.month",
              "document.gatherings.conversions.day",
              "document.gatherings.higherGeography",
              "document.gatherings.province",
              "document.gatherings.coordinatesVerbatim",
              "document.gatherings.interpretations.coordinateAccuracy",
              "document.gatherings.interpretations.coordinates.latMin",
              "document.gatherings.interpretations.coordinates.latMax",
              "document.gatherings.interpretations.coordinates.lonMin",
              "document.gatherings.interpretations.coordinates.lonMax",
              "document.gatherings.interpretations.coordinates.type",
              "document.gatherings.municipality")

#Filter out potentially interesting properties ("values") that are named by other variables ("facts"). 
#jam will contain the property names, vjam their values

jam=fulltd[,grepl("document.gatherings.units.facts.fact*",colnames(fulltd))==T|
             grepl("document.gatherings.facts.fact*",colnames(fulltd))==T|
             grepl("document.facts.fact*",colnames(fulltd))==T]
vjam=fulltd[,grepl("document.gatherings.units.facts.value*",colnames(fulltd))==T|
              grepl("document.gatherings.facts.value*",colnames(fulltd))==T|
              grepl("document.facts.value*",colnames(fulltd))==T]

#Remove URL's and other identifiers, from both values and variable names
#make a list of all possible properties, to make a new data frame ljam which has one column 
#for each such property
jam2=as.data.frame(sapply(seq(1:length(jam)),function(x) gsub("http://tun.fi/","",jam[,x])))
vjam2=as.data.frame(sapply(seq(1:length(vjam)),function(x) gsub("http://tun.fi/","",vjam[,x])))
jamterms=unlist(jam2)
jamterms=unique(jamterms)
jamterms2=gsub("MY.","",jamterms)
jamterms2=gsub("MZ.","",jamterms2)
jamterms2=gsub("MY .","",jamterms2)
jamterms2=jamterms2[is.na(jamterms2)==F]
ljam=data.frame(matrix(nrow=length(jam),ncol=length(jamterms2)))
colnames(ljam)=jamterms2
jam2=as.data.frame(sapply(seq(1:length(jam)),function(x) gsub("MZ.","",gsub("MY .","",gsub("MY.","",jam2[,x])))))

#fill in this data frame with the proper values, then combine with the properties extracted before
for (i in 1:200) {
  for (j in 1:length(jam)) {
    if (is.na(jam2[i,j])==F) {
      ljam[i,match(jam2[i,j],colnames(ljam))]=as.character(vjam2[i,j])
    }
  }
}
newtd=cbind(lesstd,ljam)

#drop some properties we won't use
newtd.s=select(newtd,
               -document.gatherings.units.media.licenseId,
               -document.createdDate,
               -document.linkings.editors.fullName,
               -microscopy,
               -speciesQualifier,
               -georeferenceSource,
               -inMustikka,
               -editor,
               -owner,
               -documentLocation,
               -dateCreated,
               -entered,
               -additionalIDs,
               -creator,
               -status,
               -acquiredFrom,
               -editNotes,
               -exsiccatum,
               -duplicatesIn,
               -unreliableFields,
               -typifDate,
               -typif,
               -typeSpecies,
               -typeBasionymePubl,
               -typePubl)

#Basic remappings to pseudo-DWC
newtd.m=rename(newtd.s,
               "collectionID1"="document.collectionId",
               "occurrenceID"="document.documentId",
               "recordedBy"="document.gatherings.team",
               "canonicalName"="document.gatherings.units.taxonVerbatim",
               "scientificNameAuthorship"="document.gatherings.units.author",
               "identifiedBy"="document.gatherings.units.det",
               "basisOfRecord"="document.gatherings.units.recordBasis",
               "identifier"="document.gatherings.units.unitId",
               "isType"="document.gatherings.units.typeSpecimen",
               "rightsHolder"="document.gatherings.units.media.copyrightOwner",
               "eventDate1"="document.gatherings.eventDate.begin",
               "eventDate2"="document.gatherings.eventDate.end",
               "year"="document.gatherings.conversions.year",
               "locality1"="document.gatherings.locality",
               "country1"="document.gatherings.country",
               "modified"="document.modifiedDate",
               "organismRemarks"="document.notes",
               "month"="document.gatherings.conversions.month",
               "day"="document.gatherings.conversions.day",
               "continent"="document.gatherings.higherGeography",
               "stateProvince"="document.gatherings.province",
               "verbatimCoordinates"="document.gatherings.coordinatesVerbatim",
               "coordinateUncertaintyInMeters"="document.gatherings.interpretations.coordinateAccuracy",
               "decimalLatitudeMin"="document.gatherings.interpretations.coordinates.latMin",
               "decimalLatitudeMax"="document.gatherings.interpretations.coordinates.latMax",
               "decimalLongitudeMin"="document.gatherings.interpretations.coordinates.lonMin",
               "decimalLongitudeMax"="document.gatherings.interpretations.coordinates.lonMax",
               "geodeticDatum"="document.gatherings.interpretations.coordinates.type",
               "county1"="document.gatherings.municipality",
               "basisOfRecord2"="recordBasis",
               "canonicalName2"="taxonVerbatim",
               "establishmentMeans"="provenance",
               "relationshipRemarks"="macroscopy",
               "measurementRemarks"="chemistry",
               "canonicalName3"="taxon",
               "identifiedBy2"="det",
               "dateIdentified"="detDate",
               "identificationRemarks"="preferredIdentification",
               "taxonRank2"="infraRank",
               "identificationRemarks2"="identificationNotes",
               "scientificNameAuthorship2"="infraAuthor",
               "verbatimIdentifiedBy"="detVerbatim",
               "verbatimEventDate"="dateVerbatim",
               "verbatimLocality"="localityVerbatim",
               "eventDate4"="dateEnd",
               "verbatimRecordedBy"="legVerbatim",
               "continent2"="higherGeography",
               "stateProvince2"="administrativeProvince",
               "recordedBy2"="leg",
               "county2"="municipality",
               "eventDate3"="dateBegin",
               "country2"="country",
               "verbatimLongitude"="longitude",
               "verbatimCoordinates2"="coordinatesVerbatim",
               "locality2"="locality",
               "verbatimLocality2"="localityDescription",
               "verbatimLatitude"="latitude",
               "fieldNotes"="notes",
               "county3"="county",
               "habitat"="habitatDescription",
               "decimalLatitude"="wgs84Latitude",
               "verbatimCoordinateSystem2"="coordinateSystem",
               "decimalLongitude"="wgs84Longitude",
               "habitat2"="biologicalProvince",
               "recordedBy3"="collectingEventName",
               "verbatimElevation"="alt",
               "coordinateUncertaintyInMeters2"="coordinateRadius",
               "georeferenceRemarks"="coordinateNotes",
               "collectionID2"="collectionID",
               "type2"="datatype",
               "preparations"="preservation",
               "modified2"="dateEdited",
               "otherCatalogNumbers"="originalSpecimenID",
               "recordNumber"="legID",
               "infraSpecificEpithet"="infraEpithet",
               "scientificNameAuthorship3"="author")

##Fix duplicated properties
#least informative are removed, sometimes remapping to other DwC properties is needed
#!!!!! The following steps were performed based on the data as it was available at the time. They might be
#less optimal if the data in the platform is updated or otherwise changed. Some steps may even produce errors
#if new use cases of the data fields come up !!!!!!

#collectionID
newtd.m$collectionID=newtd.m$collectionID1
hd=select(newtd.m,-collectionID1,-collectionID2)

#recordedBy
hd$recordedBy=hd$recordedBy2
hd$eventRemarks=hd$recordedBy3
hd=select(hd,-recordedBy2,-recordedBy3)

#scientificName: we combine the different ones to derive a canonicalName and a full scientificName 
#with as much information as available (and apparently in hideous code)
hd=select(hd,-scientificNameAuthorship3,-canonicalName3)
hd$canonicalNameF=hd$canonicalName

#if name2 is a full scientific name and not identical to name1, it was a longer name
hd$canonicalNameF[hd$canonicalName!=hd$canonicalName2&
                    is.na(hd$canonicalName2)==F&
                    toupper(substr(hd$canonicalName2,1,1))==substr(hd$canonicalName2,1,1)]=
  hd$canonicalName2[hd$canonicalName!=hd$canonicalName2&
                      is.na(hd$canonicalName2)==F&
                      toupper(substr(hd$canonicalName2,1,1))==substr(hd$canonicalName2,1,1)]

#if name2 is not a full name, it provides some extra info to past onto it
hd$canonicalNameF[is.na(hd$canonicalName2)==F&
                    toupper(substr(hd$canonicalName2,1,1))!=substr(hd$canonicalName2,1,1)]=
  paste(hd$canonicalNameF[is.na(hd$canonicalName2)==F&
                            toupper(substr(hd$canonicalName2,1,1))!=substr(hd$canonicalName2,1,1)],
        hd$canonicalName2[is.na(hd$canonicalName2)==F&
                            toupper(substr(hd$canonicalName2,1,1))!=substr(hd$canonicalName2,1,1)])

#also an ad hoc fix:
hd$canonicalNameF=gsub("Häkkinen","",hd$canonicalNameF)
hd$canonicalNameF=gsub("Räsänen","",hd$canonicalNameF)

#Construct scientificName based on how many words present
sciname=strsplit(hd$canonicalNameF,split=" ")
hd$scientificName=NA
for (i in 1:200) {
  if (length(sciname[[i]])==1) {
    hd$scientificName[i]=paste(sciname[[i]][1],hd$scientificNameAuthorship[i])  
  }
  if (length(sciname[[i]])==2) {
    hd$scientificName[i]=paste(sciname[[i]][1],sciname[[i]][2],hd$scientificNameAuthorship[i])  
  }
  if (length(sciname[[i]])==4) {
    hd$scientificName[i]=paste(sciname[[i]][1],sciname[[i]][2],hd$scientificNameAuthorship[i],
                               sciname[[i]][3],sciname[[i]][4],hd$scientificNameAuthorship2[i])  
  }
}

#taxonRank
hd$scientificName=gsub("NA","",hd$scientificName)
hd$taxonRank[is.na(hd$taxonRank2)==F]=hd$taxonRank2[is.na(hd$taxonRank2)==F]
hd$taxonRank=gsub("MX.","",hd$taxonRank)
hd$taxonRank=gsub("MY.","",hd$taxonRank)
hd$taxonRank=gsub("infraRank","",hd$taxonRank)
hd$canonicalName=hd$canonicalNameF
hd=select(hd,-canonicalName2,-canonicalNameF,-taxonRank2,-scientificNameAuthorship2)

#genus, specificEpithet
hd$genus=sapply(seq(1,200),function(x) sciname[[x]][1])
hd$specificEpithet=sapply(seq(1,200),function(x) sciname[[x]][2])

#identification
hd$identificationRemarks=hd$identificationRemarks2
hd=select(hd,-identifiedBy2,-verbatimIdentifiedBy,-identificationRemarks2)

#basisOfRecord
hd=select(hd,-basisOfRecord2,-type2,-preparations)

#dates
hd$modified=hd$modified2
hd=select(hd,-modified2,-eventDate3,-eventDate4)
hd$eventDate=ifelse(hd$eventDate1!=hd$eventDate2&is.na(hd$eventDate1)==F,
                    paste(hd$eventDate1,hd$eventDate2,sep="/"),
                    hd$eventDate1)
hd=select(hd,-eventDate1,-eventDate2)

#Location, coordinates
hd=select(hd,-stateProvince2,-county2,-country2)
hd=select(hd,-continent2,-habitat2)
hd$country=hd$country1
hd=select(hd,-country1)
hd$county=ifelse(is.na(hd$county3)==F&hd$county1!=hd$county3,
                 paste(hd$county1,hd$county3,sep=", "),
                 hd$county1)
hd=select(hd,-county1,-county3)

hd$geodeticDatum[hd$verbatimCoordinateSystem2=="MY.coordinateSystemDms"]="Dms"
hd=select(hd,-verbatimCoordinateSystem2)
hd.c=select(hd,
            verbatimCoordinates,
            decimalLatitudeMin,
            decimalLatitudeMax,
            decimalLongitudeMin,
            decimalLongitudeMax,
            verbatimLongitude,
            verbatimCoordinates2,
            verbatimLatiude,decimalLatitude,
            decimalLongitude)
hd$verbatimCoordinates2[is.na(hd$verbatimCoordinates2)==T]=
  hd$verbatimCoordinates[is.na(hd$verbatimCoordinates2)==T]
hd$verbatimCoordinates=hd$verbatimCoordinates2
hd=select(hd,-verbatimCoordinates2)
hd=select(hd,
          -decimalLatitudeMax,
          -decimalLatitudeMin,
          -decimalLongitudeMax,
          -decimalLongitudeMin)

hd$verbatimCoordinates[is.na(hd$verbatimLatitude)==F&
                         is.na(hd$verbatimLongitude)==T]=
  hd$verbatimLatitude[is.na(hd$verbatimLatitude)==F&
                        is.na(hd$verbatimLongitude)==T]

hd$verbatimLatitude[is.na(hd$verbatimLatitude)==F&
                      is.na(hd$verbatimLongitude)==T]=NA

hd=select(hd,-coordinateUncertaintyInMeters2)
hd$verbatimLocality2[is.na(hd$verbatimLocality2)==T]=
  hd$locality1[is.na(hd$verbatimLocality2)==T]

hd=rename(hd,"locality"="verbatimLocality2")
hd=select(hd,-locality1,-locality2)

#organismRemarks
hd$organismRemarks=hd$fieldNotes
#typeStatus
hd=select(hd,-isType)
hd$typeStatus=gsub("MY.typeStatus","",hd$typeStatus)
hd=select(hd,-fieldNotes)

#fix countryCode into something consistent
hd$countryCode=hd$country
hd$countryCode=gsub("Australia","AU",hd$countryCode)
hd$countryCode=gsub("BOLIVIA","BO",hd$countryCode)
hd$countryCode=gsub("BRAZIL","BR",hd$countryCode)
hd$countryCode=gsub("Brazil","BR",hd$countryCode)
hd$countryCode=gsub("China","CN",hd$countryCode)
hd$countryCode=gsub("Ethiopia","ET",hd$countryCode)
hd$countryCode=gsub("Finland","FI",hd$countryCode)
hd$countryCode=gsub("Indonesia","ID",hd$countryCode)
hd$countryCode=gsub("MALAYSIA","MY",hd$countryCode)
hd$countryCode=gsub("United States","US",hd$countryCode)
helsinki.data=hd