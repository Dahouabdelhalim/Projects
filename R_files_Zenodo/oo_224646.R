setwd("YOUR WORKING DIRECTORY")

library(tidyverse)
library(rgbif)
library(rworldmap)

###Compile the dataset
##All but H and B are taken from GBIF
edinburgh=read_tsv("sample set Edinburgh.csv")
kew=read_tsv("sample set Kew.csv")
london=read_csv2("sample set London.csv")
meise=read_csv2("sample set Meise.csv")
naturalis=read_csv2("sample set Naturalis.csv")
paris=read_csv2("sample set Paris.csv")
tartu=read_csv2("sample set Tartu.csv")

##H and B taken from data files
helsinki.data=read_csv("helsinkidata api.csv")
berlin.data=read_csv("berlin data.txt")

#Extract catalogNumbers for P
paris$catalogNumber=substr(paris$SPECIMEN,43,100)
paris$catalogNumber[184:200]=substr(paris$SPECIMEN[184:200],44,100)

#Retrieve occurrence records from GBIF
paris.all=occ_search(catalogNumber = paris$catalogNumber)
edinburgh.all=occ_search(catalogNumber = edinburgh$catalognumber)
kew.all=occ_search(catalogNumber = kew$catalognumber)
london.all=occ_search(catalogNumber = london$barcode)
tartu.all=occ_search(catalogNumber = tartu$catalognumber)
naturalis.all=occ_search(catalogNumber = naturalis$catalognumber)
meise.all=occ_search(catalogNumber = meise$catalogNumber)

#Turn them into a tibble
paris.data=paris.all[[1]]$data
for (i in 2:200) {
  paris.data=merge(paris.data,paris.all[[i]]$data,all.x=T,all.y=T)
}
edinburgh.data=edinburgh.all[[1]]$data
for (i in 2:200) {
  edinburgh.data=merge(edinburgh.data,edinburgh.all[[i]]$data,all.x=T,all.y=T)
}
kew.data=kew.all[[1]]$data
for (i in 2:200) {
  kew.data=merge(kew.data,kew.all[[i]]$data,all.x=T,all.y=T)
}
london.data=london.all[[1]]$data
for (i in 2:200) {
  london.data=merge(london.data,london.all[[i]]$data,all.x=T,all.y=T)
}
tartu.data=tartu.all[[1]]$data
for (i in 2:200) {
  tartu.data=merge(tartu.data,tartu.all[[i]]$data,all.x=T,all.y=T)
}
meise.data=meise.all[[1]]$data
for (i in 2:200) {
  if (is.null(meise.all[[i]]$data)==F) {
    meise.data=merge(meise.data,meise.all[[i]]$data,all.x=T,all.y=T)
  }
}
naturalis.data=naturalis.all[[1]]$data
for (i in 2:200) {
  naturalis.data=merge(naturalis.data,naturalis.all[[i]]$data,all.x=T,all.y=T)
}

#Add a setID variable to easily distinguish the origin
#Also remove older duplicates in the Meise due to obsolete but not yet removed GBIF datasets
paris.data$setID="P"
kew.data$setID="K"
london.data$setID="BM"
tartu.data$setID="TU"
meise.data$setID="BR"
meise.data=filter(meise.data,is.na(preparations)==F)
berlin.data$setID="B"
edinburgh.data$setID="E"
naturalis.data$setID="L"
helsinki.data$setID="H"

#Extract year for the graph. done here already to simplify the formula
berlin.data$year=ifelse(is.na(berlin.data$verbatimEventDate)==T,
                         NA,
                         ifelse(substr(berlin.data$verbatimEventDate,3,3)!="."&
                                  substr(berlin.data$verbatimEventDate,2,2)!=".",
                                substr(berlin.data$verbatimEventDate,1,4),
                                substr(berlin.data$verbatimEventDate,
                                       nchar(berlin.data$verbatimEventDate)-3,
                                       nchar(berlin.data$verbatimEventDate))))

edinburgh.data$year=ifelse(is.na(edinburgh.data$verbatimEventDate)==F&
                             is.na(edinburgh.data$eventDate)==T,
                           substr(edinburgh.data$verbatimEventDate,
                                  nchar(edinburgh.data$verbatimEventDate)-3,
                                  nchar(edinburgh.data$verbatimEventDate)),
                           edinburgh.data$year)

helsinki.data$catalogNumber=gsub("http://id.luomus.fi/","",helsinki.data$occurrenceID)

##Merge all into one dataset
alles=merge(paris.data,meise.data,all.x=T,all.y=T)
alles=merge(alles,london.data,all.x=T,all.y=T)
alles=merge(alles,kew.data,all.x=T,all.y=T)
alles=merge(alles,tartu.data,all.x=T,all.y=T)
alles=merge(alles,edinburgh.data,all.x=T,all.y=T)
alles=merge(alles,berlin.data,all.x=T,all.y=T)
alles=merge(alles,naturalis.data,all.x=T,all.y=T)
alles=merge(alles,helsinki.data,all.x=T,all.y=T)

#Filter out the extensions
alles2=alles[,-grep("unknown.org",colnames(alles))]

#Write to file
write_csv(alles2,"all data.csv",na="")

################################################################

###Graph of years
#sort by institute, then give each specimen an index of 1-200 per institute
alles.y=arrange(alles,setID)
alles.y$index=as.numeric(rownames(alles.y))-((as.numeric(rownames(alles.y))-1)%/%200)*200
#extract year from event date if not done already
alles.y$year[is.na(alles.y$year)==T&
               is.na(alles.y$eventDate)==F]=
  substr(alles.y$eventDate[is.na(alles.y$year)==T&
                             is.na(alles.y$eventDate)==F],1,4)
#extract year from Naturalis verbatimEventDate
alles.y$year[alles.y$setID=="L"&
               is.na(alles.y$verbatimEventDate)==F]=
  as.integer(substr(alles.y$verbatimEventDate[alles.y$setID=="L"&
                                                is.na(alles.y$verbatimEventDate)==F],1,4))
alles.y=filter(alles.y,is.na(year)==F)
alles.y$year=as.integer(alles.y$year)

#Calculate 10 year counts
breaks=seq(1730,2020,10)
X1=seq(1,29)
X2=X1
hists=tibble(X1)
hists$B=hist(as.numeric(alles.y$year[is.na(alles.y$year)==F&
                                       alles.y$setID=="B"]),breaks=breaks)$counts
hists$BM=hist(as.numeric(alles.y$year[is.na(alles.y$year)==F&
                                        alles.y$setID=="BM"]),breaks=breaks)$counts
hists$BR=hist(as.numeric(alles.y$year[is.na(alles.y$year)==F&
                                        alles.y$setID=="BR"]),breaks=breaks)$counts
hists$E=hist(as.numeric(alles.y$year[is.na(alles.y$year)==F&
                                       alles.y$setID=="E"]),breaks=breaks)$counts
hists$K=hist(as.numeric(alles.y$year[is.na(alles.y$year)==F&
                                       alles.y$setID=="K"]),breaks=breaks)$counts
hists$L=hist(as.numeric(alles.y$year[is.na(alles.y$year)==F&
                                       alles.y$setID=="L"]),breaks=breaks)$counts
hists$P=hist(as.numeric(alles.y$year[is.na(alles.y$year)==F&
                                       alles.y$setID=="P"]),breaks=breaks)$counts
hists$TU=hist(as.numeric(alles.y$year[is.na(alles.y$year)==F&
                                        alles.y$setID=="TU"]),breaks=breaks)$counts
hists$H=hist(as.numeric(alles.y$year[is.na(alles.y$year)==F&
                                       alles.y$setID=="H"]),breaks=breaks)$counts
hists$mids=hist(as.numeric(alles.y$year[is.na(alles.y$year)==F&
                                          alles.y$setID=="TU"]),breaks=breaks)$mids
#stackk all in one variable
ydata=stack(select(hists,-mids,-X1))
ydata$index=rep(hists$mids,9)
ydata=rename(ydata,"setID"="ind")

#derive the number of years per institute
ydata$sums=ydata$values
for (i in 3:11) {
  ydata$sums[ydata$setID==colnames(hists)[i]]=sum(hists[,i])
}

for (i in 1:length(alles.y$basisOfRecord)) {
  alles.y$nlabel[i]=length(filter(alles.y,setID==alles.y$setID[i])$setID)
}
alles.y$nlabel=paste0("n=",alles.y$nlabel)

#make the plot
ggplot() + 
  geom_bar(data=ydata,
           mapping=aes(x=index,
                       y=200,
                       fill=values,
                       color=values),
           width=10,
           stat="identity",
           position=position_dodge(width=0)) +
  scale_fill_gradient("Counts",high="red",low="steelblue") +
  scale_color_gradient("Counts",high="red",low="steelblue") +
  geom_point(data=alles.y,aes(x=year,y=index),alpha=0.3) +
  facet_grid(setID~.) +
  scale_x_continuous(name="Year",breaks=c(1750,1800,1850,1900,1950,2000)) +
  ylab("") +
  theme_classic() +
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=16),
        strip.text.y=element_text(size=16)) +
  geom_text(data=alles.y,mapping=aes(x=1720,y=100,label=nlabel)) 

###############################################################

###World graph
#set up coordinates
alles2.c=filter(alles,is.na(decimalLatitude)==F)
alles2.c$decimalLatitude=as.numeric(alles2.c$decimalLatitude)
alles2.c$decimalLongitude=as.numeric(alles2.c$decimalLongitude)

#set up countries
map.world=map_data(map="world")
map.world$countryCode=iso.alpha(map.world$region)
cccount=count(alles,countryCode)
cccount=filter(cccount,countryCode!="none",countryCode!="ZZ")
map.world$count=0
for (i in 1:length(cccount$countryCode)) {
  map.world$count[map.world$countryCode%in%cccount$countryCode[i]]=cccount$n[i]
}
ggplot() + 
  geom_polygon(data=map.world,aes(x=long,y=lat,group=group,fill=count),color="black") +
  scale_fill_gradient(low="white",high="gray19", name="# of specimens") +
  geom_point(data=alles2.c,aes(x=decimalLongitude,y=decimalLatitude),color="black") +
  coord_map("mollweide",xlim=c(-180,180)) +
  scale_x_continuous(breaks=c(-100,100)) +
  theme_grey() +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=16),
        strip.text.y=element_text(size=16))

################################################################

###Taxo graph
alles.t=alles

#set family to lower case (first letter capitalized)
alles.t$family[is.na(alles.t$family)==F]=
  paste0(substr(alles.t$family[is.na(alles.t$family)==F],1,1),
         tolower(substr(alles.t$family[is.na(alles.t$family)==F],
                        2,
                        nchar(alles.t$family[is.na(alles.t$family)==F]))))

#do away with combined family names
alles.t$family=gsub("Leguminosae-caesalpinioideae","Leguminosae",alles.t$family)
alles.t$family=gsub("Leguminosae-papilionoideae","Leguminosae",alles.t$family)

#get higher classification info from GBIF backbone
for (i in 1:length(alles.t$class)) {
  if (is.na(alles.t$phylum[i])==T&is.na(alles.t$kingdom[i])==F) {
    k=name_backbone(name = alles.t$family[i])
    if (is.null(k$kingdom)==F) {
      alles.t$class[i]=k$class
      alles.t$order[i]=k$order
      alles.t$phylum[i]=k$phylum
    }
  }
}

#no family for Helsinki data, use genus
for (i in 1:length(alles.t$basisOfRecord)) {
  if (alles.t$setID[i]=="H") {
    k=name_backbone(alles.t$genus[i])
    if (is.null(k$kingdom)==F) {
      alles.t$kingdom[i]=k$kingdom
      alles.t$family[i]=k$family
      alles.t$class[i]=k$class
      alles.t$order[i]=k$order
      alles.t$phylum[i]=k$phylum
    }
  }
}
#tally them for Krona graph, to be imported in the excel macro
list=count(alles.t,phylum,order,class,family)
list2=filter(list,is.na(phylum)==F)
write_csv(list2,"krona dat all.csv",na="")

##########################################################

###Graph of languages
#retrieve the sorted data
lang=read_csv2("icedig lang.csv")

#reformat them
langp=stack(lang[,2:10])
langp$l=rep(lang$language,9)
langp$tot=rep(lang$total,9)
langp$l=factor(langp$l)
#get percentages and sort from high to low
langp$tot=langp$tot*-1
langp$values2=100*langp$values/1800
langp$tot[langp$l=="ZZ"]=0
langp$l=reorder(langp$l,langp$tot)

#make graph
ggplot(langp) + 
  geom_bar(aes(x=l,y=values2,fill=ind),position="stack",stat="identity") + 
  theme_classic() +   
  scale_fill_manual(values=c("darkgreen","seagreen4","springgreen3",
                             "palegreen3","chartreuse4","olivedrab4",
                             "lightseagreen","deepskyblue4","royalblue3")) + 
  labs(title="",x="Language code",y="%",legend.title="") +
  theme(axis.text.x=element_text(size=16),
        axis.title.x=element_text(size=20),
        title=element_text(size=20),
        legend.text=element_text(size=16),legend.title=element_blank())

###################################################

###Export metadata
#We omit some variables less important to label text.
alles.s=select(alles,
               -name,
               -key,
               -issues,
               -datasetKey,
               -publishingOrgKey,
               -networkKeys,
               -installationKey,
               -publishingCountry,
               -protocol,
               -lastCrawled,
               -lastParsed,
               -crawlId,
               -extensions,
               -taxonKey,
               -lastInterpreted,
               -license,
               -identifiers,
               -facts,
               -relations,
               -gbifID,
               -kingdomKey,
               -phylumKey,
               -classKey,
               -orderKey,
               -familyKey,
               -genericName,
               -genusKey,
               -speciesKey,
               -species,
               -typifiedName,
               -dynamicProperties,
               -language,
               -subgenus)
#remove the DwC extension variables
alles.s=alles.s[,grepl("unknown.org",colnames(alles.s))==F]

#add the language as classified (and languageRemarks if language was unclear)
lang=read_csv("languageperspecimen.csv")
alles.s=left_join(alles.s,lang,by="catalogNumber")

#map terms to their overarching category
terms=read_csv2("mapping metadata.csv")
terms$cat[is.na(terms$cat)==T]=""

#remove meaningless time values from dates
timez=strsplit(alles.s$eventDate,split="T")
alles.s$eventDate=as.character(sapply(seq(1,length(alles.s$eventDate)),function(x) timez[[x]][1]))

timez=strsplit(alles.s$modified,split="T")
alles.s$modified=as.character(sapply(seq(1,length(alles.s$eventDate)),function(x) timez[[x]][1]))

timez=strsplit(alles.s$dateIdentified,split="T")
alles.s$dateIdentified=as.character(sapply(seq(1,length(alles.s$dateIdentified)),function(x) timez[[x]][1]))

#language codes to ISO-2, for Zenodo

library(plyr)
alles.s$language=revalue(alles.s$language,c("DE"="deu",
                                            "EN"="eng",
                                            "ES"="spa",
                                            "ET"="est",
                                            "FI"="fin",
                                            "FR"="fra",
                                            "IT"="ita",
                                            "LA"="lat",
                                            "NL"="nld",
                                            "PT"="por",
                                            "RU"="rus",
                                            "SV"="swe",
                                            "ZZ"="und"))

#shuttle all data into a list to prepare for json structure
#don't include empty values (NA)
test=list()
for (i in 1:length(alles.s$eventDate)) {
  test[[i]]=select_if(select(alles.s[i,],terms$terms[terms$cat==""]),~ !any(is.na(.)))
  names(test)[i]=alles.s$catalogNumber[i]
  test[[i]]$Occurrence=select_if(select(alles.s[i,],terms$terms[terms$cat=="Occurrence"]),~ !any(is.na(.)))
  test[[i]]$Event=select_if(select(alles.s[i,],terms$terms[terms$cat=="Event"]),~ !any(is.na(.)))
  test[[i]]$Location=select_if(select(alles.s[i,],terms$terms[terms$cat=="Location"]),~ !any(is.na(.)))
  test[[i]]$Identification=select_if(select(alles.s[i,],terms$terms[terms$cat=="Identification"]),~ !any(is.na(.)))
  test[[i]]$Taxon=select_if(select(alles.s[i,],terms$terms[terms$cat=="Taxon"]),~ !any(is.na(.)))
  test[[i]]$Organism=select_if(select(alles.s[i,],terms$terms[terms$cat=="Organism"]),~ !any(is.na(.)))
  test[[i]]$ResourceRelationship=select_if(select(alles.s[i,],terms$terms[terms$cat=="ResourceRelationship"]),~ !any(is.na(.)))
  test[[i]]$Custom=select_if(select(alles.s[i,],terms$terms[terms$cat=="Custom"]),~ !any(is.na(.)))
  test[[i]]$MeasurementOrFact=select_if(select(alles.s[i,],terms$terms[terms$cat=="MeasurementOrFact"]),~ !any(is.na(.)))
}
#export json files
library(jsonlite)
setwd("SOME DIR")
for (i in 1:length(alles.s$basisOfRecord)) {
  k=toJSON(test[[i]])
  k=substr(k,2,nchar(k)-1) #remove [] brackets, may cause issues in python
  k2=prettify(k,indent=4)
  write(k2,paste0(names(test)[i],".json"))
}

#export csv table
write_csv(alles.s,"full metadata.csv",na="")