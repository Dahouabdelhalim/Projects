
### Replace NA values in data files with values from the previous row where needed
no.replace.na <- c(17,19,20,21,22,24,25,26) # identify columns in step 4 files whose NA values should not be replaced with values from the previous row

att3 <- rbind(att3.wos,att3.no.wos)
att4 <- rbind(att4.wos,att4.no.wos)

# Identify, within each study, which drivers actually have ranks
study.drivers <- data.frame(Study = att3$UI)
study.drivers$CC <- ifelse(att3$Assessment.of.climate.change.impact == "Yes", TRUE, FALSE)
study.drivers$DE <- ifelse(att3$Assessment.of.resource.extraction.impact == "Yes", TRUE, FALSE)
study.drivers$IAS <- ifelse(att3$Assessment.of.invasive.alien.species.impact == "Yes", TRUE, FALSE)
study.drivers$LU <- ifelse(att3$Assessment.of.land.sea.use.change.impact == "Yes", TRUE, FALSE)
study.drivers$PO <- ifelse(att3$Assessment.of.pollution.impact == "Yes", TRUE, FALSE)
study.drivers$OT <- ifelse(att3$Assessment.of.other.driver.s..impact == "Yes", TRUE, FALSE)
if (exclude.other.drivers == "Yes") study.drivers$OT <- FALSE # We are not considering other drivers
saveRDS(study.drivers, "output/study.drivers.Rds")

levels(att4$Taxonomic.group..2.) <- c(levels(att4$Taxonomic.group..2.),"Not applicable")
levels(att4$Original.indicator.name) <- c(levels(att4$Original.indicator.name),"Not applicable")
levels(att4$Relevance.of.direct.driver.magnitude.assessment) <- c(levels(att4$Relevance.of.direct.driver.magnitude.assessment),"Not applicable")

for (i in 1:nrow(att4)) {  
  if (!is.na(att4$Indicator.name[i])) {
    if (is.na(att4$Original.indicator.name[i])) {
      att4$Original.indicator.name[i] <- "Not applicable"
    }
    if (is.na(att4$Taxonomic.group..1.[i])) {
      att4$Taxonomic.group..1.[i] <- "Not specified"
    }
    if (is.na(att4$Relevance.of.direct.driver.magnitude.assessment[i])) {
      att4$Relevance.of.direct.driver.magnitude.assessment[i] <- "Not applicable"
    }
}
  if (!is.na(att4$Taxonomic.group..1.[i])) {
    if (is.na(att4$Taxonomic.group..2.[i])) {
      att4$Taxonomic.group..2.[i] <- "Not applicable"
    }
  }
  col <- which(is.na(att4[i,]) | att4[i,]=='')
  col <- col[-which(col %in% no.replace.na)]
  if (is.na(att4$Unique.identifier..UI.[i])){ 
    att4[i,][col] <- att4[i-1,][col]   
  } else {
    att4[i,][col] <- att4[i,][col]
  }
}

att <- merge(att3,att4,by.x=c("UI"),by.y=c("Unique.identifier..UI."),all.x=F,all.y=F)
att$IPBES.direct.driver <- as.factor(att$IPBES.direct.driver) #Added 2021-06-02
### Rename levels of some factors

# Drivers
levels(att$IPBES.direct.driver) <- c(levels(att$IPBES.direct.driver),"Direct exploitation") 
att$IPBES.direct.driver[att$IPBES.direct.driver=="Resource extraction"]  <- "Direct exploitation"
old.lvl <- levels(att$IPBES.direct.driver)

if (exclude.other.drivers == TRUE){
    # This condition and the next line added 2020-09-02 to tidy up conditional removal of Other as a direct driver; previously the 'else' condition was unconditional
    att$IPBES.direct.driver <- factor(att$IPBES.direct.driver,levels=sort(old.lvl[old.lvl!="Other"],decreasing=F))
  }else{
    att$IPBES.direct.driver <- factor(att$IPBES.direct.driver,levels=c(sort(old.lvl[old.lvl!="Other"],decreasing=F),"Other"))
  }
att$IPBES.direct.driver <- factor(att$IPBES.direct.driver,levels=levels(factor(att$IPBES.direct.driver)))

# EBV classes
levels(att$EBV.class) <- c(levels(att$EBV.class),"Genetic\\ncomposition","Community\\ncomposition","Ecosystem\\nfunction","Ecosystem\\nstructure","Species\\npopulations","Species\\ntraits") 
att$EBV.class[att$EBV.class=="Genetic composition"]  <- "Genetic\\ncomposition"
att$EBV.class[att$EBV.class=="Community composition"]  <- "Community\\ncomposition"
att$EBV.class[att$EBV.class=="Ecosystem function"]  <- "Ecosystem\\nfunction"
att$EBV.class[att$EBV.class=="Ecosystem structure"]  <- "Ecosystem\\nstructure"
att$EBV.class[att$EBV.class=="Species populations"]  <- "Species\\npopulations"
att$EBV.class[att$EBV.class=="Species traits"]  <- "Species\\ntraits"
old.lvl <- levels(att$EBV.class)
att$EBV.class <- factor(att$EBV.class,levels=c("Genetic\\ncomposition","Species\\npopulations","Species\\ntraits","Community\\ncomposition","Ecosystem\\nstructure","Ecosystem\\nfunction"))

### Selection of subsets of data to be included in the analysis depending on options specified above by the user

# Realm selection
if(realm.selection=="Yes"){
  selected.realm.vector <- att[grep(selected.realm,att$Realm.s),]
  selected.realm.vector <- as.vector(unique(selected.realm.vector$Realm.s))
  if(include.all.realms=="Yes"){
    selected.realm.vector <- c(selected.realm.vector,"All realms")
  } else {
    selected.realm.vector <- selected.realm.vector
  }
  att <- att[is.element(att$Realm.s,selected.realm.vector),]
  print(paste("Selected realms:", paste(selected.realm.vector, collapse=", ")))
}

# Region selection
if(region.selection=="Yes"){
  selected.region.vector <- att[grep(selected.region,att$IPBES.Region.s),]
  selected.region.vector <- as.vector(unique(selected.region.vector$IPBES.Region.s))
  if(include.all.regions=="Yes"){
    selected.region.vector <- c(selected.region.vector,"All regions")
  } else {
    selected.region.vector <- selected.region.vector
  }
  att <- att[is.element(att$IPBES.Region.s,selected.region.vector),]
  print(paste("Set of regions included:", paste(selected.region.vector, collapse=", ")))
  
}

# Domain selection
if(domain.selection=="Yes"){
  selected.domain.vector <- att[grep(selected.domain,att$Climate.domain.s),]
  selected.domain.vector <- as.vector(unique(selected.domain.vector$Climate.domain.s))
  if(include.all.domains=="Yes"){
    selected.domain.vector <- c(selected.domain.vector,"All domains")
  } else {
    selected.domain.vector <- selected.domain.vector
  }
  att <- att[is.element(att$Climate.domain.s,selected.domain.vector),]
  print(paste("Selected domains:", paste(selected.domain.vector, collapse=", ")))
  
}

# Manual selection
if(manual.selection=="Yes"){
  if(region.selection=="Yes" & realm.selection=="Yes"){
    att <- att[is.element(att$Include,c("All","Realm","Region","Region/Realm")),]
  } else if(region.selection=="Yes" & realm.selection=="No") {
    att <- att[is.element(att$Include,c("All","Region","Region/Realm")),]
  } else if(region.selection=="No" & realm.selection=="Yes") {
    att <- att[is.element(att$Include,c("All","Realm","Region/Realm")),]
  } else if(region.selection=="No" & realm.selection=="No") {
    att <- att[is.element(att$Include,c("All")),]
  }
}

# Indicator selection
if(indicator.selection=="Yes"){
  att <- att[is.element(att$Indicator.name,selected.indicator),]
}
if(indicator.exclusion=="Yes"){
  att <- att[!is.element(att$Indicator.name,excluded.indicator),]
}

# Scale selection
if(scale.selection=="Yes"){
  att <- att[is.element(att$Spatial.coverage,selected.scale),]
}

# Number of drivers assessed
if(n.drivers.selection=="Yes"){
  att <- att[is.element(att$Number.of.driver.s..analysed.assessed,selected.n.drivers),]
}

### Calculate statistics related to the numbers of studies and assessments* included in the analysis
# * each study is identified with a unique identifier (UI) and includes one or several assessments (Levels.of.analysis)

n.assessments.per.indicator.x.driver <- att %>% group_by(Indicator.name,IPBES.direct.driver) %>% summarise(N=n())

n.studies.per.indicator <- att %>% group_by(UI,Indicator.name) %>% summarise(N=n())
n.studies.per.indicator <- n.studies.per.indicator %>% group_by(Indicator.name) %>% summarise(N=n())
n.studies.per.indicator.x.driver <- att %>% group_by(UI,Indicator.name,IPBES.direct.driver) %>% summarise(N=n())
n.studies.per.indicator.x.driver <- n.studies.per.indicator.x.driver %>% group_by(Indicator.name,IPBES.direct.driver) %>% summarise(N=n())

n.studies.per.indicator.x.ebv <- att %>% group_by(UI,Indicator.name,EBV.class) %>% summarise(N=n())
n.studies.per.indicator.x.ebv <- n.studies.per.indicator.x.ebv %>% group_by(Indicator.name,EBV.class) %>% summarise(N=n())
n.studies.per.indicator.x.ebv.x.driver <- att %>% group_by(UI,Indicator.name,EBV.class,IPBES.direct.driver) %>% summarise(N=n())
n.studies.per.indicator.x.ebv.x.driver <- n.studies.per.indicator.x.ebv.x.driver %>% group_by(Indicator.name,EBV.class,IPBES.direct.driver) %>% summarise(N=n())

n.studies.per.ebv <- att %>% group_by(UI,EBV.class) %>% summarise(N=n())
n.studies.per.ebv <- n.studies.per.ebv %>% group_by(EBV.class) %>% summarise(N=n())
n.studies.per.ebv.x.driver <- att %>% group_by(UI,EBV.class,IPBES.direct.driver) %>% summarise(N=n())
n.studies.per.ebv.x.driver <- n.studies.per.ebv.x.driver %>% group_by(EBV.class,IPBES.direct.driver) %>% summarise(N=n())

n.studies <- length(unique(att$UI)) # total number of studies included in the analysis
att$Unique.identifier.level <- paste(att$UI,att$Level.of.analysis,sep="-")
n.assessments <- length(unique(att$Unique.identifier.level)) # total number of assessments included in the analysis

### Add rows for drivers whose impact is not assessed for each assessment within each study
stats.ranks <- att %>% group_by(Unique.identifier.level) %>% summarise(Min.rank=min(Direct.driver.rank,na.rm=TRUE),Max.rank=max(Direct.driver.rank,na.rm=T),Avg.rank=mean(Direct.driver.rank,na.rm=T),N.direct.drivers=n())
att.other <- att[is.element(att$IPBES.direct.driver,c("Other")),]
assessments.with.other <- att.other %>% group_by(Unique.identifier.level) %>% summarise(Other.direct.driver=n())

ui.level.list <- as.vector(unique(att$Unique.identifier.level))
driver.list <- as.vector(unique(att$IPBES.direct.driver))
if (exclude.other.drivers == "Yes") driver.list <- driver.list[driver.list != "Other"] # Remove Other from the set
ui.level.x.driver.list <- expand.grid(driver.list=driver.list,ui.level.list=ui.level.list,stringsAsFactors=TRUE)
colnames(ui.level.x.driver.list) <- c("IPBES.direct.driver","Unique.identifier.level")
att <- att[with(att,order(Unique.identifier.level,IPBES.direct.driver)),]
ui.level.x.driver.list <- ui.level.x.driver.list[with(ui.level.x.driver.list,order(Unique.identifier.level,IPBES.direct.driver)),]
att <- merge(att,ui.level.x.driver.list,by=c("Unique.identifier.level","IPBES.direct.driver"),all=TRUE)

