
att <- merge(att,stats.ranks,by.x=c("Unique.identifier.level"),by.y=c("Unique.identifier.level"),all.x=T,all.y=T)
att <- merge(att,assessments.with.other,by.x=c("Unique.identifier.level"),by.y=c("Unique.identifier.level"),all.x=T,all.y=T)

for (i in 1:nrow(att)) {
  if (is.na(att$Other.direct.driver[i])){ 
    att$Other.direct.driver[i] <- 0  
  }
}

if (exclude.other.drivers == "Yes"){
  att <- subset(att, IPBES.direct.driver != "Other")
}

att$N.IPBES.direct.drivers <- att$N.direct.drivers

for (i in 1:nrow(att)) {
  if(include.non.assessed.drivers=="No"){
    att$Direct.driver.rank.reversed[i] <- att$Max.rank[i]-att$Direct.driver.rank[i]+1
    if (is.na(att$Direct.driver.rank[i])){
      att$Direct.driver.rank.reversed[i] <- 0 # non assessed drivers have NO impact 
    }
  }   
  if(include.non.assessed.drivers=="Less"){
    if (is.na(att$Direct.driver.rank[i])){
      att$Direct.driver.rank[i] <- att$Max.rank[i]+1 # non assessed drivers have LESS impact than any assessed driver
    }
    if (att$N.direct.drivers[i]==6){ 
      att$Direct.driver.rank.reversed[i] <- att$Max.rank[i]-att$Direct.driver.rank[i]+1  
    } else { 
      att$Direct.driver.rank.reversed[i] <- att$Max.rank[i]+1-att$Direct.driver.rank[i]+1  
    }
  }
  if(include.non.assessed.drivers=="Equal"){
    if (is.na(att$Direct.driver.rank[i])){
      att$Direct.driver.rank[i] <- att$Max.rank[i] # non assessed drivers have SAME impact as lowest-ranked assessed driver
    }
    att$Direct.driver.rank.reversed[i] <- att$Max.rank[i]-att$Direct.driver.rank[i]+1
  } 
  if(include.non.assessed.drivers=="Unknown"){
    if (is.na(att$Direct.driver.rank[i])){
      att$Direct.driver.rank[i] <- att$Avg.rank[i] # non assessed drivers have UNKNOWN impact
    }
    att$Direct.driver.rank.reversed[i] <- att$Max.rank[i]-att$Direct.driver.rank[i]+1
  } 
}  
  
att <- att[with(att,order(Unique.identifier.level,UI)),]

for (i in 1:nrow(att)) {  
  col <- which(is.na(att[i,]))
  col <- col[-which(col %in% c(62,63,64,65,67,68))]
  att[i,][col] <- att[i-1,][col]   
}

### Rescale the ranking of drivers (--> sum of ranks = 1 for each assessment within each study)

stats.ranks <- att %>% group_by(UI,Level.of.analysis) %>% summarise(Sum.rank=sum(Direct.driver.rank.reversed,na.rm=T))
att <- merge(att,stats.ranks,by.x=c("UI","Level.of.analysis"),by.y=c("UI","Level.of.analysis"),all.x=T,all.y=T)
att$Direct.driver.rank.rescaled <- att$Direct.driver.rank.reversed/att$Sum.rank
stats.ranks.rescaled <- att %>% group_by(UI,Level.of.analysis) %>% summarise(Sum.rank.rescaled=sum(Direct.driver.rank.rescaled,na.rm=T))
att <- merge(att,stats.ranks.rescaled,by.x=c("UI","Level.of.analysis"),by.y=c("UI","Level.of.analysis"),all.x=T,all.y=T)

### Rescale the magnitude of drivers (-->sum of magnitudes = 1 for each assessment within each study)

stats.magnitude <- att %>% group_by(UI,Level.of.analysis) %>% summarise(Sum.magnitude=sum(Direct.driver.magnitude,na.rm=T))
att <- merge(att,stats.magnitude,by.x=c("UI","Level.of.analysis"),by.y=c("UI","Level.of.analysis"),all.x=T,all.y=T)
att$Direct.driver.magnitude.rescaled <- att$Direct.driver.magnitude/att$Sum.magnitude
stats.magnitude.rescaled <- att %>% group_by(UI,Level.of.analysis) %>% summarise(Sum.magnitude.rescaled=sum(Direct.driver.magnitude.rescaled,na.rm=T))
att <- merge(att,stats.magnitude.rescaled,by.x=c("UI","Level.of.analysis"),by.y=c("UI","Level.of.analysis"),all.x=T,all.y=T)


for (i in 1:nrow(att)) {   
  if (is.na(att$Direct.driver.magnitude.rescaled[i])){
    if (att$Sum.magnitude.rescaled[i]==0){
      att$Direct.driver.importance.rescaled[i] <- att$Direct.driver.rank.rescaled[i]
    } else {
      att$Direct.driver.magnitude.rescaled[i] <- 0
      if(combine.rank.magnitude=="Yes"){
        if(att$Relevance.of.direct.driver.magnitude.assessment[i]=="Yes"){
          att$Direct.driver.importance.rescaled[i] <- att$Direct.driver.magnitude.rescaled[i]
        } else {
          att$Direct.driver.importance.rescaled[i] <- att$Direct.driver.rank.rescaled[i]
        } 
      } else {
        att$Direct.driver.importance.rescaled[i] <- att$Direct.driver.rank.rescaled[i]
      }
    }
  } else {
    if(combine.rank.magnitude=="Yes"){
      if(att$Relevance.of.direct.driver.magnitude.assessment[i]=="Yes"){
        att$Direct.driver.importance.rescaled[i] <- att$Direct.driver.magnitude.rescaled[i]
      } else {
        att$Direct.driver.importance.rescaled[i] <- att$Direct.driver.rank.rescaled[i]
      } 
    } else {
      att$Direct.driver.importance.rescaled[i] <- att$Direct.driver.rank.rescaled[i]
    }
  }
}
stats.importance.rescaled <- att %>% group_by(UI,Level.of.analysis) %>% summarise(Sum.importance.rescaled=sum(Direct.driver.importance.rescaled,na.rm=T))
att <- merge(att,stats.importance.rescaled,by.x=c("UI","Level.of.analysis"),by.y=c("UI","Level.of.analysis"),all.x=T,all.y=T)

### Remove redundant information among different assessments within each study (based on column "Sublevel of analysis")
#   --> assessments within studies may result from the same analysis carried out at varying levels (e.g. for different taxa AND overall across all studied taxa)
#   --> only the results from the analysis carried out at the highest level is retained as non-redundant source of information

sublevel <- aggregate(data=att,Sublevel.of.analysis~UI, function(x) length(unique(x)))
colnames(sublevel) <- c("UI","Types.sublevel.of.analysis")
att <- merge(att,sublevel,by="UI",all=T)
stats.sublevel <- sublevel %>% summarise(Min.sublevel=min(Types.sublevel.of.analysis),Max.sublevel=max(Types.sublevel.of.analysis))
if(stats.sublevel$Max.sublevel==2){
  att <- att[-which(att$Types.sublevel.of.analysis==2 & att$Sublevel.of.analysis=="Yes"),]
}

### Relative impacts of drivers aggregated (median) across the different assessments WITHIN each study included in the analysis
#   --> for each indicator (option "analysis.level" = indicators or overall)
#   --> for each indicator and ebv class (option "analysis.level" = ebv or indicators.x.ebv) 

if (analysis.level=="indicators" || analysis.level=="overall") { 
  attplot <- aggregate(Direct.driver.importance.rescaled~UI+IPBES.Region.s+Realm.s+Taxonomic.group..1.+Indicator.name+IPBES.direct.driver+N.IPBES.direct.drivers+Spatial.coverage,data=att,median)
  attplot <- aggregate(Direct.driver.importance.rescaled~UI+Indicator.name+IPBES.direct.driver+N.IPBES.direct.drivers+Spatial.coverage,data=attplot,median)
  attplot <- attplot[with(attplot,order(UI,Indicator.name,IPBES.direct.driver)),]
} else {
  attplot <- aggregate(Direct.driver.importance.rescaled~UI+IPBES.Region.s+Realm.s+Taxonomic.group..1.+Indicator.name+EBV.class+IPBES.direct.driver+N.IPBES.direct.drivers+Spatial.coverage,data=att,median)
  attplot <- aggregate(Direct.driver.importance.rescaled~UI+Indicator.name+EBV.class+IPBES.direct.driver+N.IPBES.direct.drivers+Spatial.coverage,data=attplot,median)
  attplot <- attplot[with(attplot,order(UI,Indicator.name,EBV.class,IPBES.direct.driver)),]
}

# Set weights using options chosen earlier
attplot$Scale.weight <- as.numeric(scale.weights[as.character(attplot$Spatial.coverage)]) #ERROR fixed on 2020-01-29 - was not ordering them correctly!
attplot$Ndrivers.weight <- ndriver.weights[attplot$N.IPBES.direct.drivers]

#Set OLD indicator reweights, intended to equalise the contributions of the indicators present in attplot to the David's score analysis
ind.in.data <- tapply(attplot$Scale.weight, attplot$Indicator.name, "sum")
ind.in.data <- ind.in.data/6 #The total indicator scale weight in the data
ave.ind <- mean(ind.in.data, na.rm=TRUE)

attplot$Old.indicator.reweight <- ave.ind/ind.in.data[attplot$Indicator.name]

if (analysis.level=="indicators" || analysis.level=="overall") { 
  attplot.long <- attplot %>% select(UI,Indicator.name,IPBES.direct.driver,Direct.driver.importance.rescaled, Scale.weight, Ndrivers.weight, Old.indicator.reweight)
  attplot.wide <- reshape(attplot.long,v.names="Direct.driver.importance.rescaled",
                          timevar="IPBES.direct.driver",idvar=c("UI","Indicator.name","Scale.weight", "Ndrivers.weight", "Old.indicator.reweight"),direction="wide")
} else {
  attplot.long <- attplot %>% select(UI,Indicator.name,EBV.class,IPBES.direct.driver,Direct.driver.importance.rescaled, Scale.weight, Ndrivers.weight, Old.indicator.reweight)
  attplot.wide <- reshape(attplot.long,v.names="Direct.driver.importance.rescaled",
                          timevar="IPBES.direct.driver",idvar=c("UI","Indicator.name","EBV.class","Scale.weight", "Ndrivers.weight", "Old.indicator.reweight"),direction="wide")
}

