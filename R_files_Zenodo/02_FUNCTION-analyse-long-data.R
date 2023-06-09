analyse.long.data <- function(longdata = attplot.long){

	### Relative impacts of drivers aggregated (mean or weighted mean) ACROSS the different studies included in the analysis 
	#   --> for each indicator (option "analysis.level" = indicators)
	#   --> for each indicator and ebv class (option "analysis.level" = indicators.x.ebv)
	#   --> for each indicator and then for each EBV class (option "analysis.level" = ebv)
	#   --> across all indicators irrespective of ebv class (option "analysis.level" = overall)
	
	number.of.rows <- nrow(longdata)
	number.of.different.studies <- length(unique(longdata[,1]))
	
	if(weighted.mean.procedure=="No"){ # arithmetic mean (irrespective of number of drivers assessed or spatial coverage)
	  if(analysis.level=="indicators" || analysis.level=="overall"){
		longdata <- aggregate(Direct.driver.importance.rescaled~Indicator.name+IPBES.direct.driver,data=longdata,mean)
		if(analysis.level=="overall"){
		  longdata <- aggregate(Direct.driver.importance.rescaled~IPBES.direct.driver,data=longdata,mean)      
		}
	  } else {
		longdata <- aggregate(Direct.driver.importance.rescaled~Indicator.name+EBV.class+IPBES.direct.driver,data=longdata,mean)
		if(analysis.level=="ebv"){
		  longdata <- aggregate(Direct.driver.importance.rescaled~EBV.class+IPBES.direct.driver,data=longdata,mean)      
		}
	  }
	} else { # weighted mean (importance of each study varies according to the number of drivers assessed and the spatial coverage)
	  longdata$Weight <- (longdata$Ndrivers.weight+longdata$Scale.weight)/2 # weight is the average of scale and ndriver weights, in this analysis
	  if(analysis.level=="indicators" || analysis.level=="overall"){  
		longdata <- ldply(split(longdata,list(longdata$Indicator.name,longdata$IPBES.direct.driver)), function(longdata) weighted.mean(longdata$Direct.driver.importance.rescaled,w=longdata$Weight))
		longdata[,3:4] <- data.frame(do.call('rbind',strsplit(as.character(longdata$.id),'.',fixed=TRUE)))
		longdata <- cbind(longdata[,3:4],longdata[,2])
		colnames(longdata) <- c("Indicator.name","IPBES.direct.driver","Direct.driver.importance.rescaled")
		if(analysis.level=="overall"){
		  longdata <- aggregate(Direct.driver.importance.rescaled~IPBES.direct.driver,data=longdata,mean)  # each indicator is assumed to have the same importance for the overall pattern
		}
	  } else {
		longdata <- ldply(split(longdata,list(longdata$Indicator.name,longdata$EBV.class,longdata$IPBES.direct.driver)), function(longdata) weighted.mean(longdata$Direct.driver.importance.rescaled,w=longdata$Weight))
		longdata[,3:5] <- data.frame(do.call('rbind',strsplit(as.character(longdata$.id),'.',fixed=TRUE)))
		longdata <- cbind(longdata[,3:5],longdata[,2])
		colnames(longdata) <- c("Indicator.name","EBV.class","IPBES.direct.driver","Direct.driver.importance.rescaled")
		if(analysis.level=="ebv"){
		  longdata <- aggregate(Direct.driver.importance.rescaled~EBV.class+IPBES.direct.driver,data=longdata,mean)  
		}  
	  }
	}

	# Rescale again the ranking of drivers (--> sum of ranks = 1) (needed after the averaging procedure applied above)

	if(analysis.level=="indicators"){
	  stats.importance <- longdata %>% group_by(Indicator.name) %>% summarise(Sum.importance=sum(Direct.driver.importance.rescaled,na.rm=T))
	  longdata <- merge(longdata,stats.importance,by.x=c("Indicator.name"),by.y=c("Indicator.name"),all.x=T,all.y=T)
	  longdata$Direct.driver.importance.rescaled <- longdata$Direct.driver.importance.rescaled/longdata$Sum.importance
	  longdata <- longdata[,-ncol(longdata)]
	  longdata <- merge(longdata,n.studies.per.indicator,by.x=c("Indicator.name"),by.y=c("Indicator.name"),all.x=T,all.y=T)
	} else if(analysis.level=="overall"){
	  stats.importance <- longdata %>% group_by() %>% summarise(Sum.importance=sum(Direct.driver.importance.rescaled,na.rm=T))
	  stats.importance <- as.numeric(stats.importance)
	  longdata$Direct.driver.importance.rescaled <- longdata$Direct.driver.importance.rescaled/stats.importance
	} else if(analysis.level=="indicators.x.ebv"){
	  stats.importance <- longdata %>% group_by(Indicator.name,EBV.class) %>% summarise(Sum.importance=sum(Direct.driver.importance.rescaled,na.rm=T))
	  longdata <- merge(longdata,stats.importance,by.x=c("Indicator.name","EBV.class"),by.y=c("Indicator.name","EBV.class"),all.x=T,all.y=T)
	  longdata$Direct.driver.importance.rescaled <- longdata$Direct.driver.importance.rescaled/longdata$Sum.importance
	  longdata <- longdata[,-ncol(longdata)]
	  longdata <- merge(longdata,n.studies.per.indicator.x.ebv,by.x=c("Indicator.name","EBV.class"),by.y=c("Indicator.name","EBV.class"),all.x=T,all.y=T)
	} else {
	  stats.importance <- longdata %>% group_by(EBV.class) %>% summarise(Sum.importance=sum(Direct.driver.importance.rescaled,na.rm=T))
	  longdata <- merge(longdata,stats.importance,by.x=c("EBV.class"),by.y=c("EBV.class"),all.x=T,all.y=T)
	  longdata$Direct.driver.importance.rescaled <- longdata$Direct.driver.importance.rescaled/longdata$Sum.importance
	  longdata <- longdata[,-ncol(longdata)]
	  longdata <- merge(longdata,n.studies.per.ebv,by.x=c("EBV.class"),by.y=c("EBV.class"),all.x=T,all.y=T)
	}

	longdata <- longdata[!is.na(longdata$Direct.driver.importance.rescaled),]
	# longdata: object 'longdata' documents the relative impacts of each driver 
	# --> on each indicator (option "analysis.level" = indicators)
	# --> on each indicator within each ebv class (option "analysis.level" = indicators.x.ebv)
	# --> on each EBV class (option "analysis.level" = ebv)
	# --> across all indicators and ebv classes (option "analysis.level" = overall)

	### Replace indicator acronyms with full names of indicators

	for (i in 1:nrow(indicators)) {  
	  if (is.na(indicators$Type.of.indicator[i])){ 
		indicators$Type.of.indicator[i] <- indicators$Type.of.indicator[i-1] 
	  } else {
		indicators$Type.of.indicator[i] <- indicators$Type.of.indicator[i]
	  }
	}

	levels(indicators$Simplified.name.of.indicator) <- c(levels(indicators$Simplified.name.of.indicator),"IUCN red list of\\nthreatened species","Proportion of fish stocks\\nwithin sustainable levels","Percentage of live\\ncoral cover","Mangrove forest cover","Extent of intact forest\\nlandscapes") 
	indicators$Simplified.name.of.indicator[indicators$Simplified.name.of.indicator=="IUCN red list of threatened species"]  <- "IUCN red list of\\nthreatened species"
	indicators$Simplified.name.of.indicator[indicators$Simplified.name.of.indicator=="Proportion of fish stocks within biologically sustainable levels"]  <- "Proportion of fish stocks\\nwithin sustainable levels"
	indicators$Simplified.name.of.indicator[indicators$Simplified.name.of.indicator=="Percentage of live coral cover"]  <- "Percentage of live\\ncoral cover"
	indicators$Simplified.name.of.indicator[indicators$Simplified.name.of.indicator=="Area of mangrove forest cover"]  <- "Mangrove forest cover"
	indicators$Simplified.name.of.indicator[indicators$Simplified.name.of.indicator=="Extent of intact forest landscapes"]  <- "Extent of intact forest\\nlandscapes"
	old.lvl <- levels(indicators$Simplified.name.of.indicator)
	indicators$Simplified.name.of.indicator <- factor(indicators$Simplified.name.of.indicator,levels=sort(levels(factor(indicators$Simplified.name.of.indicator)),decreasing=F))

	if(analysis.level=="indicators" || analysis.level=="indicators.x.ebv"){
	  longdata <- merge(longdata,indicators,by.x=c("Indicator.name"),by.y=c("Acronym.of.indicator"),all.x=T,all.y=F)
	  if(analysis.level=="indicators.x.ebv"){
		longdata$EBV.class <- longdata$EBV.class.x
		longdata <- subset(longdata,select=-c(EBV.class.x,EBV.class.y))    
	  }
	  longdata$Indicator.acronym <- longdata$Indicator.name
	  longdata$Indicator.name <- longdata$Simplified.name.of.indicator
	}

	### Export results with estimation of relative impacts of drivers

	fig.region <- ifelse(region.selection=="No","All regions",ifelse(selected.region=="All regions","All regions (cross-region assessments)",selected.region))
	fig.realm <- ifelse(realm.selection=="No","All realms",ifelse(selected.realm=="All realms","All realms (cross-realm assessments)",paste(selected.realm,"realm",sep=" ")))
	fig.level <- ifelse(analysis.level=="indicators","Indicators",ifelse(analysis.level=="ebv","EBV","Overall")) # Modified 2010-12-03 correcting export.level to analysis.level (AP)
	fig.method <- include.non.assessed.drivers
	fig.title <- paste(paste(fig.region,fig.realm,fig.level,fig.method,sep=" - "),sep=" ")

	if(analysis.level=="indicators"){
	  longdata$Indicator.name <- factor(longdata$Indicator.name,levels=levels(factor(longdata$Indicator.name)))
	} else if(analysis.level=="indicators.x.ebv"){
	  longdata$Indicator.name <- factor(longdata$Indicator.name,levels=levels(factor(longdata$Indicator.name)))
	  longdata$EBV.class <- factor(longdata$EBV.class,levels=levels(factor(att$EBV.class)))
	} else if(analysis.level=="ebv"){
	  longdata$EBV.class <- factor(longdata$EBV.class,levels=levels(factor(att$EBV.class)))
	}
	longdata$IPBES.direct.driver <- factor(longdata$IPBES.direct.driver,levels=levels(factor(att$IPBES.direct.driver)))
	
	
	#Prepare data frame with just the required columns
	if(analysis.level=="overall"){
	  longdata <- longdata[with(longdata,order(IPBES.direct.driver)),]
	  if(export.tables=="Yes"){
		to.return <- longdata
	  }
	} else if(analysis.level=="indicators"){
	  longdata <- longdata[with(longdata,order(Indicator.name,IPBES.direct.driver)),]
	  if(export.tables=="Yes"){
		to.return <- longdata[,c(1:4)]
	  }  
	} else if(analysis.level=="ebv"){
	  longdata <- longdata[with(longdata,order(IPBES.direct.driver,EBV.class)),]
	  if(export.tables=="Yes"){
		to.return <- longdata[,c(1:4)]
	  }        
	} else if(analysis.level=="indicators.x.ebv"){
	  longdata <- longdata[with(longdata,order(EBV.class,Indicator.name,IPBES.direct.driver)),]
	  if(export.tables=="Yes"){
		to.return <- longdata[,c(8,1:4)]
	  }        
	}

	if(export.tables=="Yes"){
		#Add more attributes to make sure useful metadata are associated with the data frame - analysis options and size of dataset
		attr(to.return, "analysis") <- "Average relative impact"
		attr(to.return, "number_of_rows") <- number.of.rows
		attr(to.return, "number_of_different_studies") <- number.of.different.studies
		attr(to.return, "include_non_assessed_drivers") <- include.non.assessed.drivers
		attr(to.return, "combine_rank_magnitude") <- combine.rank.magnitude
		attr(to.return, "analysis_level") <- analysis.level
		attr(to.return, "selected_indicators") <- ifelse(indicator.selection == "Yes", toString(selected.indicator), "all")
		attr(to.return, "excluded_indicators") <- ifelse(indicator.exclusion == "Yes", toString(excluded.indicator), "none")
		attr(to.return, "selected_realms") <- ifelse(realm.selection == "Yes", selected.realm, "All realms")
		attr(to.return, "selected_regions") <- ifelse(region.selection == "Yes", selected.region, "All regions")
		attr(to.return, "selected_domains") <- ifelse(domain.selection == "Yes", selected.domain, "All domains")
		attr(to.return, "selected_scales") <- ifelse(scale.selection == "Yes", selected.scale, "All scales")
		attr(to.return, "selected_no_of_drivers") <- ifelse(n.drivers.selection == "Yes", toString(selected.n.drivers), toString(c(2:6)))
		attr(to.return, "manual_selection") <- manual.selection
		attr(to.return, "weighted_mean_procedure") <- weighted.mean.procedure
		return(to.return)
		}
}