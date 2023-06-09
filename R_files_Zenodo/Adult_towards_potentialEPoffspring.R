# feeding data
data <- read.csv("dataQ.csv",stringsAsFactors=FALSE, sep=",")


# chick info
chicks <- read.csv("dataChickID_toPublish20200223.csv",stringsAsFactors=FALSE,sep=",")


# parent info
parents <- read.csv("dataParentsID.csv",stringsAsFactors=FALSE,sep=",")


# FLEDGE DAY (TO FIX)
fledge_day <- 17

# only feeding data
#data <- data[data$fed == "y",]

# necessary feeding data
data <- data[,which(colnames(data) %in% c("Day", "Aviary", "AdultID", "fed", "Who", "Count","EventID", "IDchicks_AnyTime", "Num_NQ.anytime"))]
data$Day <- as.numeric(gsub("-","",data$Day))

# cohorts (1) 20180810 to 20180829, (2) 20180926 - 20181020
data$Cohort <- NA
data$Cohort[which(data$Day >= 20180810 & data$Day <= 20180829)] <- 1
data$Cohort[which(data$Day >= 20180926 & data$Day <= 20181020)] <- 2

# remove non-cohort data
data <- data[which(data$Cohort %in% c(1,2)), ]

# necessary chick data
#chicks <- chicks[which(is.na(chicks$deathdate)),]
chicks <- chicks[,which(colnames(chicks) %in% c("fosteraviary","Cohort","infidelity","geneticmother", "geneticfather", "Idfosterbrood", "backpack", "datebackpack", "hatchdate", "deathdate"))]
chicks$datebackpack <- as.numeric(gsub("-","",chicks$datebackpack))
chicks$hatchdate <- as.numeric(gsub("-","",chicks$hatchdate))
chicks$deathdate <- as.numeric(gsub("-","",chicks$deathdate))
chicks$fledgedate <- as.character(as.Date(as.character(chicks$hatchdate), format="%Y%m%d") + fledge_day)
chicks$fledgedate <- as.numeric(gsub("-","",chicks$fledgedate))

# add a date that has either the backpacking (if occurred) OR death (if not backpacked)
chicks$lastNQ <- chicks$datebackpack
chicks$lastNQ[is.na(chicks$lastNQ)] <- chicks$deathdate[is.na(chicks$lastNQ)]

# necessary adult data
parents <- parents[,which(colnames(parents) %in% c("ind_ID", "aviary", "qr_code","sex"))]

# only males
parents <- parents[which(parents$sex == "m"),]

# find offspring of parents
parents.offspring_c1 <- list()
parents.offspring_c2 <- list()
parents$first_backpack_c1 <- NA
parents$first_backpack_c2 <- NA
parents$all_backpack_c1 <- NA
parents$all_backpack_c2 <- NA
parents$first_c2_chick <- NA
for (i in 1:nrow(parents)) {
	parents.offspring_c1[[i]] <- as.numeric(chicks$backpack[which(chicks$geneticfather == parents$ind_ID[i] & chicks$Cohort == 1 & chicks$infidelity == 1)])
	parents.offspring_c2[[i]] <- as.numeric(chicks$backpack[which(chicks$geneticfather == parents$ind_ID[i] & chicks$Cohort == 2 & chicks$infidelity == 1)])
	if (sum(!is.na(parents.offspring_c1[[i]])) > 0) {
		# first add the other chicks from the same brood
		parents.offspring_c1[[i]] <- as.numeric(chicks$backpack[which(chicks$Idfosterbrood %in% chicks$Idfosterbrood[which(chicks$backpack %in% parents.offspring_c1[[i]])])])

		parents$first_backpack_c1[i] <- min(chicks$datebackpack[which(chicks$backpack %in% parents.offspring_c1[[i]])], na.rm=T)
		parents$all_backpack_c1[i] <- max(chicks$lastNQ[which(chicks$backpack %in% parents.offspring_c1[[i]])])
	}
	if (sum(!is.na(parents.offspring_c2[[i]])) > 0) {
		# first add the other chicks from the same brood
		parents.offspring_c2[[i]] <- as.numeric(chicks$backpack[which(chicks$Idfosterbrood %in% chicks$Idfosterbrood[which(chicks$backpack %in% parents.offspring_c2[[i]])])])

		parents$first_backpack_c2[i] <- min(chicks$datebackpack[which(chicks$backpack %in% parents.offspring_c2[[i]])], na.rm=T)
		parents$all_backpack_c2[i] <- max(chicks$lastNQ[which(chicks$backpack %in% parents.offspring_c2[[i]])])
		parents$first_c2_chick[i] <- min(chicks$fledgedate[which(chicks$backpack %in% parents.offspring_c2[[i]])], na.rm=T)
	}
}


# add metadata to feeding data
data$allbackpacked <- FALSE
data$anybackpacked <- FALSE
data$noc2chicks <- FALSE

data$offspring_present <- ""
data$non_offspring_present <- ""

data$offspring_fed <- ""
data$non_offspring_fed <- ""

data$n_offspring_present <- 0
data$n_non_offspring_present <- 0
data$n_unknown_present <- 0

data$offspring_fed_count <- 0
data$nonoffspring_fed_count <- 0
data$unknown_fed_count <- 0
data$offspring_unique_fed <- 0
data$nonoffspring_unique_fed <- 0
data$unknown_unique_fed <- 0

# remove data from females
data <- data[which(data$AdultID %in% parents$qr_code),]

# loop through non-empty feeds
for (i in 1:nrow(data)) {

	id <- which(parents$qr_code == data$AdultID[i])
	ch_fed <- as.numeric(strsplit(data$Who[i], " ")[[1]])
	counts <- as.numeric(strsplit(data$Count[i], " ")[[1]])
	ids_present <- as.numeric(strsplit(data$IDchicks_AnyTime[i], " ")[[1]])

	if (data$Cohort[i] == 1) {
		off <- parents.offspring_c1[[id]]
		if (!is.na(parents$all_backpack_c1[id])) {
		if (data$Day[i] >= parents$all_backpack_c1[id]) {
			data$allbackpacked[i] <- TRUE
		}}
		if (!is.na(parents$first_backpack_c1[id])) {
		if (data$Day[i] >= parents$first_backpack_c1[id]) {
			data$anybackpacked[i] <- TRUE
		}}
		if (!is.na(parents$first_c2_chick[id])) {
		#if (data$Day[i] >= parents$first_c2_chick[id]) {  # check if they ever had c2 chicks
			data$noc2chicks[i] <- TRUE
		}#}
	} else {
		off <- parents.offspring_c2[[id]]
		if (!is.na(parents$all_backpack_c2[id])) {
		if (data$Day[i] >= parents$all_backpack_c2[id]) {
			data$allbackpacked[i] <- TRUE
		}}
		if (!is.na(parents$first_backpack_c2[id])) {
		if (data$Day[i] >= parents$first_backpack_c2[id]) {
			data$anybackpacked[i] <- TRUE
		}}
	}

	if (length(ch_fed) > 0) {
	for (j in 1:length(ch_fed)) {
		if (is.na(ch_fed[j])) {
			data$unknown_fed_count[i] <- data$unknown_fed_count[i] + counts[j]
			data$unknown_unique_fed[i] <- data$unknown_unique_fed[i] + 1
		} else if (ch_fed[j] %in% off) {
			data$offspring_fed_count[i] <- data$offspring_fed_count[i] + counts[j]
			data$offspring_unique_fed[i] <- data$offspring_unique_fed[i] + 1
			data$offspring_fed[i] <- trimws(paste(data$offspring_fed[i],ch_fed[j],sep=" "))
		} else {
			data$nonoffspring_fed_count[i] <- data$nonoffspring_fed_count[i] + counts[j]
			data$nonoffspring_unique_fed[i] <- data$nonoffspring_unique_fed[i] + 1
			data$non_offspring_fed[i] <- trimws(paste(data$non_offspring_fed[i],ch_fed[j],sep=" "))
		}
	}
	}

	if (length(ids_present) > 0) {
	for (j in 1:length(ids_present)) {
		if (is.na(ids_present[j])) {
			data$n_unknown_present[i] <- data$n_unknown_present[i] + 1
		} else if (ids_present[j] %in% off) {
			data$n_offspring_present[i] <- data$n_offspring_present[i] + 1
			data$offspring_present[i] <- trimws(paste(data$offspring_present[i],ids_present[j],sep=" "))
		} else {
			data$n_non_offspring_present[i] <- data$n_non_offspring_present[i] + 1
			data$non_offspring_present[i] <- trimws(paste(data$non_offspring_present[i],ids_present[j],sep=" "))
		}
	}
	}

	data$n_unknown_present[i] <- data$n_unknown_present[i] + data$Num_NQ.anytime[i]


}



#Weird individuals
#data <- data[which(!(data$AdultID %in% c("1","3","14","144"))),]


# Summarise all data
nrow(data)#1917 (without weirdos) #2135 (with weirdos)
data <- data[which(data$EventID!="3334"),]#because $who is empty but $fed=="y"
nrow(data[which(data$n_offspring_present>0),])#51 (with weirdos)
unique(length(data$EventID[which(data$n_offspring_present>0)]))
nrow(data[which(data$n_offspring_present>0 & (data$offspring_unique_fed + data$nonoffspring_unique_fed+ data$unknown_unique_fed)>0),])#14 (with weirdos)
nrow(data[which(data$n_offspring_present>0 & (data$offspring_unique_fed + data$nonoffspring_unique_fed+ data$unknown_unique_fed)>0 & data$offspring_unique_fed>0),])#9 (with weirdos)



# number of male adults who had at least one backpacked potential extra pair offspring
nrow(parents) # number of all males #56
nrow(parents[!is.na(parents$first_backpack_c1) | !is.na(parents$first_backpack_c2),])#15


