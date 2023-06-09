# Analysis of Summer 2013 EID data
# William E. Stutz
# January 13, 2016
# R version 3.0.2 "Frisbee Sailing"


###############
# Preliminaries
###############
{
	## upload the full raw data
	dat <- read.csv("data/raw/coinfection_data.csv")

	## libraries
	library(ggplot2) 	# for plotting
	library(Cairo)		# for plotting
	library(reshape)	# adjusting data frames
	library(plyr)		# apply functions to data frames
	library(MCMCglmm)	# mixed models
	library(coda)		# calculate HPD intervals for MCMC simulations
	library(gridExtra)	# put multiple plots in one figure
	library(grid)

	## options
	options(scipen = 9) # remove scientific notation

	## functions
	source("scripts/R_functions/analysis_functions.R")

	## set Cairo Fonts
	CairoFonts(
		regular = "Linux Biolinum O",
		italic = "Linux Biolinum O:style=Italic",
		bold = "Linux Biolinum O:style=Bold",
		bolditalic = "Linux Biolinum O"
	)

	## open plotting window
	#Cairo(type = "X11")
}


###########################################################
# create full individual level data set to use for modeling
###########################################################
{

	## deal with Ranaviral data first
	{
		# multiply ranaviral loads by correction factor (see metadata)
		dat$RVLOAD <- dat$RVLOAD*10.52

		# look at ranavirus outliers
		tail(sort(dat$RVLOAD))

		# remove two extreme ranavirus outliers
		dat[dat$RVLOAD > 10000 & is.na(dat$RVLOAD) == FALSE,c("RVLOAD")] <- NA
	}

	## reorder species levels
	dat$HOSTSPECIES <- factor(dat$HOSTSPECIES, levels = c("PSRE","BUBO","RACA","TATO","TAGR"))

	## remove adults and sub-adults
	dat <- dat[dat$STAGE != "ADULT" & dat$STAGE != "Subadult" & dat$STAGE != "Adult" & dat$STAGE != "0" & dat$STAGE != " ",]
	dat$STAGE <- factor(dat$STAGE)


	## rename Rib load
	dat$RIBLOAD <- dat$RIB

	## rename BDZE to BDLOAD
	dat$BDLOAD <- dat$BDZE

	## rename larval trematode abundance
	dat$TREMLOAD <- dat$LARVTREMTOT

	## rename echinoload
	dat$ECHINOLOAD <- dat$ECHINO

	## add postive status
	dat$RIBPOS <- as.numeric(as.logical(dat$RIB > 0))
	dat$TREMPOS <- as.numeric(as.logical(dat$TREMLOAD > 0))
	dat$ECHINOPOS <- as.numeric(as.logical(dat$ECHINOLOAD > 0))

	## round loads to nearest whole number
	{
		## round BDLOAD
		dat[dat$BDLOAD < 1 & is.na(dat$BDLOAD) == FALSE,"BDLOAD"] <- ceiling(dat[dat$BDLOAD < 1 & is.na(dat$BDLOAD) == FALSE,"BDLOAD"])
		dat$BDLOAD <- round(dat$BDLOAD)

		## round RVLOAD
		dat[dat$RVLOAD < 1 & is.na(dat$RVLOAD) == FALSE,"RVLOAD"] <-ceiling(dat[dat$RVLOAD < 1 & is.na(dat$RVLOAD) == FALSE,"RVLOAD"])
		dat$RVLOAD <- round(dat$RVLOAD)

	}

	## standardize within host species
	dat <- ddply(dat,.(HOSTSPECIES),scale.species)


}

###################################
# create site level data to use
###################################
{

	## create species specific data frames
	psre <- dat[dat$HOSTSPECIES == "PSRE",]
	bubo <- dat[dat$HOSTSPECIES == "BUBO",]
	raca <- dat[dat$HOSTSPECIES == "RACA",]
	tato <- dat[dat$HOSTSPECIES == "TATO",]
	tagr <- dat[dat$HOSTSPECIES == "TAGR",]

	## function to calculate number of animals sampled for each parasite
	length.na.rm <- function(x){x <- x[is.na(x) == FALSE]; length(x)}

	## make site data frame
	{
		site <- data.frame(SiteCode = levels(dat$SiteCode),
			RIBPREZ = as.numeric(as.logical(tapply(dat$RIBPOS,dat$SiteCode,mean, na.rm = TRUE))),
			BDPREZ = as.numeric(as.logical(tapply(dat$BDPOS,dat$SiteCode,mean, na.rm = TRUE))),
			RVPREZ = as.numeric(as.logical(tapply(dat$RVPOS,dat$SiteCode,mean, na.rm = TRUE))),
			TREMPREZ = as.numeric(as.logical(tapply(dat$TREMPOS,dat$SiteCode,mean, na.rm = TRUE))),
			ECHINOPREZ = as.numeric(as.logical(tapply(dat$ECHINOPOS,dat$SiteCode,mean, na.rm = TRUE))),
			Ribtotal.count = tapply(dat$RIBPOS,dat$SiteCode,length.na.rm),
			BDtotal.count = tapply(dat$BDPOS,dat$SiteCode,length.na.rm),
			RVtotal.count = tapply(dat$RVPOS,dat$SiteCode,length.na.rm),
			Tremtotal.count = tapply(dat$TREMPOS,dat$SiteCode,length.na.rm),
			Echinototal.count = tapply(dat$ECHINOPOS,dat$SiteCode,length.na.rm),
			Ribtotal.ninf = table(dat$SiteCode,dat$RIBPOS)[,2],
			Bdtotal.ninf = table(dat$SiteCode,dat$BDPOS)[,2],
			RVtotal.ninf = table(dat$SiteCode,dat$RVPOS)[,2],
			tremtotal.ninf = table(dat$SiteCode,dat$TREMPOS)[,2],
			echinototal.ninf = table(dat$SiteCode,dat$ECHINOPOS)[,2],
			PSRE.ribprez = as.numeric(as.logical(tapply(psre$RIBPOS,psre$SiteCode,mean, na.rm = TRUE))),
			PSRE.Bdprez = as.numeric(as.logical(tapply(psre$BDPOS,psre$SiteCode,mean, na.rm = TRUE))),
			PSRE.Rvprez = as.numeric(as.logical(tapply(psre$RVPOS,psre$SiteCode,mean, na.rm = TRUE))),
			PSRE.tremprez = as.numeric(as.logical(tapply(psre$TREMPOS,psre$SiteCode,mean, na.rm = TRUE))),
			PSRE.echinoprez = as.numeric(as.logical(tapply(psre$ECHINOPOS,psre$SiteCode,mean, na.rm = TRUE))),
			PSRE.ribcount = tapply(psre$RIBPOS,psre$SiteCode,length.na.rm),
			PSRE.bdcount = tapply(psre$BDPOS,psre$SiteCode,length.na.rm),
			PSRE.rvcount = tapply(psre$RVPOS,psre$SiteCode,length.na.rm),
			PSRE.tremcount = tapply(psre$TREMPOS,psre$SiteCode,length.na.rm),
			PSRE.echinocount = tapply(psre$ECHINOPOS,psre$SiteCode,length.na.rm),
			BUBO.ribprez = as.numeric(as.logical(tapply(bubo$RIBPOS,bubo$SiteCode,mean, na.rm = TRUE))),
			BUBO.Bdprez = as.numeric(as.logical(tapply(bubo$BDPOS,bubo$SiteCode,mean, na.rm = TRUE))),
			BUBO.Rvprez = as.numeric(as.logical(tapply(bubo$RVPOS,bubo$SiteCode,mean, na.rm = TRUE))),
			BUBO.tremprez = as.numeric(as.logical(tapply(bubo$TREMPOS,bubo$SiteCode,mean, na.rm = TRUE))),
			BUBO.echinoprez = as.numeric(as.logical(tapply(psre$ECHINOPOS,psre$SiteCode,mean, na.rm = TRUE))),
			BUBO.ribcount = tapply(bubo$RIBPOS,bubo$SiteCode,length.na.rm),
			BUBO.bdcount = tapply(bubo$BDPOS,bubo$SiteCode,length.na.rm),
			BUBO.rvcount = tapply(bubo$RVPOS,bubo$SiteCode,length.na.rm),
			BUBO.tremcount = tapply(bubo$TREMPOS,bubo$SiteCode,length.na.rm),
			BUBO.echinocount = tapply(bubo$ECHINOPOS,bubo$SiteCode,length.na.rm),
			RACA.ribprez = as.numeric(as.logical(tapply(raca$RIBPOS,raca$SiteCode,mean, na.rm = TRUE))),
			RACA.Bdprez = as.numeric(as.logical(tapply(raca$BDPOS,raca$SiteCode,mean, na.rm = TRUE))),
			RACA.Rvprez = as.numeric(as.logical(tapply(raca$RVPOS,raca$SiteCode,mean, na.rm = TRUE))),
			RACA.tremprez = as.numeric(as.logical(tapply(raca$TREMPOS,raca$SiteCode,mean, na.rm = TRUE))),
			RACA.echinoprez = as.numeric(as.logical(tapply(raca$ECHINOPOS,raca$SiteCode,mean, na.rm = TRUE))),
			RACA.ribcount = tapply(raca$RIBPOS,raca$SiteCode,length.na.rm),
			RACA.bdcount = tapply(raca$BDPOS,raca$SiteCode,length.na.rm),
			RACA.rvcount = tapply(raca$RVPOS,raca$SiteCode,length.na.rm),
			RACA.tremcount = tapply(raca$TREMPOS,raca$SiteCode,length.na.rm),
			RACA.echinocount = tapply(raca$ECHINOPOS,raca$SiteCode,length.na.rm),
			TATO.ribprez = as.numeric(as.logical(tapply(tato$RIBPOS,tato$SiteCode,mean, na.rm = TRUE))),
			TATO.Bdprez = as.numeric(as.logical(tapply(tato$BDPOS,tato$SiteCode,mean, na.rm = TRUE))),
			TATO.Rvprez = as.numeric(as.logical(tapply(tato$RVPOS,tato$SiteCode,mean, na.rm = TRUE))),
			TATO.tremprez = as.numeric(as.logical(tapply(tato$TREMPOS,tato$SiteCode,mean, na.rm = TRUE))),
			TATO.echinoprez = as.numeric(as.logical(tapply(tato$ECHINOPOS,tato$SiteCode,mean, na.rm = TRUE))),
			TATO.ribcount = tapply(tato$RIBPOS,tato$SiteCode,length.na.rm),
			TATO.bdcount = tapply(tato$BDPOS,tato$SiteCode,length.na.rm),
			TATO.rvcount = tapply(tato$RVPOS,tato$SiteCode,length.na.rm),
			TATO.tremcount = tapply(tato$TREMPOS,tato$SiteCode,length.na.rm),
			TATO.echinocount = tapply(tato$ECHINOPOS,tato$SiteCode,length.na.rm),
			TAGR.ribprez = as.numeric(as.logical(tapply(tagr$RIBPOS,tagr$SiteCode,mean, na.rm = TRUE))),
			TAGR.Bdprez = as.numeric(as.logical(tapply(tagr$BDPOS,tagr$SiteCode,mean, na.rm = TRUE))),
			TAGR.Rvprez = as.numeric(as.logical(tapply(tagr$RVPOS,tagr$SiteCode,mean, na.rm = TRUE))),
			TAGR.tremprez = as.numeric(as.logical(tapply(tagr$TREMPOS,tagr$SiteCode,mean, na.rm = TRUE))),
			TAGR.echinoprez = as.numeric(as.logical(tapply(tagr$ECHINOPOS,tagr$SiteCode,mean, na.rm = TRUE))),
			TAGR.ribcount = tapply(tagr$RIBPOS,tagr$SiteCode,length.na.rm),
			TAGR.bdcount = tapply(tagr$BDPOS,tagr$SiteCode,length.na.rm),
			TAGR.rvcount = tapply(tagr$RVPOS,tagr$SiteCode,length.na.rm),
			TAGR.tremcount = tapply(tagr$TREMPOS,tagr$SiteCode,length.na.rm),
			TAGR.echinocount = tapply(tagr$ECHINOPOS,tagr$SiteCode,length.na.rm)
		)

	}

	## merge in metacommunity
	site <- merge(site, dat[,c("SiteCode","meta")])
	site <- unique(site,margin = 1)

	## change metacommunity column name
	colnames(site)[colnames(site) == "meta"] <- "Meta"

}

#################################
# create a table of sampled sites
#################################
{

	## upload data file with latitude and longitude data
	geo_data <- read.csv("data/raw/master_site_list.csv")

	## change site name to a character rather than factor
	geo_data$SiteCode_1 <- as.character(geo_data$SiteCode_1)

	## alter site codes in geo dataa to match our data frame
	{
		geo_data[geo_data$SiteCode == "CA-SHROM2","SiteCode_1"] <- "Ca-Step"
		geo_data[geo_data$SiteCode == "BassGCP","SiteCode_1"] <- "BASSGCP"
		geo_data[geo_data$SiteCode == "CA-SF30","SiteCode_1"] <- "Ca-SF30"
		geo_data[geo_data$SiteCode == "SF-40","SiteCode_1"] <- "SF40"
		geo_data[geo_data$SiteCode == "SF-41","SiteCode_1"] <- "SF41"
		geo_data[geo_data$SiteCode == "SF-42","SiteCode_1"] <- "SF42"
		geo_data[geo_data$SiteCode == "SF27","SiteCode_1"] <- "CA-SF27"
		geo_data[geo_data$SiteCode == "SF-DAM","SiteCode_1"] <- "SFDAM"
		geo_data[geo_data$SiteCode == "CA-SF85a","SiteCode_1"] <- "CA-SF85A"
		geo_data[geo_data$SiteCode == "HeronGCP","SiteCode_1"] <- "HERONGCP"
		geo_data[geo_data$SiteCode == "Barn","SiteCode_1"] <- "BARN"
		geo_data[geo_data$SiteCode == "EagleGCP","SiteCode_1"] <- "EAGLEGCP"
		geo_data[geo_data$SiteCode == "Frog","SiteCode_1"] <- "FROG"
		geo_data[geo_data$SiteCode == "Hidden","SiteCode_1"] <- "HIDDEN"
		geo_data[geo_data$SiteCode == "WestWing","SiteCode_1"] <- "WESTWING"
	}

	## keep only the data I want
	geo_data <- geo_data[,c("SiteCode_1","Latitude_1","Longitude_1")]
	names(geo_data) <- c("SiteCode","Latitude","Longitude")
	geo_data <- unique(geo_data)

	## merge with sampled site data
	site <- merge(site,geo_data, by = "SiteCode", all.x = "TRUE")

	## add host species totals
	species.table <- table(dat$SiteCode,dat$HOSTSPECIES)
	site <- data.frame(site, "PSRE" = species.table[,1], "BUBO" = species.table[,2], "RACA" = species.table[,3], "TATO" = species.table[,4], "TAGR" = species.table[,5])

	## table to output
	output <- site[,c("SiteCode","Longitude","Latitude","PSRE","BUBO","RACA","TATO","TAGR")]

	## write to file
	write.csv(output, file = "docs/tables/site_table.csv", row.names = FALSE)

}

########################################
# create a table of sampled host species
########################################
{


	## table to output
	output <- table(dat$HOSTSPECIES)

	## write to file
	write.csv(output, file = "docs/tables/host_species_table.csv", row.names = FALSE)

}

####################################
# Basic Statistics on Site Sampling
####################################
{

	## how many sites overall?
	use <- occ.data(site,
		nrib = 0, 	# pulls all sites with this minimum number of animals dissected for rib
		ntrem = 0,  # pulls all sites with this minimum number of animals dissected for trematodes
		nbd = 0, 	# pulls all sites with this minimum number of animals scanned for bd
		nrv = 0,		# pulls all sites with this minimum number of animals scanned for ranavirus
		nechino = 0, # pulls all sites with this minimum number of animals scanned for echinostomes
		silveroaks = "yes",	# include the three silver oaks sites?
		overall = "yes"		# data at site level?
	); nrow(use) # 93 sites

	## how many total individuals with at least one parasite sampled
	nrow(dat[is.na(dat$RIBPOS) == FALSE | is.na(dat$BDPOS) == FALSE | is.na(dat$RVPOS) == FALSE | is.na(dat$ECHINOPOS) == FALSE,]) #2152

	## how many total individuals with RIB,BD, and Ranavirus sampled
	nrow(dat[is.na(dat$RIBPOS) == FALSE & is.na(dat$BDPOS) == FALSE & is.na(dat$RVPOS) == FALSE,]) #1521

	## how many total individuals with all four sampled
	nrow(dat[is.na(dat$RIBPOS) == FALSE & is.na(dat$BDPOS) == FALSE & is.na(dat$RVPOS) == FALSE & is.na(dat$ECHINOPOS) == FALSE,]) #1521

	## how many NOT dissected for each parasite
	summary(dat$RIBPOS)	# 384
	summary(dat$BDPOS)	# 61
	summary(dat$RVPOS)	# 594
	summary(dat$ECHINOPOS)	# 384

	## how many NOT dissected for each parasite by species
	table(dat$HOSTSPECIES,dat$RIBPOS, exclude = NULL)
	table(dat$HOSTSPECIES,dat$BDPOS, exclude = NULL)
	table(dat$HOSTSPECIES,dat$RVPOS, exclude = NULL)
	table(dat$HOSTSPECIES,dat$ECHINOPOS, exclude = NULL)

	## how many collections overall?
	use <- occ.data(site,
		nrib = 0, 	# pulls all sites with this minimum number of animals dissected for rib
		ntrem = 0,  # pulls all sites with this minimum number of animals dissected for trematodes
		nbd = 0, 	# pulls all sites with this minimum number of animals scanned for bd
		nrv = 0,		# pulls all sites with this minimum number of animals scanned for ranavirus
		nechino = 0, # pulls all sites with this minimum number of animals scanned for echinostomes
		silveroaks = "yes",	# include the three silver oaks sites?
		overall = "no"		# data at site level?
	); nrow(use) # 465 collections


	## how many sites had at least one individual sampled for rib
	use <- occ.data(site,
		nrib = 1, 	# pulls all sites with this minimum number of animals dissected for rib
		ntrem = 0,  # pulls all sites with this minimum number of animals dissected for trematodes
		nbd = 0, 	# pulls all sites with this minimum number of animals scanned for bd
		nrv = 0,		# pulls all sites with this minimum number of animals scanned for ranavirus
		nechino = 0, # pulls all sites with this minimum number of animals scanned for echinostomes
		silveroaks = "yes",	# include the three silver oaks sites?
		overall = "yes"		# data at site level?
	); nrow(use) # 91 sites

	## how many sites were sampled for bd
	use <- occ.data(site,
		nrib = 0, 	# pulls all sites with this minimum number of animals dissected for rib
		ntrem = 0,  # pulls all sites with this minimum number of animals dissected for trematodes
		nbd = 1, 	# pulls all sites with this minimum number of animals scanned for bd
		nrv = 0,		# pulls all sites with this minimum number of animals scanned for ranavirus
		nechino = 0, # pulls all sites with this minimum number of animals scanned for echinostomes
		silveroaks = "yes",	# include the three silver oaks sites?
		overall = "yes"		# data at site level?
	); nrow(use) # 90 sites


	## how many sites were sampled for ranavirus
	use <- occ.data(site,
		nrib = 0, 	# pulls all sites with this minimum number of animals dissected for rib
		ntrem = 0,  # pulls all sites with this minimum number of animals dissected for trematodes
		nbd = 0, 	# pulls all sites with this minimum number of animals scanned for bd
		nrv = 1,		# pulls all sites with this minimum number of animals scanned for ranavirus
		nechino = 0, # pulls all sites with this minimum number of animals scanned for echinostomes
		silveroaks = "yes",	# include the three silver oaks sites?
		overall = "yes"		# data at site level?
	); nrow(use) # 91 sites

	## how many sites were sampled for ranavirus
	use <- occ.data(site,
		nrib = 0, 	# pulls all sites with this minimum number of animals dissected for rib
		ntrem = 0,  # pulls all sites with this minimum number of animals dissected for trematodes
		nbd = 0, 	# pulls all sites with this minimum number of animals scanned for bd
		nrv = 0,		# pulls all sites with this minimum number of animals scanned for ranavirus
		nechino = 1, # pulls all sites with this minimum number of animals scanned for echinostomes
		silveroaks = "yes",	# include the three silver oaks sites?
		overall = "yes"		# data at site level?
	); nrow(use) # 91 sites



	## how many sites were sampled for all four
	use <- occ.data(site,
		nrib = 1, 	# pulls all sites with this minimum number of animals dissected for rib
		ntrem = 0,  # pulls all sites with this minimum number of animals dissected for trematodes
		nbd = 1, 	# pulls all sites with this minimum number of animals scanned for bd
		nrv = 1,		# pulls all sites with this minimum number of animals scanned for ranavirus
		nechino = 1, # pulls all sites with this minimum number of animals scanned for echinostomes
		silveroaks = "yes",	# include the three silver oaks sites?
		overall = "yes"		# data at site level?
	); nrow(use) # 88 sites


	## what stages were collected
	table(dat$STAGE,dat$HOSTSPECIES)

	## how many of each species at each site
	table(dat$SiteCode,dat$HOSTSPECIES)




}

#########################################
# Basic Statistics on Individual Sampling
#########################################
{
	## extract individuals where we checked for all three parasites
	{ # 93 sites, 2152 individuals
		use <- coinf.data(dat,site,
			check.rib = 0,			# sites where we checked for Rib in at least this many animals
			found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
			check.bd = 0,			# sites where we checked for Bd in at least this many animals
			found.bd = 0,			# sites where we found Bd in at least this many animals
			check.rv = 0,		,	# sites where we checked for Ranavirus in at least this many animals
			found.rv = 0,			# sites where we found Ranavirus in at least this many animals
			check.echino = 0,		# sites where we checked for echino in at least this many animals
			found.echino = 0,		# sites where we found echino in at least this many animals
			use.NAs = "yes"		# keep individuals if some parasites were not checked
		)
	}

	# calculate observed levels of coinfection at the individual level
	{

		# real data
		obs.dat <- use[,c("HOSTSPECIES","SiteCode","RIBPOS","BDPOS","RVPOS","ECHINOPOS")]


		# combine test statistics into a single table
		{
			real <- data.frame(
				rib.echino.ind = nrow(obs.dat[obs.dat$RIBPOS == 1 & obs.dat$ECHINOPOS == 1 & is.na(obs.dat$RIBPOS) == FALSE & is.na(obs.dat$ECHINOPOS) == FALSE,]),
				rib.bd.ind = nrow(obs.dat[obs.dat$RIBPOS == 1 & obs.dat$BDPOS == 1 & is.na(obs.dat$RIBPOS) == FALSE & is.na(obs.dat$BDPOS) == FALSE,]),
				rib.rv.ind = nrow(obs.dat[obs.dat$RVPOS == 1 & obs.dat$RIBPOS == 1 & is.na(obs.dat$RIBPOS) == FALSE & is.na(obs.dat$RVPOS) == FALSE,]),
				bd.echino.ind = nrow(obs.dat[obs.dat$BDPOS == 1 & obs.dat$ECHINOPOS == 1 & is.na(obs.dat$BDPOS) == FALSE & is.na(obs.dat$ECHINOPOS) == FALSE,]),
				rv.echino.ind = nrow(obs.dat[obs.dat$RVPOS == 1 & obs.dat$ECHINOPOS == 1 & is.na(obs.dat$RVPOS) == FALSE & is.na(obs.dat$ECHINOPOS) == FALSE,]),
				bd.rv.ind = nrow(obs.dat[obs.dat$BDPOS == 1 & obs.dat$RVPOS == 1 & is.na(obs.dat$BDPOS) == FALSE & is.na(obs.dat$RVPOS) == FALSE,]),
				rib.echino.total = nrow(obs.dat[is.na(obs.dat$RIBPOS) == FALSE & is.na(obs.dat$ECHINOPOS) == FALSE,]),
				rib.bd.total = nrow(obs.dat[is.na(obs.dat$RIBPOS) == FALSE & is.na(obs.dat$BDPOS) == FALSE,]),
				rib.rv.total = nrow(obs.dat[is.na(obs.dat$RIBPOS) == FALSE & is.na(obs.dat$RVPOS) == FALSE,]),
				bd.rv.total = nrow(obs.dat[is.na(obs.dat$BDPOS) == FALSE & is.na(obs.dat$RVPOS) == FALSE,]),
				bd.echino.total = nrow(obs.dat[is.na(obs.dat$BDPOS) == FALSE & is.na(obs.dat$ECHINOPOS) == FALSE,]),
				rv.echino.total = nrow(obs.dat[is.na(obs.dat$RVPOS) == FALSE & is.na(obs.dat$ECHINOPOS) == FALSE,]))
		}
	}
}

############################################
# relationships between length and abundance
############################################
{

	## extract individuals where we checked for all three parasites
	{ # 93 sites, 2152 individuals
		use <- coinf.data(dat,site,
			check.rib = 0,			# sites where we checked for Rib in at least this many animals
			found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
			check.bd = 0,			# sites where we checked for Bd in at least this many animals
			found.bd = 0,			# sites where we found Bd in at least this many animals
			check.rv = 0,		,	# sites where we checked for Ranavirus in at least this many animals
			found.rv = 0,			# sites where we found Ranavirus in at least this many animals
			check.echino = 0,		# sites where we checked for echino in at least this many animals
			found.echino = 0,		# sites where we found echino in at least this many animals
			use.NAs = "yes"		# keep individuals if some parasites were not checked
		)
	}

	## pull data
	length.data <- use[,c("HOSTSPECIES","RIBLOAD","BDLOAD","RVLOAD","ECHINOLOAD","svl_adj")]

	## melt
	length.data <- melt(length.data,measure.vars = c("RIBLOAD","BDLOAD","RVLOAD","ECHINOLOAD"))

	## plot
	Cairo(width = 9, height = 6, type = "svg", units = "in", file = "output/length")
	ggplot(data = length.data, aes(svl_adj,log(value + 1))) +
		geom_point(pch = 21) +
		stat_smooth(method = "glm", family = "gaussian") +
		labs(x = "\\nstandardized snout-vent length",y = "log(parasite load + 1)\\n") +
		facet_grid(variable~HOSTSPECIES, scales = "free") +
		theme(strip.background = element_blank(),
					axis.title = element_text(size = 20),
					axis.text = element_text(size = 14, color = "black"),
					panel.background = element_rect(color = "black",fill = "white"),
					panel.grid = element_blank())
	dev.off()

}

##################################
# single-response abundance models
##################################
{

	## Fit Ribeiroia distributions
	{

		## extract data from sites where we found
		{ # 58 sites, 1303 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 1,			# sites where we checked for Rib in at least this many animals
				found.rib = 1,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 0,			# sites where we checked for Bd in at least this many animals
				found.bd = 0,			# sites where we found Bd in at least this many animals
				check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 0,			# sites where we found Ranavirus in at least this many animals
				check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
				found.echino = 0		# sites where we found Echinostomes in at least this many animals

			)
		}

		### fit poisson, zero-inflated, and hurdele distributions to entire data set using site intercepts
		dist.test(use, parasite = "RIBLOAD", nitt = c(2000000,20000000,10), burnin = c(1000000,10000000,1), thin = c(1000,10000,1), snails = "no")

	 	### upload models
		for(i in c("poismod","zimod","humod")){
			load(paste0("data/results/models/abundance/RIBLOAD_",i,".RData"), verbose = TRUE)
		}

		### compare autocorrelation
		acorr(pois.mod,summ.zp(pois.mod)$cstats[3])
		acorr(zi.mod,summ.zp(zi.mod)$cstats[3])

		### what is the expected number of zeros?
		{
			ppcheck.zeros(use,pois.mod,parasite = "RIBLOAD")
			ppcheck.zeros(use,zi.mod,parasite = "RIBLOAD",prob.model = "zero-inflated")
		}

		### simulate posteriors of predicted species and site means
		sim_rib_means <- simulate_abundances(sub = use, mod = pois.mod); save(sim_rib_means, file = "output/sim_rib_means.RData")

		### pull posteriors of site and species effects
		rib_effects <- simulate_abundance_effects(sub = use, mod = pois.mod); save(rib_effects, file = "output/rib_effects.RData")

	}

	## Fit Bd distributions
	{
		## extract data from sites where we found
		{ # 52 sites, 1448 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 0,			# sites where we checked for Rib in at least this many animals
				found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 1,			# sites where we checked for Bd in at least this many animals
				found.bd = 1,			# sites where we found Bd in at least this many animals
				check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 0,			# sites where we found Ranavirus in at least this many animals
				check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
				found.echino = 0		# sites where we found Echinostomes in at least this many animals

			)
		}

		### fit poisson, zero-inflated, and hurdele distributions to entire data set using site intercepts
		#dist.test(use, parasite = "BDLOAD", nitt = c(2000000,20000000,10), burnin = c(1000000,10000000,1), thin = c(1000,10000,1), snails = "no")

		### upload models
		for(i in c("poismod","zimod","humod")){
				load(paste0("data/results/models/abundance/BDLOAD_",i,".RData"), verbose = TRUE)
		}


		### compare autocorrelations
		acorr(pois.mod,summ.zp(pois.mod)$cstats[3])
		acorr(zi.mod,summ.zp(zi.mod)$cstats[3])

		### what is the expected number of zeros?
		{
			ppcheck.zeros(use,pois.mod,parasite = "BDLOAD")
			ppcheck.zeros(use,zi.mod,parasite = "BDLOAD",prob.model = "zero-inflated")
		}

		### simulate posteriors of predicted species and site means
		sim_bd_means <- simulate_abundances(sub = use, mod = pois.mod);  save(sim_bd_means, file = "output/sim_bd_means.RData"); sum(is.na(rowSums(sim_bd_means)))

		### pull posteriors of site and species effects
		bd_effects <- simulate_abundance_effects(sub = use, mod = pois.mod); save(bd_effects, file = "output/bd_effects.RData")

	}

	## Fit Ranavirus distribution
	{

		## extract data from sites where we found
		{ # 60 sites, 1206 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 0,			# sites where we checked for Rib in at least this many animals
				found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 0,			# sites where we checked for Bd in at least this many animals
				found.bd = 0,			# sites where we found Bd in at least this many animals
				check.rv = 1,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 1,			# sites where we found Ranavirus in at least this many animals
				check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
				found.echino = 0		# sites where we found Echinostomes in at least this many animals
			)
		}

		### fit poisson, zero-inflated, and hurdele distributions to entire data set using site intercepts
		dist.test(use, parasite = "RVLOAD", nitt = c(2000000,40000000,10), burnin = c(1000000,20000000,1), thin = c(1000,20000,1), snails = "no")

		### upload models
		for(i in c("poismod","zimod","humod")){
				load(paste0("data/results/models/abundance/RVLOAD_",i,".RData"), verbose = TRUE)
		}

		### compare autocorrelations
		acorr(pois.mod,summ.zp(pois.mod)$cstats[3])
		acorr(zi.mod,summ.zp(zi.mod)$cstats[3])

		### what is the expected number of zeros?
		{
			ppcheck.zeros(use,pois.mod,parasite = "RVLOAD")
			ppcheck.zeros(use,zi.mod,parasite = "RVLOAD", prob.model = "zero-inflated")
		}

		### simulate abundances
		sim_rv_means <- simulate_abundances(sub = use, mod = pois.mod);  save(sim_rv_means, file = "output/sim_rv_means.RData")

		### pull posteriors of site and species effects
		rv_effects <- simulate_abundance_effects(sub = use, mod = pois.mod); save(rv_effects, file = "output/rv_effects.RData")

	}

	## Fit Trematode distributions
	{

		## extract data from sites where we found
		{ # 82 sites, 1680 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 0,			# sites where we checked for Rib in at least this many animals
				found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 0,			# sites where we checked for Bd in at least this many animals
				found.bd = 0,			# sites where we found Bd in at least this many animals
				check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 0,			# sites where we found Ranavirus in at least this many animals
				check.trem = 1,
				found.trem = 1
			)
		}

		### fit poisson, zero-inflated, and hurdele distributions to entire data set using site intercepts
		dist.test(use, parasite = "TREMLOAD", nitt = c(2000000,20000000,10), burnin = c(1000000,10000000,1), thin = c(1000,10000,1), snails = "no")

		### upload models
		for(i in c("poismod","zimod","humod")){
				load(paste0("data/results/models/abundance/TREMLOAD_",i,".RData"), verbose = TRUE)
		}

		### compare autocorrelations
		acorr(pois.mod,summ.zp(pois.mod)$cstats[3])
		acorr(zi.mod,summ.zp(zi.mod)$cstats[3])

		### what is the expected number of zeros?
		{
			ppcheck.zeros(use,pois.mod,parasite = "TREMLOAD")
			ppcheck.zeros(use,zi.mod,parasite = "TREMLOAD", prob.model = "zero-inflated")
		}

		### simulate abundances
		sim_trem_means <- simulate_abundances(sub = use, mod = pois.mod);  save(sim_trem_means, file = "output/sim_trem_means.RData")

		### pull posteriors of site and species effects
		trem_effects <- simulate_abundance_effects(sub = use, mod = pois.mod); save(trem_effects, file = "output/trem_effects.RData")



	}

	## Fit Echino distributions
	{

		## extract data from sites where we found
		{ # 78 sites, 1639 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 0,			# sites where we checked for Rib in at least this many animals
				found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 0,			# sites where we checked for Bd in at least this many animals
				found.bd = 0,			# sites where we found Bd in at least this many animals
				check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 0,			# sites where we found Ranavirus in at least this many animals
				check.trem = 0,
				found.trem = 0,
				check.echino = 1,
				found.echino = 1

			)
		}

		### fit poisson, zero-inflated, and hurdele distributions to entire data set using site intercepts
		dist.test(use, parasite = "ECHINOLOAD", nitt = c(2000000,20000000,10), burnin = c(1000000,10000000,1), thin = c(1000,10000,1), snails = "no")

		### upload models
		for(i in c("poismod","zimod")){
				load(paste0("data/results/models/abundance/ECHINOLOAD_",i,".RData"), verbose = TRUE)
		}

		### compare autocorrelations
		acorr(pois.mod,summ.zp(pois.mod)$cstats[3])
		acorr(zi.mod,summ.zp(zi.mod)$cstats[3])

		### what is the expected number of zeros?
		{
			ppcheck.zeros(use,pois.mod,parasite = "ECHINOLOAD")
			ppcheck.zeros(use,zi.mod,parasite = "ECHINOLOAD", prob.model = "zero-inflated")
		}

		### simulate abundances
		sim_echino_means <- simulate_abundances(sub = use, mod = pois.mod);  save(sim_echino_means, file = "output/sim_echino_means.RData")

		### pull posteriors of site and species effects
		echino_effects <- simulate_abundance_effects(sub = use, mod = pois.mod); save(echino_effects, file = "output/echino_effects.RData")



	}

	### plot predicted species means
	{

		#### upload model results
		{
			load("output/sim_rib_means.RData")
			load("output/sim_bd_means.RData")
			load("output/sim_rv_means.RData")
			load("output/sim_echino_means.RData")
		}

		#### extract posteriors and reformat
		{
			rib.means <- sim_rib_means$species_means
			bd.means <- sim_bd_means$species_means
			rv.means <- sim_rv_means$species_means
			echino.means <- sim_echino_means$species_means

		}

		### create summary tables
		{

			rib.summ <- species_abundance_summ(sim_rib_means,"rib")
			bd.summ <- species_abundance_summ(sim_bd_means,"bd")
			rv.summ <- species_abundance_summ(sim_rv_means,"rv")
			echino.summ <- species_abundance_summ(sim_echino_means,"echino")
			sp.summ <- rbind(rib.summ,bd.summ,rv.summ,echino.summ)
		}

		# change levels
		levels(sp.summ$species) <- c("A. boreas","P. regilla","L. catesbeianus","T. granulosa","T. torosa")
		levels(sp.summ$parasite) <- c("R. ondatrae","B. dendrobatidis","ranavirus","Echinostoma")

		### plot
		Cairo(file = "output/species_mean_abundances", type = "svg", units = "in", width = 7, height = 4, family = "regular", pointsize = 12)
		spec.ma.plot <- ggplot(data = sp.summ, aes(species,log10(mean))) +
			theme_bw() +
			scale_fill_brewer(guide = guide_legend(title = "host species"),type = "qual",palette = "Set3") +
			labs(y = "log10(mean abundnace)") +
			geom_errorbar(aes(ymin = log10(low95),ymax = log10(hi95)), width = 0.1) +
			geom_point(size = 3, pch = 21,aes(fill = species, group = species)) +
			facet_grid(.~parasite) +
			theme(axis.text.x = element_blank(),
				axis.title.x = element_blank(),
				legend.text = element_text(face = "italic"),
				strip.text = element_text(face = "italic"),
				strip.background = element_blank(),
				panel.grid.major.x = element_blank())
			print(spec.ma.plot)
		dev.off()



	}

	### plot predicted species effects
	{

		#### upload models
		{
			load("output/rib_effects.RData")
			load("output/bd_effects.RData")
			load("output/rv_effects.RData")
			load("output/echino_effects.RData")
		}

		### create summary tables
		{
			rib.summ <- species_effects_summ(rib_effects,"rib")
			bd.summ <- species_effects_summ(bd_effects,"bd")
			rv.summ <- species_effects_summ(rv_effects,"rv")
			echino.summ <- species_effects_summ(echino_effects,"echino")
			sp.summ <- rbind(rib.summ,bd.summ,rv.summ,echino.summ)
		}

		# change levels
		levels(sp.summ$species) <- c("A. boreas","P. regilla","L. catesbeianus","T. granulosa","T. torosa")
		levels(sp.summ$parasite) <- c("R. ondatrae","B. dendrobatidis","ranavirus","Echinostoma")

		### plot
		Cairo(file = "docs/figs/species_mean_effects", type = "svg", units = "in", width = 7, height = 4, family = "regular", pointsize = 12)
		spec.ma.plot <- ggplot(data = sp.summ, aes(species,(mean))) +
			theme_bw() +
			scale_fill_brewer(guide = guide_legend(title = "host species"),type = "qual",palette = "Set3") +
			labs(y = expression(paste("species effects (", alpha[k],")", ))) +
			geom_errorbar(aes(ymin = (low95),ymax = (hi95)), width = 0.1) +
			geom_point(size = 3, pch = 21,aes(fill = species, group = species)) +
			facet_grid(.~parasite) +
			theme(axis.text.x = element_blank(),
				axis.title.x = element_blank(),
				legend.text = element_text(face = "italic"),
				strip.text = element_text(face = "italic"),
				strip.background = element_blank(),
				panel.grid.major.x = element_blank())
		print(spec.ma.plot)
		dev.off()



}

#############################################################
# multi-response models to look for correlations in abundance
#############################################################
{



	### Correlations between Rib and Bd Loads (RBmod)
	{

		## extract data from sites where we found
		{ # 38 sites, 963 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 1,			# sites where we checked for Rib in at least this many animals
				found.rib = 1,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 1,			# sites where we checked for Bd in at least this many animals
				found.bd = 1,			# sites where we found Bd in at least this many animals
				check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 0,			# sites where we found Ranavirus in at least this many animals
				check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
				found.echino = 0		# sites where we found Echinostomes in at least this many animals

			)
		}


		## fit models for various variance/covariance structures
		{

				### with individual covariance estimated
				RBmod1 <- multi.test(use,
					predictor1 = "RIBLOAD",		# predictor variable 1
					model1 = "poisson",			# model structure for variable 1
					predictor2 = "BDLOAD",		# predictor variable 2
					model2 = "poisson",			# model structure for variable 2
					site.covar = "yes",			# fit covariance between parasites among sites?
					species.covar = "yes",		# fit covariance between parasites among species?
					ind.covar = "yes",			# fit covariance between parasites among individuals?
					fix = "no",
					nitt = 2000000,
					burnin = 1000000,
					thin = 1000
				)
				save(RBmod1, file = "output/RBmod1.RData", compress = "gzip")

		}

		## upload  models
		for(i in c(1)){
			load(paste0("data/results/models/corr.abundance/RBmod",i,".RData"), verbose = TRUE)
		}

		### check auto correlation
		acorr(RBmod1,summ.zp(RBmod1)$cstats[3])

		### add residuals
		{
			RBmod1 <- add.resid(RBmod1, response = "multi"); save(RBmod1, file = "output/RBmod1.RData", compress = "gzip")
		}

		### plot MCMC
		{
			Cairo(file = "output/RBmod1", type = "pdf", units = "in", width = 5, height = 7, family = "regular", pointsize = 7, onefile = TRUE)
			plot(RBmod1)
			dev.off()
		}

		## generate model summaries
		compare <- mod.summ(RBmod1, predictor1 = "RIBLOAD", predictor2 = "BDLOAD"); compare





	}

	### Correlations between Rib and Ranavirus Loads (RRmod)
	{

		## extract data from sites where we found
		{ # 41 sites, 927 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 1,			# sites where we checked for Rib in at least this many animals
				found.rib = 1,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 0,			# sites where we checked for Bd in at least this many animals
				found.bd = 0,			# sites where we found Bd in at least this many animals
				check.rv = 1,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 1,			# sites where we found Ranavirus in at least this many animals
				check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
				found.echino = 0		# sites where we found Echinostomes in at least this many animals
			)
		}


		## fit model fits for various variance/covariance structures
		{


				### with individual covariance estimated
				RRmod1 <- multi.test(use,
					predictor1 = "RIBLOAD",		# predictor variable 1
					model1 = "poisson",			# model structure for variable 1
					predictor2 = "RVLOAD",		# predictor variable 2
					model2 = "poisson",			# model structure for variable 2
					site.covar = "yes",			# fit covariance between parasites among sites?
					species.covar = "yes",		# fit covariance between parasites among species?
					ind.covar = "yes",			# fit covariance between parasites among individuals?
					fix = "no",
					nitt = 8000000,
					burnin = 4000000,
					thin = 4000
			)
			save(RRmod1, file = "output/RRmod1.RData", compress = "gzip")


		}

		## upload poisson/poisson models
		for(i in c(1)){
			load(paste0("data/results/models/corr.abundance/RRmod",i,".RData"), verbose = TRUE)
		}

		### check auto correlation
		acorr(RRmod1,summ.zp(RRmod1)$cstats[3])

		### add residuals
		{
			RRmod1 <- add.resid(RRmod1, response = "multi"); save(RRmod1, file = "output/RRmod1.RData", compress = "gzip")
		}

		### plot MCMC
		{
			Cairo(file = "output/RRmod1", type = "pdf", units = "in", width = 5, height = 7, family = "regular", pointsize = 7, onefile = TRUE)
			plot(RRmod1)
			dev.off()
		}

		## generate model summaries
		compare <- mod.summ(RRmod1, predictor1 = "RIBLOAD", predictor2 = "RVLOAD");compare

	}

	### Correlations between Bd and Ranavirus Loads (BRmod)
	{

		## extract data from sites where we found
		{ # 37 sites, 872 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 0,			# sites where we checked for Rib in at least this many animals
				found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 1,			# sites where we checked for Bd in at least this many animals
				found.bd = 1,			# sites where we found Bd in at least this many animals
				check.rv = 1,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 1,			# sites where we found Ranavirus in at least this many animals
				check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
				found.echino = 0		# sites where we found Echinostomes in at least this many animals
			)
		}


		## fit model fits for various variance/covariance structures
		{
				### with covariance
				BRmod1 <- multi.test(use,
					predictor1 = "BDLOAD",		# predictor variable 1
					model1 = "poisson",			# model structure for variable 1
					predictor2 = "RVLOAD",		# predictor variable 2
					model2 = "poisson",			# model structure for variable 2
					site.covar = "yes",			# fit covariance between parasites among sites?
					species.covar = "yes",		# fit covariance between parasites among species?
					ind.covar = "yes",				# fit covariance between parasites among individuals?
					fix = "no",
					nitt = 6000000,
					burnin = 3000000,
					thin = 3000
				)
				save(BRmod1, file = "output/BRmod1.RData", compress = "gzip")
		}

		## upload poisson/poisson models
		for(i in c(1)){
			load(paste0("data/results/models/corr.abundance/BRmod",i,".RData"), verbose = TRUE)
		}

		### check auto correlation
		acorr(BRmod1,summ.zp(BRmod1)$cstats[3])

		### add residuals
		{
			BRmod1 <- add.resid(BRmod1, response = "multi"); save(BRmod1, file = "output/BRmod1.RData", compress = "gzip")
		}

		### plot MCMC
		{
			Cairo(file = "output/BRmod1", type = "pdf", units = "in", width = 5, height = 7, family = "regular", pointsize = 7, onefile = TRUE)
			plot(BRmod1)
			dev.off()
		}

		## generate model summaries
		compare <- mod.summ(BRmod1, predictor1 = "BDLOAD", predictor2 = "RVLOAD");compare

	}

	### Correlations between Rib and Echino Loads (REmod)
	{

		## extract data from sites where we found
		{ # 55 sites, 1270 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 1,			# sites where we checked for Rib in at least this many animals
				found.rib = 1,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 0,			# sites where we checked for Bd in at least this many animals
				found.bd = 0,			# sites where we found Bd in at least this many animals
				check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 0,			# sites where we found Ranavirus in at least this many animals
				check.trem = 0,		# sites where we checked for Rib in at least this many animals
				found.trem = 0,			# sites where we found Ribeiroia in at least this many animals
				check.echino = 1,		# sites where we checked for Rib in at least this many animals
				found.echino = 1			# sites where we found Ribeiroia in at least this many animals

			)
		}


		## fit model fits for various variance/covariance structures
		{

				### with individual covariance estimated
				REmod1 <- multi.test(use,
					predictor1 = "RIBLOAD",	# predictor variable 1
					model1 = "poisson",			# model structure for variable 1
					predictor2 = "ECHINOLOAD",		# predictor variable 2
					model2 = "poisson",			# model structure for variable 2
					site.covar = "yes",			# fit covariance between parasites among sites?
					species.covar = "yes",		# fit covariance between parasites among species?
					ind.covar = "yes",			# fit covariance between parasites among individuals?
					fix = "no",
					nitt = 6000000,
					burnin = 3000000,
					thin = 3000
				)
				save(REmod1, file = "output/REmod1.RData", compress = "gzip")


		}

		## upload poisson/poisson models
		for(i in c(1)){
			load(paste0("data/results/models/corr.abundance/REmod",i,".RData"), verbose = TRUE)
		}

		### check auto correlation
		acorr(REmod1,summ.zp(REmod1)$cstats[3])

		### add residuals
		{
			REmod1 <- add.resid(REmod1, response = "multi"); save(REmod1, file = "output/REmod1.RData", compress = "gzip")
		}

		### plot MCMC
		{
			Cairor(file = "output/REmod1", type = "pdf", units = "in", width = 5, height = 7, family = "regular", pointsize = 7, onefile = TRUE)
			plot(REmod1)
			dev.off()
		}

		## generate model summaries
		compare <- mod.summ(REmod1, predictor1 = "RIBLOAD", predictor2 = "ECHINOLOAD"); compare

	}

	### Correlations between Echino and Bd Loads (EBmod)
	{

		## extract data from sites where we found
		{ # 50 sites, 1145 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 0,			# sites where we checked for Rib in at least this many animals
				found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 1,			# sites where we checked for Bd in at least this many animals
				found.bd = 1,			# sites where we found Bd in at least this many animals
				check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 0,			# sites where we found Ranavirus in at least this many animals
				check.trem = 0,		# sites where we checked for Rib in at least this many animals
				found.trem = 0,			# sites where we found Ribeiroia in at least this many animals
				check.echino = 1,		# sites where we checked for Rib in at least this many animals
				found.echino = 1			# sites where we found Ribeiroia in at least this many animals

			)
		}


		## fit model fits for various variance/covariance structures
		{

				### with individual covariance estimated
				EBmod1 <- multi.test(use,
					predictor1 = "ECHINOLOAD",	# predictor variable 1
					model1 = "poisson",			# model structure for variable 1
					predictor2 = "BDLOAD",		# predictor variable 2
					model2 = "poisson",			# model structure for variable 2
					site.covar = "yes",			# fit covariance between parasites among sites?
					species.covar = "yes",		# fit covariance between parasites among species?
					ind.covar = "yes",			# fit covariance between parasites among individuals?
					fix = "no",
					nitt = 4000000,
					burnin = 2000000,
					thin = 2000
				)
				save(EBmod1, file = "output/EBmod1.RData", compress = "gzip")


		}

		## upload poisson/poisson models
		for(i in c(1)){
			load(paste0("data/results/models/corr.abundance/EBmod",i,".RData"), verbose = TRUE)
		}

		### check auto correlation
		acorr(EBmod1,summ.zp(EBmod1)$cstats[3])

		### add residuals
		{
			EBmod1 <- add.resid(EBmod1, response = "multi"); save(EBmod1, file = "output/EBmod1.RData", compress = "gzip")
		}

		### plot MCMC
		{
			Cairo(file = "output/EBmod1", type = "pdf", units = "in", width = 5, height = 7, family = "regular", pointsize = 7, onefile = TRUE)
			plot(EBmod1)
			dev.off()
		}

		## generate model summaries
		compare <- mod.summ(EBmod1, predictor1 = "ECHINOLOAD", predictor2 = "BDLOAD"); compare


	}

	### Correlations between Echino and Ranavirus Loads (ERmod)
	{

		## extract data from sites where we found
		{ # 55 sites, 1168 individuals
			use <- coinf.data(dat,site,use.NAs = "no",
				check.rib = 0,			# sites where we checked for Rib in at least this many animals
				found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
				check.bd = 0,			# sites where we checked for Bd in at least this many animals
				found.bd = 0,			# sites where we found Bd in at least this many animals
				check.rv = 1,			# sites where we checked for Ranavirus in at least this many animals
				found.rv = 1,			# sites where we found Ranavirus in at least this many animals
				check.echino = 1,		# sites where we checked for Rib in at least this many animals
				found.echino = 1			# sites where we found Ribeiroia in at least this many animals
			)
		}


		## fit model fits for various variance/covariance structures
		{

				### with individual covariance estimated
				ERmod1 <- multi.test(use,
					predictor1 = "ECHINOLOAD",	# predictor variable 1
					model1 = "poisson",			# model structure for variable 1
					predictor2 = "RVLOAD",		# predictor variable 2
					model2 = "poisson",			# model structure for variable 2
					site.covar = "yes",			# fit covariance between parasites among sites?
					species.covar = "yes",		# fit covariance between parasites among species?
					ind.covar = "yes",			# fit covariance between parasites among individuals?
					fix = "no",
					nitt = 4000000,
					burnin = 2000000,
					thin = 2000
				)
				save(ERmod1, file = "output/ERmod1.RData", compress = "gzip")

		}

		## upload poisson/poisson models
		for(i in c(1)){
			load(paste0("data/results/models/corr.abundance/ERmod",i,".RData"), verbose = TRUE)
		}

		### check auto correlation
		acorr(ERmod1,summ.zp(ERmod1)$cstats[3])

		### add residuals
		{
			ERmod1 <- add.resid(ERmod1, response = "multi"); save(ERmod1, file = "output/ERmod1.RData", compress = "gzip")
		}

		### plot MCMC
		{
			Cairo(file = "output/ERmod1", type = "pdf", units = "in", width = 5, height = 7, family = "regular", pointsize = 7, onefile = TRUE)
			plot(ERmod1)
			dev.off()
		}

		## generate model summaries
		compare <- mod.summ(ERmod1, predictor1 = "ECHINOLOAD", predictor2 = "RVLOAD"); compare


	}

	### results
	{

			### upload models
			for(i in c(1)){
				load(paste0("data/results/models/corr.abundance/RBmod",i,".RData"), verbose = TRUE)
				load(paste0("data/results/models/corr.abundance/RRmod",i,".RData"), verbose = TRUE)
				load(paste0("data/results/models/corr.abundance/BRmod",i,".RData"), verbose = TRUE)
				load(paste0("data/results/models/corr.abundance/EBmod",i,".RData"), verbose = TRUE)
				load(paste0("data/results/models/corr.abundance/ERmod",i,".RData"), verbose = TRUE)
				load(paste0("data/results/models/corr.abundance/REmod",i,".RData"), verbose = TRUE)
			}

			## generate model summaries
			{
				compare <- rbind(mod.summ(RBmod1, predictor1 = "RIBLOAD", predictor2 = "BDLOAD"),
				mod.summ(RRmod1, predictor1 = "RIBLOAD", predictor2 = "RVLOAD"),
				mod.summ(BRmod1, predictor1 = "BDLOAD", predictor2 = "RVLOAD"),
				mod.summ(REmod1, predictor1 = "RIBLOAD", predictor2 = "ECHINOLOAD"),
				mod.summ(EBmod1, predictor1 = "BDLOAD", predictor2 = "ECHINOLOAD"),
				mod.summ(ERmod1, predictor1 = "RVLOAD", predictor2 = "ECHINOLOAD"))

			}

			### species level correlations
			plot.species.corr(RBmod1,RRmod1,BRmod1,REmod1,EBmod1,ERmod1,compare)

			### site level correlations
			plot.site.corr(RBmod1,RRmod1,BRmod1,REmod1,EBmod1,ERmod1,compare)

			### individual level correlations
			plot.ind.corr(RBmod1,RRmod1,BRmod1,REmod1,EBmod1, ERmod1, compare)

			### dot and line plot of correlation posteriors
			{


				## generate data
				{
					## generate model summaries (takes a long time....)
					compare <- rbind(mod.summ(RBmod1,predictor1 = "RIBLOAD", predictor2 = "BDLOAD"),
											mod.summ(RRmod1,predictor1 = "RIBLOAD", predictor2 = "RVLOAD"),
											mod.summ(BRmod1,predictor1 = "BDLOAD", predictor2 = "RVLOAD")
					)

					## add combination column
					compare$combo <- paste(compare$predictor1, compare$predictor2, sep = "_")
					compare$combo <- factor(compare$combo)
					levels(compare$combo) <- c("B. dendrobatidis vs.\\nRanavirus", "R. ondatrae vs.\\nB. dendrobatidis","R. ondatrae vs.\\nRanavirus")

					## melt data farme
					compare.melted <- rbind() # new data frame
					for(i in c(1:nrow(compare))){
						new_i <- data.frame(combo = rep(compare[i,"combo"],3),
							level = c("site","species","individual"),
							corr = c(compare[i,"site.cor"],compare[i,"species.cor"],compare[i,"ind.cor"]),
							low = c(compare[i,"site.cred.low"],compare[i,"species.cred.low"],compare[i,"ind.cred.low"]),
							hi = c(compare[i,"site.cred.hi"],compare[i,"species.cred.hi"],compare[i,"ind.cred.hi"]),
							strong = c(compare[i,"site.strong"],compare[i,"species.strong"],compare[i,"ind.strong"]))
						compare.melted <- rbind(compare.melted,new_i)
					}

					compare.melted$combo <- factor(compare.melted$combo, levels = c("R. ondatrae vs.\\nB. dendrobatidis","R. ondatrae vs.\\nRanavirus","B. dendrobatidis vs.\\nRanavirus"))
				}

				## plot
				Cairo(file = "output/all_level_correlations", type = "svg", units = "in", width = 10, height = 4, family = "regular", pointsize = 12)
				pl <- ggplot(data = compare.melted, aes(level, corr)) +
						theme_bw() +
						coord_flip() +
						#scale_fill_manual(values = c("white","black")) +
						geom_hline(y = 0, lty = 3) +
						labs(y = expression(paste("correlation (",rho,")"))) +
						geom_errorbar(aes(ymin = low,ymax = hi), width = 0.1) +
						geom_point(size = 4, pch = 21, fill = "black") +
						facet_grid(.~combo) +
						theme(legend.text = element_text(face = "italic"),
								axis.title.y = element_blank(),
							strip.text = element_text(face = "italic", size = 17),
							axis.title.x = element_text(size = 17),
							axis.text = element_text(size = 16),

							strip.background = element_blank())
				print(pl)
				dev.off()

			}

			### plot of all posteriors
			{
				## make plots
				{
					## Rib vs. Bd
					{

						# remove excess covariances
						posts <- RBmod1$VCV[,-c(3,7,15)]

						# changle names
						colnames(posts) <- c("A","B","C","D","E","F","G","H","I","J","K","L","M")

						# summary
						{
							plot.data <- data.frame(parameter = colnames(posts),
													mean = summary(posts)[1]$statistics[,1],
													low95 =  summary(posts)[2]$quantiles[,1],
													hi95 =  summary(posts)[2]$quantiles[,5])
						}

						# re factor
						plot.data$parameter <- factor(plot.data$parameter, levels = colnames(posts))


						# plot
						{
							rbplot <- ggplot(data = plot.data, aes(parameter,mean, ymin = low95, ymax = hi95), parse = TRUE) +
								coord_flip() +
								geom_hline(y = 0, color = "grey50", lty = 2, size = 0.5) +
								geom_pointrange(size = 0.5) +

								scale_x_discrete(labels = c('A' = expression({sigma^2}[eta[paste(i," Rib")]]),
																		'B' = expression(paste({rho},{sigma[eta[paste(i,", Rib")]]},{sigma[eta[paste(i,", Bd")]]})),
																		'C' = expression({sigma^2}[eta[paste(i,", Bd")]]),
																		'D' = expression({sigma^2}[theta[paste(k,", Rib")]]),
																		'E' = expression(paste({rho},{sigma[theta[paste(k,", Rib")]]},{sigma[theta[paste(k,", Bd")]]})),
																		'F' = expression({sigma^2}[theta[paste(k,", Bd")]]),
																		'G' = expression({sigma^2}[beta[paste(k,", Rib")]]),
																		'H' = expression({sigma^2}[beta[paste(k,", Bd")]]),
																		'I' = expression({sigma^2}[gamma[paste(l,", Rib")]]),
																		'J' = expression({sigma^2}[gamma[paste(l,", Bd")]]),
																		'K' = expression({sigma^2}[epsilon[paste(i," Rib")]]),
																		'L' = expression(paste({rho},{sigma[epsilon[paste(i,", Rib")]]},{sigma[epsilon[paste(i,", Bd")]]})),
																		'M' = expression({sigma^2}[epsilon[paste(i,", Bd")]]))) +
								theme(text = element_text(size = 15),
										plot.title = element_text(size = 10, face = "italic"),
										plot.margin = unit(c(0,0,0,0),"mm"),
										axis.text = element_text(color = "black"),
										axis.ticks = element_line(color = "black"),
										axis.title = element_blank(),
										panel.background = element_rect(fill = NA, size = 0.7, color = "black"),
										panel.grid = element_blank())
						}
					}

					## Rib vs. Rv
					{

						# remove excess covariances
						posts <- RRmod1$VCV[,-c(3,7,15)]

						# changle names
						colnames(posts) <- c("A","B","C","D","E","F","G","H","I","J","K","L","M")

						# summary
						{
							plot.data <- data.frame(parameter = colnames(posts),
													mean = summary(posts)[1]$statistics[,1],
													low95 =  summary(posts)[2]$quantiles[,1],
													hi95 =  summary(posts)[2]$quantiles[,5])
						}

						# re factor
						plot.data$parameter <- factor(plot.data$parameter, levels = colnames(posts))


						# plot
						{
							rrplot <- ggplot(data = plot.data, aes(parameter,mean, ymin = low95, ymax = hi95), parse = TRUE) +
								coord_flip() +
								geom_hline(y = 0, color = "grey50", lty = 2, size = 0.5) +
								geom_pointrange(size = 0.5) +

								scale_x_discrete(labels = c('A' = expression({sigma^2}[eta[paste(i," Rib")]]),
																		'B' = expression(paste({rho},{sigma[eta[paste(i,", Rib")]]},{sigma[eta[paste(i,", Rv")]]})),
																		'C' = expression({sigma^2}[eta[paste(i,", Rv")]]),
																		'D' = expression({sigma^2}[theta[paste(k,", Rib")]]),
																		'E' = expression(paste({rho},{sigma[theta[paste(k,", Rib")]]},{sigma[theta[paste(k,", Rv")]]})),
																		'F' = expression({sigma^2}[theta[paste(k,", Rv")]]),
																		'G' = expression({sigma^2}[beta[paste(k,", Rib")]]),
																		'H' = expression({sigma^2}[beta[paste(k,", Rv")]]),
																		'I' = expression({sigma^2}[gamma[paste(l,", Rib")]]),
																		'J' = expression({sigma^2}[gamma[paste(l,", Rv")]]),
																		'K' = expression({sigma^2}[epsilon[paste(i," Rib")]]),
																		'L' = expression(paste({rho},{sigma[epsilon[paste(i,", Rib")]]},{sigma[epsilon[paste(i,", Rv")]]})),
																		'M' = expression({sigma^2}[epsilon[paste(i,", Rv")]]))) +
								theme(text = element_text(size = 15),
										plot.title = element_text(size = 10, face = "italic"),
										plot.margin = unit(c(0,0,0,0),"mm"),
										axis.text = element_text(color = "black"),
										axis.ticks = element_line(color = "black"),
										axis.title = element_blank(),
										panel.background = element_rect(fill = NA, size = 0.7, color = "black"),
										panel.grid = element_blank())
						}
					}

					## Bd vs. Rv
					{

						# remove excess covariances
						posts <- BRmod1$VCV[,-c(3,7,15)]

						# changle names
						colnames(posts) <- c("A","B","C","D","E","F","G","H","I","J","K","L","M")

						# summary
						{
							plot.data <- data.frame(parameter = colnames(posts),
													mean = summary(posts)[1]$statistics[,1],
													low95 =  summary(posts)[2]$quantiles[,1],
													hi95 =  summary(posts)[2]$quantiles[,5])
						}

						# re factor
						plot.data$parameter <- factor(plot.data$parameter, levels = colnames(posts))


						# plot
						{
							brplot <- ggplot(data = plot.data, aes(parameter,mean, ymin = low95, ymax = hi95), parse = TRUE) +
								coord_flip() +
								geom_hline(y = 0, color = "grey50", lty = 2, size = 0.5) +
								geom_pointrange(size = 0.5) +

								scale_x_discrete(labels = c('A' = expression({sigma^2}[eta[paste(i," Bd")]]),
																		'B' = expression(paste({rho},{sigma[eta[paste(i,", Bd")]]},{sigma[eta[paste(i,", Rv")]]})),
																		'C' = expression({sigma^2}[eta[paste(i,", Rv")]]),
																		'D' = expression({sigma^2}[theta[paste(k,", Bd")]]),
																		'E' = expression(paste({rho},{sigma[theta[paste(k,", Bd")]]},{sigma[theta[paste(k,", Rv")]]})),
																		'F' = expression({sigma^2}[theta[paste(k,", Rv")]]),
																		'G' = expression({sigma^2}[beta[paste(k,", Bd")]]),
																		'H' = expression({sigma^2}[beta[paste(k,", Rv")]]),
																		'I' = expression({sigma^2}[gamma[paste(l,", Bd")]]),
																		'J' = expression({sigma^2}[gamma[paste(l,", Rv")]]),
																		'K' = expression({sigma^2}[epsilon[paste(i," Bd")]]),
																		'L' = expression(paste({rho},{sigma[epsilon[paste(i,", Bd")]]},{sigma[epsilon[paste(i,", Rv")]]})),
																		'M' = expression({sigma^2}[epsilon[paste(i,", Rv")]]))) +
								theme(text = element_text(size = 15),
										plot.title = element_text(size = 10, face = "italic"),
										plot.margin = unit(c(0,0,0,0),"mm"),
										axis.text = element_text(color = "black"),
										axis.ticks = element_line(color = "black"),
										axis.title = element_blank(),
										panel.background = element_rect(fill = NA, size = 0.7, color = "black"),
										panel.grid = element_blank())
						}
					}

					## Rib vs. Echino
					{

						# remove excess covariances
						posts <- REmod1$VCV[,-c(3,7,15)]

						# changle names
						colnames(posts) <- c("A","B","C","D","E","F","G","H","I","J","K","L","M")

						# summary
						{
							plot.data <- data.frame(parameter = colnames(posts),
													mean = summary(posts)[1]$statistics[,1],
													low95 =  summary(posts)[2]$quantiles[,1],
													hi95 =  summary(posts)[2]$quantiles[,5])
						}

						# re factor
						plot.data$parameter <- factor(plot.data$parameter, levels = colnames(posts))


						# plot
						{
							replot <- ggplot(data = plot.data, aes(parameter,mean, ymin = low95, ymax = hi95), parse = TRUE) +
								coord_flip() +
								geom_hline(y = 0, color = "grey50", lty = 2, size = 0.5) +
								geom_pointrange(size = 0.5) +

								scale_x_discrete(labels = c('A' = expression({sigma^2}[eta[paste(i," Rib")]]),
																		'B' = expression(paste({rho},{sigma[eta[paste(i,", Rib")]]},{sigma[eta[paste(i,", Ech")]]})),
																		'C' = expression({sigma^2}[eta[paste(i,", Ech")]]),
																		'D' = expression({sigma^2}[theta[paste(k,", Rib")]]),
																		'E' = expression(paste({rho},{sigma[theta[paste(k,", Rib")]]},{sigma[theta[paste(k,", Ech")]]})),
																		'F' = expression({sigma^2}[theta[paste(k,", Ech")]]),
																		'G' = expression({sigma^2}[beta[paste(k,", Rib")]]),
																		'H' = expression({sigma^2}[beta[paste(k,", Ech")]]),
																		'I' = expression({sigma^2}[gamma[paste(l,", Rib")]]),
																		'J' = expression({sigma^2}[gamma[paste(l,", Ech")]]),
																		'K' = expression({sigma^2}[epsilon[paste(i," Rib")]]),
																		'L' = expression(paste({rho},{sigma[epsilon[paste(i,", Rib")]]},{sigma[epsilon[paste(i,", Ech")]]})),
																		'M' = expression({sigma^2}[epsilon[paste(i,", Ech")]]))) +
								theme(text = element_text(size = 15),
										plot.title = element_text(size = 10, face = "italic"),
										plot.margin = unit(c(0,0,0,0),"mm"),
										axis.text = element_text(color = "black"),
										axis.ticks = element_line(color = "black"),
										axis.title = element_blank(),
										panel.background = element_rect(fill = NA, size = 0.7, color = "black"),
										panel.grid = element_blank())
						}
					}

					## Echimo vs. Bd
					{

						# remove excess covariances
						posts <- EBmod1$VCV[,-c(3,7,15)]

						# changle names
						colnames(posts) <- c("A","B","C","D","E","F","G","H","I","J","K","L","M")

						# summary
						{
							plot.data <- data.frame(parameter = colnames(posts),
													mean = summary(posts)[1]$statistics[,1],
													low95 =  summary(posts)[2]$quantiles[,1],
													hi95 =  summary(posts)[2]$quantiles[,5])
						}

						# re factor
						plot.data$parameter <- factor(plot.data$parameter, levels = colnames(posts))


						# plot
						{
							ebplot <- ggplot(data = plot.data, aes(parameter,mean, ymin = low95, ymax = hi95), parse = TRUE) +
								coord_flip() +
								geom_hline(y = 0, color = "grey50", lty = 2, size = 0.5) +
								geom_pointrange(size = 0.5) +

								scale_x_discrete(labels = c('A' = expression({sigma^2}[eta[paste(i," Ech")]]),
																		'B' = expression(paste({rho},{sigma[eta[paste(i,", Ech")]]},{sigma[eta[paste(i,", Bd")]]})),
																		'C' = expression({sigma^2}[eta[paste(i,", Bd")]]),
																		'D' = expression({sigma^2}[theta[paste(k,", Ech")]]),
																		'E' = expression(paste({rho},{sigma[theta[paste(k,", Ech")]]},{sigma[theta[paste(k,", Bd")]]})),
																		'F' = expression({sigma^2}[theta[paste(k,", Bd")]]),
																		'G' = expression({sigma^2}[beta[paste(k,", Ech")]]),
																		'H' = expression({sigma^2}[beta[paste(k,", Bd")]]),
																		'I' = expression({sigma^2}[gamma[paste(l,", Ech")]]),
																		'J' = expression({sigma^2}[gamma[paste(l,", Bd")]]),
																		'K' = expression({sigma^2}[epsilon[paste(i," Ech")]]),
																		'L' = expression(paste({rho},{sigma[epsilon[paste(i,", Ech")]]},{sigma[epsilon[paste(i,", Bd")]]})),
																		'M' = expression({sigma^2}[epsilon[paste(i,", Bd")]]))) +
								theme(text = element_text(size = 15),
										plot.title = element_text(size = 10, face = "italic"),
										plot.margin = unit(c(0,0,0,0),"mm"),
										axis.text = element_text(color = "black"),
										axis.ticks = element_line(color = "black"),
										axis.title = element_blank(),
										panel.background = element_rect(fill = NA, size = 0.7, color = "black"),
										panel.grid = element_blank())
						}
					}

					## Echimo vs. Bd
					{

						# remove excess covariances
						posts <- ERmod1$VCV[,-c(3,7,15)]

						# changle names
						colnames(posts) <- c("A","B","C","D","E","F","G","H","I","J","K","L","M")

						# summary
						{
							plot.data <- data.frame(parameter = colnames(posts),
													mean = summary(posts)[1]$statistics[,1],
													low95 =  summary(posts)[2]$quantiles[,1],
													hi95 =  summary(posts)[2]$quantiles[,5])
						}

						# re factor
						plot.data$parameter <- factor(plot.data$parameter, levels = colnames(posts))


						# plot
						{
							erplot <- ggplot(data = plot.data, aes(parameter,mean, ymin = low95, ymax = hi95), parse = TRUE) +
								coord_flip() +
								geom_hline(y = 0, color = "grey50", lty = 2, size = 0.5) +
								geom_pointrange(size = 0.5) +

								scale_x_discrete(labels = c('A' = expression({sigma^2}[eta[paste(i," Ech")]]),
																		'B' = expression(paste({rho},{sigma[eta[paste(i,", Ech")]]},{sigma[eta[paste(i,", Rv")]]})),
																		'C' = expression({sigma^2}[eta[paste(i,", Rv")]]),
																		'D' = expression({sigma^2}[theta[paste(k,", Ech")]]),
																		'E' = expression(paste({rho},{sigma[theta[paste(k,", Ech")]]},{sigma[theta[paste(k,", Rv")]]})),
																		'F' = expression({sigma^2}[theta[paste(k,", Rv")]]),
																		'G' = expression({sigma^2}[beta[paste(k,", Ech")]]),
																		'H' = expression({sigma^2}[beta[paste(k,", Rv")]]),
																		'I' = expression({sigma^2}[gamma[paste(l,", Ech")]]),
																		'J' = expression({sigma^2}[gamma[paste(l,", Rv")]]),
																		'K' = expression({sigma^2}[epsilon[paste(i," Ech")]]),
																		'L' = expression(paste({rho},{sigma[epsilon[paste(i,", Ech")]]},{sigma[epsilon[paste(i,", Rv")]]})),
																		'M' = expression({sigma^2}[epsilon[paste(i,", Rv")]]))) +
								theme(text = element_text(size = 15),
										plot.title = element_text(size = 10, face = "italic"),
										plot.margin = unit(c(0,0,0,0),"mm"),
										axis.text = element_text(color = "black"),
										axis.ticks = element_line(color = "black"),
										axis.title = element_blank(),
										panel.background = element_rect(fill = NA, size = 0.7, color = "black"),
										panel.grid = element_blank())
						}
					}
				}

				## blank panel
				blank <- grid.rect(gp=gpar(col="white"),draw = FALSE,name = "stuff")

				## make figure
				{


					Cairo(file = "docs/figs/posteriors.svg", type = "svg", canvas = "white", units = "in",
						 width = 6.5, height = 9)


					pl <- grid.arrange(rbplot,rrplot,brplot,replot,ebplot,erplot,
						nrow = 3, ncol = 2)

					dev.off()

				}


			}

			## tables of posteriors
			{

					output.posteriors(RBmod1,"RBmod1")
					output.posteriors(RRmod1,"RRmod1")
					output.posteriors(BRmod1,"BRmod1")
					output.posteriors(REmod1,"REmod1")
					output.posteriors(EBmod1,"EBmod1")
					output.posteriors(ERmod1,"ERmod1")

			}
	}


}

#######################################################################################
# multi-response models with separate correlations with separate VCV's for each species
#######################################################################################
{



		### Correlations between Rib and Bd Loads (RBspecmod)
		{

			## extract data from sites where we found
			{ # 38 sites, 963 individuals
				RB.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 1,			# sites where we checked for Rib in at least this many animals
					found.rib = 1,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 1,			# sites where we checked for Bd in at least this many animals
					found.bd = 1,			# sites where we found Bd in at least this many animals
					check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 0,			# sites where we found Ranavirus in at least this many animals
					check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
					found.echino = 0		# sites where we found Echinostomes in at least this many animals

				)
			}


			## fit model for various variance/covariance structures
			{


					### with individual covariance estimated
					RBspecmod <- multi.test.by.species(RB.use,
						predictor1 = "RIBLOAD",		# predictor variable 1
						model1 = "poisson",			# model structure for variable 1
						predictor2 = "BDLOAD",		# predictor variable 2
						model2 = "poisson",			# model structure for variable 2
						site.covar = "yes",			# fit covariance between parasites among sites?
						species.covar = "yes",		# fit covariance between parasites among species?
						ind.covar = "yes",			# fit covariance between parasites among individuals?
						fix = "no",
						nitt = 21000000,
						burnin = 10000000,
						thin = 11000
					)
					save(RBspecmod, file = "output/RBspecmod.RData", compress = "gzip")



			}

			## upload  models
			for(i in c(1)){
				load(paste0("data/results/models/species.corr.abundance/RBspecmod.RData"), verbose = TRUE)
			}

			### check auto correlation
			acorr(RBspecmod,summ.zp(RBspecmod)$cstats[3])

			### add residuals
			{
				RBspecmod <- add.resid(RBspecmod, response = "multi"); save(RBspecmod, file = "output/RBspecmod.RData", compress = "gzip")
			}

			### plot MCMC
			{
				Cairo(file = "output/RBspecmod", type = "pdf", units = "in", width = 5, height = 7, family = "regular", pointsize = 7, onefile = TRUE)
				plot(RBspecmod)
				dev.off()
			}

			## generate model summaries
			compare <- mod.summ(RBspecmod, predictor1 = "RIBLOAD", predictor2 = "BDLOAD"); compare


	}

		### Correlations between Rib and Ranavirus Loads (RRspecmod)
		{

			## extract data from sites where we found
			{ # 41 sites, 927 individuals
				RR.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 1,			# sites where we checked for Rib in at least this many animals
					found.rib = 1,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 0,			# sites where we checked for Bd in at least this many animals
					found.bd = 0,			# sites where we found Bd in at least this many animals
					check.rv = 1,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 1,			# sites where we found Ranavirus in at least this many animals
					check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
					found.echino = 0		# sites where we found Echinostomes in at least this many animals

				)
			}


			## fit model fits for various variance/covariance structures
			{


				### with individual covariance estimated
				RRspecmod <- multi.test.by.species(RR.use,
					predictor1 = "RIBLOAD",		# predictor variable 1
					model1 = "poisson",			# model structure for variable 1
					predictor2 = "RVLOAD",		# predictor variable 2
					model2 = "poisson",			# model structure for variable 2
					site.covar = "yes",			# fit covariance between parasites among sites?
					species.covar = "yes",		# fit covariance between parasites among species?
					ind.covar = "yes",			# fit covariance between parasites among individuals?
					fix = "no",
					nitt = 10000000,
					burnin = 5000000,
					thin = 5000

				)
				save(RRspecmod, file = "output/RRspecmod.RData", compress = "gzip")


			}

			## upload poisson/poisson models
			for(i in c(1)){
				load(paste0("data/results/models/species.corr.abundance/RRspecmod.RData"), verbose = TRUE)
			}

			### check auto correlation
			acorr(RRspecmod,summ.zp(RRspecmod)$cstats[3])

			### add residuals
			{
				RRspecmod <- add.resid(RRspecmod, response = "multi"); save(RRmod1, file = "output/RRmod1.RData", compress = "gzip")
			}

			##generate model summaries
			compare <- mod.summ.species(RRspecmod, predictor1 = "RIBLOAD", predictor2 = "RVLOAD"); compare

			## plot model results
			plot.model(mod8,predictor1 = "RIBLOAD", predictor2 = "RVLOAD")


	}

		### Correlations between Bd and Ranavirus Loads (BRspecmod)
		{

			## extract data from sites where we found
			{ # 37 sites, 1167 individuals
				BR.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 0,			# sites where we checked for Rib in at least this many animals
					found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 1,			# sites where we checked for Bd in at least this many animals
					found.bd = 1,			# sites where we found Bd in at least this many animals
					check.rv = 1,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 1,			# sites where we found Ranavirus in at least this many animals
					check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
					found.echino = 0		# sites where we found Echinostomes in at least this many animals

				)
			}


			## fit model fits for various variance/covariance structures
			{

				### with covariance
				BRspecmod <- multi.test.by.species(use,
					predictor1 = "BDLOAD",		# predictor variable 1
					model1 = "poisson",			# model structure for variable 1
					predictor2 = "RVLOAD",		# predictor variable 2
					model2 = "poisson",			# model structure for variable 2
					site.covar = "yes",			# fit covariance between parasites among sites?
					species.covar = "yes",		# fit covariance between parasites among species?
					ind.covar = "yes",				# fit covariance between parasites among individuals?
					fix = "no",
					itt = 10000000,
					burnin = 5000000,
					thin = 5000
				)
				save(BRspecmod, file = "output/BRspecmod.RData", compress = "gzip")

			}

			## upload poisson/poisson models
			for(i in c(1)){
				load(paste0("data/results/models/species.corr.abundance/BRspecmod.RData"), verbose = TRUE)
			}

			### check auto correlation
			acorr(BRspecmod,summ.zp(BRspecmod)$cstats[3])

			### add residuals
			{
				BRmod1 <- add.resid(BRmod1, response = "multi"); save(BRmod1, file = "output/BRmod1.RData", compress = "gzip")
			}

			##generate model summaries
			compare <- mod.summ.species(BRspecmod, predictor1 = "BDLOAD", predictor2 = "RVLOAD"); compare

			## plot model results
			plot.model(mod8,predictor1 = "RIBLOAD", predictor2 = "RVLOAD")

			}

		### Correlations between Rib and Echino Loads (REspecmod)
		{

			## extract data from sites where we found
			{ # 55 sites, 1270 individuals
				RE.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 1,			# sites where we checked for Rib in at least this many animals
					found.rib = 1,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 0,			# sites where we checked for Bd in at least this many animals
					found.bd = 0,			# sites where we found Bd in at least this many animals
					check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 0,		# sites where we found Ranavirus in at least this many animals
					check.echino = 1,
					found.echino = 1
				)
			}


			## fit model for various variance/covariance structures
			{


					### with individual covariance estimated
					REspecmod <- multi.test.by.species(RE.use,
						predictor1 = "RIBLOAD",		# predictor variable 1
						model1 = "poisson",			# model structure for variable 1
						predictor2 = "ECHINOLOAD",		# predictor variable 2
						model2 = "poisson",			# model structure for variable 2
						site.covar = "yes",			# fit covariance between parasites among sites?
						species.covar = "yes",		# fit covariance between parasites among species?
						ind.covar = "yes",			# fit covariance between parasites among individuals?
						fix = "no",
						nitt = 20000000,
						burnin = 10000000,
						thin = 10000
					)
					save(REspecmod, file = "output/REspecmod.RData", compress = "gzip")



			}

			## upload  models
			for(i in c(1)){
				load(paste0("data/results/models/species.corr.abundance/REspecmod.RData"), verbose = TRUE)
			}

			### check auto correlation
			acorr(REspecmod,summ.zp(REspecmod)$cstats[3])

			### add residuals
			{
				REspecmod <- add.resid(REspecmod, response = "multi"); save(REspecmod, file = "output/REspecmod.RData", compress = "gzip")
			}

			### plot MCMC
			{
				Cairo(file = "output/REspecmod", type = "pdf", units = "in", width = 5, height = 7, family = "regular", pointsize = 7, onefile = TRUE)
				plot(REspecmod)
				dev.off()
			}

			## generate model summaries
			compare <- mod.summ(REspecmod, predictor1 = "RIBLOAD", predictor2 = "ECHINOLOAD"); compare


			}

		### Correlations between Echino and Bd Loads (EBspecmod)
		{

			## extract data from sites where we found
			{ # 50 sites, 1411 individuals
				EB.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 0,			# sites where we checked for Rib in at least this many animals
					found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 1,			# sites where we checked for Bd in at least this many animals
					found.bd = 1,			# sites where we found Bd in at least this many animals
					check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 0,			# sites where we found Ranavirus in at least this many animals
					check.echino = 1,		# sites where we checked for Rib in at least this many animals
					found.echino = 1			# sites where we found Ribeiroia in at least this many animals
				)
			}


			## fit model fits for various variance/covariance structures
			{

				### Both parasites as poisson
				{

					### with individual covariance estimated
					EBspecmod <- multi.test.by.species(EB.use,
						predictor1 = "BDLOAD",	# predictor variable 1
						model1 = "poisson",			# model structure for variable 1
						predictor2 = "ECHINOLOAD",		# predictor variable 2
						model2 = "poisson",			# model structure for variable 2
						site.covar = "yes",			# fit covariance between parasites among sites?
						species.covar = "yes",		# fit covariance between parasites among species?
						ind.covar = "yes",			# fit covariance between parasites among individuals?
						fix = "no",
						nitt = 60000000,
						burnin = 30000000,
						thin = 30000
				)
					save(EBspecmod, file = "output/EBspecmod.RData", compress = "gzip")


				}

			}

			## upload poisson/poisson models
			for(i in c(1)){
				load(paste0("data/results/models/species.corr.abundance/EBspecmod.RData"), verbose = TRUE)
			}

			### check auto correlation
			acorr(EBspecmod,summ.zp(EBspecmod)$cstats[3])

			### add residuals
			{
				EBspecmod <- add.resid(EBspecmod, response = "multi"); save(EBspecmod, file = "output/EBspecmod.RData", compress = "gzip")
			}

			## generate model summaries (takes a long time....)
			compare <- mod.summ(EBspecmod, predictor1 = "BDLOAD", predictor2 = "ECHINOLOAD")

			## plot model results
			plot.model(EBspecmod,predictor1 = "BDLOAD", predictor2 = "ECHINOLOAD",use = EB.use)


			}

		### Correlations between Echino and Ranavirus Loads (ERspecmod)
		{

			## extract data from sites where we found
			{ # 50 sites, 1411 individuals
				ER.use <- coinf.data(dat,site, use.NA = "no",
					check.rib = 0,			# sites where we checked for Rib in at least this many animals
					found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 0,			# sites where we checked for Bd in at least this many animals
					found.bd = 0,			# sites where we found Bd in at least this many animals
					check.rv = 1,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 1,			# sites where we found Ranavirus in at least this many animals
					check.echino = 1,		# sites where we checked for Rib in at least this many animals
					found.echino = 1			# sites where we found Ribeiroia in at least this many animals
				)
			}


			## fit model fits for various variance/covariance structures
			{

				### Both parasites as poisson
				{

					### with individual covariance estimated
					ERspecmod <- multi.test.by.species(ER.use,
						predictor1 = "RVLOAD",	# predictor variable 1
						model1 = "poisson",			# model structure for variable 1
						predictor2 = "ECHINOLOAD",		# predictor variable 2
						model2 = "poisson",			# model structure for variable 2
						site.covar = "yes",			# fit covariance between parasites among sites?
						species.covar = "yes",		# fit covariance between parasites among species?
						ind.covar = "yes",			# fit covariance between parasites among individuals?
						fix = "no",
						nitt = 20000000,
						burnin = 10000000,
						thin = 10000
					)
					save(ERspecmod, file = "output/ERspecmod.RData", compress = "gzip")


				}

			}

			## upload poisson/poisson models
			for(i in c(1)){
				load(paste0("data/results/models/species.corr.abundance/ERspecmod.RData"), verbose = TRUE)
			}

			### check auto correlation
			acorr(ERspecmod,summ.zp(ERspecmod)$cstats[3])

			### add residuals
			{
				ERspecmod <- add.resid(ERspecmod, response = "multi"); save(ERspecmod, file = "output/ERspecmod.RData", compress = "gzip")
			}

			## generate model summaries (takes a long time....)
			compare <- mod.summ.species(ERspecmod, predictor1 = "RVLOAD", predictor2 = "ECHINOLOAD")


			## plot model results
			plot.model(ERspecmod,predictor1 = "RVLOAD", predictor2 = "ECHINOLOAD")


			}

		### plot the results
		{

			## extract data from sites where we found
			{
				RB.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 1,			# sites where we checked for Rib in at least this many animals
					found.rib = 1,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 1,			# sites where we checked for Bd in at least this many animals
					found.bd = 1,			# sites where we found Bd in at least this many animals
					check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 0,			# sites where we found Ranavirus in at least this many animals
					check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
					found.echino = 0		# sites where we found Echinostomes in at least this many animals

				)

				RR.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 1,			# sites where we checked for Rib in at least this many animals
					found.rib = 1,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 0,			# sites where we checked for Bd in at least this many animals
					found.bd = 0,			# sites where we found Bd in at least this many animals
					check.rv = 1,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 1,			# sites where we found Ranavirus in at least this many animals
					check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
					found.echino = 0		# sites where we found Echinostomes in at least this many animals

				)

				BR.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 0,			# sites where we checked for Rib in at least this many animals
					found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 1,			# sites where we checked for Bd in at least this many animals
					found.bd = 1,			# sites where we found Bd in at least this many animals
					check.rv = 1,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 1,			# sites where we found Ranavirus in at least this many animals
					check.echino = 0,		# sites where we checked for Echinostomes in at least this many animals
					found.echino = 0		# sites where we found Echinostomes in at least this many animals

				)

				RE.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 1,			# sites where we checked for Rib in at least this many animals
					found.rib = 1,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 0,			# sites where we checked for Bd in at least this many animals
					found.bd = 0,			# sites where we found Bd in at least this many animals
					check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 0,			# sites where we found Ranavirus in at least this many animals
					check.echino = 1,		# sites where we checked for Echinostomes in at least this many animals
					found.echino = 1		# sites where we found Echinostomes in at least this many animals

				)

				EB.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 0,			# sites where we checked for Rib in at least this many animals
					found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 1,			# sites where we checked for Bd in at least this many animals
					found.bd = 1,			# sites where we found Bd in at least this many animals
					check.rv = 0,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 0,			# sites where we found Ranavirus in at least this many animals
					check.echino = 1,		# sites where we checked for Echinostomes in at least this many animals
					found.echino = 1		# sites where we found Echinostomes in at least this many animals

				)

				ER.use <- coinf.data(dat,site,use.NAs = "no",
					check.rib = 0,			# sites where we checked for Rib in at least this many animals
					found.rib = 0,			# sites where we found Ribeiroia in at least this many animals
					check.bd = 0,			# sites where we checked for Bd in at least this many animals
					found.bd = 0,			# sites where we found Bd in at least this many animals
					check.rv = 1,			# sites where we checked for Ranavirus in at least this many animals
					found.rv = 1,			# sites where we found Ranavirus in at least this many animals
					check.echino = 1,		# sites where we checked for Echinostomes in at least this many animals
					found.echino = 1		# sites where we found Echinostomes in at least this many animals

				)


			}


			### upload models
			{
				load(paste0("data/results/models/species.corr.abundance/RBspecmod.RData"), verbose = TRUE)
				load(paste0("data/results/models/species.corr.abundance/RRspecmod.RData"), verbose = TRUE)
				load(paste0("data/results/models/species.corr.abundance/BRspecmod.RData"), verbose = TRUE)
				load(paste0("data/results/models/species.corr.abundance/REspecmod.RData"), verbose = TRUE)
				load(paste0("data/results/models/species.corr.abundance/EBspecmod.RData"), verbose = TRUE)
				load(paste0("data/results/models/species.corr.abundance/ERspecmod.RData"), verbose = TRUE)
			}


			## generate model summaries
			compare <- rbind(mod.summ.species(RBspecmod, predictor1 = "RIBLOAD", predictor2 = "BDLOAD"),
				mod.summ.species(RRspecmod, predictor1 = "RIBLOAD", predictor2 = "RVLOAD"),
				mod.summ.species(BRspecmod, predictor1 = "BDLOAD", predictor2 = "RVLOAD"),
				mod.summ.species(REspecmod, predictor1 = "RIBLOAD", predictor2 = "ECHINOLOAD"),
				mod.summ.species(EBspecmod, predictor1 = "BDLOAD", predictor2 = "ECHINOLOAD"),
				mod.summ.species(ERspecmod, predictor1 = "RVLOAD", predictor2 = "ECHINOLOAD"))	;compare[order(compare$species),]

			### individual level correlations
			plot.ind.corr.species(RBspecmod,RRspecmod,BRspecmod,REspecmod,EBspecmod,ERspecmod, RB.use,RR.use,BR.use, RE.use, EB.use, ER.use, sample.size = 1000, remove = c("",""))

			### dot and line plot of correlation posteriors
			{


				## generate data
				{
					## generate model summaries (takes a long time....)
					compare <- rbind(mod.summ(RBmod1,predictor1 = "RIBLOAD", predictor2 = "BDLOAD"),
											mod.summ(RRmod1,predictor1 = "RIBLOAD", predictor2 = "RVLOAD"),
											mod.summ(BRmod1,predictor1 = "BDLOAD", predictor2 = "RVLOAD")
					)

					## add combination column
					compare$combo <- paste(compare$predictor1, compare$predictor2, sep = "_")
					compare$combo <- factor(compare$combo)
					levels(compare$combo) <- c("R. ondatrae vs.\\nB. dendrobatidis","R. ondatrae vs.\\nRanavirus", "B. dendrobatidis vs.\\nRanavirus")

					## melt data farme
					compare.melted <- rbind() # new data frame
					for(i in c(1:nrow(compare))){
						new_i <- data.frame(combo = rep(compare[i,"combo"],3),
							level = c("site","species","individual"),
							corr = c(compare[i,"site.cor"],compare[i,"species.cor"],compare[i,"ind.cor"]),
							low = c(compare[i,"site.cred.low"],compare[i,"species.cred.low"],compare[i,"ind.cred.low"]),
							hi = c(compare[i,"site.cred.hi"],compare[i,"species.cred.hi"],compare[i,"ind.cred.hi"]),
							strong = c(compare[i,"site.strong"],compare[i,"species.strong"],compare[i,"ind.strong"]))
						compare.melted <- rbind(compare.melted,new_i)
					}
				}

				## plot
				Cairo(file = "output/all_level_correlations", type = "svg", units = "in", width = 10, height = 4, family = "regular", pointsize = 12)
				pl <- ggplot(data = compare.melted, aes(level, corr)) +
						theme_bw() +
						coord_flip() +
						#scale_fill_manual(values = c("white","black")) +
						geom_hline(y = 0, lty = 3) +
						labs(y = expression(paste("correlation (",rho,")"))) +
						geom_errorbar(aes(ymin = low,ymax = hi), width = 0.1) +
						geom_point(size = 4, pch = 21, fill = "black") +
						facet_grid(.~combo) +
						theme(legend.text = element_text(face = "italic"),
								axis.title.y = element_blank(),
							strip.text = element_text(face = "italic", size = 17),
							axis.title.x = element_text(size = 17),
							axis.text = element_text(size = 16),

							strip.background = element_blank())
				print(pl)
				dev.off()

			}

			## tables of posteriors
			{

					output.posteriors(RBspecmod,"RBspecmod")
					output.posteriors(RRspecmod,"RRspecmod")
					output.posteriors(BRspecmod,"BRspecmod")
					output.posteriors(REspecmod,"REspecmod")
					output.posteriors(EBspecmod,"EBspecmod")
					output.posteriors(ERspecmod,"ERspecmod")

			}
		}
}
