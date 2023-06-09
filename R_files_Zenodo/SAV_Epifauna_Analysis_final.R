#### CHESAPEAKE BAY SAV EPIFAUNAL ANALYSIS ####

# Authors: Claire Murphy, Jon Lefcheck 
# Last updated: 26 August 2020
# Contact: cemurphy@ucdavis.edu


#####-----------------------------------------------------------

# Load required libraries
library(cowplot)
library(ecodist)
library(ggmap)
library(ggthemes)
library(lme4)
library(tidyverse)
library(vegan)
library(nlme)
library(colorspace)
library(indicspecies)
library(OneR)

# Get Edgar equations
source("https://gist.githubusercontent.com/jslefche/c480eeb2ad67ca7c9a46/raw/177e47862e4053337cca66e94f5f0b91d5e7086f/edgarMethod.R")

#####-----------------------------------------------------------

# Read in data
bugs <- read.csv("epifaunal data 6-30-19.csv")

macro <- read.csv("macrophytes 6-13-19.csv")

temp <- read.csv("temperature.csv")

coords <- read.csv("coordinates.csv")

cover <- read.csv("cover.csv")

distances <- read.csv("transect distances_new.csv")

colnames(distances) <- gsub("\\\\.", " ", colnames(distances))

salinity <- read.csv("salinity.csv")

#####-----------------------------------------------------------

#rename "bugs" columns 
colnames(bugs) <- c(colnames(bugs)[1:8], "5.6", "4.0", "2.8", "2.0", "1.4", "1.0", "0.71", "0.5","notes")

# Curate bug data and scale by macrophyte biomass

# Collapse Palaemonetes
bugs <- bugs %>% mutate(species.name.revised = recode(species.name.revised, "Palaemonetes vulgaris" = "Palaemonetes sp.", "Palaemonetes intermedius" = "Palaemonetes sp.","Palaemonetes pugio" = "Palaemonetes sp."))

# Collapse unknown caprellid
bugs <- bugs %>% mutate(species.name.revised = recode(species.name.revised, "unknown caprellid" = "Caprella penantis"))

# Fix formatting of "Genus sp."s
bugs <- bugs %>% mutate(species.name.revised = recode(species.name.revised,"Batea sp" = "Batea sp.", "Corophium sp" = "Corophium sp.", "Nereis sp" = "Nereis sp."))

# Change "Nereis sp." to "Alitta sp."
bugs <- bugs %>% mutate(species.name.revised = recode(species.name.revised, "Nereis sp." = "Alitta sp."))

# Change "Alitta sp." to "Alitta succinea"
bugs <- bugs %>% mutate(species.name.revised = recode(species.name.revised, "Nereis sp." = "Alitta succinea"))

# Remove unknowns and non-target species 
remove <- c("unknown amphipod", "unknown mollusc", "unknown caprellid", "Mysida sp", "Chthamalus fragilis", "megalopa", "Ostracoda sp", "unknown shrimp", "Argopecten irradians", "Polychaeta sp", "Alpheus sp", "unknown gastropod", "copepod")

bugs <- bugs[!bugs$species.name.revised %in% remove, ] 

# Collapse species IDs
bugs <- bugs %>% group_by(sample.id, date, location, replicate, species.name.revised) %>% summarise_at(vars(`5.6`:`0.5`), sum)
            
# Assign group
bugs$group <- ifelse(bugs$species.name.revised %in% c("Bittolium varium", "Crepidula fornicata", "Haminoea solitaria", "Mitrella lunata", "Odostomia bisaturalis", "Polycerella conyma", "Anachis avara"),"Mollusc", 
       ifelse(bugs$species.name.revised %in% c("Paracaprella tenuis", "Caprella penantis", "Caprella equilibra"), "Caprellid", 
       ifelse(bugs$species.name.revised %in% c("Alitta succinea"), "Polychaete", "Crustacean")))

# Compute biomass
bugs <- edgarMethod(bugs, group.colname = "group")

# Average tin weight for epiphytes
avg.tin <- macro %>% filter(species.name.revised == "epiphytes") %>% summarize(mean(tin, na.rm = T))

macro$tin <- unlist(ifelse(is.na(macro$tin), avg.tin, macro$tin))

# any(is.na(macro$tin))

# get rid of letters on the site codes
macro$sample.id <- as.factor(gsub("(.*)[a-z]", "\\\\1", macro$sample.id))

# Scale the abundance and biomass by the total macrophyte biomass
macro <- macro %>% 
  
  # remove roots
  filter(!species.name.revised %in% c("root")) %>%
  
  # compute dry mass
  transform(dry.mass = tin.and.dry - tin) %>%
  
  group_by(sample.id) %>%
  
  # summarize masses
  summarize(
    seagrass.g = sum(dry.mass[species.name.revised %in% c("Zostera marina", "Ruppia maritima")], na.rm = T),
    epiphytes.g = sum(dry.mass[species.name.revised %in% c("epiphytes", "Hydroid sp.")], na.rm = T),
    macro.g = sum(dry.mass[!species.name.revised %in% c("Zostera marina", "Ruppia maritima", "epiphytes", "Hydroid sp.")], na.rm = T)
    ) 

# merge with the bug dataset
bugs <- left_join(bugs, macro, by = "sample.id")

# Remove any observation missing seagrass biomass
bugs <- subset(bugs, seagrass.g != 0)

# bring in lat/longs
names(coords)[2] <- "location"

bugs <- left_join(bugs, coords, by = "location")

# remove sites that didn't have transects (Coastal Bays + Goose & Goodwin Islands)
remove.loc <- c("South Bay", "Cobb Bay", "Goose Island", "Goodwin")

bugs <- bugs %>% filter(!location %in% remove.loc) 

# remove any species with <= 5 individuals across the whole dataset
remove <- bugs %>% group_by(species.name.revised) %>% summarize(abund = sum(total.abund)) %>%
  filter(abund <= 5) %>% .$species.name.revised

bugs <- bugs %>% filter(!species.name.revised %in% remove)

# scale abundance and biomass by total seagrass
bugs[, grepl("abund|biomass", colnames(bugs))] <- bugs[, grepl("abund|biomass", colnames(bugs))] / bugs[, c("seagrass.g")]

#####-----------------------------------------------------------

# Create community composition dissimilarity matrices and further curate dataset

# summarize by location
bugs.summary <- bugs %>%
  
  group_by(location, Lat, Long) %>%
  
  summarize(
    abund = mean(total.abund),
    biomass = mean(total.biomass),
    rich = length(unique(species.name.revised))
  )

### Create sample-by-sample distance matrix - uncomment the type of community composition index you would like to use for the analysis and comment out the other 2


# #using biomass and Bray Curtis distances
#
# # Cast site by species biomass matrix
bugs.mat <- bugs %>%

  group_by(location, species.name.revised) %>%

  summarise(total.biomass = mean(total.biomass)) %>%

  spread(species.name.revised, total.biomass, fill = 0)


bugs.dist <- vegdist(bugs.mat[,2:ncol(bugs.mat)])

# #using abundance and Bray Curtis distances
# bugs.mat <- bugs %>%
# 
#   group_by(location, species.name.revised) %>%
# 
#   summarise(total.abund = mean(total.abund)) %>%
# 
#   spread(species.name.revised, total.abund, fill = 0)
# 
# bugs.dist <- vegdist(bugs.mat[,2:ncol(bugs.mat)])

#using Jaccard and presence/absesnce
# bugs.dist <- vegdist(bugs.mat[,2:ncol(bugs.mat)], binary = TRUE, method = "jaccard")


# compute average dissimilarity
bugs.summary <- left_join(bugs.summary, data.frame(location = bugs.mat$location, dissim = colMeans(as.matrix(bugs.dist))))

#####-----------------------------------------------------------

# Run NMDS

bugs.nmds <- metaMDS(bugs.mat[, 2:ncol(bugs.mat)])

# Extract output and bind with metadata
bugs.nmds.df <- cbind.data.frame(location = bugs.mat$location, bugs.nmds$points)

# assign regions
bugs.nmds.df$Region <- with(bugs.nmds.df, 
    ifelse(location %in% c("Hungars Creek 1", "Hungars Creek 2"), "Eastern Shore\\n",
            ifelse(location %in% c("Poquoson 1", "Poquoson 2", "Back River"), "Poquoson\\n",
                     ifelse(location %in% c("Allens Island", "Guinea"), "York\\n",
                     ifelse(location %in% c("Browns Bay 1", "Browns Bay 2"), "Browns Bay\\n",
                     ifelse(location %in% c("East River 1", "East River 2"), "East River\\n",
                     ifelse(location %in% c("Four Points South", "Ware River"), "Four Points\\n",
                     ifelse(location == "Dameron Marsh", "Dameron\\n", "Mobjack\\n"))))))))

# Run NMDS with each individual sample (not collapsed by site first)
bugs.sample.id.mat <- bugs %>% 
  
  group_by(sample.id, location, species.name.revised) %>%
  
  summarise(total.biomass = mean(total.biomass)) %>%
  
  spread(species.name.revised, total.biomass, fill = 0)


bugs.nmds.id <- metaMDS(bugs.sample.id.mat[, 3:ncol(bugs.sample.id.mat)])



# Extract output and bind with metadata
bugs.nmds.df.id <- cbind.data.frame(location = bugs.sample.id.mat$location, bugs.nmds.id$points)

# assign regions
bugs.nmds.df.id$Region <- with(bugs.nmds.df.id, 
                               ifelse(location %in% c("Hungars Creek 1", "Hungars Creek 2"), "Eastern Shore\\n",
                                      ifelse(location %in% c("Poquoson 1", "Poquoson 2", "Back River"), "Poquoson\\n",
                                             ifelse(location %in% c("Allens Island", "Guinea"), "York\\n",
                                                    ifelse(location %in% c("Browns Bay 1", "Browns Bay 2"), "Browns Bay\\n",
                                                           ifelse(location %in% c("East River 1", "East River 2"), "East River\\n",
                                                                  ifelse(location %in% c("Four Points South", "Ware River"), "Four Points\\n",
                                                                         ifelse(location == "Dameron Marsh", "Dameron\\n", "Mobjack\\n"))))))))

#Spider plot by site
cent <- bugs.nmds.df.id %>%
  
  group_by(location, Region) %>%
  
  summarise(centMDS1 = as.numeric(mean(MDS1)), centMDS2 = as.numeric(mean(MDS2)))

#add the centroids back to the individual NMDS points
bugs.nmds.df.id.Sit <- merge(bugs.nmds.df.id, cent,
                             by = c("location", "Region"), sort = FALSE)

# Plot results
ggplot(bugs.nmds.df.id.Sit, aes(color = Region)) +
  scale_colour_manual(values = qualitative_hcl(8, "dark3")) +
  geom_point(aes(x = centMDS1, y = centMDS2,  shape = Region), size = 5, stroke = 1.5) +
  scale_shape_manual(values=c(0:3, 5, 6, 8, 13)) +
  geom_point(aes(x = MDS1, y = MDS2), size = 3, alpha = 0.25, show.legend = FALSE) +
  geom_segment(aes(x = MDS1, y = MDS2, xend = centMDS1, yend = centMDS2), 
               alpha = 0.25,show.legend = FALSE) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "NMDS1", y = "NMDS2") +
  annotate(geom = "text", label = paste("Stress =", round(bugs.nmds$stress, 3)), 
           x = -1.0, y = -0.97, size = 4.5) 


#####-----------------------------------------------------------

# Run MRMs to explain community dissimilarity matrices

# get macro matrix
macro.mat <- bugs %>% 
  
  group_by(location) %>%
  
  summarise(total.biomass = sum(epiphytes.g + macro.g)/length(unique(sample.id)))

macro.mat <- macro.mat[order(macro.mat$location), ]

macro.dist <- dist(macro.mat$total.biomass)

# get distances matrix
distances.sub <- distances[, colnames(distances) %in% c("Transect", unique(bugs$location))]

distances.sub <- subset(distances.sub, Transect %in% unique(bugs$location))

distances.dist <- as.dist((distances.sub[, -1]))


# get temperature distance matrix
temp.sub <- subset(temp, location %in% unique(bugs$location))

temp.dist <- dist(temp.sub[order(temp.sub$location), ]$temp)

# get cover distance matrix
cover.sub <- cover %>% 
  
  filter(location %in% unique(bugs$location)) %>% 
  
  group_by(location) %>% 
  
  summarize(cover = median(Cover), Dist = max(Dist, na.rm = T)) %>%
  
  filter(cover != 0)

cover.dist <- dist(cover.sub$cover)

# get date distance matrix
date.sub <- bugs %>%
  
  group_by(location) %>%
  
  subset(select = c(location, date)) %>%
  
  summarise(date = mean(as.Date(date, "%m/%d/%y")))

date.sub <- date.sub[order(date.sub$location), ]

date.dist <- dist(date.sub$date)

#get salinity matrix
salinity.sub.mat <- salinity %>% 
  
  filter(location %in% unique(bugs$location)) %>%
  
  filter(wq.depth == 1) %>% #include only shallow water measurements
  
  group_by(location) %>%

  summarise(salinity = mean(salinity))

salinity.sub.mat <- salinity.sub.mat[order(salinity.sub.mat$location), ]

salinity.dist <- dist(salinity.sub.mat$salinity)

# get evenness matrix (using eq 5 of Jost 2010)
even.mat <- as.matrix(log(exp(diversity(bugs.mat[,2:ncol(bugs.mat)])))/
                     log(specnumber(bugs.mat[,2:ncol(bugs.mat)])), 
                     ncol = 1, nrow = 16)


even.mat <- cbind(bugs.summary$location, even.mat)

even.mat <- even.mat[order(even.mat[,1]), ]

even.dist <- dist(even.mat[,2])

# get dissimilarity matrices for abundance, biomass, and richness
bugs.summary <- bugs.summary[order(bugs.summary$location), ]

abund.dist <- dist(bugs.summary$abund)

biom.dist <- dist(bugs.summary$biomass)

rich.dist <- dist(bugs.summary$rich)

# check dims
dim(as.matrix(abund.dist))
dim(as.matrix(biom.dist))
dim(as.matrix(rich.dist))
dim(as.matrix(even.dist))
dim(as.matrix(bugs.dist)) 
dim(as.matrix(distances.dist)) 
dim(as.matrix(temp.dist))
dim(as.matrix(cover.dist))
dim(as.matrix(macro.dist))
dim(as.matrix(date.dist))
dim(as.matrix(salinity.dist))

# regress against other distance matrices
abund.MRM <- MRM(abund.dist ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

biom.MRM <- MRM(biom.dist ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

rich.MRM <- MRM(rich.dist ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

even.MRM <- MRM(even.dist ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

comm.MRM <- MRM(bugs.dist ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

abund.MRM
biom.MRM
rich.MRM
comm.MRM
even.MRM


#####-----------------------------------------------------------
# Look at role of development type 

#load data
dev.type <- read.csv("DevType.csv")
#add to all bug data
bugs.withDev <- left_join(bugs, dev.type, by = "species.name.revised")

#separate out 2 dev types
direct <- bugs.withDev[bugs.withDev$Development.type == "Direct",]
plankton <- bugs.withDev[bugs.withDev$Development.type == "Planktonic",]

#### Rerun MRMs with just direct developers
direct.summary <- direct %>%
  
  group_by(location, Lat, Long) %>%
  
  summarize(
    abund = mean(total.abund),
    biomass = mean(total.biomass),
    rich = length(unique(species.name.revised))
  )

direct.mat <- direct %>% 
  
  group_by(location, species.name.revised) %>%
  
  summarise(total.biomass = mean(total.biomass)) %>%
  
  spread(species.name.revised, total.biomass, fill = 0)


# get evenness matrix (using eq 5 of Jost 2010)
even.direct <- as.matrix(log(exp(diversity(direct.mat[,2:ncol(direct.mat)])))/
                        log(specnumber(direct.mat[,2:ncol(direct.mat)])), 
                      ncol = 1, nrow = 16)


even.direct <- cbind(direct.summary$location, even.direct)

even.direct <- even.direct[order(even.direct[,1]), ]

even.direct.dist <- dist(even.direct[,2])

# get dissimilarity matrices for abundance, biomass, and richness
direct.summary <- direct.summary[order(direct.summary$location), ]

abund.direct <- dist(direct.summary$abund)

biom.direct <- dist(direct.summary$biomass)

rich.direct <- dist(direct.summary$rich)

direct.dist <- vegdist(direct.mat[,2:ncol(direct.mat)])



# regress against other distance matrices
abund.direct.MRM <- MRM(abund.direct ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

biom.direct.MRM <- MRM(biom.direct ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

rich.direct.MRM <- MRM(rich.direct ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

even.direct.MRM <- MRM(even.direct.dist ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

comm.direct.MRM <- MRM(direct.dist ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

abund.direct.MRM
biom.direct.MRM
rich.direct.MRM
comm.direct.MRM
even.direct.MRM


#### Rerun MRMs with just planktonic developers

plankton.summary <- plankton %>%
  
  group_by(location, Lat, Long) %>%
  
  summarize(
    abund = mean(total.abund),
    biomass = mean(total.biomass),
    rich = length(unique(species.name.revised))
  )

plankton.mat <- plankton %>% 
  
  group_by(location, species.name.revised) %>%
  
  summarise(total.biomass = mean(total.biomass)) %>%
  
  spread(species.name.revised, total.biomass, fill = 0)


# get evenness matrix (using eq 5 of Jost 2010)
even.plankton <- as.matrix(log(exp(diversity(plankton.mat[,2:ncol(plankton.mat)])))/
                           log(specnumber(plankton.mat[,2:ncol(plankton.mat)])), 
                         ncol = 1, nrow = 16)


even.plankton <- cbind(plankton.summary$location, even.plankton)

even.plankton <- even.plankton[order(even.plankton[,1]), ]

even.plankton.dist <- dist(even.plankton[,2])

# get dissimilarity matrices for abundance, biomass, and richness
plankton.summary <- plankton.summary[order(plankton.summary$location), ]

abund.plankton <- dist(plankton.summary$abund)

biom.plankton <- dist(plankton.summary$biomass)

rich.plankton <- dist(plankton.summary$rich)

plankton.dist <- vegdist(plankton.mat[,2:ncol(plankton.mat)])



# regress against other distance matrices
abund.plankton.MRM <- MRM(abund.plankton ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

biom.plankton.MRM <- MRM(biom.plankton ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

rich.plankton.MRM <- MRM(rich.plankton ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

even.plankton.MRM <- MRM(even.plankton.dist ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

comm.plankton.MRM <- MRM(plankton.dist ~ distances.dist + temp.dist + cover.dist + macro.dist + salinity.dist + date.dist, nperm = 1e5)

abund.plankton.MRM
biom.plankton.MRM
rich.plankton.MRM
comm.plankton.MRM
even.plankton.MRM


#####-----------------------------------------------------------

# Explore which species are driving cover pattern


# see what the range of cover values is
unique(cover.sub$cover)

# check species names
unique(bugs$species.name.revised)

# make matrix that has a row for every species every time they occur at a site, with the biomass and cover
cov.bio <- bugs %>%
  
  group_by(location, species.name.revised) %>% 
  
  subset(select = c(location, species.name.revised, total.biomass)) %>%
  
  summarise(total.biomass = mean(total.biomass))

cov.bio <- merge(cov.bio, cover.sub, by = "location")

## Using indicspecies

# Create cover bins using the Braun-Blanquet scale
# 1 = < 5 %, 2 = 5-25%, 3 = 26-50%, 4 = 51-75%, 5 = 76-100%
cov.bins <- c(3, 3, 1, 2, 3, 2, 4, 2, 4, 3, 4, 4, 4, 4, 4, 3)

indval <- multipatt(bugs.mat[,2:ncol(bugs.mat)], cov.bins, max.order = 1, control = how(nperm = 999))

summary(indval)


#####-----------------------------------------------------------
# Run MRMs excluding cover associated species (Erichsonella, Elasmopus, Paracaprella)

# get bug matrix
bugs.cov <- bugs %>% 
  
  filter(!species.name.revised %in% c("Elasmopus levis", "Erichsonella attenuata", "Paracaprella tenuis")) %>%
  
  group_by(location, species.name.revised) %>%
  
  summarise(total.biomass = mean(total.biomass)) %>%
  
  spread(species.name.revised, total.biomass, fill = 0)

bugs.cov <- bugs.cov[order(bugs.cov$location), ]

bugs.cov.dist <- vegdist(bugs.cov[,2:ncol(bugs.cov)])

# get evenness matrix (using eq 5 of Jost 2010)
even.cov.mat <- as.matrix(log(exp(diversity(bugs.cov[,2:ncol(bugs.cov)])))/
                        log(specnumber(bugs.cov[,2:ncol(bugs.cov)])), 
                      ncol = 1, nrow = 16)

even.cov.mat <- cbind(bugs.summary$location, even.cov.mat)

even.cov.mat <- even.cov.mat[order(even.cov.mat[,1]), ]

even.cov <- dist(even.cov.mat[,2])


# summarize by location
bugs.cov.sum <- bugs %>%
  
  filter(!species.name.revised %in% c("Erichsonella attenuata", "Elasmopus levis", "Paracaprella tenuis")) %>%
  
  group_by(location, Lat, Long) %>%
  
  summarize(
    abund = mean(total.abund),
    biomass = mean(total.biomass),
    rich = length(unique(species.name.revised))
  ) 

bugs.cov.sum <- bugs.cov.sum[order(bugs.cov.sum$location), ]

# get dissimilarity matrices for abundance, biomass, and richness
abund.cov <- dist(bugs.cov.sum$abund)

biom.cov <- dist(bugs.cov.sum$biomass)

rich.cov <- dist(bugs.cov.sum$rich)

# run new MRMs
abund.MRM.cov <- MRM(abund.cov ~ distances.dist + temp.dist + cover.dist + macro.dist + date.dist + salinity.dist, nperm = 1e5)

biom.MRM.cov <- MRM(biom.cov ~ distances.dist + temp.dist + cover.dist + macro.dist + date.dist + salinity.dist, nperm = 1e5)

rich.MRM.cov <- MRM(rich.cov ~ distances.dist + temp.dist + cover.dist + macro.dist + date.dist + salinity.dist, nperm = 1e5)

comm.MRM.cov <- MRM(bugs.cov.dist ~ distances.dist + temp.dist + cover.dist + macro.dist + date.dist + salinity.dist, nperm = 1e5)

even.MRM.cov <- MRM(even.cov ~ distances.dist + temp.dist + cover.dist + macro.dist + date.dist + salinity.dist, nperm = 1e5)


abund.MRM.cov
biom.MRM.cov
rich.MRM.cov
comm.MRM.cov
even.MRM.cov

#####-----------------------------------------------------------
## Make a map


#make a new formatted table to plot each location with the same color and symbol as in the NMDS
map.coords <- bugs.summary[,1:3]
# assign regions
map.coords$Region <- with(map.coords, 
                          ifelse(location %in% c("Hungars Creek 1", "Hungars Creek 2"), "Eastern Shore\\n",
                                 ifelse(location %in% c("Poquoson 1", "Poquoson 2", "Back River"), "Poquoson\\n",
                                        ifelse(location %in% c("Allens Island", "Guinea"), "York\\n",
                                               ifelse(location %in% c("Browns Bay 1", "Browns Bay 2"), "Browns Bay\\n",
                                                      ifelse(location %in% c("East River 1", "East River 2"), "East River\\n",
                                                             ifelse(location %in% c("Four Points South", "Ware River"), "Four Points\\n",
                                                                    ifelse(location == "Dameron Marsh", "Dameron\\n", "Mobjack\\n"))))))))

#get a map of the lower Chesapeake Bay region to use as a background 

cb_map <- ggmap(get_stamenmap(bbox = c(-76.6, 37, -75.8, 37.9), maptype = "terrain-background"))


#map it all
cb_map + 
  geom_point(data = map.coords, aes(x = Long, y = Lat, color = Region, shape = Region),
             stroke = 2, size = 3) +
  scale_colour_manual(values = qualitative_hcl(8, "dark3")) +
  scale_shape_manual(values=c(0:3, 5, 6, 8, 13)) +
  theme_bw(base_size = 15) +
  labs(x = "Longitude", y = "Latitude") 


## Get map of Mid Atlantic for the inset
cb <- get_stamenmap(bbox = c(-80, 32.5, -73, 42.5), maptype = "terrain")

#filter out only the lats and longs we're interested in
states <- map_data("state")

states <- states[(states$long > -80),]
states <- states[(states$long < -73),]
states <- states[(states$lat > 32.5),]
states <- states[(states$lat < 42.5),]



#plot mid-atlantic with state outlines
ggmap(cb) +
  geom_polygon( data=states, aes(x=long, y=lat, group = group),
                color="black", fill="NA" )