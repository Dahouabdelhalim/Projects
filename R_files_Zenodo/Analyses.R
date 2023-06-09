###===========================================================================
###===========================================================================
### Replication code to reproduce analyses in main text
###===========================================================================
###===========================================================================

### Loading packages
library("ggplot2") # package for plotting
library("readr")   # package for reading tab/comma/arbitrarily-delimited data ("spreadsheet data")
library("sp")      # package for interacting with spatial data in R
library("sf")      # package for interacting with spatial data in a tidier format in R
library("lme4")    # package for performing mixed effect regression models
library("lmerTest")# package that extends the functionality of lme4

###===================================
### 1: Reading in replication datasets
###===================================

issueIntensities <- readr::read_tsv("IssueIntensities.tsv")
  # dataset containing mean issue intensity by country

CONUSdf <- sf::st_read("CONUSprominence.geojson")
  # dataset containing issue prominence by state with state geography in simple features format

countyDF <- readr::read_tsv("countyProminence.tsv")
  # dataset with county-level information

###===================================
### 2: Plotting intensities of
### issues across countries
###===================================
### Inputs for plot
limits <- aes(ymin=MinSE,ymax=MaxSE) # specify the limits for the error bar using ggplot2's "aes" convention
issue_labeller <- c("Agriculture","Animal Welfare","Climate action",
                    "Climate belief","Climate policy","Renewable energy",
                    "Diets & consumer goods","CSR","Fossil fuels",
                    "Habitats & species","Marine","Birdwatching",
                    "Gardening","Hiking","Hunting & angling",
                    "Environmental policy","Pollution","Public lands",
                    "Transit infrastructure","Freshwater","Extreme weather")
names(issue_labeller) <- sort( unique( issueIntensities$variable ) )
xlabels <- issue_labeller[unique( issueIntensities$variable )]

### Jittering points
set.seed(8)
pos <- position_dodge(width=0.9, preserve="total")

### Generating plot
p_intensity <-  ggplot(issueIntensities,aes(x=varnum,y=mean,color=Country,shape=Country)) # initiate a plot called p_views
p_intensity <-  p_intensity + geom_point(size=1.5,position=pos) +
  geom_errorbar(limits,width=0.5,size=0.4,position=pos)
p_intensity <-  p_intensity + labs(x="",y="Intensity")
p_intensity <- p_intensity + scale_x_continuous(breaks=rev(c(1:21)),labels=xlabels)
p_intensity <- p_intensity + coord_flip() + theme_classic()
p_intensity # display plot

###===================================
### 3: Visualizing US geographic
### distribution
###===================================

p_US <- ggplot(data = CONUSdf,aes(fill=Prominence,group=variable)) # initiate a plot called p_US
p_US <- p_US + geom_sf(aes(fill=Prominence,geometry=geometry),color="gray47",size=0.15) +
  facet_wrap(vars(variable),ncol = 3,labeller=labeller(variable=issue_labeller))
p_US <- p_US + theme_void()
p_US

###===================================
### 4: Constructing regression model
### for issue prominence at county-scale
###===================================

### Fitting models
ag_model <- lmerTest::lmer(ag~polIdeo + broadband + rurban + region+ (1|state), data=countyDF) # agriculture
clim_act_model <- lmerTest::lmer(clim_act~polIdeo + broadband + rurban + region+ (1|state), data=countyDF) # climate action / public mobilization for decarbonization
hunt_model <- lmerTest::lmer(out_hunt_ang~polIdeo + broadband + rurban + region+ (1|state), data=countyDF) # hunting and angling
pub_model <- lmerTest::lmer(pub~polIdeo + broadband + rurban + region+ (1|state), data=countyDF) # public lands

### Model outputs
summary( ag_model ) 
summary( clim_act_model ) 
summary( hunt_model ) 
summary( pub_model ) 
