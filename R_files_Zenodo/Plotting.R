###===========================================================================
###===========================================================================
### Replication code to reproduce plots and analyses in main text
###===========================================================================
###===========================================================================

### Loading packages
library("ggplot2") # package for plotting
library("readr")   # package for reading tab/comma/arbitrarily-delimited data ("spreadsheet data")
library("sp")        # package for interacting with spatial data in R
library("sf")        # package for interacting with spatial data in a tidier format in R
library("albersusa") # package that includes a dataset projecting the US into a cleaner plotting window (projects Alaska and Hawaii into plot range)

###===================================
### 1: Reading in replication datasets
###===================================

mean_expressions_DF <- readr::read_tsv("Mean21Expressions.tsv")
  # dataset showing mean viewpoint with standard error of the mean (SEM) for all 21 environmental issues and 6 personas
US_geography_DF <- readr::read_tsv("US_geography.tsv")
  # dataset with simple features geometry and state-level ranking for each persona
political_ideology_DF <- readr::read_tsv("Persona_PoliticalIdeology.tsv")
  # dataset with mean and SEM political ideology for each persona

# Note that in these datasets, "Clust" is shorthand for "Cluster", which is equivalent to the term we use in the main text, "Persona"
# (That is, each cluster constituted its own unique persona.)

###===================================
### 2: Plotting distribution of
### expressions
###===================================
limits = aes(ymin=MinSE,ymax=MaxSE) # specify the limits for the error bar using ggplot2's "aes" convention

p_views = ggplot(mean_expressions_DF,aes(x=Clust,y=mean,color=Clust)) # initiate a plot called p_views
p_views = p_views + geom_point(size=1.2) +
    geom_errorbar(limits,width=0.5,size=0.4) +
    facet_wrap(.~variable,ncol=3,scales="free")
p_views = p_views + labs(x="",y="Mean issue viewpoint")
p_views # display plot

###===================================
### 3: Visualizing US geographic
### distribution
###===================================

p_US <- ggplot(data = US_geography_DF,aes(fill=Rank,group=Persona)) # initiate a plot called p_US
p_US <- p_US + geom_sf(aes(fill=Rank),color="gray47",size=0.15) +
      facet_wrap(~.Persona,ncol = 2)
p_US

###===================================
### 4: Plotting political ideology
###===================================

plot(mean~Clust,political_ideology_DF,pch=19)
