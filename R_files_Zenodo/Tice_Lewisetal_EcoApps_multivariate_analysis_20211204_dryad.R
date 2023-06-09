library(lubridate)
library(ggplot2)
library(ggpmisc)
library(scales)
library(gridExtra)
library(plyr)
library(dplyr)
library(scales)
library(readr)
library(zoo)
library(tidyr) 
library(data.table)
library(dtplyr)
library(reshape2)
require(grid)
#cowplot theme
require(cowplot) 

##### import csv ####
#setwd
setwd("~/Dropbox/Lab Projects/Oyster Projects/oyster resilence 2013 - 2015/New Versions/CRFL-Wells Comparison")


###Load crfl.wells from csv ####
crfl.wells <- #Load .csv file, "
crfl.wells <- crfl.wells[ ,-c(1)]

library(vegan)
##drop sites w/out salinity data (mean by season/year) 
crfl.wells <- crfl.wells %>% drop_na(Mean)
crfl.wells.env <- crfl.wells %>% select(1:35)
crfl.wells.fam <- crfl.wells %>% select(36:92)

## separate data based on study for analyses on individual studies
crfl.fam <- crfl.wells.fam %>% slice(c(36:81))
crfl.env <- crfl.wells.env %>% slice(c(36:81))
wells.fam <- crfl.wells.fam %>% slice(c(1:35))
wells.env <- crfl.wells.env %>% slice(c(1:35))

###Create a dissimilary matric based on B-C
cwf.dist <- vegdist(crfl.wells.fam, method = "jaccard", binary=T)
print(cwf.dist)
summary(cwf.dist)
plot(cwf.dist)

##Generate jaccard distance matrix of all faunal data
jaccard.fam <- vegdist(crfl.wells.fam, method = "jaccard")

##generate a distance matrix of all salinity differences 
saldiff.cw <- vegdist(crfl.wells.env$Mean, method = "euclidean",na.rm = TRUE)

##perform a mantel test to test the correlation between salinity differences and jaccard distances
mantel(saldiff.cw, jaccard.fam) 

##Generate jaccard distance matrix of crfl data
jaccard.crfl <- vegdist(crfl.fam)

##generate jaccard distance matrix of wells data
jaccard.wells <- vegdist(wells.fam)



## Run function to make a list of community dissimilarities with columns that will be used for B-diversity comparisons between studies 
## show comparison being made. From: https://stackoverflow.com/questions/23474729/convert-object-of-class-dist-into-data-frame-in-r
dist.df <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
}

##Create column lists of dissilimilarity measures with corresponding sample comparisons
dist.list.fam <- dist.df(jaccard.fam)
dist.list.crfl <- dist.df(jaccard.crfl)
dist.list.wells <- dist.df(jaccard.wells)

##inlet distance
inlet.crfl <- crfl.env[ ,c(11)]
inlet.dist <- vegdist(inlet.crfl, method = "euclidean")
inlet.dist.crfl <- dist.df(inlet.dist)
dist.env <- bind_cols(inlet.dist.crfl, dist.list.crfl)
names(dist.env)[3] <- "inlet.dist"
dist.env <- dist.env[ ,-c(4,5)]
names(dist.env)[4] <- "jacc.diss"

##seasonal mean salinity
salseas.crfl <- crfl.env[ ,c(32)]
salseas.dist <- vegdist(salseas.crfl, method = "euclidean")#make euclidean distance matrix of refractometer differences          
salseas.dist.crfl <- dist.df(salseas.dist)
dist.env <- bind_cols(salseas.dist.crfl, dist.env)
dist.env <- dist.env[ ,-c(4,5)]
names(dist.env)[3] <- "seassal.diff"


##create factors from distances
dist.env$site.comp <- NA 
dist.env$site.comp[dist.env$inlet.dist==0] <- NA
dist.env$site.comp[dist.env$inlet.dist>0 & dist.env$inlet.dist<2] <- "CR.WRR"
dist.env$site.comp[dist.env$inlet.dist>2 & dist.env$inlet.dist <10] <- "WRR.PI"
dist.env$site.comp[dist.env$inlet.dist>10] <- "CR.PI"
dist.env$site.comp[dist.env$row > 0 & dist.env$row < 17 & dist.env$col > 0 & dist.env$col < 17] <- "CR.CR"
dist.env$site.comp[dist.env$row > 16 & dist.env$row < 33 & dist.env$col > 16 & dist.env$col < 33] <- "PI.PI"
dist.env$site.comp[dist.env$row > 32 & dist.env$row < 47 & dist.env$col > 32 & dist.env$col < 47] <- "WRR.WRR"



##__________________________WELLS____________________##

##refractometer readings
refract.wells <- wells.env[ ,c(12)]   
refract.dist.wells <- vegdist(refract.wells, method = "euclidean")#make euclidean distance matrix of refractometer differences          
refract.dist.wells <- dist.df(refract.dist.wells)
dist.env.wells <- bind_cols(refract.dist.wells, dist.list.wells)
dist.env.wells <- dist.env.wells[ ,-c(4,5)]
names(dist.env.wells)[3] <- "ssal.diff"
names(dist.env.wells)[4] <- "jacc.diss"

##inlet distance
inlet.wells <- wells.env[ ,c(11)]
inlet.dist.wells <- vegdist(inlet.wells, method = "euclidean")
inlet.dist.wells <- dist.df(inlet.dist.wells)
dist.env.wells <- bind_cols(inlet.dist.wells, dist.env.wells)
names(dist.env.wells)[3] <- "inlet.dist"
dist.env.wells <- dist.env.wells[ ,-c(4,5)]

##seasonal mean salinity
salseas.wells <- wells.env[ ,c(32)]
salseas.dist.wells <- vegdist(salseas.wells, method = "euclidean")#make euclidean distance matrix of refractometer differences          
salseas.dist.wells <- dist.df(salseas.dist.wells)
dist.env.wells <- bind_cols(salseas.dist.wells, dist.env.wells)
dist.env.wells <- dist.env.wells[ ,-c(4,5)]
names(dist.env.wells)[3] <- "seassal.diff"

##create factors from distances
dist.env.wells$site.comp <- NA 
dist.env.wells$site.comp[dist.env.wells$inlet.dist==0] <- NA
dist.env.wells$site.comp[dist.env.wells$inlet.dist>0 & dist.env.wells$inlet.dist<2] <- "CR.WRR"
dist.env.wells$site.comp[dist.env.wells$inlet.dist>2 & dist.env.wells$inlet.dist <10] <- "WRR.PI"
dist.env.wells$site.comp[dist.env.wells$inlet.dist>10] <- "CR.PI"
dist.env.wells$site.comp[dist.env.wells$row > 0 & dist.env.wells$row < 16 & dist.env.wells$col > 0 & dist.env.wells$col < 16] <- "CR.CR"
dist.env.wells$site.comp[dist.env.wells$row > 15 & dist.env.wells$row < 30 & dist.env.wells$col > 15 & dist.env.wells$col < 30] <- "PI.PI"
dist.env.wells$site.comp[dist.env.wells$row > 29 & dist.env.wells$row < 36 & dist.env.wells$col > 29 & dist.env.wells$col < 36] <- "WRR.WRR"

##aggregate distance-difference df for mean and sd of differences 
dist.env.wells$study <- "W"
dist.env.wells$site.comp <- as.character(dist.env.wells$site.comp)
dist.env$study <- "C"
dist.wells.mean <- dist.env.wells %>% group_by(site.comp) %>% summarise_each_(vars(mean)) #aggregate wells distances
dist.wells.mean <- dist.env.wells %>% 
                  dplyr::group_by(study,site.comp) %>% 
                  dplyr::summarise(
                            N=n(),
                            seassal.diff = mean(seassal.diff),
                            inlet.dist = mean(inlet.dist),
                            ssal.diff = mean(ssal.diff),
                            jacc.diss = mean(jacc.diss)) #aggregate wells distances

dist.wells.sd <- dist.env.wells %>% 
  dplyr::group_by(study,site.comp) %>% 
  dplyr::summarise(N=n(),
    seassal.diff = sd(seassal.diff),
                   inlet.dist = sd(inlet.dist),
                   ssal.diff = sd(ssal.diff),
                   jacc.diss = sd(jacc.diss)) #aggregate wells distances

dist.crfl.mean <- dist.env %>% 
  dplyr::group_by(study,site.comp) %>% 
  dplyr::summarise(N=n(),
    seassal.diff = mean(seassal.diff),
                   inlet.dist = mean(inlet.dist),
                   jacc.diss = mean(jacc.diss)) #aggregate wells distances

dist.crfl.sd <- dist.env %>% 
  dplyr::group_by(study,site.comp) %>% 
  dplyr::summarise(N=n(),
    seassal.diff = sd(seassal.diff),
                   inlet.dist = sd(inlet.dist),
                   jacc.diss = sd(jacc.diss)) #aggregate wells distances

dist.crfl.sd$study <- "C" 

##join mean and sd for crfl
dist.crfl.agg<-full_join(dist.crfl.mean, dist.crfl.sd, by = "site.comp")
dist.wells.agg<-full_join(dist.wells.mean, dist.wells.sd, by = "site.comp")
dist.agg.comb <- bind_rows(dist.crfl.agg, dist.wells.agg)

##Subset dist.agg.comb by between site comparisons
dist.agg.comb.1 <- dist.agg.comb %>%
  filter(inlet.dist.x>0) ##between site
dist.agg.comb.1$site.ord <- "na"
dist.agg.comb.1$site.ord[dist.agg.comb.1$site.comp=="CR.PI"] <- "PI-CR"
dist.agg.comb.1$site.ord[dist.agg.comb.1$site.comp=="CR.WRR"] <- "WR/WRR-CR"
dist.agg.comb.1$site.ord[dist.agg.comb.1$site.comp=="WRR.PI"] <- "PI-WR/WRR"
#order factors for graph
dist.agg.comb.1$site.ord <- as.factor(dist.agg.comb.1$site.ord)
dist.agg.comb.1$site.ord <- factor(dist.agg.comb.1$site.ord, levels=c('PI-WR/WRR', 'PI-CR', 'WR/WRR-CR'))
#order factors for graph study
dist.agg.comb.1$study.x <- as.factor(dist.agg.comb.1$study.x)
dist.agg.comb.1$study.x <- factor(dist.agg.comb.1$study.x, levels=c('W','C'))

##################################################This plot to be used in pubs###############################
##Plot all jaccard diss. between sites
jacc.pairs.final<-
ggplot(dist.agg.comb.1, aes(as.factor(site.ord), jacc.diss.x, fill = study.x)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) + 
  geom_bar(stat = "identity", position = position_dodge(0.9), color="black") + 
  #geom_errorbar(aes(ymax=jacc.diss.x+jacc.diss.y, ymin=jacc.diss.x-jacc.diss.y), 
                #position = position_dodge(0.9), width = .5, colour = 'black') +
  geom_errorbar(aes(ymax=jacc.diss.x+jacc.diss.y, ymin=jacc.diss.x), 
                position = position_dodge(0.9), width = .2, colour = 'black') +
  xlab("Site Comparison")+ylab("Jaccard distance") + 
  theme_classic() +
  scale_fill_manual(name = "Study", label=c("1955-1956 (Wells 1961)", "2013-2015 (This study)"),
                  values = c("black", "white"),
                  na.value = "red") +
  theme(axis.title.y=element_text(size=15), 
        axis.text.y = element_text(size=15,color="#000000"), 
        axis.title.x=element_text(size=15), 
        axis.text.x = element_text(size=15, color="#000000"),
        #legend.title = element_text(color = "#000000", size = 14),
        legend.title = element_blank(),
        legend.text = element_text(color = "#000000", size = 15),
        legend.key.height = unit(1.5,'lines'),
        legend.position = "bottom")
jacc.pairs.final
#ggsave(plot=jacc.pairs.final, file='Figures/Wells + Filtered CRFL/FINAL FIGS/jacc.pairs.final_20210522.png', width=20, height=15, units="cm", dpi=300)



##Subset dist.agg.comb by between site comparisons
dist.agg.comb.2 <- dist.agg.comb %>%
  filter(inlet.dist.x==0) ##within site
dist.agg.comb.2$site.ord <- "na"
dist.agg.comb.2$site.ord[dist.agg.comb.2$site.comp=="PI.PI"] <- "3-PI-PI"
dist.agg.comb.2$site.ord[dist.agg.comb.2$site.comp=="WRR.WRR"] <- "2-WRR-WRR"
dist.agg.comb.2$site.ord[dist.agg.comb.2$site.comp=="CR.CR"] <- "1-CR-CR"
dist.agg.comb.2<- dist.agg.comb.2 %>% arrange(site.ord)


##combine all difference/distance data
dist.env.cw <- bind_rows(dist.env.wells, dist.env)
dist.env.cw$comp.study <- NA
dist.env.cw$comp.study<- paste0(dist.env.cw$site.comp,".",dist.env.cw$study)
dist.env.cw$study<-as.factor(dist.env.cw$study)
dist.env.cw$study <- factor(dist.env.cw$study, levels=c('W','C'))

####_______________________________________________________________________________________________________________________________##

#####################################################PAIRWISE ANALYSIS####################################################################################_###############################
dist.env.cw.agg <- bind_rows(dist.wells.agg, dist.crfl.agg)
dist.env.cw.agg$jacc.diss.var <- NA
dist.env.cw.agg$jacc.diss.var <- dist.env.cw.agg$jacc.diss.y * dist.env.cw.agg$jacc.diss.y

#Kruskal Wallis test of jaccard differences and salinity differences: 

dist.env.cw$study <- as.factor(dist.env.cw$study)
dist.env.cw$comp.study <- as.factor(dist.env.cw$comp.study)
##Run Shapiro-Wilks test to test for normality: 
#shapiro.test(dist.env.cw$jacc.diss)
histogram(dist.env.cw$jacc.diss[dist.env.cw$comp.study == "CR.PI.W"], dist.env.cw)
histogram(dist.env.cw$jacc.diss[dist.env.cw$comp.study == "CR.PI.C"], dist.env.cw)
histogram(dist.env.cw$jacc.diss[dist.env.cw$comp.study == "CR.WRR.W"], dist.env.cw)
histogram(dist.env.cw$jacc.diss[dist.env.cw$comp.study == "CR.WRR.C"], dist.env.cw)
histogram(dist.env.cw$jacc.diss[dist.env.cw$comp.study == "WRR.PI.W"], dist.env.cw)
histogram(dist.env.cw$jacc.diss[dist.env.cw$comp.study == "WRR.PI.C"], dist.env.cw)

#Pairwise tests Nonparametric

##CR_PI Mann Whitney U
mw.crpi.jacc <- wilcox.test(dist.env.cw$jacc.diss[dist.env.cw$site.comp =="CR.PI"] ~ dist.env.cw$study[dist.env.cw$site.comp =="CR.PI"], data=dist.env.cw)
print(mw.crpi.jacc)
length(dist.env.cw$jacc.diss[dist.env.cw$site.comp =="CR.PI"])

##CR_WRR Mann Whitney U
mw.crwrr.jacc <- wilcox.test(dist.env.cw$jacc.diss[dist.env.cw$site.comp =="CR.WRR"] ~ dist.env.cw$study[dist.env.cw$site.comp =="CR.WRR"], data=dist.env.cw)
print(mw.crwrr.jacc)
length(dist.env.cw$jacc.diss[dist.env.cw$site.comp =="CR.WRR"])

##WRR-PI Mann Whitney U
mw.wrrpi.jacc <- wilcox.test(dist.env.cw$jacc.diss[dist.env.cw$site.comp =="WRR.PI"] ~ dist.env.cw$study[dist.env.cw$site.comp =="WRR.PI"], data=dist.env.cw)
print(mw.wrrpi.jacc)
length(dist.env.cw$jacc.diss[dist.env.cw$site.comp =="WRR.PI"])

#pairwise wilcoxon tests: 
wilcox <- wilcox.test(dist.env.cw$jacc.diss[dist.env.cw$comp.study])
print(wilcox)


####______________Start nmds plotting process_________________###### Based on nmds ggplot example that uses dune. #####

###Create an MDS plot of salbin.study
meta.nmds.fam <- metaMDS(crfl.wells.fam, distance = "jaccard",  trymax = 1000)
###stress value for plot
str(meta.nmds.fam)
##
####create a shepard diagram to examine relationship between dissimilarity and ord distances **needs to be done immediately after above code
stressplot(meta.nmds.fam)

##take out crfl.wells.env columns that do are not relevant to analysis
crfl.wells.env <- crfl.wells.env[ ,-c(7,8,13:22,28,29,30)]
names(crfl.wells.env)[17] <- "Salinity"
names(crfl.wells.env)[7] <- "Site Salinity"
names(crfl.wells.env)[8] <- "salinity"
names(crfl.wells.env)[9] <- "Inlet Distance"
names(crfl.wells.env)[10] <- "Sample Salinity"
names(crfl.wells.env)[11] <- "Julian Day"

crfl.wells.env$Site <- NA
crfl.wells.env$Site[crfl.wells.env$site.x == "npcr"] <- "CR"
crfl.wells.env$Site[crfl.wells.env$site.x == "npwrr" & crfl.wells.env$study == "w"] <- "WR"
crfl.wells.env$Site[crfl.wells.env$site.x == "nppi"] <- "PI"
crfl.wells.env$Site[crfl.wells.env$site.x == "npwrr" & crfl.wells.env$study == "c"] <- "WRR"

crfl.wells.env$Study <- NA
crfl.wells.env$Study[crfl.wells.env$study == "w"] <- "1955-1956 (Wells 1961)"
crfl.wells.env$Study[crfl.wells.env$study == "c"] <- "2013-2015 (This study)"

crfl.wells.env$site.study <- NA
crfl.wells.env$site.study <- paste0(crfl.wells.env$Site," ",crfl.wells.env$Study)
crfl.wells.env$site.study <- as.factor(crfl.wells.env$site.study)
crfl.wells.env$site.study <- factor(crfl.wells.env$site.study, c("PI 1955-1956 (Wells 1961)", "WR 1955-1956 (Wells 1961)", "CR 1955-1956 (Wells 1961)","PI 2013-2015 (This study)", "WRR 2013-2015 (This study)", "CR 2013-2015 (This study)"))
names(crfl.wells.env)[6] <- "Year"


#envfit
cwf.envfit <- envfit(meta.nmds.fam, env = crfl.wells.env[ ,c(17,6,9,11)], perm = 999, na.rm = TRUE) #standard envfit
cwf.envfit

#data for plotting 
##NMDS points
cwf.nmds.data<-crfl.wells.env #there are other ways of doing this. But this is the way I do it for ease of plotting
cwf.nmds.data$NMDS1<-meta.nmds.fam$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
cwf.nmds.data$NMDS2<-meta.nmds.fam$points[ ,2] 

##family data
fam.occur<-colSums(crfl.wells.fam) #total abundances for each species
fams <- data.frame(scores(meta.nmds.fam, display = "species")) #dataframe of species scoes for plotting
fams$families <- row.names(fams) # making a column with species names
fams$colsums <- fam.occur #adding the colSums from above
fams<-fams[!is.na(fams$NMDS1) & !is.na(fams$NMDS2),] #removes NAs
fams.colmedian <- median(fams$colsums) #create an object that is the median of the abundance of the measured species
fams.colmean <- mean(fams$colsums) #creates a mean instead if you wish to use
fams2 <- subset(fams,fams$colsums > fams.colmean) #select the most abundant species. Could discard fewer by going something like - spps$colsums>(spps.colmedian/2) instead
fams2$families <- factor(fams2$families) #otherwise factor doesn't drop unused levels and it will throw an error


# data for the envfit arrows
env.scores <- as.data.frame(scores(cwf.envfit, display = "vectors")) #extracts relevant scores from envifit
env.scores <- cbind(env.scores, env.variables = rownames(env.scores)) #and then gives them their names
env.scores <- env.scores[-c(3,4), ]

###########NMDS PLOTTING##########################################################################

# function for ellipsess - just run this, is used later
#taken from the excellent stackoverflow Q+A: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the salbin.study

ellipse.site.study <- data.frame() #sets up a data frame before running the function.
for(g in levels(cwf.nmds.data$site.study)) {
  ellipse.site.study <- rbind(ellipse.site.study, 
                              cbind(as.data.frame(with(cwf.nmds.data [crfl.wells.env$site.study==g,], 
                                                       veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),
                                                                                                        length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                    ,site.study=g))
}

# data for labelling the ellipse
NMDS.mean=aggregate(cwf.nmds.data[ ,c("NMDS1", "NMDS2")], 
                    list(group = cwf.nmds.data$site.study), mean)
library(stringr)
##study.salbin group nmds#### 
## nmds plotting - group by study.salbin
mult <- 2 #multiplier for the arrows and text for envfit below. You can change this and then rerun the plot command.
#tiff("nmds.cwfams.wellsfilt.tiff", width=7.5, height=6.5, units="in",pointsize=8, res=600) 
nMDS.jacc.final<-
ggplot(data = cwf.nmds.data, aes(x = NMDS1, y = NMDS2))+ #sets up the plot. brackets around the entire thing to make it draw automatically
  geom_path(data = ellipse.site.study, aes(x = NMDS1, y = NMDS2, group = site.study, alpha=site.study, linetype=site.study))+ #this is the ellipse, seperate ones by Site. If you didn't change the "alpha" (the shade) then you need to keep the "group 
  scale_alpha_manual(guide = FALSE,values=c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9))+ #sets the shade for the ellipse
  #scale_linetype_manual(values=c("solid", "solid", "solid", "dashed", "dashed", "dashed")) +
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "solid", "solid")) +
  #scale_color_manual(values=c("orange2", "orangered2", "red4","steelblue1", "blue2", "navy")) +
  geom_point(aes(shape = site.study), size = 3) + #puts the site points in from the ordination, shape determined by site, size refers to size of point
  #geom_text(data=fams2, aes(x=fams2$NMDS1, y=fams2$NMDS2, label=families), size = 3.3, hjust=1.1)+ #labelling the species. hjust used to shift them slightly from their points
  #annotate("text",x = NMDS.mean$NMDS1,y = NMDS.mean$NMDS2,label=NMDS.mean$group) + #labels for the centroids - I haven't used this since we have a legend. but you could also dithc the legend, but plot will get v messy
  geom_segment(data = env.scores,aes(x = 0, xend = mult*NMDS1, y = 0, yend = mult*NMDS2),arrow = arrow(length = unit(0.25, "cm")), colour = "grey35") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores, #labels the environmental variable arrows * "mult" as for the arrows
            aes(x = NMDS1, y = NMDS2, label=env.variables), size = 7, hjust = .9, vjust = -5)+
  #geom_point(data=fams2, alpha = .6, shape = 4)+ #these are the species points, made lighter and a specific shape
  scale_shape_manual(values = c(17,16,15,2,1,0))+ #sets the shape of the plot points instead of using whatever ggplot2 automatically provides
  #scale_color_manual(values=c("orange2", "orangered2", "red4","steelblue1", "blue2", "navy")) +
  #coord_cartesian(xlim = c(-1.5,1.0)) +  ## NB this changes the visible area of the plot only (this is a good thing, apparently). Can also specify ylim. Here in case you want to set xaxis manually.
  xlab("nMDS1")+ ylab("nMDS2") + 
  theme(
    #text=element_text(family="Times New Roman"),
    panel.grid.major = element_blank(), # removes major gridlines
    panel.grid.minor = element_blank(), # removes minor gridlines
    panel.background = element_blank(), # removes grey background
    axis.line = element_line(size = .5,colour = "black"), # keeps axis lines
    strip.background = element_blank(),
    legend.key = element_blank(),  # ledgend specifics
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    #legend.key.size = unit(2,"line"),
    #legend.position = "bottom",
    #legend.spacing.y=unit(3, "cm"),
    axis.title.x = element_text(size=20), # sets size of x axis label
    axis.title.y = element_text(size=20),
    axis.text.y = element_text(size=20, color="black"),
    axis.text.x = element_text(size=20, color="black"),
    legend.position="bottom") +
    guides(linetype = guide_legend(ncol = 2,bycol=TRUE))
nMDS.jacc.final
#ggsave(plot=nMDS.jacc.final, file='Figures/Wells + Filtered CRFL/FINAL FIGS/nMDS.jacc.final_20201025.png', width=28, height=20, units="cm", dpi=300)
#ggsave(plot=nMDS.jacc.final, file='Figures/Wells + Filtered CRFL/FINAL FIGS/nMDS.jacc.final_20210522.png', width=28, height=20, units="cm", dpi=300)

## separate data based on study for individual paranova tests 
wells.fam <- crfl.wells.fam[c(1:35), ]
wells.env <- crfl.wells.env[c(1:35), ]
crfl.fam <- crfl.wells.fam[c(36:81), ]
crfl.env <- crfl.wells.env[c(36:81), ]


##adonis study, site, crossed
adonis(crfl.wells.fam ~ study + use.y + study*use.y , data = crfl.wells.env, method = "jaccard", permutation = 9999 )
##adonis site.study
adonis(crfl.wells.fam ~ site.study, data = crfl.wells.env, method = "jaccard", permutation = 9999 )
##adonis study
adonis(crfl.wells.fam ~ study, data = crfl.wells.env, method = "jaccard", permutation = 9999 )


### Homog. dispersions study.salbn

permdisp.fam<- betadisper(jaccard.fam, group = crfl.wells.env$site.study, type = "centroid")
print(permdisp.fam)
plot(permdisp.fam)
boxplot(permdisp.fam)
permutest(permdisp.fam, pairwise = TRUE)


## Homog. dispersion salbin for crfl and wells matrices individually
jaccard.crfl <- vegdist(crfl.fam, method="jaccard")
jaccard.wells <- vegdist(wells.fam, method="jaccard")

permdisp.crfl<- betadisper(jaccard.crfl, group = crfl.env$site.study, type = "centroid")
permutest(permdisp.crfl, pairwise = TRUE)
boxplot(permdisp.crfl)

permdisp.wells<- betadisper(jaccard.wells, group = wells.env$site.study, type = "centroid")
permutest(permdisp.wells, pairwise = TRUE)
boxplot(permdisp.wells)

### Simper Analysis study
simper.fam <- simper(crfl.wells.fam, group=crfl.wells.env$study)
summary(simper.fam)
print(simper.fam)
lapply(simper.fam, FUN=function(x){x$overall})
### Simper Analysis study.site
simper.study.site <- simper(crfl.wells.fam, group=crfl.wells.env$site.study)
summary(simper.study.site)

##overall dissimilarity 
lapply(simper.study.site, FUN=function(x){x$overall})

###Permanova based on study 
adonis(crfl.wells.fam ~ crfl.wells.env$study, data = crfl.wells.env, method = "jaccard", permutation = 9999 )


#################################################################################################################################
##beta diversity 
betadiv.jacc <- betadiver(crfl.wells.fam, "j")
betadiv.jacc<- data.frame(scores(betadiv.jacc, display = "species"))

### Homog. dispersions study
permdisp.fam.study<- betadisper(jaccard.fam, group = crfl.wells.env$study, type = "centroid")
print(permdisp.fam.study)
plot(permdisp.fam.study)
boxplot(permdisp.fam.study)
permutest(permdisp.fam.study, pairwise = TRUE)
scores(permdisp.fam.study, display = c("sites", "centroids"),
       choices = c(1,2))
jacc.disp <- data.frame(scores(permdisp.fam.study, display = c("sites"),
                               choices = c(1,2)))
jacc.disp2 <- data.frame(scores(permdisp.fam.study, display = c("centroid"),
                               choices = c(1,2)))
jacc.disp$centroidx <- NA
jacc.disp$centroidy <- NA
jacc.disp$centroidx

#################################################################################################################################
##Post-hoc fxn 
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 
pairwise.adonis(crfl.wells.fam, crfl.wells.env$site.study)
pairwise.adonis(crfl.fam, crfl.env$site.x)
pairwise.adonis(wells.fam, wells.env$site.x)



