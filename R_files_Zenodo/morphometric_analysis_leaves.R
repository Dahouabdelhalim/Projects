#load all the libraries we will need
library(Momocs)   #for shape analysis
library(ggplot2)  #for graphing
library(dplyr)    #for dataframe manipulation
library(tidyr)    #for dataframe manipulation

#********************************************************************
#GET LEAF SIZE AND SHAPE DATA FOR HEATHER'S LEAVES (2017) -----------------------
#********************************************************************

#set the location to the folder where your leaf scans and Groups file live on your computer
setwd("Leaf_scans-Processed/Final_Images_JPEG-rotated/")

#inputting leaf traces----
#get the library that we will use for morphometric analysis
#( NOTE - the first time you use a new library run: e.g. > install.packages("Momocs") )

#read all of the leaf scans in (this reads in all .jpg files in folder we specified above)
jpg.list_leaves <- list.files(pattern = c( ".jpg")) 
#look at the first few lines of the list of jpg files - do these look like the scans you want?
head(jpg.list_leaves)
#get a trace of each leaf
returns_leaves <- import_jpg(jpg.list_leaves, auto.notcentered = T) 
leaves <- Out(returns_leaves)
#getting leaf measurements ----
names(leaves)  ### look at names of all leaves

par(mfrow=c(1,2))

#get object that is xy coords for each leaf ------------
lf <- leaves$coo

#get length and width for all leaves (can input whole Out object into this function)
lw <- coo_lw(leaves) %>% t() %>% as.data.frame() %>% mutate(sampid=row.names(.)) %>% 
  mutate(length_cm=((V1/300)*2.54)) %>% mutate(width_cm=((V2/300)*2.54)) %>% dplyr::select(sampid,length_cm,width_cm) %>% 
  mutate(ratio_LoverW_cm=length_cm/width_cm)


#get solidity measure for each leaf (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.solid = data.frame(sampid = names(lf), solidity_score = NA)
head(lf.solid)
dim(lf.solid)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.solid$solidity_score[i] = coo_solidity(lf[lf.solid$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.solid)


#get area measure (can input whole Out object into this function)
area <- coo_area(leaves) %>% as.data.frame()
colnames(area)[colnames(area) == "."] <- "Area"
area <- area %>% mutate(sampid=row.names(.)) %>% mutate(area_cm2=((Area/(300*300))*(2.54^2))) %>% 
  dplyr::select(sampid,area_cm2)


#get perimeter for each leaf (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.perim = data.frame(sampid = names(lf), perim = NA)
head(lf.perim)
dim(lf.perim)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.perim$perim[i] = coo_perim(lf[lf.perim$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.perim)
#convert to cm
lf.perim <- lf.perim %>% mutate(perim_cm=(perim/300)*2.54) %>% dplyr::select(sampid, perim_cm)


#get calliper measurement for each leaf - longest distance between points (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.calip = data.frame(sampid = names(lf), longest_dist = NA)
head(lf.calip)
dim(lf.calip)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.calip$longest_dist[i] = coo_calliper(lf[lf.calip$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.calip)
#convert to cm
lf.calip <- lf.calip %>% mutate(longest_dist_cm=(longest_dist/300)*2.54) %>% dplyr::select(sampid, longest_dist_cm)

#merge all measurements together into one dataframe--------
mdf <- merge(lw,lf.calip, by="sampid")
mdf <- merge(mdf,lf.perim, by="sampid")
mdf <- merge(mdf,area, by="sampid")
mdf <- merge(mdf,lf.solid, by="sampid")


#getting leaves ready for morphometric analysis----
#read in the Groups file - this has the detailed information about each leaf
groups_leaves <- read.csv("Groups_file.csv", header = T)
#IMPORTANT: the numbers at the end of the next line represent how many columns are in your Groups file 
#e.g. 1:4 means there are 4 columns - change this to match number of columns in your Groups file
groups_leaves <-data.frame(groups_leaves[,1:7])
#look at how the variables are coded - grouping varibales should be "Factor"
str(groups_leaves)
#change varibales to factors
groups_leaves[,'leaf'] <- as.factor(groups_leaves[,'leaf'])
groups_leaves[,'envelope'] <- as.factor(groups_leaves[,'envelope'])
#check if all variables that you want to be Factors are now
str(groups_leaves)
#attach the leaf info to the leaf traces
leaves$fac <- groups_leaves

#look at a panel of all of the leaf outlines - do any look weird?
panel(leaves, col=0, border=1)
#if you need to see the names to find a specific leaf, display the names
#change the size of the print with "cex.names"
panel(leaves, col=0, border=1, names = T, cex.names = 0.1)
      
#interpolate n coordinates along the perimeter of the trace
leaves$coo <- lapply(leaves$coo, coo_interpolate, n=500) #coo object, stack(leaves) gives all orig. traces, not centered or scaled 
#center all of the traces on the origin/on top of each other 
leaves <- coo_center(leaves) #coo object, stack(leaves) gives all orig. traces, centered not scaled 
leaves.plotting <- leaves
#scale the traces
leaves$coo <- lapply(leaves$coo, FUN=function(x)(coo_scale(x, max(x[,2]/1000)))) 
leaves <- coo_scale(leaves, 1000)

#put landmarks on each leaf, so they all end up in the correct final orientation in the analysis
ldks2 <- vector("list", length=length(leaves$coo))
names(ldks2) <- names(leaves$coo)
centroids <- coo_centpos(leaves) 
for (i in 1:length(leaves)) {
      ldk <- numeric(4)
      leaves$coo[[i]] -> x
      which.min(x[,1]) -> min.x
      which.max(x[,1]) -> max.x
      which.min(x[,2]) -> min.y
      which.max(x[,2]) -> max.y
      p <- c(centroids[i,1], min(x[,2]))
      l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
      ldk[1] <- which.min(l)
      p <- c(centroids[i,1], max(x[,2]))
      l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
      ldk[2] <- which.min(l)
      p <- c(min(x[,1]), centroids[i,2])
      l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
      ldk[3] <- which.min(l)
      p <- c(max(x[,1]), centroids[i,2])
      l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
      ldk[4] <- which.min(l)
      #plot(x, asp=1, axes=F, xlab="", ylab="")
      #points(x[ldk[1],1], x[ldk[1],2], pch=20, col="red", cex=4)
      #points(x[ldk[2],1], x[ldk[2],2], pch=20, col="red", cex=4)
      #points(x[ldk[3],1], x[ldk[3],2], pch=20, col="red", cex=4)
      #points(x[ldk[4],1], x[ldk[4],2], pch=20, col="red", cex=4)
      ldk -> ldks2[[i]]
    }
leaves$ldk <- ldks2
leaves <- coo_slide(leaves, ldk=1)
leaves <- fgProcrustes(leaves)

#visualize where the landmarks are
plot(x, asp=1, axes=F, xlab="", ylab="")
points(x[ldk[1],1], x[ldk[1],2], pch=20, col="red", cex=2)
points(x[ldk[2],1], x[ldk[2],2], pch=20, col="red", cex=2)
points(x[ldk[3],1], x[ldk[3],2], pch=20, col="red", cex=2)
points(x[ldk[4],1], x[ldk[4],2], pch=20, col="red", cex=2)

#Do the elliptic Fourier analyses and PCA of the results
efou_leaves <- efourier(leaves, norm=F, smooth.it = 0, verbose = T, nb.h = 20)
#Summarize the data by running a principal components analysis (PCA) on it
pca_efou_leaves <- PCA(efou_leaves, center = T) 
#look at the output from the PCA - how much variance do the first few principle components explain?
#e.g. 0.68 under PC1, Cumulative Proportion means the first PC explains/captures 68% of the variation in the data
summary(pca_efou_leaves)

#get PC scores into dataframe-----
pcscores <- pca_efou_leaves$x
pcscores <- as.data.frame(pcscores)
pcscores$sampid <- rownames(pcscores) 
pcscores.sub <- pcscores[,c(1:10)]
pcscores.sub <- as.data.frame(pcscores.sub)
pcscores.sub$sampid <- rownames(pcscores.sub)

dim(mdf)
dim(pcscores.sub)
masterdf <- merge(mdf,pcscores.sub,by="sampid")
dim(masterdf)

groups_leaves$sampid <- gsub(".jpg","",groups_leaves$image, fixed = T)
masterdf <- merge(groups_leaves,masterdf,by = "sampid")
colors <- read.csv("../../leaf_trait_data_files/RGB_histogram_values_for_subset.csv")
colors$sampid <- gsub(".tif","",colors$image_name, fixed = T)
masterdf <- merge(masterdf,colors,by = "sampid", all = T)
dim(masterdf)

#write out data ------------
#write data
write.csv(masterdf, "../../leaf_trait_data_files/full_data_10pcs_20harmonics_justheathersleaves.csv")
#write pca R object
save(pca_efou_leaves, 
     file = "../../output_results/pca_efou_leaves-10pcs_20harmonics_justheathersleaves.RData")
#write out efourier coe object
save(efou_leaves, 
     file = "../../output_results/efou_leaves-10pcs_20harmonics_justheathersleaves.RData")
#write out original traces (centered but not scaled)
save(leaves.plotting, 
     file = "../../output_results/coo_leaves-centeredtraces-justheathersleaves.RData")



#write out centered (but NOT scaled or run through efourier analysis) outlines for all of Heather's leaves -------------
#load shape outlines
load("../../output_results/coo_leaves-centeredtraces-justheathersleaves.RData")
#turn list of coords into long dataframe
df <- leaves.plotting$coo
df.coo <- do.call(rbind.data.frame, df) %>% mutate(id.temp = rownames(.)) %>% 
  mutate(x = V1) %>% mutate(y = V2) %>% dplyr::select(-V1, -V2) %>% 
  separate(.,col = "id.temp", into = c("id.temp","leaf","point"), sep = "\\\\.") %>% 
  separate(.,col = "id.temp", into = c("species","site","area","plant"), sep = "_")

rm(df)

write.csv(df.coo,"../../output_results/Coo_leaf_outlines_centerednotscaled_justheathersleaves.csv")


#write out avg. shape outlines at different factor levels heather's leaves only ------------
#load shape outlines
load("../../output_results/efou_leaves-10pcs_20harmonics_justheathersleaves.RData")
#get mean shape coordinates for each population
ms.s <- MSHAPES(efou_leaves, fac = 'species', FUN = mean, nb.pts = 10000)
ms.s <- ms.s$shp
#turn list of coords into long dataframe
df.s <- do.call(rbind.data.frame, ms.s) %>% mutate(id.temp = rownames(.)) %>% mutate(level = "species") %>% 
  separate(.,col = "id.temp", into = c("id","trash"), sep = "\\\\.") %>% dplyr::select(-trash)
#repeat
#for species_leaf
ms.sl <- MSHAPES(efou_leaves, fac = 'species_leaf', FUN = mean, nb.pts = 10000)
ms.sl <- ms.sl$shp
df.sl <- do.call(rbind.data.frame, ms.sl) %>% mutate(id.temp = rownames(.)) %>% mutate(level = "species_leaf") %>% 
  separate(.,col = "id.temp", into = c("id","trash"), sep = "\\\\.") %>% dplyr::select(-trash)

df.f <- rbind(df.s,df.sl) #all outline data

write.csv(df.f,"../../output_results/MSHAPES_efourier_avg_leaf_outlines_justheathersleaves.csv")


#see the average leaf - shape and size ("original" outlines) - *FIG. 3* ------------
df.outline <- read.csv("../../output_results/Coo_leaf_outlines_centerednotscaled_justheathersleaves.csv", header = T) %>% dplyr::select(-X)
df.outline$plant <- as.factor(df.outline$plant)
df.outline$leaf <- as.factor(df.outline$leaf)
df.outline$point <- as.factor(df.outline$point)
str(df.outline)
df.outline <- df.outline %>% group_by(point,species,leaf) %>% summarise(x = mean(x), y = mean(y)) %>% as.data.frame(.)
df.outline$species <- as.character(df.outline$species) 
df.outline$species[df.outline$species == "Pallida"] = "apallida"
df.outline$species <- as.factor(df.outline$species)

#### adding scale bar
#get height of pallida leaf 1 in the graph ( = 10.913 cm)
ma <- df.outline %>% filter(species == "apallida") %>% filter(leaf == "1") 
maxht <- max(ma$y, na.rm = TRUE) - min(ma$y, na.rm = TRUE) 
bar <- (maxht*2)/10.913 #get graph length of a 2cm scale bar (cross multiply/do a proportion)
#find center value of pallida leaf 3 in the graph (so we can center scale bar)
ma <- df.outline %>% filter(species == "apallida") %>% filter(leaf == "3")
md <- median(ma$x, na.rm = TRUE)
#calc x and xend values so scale bar is 5cm and centered on (pallida) leaf 3
xxend <- (bar/2) + md
xx <- md - (bar/2)

scalebar<-data.frame(x=xx+85,y=-490,xend=xxend+85,yend=-490,species="Capensis",leaf="3") #build scale bar
ann_text <- data.frame(x = 135, y = -525, species="Capensis", leaf="3") #add label


#per species and leaf
df.outline %>%
  ggplot(aes(x = x, y = y)) +
  geom_polygon(aes(fill = species, group = species, colour = species), alpha = 0.75, size = 0.7) +
  geom_segment(data=scalebar,aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, color = "black") +
  geom_text(data = ann_text, label = "2 cm", size = 3.25) +
  scale_fill_manual(values = c("#2D8D53",NA),
                    name = "Species",
                    labels= c("I. capensis","I. pallida"),
                    breaks = c("Capensis","apallida")) +
  scale_colour_manual(values = c(NA,"black"),
                      name = "Species",
                      labels= c("I. capensis","I. pallida"),
                      breaks = c("Capensis","apallida")) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_blank(),
        axis.title = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(face = "italic", size = 9),
        strip.text = element_text(size = 9, hjust=0.52, vjust=-1.6),
        legend.box.spacing = unit(-0.5,"line"),
        legend.title = element_text(size = 10),
        panel.spacing = unit(0.02, "lines"),
        plot.margin=unit(c(0.25,0.25,0.25,0.25),"cm")) +
  coord_fixed(expand = F) +
  facet_wrap(~ leaf, strip.position = c("bottom"))

ggsave(paste("../../figures_and_tables/avgleaf_sizeandshape_preefourier_justheathersleaves","pdf", sep = ".") , 
       width = 13.8, height = 8, dpi = 600, units = c("cm"))


#see the average shapes (post efourier analysis - scaled to standard size) ------------
df <- read.csv("../../output_results/MSHAPES_efourier_avg_leaf_outlines_justheathersleaves.csv") %>% dplyr::select(-X)
str(df)
#per species and leaf
df %>% filter(level == "species_leaf") %>% separate(.,col = "id", into = c("species","leaf"), sep = "_") %>% 
  ggplot(aes(x = x, y = y)) +
  geom_polygon(aes(fill = species, group = species, colour = species), alpha = 0.5) +
  scale_fill_manual(values = c(NA,"forestgreen")) +
  scale_colour_manual(values = c("black",NA)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text = element_blank(),
                     axis.title = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(),
                     axis.ticks = element_blank(),
                     strip.background = element_blank()) +
  coord_fixed() +
  facet_wrap(~ leaf, strip.position = c("bottom"))

ggsave(paste("../../figures_and_tables/avgshape_postefourier_justheathersleaves","pdf", sep = ".") , 
       width = 6.7, height = 6, dpi = 600, units = c("cm"))





#Morphometric analysis VISUALIZATION -------------------

#load object needed for plotting previously saved to avoid rerunning all of the code
load("../../output_results/pca_efou_leaves-10pcs_20harmonics_justheathersleaves.RData")
load("../../output_results/efou_leaves-10pcs_20harmonics_justheathersleaves.RData")
df <- read.csv("../../leaf_trait_data_files/full_data_10pcs_20harmonics_justheathersleaves.csv")


#Load the library ggplot2 for graphing/viewing our results
library(ggplot2)
library(ggridges)

#see average shapes for X number of PCs
PCcontrib(pca_efou_leaves, nax=c(1:10), sd.r = c(-2, -1, 0, 1, 2), cex.labels=10) 

#see the average shape for each species
ms <- MSHAPES(leaves.plotting)

ms <- MSHAPES(efou_leaves, fac = 'species', FUN = mean)
ms <- ms$shp
coo_plot(ms$Capensis)
coo_draw(ms$Pallida, border = 'purple')

#just leaf3
ms3 <- efou_leaves %>% Momocs::filter(leaf == 3)
ms <- MSHAPES(ms3, fac = 'species', FUN = mean)
ms$Coe
class(ms$Coe)
ms <- ms$shp
coo_plot(ms$Capensis)
coo_draw(ms$Pallida, border = 'purple')

#just leaf2
ms2 <- efou_leaves %>% Momocs::filter(leaf == 2)
ms <- MSHAPES(ms2, fac = 'species', FUN = mean)
ms$Coe
class(ms$Coe)
ms <- ms$shp
coo_plot(ms$Capensis)
coo_draw(ms$Pallida, border = 'purple')

#just leaf1
ms1 <- efou_leaves %>% Momocs::filter(leaf == 1)
ms <- MSHAPES(ms1, fac = 'species', FUN = mean)
ms$Coe
class(ms$Coe)
ms <- ms$shp
coo_plot(ms$Capensis)
coo_draw(ms$Pallida, border = 'purple')



#********************************************************************
#GET LEAF SIZE AND SHAPE DATA FOR RACHEL'S LEAVES (2014) -----------------------
#********************************************************************
#collected about 80 I.pallida leaves randomly around Madison to test this idea
#added a random subset of ~80 more I.capensis leaves from Rachel's field sites

setwd("Leaf_scans-Processed/Final_Images_JPEG-rotated-RT/")

#load all the libraries we will need
library(Momocs)   #for shape analysis
library(ggplot2)  #for graphing
library(dplyr)    #for dataframe manipulation
library(tidyr)    #for dataframe manipulation

#inputting leaf traces----
#read all of the leaf scans in (this reads in all .jpg files in folder we specified above)
jpg.list_leaves <- list.files(pattern = c( ".jpg")) 
#look at the first few lines of the list of jpg files - do these look like the scans you want?
head(jpg.list_leaves)
#get a trace of each leaf
returns_leaves <- import_jpg(jpg.list_leaves, auto.notcentered = T) 
leaves <- Out(returns_leaves)
#getting leaf measurements ----
names(leaves)  ### look at names of all leaves

par(mfrow=c(1,2))

#get object that is xy coords for each leaf
lf <- leaves$coo

#get length and width for all leaves (can input whole Out object into this function)
lw <- coo_lw(leaves) %>% t() %>% as.data.frame() %>% mutate(sampid=row.names(.)) %>% 
  mutate(length_cm=((V1/300)*2.54)) %>% mutate(width_cm=((V2/300)*2.54)) %>% dplyr::select(sampid, length_cm, width_cm) %>% 
  mutate(ratio_LoverW_cm=length_cm/width_cm)


#get solidity measure for each leaf (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.solid = data.frame(sampid = names(lf), solidity_score = NA)
head(lf.solid)
dim(lf.solid)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.solid$solidity_score[i] = coo_solidity(lf[lf.solid$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.solid)


#get area measure (can input whole Out object into this function)
area <- coo_area(leaves) %>% as.data.frame()
colnames(area)[colnames(area) == "."] <- "Area"
area <- area %>% mutate(sampid=row.names(.)) %>% mutate(area_cm2=((Area/(300*300))*(2.54^2))) %>% 
  dplyr::select(sampid,area_cm2)


#get perimeter for each leaf (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.perim = data.frame(sampid = names(lf), perim = NA)
head(lf.perim)
dim(lf.perim)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.perim$perim[i] = coo_perim(lf[lf.perim$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.perim)
#convert to cm
lf.perim <- lf.perim %>% mutate(perim_cm=(perim/300)*2.54) %>% dplyr::select(sampid, perim_cm)


#get calliper measurement for each leaf - longest distance between points (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.calip = data.frame(sampid = names(lf), longest_dist = NA)
head(lf.calip)
dim(lf.calip)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.calip$longest_dist[i] = coo_calliper(lf[lf.calip$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.calip)
#convert to cm
lf.calip <- lf.calip %>% mutate(longest_dist_cm=(longest_dist/300)*2.54) %>% dplyr::select(sampid, longest_dist_cm)

#merge all measurements together into one dataframe--------
mdf <- merge(lw,lf.calip, by="sampid")
mdf <- merge(mdf,lf.perim, by="sampid")
mdf <- merge(mdf,area, by="sampid")
mdf <- merge(mdf,lf.solid, by="sampid")

#getting leaves ready for morphometric analysis----
#read in the Groups file - this has the detailed information about each leaf
groups_leaves <- read.csv("Groups_file.csv", header = T)
#IMPORTANT: the numbers at the end of the next line represent how many columns are in your Groups file 
#e.g. 1:4 means there are 4 columns - change this to match number of columns in your Groups file
groups_leaves <-data.frame(groups_leaves[,1:2])
#look at how the variables are coded - grouping varibales should be "Factor"
str(groups_leaves)
#attach the leaf info to the leaf traces
leaves$fac <- groups_leaves

#look at a panel of all of the leaf outlines - do any look weird?
panel(leaves, col=0, border=1)
#if you need to see the names to find a specific leaf, display the names
#change the size of the print with "cex.names"
panel(leaves, col=0, border=1, names = T, cex.names = 0.1)

#interpolate n coordinates along the perimeter of the trace
leaves$coo <- lapply(leaves$coo, coo_interpolate, n=500)
#center all of the traces on the origin/on top of each other 
leaves <- coo_center(leaves)
#scale the traces
leaves$coo <- lapply(leaves$coo, FUN=function(x)(coo_scale(x, max(x[,2]/1000)))) 
leaves <- coo_scale(leaves, 1000)

#put landmarks on each leaf, so they all end up in the correct final orientation in the analysis
ldks2 <- vector("list", length=length(leaves$coo))
names(ldks2) <- names(leaves$coo)
centroids <- coo_centpos(leaves) 
for (i in 1:length(leaves)) {
  ldk <- numeric(4)
  leaves$coo[[i]] -> x
  which.min(x[,1]) -> min.x
  which.max(x[,1]) -> max.x
  which.min(x[,2]) -> min.y
  which.max(x[,2]) -> max.y
  p <- c(centroids[i,1], min(x[,2]))
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[1] <- which.min(l)
  p <- c(centroids[i,1], max(x[,2]))
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[2] <- which.min(l)
  p <- c(min(x[,1]), centroids[i,2])
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[3] <- which.min(l)
  p <- c(max(x[,1]), centroids[i,2])
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[4] <- which.min(l)
  #plot(x, asp=1, axes=F, xlab="", ylab="")
  #points(x[ldk[1],1], x[ldk[1],2], pch=20, col="red", cex=4)
  #points(x[ldk[2],1], x[ldk[2],2], pch=20, col="red", cex=4)
  #points(x[ldk[3],1], x[ldk[3],2], pch=20, col="red", cex=4)
  #points(x[ldk[4],1], x[ldk[4],2], pch=20, col="red", cex=4)
  ldk -> ldks2[[i]]
}
leaves$ldk <- ldks2
leaves <- coo_slide(leaves, ldk=1)
leaves <- fgProcrustes(leaves)

#visualize where the landmarks are
plot(x, asp=1, axes=F, xlab="", ylab="")
points(x[ldk[1],1], x[ldk[1],2], pch=20, col="red", cex=2)
points(x[ldk[2],1], x[ldk[2],2], pch=20, col="red", cex=2)
points(x[ldk[3],1], x[ldk[3],2], pch=20, col="red", cex=2)
points(x[ldk[4],1], x[ldk[4],2], pch=20, col="red", cex=2)

#Do the elliptic Fourier analyses and PCA of the results
efou_leaves <- efourier(leaves, norm=F, smooth.it = 0, verbose = T, nb.h = 20)
#Summarize the data by running a principal components analysis (PCA) on it
pca_efou_leaves <- PCA(efou_leaves, center = T) 
#look at the output from the PCA - how much variance do the first few principle components explain?
#e.g. 0.68 under PC1, Cumulative Proportion means the first PC explains/captures 68% of the variation in the data
summary(pca_efou_leaves)

PCcontrib(pca_efou_leaves, nax=c(1:10), sd.r = c(-2, -1, 0, 1, 2), cex.labels=10) 

#get PC scores into dataframe-----
pcscores <- pca_efou_leaves$x
pcscores <- as.data.frame(pcscores)
pcscores$sampid <- rownames(pcscores) 
pcscores.sub <- pcscores[,c(1:10)]
pcscores.sub <- as.data.frame(pcscores.sub)
pcscores.sub$sampid <- rownames(pcscores.sub)

dim(mdf)
dim(pcscores.sub)
masterdf <- merge(mdf,pcscores.sub,by="sampid")
dim(masterdf)

groups_leaves$sampid <- gsub(".jpg","",groups_leaves$image, fixed = T)
masterdf <- merge(groups_leaves,masterdf,by = "sampid")

#write out data ---------------
#write out data
write.csv(masterdf, "../../leaf_trait_data_files/full_data_10pcs_20harmonics_RTleaves.csv")
#write out pca R object
save(pca_efou_leaves, 
     file = "../../output_results/pca_efou_leaves-10pcs_20harmonics_RTleaves.RData")


#Morphometric analysis VISUALIZATION -------------------

#load object needed for plotting previously saved to avoid rerunning all of the code
load("../../output_results/pca_efou_leaves-10pcs_20harmonics_RTleaves.RData")

#see average shapes for X number of PCs
PCcontrib(pca_efou_leaves, nax=c(1:10), sd.r = c(-2, -1, 0, 1, 2), cex.labels=10) 




#********************************************************************
#GET LEAF SIZE AND SHAPE DATA FOR HEATHER AND RACHEL'S LEAVES JOINTLY (2017 and 2014)-----------------------
#********************************************************************
setwd("Leaf_scans-Processed/Final_Images_JPEG-rotated-RTandHW/")

#load all the libraries we will need
library(Momocs)   #for shape analysis
library(ggplot2)  #for graphing
library(dplyr)    #for dataframe manipulation
library(tidyr)    #for dataframe manipulation

#inputting leaf traces----
#read all of the leaf scans in (this reads in all .jpg files in folder we specified above)
jpg.list_leaves <- list.files(pattern = c( ".jpg")) 
#look at the first few lines of the list of jpg files - do these look like the scans you want?
head(jpg.list_leaves)
#get a trace of each leaf
returns_leaves <- import_jpg(jpg.list_leaves, auto.notcentered = T) 
leaves <- Out(returns_leaves)
#getting leaf measurements ----
names(leaves)  ### look at names of all leaves

par(mfrow=c(1,2))

#get object that is xy coords for each leaf
lf <- leaves$coo

#get length and width for all leaves (can input whole Out object into this function)
lw <- coo_lw(leaves) %>% t() %>% as.data.frame() %>% mutate(sampid=row.names(.)) %>% 
  mutate(length_cm=((V1/300)*2.54)) %>% mutate(width_cm=((V2/300)*2.54)) %>% dplyr::select(sampid, length_cm, width_cm) %>% 
  mutate(ratio_LoverW_cm=length_cm/width_cm)


#get solidity measure for each leaf (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.solid = data.frame(sampid = names(lf), solidity_score = NA)
head(lf.solid)
dim(lf.solid)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.solid$solidity_score[i] = coo_solidity(lf[lf.solid$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.solid)


#get area measure (can input whole Out object into this function)
area <- coo_area(leaves) %>% as.data.frame()
colnames(area)[colnames(area) == "."] <- "Area"
area <- area %>% mutate(sampid=row.names(.)) %>% mutate(area_cm2=((Area/(300*300))*(2.54^2))) %>% 
  dplyr::select(sampid,area_cm2)


#get perimeter for each leaf (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.perim = data.frame(sampid = names(lf), perim = NA)
head(lf.perim)
dim(lf.perim)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.perim$perim[i] = coo_perim(lf[lf.perim$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.perim)
#convert to cm
lf.perim <- lf.perim %>% mutate(perim_cm=(perim/300)*2.54) %>% dplyr::select(sampid, perim_cm)


#get calliper measurement for each leaf - longest distance between points (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.calip = data.frame(sampid = names(lf), longest_dist = NA)
head(lf.calip)
dim(lf.calip)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.calip$longest_dist[i] = coo_calliper(lf[lf.calip$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.calip)
#convert to cm
lf.calip <- lf.calip %>% mutate(longest_dist_cm=(longest_dist/300)*2.54) %>% dplyr::select(sampid, longest_dist_cm)

#merge all measurements together into one dataframe--------
mdf <- merge(lw,lf.calip, by="sampid")
mdf <- merge(mdf,lf.perim, by="sampid")
mdf <- merge(mdf,area, by="sampid")
mdf <- merge(mdf,lf.solid, by="sampid")

#getting leaves ready for morphometric analysis----
#read in the Groups file - this has the detailed information about each leaf
groups_leaves <- read.csv("Groups_file.csv", header = T)
#IMPORTANT: the numbers at the end of the next line represent how many columns are in your Groups file 
#e.g. 1:4 means there are 4 columns - change this to match number of columns in your Groups file
groups_leaves <-data.frame(groups_leaves[,1:3])
#look at how the variables are coded - grouping varibales should be "Factor"
str(groups_leaves)
#attach the leaf info to the leaf traces
leaves$fac <- groups_leaves

#look at a panel of all of the leaf outlines - do any look weird?
panel(leaves, col=0, border=1)
#if you need to see the names to find a specific leaf, display the names
#change the size of the print with "cex.names"
panel(leaves, col=0, border=1, names = T, cex.names = 0.1)

#interpolate n coordinates along the perimeter of the trace
leaves$coo <- lapply(leaves$coo, coo_interpolate, n=500)
#center all of the traces on the origin/on top of each other 
leaves <- coo_center(leaves)
#scale the traces
leaves$coo <- lapply(leaves$coo, FUN=function(x)(coo_scale(x, max(x[,2]/1000)))) 
leaves <- coo_scale(leaves, 1000)

#put landmarks on each leaf, so they all end up in the correct final orientation in the analysis
ldks2 <- vector("list", length=length(leaves$coo))
names(ldks2) <- names(leaves$coo)
centroids <- coo_centpos(leaves) 
for (i in 1:length(leaves)) {
  ldk <- numeric(4)
  leaves$coo[[i]] -> x
  which.min(x[,1]) -> min.x
  which.max(x[,1]) -> max.x
  which.min(x[,2]) -> min.y
  which.max(x[,2]) -> max.y
  p <- c(centroids[i,1], min(x[,2]))
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[1] <- which.min(l)
  p <- c(centroids[i,1], max(x[,2]))
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[2] <- which.min(l)
  p <- c(min(x[,1]), centroids[i,2])
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[3] <- which.min(l)
  p <- c(max(x[,1]), centroids[i,2])
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[4] <- which.min(l)
  #plot(x, asp=1, axes=F, xlab="", ylab="")
  #points(x[ldk[1],1], x[ldk[1],2], pch=20, col="red", cex=4)
  #points(x[ldk[2],1], x[ldk[2],2], pch=20, col="red", cex=4)
  #points(x[ldk[3],1], x[ldk[3],2], pch=20, col="red", cex=4)
  #points(x[ldk[4],1], x[ldk[4],2], pch=20, col="red", cex=4)
  ldk -> ldks2[[i]]
}
leaves$ldk <- ldks2
leaves <- coo_slide(leaves, ldk=1)
leaves <- fgProcrustes(leaves)

#visualize where the landmarks are
plot(x, asp=1, axes=F, xlab="", ylab="")
points(x[ldk[1],1], x[ldk[1],2], pch=20, col="red", cex=2)
points(x[ldk[2],1], x[ldk[2],2], pch=20, col="red", cex=2)
points(x[ldk[3],1], x[ldk[3],2], pch=20, col="red", cex=2)
points(x[ldk[4],1], x[ldk[4],2], pch=20, col="red", cex=2)

#Do the elliptic Fourier analyses and PCA of the results
efou_leaves <- efourier(leaves, norm=F, smooth.it = 0, verbose = T, nb.h = 20)
#Summarize the data by running a principal components analysis (PCA) on it
pca_efou_leaves <- PCA(efou_leaves, center = T) 
#look at the output from the PCA - how much variance do the first few principle components explain?
#e.g. 0.68 under PC1, Cumulative Proportion means the first PC explains/captures 68% of the variation in the data
summary(pca_efou_leaves)

PCcontrib(pca_efou_leaves, nax=c(1:10), sd.r = c(-2, -1.5, -1, 0, 1, 1.5, 2), cex.labels=10) 

#get PC scores into dataframe-----
pcscores <- pca_efou_leaves$x
pcscores <- as.data.frame(pcscores)
pcscores$sampid <- rownames(pcscores) 
pcscores.sub <- pcscores[,c(1:10)]
pcscores.sub <- as.data.frame(pcscores.sub)
pcscores.sub$sampid <- rownames(pcscores.sub)

dim(mdf)
dim(pcscores.sub)
masterdf <- merge(mdf,pcscores.sub,by="sampid")
dim(masterdf)

groups_leaves$sampid <- gsub(".jpg","",groups_leaves$image, fixed = T)
masterdf <- merge(groups_leaves,masterdf,by = "sampid")

#write out data --------
#write out data
write.csv(masterdf, "../../leaf_trait_data_files/full_data_10pcs_20harmonics_RTandHWleaves.csv")
#write out pca R object
save(pca_efou_leaves, 
     file = "../../output_results/pca_efou_leaves-10pcs_20harmonics_RTandHWleaves.RData")


#Morphometric analysis VISUALIZATION -------------------

#load object needed for plotting previously saved to avoid rerunning all of the code
load("../../output_results/pca_efou_leaves-10pcs_20harmonics_RTandHWleaves.RData")

#Load the library ggplot2 for graphing/viewing our results
library(ggplot2)

#see average shapes for X number of PCs
PCcontrib(pca_efou_leaves, nax=c(1:10), sd.r = c(-2, -1, 0, 1, 2), cex.labels=10) 



#********************************************************************
#GET LEAF SIZE AND SHAPE DATA FOR HEATHER'S LEAVES ONLY LEAF 3 (2017) -----------------------
#********************************************************************

#set the location to the folder where your leaf scans and Groups file live on your computer
setwd("Leaf_scans-Processed/Final_Images_JPEG-rotated-onlyleaf3/")

#load all the libraries we will need
library(Momocs)   #for shape analysis
library(ggplot2)  #for graphing
library(dplyr)    #for dataframe manipulation
library(tidyr)    #for dataframe manipulation

#inputting leaf traces----
#get the library that we will use for morphometric analysis
#( NOTE - the first time you use a new library run: e.g. > install.packages("Momocs") )
library(Momocs)

#read all of the leaf scans in (this reads in all .jpg files in folder we specified above)
jpg.list_leaves <- list.files(pattern = c( ".jpg")) 
#look at the first few lines of the list of jpg files - do these look like the scans you want?
head(jpg.list_leaves)
#get a trace of each leaf
returns_leaves <- import_jpg(jpg.list_leaves, auto.notcentered = T) 
leaves <- Out(returns_leaves)
#getting leaf measurements ----
names(leaves)  ### look at names of all leaves

par(mfrow=c(1,2))

#get object that is xy coords for each leaf
lf <- leaves$coo

#get length and width for all leaves (can input whole Out object into this function)
lw <- coo_lw(leaves) %>% t() %>% as.data.frame() %>% mutate(sampid=row.names(.)) %>% 
  mutate(length_cm=((V1/300)*2.54)) %>% mutate(width_cm=((V2/300)*2.54)) %>% select(sampid,length_cm,width_cm) %>% 
  mutate(ratio_LoverW_cm=length_cm/width_cm)


#get solidity measure for each leaf (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.solid = data.frame(sampid = names(lf), solidity_score = NA)
head(lf.solid)
dim(lf.solid)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.solid$solidity_score[i] = coo_solidity(lf[lf.solid$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.solid)


#get area measure (can input whole Out object into this function)
area <- coo_area(leaves) %>% as.data.frame()
colnames(area)[colnames(area) == "."] <- "Area"
area <- area %>% mutate(sampid=row.names(.)) %>% mutate(area_cm2=((Area/(300*300))*(2.54^2))) %>% 
  select(sampid,area_cm2)


#get perimeter for each leaf (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.perim = data.frame(sampid = names(lf), perim = NA)
head(lf.perim)
dim(lf.perim)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.perim$perim[i] = coo_perim(lf[lf.perim$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.perim)
#convert to cm
lf.perim <- lf.perim %>% mutate(perim_cm=(perim/300)*2.54) %>% select(sampid, perim_cm)


#get calliper measurement for each leaf - longest distance between points (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.calip = data.frame(sampid = names(lf), longest_dist = NA)
head(lf.calip)
dim(lf.calip)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.calip$longest_dist[i] = coo_calliper(lf[lf.calip$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.calip)
#convert to cm
lf.calip <- lf.calip %>% mutate(longest_dist_cm=(longest_dist/300)*2.54) %>% select(sampid, longest_dist_cm)

#merge all measurements together into one dataframe--------
mdf <- merge(lw,lf.calip, by="sampid")
mdf <- merge(mdf,lf.perim, by="sampid")
mdf <- merge(mdf,area, by="sampid")
mdf <- merge(mdf,lf.solid, by="sampid")

#getting leaves ready for morphometric analysis----
#read in the Groups file - this has the detailed information about each leaf
groups_leaves <- read.csv("Groups_file.csv", header = T)
#IMPORTANT: the numbers at the end of the next line represent how many columns are in your Groups file 
#e.g. 1:4 means there are 4 columns - change this to match number of columns in your Groups file
groups_leaves <-data.frame(groups_leaves)
#look at how the variables are coded - grouping varibales should be "Factor"
str(groups_leaves)
#change varibales to factors
groups_leaves[,'leaf'] <- as.factor(groups_leaves[,'leaf'])
groups_leaves[,'envelope'] <- as.factor(groups_leaves[,'envelope'])
#check if all variables that you want to be Factors are now
str(groups_leaves)
#attach the leaf info to the leaf traces
leaves$fac <- groups_leaves

#look at a panel of all of the leaf outlines - do any look weird?
panel(leaves, col=0, border=1)
#if you need to see the names to find a specific leaf, display the names
#change the size of the print with "cex.names"
panel(leaves, col=0, border=1, names = T, cex.names = 0.1)

#interpolate n coordinates along the perimeter of the trace
leaves$coo <- lapply(leaves$coo, coo_interpolate, n=500)
#center all of the traces on the origin/on top of each other 
leaves <- coo_center(leaves)
#scale the traces
leaves$coo <- lapply(leaves$coo, FUN=function(x)(coo_scale(x, max(x[,2]/1000)))) 
leaves <- coo_scale(leaves, 1000)

#put landmarks on each leaf, so they all end up in the correct final orientation in the analysis
ldks2 <- vector("list", length=length(leaves$coo))
names(ldks2) <- names(leaves$coo)
centroids <- coo_centpos(leaves) 
for (i in 1:length(leaves)) {
  ldk <- numeric(4)
  leaves$coo[[i]] -> x
  which.min(x[,1]) -> min.x
  which.max(x[,1]) -> max.x
  which.min(x[,2]) -> min.y
  which.max(x[,2]) -> max.y
  p <- c(centroids[i,1], min(x[,2]))
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[1] <- which.min(l)
  p <- c(centroids[i,1], max(x[,2]))
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[2] <- which.min(l)
  p <- c(min(x[,1]), centroids[i,2])
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[3] <- which.min(l)
  p <- c(max(x[,1]), centroids[i,2])
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[4] <- which.min(l)
  #plot(x, asp=1, axes=F, xlab="", ylab="")
  #points(x[ldk[1],1], x[ldk[1],2], pch=20, col="red", cex=4)
  #points(x[ldk[2],1], x[ldk[2],2], pch=20, col="red", cex=4)
  #points(x[ldk[3],1], x[ldk[3],2], pch=20, col="red", cex=4)
  #points(x[ldk[4],1], x[ldk[4],2], pch=20, col="red", cex=4)
  ldk -> ldks2[[i]]
}
leaves$ldk <- ldks2
leaves <- coo_slide(leaves, ldk=1)
leaves <- fgProcrustes(leaves)

#visualize where the landmarks are
plot(x, asp=1, axes=F, xlab="", ylab="")
points(x[ldk[1],1], x[ldk[1],2], pch=20, col="red", cex=2)
points(x[ldk[2],1], x[ldk[2],2], pch=20, col="red", cex=2)
points(x[ldk[3],1], x[ldk[3],2], pch=20, col="red", cex=2)
points(x[ldk[4],1], x[ldk[4],2], pch=20, col="red", cex=2)

#Do the elliptic Fourier analyses and PCA of the results
efou_leaves <- efourier(leaves, norm=F, smooth.it = 0, verbose = T, nb.h = 20)
#Summarize the data by running a principal components analysis (PCA) on it
pca_efou_leaves <- PCA(efou_leaves, center = T) 
#look at the output from the PCA - how much variance do the first few principle components explain?
#e.g. 0.68 under PC1, Cumulative Proportion means the first PC explains/captures 68% of the variation in the data
summary(pca_efou_leaves)

#get PC scores into dataframe-----
pcscores <- pca_efou_leaves$x
pcscores <- as.data.frame(pcscores)
pcscores$sampid <- rownames(pcscores) 
pcscores.sub <- pcscores[,c(1:10)]
pcscores.sub <- as.data.frame(pcscores.sub)
pcscores.sub$sampid <- rownames(pcscores.sub)

dim(mdf)
dim(pcscores.sub)
masterdf <- merge(mdf,pcscores.sub,by="sampid")
dim(masterdf)

groups_leaves$sampid <- gsub(".jpg","",groups_leaves$image, fixed = T)
masterdf <- merge(groups_leaves,masterdf,by = "sampid")
dim(masterdf)

#write out data ------------
#write data
write.csv(masterdf, "../../leaf_trait_data_files/full_data_10pcs_20harmonics_justheathersleavesonlyleaf3.csv")
#write pca R object
save(pca_efou_leaves, 
     file = "../../output_results/pca_efou_leaves-10pcs_20harmonics_justheathersleavesonlyleaf3.RData")


#Morphometric analysis VISUALIZATION -------------------
load("../../output_results/pca_efou_leaves-10pcs_20harmonics_justheathersleavesonlyleaf3.RData")

#Load the library ggplot2 for graphing/viewing our results
library(ggplot2)

#see average shapes for X number of PCs
PCcontrib(pca_efou_leaves, nax=c(1:10), sd.r = c(-2, -1, 0, 1, 2), cex.labels=10) 

#for *PUBLICATION* - SUPP. FIG. 3 - make PC average shapes pretty -----------------
shps <- PCcontrib(pca_efou_leaves, nax=c(1:10), sd.r = c(-2, -1, 0, 1, 2), cex.labels=10) 
shps <- shps$shp

shps <- data.frame(do.call(rbind.data.frame, shps)) 
shps <- shps %>% mutate(sd = shp1, pc = nax) %>% dplyr::select(x,y,pc,sd)
shps$sd[shps$sd == 1] = -2
shps$sd[shps$sd == 2] = -1
shps$sd[shps$sd == 3] = 0
shps$sd[shps$sd == 4] = 1
shps$sd[shps$sd == 5] = 2
shps$pc[shps$pc == 1] = "1\\n(40.4%)"
shps$pc[shps$pc == 2] = "2\\n(23.9%)"
shps$pc[shps$pc == 3] = "3\\n(9.0%)"
shps$pc[shps$pc == 4] = "4\\n(7.2%)"
shps$pc[shps$pc == 5] = "5\\n(4.6%)"
shps$pc[shps$pc == 6] = "6\\n(3.0%)"
shps$pc[shps$pc == 7] = "7\\n(2.0%)"
shps$pc[shps$pc == 8] = "8\\n(1.5%)"
shps$pc[shps$pc == 9] = "9\\n(1.2%)"
shps$pc[shps$pc == 10] = "10\\n(1.0%)"
shps$pc_ordered = factor(shps$pc, levels=c("1\\n(40.4%)","2\\n(23.9%)","3\\n(9.0%)","4\\n(7.2%)","5\\n(4.6%)",
                                           "6\\n(3.0%)","7\\n(2.0%)","8\\n(1.5%)","9\\n(1.2%)","10\\n(1.0%)"))

shps %>%
  ggplot(aes(x = x, y = y)) +
  geom_polygon(fill = "#8BC34A", alpha = 0.9, colour = "black", size = 0.5)  +
  labs(x = "St. devs. from mean shape", y = "Principal component") +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text = element_blank(),
        axis.title.y = element_text(size = 10, margin = margin(r = -3)),
        axis.title.x = element_text(size = 10, margin = margin(t = -3)), 
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(colour=NA, fill=NA),
        strip.text.y = element_text(angle = 180, size = 7),
        strip.text.x = element_text(size = 7)) +
  coord_fixed() +
  facet_grid(pc_ordered ~ sd, switch = "both")

ggsave(paste("../../figures_and_tables/avgleaves_pca_heathersleaf3","pdf", sep = ".") , 
       width = 6.7, height = 16.5, dpi = 600, units = c("cm"))



#********************************************************************
#GET LEAF SIZE AND SHAPE DATA FOR HEATHER'S LEAVES ONLY LEAF 2 (2017) -----------------------
#********************************************************************

#set the location to the folder where your leaf scans and Groups file live on your computer
setwd("Leaf_scans-Processed/Final_Images_JPEG-rotated-onlyleaf2/")

#load all the libraries we will need
library(Momocs)   #for shape analysis
library(ggplot2)  #for graphing
library(dplyr)    #for dataframe manipulation
library(tidyr)    #for dataframe manipulation

#inputting leaf traces----
#get the library that we will use for morphometric analysis
#( NOTE - the first time you use a new library run: e.g. > install.packages("Momocs") )
library(Momocs)

#read all of the leaf scans in (this reads in all .jpg files in folder we specified above)
jpg.list_leaves <- list.files(pattern = c( ".jpg")) 
#look at the first few lines of the list of jpg files - do these look like the scans you want?
head(jpg.list_leaves)
#get a trace of each leaf
returns_leaves <- import_jpg(jpg.list_leaves, auto.notcentered = T) 
leaves <- Out(returns_leaves)
#getting leaf measurements ----
names(leaves)  ### look at names of all leaves

par(mfrow=c(1,2))

#get object that is xy coords for each leaf
lf <- leaves$coo

#get length and width for all leaves (can input whole Out object into this function)
lw <- coo_lw(leaves) %>% t() %>% as.data.frame() %>% mutate(sampid=row.names(.)) %>% 
  mutate(length_cm=((V1/300)*2.54)) %>% mutate(width_cm=((V2/300)*2.54)) %>% select(sampid,length_cm,width_cm) %>% 
  mutate(ratio_LoverW_cm=length_cm/width_cm)


#get solidity measure for each leaf (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.solid = data.frame(sampid = names(lf), solidity_score = NA)
head(lf.solid)
dim(lf.solid)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.solid$solidity_score[i] = coo_solidity(lf[lf.solid$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.solid)


#get area measure (can input whole Out object into this function)
area <- coo_area(leaves) %>% as.data.frame()
colnames(area)[colnames(area) == "."] <- "Area"
area <- area %>% mutate(sampid=row.names(.)) %>% mutate(area_cm2=((Area/(300*300))*(2.54^2))) %>% 
  select(sampid,area_cm2)


#get perimeter for each leaf (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.perim = data.frame(sampid = names(lf), perim = NA)
head(lf.perim)
dim(lf.perim)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.perim$perim[i] = coo_perim(lf[lf.perim$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.perim)
#convert to cm
lf.perim <- lf.perim %>% mutate(perim_cm=(perim/300)*2.54) %>% select(sampid, perim_cm)


#get calliper measurement for each leaf - longest distance between points (cannot input whole Out object so loop thru each leaf in lf)
# make blank data frame to write output to with unique leaf ids
lf.calip = data.frame(sampid = names(lf), longest_dist = NA)
head(lf.calip)
dim(lf.calip)
# loop through and get solidity measure for each leaf in lf
for(i in 1:length(lf)) {
  lf.calip$longest_dist[i] = coo_calliper(lf[lf.calip$sampid[i]])
  message(paste("Iteration", i, "complete"))
}
#check that it worked
head(lf.calip)
#convert to cm
lf.calip <- lf.calip %>% mutate(longest_dist_cm=(longest_dist/300)*2.54) %>% select(sampid, longest_dist_cm)

#merge all measurements together into one dataframe--------
mdf <- merge(lw,lf.calip, by="sampid")
mdf <- merge(mdf,lf.perim, by="sampid")
mdf <- merge(mdf,area, by="sampid")
mdf <- merge(mdf,lf.solid, by="sampid")

#getting leaves ready for morphometric analysis----
#read in the Groups file - this has the detailed information about each leaf
groups_leaves <- read.csv("Groups_file.csv", header = T)
groups_leaves <-data.frame(groups_leaves)
#look at how the variables are coded - grouping varibales should be "Factor"
str(groups_leaves)
#change varibales to factors
groups_leaves[,'leaf'] <- as.factor(groups_leaves[,'leaf'])
groups_leaves[,'envelope'] <- as.factor(groups_leaves[,'envelope'])
#check if all variables that you want to be Factors are now
str(groups_leaves)
#attach the leaf info to the leaf traces
leaves$fac <- groups_leaves

#look at a panel of all of the leaf outlines - do any look weird?
panel(leaves, col=0, border=1)
#if you need to see the names to find a specific leaf, display the names
#change the size of the print with "cex.names"
panel(leaves, col=0, border=1, names = T, cex.names = 0.1)

#interpolate n coordinates along the perimeter of the trace
leaves$coo <- lapply(leaves$coo, coo_interpolate, n=500)
#center all of the traces on the origin/on top of each other 
leaves <- coo_center(leaves)
#scale the traces
leaves$coo <- lapply(leaves$coo, FUN=function(x)(coo_scale(x, max(x[,2]/1000)))) 
leaves <- coo_scale(leaves, 1000)

#put landmarks on each leaf, so they all end up in the correct final orientation in the analysis
ldks2 <- vector("list", length=length(leaves$coo))
names(ldks2) <- names(leaves$coo)
centroids <- coo_centpos(leaves) 
for (i in 1:length(leaves)) {
  ldk <- numeric(4)
  leaves$coo[[i]] -> x
  which.min(x[,1]) -> min.x
  which.max(x[,1]) -> max.x
  which.min(x[,2]) -> min.y
  which.max(x[,2]) -> max.y
  p <- c(centroids[i,1], min(x[,2]))
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[1] <- which.min(l)
  p <- c(centroids[i,1], max(x[,2]))
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[2] <- which.min(l)
  p <- c(min(x[,1]), centroids[i,2])
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[3] <- which.min(l)
  p <- c(max(x[,1]), centroids[i,2])
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[4] <- which.min(l)
  #plot(x, asp=1, axes=F, xlab="", ylab="")
  #points(x[ldk[1],1], x[ldk[1],2], pch=20, col="red", cex=4)
  #points(x[ldk[2],1], x[ldk[2],2], pch=20, col="red", cex=4)
  #points(x[ldk[3],1], x[ldk[3],2], pch=20, col="red", cex=4)
  #points(x[ldk[4],1], x[ldk[4],2], pch=20, col="red", cex=4)
  ldk -> ldks2[[i]]
}
leaves$ldk <- ldks2
leaves <- coo_slide(leaves, ldk=1)
leaves <- fgProcrustes(leaves)

#visualize where the landmarks are
plot(x, asp=1, axes=F, xlab="", ylab="")
points(x[ldk[1],1], x[ldk[1],2], pch=20, col="red", cex=2)
points(x[ldk[2],1], x[ldk[2],2], pch=20, col="red", cex=2)
points(x[ldk[3],1], x[ldk[3],2], pch=20, col="red", cex=2)
points(x[ldk[4],1], x[ldk[4],2], pch=20, col="red", cex=2)

#Do the elliptic Fourier analyses and PCA of the results
efou_leaves <- efourier(leaves, norm=F, smooth.it = 0, verbose = T, nb.h = 20)
#Summarize the data by running a principal components analysis (PCA) on it
pca_efou_leaves <- PCA(efou_leaves, center = T) 
#look at the output from the PCA - how much variance do the first few principle components explain?
#e.g. 0.68 under PC1, Cumulative Proportion means the first PC explains/captures 68% of the variation in the data
summary(pca_efou_leaves)

#get PC scores into dataframe-----
pcscores <- pca_efou_leaves$x
pcscores <- as.data.frame(pcscores)
pcscores$sampid <- rownames(pcscores) 
pcscores.sub <- pcscores[,c(1:10)]
pcscores.sub <- as.data.frame(pcscores.sub)
pcscores.sub$sampid <- rownames(pcscores.sub)

dim(mdf)
dim(pcscores.sub)
masterdf <- merge(mdf,pcscores.sub,by="sampid")
dim(masterdf)

groups_leaves$sampid <- gsub(".jpg","",groups_leaves$image, fixed = T)
masterdf <- merge(groups_leaves,masterdf,by = "sampid")
dim(masterdf)

#write out data ------------
#write data
write.csv(masterdf, "../../leaf_trait_data_files/full_data_10pcs_20harmonics_justheathersleavesonlyleaf2.csv")
#write pca R object
save(pca_efou_leaves, 
     file = "../../output_results/pca_efou_leaves-10pcs_20harmonics_justheathersleavesonlyleaf2.RData")


#Morphometric analysis VISUALIZATION -------------------

#load object needed for plotting previously saved to avoid rerunning all of the code
load("../../output_results/pca_efou_leaves-10pcs_20harmonics_justheathersleavesonlyleaf3.RData")

#Load the library ggplot2 for graphing/viewing our results
library(ggplot2)

#see average shapes for X number of PCs
PCcontrib(pca_efou_leaves, nax=c(1:10), sd.r = c(-2, -1, 0, 1, 2), cex.labels=10) 

