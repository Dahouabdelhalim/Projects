
#########################Filename: HouseDustMyco_statistics_script_environmental_variables.R

#Publication: Analyzing indoor mycobiomes through a large-scale citizen science study in Norway
#Authors: Pedro M. Martin-Sanchez, Eva Lena F. Estensmo, Luis N. Morgado, Sundy Maurice, Ingeborg B. Engh, Inger Skrede and Håvard Kauserud

##Related to the sections: 
        #Materials and Methods- 2.2. Environmental data; 2.5. Statistical analyses
        #Supplemental Information

## Goals: 
#(1) To test the contribution of the continuous variables by PCA
#(2) To illustrate the geographical variation of some selected climatic variables  

#Inputs: Table of the 52 environmental continuous variables explored ("52_continuous_variables.xlsx")
        #Notes: In order to protect the personal data of citizen scientists, the geographic coordinates of their houses have intentionally been omitted from this input file.
        #The geographic coordinates are needed to run the second part of the script (goal 2). 
#Outputs: 
        #Plots for the Figure S2 (PCA analysis for data from environmental continuous variables)
        #Plots for the Figure 1 and Figure S2 (maps showing the geographical variation of some climatic variables)


#########################
#PREPARATIONS
#########################

#Load packages

library(openxlsx)
library(tidyverse)
library(ade4)
library(factoextra)
library(fhidata) 
library(gridExtra)

#Set your path
setwd("C:/WP1_data_analyses/script_preparation/")

#Load the data of the 52 continuous environmental variables that we want to test by PCA
continuous_variables <-read.xlsx("52_continuous_variables.xlsx")
#Remark: This includes 46 explanatory variables from a recent study modeling the vegetation types in Norway (Horvath et al. 2019),
        #and 6 climatic variables previously extracted from WorldClim 2 (BIO1, BIO4, BIO9, BIO10, BIO11 and BIO12)
        #Data for 267 study houses; 4 houses were removed considering the lack of some data (like the 2 houses from Svalbard)

#Check the data
dim(variables)
str(variables) 


#########################
#Principal Component Analyses (PCA) using the R package ade4 
#########################

#Move the column "house" to the rownames
variables <-continuous_variables %>%
  remove_rownames() %>% 
  column_to_rownames("house")

# PCA using dudi
res.pca <- dudi.pca(variables,
                    scannf = FALSE,   # Hide scree plot
                    nf = 5            # Number of components kept in the results
)

#Scree plot
fviz_eig(res.pca)


#Plotting the individuals. Individuals with a similar profile are grouped together 
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#Graph of variables. Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
######Remark: This plot corresponds to the Supplemental Figure S3a  

#Biplot of individuals and variables
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)


#Visualise using ade4
screeplot(res.pca, main = "Screeplot - Eigenvalues")

#Correlation circle of variables
s.corcircle(res.pca$co)


#s.corcircle(res.pca$co)
s.label(res.pca$li, 
        xax = 1,     # Dimension 1
        yax = 2)     # Dimension 2

# Biplot of individuals and variables
scatter(res.pca,
        posieig = "none", # Hide the scree plot
        clab.row = 0      # Hide row labels
)

# Check the Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation


##### Based on these results we selected 10 variables
selected_variables <- variables %>% select (Growing_season_length, swe_4, sca_2, Total_insolation,BIO1, BIO4, BIO9, BIO10, BIO11, BIO12)

#Check the selected data
dim(selected_variables) 
str(selected_variables) 


# Compute PCA using dudi
res.pca <- dudi.pca(selected_variables,
                    scannf = FALSE,   # Hide scree plot
                    nf = 5            # Number of components kept in the results
)

fviz_eig(res.pca)

#Graph of individuals. Individuals with a similar profile are grouped together 
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#Graph of the selected variables. Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
######Remark: This plot corresponds to the Supplemental Figure S3b

#Biplot of individuals and variables
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

#Visualise using ade4
screeplot(res.pca, main = "Screeplot - Eigenvalues")

#Correlation circle of variables
s.corcircle(res.pca$co)


#s.corcircle(res.pca$co)
s.label(res.pca$li, 
        xax = 1,     # Dimension 1
        yax = 2)     # Dimension 2

# Biplot of individuals and variables
scatter(res.pca,
        posieig = "none", # Hide the scree plot
        clab.row = 0      # Hide row labels
)

# Eigenvalues
eig.val <- get_eigenvalue(res.pca)
eig.val

# Results for Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation


#########################
#Norway maps illustrating the geographical variation of some selected climatic variables 
#########################

### YOU NEED TO ADD THE GEOGRAPHIC CORDINATES (LONGITUDE AND LATITUDE) OF HOUSES IN THE TABLE OF THE VARIABLES "continuous_variables"
### We have intentionally omitted these two columns ("longitude" and "latitude") to protect the personal data of the citizen scientists

#Extract the Norway map for plots
test_maps_norway<-data.frame(fhidata::norway_map_counties$long,fhidata::norway_map_counties$lat, group=fhidata::norway_map_counties$group) %>% 
  rename(long=fhidata..norway_map_counties.long, lat=fhidata..norway_map_counties.lat)

#Examples of plots for the 6 selected Bioclim variables from WorldClim 2 (BIO1, BIO4, BIO9, BIO10, BIO11 and BIO12):    
BIO1_map<-ggplot()+
  geom_polygon(aes(x = long, y = lat, group=group),col="grey", fill=NA,data = test_maps_norway)+
  geom_point(aes(x=longitude,y=latitude,col=BIO1/10),size=2,data = continuous_variables)+
  theme_void()+
  coord_quickmap()+
  scale_color_viridis_c(option = "inferno")+
  theme(legend.position = c(0.6,0.3),
        plot.title = element_text(size = 20))+
  labs(col="Degree Celsius")+
  ggtitle("(a) BIO1")

BIO4_map<-ggplot()+
  geom_polygon(aes(x = long, y = lat, group=group),col="grey", fill=NA,data = test_maps_norway)+
  geom_point(aes(x=longitude,y=latitude,col=BIO4/1000),size=2,data = continuous_variables)+
  theme_void()+
  coord_quickmap()+
  scale_color_viridis_c(option = "inferno")+
  theme(legend.position = c(0.6,0.3),
        plot.title = element_text(size = 20))+
  labs(col="Degree Celsius")+
  ggtitle("(b) BIO4/100")

BIO9_map<-ggplot()+
  geom_polygon(aes(x = long, y = lat, group=group),col="grey", fill=NA,data = test_maps_norway)+
  geom_point(aes(x=longitude,y=latitude,col=BIO9/10),size=2,data = continuous_variables)+
  theme_void()+
  coord_quickmap()+
  scale_color_viridis_c(option = "inferno")+
  theme(legend.position = c(0.6,0.3),
        plot.title = element_text(size = 20))+
  labs(col="Degree Celsius")+
  ggtitle("(c) BIO9")

BIO10_map<-ggplot()+
  geom_polygon(aes(x = long, y = lat, group=group),col="grey", fill=NA,data = test_maps_norway)+
  geom_point(aes(x=longitude,y=latitude,col=BIO10/10),size=2,data = continuous_variables)+
  theme_void()+
  coord_quickmap()+
  scale_color_viridis_c(option = "inferno")+
  theme(legend.position = c(0.6,0.3),
        plot.title = element_text(size = 20))+
  labs(col="Degree Celsius")+
  ggtitle("(d) BIO10")

BIO11_map<-ggplot()+
  geom_polygon(aes(x = long, y = lat, group=group),col="grey", fill=NA,data = test_maps_norway)+
  geom_point(aes(x=longitude,y=latitude,col=BIO11/10),size=2,data = continuous_variables)+
  theme_void()+
  coord_quickmap()+
  scale_color_viridis_c(option = "inferno")+
  theme(legend.position = c(0.6,0.3),
        plot.title = element_text(size = 20))+
  labs(col="Degree Celsius")+
  ggtitle("(e) BIO11")

BIO12_map<-ggplot()+
  geom_polygon(aes(x = long, y = lat, group=group),col="grey", fill=NA,data = test_maps_norway)+
  geom_point(aes(x=longitude,y=latitude,col=BIO12),size=2,data = continuous_variables)+
  theme_void()+
  coord_quickmap()+
  scale_color_viridis_c(option = "inferno")+
  theme(legend.position = c(0.6,0.3),
        plot.title = element_text(size = 20))+
  labs(col="mm")+
  ggtitle("(f) BIO12")

#Combined figure with the 6 maps 
grid.arrange(BIO1_map, BIO4_map, BIO9_map,  
             BIO10_map,BIO11_map,BIO12_map,ncol=3,nrow=2)

######Remark: These plots correspond to the Figure 1b and Supplemental Figure S2

