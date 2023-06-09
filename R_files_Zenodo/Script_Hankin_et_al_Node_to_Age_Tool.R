# Script_Hankin_et_al_Node_to_Age_Tool.R
#
# This script will reproduce the tree aging tool presented in 
# Hankin et al. (2018).
#
# Hankin, L.E., Higuera, P.E., Davis, K.T., Dobrowski, S.Z. 2018. Accuracy
# of node and bud-scar counts for aging two dominant conifers in western
# North America. Forest Ecology and Management. In Press.
#
# Forest Ecology and Management DOI: https://doi.org/10.1016/j.foreco.2018.06.001
# Dryad Repository DOI: https://doi.org/10.5061/dryad.63m17n5
#
# File Requirements: 
# (1) Data_AgeTool.csv - data for running the tool.
# (2) Template_NodeData.csv - an input data file with the samples for which age 
#                             estimated will be made. See "Example..." section 
#                             below. 
# *Both of these file need to be in the working directory for the sctip to work.
#
# Created by: L.E. Hankin
# Created: June 2018
# University of Montana, PaleoEcology and Fire Ecology Lab
#
# Contact:
# Lacey Hankin: lacey.hankin@gmail.com
# Philip Higuera: phiguera@umontana.edu
#
# Set working directory ---------------------------------------------------

# The directory below should be changed to reflect your data file location:
setwd("L:\\\\4_archivedData\\\\Hankin_et_al_2018_FEM") 

# Load required package ---------------------------------------------------

# If package has not been previously installed, uncomment the following code:
# install.packages("plyr")
require(plyr)

# Tree Aging Tool ---------------------------------------------------------

# Function for tree age predictions and 95% prediction intervals from node 
# counts. Import data file of node counts, species, and region from your 
# data (see example).
# Regions = "All", "CA" (California), "NR" (Northern Rockies), or 
#           # "SW" (Southwest)
# Species = "PIPO" (ponderosa pine); "PSME" (Douglas-fir)

tree_age <- function(x){
  y <- data.frame()
  x <- x[x$Nodes<23,]
  for(i in 1:nrow(x)){
    dat <- read.csv("Data_AgeTool.csv")
    Nodes <- x[i,1]
    Species <- x[i,2]
    Species <- as.character(Species)
    Region <- x[i,3]
    Region <- as.character(Region)
    dat <- dat[dat$Species==Species,]
    dat <- dat[dat$Region==Region,]
    dat <- dat[dat$Nodes==Nodes,]
    EstTreeAge <- dat$pred
    UprPredInt <- dat$upr
    LwrPredInt <- dat$lwr
    age <- data.frame(Nodes,Species,Region,EstTreeAge,UprPredInt,LwrPredInt)
    y <- rbind(y, age)
  }
  return(y)
}

# Example -----------------------------------------------------------------
#
# Read in the dataset from a .csv file, which includes at least three
# columns to identify the node count, species, and region of each sample (row). 
# These columns must be labeled "Nodes", "Species", and "Region". Other columns
# will not affect running the tool. The age estimates will be appended to the 
# end of this input file and saved. The input file should be closed when running 
# the tool.
# 
# Follow the code below and read in node data in a data file

fileName = "Template_NodeData.csv"  # Replace file name with your file name. 

InputData <- read.csv(fileName) # Read in data
NodeData <- subset(InputData,select=c("Nodes", "Species","Region"))
head(NodeData) # Check that 1st 3 columns are: Nodes, Species, Region

EstTreeAge <- tree_age(NodeData) # Use tree_age function to estimate tree age

TreeAge <- cbind(InputData,EstTreeAge[,4:6])
#TreeAge <- join(InputData,EstTreeAge[,4:6],by=NULL,type ="left",match="first") 
                                            # Merge together data in original 
                                            # data file and data output file

# Save tool output: By default, the results are appended after the last column 
# of data in the input file, and the file is resaved. If you want to save the 
# output data in a new file, then change "fileName" below to a new file name 
# (in quotes).  

write.csv(TreeAge,fileName,row.names=FALSE)
