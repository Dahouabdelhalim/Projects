
##############################################################################
### REFERENCES
### https://cran.r-project.org/web/packages/factoextra/factoextra.pdf
##############################################################################



##############################################################################
########### PCA - Female morphometric dataset ################################
##############################################################################
##############################################################################

### Load dataset. 
# file = : designate the location (path) of the csv-formatted file in your computer exactly; for example, in here the directory of the morphometric file is in "C:/Users/Desktop/Morphometry/" 
# row.names = 1 designate 1st column as specimen ID.
d <- read.csv(file = "C:/Users/Desktop/Morphometry/Revised Supplementary 2. female morphometric dataset.csv", sep = ",", dec = ".", header = TRUE, row.names = 1) 

# Display dataset
d

# Create a new dataset (d2) that excludes the taxon name column. When counting columns, ignore the 1st column in the file “example”, because it is already used as “row.name” (ex. the 27th column in the file “example” must be designated as “26”); d[,-(a:b)]: create a new dataset without the columns from a to b, a and b are numerous.
d2 <- d[,-(27:28)]

# Create a new dataset "taxon_names"
taxon_names <- d[,27]

# Display the new dataset d2
d2

# Display the field "taxon_names"
taxon_names

# Display the number of rows and columns (field of variables) of d2
dim(d2)

#### Principal component analysis (PCA)
# Run the required package(s)
library(factoextra)

res.pca <- prcomp(d2, scale = TRUE)

## Graph PCA 2D-plot
# PCA 2D-plot. Use 'point' to indicate individuals.
fviz_pca_ind(res.pca, geom.ind = "point", pointshape = 21, pointsize = 2, fill.ind = taxon_names, col.ind = "black", 
             palette = c("#de84a7ff", "#5f8dd3ff", "#5fbcd3ff", "#ff3535ff", "#ffd42aff", "#5aa02cff", "#2986cc", "#16537e", "#6a329f", "#c90076", "#cc0000", "#f3e400"), 
             addEllipses = TRUE, label = "var", col.var = "black", repel = TRUE, legend.title = "taxon") +
  ggtitle("2D PCA-plot from morphological dataset") +
  theme(plot.title = element_text(hjust = 0.5))


## Graph of variables (morphometric characters)
# Plot graph of variables which illustrates variables' contributions by colors
fviz_pca_var(res.pca, axes = c(1, 2), col.var = "contrib", gradient.cols = c("white", "blue", "red"), ggtheme = theme_minimal())

### Select the top 10 contributing variables
res.pca1 <- facto_summarize(res.pca, "var", axes = 1:5, select = list(contrib = 15))[,-1]

# Write “res.pca1” as a CSV, file = : designate the location (path) where you want to save the file.
write.csv(res.pca1, file = "C:/Users/Desktop/Morphometry/female_top_10_contributing_variables.csv")



##############################################################################
########### PCA - Male morphometric dataset ##################################
##############################################################################
##############################################################################

### Load dataset. 
# file = : designate the location (path) of the csv-formatted file in your computer exactly; for example, in here the directory of the morphometric file is in "C:/Users/Desktop/Morphometry/" 
# row.names = 1 designate 1st column as specimen ID.
d <- read.csv(file = "C:/Users/Desktop/Morphometry/Revised Supplementary 3. male morphometric dataset.csv", sep = ",", dec = ".", header = TRUE, row.names = 1) 

# Display dataset
d

# Create a new dataset (d2) that excludes the taxon name column. When counting columns, ignore the 1st column in the file “example”, because it is already used as “row.name” (ex. the 27th column in the file “example” must be designated as “26”); d[,-(a:b)]: create a new dataset without the columns from a to b, a and b are numerous.
d2 <- d[,-(28:29)]

# Create a new dataset "taxon_names"
taxon_names <- d[,28]

# Display the new dataset d2
d2

# Display the field "taxon_names"
taxon_names

# Display the number of rows and columns (field of variables) of d2
dim(d2)

#### Principal component analysis (PCA)
# Run the required package(s)
library(factoextra)

res.pca <- prcomp(d2, scale = TRUE)

## Graph PCA 2D-plot
# PCA 2D-plot. Use 'point' to indicate individuals.
fviz_pca_ind(res.pca, geom.ind = "point", pointshape = 21, pointsize = 2, fill.ind = taxon_names, col.ind = "black", 
             palette = c("#de84a7ff", "#5f8dd3ff", "#5fbcd3ff", "#ff3535ff", "#ffd42aff", "#5aa02cff", "#2986cc", "#16537e", "#6a329f", "#c90076", "#cc0000", "#f3e400"), 
             addEllipses = TRUE, label = "var", col.var = "black", repel = TRUE, legend.title = "taxon") +
  ggtitle("2D PCA-plot from morphological dataset") +
  theme(plot.title = element_text(hjust = 0.5))


## Graph of variables (morphometric characters)
# Plot graph of variables which illustrates variables' contributions by colors
fviz_pca_var(res.pca, axes = c(1, 2), col.var = "contrib", gradient.cols = c("white", "blue", "red"), ggtheme = theme_minimal())

### Select the top 10 contributing variables
res.pca1 <- facto_summarize(res.pca, "var", axes = 1:5, select = list(contrib = 15))[,-1]

# Write “res.pca1” as a CSV. file = : designate the location (path) where you want to save the file.
write.csv(res.pca1, file = "C:/Users/Desktop/Morphometry/male_top_10_contributing_variables.csv")


