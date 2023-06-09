library(dplyr)
library(reshape2)
library(tools)
library(ggplot2)

# This script generates Figure 3 (Median cluster sizes of evolved strains at the time-point where strains reached maximum size)
# Run in R version 3.5.1
# Set Master Folder Directory of raw data and Output folder
# MODIFY TO MATCH YOUR DIRECTORY STRUCTURE
master_dir <- ("/Users/matt/Documents/Dropbox/LifeCycleRevision/ScientificReports/Dryad/Figure3Data")
output <- setwd("/Users/matt/Documents/Dropbox/LifeCycleRevision/ScientificReports/Dryad//Figure3Output/")

# List the files (with full path) within master directory
file_list <- list.files(path = master_dir, full.names=T, recursive=F)

# This loop does all the heavy lifting to modify and save each .csv  
    for (j in 1:length(file_list)) {
    # save the csv filename without path or extention
    #  this is important in writing the Strain and naming the output file
    csv_name <- file_path_sans_ext(basename(file_list[j]))
    print(csv_name)
    table <- read.csv(file_list[j], header = T)
    # After reading in the csv, add a last column that divides
    #  the previously last column by 255 to get cells per cluster
    table[,ncol(table)+1] <- (table[,ncol(table)])/255
    colnames(table)[ncol(table)] <- "CellsperCluster"
    # Remove all but the last two columns of the data frame
    #  the last two are pixel intensity and cells per cluster
    table <- table[,-1:-(ncol(table)-2)]
    # Convert all zeroes to NAs and remove the NAs
    table[table==0] <- NA
    table <- na.omit(table)
    # Add a column that is the technical replicate (sampled from the file name)
    # note, starting at string position 8 and ending at string position 8
    table[,ncol(table)+1] <- substr(csv_name,8,8)
    colnames(table)[ncol(table)] <- "TechRep"
    # Add a column that is the Strain name.
    # note, starting at string position 1 and ending at string position 5
    table[,ncol(table)+1] <- substr(csv_name,1,5)
    colnames(table)[ncol(table)] <- "Strain"
    # Write the CSV to the output folder
    write.csv(x= table, paste(output, paste(csv_name, "csv", sep="."), sep="/"))
    }


# The following appends all csv in the output file with the same character prefix
#  In this case, it appends all csv of the csv with a '.' in position 3 
#  (which is all of them)
file.prefixes <- unique(sapply(list.files(path=output, pattern="*.csv"),
                               substr, 3,3))
file.list <- lapply(file.prefixes, function(x)list.files(pattern=paste("^",x,".*.csv",sep=""), path=output))
names(file.list) <- file.prefixes                         
tables <- lapply(file.list, function(filenames)lapply(filenames, function(file)read.csv(file, header=TRUE)))
joined.tables <- lapply(tables, function(t) do.call(rbind, t))
# Save table as CSV                                                                
write.csv(joined.tables, paste("Concat_CleanedData_atMax", ".csv", sep=""))

# Now going to see whether calculations from these raw data files matches
# what I found from calculations from Jill's concatenated data set

# MODIFY TO MATCH YOUR DIRECTORY STRUCTURE
Concat_data <- read.csv(file = "/Users/matt/Documents/Dropbox/LifeCycleRevision/ScientificReports/Dryad//Figure3Output/Concat_CleanedData_atMax.csv", header = T)
Concat_data <- Concat_data[,-1:-3]
colnames(Concat_data) <- c("CellsPerCluster","TechRep","Strain")

Concat_data <- subset(Concat_data, Concat_data[,1] < 256)

# This computes the median cluster size for each image and arranges it
data_median <- aggregate(CellsPerCluster ~ Strain + TechRep, 
                         FUN = median, data = Concat_data)
data_median <- arrange(data_median, Strain, TechRep)

# This computes the mean of the medians (collapsing on tech rep)
data_meanofmedian <- aggregate(CellsPerCluster ~ Strain, 
                               FUN = mean, data = data_median)
data_meanofmedian <- arrange(data_meanofmedian, Strain)

data_meanofmedian
# Great, looks like the data matches the results from Jill's concatenated data set

#==================================================================
# Plotting of the data as previously done in analyses

# Double-check that max is 256 cells
subset(Concat_data, Concat_data[,1] == max(Concat_data$CellsPerCluster))
# Great, largest is B2.11 tech rep 5 with 247 cells in a cluster

# Modify the string in character column for desired format (e.g. B2-01)
Concat_data$Strain <- as.character(Concat_data$Strain)
Concat_data$Strain <- gsub("[.]", "-", Concat_data$Strain)

# Compute the median cells per cluster for each strain at this max timepoint
data_medians <- aggregate(CellsPerCluster ~ Strain + TechRep,
                          FUN = median, data = Concat_data)
data_medians <- arrange(data_medians, Strain, TechRep)
# Now I have the median cluster size for each technical replicate of each strain
#  for the cleaned ROI
#   note that ROI were only cleaned for CSVs at the original timepoint of max size
#    but I think we can safely assume since there was little change (2 increased in size)
#    that these are still the timepoints of max size

# For plotting order below...
# Calculate the mean of medians at max timepoint
data_meanofmedians <- aggregate(CellsPerCluster ~ Strain,
                                FUN = mean, data = data_medians)
data_meanofmedians <- arrange(data_meanofmedians, Strain)
colnames(data_meanofmedians)[2] <- "Means"
# And join into data
data_medians <- inner_join(data_medians, data_meanofmedians, by="Strain")

# These steps are necessary for plot formatting so the technical
#  replicates are plotted in order of mean but separately and not combined
#  as a single factor
data_medians <- data_medians[
  order(data_medians[,1],
        data_medians[,3]),]
data_medians$Label2 <- rep(seq(1:45))

pos <- position_dodge(width=0.4)

is.factor(data_medians$Strain)
is.character(data_medians$Strain)
levels(data_medians$Strain)

# Plot order by propagule size from timelapse plot
data_medians$Strain <- factor(data_medians$Strain,
                              levels=c('K1-01','K1-06','B5-06',
                                       'B2-10','B2-01','B5-05',
                                       'B2-11','B2-03','B2-04'))
levels(data_medians$Strain)

median_plot <- ggplot(data = data_medians,
                      aes(x=Strain, y=CellsPerCluster, color=factor(Label2))) +
  geom_point(position = pos, shape = 21, size = 2, stroke =.7) +
  scale_colour_manual(values=rep('black', times = 45)) +
  labs(x = "Strain", y="Median Cells per Cluster at \\nTime of Max Size", color = 'Replicate') + 
  scale_y_continuous(limits=c(0,17), expand = c(0,0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.line = element_line(colour = 'black'),
        text = element_text(size=20),
        axis.text.x = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black')) +
  theme(legend.key = element_blank(),
        legend.position = 'none') +
  annotate("rect",
           xmin=c(1.5,3.5,5.5,7.5),
           xmax=c(2.5,4.5,6.5,8.5),
           ymin=0, ymax=Inf, alpha=0.3, fill='grey') +
  geom_point(position = pos, shape = 21, size = 2, stroke =.7)

# MODIFY TO MATCH YOUR DIRECTORY STRUCTURE
ggsave(filename=paste("/Users/matt/Documents/Dropbox/LifeCycleRevision/ScientificReports/Dryad//Figure3Output/",
                      paste("MedianAtMaxTimepoint_xbyProgeny_transparent", ".png", sep=""), sep=""), plot=median_plot,
       width = 10, height = 7, bg = "transparent")


