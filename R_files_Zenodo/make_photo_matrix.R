# File:     make_photo_matrix.R
# Author:   Sylvia Tan, James Jackson Potter
# Date:     31 May 2021
#
# Notes:    This script takes a set of photographs taken by the Robotany
#           portable field plot imager and puts them into one big image.
#
#           This script assumes that your WORKING DIRECTORY is set to
#           the location of make_photo_matrix.R, and that there is a set of
#           25 photos (with color-checker card first)
#           or 24 photos (without color checker)
#           in a subdirectory called "data"
#
#           This script assumes that photo file names are sequential, with
#           earlier photographs numerically/alphabetically "higher" on the
#           file list than later photographs.

# Install and load the required package "magick"
install.packages("magick")
library(magick)

# Get photo files
allFiles = list.files(path="data", pattern="*.JPG", all.files=TRUE, full.names=TRUE) # get all file names

# The trolley moves in a snake-like pattern, so the photos must be arranged
# in a strange order... it's ugly, but hardcoding made the most sense, here
if (length(allFiles) == 25){
  col1 <- c(image_read(allFiles[7]),image_read(allFiles[6]),image_read(allFiles[5]),image_read(allFiles[4]),image_read(allFiles[3]),image_read(allFiles[2]))
  col2 <- c(image_read(allFiles[8]),image_read(allFiles[9]),image_read(allFiles[10]),image_read(allFiles[11]),image_read(allFiles[12]),image_read(allFiles[13]))
  col3 <- c(image_read(allFiles[19]),image_read(allFiles[18]),image_read(allFiles[17]),image_read(allFiles[16]),image_read(allFiles[15]),image_read(allFiles[14]))
  col4 <- c(image_read(allFiles[20]),image_read(allFiles[21]),image_read(allFiles[22]),image_read(allFiles[23]),image_read(allFiles[24]),image_read(allFiles[25]))
}else if (length(allFiles) == 24){
  col1 <- c(image_read(allFiles[6]),image_read(allFiles[5]),image_read(allFiles[4]),image_read(allFiles[3]),image_read(allFiles[2]),image_read(allFiles[1]))
  col2 <- c(image_read(allFiles[7]),image_read(allFiles[8]),image_read(allFiles[9]),image_read(allFiles[10]),image_read(allFiles[11]),image_read(allFiles[12]))
  col3 <- c(image_read(allFiles[18]),image_read(allFiles[17]),image_read(allFiles[16]),image_read(allFiles[15]),image_read(allFiles[14]),image_read(allFiles[13]))
  col4 <- c(image_read(allFiles[19]),image_read(allFiles[20]),image_read(allFiles[21]),image_read(allFiles[22]),image_read(allFiles[23]),image_read(allFiles[24]))
}else {
  stop('Unexpected number of photos in data/ directory (there should be 24 or 25). Please investigate.')
}

# Reduce image resolution
# (comment this out to retain the original resolution)
col1 <- image_scale(col1, "300")
col2 <- image_scale(col2, "300")
col3 <- image_scale(col3, "300")
col4 <- image_scale(col4, "300")
  
# Append images in the columns
col1 <- image_append(col1,stack=TRUE)
col2 <- image_append(col2,stack=TRUE)
col3 <- image_append(col3,stack=TRUE)
col4 <- image_append(col4,stack=TRUE)

# Put columns together
photo_matrix <- c(col1,col2,col3,col4)
photo_matrix <- image_append(photo_matrix)
photo_matrix # show results

# Save results
image_write(photo_matrix, "photo_matrix.jpg")