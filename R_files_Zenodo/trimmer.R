#randomly trims the group to the size requested
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

require(stringr)
require(dplyr)

#reads a .csv (with ";"  as separators) file as data input
#only one .csv file should be present in the folder
#the number of rows 

final.number.of.rows = 8

input.file <- read.csv2(list.files()[which(str_detect(list.files(), ".csv"))])

input.file <- filter(input.file, is.finite(Mouse))

if (nrow(input.file) >= final.number.of.rows){
  trimmed <- input.file[c(sample(nrow(input.file), final.number.of.rows, replace = FALSE)), ]}else{
    cat("Number of rows requested is greater or equal to the number of rows in the file")
  }

write.csv2(trimmed, file = "trimmed.csv")