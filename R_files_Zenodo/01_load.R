#01_load
#load packages
library(tidyverse)
library(readxl)
library(ggrepel)
library(cowplot)
library(readxl)
library(patchwork)
library(scales)

#load data

D <- read_excel("MASTER copy survey2020 covid19.xlsx", 
                sheet = "Quantitative coded data")
# NOTE: The raw data generated for this study cannot be publicly available due to 
# restrictions imposed by ethics requirements and agreements with research 
# participants [providing such data would compromise the participants]. 
# Selected data that does not breach any ethics requirements are available 
# on request from the corresponding author upon reasonable request

#Convert character to numeric
D<-mutate_if(D, is.character,as.numeric)

#Add some small functions
reorder_size <- function(x) {
  factor(x, levels = names(sort(table(x))))} #useful for the ordering