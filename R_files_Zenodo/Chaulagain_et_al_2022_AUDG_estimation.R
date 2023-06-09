
###############################################################
###### AUDG ESTIMATION #######
############ BY Bhim Chaulagain #################
############### Oregon State University ######################

# install required packages
install.packages(c("knitr", "dplyr", "tidyr", "agricolae", "ggplot2", "magrittr", "ggpubr", "readxl", "reshape", "tidyverse"))

# Load required packages
library('knitr')     
library('dplyr')     
library('tidyr')     
library('agricolae') 
library('ggplot2')   
library('magrittr')  
library("ggpubr")
library("readxl")
library("reshape")
library("xlsx")
library("tidyverse")

# read in data
data <- read_excel("Chaulagain_et_al_2022_Data_AUDG_estimation.xlsx")
data <- as.data.frame(data)
str(data)

#reshape data
test_data_long <- data %>% 
  filter(Year == "2020", Location == "Madras", Type == "Field")%>% 
  melt(id=c("Distance_ft", "Treatment", "Treatment_description", "Year", "Location", "Type"))

##we excluded disease from the center of each plot equivalent to the largest cull     
#size treatment when computing the AUDG values for each treatment.For this, AUDG for each side of the     
#disease gradient (North and South) from the focus after excluding the given area was computed and combined together

############South of plot
test_data_long_S <- test_data_long %>% filter(Distance_ft < 50)
head(test_data_long_S, 10)

# AUDG calculation
audg_data <- test_data_long_S %>%
  group_by(Treatment, variable) %>% # grouping
  summarize(AUDG = audpc(value, Distance_ft))

openxlsx::write.xlsx(audg_data, "AUDG_S_2020_M_F.xlsx", sheetName = "South")

audg_data %>% 
  spread(variable, AUDG) %>% # create a table 
  arrange(Treatment) %>%     # Sort the variable to have in table
  kable() # display as a table 

##############North of plot
test_data_long_N <- test_data_long %>% filter(Distance_ft > 60)

head(test_data_long_N, 10)
# AUDPC calculation
audg_data <- test_data_long_N %>%
  group_by(Treatment, variable) %>% # grouping
  summarize(AUDG = audpc(value, Distance_ft))

openxlsx::write.xlsx(audg_data, "AUDG_N_2020_M_F.xlsx", sheetName = "North")

audg_data %>% 
  spread(variable, AUDG) %>% # create a table 
  arrange(Treatment) %>%      # Sort the variable to have in table
  kable() # display as a table 

#Combine AUDG from both side of the disease gradient for the final AUDG values for each treatment and replications

