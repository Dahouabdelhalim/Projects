#--------------------------------------------------------------
# Chapter 2: Data science in biomedicine
#--------------------------------------------------------------

# Set working directory to source file location


# Define the countries
countries <- c("USA", "UK", "JAPAN", "GERMANY", "AUSTRALIA",
               "SPAIN", "ITALY", "INDIA", "CHINA")
n_countries <- length(countries)
years <- 2004:2019

# Creating empty datasets
data_bg_countries <- matrix(0, ncol = n_countries, 
                            nrow = length(years))
data_ds_countries <- matrix(0, ncol = n_countries, 
                            nrow = length(years))
data_cc_countries <- matrix(0, ncol = n_countries, 
                            nrow = length(years))

colnames(data_bg_countries) <- countries
colnames(data_ds_countries) <- countries
colnames(data_cc_countries) <- countries

#--------------------------------------------------------------
# Read "DATA SCIENCE" databases
#--------------------------------------------------------------
ds_UK <- read.table("DataScience_UK.txt", header = TRUE)
ds_USA <- read.table("DataScience_USA.txt", header = TRUE)
ds_ITALY <- read.table("DataScience_Italy.txt", header = TRUE)
ds_SPAIN <- read.table("DataScience_Spain.txt", header = TRUE)
ds_JAPAN <- read.table("DataScience_Japan.txt", header = TRUE)
ds_CHINA <- read.table("DataScience_China.txt", header = TRUE)
ds_GER <- read.table("DataScience_Germany.txt", header = TRUE)
ds_INDIA <- read.table("DataScience_India.txt", header = TRUE)
ds_AUS <- read.table("DataScience_Australia.txt", header = TRUE)

data_ds_countries <- cbind(ds_USA[ds_USA[, 1]%in%years, 2],
                           ds_UK[ds_UK[, 1]%in%years, 2],
                           ds_JAPAN[ds_JAPAN[, 1]%in%years, 2],
                           ds_GER[ds_GER[, 1]%in%years, 2],
                           ds_AUS[ds_AUS[, 1]%in%years, 2], 
                           ds_SPAIN[ds_SPAIN[, 1]%in%years, 2],
                           ds_ITALY[ds_ITALY[, 1]%in%years, 2],
                           ds_INDIA[ds_INDIA[, 1]%in%years, 2],
                           ds_CHINA[ds_CHINA[, 1]%in%years, 2])
colnames(data_ds_countries) <- countries
data_ds_countries <- data_ds_countries[length(years):1, ]
rownames(data_ds_countries) <- years
data_ds_countries

#--------------------------------------------------------------
# Read "BIG DATA" databases
#--------------------------------------------------------------

bg_UK <- read.table("BigData_UK.txt", header = TRUE)
bg_USA <- read.table("BigData_USA.txt", header = TRUE)
bg_ITALY <- read.table("BigData_Italy.txt", header = TRUE)
bg_SPAIN <- read.table("BigData_Spain.txt", header = TRUE)
bg_JAPAN <- read.table("BigData_Japan.txt", header = TRUE)
bg_CHINA <- read.table("BigData_China.txt", header = TRUE)
bg_GER <- read.table("BigData_Germany.txt", header = TRUE)
bg_INDIA <- read.table("BigData_India.txt", header = TRUE)
bg_AUS <- read.table("BigData_Australia.txt", header = TRUE)


data_bg_countries <- cbind(bg_USA[bg_USA[, 1]%in%years, 2],
                           bg_UK[bg_UK[, 1]%in%years, 2],
                           bg_JAPAN[bg_JAPAN[, 1]%in%years, 2],
                           bg_GER[bg_GER[, 1]%in%years, 2],
                           bg_AUS[bg_AUS[, 1]%in%years, 2], 
                           bg_SPAIN[bg_SPAIN[, 1]%in%years, 2],
                           bg_ITALY[bg_ITALY[, 1]%in%years, 2],
                           bg_INDIA[bg_INDIA[, 1]%in%years, 2],
                           bg_CHINA[bg_CHINA[, 1]%in%years, 2])
colnames(data_bg_countries) <- countries
data_bg_countries <- data_bg_countries[length(years):1, ]
rownames(data_bg_countries) <- years
data_bg_countries

#--------------------------------------------------------------
# Read "CLOUD COMPUTING" databases
#--------------------------------------------------------------

cc_UK <- read.table("CloudComputing_UK.txt", header = TRUE)
cc_USA <- read.table("CloudComputing_USA.txt", header = TRUE)
cc_ITALY <- read.table("CloudComputing_Italy.txt", header = TRUE)
cc_SPAIN <- read.table("CloudComputing_Spain.txt", header = TRUE)
cc_JAPAN <- read.table("CloudComputing_Japan.txt", header = TRUE)
cc_CHINA <- read.table("CloudComputing_China.txt", header = TRUE)
cc_GER <- read.table("CloudComputing_Germany.txt", header = TRUE)
cc_INDIA <- read.table("CloudComputing_India.txt", header = TRUE)
cc_AUS <- read.table("CloudComputing_Australia.txt", header = TRUE)

data_cc_countries <- cbind(cc_USA[cc_USA[, 1]%in%years, 2],
                           cc_UK[cc_UK[, 1]%in%years, 2],
                           cc_JAPAN[cc_JAPAN[, 1]%in%years, 2],
                           cc_GER[cc_GER[, 1]%in%years, 2],
                           cc_AUS[cc_AUS[, 1]%in%years, 2], 
                           cc_SPAIN[cc_SPAIN[, 1]%in%years, 2],
                           cc_ITALY[cc_ITALY[, 1]%in%years, 2],
                           cc_INDIA[cc_INDIA[, 1]%in%years, 2],
                           cc_CHINA[cc_CHINA[, 1]%in%years, 2])
colnames(data_cc_countries) <- countries
data_cc_countries <- data_cc_countries[length(years):1, ]
rownames(data_cc_countries) <- years
data_cc_countries

#--------------------------------------------------------------
# Complete database and summary
#--------------------------------------------------------------

data_countries_total <- matrix(NA, ncol = length(countries),
                               nrow = length(years) * 3)
colnames(data_countries_total) <- countries
rownames(data_countries_total) <- 1:(length(years) * 3)
data_countries_total

#--

for(i in 1:length(years)){
 rownames(data_countries_total)[1 + (i - 1) * 3] <- paste("ds", 
                                                          years[i],
                                                          sep = "_")
 data_countries_total[1 + (i - 1) * 3, ] <- data_ds_countries[i, ]
}
data_countries_total

