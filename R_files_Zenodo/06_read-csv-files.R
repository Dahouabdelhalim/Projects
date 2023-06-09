### Read in the data from .csv files 

data.folder <- "data" #This is the folder to hold the raw data

### Read data
na.strings <- c("", "NA") 

# Studies automatically extracted from the literature
att3.wos <- read.csv(paste(data.folder, "Steps_2-3_compilation.csv", sep="/"), header=T, na.strings=na.strings) # file with results from screening (step 2) and eligibility assessment (step 3) - 166 rows
att4.wos <- read.csv(paste(data.folder, "Step4_compilation_revised.csv", sep="/"), header=T, na.strings=na.strings) # file with extracted information from eligible studies (step 4) - 486 rows

# Studies manually included from other sources
att3.no.wos <- read.csv(paste(data.folder, "Steps_2-3_manual-inclusion.csv", sep="/"), header=T, na.strings=na.strings) # file with results from screening (step 2) and eligibility assessment (step 3) - 138 rows
att4.no.wos <- read.csv(paste(data.folder, "Step4_manual-inclusion_revised.csv", sep="/"), header=T, na.strings=na.strings) # file with extracted information from eligible studies (step 4) - 611 rows

# List of indicators
indicators <- read.csv(paste(data.folder, "Indicators.csv", sep="/"), header=T, na.strings=na.strings) #Not copying the original file name, for clarity and to avoid confusion - 33 rows
