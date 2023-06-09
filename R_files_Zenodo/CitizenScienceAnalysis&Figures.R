# The purpose of this script is to use data collected by citizen scientists to determine 
# whether there are significant differences in squash bee abundance between farms/gardens
# using different management practices.

rm(list=ls())
library(tidyverse)
library(readxl)
library(dunn.test)
library(glmmADMB)
library(emmeans)
library(plyr)
library(bbmle)
library(ggplot2)

# Set working directory based on where the associated data is stored
setwd()

# Assign citizen science data to a dataframe object. Open the spreadsheet 
# ("CitizenScienceSquashBeeSurvey.xlsx") and move the "Form Responses 1" tab 
# into the first tab position before executing the following code. 
Resp <- data.frame(read_excel("CitizenScienceSquashBeeSurvey.xlsx"))

# Note: This R script is folded into sections. To open sections and view code
# click on the arrows in the margin next to the section names below, or the green
# double-sided arrow icons next to section names. To refold sections, click once 
# more on the arrow in the margin. 



###### Data Wrangling ######

# Change NA's to 0's
Resp[is.na(Resp)] <- 0

# Group mulch categories into None, Plant Material, Plastic, and Plant Material + Plastic
Resp[Resp$Mulch=="Plastic , Plant Material", "Mulch"] <- "Plastic + Plant Material"
Resp[Resp$Mulch=="Plant Material , Plastic", "Mulch"] <- "Plastic + Plant Material"
Resp[Resp$Mulch=="Plastic , Plant Material , None", "Mulch"] <- "Plastic + Plant Material"
Resp[Resp$Mulch=="None , Plant Material , Plastic", "Mulch"] <- "Plastic + Plant Material"
Resp[Resp$Mulch=="Plant Material , None", "Mulch"] <- "Plant Material"
Resp[Resp$Mulch=="Plastic , None", "Mulch"] <- "Plastic"

# Group tillage categories
Resp[Resp$TillageType=="Full Tillage (100% of soil cultivated)", "TillageType"] <- "Full Tillage"
Resp[Resp$TillageType=="No Tillage (no cultivation)", "TillageType"] <- "No Tillage"
Resp[Resp$TillageType=="No Tillage (no cultivation) , Reduced Tillage (partial cultivation)", "TillageType"] <- "No Tillage"
Resp[Resp$TillageType=="Reduced Tillage (partial cultivation)", "TillageType"] <- "Reduced Tillage"
Resp[Resp$TillageType=="Reduced Tillage (partial cultivation) , No Tillage (no cultivation)", "TillageType"] <- "Reduced Tillage"

# Change tillage depth to approximate cm eqivalent
Resp[Resp$TillageDepth=="0", "TillageDepth"] <- "0"
Resp[Resp$TillageDepth=="1-5", "TillageDepth"] <- "3-14"
Resp[Resp$TillageDepth=="6-10", "TillageDepth"] <- "15-25"
Resp[Resp$TillageDepth=="11-20", "TillageDepth"] <- "25-50"

# Group the < 0.5 acre and 0.5-1 acre cucurbit area
Resp[Resp$CucurbitArea=="< 0.5 acre", "CucurbitArea"] <- "<1 acre"
Resp[Resp$CucurbitArea=="0.5-1 acre", "CucurbitArea"] <- "<1 acre"

# Standardize nearest town names
Resp[Resp$NearestTown=="Grandville, MI", "NearestTown"] <- "Grandville, MI"
Resp[Resp$NearestTown=="Alpena, MI", "NearestTown"] <- "Alpena, MI"
Resp[Resp$NearestTown=="Empire, MI", "NearestTown"] <- "Empire, MI"
Resp[Resp$NearestTown=="New Boston, MI", "NearestTown"] <- "New Boston, MI"
Resp[Resp$NearestTown=="Escanaba", "NearestTown"] <- "Escanaba, MI"
Resp[Resp$NearestTown=="Plainwell MI", "NearestTown"] <- "Plainwell, MI"
Resp[Resp$NearestTown=="Ortonville", "NearestTown"] <- "Ortonville, MI"
Resp[Resp$NearestTown=="Milford", "NearestTown"] <- "Milford, MI"
Resp[Resp$NearestTown=="Allendale, mi", "NearestTown"] <- "Allendale, MI"
Resp[Resp$NearestTown=="Allendale, Mi", "NearestTown"] <- "Allendale, MI"
Resp[Resp$NearestTown=="ALLENDALE, MI", "NearestTown"] <- "Allendale, MI"
Resp[Resp$NearestTown=="Maple city, MI", "NearestTown"] <- "Maple City, MI"
Resp[Resp$NearestTown=="Maple city", "NearestTown"] <- "Maple City, MI"
Resp[Resp$NearestTown=="Maple City", "NearestTown"] <- "Maple City, MI"
Resp[Resp$NearestTown=="Maple city mi", "NearestTown"] <- "Maple City, MI"
Resp[Resp$NearestTown=="Maple city MI", "NearestTown"] <- "Maple City, MI"
Resp[Resp$NearestTown=="Maple City, mi", "NearestTown"] <- "Maple City, MI"
Resp[Resp$NearestTown=="Clio", "NearestTown"] <- "Clio, MI"
Resp[Resp$NearestTown=="Clio, mi", "NearestTown"] <- "Clio, MI"
Resp[Resp$NearestTown=="Clio, Michigan", "NearestTown"] <- "Clio, MI"
Resp[Resp$NearestTown=="Clio, MI", "NearestTown"] <- "Clio, MI"
Resp[Resp$NearestTown=="Clio,  MI", "NearestTown"] <- "Clio, MI"
Resp[Resp$NearestTown=="Holland MI", "NearestTown"] <- "Holland, MI"
Resp[Resp$NearestTown=="Holland", "NearestTown"] <- "Holland, MI"
Resp[Resp$NearestTown=="49423", "NearestTown"] <- "Holland, MI"
Resp[Resp$NearestTown=="Holland Mi", "NearestTown"] <- "Holland, MI"
Resp[Resp$NearestTown=="Holland, mi", "NearestTown"] <- "Holland, MI"
Resp[Resp$NearestTown=="Grand Rapids Michigan", "NearestTown"] <- "Grand Rapids, MI"
Resp[Resp$NearestTown=="Grand Rapids MI", "NearestTown"] <- "Grand Rapids, MI"
Resp[Resp$NearestTown=="Grand Rapids", "NearestTown"] <- "Grand Rapids, MI"
Resp[Resp$NearestTown=="Georgetown, Indiana", "NearestTown"] <- "Georgetown, IN"
Resp[Resp$NearestTown=="Georgetown, In", "NearestTown"] <- "Georgetown, IN"
Resp[Resp$NearestTown=="Detroit", "NearestTown"] <- "Detroit, MI"
Resp[Resp$NearestTown=="Novi, MI", "NearestTown"] <- "Novi, MI"
Resp[Resp$NearestTown=="Novi", "NearestTown"] <- "Novi, MI"
Resp[Resp$NearestTown=="novi mi", "NearestTown"] <- "Novi, MI"
Resp[Resp$NearestTown=="novi, mi", "NearestTown"] <- "Novi, MI"
Resp[Resp$NearestTown=="Lewiston", "NearestTown"] <- "Lewiston, MI"
Resp[Resp$NearestTown=="Northville", "NearestTown"] <- "Northville, MI"
Resp[Resp$NearestTown=="Mt Pleasant", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="mt pleasant", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="mt pleasant, mi", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="Mt Pleasant, Mi", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="Mt Pleasant,Mi", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="Mt pleasant mi", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="Mt pleasant, mi", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="Mt Pleasant mi", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="Mt pleasant mo", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="Mt Pleasant Mi", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="Mt pleasant", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="Mt Pleasant,  Mi", "NearestTown"] <- "Mt Pleasant, MI"
Resp[Resp$NearestTown=="Au Gres,MI", "NearestTown"] <- "Au Gres, MI"
Resp[Resp$NearestTown=="Au Gres", "NearestTown"] <- "Au Gres, MI"
Resp[Resp$NearestTown=="Sutton's bay", "NearestTown"] <- "Sutton's Bay, MI"
Resp[Resp$NearestTown=="Suttons Bay", "NearestTown"] <- "Sutton's Bay, MI"
Resp[Resp$NearestTown=="Suttons Bay, MI", "NearestTown"] <- "Sutton's Bay, MI"
Resp[Resp$NearestTown=="Traverse city", "NearestTown"] <- "Traverse City, MI"
Resp[Resp$NearestTown=="Traverse city, Mi", "NearestTown"] <- "Traverse City, MI"
Resp[Resp$NearestTown=="Traverse City, Mi", "NearestTown"] <- "Traverse City, MI"
Resp[Resp$NearestTown=="Blair, Traverse City", "NearestTown"] <- "Traverse City, MI"
Resp[Resp$NearestTown=="Blair, Traverse City, MI", "NearestTown"] <- "Traverse City, MI"
Resp[Resp$NearestTown=="Laingsbug", "NearestTown"] <- "Laingsburg, MI"
Resp[Resp$NearestTown=="West Branch, Michigan", "NearestTown"] <- "West Branch, MI"
Resp[Resp$NearestTown=="Berkley", "NearestTown"] <- "Berkley, MI"
Resp[Resp$NearestTown=="Royal Oak", "NearestTown"] <- "Royal Oak, MI"
Resp[Resp$NearestTown=="Bruce Township/ Romeo Michigan", "NearestTown"] <- "Romeo, MI"
Resp[Resp$NearestTown=="Flat Rock", "NearestTown"] <- "Flat Rock, MI"
Resp[Resp$NearestTown=="Gobles MI", "NearestTown"] <- "Gobles, MI"
Resp[Resp$NearestTown=="houghton", "NearestTown"] <- "Houghton, MI"
Resp[Resp$NearestTown=="houghton mi", "NearestTown"] <- "Houghton, MI"
Resp[Resp$NearestTown=="Howell Mi", "NearestTown"] <- "Howell, MI"
Resp[Resp$NearestTown=="Howell MI ( Marion TWP, Livingston county)", "NearestTown"] <- "Howell, MI"
Resp[Resp$NearestTown=="Howell MI, Livingston county, Marion twp", "NearestTown"] <- "Howell, MI"
Resp[Resp$NearestTown=="Howell MI", "NearestTown"] <- "Howell, MI"
Resp[Resp$NearestTown=="Ishpeming michigan", "NearestTown"] <- "Ishpeming, MI"
Resp[Resp$NearestTown=="Lapeer", "NearestTown"] <- "Lapeer, MI"
Resp[Resp$NearestTown=="Waaren, MI", "NearestTown"] <- "Warren, MI"


# Create a new column for counties
Resp$County <- as.character(Resp$NearestTown)
Resp[Resp$County=="Grandville, MI", "County"] <- "Kent"
Resp[Resp$County=="Alpena, MI", "County"] <- "Alpena"
Resp[Resp$County=="Empire, MI", "County"] <- "Leelanau"
Resp[Resp$County=="New Boston, MI", "County"] <- "Wayne"
Resp[Resp$County=="Escanaba, MI", "County"] <- "Delta"
Resp[Resp$County=="Plainwell, MI", "County"] <- "Allegan"
Resp[Resp$County=="Ortonville, MI", "County"] <- "Oakland"
Resp[Resp$County=="Milford, MI", "County"] <- "Oakland"
Resp[Resp$County=="Allendale, MI", "County"] <- "Ottawa"
Resp[Resp$County=="Au Gres, MI", "County"] <- "Arenac"
Resp[Resp$County=="Berkley, MI", "County"] <- "Oakland"
Resp[Resp$County=="Clio, MI", "County"] <- "Genesee"
Resp[Resp$County=="Detroit, MI", "County"] <- "Wayne"
Resp[Resp$County=="Flat Rock, MI", "County"] <- "Wayne"
Resp[Resp$County=="Georgetown, IN", "County"] <- "Floyd"
Resp[Resp$County=="Gobles, MI", "County"] <- "Van Buren"
Resp[Resp$County=="Grand Rapids, MI", "County"] <- "Kent"
Resp[Resp$County=="Holland, MI", "County"] <- "Ottawa"
Resp[Resp$County=="Houghton, MI", "County"] <- "Houghton"
Resp[Resp$County=="Howell, MI", "County"] <- "Livingston"
Resp[Resp$County=="Ishpeming, MI", "County"] <- "Marquette"
Resp[Resp$County=="Laingsburg, MI", "County"] <- "Shiawassee"
Resp[Resp$County=="Lapeer, MI", "County"] <- "Lapeer"
Resp[Resp$County=="Lewiston, MI", "County"] <- "Montmorency"
Resp[Resp$County=="Livonia, MI", "County"] <- "Wayne"
Resp[Resp$County=="Maple City, MI", "County"] <- "Leelanau"
Resp[Resp$County=="Mt Pleasant, MI", "County"] <- "Isabella"
Resp[Resp$County=="Northville, MI", "County"] <- "Oakland"
Resp[Resp$County=="Novi, MI", "County"] <- "Oakland"
Resp[Resp$County=="Romeo, MI", "County"] <- "Macomb"
Resp[Resp$County=="Redford, MI", "County"] <- "Wayne"
Resp[Resp$County=="Royal Oak, MI", "County"] <- "Oakland"
Resp[Resp$County=="Southfield, MI", "County"] <- "Oakland"
Resp[Resp$County=="Sutton's Bay, MI", "County"] <- "Leelanau"
Resp[Resp$County=="Traverse City, MI", "County"] <- "Grand Traverse"
Resp[Resp$County=="Twin Lake,  MI", "County"] <- "Muskegon"
Resp[Resp$County=="Warren, MI", "County"] <- "Macomb"
Resp[Resp$County=="Wayland, MI", "County"] <- "Allegan"
Resp[Resp$County=="West Branch, MI", "County"] <- "Ogemaw"


# Shorten names of cucurbit types grown
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow , Zucchini , Pattypan , Crookneck , etc.) , Winter Squash (Butternut , Acorn , Buttercup , Delicata , Hubbard , Kabocha , Pumpkin , Spaghetti , etc.) , Melon (Cantaloupe , Honeydew , Muskmelon , Watermelon , etc.)", "VineCropsGrown"] <- "Summer Squash, Winter Squash, Melon"
Resp[Resp$VineCropsGrown=="Cucumber (Salad , Pickling , Slicing , etc.)", "VineCropsGrown"] <- "Cucumber"
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow , Zucchini , Pattypan , Crookneck , etc.) , Winter Squash (Butternut , Acorn , Buttercup , Delicata , Hubbard , Kabocha , Pumpkin , Spaghetti , etc.) , Cucumber (Salad , Pickling , Slicing , etc.)", "VineCropsGrown"] <- "Summer Squash, Winter Squash, Cucumber"
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow , Zucchini , Pattypan , Crookneck , etc.) , Winter Squash (Butternut , Acorn , Buttercup , Delicata , Hubbard , Kabocha , Pumpkin , Spaghetti , etc.) , Cucumber (Salad , Pickling , Slicing , etc.) , Melon (Cantaloupe , Honeydew , Muskmelon , Watermelon , etc.)", "VineCropsGrown"] <- "Summer Squash, Winter Squash, Cucumber, Melon"
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow , Zucchini , Pattypan , Crookneck , etc.)", "VineCropsGrown"] <- "Summer Squash"
Resp[Resp$VineCropsGrown=="Winter Squash (Butternut , Acorn , Buttercup , Delicata , Hubbard , Kabocha , Pumpkin , Spaghetti , etc.) , Cucumber (Salad , Pickling , Slicing , etc.)", "VineCropsGrown"] <- "Winter Squash, Cucumber"
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow , Zucchini , Pattypan , Crookneck , etc.) , Cucumber (Salad , Pickling , Slicing , etc.)", "VineCropsGrown"] <- "Summer Squash, Cucumber"
Resp[Resp$VineCropsGrown=="Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Cucumber (Salad, Pickling, Slicing, etc.)", "VineCropsGrown"] <- "Winter Squash, Cucumber"
Resp[Resp$VineCropsGrown=="Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.) , Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.)", "VineCropsGrown"] <- "Summer Squash, Melon"
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Cucumber (Salad, Pickling, Slicing, etc.) , Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.)","VineCropsGrown"] <- "Summer Squash, Winter Squash, Cucumber, Melon"
Resp[Resp$VineCropsGrown=="Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Cucumber (Salad, Pickling, Slicing, etc.) , Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.)","VineCropsGrown"] <- "Summer Squash, Winter Squash, Cucumber, Melon"
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.)","VineCropsGrown"] <- "Summer Squash, Melon"
Resp[Resp$VineCropsGrown=="Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.) , Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Cucumber (Salad, Pickling, Slicing, etc.)","VineCropsGrown"] <- "Summer Squash, Winter Squash, Cucumber, Melon"
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Cucumber (Salad, Pickling, Slicing, etc.)","VineCropsGrown"] <- "Summer Squash, Cucumber"
Resp[Resp$VineCropsGrown=="Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Cucumber (Salad, Pickling, Slicing, etc.)","VineCropsGrown"] <- "Summer Squash, Winter Squash, Cucumber"
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Cucumber (Salad, Pickling, Slicing, etc.)","VineCropsGrown"] <- "Summer Squash, Winter Squash, Cucumber"
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.)","VineCropsGrown"] <- "Summer Squash, Winter Squash"
Resp[Resp$VineCropsGrown=="Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Cucumber (Salad, Pickling, Slicing, etc.) , Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.)","VineCropsGrown"] <- "Winter Squash, Cucumber, Melon"
Resp[Resp$VineCropsGrown=="Cucumber (Salad, Pickling, Slicing, etc.) , Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.)","VineCropsGrown"] <- "Summer Squash, Winter Squash, Cucumber"
Resp[Resp$VineCropsGrown=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.)","VineCropsGrown"] <- "Summer Squash"
Resp[Resp$VineCropsGrown=="Cucumber (Salad, Pickling, Slicing, etc.) , Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.)","VineCropsGrown"] <- "Cucumber, Melon"
Resp[Resp$VineCropsGrown=="Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.)","VineCropsGrown"] <- "Winter Squash"
Resp[Resp$VineCropsGrown=="Cucumber (Salad, Pickling, Slicing, etc.)","VineCropsGrown"] <- "Cucumber"


# Shorten names of cucurbit types observed
Resp[Resp$VineCropObserved=="Cucumber (Salad, Pickling, Slicing, etc.)", "VineCropObserved"] <- "Cucumber"
Resp[Resp$VineCropObserved=="Cucumber (Salad , Pickling , Slicing , etc.)", "VineCropObserved"] <- "Cucumber"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.)", "VineCropObserved"] <- "Summer Squash"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow , Zucchini , Pattypan , Crookneck , etc.)", "VineCropObserved"] <- "Summer Squash"
Resp[Resp$VineCropObserved=="Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.)", "VineCropObserved"] <- "Winter Squash"
Resp[Resp$VineCropObserved=="Winter Squash (Butternut , Acorn , Buttercup , Delicata , Hubbard , Kabocha , Pumpkin , Spaghetti , etc.)", "VineCropObserved"] <- "Winter Squash"
Resp[Resp$VineCropObserved=="Melon (Cantaloupe , Honeydew , Muskmelon , Watermelon , etc.)", "VineCropObserved"] <- "Melon"
Resp[Resp$VineCropObserved=="Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.)", "VineCropObserved"] <- "Melon"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.)", "VineCropObserved"] <- "Summer Squash, Melon"
Resp[Resp$VineCropObserved=="Cucumber (Salad, Pickling, Slicing, etc.) , Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.)","VineCropObserved"] <- "Summer Squash, Cucumber"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Cucumber (Salad, Pickling, Slicing, etc.)","VineCropObserved"] <- "Summer Squash, Cucumber"
Resp[Resp$VineCropObserved=="Cucumber (Salad , Pickling , Slicing , etc.) , Summer Squash (Yellow , Zucchini , Pattypan , Crookneck , etc.)","VineCropObserved"] <- "Summer Squash, Cucumber"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow , Zucchini , Pattypan , Crookneck , etc.) , Cucumber (Salad , Pickling , Slicing , etc.)","VineCropObserved"] <- "Summer Squash, Cucumber"
Resp[Resp$VineCropObserved=="Winter Squash (Butternut , Acorn , Buttercup , Delicata , Hubbard , Kabocha , Pumpkin , Spaghetti , etc.) , Cucumber (Salad , Pickling , Slicing , etc.)","VineCropObserved"] <- "Winter Squash, Cucumber"
Resp[Resp$VineCropObserved=="Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Cucumber (Salad, Pickling, Slicing, etc.)","VineCropObserved"] <- "Winter Squash, Cucumber"
Resp[Resp$VineCropObserved=="Cucumber (Salad, Pickling, Slicing, etc.) , Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.)","VineCropObserved"] <- "Winter Squash, Cucumber"
Resp[Resp$VineCropObserved=="Cucumber (Salad , Pickling , Slicing , etc.) , Winter Squash (Butternut , Acorn , Buttercup , Delicata , Hubbard , Kabocha , Pumpkin , Spaghetti , etc.)","VineCropObserved"] <- "Winter Squash, Cucumber"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow , Zucchini , Pattypan , Crookneck , etc.) , Winter Squash (Butternut , Acorn , Buttercup , Delicata , Hubbard , Kabocha , Pumpkin , Spaghetti , etc.) , Melon (Cantaloupe , Honeydew , Muskmelon , Watermelon , etc.)","VineCropObserved"] <- "Summer Squash, Winter Squash, Melon"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.)","VineCropObserved"] <- "Summer Squash, Winter Squash"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Cucumber (Salad, Pickling, Slicing, etc.) , Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.)","VineCropObserved"] <- "Summer Squash, Winter Squash, Cucumber, Melon"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Cucumber (Salad, Pickling, Slicing, etc.)","VineCropObserved"] <- "Summer Squash, Winter Squash, Cucumber"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Cucumber (Salad, Pickling, Slicing, etc.) , Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.)","VineCropObserved"] <- "Summer Squash, Cucumber, Melon"
Resp[Resp$VineCropObserved=="Summer Squash (Yellow, Zucchini, Pattypan, Crookneck, etc.) , Winter Squash (Butternut, Acorn, Buttercup, Delicata, Hubbard, Kabocha, Pumpkin, Spaghetti, etc.) , Melon (Cantaloupe, Honeydew, Muskmelon, Watermelon, etc.)","VineCropObserved"] <- "Summer Squash, Winter Squash, Melon"


# Create a new column indicating whether or not summer and/or winter squash were at least one of the vine
# crops observed during surveys. 
Resp$SumWin <- Resp$VineCropObserved
Resp[Resp$SumWin=="Summer Squash","SumWin"] <- "Summer"
Resp[Resp$SumWin=="Winter Squash","SumWin"] <- "Winter"
Resp[Resp$SumWin=="Summer Squash, Cucumber","SumWin"] <- "Mixed"
Resp[Resp$SumWin=="Summer Squash, Melon","SumWin"] <- "Mixed"
Resp[Resp$SumWin=="Summer Squash, Winter Squash","SumWin"] <- "Summer & Winter"
Resp[Resp$SumWin=="Summer Squash, Winter Squash, Cucumber","SumWin"] <- "Mixed"
Resp[Resp$SumWin=="Summer Squash, Winter Squash, Cucumber, Melon","SumWin"] <- "Mixed"
Resp[Resp$SumWin=="Winter Squash, Cucumber","SumWin"] <- "Mixed"
Resp[Resp$SumWin=="Cucumber","SumWin"] <- "Cucumber"
Resp[Resp$SumWin=="Melon","SumWin"] <- "Melon"
Resp[Resp$SumWin=="Summer Squash, Cucumber, Melon","SumWin"] <- "Mixed"
Resp[Resp$SumWin=="Summer Squash, Winter Squash, Melon","SumWin"] <- "Mixed"


# Shorten names for irrigation categories
Resp[Resp$Irrigation=="Overhead/Sprinkler/By hand", "Irrigation"] <- "Overhead"
Resp[Resp$Irrigation=="Overhead/Sprinkler/By hand , Trickle/Drip", "Irrigation"] <- "Overhead + Drip"
Resp[Resp$Irrigation=="Trickle/Drip , Overhead/Sprinkler/By hand", "Irrigation"] <- "Overhead + Drip"
Resp[Resp$Irrigation=="Trickle/Drip", "Irrigation"] <- "Drip"
Resp[Resp$Irrigation=="Overhead/Sprinkler/By hand , None", "Irrigation"] <- "Overhead"


# Group/shorten names for insecticide categories
Resp[Resp$Insecticides=="Approved for organic use, Biopesticides", "Insecticides"] <- "Organic"
Resp[Resp$Insecticides=="Approved for organic use , Biopesticides", "Insecticides"] <- "Organic"
Resp[Resp$Insecticides=="Approved for organic use , Biopesticides , Restricted use or conventional synthetic pesticides", "Insecticides"] <- "Organic + Conventional"
Resp[Resp$Insecticides=="Approved for organic use, Biopesticides , Restricted use or conventional synthetic pesticides", "Insecticides"] <- "Organic + Conventional"
Resp[Resp$Insecticides=="None", "Insecticides"] <- "None"
Resp[Resp$Insecticides=="Restricted use or conventional synthetic pesticides", "Insecticides"] <- "Conventional"
Resp[Resp$Insecticides=="Restricted use or conventional synthetic pesticides , Approved for organic use, Biopesticides", "Insecticides"] <- "Organic + Conventional"


# Sort the data set by date; pdf surveys that were manually entered into the spreadsheet
# will appear out of order otherwise. 
Resp$Date <- as.character(Resp$Date)
as.Date(Resp$Date)
Resp <- arrange(Resp, Date)


# Make a new column to group by month
Resp$Month <- Resp$Date

# Change dates in month column to months
Resp[1:2, "Month"] <- "June"
Resp[144:153, "Month"] <- "June"
Resp[3:21, "Month"] <- "July"
Resp[60:99, "Month"] <- "July"
Resp[154:204, "Month"] <- "July"
Resp[22:57, "Month"] <- "August"
Resp[100:143, "Month"] <- "August"
Resp[205:286, "Month"] <- "August"
Resp[58, "Month"] <- "September"
Resp[287:291, "Month"] <- "September"
Resp[59, "Month"] <- "October"

# Create a year column
Resp$Year <- Resp$Date
Resp[1:59, "Year"] <- "2017"
Resp[60:143, "Year"] <- "2018"
Resp[144:291, "Year"] <- "2019"


# Remove surveys with incorrect identifications
Resp <- Resp[Resp$PhotographAccurate!="No",]


# Remove surveys from October
Resp <- Resp[Resp$Month!="October",]


# Determine numeric outliers for each bee group (bee counts greater than the third quartile
# plus 1.5 times the interquartile range of each bee group)

# Squash Bees

# Mean number of squash bees overall
mean(Resp$SquashBees)

# Below are the squash bee quartiles, median, interquartile range, and outlier threshold
# for squash bees
quantile(Resp$SquashBees)

Q1.sb <- 0
median.sb <- 2
Q3.sb <- 4
IQR.sb <- Q3.sb - Q1.sb
Outliers.sb = Q3.sb + 1.5*IQR.sb
Outliers.sb

# Any squash bee count above 10 is an outlier; create a dataframe for squash bees excluding 
# outliers. 
Resp.SB <- Resp[Resp$SquashBees<=10,]


# Honey Bees

# Mean number of honey bees overall
mean(Resp$HoneyBees)

# Below are the honey bee quartiles, median, interquartile range, and outlier threshold
# for squash bees
quantile(Resp$HoneyBees)

Q1.hb <- 0
median.hb <- 0
Q3.hb <- 1
IQR.hb <- Q3.hb - Q1.hb
Outliers.hb = Q3.hb + 1.5*IQR.hb
Outliers.hb

# Any honey bee count above 2.5 is an outlier; create a dataframe for honey bees excluding 
# outliers.
Resp.HB <- Resp[Resp$HoneyBees<=2.5,]


# Bumble Bees

# Mean number of bumble bees overall
mean(Resp$BumbleBees)

# Below are the bumble bee quartiles, median, interquartile range, and outlier threshold
# for squash bees
quantile(Resp$BumbleBees)

Q1.bb <- 0
median.bb <- 0
Q3.bb <- 1
IQR.bb <- Q3.bb - Q1.bb
Outliers.bb = Q3.bb + 1.5*IQR.bb
Outliers.bb

# Any bumble bee count above 2.5 is an outlier; create a dataframe for bumble bees excluding 
# outliers.
Resp.BB <- Resp[Resp$BumbleBees<=2.5,]


# Other Bees

# Mean number of other bees overall
mean(Resp$OtherBees)

# Below are the other bee quartiles, median, interquartile range, and outlier threshold
# for squash bees
quantile(Resp$OtherBees)

Q1.ob <- 0
median.ob <- 1
Q3.ob <- 2
IQR.ob <- Q3.ob - Q1.ob
Outliers.ob = Q3.ob + 1.5*IQR.ob
Outliers.ob

# Any other bee count above 5 is an outlier; create a dataframe for other bees excluding 
# outliers.
Resp.OB <- Resp[Resp$OtherBees<=5,]


# Exclude the following categories within factors due to low sample sizes (lower than 5 for at least
# of the bee groups)

# Exclude 25-50 tillage depth
Resp.SB <- Resp.SB[Resp.SB$TillageDepth!="25-50",]
Resp.HB <- Resp.HB[Resp.HB$TillageDepth!="25-50",]
Resp.BB <- Resp.BB[Resp.BB$TillageDepth!="25-50",]
Resp.OB <- Resp.OB[Resp.OB$TillageDepth!="25-50",]


# Exclude None and Trickle/Drip , None irrigation
Resp.SB <- Resp.SB[Resp.SB$Irrigation!="None",]
Resp.SB <- Resp.SB[Resp.SB$Irrigation!="Trickle/Drip , None",]

Resp.HB <- Resp.HB[Resp.HB$Irrigation!="None",]
Resp.HB <- Resp.HB[Resp.HB$Irrigation!="Trickle/Drip , None",]

Resp.BB <- Resp.BB[Resp.BB$Irrigation!="None",]
Resp.BB <- Resp.BB[Resp.BB$Irrigation!="Trickle/Drip , None",]

Resp.OB <- Resp.OB[Resp.OB$Irrigation!="None",]
Resp.OB <- Resp.OB[Resp.OB$Irrigation!="Trickle/Drip , None",]


# Remove Organic + Conventional insecticides
Resp.SB <- Resp.SB[Resp.SB$Insecticides!="Organic + Conventional",]
Resp.HB <- Resp.HB[Resp.HB$Insecticides!="Organic + Conventional",]
Resp.BB <- Resp.BB[Resp.BB$Insecticides!="Organic + Conventional",]
Resp.OB <- Resp.OB[Resp.OB$Insecticides!="Organic + Conventional",]


# Exclude Cucumber and Melon in vine crop observed (when the only type observed) 
# since squash bees do not visit these plants; sample sizes also low
Resp.SB <- Resp.SB[Resp.SB$SumWin!="Cucumber",]
Resp.SB <- Resp.SB[Resp.SB$SumWin!="Melon",]

Resp.HB <- Resp.HB[Resp.HB$SumWin!="Cucumber",]
Resp.HB <- Resp.HB[Resp.HB$SumWin!="Melon",]

Resp.BB <- Resp.BB[Resp.BB$SumWin!="Cucumber",]
Resp.BB <- Resp.BB[Resp.BB$SumWin!="Melon",]

Resp.OB <- Resp.OB[Resp.OB$SumWin!="Cucumber",]
Resp.OB <- Resp.OB[Resp.OB$SumWin!="Melon",]


# View the dataframes
View(Resp.SB)
View(Resp.HB)
View(Resp.BB)
View(Resp.OB)






#######################################################################################
###### Squash Bee Analysis and Figures ######

# Note:Glmm's were are used only for squash bees. Due to consistent model convergence issues 
# Kruskal-Wallis tests were used for the other 3 bee groups. 
# Not all figures seen below were presented in the submitted manuscript. 

# Tillage Type
SB.tt <- ddply(Resp.SB, "TillageType", summarise,
               N = length(SquashBees),
               mean = mean(SquashBees),
               median = median(SquashBees),
               sum = sum(SquashBees),
               sd = sd(SquashBees), 
               se = sd/sqrt(N)
               
)               
SB.tt


# Bar plot
sb.tt<-ggplot(SB.tt, aes(x=TillageType, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("No Tillage","Reduced Tillage","Full Tillage")) +
  xlab(NULL)
sb.tt

sb.tt<-sb.tt + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
sb.tt<-sb.tt + scale_y_continuous("Mean ± SEM squash bees per survey", expand=c(0,0), limits = c(0,3.5))
sb.tt


# Tillage Depth
SB.td <- ddply(Resp.SB, "TillageDepth", summarise,
               N = length(SquashBees),
               mean = mean(SquashBees),
               median = median(SquashBees),
               sum = sum(SquashBees),
               sd = sd(SquashBees), 
               se = sd/sqrt(N)
               
)               
SB.td


sb.td<-ggplot(SB.td, aes(x=TillageDepth, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("0", "3-14", "15-25")) +
  xlab(NULL)
sb.td

sb.td<-sb.td + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
sb.td<-sb.td + scale_y_continuous("Mean ± SEM squash bees per survey", expand=c(0,0), limits = c(0,4))
sb.td


# Mulch
SB.m <- ddply(Resp.SB, "Mulch", summarise,
               N = length(SquashBees),
               mean = mean(SquashBees),
               median = median(SquashBees),
               sum = sum(SquashBees),
               sd = sd(SquashBees), 
               se = sd/sqrt(N)
               
)               
SB.m


sb.m<-ggplot(SB.m, aes(x=Mulch, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey30", "grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("None", "Plant Material", "Plastic", "Plastic + Plant Material")) +
  xlab(NULL)
sb.m

sb.m<-sb.m + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
sb.m<-sb.m + scale_y_continuous("Mean ± SEM squash bees per survey", expand=c(0,0), limits = c(0,6.5))
sb.m


# Insecticides
SB.in <- ddply(Resp.SB, "Insecticides", summarise,
              N = length(SquashBees),
              mean = mean(SquashBees),
              median = median(SquashBees),
              sum = sum(SquashBees),
              sd = sd(SquashBees), 
              se = sd/sqrt(N)
              
)               
SB.in


sb.in<-ggplot(SB.in, aes(x=Insecticides, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("None", "Organic", "Conventional")) +
  xlab(NULL)
sb.in

sb.in<-sb.in + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
sb.in<-sb.in + scale_y_continuous("Mean ± SEM squash bees per survey", expand=c(0,0), limits = c(0,4))
sb.in


# Irrigation
SB.ir <- ddply(Resp.SB, "Irrigation", summarise,
               N = length(SquashBees),
               mean = mean(SquashBees),
               median = median(SquashBees),
               sum = sum(SquashBees),
               sd = sd(SquashBees), 
               se = sd/sqrt(N)
               
)               
SB.ir


sb.ir<-ggplot(SB.ir, aes(x=Irrigation, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("Drip", "Overhead", "Overhead + Drip")) +
  xlab(NULL)
sb.ir

sb.ir<-sb.ir + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
sb.ir<-sb.ir + scale_y_continuous("Mean ± SEM squash bees per survey", expand=c(0,0), limits = c(0,8))
sb.ir


# Cucurbit Area
SB.ca <- ddply(Resp.SB, "CucurbitArea", summarise,
               N = length(SquashBees),
               mean = mean(SquashBees),
               median = median(SquashBees),
               sum = sum(SquashBees),
               sd = sd(SquashBees), 
               se = sd/sqrt(N)
               
)               
SB.ca


sb.ca<-ggplot(SB.ca, aes(x=CucurbitArea, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("<1 acre", ">1 acre")) +
  xlab(NULL)
sb.ca

sb.ca<-sb.ca + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
sb.ca<-sb.ca + scale_y_continuous("Mean ± SEM squash bees per survey", expand=c(0,0), limits = c(0,8))
sb.ca


# Vine crop observed
SB.sw <- ddply(Resp.SB, "SumWin", summarise,
               N = length(SquashBees),
               mean = mean(SquashBees),
               median = median(SquashBees),
               sum = sum(SquashBees),
               sd = sd(SquashBees), 
               se = sd/sqrt(N)
               
)               
SB.sw


sb.sw<-ggplot(SB.sw, aes(x=SumWin, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey30","grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("Mixed", "Summer", "Summer & Winter", "Winter")) +
  xlab(NULL)
sb.sw

sb.sw<-sb.sw + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
sb.sw<-sb.sw + scale_y_continuous("Mean ± SEM squash bees per survey", expand=c(0,0), limits = c(0,4))
sb.sw


### Squash bee models


# Model Factors
date.sb <- as.factor(Resp.SB$Date)
year.sb <- as.factor(Resp.SB$Year)
mulch.sb <- as.factor(Resp.SB$Mulch)
tilltype.sb <- as.factor(Resp.SB$TillageType)
tilldep.sb <- as.factor(Resp.SB$TillageDepth)
insecticide.sb <- as.factor(Resp.SB$Insecticides)
cuc.area.sb <- as.factor(Resp.SB$CucurbitArea)
vc.obs.sb <- as.factor(Resp.SB$VineCropObserved)
irrigation.sb <- as.factor(Resp.SB$Irrigation)
sum.win.sb <- as.factor(Resp.SB$SumWin)
county.sb <- as.factor(Resp.SB$County)
year.sb <- as.factor(Resp.SB$Year)

# Models

msb.0 <- glmmadmb (SquashBees ~ 1 + (1|date.sb/county.sb), data = Resp.SB, family="nbinom2")

msb.1 <- glmmadmb (SquashBees ~ tilltype.sb + (1|date.sb/county.sb), data = Resp.SB, family="nbinom2")

msb.2 <- glmmadmb (SquashBees ~ tilldep.sb + (1|date.sb/county.sb), data = Resp.SB, family="nbinom2")

msb.3 <- glmmadmb (SquashBees ~ mulch.sb + (1|date.sb/county.sb), data = Resp.SB, family="nbinom2")

msb.4 <- glmmadmb (SquashBees ~ insecticide.sb + (1|date.sb/county.sb), data = Resp.SB, family="nbinom2")

msb.5 <- glmmadmb (SquashBees ~ irrigation.sb + (1|date.sb/county.sb), data = Resp.SB, family="nbinom2")

msb.6 <- glmmadmb (SquashBees ~ cuc.area.sb + (1|date.sb/county.sb), data = Resp.SB, family="nbinom2")

msb.7 <- glmmadmb (SquashBees ~ sum.win.sb + (1|date.sb/county.sb), data = Resp.SB, family="nbinom2")

# Use AIC to find the best model:
AICctab(msb.0, msb.1, msb.2, msb.3, msb.4, msb.5, msb.6, msb.7, weights=TRUE)
# msb.1 is the only model with and AIC weight greater than the null model
summary(msb.1)
summary(msb.2)
summary(msb.3)
summary(msb.4)
summary(msb.5)
summary(msb.6)
summary(msb.7)


anova(msb.0, msb.1)
# tillage type significant
anova(msb.0, msb.2)
# tillage depth not significant
anova(msb.0, msb.3)
# mulch not significant
anova(msb.0, msb.4)
# incsecticides not significant
anova(msb.0, msb.5)
# irrigation not significant
anova(msb.0, msb.6)
# cucurbit area significant
anova(msb.0, msb.7)
# vine crop observed not significant


tt.emm.sb <- emmeans(msb.1, pairwise ~ tilltype.sb, adjust="fdr")
tt.emm.sb
# Significant difference between Full Tillage and No Tillage,
# and between Full Tillage and Reduced Tillage. 





#######################################################################################
###### Honey Bee Analysis and Figures ######

# Tillage Type
HB.tt <- ddply(Resp.HB, "TillageType", summarise,
               N = length(HoneyBees),
               mean = mean(HoneyBees),
               median = median(HoneyBees),
               sum = sum(HoneyBees),
               sd = sd(HoneyBees), 
               se = sd/sqrt(N)
               
)               
HB.tt


# Bar plot
hb.tt<-ggplot(HB.tt, aes(x=TillageType, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("No Tillage","Reduced Tillage","Full Tillage")) +
  xlab(NULL)
hb.tt

hb.tt<-hb.tt + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
hb.tt<-hb.tt + scale_y_continuous("Mean ± SEM honey bees per survey", expand=c(0,0), limits = c(0,0.75))
hb.tt


# Tillage Depth
HB.td <- ddply(Resp.HB, "TillageDepth", summarise,
               N = length(HoneyBees),
               mean = mean(HoneyBees),
               median = median(HoneyBees),
               sum = sum(HoneyBees),
               sd = sd(HoneyBees), 
               se = sd/sqrt(N)
               
)               
HB.td


hb.td<-ggplot(HB.td, aes(x=TillageDepth, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("0", "3-14", "15-25")) +
  xlab(NULL)
hb.td

hb.td<-hb.td + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
hb.td<-hb.td + scale_y_continuous("Mean ± SEM honey bees per survey", expand=c(0,0), limits = c(0,0.75))
hb.td


# Mulch
HB.m <- ddply(Resp.HB, "Mulch", summarise,
              N = length(HoneyBees),
              mean = mean(HoneyBees),
              median = median(HoneyBees),
              sum = sum(HoneyBees),
              sd = sd(HoneyBees), 
              se = sd/sqrt(N)
              
)               
HB.m


hb.m<-ggplot(HB.m, aes(x=Mulch, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey30", "grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("None", "Plant Material", "Plastic", "Plastic + Plant Material")) +
  xlab(NULL)
hb.m

hb.m<-hb.m + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
hb.m<-hb.m + scale_y_continuous("Mean ± SEM honey bees per survey", expand=c(0,0), limits = c(0,1))
hb.m


# Insecticides
HB.in <- ddply(Resp.HB, "Insecticides", summarise,
               N = length(HoneyBees),
               mean = mean(HoneyBees),
               median = median(HoneyBees),
               sum = sum(HoneyBees),
               sd = sd(HoneyBees), 
               se = sd/sqrt(N)
               
)               
HB.in


hb.in<-ggplot(HB.in, aes(x=Insecticides, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("None", "Organic", "Conventional")) +
  xlab(NULL)
hb.in

hb.in<-hb.in + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
hb.in<-hb.in + scale_y_continuous("Mean ± SEM honey bees per survey", expand=c(0,0), limits = c(0,0.75))
hb.in


# Irrigation
HB.ir <- ddply(Resp.HB, "Irrigation", summarise,
               N = length(HoneyBees),
               mean = mean(HoneyBees),
               median = median(HoneyBees),
               sum = sum(HoneyBees),
               sd = sd(HoneyBees), 
               se = sd/sqrt(N)
               
)               
HB.ir


hb.ir<-ggplot(HB.ir, aes(x=Irrigation, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("Drip", "Overhead", "Overhead + Drip")) +
  xlab(NULL)
hb.ir

hb.ir<-hb.ir + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
hb.ir<-hb.ir + scale_y_continuous("Mean ± SEM honey bees per survey", expand=c(0,0), limits = c(0,0.75))
hb.ir


# Cucurbit Area
HB.ca <- ddply(Resp.HB, "CucurbitArea", summarise,
               N = length(HoneyBees),
               mean = mean(HoneyBees),
               median = median(HoneyBees),
               sum = sum(HoneyBees),
               sd = sd(HoneyBees), 
               se = sd/sqrt(N)
               
)               
HB.ca


hb.ca<-ggplot(HB.ca, aes(x=CucurbitArea, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("<1 acre", ">1 acre")) +
  xlab(NULL)
hb.ca

hb.ca<-hb.ca + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
hb.ca<-hb.ca + scale_y_continuous("Mean ± SEM honey bees per survey", expand=c(0,0), limits = c(0,1))
hb.ca


# Vine crop observed
HB.sw <- ddply(Resp.HB, "SumWin", summarise,
               N = length(HoneyBees),
               mean = mean(HoneyBees),
               median = median(HoneyBees),
               sum = sum(HoneyBees),
               sd = sd(HoneyBees), 
               se = sd/sqrt(N)
               
)               
HB.sw


hb.sw<-ggplot(HB.sw, aes(x=SumWin, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey30","grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("Mixed", "Summer", "Summer & Winter", "Winter")) +
  xlab(NULL)
hb.sw

hb.sw<-hb.sw + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
hb.sw<-hb.sw + scale_y_continuous("Mean ± SEM honey bees per survey", expand=c(0,0), limits = c(0,0.75))
hb.sw


### Honey bee Kruskal-Wallis tests

# Tillage Type
kruskal.test(HoneyBees ~ as.factor(TillageType), data=Resp.HB)
# Not significant


# Tillage Depth
kruskal.test(HoneyBees ~ as.factor(TillageDepth), data=Resp.HB)
# Not significant


# Mulch
kruskal.test(HoneyBees ~ as.factor(Mulch), data=Resp.HB)
# Significant; use Dunn's test to examine pairwise differences
dunn.test(Resp.HB$HoneyBees, Resp.HB$Mulch, altp=T, method="holm")
# Significant difference between None and Plant Material categories


# Insecticides
kruskal.test(HoneyBees ~ as.factor(Insecticides), data=Resp.HB)
# Not significant


# Irrigation
kruskal.test(HoneyBees ~ as.factor(Irrigation), data=Resp.HB)
# Not significant


# Vine Crop Observed
kruskal.test(HoneyBees ~ as.factor(SumWin), data=Resp.HB)
# Not significant





#######################################################################################
###### Bumble Bee Analysis and Figures ######


# Tillage Type
BB.tt <- ddply(Resp.BB, "TillageType", summarise,
               N = length(BumbleBees),
               mean = mean(BumbleBees),
               median = median(BumbleBees),
               sum = sum(BumbleBees),
               sd = sd(BumbleBees), 
               se = sd/sqrt(N)
               
)               
BB.tt


# Bar plot
bb.tt<-ggplot(BB.tt, aes(x=TillageType, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("No Tillage","Reduced Tillage","Full Tillage")) +
  xlab(NULL)
bb.tt

bb.tt<-bb.tt + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
bb.tt<-bb.tt + scale_y_continuous("Mean ± SEM bumble bees per survey", expand=c(0,0), limits = c(0,1))
bb.tt


# Tillage Depth
BB.td <- ddply(Resp.BB, "TillageDepth", summarise,
               N = length(BumbleBees),
               mean = mean(BumbleBees),
               median = median(BumbleBees),
               sum = sum(BumbleBees),
               sd = sd(BumbleBees), 
               se = sd/sqrt(N)
               
)               
BB.td


bb.td<-ggplot(BB.td, aes(x=TillageDepth, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("0", "3-14", "15-25")) +
  xlab(NULL)
bb.td

bb.td<-bb.td + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
bb.td<-bb.td + scale_y_continuous("Mean ± SEM bumble bees per survey", expand=c(0,0), limits = c(0,1))
bb.td


# Mulch
BB.m <- ddply(Resp.BB, "Mulch", summarise,
              N = length(BumbleBees),
              mean = mean(BumbleBees),
              median = median(BumbleBees),
              sum = sum(BumbleBees),
              sd = sd(BumbleBees), 
              se = sd/sqrt(N)
              
)               
BB.m


bb.m<-ggplot(BB.m, aes(x=Mulch, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey30", "grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("None", "Plant Material", "Plastic", "Plastic + Plant Material")) +
  xlab(NULL)
bb.m

bb.m<-bb.m + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
bb.m<-bb.m + scale_y_continuous("Mean ± SEM bumble bees per survey", expand=c(0,0), limits = c(0,1))
bb.m


# Insecticides
BB.in <- ddply(Resp.BB, "Insecticides", summarise,
               N = length(BumbleBees),
               mean = mean(BumbleBees),
               median = median(BumbleBees),
               sum = sum(BumbleBees),
               sd = sd(BumbleBees), 
               se = sd/sqrt(N)
               
)               
BB.in


bb.in<-ggplot(BB.in, aes(x=Insecticides, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("None", "Organic", "Conventional")) +
  xlab(NULL)
bb.in

bb.in<-bb.in + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
bb.in<-bb.in + scale_y_continuous("Mean ± SEM bumble bees per survey", expand=c(0,0), limits = c(0,1))
bb.in


# Irrigation
BB.ir <- ddply(Resp.BB, "Irrigation", summarise,
               N = length(BumbleBees),
               mean = mean(BumbleBees),
               median = median(BumbleBees),
               sum = sum(BumbleBees),
               sd = sd(BumbleBees), 
               se = sd/sqrt(N)
               
)               
BB.ir


bb.ir<-ggplot(BB.ir, aes(x=Irrigation, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("Drip", "Overhead", "Overhead + Drip")) +
  xlab(NULL)
bb.ir

bb.ir<-bb.ir + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
bb.ir<-bb.ir + scale_y_continuous("Mean ± SEM bumble bees per survey", expand=c(0,0), limits = c(0,0.75))
bb.ir


# Cucurbit Area
BB.ca <- ddply(Resp.BB, "CucurbitArea", summarise,
               N = length(BumbleBees),
               mean = mean(BumbleBees),
               median = median(BumbleBees),
               sum = sum(BumbleBees),
               sd = sd(BumbleBees), 
               se = sd/sqrt(N)
               
)               
BB.ca


bb.ca<-ggplot(BB.ca, aes(x=CucurbitArea, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("<1 acre", ">1 acre")) +
  xlab(NULL)
bb.ca

bb.ca<-bb.ca + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
bb.ca<-bb.ca + scale_y_continuous("Mean ± SEM bumble bees per survey", expand=c(0,0), limits = c(0,1))
bb.ca


# Vine crop observed
BB.sw <- ddply(Resp.BB, "SumWin", summarise,
               N = length(BumbleBees),
               mean = mean(BumbleBees),
               median = median(BumbleBees),
               sum = sum(BumbleBees),
               sd = sd(BumbleBees), 
               se = sd/sqrt(N)
               
)               
BB.sw


bb.sw<-ggplot(BB.sw, aes(x=SumWin, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey30","grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("Mixed", "Summer", "Summer & Winter", "Winter")) +
  xlab(NULL)
bb.sw

bb.sw<-bb.sw + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
bb.sw<-bb.sw + scale_y_continuous("Mean ± SEM bumble bees per survey", expand=c(0,0), limits = c(0,1.25))
bb.sw


### Bumble bee Kruskal-Wallis tests

# Tillage Type
kruskal.test(BumbleBees ~ as.factor(TillageType), data=Resp.BB)
# Not significant


# Tillage Depth
kruskal.test(BumbleBees ~ as.factor(TillageDepth), data=Resp.BB)
# Not significant


# Mulch
kruskal.test(BumbleBees ~ as.factor(Mulch), data=Resp.BB)
# Not significant


# Insecticides
kruskal.test(BumbleBees ~ as.factor(Insecticides), data=Resp.BB)
# Not significant


# Irrigation
kruskal.test(BumbleBees ~ as.factor(Irrigation), data=Resp.BB)
# Not significant


# Vine Crop Observed
kruskal.test(BumbleBees ~ as.factor(SumWin), data=Resp.BB)
# Significant; use Dunn's test to examine pairwise differences. 
dunn.test(Resp.BB$BumbleBees, Resp.BB$SumWin, altp=T, method="holm")
# Significant difference between Mixes and Winter, and between Summer and Winter





#######################################################################################
###### Other Bee Analysis and Figures ######


# Tillage Type
OB.tt <- ddply(Resp.OB, "TillageType", summarise,
               N = length(OtherBees),
               mean = mean(OtherBees),
               median = median(OtherBees),
               sum = sum(OtherBees),
               sd = sd(OtherBees), 
               se = sd/sqrt(N)
               
)               
OB.tt


# Bar plot
ob.tt<-ggplot(OB.tt, aes(x=TillageType, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("No Tillage","Reduced Tillage","Full Tillage")) +
  xlab(NULL)
ob.tt

ob.tt<-ob.tt + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
ob.tt<-ob.tt + scale_y_continuous("Mean ± SEM other bees per survey", expand=c(0,0), limits = c(0,2))
ob.tt


# Tillage Depth
OB.td <- ddply(Resp.OB, "TillageDepth", summarise,
               N = length(OtherBees),
               mean = mean(OtherBees),
               median = median(OtherBees),
               sum = sum(OtherBees),
               sd = sd(OtherBees), 
               se = sd/sqrt(N)
               
)               
OB.td


ob.td<-ggplot(OB.td, aes(x=TillageDepth, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("0", "3-14", "15-25")) +
  xlab(NULL)
ob.td

ob.td<-ob.td + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
ob.td<-ob.td + scale_y_continuous("Mean ± SEM other bees per survey", expand=c(0,0), limits = c(0,1.5))
ob.td


# Mulch
OB.m <- ddply(Resp.OB, "Mulch", summarise,
              N = length(OtherBees),
              mean = mean(OtherBees),
              median = median(OtherBees),
              sum = sum(OtherBees),
              sd = sd(OtherBees), 
              se = sd/sqrt(N)
              
)               
OB.m


ob.m<-ggplot(OB.m, aes(x=Mulch, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey30", "grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("None", "Plant Material", "Plastic", "Plastic + Plant Material")) +
  xlab(NULL)
ob.m

ob.m<-ob.m + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
ob.m<-ob.m + scale_y_continuous("Mean ± SEM other bees per survey", expand=c(0,0), limits = c(0,3))
ob.m


# Insecticides
OB.in <- ddply(Resp.OB, "Insecticides", summarise,
               N = length(OtherBees),
               mean = mean(OtherBees),
               median = median(OtherBees),
               sum = sum(OtherBees),
               sd = sd(OtherBees), 
               se = sd/sqrt(N)
               
)               
OB.in


ob.in<-ggplot(OB.in, aes(x=Insecticides, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("None", "Organic", "Conventional")) +
  xlab(NULL)
ob.in

ob.in<-ob.in + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
ob.in<-ob.in + scale_y_continuous("Mean ± SEM other bees per survey", expand=c(0,0), limits = c(0,1.5))
ob.in


# Irrigation
OB.ir <- ddply(Resp.OB, "Irrigation", summarise,
               N = length(OtherBees),
               mean = mean(OtherBees),
               median = median(OtherBees),
               sum = sum(OtherBees),
               sd = sd(OtherBees), 
               se = sd/sqrt(N)
               
)               
OB.ir


ob.ir<-ggplot(OB.ir, aes(x=Irrigation, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("Drip", "Overhead", "Overhead + Drip")) +
  xlab(NULL)
ob.ir

ob.ir<-ob.ir + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
ob.ir<-ob.ir + scale_y_continuous("Mean ± SEM other bees per survey", expand=c(0,0), limits = c(0,1.5))
ob.ir


# Cucurbit Area
OB.ca <- ddply(Resp.OB, "CucurbitArea", summarise,
               N = length(OtherBees),
               mean = mean(OtherBees),
               median = median(OtherBees),
               sum = sum(OtherBees),
               sd = sd(OtherBees), 
               se = sd/sqrt(N)
               
)               
OB.ca


ob.ca<-ggplot(OB.ca, aes(x=CucurbitArea, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey50","grey80")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("<1 acre", ">1 acre")) +
  xlab(NULL)
ob.ca

ob.ca<-ob.ca + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
ob.ca<-ob.ca + scale_y_continuous("Mean ± SEM other bees per survey", expand=c(0,0), limits = c(0,1.25))
ob.ca


# Vine crop observed
OB.sw <- ddply(Resp.OB, "SumWin", summarise,
               N = length(OtherBees),
               mean = mean(OtherBees),
               median = median(OtherBees),
               sum = sum(OtherBees),
               sd = sd(OtherBees), 
               se = sd/sqrt(N)
               
)               
OB.sw


ob.sw<-ggplot(OB.sw, aes(x=SumWin, y=mean)) + 
  geom_bar(color="black", stat="identity", fill=c("grey30","grey50","grey80","grey100")) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=.2,                 
                position=position_dodge(.9)) + theme(panel.background = element_rect(fill="white")) +
  scale_x_discrete(limits=c("Mixed", "Summer", "Summer & Winter", "Winter")) +
  xlab(NULL)
ob.sw

ob.sw<-ob.sw + theme(axis.line = element_line(size=0.8, color = "black", linetype = "solid"), axis.text = element_text(color="black", size = 14), axis.ticks = element_line(color="black", size = 1))
ob.sw<-ob.sw + scale_y_continuous("Mean ± SEM other bees per survey", expand=c(0,0), limits = c(0,1.5))
ob.sw


### Bumble bee Kruskal-Wallis tests

# Tillage Type
kruskal.test(OtherBees ~ as.factor(TillageType), data=Resp.OB)
# Significant; use Dunn's test to examine pairwise differences. 
dunn.test(Resp.OB$OtherBees, Resp.OB$TillageType, altp=T, method="holm")
# Significant difference between No Tillage and Reduced Tillage


# Tillage Depth
kruskal.test(OtherBees ~ as.factor(TillageDepth), data=Resp.OB)
# Significant; use Dunn's test to examine pairwise differences. 
dunn.test(Resp.OB$OtherBees, Resp.OB$TillageDepth, altp=T, method="holm")
# Significant difference between 0 and 3-14 cm


# Mulch
kruskal.test(OtherBees ~ as.factor(Mulch), data=Resp.OB)
# Not significant


# Insecticides
kruskal.test(OtherBees ~ as.factor(Insecticides), data=Resp.OB)
# Not significant


# Irrigation
kruskal.test(OtherBees ~ as.factor(Irrigation), data=Resp.OB)
# Not significant


# Vine Crop Observed
kruskal.test(OtherBees ~ as.factor(SumWin), data=Resp.OB)
# Not significant






