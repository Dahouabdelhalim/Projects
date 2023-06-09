library(dplyr)
library(readxl)
library(tidyr)
library(stringr)

'''
Note! Three input files are needed. Please change the file names accordingly. 

1) The file from deduplication ("images_matches.csv") - line 18
2) The file from Pykogition or annotation - line 49
3) The file from the Facebook Ad Library -79
'''

### The first step is to load in the images_matches.csv from the 'Images as Data' Deduplication Script 
### See: https://github.com/norawebbwilliams/images_as_data/tree/master/code

IAD_data <- read.csv("images_matches.csv")
IAD_data <- as_tibble(IAD_data)

### We need to wrangle the IAD_data and prep it for joining with the Pykogniton data
### Add a group_id for each group of unique images in IAD_Data

IAD_data$group_id <- seq.int(nrow(IAD_data))

## Create a new column that combines all images within each group

IAD_data$matches_combined <-  paste(IAD_data$id, IAD_data$matches, sep=", ")

## Unnest the groups, so that each image within the group is its own row but labelled with its group_id

IAD_data <- as_tibble(IAD_data) %>% 
  select(matches_combined, group_id) %>% 
  mutate(matches_combined = strsplit(as.character(matches_combined), ",")) %>% 
  unnest(matches_combined) %>% 
  filter(!grepl("\\\\[]", matches_combined)) 

## Clean up leftover whitespace

IAD_data$matches_combined <- str_replace_all(IAD_data$matches_combined, " ", "")

## To prep IAD_data for hydration with the Rekognition API data, rename matches_combined to imageName

names(IAD_data)[1] <- "imageName"

## Now, load in the output from Pykognition with or without your annotations
## We name this data AWS_data for "Amazon Web Service"

AWS_data <- read.csv("NAME OF YOUR PYKOGNITION OR ANNOTATED FILE.csv")

## Optional: If you've added cluster numbers to adlib_id/imageName, you need to remove them before joining
## In this case, we concatenated cluster numbers with imageName using "_" as a separator

AWS_data$imageName <-str_replace_all(AWS_data$imageName,".*_", "")

## Join the AWS_Data with IAD_data
labelled_data <- full_join(IAD_data, AWS_data, by = "imageName")

## We use fill to extend our labels from core images to duplicates 

labelled_data <- labelled_data %>% 
  group_by(group_id) %>% 
  fill(Page:Type, .direction = "down")

## At this stage, you may see many NAs in your data
## You need to filter out images that are not labelled
## In our case, we can do this using the Emotion column, but the exact column you need may vary based on your data

labelled_data <- labelled_data %>% 
  filter(!is.na(Emotion))
  
## To prep your labeled data for hydration with the Ad Library API:
## You need to remove the ".jpg" tag from images and rename imageName to adlib_id
## Here we add a new column called ad_lib instead of renaming, while removing ".jpg"

labelled_data$adlib_id <- str_remove_all(labelled_data$imageName,"\\\\.jpg")

## Now you can load in your Ad Library API file (at the time of writing only .xlsx is supported)

ad_lib_data <- read_excel("NAME OF YOUR AD LIBRARY DATA FILE.xlsx")
ad_lib_data <- as_tibble(ad_lib_data)

## Merge the two finals to complete the hydration

hydrated_data <- merge(ad_lib_data, labelled_data, by = "adlib_id")

## Save the hydrated data

write.csv(hydrated_data, "Hydrated and Annotated Ad Library API Data.csv")

## Thanks for using our tools!
## Please don't forget to cite if you use them! 
