# downloading lilac phenometrics data for comparing with Constant Contact message open stats
# T Crimmins (theresa@usanpn.org) updated 1-6-21

library(rnpn)
library(dplyr)

lilacdata <- npn_download_individual_phenometrics(
                                     request_source = "TCrimmins",
                                     years=c("2019"),
                                     species_ids = c(35,36), 
                                     six_leaf_layer = TRUE, 
                                     six_sub_model = "lilac",
                                     )


# subset to only LEAF or BLOOM 
# unique(lilacdata$phenophase_description)

lilacbloom <- subset(lilacdata, phenophase_description=="Open flowers (lilac)")
#OR
lilacleaf <- subset(lilacdata, phenophase_description=="Breaking leaf buds (lilac/honeysuckle)")

# thin out multiple first y's -- only one individual_id per year for a single phenophase
lilacleafsub <- lilacleaf %>% group_by(individual_id,first_yes_year) %>% filter(first_yes_doy==min(as.numeric(as.character(first_yes_doy))))


# load 20XX LEAF or BLOOM message open file
leaf2019msgstats <- read.csv(file = "c:/NPN/Manuscripts/Working/Springcasting_evaluation/data_files/2019_Leaf_msg_stats.csv")

# fix column name, if necessary
colnames(leaf2019msgstats)[1]<-"Station_ID"

# join w/ConCon message open info 
###I THINK ON THIS STEP I WANT TO RETAIN ALL ROWS - even situations where folks didn't log obs - if what this option, keep the "all.y = TRUE"
leaf2019merge <- merge(x = lilacleafsub, y = leaf2019msgstats, by.x = "site_id", by.y = "Station_ID", all.x = TRUE, all.y = TRUE)

# Look at number of rows in MERGED df - if greater than Pheno data, need to do next steps:
# in 2018, had duplicates - multiple people received emails for same plants

# assign codes for email open codes:
# Opened = 1, Did not open = 2, Did not receive = 3 - Ensure the order is correct here!!
leaf2019merge$Code2 <- c("2", "3", "1")[leaf2019merge$LeafCode]

# generate stats on open rate & data logging rate
table <- leaf2019merge %>%
  group_by(first_yes_year, LeafCode) %>%
  tally()



# CHECK WHAT's REMOVED IN NEXT STEPS.....
# remove duplicates based on code
leaf2019trim <- leaf2019merge %>% group_by(individual_id) %>% filter(Code2==min(Code2))

# remove additional duplicates - this preps file for other analyses
leaf2019trim2 <- leaf2019trim[!duplicated(leaf2019trim$individual_id), ]

# save data file...fix filename here
write.csv(leaf2019trim, "c:/NPN/Manuscripts/Working/Springcasting_evaluation/merged_files/2019_lilac_leaf-submission_rate.csv")



