# Extract displacement events from feeder visitation data

# Required libraries
  library(lubridate)
  library(EloRating)

# Load in detectdisplacement function
source("eth12720-sup-0001-Code.R")

  #=============================================================================================
# Load feeder visitation data
sf21 = read.csv("SingleFeederVisits_20-21.csv")
sf20 = read.csv("SingleFeederVisits_19-20.csv")

# Get time in seconds
sf21$startTime = as.numeric(ymd_hms(sf21$startTime))
sf21$endTime = as.numeric(ymd_hms(sf21$endTime))
  
sf20$startTime = as.numeric(ymd_hms(sf20$startTime))
sf20$endTime = as.numeric(ymd_hms(sf20$endTime))
  
# Create locationDay variable 
sf21$locationDay = paste(sf21$location, sf21$date, sep = "_")
sf20$locationDay = paste(sf20$location, sf20$date, sep = "_")


# Detect displacement events, get dominance matrix
sf21_dom = detectdisplacement(arrivetime = sf21$startTime,
                                 departtime = sf21$endTime,
                                 ids = sf21$tag,
                                 sites = sf21$location,
                                 divisions = sf21$locationDay,
                                 displacetime = 1,
                                 remaintime = 2)


sf20_dom = detectdisplacement(arrivetime = sf21$startTime,
                              departtime = sf21$endTime,
                              ids = sf21$tag,
                              sites = sf21$location,
                              divisions = sf21$locationDay,
                              displacetime = 1,
                              remaintime = 2)

# Pull displacement interactions 
disp_21 = sf21_dom$interactions
disp_20 = sf20_dom$interactions


# Get groups 
groups_21 = read.csv("Communities_20-21.csv")
groups_20 = read.csv("Communities_19-20.csv")


# Match winner and loser IDs to community IDs
# Only use interactions where both winner and loser are members of the same community

disp_21$wgroup = groups_21$finalcommunity[match(unlist(disp_21$win), groups_21$tag)]
disp_21$lgroup = groups_21$finalcommunity[match(unlist(disp_21$loss), groups_21$tag)]
disp_21$group = ifelse(disp_21$wgroup == disp_21$lgroup, disp_21$wgroup, NA)
disp_21 = disp_21[!is.na(disp_21$group),]
disp_21$wgroup = NULL
disp_21$lgroup = NULL

disp_20$wgroup = groups_20$finalcommunity[match(unlist(disp_20$win), groups_20$tag)]
disp_20$lgroup = groups_20$finalcommunity[match(unlist(disp_20$loss), groups_20$tag)]
disp_20$group = ifelse(disp_20$wgroup == disp_20$lgroup, disp_20$wgroup, NA)
disp_20 = disp_20[!is.na(disp_20$group),]
disp_20 = disp_20 %>% select(-c(wgroup, lgroup))

# Combine seasons
disp_21$season = "S20-21"
disp_20$season = "S19-20"

disp_21$group = paste0(disp_21$group, "_20-21")
disp_20$group = paste0(disp_20$group, "_19-20")

disp_both = rbind(disp_20, disp_21)

write.csv(disp_both, file="DisplacementEvents.csv", row.names=F)

#==================================================================================
# Get within-group rankings

# Function to subset by group and calculate David's scores w/in that group
scores_bygroup = function(dat, groupname) {
  
  dat = dat %>% 
    filter(group == groupname) %>%
    arrange(dt)
  
  seq = elo.seq(winner = dat$win, 
                loser=dat$loss, 
                Date=dat$date, 
                progressbar=F)
  mat = creatematrix(seq)
  scores = DS(mat)
  scores$group = groupname
  scores$season = dat[1,7]
  scores$rank = rank(scores[,2])
  
  return(scores)
}

# Loop through groups and get rankings
superlist = vector("list", length=length(unique(disp_both$group)))  

for (i in 1:length(unique(disp_both$group))) {
  
  scores = scores_bygroup(disp_both, unique(disp_both$group)[[i]])
  superlist[[i]] = scores
}

# Combine scores
scores = do.call("rbind", superlist)

# These scores were then combined with cognitive testing data, etc 
# The final output is saved as "Ranks_CogScores.csv"

#====================================================================================

