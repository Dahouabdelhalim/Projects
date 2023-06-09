# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 1
# ------------------------------------------------

# Load Libraries
library(dplyr)
library(readr)

# Set working directory
setwd("~/Data Repository/eBird Data")

# Connecticut
Connecticut <- read_csv("Routine 1/Connecticut_1.csv")
length(unique(Connecticut$`COMMON NAME`))
Connecticut <- Connecticut %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Connecticut")
write_csv(Connecticut, "Routine 1A/Connecticut_1_summary.csv")
gc()

# Delaware
Delaware <- read_csv("Routine 1/Delaware_1.csv")
length(unique(Delaware$`COMMON NAME`))
Delaware <- Delaware %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Delaware")
write_csv(Delaware, "Routine 1A/Delaware_1_summary.csv")
gc()


# Illinois
Illinois <- read_csv("Routine 1/Illinois_1.csv")
length(unique(Illinois$`COMMON NAME`))
Illinois <- Illinois %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Illinois")
write_csv(Illinois, "Routine 1A/Illinois_1_summary.csv")
gc()


# Indiana
Indiana <- read_csv("Routine 1/Indiana_1.csv")
length(unique(Indiana$`COMMON NAME`))
Indiana <- Indiana %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Indiana")
write_csv(Indiana, "Routine 1A/Indiana_1_summary.csv")
gc()


# Kentucky
Kentucky <- read_csv("Routine 1/Kentucky_1.csv")
length(unique(Kentucky$`COMMON NAME`))
Kentucky <- Kentucky %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Kentucky")
write_csv(Kentucky, "Routine 1A/Kentucky_1_summary.csv")
gc()

# Maine
Maine <- read_csv("Routine 1/Maine_1.csv")
length(unique(Maine$`COMMON NAME`))
Maine <- Maine %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Maine")
write_csv(Maine, "Routine 1A/Maine_1_summary.csv")
gc()


# Maryland
Maryland <- read_csv("Routine 1/Maryland_1.csv")
length(unique(Maryland$`COMMON NAME`))
Maryland <- Maryland %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Maryland")
write_csv(Maryland, "Routine 1A/Maryland_1_summary.csv")
gc()


# Massachusetts
Massachusetts <- read_csv("Routine 1/Massachusetts_1.csv")
length(unique(Massachusetts$`COMMON NAME`))
Massachusetts <- Massachusetts %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Massachusetts")
write_csv(Massachusetts, "Routine 1A/Massachusetts_1_summary.csv")
gc()


# Michigan
Michigan <- read_csv("Routine 1/Michigan_1.csv")
length(unique(Michigan$`COMMON NAME`))
Michigan <- Michigan %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Michigan")
write_csv(Michigan, "Routine 1A/Michigan_1_summary.csv")
gc()


# New_Brunswick
New_Brunswick <- read_csv("Routine 1/New_Brunswick_1.csv")
length(unique(New_Brunswick$`COMMON NAME`))
New_Brunswick <- New_Brunswick %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_New_Brunswick")
write_csv(New_Brunswick, "Routine 1A/New_Brunswick_1_summary.csv")
gc()


# New_Hampshire
New_Hampshire <- read_csv("Routine 1/New_Hampshire_1.csv")
length(unique(New_Hampshire$`COMMON NAME`))
New_Hampshire <- New_Hampshire %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_New_Hampshire")
write_csv(New_Hampshire, "Routine 1A/New_Hampshire_1_summary.csv")
gc()


# New_Jersey
New_Jersey <- read_csv("Routine 1/New_Jersey_1.csv")
length(unique(New_Jersey$`COMMON NAME`))
New_Jersey <- New_Jersey %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_New_Jersey")
write_csv(New_Jersey, "Routine 1A/New_Jersey_1_summary.csv")
gc()


# New_York
New_York <- read_csv("Routine 1/New_York_1.csv")
length(unique(New_York$`COMMON NAME`))
New_York <- New_York %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_New_York")
write_csv(New_York, "Routine 1A/New_York_1_summary.csv")
gc()


# Nova_Scotia
Nova_Scotia <- read_csv("Routine 1/Nova_Scotia_1.csv")
length(unique(Nova_Scotia$`COMMON NAME`))
Nova_Scotia <- Nova_Scotia %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Nova_Scotia")
write_csv(Nova_Scotia, "Routine 1A/Nova_Scotia_1_summary.csv")
gc()


# Ohio
Ohio <- read_csv("Routine 1/Ohio_1.csv")
length(unique(Ohio$`COMMON NAME`))
Ohio <- Ohio %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Ohio")
write_csv(Ohio, "Routine 1A/Ohio_1_summary.csv")
gc()


# Ontario
Ontario <- read_csv("Routine 1/Ontario_1.csv")
length(unique(Ontario$`COMMON NAME`))
Ontario <- Ontario %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Ontario")
write_csv(Ontario, "Routine 1A/Ontario_1_summary.csv")
gc()


# Pennsylvania
Pennsylvania <- read_csv("Routine 1/Pennsylvania_1.csv")
length(unique(Pennsylvania$`COMMON NAME`))
Pennsylvania <- Pennsylvania %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Pennsylvania")
write_csv(Pennsylvania, "Routine 1A/Pennsylvania_1_summary.csv")
gc()


# Quebec
Quebec <- read_csv("Routine 1/Quebec_1.csv")
length(unique(Quebec$`COMMON NAME`))
Quebec <- Quebec %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Quebec")
write_csv(Quebec, "Routine 1A/Quebec_1_summary.csv")
gc()


# Rhode_Island
Rhode_Island <- read_csv("Routine 1/Rhode_Island_1.csv")
length(unique(Rhode_Island$`COMMON NAME`))
Rhode_Island <- Rhode_Island %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Rhode_Island")
write_csv(Rhode_Island, "Routine 1A/Rhode_Island_1_summary.csv")
gc()


# Vermont
Vermont <- read_csv("Routine 1/Vermont_1.csv")
length(unique(Vermont$`COMMON NAME`))
Vermont <- Vermont %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Vermont")
write_csv(Vermont, "Routine 1A/Vermont_1_summary.csv")
gc()


# Virginia
Virginia <- read_csv("Routine 1/Virginia_1.csv")
length(unique(Virginia$`COMMON NAME`))
Virginia <- Virginia %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Virginia")
write_csv(Virginia, "Routine 1A/Virginia_1_summary.csv")
gc()


# West_Virginia
West_Virginia <- read_csv("Routine 1/West_Virginia_1.csv")
length(unique(West_Virginia$`COMMON NAME`))
West_Virginia <- West_Virginia %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_West_Virginia")
write_csv(West_Virginia, "Routine 1A/West_Virginia_1_summary.csv")
gc()


# Wisconsin
Wisconsin <- read_csv("Routine 1/Wisconsin_1.csv")
length(unique(Wisconsin$`COMMON NAME`))
Wisconsin <- Wisconsin %>% 
  group_by(`COMMON NAME`) %>% 
  count(name = "n_Wisconsin")
write_csv(Wisconsin, "Routine 1A/Wisconsin_1_summary.csv")
gc()



# Join together
Birds <- left_join(Connecticut, Delaware) %>% 
  left_join(Illinois) %>% 
  left_join(Indiana) %>% 
  left_join(Kentucky) %>% 
  left_join(Maine) %>% 
  left_join(Maryland) %>% 
  left_join(Massachusetts) %>% 
  left_join(Michigan) %>% 
  left_join(New_Brunswick) %>% 
  left_join(New_Hampshire) %>% 
  left_join(New_Jersey) %>% 
  left_join(New_York) %>% 
  left_join(Nova_Scotia) %>% 
  left_join(Ohio) %>% 
  left_join(Ontario) %>% 
  left_join(Pennsylvania) %>% 
  left_join(Quebec) %>% 
  left_join(Rhode_Island) %>% 
  left_join(Vermont) %>% 
  left_join(Virginia) %>% 
  left_join(West_Virginia) %>% 
  left_join(Wisconsin)
# Create new column
Birds <- Birds %>% 
  mutate(n_all = rowSums(across(where(is.numeric)), na.rm = TRUE))

write_csv(Birds, "Routine 1A/Birds.csv")
