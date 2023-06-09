source ("Library.R")

#Original data -----------------------------------------------------------------

Data = readr::read_csv("Coordinate_second.csv") %>%
  dplyr::select(-1) %>%
  filter(Time == 1200)

Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")

data_all <- c()
for (TypeID in 1:length(Type_Name)) {
  data_type = Data %>%
    filter(Type == Type_Name[TypeID])
  
  for (MovieID in 1:max(data_type$Movie_number)) {
    fly_all = data_type %>%
      filter(Movie_number == MovieID)
    
    for (FlyID in 0:19) {
      fly_one = fly_all %>%
        filter(Fly_number == FlyID)
      
      Difference_x <- (fly_all$x - fly_one$x)
      Difference_y <- (fly_all$y - fly_one$y)
      
      data_difference <- data.frame(fly_all[, 1:3], 
                                    Other_number = as.character(FlyID),
                                    Difference_x,
                                    Difference_y) %>%
        filter(Fly_number != Other_number) %>%
        mutate(Distance = sqrt(Difference_x^2 + Difference_y^2))
      
      data_all <- bind_rows(data_all, data_difference)
    }
  }
}

write.csv(data_all, "Accumulative_distribution.csv")

#Random data -------------------------------------------------------------------

source("Random_data_creation.R")

Data = Random_data

Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")

data_all <- c()
for (TypeID in 1:length(Type_Name)) {
  data_type = Data %>%
    filter(Type == Type_Name[TypeID])
  
  for (Trial in 1:max(data_type$Trial_number)) {
    fly_all = data_type %>%
      filter(Trial_number == Trial)
    
    for (FlyID in 0:19) {
      fly_one = fly_all %>%
        filter(Fly_number == FlyID)
      
      Difference_x <- (fly_all$x - fly_one$x)
      Difference_y <- (fly_all$y - fly_one$y)
      
      data_difference <- data.frame(fly_all[, 1:3], 
                                    Other_number = as.character(FlyID),
                                    Difference_x,
                                    Difference_y) %>%
        filter(Fly_number != Other_number) %>%
        mutate(Distance = sqrt(Difference_x^2 + Difference_y^2))
      
      data_all <- bind_rows(data_all, data_difference)
    }
  }
}

write.csv(data_all, "Accumulative_distribution_random.csv")