source ("Library.R")

adjust_Dmel = 127/838
adjust_Calo = 127/823
adjust_IR = 127/419
Adjust = c(adjust_Dmel, adjust_Calo, adjust_IR)

Center = readr::read_csv("Center_coordinate.csv")
Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")

Center_adjust = c()
for (TypeID in 1:length(Type_Name)) {
  Center_type = Center %>%
    filter(Type == Type_Name[TypeID]) %>%
    mutate(x = x * Adjust[TypeID], y = y * Adjust[TypeID])
  Center_adjust = bind_rows(Center_adjust, Center_type)
}

Data = readr::read_csv("Coordinate_second.csv") %>%
  dplyr::select(-1)

data_Adjust <- c()
for (TypeID in 1:length(Type_Name)) {
  data_type = Data %>%
    filter(Type == Type_Name[TypeID])
  
  for (MovieID in 1:max(data_type$Movie_number)) {
    data_movie = data_type %>%
      filter(Movie_number == MovieID)
    center_movie = Center_adjust %>%
      filter(Type == Type_Name[TypeID], Movie_number == MovieID)
    
    Adjust_x <- (data_movie$x - center_movie$x)
    Adjust_y <- (data_movie$y - center_movie$y)
    
    df <- data.frame(data_movie[, 1:4],
                     x = Adjust_x,
                     y = Adjust_y)
    
    data_Adjust = bind_rows(data_Adjust, df)
  }
}


Random_data <- c()
for (TypeID in 1:length(Type_Name)) {
  data_Adjust_type = data_Adjust %>%
    filter(Type == Type_Name[TypeID])
  
  for (Trial_num in 1:max(data_Adjust_type$Movie_number)) {
    print(Trial_num)
    random_one <- data.frame(sample_n(tbl = data_Adjust_type[, -c(2:4)], size = 20), 
                             Trial_number = as.numeric(Trial_num),
                             Fly_number = as.numeric(0:19)) %>%
      dplyr::select(1, 4, 5, 2, 3)
    Random_data <- bind_rows(Random_data, random_one)
  }
}