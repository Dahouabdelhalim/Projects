source ("Library.R")
source ("xmeans.R")

Size_mean = readr::read_csv("Size_mean_sd.csv") %>%
  dplyr::select(-1)

size_Calo = Size_mean$C.alo_mean
size_Dmel = Size_mean$D.mel_mean
size_IR = Size_mean$C.alo_mean

Size <- c(size_Dmel, size_Calo, size_IR)


Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")

Group = readr::read_csv("Group_judgment_correction.csv") %>%
  dplyr::select(-1)

Coordinate = readr::read_csv("Coordinate_second.csv") %>%
  dplyr::select(-1)

Data = as.data.frame(left_join(Group, Coordinate))

data_all <- c()
for (TypeID in 1:length(Type_Name)) {
  data_type = Data %>%
    filter(Type == Type_Name[TypeID])
  
  for (MovieID in 1:max(data_type$Movie_number)) {
    group_in = data_type %>%
      filter(Movie_number == MovieID,
             Group_in == 1) %>%
      arrange(Time)
    
    Time_group_in <- as.vector(unique(group_in$Time))
    group_number_all <- c()
    for (t in Time_group_in) {
      print(TypeID)
      print(MovieID)
      group_in_one = group_in %>%
        filter(Time == t)
      
      if(nrow(group_in_one) == 1) {
        group_number <- 1
      }
      
      if(nrow(group_in_one) == 2){
        data_distance <- sqrt((group_in_one[1, 8] - group_in_one[2, 8])^2 + (group_in_one[1, 9] - group_in_one[2, 9])^2)

        Threshold_Distance = Size[TypeID] * 2
        if(data_distance < Threshold_Distance) {
          group_number <- c(1, 1)
        }
        if(data_distance >= Threshold_Distance){
          group_number <- c(1, 2)
        }
      }
      
      if(nrow(group_in_one) >= 3){
        xkm  <- xmeans$new()
        xkm$ishioka_xmeans(group_in_one[8:9])
        group_number <- xkm$Cluster
      }
      
      group_number_all <- append(group_number_all, group_number)
    }
    data_result <- data.frame(group_in, Group_number = group_number_all) %>%
      dplyr::select(-c(Proximate_population, Group_in, x, y))
    data_all = bind_rows(data_all, data_result)
  }
}

write.csv(data_all, "Group_identification.csv")