source ("Library.R")

Size_mean = readr::read_csv("Size_mean_sd.csv") %>%
  dplyr::select(-1)

size_Calo = Size_mean$C.alo_mean
size_Dmel = Size_mean$D.mel_mean
size_IR = Size_mean$C.alo_mean

Size <- c(size_Dmel, size_Calo, size_IR)

Coordinate = readr::read_csv("Coordinate_second.csv") %>%
  dplyr::select(-1)

Data <- as.data.frame(Coordinate)

Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")


data_all = c()
for (TypeID in 1:length(Type_Name)) {
  data_type = Data %>%
    filter(Type == Type_Name[TypeID])
  
  Threshold_Distance = Size[TypeID] * 2
  Def_Distance = function(Distance)
  {
    if(Distance < Threshold_Distance){
      Distance <- 1
    }
    else{
      Distance <- 0
    }
  }
  
  for (MovieID in 1:max(data_type$Movie_number)) {
    print(MovieID)
    
    for (FlyID in 0:19) {
      me = data_type %>%
        filter(Movie_number == MovieID, Fly_number == FlyID)
      
      other = data_type %>%
        filter(Movie_number == MovieID) %>%
        rename(Other_number = Fly_number,
               Other_Time = Time,
               Other_x = x,
               Other_y = y) %>%
        arrange(Other_number)
      
      data_bind <- data.frame(me, other[3:6])
      
      data_distance <- sqrt((data_bind$x - data_bind$Other_x)^2 + (data_bind$y - data_bind$Other_y)^2)
      
      judge_distance <- mapply(Def_Distance, data_distance)
      
      data_judge <- data.frame(me[1:4], other[3], Judge = judge_distance)
      
      data_result = data_judge %>%
        group_by(Type, Movie_number, Fly_number, Time) %>%
        summarise(Proximate_population = sum(Judge)) %>%
        ungroup
      
      data_all <- bind_rows(data_all, data_result)
    }
  }
}


Threshold_Population = 3
Def_Population = function(Population)
{
  if(Population < Threshold_Population){
    Population <- 0
  }
  else{
    Population <- 1
  }
}

judge_population <- mapply(Def_Population, data_all$Proximate_population)

Result <- data.frame(data_all, Group_in = judge_population)

write.csv(Result, "Group_judgment.csv")



#-------------------------------------------------------------------------------------------------------------

Data = readr::read_csv("Group_judgment.csv") %>%
  dplyr::select(-1)

Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")

Out_Time = 5

Result_correction = c()
for (TypeID in 1:length(Type_Name)) {
  data_type = Data %>%
    filter(Type == Type_Name[TypeID])
  
  for (MovieID in 1:max(data_type$Movie_number)) {
    print(MovieID)
    for (FlyID in 0:19) {
      data_one = data_type %>%
        filter(Movie_number == MovieID, Fly_number == FlyID) %>%
        mutate(Lag = lag(Group_in)) %>%
        mutate(State = Group_in - Lag) %>%
        dplyr::select(-c(Lag))
      
      df_time = data_one %>%
        filter(State == 1 | State == -1) %>%
        mutate(Lag=lag(Time)) %>%
        mutate(State_Time = Time - Lag) %>%
        dplyr::select(-c(Lag)) %>%
        drop_na() %>%
        filter(State == 1, State_Time <= Out_Time) %>%
        mutate(Out_Start = Time - State_Time, Out_End = Time - 1) %>%
        dplyr::select(Out_Start, Out_End)

      if(nrow(df_time)>0) {
        Start_Time <- as.vector(df_time$Out_Start)
        End_Time <- as.vector(df_time$Out_End)
        
        for (t in 1:nrow(df_time)) {
          data_one[((Start_Time[t] + 1):(End_Time[t] + 1)), 6] <- 1
        }
      }
      
      data_correction = data_one %>%
        dplyr::select(-State) %>%
        mutate(Lag=lag(Group_in)) %>%
        mutate(State = Group_in - Lag) %>%
        dplyr::select(-c(Lag)) %>%
        drop_na()
      
      Result_correction = bind_rows(Result_correction, data_correction)
    }
  }
}

write.csv(Result_correction, "Group_judgment_correction.csv")