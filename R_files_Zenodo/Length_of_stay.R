source ("Library.R")

Data = readr::read_csv("Group_judgment_correction.csv") %>%
  dplyr::select(-1)

Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")


data_all <- c()
Stay_time = 10
for (TypeID in 1:length(Type_Name)) {
  data_type = Data %>%
    filter(Type == Type_Name[TypeID])
  
  for (MovieID in 1:max(data_type$Movie_number)) {
    print(MovieID)
    
    for (FlyID in 0:19) {
      data_one = data_type %>%
        filter(Movie_number == MovieID, Fly_number == FlyID)
      
      data_state = data_one %>%
        filter(State != 0)
      
      if(nrow(data_state) > 0){
        
        if(data_one[1200, 6] == 1){
          data_last = data_one[1200, ]
          data_last$State <- -1
          data_state = bind_rows(data_state, data_last)
        }
        
        data_belong = data_state %>%
          mutate(Lag = lag(Time)) %>%
          mutate(Belong_time = Time - Lag - 1) %>%
          filter(State == -1, Belong_time >= Stay_time) %>%
          mutate(Total_duration = sum(Belong_time)/60,
                 Mean_duration = mean(Belong_time)/60) %>%
          add_tally() %>%
          rename(Number_of_stay_events = n) %>%
          slice(1) %>%
          dplyr::select(Type, 
                        Movie_number, 
                        Fly_number, 
                        Total_duration, 
                        Mean_duration, 
                        Number_of_stay_events)
        
        data_all = bind_rows(data_all, data_belong)
      }
    }
  }
}

write.csv(data_all, "Length_of_stay.csv")