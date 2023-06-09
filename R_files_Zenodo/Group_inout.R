source ("Library.R")

Data = readr::read_csv("Group_judgment_correction.csv") %>%
  dplyr::select(-1)

Speed = readr::read_csv("Speed.csv") %>%
  dplyr::select(-1)

Group_number = readr::read_csv("Group_identification.csv") %>%
  dplyr::select(-1)

Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")

After_time = 10
Before_time = 5

Result <- c()
for (a in 1:length(Type_Name)) {
  d0 = Data %>%
    filter(Type == Type_Name[a])

  for (b in 1:max(d0$Movie_number)) {
    d1 = d0 %>%
      filter(Movie_number == b)

    for (c in 0:19) {
      d2 = d1 %>%
        filter(Fly_number == c)
      
      d3 = d2 %>%
        filter(Time > Before_time , Time <= (1200 - After_time)) %>%
        filter(State == 1)
      
      if(nrow(d3) > 0){
        group_in_time <- as.vector(d3$Time)

        d5 = NULL
        for (d in group_in_time) {
          d4 = d2 %>%
            filter(Time >= d - Before_time , Time < d) %>%
            summarise(Group_sum = sum(Group_in))
          d5 = bind_rows(d5, d4)
        }
        d6 = as.data.frame(bind_cols(d3, d5)) %>%
          filter(Group_sum == 0) %>%
          dplyr::select(Type, Movie_number, Fly_number, Time)
        
        if(nrow(d6) > 0){
          s1 = Speed %>%
            filter(Type == Type_Name[a], Movie_number == b, Fly_number == c)
          
          s2 = NULL
          for (e in 1:nrow(d6)) {
            s3 = s1 %>%
              filter(Time > (d6[e, 4] - Before_time), Time <= d6[e, 4]) %>%
              dplyr::select(State) %>%
              summarise(Move = sum(State))
            s2 <- bind_rows(s2, s3)
          }
          d7 = bind_cols(d6, s2) %>%
            filter(Move == 5)
          
          if(nrow(d7) > 0){
            g1 = Group_number %>%
              filter(Type == Type_Name[a], Movie_number == b)
            g2 = left_join(d7, g1)
            for (f in 1:nrow(g2)) {
              d8 = d2 %>%
                filter(Time >= g2[f, 4], Time < g2[f, 4] + After_time) %>%
                summarise(Continuous_time = sum(Group_in))
              g3 = g1 %>%
                filter(Time == g2[f, 4], Group_number == g2[f, 7]) %>%
                group_by(Type, Movie_number, Time) %>%
                tally() %>%
                rename(Group_size = n) %>%
                ungroup() %>%
                dplyr::select(Group_size)
              d9 = bind_cols(g2[f, 1:4], g3, d8)
              Result = bind_rows(Result, d9)
            }
          }
        }
      }
    }
  }
}

Threshold_Time = 10
Def_Move_inout = function(Belong_time)
{
  if(Belong_time < Threshold_Time){
    Belong_time <- "Move out"
  }
  else{
    Belong_time <- "Move in"
  }
}

Move_inout <- mapply(Def_Move_inout, Result$Continuous_time)
data_move_inout = data.frame(Result, Move_inout)

Result_move_inout = data_move_inout %>%
  dplyr::select(Type, Movie_number, Fly_number, Move_inout) %>%
  group_by(Type, Movie_number, Fly_number, Move_inout) %>%
  tally() %>%
  ungroup() %>% 
  complete(nesting(Type, Movie_number, Fly_number), Move_inout, fill = list(n = 0)) %>%
  group_by(Type, Movie_number, Fly_number) %>%
  rename(Number_of_events = n) %>%
  mutate(Total_number_of_events = sum(Number_of_events)) %>%
  mutate(Ratio = (Number_of_events/Total_number_of_events) * 100) %>%
  filter(Move_inout == "Move in")

write.csv(Result_move_inout, "Group_inout.csv")

#-------------------------------------------------------------------------------

Result_move_inout_size = data_move_inout %>%
  dplyr::select(Type, Movie_number, Fly_number, Group_size, Move_inout) %>%
  group_by(Type, Movie_number, Fly_number, Group_size, Move_inout) %>%
  tally() %>%
  ungroup() %>% 
  complete(nesting(Type, Movie_number, Fly_number, Group_size), Move_inout, fill = list(n = 0)) %>%
  group_by(Type, Movie_number, Fly_number, Group_size) %>%
  rename(Number_of_events = n) %>%
  mutate(Total_number_of_events = sum(Number_of_events)) %>%
  mutate(Ratio = (Number_of_events/Total_number_of_events) * 100) %>%
  ungroup() %>%
  filter(Group_size >= 3) %>%
  mutate(Group_size = Group_size - 1) %>%
  filter(Move_inout == "Move in") %>%
  group_by(Type, Group_size) %>%
  add_tally() %>%
  rename(Number_of_events_per_group_size = n) %>%
  ungroup()

write.csv(Result_move_inout_size, "Group_inout_size.csv")