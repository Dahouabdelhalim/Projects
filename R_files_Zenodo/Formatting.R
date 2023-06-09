source ("Library.R")

Minute = 20
Frame = Minute * 60 * 30 #fps:30

InitialCut <- 300;

adjust_Dmel = 127/838
adjust_Calo = 127/823
adjust_IR = 127/419
Adjust = c(adjust_Dmel, adjust_Calo, adjust_IR)

coln <- as.integer(1+6*(0:19));
colx <- as.integer(2+6*(0:19));
coly <- as.integer(3+6*(0:19));

File_Name = c("Dmel_f", "Calo_f", "IR_f") 
Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")

data_all <- c()
for (Type in 1:length(File_Name)) {
  source_csv_files = list.files(pattern = File_Name[Type])
  
  for (FileID in 1:length(source_csv_files)) {
    current_file <- source_csv_files[FileID];
    sample <- read_csv(current_file, col_names = F);
    samplecut <- sample[(1 + InitialCut):(1 + InitialCut + Frame),] #20分間のデータにカット
    sampleadjust <- samplecut %>%
      mutate(across(c(colx, coly), ~(.x * Adjust[Type]))) #ピクセルからmmに変換
    data <- as.data.frame(sampleadjust)
    
    for(FlyID in 1:20){
      Frame_num <- c(0:Frame)
      Fly_num <- data[, coln[FlyID]]
      X <- data[, colx[FlyID]]
      Y <- data[, coly[FlyID]]
      data_one <- data.frame(Type = Type_Name[Type], 
                             Movie_number = as.character(FileID),
                             Fly_number = Fly_num,
                             Frame_number = Frame_num,
                             x = X,
                             y = Y)
      
      data_all <- bind_rows(data_all, data_one)
    }
  }
}

#---------------------------------------------------------------------------#

data_second = c()
cut = 36000
for (cut in seq(0, 36000, 30)) {
  print(cut)
  data_cut = data_all %>%
    filter(Frame_number == cut) %>%
    mutate(Time = Frame_number/30) %>%
    dplyr::select(1, 2, 3, 7, 5, 6)
  data_second = bind_rows(data_second, data_cut)
}

write.csv(data_second, "Coordinate_second.csv")

#---------------------------------------------------------------------------#

Data = data_all %>%
  filter(Frame_number <= 36000)

Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")

State_value = c(1.0, 1.0, 2.0)

data_all = c()
for (TypeID in 1:length(Type_Name)) {
  data_type = Data %>%
    filter(Type == Type_Name[TypeID])
  
  Threshold_Speed = State_value[TypeID]
  Def_State = function(Speed)
  {
    if(Speed < Threshold_Speed){
      Speed <- 0
    }
    else{
      Speed <- 1
    }
  }
  
  for (MovieID in 1:max(data_type$Movie_number)) {
    print(MovieID)
    
    for (FlyID in 0:19) {
      data_one = data_type %>%
        filter(Movie_number == MovieID, Fly_number == FlyID) %>%
        mutate(xx = lag(x), yy = lag(y))
      
      data_displacement <- sqrt((data_one$x - data_one$xx)^2 + (data_one$y - data_one$yy)^2)
      df <- data.frame(data_one[1:4], Displacement = data_displacement) %>%
        drop_na()
      
      data_speed <- c()
      for (cut in 1:1200) {
        df_cut = df[c(1:30) + 30 * (cut - 1), ]
        sum_displacement <- sum(df_cut$Displacement)
        data_speed <- append(data_speed, sum_displacement)
      }
      
      data_state <- mapply(Def_State, data_speed)
      
      data_result <- data.frame(Type = Type_Name[TypeID], 
                                Movie_number = as.character(MovieID),
                                Fly_number = as.character(FlyID),
                                Time = c(1:1200),
                                Speed = data_speed,
                                State = data_state)
      
      data_all <- bind_rows(data_all, data_result)
    }
  }
}

write.csv(data_all, "Speed.csv")

#---------------------------------------------------------------------------#

ad_Dmel = 127/838
ad_Calo = 127/823

mean_sd <- list(
  mean = ~mean(.x, na.rm = TRUE),
  sd = ~sd(.x, na.rm = TRUE))

Size = readr::read_csv("Size_individual.csv") %>%
  summarise('C.alo' = Calocasiae * ad_Calo, 'D.mel' = Dmelanogaster * ad_Dmel) %>%
  summarise(across(where(is.numeric), mean_sd))

write.csv(Size, "Size_mean_sd.csv")