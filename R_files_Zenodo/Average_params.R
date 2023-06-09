source ("Library.R")

Minute = 20
Frame = Minute * 60 * 30 #fps:30
InitialCut <- 300;
interval <- 2^(7:7);

ad_Dmel = 127/838
ad_Calo = 127/823
ad_IR = 127/419
Adjust = c(ad_Dmel, ad_Calo, ad_IR)

colx <- as.integer(2+6*(0:19));
coly <- as.integer(3+6*(0:19));

File_Name = c("Dmel_f", "Calo_f", "IR_f") 
Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")

df_params_all <- c();
for (Type in 1:length(File_Name)) {
  source_csv_files <- list.files(pattern = File_Name[Type])
  
  for(fileID in 1:length(source_csv_files)){
    current_file <- source_csv_files[fileID];
    sample <- read_csv(current_file, col_names = F);
    samplecut <- sample[1:(1 + InitialCut + Frame),]
    samplead <- samplecut %>%
      mutate(across(c(colx, coly), ~(.x * Adjust[Type])))
    data <- as.data.frame(samplead)
    
    dt = interval
    print(dt);
    dt_in_sec <- dt/30.0;
    num_frame <- (nrow(data) - InitialCut -1) %/% dt;
    each_frame <- c(1:num_frame);
    
    aveV <- c();
    aveD <-c();
    avesqD <-c();
    
    for(t in each_frame){
      print(t);
      now <- 1 + InitialCut + t*dt;
      last <- now - dt;
      
      sumVx <- 0.0;
      sumVy <- 0.0;    
      sumVnorm <- 0.0;
      
      for(fid in 1:20){
        Vx <- data[now, colx[fid]]-data[last, colx[fid]];
        Vy <- data[now, coly[fid]]-data[last, coly[fid]];
        sumVx <- sumVx + Vx;
        sumVy <- sumVy + Vy;        
        sumVnorm <- sumVnorm + sqrt(Vx^2 + Vy^2);
      }
      aveV <- append(aveV, sumVnorm/20.0/dt_in_sec)
      
      sumD <- 0.0;
      sumsqD <- 0.0;
      
      for(me in 1:19){
        others <- (me + 1):20
        
        for(you in others){
          tmpsq = (data[now, colx[me]] - data[now, colx[you]])^2 + (data[now, coly[me]] - data[now, coly[you]])^2
          sumD <- sumD + sqrt(tmpsq)
          sumsqD <- sumsqD + tmpsq
        }
      }
      aveD <- append(aveD, sumD/20.0/19.0*2)
      avesqD <- append(avesqD, sumsqD/20.0/19.0*2)
    }
    
    dt_in_min = dt/(30*60)
    df_params <- data.frame(Type = Type_Name[Type],
                            Movie_number = as.character(fileID),
                            Time = seq(dt_in_min, dt_in_min * num_frame, length = num_frame), 
                            Speed = aveV, 
                            Distance = aveD)
    
    df_params_all = bind_rows(df_params_all, df_params)
  }
}

write.csv(df_params_all, "Average_params.csv")