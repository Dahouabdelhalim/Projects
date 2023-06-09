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
cola <- as.integer(6+6*(0:19));

Size_mean = readr::read_csv("Size_mean_sd.csv") %>%
  dplyr::select(-1)
size_Calo = Size_mean$C.alo_mean
size_Dmel = Size_mean$D.mel_mean
size_IR = Size_mean$C.alo_mean
Size <- c(size_Dmel, size_Calo, size_IR)

SpeedThreshold = c(1.0, 1.0, 2.0) * 128/30

File_Name = c("Dmel_f", "Calo_f", "IR_f") 
Type_Name = c("Dmelanogaster", "Calocasiae", "Calocasiae(IR)")

numNN_all <- c()
lambda_all <- c()
for (Type in 1:length(File_Name)){
  source_csv_files <- list.files(pattern = File_Name[Type])
  
  for(N_times in 4:6){
    print(N_times)
    DefDistance = Size[Type] * N_times
    
    for(fileID in 1:length(source_csv_files)){
      current_file <- source_csv_files[fileID];
      sample <- read_csv(current_file, col_names = F);
      samplecut <- sample[1:(1 + InitialCut + Frame),]
      samplead <- samplecut %>%
        mutate(across(c(colx, coly), ~(.x * Adjust[Type])))
      data <- as.data.frame(samplead)
      
      for(dt in interval){
        dt_in_sec <- dt/30.0;
        num_frame <- (nrow(data) - InitialCut -1) %/% dt;
        each_frame <- c(1:num_frame);
        
        ave_numNN <-c();
        normVs <- c();
        normVs_active <- c();
        normVs_rest <- c();
        disToNN <- c();
        disToNN_active <- c();
        disToNN_rest <- c();
        min_disToNN <- c();
        numNN <- c();
        numNN_active <- c();  
        numNN_rest <- c();
        
        for(t in each_frame){
          print(t);
          now <- 1 + InitialCut + t*dt;
          last <- now - dt;
          
          sum_numNN <-0.0;
          min_sq_min_dis <- 999999999;
          
          for(me in 1:20){
            Vx <- (data[now, colx[me]]-data[last, colx[me]])/dt_in_sec;
            Vy <- (data[now, coly[me]]-data[last, coly[me]])/dt_in_sec;
            Vsq <- Vx^2 + Vy^2;
            Vnorm <- sqrt(Vsq);
            
            others <- 1:20;
            sq_min_dis <- 999999999;
            his_numNN <- 0.0;
            
            for(you in others){
              tmpsq = (data[now, colx[me]] - data[now, colx[you]])^2 + (data[now, coly[me]] - data[now, coly[you]])^2;

              if(tmpsq > 1.0){
                if(tmpsq < DefDistance){his_numNN <- his_numNN + 1};
                if(sq_min_dis > tmpsq){sq_min_dis <- tmpsq};
              }
            }
            sum_numNN <- sum_numNN + his_numNN;
            if(min_sq_min_dis > sq_min_dis){min_sq_min_dis <- sq_min_dis};
            
            normVs <- append(normVs, Vnorm);
            numNN <- append(numNN, his_numNN);
            disToNN <- append(disToNN, sqrt(sq_min_dis));
            if(Vnorm > SpeedThreshold[Type]){
              normVs_active <- append(normVs_active, Vnorm);
              numNN_active <- append(numNN_active, his_numNN);
              disToNN_active <- append(disToNN_active, sqrt(sq_min_dis));
            }else{
              normVs_rest <- append(normVs_rest, Vnorm);          
              numNN_rest <- append(numNN_rest, his_numNN);
              disToNN_rest <- append(disToNN_rest, sqrt(sq_min_dis));
            }
          }
          min_disToNN <- append(min_disToNN, sqrt(min_sq_min_dis));
          ave_numNN <- append(ave_numNN, sum_numNN/20.0);
        }
        
        fit_active <- fitdist(numNN_active,"pois")
        lambda_active  = data.frame(Type = Type_Name[Type], 
                                    Definition = as.character(N_times), 
                                    Movie_number = as.character(fileID), 
                                    Status = "active", 
                                    Lambda = fit_active[[1]][[1]], 
                                    Max_numNN = max(numNN_active))
        
        fit_rest <- fitdist(numNN_rest,"pois")
        lambda_rest  = data.frame(Type = Type_Name[Type], 
                                  Definition = as.character(N_times), 
                                  Movie_number = as.character(fileID), 
                                  Status = "rest", 
                                  Lambda = fit_rest[[1]][[1]], 
                                  Max_numNN = max(numNN_rest))
        
        
        df_numNN_active = data.frame(Type = Type_Name[Type], 
                                     Definition = as.character(N_times), 
                                     Movie_number = as.character(fileID), 
                                     active = numNN_active) %>%
          pivot_longer(active, names_to = "Status", values_to = "numNN")
        
        df_numNN_rest = data.frame(Type = Type_Name[Type], 
                                   Definition = as.character(N_times), 
                                   Movie_number = as.character(fileID), 
                                   rest = numNN_rest) %>%
          pivot_longer(rest, names_to = "Status", values_to = "numNN")
        
        lambda_all = bind_rows(lambda_all, lambda_active, lambda_rest)
        numNN_all = bind_rows(numNN_all, df_numNN_active, df_numNN_rest)
      }
    }
  }
}

write.csv(numNN_all, "Nearest_neighbor_population.csv")
write.csv(lambda_all, "Nearest_neighbor_population_lambda.csv")