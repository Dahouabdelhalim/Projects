library(jsonlite)            	
library(readr)
library(dplyr)			

substrLeft = function(text, num_char) {
  substr(text, 1, num_char)
}


data_directory <- "D:\\\\WoS\\\\"
FOS_file <- "D:\\\\FOS_WoS.csv"
MyFOS <- read.csv(FOS_file, header=TRUE, sep=",")
discipline_order <- c(1,4,5,3,2)

current_file <- "WoSwFOS_Random"
MyWoS <- read.csv(paste(data_directory,"\\\\",current_file,".csv",sep=""),header=TRUE, sep=",")

Current_FOS_list <- list() 

n_row <- nrow(MyWoS)

MyPaper <- list()
for (list_count in 1:5)
{MyPaper[list_count] <- NULL}
FOS_found <- ""
assignment_found <- TRUE

for (i in 1:n_row) 
{
  if ((assignment_found)&(MyWoS[i,3]>1))
  {
  }else
  {
    array_found <- seq(from = 0, to = 0, length.out = 5)
    
    for (j in 1:length(MyFOS))
    {
      Current_FOS_list <- MyFOS[,j] 
      Current_FOS_list <- Current_FOS_list[Current_FOS_list!=""] 
      current_record <- MyWoS[i,2]
      position_found <- which(current_record== Current_FOS_list)
      
      if (length(position_found)==0)
      {
        comma_position <- gregexpr(pattern =',',current_record)
        if (comma_position[[1]][1]>0)
        {
          left_phrase <- substrLeft(current_record, comma_position[[1]][1]-1)
          position_found <- which(left_phrase== Current_FOS_list) 
        }
      }else
      {
        list_phrase <- strsplit(current_record, ",")
        if (length(list_phrase[[1]])>1)
        {
          partly_found<-FALSE
          for (k in 1: length(list_phrase[[1]]))
          {
            partly_found <- partly_found | length(which(gsub(" ", "", list_phrase[[1]][k])== Current_FOS_list))>0
          }
          if (partly_found == TRUE)
          {
            result_file<-paste(data_directory, "\\\\", current_file,  "_", as.character(discipline), "-additional.csv", sep="")
            write.table(MyWoS[i,], result_file, append = TRUE,sep = ",", row.names = FALSE, col.names = FALSE)
          }
        }
      }
      
      
      if (length(position_found)>0)
      {
        if (sum(array_found)==0) 
        {
          array_found[j] <- array_found[j]+1
        }
      }
    }
    
    if (max(array_found)>0)
    {
      assignment_found <- TRUE
      FOS_found <- which(max(array_found)==array_found)[1]
      discipline <- discipline_order[FOS_found]
      result_file<-paste(data_directory, "\\\\", current_file,  "_", as.character(discipline), ".csv", sep="")
      
      write.table(MyWoS[i,], result_file, append = TRUE,sep = ",", row.names = FALSE, col.names = FALSE)
    }
    else
    {
      assignment_found <- FALSE
    }
    
  }
  
}

