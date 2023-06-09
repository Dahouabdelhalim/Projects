library(jsonlite)            	
library(readr)
library(jsonlite)            	
library(dplyr)			

data_directory <- "D:\\\\MAG\\\\unzipped-all\\\\"
FOS_file <- "FOS_name_list_top50from10000.csv"
MyFOS <- read.csv(FOS_file, header=TRUE, sep=",")

Current_FOS_list <- list() 

n_row <- nrow(MyFOS)

file_list <-list.files(data_directory)

for (current_file in file_list) 
{
  print(paste("FILE: ", current_file))
  
  MyPaper <- list()
  for (list_count in 1:22)
  {MyPaper[list_count] <- NULL}
  
  con = file(paste(data_directory,current_file,sep=""), "r")
  
  while ( TRUE ) 
  {
    an.error.occured <- FALSE
    tryCatch( { line = readLines(con, n = 1)}
              , error = function(e) {an.error.occured <- TRUE})
    if (an.error.occured==TRUE)
    {
      print("Error in reading line")
      break
    }
    else
    {
      if ( length(line) == 0 ) {
        break
      }
      
      paper <- jsonlite::fromJSON(line)
      
      pfos <- paper$fos
      if (length(pfos)>0)
      {
        array_found <- seq(from = 0, to = 0, length.out = 22)
        
        for (i in 1:length(pfos))
        {
          for (j in 1:length(MyFOS))
          {
            Current_FOS_list <- MyFOS[,j] 
            Current_FOS_list <- Current_FOS_list[Current_FOS_list!=""] 
            position_found <- which(pfos[i]== Current_FOS_list) 
            
            if (length(position_found)>0)
            {
              array_found[j] <- array_found[j]+1
            }
          }
        }
        
        if (max(array_found)>0)
        {
          FOS_found <- which(max(array_found)==array_found)[1]
          result_file<-paste(as.character(FOS_found), "_", current_file, sep="")
          paper %>% toJSON() %>% write_lines(result_file, append=TRUE)
        }
      }
    }
  }
  close(con)
  break
}
