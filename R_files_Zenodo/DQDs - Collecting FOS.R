library(jsonlite)            
library(dplyr)			

data_directory <- "D:\\\\MAG\\\\unzipped-all\\\\"
file_name <- "FOS-from10000.csv"
Title_file <- "Titles-10000.csv"
MyData <- read.csv(Title_file, header=TRUE, sep=",")

n_row <- nrow(MyData)

file_list <-list.files(data_directory)

for (current_file in file_list) {
    k <- 1 
    
    print(paste("FILE: ", current_file))

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
            
            ptitle <- paper$title
            pyear <- paper$year
            
            if(length(which(MyData$Title == ptitle))>0) 
            {
              location <- which(MyData$Title == ptitle)
              discipline <- MyData$Discipline[location]
              if (length(paper$fos)>0)
              {
                file <- data.frame()
                
                for (fos_count in 1: length(paper$fos))
                {
                  pfos <- paper$fos[fos_count]
                  file <- rbind(file, c(current_file, discipline, ptitle, pfos))
                }
                write.table(file,file=file_name,col.names=F , sep = ",",row.names = F, append = TRUE)  
              }
            } 
          }
        }
      
        close(con)
      
      }

