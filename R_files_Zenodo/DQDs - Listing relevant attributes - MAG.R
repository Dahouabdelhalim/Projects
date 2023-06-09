library(jsonlite)            	# for opening json files
library(stringr)            	# for processing strings
library(dplyr)			          # for grouping and summarizing
library(RCurl)			          # for composing general HTTP requests

f1_min <- 1
f1_max <- 166
f2_min <- 1
f2_max <- 5

f1 <- f1_min
f2 <- f2_min 

file <- data.frame()
temp <- NULL

for (s in f1_min:f1_max) 
{
  for (t in f2_min:f2_max)
  { 
    current_file <- paste("MAG_filtered_",toString(f1),"_",toString(f2),".txt",sep="")
    print(current_file)
    j <- 1
    
    Sys.time()
    
    con = file(current_file, "r")
    while ( TRUE ) 
    {
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }

      paper <- jsonlite::fromJSON(line)
      
      print(paper$authors$name)
      print(paper$authors$org)
      
      if ((length(paper$authors$name)<1) & (length(paper$authors$org)<1))
      {
        print("Error!")
        print(length(paper$authors$name))
        print(paper$authors$name)
        print(length(paper$authors$org))
        print(paper$authors$org)
      }
      else
      {
        
        if (length(paper$authors$name)==0)
        {
          for (i in 1:length(paper$authors$org))
          {
            temp$year <- paper$year
            temp$name <- "Empty value"
            temp$org <- paper$authors$org[i]
            file <- rbind(file, c(temp$year,temp$name,temp$org))  
            
            print(c(temp$year,temp$name,temp$org))
            readline()
            
          }
        }
        else
        {
          if (length(paper$authors$org)==0)
          {
            for (i in 1:length(paper$authors$name))
            {
              temp$year <- paper$year
              temp$name <- paper$authors$name[i]
              temp$org <- "Empty value"
              file <- rbind(file, c(temp$year,temp$name,temp$org))     
              
              print(c(temp$year,temp$name,temp$org))
              readline()
              
            }
          }
          else 
          {
            min_count <- min(length(paper$authors$name),length(paper$authors$org))
            max_count <- max(length(paper$authors$name),length(paper$authors$org))
            
            for (i in 1:min_count)
            {
              temp$year <- paper$year
              temp$name <- paper$authors$name[i]
              temp$org <- paper$authors$org[i]
              file <- rbind(file, c(temp$year,temp$name,temp$org))    
              
              print(c(temp$year,temp$name,temp$org))
              readline()
              
            }
            
            if (length(paper$authors$name)>length(paper$authors$org))
            {
              for (i in (min_count+1):max_count)
              {
                temp$year <- paper$year
                temp$name <- paper$authors$name[i]
                temp$org <- "Empty value"
                file <- rbind(file, c(temp$year,temp$name,temp$org))    
                
                print(c(temp$year,temp$name,temp$org))
                readline()
                
              }
            }
            if (length(paper$authors$name)<length(paper$authors$org))
            {
              for (i in (min_count+1):max_count)
              {
                temp$year <- paper$year
                temp$name <- "Empty value"
                temp$org <- paper$authors$org[i]
                file <- rbind(file, c(temp$year,temp$name,temp$org))    
                
                print(c(temp$year,temp$name,temp$org))
                readline()
                
              }
            }    
          }
        }
      }
    }

    close(con)
    
    f2 <- f2+1
    if (f2>f2_max) 
    {
      f2 <- 1
      f1 <- f1+1
      if (f1>f1_max) {break}
    }
    if ((f1==166) &(f2==2)) {break}
  }
}

