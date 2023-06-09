library(jsonlite)            	
data_directory <- "/Dimensions/"
temp <- NULL
count_errors <- 0

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))}

substrLeft = function(text, num_char) {
  substr(text, 1, num_char)
}

trim_org = function(text) {
  left_part <- text
  last_char <- substrRight(text,1)
  while (last_char %in% c(","," "))
  {
    left_part <- substrLeft(left_part, nchar(left_part)-1)
    last_char <- substrRight(left_part,1)
  }
  return(left_part)
}

get_from_affiliation = function(affiliation) {

  temp_org <- NULL
  if ((!anyNA(affiliation$name)) & (!is.null(affiliation$name)))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$name} 
    else {temp_org <- paste(temp_org, affiliation$name, sep = ", ")}
  }
  if (!anyNA(affiliation$city)&!is.null(affiliation$city))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$city} 
    else {temp_org <- paste(temp_org, affiliation$city, sep = ", ")}
  }
  if (!anyNA(affiliation$city_id)&!is.null(affiliation$city_id))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$city_id} 
    else {temp_org <- paste(temp_org, affiliation$city_id, sep = ", ")}
  }
  if (!anyNA(affiliation$country)&!is.null(affiliation$country))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$country} 
    else {temp_org <- paste(temp_org, affiliation$country, sep = ", ")}
  }
  if (!anyNA(affiliation$country_code)&!is.null(affiliation$country_code))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$country_code} 
    else {temp_org <- paste(temp_org, affiliation$country_code, sep = ", ")}
  }
  if (!anyNA(affiliation$state)&!is.null(affiliation$state))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$state} 
    else {temp_org <- paste(temp_org, affiliation$state, sep = ", ")}
  }
  if (!anyNA(affiliation$state_code)&!is.null(affiliation$state_code))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$state_code} 
    else {temp_org <- paste(temp_org, affiliation$state_code, sep = ", ")}
  }
  
  if (length(temp_org)==0){temp_org <- "Empty value"}
  else {temp_org <- trim_org(temp_org)}
  
  return(temp_org)
}

folder_list <- list.dirs(data_directory, full.names = FALSE, recursive = FALSE) 

for (current_folder in folder_list) 
{
  print(paste('Current folder: ', current_folder))
  
  file_list <-list.files(paste(data_directory,current_folder,sep="")) 
  
  file <- data.frame() 
  
  for(current_file in file_list)
  {
    print(paste("Curent number of Errors: ", as.character(count_errors)))
    print(paste('Current file: ', current_file))
    con = file(paste(data_directory,current_folder,"/",current_file,sep=""), "r")
    
    number_of_record <- 0
    my_block <- 100
  
    while ( TRUE ) 
    {
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }
     
      mod2 <- try({
        paper <- jsonlite::fromJSON(line)
      }, TRUE)
      
      if (isTRUE(class(mod2)=="try-error"))
      {
        count_errors <- count_errors +1   
      }
      else
      {
        if (length(paper$author_affiliations)>0)
        {
          for (i in 1: length(paper$author_affiliations))
          { 
            if ( number_of_record >= my_block ) {
              break
            }
            if ((length(paper$author_affiliations[[i]][[1]]$first_name)<1) & (length(paper$author_affiliations[[i]][[1]]$affiliations)<1)) 
            {
              print("Error!")
              print(length(paper$author_affiliations[[i]][[1]]$first_name))
              print(paper$author_affiliations[[i]][[1]]$first_name)
              print(length(paper$author_affiliations[[i]][[1]]$affiliations))
              print(paper$author_affiliations[[i]][[1]]$affiliations)
            }
            else
            {
              if (length(paper$author_affiliations[[i]][[1]]$first_name)==0)
              {
                for (j in 1:length(paper$author_affiliations[[i]][[1]]$affiliations))
                {
                  if ( number_of_record >= my_block ) {
                    break
                  }
                  temp$year <- paper$year[i]
                  temp$first_name <- "Empty value"
                  if (length(paper$author_affiliations[[i]][[1]]$affiliations[j])==0)
                  {
                    temp$org <- "Empty value"
                    temp$country <- NULL 
                  }
                  else 
                  {
                    list_returned <- get_from_affiliation(paper$author_affiliations[[i]][[1]]$affiliations[[j]])
                    temp$org <- list_returned$org
                    temp$country <- list_returned$country
                  }                  
                  file <- rbind(file, c(temp$year, temp$first_name, temp$org, temp$country)) 
                  number_of_record <- number_of_record +1
                }
              }
              else
              {
                if (length(paper$author_affiliations[[i]][[1]]$affiliations)==0)
                {
                  for (j in 1:length(paper$author_affiliations[[i]][[1]]$first_name))
                  {
                    if ( number_of_record >= my_block ) {
                      break
                    }
                    temp$year <- paper$year[i]
                    temp$org <- "Empty value"
                    temp$country <- NULL 
                    
                    if (length(paper$author_affiliations[[i]][[1]]$first_name[[j]])==0)
                    {
                      temp$first_name <- "Empty value"
                    }
                    else 
                    {
                      temp$first_name <- paper$author_affiliations[[i]][[1]]$first_name[[j]]
                    }
                    file <- rbind(file, c(temp$year, temp$first_name, temp$org, temp$country)) 
                    number_of_record <- number_of_record +1
                  }
                }
                else
                {
                  {
                    min_count <- min(length(paper$author_affiliations[[i]][[1]]$first_name),length(paper$author_affiliations[[i]][[1]]$affiliations))
                    max_count <- max(length(paper$author_affiliations[[i]][[1]]$first_name),length(paper$author_affiliations[[i]][[1]]$affiliations))
                    
                    for (j in 1:min_count)
                    {
                      if ( number_of_record >= my_block ) {
                        break
                      }
                      temp$year <- paper$year[i]
                      if (length(paper$author_affiliations[[i]][[1]]$first_name[[j]])==0)
                      {
                        temp$first_name <- "Empty value"
                      }
                      else 
                      {
                        temp$first_name <- paper$author_affiliations[[i]][[1]]$first_name[[j]]
                      }
                      if (length(paper$author_affiliations[[i]][[1]]$affiliations[j])==0)
                      {
                        temp$org <- "Empty value"
                        temp$country <- NULL 
                      }
                      else 
                      {
                        list_returned <- get_from_affiliation(paper$author_affiliations[[i]][[1]]$affiliations[[j]])
                        temp$org <- list_returned$org
                        temp$country <- list_returned$country
                      }
                      file <- rbind(file, c(temp$year, temp$first_name, temp$org, temp$country)) 
                      number_of_record <- number_of_record +1
                    }
                    
                    if (length(paper$author_affiliations[[i]][[1]]$first_name)>length(paper$author_affiliations[[i]][[1]]$affiliations))
                    {
                      for (j in (min_count+1):max_count)
                      {
                        if ( number_of_record >= my_block ) {
                          break
                        }
                        temp$year <- paper$year[i]
                        temp$first_name <- paper$author_affiliations[[i]][[1]]$first_name[[j]]
                        temp$org <- "Empty value"
                        temp$country <- NULL 
                        
                        file <- rbind(file, c(temp$year, temp$first_name, temp$org, temp$country)) 
                        number_of_record <- number_of_record +1
                      }
                    }
                    if (length(paper$author_affiliations[[i]][[1]]$first_name)<length(paper$author_affiliations[[i]][[1]]$affiliations))
                    {
                      for (j in (min_count+1):max_count)
                      {
                        if ( number_of_record >= my_block ) {
                          break
                        }
                        temp$year <- paper$year[i]
                        temp$first_name <- "Empty value"
                        if (length(paper$author_affiliations[[i]][[1]]$affiliations[j])==0)
                        {
                          temp$org <- "Empty value"
                          temp$country <- NULL 
                        }
                        else 
                        {
                          list_returned <- get_from_affiliation(paper$author_affiliations[[i]][[1]]$affiliations[[j]])
                          temp$org <- list_returned$org
                          temp$country <- list_returned$country
                        }                        
                        file <- rbind(file, c(temp$year, temp$first_name, temp$org, temp$country)) 
                        number_of_record <- number_of_record +1
                      }
                    }    
                  }
                }
              }
            }
          }
        }      
      }
  
    }
    
    close(con)
    
    file_name <- paste("DQDs - Listed relevent attributes - Dimensions - ",current_folder,".csv", sep="") 
    write.table(file, file = file_name, 
                sep=",",append = T,row.names = F, col.names = F)   
    
  }

}