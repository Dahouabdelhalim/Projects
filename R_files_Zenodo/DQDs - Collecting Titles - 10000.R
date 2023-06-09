file_name <- "D:\\\\Titles-10000.csv"

my_block <- 1000
my_block_file <- 10
my_sample_iteration <- 10

data_directory <- "D:\\\\Dimensions\\\\" 

folder_list <- list.dirs(data_directory, full.names = FALSE, recursive = FALSE) 

for (current_folder in folder_list) 
{
  file_list <-list.files(paste(data_directory,current_folder,sep=""))
  if (length(file_list)>0)
  {
    print(paste('CURRENT FOLDER: ', current_folder))

    block_count <- length(file_list)%/%my_block_file 
    
    if (block_count<my_sample_iteration) {block_count<- length(file_list)-my_block_file} 
    
    dataset_index <- c(1:block_count)
    sample_index <- sample(dataset_index, my_sample_iteration, FALSE)
    
    print(sample_index)
    
    file <- data.frame() 
    
    for (index_count in 1:length(sample_index))
    {
      print(paste("Sample number: ", index_count, ", starting from the file: ", file_list[sample_index[index_count]]))
      print(paste('Starting time: ', Sys.time()))
      
      number_of_record <- 0 
      
      for(number_file_count in 1:my_block_file) 
      {
        if ( number_of_record >= my_block ) {
          break
        }
        
        current_file <- file_list[(sample_index[index_count]-1)*my_block_file+number_file_count] 
        if (block_count == length(file_list)-my_block_file){current_file <- file_list[sample_index[index_count]]} 
        
        print(paste('Current file: ', current_file))
        con = file(paste(data_directory,current_folder,"\\\\",current_file,sep=""), "r")
        
        while ( TRUE ) 
        {
          if ( number_of_record >= my_block ) {
            break
          }
          
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
            if (!is.null(nrow(paper)))
            {
              if (nrow(paper)>0)
              {
                for (i in 1:nrow(paper))
                {
                  file <- rbind(file, c(current_folder,paper$title[i], paper$year[i]))  
                  number_of_record <- number_of_record +1
                  if ( number_of_record >= my_block ) {
                    break
                  }
                }
              }
            }
            
          }
        }
        
        close(con)
        
        print(paste('Ending time: ', Sys.time()))
      }
      
    }  
    
    write.table(file, file = file_name, 
                sep=",",append = T,row.names = F, col.names = F)   
  }  
}
