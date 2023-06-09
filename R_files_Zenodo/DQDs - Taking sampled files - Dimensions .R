index_list <- c(1,3,4,1,1,5,1,2,4,1,11,1,1,1,1,1,1,1,1,1,1,1)

data_directory <- "D:\\\\Dimensions" 
folder_list <- list.files(data_directory,recursive = FALSE)

new_folder <- "D:\\\\Dimensions\\\\sampled\\\\"

for (current_folder in folder_list)
{
  
  file_list <-list.files(paste(data_directory,current_folder,sep=""))
  
  dataset_index <- c(1:length(file_list))
  current_folder_index <- as.numeric(current_folder)
  my_sample_iteration <- index_list[current_folder_index]
  sample_index <- sample(dataset_index, my_sample_iteration*50, FALSE)
  for (i in 1:length(sample_index))
  {
    file.copy(paste(data_directory,current_folder,"\\\\",file_list[sample_index[i]], sep=""), new_folder)
  }
}
