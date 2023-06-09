FOS_name_list <- read.csv("FOS_name_list.csv", header = TRUE, sep=",")
FOS_count_list <- read.csv("FOS_count_list.csv", header = TRUE, sep=",")
length(FOS_name_list)

for (i in 1:(length(FOS_name_list)-1))
{
  for (j in (i+1):length(FOS_name_list))
  {
    for (k in 1:length(FOS_name_list[,i]))
    {
      for (l in 1:length(FOS_name_list[,j]))
      {
        if (FOS_name_list[k,i]==FOS_name_list[l,j])
        {
          print(c(i,k,j,l,FOS_name_list[k,i],FOS_count_list[k,i], FOS_name_list[l,j], FOS_count_list[l,j]))
          if (k>l) 
          {
          FOS_name_list[k,i] <- "" 
          }
          else 
          {
          FOS_name_list[l,j] <-""
          }
        }
      }
    }
  }
}

write.csv(FOS_name_list, file = "FOS_name_list_updated.csv", row.names = FALSE)

