#load the required library 
library(tidyverse)

#Z-transform FST values for each population pair, select the 50 highest ZFST values, 
#and print the average of these 50 ZFST estimates for each population pair

#ZFST for Florida population pairs
max_ZFSTs_Florida <- tibble(population_pair=character(), max_ZFST=numeric())
files<-read.csv(file="unique_Florida_population_combinations.txt", header = FALSE)
poplist<-dplyr::pull(files, V1)

for (i in poplist){
  name<-paste ("pairwise_FST_unique_pop_pairs/",i,".windowed.weir.fst",sep="")
  input<-read.csv(file=name, header = TRUE, sep='\\t')
  input$Z_WEIGHTED_FST<-(input$WEIGHTED_FST-mean(input$WEIGHTED_FST))/sd(input$WEIGHTED_FST)
  top50<-top_n(input, 50, Z_WEIGHTED_FST)
  top50_zfsts<-mean(top50$Z_WEIGHTED_FST)
  max_ZFSTs_Florida %<>% summarise(population_pair = i, max_ZFST = top50_zfsts) %>% bind_rows(max_ZFSTs_Florida)
}
write.table(max_ZFSTs_Florida, file="Florida_max_ZFSTs.csv",row.names = FALSE, col.names = FALSE, quote = FALSE)


#ZFST for Cuba population pairs
max_ZFSTs_Cuba <- tibble(population_pair=character(), max_ZFST=numeric())
files<-read.csv(file="unique_Cuba_population_combinations.txt", header = FALSE)
poplist<-dplyr::pull(files, V1)

for (i in poplist){
  name<-paste ("pairwise_FST_unique_pop_pairs/",i,".windowed.weir.fst",sep="")
  input<-read.csv(file=name, header = TRUE, sep='\\t')
  input$Z_WEIGHTED_FST<-(input$WEIGHTED_FST-mean(input$WEIGHTED_FST))/sd(input$WEIGHTED_FST)
  top50<-top_n(input, 50, Z_WEIGHTED_FST)
  top50_zfsts<-mean(top50$Z_WEIGHTED_FST)
  max_ZFSTs_Cuba %<>% summarise(population_pair = i, max_ZFST = top50_zfsts) %>% bind_rows(max_ZFSTs_Cuba)
}
write.table(max_ZFSTs_Cuba, file="Cuba_max_ZFSTs.csv",row.names = FALSE, col.names = FALSE, quote = FALSE)
