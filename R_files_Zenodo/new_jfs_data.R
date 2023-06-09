# Create a new sfs_data_table 

library(dplyr)
library(tidyr)

#setwd("../Desktop/February_2021/sfs_data/")
#indv_data <- read.table("all_individual_blocks_012121b.txt", header=FALSE, sep="\\t")
#indv_data_test <- indv_data[1:300,]


args = commandArgs(trailingOnly=TRUE)

indv_data_initial <- read.table(args[1], header=FALSE, sep="\\t")
indv_data <- unique(indv_data_initial)


blocks <- separate_rows(indv_data,V8,sep=" ", convert = TRUE)
colnames(blocks) <- c("INFO","GEN","MIG","DEME","cluster","process","indv","V8")

blocks_to_add1 <- blocks %>% filter(grepl(",0,2",V8))
blocks_to_add2 <- blocks %>% filter(grepl(",2,0",V8))

blocks$V8 <- gsub(",0,2",",0,1",blocks$V8)
blocks_to_add1$V8 <- gsub(",0,2",",0,1",blocks_to_add1$V8)
blocks_to_add2$V8 <- gsub(",0,2",",0,1",blocks_to_add2$V8)

blocks$V8 <- gsub(",2,0",",1,0",blocks$V8)
blocks_to_add1$V8 <- gsub(",2,0",",1,0",blocks_to_add1$V8)
blocks_to_add2$V8 <- gsub(",2,0",",1,0",blocks_to_add2$V8)

blocks_table1 <- rbind(blocks,blocks_to_add1,blocks_to_add2)

blocks_table2 <- blocks_table1 %>% select(-INFO,-indv) %>% filter(!grepl("NA",V8)) %>% group_by(GEN,MIG,DEME,cluster,process,V8) %>% tally()

blocks_table3 <- blocks_table2 %>% group_by(GEN,MIG,DEME,cluster,process,n) %>% tally()

blocks_table4 <- blocks_table3 %>% group_by(GEN,MIG,DEME,n) %>% summarise(freq_total=sum(nn)) 

colnames(blocks_table4)[4] <- "category"

blocks_table5 <- blocks_table4 %>% add_tally(freq_total)

blocks_table5$prop <- blocks_table5$freq_total/blocks_table5$n

colnames(blocks_table5) <- c("GEN","MIG","DEME","freq_category","freq_total","total_junctions","prop")

save(blocks_table5,file="sfs_data.Rdata")
#load("sfs_data.Rdata")


#blocks_table2 <- separate(blocks_table3,V8,c("SNP","From","To"),sep=",") 
#blocks_table <- blocks_table2 %>% filter(From!="NA") %>% filter(To!="NA")
