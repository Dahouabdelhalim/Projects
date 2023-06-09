table1 <- read.csv("./table1_all_reports/all_parks.csv", stringsAsFactors = F)
table1[,"Invasives"] <- c("")

inv_list <- read.csv("https://raw.githubusercontent.com/KelseyDCooper/USGS-NPS-App/master/invasives_list.csv")
inv_list <- inv_list[,"Scientific.Name"]
inv_list <- as.character(inv_list)
inv_list2 <- read.csv("./table1_all_reports/InvasivePlantsList.csv", stringsAsFactors = F)
inv_list2 <- inv_list2[,"ScientificName2"]
inv_list2 <- as.character(inv_list2)

inv_list3 <- c(inv_list, inv_list2, "Plantago coronopus", "Senna artemisioides", "Soliva sessilis", "Veronica persica")


for (i in 1:length(table1$ï..Scientific_name)){
  if (table1$ï..Scientific_name[i] %in% inv_list3){
    table1$Invasives[i] <- c("Yes")
  } else {
    table1$Invasives[i] <- c("No")
  }
}

write.csv(table1, file="20210128_invasives_all_RG_entries.csv")
