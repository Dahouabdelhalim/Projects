library(rioja)
library(vegan)
library(dplyr)
library(ggplot2)

#Pollen data processing, numerical analyses and diagram plotting
##################################################################
##################################################################
#Creating tables and readable datasets for pollen
pollen<-read.table("pollen.csv", header=TRUE, sep=",", stringsAsFactors
                   = FALSE)
dictionary<-read.table("dictionary.csv", header=TRUE, sep=",",
                       stringsAsFactors = FALSE)
#defining time, depth and spike objects
depth<-pollen [, "depth.cm"]
age<-pollen [, "AgeCal"]
lycopodium<-pollen [,"Lycopodium"]
undet <- pollen[,"Undet"]
#percentages for pollen sum taxa
################################################
pollen.types<-dictionary[(dictionary$functional!="fern" &
                            dictionary$functional!="hydrophyte" & 
                            dictionary$functional!= "hygrophyte"),"variables"]
pollen.sum<-pollen[, pollen.types]
percentages.pollen.sum<-pollen.sum/rowSums(pollen.sum) * 100#percentages for non upland, non terrestrial pollen sum taxa
##########################################################
non.pollen.types <-dictionary[dictionary$functional=="fern" |
                                dictionary$functional=="hydrophyte" |
                                dictionary$functional=="hygrophyte","variables"]
non.pollen.sum <-pollen[, non.pollen.types]
#percentages of these taxa are based on the whole pollen and spores
assemblage
#unlike the pollen sum taxa, whose percentages are based on just
#terrestrial pollen (woody and herbs)
pollen.total <- dictionary[,"variables"]
pollen.total.sum <-pollen [,pollen.total]
percentages.non.pollen.sum <-non.pollen.sum/rowSums(pollen.total.sum) *
  100
##percentages for functional groups
################################################
#getting the functional groups
functional.groups<-unique(dictionary$functional)
#functional groups non fern or aquatic
functional.groups.non.aquatic <-functional.groups[functional.groups!="fern"
                                                & functional.groups!="hydrophyte"
                                                & functional.groups
                                                !="hygrophyte" & functional.groups!="pine"]
                                              
#create empty dataframe for the non aquatic functional groups
functional.groups.non.aquatic.df<-data.frame(matrix(NA,
                                                    ncol=length(functional.groups.non.aquatic), nrow=nrow(pollen)))
names(functional.groups.non.aquatic.df)<-functional.groups.non.aquatic
for (functional.group in functional.groups.non.aquatic){
  #selecting the columns in pollen belonging to functional.group
  pollen.functional.group<-pollen[,
                                  dictionary[dictionary$functional==functional.group ,"variables"]]
  #summing the rows (if more than one column, we use rowSums, otherwise
  #we get the column itself)
if(ncol(pollen.functional.group)>1){
  functional.groups.non.aquatic.df[,functional.group]<-
    rowSums(pollen.functional.group)
} else {
  functional.groups.non.aquatic.df[,functional.group]<-
    pollen.functional.group
}
}
#computing percentage functional groups non aquatic
functional.groups.non.aquatic.df<-functional.groups.non.aquatic.df/
  rowSums(functional.groups.non.aquatic.df) * 100
#functional groups for fern or aquatic
functional.groups.aquatic <-functional.groups[functional.groups
                                              !="woody" & functional.groups!="herbs"]
#create empty dataframe for the aquatic functional groups
functional.groups.aquatic.df<-data.frame(matrix(NA,
                                                ncol=length(functional.groups.aquatic), nrow=nrow(pollen)))
names(functional.groups.aquatic.df)<-functional.groups.aquatic
for (functional.group in functional.groups.aquatic){
  #selecting the columns in pollen belonging to functional.group
  pollen.functional.group<-pollen[,
                                  dictionary[dictionary$functional==functional.group ,"variables"]]
  #summing the rows (if more than one column, we use rowSums, otherwise we get the column itself)
if(ncol(pollen.functional.group)>1){
  functional.groups.aquatic.df[,functional.group]<-
    rowSums(pollen.functional.group)
} else {
  functional.groups.aquatic.df[,functional.group]<-
    pollen.functional.group
}
}
#computing percentage functional groups aquatic
functional.groups.aquatic.df<-functional.groups.aquatic.df/
  rowSums(functional.groups.aquatic.df) * 100
#Binding all objects and data frames together
##############################################
pollen_p<-cbind(depth, age, percentages.pollen.sum,percentages.non.pollen.sum,
functional.groups.non.aquatic.df,functional.groups.aquatic.df)
pollen_raw_p <-cbind(pollen, percentages.pollen.sum,
                     percentages.non.pollen.sum,
                     functional.groups.non.aquatic.df,functional.groups.aquatic.df)
#writing result
save(pollen_p, file = "pollen_p.RData")
save(pollen_raw_p, file = "pollen.RData")
write.table(pollen_raw_p, file = "pollen_raw_p.csv", row.names = FALSE,
            col.names = TRUE, sep=",")
write.table(pollen_p, file = "pollen_p.csv", row.names = FALSE,
            col.names = TRUE, sep=",")

##Plotting pollen diagram
#######################################################################

age <- pollen_raw_p$AgeCal
depth <-pollen_raw_p$depth.cm
pollen.p<-read.table("pollen_p.csv", header=TRUE, sep=",",
                     stringsAsFactors = FALSE)
pollen.plot <- pollen.p[-1:-2]
write.table(pollen.plot, file="pollen.plot.csv", sep= ",")
## Change any NA it may have changed from 0 back to 0
pollen.plot[is.na(pollen.plot)] = 0
##Filtering those with <2% 
abundancemx <- apply(pollen.plot, 2, max)
# Remove taxa where total abundance isless than 2% try both
pollen.plot.p2 <- pollen.plot[, mx>2]
##Setting the axis limits
ylim <- range(0,135000)
ymin <-0
ymax <-135000
y.tks <-seq(ymin, 135000,5000)
##Plot diagram
dev.off()
#partial diagrams by making new data frames for each plot
names(pollen.plot.p2)
all.part.1 <-pollen.plot.p2 [1:18]
all.part.2 <- pollen.plot.p2 [19:36]
all.part.3 <- pollen.plot.p2 [37:54]
all.part.4 <- pollen.plot.p2 [55:66]
names(all.part.1)
names(all.part.2)
names(all.part.3)
names(all.part.4)
all.1 <-strat.plot(d=all.part.1, yvar=age, y.rev=TRUE,
                   ylim=ylim, y.tks=y.tks, plot.poly=TRUE, plot.line =
                     TRUE, col.poly ="forestgreen", col.line ="forestgreen",
                   title="Villarquemado (1050m asl)")all.2 <- strat.plot(d=all.part.2, yvar=age, y.rev=TRUE,
                                                                         ylim=ylim, y.tks=y.tks, plot.poly=TRUE, plot.line =
                                                                           TRUE, col.poly="forestgreen", col.line = "forestgreen",
                                                                         title="Villarquemado (1050m asl)")
all.3 <- strat.plot(d=all.part.3, yvar=age, y.rev=TRUE,
                    ylim=ylim, y.tks=y.tks, plot.poly=TRUE, plot.line =
                      TRUE, col.poly="forestgreen", col.line = "forestgreen",
                    title="Villarquemado (1050m asl)")
all.4 <- strat.plot(d=all.part.4, yvar=age, y.rev=TRUE,
                    ylim=ylim, y.tks=y.tks, plot.poly=TRUE, plot.line =
                      TRUE, col.poly="forestgreen", col.line = "forestgreen",
                    title="Villarquemado (1050m asl)")