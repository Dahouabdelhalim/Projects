#This script converts the networks of nests into measures of nest properties for each nest

library("igraph")
library("moments")
library("statnet.common")


##Initialize empty lists for data to be written into
cycles.list<-c()  #make an empty list for number of cycles in a network
mean.dist.all<-c()  #Mean distance for whole nest
num.nodes<-c()

filename.list<-c()   #make an empty list for the file names


# Run analyses:   ----
#Nests starting with E1 - entrance node is a hole that starts as a tunnel
#E1----
setwd("....") #set working directory
input_files<-list.files("....Data Files\\\\E1\\\\") #Make a list of the names of all the files in that folder, now input_files is a list of csv's
data_folder<-("....Data Files\\\\E1\\\\")


for (i in 1:length(input_files)){ #For each csv file in the folder...
  
  
  
  #Prepare data----
  #_____________________
  data<-paste(data_folder,input_files[i],sep="") #Select the ith file in the Data Files folder; this yields a path, not the actual file, data = path
  df<-read.csv(data, header=FALSE) #read this file and convert into a data frame
  subset.data<-df[c(3,4)] #in Ant Nest Excel Sheet, the 2-column matrix is in colmns 3 & 4
  edgelist<-as.matrix(subset.data) #convert data to a matrix - 2 columns
  graph<-graph_from_edgelist(edgelist,directed=FALSE) #make an igraph object
  plot(graph,main=c(input_files[i]))

  
  #Cycles----
  #__________________
  n<-vcount(graph) #number of nodes or vertices
  m<-ecount(graph) #number of edges
  num.cycles<-m-n+1
  cycles.list<-c(cycles.list,num.cycles)
  num.nodes<-c(num.nodes,n)
  
  
  #Mean Distance ----
  #_______________________________________
  #Mean Distance for whole nest
  mean.dist.all<-c(mean.dist.all,mean_distance(graph))
  
  
  filename.list<-c(filename.list,input_files[i])
  
}

# Run analyses: 
#Nests starting with C1 - entrance chamber is when a hole that leads to a chamber
#C1----
setwd("....") #set working directory
input_files<-list.files("....\\\\Data Files\\\\C1\\\\") #Make a list of the names of all the files in that folder, now input_files is a list of csv's
data_folder<-("....\\\\Data Files\\\\C1\\\\")

for (i in 1:length(input_files)){ #For each csv file in the folder...
  
  #Prepare data ----
  #_______________
  data<-paste(data_folder,input_files[i],sep="") #Select the ith file in the Data Files folder; this yields a path, not the actual file, data = path
  df<-read.csv(data, header=FALSE) #read this file and convert into a data frame
  subset.data<-df[c(3,4)] #in Ant Nest Excel Sheet, the 2-column matrix is in colmns 3 & 4
  edgelist<-as.matrix(subset.data) #convert data to a matrix - 2 columns
  graph<-graph_from_edgelist(edgelist,directed=FALSE) #make an igraph object
  plot(graph)
 
  #Cycles----
  #__________________
  n<-vcount(graph) #number of nodes or vertices
  m<-ecount(graph) #number of edges
  num.cycles<-m-n+1
  cycles.list<-c(cycles.list,num.cycles)
  num.nodes<-c(num.nodes,n)
  
  #Mean Distance ----
  #_______________________________________
  #Mean Distance for whole nest
  mean.dist.all<-c(mean.dist.all,mean_distance(graph))
  
  filename.list<-c(filename.list,input_files[i])
  
}


#MULTIPLE ENTRIES----
# Run analyses: Nests with multiple entries
setwd("....") #set working directory
input_files<-list.files("....\\\\Data Files\\\\Multiple Entries\\\\") #Make a list of the names of all the files in that folder, now input_files is a list of csv's
data_folder<-("....\\\\Data Files\\\\Multiple Entries\\\\")

for (i in 1:length(input_files)){ #For each csv file in the folder...
  
  #Prepare data ----
  #_______________
  data<-paste(data_folder,input_files[i],sep="") #Select the ith file in the Data Files folder; this yields a path, not the actual file, data = path
  df<-read.csv(data, header=FALSE) #read this file and convert into a data frame
  subset.data<-df[c(3,4)] #in Ant Nest Excel Sheet, the 2-column matrix is in colmns 3 & 4
  edgelist<-as.matrix(subset.data) #convert data to a matrix - 2 columns
  graph<-graph_from_edgelist(edgelist,directed=FALSE) #make an igraph object
  plot(graph)
  

    #Cycles----
    #__________________
    n<-vcount(graph) #number of nodes or vertices
    m<-ecount(graph) #number of edges
    num.cycles<-m-n+1
    cycles.list<-c(cycles.list,num.cycles)
    num.nodes<-c(num.nodes,n)
    
    
    #Mean Distance ----
    #_______________________________________
    #Mean Distance for whole nest
    mean.dist.all<-c(mean.dist.all,mean_distance(graph))
    
    
    filename.list<-c(filename.list,input_files[i])

}



#AFTER the above for loops have run with each file type, E1, C1 and multiple entries, run this
result<-data.frame(filename.list,cycles.list, mean.dist.all,num.nodes)

setwd("....")
write.csv(result,"Network_Analyses_WholeNest_2021.csv")


#Add species name to results spreadsheet
NameKey<-read.csv("....filename to species name key.csv")
TraitData<-merge(result,NameKey,by="filename.list",all.x=TRUE)
write.csv(TraitData,"Trait_Data_WholeNest_2021.csv")

TraitData<-read.csv("Trait_Data_WholeNest_2021.csv")

#Combine network data with colony size data
ColSize<-read.csv("....Colony size.csv")
dall<-merge(TraitData,ColSize,by="Name",all.x=TRUE)

#Write merged data to new csv file
setwd("....")
write.csv(dall,"dall_WholeNest.csv")

#Now use 'dall' in "PGLS with caper.r" script file

