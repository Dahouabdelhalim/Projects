#Relatedness script
#Testing whether relatedness among dyads within groups and commmunities is more than expected by chance


#Load relatedness data
rel = read.csv(here::here("Input files","RelatednessEstimates_8xdenovo80_14to19.csv"))
rel = rel %>% select(-Pair)

#Load network data
birdlist16noKsubc = read.csv(here::here("Output files","birdlist16noKsubc.csv"))
birdlist17noKsubc = read.csv(here::here("Output files","birdlist17noKsubc.csv"))
birdlist18noKsubc = read.csv(here::here("Output files","birdlist18noKsubc.csv"))
birdlist19noKsubc = read.csv(here::here("Output files","birdlist19noKsubc.csv"))

library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(plotrix)

####Relatedness within groups yearly function####
related_groups = function(rel, birdlist) {
  
  #Get all possible dyads
  dyads = expand.grid(birdlist$FWnumber,birdlist$FWnumber)
  dyads = dyads[which(dyads$Var1!=dyads$Var2),]
  dyads = data.frame(t(apply(dyads, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
  colnames(dyads) = c("Bird1","Bird2")
  dyads = dyads %>% distinct(Bird1,Bird2) #Remove duplicate pairs
  
  #Merge relatedness
  dyads = merge(dyads,rel,by.x=c("Bird1","Bird2"),by.y=c("ind1","ind2"))
  
  #Bring birdlist info
  birdlist = birdlist %>% select(FWnumber, Social.Group,Community,Sex,class)
  dyads = merge(dyads,birdlist,by.x="Bird1",by.y="FWnumber")
  dyads = merge(dyads,birdlist,by.x="Bird2",by.y="FWnumber")
  
  #Grab dyads that are in the same social group
  dyads.same = dyads[which(dyads$Social.Group.x==dyads$Social.Group.y),]
  
  #Mean relatedness within social groups
  mean.relatedness = mean(dyads.same$wang)
  
  #Look at relatedness among social classes
  classes = c("Resident Male","Resident Female","Auxilliary Male","1 Year-old Male","Dispersing Female")
  classes = expand.grid(classes,classes)
  classes = data.frame(t(apply(classes, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
  colnames(classes) = c("class.x","class.y")
  classes = classes %>% distinct(class.x,class.y)
  classes$wang.mean = NA
  classes$wang.min = NA
  classes$wang.max = NA
  
  #Order classes within rows so they match the class list order within rows later 
  #This now mismatches to the bird IDs but don't need them anymore, so just grab relatedness values
  dyads.same.classes = dyads.same %>% select(class.x,class.y)
  dyads.same.classes = data.frame(t(apply(dyads.same.classes, 1, sort)),dyads.same$wang)
  colnames(dyads.same.classes) = c("class.x","class.y","wang") 
  
  #Calculate the relatedness for each dyad of social classes 
    for (i in 1:nrow(classes)) {
    cl = dyads.same.classes[which(dyads.same.classes$class.x %in% classes$class.x[i] & 
                                       dyads.same.classes$class.y %in% classes$class.y[i]),]
    if(nrow(cl)==0) {
      classes$wang.mean[i] = NA
      classes$wang.min[i] = NA
      classes$wang.max[i] = NA
    } else {
      classes$wang.mean[i] = mean(cl$wang)
      classes$wang.min[i] = min(cl$wang)
      classes$wang.max[i] = max(cl$wang) }
  }
  
  #Get counts of number of dyads
  classes$count = NA
  for (i in 1:nrow(classes)) {
    cl = dyads.same.classes[which(dyads.same.classes$class.x %in% classes$class.x[i] & 
                                    dyads.same.classes$class.y %in% classes$class.y[i]),]
    classes$count[i] = nrow(cl)
  }
  
  #Get arrays
  group.mean.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.mean",fun.aggregate = sum,fill=NA_real_))
  group.min.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.min",fun.aggregate = sum,fill=NA_real_))
  group.max.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.max",fun.aggregate = sum,fill=NA_real_))
  group.count.array = as.matrix(acast(classes,class.x~class.y,value.var="count",fun.aggregate = sum,fill=NA_real_))
  
  #Print arrays
  print(mean.relatedness)
  print(group.mean.array)
  print(group.min.array)
  print(group.max.array)
  print(group.count.array)
  
  #Don't get multi-level labels for each factor here - need dyads.same.classes output to match single level factors in other lists later
  #Not going to publish individual year graphs so that doens't matter here. 
  
  #Get plots of the arrays
  dyads.same.classes = dyads.same.classes[which(dyads.same.classes$class.y!="Unknown"),]
  plot = ggplot(dyads.same.classes,aes(x=wang)) + geom_histogram(binwidth = 0.02) + facet_grid(class.x~class.y) + theme_bw() +
    theme(strip.background = element_rect(fill="white"),strip.text = element_text(size=12)) + 
    scale_x_continuous(limits=c(-0.25,1),breaks=c(-0.25,0,0.25,0.5,0.75,1)) + xlab("Relatedness") + ylab("Count") 
  
  print(plot)
  return(dyads.same.classes)
}

#Run for each year
rel_groups16 = related_groups(rel=rel, birdlist = birdlist16noKsubc)
rel_groups17 = related_groups(rel=rel, birdlist = birdlist17noKsubc)
rel_groups18 = related_groups(rel=rel, birdlist = birdlist18noKsubc)
rel_groups19 = related_groups(rel=rel, birdlist = birdlist19noKsubc)





####Add all years together and get relatedness within groups plot####
rel_groupsall = rbind(rel_groups16,rel_groups17,rel_groups18,rel_groups19)

related_groups_all = function(dyads.same.classes) {
  #Look at relatedness among social classes
  classes = c("Resident Male","Resident Female","Auxilliary Male","1 Year-old Male","Dispersing Female")
  classes = expand.grid(classes,classes)
  classes = data.frame(t(apply(classes, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
  colnames(classes) = c("class.x","class.y")
  classes = classes %>% distinct(class.x,class.y)
  classes$wang.mean = NA
  classes$wang.min = NA
  classes$wang.max = NA
  classes$std.error = NA
  
  #Calculate the relatedness for each dyad of social classes 
  for (i in 1:nrow(classes)) {
    cl = dyads.same.classes[which(dyads.same.classes$class.x %in% classes$class.x[i] & 
                                    dyads.same.classes$class.y %in% classes$class.y[i]),]
    if(nrow(cl)==0) {
      classes$wang.mean[i] = NA
      classes$wang.min[i] = NA
      classes$wang.max[i] = NA
      classes$std.error[i] = NA
    } else {
      classes$wang.mean[i] = mean(cl$wang)
      classes$wang.min[i] = min(cl$wang)
      classes$wang.max[i] = max(cl$wang) 
      classes$std.error[i] = std.error(cl$wang)}
  }
  
  #Get counts of number of dyads
  classes$count = NA
  for (i in 1:nrow(classes)) {
    cl = dyads.same.classes[which(dyads.same.classes$class.x %in% classes$class.x[i] & 
                                    dyads.same.classes$class.y %in% classes$class.y[i]),]
    classes$count[i] = nrow(cl)
  }
  
  #Put into array
  group.mean.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.mean",fun.aggregate = sum,fill=NA_real_))
  group.min.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.min",fun.aggregate = sum,fill=NA_real_))
  group.max.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.max",fun.aggregate = sum,fill=NA_real_))
  group.count.array = as.matrix(acast(classes,class.x~class.y,value.var="count",fun.aggregate = sum,fill=NA_real_))
  group.se.array = as.matrix(acast(classes,class.x~class.y,value.var="std.error",fun.aggregate = sum,fill=NA_real_))
  
  #Print arrays
  print(group.mean.array)
  print(group.min.array)
  print(group.max.array)
  print(group.count.array)
  print(group.se.array)
  
  #Get multi-level labels for each factor
  levels(dyads.same.classes$class.x) <- gsub(" ", "\\n", levels(dyads.same.classes$class.x))
  levels(dyads.same.classes$class.y) <- gsub(" ", "\\n", levels(dyads.same.classes$class.y))
  
  #Plot
  plot = ggplot(dyads.same.classes,aes(x=wang)) + geom_histogram(binwidth = 0.04) + facet_grid(class.x~class.y) + theme_bw() +
    theme(strip.background = element_rect(fill="white"),strip.text = element_text(size=12)) + 
    scale_x_continuous(limits=c(-0.25,1),breaks=c(-0.25,0,0.25,0.5,0.75,1)) + xlab("Relatedness") + ylab("Count")  +
    theme(panel.grid.minor = element_blank()) + scale_y_continuous(breaks=c(0,5,10,15),limits=c(0,15)) + 
    theme(panel.spacing = unit(0.5, "lines"))
  print(plot)
  return(group.mean.array)
}

#Get relatedness among dyads for all years combined
groupsall_obs = related_groups_all(dyads.same.classes = rel_groupsall) 


##Get overall mean within group relatedness 
mean(rel_groupsall$wang)
std.error(rel_groupsall$wang)





####Compare observed relatedness within groups to random####

#First a function to get n # of randomized dyad.same.classes's for a year. Only randomize within social classes. 
#Randomly assign a new social group to each individual. 

dyads_same_rand = function(rel, birdlist,permutations) {
  
  #Create an empty list
  dyads_list <- vector(mode = "list", length = permutations)
  
  for (i in 1:permutations) {
    #Order by class
    birdlist.r = birdlist[order(birdlist$class),]
    #Randomly sample within group without replacement - this randomly orders within class, but keeps class order the same
    birdlist.r2 = birdlist.r %>% group_by(class) %>% sample_frac(replace=F) %>% select(class,Social.Group)
    #Then add the random social group to the birdlist dataframe
    birdlist.r$Social.Group = birdlist.r2$Social.Group
  
    #Get all possible dyads
    dyads = expand.grid(birdlist.r$FWnumber,birdlist.r$FWnumber)
    dyads = dyads[which(dyads$Var1!=dyads$Var2),]
    dyads = data.frame(t(apply(dyads, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
    colnames(dyads) = c("Bird1","Bird2")
    dyads = dyads %>% distinct(Bird1,Bird2) #Remove duplicate pairs
    
    #Merge relatedness
    dyads = merge(dyads,rel,by.x=c("Bird1","Bird2"),by.y=c("ind1","ind2"))
    
    #Bring birdlist.r info
    birdlist.r = birdlist.r %>% select(FWnumber, Social.Group,Community,Sex,class)
    dyads = merge(dyads,birdlist.r,by.x="Bird1",by.y="FWnumber")
    dyads = merge(dyads,birdlist.r,by.x="Bird2",by.y="FWnumber")
    
    #Grab dyads that are in the same social group
    dyads.same = dyads[which(dyads$Social.Group.x==dyads$Social.Group.y),] 
    
    #Order classes within rows so they match the class list order within rows later 
    #This now mismatches to the bird IDs but don't need them anymore, so just grab relatedness values
    dyads.same.classes = dyads.same %>% select(class.x,class.y)
    dyads.same.classes = data.frame(t(apply(dyads.same.classes, 1, sort)),dyads.same$wang)
    colnames(dyads.same.classes) = c("class.x","class.y","wang")
    
    #Add to list
    dyads_list[[i]] = dyads.same.classes
  }
  return(dyads_list)
}

#Get random relatedness for dyads for each year - dataframe lengths will be slightly different because not all birds have 
#genetic data so depending on what group they get placed into the number of dyads will be slightly different. 
#dyads_rand16 = dyads_same_rand(rel = rel, birdlist = birdlist16noKsubc, permutations=1000)
#dyads_rand17 = dyads_same_rand(rel = rel, birdlist = birdlist17noKsubc, permutations=1000)
#dyads_rand18 = dyads_same_rand(rel = rel, birdlist = birdlist18noKsubc, permutations=1000)
#dyads_rand19 = dyads_same_rand(rel = rel, birdlist = birdlist19noKsubc, permutations=1000)


#Now a function to combine the random lists and calculate average relatedness values for each set of random lists
dyads_same_rand_combine = function(dyads_rand1,dyads_rand2,dyads_rand3,dyads_rand4,groupsall_obs) {
  
  #Make a list to hold the randomized values from each permutation
  dyads_list.all <- vector(mode = "list", length = length(dyads_rand1))
  
  #Combine the four for every i
  for (i in 1:length(dyads_rand1)) {
    dyads_rand1.i = dyads_rand1[[i]]
    dyads_rand2.i = dyads_rand2[[i]]
    dyads_rand3.i = dyads_rand3[[i]]
    dyads_rand4.i = dyads_rand4[[i]]
    dyads_rand.all = rbind(dyads_rand1.i,dyads_rand2.i,dyads_rand3.i,dyads_rand4.i)
    
    #Code copied from rel_groups_all
    #Look at relatedness among social classes
    classes = c("Resident Male","Resident Female","Auxilliary Male","1 Year-old Male","Dispersing Female")
    classes = expand.grid(classes,classes)
    classes = data.frame(t(apply(classes, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
    colnames(classes) = c("class.x","class.y")
    classes = classes %>% distinct(class.x,class.y)
    classes$wang.mean = NA
    
    #Calculate the relatedness for each dyad of social classes 
    for (h in 1:nrow(classes)) {
      cl = dyads_rand.all[which(dyads_rand.all$class.x %in% classes$class.x[h] & 
                                      dyads_rand.all$class.y %in% classes$class.y[h]),]
      if(nrow(cl)==0) {
        classes$wang.mean[h] = NA
        } else {
        classes$wang.mean[h] = mean(cl$wang)}
    }
    
    #Get array
    group.mean.array.rand = as.matrix(acast(classes,class.x~class.y,value.var="wang.mean",fun.aggregate = sum,fill=NA_real_))
    
    #Put array into list
    dyads_list.all[[i]] = group.mean.array.rand
    }
  
    #Now melt the list of arrays and melt the observed list and combine in a graph
  dyads_list.all.m = melt(dyads_list.all)
  groupsall_obs.m = melt(groupsall_obs)
  
  #Get multi-level labels for each factor
  levels(dyads_list.all.m$Var1) <- gsub(" ", "\\n", levels(dyads_list.all.m$Var1))
  levels(dyads_list.all.m$Var2) <- gsub(" ", "\\n", levels(dyads_list.all.m$Var2))
  levels(groupsall_obs.m$Var1) <- gsub(" ", "\\n", levels(groupsall_obs.m$Var1))
  levels(groupsall_obs.m$Var2) <- gsub(" ", "\\n", levels(groupsall_obs.m$Var2))
  
  plot = ggplot(dyads_list.all.m,aes(x=value)) + geom_histogram(binwidth = 0.02) + facet_grid(Var1~Var2) +
    geom_vline(data=groupsall_obs.m,aes(xintercept=value),color="red",size=0.75) + theme_bw() +
    theme(strip.background = element_rect(fill="white"),strip.text = element_text(size=12)) + 
    scale_x_continuous(limits=c(-0.25,1),breaks=c(-0.25,0,0.25,0.5,0.75,1)) + xlab("Relatedness") + ylab("Count") 
  print(plot)
  
  #Get p-values for each comparison of obeserved to random
  groupsall_obs.m$p.value = NA
  for (i in 1:nrow(groupsall_obs.m)) {
    rand.m = dyads_list.all.m[which(dyads_list.all.m$Var1 %in% groupsall_obs.m$Var1[i] &
                                    dyads_list.all.m$Var2 %in% groupsall_obs.m$Var2[i]),]
    groupsall_obs.m$p.value[i] = 1 - sum(groupsall_obs.m$value[i] > rand.m$value)/length(dyads_rand1)
  }
  
  #Get array of p-values and print
  groupsall_obs.m.array = as.matrix(acast(groupsall_obs.m,Var1~Var2,value.var="p.value",fun.aggregate = sum,fill=NA_real_))
  print(groupsall_obs.m.array)
}

#Run function to compare observed mean relatedness values for each class dyad to random
#dyads_same_rand_combine(dyads_rand16,dyads_rand17,dyads_rand18,dyads_rand19,groupsall_obs)







####Communities####

#Important to note - changed this function so only looking at dyads not in the same social group

related_coms = function(rel, birdlist) {
  
  #Get all possible dyads
  dyads = expand.grid(birdlist$FWnumber,birdlist$FWnumber)
  dyads = dyads[which(dyads$Var1!=dyads$Var2),]
  dyads = data.frame(t(apply(dyads, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
  colnames(dyads) = c("Bird1","Bird2")
  dyads = dyads %>% distinct(Bird1,Bird2) #Remove duplicate pairs
  
  #Merge relatedness
  dyads = merge(dyads,rel,by.x=c("Bird1","Bird2"),by.y=c("ind1","ind2"))
  
  #Bring birdlist info
  birdlist = birdlist %>% select(FWnumber, Social.Group,Community,Sex,class)
  dyads = merge(dyads,birdlist,by.x="Bird1",by.y="FWnumber")
  dyads = merge(dyads,birdlist,by.x="Bird2",by.y="FWnumber")
  
  #Grab dyads that are in the same community
  dyads.same = dyads[which(dyads$Community.x==dyads$Community.y),]
  
  #Get dyads that are in different social groups
  dyads.same = dyads.same[which(dyads.same$Social.Group.x!=dyads.same$Social.Group.y),]
  
  #Mean relatedness within communities
  mean.relatedness = mean(dyads.same$wang)
  
  #Look at relatedness among social classes
  classes = c("Resident Male","Resident Female","Auxilliary Male","1 Year-old Male","Dispersing Female")
  classes = expand.grid(classes,classes)
  classes = data.frame(t(apply(classes, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
  colnames(classes) = c("class.x","class.y")
  classes = classes %>% distinct(class.x,class.y)
  classes$wang.mean = NA
  classes$wang.min = NA
  classes$wang.max = NA
  
  #Order classes within rows so they match the class list order within rows later 
  #This now mismatches to the bird IDs but don't need them anymore, so just grab relatedness values
  dyads.same.classes = dyads.same %>% select(class.x,class.y)
  dyads.same.classes = data.frame(t(apply(dyads.same.classes, 1, sort)),dyads.same$wang)
  colnames(dyads.same.classes) = c("class.x","class.y","wang") 
  
  #Calculate the relatedness for each dyad of social classes 
  for (i in 1:nrow(classes)) {
    cl = dyads.same.classes[which(dyads.same.classes$class.x %in% classes$class.x[i] & 
                                    dyads.same.classes$class.y %in% classes$class.y[i]),]
    if(nrow(cl)==0) {
      classes$wang.mean[i] = NA
      classes$wang.min[i] = NA
      classes$wang.max[i] = NA
    } else {
      classes$wang.mean[i] = mean(cl$wang)
      classes$wang.min[i] = min(cl$wang)
      classes$wang.max[i] = max(cl$wang) }
  }
  
  #Get counts of number of dyads
  classes$count = NA
  for (i in 1:nrow(classes)) {
    cl = dyads.same.classes[which(dyads.same.classes$class.x %in% classes$class.x[i] & 
                                    dyads.same.classes$class.y %in% classes$class.y[i]),]
    classes$count[i] = nrow(cl)
  }
  
  #Get arrays
  group.mean.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.mean",fun.aggregate = sum,fill=NA_real_))
  group.min.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.min",fun.aggregate = sum,fill=NA_real_))
  group.max.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.max",fun.aggregate = sum,fill=NA_real_))
  group.count.array = as.matrix(acast(classes,class.x~class.y,value.var="count",fun.aggregate = sum,fill=NA_real_))
  
  #Print arrays
  print(group.mean.array)
  print(group.min.array)
  print(group.max.array)
  print(group.count.array)
  
  #Get plots of the arrays
  dyads.same.classes = dyads.same.classes[which(dyads.same.classes$class.y!="Unknown"),]
  plot = ggplot(dyads.same.classes,aes(x=wang)) + geom_histogram(binwidth = 0.02) + facet_grid(class.x~class.y) + theme_bw() +
    theme(strip.background = element_rect(fill="white"),strip.text = element_text(size=12)) + 
    scale_x_continuous(limits=c(-0.25,1),breaks=c(-0.25,0,0.25,0.5,0.75,1)) + xlab("Relatedness") + ylab("Count") 
  
  print(plot)
  return(dyads.same.classes)
}

#Run for each year
rel_coms16 = related_coms(rel=rel, birdlist = birdlist16noKsubc)
rel_coms17 = related_coms(rel=rel, birdlist = birdlist17noKsubc)
rel_coms18 = related_coms(rel=rel, birdlist = birdlist18noKsubc)
rel_coms19 = related_coms(rel=rel, birdlist = birdlist19noKsubc)





####Add all years together and get relatedness within coms plot####
rel_comsall = rbind(rel_coms16,rel_coms17,rel_coms18,rel_coms19)

related_coms_all = function(dyads.same.classes) {
  #Look at relatedness among social classes
  classes = c("Resident Male","Resident Female","Auxilliary Male","1 Year-old Male","Dispersing Female")
  classes = expand.grid(classes,classes)
  classes = data.frame(t(apply(classes, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
  colnames(classes) = c("class.x","class.y")
  classes = classes %>% distinct(class.x,class.y)
  classes$wang.mean = NA
  classes$wang.min = NA
  classes$wang.max = NA
  classes$std.error = NA
  
  #Calculate the relatedness for each dyad of social classes 
  for (i in 1:nrow(classes)) {
    cl = dyads.same.classes[which(dyads.same.classes$class.x %in% classes$class.x[i] & 
                                    dyads.same.classes$class.y %in% classes$class.y[i]),]
    if(nrow(cl)==0) {
      classes$wang.mean[i] = NA
      classes$wang.min[i] = NA
      classes$wang.max[i] = NA
    } else {
      classes$wang.mean[i] = mean(cl$wang)
      classes$wang.min[i] = min(cl$wang)
      classes$wang.max[i] = max(cl$wang) 
      classes$std.error[i] = std.error(cl$wang)}
  }
  
  #Get counts of number of dyads
  classes$count = NA
  for (i in 1:nrow(classes)) {
    cl = dyads.same.classes[which(dyads.same.classes$class.x %in% classes$class.x[i] & 
                                    dyads.same.classes$class.y %in% classes$class.y[i]),]
    classes$count[i] = nrow(cl)
  }
  
  #Get arrays
  group.mean.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.mean",fun.aggregate = sum,fill=NA_real_))
  group.min.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.min",fun.aggregate = sum,fill=NA_real_))
  group.max.array = as.matrix(acast(classes,class.x~class.y,value.var="wang.max",fun.aggregate = sum,fill=NA_real_))
  group.count.array = as.matrix(acast(classes,class.x~class.y,value.var="count",fun.aggregate = sum,fill=NA_real_))
  group.se.array = as.matrix(acast(classes,class.x~class.y,value.var="std.error",fun.aggregate = sum,fill=NA_real_))
  
  #Print arrays
  print(group.mean.array)
  print(group.min.array)
  print(group.max.array)
  print(group.count.array)
  print(group.se.array)
  
  #Get multi-level labels for each factor
  levels(dyads.same.classes$class.x) <- gsub(" ", "\\n", levels(dyads.same.classes$class.x))
  levels(dyads.same.classes$class.y) <- gsub(" ", "\\n", levels(dyads.same.classes$class.y))
  
  #Plot
  plot = ggplot(dyads.same.classes,aes(x=wang)) + geom_histogram(binwidth = 0.04) + facet_grid(class.x~class.y) + theme_bw() +
    theme(strip.background = element_rect(fill="white"),strip.text = element_text(size=12)) + 
    scale_x_continuous(limits=c(-0.25,1),breaks=c(-0.25,0,0.25,0.5,0.75,1)) + xlab("Relatedness") + ylab("Count") +
    theme(panel.grid.minor = element_blank()) + #scale_y_continuous(breaks=c(0,5,10,15)) + 
    theme(panel.spacing = unit(0.5, "lines"))
  print(plot)
  return(group.mean.array)
}

comsall_obs = related_coms_all(dyads.same.classes = rel_comsall) 

##Get overall mean within community relatedness - still does not include within-group connections
mean(rel_comsall$wang)
std.error(rel_comsall$wang)




####Compare observed relatedness within communities to random####

#First a function to get n # of randomized dyad.same.classes's for a year. Only randomize within social classes. 
#Randomly assign a new social group to each individual. 
#Do not randomize within social groups here, but do take out dyads in the same social group as above. 
#The randomization just mixes up which community each individual belongs to. I want to essentially work 
#with the same pool of dyads as in the observed test, so don't want to compare two birds that are actually
#in the same social group. If the randomization puts two individuals in the same social group in the same community
#do not include them, as is done above. 

dyads_same_rand_coms = function(rel, birdlist,permutations) {
  
  #Create an empty list
  dyads_list <- vector(mode = "list", length = permutations)
  
  for (i in 1:permutations) {
    #Order by class
    birdlist.r = birdlist[order(birdlist$class),]
    #Randomly sample within group without replacement - this randomly orders within class, but keeps class order the same
    birdlist.r2 = birdlist.r %>% group_by(class) %>% sample_frac(replace=F) %>% select(class,Community)
    #Then add the random community to the birdlist dataframe
    birdlist.r$Community = birdlist.r2$Community

    #Get all possible dyads
    dyads = expand.grid(birdlist.r$FWnumber,birdlist.r$FWnumber)
    dyads = dyads[which(dyads$Var1!=dyads$Var2),]
    dyads = data.frame(t(apply(dyads, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
    colnames(dyads) = c("Bird1","Bird2")
    dyads = dyads %>% distinct(Bird1,Bird2) #Remove duplicate pairs
    
    #Merge relatedness
    dyads = merge(dyads,rel,by.x=c("Bird1","Bird2"),by.y=c("ind1","ind2"))
    
    #Bring birdlist.r info
    birdlist.r = birdlist.r %>% select(FWnumber, Social.Group,Community,Sex,class)
    dyads = merge(dyads,birdlist.r,by.x="Bird1",by.y="FWnumber")
    dyads = merge(dyads,birdlist.r,by.x="Bird2",by.y="FWnumber")
    
    #Grab dyads that are in the same community
    dyads.same = dyads[which(dyads$Community.x==dyads$Community.y),] 
    
    #Get dyads that are in different social groups
    dyads.same = dyads.same[which(dyads.same$Social.Group.x!=dyads.same$Social.Group.y),]
    
    #Order classes within rows so they match the class list order within rows later 
    #This now mismatches to the bird IDs but don't need them anymore, so just grab relatedness values
    dyads.same.classes = dyads.same %>% select(class.x,class.y)
    dyads.same.classes = data.frame(t(apply(dyads.same.classes, 1, sort)),dyads.same$wang)
    colnames(dyads.same.classes) = c("class.x","class.y","wang")
    
    #Add to list
    dyads_list[[i]] = dyads.same.classes
  }
  return(dyads_list)
}

##Get random relatedness for dyads for each year - dataframe lengths will be slightly different because not all birds have 
##genetic data so depending on what group they get placed into the number of dyads will be slightly different. 
#dyads_rand16.coms = dyads_same_rand_coms(rel = rel, birdlist = birdlist16noKsubc, permutations=1000)
#dyads_rand17.coms = dyads_same_rand_coms(rel = rel, birdlist = birdlist17noKsubc, permutations=1000)
#dyads_rand18.coms = dyads_same_rand_coms(rel = rel, birdlist = birdlist18noKsubc, permutations=1000)
#dyads_rand19.coms = dyads_same_rand_coms(rel = rel, birdlist = birdlist19noKsubc, permutations=1000)


#Now a function to combine the random lists and calculate average relatedness values for each set of random lists
dyads_same_rand_combine_coms = function(dyads_rand1,dyads_rand2,dyads_rand3,dyads_rand4,comsall_obs) {
  
  #Make a list to hold the randomized values from each permutation
  dyads_list.all <- vector(mode = "list", length = length(dyads_rand1))
  
  #Combine the four for every i
  for (i in 1:length(dyads_rand1)) {
    dyads_rand1.i = dyads_rand1[[i]]
    dyads_rand2.i = dyads_rand2[[i]]
    dyads_rand3.i = dyads_rand3[[i]]
    dyads_rand4.i = dyads_rand4[[i]]
    dyads_rand.all = rbind(dyads_rand1.i,dyads_rand2.i,dyads_rand3.i,dyads_rand4.i)
    
    #Code copied from rel_groups_all
    #Look at relatedness among social classes
    classes = c("Resident Male","Resident Female","Auxilliary Male","1 Year-old Male","Dispersing Female")
    classes = expand.grid(classes,classes)
    classes = data.frame(t(apply(classes, 1, sort))) #Sort within rows - otherwise get duplicates of each pair
    colnames(classes) = c("class.x","class.y")
    classes = classes %>% distinct(class.x,class.y)
    classes$wang.mean = NA
    
    #Calculate the relatedness for each dyad of social classes 
    for (h in 1:nrow(classes)) {
      cl = dyads_rand.all[which(dyads_rand.all$class.x %in% classes$class.x[h] & 
                                  dyads_rand.all$class.y %in% classes$class.y[h]),]
      if(nrow(cl)==0) {
        classes$wang.mean[h] = NA
      } else {
        classes$wang.mean[h] = mean(cl$wang)}
    }
    
    #Get array
    group.mean.array.rand = as.matrix(acast(classes,class.x~class.y,value.var="wang.mean",fun.aggregate = sum,fill=NA_real_))
    
    #Put array into list
    dyads_list.all[[i]] = group.mean.array.rand
  }
  
  #Now melt the list of arrays and melt the observed list and combine in a graph
  dyads_list.all.m = melt(dyads_list.all)
  comsall_obs.m = melt(comsall_obs)
  
  #Get multi-level labels for each factor
  levels(dyads_list.all.m$Var1) <- gsub(" ", "\\n", levels(dyads_list.all.m$Var1))
  levels(dyads_list.all.m$Var2) <- gsub(" ", "\\n", levels(dyads_list.all.m$Var2))
  levels(comsall_obs.m$Var1) <- gsub(" ", "\\n", levels(comsall_obs.m$Var1))
  levels(comsall_obs.m$Var2) <- gsub(" ", "\\n", levels(comsall_obs.m$Var2))
  
  plot = ggplot(dyads_list.all.m,aes(x=value)) + geom_histogram(binwidth = 0.02) + facet_grid(Var1~Var2) +
    geom_vline(data=comsall_obs.m,aes(xintercept=value),color="red") + theme_bw() +
    theme(strip.background = element_rect(fill="white"),strip.text = element_text(size=12)) + 
    scale_x_continuous(limits=c(-0.25,1),breaks=c(-0.25,0,0.25,0.5,0.75,1)) + xlab("Relatedness") + ylab("Count") +
    theme(panel.grid.minor = element_blank()) + #scale_y_continuous(breaks=c(0,5,10,15)) + 
    theme(panel.spacing = unit(0.5, "lines"))
  print(plot)
  
  #Get p-values for each comparison of obeserved to random
  comsall_obs.m$p.value = NA
  for (i in 1:nrow(comsall_obs.m)) {
    rand.m = dyads_list.all.m[which(dyads_list.all.m$Var1 %in% comsall_obs.m$Var1[i] &
                                      dyads_list.all.m$Var2 %in% comsall_obs.m$Var2[i]),]
    comsall_obs.m$p.value[i] = 1 - sum(comsall_obs.m$value[i] > rand.m$value)/length(dyads_rand1)
  }
  
  #Get array of p-values and print
  comsall_obs.m.array = as.matrix(acast(comsall_obs.m,Var1~Var2,value.var="p.value",fun.aggregate = sum,fill=NA_real_))
  print(comsall_obs.m.array)
}

#Run function to compare observed mean relatedness values for each class dyad to random
#dyads_same_rand_combine_coms(dyads_rand16.coms,dyads_rand17.coms,dyads_rand18.coms,dyads_rand19.coms,comsall_obs)


