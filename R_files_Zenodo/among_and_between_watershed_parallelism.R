#######################################################################################################################################################
## Calculate among- and between-watershed QTL parallelism #############################################################################################
## First, calculate the observed parallelisms using the real QTL data #################################################################################
## Then, calculate parallelisms in randomly generated data sets that have the same structure as the real data set #####################################
## Finally, calculate statistical significance of the empirically observed parallelisms by comparing them to the parallelisms in the random data sets #
#######################################################################################################################################################

rm (list = ls(all=TRUE))
setwd("") # Set the working directory
d<-read.csv("watershed_family_trait_LG.csv", header=T) # The empirical QTL data set
iter<-50 # Number of randomly generated data sets; set manually


#########################################################
### OBSERVED BETWEEN- AND AMONG-WATERSHED PARALLELISM ###
#########################################################

# How many traits are there?
unique.traits<-unique(d$Trait)
length(unique.traits)

data<-NULL # Collect the statistics for among-watershed parallelism
data.pairwise<-NULL # Collect the statistics for between-watershed parallelism

for(i in 1:length(unique.traits)){ # loop opens for trait
    
    out<-NULL # Just making sure 'out' is empty
    dd<-d[which(d$Trait==unique.traits[i]),]
	trait<-dd[1,3] # The actual trait
    
    # A trait must have been mapped in at least 2 different watersheds to be able to calculate among- or between-watershed parallelism
    
    # A trait has mapped QTL in only one watershed
    if(length(unique(dd$Watershed))==1){
        N.qtl.tot.mapped<-length(dd[,1])
        N.watersheds<-1
        out<-cbind.data.frame(trait,N.watersheds,N.qtl.tot.mapped,N.qtl.tot.mapped,N.qtl.tot.mapped,"NA")
        }
    
    # A trait has mapped QTL in multiple watersheds
    if(length(unique(dd$Watershed))>1){
        
        N.watersheds<-length(unique(dd$Watershed)) # In how many watersheds was at least one QTL found?
        N.qtl.tot.mapped<-length(dd[,1]) # Total number of mapped QTL across all families and watersheds for this trait
        
        # Remove QTL that were repeatedly mapped *within* a watershed
        ddd<-dd[!duplicated(paste0(dd$Watershed,"_",dd$LinkageGroup)),]
        N.qtl.among.watershed<-length(ddd[,1]) # Sum of unique within-watershed LG from all watersheds
        N.qtl.unique<-length(unique(ddd$LinkageGroup)) # N of unqiue LG mapped from across all watersheds
        
        # Note that traits that mapped in only one watershed are not considered here
        pairwise.watershed.comp<-data.frame(combn(c("Boot","Misty","Pye","Roberts"),2))
        for(p in 1:length(pairwise.watershed.comp[1,])){
                
                # For the case where each of two watersheds from a specific pairwise comparison has mapped QTL for a trait
                # -> in this case, between-watershed parallelism can be calculated
                if(length(grep(pairwise.watershed.comp[1,p],ddd$Watershed))!=0 &&
                   length(grep(pairwise.watershed.comp[2,p],ddd$Watershed))!=0){
                       ddd.a<-ddd[grep(pairwise.watershed.comp[1,p],ddd$Watershed),]
                       ddd.b<-ddd[grep(pairwise.watershed.comp[2,p],ddd$Watershed),]
                       d4<-rbind(ddd.a,ddd.b)
                       N.qtl.watershed.A<-length(ddd.a[,1])
                       N.qtl.watershed.B<-length(ddd.b[,1])
                       N.qtl.between.pairwise.watersheds<-length(d4[,1]) # Sum of unique LG
                       N.qtl.unique.between.pairwise.watersheds<-length(unique(d4$LinkageGroup)) # N of unique LG
                       prop.pairwise.parallel<-round((N.qtl.between.pairwise.watersheds-N.qtl.unique.between.pairwise.watersheds)/N.qtl.unique.between.pairwise.watersheds,3)
                       out.pairwise<-cbind.data.frame(trait,pairwise.watershed.comp[1,p],pairwise.watershed.comp[2,p],N.qtl.watershed.A,N.qtl.watershed.B,N.qtl.between.pairwise.watersheds,N.qtl.unique.between.pairwise.watersheds,prop.pairwise.parallel)
                       names(out.pairwise)[2]<-"watershed.A"
                       names(out.pairwise)[3]<-"watershed.B"
                       data.pairwise<-rbind.data.frame(data.pairwise,out.pairwise)
                       }
                   # For the case where one of the two watersheds from a specific pairwise comparison has no mapped QTL
                   # -> in this case, between-watershed parallelism cannot be calculated
                   if(length(grep(pairwise.watershed.comp[1,p],ddd$Watershed))==0 |
                      length(grep(pairwise.watershed.comp[2,p],ddd$Watershed))==0){
                          ddd.a<-ddd[grep(pairwise.watershed.comp[1,p],ddd$Watershed),]
                          ddd.b<-ddd[grep(pairwise.watershed.comp[2,p],ddd$Watershed),]
                          N.qtl.watershed.A<-length(ddd.a[,1])
                          N.qtl.watershed.B<-length(ddd.b[,1])
                          N.qtl.between.pairwise.watersheds<-"NA"
                          N.qtl.unique.between.pairwise.watersheds<-"NA"
                          prop.pairwise.parallel<-"NA"
                          out.pairwise<-cbind.data.frame(trait,pairwise.watershed.comp[1,p],pairwise.watershed.comp[2,p],N.qtl.watershed.A,N.qtl.watershed.B,N.qtl.between.pairwise.watersheds,N.qtl.unique.between.pairwise.watersheds,prop.pairwise.parallel)
                          names(out.pairwise)[2]<-"watershed.A"
                          names(out.pairwise)[3]<-"watershed.B"
                          data.pairwise<-rbind.data.frame(data.pairwise,out.pairwise)
                          }
        }
        
        # Check for the rare case in which for the same trait, the same LG was mapped in three
        # watersheds. In the real data, this only happened once for "caudal peduncle depth"
        rep.measured<-as.vector(table(ddd$LinkageGroup))
        if(length(rep.measured[rep.measured>2])>0 && rep.measured[rep.measured>2]==3){
                prop.parallel<-round(((N.qtl.among.watershed-1)-N.qtl.unique)/N.qtl.unique,3)}
        
        # For the same trait, the same LG was normally not mapped more than twice among all watersheds
        if(length(rep.measured[rep.measured>2])==0){
                prop.parallel<-round((N.qtl.among.watershed-N.qtl.unique)/N.qtl.unique,3)}
        
        # Write out data:
        out<-cbind.data.frame(trait,N.watersheds,N.qtl.tot.mapped,N.qtl.among.watershed,N.qtl.unique,prop.parallel)
        }
    
    # Rename columns to for more intuitive names and make them consistent for attaching 'out' to 'data'
    names(out)[1]<-"trait"
    names(out)[2]<-"N.watersheds.with.qtl" # In how many watersheds was at least one QTL found for this trait?
    names(out)[3]<-"N.qtl.tot.mapped"      # Sum of total (not unique!) QTL mapped across all families and watersheds
    names(out)[4]<-"N.qtl.among.watershed" # Sum of unique LG across all watersheds (i.e., repeatedly mapped LG among families from the same watershed were counted only once
    names(out)[5]<-"N.qtl.unique"          # How many unique LG are there across all families and watersheds?
    names(out)[6]<-"prop.parallel"         # Proportion of parallel QTL among watersheds
    data<-rbind.data.frame(data,out)
}

### (i) Among-watershed parallelism ###
head(data) # Among-watershed parallelism for every trait
mean.observed.parallelism<-suppressWarnings(mean(as.numeric(data$prop.parallel),na.rm=T))
round(mean.observed.parallelism,4) # Among-watershed parallelism across all traits

### (ii) Between-watershed parallelism ###
data.pairwise$pop.pair<-paste0(data.pairwise$watershed.A,"_",data.pairwise$watershed.B)
unique.pop.pairs<-unique(data.pairwise$pop.pair)
stats.pairwise.parallelism<-NULL
for(u in 1:length(unique.pop.pairs)){
    sub<-subset(data.pairwise,data.pairwise$pop.pair==unique.pop.pairs[u])
    average.pairwise.parallelism<-suppressWarnings(mean(as.numeric(sub$prop.pairwise.parallel),na.rm=T))
    sub.out<-cbind.data.frame(sub[1,2],sub[1,3],average.pairwise.parallelism)
    names(sub.out)[1]<-"watershed.A"
    names(sub.out)[2]<-"watershed.B"
    stats.pairwise.parallelism<-rbind.data.frame(stats.pairwise.parallelism,sub.out)
}
head(data.pairwise) # Between-watershed parallelism for every trait
stats.pairwise.parallelism # Between-watershed parallelisms across all traits


#########################################################################################
### OBSERVED BETWEEN- AND AMONG-WATERSHED PARALLELISM IN RANDOMLY GENERATED DATA SETS ###
#########################################################################################
# -> In essence, we now run the script from above but on randomly generated data sets that follow the same structure as the real (empirical) data set

# Note that we'll use chromosome lengths (from Jones et al. 2012, Nature) to weigh the probability of a chromosome to be picked when generating random data sets
chrL<-structure(list(LG = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L,
11L, 12L, 13L, 14L, 15L, 16L, 17L, 18L, 20L, 21L), length = c(28185914L,
23295652L, 16798506L, 32632948L, 12251397L, 17083675L, 27937443L,
19368704L, 20249479L, 15657440L, 16706052L, 18401067L, 20083130L,
15246461L, 16198764L, 18115788L, 14603141L, 16282716L, 19732071L,
11717487L), genes = c(1257L, 860L, 932L, 1323L, 732L, 721L, 1320L,
881L, 1012L, 815L, 1058L, 1003L, 970L, 736L, 778L, 801L, 702L,
762L, 931L, 463L)), class = "data.frame", row.names = c(NA, -20L
))
chrL$prob<-(chrL$length)/sum(chrL$length) # Calculate probabilities based on relative chromosome lengths
collect.picked.chrs<-NULL # This object will collect the picked chromosomes; can be used to check whether
                          # longer chromosomes are indeed picked more often, as expected given the above weighing

random.parallelisms<-NULL # Object to collect stats for among-watershed parallelism
random.data.pairwise<-NULL # Object to collect stats for between-watershed parallelism

for(t in 1:iter){ # Iter-loop opens
    
    # Collect the exact same information as for the real data, but using randomly-generated data sets
    data.neutral<-NULL

    for(i in 1:length(unique.traits)){ # Loop that cycles through each trait
        out<-NULL # Just making sure 'out' is empty
        dd<-d[which(d$Trait==unique.traits[i]),]
        trait<-dd[1,3] # The actual trait
        
        dd.temp<-NULL
        dd$family.watershed<-paste0(dd$Family,"_",dd$Watershed)
        unique.family.watershed<-unique(dd$family.watershed)
        for(w in 1:length(unique.family.watershed)){
            w.temp<-subset(dd,dd$family.watershed==unique.family.watershed[w])
            
            # Drawing random chromosomes, weighted by chromosome length
            w.temp$LinkageGroup<-sample(c(chrL$LG),length(w.temp$LinkageGroup), prob=c(chrL$prob),replace = F)
            collect.picked.chrs<-c(collect.picked.chrs,w.temp$LinkageGroup)
            dd.temp<-rbind(dd.temp,w.temp)
        }
        dd<-dd.temp
    
        # A trait must have been mapped in at least two watersheds to be able to calculate among- or between-watershed parallelism
    
        # A trait has mapped QTL in only one watershed
        if(length(unique(dd$Watershed))==1){
            N.qtl.tot.mapped<-length(dd[,1])
            N.watersheds<-1
            out<-cbind.data.frame(trait,N.watersheds,N.qtl.tot.mapped,N.qtl.tot.mapped,N.qtl.tot.mapped,"NA")}
    
        # A trait has mapped QTL in multiple watersheds
        if(length(unique(dd$Watershed))>1){
        
            N.watersheds<-length(unique(dd$Watershed)) # In how many watersheds was at least one QTL found?
            N.qtl.tot.mapped<-length(dd[,1]) # Total number of mapped QTL across all families and watersheds for this trait
        
            ## Remove QTL that were repeatedly mapped *within* a watershed
            ddd<-dd[!duplicated(paste0(dd$Watershed,"_",dd$LinkageGroup)),]
            N.qtl.among.watershed<-length(ddd[,1]) # Total number of mapped LG
            N.qtl.unique<-length(unique(ddd$LinkageGroup)) # Number of uniquely mapped LG
            
            # For the case where each of two watersheds from a specific pairwise comparison has mapped QTL for a trait
            # -> in this case, between-watershed parallelism can be calculated
            pairwise.watershed.comp<-data.frame(combn(c("Boot","Misty","Pye","Roberts"),2))
            for(p in 1:length(pairwise.watershed.comp[1,])){
                    if(length(grep(pairwise.watershed.comp[1,p],ddd$Watershed))!=0 &&
                       length(grep(pairwise.watershed.comp[2,p],ddd$Watershed))!=0){
                           ddd.a<-ddd[grep(pairwise.watershed.comp[1,p],ddd$Watershed),]
                           ddd.b<-ddd[grep(pairwise.watershed.comp[2,p],ddd$Watershed),]
                           d4<-rbind(ddd.a,ddd.b)
                           N.qtl.watershed.A<-length(ddd.a[,1])
                           N.qtl.watershed.B<-length(ddd.b[,1])
                           N.qtl.between.pairwise.watersheds<-length(d4[,1])
                           N.qtl.unique.between.pairwise.watersheds<-length(unique(d4$LinkageGroup))
                           prop.pairwise.parallel<-round((N.qtl.between.pairwise.watersheds-N.qtl.unique.between.pairwise.watersheds)/N.qtl.unique.between.pairwise.watersheds,3)
                           out.pairwise<-cbind.data.frame(t,trait,pairwise.watershed.comp[1,p],pairwise.watershed.comp[2,p],N.qtl.watershed.A,N.qtl.watershed.B,N.qtl.between.pairwise.watersheds,N.qtl.unique.between.pairwise.watersheds,prop.pairwise.parallel)
                           names(out.pairwise)[1]<-"random.iteration"
                           names(out.pairwise)[3]<-"watershed.A"
                           names(out.pairwise)[4]<-"watershed.B"
                           random.data.pairwise<-rbind.data.frame(random.data.pairwise,out.pairwise)
                           }
                       
            }
    
            # Check for the rare case in which for a trait, the same LG was mapped in three watersheds
            rep.measured<-as.vector(table(ddd$LinkageGroup))
            if(length(rep.measured[rep.measured>2])>0 && rep.measured[rep.measured>2]==3){ # C-loop opens
                prop.parallel<-round(((N.qtl.among.watershed-1)-N.qtl.unique)/N.qtl.unique,3)}
            
            # For the same trait, the same LG was normally not mapped more than twice among all watersheds
            if(length(rep.measured[rep.measured>2])==0){
                prop.parallel<-round((N.qtl.among.watershed-N.qtl.unique)/N.qtl.unique,3)}
            
            # Write out data
            out<-cbind.data.frame(trait,N.watersheds,N.qtl.tot.mapped,N.qtl.among.watershed,N.qtl.unique,prop.parallel)
            }
        
        # Rename columns to for more intuitive names and make them consistant for attaching 'out' to 'data'
        names(out)[1]<-"trait"
        names(out)[2]<-"N.watersheds.with.qtl"  # In how many watersheds was at least one QTL found for this trait?
        names(out)[3]<-"N.qtl.tot.mapped"       # Sum of total (not unique!) QTL mapped across all families and watersheds
        names(out)[4]<-"N.qtl.among.watershed"  # Sum of unique LG across all watersheds (i.e., repeatedly mapped LG among families from the same watershed were counted only once
        names(out)[5]<-"N.qtl.unique"           # How many unique LG are there across all families and watersheds?
        names(out)[6]<-"prop.parallel"          # Proportion of parallel QTL among watersheds
        data.neutral<-rbind(data.neutral,out)
        
    }
    # Among-watershed parallelism
    mean.parallelism.random<-suppressWarnings(mean(as.numeric(data.neutral$prop.parallel),na.rm=T))
    random.parallelisms<-c(random.parallelisms,mean.parallelism.random)
}

#######################################
### (i) Among-watershed parallelism ###

# First, check whether the 'random' chromosomes picked weigh by length. Expected is a positive correlation
# between relative chromosome length and the number of times a chromosome was picked:
plot(chrL$prob,table(collect.picked.chrs))
# looking good!

# Second, how many iterations produced greater parallelism than observed in the real data?
cnt<-length(random.parallelisms[random.parallelisms>=mean.observed.parallelism])
# Total number of iterations:
iter<-length(random.parallelisms)
# Calculate P-value accordingly:
p.value<-round(cnt/iter, digits=4)
p.value

# Plot these results as histogram
hist(random.parallelisms,breaks=20,las=1,xlab = "Neutral expectation for among-watershed parallelism",xlim=c(0,0.25),main=paste0("p.value = ",p.value))
segments(mean.observed.parallelism,0,mean.observed.parallelism,4000,col="red",lwd=2)
segments(mean(random.parallelisms),0,mean(random.parallelisms),4000,col="blue",lwd=2)
text(0.17, y = 1700, labels = round(mean.observed.parallelism,3),cex=0.85,col="red")
text(0.11, y = 1700, labels = round(mean(random.parallelisms),3),cex=0.85,col="blue")


##########################################
### (ii) Between-watershed parallelism ###

N.iter<-length(unique(random.data.pairwise$random.iteration))
iteration<-unique(random.data.pairwise$random.iteration)
stats.random.pairwise.parallelism<-NULL
random.data.pairwise$pop.pair<-paste0(random.data.pairwise$watershed.A,"_",random.data.pairwise$watershed.B)
for(i in 1 : N.iter){
    pair.sub<-subset(random.data.pairwise,random.data.pairwise$random.iteration==iteration[i])
    ## Pairwise between-watershed parallelism
    unique.pop.pairs<-unique(random.data.pairwise$pop.pair)
    for(u in 1:length(unique.pop.pairs)){
        sub<-subset(pair.sub,pair.sub$pop.pair==unique.pop.pairs[u])
        average.pairwise.parallelism<-suppressWarnings(mean(as.numeric(sub$prop.pairwise.parallel),na.rm=T))
        sub.out<-cbind.data.frame(iteration[i],sub[1,3],sub[1,4],average.pairwise.parallelism)
        names(sub.out)[1]<-"iter"
        names(sub.out)[2]<-"watershed.A"
        names(sub.out)[3]<-"watershed.B"
        stats.random.pairwise.parallelism<-rbind(stats.random.pairwise.parallelism,sub.out)
    }
}
# "stats.random.pairwise.parallelism" contains the average pairwise parallelism for each random iteration
# Now, compare how many of these random iterations have a higher parallelism than the observed parallelism in the real data
# These real data are in "stats.pairwise.parallelism"
stats.pairwise.parallelism$pop.pair<-paste0(stats.pairwise.parallelism$watershed.A,"_",stats.pairwise.parallelism$watershed.B)
stats.random.pairwise.parallelism$pop.pair<-paste0(stats.random.pairwise.parallelism$watershed.A,"_",stats.random.pairwise.parallelism$watershed.B)
unique.pop.pairs<-unique(stats.pairwise.parallelism$pop.pair)

# Calculate stats for each pairwise comparison
stats.pairwise<-NULL
for(c in 1:length(unique.pop.pairs)){
    random.parallelisms<-subset(stats.random.pairwise.parallelism,stats.random.pairwise.parallelism$pop.pair==unique.pop.pairs[c])
    observed.parallelism<-subset(stats.pairwise.parallelism,stats.pairwise.parallelism$pop.pair==unique.pop.pairs[c])$average.pairwise.parallelism
    # First, how many iterations produced greater parallelism than observed in the real data?
    cnt<-length(which(random.parallelisms$average.pairwise.parallelism>=observed.parallelism))
    # Total number of iterations:
    random.iter<-length(random.parallelisms[,1])
    # Calculate P-value
    p.value<-round(cnt/random.iter, digits=4)
    average.random.parallelism<-round(mean(random.parallelisms$average.pairwise.parallelism),4)
    stats.pairwise<-rbind.data.frame(stats.pairwise,cbind.data.frame(unique.pop.pairs[c],observed.parallelism,random.iter,average.random.parallelism,p.value))
}
names(stats.pairwise)[1]<-"pop.pair"
names(stats.pairwise)[3]<-"totN.random.iter"
stats.pairwise
