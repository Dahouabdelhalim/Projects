#######################################################################################################################################################
## Calculate within-watershed QTL parallelism (i.e., among different F1 families from the same watershed) #############################################
## First, calculate the observed within-watershed parallelism for each watershed using the real QTL data ##############################################
## Then, calculate within-watershed parallelism in randomly generated data sets that have the same structure as the real data set #####################
## Finally, calculate statistical significance of the empirically observed parallelisms by comparing them to the parallelisms in the random data sets #
#######################################################################################################################################################

rm (list = ls(all=TRUE))
setwd("") # Set the working directory
d<-read.csv("watershed_family_trait_LG.csv", header=T) # The empirical QTL data set
iter<-1000 # Number of randomly generated data sets

#############################################
### OBSERVED WITHIN-WATERSHED PARALLELISM ###
#############################################

# How many traits are there?
unique.traits<-unique(d$Trait)
length(unique.traits)

data<-NULL # Collect the statistics
for(i in 1:length(unique.traits)){ # loop opens for trait
    dd<-d[which(d$Trait==unique.traits[i]),]
	trait<-dd[1,3] # The actual trait
    watersheds<-unique(dd$Watershed)

    out<-NULL # Just making sure 'out' is empty
    for(w in 1:length(watersheds)){ # loop opens for watershed
        ddd<-dd[which(dd$Watershed==watersheds[w]),]
        
        # If there are no mapped QTL for this watershed
        if(length(unique(ddd$Family))==0){
            N.families<-0
            N.qtl.tot.mapped<-"NA"
            out<-cbind.data.frame(t,trait,watersheds[w],N.families,N.qtl.tot.mapped,N.qtl.tot.mapped,N.qtl.tot.mapped,as.numeric("NA"))
            }
        
        # If there is only one F1 cross-family within a certain watershed with mapped QTL (if so, no within-watershed parallelism can be calculated)
        if(length(unique(ddd$Family))==1){
            N.families<-1
            N.qtl.tot.mapped<-length(ddd[,1])
            out<-cbind.data.frame(trait,watersheds[w],N.families,N.qtl.tot.mapped,N.qtl.tot.mapped,N.qtl.tot.mapped,as.numeric("NA"))
            }
    
        # If there are multiple F1 families within this watershed with mapped QTL
        if(length(unique(ddd$Family))>1){
        
            N.families<-length(unique(ddd$Family)) # In how many families was at least one QTL found?
            N.qtl.tot.mapped<-length(ddd[,1]) # Total number of mapped QTL across all families and watersheds for this trait
        
            # Remove QTL that were repeatedly mapped *within* a watershed
            d4<-ddd[!duplicated(paste0(ddd$Family,"_",ddd$LinkageGroup)),]
            N.qtl.among.families<-length(d4[,1]) # Sum of unique within-watershed LG from all families
            N.qtl.unique<-length(unique(d4$LinkageGroup)) # N of unique LG mapped across all families
        
            # Check for the rare case in which the same LG was mapped three times for the same trait within a watershed (not sure this occurs at all)
            rep.measured<-as.vector(table(d4$LinkageGroup))
            if(length(rep.measured[rep.measured>2])>0 && rep.measured[rep.measured>2]==3){
                prop.parallel<-round(((N.qtl.among.families-1)-N.qtl.unique)/N.qtl.unique,3)}
        
            # For the same trait, the same LG was normally not mapped more than twice across the different F1 families from the same watershed
            if(length(rep.measured[rep.measured>2])==0){
                prop.parallel<-round((N.qtl.among.families-N.qtl.unique)/N.qtl.unique,3)}
        
            # Write out data statistics
            out<-cbind.data.frame(trait,watersheds[w],N.families,N.qtl.tot.mapped,N.qtl.among.families,N.qtl.unique,as.numeric(prop.parallel))
            }
    
        # Rename columns to for more intuitive names and make them consistent for attaching 'out' to 'data'
        names(out)[1]<-"trait"
        names(out)[2]<-"watershed"
        names(out)[3]<-"N.families.with.qtl"   # In how many families was at least one QTL found for this trait?
        names(out)[4]<-"N.qtl.tot.mapped"      # Sum of total (not unique!) QTL mapped across all families across all watersheds
        names(out)[5]<-"N.qtl.among.families"  # Sum of unique within-watershed LG from all watersheds (i.e., repeatedly
                                               # mapped LG among families from the same watershed were counted only once) -> should be the same as "N.qtl.tot.mapped"
        names(out)[6]<-"N.qtl.unique"          # How many unique LG are there across all families and watersheds?
        names(out)[7]<-"prop.parallel"         # Proportion of parallel QTL within a watershed for each trait
        data<-rbind.data.frame(data,out)
    } # watershed-loop closes
} # trait-loop closes
# -> warnings are no problem here, ignore

# Within-watershed parallelism
unique.watersheds<-unique(data$watershed)
stats.withinWatershed.parallelism<-NULL
for(u in 1:length(unique.watersheds)){
    sub<-subset(data,data$watershed==unique.watersheds[u])
    average.withinWatershed.parallelism<-suppressWarnings(mean(as.numeric(sub$prop.parallel),na.rm=T))
    sub.out<-cbind.data.frame(unique.watersheds[u],average.withinWatershed.parallelism)
    names(sub.out)[1]<-"watershed"
    stats.withinWatershed.parallelism<-rbind.data.frame(stats.withinWatershed.parallelism,sub.out)
}
stats.withinWatershed.parallelism

####################################################################
### WITHIN-WATERSHED PARALLELISM IN RANDOMLY GENERATED DATA SETS ###
####################################################################
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

# Now calculate the exact same information as for the real data, but using randomly generated data sets

data.neutral.parallelism<-NULL # Object to collect data / statistics

for(t in 1:iter){ # iter-loop opens
    
    # First, generate a random data set
    unique.traits<-unique(d$Trait)
    
    d.randomized<-NULL
    for(i in 1:length(unique.traits)){ # Loop that cycles through each trait
            a<-subset(d,d$Trait==unique.traits[i])
            unique.watersheds<-unique(a$Watershed)
            
            for(f in 1:length(unique.watersheds)){
                b<-subset(a,a$Watershed==unique.watersheds[f])
                
                unique.families<-unique(b$Family)
                for(x in 1:length(unique.families)){
                    c<-subset(b,b$Family==unique.families[x])
                    r.temp<-c
                    
                    # Drawing random chromosomes, weighted by chromosome length
                    r.temp$LinkageGroup<-sample(c(chrL$LG),length(r.temp$LinkageGroup), prob=c(chrL$prob),replace = F)
                    collect.picked.chrs<-c(collect.picked.chrs,r.temp$LinkageGroup)
                    d.randomized<-rbind(d.randomized,r.temp)
                    }
                }
            }
    # Random data generated!
    
    # Now, calculate within-watershed parallelism in the random data set
    watersheds<-unique(d.randomized$Watershed) # different watersheds

    for(i in 1:length(unique.traits)){ # loop opens for trait
        dd<-d.randomized[which(d.randomized$Trait==unique.traits[i]),]
        trait<-dd[1,3] # The actual trait

        out<-NULL # Just making sure 'out' is empty
        for(w in 1:length(watersheds)){ # loop opens for watersheds
            ddd<-dd[which(dd$Watershed==watersheds[w]),]
            
            # If there are no mapped QTL for this watershed
            if(length(unique(ddd$Family))==0){
                N.families<-0
                N.qtl.tot.mapped<-"NA"
                out<-cbind.data.frame(t,trait,watersheds[w],N.families,N.qtl.tot.mapped,N.qtl.tot.mapped,N.qtl.tot.mapped,as.numeric("NA"))
                }
            
            # If there is only one F1 cross-family within a certain watershed with mapped QTL (if so, no within-watershed parallelism can be calculated)
            if(length(unique(ddd$Family))==1){
                N.families<-1
                N.qtl.tot.mapped<-length(ddd[,1])
                out<-cbind.data.frame(t,trait,watersheds[w],N.families,N.qtl.tot.mapped,N.qtl.tot.mapped,N.qtl.tot.mapped,as.numeric("NA"))
                }
        
            # If there are multiple F1 families within this watershed with mapped QTL
            if(length(unique(ddd$Family))>1){
            
                N.families<-length(unique(ddd$Family)) # In how many families was at least one QTL found?
                N.qtl.tot.mapped<-length(ddd[,1]) # Total number of mapped QTL across all families and watersheds for this trait
            
                # Remove QTL that were repeatedly mapped *within* a watershed
                d4<-ddd[!duplicated(paste0(ddd$Family,"_",ddd$LinkageGroup)),]
                N.qtl.among.families<-length(d4[,1]) # Sum of unique within-watershed LG from all families
                N.qtl.unique<-length(unique(d4$LinkageGroup)) # N of unique LG mapped from across all families
            
                # Check for the rare case in which the same LG was mapped three times for the same trait within a watershed (not sure this occurs at all)
                rep.measured<-as.vector(table(d4$LinkageGroup))
                if(length(rep.measured[rep.measured>2])>0 && rep.measured[rep.measured>2]==3){
                    prop.parallel<-round(((N.qtl.among.families-1)-N.qtl.unique)/N.qtl.unique,3)}
            
                # For the same trait, the same LG was normally not mapped more than twice across the different F1 families from the same watershed
                if(length(rep.measured[rep.measured>2])==0){
                    prop.parallel<-round((N.qtl.among.families-N.qtl.unique)/N.qtl.unique,3)}
            
                # Write out data:
                out<-cbind.data.frame(t,trait,watersheds[w],N.families,N.qtl.tot.mapped,N.qtl.among.families,N.qtl.unique,as.numeric(prop.parallel))
                }
        
            # Rename columns to for more intuitive names and make them consistant for attaching 'out' to 'data'
            names(out)[1]<-"random.iteration"
            names(out)[2]<-"trait"
            names(out)[3]<-"watershed"
            names(out)[4]<-"N.families.with.qtl"
                out[,4]<-as.numeric(as.character(out[,4]))
            names(out)[5]<-"N.qtl.tot.mapped"      # Sum of total (not unique!) QTL mapped across all families and watersheds.
                out[,5]<-as.numeric(as.character(out[,5]))
            names(out)[6]<-"N.qtl.among.families"  # Sum of unique within-watershed LG from all watersheds (i.e., repeatedly
                                                   # mapped LG among families from the same watershed were counted only once)
                out[,6]<-as.numeric(as.character(out[,6]))
            names(out)[7]<-"N.qtl.unique"          # How many unique LG are there across all families and watersheds?
                out[,7]<-as.numeric(as.character(out[,7]))
            names(out)[8]<-"prop.parallel"         # Proportion of parallel QTL among watersheds
                out[,8]<-as.numeric(as.character(out[,8]))
            data.neutral.parallelism<-rbind.data.frame(data.neutral.parallelism,out)
        } # watershed-loop closes
    } # trait-loop closes
} # iteration loop closes
# -> warnings are no problem, ignore
#write.csv(data.neutral.parallelism,"withinWatershed_data.neutral.parallelism.csv", quote=F, row.names=F) # Write out data, if you want


###################################################################
### WITHIN-WATRSHED PARALLELISM IN RANDOMLY GENERATED DATA SETS ###
###################################################################

#data.neutral.parallelism<-read.csv("withinWatershed_data.neutral.parallelism.csv",stringsAsFactor=F) # Read-in data if needed

# Stats for within-watershed parallelism
unique.watersheds<-unique(data.neutral.parallelism$watershed)
withinWatershed.parallelism.significance<-NULL

for(u in 1:length(unique.watersheds)){ # watershed-loop opens
    sub<-subset(data.neutral.parallelism,data.neutral.parallelism$watershed==unique.watersheds[u])
    iteration.averages<-NULL
    iter<-unique(sub$random.iteration)
    for(i in 1:length(iter)){
        temp<-subset(sub,sub$random.iteration==iter[i])
        iteration.averages<-c(iteration.averages,suppressWarnings(mean(as.numeric(temp$prop.parallel),na.rm=T)))
        }
    average.random.parallelism<-mean(iteration.averages)
    
    # Retrieve the observed within-watershed parallelism for this watershed
    observed.parallelism<-subset(stats.withinWatershed.parallelism,stats.withinWatershed.parallelism$watershed==unique.watersheds[u])[1,2]

    # How many random data sets produced greater within-watershed parallelism than observed in the real data?
    cnt<-length(which(iteration.averages>=observed.parallelism))
    totN.random.iter<-length(iter)
    
    # Calculate resampling P-value
    p.value<-round(cnt/totN.random.iter, digits=4)
    withinWatershed.parallelism.significance<-rbind.data.frame(withinWatershed.parallelism.significance,cbind.data.frame(unique.watersheds[u],observed.parallelism,totN.random.iter,average.random.parallelism,p.value))
} # watershed-loop closes

# Show/write out significance statistics
head(data.neutral.parallelism)
#write.csv(data.neutral.parallelism,"withinWatershed_neutralExpectation_parallelism.per.watershed.and.iteration.csv",quote=F, row.names=F)
names(withinWatershed.parallelism.significance)[1]<-"watershed"
withinWatershed.parallelism.significance
#write.csv(withinWatershed.parallelism.significance,"withinWatershed_neutralExpectation_overall.parallelism&stats.per.comparison.csv",quote=F, row.names=F)
