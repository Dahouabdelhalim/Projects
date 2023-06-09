### Disentangling the causes of age-assortative mating in bird populations with 
  # contrasting life-history strategies - ANALYSIS

### DATA
# XXBreedingData        ===>  Dataframe (either "GTBreedingData" or "MSBreedingData") with columns 
#                             'Year' (character), 'Male' (character), 'Female' (character), 'Male.age' 
#                             (numeric) & 'Female.age' (numeric). Each row represents a breeding attempt 
#                             where at least one egg was laid in the given year.
# MSRoundup_X           ===>  Dataframe (either "MSRoundup_Males" or "MSRoundup_Females) with columns 
#                             'Year' (character), 'Male' (character), 'Female' (character), 'Male.age' 
#                             (numeric) & 'Female.age' (numeric). Each row represents an individual of the 
#                             focal sex (depending on the dataframe), which has either been caught breeding 
#                             (in which case, the other sex's column includes ring number of its partner), or 
#                             was caught as a non-breeding individual in the roundup (in which case, the other 
#                             sex's column is NA).
# MSRoundup_ExpProbX    ===>  Dataframe (either "MSRoundup_ExpProbMales" or "MSRoundup_ExpProbFemales) with 
#                             columns 'Year' (character), 'Male' (character), 'Female' (character), 'Male.age' 
#                             (numeric) & 'Female.age' (numeric). Each row represents an individual of the focal
#                             sex (depending on the dataframe), which has either been caught breeding (in which 
#                             case, the other sex's column includes ring number of its partner), or was caught 
#                             as a non-breeding individual in the roundup (in which case, the other sex's column
#                             is NA). The dataframe is restricted such that the number of individuals of the focal
#                             sex caught at the roundup each year is determined by age-specific breeding probabilities 
#                             calculated from the longterm data (i.e. the number of individuals of a given age from
#                             the roundup included each year = the number of individuals of the given age in the 
#                             roundup * the age-specific breeding probability calculated from the longterm data).
# XXAnnualProportions ===>    Dataframe (either "GTAnnualProportions" or "MSAnnualProportions") with columns
#                             'Year' (character), 'Prop.assorted' (= proportion of breeding pairs in the
#                             population that are assorted by age; numeric), 'Prop.first' (proportion of
#                             individuals that are breeding for the first-time in their life-history;
#                             numeric), 'Prop.young' (= proportion of individuals that are 'young breeders' 
#                             (yearling great tits; or 2-4-year-old mute swans); numeric), 'Prop.new.exp'
#                             (= proportion of breeding pairs where both partners have attempted to breed in a 
#                             previous year, but with a different individual to the current year, therefore 
#                             representing pairs that did not form at either partnersâ€™ first breeding; numeric).

### Load packages
library(dplyr)
library(tidyverse)

### ANALYSIS

## Spearman's rank correlation in partner age ####

# Across all time
cor.test(XXBreedingData$Male.age,
         XXBreedingData$Female.age,  
         method = "spearman", exact = FALSE)

## END ####

## Permutation analysis to simulate random pairing ####

# Create a vector of assessed years
Years <- XXBreedingData$Year %>%        # Or substitute "XXBreedingData" with
                                        # "MuteSwanRoundUp_XX" to assess years with
                                        # non-breeding individuals from the roundup
  unique %>% sort

# Create list from breeding data, split by year
Breeding.list<-lapply(Years,
                 function(a)XXBreedingData[XXBreedingData$Year %in%a,]) # Or substitute "XXBreedingData" 
                                                                        # with "MuteSwanRoundUp_XX" to 
                                                                        # include non-breeding individuals 
                                                                        # from the roundup
names(Breeding.list)<-Years

# Calculate age-difference in 1000 permutations of random pairing
nperm<-1000           # set as 1000 permutations
Perm.age.diff<-NULL   # create a blank variable to store into
Perm.number<-NULL     # create a blank variable to store into

for(i in 1:nperm){
  print(i)            # Track permutation progression
  Age.difference<-unlist(lapply(Breeding.list,function(a){
    abs(sample(a$Male.age)-
              (a$Female.age))}))
  Perm.age.diff<-append(Perm.age.diff,Age.difference)
  Perm.number<-append(Perm.number,rep(i,
                                      length(Age.difference)))
}

# Combine into dataframe (and add in a column listing the number of pairs
  # each year)
Age.difference<-as.data.frame(
  cbind(Perm.number,Perm.age.diff))

Age.difference$Check<-c(1:nrow(Age.difference))

# Turn into wide format
Age.difference<-Age.difference%>%
  spread(Perm.number,Perm.age.diff)%>%
  select(2:1001)

# Frequency table for each permutation
Age.difference.freq.table<-Age.difference%>%
  map(table)

# Make dataset for freq of age-assorted pairs
Freq.age.assorted<-
  as.data.frame(lapply(Age.difference.freq.table, `[`, '0'))

# How many age-assorted pairs compared to null?
Freq.age.assorted<-as.data.frame(Freq.age.assorted)
Freq.age.assorted<-
  gather(Freq.age.assorted, key="Perm",
         value="Freq", 1:1000)

# Quantiles...
quantile(Freq.age.assorted$Freq,probs=c(0.025,0.975))

# ...Compared to observed
table(XXBreedingData$Male.age==XXBreedingData$Female.age)

## END ####

## Spearman's rank in annual breeding population characteristics ####

# Relationship between within-year proportion of first-time breeders and 
  # proportion of age-assorted pairs 
cor.test(XXAnnualProportions$Prop.assorted,
         XXAnnualProportions$Prop.first,
         method = "spearman", exact = FALSE)

# Relationship between within-year proportion of young breeders and 
  # proportion of age-assorted pairs
cor.test(XXAnnualProportions$Prop.assorted,
         XXAnnualProportions$Prop.young,
         method = "spearman", exact = FALSE)

# Relationship between within-year proportion of newly-formed experienced pairs and 
  # proportion of age-assorted pairs
cor.test(XXAnnualProportions$Prop.assorted,
         XXAnnualProportions$Prop.new.exp,
         method = "spearman", exact = FALSE)


## END ####

## Permutation analysis to simulate the relationship between population age-structure
  # and age-assortative mating consistent with random pairing ####
  
# Create a vector of assessed year
Years <- XXBreedingData$Year %>%        # Or substitute "XXBreedingData" with
                                        # "MuteSwanRoundUp_XX" to assess years with
                                        # non-breeding individuals from the roundup
  unique %>% sort

# Create list from breeding data, split by year
Breeding.list<-lapply(Years,
                      function(a)XXBreedingData[XXBreedingData$Year %in%a,])

names(Breeding.list)<-Years

# Caclulate 1000 correlation coefficients from within-year permutations of random pairing

nperm<-1000     # set as 1000 permutations
perm.cor<-NA    # create a blank variable to store into
for(i in 1:nperm){
  print(i)            # Track permutation progression
  
  # Calculate the proportion age-assorted (Prop.assorted.perm) for each year from random pairing:
  Prop.assorted.perm<-unlist(lapply(Breeding.list,function(a){
    mean(sample(a$Female.age)==
        (a$Male.age))}
  ))
  
  # Check the correlation between the proportion of age-assorted pairs in these permuted 
  # datasets and proportion of 'young breeders' (which remains constant across permutations)
  
  Prop.young_Prop.assorted.PERM <- cor.test(XXAnnualProportions$Prop.young,
                                            Prop.assorted.perm[as.character(XXAnnualProportions$Year)],
                                            method = "spearman", exact = FALSE) 
  
  perm.cor[i]<-Prop.young_Prop.assorted.PERM$estimate # Store the estimate
}

# How does the real correlation compare to the permuted correlations?
# Quantiles...
quantile(perm.cor,probs=c(0.025,0.975))

# ...Compared to observed
(cor.test(XXAnnualProportions$Prop.assorted,
         XXAnnualProportions$Prop.young,
         method = "spearman", exact = FALSE))$estimate

## END ####

