library(rlang)
library(metafor)
library(effsize)
library(compute.es)
library(esc)
library(meta)
library(forestplot)
library(MAd)
library(tidyr)
library(plyr)
library(dplyr)
library(splitstackshape)
library(schoolmath)
library(Gmisc)


# first set working directory
setwd("")
# - where you stored the file "Data_Meta-analysis_Schafft_etal.csv"

####load data set ####
recrfull <- read.csv("Data_Meta-analysis_Schafft_etal.csv", sep=";", header=T, dec=".")

#include only those that are definately freshwater studies
recrfull <- recrfull[(recrfull$Freshwater== 2),]

# convert variables to numeric
recrfull$Longitud <-as.numeric(gsub(",","",recrfull$Longitud,fixed=TRUE))
recrfull$C.Mean <-as.numeric(gsub(",","",recrfull$C.Mean,fixed=TRUE))
recrfull$C.neg <-as.numeric(as.character(recrfull$C.neg))
recrfull$C.pos <-as.numeric(as.character(recrfull$C.pos))
recrfull$C.N<-as.numeric(as.character(recrfull$C.N))
recrfull$T.Mean <-as.numeric(gsub(",","",recrfull$T.Mean,fixed=TRUE))
recrfull$F.Test <-as.numeric(gsub(",","",recrfull$F.Test,fixed=TRUE))
recrfull$T.neg <-as.numeric(as.character(recrfull$T.neg))
recrfull$T.pos <-as.numeric(as.character(recrfull$T.pos))
recrfull$T.positive... <-as.numeric(as.character(recrfull$T.positive...))
recrfull$T.N<-as.numeric(as.character(recrfull$T.N))
recrfull$C.SD <-as.numeric(as.character(recrfull$C.SD))
recrfull$C.SE <-as.numeric(as.character(recrfull$C.SE))
recrfull$T.SD <-as.numeric(as.character(recrfull$T.SD))
recrfull$T.SE <-as.numeric(as.character(recrfull$T.SE))
recrfull$tvalue <-as.numeric(as.character(recrfull$t.value))
recrfull$df <-as.numeric(as.character(recrfull$df))
recrfull$p.value <-as.numeric(gsub(",","",recrfull$p.value,fixed=TRUE))
recrfull$C.SE <-as.numeric(as.character(recrfull$C.SE))

rec <- recrfull
names(rec)

# exclude studies that compare different recreational activties 
#instead of comparing to control or different intensities 

rec$Real.control.............[is.na(rec$Real.control.............)]<-'999'
rec<-subset(rec, Real.control............. %in% c('1', '2', '3', '5','999'))

#### correcting direction of values ####

#first get rid of NA
rec$Meaning[is.na(rec$Meaning)]<-'1'

for(r in 1:nrow(rec)) {
  if (rec[r,43]=='2') {
    rec[r,49]<-rec[r,49]*-1;
    rec[r,53]<-rec[r,53]*-1
  } else 
    if (rec[r,43]=='3') {
      if (rec[r,42]=="+"){
        rec[r,49]<-rec[r,49]*-1;
        rec[r,53]<-rec[r,53]*-1
      }
    }
}

####calculating effect size ####

rec.sd <- rec[!is.na(rec$T.SD),]

ESg<-esc_mean_sd(grp1m = as.numeric(rec.sd$T.Mean), 
                 grp1sd = rec.sd$T.SD, 
                 grp1n = as.numeric(rec.sd$T.N),
                 grp2m =as.numeric(rec.sd$C.Mean), 
                 grp2sd = rec.sd$C.SD, 
                 grp2n = as.numeric(rec.sd$C.N),
                 es.type = "g")

rec.sd$es<-ESg$es
rec.sd$se<-ESg$var

### create colomn es and se 
### to have identical number of rows in rec and rec.sd#

rec$es<-NA
rec$se<-NA

#create the complementary to rec.sd# to be able to combine the data again 

rec.sd.t.na <- rec[is.na(rec$T.SD),]
rec<-rbind(rec.sd,rec.sd.t.na) #combine rows to have have full data set again

rec$Cohen.s.d<-NA
rec$Cohen.s.d.1<- NA

# calculate Hedges g with t-value / p-value

rec.p<- rec[!is.na(rec$p.value),]

EStg<-esc_t(t = as.numeric(rec.p$t.value), 
            p = as.numeric(rec.p$p.value), 
            grp1n = as.numeric(rec.p$T.N), 
            grp2n = as.numeric(rec.p$C.N), 
            es.type = "g")


rec.p$Cohen.s.d<-EStg$es
rec.p$Cohen.s.d.1<- EStg$var

###combine data to full data set ###
rec.p.na<- rec[is.na(rec$p.value),] #rec without rec.p
rec<-rbind(rec.p,rec.p.na) #combine rows to have have full data set again


#put values in right columns
for(r in 1:nrow(rec)) {if (is.na(rec[r,96])) {
  rec[r,96]<-rec[r,60];
  rec[r,97]<-rec[r,84] 
}
}

### calculate Odds Ratio ####

ES.OR<-escalc(measure = "OR", ai = T.pos, bi = T.neg, 
              ci = C.pos, di = C.neg, data = rec)

rec$Odds.ratio<-ES.OR$yi
rec$Risk.ratio<-ES.OR$vi


### Converting Effect Sizes ####

for(r in 1:nrow(rec)) {if (!is.na(rec[r,68])) 
  {rec[r,96]<-rec[r,68]*sqrt(3)/pi;rec[r,97]<-rec[r,69]*sqrt(3)/pi^2 }}


#### Effect size of correlations ####
r.es <- res(Correlation.coefficient......R.., var.r = NULL, sample.size, id = rec$Study_ID ,data = rec)

rec$Fishers.z <- r.es$d
rec$Variance <- r.es$var.d

# put values in right colums
for(r in 1:nrow(rec)) {if (is.na(rec[r,96])) {
  rec[r,96]<-rec[r,91];
  rec[r,97]<-rec[r,89] 
}
}

#### Effect size with F statistics ####

f.es<-fes(F.Test, T.N, C.N, 
          level = 95, 
          cer = 0.2, dig = 2, 
          verbose = TRUE, id=Study_ID, data=rec)

rec$Cohen.s.d <- f.es$g
rec$Cohen.s.d.1 <- f.es$var.g

#put values in right columns
for(r in 1:nrow(rec)) {if (is.na(rec[r,96])) {
  rec[r,96]<-rec[r,60];
  rec[r,97]<-rec[r,84] 
}
}

####correct positive effect size where negative effect was seen####
rec<-rec[!is.na(rec$es), ]
rec<-rec[!is.na(rec$se), ]
rec<-rec[!is.infinite(rec$se), ]

recneg<-rec[is.negative(rec$es),]
recpos<-rec[!is.negative(rec$es),]
#negative
for(r in 1:nrow(recneg)) {
  if (recneg[r,43]=='1') {
    if (recneg[r,42]=='+'){
      recneg[r,96]<-recneg[r,96]*-1
    }
  } else 
    if (recneg[r,43]=='2') {
      if (recneg[r,42]=='-'){
        recneg[r,96]<-recneg[r,96]*-1
      }
    }
}

#positive
for(r in 1:nrow(recpos)) {
  if (recpos[r,43]=='1') {
    if (recpos[r,42]=='-'){
      recpos[r,96]<-recpos[r,96]*-1
    }
  }else 
    if (recpos[r,43]=='3') {
      if (recpos[r,42]=='-'){
        recpos[r,96]<-recpos[r,96]*-1
      }
    }else 
      if (recpos[r,43]=='3') {
        if (recpos[r,42]=='+'){
          recpos[r,96]<-recpos[r,96]*-1
        }
      }else 
        if (recpos[r,43]=='3') {
          if (recpos[r,42]=='?'){
            recpos[r,96]<-recpos[r,96]*-1
          }
        }else 
          if (recpos[r,43]=='3') {
            if (recpos[r,42]=='0'){
              recpos[r,96]<-recpos[r,96]*-1
            }
          }else 
            if (recpos[r,43]=='2') {
              if (recpos[r,42]=='+'){
                recpos[r,96]<-recpos[r,96]*-1
              }
            }
}

rec<-rbind(recpos,recneg) #combine rows to have have full data set again

#check positive avoidance effect sizes
recneg<-rec[is.negative(rec$es),]
recpos<-rec[!is.negative(rec$es),]

rec.avoidpos <- recpos[(recpos$Response.measured =='avoidance'),]
rec.nonavoidpos <- recpos[!(recpos$Response.measured =='avoidance'),]
rec.avoidpos$Response.specification

rec.avoidpos$Response.specification[is.na(rec.avoidpos$Response.specification)]<-'flight / birds flying up'

for(r in 1:nrow(rec.avoidpos)) {
  if (rec.avoidpos[r,40]=='reaction distance (short vigilance)') {
    rec.avoidpos[r,96]<-rec.avoidpos[r,96]*-1
  }else if (rec.avoidpos[r,40]=='reaction distance (intense vigilance)') {
    rec.avoidpos[r,96]<-rec.avoidpos[r,96]*-1
  }else if (rec.avoidpos[r,40]=='flush distance') {
    rec.avoidpos[r,96]<-rec.avoidpos[r,96]*-1
  }else 
    if (rec.avoidpos[r,40]=='flight / birds flying up') {
      rec.avoidpos[r,96]<-rec.avoidpos[r,96]*-1
    }else 
      if (rec.avoidpos[r,40]=='Change in frequency of alert behavior of pairs of dabchicks') {
        rec.avoidpos[r,96]<-rec.avoidpos[r,96]*-1
      }else 
        if (rec.avoidpos[r,40]=='Change in frequency of pairs of dabchicks being out of sight') {
          rec.avoidpos[r,96]<-rec.avoidpos[r,96]*-1
        }else 
          if (rec.avoidpos[r,40]=='flush response of eagles perched in trees') {
            rec.avoidpos[r,96]<-rec.avoidpos[r,96]*-1
          }else 
            if (rec.avoidpos[r,40]=='flush response of eagles standing/feeding on ground') {
              rec.avoidpos[r,96]<-rec.avoidpos[r,96]*-1
            }else 
              if (rec.avoidpos[r,40]=='Disturbance flight time') {
                rec.avoidpos[r,96]<-rec.avoidpos[r,96]*-1
              }else 
                if (rec.avoidpos[r,40]=='tolerance distance') {
                  rec.avoidpos[r,96]<-rec.avoidpos[r,96]*-1
                  
                }}
recpos<-rbind(rec.avoidpos, rec.nonavoidpos)
rec<-rbind(recpos,recneg) #combine rows to have have full data set again

#### decode recreational activities ####

#create new coloumnsfor decoded variables 
rec$Rec.act<- NA #this was our column for german names for recreation categories. 
#We left this code inside (even though it will stay empty) 
#to have the same number and order of variables, 
#because otherwise we would have to change all column numbers in the following code. 

rec$Rec.a<- NA #this is the column, that we will fill now

#recreational activities by shore-open water gradient
for(r in 1:nrow(rec)) {
  if (rec[r,32]=='1') {rec[r,99]<-'Angling'} 
  else if (rec[r,32]=='1.2') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='1,2') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='2') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='2.27') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='3') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='4') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='5') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='6') {rec[r,99]<-"Plane" }
  else if (rec[r,32]=='7') {rec[r,99]<-"Swimming"}
  else if (rec[r,32]=='8') {rec[r,99]<-"Swimming"}
  else if (rec[r,32]=='9') {rec[r,99]<-"Swimming"}
  else if (rec[r,32]=='10') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='11') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='12') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='13') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='14') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='14.1') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='14,1') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='15') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='16') {rec[r,99]<-"Boating" }
  else if (rec[r,32]=='17') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='18') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='19') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='20') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='21') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='22') {rec[r,99]<-"Boating"}
  else if (rec[r,32]=='23') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='24') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='25') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='26') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='27') {rec[r,99]<-"Multiple use"}
  else if (rec[r,32]=='28') {rec[r,99]<-"Shore use" }
  else if (rec[r,32]=='29') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='30') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='31') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='32') {rec[r,99]<-"Shore use"}
  else if (rec[r,32]=='999') {rec[r,99]<-"Boating"}
}

#### decode Taxonomic groups ####
rec$taxa<-NA
rec$Group.[is.na(rec$Group.)]<-'999'

for(r in 1:nrow(rec)) {
  if (rec[r,30]=='1') {rec[r,100]<-'Invertebrates'} 
  else if (rec[r,30]=='2') {rec[r,100]<-"Fish"}
  else if (rec[r,30]=='3') {rec[r,100]<-"Amphibians"}
  else if (rec[r,30]=='4') {rec[r,100]<-"Reptiles"}
  else if (rec[r,30]=='5') {rec[r,100]<-"Birds"}
  else if (rec[r,30]=='6') {rec[r,100]<-"Mammals" }
  else if (rec[r,30]=='7') {rec[r,100]<-"Plants"}
  else if (rec[r,30]=='8') {rec[r,100]<-"Other"}
  else if (rec[r,30]=='9') {rec[r,100]<-"Other"}
  else if (rec[r,30]=='999') {rec[r,100]<-"Other"}
}

#Unify Study design column
rec$Study.type. [is.na(rec$Study.type. )]<-'nodata'

for(r in 1:nrow(rec)) {
  if (rec[r,11]=='A ') {rec[r,11]<-'A'} 
  else if (rec[r,11]=='CI ') {rec[r,11]<-"CI"}
  else if (rec[r,11]=='G ') {rec[r,11]<-"G"}
  else if (rec[r,11]=='BACI ') {rec[r,11]<-"BACI"}
  else if (rec[r,11]=='BA ') {rec[r,11]<-"BA"}
  
}

#### coloum with study, recreational activity,level and taxon ####
recu<-unite(rec, Author.s.,Study.Year, Group. ,taxa ,Level.of.biological.organization, Rec.a, col= ATR, remove = FALSE)

rec<-recu

#### get rid of NA ####
rec<-rec[!is.na(rec$es), ]
rec<-rec[!is.na(rec$se), ]
rec<-rec[!is.infinite(rec$se), ]
rec<-rec[!is.infinite(rec$se), ]
rec<-rec[!is.negative(rec$se),]


####  exclude studies ####
rec <- rec[(rec$Freshwater== 2),] #exclude marine studies # already done in the beginning
rec <- rec[!(rec$Recreational.activity== 27),] # exclude multiple uses

#exculde outlier
rec <- rec[(rec$es> -25),]
rec <- rec[(rec$se >0),]
#rename variables for modeling
rec$yi<-rec$es
rec$vi<-rec$se

#create subsample without impacts of angling on fish 
rec.angling <- rec[(rec$Rec.a=='Angling'),]#only angling studies
rec.NOTangling <- rec[!(rec$Rec.a=='Angling'),]#all except angling

rec.Fishangling <- rec.angling[(rec.angling$Group.=='2'),]#only angling+fish studies
rec.nofish.angling <- rec.angling[!(rec.angling$Group.=='2'),]#only angling+fish studies

rec1<-rbind(rec.NOTangling,rec.nofish.angling)
rec<-rec1

#### Give weights according to study quality ####

#create new variable
rec$weights<-NA

# weights for study design #
for(r in 1:nrow(rec)) {
  if (rec[r,11]=='A') {
    rec[r,104]<-'0.0904'
  } else 
    if (rec[r,11]=='A ') {
      rec[r,104]<-'0.0904'
    } else 
      if (rec[r,11]=='CI') {
        rec[r,104]<-'0.206'
      } else 
        if (rec[r,11]=='CI ') {
          rec[r,104]<-'0.206'
        }else 
          if (rec[r,11]=='BA') {
            rec[r,104]<-'0.226'
          }else 
            if (rec[r,11]=='BA ') {
              rec[r,104]<-'0.226'
            }else 
              if (rec[r,11]=='BACI') {
                rec[r,104]<-'0.4'
              }else 
                if (rec[r,11]=='BACI ') {
                  rec[r,104]<-'0.4'
                }else 
                  if (rec[r,11]=='G') {
                    rec[r,104]<-'0.3'
                  }else 
                    if (rec[r,11]=='G ') {
                      rec[r,104]<-'0.3'
                    }
}

#make variable numeric to make calculations possible

rec$weights <-as.numeric(gsub(",","",rec$weights,fixed=TRUE))

# weights for sample size of control units
for(r in 1:nrow(rec)) {
  if (rec[r,21]=='1') {
    rec[r,104]<-rec[r,104]+0.2
  } else 
    if (rec[r,21]=='1 ') {
      rec[r,104]<-rec[r,104]+0.2
    } else 
      if (rec[r,21]=='>1') {
        rec[r,104]<-rec[r,104]+0.3
      } else 
        if (rec[r,21]=='>1 ') {
          rec[r,104]<-rec[r,104]+0.3
        }
}

#weights for sample size of impact units
for(r in 1:nrow(rec)) {
  if (rec[r,22]=='1') {
    rec[r,104]<-rec[r,104]+0.1
  } else 
    if (rec[r,22]=='1 ') {
      rec[r,104]<-rec[r,104]+0.1
    } else 
      if (rec[r,22]=='>1') {
        rec[r,104]<-rec[r,104]+0.2
      } else 
        if (rec[r,22]=='>1 ') {
          rec[r,104]<-rec[r,104]+0.2
          
        }
} 

#weights for sample size of gradient response
for(r in 1:nrow(rec)) {
  if (rec[r,23]=='4') {
    rec[r,104]<-rec[r,104]+0.1
  } else 
    if (rec[r,23]=='4 ') {
      rec[r,104]<-rec[r,104]+0.1
    } else 
      if (rec[r,23]=='5') {
        rec[r,104]<-rec[r,104]+0.2
      } else 
        if (rec[r,23]=='>5') {
          rec[r,104]<-rec[r,104]+0.3
          
        }
} 

#controlled for confounding factors
#first get rid of NA
rec$confounding.factors[is.na(rec$confounding.factors)]<-'1'
rec$confounding.factors <-as.numeric(gsub(",","",rec$confounding.factors,fixed=TRUE))

#controlled for confounding factors
for(r in 1:nrow(rec)) {
  if (rec[r,24]=='2') {
    rec[r,104]<-rec[r,104]+0.1
  }
} 

##weights for real control
for(r in 1:nrow(rec)) {
  if (rec[r,47]=='2') {
    rec[r,104]<-rec[r,104]/2
  }
} 

#weights for waterbody vs zone
rec$Scale[is.na(rec$Scale)]<-'1'

for(r in 1:nrow(rec)) {
  if (rec[r,26]=='2') {
    rec[r,104]<-rec[r,104]/2
  }
}


#### study quality ####

rec$quality<-NA

# low quality #
for(r in 1:nrow(rec)) {
  if (rec[r,11]=='A') {
    rec[r,105]<-'low'
  } else 
    if (rec[r,11]=='A ') {
      rec[r,105]<-'low'
    } else 
      if (rec[r,21]=='0') { #control sampling units
        rec[r,105]<-'low'
      } else 
        if (rec[r,21]=='0 ') {
          rec[r,105]<-'low'
        } else 
          if (rec[r,22]=='0') { #impact sampling units
            rec[r,105]<-'low'
          } else 
            if (rec[r,22]=='0 ') {
              rec[r,105]<-'low'
            } else 
              if (rec[r,23]=='<4') { # replication of gradient response model
                rec[r,105]<-'low'
              } else 
                if (rec[r,23]=='<4 ') {
                  rec[r,105]<-'low'
                } else 
                  if (rec[r,24]=='1') { #confounding factors
                    rec[r,105]<-'low'
                  } else 
                    if (rec[r,24]=='1 ') {
                      rec[r,105]<-'low'
                      #} else 
                      # if (rec[r,20]=='4') { #randomization
                      #  rec[r,105]<-'low'
                      # } else 
                      #  if (rec[r,20]=='5') {
                      #   rec[r,105]<-'low'
                      # } else 
                      #  if (rec[r,20]=='999'){ 
                      #  rec[r,105]<-'low'
                    } else 
                      if (rec[r,47]=='2') { #real control
                        rec[r,105]<-'low'
                      } else 
                        if (rec[r,26]=='2'){ # scale - water body or zone
                          rec[r,105]<-'low'
                        }
}

#high quality  

for(r in 1:nrow(rec)) {
  if (rec[r,11]=='BACI') {
    if (rec[r,21]=='>1') { #control sampling units
      if (rec[r,22]=='>1') { #impact sampling units
        #if (rec[r,23]=='>5') { # replication of gradient response model
        if (rec[r,24]=='2') { #confounding factors
          if (rec[r,20]=='1') { #randomization
            # if (rec[r,47]=='2') { #real control
            if (rec[r,26]=='1'){ # scale - water body or zone
              rec[r,105]<-'high'}
          }
        }
      }
    }
  }
}

for(r in 1:nrow(rec)) {
  if (rec[r,11]=='BACI') {
    if (rec[r,21]=='>1') { #control sampling units
      if (rec[r,22]=='>1') { #impact sampling units
        #if (rec[r,23]=='>5') { # replication of gradient response model
        if (rec[r,24]=='2') { #confounding factors
          if (rec[r,20]=='2') { #randomization
            # if (rec[r,47]=='2') { #real control
            if (rec[r,26]=='1'){ # scale - water body or zone
              rec[r,105]<-'high'}
          }
        }
      }
    }
  }
}
#medium quality
rec$quality[is.na(rec$quality)]<-'medium'

#study quality
recqmedium <- rec[(rec$quality== 'medium'),]
recqlow<- rec[(rec$quality== 'low'),]


reclow<-rec[rec$weights <0.4,]
recmedium<-rec[!rec$weights<0.4,]

#### real controll vs high vs low ####

### data only high vs low ###
rec2 <- rec[(rec$Real.control.............== 2),]
### data only not real BA ###
rec3 <- rec[(rec$Real.control.............== 3),]
### data only with real control ###
rec1 <- rec[(rec$Real.control.............== 1),]
rec13<-rbind(rec1,rec3)
rec123<-rbind(rec2,rec13)


for(r in 1:nrow(rec123)) {
  if (rec123[r,47]=='3') {rec123[r,47]<-'1'} }

#first rename variables for modeling
rec$yi<-rec$es
rec$vi<-rec$se

res <- rma.mv(yi, V=vi, 
              W=weights,
              random = ~ 1 | Study_ID/taxa/Response.measured, #method = "ML",
              data=rec)
res

#########################################
#############################################
#######################################

#### forest plot for recreational activities ####
# 4 shore-open-gradient categories

#categories
rec.shore <- rec[(rec$Rec.a =='Shore use'),]
rec.angl <- rec[(rec$Rec.a =='Angling'),]
rec.swim <- rec[(rec$Rec.a =='Swimming'),]
rec.boat <- rec[(rec$Rec.a =='Boating'),]

#model for each activity
m.rec.shore<- rma.mv(yi, vi, W=weights, 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     data=rec.shore)

m.rec.angl<- rma.mv(yi, vi, W=weights, 
                    random = ~ 1 | Study_ID/taxa/Response.measured, 
                    data=rec.angl)

m.rec.swim<- rma.mv(yi, vi, W=weights, 
                    random = ~ 1 | Study_ID/taxa/Response.measured, 
                    data=rec.swim)

m.rec.boat<- rma.mv(yi, vi, W=weights, 
                    random = ~ 1 | Study_ID/taxa/Response.measured, 
                    data=rec.boat)

#Get data#

Mean<-c(           
  m.rec.shore[["beta"]],
  m.rec.angl[["beta"]],
  m.rec.swim[["beta"]],
  m.rec.boat[["beta"]]
)
Lower <-c(
  m.rec.shore[["ci.lb"]],
  m.rec.angl[["ci.lb"]],
  m.rec.swim[["ci.lb"]],
  m.rec.boat[["ci.lb"]]
)
Upper <-c(
  m.rec.shore[["ci.ub"]],
  m.rec.angl[["ci.ub"]],
  m.rec.swim[["ci.ub"]],
  m.rec.boat[["ci.ub"]]
)
N<-c(
  m.rec.shore[["k.eff"]],
  m.rec.angl[["k.eff"]],
  m.rec.swim[["k.eff"]],
  m.rec.boat[["k.eff"]]
)
Nst<-c(
  m.rec.shore$s.nlevels[1],
  m.rec.angl$s.nlevels[1],
  m.rec.swim$s.nlevels[1],
  m.rec.boat$s.nlevels[1]
)

# categories
Activity<-c(
  "Shore use",
  "Angling",
  "Swimming",
  "Boating"
)

values<- cbind(Mean, Lower,Upper)
values<-format(round(values, 2), nsmall = 2)
CI<- cbind(Lower,Upper)
CI<-format(round(CI, 2), nsmall = 2)
CI <- apply(CI,1,function(x){
  paste0("[",paste(x, collapse="; "),"]")
})

Nk <- cbind(Nst,N)
Nk <- apply(Nk,1,function(x){
  paste0(paste(x, collapse="("),")")
})

Mean<-format(round(Mean, 2), nsmall = 2)


tabletext<-cbind(Activity, 
                 Mean, CI,
                 Nk)


### forest plot ###
library(lattice)
trellis.device(device="windows", height = 3, width = 8, color=TRUE)

forestplot(rbind(c("Activity", "Mean", "95% CI", "N(k)"),tabletext),
           graph.pos = 2,
           mean = c(NA,Mean),
           lower = c(NA,Lower),
           upper = c(NA,Upper),
           boxsize=0.3,
           #col=fpColors(box=c("brown","darkgreen","blue","darkblue"),
           #            line=c("brown","darkgreen","blue","darkblue")),
           xlab="Hedges G",
           #clip=c(-3.5,2),
           xticks = c(  -1, -0.5,0, 0.5),
           new_page=TRUE,
           #legend = c("Shore use","Angling","Swimming","Boating"),
           #legend_args = fpLegend(pos = list(x=1.05, y=0.625), 
           #gp=gpar(col="black", fill="#F9F9F9")),
           hrzl_lines = list("2" = gpar(lwd=2, col = "black"),
                             "6" = gpar(lwd=2, col = "Black")),
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=1.5),
                            ticks=gpar(cex=1.2), 
                            xlab=gpar(cex=1.2)),
           fn.ci_norm = fpDrawCircleCI)

#########################################
#############################################
#######################################

#### Real control forestplot ####

rec1 <- rec[(rec$Real.control.............== '1'),]
### data only high vs low ###
rec2 <- rec[(rec$Real.control.............== '2'),]
### data only not real BA ###
rec3 <- rec[(rec$Real.control.............== '3'),]

rec.real<-rbind(rec1,rec3)
rec.lvh <- rec2

#categories
rec.real.shore<- rec.real[(rec.real$Rec.a =='Shore use'),]
rec.real.angl <- rec.real[(rec.real$Rec.a =='Angling'),]
rec.real.swim <- rec.real[(rec.real$Rec.a =='Swimming'),]
rec.real.boat <- rec.real[(rec.real$Rec.a =='Boating'),]

#categories
rec.lvh.shore<- rec.lvh[(rec.lvh$Rec.a =='Shore use'),]
rec.lvh.angl <- rec.lvh[(rec.lvh$Rec.a =='Angling'),]
rec.lvh.swim <- rec.lvh[(rec.lvh$Rec.a =='Swimming'),]
rec.lvh.boat <- rec.lvh[(rec.lvh$Rec.a =='Boating'),]

#seperate models

m.rec.real.shore<- rma.mv(yi, vi,W=weights,
                          random = ~ 1 | Study_ID/taxa/Response.measured, 
                          data=rec.real.shore)
m.rec.lvh.shore<- rma.mv(yi, vi, W=weights,
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         data=rec.lvh.shore)
m.rec.real.angl<- rma.mv(yi, vi,W=weights,
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         data=rec.real.angl)
m.rec.lvh.angl<- rma.mv(yi, vi, W=weights,
                        random = ~ 1 | Study_ID/taxa/Response.measured, 
                        data=rec.lvh.angl)
m.rec.real.swim<- rma.mv(yi, vi,W=weights,
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         data=rec.real.swim)
m.rec.lvh.swim<- rma.mv(yi, vi, W=weights,
                        random = ~ 1 | Study_ID/taxa/Response.measured, 
                        data=rec.lvh.swim)
m.rec.real.boat<- rma.mv(yi, vi,W=weights,
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         data=rec.real.boat)
m.rec.lvh.boat<- rma.mv(yi, vi, W=weights,
                        random = ~ 1 | Study_ID/taxa/Response.measured, 
                        data=rec.lvh.boat)

### get data from models ####
#and put them in table

Act<-c("Shore use", "",
       "Angling", "",
       "Swimming", "",
       "Boating", "")

Comparison<-c(
  "Control vs Impact",
  "Low vs high Impact",
  "Control vs Impact",
  "Low vs high Impact",
  "Control vs Impact",
  "Low vs high Impact",
  "Control vs Impact",
  "Low vs high Impact")
Mean <- c(
  m.rec.real.shore[["beta"]],
  m.rec.lvh.shore[["beta"]],
  m.rec.real.angl[["beta"]],
  m.rec.lvh.angl[["beta"]],
  m.rec.real.swim[["beta"]],
  m.rec.lvh.swim[["beta"]],
  m.rec.real.boat[["beta"]],
  m.rec.lvh.boat[["beta"]])

Lower <-c(
  m.rec.real.shore[["ci.lb"]],
  m.rec.lvh.shore[["ci.lb"]],
  m.rec.real.angl[["ci.lb"]],
  m.rec.lvh.angl[["ci.lb"]],
  m.rec.real.swim[["ci.lb"]],
  m.rec.lvh.swim[["ci.lb"]],
  m.rec.real.boat[["ci.lb"]],
  m.rec.lvh.boat[["ci.lb"]])
Upper<-c(
  m.rec.real.shore[["ci.ub"]],
  m.rec.lvh.shore[["ci.ub"]],
  m.rec.real.angl[["ci.ub"]],
  m.rec.lvh.angl[["ci.ub"]],
  m.rec.real.swim[["ci.ub"]],
  m.rec.lvh.swim[["ci.ub"]],
  m.rec.real.boat[["ci.ub"]],
  m.rec.lvh.boat[["ci.ub"]])
N<-c(
  m.rec.real.shore[["k.eff"]],
  m.rec.lvh.shore[["k.eff"]],
  m.rec.real.angl[["k.eff"]],
  m.rec.lvh.angl[["k.eff"]],
  m.rec.real.swim[["k.eff"]],
  m.rec.lvh.swim[["k.eff"]],
  m.rec.real.boat[["k.eff"]],
  m.rec.lvh.boat[["k.eff"]])
Nst<-c(
  m.rec.real.shore$s.nlevels[1],
  m.rec.lvh.shore$s.nlevels[1],
  m.rec.real.angl$s.nlevels[1],
  m.rec.lvh.angl$s.nlevels[1],
  m.rec.real.swim$s.nlevels[1],
  m.rec.lvh.swim$s.nlevels[1],
  m.rec.real.boat$s.nlevels[1],
  m.rec.lvh.boat$s.nlevels[1])
#separate
Mean.real <- c(
  m.rec.real.shore[["beta"]],
  m.rec.real.angl[["beta"]],
  m.rec.real.swim[["beta"]],
  m.rec.real.boat[["beta"]])
Mean.lvh<- c(
  m.rec.lvh.shore[["beta"]],
  m.rec.lvh.angl[["beta"]],
  m.rec.lvh.swim[["beta"]],
  m.rec.lvh.boat[["beta"]])
Lower.real <-c(
  m.rec.real.shore[["ci.lb"]],
  m.rec.real.angl[["ci.lb"]],
  m.rec.real.swim[["ci.lb"]],
  m.rec.real.boat[["ci.lb"]])
Lower.lvh <-c(
  m.rec.lvh.shore[["ci.lb"]],
  m.rec.lvh.angl[["ci.lb"]],
  m.rec.lvh.swim[["ci.lb"]],
  m.rec.lvh.boat[["ci.lb"]])
Upper.real<-c(
  m.rec.real.shore[["ci.ub"]],
  m.rec.real.angl[["ci.ub"]],
  m.rec.real.swim[["ci.ub"]],
  m.rec.real.boat[["ci.ub"]])
Upper.lvh<-c(
  m.rec.lvh.shore[["ci.ub"]],
  m.rec.lvh.angl[["ci.ub"]],
  m.rec.lvh.swim[["ci.ub"]],
  m.rec.lvh.boat[["ci.ub"]])

Activ<-c("Shore use","Angling","Swimming","Boating")
Empty1<-c(
  "                    ",
  "                    ",
  "                    ",
  "                    "
)
Empty2<-c(
  "                   ",
  "                   ",
  "                   ",
  "                   "
)

Activity<-cbind(Activ,Empty1,Empty2)

CI<- cbind(Lower,Upper)
CI<-format(round(CI, 2), nsmall = 2)
CI <- apply(CI,1,function(x){
  paste0("[",paste(x, collapse="; "),"]")
})
Nk <- cbind(Nst,N)
Nk <- apply(Nk,1,function(x){
  paste0(paste(x, collapse="("),")")
})

Mean<-format(round(Mean, 2), nsmall = 2)

tabletext<-cbind(Act, Comparison, 
                 Mean, CI,
                 Nk)
tabletext2<-cbind(Act, 
                  Mean, CI,
                  Nk)

#### Multi-line real control forestplot ####
trellis.device(device="windows", height = 4, width = 9, color=TRUE)

forestplot(Activity,
           graph.pos = 2,
           mean = cbind(Mean.real, Mean.lvh),
           lower = cbind(Lower.real,Lower.lvh),
           upper = cbind(Upper.real,Upper.lvh),
           boxsize=0.17,
           col=fpColors(box=c("black","grey50"),
                        line=c("black","grey50")),
           xlab="Hedges G",
           #clip=c(-3.5,2),
           xticks = c( -1, -0.5,0, 0.5),
           new_page=TRUE,
           legend = c("Control vs. Impact",
                      "Low vs. high Impact"),
           legend_args = fpLegend(pos = list(x=1.7, y=-0.13)#, 
                                  #gp=gpar(col="black", fill="#F9F9F9")
           ),
           hrzl_lines = list("1" = gpar(lwd=2, col = "black"),
                             #"2" = gpar(lwd=1, col = "grey"),
                             #"3" = gpar(lwd=1, col = "grey"),
                             #"4" = gpar(lwd=1, col = "grey"),
                             "5" = gpar(lwd=2, col = "Black")),
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=1.5),
                            ticks=gpar(cex=1.2), 
                            xlab=gpar(cex=1.2)),
           #fn.ci_norm = fpDrawCircleCI,
           fn.ci_norm=list(fpDrawCircleCI,fpDrawDiamondCI))

#########################################
#############################################
#######################################

#### levels forest plot ####

#### Models ####

### subsample only level 1 of biological organization ###
rec.ind <- rec[(rec$Level.of.biological.organization==1),]#level 1 Individual

m.rec.ind.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.ind)
m.rec.ind.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                       random = ~ 1 | Study_ID/taxa/Response.measured, 
                       W=weights, data=rec.ind)
m.rec.ind.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.ind)
m.rec.ind.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.ind)

### subsample only level 2 of biological organization ###
rec.pop <- rec[(rec$Level.of.biological.organization==2),]#level 2 Population

m.rec.pop.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.pop)
m.rec.pop.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                       random = ~ 1 | Study_ID/taxa/Response.measured, 
                       W=weights, data=rec.pop)
m.rec.pop.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.pop)
m.rec.pop.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.pop)

### subsample only level 3 of biological organization ###
rec.bio <- rec[(rec$Level.of.biological.organization==3),]#level 3 Community

m.rec.bio.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.bio)
m.rec.bio.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                       random = ~ 1 | Study_ID/taxa/Response.measured, 
                       W=weights, data=rec.bio)
m.rec.bio.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.bio)
m.rec.bio.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.bio)

### subsample only level 4 of biological organization ###
rec.eco <- rec[(rec$Level.of.biological.organization==4),]#level 4 Ecosystem

m.rec.eco.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.eco)
m.rec.eco.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                       random = ~ 1 | Study_ID/taxa/Response.measured, 
                       W=weights, data=rec.eco)
m.rec.eco.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.eco)
m.rec.eco.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.eco)

## get data from models ####
#and put them together in a table

Level<-c(
  "Individual",
  "Population",
  "Community",
  "Environment"
)

Empty1<-c(
  "   ",
  "   ",
  "   ",
  "   "
)
Empty2<-c(
  "         ",
  "         ",
  "         ",
  "         "
)

Level<-cbind(Level,Empty2,Empty2, Empty2)

Mean.a <- c(
  m.rec.ind.a[["beta"]],
  m.rec.pop.a[["beta"]],
  m.rec.bio.a[["beta"]],
  m.rec.eco.a[["beta"]])
Mean.ab<-c(
  m.rec.ind.ab[["beta"]],
  m.rec.pop.ab[["beta"]],
  m.rec.bio.ab[["beta"]],
  m.rec.eco.ab[["beta"]])
Mean.b<-c(
  m.rec.ind.b[["beta"]],
  m.rec.pop.b[["beta"]],
  m.rec.bio.b[["beta"]],
  m.rec.eco.b[["beta"]]
)
Mean.h<-c(           
  m.rec.ind.h[["beta"]],
  m.rec.pop.h[["beta"]],
  m.rec.bio.h[["beta"]],
  m.rec.eco.h[["beta"]]
)
Lower.a <-c(
  m.rec.ind.a[["ci.lb"]],
  m.rec.pop.a[["ci.lb"]],
  m.rec.bio.a[["ci.lb"]],
  m.rec.eco.a[["ci.lb"]])
Lower.ab <-c(           
  m.rec.ind.ab[["ci.lb"]],
  m.rec.pop.ab[["ci.lb"]],
  m.rec.bio.ab[["ci.lb"]],
  m.rec.eco.ab[["ci.lb"]])
Lower.b <-c(
  m.rec.ind.b[["ci.lb"]],
  m.rec.pop.b[["ci.lb"]],
  m.rec.bio.b[["ci.lb"]],
  m.rec.eco.b[["ci.lb"]]
)
Lower.h <-c(
  m.rec.ind.h[["ci.lb"]],
  m.rec.pop.h[["ci.lb"]],
  m.rec.bio.h[["ci.lb"]],
  m.rec.eco.h[["ci.lb"]]
)
Upper.a <-c(
  m.rec.ind.a[["ci.ub"]],
  m.rec.pop.a[["ci.ub"]],
  m.rec.bio.a[["ci.ub"]],
  m.rec.eco.a[["ci.ub"]]
)
Upper.ab <-c(
  m.rec.ind.ab[["ci.ub"]],
  m.rec.pop.ab[["ci.ub"]],
  m.rec.bio.ab[["ci.ub"]],
  m.rec.eco.ab[["ci.ub"]])
Upper.b <-c(
  m.rec.ind.b[["ci.ub"]],
  m.rec.pop.b[["ci.ub"]],
  m.rec.bio.b[["ci.ub"]],
  m.rec.eco.b[["ci.ub"]]
)
Upper.h <-c(
  m.rec.ind.h[["ci.ub"]],
  m.rec.pop.h[["ci.ub"]],
  m.rec.bio.h[["ci.ub"]],
  m.rec.eco.h[["ci.ub"]]
)

N.a<-c(
  m.rec.ind.a[["k.eff"]],
  m.rec.pop.a[["k.eff"]],
  m.rec.bio.a[["k.eff"]],
  m.rec.eco.a[["k.eff"]]
)
N.ab<-c(
  m.rec.ind.ab[["k.eff"]],
  m.rec.pop.ab[["k.eff"]],
  m.rec.bio.ab[["k.eff"]],
  m.rec.eco.ab[["k.eff"]]
)
N.b<-c(
  m.rec.ind.b[["k.eff"]],
  m.rec.pop.b[["k.eff"]],
  m.rec.bio.b[["k.eff"]],
  m.rec.eco.b[["k.eff"]]
)
N.h<-c(
  m.rec.ind.h[["k.eff"]],
  m.rec.pop.h[["k.eff"]],
  m.rec.bio.h[["k.eff"]],
  m.rec.eco.h[["k.eff"]]
)


Nst.a<-c(
  m.rec.ind.a$s.nlevels[1],
  m.rec.pop.a$s.nlevels[1],
  m.rec.bio.a$s.nlevels[1],
  m.rec.eco.a$s.nlevels[1]
)
Nst.ab<-c(
  m.rec.ind.ab$s.nlevels[1],
  m.rec.pop.ab$s.nlevels[1],
  m.rec.bio.ab$s.nlevels[1],
  m.rec.eco.ab$s.nlevels[1]
)
Nst.b<-c(
  m.rec.ind.b$s.nlevels[1],
  m.rec.pop.b$s.nlevels[1],
  m.rec.bio.b$s.nlevels[1],
  m.rec.eco.b$s.nlevels[1]
)
Nst.h<-c(
  m.rec.ind.h$s.nlevels[1],
  m.rec.pop.h$s.nlevels[1],
  m.rec.bio.h$s.nlevels[1],
  m.rec.eco.h$s.nlevels[1]
)
Levels<-c("Individual","","","",
          "Population","","","",
          "Community","","","",
          "Environment","","","")
Activity<-c("Shore use","Angling","Swimming","Boating",
            "Shore use","Angling","Swimming","Boating",
            "Shore use","Angling","Swimming","Boating",
            "Shore use","Angling","Swimming","Boating")
Mean<-c(Mean.a, Mean.ab, Mean.b, Mean.h)
Lower<-c(Lower.a, Lower.ab, Lower.b, Lower.h)
Upper<-c(Upper.a, Upper.ab, Upper.b, Upper.h)
N<-c(N.a, N.ab, N.b, N.h)
Nst<-c(Nst.a, Nst.ab, Nst.b, Nst.h)

CI<- cbind(Lower,Upper)
CI<-format(round(CI, 2), nsmall = 2)
CI <- apply(CI,1,function(x){
  paste0("[",paste(x, collapse="; "),"]")
})
Nk <- cbind(Nst,N)
Nk <- apply(Nk,1,function(x){
  paste0(paste(x, collapse="("),")")
})

Mean<-format(round(Mean, 2), nsmall = 2)


tabletext<-cbind(Levels, Activity, 
                 Mean, CI,
                 Nk)


tabletext2<-cbind(Level,N.a, N.ab, N.b, N.h,
                  Nst.a, Nst.ab, Nst.b, Nst.h)

#### Multi-line Level forestplot ####

library(lattice)
trellis.device(device="windows", height = 8, width = 9, color=TRUE)

forestplot(Level,
           graph.pos = 2,
           mean = cbind(Mean.a, Mean.ab, Mean.b, Mean.h),
           lower = cbind(Lower.a,Lower.ab,Lower.b,Lower.h),
           upper = cbind(Upper.a,Upper.ab,Upper.b,Upper.h),
           boxsize=0.1,
           col=fpColors(box=c("brown","darkgreen","blue","darkblue"),
                        line=c("brown","darkgreen","blue","darkblue")),
           #txt_gp = fpTxtGp(ticks=gpar(cex=4)),
           xlab="Hedges G", 
           #clip=c(-3.5,2),
           xticks = c( -3,-2, -1, 0, 1),
           new_page=TRUE,
           legend = c("Shore use","Angling","Swimming","Boating"),
           legend_args = fpLegend(pos = list(x=1.2, y=0.61), 
                                  gp=gpar(col="black", fill="white")),
           hrzl_lines = list("1" = gpar(lwd=2, col = "black"),
                             "2" = gpar(lwd=1, col = "grey"),
                             "3" = gpar(lwd=1, col = "grey"),
                             "4" = gpar(lwd=2, col = "Black"),
                             "5" = gpar(lwd=2, col = "Black")),
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=1.5),
                            ticks=gpar(cex=1.5), 
                            xlab=gpar(cex=1.5)),
           fn.ci_norm = fpDrawCircleCI)


#########################################
#############################################
#######################################

#### Taxa forest plot ####

#models
### subsample only taxa 1 ###
rec.t1 <- rec[(rec$Group.==1),]#taxa 1 Invertebrates

m.rec.t1.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t1)
m.rec.t1.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.t1)
m.rec.t1.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t1)
m.rec.t1.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t1)

### subsample only taxa 2 ###
rec.t2 <- rec[(rec$Group.==2),]#taxa 2

m.rec.t2.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t2)
m.rec.t2.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.t2)
m.rec.t2.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t2)
m.rec.t2.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t2)

### subsample only taxa 3 ###
rec.t3 <- rec[(rec$Group.==3),]#taxa 3

m.rec.t3.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t3)
m.rec.t3.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.t3)
m.rec.t3.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t3)
m.rec.t3.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t3)

### subsample taxa 4 ###
rec.t4 <- rec[(rec$Group.==4),]#taxa 4

m.rec.t4.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t4)
m.rec.t4.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                      #random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.t4)
m.rec.t4.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t4)
m.rec.t4.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t4)

### subsample taxa 5 ###
rec.t5 <- rec[(rec$Group.==5),]#taxa 5

m.rec.t5.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, digits = 2, 
                     W=weights, data=rec.t5)
m.rec.t5.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.t5)
m.rec.t5.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t5)
m.rec.t5.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t5)

### subsample only taxa 6 ###
rec.t6 <- rec[(rec$Group.==6),]#taxa 6

m.rec.t6.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t6)
m.rec.t6.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.t6)
m.rec.t6.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t6)
m.rec.t6.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t6)
### subsample only taxa 7 ###
rec.t7 <- rec[(rec$Group.==7),]#taxa 7

m.rec.t7.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t7)
m.rec.t7.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                      random = ~ 1 | Study_ID/taxa/Response.measured, 
                      W=weights, data=rec.t7)
m.rec.t7.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t7)
m.rec.t7.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     W=weights, data=rec.t7)

### subsample others ###
rec.t999<-rbind(rec.t999,rec.t9,rec.t8) #combine rows 
### subsample only Environment ###
rec.t8 <- rec[(rec$Group.==8),]#phytoplancton
rec.t9 <- rec[(rec$Group.==9),]#zooplancton
rec.t999 <- rec[(rec$Group.==999),]#other
rec.t999<-rbind(rec.t999,rec.t9,rec.t8) #combine rows

m.rec.t999.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                       random = ~ 1 | Study_ID/taxa/Response.measured, 
                       W=weights, data=rec.t999)
m.rec.t999.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                        random = ~ 1 | Study_ID/taxa/Response.measured, 
                        W=weights, data=rec.t999)
m.rec.t999.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                       random = ~ 1 | Study_ID/taxa/Response.measured, 
                       W=weights, data=rec.t999)
m.rec.t999.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                       random = ~ 1 | Study_ID/taxa/Response.measured, 
                       W=weights, data=rec.t999)

###get data ####
Taxa2<-c(
  "Invertebrates",
  "Fish",
  "Amphibians",
  "Reptiles",
  "Birds",
  "Plants"
)# tabletext for plot
Empty1<-c(
  "                    ",
  "                    ",
  "                    ",
  "                    ",
  "                    ",
  "                    "
)
Empty2<-c(
  "                   ",
  "                   ",
  "                   ",
  "                   ",
  "                   ",
  "                   "
)# empty space to add Sample sizes manually in Powerpoint
Taxa2<-cbind(Taxa2,Empty1,Empty2)


count.all.ID<- count(rec,Rec.a,taxa,Study_ID)
count.all<- count(rec,Rec.a,taxa)
count.all$N_ES <- count.all$n
count.all.studies<-count(count.all.ID, Rec.a,taxa)
count.all.studies$N_Studies<-count.all.studies$n
count.all.both<-cbind(count.all,count.all.studies)
count.all.= subset(count.all.both, select = c( Rec.a,taxa,N_Studies, N_ES) )
names(count.all.both)

count.all.[order(count.all.$Rec.a,count.all.$taxa),]


## values from models
Mean.a <- c(
  m.rec.t1.a[["beta"]],
  NA,#m.rec.t2.a[["beta"]],
  NA,#m.rec.t3.a[["beta"]],
  -6.68,#m.rec.t4.a[["beta"]],
  m.rec.t5.a[["beta"]],
  m.rec.t7.a[["beta"]])
Mean.ab<-c(
  m.rec.t1.ab[["beta"]],
  NA,#m.rec.t2.ab[["beta"]],
  m.rec.t3.ab[["beta"]],
  m.rec.t4.ab[["beta"]],
  m.rec.t5.ab[["beta"]],
  m.rec.t7.ab[["beta"]])
Mean.b<-c(
  m.rec.t1.b[["beta"]],
  m.rec.t2.b[["beta"]],
  NA,#m.rec.t3.b[["beta"]],
  NA,#m.rec.t4.b[["beta"]],
  NA,#m.rec.t5.b[["beta"]],
  NA#m.rec.t7.b[["beta"]],
)
Mean.h<-c(           
  m.rec.t1.h[["beta"]],
  m.rec.t2.h[["beta"]],
  NA,#m.rec.t3.h[["beta"]],
  m.rec.t4.h[["beta"]],
  m.rec.t5.h[["beta"]],
  m.rec.t7.h[["beta"]]
)
Lower.a <-c(
  m.rec.t1.a[["ci.lb"]],
  NA,#m.rec.t2.a[["ci.lb"]],
  -1.17,#m.rec.t3.a[["ci.lb"]],
  NA,#m.rec.t4.a[["ci.lb"]],
  m.rec.t5.a[["ci.lb"]],
  m.rec.t7.a[["ci.lb"]])
Lower.ab <-c(           
  m.rec.t1.ab[["ci.lb"]],
  NA,#m.rec.t2.ab[["ci.lb"]],
  m.rec.t3.ab[["ci.lb"]],
  m.rec.t4.ab[["ci.lb"]],
  m.rec.t5.ab[["ci.lb"]],
  m.rec.t7.ab[["ci.lb"]])
Lower.b <-c(
  m.rec.t1.b[["ci.lb"]],
  m.rec.t2.b[["ci.lb"]],
  NA,#m.rec.t3.b[["ci.lb"]],
  NA,#m.rec.t4.b[["ci.lb"]],
  NA,#m.rec.t5.b[["ci.lb"]],
  NA#m.rec.t7.b[["ci.lb"]]
)
Lower.h <-c(
  m.rec.t1.h[["ci.lb"]],
  m.rec.t2.h[["ci.lb"]],
  NA,#m.rec.t3.h[["ci.lb"]],
  m.rec.t4.h[["ci.lb"]],
  m.rec.t5.h[["ci.lb"]],
  m.rec.t7.h[["ci.lb"]])
Upper.a <-c(
  m.rec.t1.a[["ci.ub"]],
  NA,#m.rec.t2.a[["ci.ub"]],
  -0.60,#m.rec.t3.a[["ci.ub"]],
  NA,#m.rec.t4.a[["ci.ub"]],
  m.rec.t5.a[["ci.ub"]],
  m.rec.t7.a[["ci.ub"]])
Upper.ab <-c(
  m.rec.t1.ab[["ci.ub"]],
  NA,#m.rec.t2.ab[["ci.ub"]],
  m.rec.t3.ab[["ci.ub"]],
  m.rec.t4.ab[["ci.ub"]],
  m.rec.t5.ab[["ci.ub"]],
  m.rec.t7.ab[["ci.ub"]])
Upper.b <-c(
  m.rec.t1.b[["ci.ub"]],
  m.rec.t2.b[["ci.ub"]],
  NA,#m.rec.t3.b[["ci.ub"]],
  NA,#m.rec.t4.b[["ci.ub"]],
  NA,#m.rec.t5.b[["ci.ub"]],
  NA#m.rec.t7.b[["ci.ub"]]
)
Upper.h <-c(
  m.rec.t1.h[["ci.ub"]],
  m.rec.t2.h[["ci.ub"]],
  NA,#m.rec.t3.h[["ci.ub"]],
  m.rec.t4.h[["ci.ub"]],
  m.rec.t5.h[["ci.ub"]],
  m.rec.t7.h[["ci.ub"]]
)

N.a<-c(
  m.rec.t1.a[["k.eff"]],
  0,#m.rec.t2.a[["k.eff"]],
  m.rec.t3.a[["k.eff"]],
  1,#m.rec.t4.a[["k.eff"]],
  m.rec.t5.a[["k.eff"]],
  m.rec.t7.a[["k.eff"]])
N.ab<-c(
  m.rec.t1.ab[["k.eff"]],
  0,#m.rec.t2.ab[["k.eff"]],
  m.rec.t3.ab[["k.eff"]],
  m.rec.t4.ab[["k.eff"]],
  m.rec.t5.ab[["k.eff"]],
  m.rec.t7.ab[["k.eff"]])
N.b<-c(
  m.rec.t1.b[["k.eff"]],
  m.rec.t2.b[["k.eff"]],
  0,#m.rec.t3.b[["k.eff"]],
  0,#m.rec.t4.b[["k.eff"]],
  0,#m.rec.t5.b[["k.eff"]],
  0#m.rec.t7.b[["k.eff"]],
)
N.h<-c(
  m.rec.t1.h[["k.eff"]],
  m.rec.t2.h[["k.eff"]],
  0,#m.rec.t3.h[["k.eff"]],
  m.rec.t4.h[["k.eff"]],
  m.rec.t5.h[["k.eff"]],
  m.rec.t7.h[["k.eff"]]
)
Nst.a<-c(
  m.rec.t1.a$s.nlevels[1],
  0,#m.rec.t2.a$s.nlevels[1],
  m.rec.t3.a$s.nlevels[1],
  1,#m.rec.t4.a$s.nlevels[1],
  m.rec.t5.a$s.nlevels[1],
  m.rec.t7.a$s.nlevels[1])
Nst.ab<-c(
  m.rec.t1.ab$s.nlevels[1],
  0,#m.rec.t2.ab$s.nlevels[1],
  m.rec.t3.ab$s.nlevels[1],
  1,#m.rec.t4.ab$s.nlevels[1],
  m.rec.t5.ab$s.nlevels[1],
  m.rec.t7.ab$s.nlevels[1])
Nst.b<-c(
  m.rec.t1.b$s.nlevels[1],
  m.rec.t2.b$s.nlevels[1],
  0,#m.rec.t3.b$s.nlevels[1],
  0,#m.rec.t4.b$s.nlevels[1],
  0,#m.rec.t5.b$s.nlevels[1],
  0#m.rec.t7.b$s.nlevels[1],
)
Nst.h<-c(
  m.rec.t1.h$s.nlevels[1],
  m.rec.t2.h$s.nlevels[1],
  0,#m.rec.t3.h$s.nlevels[1],
  m.rec.t4.h$s.nlevels[1],
  m.rec.t5.h$s.nlevels[1],
  m.rec.t7.h$s.nlevels[1]
)

Taxa<-c(
  "Invertebrates",
  "Fish",
  "Amphibians",
  "Reptiles",
  "Birds",
  "Plants",
  "Invertebrates",
  "Fish",
  "Amphibians",
  "Reptiles",
  "Birds",
  "Plants",
  "Invertebrates",
  "Fish",
  "Amphibians",
  "Reptiles",
  "Birds",
  "Plants",
  "Invertebrates",
  "Fish",
  "Amphibians",
  "Reptiles",
  "Birds",
  "Plants")
Activity<-c("Shore use","","","","","",
            "Angling","","","","","",
            "Swimming","","","","","",
            "Boating","","","","",""
)
Mean<-c(Mean.a, Mean.ab, Mean.b, Mean.h)
Lower<-c(Lower.a, Lower.ab, Lower.b, Lower.h)
Upper<-c(Upper.a, Upper.ab, Upper.b, Upper.h)
N<-c(N.a, N.ab, N.b, N.h)
Nst<-c(Nst.a, Nst.ab, Nst.b, Nst.h)

CI<- cbind(Lower,Upper)
CI<-format(round(CI, 2), nsmall = 2)
CI <- apply(CI,1,function(x){
  paste0("[",paste(x, collapse="; "),"]")
})
Nk <- cbind(Nst,N)
Nk <- apply(Nk,1,function(x){
  paste0(paste(x, collapse="("),")")
})

Mean<-format(round(Mean, 2), nsmall = 2)


tabletext<-cbind(Taxa, Activity, 
                 Mean, CI,
                 Nk)



tabletext2<-cbind(Taxa2,N.a, N.ab, N.b, N.h,
                  Nst.a, Nst.ab, Nst.b, Nst.h)

#### Multi-line TAXA forestplot ####

library(lattice)
trellis.device(device="windows", height = 8, width = 9, color=TRUE)


forestplot(Taxa2,
           graph.pos = 2,
           mean = cbind(Mean.a, Mean.ab, Mean.b, Mean.h),
           lower = cbind(Lower.a,Lower.ab,Lower.b,Lower.h),
           upper = cbind(Upper.a,Upper.ab,Upper.b,Upper.h),
           boxsize=0.13,
           col=fpColors(box=c("brown","darkgreen","blue","darkblue"),
                        line=c("brown","darkgreen","blue","darkblue")),
           xlab="Hedges G",
           #
           clip=c(-3.7,1.05),
           xticks = c( -3,-2, -1, 0, 1),
           new_page=TRUE,
           legend = c("Shore use","Angling","Swimming","Boating"),
           legend_args = fpLegend(pos = list(x=1.2, y=0.75)),
           hrzl_lines = list("1" = gpar(lwd=2, col = "black"),
                             "2" = gpar(lwd=1, col = "grey"),
                             "3" = gpar(lwd=1, col = "grey"),
                             "4" = gpar(lwd=1, col = "grey"),
                             "5" = gpar(lwd=1, col = "grey"),
                             "6" = gpar(lwd=1, col = "grey"),
                             "7" = gpar(lwd=2, col = "black")),
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=1.5),
                            ticks=gpar(cex=1.5), 
                            xlab=gpar(cex=1.5),
                            legend = gpar(cex=1)),
           fn.ci_norm = fpDrawCircleCI)

#########################################
#############################################
#######################################


#### only birds forest plot ####

#subsample only birds
rec.t5 <- rec[(rec$Group.==5),]#taxa 5

#### Models ####

## subsample only level 1 of biological organization ###
rec.t5.ind <- rec.t5[(rec.t5$Level.of.biological.organization==1),]#level 1 Individual

m.rec.t5.ind.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.ind)
m.rec.t5.ind.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                          random = ~ 1 | Study_ID/taxa/Response.measured, 
                          W=weights, 
                          data=rec.t5.ind)
m.rec.t5.ind.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.ind)
m.rec.t5.ind.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.ind)

### subsample only level 2 of biological organization ###
rec.t5.pop <- rec.t5[(rec.t5$Level.of.biological.organization==2),]#level 2 Population

m.rec.t5.pop.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.pop)
m.rec.t5.pop.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                          random = ~ 1 | Study_ID/taxa/Response.measured, 
                          W=weights, 
                          data=rec.t5.pop)
m.rec.t5.pop.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.pop)
m.rec.t5.pop.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.pop)

### subsample only level 3 of biological organization ###
rec.t5.bio <- rec.t5[(rec.t5$Level.of.biological.organization==3),]#level 3 Community

m.rec.t5.bio.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.bio)
m.rec.t5.bio.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                          random = ~ 1 | Study_ID/taxa/Response.measured, 
                          W=weights, 
                          data=rec.t5.bio)
m.rec.t5.bio.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.bio)
m.rec.t5.bio.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.bio)

### subsample only level 4 of biological organization ###
rec.t5.eco <- rec.t5[(rec.t5$Level.of.biological.organization==4),]#level 4 Ecosystem

m.rec.t5.eco.a <- rma.mv(yi, vi,  subset = (Rec.a=="Shore use"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.eco)
m.rec.t5.eco.ab <- rma.mv(yi, vi,  subset = (Rec.a=="Angling"), 
                          random = ~ 1 | Study_ID/taxa/Response.measured, 
                          W=weights, 
                          data=rec.t5.eco)
m.rec.t5.eco.b <- rma.mv(yi, vi,  subset = (Rec.a=="Swimming"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.eco)
m.rec.t5.eco.h <- rma.mv(yi, vi,  subset = (Rec.a=="Boating"), 
                         random = ~ 1 | Study_ID/taxa/Response.measured, 
                         W=weights, 
                         data=rec.t5.eco)

## get data from models ####
Level<-c(
  "Individual birds",
  " Bird Populations",
  "Bird communities"
)

Empty1<-c(
  "                    ",
  "                    ",
  "                    "
)
Empty2<-c(
  "                   ",
  "                   ",
  "                   "
)

Level<-cbind(Level,Empty1,Empty2)

Mean.a <- c(
  m.rec.t5.ind.a[["beta"]],
  m.rec.t5.pop.a[["beta"]],
  NA#-0.47#m.rec.t5.bio.a[["beta"]]
)
Mean.ab<-c(
  m.rec.t5.ind.ab[["beta"]],
  m.rec.t5.pop.ab[["beta"]],
  m.rec.t5.bio.ab[["beta"]])
Mean.h<-c(           
  m.rec.t5.ind.h[["beta"]],
  5,#m.rec.t5.pop.h[["beta"]],
  m.rec.t5.bio.h[["beta"]]
)
Lower.a <-c(
  m.rec.t5.ind.a[["ci.lb"]],
  m.rec.t5.pop.a[["ci.lb"]],
  NA#m.rec.t5.bio.a[["ci.lb"]]
)
Lower.ab <-c(           
  m.rec.t5.ind.ab[["ci.lb"]],
  m.rec.t5.pop.ab[["ci.lb"]],
  m.rec.t5.bio.ab[["ci.lb"]])

Lower.h <-c(
  m.rec.t5.ind.h[["ci.lb"]],
  -0.3,#m.rec.t5.pop.h[["ci.lb"]],
  m.rec.t5.bio.h[["ci.lb"]]
)
Upper.a <-c(
  m.rec.t5.ind.a[["ci.ub"]],
  m.rec.t5.pop.a[["ci.ub"]],
  NA#m.rec.t5.bio.a[["ci.ub"]]
)
Upper.ab <-c(
  m.rec.t5.ind.ab[["ci.ub"]],
  m.rec.t5.pop.ab[["ci.ub"]],
  m.rec.t5.bio.ab[["ci.ub"]])

Upper.h <-c(
  m.rec.t5.ind.h[["ci.ub"]],
  -0.05,#m.rec.t5.pop.h[["ci.ub"]],
  m.rec.t5.bio.h[["ci.ub"]]
)

N.a<-c(
  m.rec.t5.ind.a[["k.eff"]],
  m.rec.t5.pop.a[["k.eff"]],
  1#m.rec.t5.bio.a[["k.eff"]]
)
N.ab<-c(
  m.rec.t5.ind.ab[["k.eff"]],
  m.rec.t5.pop.ab[["k.eff"]],
  m.rec.t5.bio.ab[["k.eff"]]
)

N.h<-c(
  m.rec.t5.ind.h[["k.eff"]],
  m.rec.t5.pop.h[["k.eff"]],
  m.rec.t5.bio.h[["k.eff"]]
)


Nst.a<-c(
  m.rec.t5.ind.a$s.nlevels[1],
  m.rec.t5.pop.a$s.nlevels[1],
  1#m.rec.t5.bio.a$s.nlevels[1]
)
Nst.ab<-c(
  m.rec.t5.ind.ab$s.nlevels[1],
  m.rec.t5.pop.ab$s.nlevels[1],
  m.rec.t5.bio.ab$s.nlevels[1]
)
Nst.h<-c(
  m.rec.t5.ind.h$s.nlevels[1],
  m.rec.t5.pop.h$s.nlevels[1],
  m.rec.t5.bio.h$s.nlevels[1]
)


Mean<-c(Mean.a, Mean.ab, Mean.h)
Lower<-c(Lower.a, Lower.ab, Lower.h)
Upper<-c(Upper.a, Upper.ab, Upper.h)
N<-c(N.a, N.ab, N.h)
Nst<-c(Nst.a, Nst.ab, Nst.h)

CI<- cbind(Lower,Upper)
CI<-format(round(CI, 2), nsmall = 2)
CI <- apply(CI,1,function(x){
  paste0("[",paste(x, collapse="; "),"]")
})
Nk <- cbind(Nst,N)
Nk <- apply(Nk,1,function(x){
  paste0(paste(x, collapse="("),")")
})

Mean<-format(round(Mean, 2), nsmall = 2)

Levels<-c(
  "Individual birds",
  " Bird Populations",
  "Bird communities",
  "Individual birds",
  " Bird Populations",
  "Bird communities",
  "Individual birds",
  " Bird Populations",
  "Bird communities"
)
Activity<-c("Shore use","","",
            "Angling","","",
            "Boating","",""
)

tabletext<-cbind(Levels, Activity, 
                 Mean, CI,
                 Nk)


#### Multi-line TAXA forestplot ####

library(lattice)
trellis.device(device="windows", height = 6, width = 9, color=TRUE)

forestplot(Level,
           graph.pos = 2,
           mean = cbind(Mean.a, Mean.ab,  Mean.h),
           lower = cbind(Lower.a,Lower.ab,Lower.h),
           upper = cbind(Upper.a,Upper.ab,Upper.h),
           boxsize=0.1,
           col=fpColors(box=c("brown","darkgreen","darkblue"),
                        line=c("brown","darkgreen","darkblue")),
           xlab="Hedges G",
           clip=c(-2.5,1),
           xticks = c( -2, -1, 0, 1),
           new_page=TRUE,
           legend = c("Shore use","Angling","Boating"),
           legend_args = fpLegend(pos = list(x=1.1, y=0.85), 
                                  gp=gpar(col="black")),
           hrzl_lines = list("1" = gpar(lwd=2, col = "black"),
                             "2" = gpar(lwd=1, col = "grey"),
                             "3" = gpar(lwd=1, col = "grey"),
                             #"4" = gpar(lwd=1, col = "grey"),
                             "4" = gpar(lwd=2, col = "Black")),
           txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=1.5),
                            ticks=gpar(cex=1.5), 
                            xlab=gpar(cex=1.5)),
           fn.ci_norm = fpDrawCircleCI)

################################
###################################
##################################

#### moderator analysis ####

### moderators ###

rec.copy<-rec

#get rid of NA
rec$Water.body[is.na(rec$Water.body)]<-'999'

#decode moderators for moderator analysis

for(r in 1:nrow(rec)) {
  if (rec[r,28]=='1') {rec[r,28]<-'marine'} 
  else if (rec[r,28]=='2') {rec[r,28]<-"freshwater"}
  else if (rec[r,28]=='21') {rec[r,28]<-"lentic"}
  else if (rec[r,28]=='211') {rec[r,28]<-"lake"}
  else if (rec[r,28]=='2111') {rec[r,28]<-"gravel pit lake"}
  else if (rec[r,28]=='212') {rec[r,28]<-"pond" }
  else if (rec[r,28]=='22') {rec[r,28]<-"lotic"}
  else if (rec[r,28]=='221') {rec[r,28]<-"stream"}
  else if (rec[r,28]=='222') {rec[r,28]<-"river"}
  else if (rec[r,28]=='999') {rec[r,28]<-"Other"}
}

count(rec, Water.body)

rec<-rec.copy
#get rid of NA
rec$Water.body[is.na(rec$Water.body)]<-'999'

for(r in 1:nrow(rec)) {
  if (rec[r,28]=='1') {rec[r,28]<-'marine'} 
  else if (rec[r,28]=='2') {rec[r,28]<-"freshwater"}
  else if (rec[r,28]=='21') {rec[r,28]<-"lentic"}
  else if (rec[r,28]=='211') {rec[r,28]<-"lentic"}
  else if (rec[r,28]=='2111') {rec[r,28]<-"lentic"}
  else if (rec[r,28]=='212') {rec[r,28]<-"lentic" }
  else if (rec[r,28]=='22') {rec[r,28]<-"lotic"}
  else if (rec[r,28]=='221') {rec[r,28]<-"lotic"}
  else if (rec[r,28]=='222') {rec[r,28]<-"lotic"}
  else if (rec[r,28]=='999') {rec[r,28]<-"Other"}
}

#subsample with only lentic lotic
rec.water.len <- rec[(rec$Water.body=="lentic"),]
rec.water.lot <- rec[(rec$Water.body=="lotic"),]

rec.water<-rbind(rec.water.len,rec.water.lot)
rec<-rec.water

#### model ####
res <- rma.mv(yi, V=vi, 
              W=weights,
              mod = ~factor(Water.body),
              random = ~ 1 | Study_ID/taxa/Response.measured, data=rec)
res

#### HABITAT ####

rec<-rec.copy
names(rec)

rec$Habitat[is.na(rec$Habitat)]<-'999'

for(r in 1:nrow(rec)) {
  if (rec[r,29]=='1') {rec[r,29]<-'beach'} 
  else if (rec[r,29]=='2') {rec[r,29]<-"coast"}
  else if (rec[r,29]=='3') {rec[r,29]<-"estuary"}
  else if (rec[r,29]=='4') {rec[r,29]<-"open water"}
  else if (rec[r,29]=='5') {rec[r,29]<-"tidelands"}
  else if (rec[r,29]=='6') {rec[r,29]<-"mudflats" }
  else if (rec[r,29]=='7') {rec[r,29]<-"sand dunes"}
  else if (rec[r,29]=='8') {rec[r,29]<-"herbaceous bank"}
  else if (rec[r,29]=='9') {rec[r,29]<-"woody bank"}
  else if (rec[r,29]=='10') {rec[r,29]<-"reeds"}
  else if (rec[r,29]=='11') {rec[r,29]<-"urban bank"}
  else if (rec[r,29]=='12') {rec[r,29]<-"undefined bank"}
  else if (rec[r,29]=='13') {rec[r,29]<-"benthos"}
  else if (rec[r,29]=='14') {rec[r,29]<-"moors"}
  else if (rec[r,29]=='999') {rec[r,29]<-"other"}
}

count(rec, Habitat)
c<- count(rec, Habitat, Study_ID)
count(c, Habitat)

rec.shore <- rec[(rec$Rec.a =='Shore use'),]
rec.angl <- rec[(rec$Rec.a =='Angling'),]
rec.swim <- rec[(rec$Rec.a =='Swimming'),]
rec.boat <- rec[(rec$Rec.a =='Boating'),]

c.shore<-count(rec.shore, Habitat, Study_ID)
c.shore<-count(c.shore, Habitat)
c.shore
c.angl<-count(rec.angl, Habitat, Study_ID)
c.angl<-count(c.angl, Habitat)
c.angl
c.swim<-count(rec.swim, Habitat, Study_ID)
c.swim<-count(c.swim, Habitat)
c.swim
c.boat<-count(rec.boat, Habitat, Study_ID)
c.boat<-count(c.boat, Habitat)
c.boat


rec<-rec.copy
rec$Habitat[is.na(rec$Habitat)]<-'999'

for(r in 1:nrow(rec)) {
  if (rec[r,29]=='1') {rec[r,29]<-'beach'} 
  else if (rec[r,29]=='2') {rec[r,29]<-"other"}
  else if (rec[r,29]=='3') {rec[r,29]<-"other"}
  else if (rec[r,29]=='4') {rec[r,29]<-"open water"}
  else if (rec[r,29]=='5') {rec[r,29]<-"other"}
  else if (rec[r,29]=='6') {rec[r,29]<-"other" }
  else if (rec[r,29]=='7') {rec[r,29]<-"beach"}
  else if (rec[r,29]=='8') {rec[r,29]<-"herbaceous bank"}
  else if (rec[r,29]=='9') {rec[r,29]<-"woody bank"}
  else if (rec[r,29]=='10') {rec[r,29]<-"reeds"}
  else if (rec[r,29]=='11') {rec[r,29]<-"undefined bank"}
  else if (rec[r,29]=='12') {rec[r,29]<-"undefined bank"}
  else if (rec[r,29]=='13') {rec[r,29]<-"benthos"}
  else if (rec[r,29]=='14') {rec[r,29]<-"other"}
  else if (rec[r,29]=='999') {rec[r,29]<-"other"}
}


rec<- rec[!(rec$Habitat=='other'),]#exclude category "other"

rec.shore <- rec[(rec$Rec.a =='Shore use'),]
rec.angl <- rec[(rec$Rec.a =='Angling'),]
rec.swim <- rec[(rec$Rec.a =='Swimming'),]
rec.boat <- rec[(rec$Rec.a =='Boating'),]

c.shore<-count(rec.shore, Habitat, Study_ID)
c.shore<-count(c.shore, Habitat)
c.shore # reeds has too small sample size
rec.shore <- rec.shore[!(rec.shore$Habitat=='reeds'),]#exclude category "reeds"

c.angl<-count(rec.angl, Habitat, Study_ID)
c.angl<-count(c.angl, Habitat)
c.angl# reeds and woody bank have only one study each
rec.angl <- rec.angl[!(rec.angl$Habitat=='reeds'),]#exclude category "reeds"
for(r in 1:nrow(rec.angl)) {
  if (rec.angl[r,29]=='woody bank') {rec.angl[r,29]<-'undefined bank'}# put woody bank into category undefined bank
}
c.swim<-count(rec.swim, Habitat, Study_ID)
c.swim<-count(c.swim, Habitat)
c.swim # only two categories and one of them only one Study , therefore no moderator analysis possible

c.boat<-count(rec.boat, Habitat, Study_ID)
c.boat<-count(c.boat, Habitat)
c.boat#beach, herbaceous bank and reeds only one study
rec.boat <- rec.boat[!(rec.boat$Habitat=='reeds'),]#exclude category "reeds"
for(r in 1:nrow(rec.boat)) {
  if (rec.boat[r,29]=='herbaceous bank') {rec.boat[r,29]<-'undefined bank'}# put woody bank into category undefined bank
  else if (rec.boat[r,29]=='beach') {rec.boat[r,29]<-"undefined bank"}
}




#### grey literature ###

rec<-rec.copy
rec$Publication.type.

rec$Publication.type.[is.na(rec$Publication.type.)]<-'999'

for(r in 1:nrow(rec)) {
  if (rec[r,8]=='1') {rec[r,8]<-'journal'} 
  else if (rec[r,8]=='2') {rec[r,8]<-"book"}
  else if (rec[r,8]=='3') {rec[r,8]<-"book chapter"}
  else if (rec[r,8]=='4') {rec[r,8]<-"report"}
  else if (rec[r,8]=='5') {rec[r,8]<-"conference proceedings"}
  else if (rec[r,8]=='6') {rec[r,8]<-"published other" }
  else if (rec[r,8]=='7') {rec[r,8]<-"unplublished thesis"}
  else if (rec[r,8]=='8') {rec[r,8]<-"unpublished other"}
}


count(rec, Publication.type.)
c<- count(rec, Publication.type., Study_ID)
count(c, Publication.type.)

rec.shore <- rec[(rec$Rec.a =='Shore use'),]
rec.angl <- rec[(rec$Rec.a =='Angling'),]
rec.swim <- rec[(rec$Rec.a =='Swimming'),]
rec.boat <- rec[(rec$Rec.a =='Boating'),]

c.shore<-count(rec.shore, Publication.type., Study_ID)
c.shore<-count(c.shore, Publication.type.)
c.shore
c.angl<-count(rec.angl, Publication.type., Study_ID)
c.angl<-count(c.angl, Publication.type.)
c.angl
c.swim<-count(rec.swim, Publication.type., Study_ID)
c.swim<-count(c.swim, Publication.type.)
c.swim
c.boat<-count(rec.boat, Publication.type., Study_ID)
c.boat<-count(c.boat, Publication.type.)
c.boat

### published vs unpublished

rec<-rec.copy
rec$Peer.reviewed

rec$Peer.reviewed[is.na(rec$Peer.reviewed)]<-'999'

for(r in 1:nrow(rec)) {
  if (rec[r,9]=='1') {rec[r,9]<-'no'} 
  else if (rec[r,9]=='2') {rec[r,9]<-"yes"}
}

rec<- rec[!(rec$Peer.reviewed=='999'),]#exclude category "other"

count(rec, Peer.reviewed)
c<- count(rec, Peer.reviewed, Study_ID)
count(c, Peer.reviewed)

rec.shore <- rec[(rec$Rec.a =='Shore use'),]
rec.angl <- rec[(rec$Rec.a =='Angling'),]
rec.swim <- rec[(rec$Rec.a =='Swimming'),]
rec.boat <- rec[(rec$Rec.a =='Boating'),]

c.shore<-count(rec.shore, Peer.reviewed, Study_ID)
c.shore<-count(c.shore, Peer.reviewed)
c.shore
c.angl<-count(rec.angl, Peer.reviewed, Study_ID)
c.angl<-count(c.angl, Peer.reviewed)
c.angl
c.swim<-count(rec.swim, Peer.reviewed, Study_ID)
c.swim<-count(c.swim, Peer.reviewed)
c.swim
c.boat<-count(rec.boat, Peer.reviewed, Study_ID)
c.boat<-count(c.boat, Peer.reviewed)
c.boat


#### scale as moderator ###

rec<-rec-copy

rec.shore <- rec[(rec$Rec.a =='Shore use'),]
rec.angl <- rec[(rec$Rec.a =='Angling'),]
rec.swim <- rec[(rec$Rec.a =='Swimming'),]
rec.boat <- rec[(rec$Rec.a =='Boating'),]

#FSN
FSN.swim <- fsn(yi, vi, data = rec.swim)
FSN.boat <- fsn(yi, vi, data = rec.boat)
FSN.shore <- fsn(yi, vi, data = rec.shore)
FSN.angl <- fsn(yi, vi, data = rec.angl)

FSN<-c(
  FSN.shore[["fsnum"]],
  FSN.angl[["fsnum"]],
  FSN.swim[["fsnum"]],
  FSN.boat[["fsnum"]]
)

###models####



m.rec.shore<- rma.mv(yi, vi, W=weights,
                     mod = ~factor(Scale),
                     random = ~ 1 | Study_ID/taxa/Response.measured, 
                     data=rec.shore)

m.rec.angl<- rma.mv(yi, vi, W=weights,
                    mod = ~factor(Scale), 
                    random = ~ 1 | Study_ID/taxa/Response.measured, 
                    data=rec.angl)

m.rec.swim<- rma.mv(yi, vi, W=weights,
                    mod = ~factor(Scale),
                    random = ~ 1 | Study_ID/taxa/Response.measured, #control=list(optimizer="optim"),
                    data=rec.swim)

m.rec.boat<- rma.mv(yi, vi, W=weights,
                    mod = ~factor(Scale),
                    random = ~ 1 | Study_ID/taxa/Response.measured, 
                    data=rec.boat)



#Get data#

QE<-c(           
  m.rec.shore[["QE"]],
  m.rec.angl[["QE"]],
  m.rec.swim[["QE"]],
  m.rec.boat[["QE"]]
)
QEp <-c(
  m.rec.shore[["QEp"]],
  m.rec.angl[["QEp"]],
  m.rec.swim[["QEp"]],
  m.rec.boat[["QEp"]]
)
QM <-c(
  m.rec.shore[["QM"]],
  m.rec.angl[["QM"]],
  m.rec.swim[["QM"]],
  m.rec.boat[["QM"]]
)
QMp<-c(
  m.rec.shore[["QMp"]],
  m.rec.angl[["QMp"]],
  m.rec.swim[["QMp"]],
  m.rec.boat[["QMp"]]
)
N<-c(
  m.rec.shore[["k.eff"]],
  m.rec.angl[["k.eff"]],
  m.rec.swim[["k.eff"]],
  m.rec.boat[["k.eff"]]
)
Nst<-c(
  m.rec.shore$s.nlevels[1],
  m.rec.angl$s.nlevels[1],
  m.rec.swim$s.nlevels[1],
  m.rec.boat$s.nlevels[1]
)

# categories
Activity<-c(
  "Shore use",
  "Angling",
  "Swimming",
  "Boating"
)

QE<-format(round(QE, 0), nsmall = 0)
QEp<-format(round(QEp, 2), nsmall = 2)
QM<-format(round(QM, 2), nsmall = 2)
QMp<-format(round(QMp, 2), nsmall = 2)

tabletext<-cbind(Activity,
                 N,Nst, 
                 QE, QEp, QM, QMp)


QE<- cbind(QE,QEp)
QE<-format(round(QE, 2), nsmall = 2)
QE <- apply(QE,1,function(x){
  paste0(paste(x, collapse=" (p = "),")")
})

QM<- cbind(QM,QMp)
QM<-format(round(QM, 2), nsmall = 2)
QM <- apply(QM,1,function(x){
  paste0(paste(x, collapse=" (p = "),")")
})

Nk <- cbind(Nst,N)
Nk <- apply(Nk,1,function(x){
  paste0(paste(x, collapse="("),")")
})


tabletext<-cbind(Activity, 
                 Nk, QE, QM, N)

tabletext
