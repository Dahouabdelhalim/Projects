# library(devtools)
# install_github("edseab/Counterfact")
library(Counterfact)
library(dplyr)

## RAW DATASETS ##
## ta <- ta_full <- read.csv("")  ## the full raw time allocation database
## census_ta_data<- read.csv("")  ## the contemporaneous census data 
codes <- read.csv("codes.csv")  ## the activity codes file
## com <- read.csv("")  ## data on communities
## r_matrix<- read.csv("") ## relatedness matrix for all communities
## reg <-read.csv("") ## full population register with data on parents
## med <- read.csv("") ## medical database with some additional data for cleaning and verification
## vis <- read.csv("") ## visit register with family ID and residence data

reg$pid<-as.character(reg$pid)
reg$father_pid<-as.character(reg$father_pid)
reg$mother_pid<-as.character(reg$mother_pid)

## Here we fun a family function which returns the relationship between any 2 individuals as a categorical variable (including affinal relationships) 

## First an "identical" function
i<-function(x,y){
  if(is.na(x)|is.na(y))return(FALSE)
  else if(x==y) return(TRUE)
  else return(FALSE)}

## Now the "family" function

family<-function (pid1,pid2, paste=F){
  
  
  pid1<-as.character(pid1)
  pid2<-as.character(pid2)
  
  a<-reg[match(pid1,reg$pid),]
  b<-reg[match(pid2,reg$pid),]
  c<-character(0)
  
  d1<-data.frame(pid=as.character(pid1))
  d2<-data.frame(pid=as.character(pid2))
  
  #parents
  d1$f<-a$father_pid
  d1$m<-a$mother_pid
  d2$f<-b$father_pid
  d2$m<-b$mother_pid
  #grandparents
  d1$pgf<-reg[match(d1$f, reg$pid),"father_pid"]
  d1$mgf<-reg[match(d1$m, reg$pid),"father_pid"]
  d1$pgm<-reg[match(d1$f, reg$pid),"mother_pid"]
  d1$mgm<-reg[match(d1$m, reg$pid),"mother_pid"]
  d2$pgf<-reg[match(d2$f, reg$pid),"father_pid"]
  d2$mgf<-reg[match(d2$m, reg$pid),"father_pid"]
  d2$pgm<-reg[match(d2$f, reg$pid),"mother_pid"]
  d2$mgm<-reg[match(d2$m, reg$pid),"mother_pid"]
  
  d1[is.na(d1)] <- "12345placeholder"
  d2[is.na(d2)] <- "67890placeholder"
  
  #siblings 
  as<-reg$pid[which((reg$father_pid %in% d1$f | reg$mother_pid %in% d1$m)& reg$pid!=pid1)]
  bs<-reg$pid[which((reg$father_pid %in% d2$f | reg$mother_pid %in% d2$m)& reg$pid!=pid2)]
  
  #joint children
  children<-intersect(reg[match(pid1,reg$father_pid),"pid"], reg[match(pid2,reg$mother_pid),"pid"])
  children<-append(children,intersect(reg[match(pid2,reg$father_pid),"pid"], reg[match(pid1,reg$mother_pid),"pid"]))
  
  
  #spouse
  
  aspouse<-rep(NA,nrow(a))
  aspouse <- ifelse(a$male==1,unique(reg$mother_pid[match(pid1,reg$father_pid)]),
                    unique(reg$father_pid[match(pid1,reg$mother_pid)]))
  
  bspouse<-rep(NA,nrow(b))
  bspouse<-ifelse(b$male==1,unique(reg$mother_pid[match(pid2,reg$father_pid)]),
                  unique(reg$father_pid[match(pid2,reg$mother_pid)]))
  
  aspouse[which(is.na(aspouse))] <- "aspouseplaceholder"	
  bspouse[which(is.na(bspouse))] <- "bspouseplaceholder"	
  
  #same person
  
  if(pid1==pid2){c<-"same person"} 
  
  if(pid1!=pid2){
    if(i(pid1,d2$f)) c<-append(c,"father")
    if(i(pid1,d2$m)) c<-append(c,"mother")
    if(i(pid2,d1$f) | i(pid2,d1$m)) c<-append(c,"child")
    if(i(pid1,d1$pgf)) c<-append(c,"paternal grandfather")
    if(i(pid1,d1$mgf)) c<-append(c,"maternal grandfather")
    if(i(pid1,d1$pgm)) c<-append(c,"paternal grandmother")
    if(i(pid1,d1$mgm)) c<-append(c,"maternal grandmother")
    if(i(pid2,d2$pgf) | i(pid2,d2$pgm))c<-append(c,"son's child")
    if(i(pid2,d2$mgf) | i(pid2,d2$mgm))c<-append(c,"daughter's child")
    if(i(d1$f,d2$f) & (i(d1$m,d2$m)==FALSE))c<-append(c,"paternal half-sibling")
    if(i(d1$m,d2$m) & (i(d1$f,d2$f)==FALSE))c<-append(c,"maternal half-sibling")
    if(i(d1$f,d2$f) & i(d1$m,d2$m)) c<-append(c,"full sibling")
    if(length(bspouse)!=0) {if(pid1 %in% reg$father_pid[match(bspouse,reg$pid)]==TRUE) c<-append(c,"father-in-law")
    if(pid1 %in% reg$mother_pid[match(bspouse,reg$pid)]==TRUE)c<-append(c,"mother-in-law")}
    if(!i(d1$f,d2$f)&!i(d1$m,d2$m) & (i(d1$pgf,d2$pgf) | i(d1$pgm,d2$pgm)))c<-append(c,"paternal parallel cousin")
    if(!i(d1$f,d2$f)&!i(d1$m,d2$m) & (i(d1$mgf,d2$mgf) | i(d1$mgm,d2$mgm)))c<-append(c,"maternal parallel cousin")
    if(!i(d1$f,d2$f)&!i(d1$m,d2$m) & (i(d1$pgf,d2$mgf) | i(d1$pgm,d2$mgm)))c<-append(c,"maternal cross cousin")
    if(!i(d1$f,d2$f)&!i(d1$m,d2$m) & (i(d1$mgf,d2$pgf) | i(d1$mgm,d2$pgm)))c<-append(c,"paternal cross cousin")
    if(length(bspouse)!=0){if(pid2 %in% reg$father_pid[match(aspouse,reg$pid)] == TRUE| pid2 %in% reg$mother_pid[match(aspouse,reg$pid)] == TRUE) c<-append(c,"offspring-in-law")
    if(pid1 %in% bspouse==TRUE) c<-append(c,"married")}
    if(length(bspouse)!=0 & length(as)!=0){if(bspouse %in% as )c<-append(c,"sibling's spouse")}
    if(length(aspouse)!=0 & length(bs)!=0){if(aspouse %in% bs )c<-append(c,"spouse's sibling")}
    if(length(aspouse)!=0 & length(d2$f)!=0){if(!i(pid1,d2$m) & d2$f %in% aspouse==TRUE)c<-append(c,"stepmother")}
    if(length(bspouse)!=0 & length(d1$f)!=0){if(!i(pid2,d1$m) & d1$f %in% bspouse==TRUE)c<-append(c,"stepchild")}
    if(length(aspouse)!=0 & length(bspouse)!=0) {if(length(na.omit(intersect(reg[which(reg$pid %in% aspouse==TRUE),"father_pid"],reg[which(reg$pid %in% bspouse==TRUE),"father_pid"])))!=0
                                                    | length(na.omit(intersect(reg[which(reg$pid %in% aspouse==TRUE),"mother_pid"],reg[which(reg$pid %in% bspouse==TRUE),"mother_pid"])))!=0)c<-append(c,"spouse's sibling's spouse")}
    if(length(bs)!=0){
      if(a$father_pid %in% bs | a$mother_pid %in% bs){
        if(is.na(a$male)){c<-append(c,"nibling")
        }else{if(a$male==1) c<-append(c,"nephew")
        if(a$male==0) c<-append(c,"niece")}
      }}
    
    
    if(length(as)!=0){
      if(b$father_pid %in% as | b$mother_pid %in% as){
        if(is.na(a$male)){c<-append(c,"auntle")
        } else {if(a$male==1) c<-append(c,"uncle")
        if(a$male==0) c<-append(c,"aunt")}}
    }
  }	
  if(length(c)==0) c<-"None"
  if(paste) c <- paste(c, collapse=";")
  return(c)
  end
  
}

## Here we clean the time allocation database
### Clean clusters, activities, family IDs ###

for(i in 1:nrow(ta)){
  
  if(is.na(ta$cluster[i])){
    if(length(unique(ta$cluster[which(ta$timeblock==ta$timeblock[i] & is.na(ta$cluster))]))==1){
      ta$cluster[i] <- unique(ta$cluster[which(ta$timeblock==ta$timeblock[i] & is.na(ta$cluster))])[1]
    }
  }
}
# Match descriptions
colnames(codes)[grep('code',colnames(codes))] <- "code"
act_codes<-subset(codes, codes$type=="activity")
ta$act1 <- toupper(ta$act1)
ta$cat1<-act_codes$category[match(ta$act1,act_codes$code)]
ta$desc1<-act_codes$desc[match(ta$act1,act_codes$code)]

ta$act2 <- toupper(ta$act2)
ta$cat2 <- act_codes$category[match(ta$act2,act_codes$code)]
ta$desc2<-act_codes$desc[match(ta$act2,act_codes$code)]

ta$cat1[is.na(ta$cat1)] <- "unspecified"
ta$cat2[is.na(ta$cat2)] <- "unspecified"


##Also clean unknown pids
ta <- ta[which(ta$pid %in% reg$pid),]
ta$timeblock <- paste(ta$fecha,ta$hora)
ta$sgroupID <- paste(ta$cluster,ta$timeblock,ta$sgroup)
ta$agroupID <- paste(ta$cluster,ta$timeblock,ta$agroup)
ta$male <- reg$male[match(ta$pid, reg$pid)]
ta$dob <- reg$date_of_birth[match(ta$pid, reg$pid)]

ta$age <- as.numeric(as.Date(ta$fecha) - as.Date(ta$dob))/365.25
ta$age[is.na(ta$age)] <- ta$AgeYr[is.na(ta$age)]
ta$age[which(ta$pid=="JALL")] <- 75
ta$age <- as.numeric(ta$age)

for(i in 1:nrow(ta)){
  if (ta$age[i]<0){
    other_dobs <- as.Date(med$FechNacim,'%m/%d/%Y')[which(med$pid==ta$pid[i])]
    ta$dob[i] <- min(other_dobs)
    ta$age[i] <- as.numeric(as.Date(ta$fecha[i]) - as.Date(ta$dob[i]))/365.25
    
  }
  progress(i, nrow(ta))
}

ta$agesexgrp <- NA
ta$agesexgrp[which(ta$age<5)] <- "Infant"
ta$agesexgrp[which(ta$age<14)] <- "Child"
ta$agesexgrp[which(ta$age>=14 & ta$male==0)] <- "Woman"
ta$agesexgrp[which(ta$age>=14 & ta$male==1)] <- "Man"
ta <- ta[!is.na(ta$age),]

### Create data frame of descriptive statistics for different activities ###
# New data frame with all we need 
activities <- data.frame(act=unique(ta$cat1))

for (i in 1:nrow(activities)){
  print(paste0(i,"/",nrow(activities)))
  act_i <- ta[which(ta$cat1==activities$act[i]),]
  act_f <- act_i[which(act_i$agesexgrp=="Woman"),] 
  act_m <- act_i[which(act_i$agesexgrp=="Man"),] 
  
  sgroup_act_i <- ta[which(ta$sgroupID %in% unique(act_i$sgroupID)),]
  sgroup_act_f <- sgroup_act_i[which(sgroup_act_i$agesexgrp=="Woman"),]
  sgroup_act_m <- sgroup_act_i[which(sgroup_act_i$agesexgrp=="Man"),]
  
  activities$grpsize_T[i] <- mean(sapply(split(sgroup_act_i, sgroup_act_i$timeblock),nrow))
  
  activities$prct_female_in_group[i] <- mean(sapply(split(sgroup_act_i, sgroup_act_i$timeblock),function(x)sum(x$male==0, na.rm=T)/nrow(x)))
  activities$prct_female_doing_act[i] <- mean(sapply(split(act_i, act_i$timeblock),function(x)sum(x$male==0, na.rm=T)/nrow(x)))
  
  activities$avg_age_sgroup_T[i] <- mean(sapply(split(sgroup_act_i$age, sgroup_act_i$pid), function(x)mean(x,na.rm=T)),na.rm=T)
  activities$avg_age_sgroup_F[i] <- mean(sapply(split(sgroup_act_f$age, sgroup_act_f$pid), mean))
  activities$avg_age_sgroup_M[i] <- mean(sapply(split(sgroup_act_m$age, sgroup_act_m$pid), mean))
  
  activities$avg_age_T[i] <- mean(sapply(split(act_i$age, act_i$pid), function(x)mean(x,na.rm=T)),na.rm=T)
  activities$avg_age_F[i] <- mean(sapply(split(act_f$age, act_f$pid), mean))
  activities$avg_age_M[i] <- mean(sapply(split(act_m$age, act_m$pid), mean))
  
  activities$prct_observations[i] <- nrow(act_i)/nrow(ta)
  activities$prct_obs_F[i] <- nrow(act_f)/length(which(ta$agesexgrp=="Woman"))
  activities$prct_obs_M[i] <- nrow(act_m)/length(which(ta$agesexgrp=="Man"))
  
}

#### CLEAN FAMILY IDS ####

# Start
ta$day <- as.Date(ta$fecha)
ta$timeblock <-paste(ta$fecha,ta$hora)
ta$hour <- substr(ta$timeblock,1,13)
ta$month <- substr(ta$day, 1,7)
ta$year <- substr(ta$day, 1,4)

# Remove visitors for now
sum(ta$pid=="VISITOR") #475
ta <- ta[which(ta$pid!="VISITOR"),]

# deal with repeated measures
ta<-ta[which(ta$pid=="VISITOR" | !duplicated(paste(ta$pid,ta$timeblock))),]

ta$contact_group <- paste(ta$sgroup, ta$timeblock)
ta <- ta[!is.na(ta$contact_group),]

# Fix family ids and community
for(i in 1:nrow(ta)){
  options <- vis[which(vis$pid==ta$pid[i] & !is.na(vis$family_id)),]
  if(nrow(options!=0)){
    ta$fid[i] <- options$family_id[which.min(abs(as.Date(options$date)-ta$day[i]))]
  }
  progress(i, nrow(ta))
}
ta$Comunidad <- vis$community[match(ta$ComID, vis$community_id)]
ta <- ta[!is.na(ta$sgroup),]

ta$fid <- ta$FamID
ta$fid[is.na(ta$fid)] <- ta$familiaID[is.na(ta$fid)]
ta$fid <- paste(ta$Comunidad,ta$fid, sep="_")


# update pids for activity objects
ta$obj1pid <- ta$pid[match(ta$obj1,ta$original_pid_time_allocation_file)]
ta$obj1pid[is.na(ta$obj1pid) & !is.na(ta$obj1)] <- ta$pid[match(ta$obj1[is.na(ta$obj1pid)& !is.na(ta$obj1)],ta$MidPId)]
ta$obj1pid[is.na(ta$obj1pid) & !is.na(ta$obj1)] <- ta_full$pid[match(ta$obj1[is.na(ta$obj1pid)& !is.na(ta$obj1)],ta_full$midpid_16)]

ta$iobj1pid <- ta$pid[match(ta$iobj1,ta$original_pid_time_allocation_file)]
ta$iobj1pid[is.na(ta$iobj1pid)& !is.na(ta$iobj1)] <- ta$pid[match(ta$iobj1[is.na(ta$iobj1pid)& !is.na(ta$iobj1)],ta$MidPId)]
ta$iobj1pid[is.na(ta$iobj1pid)& !is.na(ta$iobj1)] <- ta_full$pid[match(ta$iobj1[is.na(ta$iobj1pid)& !is.na(ta$iobj1)],ta_full$midpid_16)]

ta$obj2pid <- ta$pid[match(ta$obj2,ta$original_pid_time_allocation_file)]
ta$obj2pid[is.na(ta$obj2pid)& !is.na(ta$obj2)] <- ta$pid[match(ta$obj2[is.na(ta$obj2pid)& !is.na(ta$obj2)],ta$MidPId)]
ta$obj2pid[is.na(ta$obj2pid)& !is.na(ta$obj2)] <- ta_full$pid[match(ta$obj2[is.na(ta$obj2pid)& !is.na(ta$obj2)],ta_full$midpid_16)]

ta$iobj2pid <- ta$pid[match(ta$iobj2,ta$original_pid_time_allocation_file)]
ta$iobj2pid[is.na(ta$iobj2pid)& !is.na(ta$iobj2)] <- ta$pid[match(ta$iobj2[is.na(ta$iobj2pid)& !is.na(ta$iobj2)],ta$MidPId)]
ta$iobj2pid[is.na(ta$iobj2pid)& !is.na(ta$iobj2)] <- ta_full$pid[match(ta$iobj2[is.na(ta$iobj2pid)& !is.na(ta$iobj2)],ta_full$midpid_16)]


## Next there was some individual ID cleaning which we can't share for confidentiality purposes ##

# Now add object ages
ta$obj1age <- ta$age[match(paste(ta$obj1pid,ta$fecha),paste(ta$pid,ta$fecha))]
ta$iobj1age <- ta$age[match(paste(ta$iobj1pid,ta$fecha),paste(ta$pid,ta$fecha))]
ta$obj2age <- ta$age[match(paste(ta$obj2pid,ta$fecha),paste(ta$pid,ta$fecha))]
ta$iobj2age <- ta$age[match(paste(ta$iobj2pid,ta$fecha),paste(ta$pid,ta$fecha))]

### CHILDCARE ####

#Of the times they are observed at home, how often is their kid taken care of by somebody else?
#Unit of analysis is the kid of a particular woman
#Need to create db with timeblock,pid,carer, rel, cared for (yes/no), type of care

# Check if there are any mischaracterized childcare (eg adult tending to elderly person)

sum(ta$cat1=="childcare" & ta$age>=15 & rowSums(cbind(ta$obj1age <15, ta$iobj1age<15),na.rm=T)<1,na.rm=T)

# recode these as just care

ta$cat1[which(ta$cat1=="childcare" & ta$age>=15 & rowSums(cbind(ta$obj1age <15, ta$iobj1age<15),na.rm=T)<1)] <-"care"
ta$cat2[which(ta$cat2=="childcare" & ta$age>=15 & rowSums(cbind(ta$obj2age <15, ta$iobj2age<15),na.rm=T)<1)] <-"care"




#Create cat receiving vs giving childcare by choosing youngest as receiver of childcare
#First see there are any rows where object columns (iobj and obj) have 1 person younger and 1 person older than ego
sum((ta$age>ta$obj1age & ta$age<ta$iobjage)|(ta$age<ta$obj1age & ta$age>ta$iobjage),na.rm=T) #0 ok fine


ta$cat1[which(ta$cat1=="childcare" & (ta$age<ta$obj1age | ta$age<ta$iobj1age))] <- "receive childcare"
ta$cat2[which(ta$cat2=="childcare" & (ta$age<ta$obj2age | ta$age<ta$iobj2age))] <- "receive childcare"

#Only look at kids who are in the home cluster
ta$home_cluster <- famreg$cluster[match(ta$fid,famreg$fid)]

# parental pids
ta$m_pid <- reg$mother_pid[match(ta$pid, reg$pid)]
ta$f_pid <- reg$father_pid[match(ta$pid, reg$pid)]

## Ok now create analysis data frame for children
d <- ta[which(ta$age<=14 & !is.na(ta$home_cluster) & !is.na(ta$f_pid)),]

d$cared_for <- as.numeric(d$cat1=="receive childcare" | d$cat2=="receive childcare")

# Who are the carers?
d$n_carer <- d$carer1 <- d$carer2 <- d$carer3 <- d$carer4 <- NA
d$cared_for_social <- d$social_nonsib_allocare <- 0

for (i in 1:nrow(d)){
  sg <- ta[which(ta$sgroupID==d$sgroupID[i] & ta$age >=11),]
  
  if(nrow(sg)>0){
    social_carers <- sg$pid[which((sg$cat1=="social" & (sg$obj1pid==d$pid[i] | sg$iobj1pid==d$pid[i]))
                                  |(sg$cat2=="social" & (sg$obj2pid==d$pid[i] | sg$iobj2pid==d$pid[i])))]
    if(length(social_carers)>0){
      d$cared_for_social[i] <- 1
      social_carers_rel <- sapply(social_carers,family,d$pid[i])
      if(any(!social_carers_rel %in% c("mother","father","full sibling","maternal half-sibling"))){
        d$social_nonsib_allocare[i]<- 1
      }
    }
    
    if (d$cared_for[i]==1 | d$cared_for_social[i]==1){
      carers <- sg[which(sg$cat1=="childcare" | sg$cat2=="childcare"),]
      defcarers <- apply(carers[,c("obj1","iobj1","obj2","iobj2")],1, function(x)d$original_pid_time_allocation_file[i] %in% x)
      
      if(sum(defcarers)>0){
        all_carers <- c(carers$pid[which(defcarers)])
      }else all_carers <- c(carers$pid[!is.na(carers$pid)])
      d$n_carer[i] <- length(all_carers)
      d$carer1[i] <- all_carers[1]
      d$carer2[i] <- all_carers[2]
      d$carer3[i] <- all_carers[3]
      d$carer4[i] <- all_carers[4]
    }
  }
  progress(i,nrow(d))  
}

d$carer1rel <- mapply(function(x,y)ifelse(!is.na(x) & !is.na(y),paste(family(x,y),collapse="/"),NA),d$carer1,d$pid)
d$carer1male <- reg$male[match(d$carer1,reg$pid)]

#Repeat for other carers
d$carer2rel <- mapply(function(x,y)ifelse(!is.na(x) & !is.na(y),paste(family(x,y),collapse="/"),NA),d$carer2,d$pid)
d$carer2male <- reg$male[match(d$carer2,reg$pid)]

d$carer3rel <- mapply(function(x,y)ifelse(!is.na(x) & !is.na(y),paste(family(x,y),collapse="/"),NA),d$carer3,d$pid)
d$carer3male <- reg$male[match(d$carer3,reg$pid)]

d$carer4rel <- mapply(function(x,y)ifelse(!is.na(x) & !is.na(y),paste(family(x,y),collapse="/"),NA),d$carer4,d$pid)
d$carer4male <- reg$male[match(d$carer4,reg$pid)]


d$allocare <- apply(d[,grep("rel",colnames(d))],1,function(x)as.numeric(any(!is.na(x) & ! x %in% c("mother","father")))) # help w your kid regardless of whether mom is around
d$non_sibling_allocare<- apply(d[,grep("rel",colnames(d))],1,function(x)as.numeric(any(!is.na(x) & ! x %in% c("mother","father","full sibling","maternal half-sibling")))) # help w your kid not from your other kids

## Now calculate for each child whether they are in a social group that doesnt include their parents but includes another adult ## 
d$passive_allocare <- d$unsupervized <- d$passive_ac_sgroup <- d$passive_ac_agroup <- NA

for (i in 1:nrow(d)){
  other_sgroup_members <- ta[which(ta$sgroupID==d$sgroupID[i]),]
  other_sgroup_members$rel_to_child <- sapply(other_sgroup_members$pid,family,pid2=d$pid[i])
  
  other_agroup_members <- ta[which(ta$agroupID==d$agroupID[i]),]
  other_agroup_members$rel_to_child <- sapply(other_agroup_members$pid,family,pid2=d$pid[i])
  
  d$passive_ac_sgroup[i] <- as.numeric(d$cared_for[i]==0 & 
                                         !any(c("mother","father","full sibling","maternal half-sibling","paternal half_sibling") %in% 
                                                other_sgroup_members$rel_to_child[which(other_sgroup_members$age>=14)]) &
                                         any(other_sgroup_members$age>=14))
  
  
  d$passive_ac_agroup[i] <- as.numeric(d$cared_for[i]==0 & 
                                         !any(c("mother","father","full sibling","maternal half-sibling","paternal half_sibling") %in% 
                                                other_agroup_members$rel_to_child[which(other_agroup_members$age>=14)]) &
                                         any(other_agroup_members$age>=14))
  
  
  d$unsupervized[i] <- as.numeric(!any(other_sgroup_members$age>=14))
  
  if(d$passive_ac_sgroup[i]==1){
    d$passive_sgroup_carers[i] <- paste(other_sgroup_members$rel_to_child[which(other_sgroup_members$age>=14)], collapse="/")
  }
  
  if(d$passive_ac_agroup[i]==1){
    d$passive_agroup_carers[i] <- paste(other_agroup_members$rel_to_child[which(other_agroup_members$age>=14)], collapse="/")
  }  
  
  progress(i,nrow(d))
}

d$any_allocare <- as.numeric(d$allocare | d$passive_ac_sgroup | d$passive_ac_sgroup[i])


#### adding grand-parental ids ####
d$mgm_pid<-reg$mother_pid[match(d$m_pid, reg$pid)]
d$mgf_pid<-reg$father_pid[match(d$m_pid, reg$pid)]

d$pgm_pid<-reg$mother_pid[match(d$f_pid, reg$pid)]
d$pgf_pid<-reg$father_pid[match(d$f_pid, reg$pid)]


#### Remove kids who are missing from census ####
d <- d[which(d$pid %in% census_ta_data$pid),]

# Calculate average relatedness to members of the community and number of aunts and uncles in the community

census_ta_data$fid <- census_ta_data$origin_FamID

for(i in 1:nrow(d)){
  family_id <- census_ta_data$fid[match(d$pid[i],census_ta_data$pid)]
  coresidents <- census_ta_data$pid[which(census_ta_data$origin_com==d$Comunidad[i] & census_ta_data$fid!=family_id & census_ta_data$pid %in% colnames(r_matrix))]
  if(d$m_pid[i] %in% rownames(r_matrix))  d$r_mother[i] <- mean(r_matrix[d$m_pid[i], coresidents],na.rm=T)
  if(d$f_pid[i] %in% rownames(r_matrix))  d$r_father[i] <- mean(r_matrix[d$f_pid[i], coresidents],na.rm=T)
  if(d$pid[i] %in% rownames(r_matrix)) d$r_child[i] <- mean(r_matrix[d$pid[i], coresidents],na.rm=T) 
  
  d$mat_aunts[i] <- sum(reg$male[match(coresidents,reg$pid)]==0 & (reg$mother_pid[match(coresidents,reg$pid)]==d$mgm_pid[i] | reg$father_pid[match(coresidents,reg$pid)]==d$mgf_pid[i]),na.rm=T)
  
  d$mat_uncles[i] <- sum(reg$male[match(coresidents,reg$pid)]==1 & (reg$mother_pid[match(coresidents,reg$pid)]==d$mgm_pid[i] | reg$father_pid[match(coresidents,reg$pid)]==d$mgf_pid[i]),na.rm=T)
  
  d$pat_aunts[i] <- sum(reg$male[match(coresidents,reg$pid)]==0 & (reg$mother_pid[match(coresidents,reg$pid)]==d$pgm_pid[i] | reg$father_pid[match(coresidents,reg$pid)]==d$pgf_pid[i]),na.rm=T)
  
  d$pat_uncles[i] <- sum(reg$male[match(coresidents,reg$pid)]==1 & (reg$mother_pid[match(coresidents,reg$pid)]==d$pgm_pid[i] | reg$father_pid[match(coresidents,reg$pid)]==d$pgf_pid[i]),na.rm=T)
  
  progress(i,nrow(d))
}

## add gps data from the census ##

d$latitude <- census_ta_data$LAT[match(d$pid, census_ta_data$pid)]
d$longitude <- census_ta_data$LONG[match(d$pid, census_ta_data$pid)]

## Create database of grandparents
gp_pids <- unique(c(d$mgm_pid,d$mgf_pid,d$pgm_pid,d$pgf_pid))
gp <- data.frame(pid=gp_pids, community = census_ta_data$origin_com[match(gp_pids, census_ta_data$pid)])
gp$in_census <- as.numeric(gp$pid %in% census_ta_data$pid)
gp$in_vis <- as.numeric(gp$pid %in% vis$pid)

## extract grandparent gps from visit register, taking visit closest to time allocation visit
gp$community[is.na(gp$community)] <- sapply(gp$pid[is.na(gp$community)], function(x) return(vis$community[match(x,vis$pid)][which.min(as.Date(vis$date[match(x,vis$pid)]))]))
gp$community <- unlist(lapply(gp$community, function(x)ifelse(identical(x,character(0)), NA,x)))
gp$community[which(gp$community=="")] <- NA

for(i in 1:nrow(gp)){
  if(gp$pid[i] %in% census_ta_data$pid){
    gp$latitude[i] <-  census_ta_data$LAT[match(gp$pid[i], census_ta_data$pid)]
    gp$longitude[i] <-  census_ta_data$LONG[match(gp$pid[i], census_ta_data$pid)]
  }else {
    gp$latitude[i] <-  com$lat[match(gp$community[i], com$comunidad)]
    gp$longitude[i] <-  com$long[match(gp$community[i], com$comunidad)] 
  }
}

## create new columns for number of aunts, uncles, and grandparents ##

d$uncles <- d$pat_uncles+d$mat_uncles
d$aunts <- d$pat_aunts + d$mat_aunts

d$mat_kin <- d$mat_uncles + d$mat_aunts
d$pat_kin <- d$pat_uncles + d$pat_aunts

d$pgm_in_com <- d$Comunidad==gp$community[match(d$pgm_pid,gp$pid)]
d$pgf_in_com <- d$Comunidad==gp$community[match(d$pgf_pid,gp$pid)]
d$pgm_in_com[is.na(d$pgm_in_com)] <- d$pgf_in_com[is.na(d$pgm_in_com)] 
d$pgf_in_com[is.na(d$pgf_in_com)] <- d$pgm_in_com[is.na(d$pgf_in_com)] 

d$mgf_in_com <- d$Comunidad==gp$community[match(d$mgf_pid,gp$pid)]
d$mgm_in_com <- d$Comunidad==gp$community[match(d$mgm_pid,gp$pid)]
d$mgm_in_com[is.na(d$mgm_in_com)] <- d$mgf_in_com[is.na(d$mgm_in_com)] 
d$mgf_in_com[is.na(d$mgf_in_com)] <- d$mgm_in_com[is.na(d$mgf_in_com)] 

## create new columns indicating whether children live with at least some maternal or paternal kin in their community ##
d$any_mg_in_com_sure <- as.numeric(d$mgf_in_com | d$mgm_in_com)
d$any_mg_in_com <- d$any_mg_in_com_sure
d$any_mg_in_com[is.na(d$any_mg_in_com)] <- as.numeric(d$mat_aunts[is.na(d$any_mg_in_com)] >0 | d$mat_uncles[is.na(d$any_mg_in_com)] >0)

d$any_pg_in_com_sure <- as.numeric(d$pgf_in_com | d$pgm_in_com)
d$any_pg_in_com <- d$any_pg_in_com_sure
d$any_pg_in_com[is.na(d$any_pg_in_com)] <- as.numeric(d$pat_aunts[is.na(d$any_pg_in_com)] >0 | d$pat_uncles[is.na(d$any_pg_in_com)] >0)


## code the post-marital residence as a function of whether children live in the same community as their maternal grandparents, paternal grandparents, both, or neither ##
d$Residence_Pattern_sure <- NA
d$Residence_Pattern_sure[which(d$any_pg_in_com_sure & !d$any_mg_in_com_sure)] <- "Patrilocal"
d$Residence_Pattern_sure[which(d$any_mg_in_com_sure & !d$any_pg_in_com_sure)] <- "Matrilocal"
d$Residence_Pattern_sure[which(d$any_mg_in_com_sure & d$any_pg_in_com_sure)] <- "Bilocal"
d$Residence_Pattern_sure[which(!d$any_mg_in_com_sure & !d$any_pg_in_com_sure)] <- "Neolocal"
d$Residence_Pattern_sure <- relevel(as.factor(d$Residence_Pattern_sure),ref="Neolocal")

## Code the post-marital residence as a function of whether they live in the same community as their grandparents OR aunts or uncles (see methods section)##
d$Residence_Pattern <- NA
d$Residence_Pattern[which(d$any_pg_in_com & !d$any_mg_in_com)] <- "Patrilocal"
d$Residence_Pattern[which(d$any_mg_in_com & !d$any_pg_in_com)] <- "Matrilocal"
d$Residence_Pattern[which(d$any_mg_in_com & d$any_pg_in_com)] <- "Bilocal"
d$Residence_Pattern[which(!d$any_mg_in_com & !d$any_pg_in_com)] <- "Neolocal"
d$Residence_Pattern <- relevel(as.factor(d$Residence_Pattern),ref="Bilocal")

d$any_mat_aunts <- as.numeric(d$mat_aunts>0)
d$any_pat_aunts <- as.numeric(d$pat_aunts>0)
d$any_mat_uncles <- as.numeric(d$mat_uncles>0)
d$any_pat_uncles <- as.numeric(d$pat_uncles>0)
d$any_aunts <- as.numeric(d$aunts>0)
d$any_uncles <- as.numeric(d$uncles>0)

## calculate maternal age  ##

d$age <- as.numeric(d$age)
d$m_age<- as.numeric(as.Date(d$fecha)-as.Date(reg$date_of_birth[match(d$m_pid,reg$pid)]))/365.25


## now calculate distance to grandparents using gps data ##
#### Now also do distance to maternal and paternal grandparents

for (i in 1:nrow(d)){
  d$dist_to_pgm[i] <- distm(d[i,c("longitude","latitude")],as.numeric(gp[match(d$pgm_pid[i],gp$pid),c("longitude","latitude")]))
  d$dist_to_mgm[i] <- distm(d[i,c("longitude","latitude")],as.numeric(gp[match(d$mgm_pid[i],gp$pid),c("longitude","latitude")]))
  d$dist_to_pgf[i] <- distm(d[i,c("longitude","latitude")],as.numeric(gp[match(d$pgf_pid[i],gp$pid),c("longitude","latitude")]))
  d$dist_to_mgf[i] <- distm(d[i,c("longitude","latitude")],as.numeric(gp[match(d$mgf_pid[i],gp$pid),c("longitude","latitude")]))
  
  if((is.na(d$dist_to_pgm[i])&  is.na(d$dist_to_pgf[i])) | (is.na(d$dist_to_mgm[i])&  is.na(d$dist_to_mgf[i])) ){
    d$closer_to_maternal_fam[i] <-  NA
  } else{
    d$closer_to_maternal_fam[i] <- as.numeric(min(d$dist_to_mgm[i],d$dist_to_mgf[i],na.rm=T)<min(d$dist_to_pgm[i],d$dist_to_pgf[i],na.rm=T))
    d$equidistant[i] <- as.numeric(min(d$dist_to_mgm[i],d$dist_to_mgf[i],na.rm=T)==min(d$dist_to_pgm[i],d$dist_to_pgf[i],na.rm=T))
  }  
  progress(i,nrow(d))  
}

## code which observations were during the morning timeblocks  ##

d$Morning <- as.numeric(as.numeric(substr(as.character(d$hora),1,2))<=12)
d$Morning[is.na(d$Morning)] <- 1 # this is just for convenience as some mornings were miscoded

## code distance to grandparents  ##


d$dist_to_mat_gramps <- apply(d[,c("dist_to_mgm","dist_to_mgf")], 1, 
                                function(x) ifelse(all(is.na(x)),NA,min(x, na.rm=T)))

d$dist_to_pat_gramps <- apply(d[,c("dist_to_pgm","dist_to_pgf")], 1, 
                                function(x) ifelse(all(is.na(x)),NA,min(x, na.rm=T)))

### Now create the data frame for the mothers, dmom ###
### MOTHER DATABASE ###

dmom <- ta[which(ta$pid %in% d$m_pid & !is.na(ta$home_cluster)),]
dmom$nkids <- sapply(dmom$pid,function(x)sum(d$m_pid[!duplicated(d$pid)]==x,na.rm=T))
dmom$nkids_under7 <- sapply(dmom$pid,function(x)sum(d7$m_pid[!duplicated(d7$pid)]==x,na.rm=T))
dmom$youngest_child <- unlist(sapply(dmom$pid,function(x)d$pid[which(d$m_pid==x & !is.na(d$f_pid))][which.min(d$age[which(d$m_pid==x & !is.na(d$f_pid))])]))
dmom$Residence_Pattern <- d$Residence_Pattern[match(dmom$youngest_child,d$pid)]
dmom$Residence_Pattern_sure <- d$Residence_Pattern_sure[match(dmom$youngest_child,d$pid)]

dmom$dist_to_fam <- d$dist_to_mat_gramps[match(dmom$youngest_child,d$pid)]

dmom$dist_to_inlaws <- d$dist_to_pat_gramps[match(dmom$youngest_child,d$pid)]

dmom$community_r <- d$r_mother[match(dmom$youngest_child,d$pid)]
dmom$community_affinal_r <- d$r_father[match(dmom$youngest_child,d$pid)]
dmom$active_care <- dmom$cat1=="childcare" | dmom$cat2=="childcare"

dmom[,c("sgroupsize","help_with_labor","help_with_food","help_with_manufacturing","resource_acquisition_helpers","agroupsize","labor_partners","food_partners","manufacturing_partners","resource_partners")] <-0

helpers <- c() # reg in same social group (talking or within 3 meters)
partners <- c() # reg in same activity group (engaged in same activity)

## loop through all scans and record who is in a social or activity group with the focal woman, and calculate their relationship to her using the family function ##

for(i in 1:nrow(dmom)){
  
  husband <- unique(d14$f_pid[which(d14$pid == dmom$youngest_child[i])])
  dmom$husband[i] <- husband
  dmom$husband_rel[i] <- family(dmom$pid[i],husband,paste=T)
  if(length(husband)>1) print(paste0(i,"/ WARNING: ",dmom$pid[i]," has more han 1 baby daddy"))
  kids <- reg$pid[which(reg$mother_pid == dmom$pid[i])]
  
  # social group first (here, helpers)
  {dmom_helpers <- ta[which(ta$sgroupID==dmom$sgroupID[i] & ta$age>=14 & !ta$pid %in% c(dmom$pid[i],husband,kids)),]
    
    dmom_helpers$relation <- unlist(sapply(dmom_helpers$pid,family,pid2=dmom$pid[i],paste=T))
    dmom_helpers_relation <- gsub("mother-in-law|father-in-law|spouse's sibling's spouse|spouse's sibling|sibling's spouse","",dmom_helpers$relation)
    
    dmom_helpers$husbandrelation <- unlist(sapply(dmom_helpers$pid,family,pid2=husband[1],paste=T))
    dmom_helpers_husbandrelation <- gsub("mother-in-law|father-in-law|spouse's sibling's spouse|spouse's sibling|sibling's spouse","",dmom_helpers$husbandrelation)
    
    dmom_helpers_husbandrelation[dmom_helpers_husbandrelation %in% c("",";")] <- dmom_helpers_relation[dmom_helpers_relation %in% c("",";")] <- NA
    
    dmom_helpers$moms_fam <- unlist(sapply(dmom_helpers_relation,function(x)as.numeric(!grepl(c("law|None|married|son|daughter|child|offspring"),x) & !is.na(x))))
    dmom_helpers$dads_fam <- unlist(sapply(dmom_helpers_husbandrelation,function(x)as.numeric(!grepl(c("law|None|married|son|daughter|offspring|child"),x) & !is.na(x))))
    
    
    if(!grepl("cousin",dmom$husband_rel[i]) & any(dmom_helpers$moms_fam==1 & dmom_helpers$dads_fam==1 & !dmom_helpers$relation %in% c("nephew","niece"))){
      print(paste0(i,"/ Warning: ",dmom$pid[i], "shares a relative with her husband, ",dmom_helpers$pid[which(dmom_helpers$moms_fam==1 & dmom_helpers$dads_fam==1)]))
    }
    dmom$m_rels[i] <- sum(dmom_helpers$moms_fam)
    dmom$f_rels[i] <- sum(dmom_helpers$dads_fam)
    
    dmom$sgroupsize[i] <- sum(ta$sgroupID==dmom$sgroupID[i])
    dmom$out_sgroupsize[i] <- nrow(dmom_helpers)
    dmom$adult_out_sgroupsize[i] <- length(which(dmom_helpers$age>=18))
    
    if(dmom$pid[i] %in% colnames(r_matrix)){
      dmom$avg_r [i] <- mean(r_matrix[dmom$pid[i],colnames(r_matrix)[which(colnames(r_matrix) %in% ta$pid[which(ta$sgroupID==dmom$sgroupID[i])])]])
      dmom$avg_r_helpers [i] <- mean(r_matrix[dmom$pid[i],colnames(r_matrix)[which(colnames(r_matrix) %in% dmom_helpers$pid)]])
    }
    
    if(nrow(dmom_helpers)>0){
      dmom_helpers$labor_helpers <- as.numeric(grepl("garden labor|wage labor",dmom_helpers$desc1) | grepl("garden labor|wage labor",dmom_helpers$desc2))
      dmom$help_with_labor[i] <- sum(dmom_helpers$labor_helpers)
      
      dmom_helpers$resource_helpers <- 0
      dmom_helpers$resource_helpers[which(dmom_helpers$cat1=="resource acquisition" | dmom_helpers$cat2=="resource acquisition")] <- 1
      dmom$resource_acquisition_helpers[i] <- sum(dmom_helpers$resource_helpers)
      
      dmom_helpers$manu_helpers <-0
      dmom_helpers$manu_helpers[which(dmom_helpers$cat1=="manufacture"| dmom_helpers$desc1=="bring jatata"| dmom_helpers$cat2=="manufacture" | dmom_helpers$desc2=="bring jatata")] <- 1
      dmom$help_with_manufacturing[i] <- sum(dmom_helpers$manu_helpers)
      
      dmom_helpers$food_helpers <- 0
      dmom_helpers$food_helpers[which(dmom_helpers$cat1=="food processing" | dmom_helpers$cat2=="food processing")] <- 1
      dmom$help_with_food[i] <- sum(dmom_helpers$food_helpers)
      
      
      dmom_helpers$ego_id <- dmom$pid[i]
      helpers <- rbind(helpers,dmom_helpers)
    }
  }
  
  # Now activity group (here, partners)
  if(!is.na(dmom$agroup[i])){dmom_partners <- ta[which(ta$agroupID==dmom$agroupID[i] & ta$age>=14 & !ta$pid %in% c(dmom$pid[i],husband,kids)),]
  
  dmom_partners$relation <- unlist(sapply(dmom_partners$pid,family,pid2=dmom$pid[i],paste=T))
  dmom_partners_relation <- gsub("mother-in-law|father-in-law|spouse's sibling's spouse|spouse's sibling|sibling's spouse","",dmom_partners$relation)
  
  dmom_partners$husbandrelation <- unlist(sapply(dmom_partners$pid,family,pid2=husband[1],paste=T))
  dmom_partners_husbandrelation <- gsub("mother-in-law|father-in-law|spouse's sibling's spouse|spouse's sibling|sibling's spouse","",dmom_partners$husbandrelation)
  
  dmom_partners_husbandrelation[dmom_partners_husbandrelation %in% c("",";")] <- dmom_partners_relation[dmom_partners_relation %in% c("",";")] <- NA
  
  dmom_partners$moms_fam <- unlist(sapply(dmom_partners_relation,function(x)as.numeric(!grepl(c("law|None|married|son|daughter|child|offspring"),x) & !is.na(x))))
  dmom_partners$dads_fam <- unlist(sapply(dmom_partners_husbandrelation,function(x)as.numeric(!grepl(c("law|None|married|son|daughter|offspring|child"),x) & !is.na(x))))
  
  
  if(!grepl("cousin",dmom$husband_rel[i]) & any(dmom_partners$moms_fam==1 & dmom_partners$dads_fam==1 & !dmom_partners$relation %in% c("nephew","niece"))){
    print(paste0(i,"/ Warning: ",dmom$pid[i], "shares a relative with her husband, ",dmom_partners$pid[which(dmom_partners$moms_fam==1 & dmom_partners$dads_fam==1)]))
  }
  dmom$m_rels[i] <- sum(dmom_partners$moms_fam)
  dmom$f_rels[i] <- sum(dmom_partners$dads_fam)
  
  dmom$agroupsize[i] <- sum(ta$agroupID==dmom$agroupID[i])
  dmom$out_agroupsize[i] <- nrow(dmom_partners)
  dmom$adult_out_agroupsize[i] <- length(which(dmom_partners$age>=18))
  
  if(dmom$pid[i] %in% colnames(r_matrix)){
    dmom$avg_r [i] <- mean(r_matrix[dmom$pid[i],colnames(r_matrix)[which(colnames(r_matrix) %in% ta$pid[which(ta$agroupID==dmom$agroupID[i])])]])
    dmom$avg_r_partners [i] <- mean(r_matrix[dmom$pid[i],colnames(r_matrix)[which(colnames(r_matrix) %in% dmom_partners$pid)]])
  }
  
  if(nrow(dmom_partners)>0){
    dmom_partners$labor_partners <- as.numeric(grepl("garden labor|wage labor",dmom_partners$desc1) | grepl("garden labor|wage labor",dmom_partners$desc2))
    dmom$labor_partners[i] <- sum(dmom_partners$labor_partners)
    
    dmom_partners$resource_partners <- 0
    dmom_partners$resource_partners[which(dmom_partners$cat1=="resource acquisition" | dmom_partners$cat2=="resource acquisition")] <- 1
    dmom$resource_partners[i] <- sum(dmom_partners$resource_partners)
    
    dmom_partners$manu_partners <-0
    dmom_partners$manu_partners[which(dmom_partners$cat1=="manufacture"| dmom_partners$desc1=="bring jatata"| dmom_partners$cat2=="manufacture" | dmom_partners$desc2=="bring jatata")] <- 1
    dmom$manufacturing_partners[i] <- sum(dmom_partners$manu_partners)
    
    dmom_partners$food_partners <- 0
    dmom_partners$food_partners[which(dmom_partners$cat1=="food processing" | dmom_partners$cat2=="food processing")] <- 1
    dmom$food_partners[i] <- sum(dmom_partners$food_partners)
    
    
    dmom_partners$ego_id <- dmom$pid[i]
    partners <- rbind(partners,dmom_partners)
  }
  }
  
  progress(i,nrow(dmom))
}


dmom$Residence_Pattern <- relevel(as.factor(dmom$Residence_Pattern),ref="Patrilocal")
dmom$closer_to_own_family <- d$closer_to_maternal_fam[match(dmom$pid,d$m_pid)]
dmom$equidistant <- d$equidistant[match(dmom$pid,d$m_pid)]
dmom$age_centered <- (dmom$age-min(dmom$age,na.rm=T))/10
dmom$age_centered2 <-dmom$age_centered^2

dmom$Morning <- as.numeric(as.numeric(substr(as.character(dmom$hora),1,2))<=12)
dmom$Morning[is.na(dmom$Morning)] <- 1


## Now subset the children dataframe to those under 7 for the analyses ##
d7 <- d[which(d$age<=7),]
d7$age_centered <- d7$age - mean(d7$age,na.rm=T)
d7$age_centered2 <- d7$age_centered^2
d7$all_direct_nonsib_allocare <- as.numeric(d7$non_sibling_allocare | d7$social_nonsib_allocare)
d7$all_nonsib_activity_allocare <- as.numeric(d7$non_sibling_allocare | d7$social_nonsib_allocare | d7$passive_ac_agroup)
d7$all_nonsib_allocare <- as.numeric(d7$non_sibling_allocare | d7$social_nonsib_allocare | d7$passive_ac_agroup | d7$passive_ac_sgroup)
d7$Residence_Pattern <- relevel(as.factor(d7$Residence_Pattern),ref="Patrilocal")
d7_dist <- d7[which(d7$equidistant==0 & !is.na(d7$closer_to_maternal_fam)),]


## Switch IDs and communities ##
dmom$ID <- paste0("mom",match(dmom$pid,unique(dmom$pid)))
dmom$communityID <-paste0("Community",match(dmom$Comunidad,unique(dmom$Comunidad)))

dmom_columns_to_retain <- c(
  'Residence_Pattern',
  'communityID',
  'ID',
  'Morning',
  'out_sgroupsize',
  'out_agroupsize',
  'age_centered',
  'age_centered2',
  'community_r',
  'community_affinal_r',
  'active_care',
  'help_with_food',
  'help_with_labor',
  'help_with_manufacturing',
  'resource_acquisition_helpers',
  'food_partners',
  'labor_partners',
  'resource_partners',
  'manufacturing_partners',
  'dist_to_fam',
  'dist_to_inlaws'
  
)
dmom_public <- dmom[,dmom_columns_to_retain]

## Now do the same for the children's database ##
d7$ID <- paste0("child",match(d7$pid,unique(d7$pid)))
d7$communityID <-paste0("Community",match(d7$Comunidad,unique(d7$Comunidad)))
d7$motherID<-paste0("Mother",match(d7$m_pid,unique(d7$m_pid)))

d7_columns_to_retain <- c(
  'Residence_Pattern',
  'communityID',
  'ID',
  'motherID',
  'Morning',
  'non_sibling_allocare',
  'all_direct_nonsib_allocare',
  'age_centered',
  'age_centered2',
  'unsupervized',
  'dist_to_pat_gramps',
  'dist_to_mat_gramps'
)
d7_public <- d7[,d7_columns_to_retain]


## write.csv(d7_public, "data_for_public/final_analysis_table_children.csv",row.names=F)
##

## write.csv(dmom_public, "data_for_public/final_analysis_table_mothers.csv",row.names=F)

###### PREP FOR SOME ADDITIONAL GRAPHS  #######

{
  helpers$which <- "Related to neither"
  helpers$which[which(helpers$moms_fam==1)] <- "Related to wife"
  helpers$which[which(helpers$dads_fam==1)] <- "Related to husband"
  helpers$which[which(helpers$dads_fam==1 & helpers$moms_fam==1)] <-"Related to both" 
  
  helpers_full <- helpers
  helpers <- helpers[!duplicated(paste(helpers$ego_id, helpers$pid)),]
  
  
  helpers_plot <- cbind(prop.table(table(helpers$which[helpers$labor_helpers==1])),
                        prop.table(table(helpers$which[helpers$food_helpers==1])),
                        prop.table(table(helpers$which[helpers$resource_helpers==1])),
                        prop.table(table(helpers$which[helpers$manu_helpers==1]))
  )[c("Related to wife","Related to both","Related to husband","Related to neither"),]
  
  colnames(helpers_plot) <- c("Labor", "Food preparation", "Hunting/Fishing/Gathering","Manufacturing")
  
}
  ## write.csv(helpers_plot,'social_group_plot.csv')
  
  
  ta_acts <-ta[which(ta$cat1 %in% c('manufacture','resource acquisition','food processing','labor')),]
  sexbias <- prop.table(table(ta_acts$male,ta_acts$cat1),margin=2)[,c(2,1,4,3)]
  sexbias[2,] <- sqrt(sexbias[1,]*(1-sexbias[1,])/table(ta_acts$cat1))
  colnames(sexbias)<-colnames(helpers_plot)
  rownames(sexbias)<-c('Proportion Women','Standard Error')
 
  ## write.csv(sexbias,'activity_sexbias.csv')
  
  
  partners$which <- "Related to neither"
  partners$which[which(partners$moms_fam==1)] <- "Related to wife"
  partners$which[which(partners$dads_fam==1)] <- "Related to husband"
  partners$which[which(partners$dads_fam==1 & partners$moms_fam==1)] <-"Related to both" 
  
  partners_full <- partners
  partners <- partners[!duplicated(paste(partners$ego_id, partners$pid)),]
  partners_plot <- cbind(prop.table(table(partners$which[partners$labor_partners==1])),
                         prop.table(table(partners$which[partners$food_partners==1])),
                         prop.table(table(partners$which[partners$resource_partners==1])),
                         prop.table(table(partners$which[partners$manu_partners==1]))
  )[c("Related to wife","Related to both","Related to husband","Related to neither"),]
  colnames(partners_plot) <- c("Labor", "Food preparation", "Hunting/Fishing/Gathering","Manufacturing")
  
  ## write.csv(partners_plot,'activity_group_plot.csv')