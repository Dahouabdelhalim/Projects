
library(dplyr)
library(purrr)
library(tidyverse)
library(stringr)
require(devtools)
library(vegan)
library(SYNCSA)

source('01_extracting_matrices_env.R')

'%nin%' <- function(x,y)!('%in%'(x,y))

#just to get function to load

Broms.all <-read_csv("bromeliads_syncsa_ready.csv", 
                         col_types = cols(
                           visit_id = col_character(),
                           bromeliad_id = col_character(),
                           dataset_name = col_character(),
                           total_detritus = col_number(),
                           open.canopy = col_number(),
                           actual_water = col_number(),
                           sites = col_character(),
                           total_detritus_center= col_number(),
                           logtotal_detritus= col_number(),
                           logtotal_detritus_center= col_number(),
                           include.detritus= col_number(),
                           open.canopy_center= col_number(),
                           include.canopy= col_number(),
                           actual_water_center = col_number(),
                           logactual_water = col_number(),
                           logactual_water_center = col_number(),
                           include.water= col_number()
                         ) 
) %>% 
  as.data.frame()

####Define function For each region, we do the following analysis

##this has an put.together traits error for fuzzy trait that appeared after I inserted rao, entirely due to pattern = rao (repalcing with pattern = tdap removes issue)

single.site.syncsa<-function(sitename, scale, trait.type, env_var,  trait.name = "all", ro.method = "procrustes", N = "all", stock = "abundance")
  {
  
  
  if((trait.name!="all"))
  {
    B_axis2<-B_axis2%>%select(species_id,trait.name)
  }else{
    B_axis2<-B_axis2
  }
  
  if((stock=="abundance"))
  {
    Stock<-Abun
  }else{
    Stock<-Biom
  }


  Broms2<- Broms.all %>%
    filter(sites%in%sitename & !is.na(paste(env_var)))
  
  Broms<-Broms2[which(!is.na(Broms2[env_var])),]
 

  if(is.numeric(N))
    {
    Broms.temp<-Broms[sample(1:nrow(Broms), size = N, replace = FALSE),]
    }else{
    Broms.temp<-Broms
    }

  my_mat <- extract_matrices_env(Bromeliad = Broms.temp, Abundance = Stock, Traits = Tra, env_var)

  my_mat$B %>% dim()
  my_mat$W  %>% dim()
  my_mat$E %>% dim()

  if (scale=="log")
    {
    my_mat$E[,"logenv"]<-log(my_mat$E[,1])
    my_mat$E<-my_mat$E[-1]
    }

#Finally time to run Syncsa on fuzzy traits

if (trait.type=="fuzzy")
  {
    DATA<-organize.syncsa(comm=my_mat$W,traits=my_mat$B,envir=my_mat$E, convert.traits=TRUE, ranks=TRUE)

    traits_kept<-DATA$traits%>%colSums(na.rm =TRUE)>0
    DATA$traits<-DATA$traits[,traits_kept]
    trait_list <- split(colnames(DATA$traits), colnames(DATA$traits) %>% substring(first=1,last=2))

    RES_procrustes<-syncsa(DATA$community,DATA$traits,envir=DATA$environmental, ro.method=ro.method, ranks=TRUE, put.together = trait_list, na.rm=TRUE, permutations = 100, notification=TRUE)

    best.tcap.traits<-optimal(DATA$community, envir = DATA$environmental, DATA$traits, subset.min = 1, subset.max = 3, pattern = "tcap", ro.method= ro.method, dist = "euclidean",  method = "pearson",ranks=TRUE,na.rm=TRUE, put.together = trait_list, notification = TRUE, progressbar = TRUE)
    best.tdap.traits<-optimal(DATA$community, envir = DATA$environmental, DATA$traits, subset.min = 1, subset.max = 3, pattern = "tdap", ro.method= ro.method, dist = "euclidean",  method = "pearson",ranks=TRUE,na.rm=TRUE, put.together = trait_list, notification = TRUE, progressbar = TRUE)
    best.rao.traits<-optimal(DATA$community, envir = DATA$environmental, DATA$traits, subset.min = 1, subset.max = 3, pattern = "rao", ro.method= ro.method, dist = "euclidean",  method = "pearson",ranks=TRUE,na.rm=TRUE, put.together = trait_list, notification = TRUE, progressbar = TRUE)

    # tcap optimized
        traits_tcap<-as.list(unlist(strsplit(best.tcap.traits$optimization$Subset[1], "[[:blank:]]+"))) 
        my_mat$B2<-my_mat$B%>%select_(.dots = traits_tcap)  
            
        DATA<-organize.syncsa(comm=my_mat$W,traits=my_mat$B2,envir=my_mat$E, convert.traits=TRUE, ranks=TRUE)
    
        traits_kept<-DATA$traits%>%colSums(na.rm =TRUE)>0
        DATA$traits<-DATA$traits[,traits_kept]
        trait_list2 <- split(colnames(DATA$traits), colnames(DATA$traits) %>% substring(first=1,last=2))
        
        RES_procrustes_optim_tcap<-syncsa(DATA$community,DATA$traits,envir=DATA$environmental, ro.method=ro.method, ranks=TRUE, put.together = trait_list2, na.rm=TRUE, permutations = 100, notification=TRUE)

    #tdap optimized    
        traits_tdap<-as.list(unlist(strsplit(best.tdap.traits$optimization$Subset[1], "[[:blank:]]+")))  
        my_mat$B3<-my_mat$B%>%select_(.dots = traits_tdap)
        
        DATA<-organize.syncsa(comm=my_mat$W,traits=my_mat$B3,envir=my_mat$E, convert.traits=TRUE, ranks=TRUE)
        
        traits_kept<-DATA$traits%>%colSums(na.rm =TRUE)>0
        DATA$traits<-DATA$traits[,traits_kept]
        trait_list3<- split(colnames(DATA$traits), colnames(DATA$traits) %>% substring(first=1,last=2))
        
        RES_procrustes_optim_tdap<-syncsa(DATA$community,DATA$traits,envir=DATA$environmental, ro.method=ro.method, ranks=TRUE, put.together = trait_list3, na.rm=TRUE, permutations = 100, notification=TRUE)

    #rao optimized    
        traits_rao<-as.list(unlist(strsplit(best.rao.traits$optimization$Subset[1], "[[:blank:]]+")))  
        my_mat$B4<-my_mat$B%>%select_(.dots = traits_rao)
        
        DATA<-organize.syncsa(comm=my_mat$W,traits=my_mat$B4,envir=my_mat$E, convert.traits=TRUE, ranks=TRUE)
        
        traits_kept<-DATA$traits%>%colSums(na.rm =TRUE)>0
        DATA$traits<-DATA$traits[,traits_kept]
        trait_list4<- split(colnames(DATA$traits), colnames(DATA$traits) %>% substring(first=1,last=2))
        
        RES_procrustes_optim_rao<-syncsa(DATA$community,DATA$traits,envir=DATA$environmental, ro.method=ro.method, ranks=TRUE, put.together = trait_list4, na.rm=TRUE, permutations = 100, notification=TRUE)
        
  }
    else
  {
 
    B_axis_repeated <- rownames(my_mat$B) %>% str_split('_') %>% 
    map_df(~ data.frame(ds_id = .x[1], species_id = .x[2])) %>% 
    left_join(B_axis2, by = 'species_id') %>% 
    unite("species_id", c(ds_id, species_id), sep = "_")%>%
    select(-species_id)

    rownames(B_axis_repeated)<-rownames(my_mat$B)

    DATA<-organize.syncsa(comm=my_mat$W,traits=B_axis_repeated,envir=my_mat$E,  convert.traits=FALSE, ranks=FALSE)#note that ranks are now false as continuous

    RES_procrustes<-syncsa(DATA$community,DATA$traits,envir=DATA$environmental,  ro.method=ro.method, ranks=FALSE,  na.rm=TRUE, permutations = 100, notification=TRUE)
    
    if((trait.name=="all"))
    {
    best.tcap.traits<-optimal(DATA$community, envir = DATA$environmental, DATA$traits, subset.min = 1, subset.max = 3, pattern = "tcap",ro.method= ro.method, dist = "euclidean", ranks=FALSE, method = "pearson",na.rm=TRUE,  notification = TRUE, progressbar = TRUE)
    best.tdap.traits<-optimal(DATA$community, envir = DATA$environmental, DATA$traits, subset.min = 1, subset.max = 3, pattern = "tdap",ro.method= ro.method, dist = "euclidean", ranks=FALSE, method = "pearson",na.rm=TRUE,  notification = TRUE, progressbar = TRUE)
    best.rao.traits<-optimal(DATA$community, envir = DATA$environmental, DATA$traits, subset.min = 1, subset.max = 3, pattern = "rao",ro.method= ro.method, dist = "euclidean", ranks=FALSE, method = "pearson",na.rm=TRUE,  notification = TRUE, progressbar = TRUE)

    # tcap optimized
        traits_tcap<-as.list(unlist(strsplit(best.tcap.traits$optimization$Subset[1], "[[:blank:]]+")))  
        B_axis_optimum_tcap<-B_axis_repeated%>%select_(.dots = traits_tcap)
        
        DATA1<-organize.syncsa(comm=my_mat$W,traits=B_axis_optimum_tcap,envir=my_mat$E,  convert.traits=FALSE, ranks=FALSE)#note that ranks are now false as continuous
        
        RES_procrustes_optim_tcap<-syncsa(DATA1$community,DATA1$traits,envir=DATA1$environmental,  ro.method=ro.method, ranks=FALSE,  na.rm=TRUE, permutations = 100, notification=TRUE)
        
    # tdap optimized
        traits_tdap<-as.list(unlist(strsplit(best.tdap.traits$optimization$Subset[1], "[[:blank:]]+")))  
        B_axis_optimum_tdap<-B_axis_repeated%>%select_(.dots = traits_tdap)
        DATA2<-organize.syncsa(comm=my_mat$W,traits=B_axis_optimum_tdap,envir=my_mat$E,  convert.traits=FALSE, ranks=FALSE)#note that ranks are now false as continuous
        RES_procrustes_optim_tdap<-syncsa(DATA2$community,DATA2$traits,envir=DATA2$environmental,  ro.method=ro.method, ranks=FALSE,  na.rm=TRUE, permutations = 100, notification=TRUE)

    # rao optimized
        traits_rao<-as.list(unlist(strsplit(best.rao.traits$optimization$Subset[1], "[[:blank:]]+")))  
        B_axis_optimum_rao<-B_axis_repeated%>%select_(.dots = traits_rao)
        DATA3<-organize.syncsa(comm=my_mat$W,traits=B_axis_optimum_rao,envir=my_mat$E,  convert.traits=FALSE, ranks=FALSE)#note that ranks are now false as continuous
        RES_procrustes_optim_rao<-syncsa(DATA3$community,DATA3$traits,envir=DATA3$environmental,  ro.method=ro.method, ranks=FALSE,  na.rm=TRUE, permutations = 100, notification=TRUE)
    } else
      
    {

    }
        }

cv.env<-as.numeric(apply(DATA$environmental, 2, var))/as.numeric(apply(DATA$environmental, 2, mean))
min.env<-as.numeric(apply(DATA$environmental, 2, min))
max.env<-as.numeric(apply(DATA$environmental, 2, max))

if((trait.name=="all"))
  {
return(cbind(RES_procrustes$statistics[1,1],RES_procrustes$statistics[2,1], RES_procrustes$statistics[6,1], RES_procrustes$statistics[10,1], 
             RES_procrustes$statistics[1,2],RES_procrustes$statistics[2,2], RES_procrustes$statistics[6,2], RES_procrustes$statistics[10,2],
              best.tcap.traits$optimization$Subset[1], best.tcap.traits$optimization$ro[1], 
              best.tdap.traits$optimization$Subset[1], best.tdap.traits$optimization$ro[1], 
              best.rao.traits$optimization$Subset[1], best.rao.traits$optimization$ro[1],
              cv.env, min.env, max.env, 
              RES_procrustes_optim_tcap$statistics[1,1], RES_procrustes_optim_tcap$statistics[1,2],
              RES_procrustes_optim_tdap$statistics[6,1], RES_procrustes_optim_tdap$statistics[6,2], 
              RES_procrustes_optim_rao$statistics[10,1], RES_procrustes_optim_rao$statistics[10,2]) )
  }

else
  
  {
return(cbind(RES_procrustes$statistics[1,1],RES_procrustes$statistics[2,1], RES_procrustes$statistics[6,1], RES_procrustes$statistics[10,1], 
             RES_procrustes$statistics[1,2],RES_procrustes$statistics[2,2], RES_procrustes$statistics[6,2], RES_procrustes$statistics[10,2]))
    }

}



##===and now just to make the matrices

single.site.matrixmaker<-function(sitename, scale, trait.type, env_var, N = "all", stock="abundance")
{
  
  Broms<- Broms.all%>%filter(sites%in%sitename)
  
  if((stock=="abundance"))
  {
    Stock<-Abun
  }else{
    Stock<-Biom
  }
  
  if(is.numeric(N)){
    Broms.temp<-Broms[sample(1:nrow(Broms), size = N, replace = FALSE),]
  }else{
    Broms.temp<-Broms
  }
  
  my_mat <- extract_matrices_env(Bromeliad = Broms.temp, Abundance = Stock, Traits = Tra, env_var)
  
  my_mat$B %>% dim()
  my_mat$W  %>% dim()
  my_mat$E %>% dim()
  
  if (scale=="log")
  {
    my_mat$E[,"logenv"]<-log(my_mat$E[,1])
    my_mat$E<-my_mat$E[-1]
  }
  
  B_axis_repeated <- rownames(my_mat$B) %>% str_split('_') %>% 
    map_df(~ data.frame(ds_id = .x[1], species_id = .x[2])) %>% 
    left_join(B_axis2, by = 'species_id') %>% 
    unite("species_id", c(ds_id, species_id), sep = "_")%>%
    select(-species_id) #ds
  
  rownames(B_axis_repeated)<-rownames(my_mat$B)
  
  my_mat$B<-as.matrix(my_mat$B)
  my_mat$E<-as.matrix(my_mat$E)
  my_mat$W<-as.matrix(my_mat$W)
  
  j<-matrix(NA,nrow=dim(my_mat$W )[1],ncol=dim(my_mat$W )[2])
  for(i in 1:dim(my_mat$W )[2]){
    j[,i]<-t(as.matrix(rowSums(my_mat$W )))
  }
  
  my_mat$propW<-my_mat$W /j
  
  my_mat$CwMfuzzy<-my_mat$propW%*%my_mat$B
  my_mat$CwMaxis<-my_mat$propW%*%as.matrix(B_axis_repeated)
  
  my_mat$ECwMfuzzy<-as.data.frame(cbind(my_mat$E, my_mat$CwMfuzzy))
  my_mat$ECwMaxis<-as.data.frame(cbind(my_mat$E, my_mat$CwMaxis))
  
  BWElist<-list(
    B = my_mat$B,
    W = my_mat$W,
    E = my_mat$E,
    CwMfuzzy =my_mat$CwM,
    CwMaxis =my_mat$CwM,
    ECwMfuzzy = my_mat$ECwMfuzzy,
    ECwMaxis = my_mat$ECwMaxis,
    B_axis_repeated = B_axis_repeated)
  
  return(BWElist)
  
}



