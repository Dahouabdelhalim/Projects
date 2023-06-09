
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)


# Input the bromeliad matrix filtered to the bromeliads that you want, the full abundance and full trait 
  #data frame
#replaced eval(parse(text with as.name( on May 2018 as former stopped working (version issue?)


extract_matrices_env <- function(Bromeliad, Abundance, Traits, env_var){
  
  ## Assertions
  if(!is.data.frame(Bromeliad)){stop('Bromeliad is not a data frame')}
  if(!is.data.frame(Abundance)){stop('Abundance is not a data frame')}
  if(!is.data.frame(Traits)){stop('Traits is not a data frame')}
  if(any(c("visit_id", "bromeliad_id","actual_water", "total_detritus", "open.canopy") %in% colnames(Bromeliad))==FALSE)
  {
    stop("Bromeliad doesn't contain the right environmental columns");
  }
  
  
  # Create E matrix by taking the appropriate environmental variables. Remove
  # any sites where the collection method is NA, and also remove visit 81
  E<-eval(substitute(select(Bromeliad, visit_id, bromeliad_id,env_var), list(env_var = as.name(env_var))))%>% 
    filter(visit_id !=101 & visit_id !=81) %>% 
    dplyr::select(-visit_id)
  
  if(env_var=="open.canopy")
  {
    E$open.canopy<-as.numeric(E$open.canopy)
  }
  
  if(env_var=="elevation_m")
  {
    E$elevation_m<-as.numeric(E$elevation_m)
  }
  
  # Choose from E those bromeliads which are kept, and keep only abundances from
  # those bromeliads
  bromeliads_kept <- data.frame(bromeliad_id = E$bromeliad_id)
  W_all <- inner_join(bromeliads_kept, Abundance) %>% as.data.frame()
  
  
  E <- semi_join(E, W_all, by = "bromeliad_id" )
  
  # bromeliad IDs to rownames
  rownames(W_all) <- W_all$bromeliad_id
  
  W_all <- W_all %>% 
    select(-bromeliad_id)
  species_kept <- W_all %>% colSums(na.rm =TRUE) > 0
  W <- W_all[,species_kept]
  
  B <- colnames(W) %>% 
    str_split('_') %>% 
    purrr::map_df(~ data.frame(ds_id = .x[1], species_id = .x[2])) %>% 
    mutate(species_id = as.character(species_id)) %>% #added
    mutate(species_id = as.double(species_id)) %>% #added
    mutate(ds_id = as.character(ds_id)) %>% #added
    mutate(ds_id = as.double(ds_id)) %>% #added
    left_join(Traits, by = c("species_id" = "species_id"))%>% 
    filter(realm == 'aquatic' & micro_macro == 'macro') %>% 
    tidyr::unite("species_id", c(ds_id, species_id), sep = "_") %>%
    select(species_id, matches("^[A-Z]{2}\\\\d"))
  
  species_filtered <- B %>% select(species_id) %>% unlist()
  W <- W[species_filtered]
  
  expmaker<-function(y){y^0.5}
  
  W<-apply(W, 2, expmaker)
  
  rownames(E) <- E$bromeliad_id
  rownames(B) <- B$species_id 
  
  E <- E %>% select(-bromeliad_id) 
  B <- B %>% select(-species_id)
  
  BWE <- list('B' = B, "W" =  W, "E" = E)
  return(BWE)
}


