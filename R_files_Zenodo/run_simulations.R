#This script sets up parameter combinations and runs all simulations
#Claire Teitelbaum
#claire.teitelbaum@gmail.com

library(dplyr); library(data.table); library(stringr)
library(doParallel)
#load in functions that run model
source("move_infect_fun.R")
source("resource_sim_fun.R")
options(stringsAsFactors = FALSE)



#### 1. SET PARAMETERS ####
n_sims = 20 #number of repeated simulations for each prop_urban

#"baseline" parameter values: use these when holding values constant
pars_base = cbind(P=10, urbval = 0.5, urbvar = 0.1, 
                  mu_0 = 0.001 , mu_1 = (0.001*2),
                  beta_w = 0.15, beta_diff = 0,  gamma = 0.03, v = 0.01,
                  prop_generalist = 1, 
                  w = 7, tmax = 52, prop_infected = 0.01)
pars_base = as.data.frame(pars_base)

#set up parameter combinations to test
params_list = list(
  pars_base,
  #mortality rates
  pars_mort = expand.grid(mu_0=c(0, 0.001, 0.002) , mu_1=c(0, 0.001*2, 0.002*2),
                          beta_w = 0, gamma = 0, v = 0, urbvar = c(0, 0.1)) 
  ,
  #infection parameters
  pars_infect = expand.grid(beta_w = c(0, 0.06, 0.15, 0.17), gamma = c(0, pars_base$gamma, 0.06), v = c(0, 0.001, pars_base$v, 0.1),
                            urbvar = c(0, 0.1)) 
  ,
  #specialization
  pars_spec = expand.grid(prop_generalist = seq(0,1,0.1), beta_w = 0, gamma = 0, v = 0, urbvar = c(0, 0.1))
  ,
  #specialization with infection and urbval
  pars_spec_infect = expand.grid(prop_generalist = seq(0,1,0.1), urbval = c(0.3, pars_base$urbval), urbvar = c(0, 0.1))
  ,
  #number of patches
  pars_P = expand.grid(P = c(10,50), 
                       beta_w = c(0, pars_base$beta_w), gamma = c(0, pars_base$gamma), v = c(0, pars_base$v)) 
  ,
  #urban resource level
  pars_urbval = expand.grid(urbval = c(0.3,0.5,0.7), 
                            beta_w = c(0, pars_base$beta_w), gamma = c(0, pars_base$gamma), v = c(0, pars_base$v))
  ,
  #resource variability at urban patches
  pars_urbvar = expand.grid(urbvar = c(0,0.01,0.1,0.2,0.5,0.8),
                            beta_w = c(0, pars_base$beta_w), gamma = c(0, pars_base$gamma), v = c(0, pars_base$v)) 
  ,
  #beta_diff (difference in beta between urban and natural patches)
  pars_deltab = expand.grid(beta_diff = c(0, 0.01, 0.06, 0.09), prop_generalist = c(0,0.1,0.5,1)) 
)

params_df = bind_rows(params_list) 
params_df = sapply(names(params_df), function(x)  tidyr::replace_na(params_df[,x], pars_base[[x]])) 

prop_urban = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) #proportion of urban and natural patches

prop_urban_add = expand.grid(row = 1:nrow(params_df), prop_urban = prop_urban)

params_df = cbind(params_df[prop_urban_add$row,], prop_urban = prop_urban_add$prop_urban)
params_df = data.frame(params_df)

params_df$n_sims = ifelse(params_df$prop_urban==1 & params_df$urbvar==0, 1, n_sims)
params_df$prop_infected = with(params_df, ifelse(beta_w==0, 0, prop_infected))
params_df$beta_u = with(params_df, beta_w + beta_diff)
params_df$beta_avg = with(params_df, (beta_w+beta_u)/2)

params_df = params_df[!duplicated(params_df),]

#simulate resource availability and set initial population sizes for each parameter set
pars_list = lapply(1:nrow(params_df),function(i){
  out_list = lapply(1:params_df[i,"n_sims"], function(j){
    #extract parameter values
    out = as.list(params_df[i,])
    
    #simulate resource availability
    P_u = out$P * out$prop_urban
    if(P_u>0){
      A_u = sapply(1:P_u , function(x) sim_A_unif(out$tmax, min=out$urbval-(out$urbvar/2), max=out$urbval+(out$urbvar/2)))
    } else(A_u = NULL)
    P_w = out$P - P_u
    if(P_w>0){
      A_w = sapply(1:P_w , function(x) sim_A_unif(out$tmax,min=0,max=1))
    } else(A_w = NULL)
    #combine all site resource availability
    #rows are time-steps, columns are sites
    A_all = cbind(A_w , A_u)
    
    #save resource availability to master list
    out$A_all = A_all
    out$A_num = j
    
    #simulation population starting conditions
    #back-calculate proportion urban patches
    p_urb = out$prop_urban
    #extract which patches are urban
    urb = c(rep(F, P_w) , rep(T, P_u)) #urban patches are always appended to the end of the series
    
    
    #vector of site types and number of sites of each type
    site_types = ifelse(urb , "u" , "w")
    sites_u = length(which(urb))
    sites_w = length(which(!urb))
    
    #set initial population sizes average carrying capacity (i.e. equal to mean A over time across all patches)
    #separately for urban and natural sites
    #assign population sizes by site ID
    pops = expand.grid(site = 1:length(site_types),
                       q = c("w","u","r")) %>%
      mutate(p = site_types[site],
             N_site = ifelse(p=="u" , out$urbval , 0.5),
             N_total = ifelse(q=="r" , N_site*out$prop_generalist ,
                              ifelse(p==q , N_site-N_site*out$prop_generalist, 0))) %>%
      select(-N_site)
    
    #seed infection
    #start infection at each site type depending on j value: first site is natural and last site is urban when both types are present
    infection_start_loc = ifelse(j <= out$n_sims/2, 1 , out$P)
    out$infection_start_loc = infection_start_loc
    
    p_infect = out$prop_infected
    #proportion infected at each site of each type with infection
    N_total = sum(c(rep(0.5, sites_w), rep(out$urbval, sites_u))) #total initial starting population
    pops = mutate(pops , I = ifelse(site %in% infection_start_loc , N_total*p_infect , 0),
                  S = N_total-I)
    
    #translate to long form
    pops = select(pops , -N_total) %>%
      reshape2::melt(id.vars = c("site","q","p") , measure.vars = c("I","S"),
                     value.name = "N_0" , variable.name = "infect_status")
    
    out$pops = pops
    
    return(out)
  })
  #save output
  return(out_list)
})

pars_list = unlist(pars_list, recursive = F)


##### 2. RUN SIMULATIONS ####
#run simulations in parallel
no_cores = ifelse(detectCores() > 30 , 30 , detectCores()-1) #how many cores to use depends on what computer you're using

simFun = function(x){
  out = move_infect(pars = x, pops = x$pops, A=x$A_all)
  gc()
  return(out)
}

out_list = vector(mode = "list" , length = length(pars_list))
#run sub-lists of simulations each to save memory (sets of 10 simulations per core)
#create data frame of indices to separate blocks of simulations
start_end = data.frame(start = seq(1,length(pars_list), no_cores*10) ,
                       end = c(seq(no_cores*10, length(pars_list),no_cores*10), length(pars_list)))
startTime = Sys.time()
for(i in 1:nrow(start_end)){
  start = start_end[i,"start"] ; end = start_end[i,"end"]
  out_list[start:end] <- mclapply(pars_list[start:end] , simFun , mc.cores = no_cores)
  gc()
  cat(round(i/nrow(start_end) , 2)*100) ; cat("%: ") ; cat(Sys.time()-startTime); cat("\\n")
}


#option to run not in parallel
# out_list = vector(mode = "list" , length = length(pars_list))
# startTime = Sys.time()
# for(i in 1:length(pars_list)){
#   out_list[[i]] = simFun(pars_list[[i]])
#   cat(round(i/length(pars_list) , 3)*100) ; cat("%: ") ; cat(difftime(Sys.time(),startTime,units="hours")); cat(" hours \\n")
# }



saveRDS(list(out_list = out_list, pars = pars_list) , "sim_outputs/landscape_sim_output.Rds")


pars_df = do.call(rbind , pars_list)
pars_df = pars_df[,!colnames(pars_df) %in% c("A_all","A_num","pops")] %>% 
  as.data.frame() %>% mutate_all(unlist)
pars_df$sim_num = ifelse(!duplicated(pars_df), 1, 2)

#### 3. SUMMARIZE SIMULATION RESULTS ####
#calculate summary prevalence metrics for each sim
#this is the function to calculate them
calc_prev_metrics = function(out, time_final=52){
  #inputs: 
  #out: data from simulation
  #time: time step from which to calculate metrics
  
  #reorganize data
  pops = cbind.data.frame(out[["info"]] , out[["pops"]])
  pops = mutate(pops , p = factor(p , levels = c("u","w"))) #set factor levels so that summaries will be calculated even when no urban or natural patches are present
  N_0 = sum(pops$N_0)
  S_0 = sum(pops$N_0[pops$infect_status=="S"])
  
  pops = select(pops , 1:paste0("N_",time_final)) %>% rename("N_final"=paste0("N_",time_final))
  
  prev = sum(pops[pops$infect_status == "I","N_final"])/sum(pops[,"N_final"])
  #difference in prevalence in hosts (natural-urban)
  hostprevs = pops %>% group_by(q) %>% 
    summarize(prev = sum(N_final[infect_status=="I"])/n())
  host_diff = hostprevs$prev[hostprevs$q=="w"] - hostprevs$prev[hostprevs$q=="u"] 
  #generalist prevalence
  gen_prev = hostprevs$prev[hostprevs$q=="r"]
  
  #difference in prevalence at sites (natural-urban)
  siteprevs = pops %>% group_by(p , .drop = F) %>% 
    summarize(prev = sum(N_final[infect_status=="I"])/n()) 
  site_diff = siteprevs$prev[siteprevs$p=="w"] - siteprevs$prev[siteprevs$p=="u"]
  
  #number of sites at which prevalence is >0 (with some tolerance)
  prop_sites_positive = sum(siteprevs$prev>0.01)/length(siteprevs$prev)
  
  #varaiation in prevalence across sites
  end_prev_bysite = pops %>% group_by(site) %>% 
    summarize(prev = sum(N_final[infect_status=="I"])/n()) 
  spat_var = sd(end_prev_bysite$prev)
  
  #maximum infection prevalence
  max_prev = pops %>% select(N_1:N_final) %>% 
    apply(2 , function(x) sum(x[pops$infect_status=="I"])/sum(x)) %>% max()
  
  #time of maximum infection prevalence
  max_prev_time = pops %>% select(N_1:N_final) %>% 
    apply(2 , function(x) sum(x[pops$infect_status=="I"])/sum(x)) %>% which.max()
  
  #linear interpolation of spread rate until peak
  prev_init = sum(pops %>% filter(infect_status=="I") %>% pull(N_0))/sum(pops %>% pull(N_0))
  spread_rate = (max_prev-prev_init)/max_prev_time
  
  #translate to long form
  pops_long = select(pops , -N_final) %>%
    reshape2::melt(id.vars = c("site","p","q","infect_status"),
                   variable.name = "time",
                   value.name = "N") %>%
    mutate(time = str_split_fixed(time,"_",2)[,2] %>% as.numeric,
           N = as.numeric(N))
  #time it takes until all sites reach prevalence of 0.1 (starting prevalence of index site)
  spread_threshold = pops_long %>% group_by(site,time) %>%
    summarize(prev = sum(N[infect_status=="I"])/sum(N)) %>% 
    group_by(time) %>%
    summarize(metric=all(prev>0.1)) %>%
    filter(metric) %>%
    slice(1) %>% pull(time)
  spread_threshold = ifelse(length(spread_threshold)==0,NA,spread_threshold)
  
  #proportion of time prevalence at urban sites is greater than at wild sites
  if(any(pops$p=="u") & any(pops$p=="w")){
    u_greater = pops_long %>% group_by(time,p) %>%
      summarize(prev = sum(N[infect_status=="I"])/sum(N)) %>% 
      group_by(time) %>%
      summarize(uw = as.numeric(prev[p=="u"]>prev[p=="w"]), .groups = "drop") %>%
      pull(uw) %>% mean()
  } else{
    u_greater = NA
  }
  
  
  #overall movement rate across all sites and time (mean number moving per timestep)
  move_rate = out[["moves"]] %>% as.data.frame() %>% select_at(vars(contains("M_"))) %>% colSums() %>%
    mean()
  
  #overall movement rate across all sites and 1st 8 timesteps
  move_rate_init = out[["moves"]] %>% as.data.frame() %>% select(M_1:M_8) %>% colSums() %>%
    mean()
  
  #proportion of individuals moving per timestep
  move_rates_all = out[["moves"]][,1:time_final] %>% as.data.frame() %>% select_at(vars(contains("M_"))) %>% as.matrix() #number moving
  n_start = pops %>% select_at(vars(contains("N_"))) %>% select(-ncol(.)) %>% as.matrix() #total number initially
  n_pre_move = n_start - out[["site_dynamics"]][,1:time_final] #total number after infection and mortality
  move_rate_prop = mean(move_rates_all/n_pre_move , na.rm = T) #proportion of those surviving that move
  
  #population size
  end_pop = sum(pops$N_final)
  
  #annual survival
  annual_surv = sum(pops$N_final)/sum(pops$N_0)
  
  #difference in survival in hosts (natural-urban)
  hostsurvs = pops %>% group_by(q) %>% 
    summarize(surv = sum(N_final)/sum(N_0))
  surv_u = hostsurvs$surv[hostsurvs$q=="u"]
  surv_w = hostsurvs$surv[hostsurvs$q=="w"]
  surv_r = hostsurvs$surv[hostsurvs$q=="r"]
  
  
  return(cbind(end_prev = prev , host_diff  = host_diff, site_diff = site_diff, prop_sites_positive , annual_surv,
               spat_var , gen_prev , end_pop , max_prev , N_0 , S_0 , max_prev_time,
               spread_rate, spread_threshold , u_greater,
               move_rate , move_rate_init, move_rate_prop,
               surv_u, surv_w , surv_r))
} 

gc()

#CALCULATE in parallel
prev_metrics = mclapply(out_list , function(x) calc_prev_metrics(x,time_final=40), mc.cores = no_cores) 
prev_metrics_20w = mclapply(out_list , function(x) calc_prev_metrics(x,time_final=20), mc.cores = no_cores) 
prev_metrics = do.call(rbind,prev_metrics) 
prev_metrics_20w = do.call(rbind,prev_metrics_20w) 

#for infection start loc, also add type of site where it started
pars_df = pars_df %>% mutate(max_w = round(P*(1-prop_urban)) , #need to round to deal with R's invisible decimals
                             infection_w = as.numeric(infection_start_loc)<=max_w,
                             infection_start_type = ifelse(infection_w , "w" , "u"),
                             infection_start_type = tidyr::replace_na(infection_start_type,"all")) %>%
  select(-max_w,-infection_w)

results = cbind(pars_df, prev_metrics) 
results_20w = cbind(pars_df, prev_metrics_20w) 

#save results summary
saveRDS(results , "sim_outputs/results_summary.Rds")
saveRDS(results_20w , "sim_outputs/results_summary_20w.Rds")

#save subsets of one representative simulation from out_list (for faster download/loading)
example_sims = out_list[pars_df$sim_num==1]
saveRDS(list(pars = pars_df[pars_df$sim_num==1,] , out = example_sims),
        "sim_outputs/example_sims.Rds")
