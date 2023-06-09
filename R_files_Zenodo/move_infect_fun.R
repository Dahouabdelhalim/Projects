#This script creates functions to model infection dynamics & population sizes
#in a population of moving animals
#Claire Teitelbaum
#claire.teitelbaum@gmail.com


#FUNCTION FOR SIMULATION
move_infect = function(pars,pops,A,
                       depart_fun = depart_informed,
                       n_times = 30){
  library(data.table)
  library(dplyr)
  
  #inputs:
  #pars: parameters needed in functions below. Required parameters:
      #beta_w: contact rate at natural patches
      #beta_diff: difference in contact rate between urban and natural patches (0=identical)
      #gamma: recovery rate
      #v: maximum infection-induced mortality. see text.
      #mu_0: density-independent mortality rate
      #mu_1: density-dependent mortality rate
  #pops: initial populatinon sizes. Must include N_0, p, q, infect_status as columns
  #A: matrix of resource dynamics. rows are timesteps, columns are sites.
  #n_times: how many timesteps of output of diff eq to include? 2 will be just beginning and end (and is minimum value)
  
  #RUN SIMS
  
  #create vector of urban/natural site types in order of sites in pops
  site_types = list(pops[!duplicated(pops$site),"p"])
  
  #initiate matrices to save outputs
  moves = pops[,-which(colnames(pops)=="N_0")] #movement in each time step
  site_dynamics = pops[,-which(colnames(pops)=="N_0")] #at-site compartment changes (mortality and infection)
  
  for(t in 1:pars$tmax){
    pops$N = pops[,paste("N" , t-1 , sep = "_")]
    
    #extract resource availablity at timestep
    pops$A = A[t,pops$site]
    As = A[t,] ; names(As) = paste0("A_",1:ncol(A))
    
    #recalculate total individuals at site across all infection classes and types
    #N_i = pops %>% group_by(site) %>% summarize(N = sum(N)) %>% pull(N) #old (equivalent) way using dplyr is slower
    N_i = as.data.table(pops)
    N_i = N_i[, .(
      N = sum(N)
    ), by = site]
    N_i = N_i[,N]
    
    pops$N_i = N_i[pops$site]
    
    #departure
    M = depart_fun(A = pops$A, N = pops$N, N_i = pops$N_i)
    
    
    pops$M = M
    pops$N = pops$N - pops$M
    #sum(pops$N) + sum(pops$M)  #keep track of total population size, make sure it hasn't changed
    
    #recalculate total individuals at site across all infection classes and types
    N_i = as.data.table(pops)
    N_i = N_i[, .(
      N = sum(N)
    ), by = site]
    N_i = N_i[,N]
    
    pops$N_i = N_i[pops$site]
    
    #destinations
    M_qij = destination_fun(M = pops$M , A_all = A[t,] , pops = pops)
    
    #aggregate destinations per destination site and impose mortality cost of movement
    dt = as.data.table(M_qij)
    N_in = dt[, .(
      site = j,
      N_in = sum(N)
    ), by = .(infect_status,q,j)] %>% as.data.frame()
    
    pops = left_join(pops , N_in , by = c("site", "q", "infect_status")) %>% mutate(N = N + N_in) %>% dplyr::select(-N_in , -j)
    #sum(pops$N)
    
    #total individuals at site across all infection classes and types
    N_i = as.data.table(pops)
    N_i = N_i[, .(
      N = sum(N)
    ), by = site]
    N_i = N_i[,N]
    
    pops$N_i = N_i[pops$site]
    
    
    #infection and demography occurs at each site
    initial_pops = pops$N
    pops = mutate(pops , full_name = paste(infect_status ,q ,site , sep = "_"))
    names(initial_pops) = pops$full_name
    #run infection function through desolve
    out = deSolve::ode(y = initial_pops, 
                       times = seq(0,pars$w,length.out=n_times), 
                       func = infect_dem, 
                       parms = c(pars,As,site_types=site_types),
                       method = "lsode", mf = 10)
    
    N_new = out[nrow(out),-1]
    N_new = N_new[names(initial_pops)]
    id_result = initial_pops - N_new
    
    
    pops$N = N_new
    #sum(pops$N) #keep track of total population size
    
    
    #save the end population sizes for this timestep
    pops[,paste("N",t,sep="_")] = pops$N
    #save number moving out (for movement rate calculation)
    moves[,paste("M",t,sep="_")] = pops$M
    #save infection and demography output
    site_dynamics[,paste("ID",t,sep="_")] = id_result
    
  }
  
  #convert to matrix for saving (remove non-numeric info)
  info = select(pops , site, p , q , infect_status) %>% as.matrix() #character matrix giving info for each row in other matrices
  pops$N_final = pops[,ncol(pops)]
  pops = select(pops , -N , -A , -M , -N_i , -full_name, -p,-q,-infect_status,-site) %>% 
    as.matrix()
  
  moves = select_at(moves, vars(contains("M_"))) %>% as.matrix()
  site_dynamics = select_at(site_dynamics , vars(contains("ID_"))) %>% as.matrix()
  
  return(list(pops = pops , moves = moves , site_dynamics = site_dynamics, info = info))
}

#FUNCTIONS
#infection and mortality function in continuous time
#s-i-s
infect_dem = function(time , x , parms){
  
  class = stringr::str_split_fixed(names(x) , "_" , 3)[,1]
  types = stringr::str_split_fixed(names(x) , "_" , 3)[,2]
  sites = stringr::str_split_fixed(names(x) , "_" , 3)[,3]
  
  out = NULL
  for(i in unique(sites)){ #go through each site
    
    #generalist = r , urban = u , natural = w
    S_ri = x[class == "S" & types == "r" & sites == i]
    I_ri = x[class == "I" & types == "r" & sites == i]
    N_ri = sum(x[types == "r" & sites == i])
    
    S_wi = x[class == "S" & types == "w" & sites == i]
    I_wi = x[class == "I" & types == "w" & sites == i]
    N_wi = sum(x[types == "w" & sites == i])
    
    S_ui = x[class == "S" & types == "u" & sites == i]
    I_ui = x[class == "I" & types == "u" & sites == i]
    N_ui = sum(x[types == "u" & sites == i])
    
    I_i = sum(x[class == "I" & sites == i])
    N_i = sum(x[sites == i])
    A_i = parms[[(paste0("A_",i))]]
    
    if("beta_u" %in% names(parms)){ #assign beta value depending on habitat type. makes it possible to have beta differ across site types.
      beta = ifelse(parms$site_types[as.numeric(i)]=="u" , parms$beta_u , parms$beta_w)
    }
    
    out <- with(as.list(parms),{
      mu_1 = ifelse(N_i < A_i , 0 , mu_1)
      dSr = (-1)*beta*S_ri*I_i + gamma*I_ri - S_ri*(mu_0 + mu_1*(N_i/A_i)) 
      names(dSr) = paste("S_r",i,sep="_")
      dSw = (-1)*beta*S_wi*I_i + gamma*I_wi - S_wi*(mu_0 + mu_1*(N_i/A_i)) 
      names(dSw) = paste("S_w",i,sep="_")
      dSu = (-1)*beta*S_ui*I_i + gamma*I_ui - S_ui*(mu_0 + mu_1*(N_i/A_i)) 
      names(dSu) = paste("S_u",i,sep="_")
      
      dIr = beta*S_ri*I_i - gamma*I_ri - v*I_ri*(1-A_i) - I_ri*(mu_0 + mu_1*(N_i/A_i))
      names(dIr) = paste("I_r",i,sep="_")
      dIw = beta*S_wi*I_i - gamma*I_wi - v*I_wi*(1-A_i) - I_wi*(mu_0 + mu_1*(N_i/A_i))
      names(dIw) = paste("I_w",i,sep="_")
      dIu = beta*S_ui*I_i - gamma*I_ui - v*I_ui*(1-A_i) - I_ui*(mu_0 + mu_1*(N_i/A_i))
      names(dIu) = paste("I_u",i,sep="_")
      
      c(out,dSr,dSw,dSu,dIr,dIw,dIu)
    })
    
  }
  out = out[names(x)] #reorder output so that it is in the same order as initial conditions
  names(out) = paste0("d",names(out)) #change names
  list(out)
}



#movement (competition rule)
depart_informed = function(A , N , N_i){
  M = ifelse(N_i > A ,  #move if N>carrying capacity
             N * (1 - (A / N_i)),
             0
  )
  return(M)
}


#destination function: distribute equally to all sites, scaled by phi
destination_fun = function(M , A_all ,pops = NULL , p = NULL , q = NULL, i = NULL , class = NULL){
  
  #M: populations departing from all sites (M_qi)
  #f: preference parameter
  #A_all: A values at current timestep, same order as p
  
  #either provide these variables or they will be extracted from pops if it is not provided
  #q: individual types (n=#of site-class combos, same length as M)
  #i: site numbers
  #class: infection status
  #p: site types (n=#of sites)
  
  if(!is.null(pops)){
    p = pops$p[!duplicated(pops$site)]
    q = pops$q
    i = pops$site
    class = pops$infect_status
  }
  
  
  n_sites = length(p) - 1
  n_u = sum(p == "u")
  n_w = sum(p == "w")
  
  p_mat = matrix(p , nrow = length(q) , ncol = length(p) , byrow = T)
  q_mat = matrix(rep(q,length(p)) , nrow = length(q) , ncol = length(p))
  phis = matrix(1/n_sites , nrow = length(q) , ncol = length(p)) #set up matrix: rows are from (i), columns are to (j)
  phis[q_mat == "r"] = 1/n_sites #go equally to all other patches if generalist
  phis[q_mat == "u" & p_mat == "w"] = 0 #apply specialization parameter if specialized
  phis[q_mat == "w" & p_mat == "u"] = 0 #apply specialization parameter if specialized
  different_site = !sapply(unique(pops$site) , function(x) pops$site == x) #cannot go back to the site you left from
  if(n_u==1){ #you can go back to the site you left from if you are a specialist and there in only one site of that type
    different_site[which(pops$p=="u" & pops$q=="u"),which(p=="u")] = 1
  } else if(n_w==1){
    different_site[which(pops$p=="w" & pops$q=="w"),which(p=="w")] = 1
  }
  phis = phis*different_site 
  phis = phis/rowSums(phis) #calculate proportion going to each site
  #impute zeroes if there are no urban or natural hosts; otherwise these turn into NAs
  if(n_u==0){
    phis[q_mat=="u"] = 0
  }
  if(n_w==0){
    phis[q_mat=="w"] = 0
  }
  
  #multiply by number departing to get number going to each site
  M_ijq = M*phis
  
  #transform to long form
  colnames(M_ijq) = NULL
  out = cbind.data.frame(M_ijq , i = as.numeric(i) , q = as.character(q) , infect_status = as.character(class) , stringsAsFactors = F) %>%
    reshape2::melt(id.vars = c("i","q","infect_status") , 
                   measure.vars = 1:length(p) , variable.name = "j" , value.name = "N") %>%
    mutate(j = as.numeric(j))
  
  
  return(out)
  
}
