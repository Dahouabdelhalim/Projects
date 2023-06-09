# scenario simulations - first run
#out <- array(dim=c(50,2,100), data=NA)

#system.time(for(i in 1:100){
#out[,,i] <- popsim(nyears=50, init_pop=c(100,100), dphi_sa = .9986, pr_b= .4230, nb_dphi = .9992, a_in_night_prob = c(0.333,0.333,0.334), 
#       a_rmort = 0.0407, a_in_num_guards = c(0,0,0), guard_limit = 10, a_p_dphi = 0.9994, p_days = 14, fec = 165, lmort = 0.0523, 
#       a_out_night_prob = c(0.333,0.333,0.334), rmort2 = 0.0407, a_out_nguard = c(0,0,0), met_night_prob = c(0.333,0.333,0.334), 
#       rmort_met = 0.0407, ng_met = c(0,0,0))
#})  

# load data
source('./R/popsim.r')
source('./R/bprob.r')
source('./R/migrate.r')
source('./R/fecund.r')
source('./R/surv.r')
# fecund, migrate
load('./data/scenarios.Rdata')

# create object for simulation output
my_fun <- function(){
  list(d=array(dim=c(50,2,500), data=NA))
}
replicate(361, my_fun(), simplify=TRUE)->out

# list option
# loop through first 100 scenarios, going to do 500 iterations
system.time(for(i in 1:361){
  
   # pull in road xing guard data from scenarios object
  gds_ai <- c(scenarios$night1_adult_in_migtation[i],scenarios$night2_adult_inout_migtation[i],scenarios$night3_adult_inout_migtation[i])
  gds_ao <- c(scenarios$night2_adult_inout_migtation[i],scenarios$night3_adult_inout_migtation[i],scenarios$night4_adult_out_migtation[i])
  gd_m <- c(scenarios$night5_meta_out_migration[i], scenarios$night6_meta_out_migration[i], scenarios$night7_meta_out_migration[i])
  # loop through for 500 iterations
    for(j in 1:500){
    out[[i]][,,j] <- popsim(nyears=50, init_pop=c(100,100), dphi_sa = .9986, pr_b= .4230, nb_dphi = .9992, a_in_night_prob = c(0.333,0.333,0.334), 
            a_rmort = 0.0407, a_in_num_guards = gds_ai, guard_limit = 10, a_p_dphi = 0.9994, p_days = 14, fec = 165, lmort = 0.0523, 
            a_out_night_prob = c(0.333,0.333,0.334), rmort2 = 0.0407, a_out_nguard = gds_ao, met_night_prob = c(0.333,0.333,0.334), 
            rmort_met = 0.0407, ng_met = gd_m)
    
      } #j
    if(i%%1==0) {save(out, file='./data/outdata.Rdata')}
  }) #i

