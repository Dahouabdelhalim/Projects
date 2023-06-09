#############################################################
dataset_label <- c("all", "ff","fm")

for( dataset in 1:length(dataset_label) ){


# import data
dat <- read.csv( paste0("undirected_forage_matrix_", dataset_label[dataset],".csv" ), stringsAsFactors = F )

# get stan code
source("srm_undirected_model.R")
 
#############################################################
#### prepare data for Stan ##################################
#############################################################

## data list 
dat_list <- list(
                N = nrow(dat),
                N_i__j_kt_ID = max(dat$ikt_ID, dat$jkt_ID),
                N_ij_kt_ID = max(dat$ijkt_ID),
                y = dat$count,
                denominator = dat$exposure,
                ikt_ID = dat$ikt_ID,
                jkt_ID = dat$jkt_ID,
                ijkt_ID = dat$ijkt_ID)

n_chains = 4
n_cores = 4
n_warmup = 100
n_iter = 200
n_thin = 1

# model
m0 <- stan( model_code=model_code_0 , data = dat_list , chains=n_chains , cores=n_cores, warmup = n_warmup, iter = n_iter, thin = n_thin, control=list(adapt_delta=0.99))

# save  
save.image( file=paste0("m0_",dataset_label[dataset],"_forage.RDATA" ) )

}