

#################################################################################
###   EXTRACT INDIVIDUAL STIMATES: MEAN AND SD   ################################
#################################################################################
behavior_label <- c("groom", "support")
dataset_label <- c("all","ff", "fm")


for(i in 1:length(dataset_label) ){

  d <- read.csv( paste0("survival_dataset_incomplete_",dataset_label[i],".csv"), stringsAsFactors=FALSE )

  for( x in 1:length(behavior_label) ){

    load(paste0("m0_",dataset_label[i],"_",behavior_label[x],".RDATA"))

    # model summary
    m0p <- precis( m0, depth = 3 )
    
    # get a list of individual idx/group/year combination indexes
    idx <- sort( unique( c( dat$ikt_ID, dat$jkt_ID) ) )
    r <- data.frame( idx = idx)

    # pull out model estimated posterior giver and receiver effects
    # for grooming and support

  
    # giver posterior means
    post_mean_giver <- sapply( r$idx, function(z) 
        m0p[ which( m0p@row.names == paste0( "v_i__j_kt_ID[", z,",1]")  ) , 1]  ) 
    # giver posterior sd 
    post_sd_giver <- sapply( r$idx, function(z)
      m0p[ which( m0p@row.names == paste0( "v_i__j_kt_ID[", z,",1]")  ) , 2] )  
    
    # receiver posterior means
    post_mean_receiver <- sapply( r$idx, function(z)
      m0p[ which( m0p@row.names == paste0("v_i__j_kt_ID[", z,",2]")  ) ,1] ) 
    # receiver posterior sd
    post_sd_receiver <- sapply( r$idx, function(z)
      m0p[ which( m0p@row.names == paste0("v_i__j_kt_ID[", z,",2]")  ) ,2] )
    
    # get n_eff and Rhat4
    n_eff_giver <- sapply( r$idx, function(z) 
      m0p[ which( m0p@row.names == paste0( "v_i__j_kt_ID[", z,",1]")  ) , 5]  ) 
    Rhat4_giver <- sapply( r$idx, function(z) 
      m0p[ which( m0p@row.names == paste0( "v_i__j_kt_ID[", z,",1]")  ) , 6]  ) 
    
    n_eff_receiver <- sapply( r$idx, function(z) 
      m0p[ which( m0p@row.names == paste0( "v_i__j_kt_ID[", z,",2]")  ) , 5]  ) 
    Rhat4_receiver <- sapply( r$idx, function(z) 
      m0p[ which( m0p@row.names == paste0( "v_i__j_kt_ID[", z,",2]")  ) , 6]  ) 
    
    
    
    # make results dataset
    r$post_mean_giver <- post_mean_giver
    r$post_sd_giver <- post_sd_giver
    r$n_eff_giver = round(n_eff_giver,0)
    r$Rhat4_giver = round(Rhat4_giver,2)
    r$post_mean_receiver = post_mean_receiver
    r$post_sd_receiver = post_sd_receiver
    r$n_eff_receiver = round(n_eff_receiver,0)
    r$Rhat4_receiver = round(Rhat4_receiver,2)

  if(x == 1){
  
    d$gg_mean <- unlist(sapply( 1:length(d$female), function(z) r$post_mean_giver[ which( r$idx==d$idx_groom[z] ) ] ) )
    d$gg_sd <- unlist(sapply( 1:length(d$female), function(z) r$post_sd_giver[ which( r$idx==d$idx_groom[z] ) ] ) )
  
    d$gr_mean <- unlist(sapply( 1:length(d$female), function(z) r$post_mean_receiver[ which( r$idx==d$idx_groom[z] ) ] ) )
    d$gr_sd <- unlist(sapply( 1:length(d$female), function(z) r$post_sd_receiver[ which( r$idx==d$idx_groom[z] ) ] ) )
  
  
  }else{
  
    d$sg_mean <- unlist(sapply( 1:length(d$female), function(z) r$post_mean_giver[ which(r$idx==d$idx_support[z] ) ] ) )
    d$sg_sd <- unlist(sapply( 1:length(d$female), function(z) r$post_sd_giver[ which( r$idx==d$idx_support[z] ) ] ) )
  
    d$sr_mean <- unlist(sapply( 1:length(d$female), function(z) r$post_mean_receiver[ which( r$idx==d$idx_support[z] ) ] ) )
    d$sr_sd <- unlist(sapply( 1:length(d$female), function(z) r$post_sd_receiver[ which( r$idx==d$idx_support[z] ) ] ) )
  }

}
# save results
write.csv(d, paste0("survival_dataset_",dataset_label[i],".csv") )

}




