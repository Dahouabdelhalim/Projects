
#################################################################################
###   EXTRACT INDIVIDUAL STIMATES: MEAN AND SD   ################################
#################################################################################
dataset_label <- c("all","ff", "fm")

for(i in 1:length(dataset_label)){

  d <-read.csv(paste0("survival_dataset_",dataset_label[i],".csv"), stringsAsFactors=FALSE)
 

  # load model results  
  load( paste0("m0_",dataset_label[i],"_forage.RDATA"))

  # model summary
  m0p <- precis( m0, depth = 3 )


  # get a list of individual idx/group/year combination indices
  idx <- sort( unique( c( dat$ikt_ID, dat$jkt_ID) ) )
  r <- data.frame( idx = idx)


  # posterior mean
  post_mean <- sapply( r$idx, function(z)
      m0p[ which(m0p@row.names == paste0("v_i__j_kt_ID[", z,"]")  ) ,1] ) 

  # sd
  post_sd <- sapply( r$idx, function(z)
    m0p[ which(m0p@row.names == paste0("v_i__j_kt_ID[", z,"]")  ) ,2 ] )

  # get n_eff and Rhat4
  n_eff <- sapply( r$idx, function(z) 
    m0p[ which( m0p@row.names == paste0( "v_i__j_kt_ID[", z,"]")  ) , 5]  ) 
  Rhat4 <- sapply( r$idx, function(z) 
    m0p[ which( m0p@row.names == paste0( "v_i__j_kt_ID[", z,"]")  ) , 6]  ) 


  # make results dataset for foraging
  r$post_mean = post_mean
  r$post_sd = post_sd
  r$n_eff = round(n_eff,0)
  r$Rhat4 = round(Rhat4,2)


  # add to survival dataset
  d$f_mean <- unlist( sapply( 1:length(d$female), function(z) r$post_mean[ which( r$idx==d$idx_forage[z] ) ] ) ) 
  d$f_sd <- unlist( sapply( 1:length(d$female), function(z) r$post_sd[ which( r$idx==d$idx_forage[z] ) ]  ) )

write.csv(d, paste0( "survival_dataset_",dataset_label[i],".csv"), row.names = F)
}

