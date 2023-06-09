#Digidolph helper functions
#Vivienne Foroughirad
#vjf2@duke.edu
#Created August 8th, 2016
#Modified December 8, 2016 


#Creates unique observation IDs for Digiroo2::fAssoctable output

digu<-function(x){x$Group<-paste0(x$Permutation,"-", x$Group);return(x)}

#Convert association matrices to linear format and merge into dataframe

mergeMatrices<-function(lmat) {
  rands<-lapply(lmat, function(mat) mat<-na.omit(unmatrix(mat)))
  rands<-lapply(rands, function(x) x<-x[order(names(x))])
  rands<-as.data.frame(do.call("cbind",rands))
  return(rands)
}

#This feature adds the launch point and extra points so a sampling area can always be calculated

add_launch<-function(x) {
  launch[1]<-x[1,1]
  x<-rbind(x, launch)
  return(x)
}

#Calculates the half-weight index with a one-day sampling period
#Filters sightings for joint availabilty for each pair
#Takes an availability dataframe with ID, entry date, and depart date

hwi_filtered<-function(sightings=sightings, group_variable=group_variable, dates=dates, IDs=IDs, symmetric=TRUE, availability=availability){
  z<-sightings
  dolphins<-sort(unique(z[,IDs]))
  n<-length(dolphins)
  matnames<-list(dolphins,dolphins)
  currentmat<-matrix(c(rep(NA,n^2)),nrow=n,dimnames=matnames)
  currentdflist<-split(z, z[,IDs],drop=TRUE)
  for (i in 1:nrow(currentmat)) {
    ego<-row.names(currentmat)[i]
    #Get the list element for each dolphin
    ego_start<-availability$entry[availability$dolphin_id==ego]
    ego_end<-availability$depart[availability$dolphin_id==ego]
    
    #all_ego<-get(ego, currentdflist)
    for (j in i:ncol(currentmat)) {
      alter<-colnames(currentmat)[j]
      #Get the list element for each dolphin
      
      alter_start<-availability$entry[availability$dolphin_id==alter]
      alter_end<-availability$depart[availability$dolphin_id==alter]
      
      all_ego<-get(ego, currentdflist)[,c(dates,group_variable)]
      all_alter<-get(alter, currentdflist)[,c(dates,group_variable)]
      
      hstart<-max(ego_start, alter_start)
      hend<-min(alter_end, ego_end)
      tp<-hend-hstart
      if(tp<1){currentmat[i,j]<-NA} else {
        
        all_ego<-subset(all_ego, all_ego[,dates]>=hstart & all_ego[,dates]<=hend)
        all_alter<-subset(all_alter, all_alter[,dates]>=hstart & all_alter[,dates]<=hend)    
        
        #Take the intersection of the partycomp_dolphins to see when the dolphins were in the same group_variable
        set<-(intersect(all_ego[,group_variable], all_alter[,group_variable]))
        sample<-subset(all_ego, all_ego[,group_variable] %in% set)
        #Numerator for HWC
        X<-length(unique(sample[,dates]))
        if (X>0){
          
          X<-length(unique(all_ego[,dates][all_ego[,group_variable] %in% set]))
          #Ego without alter
          Ya<-length(setdiff(all_ego[,dates], all_alter[,dates]))
          #Alter without ego
          Yb<-length(setdiff(all_alter[,dates], all_ego[,dates]))
          #Both seen but not together
          Yab<-length(intersect(all_ego[,dates], all_alter[,dates]))-X
          #Half weight coefficient
          HWC<-(X/(X+0.5*(Ya+Yb)+Yab))
          currentmat[i,j]<-HWC}
        else{currentmat[i,j]<-0}
      }}
  }
  diag(currentmat)<-NA
  if(symmetric==TRUE){
    currentmat[lower.tri(currentmat)]=t(currentmat)[lower.tri(currentmat)]
  }
  return(currentmat)
}

#This function returns one set of simulations for the whole study period

one_sim<-function(d=d,buff_days=buff_days, udsgdf=udsgdf, 
                  dolphin_density_per_km=dolphin_density_per_km, 
                  schedule=schedule, 
                  gprox=gprox){
  
  for (i in 1:d){
    bound<-buff_days[i,]
    areakm<-area(bound)/1000000
    nd<-round(areakm*dolphin_density_per_km)
    nd<-ifelse(nd==1, nd<-2, nd) 
    dailygrid<-udsgdf
    rgrid <- raster(dailygrid)
    rgrid_msk <- mask(rgrid,bound, inverse=FALSE)
    
    #convert back to sgdf
    
    grid_ae <- as(rgrid_msk, 'SpatialGridDataFrame')
    gridded(grid_ae) <- TRUE
    grid_ae[[1]] <- as.numeric(!is.na(grid_ae[[1]])) 
    
    ## Then, you just have to multiply each column of udsgdf by the mask 
    resu <- lapply(1:ncol(dailygrid), function(i) { 
      dailygrid[[i]] * grid_ae[[1]] }) 
    resu <- as.data.frame(resu) 
    names(resu) <- names(dailygrid@data) 
    
    ## and define it as data slot for udsgdf 
    dailygrid@data <- resu 
    
    probweights<-colSums(dailygrid@data)
    
    daily_dolphins<-sample(names(probweights)[which(schedule[i,]==TRUE)],size=nd,replace=FALSE,prob=probweights[which(schedule[i,]==TRUE)])
    
    i_dd<-which(colnames(dailygrid@data) %in% daily_dolphins)
    
    rp<-fRanXY(i_dd, dailygrid)
    
    coordinates(rp) <- ~x+y
    
    dnn_digi <- dnearneigh(rp,0,gprox,row.names=as.character(rp$ID))
    
    dayAssoc<-as.data.frame(fAssoctable(dnn_digi))
    dayAssoc$Permutation<-c(rep(as.character(buff_days$id[i]),nd))
    dayAssoc$Permutation<-as.Date(dayAssoc$Permutation, format=c("%Y-%m-%d"))
    each_days_assoc[[i]]<-dayAssoc
  }
  return(each_days_assoc)
}

