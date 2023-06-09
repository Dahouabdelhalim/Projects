hwi<-function(sightings=sightings, group_variable=group_variable, dates=dates, IDs=IDs, symmetric=TRUE){
    z<-sightings
    dolphins<-sort(unique(z[,IDs]))
    n<-length(dolphins)
    matnames<-list(dolphins,dolphins)
    currentmat<-matrix(c(rep(NA,n^2)),nrow=n,dimnames=matnames)
    currentdflist<-split(z, z[,IDs],drop=TRUE)
    for (i in 1:nrow(currentmat)) {
      ego<-row.names(currentmat)[i]
      #Get the list element for each dolphin
      all_ego<-get(ego, currentdflist)
      for (j in i:ncol(currentmat)) {
        alter<-colnames(currentmat)[j]
        #Get the list element for each dolphin
        all_alter<-get(alter, currentdflist)
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
          HWC<-(X/(X+Ya+Yab))
          currentmat[i,j]<-HWC}
        else{currentmat[i,j]<-0}
      }
    }
    diag(currentmat)<-NA
    if(symmetric==TRUE){
    currentmat[lower.tri(currentmat)]=t(currentmat)[lower.tri(currentmat)]
    }
    return(currentmat)
}

digu<-function(x){x$Group<-paste0(x$Permutation,"-", x$Group);return(x)}

mergeMatrices<-function(lmat) {
  rands<-lapply(lmat, function(mat) mat<-na.omit(unmatrix(mat)))
  rands<-lapply(rands, function(x) x<-x[order(names(x))])
  rands<-as.data.frame(do.call("cbind",rands))
  return(rands)
  }


sample_sightings<-function(x,nums){
  s=sample(row.names(x), nums, replace=FALSE) 
  y=subset(x, row.names(x) %in% s, drop=TRUE) 
  return(y)
}

