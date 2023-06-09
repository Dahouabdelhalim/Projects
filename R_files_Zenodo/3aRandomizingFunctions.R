#Randomzing Functions

#Directly from Strona et al. 2014 (Nature Communications)

curve_ball<-function(m){
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  dimnames(rm)<-dimnames(m) #added retention of dimnames
  rm
}


library(SocGen) #devtools::install_github("vjf2/SocGen")

make_random<-function(sf){
  rand1 <- list()
  
  for (i in unique(sf$twomonth_interval)) {
    thesemonths <- sf[sf$twomonth_interval == i, ]
    
    infected <- unique(thesemonths$dolphin_id[which(thesemonths$infection_status == 1)])
    
    #make matrix
    m <- dat2mat(data.frame(thesemonths$observation_id, thesemonths$dolphin_id, 1))
    m[is.na(m)] <- 0
    
    m1 <- curve_ball(m)
    
    #convert back to original format
    d1 <- mat2dat(m1)
    d1$infection_status <- ifelse(d1$ID2 %in% infected, 1, 0)
    
    #add infection status
    rand1[[i]] <- d1[d1$values == 1, c(1,2,4)]
  }
  
  randsurv <- data.table::rbindlist(rand1)
  
  randsurv$date <- sf$date[match(randsurv$ID1, sf$observation_id)]
  
  names(randsurv)[1:2]<-c("observation_id", "dolphin_id")
  
  return(as.data.frame(randsurv))
  
}
