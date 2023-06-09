##### Run the two functions
##### FUNCTION 1
# The following "CHVind" function was built by Gonzales et al. 2017, Frontiers
# Function compute convex hull (verties + volume) of two set of points and their intersection in a D-dimensional space (D>=2)
# Inputs:
#  - set1: matrix(N1*D) with coordinates of the first set of points in D-dimensions
#  - set2: matrix(N2*D) with coordinates of the second set of points in D-dimensions
#  ->  N1 and N2 must be strictly higher than D
#  - NA are not allowed
# outputs:                                                                                                                   
#  - $vert_set1: matrix with coordinates of the vertices delimiting the convex hull of set1                          
#  - $vert_set2: matrix with coordinates of the vertices delimiting the convex hull of set2                          
#  - $vert_intersect: matrix with coordinates of the vertices delimiting the intersection between the two sets                       
#  - $vol: vector (length 3) with convex hull volume of set1, set2 and their intersection                            
# This code was cut and pasted from the original version

CHVind<-function(set1) {
  # checking inputs 
  # If the number of rows is less than or equal to the number of columns = STOP ncol =# + 1 points
  if (nrow(set1)<=ncol(set1))  stop(paste("error : 'set1' must contain at least ",ncol(set1)+1, " points",sep=""))
  # checking absence of NA
  if (length(which(is.na(set1)==TRUE))!=0) stop(paste("error : NA are not allowed in coordinates of 'set1'",sep=""))
  # tranforming coordinates in true rational number written as character string
  set1rep <- d2q(cbind(0, cbind(1, set1)))  
  # reduce set of points to vertices only using redundant function
  polytope1 <- redundant(set1rep, representation = "V")$output  
  # changing polytope representation: vertices to inequality constraints
  H_chset1 <- scdd(polytope1, representation = "V")$output    
  # extracting coordinates of vertices
  vert_set1 <- q2d(polytope1[ , - c(1, 2)])
  # computing convex hull volume of the two polytopes and of their intersection if it exists (no intersection if one vertex in common)
  vol<-rep(0,1) ; names(vol)<-c("vol_set1")
  vol["vol_set1"]<-convhulln(vert_set1,"FA")$vol  
  
  # results
  res<-list(vert_set1=vert_set1, vol=vol)
  return(res)    
} # end of function CHVind

##### FUNCTION 2
# The following "CHVintersect" (intersect of the hypervolumes) function was built by Gonzales et al. 2017, Frontiers
# function to compute convex hull (vertices + volume) of two set of points and their intersection in a D-dimensions space (D>=2)                                                                                                              
# inputs:                               
#  - set1: matrix(N1*D) with coordinates of the first set of points in D-dimensions                                   
#  - set2: matrix(N2*D) with coordinates of the second set of points in D-dimensions                                 
#  ->  N1 and N2 must be strictly higher than D    
#  - NA are not allowed                          
# outputs:                                      
#  - $vert_set1: matrix with coordinates of the vertices delimiting the convex hull of set1                          
#  - $vert_set2: matrix with coordinates of the vertices delimiting the convex hull of set2                          
#  - $vert_intersect: matrix with coordinates of the vertices delimiting the intersection between the two sets                                
#  - $vol: vector (length 3) with convex hull volume of set1, set2 and their intersection      
# This code was cut and pasted from the original version

CHVintersect<-function(set1,set2) {
  # checking inputs
  # number of dimensions
  if(ncol(set1)!=ncol(set2)) stop("error : different number of dimensions in the two set of points") 
  # number of points
  if (nrow(set1)<=ncol(set1))  stop(paste("error : 'set1' must contain at least ",ncol(set1)+1, " points",sep=""))
  if (nrow(set2)<=ncol(set2))  stop(paste("error : 'set2' must contain at least ",ncol(set2)+1, " points",sep=""))
  # checking absence of NA
  if (length(which(is.na(set1)==TRUE))!=0) stop(paste("error : NA are not allowed in coordinates of 'set1'",sep=""))
  if (length(which(is.na(set2)==TRUE))!=0) stop(paste("error : NA are not allowed in coordinates of 'set2'",sep=""))
  
  # tranforming coordinates in true rational number written as character string
  set1rep <- d2q(cbind(0, cbind(1, set1)))
  set2rep <- d2q(cbind(0, cbind(1, set2)))
  # reduce set of points to vertices only using redundant function
  polytope1 <- redundant(set1rep, representation = "V")$output
  polytope2 <- redundant(set2rep, representation = "V")$output
  # changing polytope representation: vertices to inequality constraints
  H_chset1 <- scdd(polytope1, representation = "V")$output
  H_chset2 <- scdd(polytope2, representation = "V")$output
  # intersection between the two polytopes
  H_inter <- rbind(H_chset1, H_chset2)
  V_inter <- scdd(H_inter, representation = "H")$output
  # extracting coordinates of vertices
  vert_set1 <- q2d(polytope1[ , - c(1, 2)])
  vert_set2 <- q2d(polytope2[ , - c(1, 2)])
  vert_inter <- q2d(V_inter[ , - c(1, 2)])
  
  # computing convex hull volume of the two polytopes and of their intersection if it exists (no intersection if one vertex in common)
  vol<-rep(0,3) ; names(vol)<-c("vol_set1","vol_set2","vol_inter")
  vol["vol_set1"]<-convhulln(vert_set1,"FA")$vol
  vol["vol_set2"]<-convhulln(vert_set2,"FA")$vol
  if (is.matrix(vert_inter)==TRUE) # vector if one vertex in common
    if( nrow(vert_inter)>ncol(vert_inter) ) vol["vol_inter"]<-convhulln(vert_inter,"FA")$vol
  
  # results
  res<-list(vert_set1=vert_set1, vert_set2=vert_set2, vert_inter=vert_inter, vol=vol)
  return(res)    
} # end of function CHVintersect

################################################################################################
# Function to calculate standard errors
se=function(x) sd(x)/sqrt(length(x))

