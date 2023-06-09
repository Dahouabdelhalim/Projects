###### CORE GROUPS MODEL ###### 
# state variables:
# g = mean group size
# G = number of groups

# parameters:
# Bf = number of breeding females
# b = birth 
# bf = number of offspring/breeding female
# d = death
# f = fission
# c = fusion
# a = aggression/attack (intraspecific mortality)
# e = group extinction rate


GrpDynCore <- function(t, x, params){
  
  ## create vectors to hold states through time
  g <- x[1]
  G <- x[2]
  
  ## model
  with(as.list(params),{
    
    ## equation 6 in text
    dg <- b*g - d*g - (f*g*g/(gf+g))*(d/(c+d)) - (e*G)*g
    ## equation 7 in text
    dG <- (c*f*g*g*G)/((gf+g)*(d+c)) - (a+e*G)*G
    
    ## output
    res <- c(dg, dG)
    list(res)
  })
}


GrpDyn <- function(t, x, params){
  
  ## create vectors to hold states through time
  g <- x[1]
  G <- x[2]
  
  ## model
  with(as.list(params),{
    
    ## equation 6 in text + breeding females specified
    dg <- bf*(g/(Bf+g)) - d*g - (f*g*g/(gf+g))*(d/(c+d)) - (a+e*G)*g
    ## equation 7 in text
    dG <- (c*f*g*g*G)/((gf+g)*(d+c)) - a*G - e*G*G
    
    ## output
    res <- c(dg, dG)
    list(res)
  })
}



###### CORE GROUPS MODEL + ALLEE ###### 
# additional parameter:
# A = allee threshold

GrpDynAllee <- function(t, x, params){
  
  ## create vectors to hold states through time
  g <- x[1]
  G <- x[2]
  
  ## model
  with(as.list(params),{
    
    ## equation 6 in text + breeding females specified + Allee term
    dg <- bf*(g/(Bf+g))*(g/(g+A)) - d*g - (f*g*g/(gf+g))*(d/(c+d)) - (a+e*G)*g
    ## equation 7 in text
    dG <- (c*f*g*g*G)/((gf+g)*(d+c)) - a*G - e*G*G
    
    ## output
    res <- c(dg, dG)
    list(res)
  })
}

